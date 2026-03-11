"""
redox_mc.py
-----------
Redox Monte Carlo module for electron transfer simulations at constant potential.

Implements the RedoxMC class for Metropolis MC sampling of electron transfer
between electrode and redox-active molecules in solution, based on the
grand-canonical ensemble at fixed electrode potential.

Reference:
  - Reed, S. K. et al. J. Chem. Phys. 126, 084704 (2007)
  - Meissner, R. et al. J. Chem. Theory Comput. 14, 6307-6319 (2018)
"""

import json
import random
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from datetime import datetime

import openmm as mm
import openmm.app as app
import openmm.unit as unit


@dataclass
class RedoxParameters:
    """
    Parameters defining redox-active molecule behavior.

    Attributes:
        name: molecular identifier (e.g., "EC", "DMC", "qbzn")
        residue_names: list of PDB residue names to identify this molecule
        charge_states: dict mapping state_int → charge (in elementary charges)
                       e.g., {0: 0.0, -1: -1.0, -2: -2.0} for sequential reduction
        electron_affinities: dict mapping (state_from, state_to) → EA (in eV)
                            positive EA means electron addition is favorable in gas phase
                            e.g., {(0, -1): 0.5, (-1, -2): -0.8}
        allowed_states: list of integer charge states accessible at this electrode
        partition_coefficient: solvation free energy difference ΔG_solv (kJ/mol)
                              relative to reference state. Used for Boltzmann weighting.
    """
    name: str
    residue_names: List[str]
    charge_states: Dict[int, float]          # state → charge (e)
    electron_affinities: Dict[Tuple[int, int], float]  # (state_from, state_to) → EA (eV)
    allowed_states: List[int]
    partition_coefficient: float = 0.0       # solvation free energy (kJ/mol)
    state_free_energy_offsets_kjmol: Optional[Dict[int, float]] = None


class RedoxMC:
    """
    Grand-canonical Monte Carlo for electron transfer at constant electrode potential.

    Energy balance for electron transfer:
        w = EA_gas(state_from → state_to)
            + E_electrode(kJ/mol)  [±1 depending on reduction/oxidation]
            + ΔPE_solv(context)     [from NonbondedForce energy change]
            + ΔE_intra              [intramolecular energy change, default 0]

    Acceptance probability:
        P_accept = min(1, exp(-w / kT))
    where kT = R*T = 8.314462e-3 * T (kJ/mol)

    At equilibrium, the Boltzmann distribution of charge states is maintained:
        ln(n_i / n_j) = -(E_i^solv - E_j^solv) / kT + (z_i - z_j) * F * E_electrode / kT

    This reproduces the electrochemical potential balance.
    """

    FARADAY_KJMOL_PER_V = 96.48533212331002  # kJ/mol per elementary charge per Volt
    KB_KJMOL_PER_K = 8.314462e-3              # Boltzmann constant (kJ/(mol·K))
    EV_TO_KJMOL = 96.4853                      # Conversion eV → kJ/mol (= F/1000)

    def __init__(
        self,
        redox_params: RedoxParameters,
        temperature_k: float,
        voltage_v: float,  # electrode potential (V vs vacuum)
        max_states: int = 5,  # max accessible charge states
    ):
        """
        Initialize redox MC sampler.

        Args:
            redox_params: RedoxParameters defining the redox-active molecule
            temperature_k: temperature in Kelvin
            voltage_v: electrode potential in Volts vs vacuum reference
            max_states: maximum number of distinct charge states to track
        """
        self.params = redox_params
        self.temperature_k = temperature_k
        self.voltage_v = voltage_v
        self.kT = self.KB_KJMOL_PER_K * temperature_k
        self.reference_state = self._select_reference_state()
        self.reference_charge = self.params.charge_states[self.reference_state]
        self.state_energy_offsets = self._build_state_energy_offsets()

        # Statistics tracking
        self.n_accepted = 0
        self.n_attempted = 0
        self.n_accepted_by_state = {s: 0 for s in redox_params.allowed_states}
        self.n_attempted_by_state = {s: 0 for s in redox_params.allowed_states}
        self.state_history = []  # track state trajectory
        self.energy_history = []

    def _select_reference_state(self) -> int:
        """Choose a default reference state, preferring the neutral state."""
        if 0 in self.params.allowed_states:
            return 0
        return min(
            self.params.allowed_states,
            key=lambda state: abs(self.params.charge_states[state]),
        )

    def _build_state_energy_offsets(self) -> Dict[int, float]:
        """
        Build intrinsic free-energy offsets for each charge state in kJ/mol.

        The reference state is assigned zero energy. For an electron-addition step
        `state_from -> state_to`, a positive electron affinity lowers the reduced
        state's free energy by EA.
        """
        if self.params.state_free_energy_offsets_kjmol:
            offsets = {
                int(state): float(offset)
                for state, offset in self.params.state_free_energy_offsets_kjmol.items()
            }
            for state in self.params.allowed_states:
                offsets.setdefault(state, 0.0)
            return offsets

        offsets = {self.reference_state: 0.0}
        pending = True
        while pending:
            pending = False
            for (state_from, state_to), ea_ev in self.params.electron_affinities.items():
                delta_kjmol = -ea_ev * self.EV_TO_KJMOL
                if state_from in offsets and state_to not in offsets:
                    offsets[state_to] = offsets[state_from] + delta_kjmol
                    pending = True
                elif state_to in offsets and state_from not in offsets:
                    offsets[state_from] = offsets[state_to] - delta_kjmol
                    pending = True

        for state in self.params.allowed_states:
            offsets.setdefault(state, 0.0)

        return offsets

    def _solvation_correction(self, state: int) -> float:
        """
        Apply the configured solvation correction relative to the reference state.
        """
        state_charge = self.params.charge_states[state]
        delta_charge = abs(state_charge - self.reference_charge)
        return delta_charge * self.params.partition_coefficient

    def compute_state_energy(
        self,
        state: int,
        pe_nonbonded: float,  # kJ/mol, from NonbondedForce
    ) -> float:
        """
        Compute total energy of molecule in given charge state.

        E_total(state) = E_electrode(state) + PE_solv(state)

        where:
        - E_electrode = -z * F * E_potential  (z = charge in e, E in V)
        - PE_solv = nonbonded potential energy

        Returns energy in kJ/mol.
        """
        if state not in self.params.allowed_states:
            return float('inf')

        z = self.params.charge_states[state]
        intrinsic_offset = self.state_energy_offsets.get(state, 0.0)

        # Electrode energy: -z * F * V.
        # More negative electrode potentials stabilize reduced (negative) states.
        e_electrode = -z * self.FARADAY_KJMOL_PER_V * self.voltage_v

        # Total
        e_total = intrinsic_offset + e_electrode + pe_nonbonded + self._solvation_correction(state)
        return e_total

    def _neighbor_states(self, current_state: int) -> List[int]:
        """Restrict proposals to adjacent redox states when possible."""
        ordered_states = sorted(
            self.params.allowed_states,
            key=lambda state: self.params.charge_states[state],
            reverse=True,
        )
        if current_state not in ordered_states:
            return [state for state in ordered_states if state != current_state]

        idx = ordered_states.index(current_state)
        neighbors = []
        if idx > 0:
            neighbors.append(ordered_states[idx - 1])
        if idx < len(ordered_states) - 1:
            neighbors.append(ordered_states[idx + 1])
        return neighbors or [state for state in ordered_states if state != current_state]

    def compute_transition_energy(
        self,
        state_from: int,
        state_to: int,
        pe_solv_from: float,  # NonbondedForce energy in state_from
        pe_solv_to: float,    # NonbondedForce energy in state_to
    ) -> float:
        """
        Compute energy change for transition from state_from to state_to.

        w = E_to - E_from

        Returns energy in kJ/mol. Negative w means favorable transition.
        """
        e_from = self.compute_state_energy(state_from, pe_solv_from)
        e_to = self.compute_state_energy(state_to, pe_solv_to)
        return e_to - e_from

    def attempt_transition(
        self,
        current_state: int,
        pe_solv_from: float,
        pe_solv_to: float,
    ) -> Tuple[bool, Optional[int]]:
        """
        Attempt one MC electron transfer move.

        Randomly select a target state from allowed_states and attempt transition
        using Metropolis criterion.

        Returns:
            (accepted: bool, new_state: int or None)
        """
        # Randomly choose a target state (different from current)
        available_states = self._neighbor_states(current_state)
        if not available_states:
            return False, None

        target_state = random.choice(available_states)

        # Compute energy change
        dE = self.compute_transition_energy(
            current_state, target_state,
            pe_solv_from, pe_solv_to
        )

        # Metropolis criterion
        self.n_attempted += 1
        self.n_attempted_by_state[current_state] = self.n_attempted_by_state.get(current_state, 0) + 1

        if dE < 0.0 or random.random() < np.exp(-dE / self.kT):
            self.n_accepted += 1
            self.n_accepted_by_state[target_state] = self.n_accepted_by_state.get(target_state, 0) + 1
            return True, target_state

        return False, None

    def run_mc_sampling(
        self,
        n_steps: int,
        initial_state: Optional[int] = None,
        pe_solv: float = 0.0,
    ) -> Dict:
        """
        Run a simple MC trajectory at fixed solvation energy.
        """
        current_state = self.reference_state if initial_state is None else initial_state
        for _ in range(n_steps):
            accepted, new_state = self.attempt_transition(
                current_state=current_state,
                pe_solv_from=pe_solv,
                pe_solv_to=pe_solv,
            )
            if accepted and new_state is not None:
                current_state = new_state

            self.state_history.append(current_state)
            self.energy_history.append(
                self.compute_state_energy(current_state, pe_nonbonded=pe_solv)
            )

        return self.summary_dict()

    def get_acceptance_rate(self) -> float:
        """Overall MC acceptance rate."""
        if self.n_attempted == 0:
            return 0.0
        return self.n_accepted / self.n_attempted

    def get_acceptance_rate_by_state(self, state: int) -> float:
        """MC acceptance rate for transitions FROM a given state."""
        n_att = self.n_attempted_by_state.get(state, 0)
        if n_att == 0:
            return 0.0
        return self.n_accepted_by_state.get(state, 0) / n_att

    def reset_statistics(self):
        """Clear MC statistics (useful for re-equilibration at new voltage)."""
        self.n_accepted = 0
        self.n_attempted = 0
        self.n_accepted_by_state = {s: 0 for s in self.params.allowed_states}
        self.n_attempted_by_state = {s: 0 for s in self.params.allowed_states}
        self.state_history = []
        self.energy_history = []

    def get_state_occupancy(self) -> Dict[int, float]:
        """
        Estimate probability of each state from trajectory history.

        Returns dict: {state: fractional_occupancy}
        """
        if not self.state_history:
            return {s: 0.0 for s in self.params.allowed_states}

        counts = {s: 0 for s in self.params.allowed_states}
        for state in self.state_history:
            counts[state] += 1

        total = len(self.state_history)
        return {s: counts[s] / total for s in self.params.allowed_states}

    def get_mean_charge(self) -> float:
        """Estimate mean charge from trajectory."""
        if not self.state_history:
            return 0.0

        total_charge = sum(
            self.params.charge_states[state]
            for state in self.state_history
        )
        return total_charge / len(self.state_history)

    def summary_dict(self) -> Dict:
        """Return summary of MC statistics as dict."""
        occupancy = self.get_state_occupancy()
        return {
            "n_attempted": self.n_attempted,
            "n_accepted": self.n_accepted,
            "acceptance_rate": self.get_acceptance_rate(),
            "temperature_k": self.temperature_k,
            "voltage_v": self.voltage_v,
            "mean_charge_e": self.get_mean_charge(),
            "state_occupancy": occupancy,
            "state_history_length": len(self.state_history),
        }


def load_redox_parameters_from_json(config_file: str) -> Dict[str, RedoxParameters]:
    """
    Load redox parameters for multiple molecules from JSON config file.

    Config format:
    {
        "molecules": {
            "EC": {
                "residue_names": ["ECA", "EC1"],
                "charge_states": {"0": 0.0, "-1": -1.0},
                "electron_affinities": {"0,-1": 0.5},
                "allowed_states": [0, -1],
                "partition_coefficient": 0.0
            },
            ...
        }
    }

    Returns:
        dict: {molecule_name: RedoxParameters, ...}
    """
    with open(config_file, 'r') as f:
        cfg = json.load(f)

    redox_params = {}
    for mol_name, mol_cfg in cfg.get("molecules", {}).items():
        # Parse charge_states: keys are strings from JSON, convert to int
        charge_states = {
            int(k): float(v)
            for k, v in mol_cfg.get("charge_states", {}).items()
        }

        # Parse electron_affinities: keys are "s1,s2" strings
        ea_dict = {}
        for key_str, ea_val in mol_cfg.get("electron_affinities", {}).items():
            s1, s2 = map(int, key_str.split(','))
            ea_dict[(s1, s2)] = float(ea_val)

        redox_params[mol_name] = RedoxParameters(
            name=mol_name,
            residue_names=mol_cfg.get("residue_names", [mol_name]),
            charge_states=charge_states,
            electron_affinities=ea_dict,
            allowed_states=list(map(int, mol_cfg.get("allowed_states", [0]))),
            partition_coefficient=float(mol_cfg.get("partition_coefficient", 0.0)),
            state_free_energy_offsets_kjmol={
                int(k): float(v)
                for k, v in mol_cfg.get("state_free_energy_offsets_kjmol", {}).items()
            } or None,
        )

    return redox_params


if __name__ == "__main__":
    # Example usage
    print("RedoxMC module loaded successfully")
    print(f"Faraday constant: {RedoxMC.FARADAY_KJMOL_PER_V} kJ/mol/V/e")
    print(f"Boltzmann constant: {RedoxMC.KB_KJMOL_PER_K} kJ/mol/K")
