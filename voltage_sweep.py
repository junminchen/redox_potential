"""
voltage_sweep.py
----------------
Voltage sweep simulation for generating CV/LSV curves from constant-potential MD.

Performs a series of constant-potential MD simulations across a range of electrode
potentials, tracking redox-active molecule charge state distribution and computing
simulated "current" proxy via d(n_redox)/dV.

Allows direct comparison with experimental CV/LSV data.
"""

import os
import json
from collections import Counter
import numpy as np
import pandas as pd
from datetime import datetime
from typing import Dict, List, Tuple, Optional

import openmm as mm
import openmm.app as app
import openmm.unit as unit

from redox_mc import RedoxMC, RedoxParameters, load_redox_parameters_from_json
from subroutines_opls import OPLSMDSimulation


class VoltageSweepSimulation:
    """
    Perform constant-potential MD simulations across a voltage range.

    At each voltage point:
    1. Run constant-potential MD with ConstantPotentialForce
    2. Sample redox-active molecule state distribution via RedoxMC
    3. Track mean charge and occupancy
    4. Estimate "current" as d(occupancy)/dV

    Results are saved as a DataFrame ready for CV analysis.
    """

    def __init__(
        self,
        pdb_file: str,
        config_file: str,
        redox_config_file: str,
        output_dir: str = "redox_sweep_results",
    ):
        """
        Initialize voltage sweep simulation.

        Args:
            pdb_file: path to PDB with electrodes + electrolyte
            config_file: OPLS MD configuration (from OPLSMDSimulation)
            redox_config_file: redox molecule parameters (EC, DMC, etc.)
            output_dir: directory for results
        """
        self.pdb_file = pdb_file
        self.config_file = config_file
        self.redox_config_file = redox_config_file
        self.output_dir = output_dir

        os.makedirs(output_dir, exist_ok=True)

        # Load redox parameters
        self.redox_params = load_redox_parameters_from_json(redox_config_file)
        print(f"[VoltageSweep] Loaded redox parameters for: {list(self.redox_params.keys())}")

        self.results = []  # Will store data as dicts

    @staticmethod
    def _voltage_label(voltage_v: float) -> str:
        """Return a filesystem-safe label for one voltage point."""
        return f"V{voltage_v:+.3f}".replace("+", "p").replace("-", "m")

    def run_sweep(
        self,
        v_start: float = -2.0,       # V vs vacuum
        v_end: float = 0.0,
        v_step: float = 0.1,
        equil_steps_per_v: int = 5000,
        sample_steps_per_v: int = 5000,
        report_interval: int = 1000,
        platform: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Run voltage sweep simulation.

        Args:
            v_start: starting voltage (V vs vacuum)
            v_end: ending voltage (V vs vacuum)
            v_step: voltage step size (V)
            equil_steps_per_v: equilibration MD steps at each voltage
            sample_steps_per_v: production sampling steps at each voltage
            report_interval: reporting interval for MD
            platform: optional OpenMM platform override ("CUDA", "OpenCL", "CPU")

        Returns:
            pandas.DataFrame with sweep results
        """
        voltages = np.arange(v_start, v_end + v_step / 2, v_step)
        print(f"[VoltageSweep] Scanning {len(voltages)} voltages: {v_start:.2f} to {v_end:.2f} V (step {v_step:.2f} V)")
        print(f"[VoltageSweep] Per voltage: {equil_steps_per_v} equil + {sample_steps_per_v} sample steps")

        for i, voltage_v in enumerate(voltages):
            print(f"\n{'='*70}")
            print(f"[VoltageSweep] Point {i+1}/{len(voltages)}: V = {voltage_v:.3f} V vs vacuum")
            print(f"{'='*70}")

            result = self._simulate_at_voltage(
                voltage_v,
                equil_steps_per_v,
                sample_steps_per_v,
                report_interval,
                platform,
            )
            self.results.append(result)

            # Checkpoint results
            self._save_checkpoint()

        # Convert to DataFrame and compute derivatives
        df = pd.DataFrame(self.results)
        df = self._compute_currents(df)

        # Save final results
        output_file = os.path.join(self.output_dir, "voltage_sweep_results.csv")
        df.to_csv(output_file, index=False)
        print(f"\n[VoltageSweep] Final results saved to {output_file}")

        return df

    def _simulate_at_voltage(
        self,
        voltage_v: float,
        equil_steps: int,
        sample_steps: int,
        report_interval: int,
        platform: Optional[str],
    ) -> Dict:
        """
        Run MD and MC sampling at a single voltage point.

        Returns:
            dict with voltage, occupancy, charges, MC statistics
        """
        # Load MD config and update voltage
        with open(self.config_file, 'r') as f:
            cfg = json.load(f)

        cfg["electrode"]["voltage_v"] = voltage_v
        if platform is not None:
            cfg["md"]["platform"] = platform
        cfg["run"]["equil_steps"] = equil_steps
        cfg["run"]["prod_steps"] = sample_steps
        cfg["run"]["report_interval"] = report_interval

        point_dir = os.path.abspath(
            os.path.join(self.output_dir, self._voltage_label(voltage_v))
        )
        os.makedirs(point_dir, exist_ok=True)
        cfg["run"]["dcd_file"] = os.path.join(point_dir, "md_opls.dcd")
        cfg["run"]["log_file"] = os.path.join(point_dir, "md_opls.log")
        cfg["run"]["charges_file"] = os.path.join(point_dir, "electrode_charges.log")
        cfg["run"]["final_pdb"] = os.path.join(point_dir, "final_opls.pdb")

        # Temporary config file for this voltage point
        temp_config = os.path.join(point_dir, "config.json")
        with open(temp_config, 'w') as f:
            json.dump(cfg, f, indent=2)

        # Run OPLS MD simulation
        print(f"  → Running OPLS MD at V={voltage_v:.3f} V...")
        sim = OPLSMDSimulation(temp_config)
        sim.run_equilibration()
        sim.run_production()

        # Extract results from simulation
        result = {
            "voltage_v": voltage_v,
            "temperature_k": cfg["md"]["temperature_k"],
            "timestamp": datetime.now().isoformat(),
        }

        # Read electrode charges from log
        charges_log = cfg["run"]["charges_file"]
        if os.path.exists(charges_log):
            charges_df = pd.read_csv(
                charges_log,
                sep=r"\s+",
                comment="#",
                names=["step", "Q_cathode_e", "Q_anode_e"],
                header=None,
            )
            # Skip equilibration, average over production steps
            prod_charges = charges_df.iloc[-sample_steps // report_interval:]
            result["q_cathode_mean"] = prod_charges["Q_cathode_e"].mean()
            result["q_anode_mean"] = prod_charges["Q_anode_e"].mean()
            result["q_cathode_std"] = prod_charges["Q_cathode_e"].std()
        else:
            result["q_cathode_mean"] = 0.0
            result["q_anode_mean"] = 0.0
            result["q_cathode_std"] = 0.0

        # MC sampling for redox molecules
        residue_counts = Counter(res.name for res in sim.pdb.topology.residues())
        mc_steps = max(200, min(2000, sample_steps // max(report_interval, 1) * 100))
        for mol_name, redox_params in self.redox_params.items():
            residue_names = set(redox_params.residue_names)
            residue_names.add(redox_params.name)
            n_molecules = sum(
                count
                for res_name, count in residue_counts.items()
                if res_name in residue_names
            )

            result[f"{mol_name}_n_molecules"] = n_molecules
            if n_molecules == 0:
                print(f"  [RedoxMC] {mol_name}: no matching residues found in topology")
                result[f"{mol_name}_n_mc_accepted"] = 0
                result[f"{mol_name}_n_mc_attempted"] = 0
                result[f"{mol_name}_fraction_reduced"] = 0.0
                result[f"{mol_name}_mean_charge"] = 0.0
                continue

            redox_mc = RedoxMC(
                redox_params,
                temperature_k=cfg["md"]["temperature_k"],
                voltage_v=voltage_v,
            )
            summary = redox_mc.run_mc_sampling(mc_steps)
            reduced_occupancy = sum(
                frac
                for state, frac in summary["state_occupancy"].items()
                if redox_params.charge_states[state] < redox_mc.reference_charge
            )

            result[f"{mol_name}_n_mc_accepted"] = summary["n_accepted"] * n_molecules
            result[f"{mol_name}_n_mc_attempted"] = summary["n_attempted"] * n_molecules
            result[f"{mol_name}_fraction_reduced"] = reduced_occupancy
            result[f"{mol_name}_mean_charge"] = summary["mean_charge_e"]
            print(
                f"  [RedoxMC] {mol_name}: n={n_molecules}, "
                f"fraction_reduced={reduced_occupancy:.3f}, "
                f"mean_charge={summary['mean_charge_e']:.3f} e"
            )

        print(f"  ✓ Completed: Q_cath={result['q_cathode_mean']:.3f} e, Q_anode={result['q_anode_mean']:.3f} e")

        return result

    def _save_checkpoint(self):
        """Save intermediate results."""
        df = pd.DataFrame(self.results)
        checkpoint_file = os.path.join(self.output_dir, "voltage_sweep_checkpoint.csv")
        df.to_csv(checkpoint_file, index=False)

    def _compute_currents(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute simulated "current" as numerical derivative of occupancy.

        I_sim(V) ∝ d(fraction_reduced) / dV

        Adds columns:
        - simulated_current_proxy_ec, simulated_current_proxy_dmc, etc.
        """
        for mol_name in self.redox_params.keys():
            if f"{mol_name}_fraction_reduced" in df.columns:
                frac_col = df[f"{mol_name}_fraction_reduced"].values
                # Numerical derivative
                if len(df) > 1:
                    voltage_col = df["voltage_v"].values
                    dI = np.gradient(frac_col, voltage_col)
                    df[f"{mol_name}_simulated_current_proxy"] = dI
                else:
                    df[f"{mol_name}_simulated_current_proxy"] = 0.0

        return df


class RedoxPotentialAnalyzer:
    """
    Analyze voltage sweep results to extract redox potentials.

    Methods:
    - Find half-wave potential (E_1/2) where fraction_reduced = 0.5
    - Compare with LSV experimental data
    - Estimate thermodynamic redox potential
    """

    @staticmethod
    def find_half_wave_potential(
        sweep_df: pd.DataFrame,
        mol_name: str,
        voltage_col: str = "voltage_v",
        fraction_col_name: Optional[str] = None,
    ) -> float:
        """
        Find voltage where fraction_reduced = 0.5 (half-wave potential).

        Args:
            sweep_df: DataFrame from VoltageSweep
            mol_name: molecule name (e.g., "EC", "DMC")
            voltage_col: column name for voltage
            fraction_col_name: optional override for fraction column name

        Returns:
            half-wave potential (V vs vacuum)
        """
        if fraction_col_name is None:
            fraction_col_name = f"{mol_name}_fraction_reduced"

        if fraction_col_name not in sweep_df.columns:
            print(f"Warning: column '{fraction_col_name}' not found")
            return np.nan

        voltages = sweep_df[voltage_col].values
        fractions = sweep_df[fraction_col_name].values

        # Linear interpolation to find E where fraction = 0.5
        if len(voltages) < 2:
            return np.nan

        for i in range(len(voltages) - 1):
            if fractions[i] == 0.5:
                return voltages[i]
            if fractions[i + 1] == 0.5:
                return voltages[i + 1]
            if (fractions[i] - 0.5) * (fractions[i+1] - 0.5) < 0:
                # Root between i and i+1
                v1, v2 = voltages[i], voltages[i+1]
                f1, f2 = fractions[i], fractions[i+1]
                e_half = v1 + (0.5 - f1) * (v2 - v1) / (f2 - f1)
                return e_half

        return np.nan

    @staticmethod
    def convert_to_li_reference(
        voltage_v_vs_vacuum: float,
        offset_ev: float = 1.4,  # typical for organic electrolyte
    ) -> float:
        """
        Convert voltage from vacuum reference to Li/Li+ reference.

        In organic electrolyte, Li/Li+ is typically about +1.4 V above the
        vacuum-referenced value used here.

        Args:
            voltage_v_vs_vacuum: voltage in V vs vacuum
            offset_ev: offset in eV (default 1.4 for org. electrolyte)

        Returns:
            voltage in V vs Li/Li+
        """
        return voltage_v_vs_vacuum + offset_ev

    @staticmethod
    def compare_with_lsv(
        sweep_df: pd.DataFrame,
        mol_name: str,
        lsv_peak_potential_v_vs_li: float,
        reference_offset_ev: float = 1.4,
    ) -> Dict:
        """
        Compare simulated redox potential with LSV experimental data.

        Args:
            sweep_df: DataFrame from VoltageSweep
            mol_name: molecule name
            lsv_peak_potential_v_vs_li: experimental LSV peak (V vs Li/Li+)
            reference_offset_ev: offset for vacuum→Li/Li+ conversion

        Returns:
            dict with comparison results
        """
        e_half_vacuum = RedoxPotentialAnalyzer.find_half_wave_potential(
            sweep_df, mol_name
        )
        e_half_li = RedoxPotentialAnalyzer.convert_to_li_reference(
            e_half_vacuum, reference_offset_ev
        )

        diff_mv = (e_half_li - lsv_peak_potential_v_vs_li) * 1000

        return {
            "molecule": mol_name,
            "e_half_wave_v_vs_vacuum": e_half_vacuum,
            "e_half_wave_v_vs_li": e_half_li,
            "lsv_peak_v_vs_li": lsv_peak_potential_v_vs_li,
            "difference_mv": diff_mv,
            "agreement_excellent": abs(diff_mv) < 50,
            "agreement_good": abs(diff_mv) < 100,
        }


def main():
    """Example usage."""
    print("VoltageSweep module loaded successfully")
    print(f"Start time: {datetime.now()}")


if __name__ == "__main__":
    main()
