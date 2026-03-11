#!/usr/bin/env python
"""
subroutines_opls.py
--------------------
OPLS fixed-charge constant-potential MD simulation using OpenMM 8.4's
built-in ConstantPotentialForce plugin.

Workflow:
  1. Load PDB + OPLS XML force field files
  2. Create OpenMM System (PME, HBonds)
  3. Collect cathode/anode atom indices by chain index (or residue name)
  4. Transfer all particle charges from NonbondedForce → ConstantPotentialForce
  5. Register electrodes with target potentials
  6. Run equilibration + production MD, logging electrode charges

Reference implementation:
  PhyNEO/example/example_constantQ_interface/Example_OPLS/
  LiPF6_EC_DMC_inert_electrode/run_openmm84_inert_electrode.py
"""

import json
import os
import sys
import numpy as np
from datetime import datetime
from math import sqrt
from pathlib import Path

import openmm as mm
import openmm.app as app
import openmm.unit as unit

# Faraday constant: kJ/mol per elementary_charge per Volt
FARADAY_KJMOL_PER_V = 96.48533212331002


# ---------------------------------------------------------------------------
# OPLS-AA geometric combination rule
# ---------------------------------------------------------------------------

def OPLS_LJ(system):
    """
    Apply OPLS-AA geometric combination rules for LJ interactions.

    OPLS-AA uses geometric mean for both sigma and epsilon:
        sigma_ij   = sqrt(sigma_i   * sigma_j)
        epsilon_ij = sqrt(epsilon_i * epsilon_j)

    OpenMM's default NonbondedForce uses Lorentz-Berthelot mixing
    (arithmetic mean for sigma, geometric mean for epsilon), which is
    incorrect for OPLS-AA. This function:
      1. Creates a CustomNonbondedForce with the correct geometric sigma rule.
      2. Moves all LJ parameters (sigma, epsilon) to the CustomNonbondedForce.
      3. Zeroes epsilon in the original NonbondedForce (retaining charges and
         sigma so that 1-4 exception terms remain correct).
      4. Sets proper 1-4 LJ exceptions using geometric sigma.

    Must be called AFTER createSystem() and BEFORE _setup_cpf().

    Reference: OpenMM OPLS-AA tutorial / ParmEd OPLS examples.
    """
    forces = {
        system.getForce(index).__class__.__name__: system.getForce(index)
        for index in range(system.getNumForces())
    }
    nonbonded_force = forces['NonbondedForce']
    lorentz = mm.CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); '
        'sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)'
    )
    # CustomNonbondedForce does not support PME/Ewald; map to CutoffPeriodic
    nb_method = nonbonded_force.getNonbondedMethod()
    if nb_method in (mm.NonbondedForce.PME,
                     mm.NonbondedForce.Ewald,
                     mm.NonbondedForce.CutoffPeriodic):
        lorentz.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    elif nb_method == mm.NonbondedForce.CutoffNonPeriodic:
        lorentz.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    else:
        lorentz.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    system.addForce(lorentz)

    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        # Store as plain floats (nm, kJ/mol) so math.sqrt can operate on them
        sig_nm  = sigma.value_in_unit(unit.nanometer)
        eps_kj  = epsilon.value_in_unit(unit.kilojoule_per_mole)
        LJset[index] = (sig_nm, eps_kj)
        lorentz.addParticle([sigma, epsilon])
        # Zero epsilon in NonbondedForce; charges and sigma stay for exceptions
        nonbonded_force.setParticleParameters(index, charge, sigma, epsilon * 0)

    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # Exclude this pair from CustomNonbondedForce (handled by NonbondedForce exceptions)
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            # 1-4 pair: recompute sigma and epsilon with geometric mean
            sig14 = sqrt(LJset[p1][0] * LJset[p2][0])  # nm
            eps14 = sqrt(LJset[p1][1] * LJset[p2][1])  # kJ/mol
            nonbonded_force.setExceptionParameters(
                i, p1, p2, q,
                sig14 * unit.nanometer,
                eps14 * unit.kilojoule_per_mole
            )

    return system


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_nonbonded_force(system):
    """Return the single NonbondedForce object from an OpenMM System."""
    nb_forces = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)]
    if len(nb_forces) != 1:
        raise ValueError(
            f"Expected exactly 1 NonbondedForce, found {len(nb_forces)}. "
            "Make sure the OPLS XML files define a standard PME system."
        )
    return nb_forces[0]


def collect_electrode_atoms_by_chain(topology, cathode_chain_idx, anode_chain_idx):
    """
    Collect electrode atom indices by chain index.
    Excludes hydrogen atoms (element H) which are typically absent in
    frozen electrode models but included for safety.

    Parameters
    ----------
    topology : openmm.app.Topology
    cathode_chain_idx : int
        Index of the chain whose atoms form the cathode (chain 0 by default).
    anode_chain_idx : int
        Index of the chain whose atoms form the anode (chain 1 by default).

    Returns
    -------
    cath_atoms : list[int]
    ano_atoms  : list[int]
    """
    cath_atoms, ano_atoms = [], []
    for chain in topology.chains():
        if chain.index == cathode_chain_idx:
            cath_atoms.extend(
                a.index for a in chain.atoms()
                if a.element is not None and a.element.symbol != "H"
            )
        elif chain.index == anode_chain_idx:
            ano_atoms.extend(
                a.index for a in chain.atoms()
                if a.element is not None and a.element.symbol != "H"
            )
    return cath_atoms, ano_atoms


def collect_electrode_atoms_by_residue(topology, cathode_resname="CAT", anode_resname="ANO"):
    """
    Alternative: collect electrode atoms by residue name.
    Useful when the PDB places electrode and electrolyte atoms in the same chain.
    """
    cath_atoms, ano_atoms = [], []
    for res in topology.residues():
        if res.name == cathode_resname:
            cath_atoms.extend(a.index for a in res.atoms())
        elif res.name == anode_resname:
            ano_atoms.extend(a.index for a in res.atoms())
    return cath_atoms, ano_atoms


# ---------------------------------------------------------------------------
# Main simulation class
# ---------------------------------------------------------------------------

class OPLSMDSimulation:
    """
    Constant-potential MD simulation using the OPLS fixed-charge force field
    and OpenMM 8.4's ConstantPotentialForce.

    Parameters
    ----------
    config_path : str
        Path to the JSON configuration file (see config_opls.json for schema).
    """

    def __init__(self, config_path):
        config_path = Path(config_path).resolve()
        self.config_path = config_path
        self.config_dir = config_path.parent
        with open(config_path, "r") as f:
            self.cfg = json.load(f)

        print(f"[OPLSMDSimulation] Loading config from: {config_path}")
        print(f"[OPLSMDSimulation] {datetime.now()}")

        # ------------------------------------------------------------------ #
        # 1. Load PDB
        # ------------------------------------------------------------------ #
        pdb_path = self._resolve_input_path(self.cfg["pdb"])
        print(f"[OPLSMDSimulation] Reading PDB: {pdb_path}")
        self.pdb = app.PDBFile(pdb_path)

        # Optionally set box vectors from config if the PDB has no CRYST1 record
        cell_cfg = self.cfg.get("cell")
        if cell_cfg is not None:
            lx = cell_cfg["box_x_nm"] * unit.nanometer
            ly = cell_cfg["box_y_nm"] * unit.nanometer
            lz = cell_cfg["box_z_nm"] * unit.nanometer
            self.pdb.topology.setUnitCellDimensions(
                mm.Vec3(lx / unit.nanometer,
                        ly / unit.nanometer,
                        lz / unit.nanometer) * unit.nanometer
            )
            print(f"[OPLSMDSimulation] Box set from config: {lx}, {ly}, {lz}")

        # ------------------------------------------------------------------ #
        # 2. Load force field XML files
        # ------------------------------------------------------------------ #
        ff_files = [self._resolve_input_path(path) for path in self.cfg["forcefield_xml"]]
        print(f"[OPLSMDSimulation] Force field files: {ff_files}")
        self.ff = app.ForceField(*ff_files)

        # ------------------------------------------------------------------ #
        # 3. Create OpenMM System
        # ------------------------------------------------------------------ #
        md_cfg = self.cfg["md"]
        cutoff_nm = md_cfg.get("nonbonded_cutoff_nm", 1.2)
        print(f"[OPLSMDSimulation] Creating system (PME, HBonds, cutoff={cutoff_nm} nm)...")
        self.system = self.ff.createSystem(
            self.pdb.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=cutoff_nm * unit.nanometer,
            constraints=app.HBonds,
            rigidWater=False,
            ewaldErrorTolerance=5e-4,
            removeCMMotion=False,
        )

        # ------------------------------------------------------------------ #
        # 4. Apply OPLS-AA geometric LJ combination rules
        #    Must be done BEFORE extracting NonbondedForce and before CPF setup,
        #    because OPLS_LJ zeroes epsilon in NonbondedForce and adds a
        #    CustomNonbondedForce with geometric sigma/epsilon mixing.
        # ------------------------------------------------------------------ #
        print("[OPLSMDSimulation] Applying OPLS-AA geometric LJ combination rule ...")
        self.system = OPLS_LJ(self.system)

        # ------------------------------------------------------------------ #
        # 5. Get NonbondedForce (charges intact; epsilon zeroed by OPLS_LJ)
        # ------------------------------------------------------------------ #
        self.nb = get_nonbonded_force(self.system)

        # ------------------------------------------------------------------ #
        # 6. Collect electrode atom indices
        # ------------------------------------------------------------------ #
        ele_cfg = self.cfg["electrode"]
        use_chain = ele_cfg.get("use_chain_index", True)

        if use_chain:
            cath_idx = ele_cfg.get("cathode_chain_index", 0)
            ano_idx  = ele_cfg.get("anode_chain_index",   1)
            self.cath_atoms, self.ano_atoms = collect_electrode_atoms_by_chain(
                self.pdb.topology, cath_idx, ano_idx
            )
        else:
            cath_res = ele_cfg.get("cathode_resname", "CAT")
            ano_res  = ele_cfg.get("anode_resname",   "ANO")
            self.cath_atoms, self.ano_atoms = collect_electrode_atoms_by_residue(
                self.pdb.topology, cath_res, ano_res
            )

        print(f"[OPLSMDSimulation] Cathode atoms: {len(self.cath_atoms)}")
        print(f"[OPLSMDSimulation] Anode   atoms: {len(self.ano_atoms)}")
        if not self.cath_atoms or not self.ano_atoms:
            raise ValueError(
                "No electrode atoms found! Check cathode_chain_index / anode_chain_index "
                "in config, or switch to use_chain_index=false and set cathode_resname/anode_resname."
            )

        # ------------------------------------------------------------------ #
        # 7. Set up ConstantPotentialForce and transfer charges from NB
        # ------------------------------------------------------------------ #
        self._setup_cpf()

        # ------------------------------------------------------------------ #
        # 8. Integrator and Simulation
        # ------------------------------------------------------------------ #
        T   = md_cfg.get("temperature_k", 300.0)
        fric = md_cfg.get("friction_ps", 1.0)
        dt  = md_cfg.get("timestep_fs", 2.0)
        self.temperature_k = T

        self.integrator = mm.LangevinMiddleIntegrator(
            T    * unit.kelvin,
            fric / unit.picosecond,
            dt   * unit.femtosecond,
        )

        platform_name = md_cfg.get("platform", "CUDA")
        try:
            platform = mm.Platform.getPlatformByName(platform_name)
            print(f"[OPLSMDSimulation] Using platform: {platform_name}")
        except Exception:
            print(f"[OPLSMDSimulation] Platform '{platform_name}' not available, falling back to CPU.")
            platform = mm.Platform.getPlatformByName("CPU")

        self.sim = app.Simulation(
            self.pdb.topology,
            self.system,
            self.integrator,
            platform,
        )
        self.sim.context.setPositions(self.pdb.positions)

        # ------------------------------------------------------------------ #
        # 9. Reporters
        # ------------------------------------------------------------------ #
        run_cfg = self.cfg.get("run", {})
        report_interval = run_cfg.get("report_interval", 1000)
        dcd_file  = run_cfg.get("dcd_file",  "md_opls.dcd")
        log_file  = run_cfg.get("log_file",  "md_opls.log")
        self.charges_file = run_cfg.get("charges_file", "electrode_charges.log")
        self.report_interval = report_interval

        self.sim.reporters.append(app.DCDReporter(dcd_file, report_interval))
        self.sim.reporters.append(
            app.StateDataReporter(
                log_file,
                report_interval,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                density=True,
                progress=True,
                remainingTime=True,
                speed=True,
                totalSteps=run_cfg.get("equil_steps", 0) + run_cfg.get("prod_steps", 0),
                separator="\t",
            )
        )
        print(f"[OPLSMDSimulation] Reporters: DCD → {dcd_file}, log → {log_file}")
        print(f"[OPLSMDSimulation] Electrode charges → {self.charges_file}")
        print(f"[OPLSMDSimulation] Initialization complete.\n")

    def _resolve_input_path(self, raw_path):
        path = Path(raw_path)
        if path.is_absolute():
            return str(path)
        candidate = (self.config_dir / path).resolve()
        if candidate.exists():
            return str(candidate)
        return str((Path.cwd() / path).resolve())

    # ---------------------------------------------------------------------- #
    # Internal: CPF setup
    # ---------------------------------------------------------------------- #

    def _setup_cpf(self):
        """
        Build ConstantPotentialForce:
          - Transfer every particle charge from NonbondedForce → CPF
          - Zero out NonbondedForce charges (LJ parameters unchanged)
          - Transfer 1-4 exception charge products to CPF
          - Register cathode and anode electrodes
          - Add CPF to system
        """
        ele_cfg = self.cfg["electrode"]

        self.cpf = mm.ConstantPotentialForce()

        # Solver settings
        self.cpf.setCutoffDistance(
            float(self.cfg["md"].get("nonbonded_cutoff_nm", 1.2))
        )
        self.cpf.setConstantPotentialMethod(mm.ConstantPotentialForce.CG)
        self.cpf.setCGErrorTolerance(float(ele_cfg.get("cg_error_tol", 1e-4)))

        use_constraint = bool(ele_cfg.get("use_charge_constraint", False))
        self.cpf.setUseChargeConstraint(use_constraint)
        if use_constraint:
            self.cpf.setChargeConstraintTarget(
                float(ele_cfg.get("charge_constraint_target_e", 0.0))
            )

        # Transfer particle charges
        print("[OPLSMDSimulation] Transferring charges: NonbondedForce → CPF ...")
        n_particles = self.system.getNumParticles()
        for i in range(n_particles):
            q, sig, eps = self.nb.getParticleParameters(i)
            q_e = q.value_in_unit(unit.elementary_charge)
            self.cpf.addParticle(float(q_e))
            self.nb.setParticleParameters(i, 0.0 * unit.elementary_charge, sig, eps)

        # Transfer exception charge products (1-4 and excluded pairs)
        print("[OPLSMDSimulation] Transferring exceptions: NonbondedForce → CPF ...")
        for j in range(self.nb.getNumExceptions()):
            p1, p2, qprod, sig, eps = self.nb.getExceptionParameters(j)
            qprod_e2 = qprod.value_in_unit(unit.elementary_charge ** 2)
            self.cpf.addException(int(p1), int(p2), float(qprod_e2))
            self.nb.setExceptionParameters(
                j, p1, p2, 0.0 * unit.elementary_charge ** 2, sig, eps
            )

        # Electrode potentials (convert Volts → kJ/mol/e).
        # The configured cell voltage is applied symmetrically across electrodes.
        voltage_v = float(ele_cfg.get("voltage_v", 2.0))
        half_cell_voltage_v = 0.5 * voltage_v
        cath_pot  =  half_cell_voltage_v * FARADAY_KJMOL_PER_V
        ano_pot   = -half_cell_voltage_v * FARADAY_KJMOL_PER_V

        gauss_width = float(ele_cfg.get("gaussian_width_nm", 0.20))
        tf_scale    = float(ele_cfg.get("thomas_fermi_scale_invnm", 5.0))

        print(f"[OPLSMDSimulation] Adding electrodes: V={voltage_v} V total, "
              f"cathode pot={cath_pot:.2f} kJ/mol/e, anode pot={ano_pot:.2f} kJ/mol/e")

        self.cpf.addElectrode(set(self.cath_atoms), cath_pot, gauss_width, tf_scale)
        self.cpf.addElectrode(set(self.ano_atoms),  ano_pot,  gauss_width, tf_scale)

        self.system.addForce(self.cpf)
        print("[OPLSMDSimulation] ConstantPotentialForce added to system.")

    # ---------------------------------------------------------------------- #
    # Public: log electrode charges for one step
    # ---------------------------------------------------------------------- #

    def _log_charges(self, step, charges_fh):
        charges = list(self.cpf.getCharges(self.sim.context))
        q_c = sum(
            charges[i].value_in_unit(unit.elementary_charge)
            for i in self.cath_atoms
        )
        q_a = sum(
            charges[i].value_in_unit(unit.elementary_charge)
            for i in self.ano_atoms
        )
        charges_fh.write(f"{step}\t{q_c:.6f}\t{q_a:.6f}\n")
        charges_fh.flush()

    # ---------------------------------------------------------------------- #
    # Public: equilibration
    # ---------------------------------------------------------------------- #

    def run_equilibration(self):
        run_cfg = self.cfg.get("run", {})
        equil_steps = run_cfg.get("equil_steps", 0)
        if equil_steps <= 0:
            print("[OPLSMDSimulation] No equilibration requested.")
            return

        print(f"[OPLSMDSimulation] Minimizing energy ...")
        self.sim.minimizeEnergy(maxIterations=1000)

        print(f"[OPLSMDSimulation] Equilibrating for {equil_steps} steps ...")
        with open(self.charges_file, "w") as fh:
            fh.write("# step\tQ_cathode_e\tQ_anode_e\n")
            for step in range(0, equil_steps, self.report_interval):
                self.sim.step(self.report_interval)
                self._log_charges(step + self.report_interval, fh)
        print(f"[OPLSMDSimulation] Equilibration done. {datetime.now()}")

    # ---------------------------------------------------------------------- #
    # Public: production
    # ---------------------------------------------------------------------- #

    def run_production(self):
        run_cfg = self.cfg.get("run", {})
        prod_steps    = run_cfg.get("prod_steps", 1000000)
        equil_steps   = run_cfg.get("equil_steps", 0)

        print(f"[OPLSMDSimulation] Production run: {prod_steps} steps ...")
        with open(self.charges_file, "a") as fh:
            for step in range(0, prod_steps, self.report_interval):
                self.sim.step(self.report_interval)
                global_step = equil_steps + step + self.report_interval
                self._log_charges(global_step, fh)

        # Save final structure
        final_pdb = run_cfg.get("final_pdb", "final_opls.pdb")
        state = self.sim.context.getState(getPositions=True)
        with open(final_pdb, "w") as fh:
            app.PDBFile.writeFile(self.sim.topology, state.getPositions(), fh)
        print(f"[OPLSMDSimulation] Production done. Final structure: {final_pdb}")
        print(f"[OPLSMDSimulation] {datetime.now()}")
