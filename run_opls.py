#!/usr/bin/env python
"""
run_opls.py
-----------
Entry point for OPLS fixed-charge constant-potential MD simulation.

Usage:
    python run_opls.py config_opls.json

The JSON config specifies:
  - pdb          : path to the input PDB (electrode chains + electrolyte)
  - forcefield_xml: list of XML force field files (OPLS salt, solvent, electrode)
  - electrode    : voltage, gaussian_width, thomas_fermi_scale, solver settings
  - md           : temperature, timestep, platform, cutoff
  - run          : equil_steps, prod_steps, report_interval, output filenames

See config_opls.json for a complete example.

Force field requirements:
  - opls_salt.xml   : LiPF6-family salts (Li+, PF6-, FSI-, TFSI-, etc.)
  - opls_solvent.xml: carbonate/ether solvents (EC, DMC, FEC, EA3F, ...)
  - electrode_residues.xml : CAT / ANO residue topology (in ffdir/)
  - electrode_ff.xml       : E-CAT / E-ANO force field parameters (in ffdir/)

Electrode PDB conventions:
  - Cathode atoms → chain 0, residue name "CAT", atom name "CA"
  - Anode   atoms → chain 1, residue name "ANO", atom name "AN"
  - Electrolyte   → remaining chains
  (can be overridden with use_chain_index=false in electrode config)
"""

import argparse
import sys
from subroutines_opls import OPLSMDSimulation


def main():
    parser = argparse.ArgumentParser(
        description="OPLS constant-potential MD with OpenMM 8.4 ConstantPotentialForce"
    )
    parser.add_argument(
        "config",
        type=str,
        help="Path to JSON configuration file (e.g. config_opls.json)",
    )
    args = parser.parse_args()

    sim = OPLSMDSimulation(args.config)
    sim.run_equilibration()
    sim.run_production()


if __name__ == "__main__":
    main()
