# References

This file centralizes the literature and documentation links already used inside this repository.

Use it for:

- slides and reports
- manuscript drafting
- tracing which paper supports which screening assumption

## 1. Core Electrolyte Screening Anchors

### R1. EC / DMC reduction thermodynamics and microkinetics

- Title: *Atomic thermodynamics and microkinetics of the reduction mechanism of electrolyte additives to facilitate the formation of solid electrolyte interphases in lithium-ion batteries*
- URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC9052788/
- Used for:
  - EC reduction target near `0.47 V vs Li/Li+`
  - DMC reduction target near `0.38 V vs Li/Li+`
  - reduction-side ordering in the effective screening model
- Used in project files:
  - `docs/provenance/REDOX_PARAMETER_PROVENANCE.md`
  - `docs/provenance/redox_parameter_provenance.csv`
  - `docs/provenance/EXPERIMENTAL_VALIDATION_AND_INFERENCE_CN.md`

### R2. EC-free / high-voltage electrolyte review

- Title: *Ethylene Carbonate-Free Electrolytes for High Voltage Lithium-Ion Batteries: Progress and Perspectives*
- URL: https://pubs.acs.org/doi/10.1021/acs.energyfuels.4c02728
- Used for:
  - EC becoming problematic in high-voltage carbonate electrolytes
  - oxidation-side ordering and pathway naming
  - qualitative support that EC is a high-voltage weak point
- Used in project files:
  - `docs/provenance/REDOX_PARAMETER_PROVENANCE.md`
  - `docs/provenance/redox_parameter_provenance.csv`
  - `docs/provenance/EXPERIMENTAL_VALIDATION_AND_INFERENCE_CN.md`

### R3. Experimental high-voltage LSV/CV context for carbonate-like electrolytes

- Title: *Enhanced high voltage performance of LiNi0.5Mn0.3Co0.2O2 cathode via the synergistic effect of LiPO2F2 and FEC in fluorinated electrolyte for lithium-ion batteries*
- URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC8695088/
- Used for:
  - experimental oxidation-response context in the `3.0-6.0 V` region
  - LSV/CV-like interpretation for high-voltage screening
  - DMC-related oxidation pathway placeholders in the current prototype
- Used in project files:
  - `docs/provenance/REDOX_PARAMETER_PROVENANCE.md`
  - `docs/provenance/redox_parameter_provenance.csv`
  - `docs/provenance/EXPERIMENTAL_VALIDATION_AND_INFERENCE_CN.md`

## 2. Example Formulation Literature

### F1. 4 m LiFSI in DME

- Title / source: *Scientific Reports* example for concentrated `LiFSI/DME`
- URL: https://www.nature.com/articles/s41598-017-16268-7
- Supporting mirror: https://www.osti.gov/biblio/1364007
- Used for:
  - Li-metal-friendly concentrated ether example
  - screening preset `4 m LiFSI in DME`
- Used in project files:
  - `docs/notes/EXAMPLE_FORMULATIONS_NOTES.md`
  - `configs/redox/config_redox_lm_4m_lifsi_dme.json`

### F2. LiFSI:EC 1:4

- Title / source: ACS Applied Energy Materials / PMC example for concentrated `LiFSI:EC`
- URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC8790720/
- Used for:
  - concentrated EC-based high-voltage example
  - screening preset `LiFSI:EC 1:4`
- Used in project files:
  - `docs/notes/EXAMPLE_FORMULATIONS_NOTES.md`
  - `configs/redox/config_redox_hce_lifsi_ec14.json`
  - `configs/oxidation/config_oxidation_hce_lifsi_ec14.json`

### F3. 1:1.1 LiFSA / LiFSI in DMC

- Title / source: *Nature Communications* / PMC example for superconcentrated `LiFSA(DMC)`-type electrolyte
- URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC4931331/
- Used for:
  - 5.2 V-class high-voltage screening example
  - screening preset `1:1.1 LiFSI in DMC`
- Used in project files:
  - `docs/notes/EXAMPLE_FORMULATIONS_NOTES.md`
  - `configs/redox/config_redox_hv_11_lifsi_dmc.json`
  - `configs/oxidation/config_oxidation_hv_11_lifsi_dmc.json`

### F4. Conventional 1 M LiPF6 in EC:EMC baseline

- Title / source: recent battery-study example using `1 M LiPF6 in EC:EMC (3:7)`
- URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC11428382/
- Used for:
  - conventional carbonate baseline context
  - screening baseline preset
- Used in project files:
  - `docs/notes/EXAMPLE_FORMULATIONS_NOTES.md`
  - `configs/redox/config_redox_baseline_lipf6_ec_emc37.json`

## 3. Software and Method Documentation

### S1. OpenMM documentation

- URL: https://openmm.org/documentation/
- Used for:
  - API reference in workflow/history docs
- Used in project files:
  - `docs/reports/TASK_3_4_5_COMPLETION.md`

### S2. OpenMM migration guide

- URL: https://openmm.org/documentation/latest/migration/
- Used for:
  - migration notes and compatibility references
- Used in project files:
  - `docs/reports/PRD.md`
  - `docs/methods/REDOX_WORKFLOW.md`

## 4. Classical References Mentioned in Docs but Without URL in Repo

These are already cited in prose in the repository, but no direct URL is currently stored in the docs.

- Reed, Madden, Papadopoulos (2007), *J. Chem. Phys.* 126, 084704
- Frenkel and Smit, *Understanding Molecular Simulation*, 2nd ed.
- Bard and Faulkner, *Electrochemical Methods*, 2nd ed.
- OpenMM main paper: Eastman et al. (2019), *PLOS Computational Biology*
- OPLS-AA reference: Jorgensen et al. (1996), *J. Am. Chem. Soc.*

## 5. How To Cite These in Project Communication

- If you need a reduction-side anchor, start from `R1`.
- If you need a high-voltage instability / EC weakness anchor, start from `R2`.
- If you need an experiment-facing LSV/CV context sentence, use `R3`.
- If you need a formulation example for comparison, use `F1-F4`.
- If you need software/method support, use `S1-S2`.

## 6. Scope Note

This file intentionally collects only links already used in the current repository.
It is not yet a complete bibliography for the field.
