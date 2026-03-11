# Redox Parameter Provenance

## Purpose

This file records where the current screening parameters come from and how they should be interpreted.
The project now uses three evidence levels:

- `literature_numeric`: directly anchored to a numeric value reported in literature.
- `literature_trend`: anchored to a qualitative or comparative literature trend.
- `semi_empirical_mapping`: an effective parameter introduced to make the screening model reproduce the intended onset/order/risk trend.

These parameters are suitable for screening and ranking. They are not all first-principles truth values.

For a Chinese, experiment-facing explanation of what is directly supported by literature versus what is still a model inference, see:

- [EXPERIMENTAL_VALIDATION_AND_INFERENCE_CN.md](/home/am3-peichenzhong-group/Documents/project/project_solv_structure/redox_potential/redox_constantV/docs/provenance/EXPERIMENTAL_VALIDATION_AND_INFERENCE_CN.md)

## Source list

1. Atomic thermodynamics and microkinetics of the reduction mechanism of electrolyte additives to facilitate the formation of solid electrolyte interphases in lithium-ion batteries
   - URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC9052788/
   - Used for: EC/DMC reduction-potential targets in Li-ion electrolyte context.

2. Ethylene Carbonate-Free Electrolytes for High Voltage Lithium-Ion Batteries: Progress and Perspectives
   - URL: https://pubs.acs.org/doi/10.1021/acs.energyfuels.4c02728
   - Used for: EC being problematic above about 4.3 V and for qualitative oxidation/ring-opening trend under high voltage.

3. Enhanced high voltage performance of LiNi0.5Mn0.3Co0.2O2 cathode via the synergistic effect of LiPO2F2 and FEC in fluorinated electrolyte for lithium-ion batteries
   - URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC8695088/
   - Used for: experimental LSV/CV context showing EC/DMC-type electrolyte oxidation response in the 3.0-6.0 V range.

## Interpretation rules

- `reduction_onset` targets for EC and DMC are the strongest anchors in the current model.
- `oxidation_onset` targets are formulation-level effective screening targets, not universal constants.
- `reaction pathway` onsets and barriers are currently semi-empirical and encode literature-consistent failure order and mechanism labels.
- Barrier values and slopes should be replaced by DFT/NEB data when available.
