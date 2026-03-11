# Electrolyte Redox Screening Report

## What This Figure Means

- Left panel: an LSV-like reduction-response proxy based on the voltage derivative of reduced-state occupancy.
- Right panel: the equilibrium reduced fraction, which is the clearest way to read reduction onset and ordering.
- Potentials are reported in the experimental convention `V vs Li/Li+`.

## Main Conclusions

- `EC`: onset(10%) `0.527 V`, E1/2 `0.470 V`, LSV-like peak `0.470 V vs Li/Li+`.
- `DMC`: onset(10%) `0.437 V`, E1/2 `0.380 V`, LSV-like peak `0.380 V vs Li/Li+`.
- Easier to reduce ranking in this model: `EC > DMC`.

## Method

1. Read the redox configuration for each molecule.
2. Convert the configured redox free-energy model into Boltzmann occupancies over a voltage grid.
3. Compute `fraction_reduced(V)` and use `-d(fraction_reduced)/dV` as an LSV-like proxy.
4. Report onset (10%), half-wave (50%), and near-complete reduction (90%) potentials.

## Principle

- The model is a constant-potential grand-canonical occupancy model.
- A more negative potential stabilizes reduced states.
- The crossing where oxidized and reduced populations are equal is the model half-wave potential.
- This is best interpreted as a thermodynamic reduction-onset screen, not a direct current prediction.

## Why Experimentalists Can Use It

- It outputs `V vs Li/Li+`, onset, and LSV-like peak positions directly.
- It makes molecular ranking obvious before expensive electrochemistry or full DFT workflows.
- It is well suited for formulation screening and sensitivity analysis.

## Cost and Throughput

- The equilibrium screening mode used here is cheap: seconds to minutes per formulation.
- The MD-coupled mode is much more expensive: minutes to hours per full sweep on one GPU.
- This supports rapid pre-screening followed by expensive refinement only for short-listed candidates.

## Accuracy and Limits

- Expected accuracy in the current effective model is qualitative to semi-quantitative.
- Ranking and onset order can be informative, but peak currents and exact LSV shapes are not predicted.
- Precision depends on how well the redox free-energy offsets reflect real solvation, ion pairing, and decomposition chemistry.
- EC/DMC reduction is often irreversible and tied to SEI-forming chemistry, so these outputs should be read as onset tendencies.

## Can It Be Used for Screening?

- Yes, for molecular ordering, identifying likely first-to-reduce species, and narrowing the voltage window of interest.
- Not yet as a direct replacement for experimental LSV/CV current curves.

## Next Optimizations

1. Couple RedoxMC to frame-resolved solvation energies from MD instead of the current effective offsets.
2. Add oxidation branches so high-voltage stability can be screened symmetrically.
3. Fit offsets against a training set of known experimental onsets.
4. Include reaction channels for irreversible solvent decomposition and SEI formation.
5. Build a formulation parser that maps composition -> molecules -> ranked redox windows.

## Future Product Direction

- Input: an electrolyte recipe (solvents, salt, additives, concentrations).
- Output: ranked molecular reduction/oxidation potentials, likely first-reduced species, and a high-voltage stability flag.
- For a target such as `5 V`, the tool would compare the oxidation side of each component against the applied window and report likely unstable species.
- The long-term value is fast formulation triage before synthesis and electrochemical testing.
