# High-Voltage Stability Workflow

## Goal

Build a workflow that can take an electrolyte formulation and answer:
- which molecules are easiest to reduce
- which molecules are easiest to oxidize
- whether a target high-voltage condition such as `5.0 V` is likely to be stable
- the ranking of redox windows for all solvents, salts, and additives

## Practical Interpretation of "5 V Stable"

For a lithium-battery electrolyte, `5 V stable` does not mean a molecule has no electron-transfer channel at all.
It means that under the target electrode potentials and interfacial environment, oxidation and reduction onsets both remain outside the operating window, and decomposition is slow enough not to dominate the electrochemistry.

## Proposed Pipeline

1. Parse formulation
- Input solvent mixture, salt, additives, concentration, and electrode material.
- Expand to a molecular list and concentrations.

2. Build molecular redox descriptors
- Reduction-side descriptors: electron affinity, anion stabilization, Li+ coordination effects.
- Oxidation-side descriptors: ionization potential / hole-removal free energy, cation-radical stabilization.
- Reaction descriptors: bond cleavage pathways after electron or hole transfer.

3. Generate effective free-energy offsets
- Fit or infer reduced and oxidized state free energies in electrolyte.
- Include salt and coordination corrections.

4. Screen equilibrium redox windows
- Predict onset, E1/2, and ranking for each molecule.
- Flag molecules whose oxidation or reduction enters the device operating window.

5. Refine with interfacial constant-potential MD
- Sample specific shortlisted molecules near the target electrode.
- Include interfacial solvation and electric field effects.

6. Escalate to reaction-path calculations when needed
- For unstable species, compute decomposition barriers and likely products.
- Distinguish reversible redox from irreversible breakdown / SEI formation.

## Output for Experimental Users

- Reduction onset ranking
- Oxidation onset ranking
- Estimated safe operating window
- A `traffic-light` flag for each formulation component:
  - green: outside window
  - yellow: near threshold
  - red: likely unstable
- A concise explanation of the limiting molecule

## What Is Needed to Predict 5 V Stability Well

The current repository can screen the reduction side with effective models.
To support 5 V high-voltage stability, the next essential additions are:
- oxidation-state modeling, not only reduction-state modeling
- formulation-aware coordination corrections
- decomposition pathways for oxidized solvents and anions
- calibration against experimental high-voltage LSV/CV data

## Recommended Development Stages

### Stage 1: Fast screening
- Use effective offsets for both reduction and oxidation.
- Rank molecules and estimate operating window.
- Best for formulation triage.

### Stage 2: MD-informed screening
- Add interfacial solvation sampling from constant-potential MD.
- Improve ranking near the stability threshold.

### Stage 3: Reaction-aware prediction
- Add decomposition barriers and product channels.
- Support claims about real stability at high voltage.

## Expected Accuracy by Stage

- Stage 1: qualitative to semi-quantitative ranking
- Stage 2: semi-quantitative onset prediction
- Stage 3: best route to predictive formulation stability claims
