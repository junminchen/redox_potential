# Example Formulations Notes

These examples are screening presets constrained by the current local force-field coverage.

Local force-field / residue coverage confirmed in this repo:
- Solvents: `EC`, `DMC`, `DEC`, `EMC`, `FEC`, `DME`, `DOL`
- Salts / anions: `LiA`, `PF6`, `FSI`, `BF4`, `BOB`, `DFO`, `TFS`

Examples added here were chosen because they are both representative in the literature and at least approximately mappable to the local residue library.

Important limitations:
- These configs are formulation-specific semi-empirical screening presets.
- Concentration effects are not inferred automatically from MD yet; they are encoded via effective redox/pathway parameters.
- Therefore, use them to compare trends and identify promising windows, not as final quantitative electrochemical truth.

Literature anchors used for the added examples:
- `4 m LiFSI in DME`
  - Scientific Reports: https://www.nature.com/articles/s41598-017-16268-7
  - OSTI summary: https://www.osti.gov/biblio/1364007
- `LiFSI:EC 1:4`
  - ACS Applied Energy Materials / PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC8790720/
- `1:1.1 LiFSI in DMC` (LiFSA/DMC in the paper; LiFSA and LiFSI refer to the same FSI anion family naming convention)
  - Nature Communications / PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC4931331/
- `1 M LiPF6 in EC:EMC (3:7)` as a conventional baseline
  - Example usage in a recent battery study: https://pmc.ncbi.nlm.nih.gov/articles/PMC11428382/
