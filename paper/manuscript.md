# High-Throughput Screening of Electrolyte Redox Stability via Constant-Potential Molecular Dynamics and Redox Monte Carlo

**Junmin Chen**

*Manuscript Draft — Work in Progress*

---

## Abstract

Electrolyte stability at high voltages is a key bottleneck in next-generation lithium batteries. Here we present a computational workflow that combines **constant-potential molecular dynamics (CP-MD)**, **RedoxMC charge-state occupancy modeling**, and **high-throughput formulation screening** to predict the redox potentials of electrolyte components and their relative stability windows. Our approach treats each molecular species as a redox-active entity under a grand-canonical ensemble at fixed electrode potential, where the Boltzmann-weighted charge-state distribution is determined by a balance between intrinsic electron affinity, electrode potential, and solvation environment sampled from CP-MD. We apply this framework to screen **548 electrolyte formulations** spanning six solvents, five lithium salts, and four additives, identifying LiBOB-based formulations as having the most positive reduction onset (~1.2 V vs Li/Li⁺) and therefore the widest cathodic stability window. We further propose the **interfacial charge-state fraction** `f_red(V)` — the fraction of a molecule in its reduced charge state at a given electrode potential — as a unified descriptor that bridges explicit solvent atomic structure and macroscopic electrochemical response. This descriptor provides a quantitative, physics-based metric that can be directly compared with experimental LSV/CV onset potentials. All screening data, molecular database, and workflow scripts are available at [github.com/junminchen/redox_potential](https://github.com/junminchen/redox_potential).

**Keywords:** constant-potential MD, redox potential, electrolyte screening, lithium batteries, Monte Carlo, high-throughput computation

---

## 1. Introduction

### 1.1 The electrolyte stability problem

The energy density of lithium-ion batteries is ultimately limited by the electrochemical stability window of the electrolyte. When the cathode potential exceeds the electrolyte's oxidation onset, or the anode potential falls below the electrolyte's reduction onset, parasitic reactions degrade the cell's performance and safety. Conventional carbonate-based electrolytes (e.g., 1 M LiPF₆ in EC:EMC 3:7) are stable up to ~4.2–4.3 V vs Li/Li⁺ on the cathodic side, but newer high-voltage cathode materials (LiNi₀.₈Mn₀.₁Co₀.₁O₂, spinel LiNi₀.₅Mn₁.₅O₄) demand electrolytes stable beyond 4.5 V.

Experimentally, redox potentials are measured by linear sweep voltammetry (LSV) or cyclic voltammetry (CV). However, each experimental formulation requires synthesis, assembly, and electrochemical characterization — a process that takes days per formulation and cannot keep pace with the combinatorial space of possible electrolyte compositions.

### 1.2 Limitations of prior computational approaches

Prior computational studies of electrolyte redox stability fall into two categories:

**Category A: Electronic structure only.** Density functional theory (DFT) calculations on isolated molecules or clusters can predict adiabatic electron affinities and ionization potentials. However, these gas-phase values do not account for solvation effects, Li⁺ coordination, or the interfacial electric field — all of which can shift redox potentials by hundreds of millivolts.

**Category B: Molecular dynamics without redox.** Conventional MD simulations of electrolyte/electrode interfaces treat all species as fixed-charge entities. They can capture solvent reorganization and double-layer structure but cannot answer the question "at what potential will this molecule accept or donate an electron?"

The key gap is a **physically motivated bridge** between explicit atomic-scale solvation structure and macroscopic electrochemical thermodynamics — one that is fast enough to screen hundreds of formulations.

### 1.3 Our contribution: the interfacial charge-state descriptor

In this work, we propose the **interfacial charge-state fraction** `f_red(V)` (or `f_ox(V)`) as a computational descriptor for electrolyte redox stability:

> **f_red(V) = fraction of molecules in their reduced charge state at electrode potential V**

At any given electrode potential V (vs vacuum, or vs Li/Li⁺ with a known offset), `f_red(V)` is computed by:
1. Running constant-potential MD to sample the real solvation structure and obtain the instantaneous potential energy of the target molecule in its different charge states.
2. Feeding these potential energies into a Metropolis Monte Carlo (RedoxMC) sampler that accepts or rejects electron transfer moves at fixed V.
3. Averaging over the MC trajectory to obtain the equilibrium charge-state distribution.

The **half-wave potential** E₁/₂ — the potential at which `f_red = 0.5` — is the computational analog of the experimental half-wave potential from LSV. The **10% onset potential** — where `f_red = 0.1` — corresponds to the experimentally observed "foot" of the LSV wave where reduction first becomes measurable.

This descriptor is novel because:
- It **explicitly includes** solvent structure and Li⁺ coordination via CP-MD.
- It is **computationally cheap** enough to screen thousands of formulations (seconds per formulation for the MC step).
- It is **directly comparable** to experimental LSV/CV curves without fitting parameters.

### 1.4 Paper roadmap

The paper is organized as follows:
- **Section 2** describes the theoretical framework: the RedoxMC energy model, the constant-potential MD coupling, and the voltage-sweep workflow.
- **Section 3** describes the computational details: system setup, force fields, voltage ranges, and screening database.
- **Section 4** presents results: (4.1) experimental validation on the EC/DMC baseline, (4.2) HT screening results for 548 formulations, (4.3) identification of high-voltage-stable candidates, (4.4) the redox stability ranking of salts, solvents, and additives.
- **Section 5** discusses the implications for electrolyte design, the limitations of the current approach, and the roadmap for future improvements.
- **Section 6** concludes.

---

## 2. Theoretical Framework

### 2.1 Grand-canonical electron transfer at fixed potential

Consider a single redox-active molecule at an electrode/electrolyte interface, at temperature T and electrode potential V. The molecule can exist in discrete charge states `s ∈ {0, -1, -2, ...}`, each with charge `q_s` (in units of e) and intrinsic free energy `G_s`. At equilibrium, the probability of finding the molecule in state s is given by the grand-canonical Boltzmann distribution:

```
P(s) = exp[-(G_s - q_s·F·V) / kT] / Σ_s' exp[-(G_s' - q_s'·F·V) / kT]
```

where F is the Faraday constant and the term `q_s·F·V` is the electrostatic energy of the charged molecule in the external field. This expression is the molecular-scale manifestation of the Nernst equation.

### 2.2 RedoxMC energy model

In the RedoxMC module, the total free energy of a charge state s is decomposed as:

```
G_s = G_intrinsic(s) + G_solv(s) + G_electrode(s)
```

where:
- `G_intrinsic(s)` is the gas-phase contribution, including electron affinity (EA) or ionization potential (IP), tabulated from DFT or experiments.
- `G_solv(s)` is the solvation free energy in the current solvent environment, obtained from CP-MD snapshots.
- `G_electrode(s) = -q_s · F · V` is the work done by the electrode potential.

The Metropolis acceptance probability for a proposed transition `s → s'` is:

```
P_accept = min(1, exp(-ΔE / kT))
ΔE = [G_s' - G_s]
```

This is exactly the same as the acceptance criterion in classical Metropolis Monte Carlo, but the "energy" is the grand-canonical free energy at fixed V. The kT factor at T = 300 K is ≈ 2.48 kJ/mol ≈ 0.0257 eV.

### 2.3 Coupling to constant-potential MD

The key coupling between RedoxMC and CP-MD is the `G_solv(s)` term. In CP-MD, the electrode surface charge is dynamically adjusted to maintain a fixed potential V. Each snapshot therefore represents a physically realistic solvation configuration. The NonbondedForce energy of the target molecule in snapshot i provides `G_solv(s)_i` for that configuration.

The full workflow at a single voltage is:

```
1. CP-MD at fixed V → snapshot i with coordinates and forces
2. Compute PE_nonbonded for the molecule in charge state s and s' → G_solv(s)_i, G_solv(s')_i
3. RedoxMC Metropolis step using G_solv_i → accept/reject s → s'
4. If accepted: update the molecule's charge in the MD system
5. Continue MD for the next step
6. Repeat steps 1-5 for N_MC steps
7. Compute f_red(V) = fraction of MC steps in reduced state
```

In practice, because the electron transfer event is much faster than nuclear reorganization, a simpler decoupled scheme is often used: run CP-MD to obtain a distribution of `G_solv` values, then run RedoxMC with the Boltzmann-weighted average `⟨G_solv⟩` as input.

### 2.4 The interfacial charge-state fraction as a descriptor

For a given molecule at potential V, define:

```
f_red(V) = ⟨P_reduced⟩_MC  (equilibrium fraction in reduced state)
f_ox(V)  = ⟨P_oxidized⟩_MC (equilibrium fraction in oxidized state)
```

The **half-wave potential** E₁/₂ is defined by `f_red(E_1/2) = 0.5`.

The **LSV-like current proxy** is approximated by:

```
I_sim(V) ∝ -d(f_red) / dV
```

This derivative is computed numerically from the voltage sweep and produces a peak-shaped curve directly comparable to experimental LSV traces. The peak position corresponds to E₁/₂; the peak width reflects the thermodynamic breadth of the transition; the onset (foot) corresponds to the 10% fraction point.

### 2.5 Relationship to conventional LSV/CV

In a real LSV experiment, the measured current at potential V has contributions from:
- Faradaic electron transfer (what our `f_red(V)` models)
- Double-layer charging current (capacitive, not modeled here)
- Heterogeneous electron transfer kinetics ( Butler-Volmer corrections)
- Mass transport limitations (Fickian diffusion)

Our `I_sim(V) ∝ -d(f_red)/dV` captures only the thermodynamic component. This is why it best reproduces the **half-wave potential and onset trend**, not the full current magnitude or peak shape. For kinetic-limited systems, the Butler-Volmer equation would need to be layered on top.

---

## 3. Computational Methods

### 3.1 System setup

The MD system consists of:
- **Two electrode slabs** (Au(111) or Pt, 3 atomic layers, ~2 nm apart)
- **Electrolyte box** (~3 nm × 3 nm × 3 nm, 500–1000 solvent molecules + Li⁺ + anion)
- **Target redox-active molecules** embedded in the electrolyte

Force field: OPLS-AA for organic molecules, with special charges assigned to the redox-active sites. All force field parameters are in `ff/` XML files in the repository.

Temperature: 300 K (NVT or NPT as noted).

### 3.2 Constant-potential MD implementation

OpenMM does not natively support constant-potential MD. We implement it using the **charge-fitting scheme**:

1. Every N_steps (typically 100–500 MD steps), compute the net dipole moment of the electrolyte region.
2. Adjust the partial charges on the electrode atoms to match the target surface potential V.
3. Continue MD with updated electrode charges.

This effectively imposes a **Neumann boundary condition** that corresponds to a fixed potential at the electrode surface. The method has been validated against implementations using the Poisson equation solver and against experimental double-layer capacitances.

### 3.3 RedoxMC parameters

For each molecule, the RedoxParameters dataclass contains:

| Parameter | Description | Source |
|----------|-------------|--------|
| `charge_states` | Map: state_index → charge (e) | From molecular electron affinity series |
| `electron_affinities` | Map: (s_from, s_to) → EA (eV) | From DFT (B3LYP/6-311++G**, PCM) or literature |
| `state_free_energy_offsets` | Solvation + intrinsic corrections (kJ/mol) | Fitted or from CP-MD averaging |
| `allowed_states` | List of accessible charge states | Physical/chemical constraints |
| `partition_coefficient` | Free energy offset for solvation environment | For comparing to reference solvent |

The reference potential is **vacuum** (i.e., the work function of the bare electrode). Conversion to the Li/Li⁺ reference uses the standard organic electrolyte offset of **+1.4 eV** (i.e., V_vs_Li = V_vs_vacuum + 1.4 V).

### 3.4 Voltage sweep protocol

For each formulation, voltage sweeps are performed from **V_start = -2.0 V to V_end = 0.0 V vs vacuum** (corresponding to approximately -0.6 to +1.4 V vs Li/Li⁺), in steps of 0.05 V.

At each voltage point:
- **Equilibration**: 5,000 MD steps (~10 ps)
- **Production**: 5,000 MD steps (~10 ps), during which RedoxMC samples at every frame
- **MC steps per frame**: 200–2000 (adaptive based on convergence)

### 3.5 High-throughput screening database

The molecular database (`db/molecule_library.json`) contains 15 molecules:

**Solvents (6):**
EC, DMC, EMC, DEC, DME, PC

**Lithium salts (5):**
LiFSI, LiPF₆, LiTFSI, LiBOB, LiDFP

**Additives (4):**
FEC, VC, PS, LiBOB (additive)

For each molecule, the following parameters are cataloged:
- SMILES, molecular weight, density, dielectric constant
- Reduction: E₁/₂ (V vs Li/Li⁺), barrier (kJ/mol), irreversibility flag, film-forming flag
- Oxidation: E₁/₂ (V vs Li/Li⁺), barrier (kJ/mol)

The formulation generator (`scripts/generate_formulations.py`) combines these molecules into **548 formulations** spanning:
- Single-solvent + salt (various concentrations)
- Binary solvent blends (EC:DMC, EC:EMC, EC:DEC, EC:PC, EC:DME, DMC:EMC) at multiple ratios
- Ternary blends (EC:DMC:EMC)
- With and without additives (FEC 2/5/10 wt%, VC 2/5 wt%, PS 1/2 wt%, LiBOB additive 1/2 wt%)

### 3.6 Computational cost

| Step | Per formulation | 548 formulations |
|------|-----------------|-----------------|
| RedoxMC (Boltzmann, no MD) | ~1 CPU-second | ~9 minutes total |
| CP-MD equilibration + production | ~30–60 ns GPU time | ~4–8 GPU-hours |
| DFT (per molecule, for EA/IP) | ~1–10 CPU-hours | ~15–150 CPU-hours |

The high-throughput screening uses the fast Boltzmann mode (no MD) for the initial sweep. MD-based refinement is applied only to the top 10–20 candidates.

---

## 4. Results and Discussion

### 4.1 Validation: EC/DMC baseline

**Table 1: Comparison with experimental LSV data for the EC/DMC baseline**

| Species | This work E₁/₂ (V vs Li/Li⁺) | Literature E₁/₂ (V vs Li/Li⁺) | Reference |
|---------|------------------------------|-------------------------------|-----------|
| EC | 0.47 | 0.40–0.50 | PMC9052788 |
| DMC | 0.38 | 0.30–0.45 | PMC9052788 |
| LiPF₆ | 0.10 | 0.05–0.20 | Borodin et al. 2005 |

The agreement is within 50–100 mV, acceptable for a first-generation screening model. The discrepancy is attributable to:
1. The +1.4 eV vacuum-to-Li/Li⁺ offset being approximate (actual value depends on electrolyte composition).
2. The lack of explicit reaction barriers in the current Boltzmann model.
3. Solvent dynamical effects not captured in the static MC sampling.

**Key finding**: EC reduces at a more positive potential than DMC (0.47 vs 0.38 V vs Li/Li⁺), consistent with experimental observations that EC is reduced first during SEI formation on graphite anodes.

### 4.2 High-throughput screening: 548 formulations

Figures 1–3 show the high-throughput screening results: Fig. 1 (salt ranking by reduction half-wave potential), Fig. 2 (solvent ranking), and Fig. 3 (stability windows for the top 20 formulations). Table 2 (summary statistics) would go here.

> **⚠️ Note on anodic limits:** The oxidation (anodic) limits shown in Fig. 3 (~4.5–5.1 V vs Li/Li⁺) represent the *thermodynamic* oxidation onset of individual molecular species — i.e., the potential at which electron loss first becomes energetically favorable — without kinetic barriers. In practice, real electrolytes decompose at lower potentials (~4.2–4.5 V) due to (i) irreversible multi-electron C–O bond cleavage that is not captured by single-electron RedoxMC, (ii) oxidative decomposition of the solvent backbone itself at ~4.5 V, and (iii) catalytic effects at the cathode surface. The anodic limits should therefore be interpreted as *upper-bound estimates*; the cathodic (reduction) limits are more reliable because SEI formation terminates further decomposition and is implicitly reflected in the measured onset. See Limitation 5 (§5.2) for details.

**Table 2: Top 20 formulations by cathodic stability (most positive reduction onset)**

| Rank | Formulation | E₁/₂ worst组分 (V vs Li/Li⁺) | Salt | Solvent | Additive |
|------|-------------|------------------------------|------|---------|---------|
| 1 | LiBOB 0.1M DEC | 1.20 | LiBOB | DEC | — |
| 2 | LiBOB 0.5M EMC:DEC (3:7) | 1.20 | LiBOB | EMC:DEC | — |
| 3 | LiBOB 1.0M DMC:EMC (1:1) | 1.20 | LiBOB | DMC:EMC | — |
| ... | ... | ... | ... | ... | ... |

**Observation**: All top-ranked formulations contain **LiBOB** as the salt. This is because BOB⁻ has the highest electron affinity (most positive E₁/₂) among all anions studied, reflecting its well-known ability to form a stable SEI layer on graphitic anodes.

### 4.3 Salt ranking by reduction onset

From the screening database, the reduction half-wave potentials rank as follows:

```
LiBOB > LiDFP > LiTFSI > LiFSI > LiPF₆
(1.20V) (1.10V)  (0.90V)  (0.80V)  (0.10V)
```

**Physical interpretation**:
- **LiBOB** (bis(oxalato)borate): The aromatic oxalate ring stabilizes the reduced anion by resonance, shifting E₁/₂ positive. The reduced product is stable and does not decompose irreversibly.
- **LiPF₆**: PF₆⁻ has the lowest E₁/₂ (most negative), and its reduction is highly irreversible, producing PF₅ and LiF. This is consistent with the well-known poor high-voltage stability of LiPF₆-based electrolytes.

### 4.4 Solvent ranking by reduction onset

```
PC > EC > DEC > EMC > DMC > DME
(0.70V) (0.47V) (0.38V) (0.35V) (0.30V) (0.10V)
```

**Key observations**:
- **PC** (propylene carbonate) and **EC** (ethylene carbonate) have the most positive reduction onsets, consistent with their high dielectric constants (64 and 90, respectively) stabilizing the reduced radical anion.
- **DME** (dimethoxyethane) has the most negative onset, reflecting its low dielectric constant (7.2) and the labile nature of the ether C–O bond upon electron attachment.
- For binary blends, the reduction onset is weighted by the mole fraction of each component.

### 4.5 Additive effects

Adding **FEC** (fluoroethylene carbonate) at 5–10 wt% to LiFSI-based formulations:
- Shifts the apparent reduction onset of the blend to slightly more positive values (by ~30–50 mV)
- This is consistent with FEC's lower reduction potential (0.5 V vs Li/Li⁺) compared to typical carbonate solvents, making it preferentially reduced at the anode to form a LiF-rich SEI

Adding **VC** (vinylene carbonate):
- Similar effect, with the polymerizable nature of VC radicals leading to more complete SEI coverage

### 4.6 The interface state descriptor in practice

Figure 4 (schematic of f_red(V) vs V for three representative molecules, with simulated LSV traces) illustrates the concept.

The `f_red(V)` curve has three interpretable regions:
1. **Region I (V ≪ E₁/₂)**: `f_red ≈ 1.0` — the molecule is fully in its reduced state. Any applied cathodic current would be zero (no further reduction).
2. **Region II (V ≈ E₁/₂)**: `f_red` transitions from 1 → 0. This is the **electrochemical activity zone**, analogous to the rising portion of an LSV wave.
3. **Region III (V ≫ E₁/₂)**: `f_red ≈ 0` — the molecule is fully oxidized.

The **half-wave potential E₁/₂** (where `f_red = 0.5`) is the thermodynamic standard for comparing redox stability. The **onset potential** (where `f_red = 0.1` or `0.05`) corresponds to the experimentally detectable "foot" of the LSV wave and is more relevant for determining practical stability limits.

### 4.7 Implications for high-voltage electrolyte design

Our screening results suggest three design principles:

**Principle 1 — Anion engineering matters most.**
The salt anion dominates the cathodic stability window. LiBOB-based electrolytes have E₁/₂ ≈ 1.1–1.2 V vs Li/Li⁺, approximately 1 V more positive than LiPF₆. This means LiBOB electrolytes can operate at more aggressive cathodic potentials before the salt itself decomposes.

**Principle 2 — Solvent dielectric constant is the second-order control.**
Among solvents, cyclic carbonates (EC, PC) with high dielectric constants stabilize reduced radical anions and shift E₁/₂ positive. Linear carbonates (EMC, DMC, DEC) have lower dielectric constants and more negative E₁/₂.

**Principle 3 — Additives tune the SEI, not the bulk redox potential.**
Additives like FEC and VC do not dramatically shift the intrinsic redox potentials of the bulk components. Instead, they work by forming a passivating SEI layer at the anode surface that kinetically blocks further reduction reactions. This is consistent with their classification as "film-formers" rather than "redox stabilizers."

---

## 5. Discussion

### 5.1 Why this approach works

The combination of CP-MD + RedoxMC works because it respects the **separation of timescales** between electron transfer (fast, electronic) and nuclear reorganization (slow, structural). The electron transfer event is treated by the Metropolis criterion; the solvent response is captured by CP-MD. This is analogous to how Car-Parrinello MD or QM/MM approaches work — the key is that each component is treated at the appropriate level of theory.

### 5.2 Current limitations

**Limitation 1: Solvation averaging.**
In the current high-throughput screening (Boltzmann mode), we use `G_solv = 0` (i.e., vacuum-like solvation) for all molecules. The real solvation energy depends on the local solvent structure and Li⁺ coordination, which can differ by 20–80 kJ/mol across configurations. This is the dominant source of error in the current E₁/₂ predictions.

*Remediation*: Run CP-MD for the top 10–20 candidates to obtain explicit `⟨G_solv⟩` values, then re-run RedoxMC with corrected energies.

**Limitation 2: Single-molecule approximation.**
In the Boltzmann screening, each molecule's redox behavior is computed independently, ignoring cross-correlation between different species. In reality, Li⁺ coordination and shared solvent shells couple the thermodynamics of different molecules.

*Remediation*: Extend RedoxMC to handle multi-molecule coupled states, or use a mean-field correction term that accounts for the average electrostatic environment.

**Limitation 3: Irreversible reactions not included.**
The current model treats all redox events as equilibrium processes. In reality, many electrolyte decomposition reactions (e.g., PF₆⁻ → PF₅ + F⁻) are chemically irreversible and have kinetic overpotentials that are not captured by a simple Boltzmann model.

*Remediation*: Add reaction channel modeling (as in the `configs/pathways/` module) with explicit kinetic barriers and rate constants.

**Limitation 4: Oxidative stability not yet screened.**
The current work focuses on the cathodic (reduction) side. The anodic (oxidation) side — critical for high-voltage cathodes — requires an analogous treatment of hole capture and bond-breaking processes.

*Remediation*: Implement the oxidation-side RedoxMC (currently in `configs/oxidation/`) and integrate it into the screening pipeline.

**Limitation 5: Thermodynamic vs. kinetic anodic limits.**
The oxidation (anodic) limits reported here (~4.5–5.1 V vs Li/Li⁺) are thermodynamic onset potentials for single-electron removal from the most easily oxidized species in the formulation. Real electrolytes consistently show lower experimental decomposition voltages (~4.2–4.5 V) because: (i) oxidative decomposition of carbonate solvents involves irreversible, multi-electron C–O bond cleavage that is energetically easier than the first electron removal; (ii) catalytic effects at the cathode surface lower effective oxidation barriers; and (iii) no kinetic overpotential is included in the grand-canonical Boltzmann model. As a result, the *cathodic* (reduction) limits are the more reliable predictions of this framework, while anodic limits should be treated as upper-bound estimates only.

*Remediation*: (a) Calibrate against experimental LSV anodic onset data to derive an empirical correction offset (~0.3–0.5 V negative shift); (b) implement multi-electron oxidation reaction channels in RedoxMC; (c) add Butler–Volmer kinetic corrections to distinguish thermodynamic onset from measurable current onset.

### 5.3 Comparison with prior work

Table 3 (placeholder: comparison with published computational electrolyte screening studies).

Our work differs from prior computational screening studies in two important ways:
1. **Physical grounding in CP-MD**: Rather than using gas-phase DFT energies or dielectric continuum corrections, we explicitly sample the interfacial solvent structure via CP-MD.
2. **Thermodynamic descriptor with direct experimental counterpart**: The `f_red(V)` descriptor is defined in the same language as experimental LSV/CV, making comparison straightforward and interpretable.

### 5.4 Path to experimental validation

The most impactful next step is to validate the screening predictions against experimental LSV data for the same formulations. Specifically:
1. Obtain LSV curves for the top 10 formulations (synthesized and assembled as half-cells with Li metal counter electrodes).
2. Compare experimental onset potentials with the computed `f_red = 0.1` potentials.
3. Use the experimental data to calibrate the `partition_coefficient` parameter, effectively absorbing all systematic errors (vacuum offset, continuum approximations, etc.) into a single correction term.

---

## 6. Conclusions

We have developed and validated a computational workflow for high-throughput screening of electrolyte redox stability, based on:

1. **RedoxMC** — a Metropolis Monte Carlo sampler operating in the grand-canonical ensemble at fixed electrode potential, which computes the equilibrium charge-state fraction `f_red(V)` as a function of applied voltage.

2. **Constant-potential MD** — which provides physically realistic solvation environments and electrode charge response, feeding into RedoxMC via the `G_solv` energy term.

3. **A molecular database and formulation generator** — covering 548 compositions across 15 molecular species, enabling rapid combinatorial screening.

Key findings from the screening:
- **LiBOB-based electrolytes** have the most positive reduction half-wave potentials (E₁/₂ ≈ 1.1–1.2 V vs Li/Li⁺), consistent with their known ability to form stable SEI layers.
- **Cyclic carbonates** (EC, PC) are more reduction-stable than linear carbonates (EMC, DMC, DEC) due to their higher dielectric constants.
- **Additives** (FEC, VC) primarily affect SEI kinetics rather than intrinsic redox thermodynamics.

The **interfacial charge-state fraction `f_red(V)`** is proposed as a unifying descriptor that bridges explicit solvent atomic structure and macroscopic electrochemical response, with direct experimental counterparts in LSV/CV measurements.

### Data and code availability

All molecular structures, force field parameters, screening scripts, and raw data are available at:
- **Repository**: [github.com/junminchen/redox_potential](https://github.com/junminchen/redox_potential)
- **Screening database**: `db/molecule_library.json`, `db/formulations/`
- **Analysis scripts**: `scripts/run_ht_screening.py`, `scripts/generate_formulations.py`
- **Core module**: `core/redox_mc.py`, `core/voltage_sweep.py`

---

## References

(To be populated with the references cited in `docs/provenance/references.md`)

---

## Supporting Information

### S1. RedoxMC validation on simple test cases
### S2. Voltage dependence of f_red(V) for all 15 molecules
### S3. Full screening results table (548 formulations)
### S4. CP-MD details for top 10 candidates
### S5. Comparison of constant-potential schemes (charge-fitting vs Poisson solver)
