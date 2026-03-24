#!/usr/bin/env python3
"""Generate figures for redox_potential manuscript."""

import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────
DATA = Path(__file__).parent.parent.parent / "results" / "ht_screening_v2"
OUT  = Path(__file__).parent
OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})

# ── colour palettes ───────────────────────────────────────────────────────────
SALT_COLORS = {
    "LiBOB":  "#E63946",
    "LiDFP":  "#F4A261",
    "LiTFSI": "#2A9D8F",
    "LiFSI":  "#457B9D",
    "LiPF6":  "#1D3557",
}
SOLVENT_COLORS = {
    "PC":  "#E63946",
    "EC":  "#F4A261",
    "DEC": "#2A9D8F",
    "EMC": "#457B9D",
    "DMC": "#8AB17D",
    "DME": "#B2857A",
}


def load_csv(name):
    p = DATA / name
    with open(p) as f:
        return list(csv.DictReader(f))


# ══════════════════════════════════════════════════════════════════════════════
# FIG 1 – Salt ranking by reduction half-wave potential (cathodic stability)
# ══════════════════════════════════════════════════════════════════════════════
def fig1_salt_ranking():
    rows = load_csv("screening_summary.csv")

    # collect salt-level data  (half-wave, onset)
    salt_hw  = {}   # salt → half-wave
    salt_onset = {} # salt → onset
    for r in rows:
        s = r["salt"].strip()
        if s not in salt_hw:
            salt_hw[s]   = []
            salt_onset[s] = []
        col = f"{s}_red_halfwave_V_vs_Li"
        if r.get(col):
            salt_hw[s].append(float(r[col]))
        col_o = f"{s}_red_onset_V_vs_Li"
        if r.get(col_o):
            salt_onset[s].append(float(r[col_o]))

    salts = ["LiBOB", "LiDFP", "LiTFSI", "LiFSI", "LiPF6"]
    # deduplicate lists (same value repeated across formulations)
    hw_mean  = [np.mean(salt_hw[s])   for s in salts]
    on_mean  = [np.mean(salt_onset[s]) for s in salts]

    x = np.arange(len(salts))
    width = 0.4

    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.bar(x - width/2, hw_mean,  width, label="E½ (half-wave)", color=[SALT_COLORS[s] for s in salts], alpha=0.9, edgecolor="white")
    ax.bar(x + width/2, on_mean,  width, label="E onset (10%)",  color=[SALT_COLORS[s] for s in salts], alpha=0.45, edgecolor="white")

    for i, s in enumerate(salts):
        ax.text(i, hw_mean[i] + 0.02, f"{hw_mean[i]:.2f}V", ha="center", va="bottom", fontsize=8.5)

    ax.set_xticks(x)
    ax.set_xticklabels(salts)
    ax.set_ylabel("Potential (V vs Li/Li⁺)")
    ax.set_ylim(0, 1.5)
    ax.set_title("Salt Ranking by Cathodic Redox Stability", fontweight="bold", pad=10)
    ax.legend(frameon=False, fontsize=9)
    ax.set_xlabel("Lithium Salt (1 M in EC)")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_salt_ranking.png", bbox_inches="tight")
    plt.close()
    print("Saved fig1_salt_ranking.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIG 2 – Solvent ranking by reduction onset
# ══════════════════════════════════════════════════════════════════════════════
def fig2_solvent_ranking():
    rows = load_csv("screening_summary.csv")

    solvents = ["PC", "EC", "DEC", "EMC", "DMC", "DME"]
    onset_mean  = []
    hw_mean     = []
    for sv in solvents:
        o_vals, h_vals = [], []
        for r in rows:
            col_o = f"{sv}_red_onset_V_vs_Li"
            col_h = f"{sv}_red_halfwave_V_vs_Li"
            if r.get(col_o) and r[col_o].strip():
                o_vals.append(float(r[col_o]))
            if r.get(col_h) and r[col_h].strip():
                h_vals.append(float(r[col_h]))
        onset_mean.append(np.mean(o_vals) if o_vals else np.nan)
        hw_mean.append(np.mean(h_vals) if h_vals else np.nan)

    x = np.arange(len(solvents))
    width = 0.38
    fig, ax = plt.subplots(figsize=(7, 4))
    bars1 = ax.bar(x - width/2, hw_mean,    width, label="E½ (half-wave)",  color=[SOLVENT_COLORS[s] for s in solvents], alpha=0.9,  edgecolor="white")
    bars2 = ax.bar(x + width/2, onset_mean,  width, label="E onset (10%)",   color=[SOLVENT_COLORS[s] for s in solvents], alpha=0.45, edgecolor="white")

    for i, (h, o) in enumerate(zip(hw_mean, onset_mean)):
        if not np.isnan(h):
            ax.text(i, h + 0.02, f"{h:.2f}", ha="center", va="bottom", fontsize=8)
        if not np.isnan(o):
            ax.text(i, o + 0.02, f"{o:.2f}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(solvents)
    ax.set_ylabel("Potential (V vs Li/Li⁺)")
    ax.set_ylim(0, 1.0)
    ax.set_title("Solvent Ranking by Reduction Stability", fontweight="bold", pad=10)
    ax.legend(frameon=False, fontsize=9)
    ax.set_xlabel("Solvent")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_solvent_ranking.png", bbox_inches="tight")
    plt.close()
    print("Saved fig2_solvent_ranking.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIG 3 – Stability window for top 20 formulations (grouped by salt)
# ══════════════════════════════════════════════════════════════════════════════
def fig3_stability_windows():
    rows = load_csv("screening_summary.csv")

    # sort by window width
    def window(row):
        w = row.get("stability_window_V", "")
        return float(w) if w.strip() else 0

    rows_sorted = sorted(rows, key=window, reverse=True)[:20]

    labels, cath, anodic, widths = [], [], [], []
    for r in rows_sorted:
        lbl = r["formulation_name"].replace("LiBOB_", "LiBOB·").replace("_", "·")
        labels.append(lbl)
        cath_val  = float(r["cathodic_limit_V_vs_Li"])  if r.get("cathodic_limit_V_vs_Li","").strip()  else np.nan
        an_val    = float(r["anodic_limit_V_vs_Li"])    if r.get("anodic_limit_V_vs_Li","").strip()    else np.nan
        cath.append(cath_val)
        anodic.append(an_val)
        widths.append(an_val - cath_val if not (np.isnan(cath_val) or np.isnan(an_val)) else np.nan)

    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.5, 7))

    salt_color = {s: c for s, c in SALT_COLORS.items()}
    for i, (lb, ca, an, wid, row) in enumerate(zip(labels, cath, anodic, widths, rows_sorted)):
        s = row["salt"].strip()
        col = salt_color.get(s, "#888")
        ax.barh(i, wid, left=ca, height=0.6, color=col, alpha=0.8, edgecolor="white")
        ax.text(ca - 0.02, i, f"{ca:.2f}", va="center", ha="right", fontsize=7.5, color=col)
        ax.text(an + 0.02, i, f"{an:.2f}", va="center", ha="left",  fontsize=7.5, color=col)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("Potential (V vs Li/Li⁺)")
    ax.set_title("Stability Windows — Top 20 Formulations", fontweight="bold", pad=10)

    patches = [mpatches.Patch(color=SALT_COLORS[s], label=s, alpha=0.8) for s in ["LiBOB","LiDFP","LiTFSI","LiFSI","LiPF6"]]
    ax.legend(handles=patches, frameon=False, fontsize=8, loc="lower right", title="Salt", title_fontsize=8)
    ax.set_xlim(-0.1, 5.5)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_stability_windows.png", bbox_inches="tight")
    plt.close()
    print("Saved fig3_stability_windows.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIG 4 – f_red(V) from real RedoxMC calculations (EC and DMC)
# Data source: experimental_redox_report / electrolyte_redox_curves.csv
# ══════════════════════════════════════════════════════════════════════════════
def fig4_fred_concept():
    REDOX_CSV = (
        Path(__file__).parent.parent.parent
        / "results" / "experimental_redox_report_20260311_234115"
        / "electrolyte_redox_curves.csv"
    )

    import csv as _csv
    data = {"EC": [], "DMC": []}
    with open(REDOX_CSV) as f:
        for row in _csv.DictReader(f):
            mol = row["molecule"].strip()
            if mol in data:
                data[mol].append({
                    "v_li":    float(row["voltage_li_v"]),
                    "f_red":   float(row["fraction_reduced"]),
                    "mean_q":  float(row["mean_charge_e"]),
                })

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    COL = {"EC": "#F4A261", "DMC": "#8AB17D"}

    ax = axes[0]
    for mol, col in COL.items():
        rows = sorted(data[mol], key=lambda r: r["v_li"])
        v  = [r["v_li"]  for r in rows]
        fr = [r["f_red"] for r in rows]
        ax.plot(v, fr, "o-", label=mol, color=col, lw=1.8, ms=4)

    ax.axhline(0.5,  color="grey", lw=0.8, ls="--", alpha=0.6)
    ax.axhline(0.1,  color="grey", lw=0.8, ls=":",  alpha=0.6)
    ax.text(ax.get_xlim()[1] - 0.02, 0.5,  "f_red = 0.5 (E½)",  va="center", ha="right", fontsize=7.5, color="grey")
    ax.text(ax.get_xlim()[1] - 0.02, 0.1,  "f_red = 0.1 (onset)", va="center", ha="right", fontsize=7.5, color="grey")
    ax.set_xlabel("Electrode potential V (V vs Li/Li⁺)")
    ax.set_ylabel("f_red(V)")
    ax.set_title("Charge-State Fraction  f_red(V)", fontweight="bold")
    ax.legend(frameon=False)
    ax.set_ylim(-0.05, 1.05)

    # Compute I_sim = -d(f_red)/dV numerically
    ax2 = axes[1]
    for mol, col in COL.items():
        rows = sorted(data[mol], key=lambda r: r["v_li"])
        v  = np.array([r["v_li"]  for r in rows])
        fr = np.array([r["f_red"] for r in rows])
        # simple first-order derivative
        dv = np.gradient(v)
        I  = -np.gradient(fr, v)   # -df_red/dV
        ax2.plot(v, I, label=mol, color=col, lw=1.8)

    ax2.set_xlabel("Electrode potential V (V vs Li/Li⁺)")
    ax2.set_ylabel("I_sim(V) ∝ −df_red/dV")
    ax2.set_title("Simulated LSV Trace (−df_red/dV)", fontweight="bold")
    ax2.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(OUT / "fig4_fred_concept.png", bbox_inches="tight")
    plt.close()
    print("Saved fig4_fred_concept.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIG 5 – Additive effect: FEC/VC on LiFSI formulations
# ══════════════════════════════════════════════════════════════════════════════
def fig5_additive_effect():
    rows = load_csv("screening_summary.csv")

    # Filter: LiFSI salt only
    lifsi = [r for r in rows if r.get("salt","").strip() == "LiFSI"]

    no_add = [r for r in lifsi if not r.get("additive","").strip() or r["additive"].strip() == "0.0"]
    fec5   = [r for r in lifsi if r.get("additive","").strip() == "FEC"]
    vc     = [r for r in lifsi if r.get("additive","").strip() == "VC"]

    def avg_window(rs):
        vals = [float(r["stability_window_V"]) for r in rs if r.get("stability_window_V","").strip()]
        return np.mean(vals) if vals else np.nan

    categories = ["No additive", "FEC (2–10 wt%)", "VC (2–5 wt%)"]
    means = [avg_window(no_add), avg_window(fec5), avg_window(vc)]
    stds  = [
        np.std([float(r["stability_window_V"]) for r in no_add if r.get("stability_window_V","").strip()]) if no_add else 0,
        np.std([float(r["stability_window_V"]) for r in fec5   if r.get("stability_window_V","").strip()]) if fec5 else 0,
        np.std([float(r["stability_window_V"]) for r in vc     if r.get("stability_window_V","").strip()]) if vc else 0,
    ]

    x = np.arange(len(categories))
    bars = ax = plt.figure(figsize=(5, 4)).add_subplot(111)
    bars = ax.bar(x, means, yerr=stds, color=["#457B9D","#E63946","#2A9D8F"],
                  alpha=0.8, edgecolor="white", capsize=5)

    for i, (m, s) in enumerate(zip(means, stds)):
        if not np.isnan(m):
            ax.text(i, m + s + 0.05, f"{m:.2f}V", ha="center", va="bottom", fontsize=10)

    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.set_ylabel("Stability Window (V)")
    ax.set_title("Effect of SEI-Forming Additives\n(LiFSI salt, EC-based solvent)", fontweight="bold")
    ax.set_ylim(0, 5.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(OUT / "fig5_additive_effect.png", bbox_inches="tight")
    plt.close()
    print("Saved fig5_additive_effect.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIG 6 – Salt concentration effect on stability window (LiFSI in EC)
# ══════════════════════════════════════════════════════════════════════════════
def fig6_concentration_effect():
    rows = load_csv("screening_summary.csv")

    # Filter: LiFSI salt, EC solvent, no additive
    target = [
        r for r in rows
        if r.get("salt","").strip() == "LiFSI"
        and "EC" in r.get("composition","")
        and (not r.get("additive","") or r["additive"].strip() in ("","0.0"))
        and "PC" not in r.get("composition","")
    ]

    conc   = sorted(set(float(r["salt_conc_M"]) for r in target if r.get("salt_conc_M","").strip()))
    hw_ec  = []
    on_ec  = []
    for c in conc:
        matching = [r for r in target if float(r["salt_conc_M"]) == c]
        hw = [float(r["EC_red_halfwave_V_vs_Li"]) for r in matching if r.get("EC_red_halfwave_V_vs_Li","").strip()]
        on = [float(r["EC_red_onset_V_vs_Li"])   for r in matching if r.get("EC_red_onset_V_vs_Li","").strip()]
        hw_ec.append(np.mean(hw) if hw else np.nan)
        on_ec.append(np.mean(on) if on else np.nan)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(conc, hw_ec, "o-", color="#457B9D", lw=2, ms=6, label="E½ (EC)")
    ax.plot(conc, on_ec, "s--", color="#F4A261", lw=1.5, ms=5, label="E onset (EC)")
    ax.set_xlabel("LiFSI concentration (M)")
    ax.set_ylabel("EC reduction potential (V vs Li/Li⁺)")
    ax.set_title("Salt Concentration Effect on EC Reduction", fontweight="bold")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(OUT / "fig6_concentration_effect.png", bbox_inches="tight")
    plt.close()
    print("Saved fig6_concentration_effect.png")


if __name__ == "__main__":
    fig1_salt_ranking()
    fig2_solvent_ranking()
    fig3_stability_windows()
    fig4_fred_concept()
    fig5_additive_effect()
    fig6_concentration_effect()
    print("\nAll figures saved to", OUT)
