#!/usr/bin/env python3
"""
Reproduce the internalizing-arm tumor-selectivity figure and per-cancer
statistics from the included extracted per-sample table (no 8 GB download).

Compares candidate internalizing-arm targets (CDCP1, MET, EphA2) against the
fast-internalizing but broadly-expressed receptors conventionally used as
eTPD/uptake effector arms (TfR1/TFRC, CXCR7/ACKR3, CI-M6PR/IGF2R, ASGPR/ASGR1,
LDLR, sortilin/SORT1, LRP1), with HER2 (ERBB2) and EGFR as validated ADC-target
controls.

Input : data/extracted/InternalizingArm_RAW_persample_log2TPM.csv
Output: results/figures/InternalizingArm_tumor_upregulation_dotsbar.(png|svg)
        data/derived/InternalizingArm_perCancer_delta.csv, _stats_summary.csv

Run:   python scripts/make_figures.py        (from repo root)   or:  make
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats

NM = "InternalizingArm"
RNG = np.random.default_rng(0)
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW  = os.path.join(HERE, "data", "extracted", f"{NM}_RAW_persample_log2TPM.csv")
FIGDIR = os.path.join(HERE, "results", "figures")
DERIVED = os.path.join(HERE, "data", "derived")

ARMS = ["CDCP1", "MET", "EPHA2"]
CTRL = ["ERBB2", "EGFR"]
EFF  = ["TFRC", "ACKR3", "IGF2R", "LDLR", "SORT1", "LRP1", "ASGR1"]
ORDER = ARMS + CTRL + EFF
ALIAS = {"TFRC": "TFRC (TfR1)", "ACKR3": "ACKR3 (CXCR7)", "IGF2R": "IGF2R (CI-M6PR)",
         "ASGR1": "ASGR1 (ASGPR)", "SORT1": "SORT1 (sortilin)", "EPHA2": "EphA2",
         "ERBB2": "ERBB2 (HER2)"}
GTEX_MATCH = {"ACC":"Adrenal Gland","BLCA":"Bladder","BRCA":"Breast","CESC":"Cervix Uteri","CHOL":"Liver",
 "COAD":"Colon","DLBC":"Blood","ESCA":"Esophagus","GBM":"Brain","KICH":"Kidney","KIRC":"Kidney","KIRP":"Kidney",
 "LAML":"Bone Marrow","LGG":"Brain","LIHC":"Liver","LUAD":"Lung","LUSC":"Lung","OV":"Ovary","PAAD":"Pancreas",
 "PCPG":"Adrenal Gland","PRAD":"Prostate","READ":"Colon","SKCM":"Skin","STAD":"Stomach","TGCT":"Testis",
 "THCA":"Thyroid","THYM":"Blood","UCEC":"Uterus","UCS":"Uterus"}

def stars(q): return "****" if q<1e-4 else "***" if q<1e-3 else "**" if q<1e-2 else "*" if q<0.05 else "ns"
def bh(p):
    m=len(p); sp=p.sort_values(); q=np.minimum.accumulate((sp.values*m/np.arange(m,0,-1))[::-1])[::-1]
    return pd.Series(np.clip(q,0,1),index=sp.index)

def main():
    df = pd.read_csv(RAW)
    cancers = sorted(GTEX_MATCH)
    delta = pd.DataFrame(index=cancers, columns=ORDER, dtype=float)
    for c in cancers:
        tum = df[(df.group=="Tumor_TCGA") & (df.cancer_abbr==c)]
        nor = df[((df.group=="Normal_TCGA") & (df.cancer_abbr==c)) |
                 ((df.group=="Normal_GTEx") & (df.tissue==GTEX_MATCH[c]))]
        if len(tum) and len(nor):
            for g in ORDER: delta.loc[c,g] = tum[g].median() - nor[g].median()
    delta.index.name = "Cancer"
    os.makedirs(DERIVED, exist_ok=True)
    delta.round(3).to_csv(os.path.join(DERIVED, f"{NM}_perCancer_delta.csv"))
    meanD = {g: delta[g].dropna().mean() for g in ORDER}
    q = bh(pd.Series({g: stats.wilcoxon(delta[g].dropna())[1] for g in ORDER}))
    pd.DataFrame({"gene": ORDER,
                  "class": ["arm" if g in ARMS else "control" if g in CTRL else "effector" for g in ORDER],
                  "mean_delta": [round(meanD[g],3) for g in ORDER],
                  "n_up": [int((delta[g].dropna()>0).sum()) for g in ORDER],
                  "wilcoxon_q_BH": [f"{q[g]:.2e}" for g in ORDER],
                  "sig": [stars(q[g]) for g in ORDER]}
                 ).sort_values("mean_delta", ascending=False
                 ).to_csv(os.path.join(DERIVED, f"{NM}_stats_summary.csv"), index=False)

    colmap = {**{g:"#c0392b" for g in ARMS}, **{g:"#e8a33d" for g in CTRL}, **{g:"#5b6b7a" for g in EFF}}
    gs = sorted(ORDER, key=lambda g: meanD[g])
    fig, ax = plt.subplots(figsize=(10.5,8.4))
    for i,g in enumerate(gs):
        d = delta[g].dropna().values
        ax.barh(i, meanD[g], color=colmap[g], alpha=.85, height=.62, zorder=2, edgecolor="black", linewidth=0.5)
        ax.scatter(d, RNG.normal(i,0.07,len(d)), s=13, color="black", alpha=.45, zorder=3, linewidths=0)
        s = stars(q[g])
        ax.text(max(meanD[g], d.max())+0.15, i, s, va="center", fontsize=10,
                fontweight="bold" if s!="ns" else "normal", color="#111" if s!="ns" else "#999")
    ax.axvline(0, color="k", lw=.8)
    ax.set_yticks(range(len(gs))); ax.set_yticklabels([ALIAS.get(g,g) for g in gs], fontsize=11)
    for t,g in zip(ax.get_yticklabels(), gs):
        t.set_color(colmap[g] if g in ARMS+CTRL else "#333"); t.set_fontweight("bold" if g in ARMS+CTRL else "normal")
    ax.set_xlabel("Tumor - Normal   delta median log2(TPM+1)   (dot = one cancer type, bar = mean)", fontsize=11)
    ax.set_title("Tumor selectivity: candidate internalizing arms vs conventional eTPD/uptake receptors",
                 fontsize=12, fontweight="bold")
    ax.legend(handles=[Patch(facecolor="#c0392b", label="Candidate internalizing arm (CDCP1, MET, EphA2)"),
                       Patch(facecolor="#e8a33d", label="Validated ADC-target control (HER2, EGFR)"),
                       Patch(facecolor="#5b6b7a", label="Conventional eTPD / uptake receptor")],
              frameon=False, fontsize=9.5, loc="lower right")
    ax.text(0.5, -0.10, "One-sample Wilcoxon of per-cancer delta vs 0, BH-FDR.  "
            "**** q<1e-4  *** q<1e-3  ** q<1e-2  * q<0.05",
            transform=ax.transAxes, ha="center", fontsize=8, color="#666")
    ax.spines[["top","right"]].set_visible(False); ax.margins(y=0.01)
    plt.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(os.path.join(FIGDIR, f"{NM}_tumor_upregulation_dotsbar.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIGDIR, f"{NM}_tumor_upregulation_dotsbar.svg"), bbox_inches="tight")
    print("Wrote figure + derived tables.")

if __name__ == "__main__":
    main()
