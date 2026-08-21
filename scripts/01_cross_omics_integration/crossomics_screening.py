from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# ============================================================
# Paths
# ============================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data" / "crossomics"
RESULT_DIR = PROJECT_ROOT / "results" / "crossomics"

RESULT_DIR.mkdir(parents=True, exist_ok=True)

workdir = DATA_DIR
outdir = RESULT_DIR
outdir.mkdir(parents=True, exist_ok=True)

scrna_file = workdir / "scRNA_Fibroblast_pseudobulk_DEG_Ischemic_vs_Control.csv"
st_file = workdir / "ST_section_level_gene_screening_Pathological_vs_Control.csv"
cross_file = workdir / "Crossomics_scRNA_fibroblast_vs_ST_screening_merged.csv"

for f in [scrna_file, st_file, cross_file]:
    if not f.exists():
        raise FileNotFoundError(f)

# ============================================================
# Style
# ============================================================
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 0.7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight"
})

COL_BG = "#D7DCE2"          # background genes
COL_UP = "#D98C3F"          # concordant/upregulated
COL_UP_LIGHT = "#E7B77A"
COL_DOWN = "#6BAED6"
COL_FSTL1 = "#D62728"
COL_KEY = "#111111"
COL_AXIS = "#333333"
COL_LINE = "#9A9A9A"

KEY_GENES = [
    "FSTL1",
    "POSTN", "COL1A1", "COL3A1", "COL1A2", "FN1",
    "SPARC", "VCAN", "FBN1", "MFAP5", "COL6A1",
    "TIMP1", "COL5A1"
]

# Fewer labels for crowded volcanoes
SCRNA_LABELS = ["FSTL1", "POSTN", "FN1", "COL1A1", "COL3A1", "SPARC", "VCAN"]
ST_LABELS = ["FSTL1", "COL1A1", "COL3A1", "COL1A2", "FN1", "SPARC", "COL6A1", "TIMP1"]
CROSS_LABELS = ["FSTL1", "POSTN", "COL1A1", "COL3A1", "COL1A2", "FN1", "SPARC", "VCAN", "COL6A1", "TIMP1", "MFAP5"]

# Manual label offsets in points: dx, dy
SCRNA_OFFSETS = {
    "FSTL1": (8, 4),
    "POSTN": (8, 6),
    "FN1": (8, 6),
    "COL1A1": (8, 2),
    "COL3A1": (8, 4),
    "SPARC": (8, -8),
    "VCAN": (8, 4)
}

ST_OFFSETS = {
    "FSTL1": (8, 4),
    "COL1A1": (8, 4),
    "COL3A1": (8, -6),
    "COL1A2": (8, 4),
    "FN1": (8, 4),
    "SPARC": (8, 4),
    "COL6A1": (8, 4),
    "TIMP1": (8, -8)
}

CROSS_OFFSETS = {
    "FSTL1": (8, -8),
    "POSTN": (8, 2),
    "COL1A1": (8, 4),
    "COL3A1": (8, 4),
    "COL1A2": (8, -8),
    "FN1": (8, 2),
    "SPARC": (8, -8),
    "VCAN": (8, 2),
    "COL6A1": (8, 2),
    "TIMP1": (8, 4),
    "MFAP5": (8, -8)
}

def save_all(fig, prefix):
    for ext in ["pdf", "png", "tiff"]:
        dpi = 600 if ext in ["png", "tiff"] else 300
        fig.savefig(outdir / f"{prefix}.{ext}", dpi=dpi)
    plt.close(fig)

def format_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(COL_AXIS)
        ax.spines[s].set_linewidth(0.7)
    ax.tick_params(axis="both", width=0.7, length=3, color=COL_AXIS)
    ax.grid(False)

def add_gene_label(ax, x, y, gene, offset=(6, 4), color=COL_KEY, size=6.8):
    ax.scatter(
        [x], [y],
        s=32 if gene == "FSTL1" else 22,
        color=COL_FSTL1 if gene == "FSTL1" else COL_KEY,
        edgecolor="white",
        linewidth=0.35,
        zorder=10
    )
    ax.annotate(
        gene,
        xy=(x, y),
        xytext=offset,
        textcoords="offset points",
        fontsize=size,
        color=COL_FSTL1 if gene == "FSTL1" else color,
        ha="left",
        va="center",
        arrowprops=dict(
            arrowstyle="-",
            color=COL_FSTL1 if gene == "FSTL1" else "#555555",
            lw=0.35,
            shrinkA=0,
            shrinkB=2
        ),
        zorder=11
    )

def add_panel_label(ax, label):
    ax.text(
        -0.18, 1.06, label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        ha="left",
        va="top"
    )

# ============================================================
# Load data
# ============================================================
scrna = pd.read_csv(scrna_file)
st = pd.read_csv(st_file)
cross = pd.read_csv(cross_file)

# Standard columns
scrna["gene"] = scrna["gene"].astype(str)
st["gene"] = st["gene"].astype(str)
cross["gene"] = cross["gene"].astype(str)

# Use nominal P for visualization; retain FDR in table/legend
scrna["neglog10_p"] = -np.log10(scrna["p_value"].clip(lower=1e-300))
scrna["neglog10_fdr"] = -np.log10(scrna["p_adj"].clip(lower=1e-300))

st["neglog10_p"] = -np.log10(st["p_value"].clip(lower=1e-300))
st["neglog10_fdr"] = -np.log10(st["p_adj"].clip(lower=1e-300))

# Categories
scrna["category"] = "Other genes"
scrna.loc[
    (scrna["scRNA_logFC_Ischemic_vs_Control"] > 0) & (scrna["p_value"] < 0.05),
    "category"
] = "Upregulated, nominal P < 0.05"
scrna.loc[
    (scrna["scRNA_logFC_Ischemic_vs_Control"] < 0) & (scrna["p_value"] < 0.05),
    "category"
] = "Downregulated, nominal P < 0.05"
scrna.loc[
    (scrna["scRNA_logFC_Ischemic_vs_Control"] > 0) & (scrna["p_adj"] < 0.05),
    "category"
] = "Upregulated, FDR < 0.05"

st["category"] = "Other genes"
st.loc[
    (st["ST_delta_Pathological_vs_Control"] > 0) & (st["p_value"] < 0.05),
    "category"
] = "Upregulated, nominal P < 0.05"
st.loc[
    (st["ST_delta_Pathological_vs_Control"] < 0) & (st["p_value"] < 0.05),
    "category"
] = "Downregulated, nominal P < 0.05"
st.loc[
    (st["ST_delta_Pathological_vs_Control"] > 0) & (st["p_adj"] < 0.10),
    "category"
] = "Upregulated, FDR < 0.10"

# ============================================================
# Figure 1: scRNA fibroblast pseudobulk DEG
# ============================================================
fig, ax = plt.subplots(figsize=(3.45, 3.1))

plot_order = [
    ("Other genes", COL_BG, 7, 0.50),
    ("Downregulated, nominal P < 0.05", COL_DOWN, 9, 0.70),
    ("Upregulated, nominal P < 0.05", COL_UP_LIGHT, 9, 0.75),
    ("Upregulated, FDR < 0.05", COL_UP, 11, 0.85),
]

for cat, color, size, alpha in plot_order:
    sub = scrna[scrna["category"] == cat]
    if len(sub) == 0:
        continue
    ax.scatter(
        sub["scRNA_logFC_Ischemic_vs_Control"],
        sub["neglog10_p"],
        s=size,
        color=color,
        alpha=alpha,
        linewidths=0,
        rasterized=True,
        zorder=2
    )

ax.axvline(0, color=COL_LINE, lw=0.6, ls="--", zorder=1)
ax.axhline(-np.log10(0.05), color=COL_LINE, lw=0.6, ls="--", zorder=1)

for g in SCRNA_LABELS:
    sub = scrna[scrna["gene"] == g]
    if len(sub) == 0:
        continue
    x = sub["scRNA_logFC_Ischemic_vs_Control"].iloc[0]
    y = sub["neglog10_p"].iloc[0]
    add_gene_label(ax, x, y, g, SCRNA_OFFSETS.get(g, (6, 4)))

ax.set_xlabel("Fibroblast pseudobulk logFC\n(Ischemic vs Control)")
ax.set_ylabel("-log10(P value)")
ax.set_title("Fibroblast pseudobulk DEG", pad=4)
format_axis(ax)

# Small method note inside plot
ax.text(
    0.02, 0.98,
    "Control n = 6; Ischemic n = 10",
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=5.8,
    color="#555555"
)

legend_handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_BG, markeredgecolor="none", markersize=4, label="Other genes"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_UP_LIGHT, markeredgecolor="none", markersize=4, label="Nominal up"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_UP, markeredgecolor="none", markersize=4, label="FDR up"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_FSTL1, markeredgecolor="white", markersize=4, label="FSTL1"),
]
ax.legend(
    handles=legend_handles,
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(0.00, -0.28),
    ncol=4,
    handletextpad=0.3,
    columnspacing=0.8
)

save_all(fig, "Figure_crossomics_1_scRNA_fibroblast_pseudobulk_volcano_beautified")

# ============================================================
# Figure 2: ST section-level screening
# ============================================================
fig, ax = plt.subplots(figsize=(3.45, 3.1))

plot_order_st = [
    ("Other genes", COL_BG, 7, 0.50),
    ("Downregulated, nominal P < 0.05", COL_DOWN, 9, 0.70),
    ("Upregulated, nominal P < 0.05", COL_UP_LIGHT, 9, 0.75),
    ("Upregulated, FDR < 0.10", COL_UP, 11, 0.85),
]

for cat, color, size, alpha in plot_order_st:
    sub = st[st["category"] == cat]
    if len(sub) == 0:
        continue
    ax.scatter(
        sub["ST_delta_Pathological_vs_Control"],
        sub["neglog10_p"],
        s=size,
        color=color,
        alpha=alpha,
        linewidths=0,
        rasterized=True,
        zorder=2
    )

ax.axvline(0, color=COL_LINE, lw=0.6, ls="--", zorder=1)
ax.axhline(-np.log10(0.05), color=COL_LINE, lw=0.6, ls="--", zorder=1)

for g in ST_LABELS:
    sub = st[st["gene"] == g]
    if len(sub) == 0:
        continue
    x = sub["ST_delta_Pathological_vs_Control"].iloc[0]
    y = sub["neglog10_p"].iloc[0]
    add_gene_label(ax, x, y, g, ST_OFFSETS.get(g, (6, 4)))

ax.set_xlabel("Section-level expression difference\n(Pathological - Control)")
ax.set_ylabel("-log10(P value)")
ax.set_title("ST section-level screening", pad=4)
format_axis(ax)

ax.text(
    0.02, 0.98,
    "Section used as statistical unit",
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=5.8,
    color="#555555"
)

legend_handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_BG, markeredgecolor="none", markersize=4, label="Other genes"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_UP_LIGHT, markeredgecolor="none", markersize=4, label="Nominal up"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_UP, markeredgecolor="none", markersize=4, label="FDR < 0.10 up"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_FSTL1, markeredgecolor="white", markersize=4, label="FSTL1"),
]
ax.legend(
    handles=legend_handles,
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(0.00, -0.28),
    ncol=4,
    handletextpad=0.3,
    columnspacing=0.8
)

save_all(fig, "Figure_crossomics_2_ST_section_level_screening_volcano_beautified")

# ============================================================
# Figure 3: Cross-omics scatter
# ============================================================
cross["concordant_up"] = cross["concordant_up"].astype(bool)

fig, ax = plt.subplots(figsize=(3.65, 3.35))

# Subtle quadrant shade for concordant-up area
xmax = np.nanpercentile(cross["scRNA_logFC_Ischemic_vs_Control"], 99.5)
ymax = np.nanpercentile(cross["ST_delta_Pathological_vs_Control"], 99.5)
xmin = np.nanpercentile(cross["scRNA_logFC_Ischemic_vs_Control"], 0.5)
ymin = np.nanpercentile(cross["ST_delta_Pathological_vs_Control"], 0.5)

# Keep axes wide enough for POSTN
xmax = max(xmax, cross.loc[cross["gene"] == "POSTN", "scRNA_logFC_Ischemic_vs_Control"].max() + 0.4)
ymax = max(ymax, cross.loc[cross["gene"] == "COL1A1", "ST_delta_Pathological_vs_Control"].max() + 0.25)

ax.add_patch(
    Rectangle(
        (0, 0),
        xmax,
        ymax,
        facecolor="#F8E7D2",
        edgecolor="none",
        alpha=0.55,
        zorder=0
    )
)

other = cross[~cross["concordant_up"]]
cu = cross[cross["concordant_up"]]

ax.scatter(
    other["scRNA_logFC_Ischemic_vs_Control"],
    other["ST_delta_Pathological_vs_Control"],
    s=7,
    color=COL_BG,
    alpha=0.48,
    linewidths=0,
    rasterized=True,
    zorder=1
)

ax.scatter(
    cu["scRNA_logFC_Ischemic_vs_Control"],
    cu["ST_delta_Pathological_vs_Control"],
    s=8,
    color=COL_UP,
    alpha=0.55,
    linewidths=0,
    rasterized=True,
    zorder=2
)

ax.axvline(0, color=COL_LINE, lw=0.65, ls="--", zorder=3)
ax.axhline(0, color=COL_LINE, lw=0.65, ls="--", zorder=3)

for g in CROSS_LABELS:
    sub = cross[cross["gene"] == g]
    if len(sub) == 0:
        continue
    x = sub["scRNA_logFC_Ischemic_vs_Control"].iloc[0]
    y = sub["ST_delta_Pathological_vs_Control"].iloc[0]
    add_gene_label(ax, x, y, g, CROSS_OFFSETS.get(g, (6, 4)))

ax.set_xlim(min(xmin, -4.4), xmax)
ax.set_ylim(min(ymin, -1.8), ymax)
ax.set_xlabel("scRNA-seq fibroblast logFC\n(Ischemic vs Control)")
ax.set_ylabel("ST section-level difference\n(Pathological - Control)")
ax.set_title("Cross-omics gene prioritization", pad=4)
format_axis(ax)

ax.text(
    0.98, 0.96,
    "Concordant up",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=6.2,
    color=COL_UP
)

legend_handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_BG, markeredgecolor="none", markersize=4, label="Other genes"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_UP, markeredgecolor="none", markersize=4, label="Concordant up"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_FSTL1, markeredgecolor="white", markersize=4, label="FSTL1"),
]
ax.legend(
    handles=legend_handles,
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(0.00, -0.27),
    ncol=3,
    handletextpad=0.3,
    columnspacing=0.9
)

save_all(fig, "Figure_crossomics_3_scRNA_ST_scatter_highlight_FSTL1_beautified")

# ============================================================
# Combined 3-panel version
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(9.8, 3.25))

# Panel A
ax = axes[0]
for cat, color, size, alpha in plot_order:
    sub = scrna[scrna["category"] == cat]
    if len(sub):
        ax.scatter(
            sub["scRNA_logFC_Ischemic_vs_Control"],
            sub["neglog10_p"],
            s=size,
            color=color,
            alpha=alpha,
            linewidths=0,
            rasterized=True
        )
ax.axvline(0, color=COL_LINE, lw=0.6, ls="--")
ax.axhline(-np.log10(0.05), color=COL_LINE, lw=0.6, ls="--")
for g in SCRNA_LABELS:
    sub = scrna[scrna["gene"] == g]
    if len(sub):
        add_gene_label(
            ax,
            sub["scRNA_logFC_Ischemic_vs_Control"].iloc[0],
            sub["neglog10_p"].iloc[0],
            g,
            SCRNA_OFFSETS.get(g, (6, 4)),
            size=6.2
        )
ax.set_xlabel("Fibroblast logFC\n(Ischemic vs Control)")
ax.set_ylabel("-log10(P value)")
ax.set_title("Fibroblast pseudobulk DEG", pad=4)
format_axis(ax)
add_panel_label(ax, "A")

# Panel B
ax = axes[1]
for cat, color, size, alpha in plot_order_st:
    sub = st[st["category"] == cat]
    if len(sub):
        ax.scatter(
            sub["ST_delta_Pathological_vs_Control"],
            sub["neglog10_p"],
            s=size,
            color=color,
            alpha=alpha,
            linewidths=0,
            rasterized=True
        )
ax.axvline(0, color=COL_LINE, lw=0.6, ls="--")
ax.axhline(-np.log10(0.05), color=COL_LINE, lw=0.6, ls="--")
for g in ST_LABELS:
    sub = st[st["gene"] == g]
    if len(sub):
        add_gene_label(
            ax,
            sub["ST_delta_Pathological_vs_Control"].iloc[0],
            sub["neglog10_p"].iloc[0],
            g,
            ST_OFFSETS.get(g, (6, 4)),
            size=6.2
        )
ax.set_xlabel("Section-level difference\n(Pathological - Control)")
ax.set_ylabel("-log10(P value)")
ax.set_title("ST section-level screening", pad=4)
format_axis(ax)
add_panel_label(ax, "B")

# Panel C
ax = axes[2]
ax.add_patch(
    Rectangle(
        (0, 0),
        xmax,
        ymax,
        facecolor="#F8E7D2",
        edgecolor="none",
        alpha=0.55,
        zorder=0
    )
)
ax.scatter(
    other["scRNA_logFC_Ischemic_vs_Control"],
    other["ST_delta_Pathological_vs_Control"],
    s=7,
    color=COL_BG,
    alpha=0.48,
    linewidths=0,
    rasterized=True
)
ax.scatter(
    cu["scRNA_logFC_Ischemic_vs_Control"],
    cu["ST_delta_Pathological_vs_Control"],
    s=8,
    color=COL_UP,
    alpha=0.55,
    linewidths=0,
    rasterized=True
)
ax.axvline(0, color=COL_LINE, lw=0.65, ls="--")
ax.axhline(0, color=COL_LINE, lw=0.65, ls="--")
for g in CROSS_LABELS:
    sub = cross[cross["gene"] == g]
    if len(sub):
        add_gene_label(
            ax,
            sub["scRNA_logFC_Ischemic_vs_Control"].iloc[0],
            sub["ST_delta_Pathological_vs_Control"].iloc[0],
            g,
            CROSS_OFFSETS.get(g, (6, 4)),
            size=6.2
        )
ax.set_xlim(min(xmin, -4.4), xmax)
ax.set_ylim(min(ymin, -1.8), ymax)
ax.set_xlabel("scRNA-seq fibroblast logFC\n(Ischemic vs Control)")
ax.set_ylabel("ST section-level difference\n(Pathological - Control)")
ax.set_title("Cross-omics prioritization", pad=4)
format_axis(ax)
add_panel_label(ax, "C")

legend_handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_BG, markeredgecolor="none", markersize=4, label="Other genes"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_UP_LIGHT, markeredgecolor="none", markersize=4, label="Nominal up"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_UP, markeredgecolor="none", markersize=4, label="FDR/concordant up"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=COL_FSTL1, markeredgecolor="white", markersize=4, label="FSTL1"),
]
fig.legend(
    handles=legend_handles,
    frameon=False,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.04),
    ncol=4,
    handletextpad=0.3,
    columnspacing=1.2
)

fig.tight_layout(w_pad=1.8)
save_all(fig, "Figure_crossomics_screening_three_panel_beautified")

print("Saved beautified figures to:", outdir)
print("Key outputs:")
print(outdir / "Figure_crossomics_1_scRNA_fibroblast_pseudobulk_volcano_beautified.pdf")
print(outdir / "Figure_crossomics_2_ST_section_level_screening_volcano_beautified.pdf")
print(outdir / "Figure_crossomics_3_scRNA_ST_scatter_highlight_FSTL1_beautified.pdf")
print(outdir / "Figure_crossomics_screening_three_panel_beautified.pdf")
