import os
import re
import gzip
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu, kruskal

warnings.filterwarnings("ignore")

# =========================
# Paths
# =========================
PROJECT_ROOT = Path(__file__).resolve().parents[2]

data_dir = PROJECT_ROOT / "data" / "spatial"

out_dir = PROJECT_ROOT / "results" / "spatial_FSTL1"

map_dir = out_dir / "individual_spatial_maps"

out_dir.mkdir(parents=True, exist_ok=True)
map_dir.mkdir(parents=True, exist_ok=True)
map_dir = out_dir / "individual_spatial_maps"
out_dir.mkdir(parents=True, exist_ok=True)
map_dir.mkdir(parents=True, exist_ok=True)

# =========================
# Figure style: Phytomedicine-like
# =========================
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 0.6,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight"
})

region_order = [
    "Control",
    "GT",
    "IZ",
    "BZ",
    "RZ",
    "FZ",
    "GT_IZ",
    "IZ_BZ",
    "RZ_BZ",
    "RZ_FZ",
    "RZ_GT",
    "FZ_GT",
    "Unknown"
]

region_palette = {
    "Control": "#4C78A8",
    "GT": "#72B7B2",
    "IZ": "#F58518",
    "BZ": "#E45756",
    "RZ": "#54A24B",
    "FZ": "#B279A2",
    "GT_IZ": "#FF9DA6",
    "IZ_BZ": "#9D755D",
    "RZ_BZ": "#59A14F",
    "RZ_FZ": "#AF7AA1",
    "RZ_GT": "#76B7B2",
    "FZ_GT": "#EDC948",
    "Unknown": "#BAB0AC"
}

genes_to_plot = [
    "FSTL1", "COL1A1", "COL1A2", "COL3A1",
    "POSTN", "DCN", "LUM", "FN1", "ACTA2"
]

ecm_genes = [
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1",
    "POSTN", "DCN", "LUM", "FN1", "THBS2", "ACTA2"
]

celltype_cols = [
    "Adipocyte", "Cardiomyocyte", "Endothelial", "Fibroblast",
    "Lymphoid", "Mast", "Myeloid", "Neuronal", "Pericyte",
    "Cycling.cells", "vSMCs"
]

# =========================
# Utilities
# =========================
def infer_region(sample_name: str) -> str:
    s = sample_name.replace("Visium_", "").replace(".h5ad", "")
    s = re.sub(r"_P\d+.*$", "", s)
    if s.lower().startswith("control"):
        return "Control"
    return s if s else "Unknown"

def infer_patient(sample_name: str) -> str:
    m = re.search(r"_P(\d+)", sample_name)
    return f"P{m.group(1)}" if m else "Unknown"

def standardize_gene_names(adata):
    adata.var_names_make_unique()

    upper_names = [str(g).upper() for g in adata.var_names]
    if "FSTL1" in upper_names:
        adata.var_names = upper_names
        adata.var_names_make_unique()
        return adata

    candidate_cols = ["features", "gene_symbols", "gene_symbol", "gene", "symbol", "name"]
    for col in candidate_cols:
        if col in adata.var.columns:
            symbols = [str(x).upper() for x in adata.var[col].values]
            if "FSTL1" in symbols:
                adata.var_names = symbols
                adata.var_names_make_unique()
                return adata

    adata.var_names = upper_names
    adata.var_names_make_unique()
    return adata

def normalize_if_needed(adata):
    ad = adata.copy()
    X = ad.X
    if hasattr(X, "data"):
        max_val = X.data.max() if X.data.size > 0 else 0
    else:
        max_val = np.max(X)

    if max_val > 50:
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)

    return ad

def get_expr(adata, gene):
    gene = gene.upper()
    if gene not in adata.var_names:
        return None
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        return x.toarray().ravel()
    return np.asarray(x).ravel()

def get_spatial_coords(adata):
    if "spatial" in adata.obsm.keys():
        return adata.obsm["spatial"]
    if "X_spatial" in adata.obsm.keys():
        return adata.obsm["X_spatial"]
    return None

def compute_module_score(adata, gene_list):
    available = [g for g in gene_list if g in adata.var_names]
    if len(available) < 3:
        return None, available

    mats = []
    for g in available:
        v = get_expr(adata, g)
        if v is None:
            continue
        sd = np.std(v)
        if sd == 0:
            z = np.zeros_like(v)
        else:
            z = (v - np.mean(v)) / sd
        mats.append(z)

    if len(mats) < 3:
        return None, available

    score = np.vstack(mats).mean(axis=0)
    return score, available

def savefig_multi(fig, basename):
    for ext in ["pdf", "png", "tiff"]:
        dpi = 600 if ext in ["tiff", "png"] else 300
        fig.savefig(out_dir / f"{basename}.{ext}", dpi=dpi)
    plt.close(fig)

def savefig_map(fig, basename):
    for ext in ["pdf", "png"]:
        dpi = 600 if ext == "png" else 300
        fig.savefig(map_dir / f"{basename}.{ext}", dpi=dpi)
    plt.close(fig)

def clean_label(x):
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(x))

def add_panel_label(ax, label):
    ax.text(
        -0.08, 1.05, label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        ha="right"
    )

def jitter(n, width=0.18):
    if n <= 1:
        return np.zeros(n)
    return np.random.uniform(-width, width, size=n)

def box_strip(ax, df, x_col, y_col, ylabel, title=None, order=None):
    if order is None:
        order = [x for x in region_order if x in df[x_col].unique()]
    data = [df.loc[df[x_col] == x, y_col].dropna().values for x in order]
    positions = np.arange(len(order))

    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=0.8),
        whiskerprops=dict(linewidth=0.7),
        capprops=dict(linewidth=0.7),
        boxprops=dict(linewidth=0.7)
    )

    for patch, cat in zip(bp["boxes"], order):
        patch.set_facecolor(region_palette.get(cat, "#CCCCCC"))
        patch.set_alpha(0.45)

    for i, cat in enumerate(order):
        vals = df.loc[df[x_col] == cat, y_col].dropna().values
        ax.scatter(
            np.full(len(vals), i) + jitter(len(vals)),
            vals,
            s=14,
            color=region_palette.get(cat, "#666666"),
            edgecolor="black",
            linewidth=0.25,
            zorder=3
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def spatial_scatter(ax, coords, values, title, cmap="viridis", vmin=None, vmax=None):
    sca = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=values,
        s=4,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0
    )
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, pad=2)
    return sca

# =========================
# Main analysis
# =========================
files = sorted(list(data_dir.glob("*.h5ad")))
if len(files) == 0:
    raise FileNotFoundError(f"No h5ad files found in: {data_dir}")

print(f"Detected {len(files)} h5ad files.")

section_rows = []
correlation_rows = []
highlow_rows = []
spot_rows = []

# 用于代表性空间图
representative_cache = {}

for path in files:
    sample = path.name.replace(".h5ad", "")
    region = infer_region(sample)
    patient = infer_patient(sample)

    print(f"\nProcessing {sample} | region={region} | patient={patient}")

    adata = sc.read_h5ad(path)
    adata = standardize_gene_names(adata)
    coords = get_spatial_coords(adata)

    if coords is None:
        print(f"[Skip] no spatial coordinate: {sample}")
        continue

    adata = normalize_if_needed(adata)
    coords = get_spatial_coords(adata)

    fstl1 = get_expr(adata, "FSTL1")
    if fstl1 is None:
        print(f"[Skip] FSTL1 missing: {sample}")
        continue

    ecm_score, ecm_available = compute_module_score(adata, ecm_genes)
    if ecm_score is None:
        print(f"[Warning] ECM score unavailable: {sample}")
        ecm_score = np.full(adata.n_obs, np.nan)

    # 主要细胞类型分数
    obs_cols = adata.obs.columns
    ct_data = {}
    for ct in celltype_cols:
        if ct in obs_cols:
            ct_data[ct] = adata.obs[ct].astype(float).values

    # FSTL1-high spots: section内top 25%
    q75 = np.quantile(fstl1, 0.75)
    high = fstl1 >= q75
    low = fstl1 < q75

    section_row = {
        "sample": sample,
        "patient": patient,
        "region": region,
        "n_spots": int(adata.n_obs),
        "fstl1_mean": float(np.mean(fstl1)),
        "fstl1_median": float(np.median(fstl1)),
        "fstl1_max": float(np.max(fstl1)),
        "fstl1_pct_positive": float(np.mean(fstl1 > 0) * 100),
        "ecm_score_mean": float(np.nanmean(ecm_score)),
        "ecm_score_median": float(np.nanmedian(ecm_score)),
        "ecm_genes_used": ";".join(ecm_available),
        "fstl1_high_threshold_q75": float(q75),
        "fstl1_high_n_spots": int(np.sum(high)),
        "fstl1_low_n_spots": int(np.sum(low)),
        "ecm_score_high_mean": float(np.nanmean(ecm_score[high])),
        "ecm_score_low_mean": float(np.nanmean(ecm_score[low])),
        "ecm_score_delta_high_minus_low": float(np.nanmean(ecm_score[high]) - np.nanmean(ecm_score[low]))
    }

    for ct, val in ct_data.items():
        section_row[f"{ct}_mean"] = float(np.mean(val))
        section_row[f"{ct}_high_mean"] = float(np.mean(val[high]))
        section_row[f"{ct}_low_mean"] = float(np.mean(val[low]))
        section_row[f"{ct}_delta_high_minus_low"] = float(np.mean(val[high]) - np.mean(val[low]))

    section_rows.append(section_row)

    # 相关性：FSTL1 vs ECM/cell type
    rho, pval = spearmanr(fstl1, ecm_score, nan_policy="omit")
    correlation_rows.append({
        "sample": sample, "patient": patient, "region": region,
        "feature": "ECM_score",
        "spearman_rho": float(rho),
        "p_value": float(pval),
        "n_spots": int(adata.n_obs)
    })

    for ct, val in ct_data.items():
        rho, pval = spearmanr(fstl1, val)
        correlation_rows.append({
            "sample": sample, "patient": patient, "region": region,
            "feature": ct,
            "spearman_rho": float(rho),
            "p_value": float(pval),
            "n_spots": int(adata.n_obs)
        })

        try:
            _, mw_p = mannwhitneyu(val[high], val[low], alternative="two-sided")
        except Exception:
            mw_p = np.nan

        highlow_rows.append({
            "sample": sample,
            "patient": patient,
            "region": region,
            "feature": ct,
            "FSTL1_high_mean": float(np.mean(val[high])),
            "FSTL1_low_mean": float(np.mean(val[low])),
            "delta_high_minus_low": float(np.mean(val[high]) - np.mean(val[low])),
            "mannwhitney_p_spot_level": float(mw_p) if not np.isnan(mw_p) else np.nan,
            "n_high_spots": int(np.sum(high)),
            "n_low_spots": int(np.sum(low))
        })

    try:
        _, mw_p = mannwhitneyu(ecm_score[high], ecm_score[low], alternative="two-sided")
    except Exception:
        mw_p = np.nan

    highlow_rows.append({
        "sample": sample,
        "patient": patient,
        "region": region,
        "feature": "ECM_score",
        "FSTL1_high_mean": float(np.nanmean(ecm_score[high])),
        "FSTL1_low_mean": float(np.nanmean(ecm_score[low])),
        "delta_high_minus_low": float(np.nanmean(ecm_score[high]) - np.nanmean(ecm_score[low])),
        "mannwhitney_p_spot_level": float(mw_p) if not np.isnan(mw_p) else np.nan,
        "n_high_spots": int(np.sum(high)),
        "n_low_spots": int(np.sum(low))
    })

    # 保存 spot-level 表：只保留关键列
    spot_df = pd.DataFrame({
        "sample": sample,
        "patient": patient,
        "region": region,
        "spot_id": adata.obs_names.astype(str),
        "x": coords[:, 0],
        "y": coords[:, 1],
        "FSTL1": fstl1,
        "ECM_score": ecm_score,
        "FSTL1_high_top25": high.astype(int)
    })
    for ct, val in ct_data.items():
        spot_df[ct] = val
    spot_rows.append(spot_df)

    # 个体空间图：FSTL1, ECM_score, Fibroblast
    for feature_name, values, cmap in [
        ("FSTL1", fstl1, "viridis"),
        ("ECM_score", ecm_score, "magma"),
        ("Fibroblast", ct_data.get("Fibroblast", None), "plasma")
    ]:
        if values is None:
            continue
        fig, ax = plt.subplots(figsize=(2.1, 2.1))
        vmax = np.nanpercentile(values, 99)
        vmin = np.nanpercentile(values, 1) if feature_name != "FSTL1" else 0
        sca = spatial_scatter(
            ax, coords, values,
            title=f"{sample}\n{feature_name}",
            cmap=cmap, vmin=vmin, vmax=vmax
        )
        cbar = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
        cbar.ax.tick_params(labelsize=5)
        fig.tight_layout(pad=0.3)
        savefig_map(fig, f"{sample}_{feature_name}_spatial")

    # 代表性缓存
    representative_cache[sample] = {
        "region": region,
        "coords": coords,
        "FSTL1": fstl1,
        "ECM_score": ecm_score,
        "Fibroblast": ct_data.get("Fibroblast", None)
    }

# =========================
# Save tables
# =========================
section_df = pd.DataFrame(section_rows)
corr_df = pd.DataFrame(correlation_rows)
highlow_df = pd.DataFrame(highlow_rows)

section_df.to_csv(out_dir / "Kuppe_all_section_level_summary.csv", index=False)
corr_df.to_csv(out_dir / "Kuppe_all_FSTL1_correlation_summary.csv", index=False)
highlow_df.to_csv(out_dir / "Kuppe_all_FSTL1_high_low_summary.csv", index=False)

if len(spot_rows) > 0:
    all_spots_df = pd.concat(spot_rows, ignore_index=True)
    all_spots_df.to_csv(out_dir / "Kuppe_all_spot_level_key_metrics.csv.gz", index=False, compression="gzip")

# =========================
# Region-level statistics
# =========================
stats_rows = []
main_metrics = [
    "fstl1_mean",
    "fstl1_pct_positive",
    "ecm_score_delta_high_minus_low",
    "Fibroblast_delta_high_minus_low",
    "Cardiomyocyte_delta_high_minus_low"
]

for metric in main_metrics:
    if metric not in section_df.columns:
        continue

    groups = []
    labels = []
    for r in [x for x in region_order if x in section_df["region"].unique()]:
        vals = section_df.loc[section_df["region"] == r, metric].dropna().values
        if len(vals) > 0:
            groups.append(vals)
            labels.append(r)

    if len(groups) >= 2:
        try:
            kw_stat, kw_p = kruskal(*groups)
        except Exception:
            kw_stat, kw_p = np.nan, np.nan
    else:
        kw_stat, kw_p = np.nan, np.nan

    stats_rows.append({
        "metric": metric,
        "test": "Kruskal-Wallis across regions",
        "comparison": "all regions",
        "statistic": kw_stat,
        "p_value": kw_p
    })

    # 与 Control 比较，section-level Mann-Whitney
    control_vals = section_df.loc[section_df["region"] == "Control", metric].dropna().values
    for r in labels:
        if r == "Control":
            continue
        vals = section_df.loc[section_df["region"] == r, metric].dropna().values
        if len(control_vals) >= 2 and len(vals) >= 2:
            try:
                stat, p = mannwhitneyu(control_vals, vals, alternative="two-sided")
            except Exception:
                stat, p = np.nan, np.nan
        else:
            stat, p = np.nan, np.nan

        stats_rows.append({
            "metric": metric,
            "test": "Mann-Whitney vs Control",
            "comparison": f"{r} vs Control",
            "statistic": stat,
            "p_value": p
        })

stats_df = pd.DataFrame(stats_rows)
stats_df.to_csv(out_dir / "Kuppe_all_section_level_statistics.csv", index=False)

# =========================
# Figure 1: representative spatial maps
# =========================
def pick_sample_contains(candidates):
    for key in candidates:
        for s in representative_cache.keys():
            if key in s:
                return s
    return None

rep_samples = [
    pick_sample_contains(["Visium_control_P1", "Visium_control_P7", "control"]),
    pick_sample_contains(["Visium_IZ_BZ_P2", "Visium_IZ_P3", "Visium_IZ"]),
    pick_sample_contains(["Visium_RZ_BZ_P3", "Visium_RZ_BZ_P12", "Visium_RZ_BZ"]),
    pick_sample_contains(["Visium_FZ_P14", "Visium_FZ_P18", "Visium_FZ"])
]
rep_samples = [s for s in rep_samples if s is not None]

if len(rep_samples) > 0:
    fig, axes = plt.subplots(
        3, len(rep_samples),
        figsize=(1.8 * len(rep_samples), 5.2),
        constrained_layout=True
    )
    if len(rep_samples) == 1:
        axes = np.array(axes).reshape(3, 1)

    all_fstl1 = np.concatenate([representative_cache[s]["FSTL1"] for s in rep_samples])
    all_ecm = np.concatenate([representative_cache[s]["ECM_score"] for s in rep_samples])
    all_fib = np.concatenate([
        representative_cache[s]["Fibroblast"]
        for s in rep_samples
        if representative_cache[s]["Fibroblast"] is not None
    ])

    vmax_fstl1 = np.nanpercentile(all_fstl1, 99)
    vmax_ecm = np.nanpercentile(all_ecm, 99)
    vmin_ecm = np.nanpercentile(all_ecm, 1)
    vmax_fib = np.nanpercentile(all_fib, 99) if len(all_fib) > 0 else None

    panel_letters = list("ABCDEFGHIJKL")
    k = 0
    for j, s in enumerate(rep_samples):
        dat = representative_cache[s]
        title_suffix = dat["region"]

        sca = spatial_scatter(
            axes[0, j], dat["coords"], dat["FSTL1"],
            title=f"{title_suffix}\nFSTL1",
            cmap="viridis", vmin=0, vmax=vmax_fstl1
        )
        add_panel_label(axes[0, j], panel_letters[k]); k += 1

        sca = spatial_scatter(
            axes[1, j], dat["coords"], dat["ECM_score"],
            title=f"{title_suffix}\nECM score",
            cmap="magma", vmin=vmin_ecm, vmax=vmax_ecm
        )
        add_panel_label(axes[1, j], panel_letters[k]); k += 1

        if dat["Fibroblast"] is not None:
            sca = spatial_scatter(
                axes[2, j], dat["coords"], dat["Fibroblast"],
                title=f"{title_suffix}\nFibroblast score",
                cmap="plasma", vmin=0, vmax=vmax_fib
            )
        else:
            axes[2, j].axis("off")
        add_panel_label(axes[2, j], panel_letters[k]); k += 1

    savefig_multi(fig, "Figure_1_Kuppe_representative_spatial_maps")

# =========================
# Figure 2: section-level quantification
# =========================
available_regions = [r for r in region_order if r in section_df["region"].unique()]
fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.3), constrained_layout=True)

box_strip(
    axes[0, 0], section_df, "region", "fstl1_mean",
    ylabel="Mean FSTL1 expression",
    title="Section-level FSTL1 expression",
    order=available_regions
)
add_panel_label(axes[0, 0], "A")

box_strip(
    axes[0, 1], section_df, "region", "fstl1_pct_positive",
    ylabel="FSTL1-positive spots (%)",
    title="FSTL1-positive spatial spots",
    order=available_regions
)
add_panel_label(axes[0, 1], "B")

# Correlation FSTL1 vs ECM_score
plot_corr_ecm = corr_df[corr_df["feature"] == "ECM_score"].copy()
box_strip(
    axes[1, 0], plot_corr_ecm, "region", "spearman_rho",
    ylabel="Spearman rho",
    title="FSTL1 vs ECM score",
    order=[r for r in available_regions if r in plot_corr_ecm["region"].unique()]
)
axes[1, 0].axhline(0, color="black", linewidth=0.6, linestyle="--")
add_panel_label(axes[1, 0], "C")

plot_corr_fib = corr_df[corr_df["feature"] == "Fibroblast"].copy()
box_strip(
    axes[1, 1], plot_corr_fib, "region", "spearman_rho",
    ylabel="Spearman rho",
    title="FSTL1 vs Fibroblast score",
    order=[r for r in available_regions if r in plot_corr_fib["region"].unique()]
)
axes[1, 1].axhline(0, color="black", linewidth=0.6, linestyle="--")
add_panel_label(axes[1, 1], "D")

savefig_multi(fig, "Figure_2_Kuppe_section_level_quantification")

# =========================
# Figure 3: FSTL1-high vs FSTL1-low enrichment
# =========================
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), constrained_layout=True)

box_strip(
    axes[0], section_df, "region", "ecm_score_delta_high_minus_low",
    ylabel="Δ ECM score\n(FSTL1-high − low)",
    title="ECM enrichment in FSTL1-high spots",
    order=available_regions
)
axes[0].axhline(0, color="black", linewidth=0.6, linestyle="--")
add_panel_label(axes[0], "A")

if "Fibroblast_delta_high_minus_low" in section_df.columns:
    box_strip(
        axes[1], section_df, "region", "Fibroblast_delta_high_minus_low",
        ylabel="Δ Fibroblast score\n(FSTL1-high − low)",
        title="Fibroblast enrichment in FSTL1-high spots",
        order=available_regions
    )
    axes[1].axhline(0, color="black", linewidth=0.6, linestyle="--")
add_panel_label(axes[1], "B")

if "Cardiomyocyte_delta_high_minus_low" in section_df.columns:
    box_strip(
        axes[2], section_df, "region", "Cardiomyocyte_delta_high_minus_low",
        ylabel="Δ Cardiomyocyte score\n(FSTL1-high − low)",
        title="Cardiomyocyte score in FSTL1-high spots",
        order=available_regions
    )
    axes[2].axhline(0, color="black", linewidth=0.6, linestyle="--")
add_panel_label(axes[2], "C")

savefig_multi(fig, "Figure_3_Kuppe_FSTL1_high_spot_enrichment")

# =========================
# Figure 4: celltype correlation heatmap
# =========================
heat_df = corr_df[corr_df["feature"].isin(celltype_cols + ["ECM_score"])].copy()
heat_mean = (
    heat_df
    .groupby(["feature", "region"], as_index=False)["spearman_rho"]
    .mean()
)

features_order = ["ECM_score", "Fibroblast", "Cardiomyocyte", "Endothelial", "Pericyte", "vSMCs", "Myeloid", "Lymphoid", "Mast", "Adipocyte"]
features_order = [f for f in features_order if f in heat_mean["feature"].unique()]
regions_heat = [r for r in available_regions if r in heat_mean["region"].unique()]

mat = np.full((len(features_order), len(regions_heat)), np.nan)
for i, f in enumerate(features_order):
    for j, r in enumerate(regions_heat):
        vals = heat_mean.loc[(heat_mean["feature"] == f) & (heat_mean["region"] == r), "spearman_rho"].values
        if len(vals) > 0:
            mat[i, j] = vals[0]

fig, ax = plt.subplots(figsize=(7.2, 3.4), constrained_layout=True)
im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
ax.set_xticks(np.arange(len(regions_heat)))
ax.set_xticklabels(regions_heat, rotation=45, ha="right")
ax.set_yticks(np.arange(len(features_order)))
ax.set_yticklabels(features_order)
ax.set_title("Mean section-level spatial correlation with FSTL1")
add_panel_label(ax, "A")

for i in range(len(features_order)):
    for j in range(len(regions_heat)):
        if not np.isnan(mat[i, j]):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=5.5)

cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
cbar.set_label("Spearman rho")
savefig_multi(fig, "Figure_4_Kuppe_FSTL1_celltype_correlation_heatmap")

print("\nAnalysis complete.")
print(f"Results saved to: {out_dir}")
print("Main files:")
print("  - Kuppe_all_section_level_summary.csv")
print("  - Kuppe_all_FSTL1_correlation_summary.csv")
print("  - Kuppe_all_FSTL1_high_low_summary.csv")
print("  - Kuppe_all_section_level_statistics.csv")
print("  - Figure_1_Kuppe_representative_spatial_maps.pdf/png/tiff")
print("  - Figure_2_Kuppe_section_level_quantification.pdf/png/tiff")
print("  - Figure_3_Kuppe_FSTL1_high_spot_enrichment.pdf/png/tiff")
print("  - Figure_4_Kuppe_FSTL1_celltype_correlation_heatmap.pdf/png/tiff")
