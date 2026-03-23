"""
Tissue enrichment and combined scoring utilities.
"""

# ============================================================
# Imports
# ============================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.stats import norm, mannwhitneyu
from statsmodels.stats.multitest import multipletests


# ============================================================
# 1) Combined score construction
# ============================================================
def make_scores(df, lfc_col="log2FoldChange", p_col="pvalue",
                use_padj=True, padj_col="padj",
                lfc_shrunk_col=None,  
                positive_only=True, z_cap=10.0):
    """
    Create a robust 'combined_score' that's stable and directional.
    """
    # choose p: prefer adjusted if available (more conservative, better across many tests)
    p = df[padj_col].to_numpy() if (use_padj and padj_col in df) else df[p_col].to_numpy()
    p = np.clip(p.astype(float), 1e-300, 1.0)

    # two-sided Z; cap extremes to avoid infinities
    z = np.minimum(norm.isf(p / 2.0), z_cap)

    # pick LFC (shrunk preferred)
    lfc = df[lfc_shrunk_col].to_numpy() if (lfc_shrunk_col and lfc_shrunk_col in df) else df[lfc_col].to_numpy()

    # directional Z (sign from LFC)
    z_dir = np.sign(lfc) * z

    # magnitude tempering: asinh keeps large LFCs from dominating while staying ~linear near 0
    mag = np.arcsinh(np.abs(lfc))

    # final score: directional Z scaled by tempered magnitude
    score = z_dir * mag

    if positive_only:
        score = np.maximum(0.0, score)  # keep only enrichment direction

    out = df.copy()
    out["combined_score"] = score
    out["z_dir"] = z_dir
    out["mag_asinh_lfc"] = mag
    return out


# ============================================================
# 2) Permutation testing helpers
# ============================================================
def _perm_pvalues_both(scores_all: np.ndarray,
                       n_set: int,
                       observed_sum: float,
                       *,
                       n_perm: int = 10000,
                       random_state: int | None = 42):
    """
    Monte-Carlo permutation p-values for the sum of scores.
    Returns (p_right, p_left, mean_null, sd_null).
    """
    U = scores_all.size
    if n_set > U:
        raise ValueError(f"Set size ({n_set}) exceeds universe ({U}).")
    rng = np.random.default_rng(random_state)

    pop_mean = scores_all.mean()
    pop_var  = scores_all.var(ddof=0)
    mean_null = n_set * pop_mean
    sd_null   = float(np.sqrt(n_set * (U - n_set) / max(U - 1, 1) * pop_var))

    ge = le = 0
    for _ in range(int(n_perm)):
        idx = rng.choice(U, size=n_set, replace=False)
        s = float(scores_all[idx].sum())
        if s >= observed_sum: ge += 1
        if s <= observed_sum: le += 1

    # add-one smoothing
    p_right = (ge + 1) / (n_perm + 1)
    p_left  = (le + 1) / (n_perm + 1)
    return float(p_right), float(p_left), float(mean_null), float(sd_null)


# ============================================================
# 3) Tissue enrichment analysis
# ============================================================
def net_tissue_enrichment_signed_two_tailed(
    inc_df: pd.DataFrame,
    dec_df: pd.DataFrame,
    tau_df: pd.DataFrame,
    *,
    id_col_hits="Transcript stable ID",
    score_col="combined_score",
    id_col_tau="Transcript stable ID",
    tissue_col="primary_tissue",
    min_set_size: int = 10,
    n_perm: int = 10000,
    random_state: int | None = 1,
):
    # per-transcript signed scores: +inc − dec
    def prep(df):
        d = df[[id_col_hits, score_col]].copy()
        d["ENST"] = d[id_col_hits].astype(str).str.replace(r"\.\d+$", "", regex=True)
        return d.groupby("ENST", as_index=False)[score_col].sum()

    inc = prep(inc_df).rename(columns={score_col: "inc"})
    dec = prep(dec_df).rename(columns={score_col: "dec"})
    uni = pd.merge(inc, dec, on="ENST", how="outer").fillna(0.0)
    uni["signed_score"] = uni["inc"] - uni["dec"]

    # map to tissues
    tau = tau_df[[id_col_tau, tissue_col]].copy()
    tau["ENST"] = tau[id_col_tau].astype(str).str.replace(r"\.\d+$", "", regex=True)

    # universe
    universe_ids = np.intersect1d(uni["ENST"].unique(), tau["ENST"].unique())
    if len(universe_ids) == 0:
        raise ValueError("No overlap between signed-score universe and tau mapping.")
    uni_u = uni.loc[uni["ENST"].isin(universe_ids), ["ENST", "signed_score"]]
    tau_u = tau.loc[tau["ENST"].isin(universe_ids), ["ENST", tissue_col]]

    scores_all = uni_u["signed_score"].to_numpy()
    total_sum = float(scores_all.sum())
    U = len(uni_u)

    rows = []
    for tissue, g in tau_u.groupby(tissue_col, sort=False):
        members = g["ENST"].unique()
        n_set = len(members)
        n_bg  = U - n_set
        if n_set < min_set_size or n_bg <= 0:
            continue

        s_in  = uni_u.loc[uni_u["ENST"].isin(members), "signed_score"].to_numpy()
        s_out = uni_u.loc[~uni_u["ENST"].isin(members), "signed_score"].to_numpy()

        # permutation on sum (both tails)
        obs_sum = float(s_in.sum())
        exp_sum = total_sum * (n_set / U) if U > 0 else np.nan
        p_right, p_left, mean_null, sd_null = _perm_pvalues_both(
            scores_all, n_set, obs_sum, n_perm=n_perm,
            random_state=None if random_state is None else (random_state + hash(tissue) % 10_000_000)
        )
        z_perm = (obs_sum - mean_null) / sd_null if sd_null > 0 else np.nan
        p_two_perm = min(1.0, 2 * min(p_right, p_left))  # two-sided

        # MWU on signed scores (both tails)
        try:
            U_greater = mannwhitneyu(s_in, s_out, alternative="greater")
            p_mw_greater = float(U_greater.pvalue)
            U_less = mannwhitneyu(s_in, s_out, alternative="less")
            p_mw_less = float(U_less.pvalue)
            p_mw_two  = min(1.0, 2 * min(p_mw_greater, p_mw_less))
            auc = float(U_greater.statistic) / (n_set * n_bg)  # AUC; <0.5 suggests depletion
        except ValueError:
            p_mw_greater = p_mw_less = p_mw_two = 1.0
            auc = np.nan

        # directional one-sided p: pick tail matching the effect sign
        sign = np.sign(obs_sum - (mean_null if np.isfinite(mean_null) else 0.0))
        p_dir_perm = p_right if sign >= 0 else p_left
        p_dir_mw   = p_mw_greater if sign >= 0 else p_mw_less

        rows.append({
            "tissue": tissue,
            "n_set": n_set,
            "universe": U,
            "obs_sum_signed": obs_sum,
            "expected_sum_signed": exp_sum,
            "fold_signed": (obs_sum / exp_sum) if exp_sum != 0 else np.nan,
            "z_perm_signed": z_perm,
            "p_perm_right": p_right,
            "p_perm_left": p_left,
            "p_perm_two": p_two_perm,
            "p_dir_perm": p_dir_perm,        # use this for signed combined scores
            "p_mw_greater": p_mw_greater,
            "p_mw_less": p_mw_less,
            "p_mw_two": p_mw_two,
            "p_dir_mw": p_dir_mw,            # use this for signed combined scores (MWU)
            "AUC_signed": auc,
            "median_in_signed": float(np.median(s_in)) if len(s_in) else np.nan,
            "median_out_signed": float(np.median(s_out)) if len(s_out) else np.nan,
        })

    res = pd.DataFrame(rows)
    # FDR corrections
    res["FDR_perm_two"] = multipletests(res["p_perm_two"].values, method="fdr_bh")[1]
    res["FDR_mw_two"]   = multipletests(res["p_mw_two"].values,   method="fdr_bh")[1]
    res["FDR_dir_perm"] = multipletests(res["p_dir_perm"].values, method="fdr_bh")[1]
    res["FDR_dir_mw"]   = multipletests(res["p_dir_mw"].values,   method="fdr_bh")[1]

    # sort by directional permutation FDR, then |fold| and |z|
    res = res.sort_values(["FDR_dir_perm", "fold_signed", "z_perm_signed"],
                          ascending=[True, False, False]).reset_index(drop=True)
    return res


# ============================================================
# 4) Plot annotation helpers
# ============================================================
def fdr_to_stars(fdr: float) -> str:
    """
    Convert FDR to significance stars.
    """
    if np.isnan(fdr):
        return ""
    if fdr < 1e-4: return "****"
    if fdr < 1e-3: return "***"
    if fdr < 1e-2: return "**"
    if fdr < 5e-2: return "*"
    return ""


# ============================================================
# 5) Tissue enrichment plotting
# ============================================================
def plot_tissue_enrichment(
    res: pd.DataFrame,
    *,
    tissue_col: str = "tissue",
    score_col: str = "combined_score",
    fdr_col: str = "FDR_mw",
    cmap: str = "viridis",          # colormap for FDR
    max_fdr_shown: float = 0.10,    # clamp color scale at this FDR
    figsize_per_bar: float = 0.35,  # height per bar (inches)
    title: str = "Tissue Enrichment by Combined Score",
    savepath: str | None = None
):
    """
    Horizontal bar plot: tissues on y, combined_score on x, colored by FDR (viridis).
    Adds significance asterisks to the right of positive bars AND to the left of negative bars.
    """
    df = res[[tissue_col, score_col, fdr_col]].dropna(subset=[score_col]).copy()

    # Order by score descending (so big positives at top)
    df = df.sort_values(score_col, ascending=False).reset_index(drop=True)

    # Colors from FDR (lower FDR -> darker color)
    fdr_vals = df[fdr_col].fillna(1.0).astype(float).to_numpy()
    fdr_clamped = np.clip(fdr_vals, 0.0, max_fdr_shown)
    norm = Normalize(vmin=0.0, vmax=max_fdr_shown)
    cmap_obj = cm.get_cmap(cmap)
    colors = cmap_obj(norm(fdr_clamped))

    # Figure sizing
    n = len(df)
    height = max(2.5, n * figsize_per_bar)
    fig, ax = plt.subplots(figsize=(8, height))

    # Bars
    y = np.arange(n)
    x = df[score_col].to_numpy()
    ax.barh(y, x, color=colors, edgecolor="black", linewidth=0.5)

    # Y labels (tissues)
    ax.set_yticks(y)
    ax.set_yticklabels(df[tissue_col].tolist())
    ax.invert_yaxis()  # largest at top

    # Labels & title
    ax.set_xlabel("Combined Score (Fold Change & q-value)", fontsize=16)
    ax.set_ylabel("Tissue (GTEx)", fontsize=16)
    ax.set_title(title, fontsize=16, pad=25)

    # Clean up spines, add 0-line for reference
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axvline(0, lw=0.8, color="k", alpha=0.3)

    # Star placement offsets & x-limits with margins so stars don't clip
    x_abs_max = float(np.nanmax(np.abs(x))) if len(x) else 1.0
    x_offset = 0.02 * x_abs_max    # star offset from bar end
    margin   = 0.14 * x_abs_max    # axis margin to fit stars

    left_lim  = min(0.0, float(np.nanmin(x))) - margin
    right_lim = max(0.0, float(np.nanmax(x))) + margin
    ax.set_xlim(left_lim, right_lim)

    # Asterisks: right of positive bars; left of negative bars; near +0 for zeros
    for i, (xi, fdr) in enumerate(zip(x, fdr_vals)):
        stars = fdr_to_stars(fdr)
        if not stars or not np.isfinite(xi):
            continue
        if xi > 0:
            ax.text(xi + x_offset, i, stars, va="center", ha="left", fontsize=11)
        elif xi < 0:
            ax.text(xi - x_offset, i, stars, va="center", ha="right", fontsize=11)
        else:
            # xi == 0: place to the right of zero
            ax.text(x_offset, i, stars, va="center", ha="left", fontsize=11)

    # Colorbar for FDR
    sm = cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(f"q-value", fontsize=14, rotation=270, labelpad=15)
    cbar.ax.invert_yaxis()  # optional: make 0 at top

    plt.tight_layout()
    if savepath:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")
    plt.show()
