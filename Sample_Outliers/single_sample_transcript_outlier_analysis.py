"""
Single-sample transcript outlier analysis utilities.
"""

# ============================================================
# Imports
# ============================================================
import numpy as np
import pandas as pd

from scipy import optimize
from scipy.special import digamma, polygamma
from scipy import stats as sps


# ============================================================
# 1) Multiple-testing correction
# ============================================================
def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR for a 1D array."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = np.arange(1, n + 1)
    q = p[order] * n / ranked
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


# ============================================================
# 2) Empirical Bayes helper
# ============================================================
def inv_trigamma(y: float) -> float:
    """
    Invert trigamma(x) = y for y>0. Returns x>0.
    """
    if y <= 0:
        raise ValueError("y must be > 0 for inv_trigamma")

    f = lambda x: polygamma(1, x) - y

    low = 1e-8
    high = max(1.0 / y, 1e-6)

    while f(high) > 0 and high < 1e8:
        high *= 2.0

    root = optimize.brentq(f, low, high, maxiter=200, xtol=1e-12)
    return root


# ============================================================
# 3) Single-sample outlier analysis
# ============================================================
def eb_robust_outliers_for_sample(
    data: pd.DataFrame,
    sample_col: str,
    id_col: str = "transcript_id",
    tpm_cutoff: float = 0.1,
    min_prop: float = 0.05,
    target_tpm_keep: float = 0.1,
    pseudo: float = 1.0,
    robust_z_thresh: float = 4.0,
    pct_hi_thresh: float = 0.99,
    pct_lo_thresh: float = 0.01,
    q_thresh: float = 0.05,
    min_abs_diff_log2: float = 1.0,
) -> pd.DataFrame:
    """
    Leave-one-out, robust + EB single-sample isoform outlier detection.

    For a target sample:
      1) Filter isoforms:
         keep if:
           - target sample TPM >= target_tpm_keep
           OR
           - TPM >= tpm_cutoff in at least ceil(min_prop * n_samples) samples
      2) Transform Y = log2(TPM + pseudo)
      3) For each isoform, compute leave-one-out baseline from OTHER samples:
           - mean_others
           - median_others
           - MAD_others
           - percentile of target vs others
      4) Compute robust z-score:
           z_robust = (target - median_others) / (1.4826 * MAD_others)
      5) Compute EB-moderated residual score:
           residual = target - mean_others
           shrink variance of OTHER samples across isoforms
      6) Compute two-sided p-values, BH FDR
      7) Return metrics and boolean outlier flags

    Notes:
      - q-values here are useful for ranking, but still should be interpreted
        as an outlier screen rather than a perfect inferential model.
      - robust_z / percentile are often the most biologically interpretable columns.
    """
    assert id_col in data.columns, f"{id_col=} not found"
    assert sample_col in data.columns, f"{sample_col=} not found"

    sample_cols = [c for c in data.columns if c != id_col]
    n = len(sample_cols)
    if n < 4:
        raise ValueError("Need at least 4 samples total for leave-one-out analysis.")

    j = sample_cols.index(sample_col)

    # TPM matrix
    tpm = data[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float64)
    ids = data[id_col].to_numpy()

    # Smarter filter:
    # keep if expressed in target sample OR expressed in enough cohort samples
    min_count = int(np.ceil(min_prop * n))
    keep = (tpm[:, j] >= target_tpm_keep) | ((tpm >= tpm_cutoff).sum(axis=1) >= min_count)

    if not np.any(keep):
        raise ValueError("No isoforms passed the expression filter; relax thresholds.")

    print(
        f"Keeping {keep.sum()} isoforms: "
        f"target TPM >= {target_tpm_keep} OR TPM >= {tpm_cutoff} in at least {min_count} samples."
    )

    tpm = tpm[keep, :]
    ids = ids[keep]

    # Log transform
    Y = np.log2(tpm + float(pseudo))
    m = Y.shape[0]

    # Split target vs others
    y_target = Y[:, j]
    Y_others = np.delete(Y, j, axis=1)  # shape (m, n-1)
    n_others = Y_others.shape[1]
    df_resid = n_others - 1

    if n_others < 3:
        raise ValueError("Need at least 3 background samples after leaving one out.")

    # Leave-one-out summary statistics
    mean_others = Y_others.mean(axis=1)
    median_others = np.median(Y_others, axis=1)
    mad_others = sps.median_abs_deviation(Y_others, axis=1, scale=1.0, nan_policy="omit")
    robust_scale = np.maximum(1e-6, 1.4826 * mad_others)

    diff_log2_mean = y_target - mean_others
    diff_log2_median = y_target - median_others

    # Robust z-score vs others
    z_robust = diff_log2_median / robust_scale

    # Percentiles: fraction of others <= target
    pct = np.mean(Y_others <= y_target[:, None], axis=1)

    # Direction
    direction = np.where(diff_log2_median >= 0, "higher_in_sample", "lower_in_sample")

    # EB variance on OTHER samples only
    R_other = Y_others - mean_others[:, None]
    s2 = np.sum(R_other**2, axis=1) / df_resid

    # limma-style shrinkage on log variances
    a = df_resid / 2.0
    log_s2 = np.log(s2 + 1e-300)
    mean_log_s2 = log_s2.mean()
    var_log_s2 = log_s2.var(ddof=1)

    s02 = np.exp(mean_log_s2 - (digamma(a) - np.log(a)))

    target_var = var_log_s2 - polygamma(1, a)
    if target_var <= 1e-12:
        df0 = np.inf
    else:
        a0 = inv_trigamma(target_var)
        df0 = 2.0 * a0

    if np.isfinite(df0):
        s2_post = (df0 * s02 + df_resid * s2) / (df0 + df_resid)
        df_post = df0 + df_resid
    else:
        s2_post = np.full_like(s2, s02)
        df_post = np.inf

    # EB moderated residual for target sample vs leave-one-out mean
    # predictive-ish SE: sqrt(s2_post * (1 + 1/n_others))
    se_target = np.sqrt(s2_post * (1.0 + 1.0 / n_others))
    t_mod = diff_log2_mean / np.maximum(se_target, 1e-12)

    # Raw p-values
    if np.isfinite(df_post):
        p_raw = 2.0 * sps.t.sf(np.abs(t_mod), df=df_post)
    else:
        p_raw = 2.0 * sps.norm.sf(np.abs(t_mod))

    # Heavy-tail correction
    kurt_excess = sps.kurtosis(t_mod, fisher=True, bias=False, nan_policy="omit")
    if np.isfinite(df_post) and np.isfinite(kurt_excess) and kurt_excess > 0:
        nu_hat = 4.0 + 6.0 / kurt_excess
        df_eff = min(df_post, nu_hat)
    else:
        df_eff = df_post

    if np.isfinite(df_eff):
        p_ttail = 2.0 * sps.t.sf(np.abs(t_mod), df=df_eff)
    else:
        p_ttail = 2.0 * sps.norm.sf(np.abs(t_mod))

    q_raw = bh_fdr(p_raw)
    q_ttail = bh_fdr(p_ttail)

    # Back to TPM for interpretability
    tpm_target = tpm[:, j]
    tpm_others = np.delete(tpm, j, axis=1)
    mean_tpm_others = tpm_others.mean(axis=1)
    median_tpm_others = np.median(tpm_others, axis=1)

    # Boolean calls
    robust_high_flag = (z_robust >= robust_z_thresh) | (pct >= pct_hi_thresh)
    robust_low_flag  = (z_robust <= -robust_z_thresh) | (pct <= pct_lo_thresh)

    eb_flag = (q_raw <= q_thresh) & (np.abs(diff_log2_median) >= min_abs_diff_log2)

    final_high_flag = robust_high_flag & eb_flag & (diff_log2_median > 0)
    final_low_flag  = robust_low_flag & eb_flag & (diff_log2_median < 0)
    final_any_flag  = final_high_flag | final_low_flag

    out = pd.DataFrame({
        "transcript_id": ids,
        "sample": sample_col,

        "TPM_sample": tpm_target,
        "mean_TPM_others": mean_tpm_others,
        "median_TPM_others": median_tpm_others,

        "log2TPM_sample": y_target,
        "mean_log2TPM_others": mean_others,
        "median_log2TPM_others": median_others,

        "diff_log2_mean": diff_log2_mean,
        "diff_log2_median": diff_log2_median,

        "mad_others": mad_others,
        "robust_scale": robust_scale,
        "z_robust": z_robust,
        "percentile_vs_others": pct,

        "t_mod": t_mod,
        "p_raw_dfpost": p_raw,
        "q_raw_dfpost": q_raw,
        "p_ttail": p_ttail,
        "q_ttail": q_ttail,
        "df_post": df_post if np.isfinite(df_post) else np.inf,
        "s2_post": s2_post,

        "direction": direction,
        "robust_high_flag": robust_high_flag,
        "robust_low_flag": robust_low_flag,
        "eb_flag": eb_flag,
        "final_high_flag": final_high_flag,
        "final_low_flag": final_low_flag,
        "final_any_flag": final_any_flag,
    })

    out = out.sort_values(
        ["final_any_flag", "q_ttail", "p_ttail", "z_robust"],
        ascending=[False, True, True, False]
    ).reset_index(drop=True)

    return out
