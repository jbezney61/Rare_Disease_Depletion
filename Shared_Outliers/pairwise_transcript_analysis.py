"""
Pairwise transcript outlier analysis utilities.
"""

# ============================================================
# Imports
# ============================================================
import numpy as np
import pandas as pd

from itertools import combinations
from joblib import Parallel, delayed
from scipy import optimize
from scipy.special import digamma, polygamma
from scipy import stats as sps


# ============================================================
# 1) Preprocessing
# ============================================================
def preprocess_condition_matrix(
    data: pd.DataFrame,
    id_col: str = "transcript_id",
    pseudo: float = 1.0,
    global_tpm_cutoff: float = 0.1,
    global_min_prop: float = 0.02,
):
    sample_cols = [c for c in data.columns if c != id_col]
    ids = data[id_col].to_numpy()

    tpm = (
        data[sample_cols]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
        .to_numpy(dtype=np.float64)
    )

    n = tpm.shape[1]
    min_count = int(np.ceil(global_min_prop * n))
    keep = (tpm >= global_tpm_cutoff).sum(axis=1) >= min_count

    tpm = tpm[keep]
    ids = ids[keep]
    Y = np.log2(tpm + pseudo)

    return {
        "ids": ids,
        "sample_cols": sample_cols,
        "tpm": tpm,
        "Y": Y,
    }


# ============================================================
# 2) Empirical Bayes helper functions
# ============================================================
def inv_trigamma(y: float) -> float:
    f = lambda x: polygamma(1, x) - y
    low = 1e-8
    high = max(1.0 / y, 1e-6)
    while f(high) > 0 and high < 1e8:
        high *= 2.0
    return optimize.brentq(f, low, high, maxiter=200, xtol=1e-12)


def estimate_global_eb_prior(Y: np.ndarray):
    """
    Estimate EB prior from all samples together once.
    """
    mean_all = Y.mean(axis=1, keepdims=True)
    R = Y - mean_all
    df_resid = Y.shape[1] - 1
    s2 = np.sum(R**2, axis=1) / df_resid

    a = df_resid / 2.0
    log_s2 = np.log(s2 + 1e-300)
    mean_log_s2 = log_s2.mean()
    var_log_s2 = log_s2.var(ddof=1)

    s02 = np.exp(mean_log_s2 - (digamma(a) - np.log(a)))

    target_var = var_log_s2 - polygamma(1, a)
    if target_var <= 1e-12:
        df0 = np.inf
    else:
        df0 = 2.0 * inv_trigamma(target_var)

    return s02, df0


def bh_fdr(pvals):
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
# 3) Pairwise worker
# ============================================================
def fast_pair_worker(
    j1,
    j2,
    ids,
    sample_cols,
    tpm,
    Y,
    s02,
    df0,
    robust_z_thresh=4.0,
    pct_hi_thresh=0.99,
    pct_lo_thresh=0.01,
    q_thresh=0.05,
    min_abs_diff_log2=1.0,
    use_q="raw",
    keep_only_significant=False,
):
    n = Y.shape[1]
    other_idx = np.array([k for k in range(n) if k not in (j1, j2)])

    Y_target = Y[:, [j1, j2]]
    Y_others = Y[:, other_idx]
    tpm_target = tpm[:, [j1, j2]]
    tpm_others = tpm[:, other_idx]

    n_other = Y_others.shape[1]

    mean_others = Y_others.mean(axis=1)
    median_others = np.median(Y_others, axis=1)
    mad_others = sps.median_abs_deviation(Y_others, axis=1, scale=1.0, nan_policy="omit")
    robust_scale = np.maximum(1e-6, 1.4826 * mad_others)

    y1 = Y_target[:, 0]
    y2 = Y_target[:, 1]
    z1 = (y1 - median_others) / robust_scale
    z2 = (y2 - median_others) / robust_scale
    pct1 = np.mean(Y_others <= y1[:, None], axis=1)
    pct2 = np.mean(Y_others <= y2[:, None], axis=1)

    y_pair_mean = Y_target.mean(axis=1)
    diff_log2_mean = y_pair_mean - mean_others
    diff_log2_median = y_pair_mean - median_others

    # reuse global EB prior
    R_other = Y_others - mean_others[:, None]
    df_resid = n_other - 1
    s2 = np.sum(R_other**2, axis=1) / df_resid

    if np.isfinite(df0):
        s2_post = (df0 * s02 + df_resid * s2) / (df0 + df_resid)
        df_post = df0 + df_resid
    else:
        s2_post = np.full_like(s2, s02)
        df_post = np.inf

    se_pair = np.sqrt(s2_post * (1/2 + 1/n_other))
    t_mod = diff_log2_mean / np.maximum(se_pair, 1e-12)

    if np.isfinite(df_post):
        p_raw = 2.0 * sps.t.sf(np.abs(t_mod), df=df_post)
    else:
        p_raw = 2.0 * sps.norm.sf(np.abs(t_mod))

    q_raw = bh_fdr(p_raw)

    robust_high_flag = ((z1 >= robust_z_thresh) & (z2 >= robust_z_thresh)) | (
        (pct1 >= pct_hi_thresh) & (pct2 >= pct_hi_thresh)
    )
    robust_low_flag = ((z1 <= -robust_z_thresh) & (z2 <= -robust_z_thresh)) | (
        (pct1 <= pct_lo_thresh) & (pct2 <= pct_lo_thresh)
    )

    q_use = q_raw
    eb_flag = (q_use <= q_thresh) & (np.abs(diff_log2_median) >= min_abs_diff_log2)

    final_high_flag = robust_high_flag & eb_flag & (diff_log2_median > 0)
    final_low_flag = robust_low_flag & eb_flag & (diff_log2_median < 0)
    final_any_flag = final_high_flag | final_low_flag

    out = pd.DataFrame({
        "transcript_id": ids,
        "pair_label": "__".join(sorted([
            sample_cols[j1].replace("-WM", "").replace("-D", ""),
            sample_cols[j2].replace("-WM", "").replace("-D", "")
        ])),
        "sample1": sample_cols[j1],
        "sample2": sample_cols[j2],
        "z1": z1,
        "z2": z2,
        "diff_log2_median": diff_log2_median,
        "q_raw_dfpost": q_raw,
        "final_high_flag": final_high_flag,
        "final_low_flag": final_low_flag,
        "final_any_flag": final_any_flag,
    })

    if keep_only_significant:
        out = out[out["final_any_flag"]].copy()

    return out


# ============================================================
# 4) Run all sample pairs
# ============================================================
def run_all_pairs_fast(
    prep,
    condition_name,
    n_jobs=8,
    keep_only_significant=False,
    robust_z_thresh=4.0,
    pct_hi_thresh=0.99,
    pct_lo_thresh=0.01,
    q_thresh=0.05,
    min_abs_diff_log2=1.0,
):
    ids = prep["ids"]
    sample_cols = prep["sample_cols"]
    tpm = prep["tpm"]
    Y = prep["Y"]

    s02, df0 = estimate_global_eb_prior(Y)

    pair_idx = list(combinations(range(len(sample_cols)), 2))

    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(fast_pair_worker)(
            j1, j2,
            ids=ids,
            sample_cols=sample_cols,
            tpm=tpm,
            Y=Y,
            s02=s02,
            df0=df0,
            robust_z_thresh=robust_z_thresh,
            pct_hi_thresh=pct_hi_thresh,
            pct_lo_thresh=pct_lo_thresh,
            q_thresh=q_thresh,
            min_abs_diff_log2=min_abs_diff_log2,
            keep_only_significant=keep_only_significant,
        )
        for j1, j2 in pair_idx
    )

    out = pd.concat(results, ignore_index=True)
    out["condition"] = condition_name
    return out


# ============================================================
# 5) Build pair × transcript matrix
# ============================================================
def build_pair_transcript_matrix(
    pair_results: pd.DataFrame,
    transcript_col: str = "transcript_id",
    pair_col: str = "pair_label",
    min_pairs: int = 2,
    score_mode: str = "avg_z",   # "avg_z", "diff_log2_median", "signed_binary"
) -> pd.DataFrame:
    """
    Build a pair × transcript matrix from long pair results.

    Entries:
      - avg_z: mean(z1, z2) if final_any_flag else 0
      - diff_log2_median: pair diff_log2_median if final_any_flag else 0
      - signed_binary: +1 for final_high_flag, -1 for final_low_flag, else 0

    Filters transcripts to those with final_any_flag == True in at least min_pairs pairs.
    """
    req = {
        transcript_col, pair_col, "final_any_flag",
        "final_high_flag", "final_low_flag",
        "z1", "z2", "diff_log2_median"
    }
    missing = req - set(pair_results.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = pair_results.copy()

    recurrent = (
        df.loc[df["final_any_flag"]]
        .groupby(transcript_col)[pair_col]
        .nunique()
    )
    keep_transcripts = recurrent[recurrent >= min_pairs].index
    df = df[df[transcript_col].isin(keep_transcripts)].copy()

    if score_mode == "avg_z":
        df["signed_score"] = np.where(
            df["final_any_flag"],
            (df["z1"] + df["z2"]) / 2.0,
            0.0
        )
    elif score_mode == "diff_log2_median":
        df["signed_score"] = np.where(
            df["final_any_flag"],
            df["diff_log2_median"],
            0.0
        )
    elif score_mode == "signed_binary":
        df["signed_score"] = np.select(
            [df["final_high_flag"], df["final_low_flag"]],
            [1.0, -1.0],
            default=0.0
        )
    else:
        raise ValueError("score_mode must be one of: avg_z, diff_log2_median, signed_binary")

    # Keep one row per pair/transcript; if duplicates somehow exist, keep most extreme
    df = (
        df.assign(abs_score=df["signed_score"].abs())
          .sort_values("abs_score", ascending=False)
          .drop_duplicates(subset=[pair_col, transcript_col], keep="first")
          .drop(columns="abs_score")
    )

    mat = (
        df.pivot(index=pair_col, columns=transcript_col, values="signed_score")
          .fillna(0.0)
          .sort_index()
    )
    return mat


# ============================================================
# 6) Save matrix
# ============================================================
def save_pair_matrix(mat: pd.DataFrame, out_tsv: str) -> None:
    mat.to_csv(out_tsv, sep="\t")
