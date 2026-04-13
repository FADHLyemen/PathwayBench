#!/usr/bin/env python3
"""
PathwayBench Advisor - Streamlit App
Recommends the best pathway scoring method for user data based on
6 robustness criteria from the PathwayBench benchmark.
"""

from __future__ import annotations

import io
import os
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy import stats

APP_VERSION = "2.1.0"
ZENODO_DOI = "10.5281/zenodo.19503595"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

APP_DIR = Path(__file__).resolve().parent
EXAMPLE_DIR = APP_DIR / "example_data"

METHOD_COLORS: dict[str, str] = {
    "ssGSEA": "#2166AC",
    "GSVA": "#4393C3",
    "zscore": "#D6604D",
    "AUCell": "#4DAF4A",
    "UCell": "#984EA3",
}

CRITERIA_LABELS: dict[str, str] = {
    "bio_relevance": "Biological\nrelevance",
    "calc_method_cor": "Aggregation\nstability",
    "outlier_robustness": "Outlier\nrobustness",
    "norm_stability": "Normalization\nstability",
    "sample_stability": "Sample-size\nstability",
}

# Short labels for chart axis (no newlines)
CRITERIA_SHORT: dict[str, str] = {
    "bio_relevance": "Biological relevance",
    "calc_method_cor": "Aggregation stability",
    "outlier_robustness": "Outlier robustness",
    "norm_stability": "Normalization stability",
    "sample_stability": "Sample-size stability",
}

DEFAULT_WEIGHTS: dict[str, float] = {
    "bio_relevance": 0.30,
    "calc_method_cor": 0.15,
    "outlier_robustness": 0.20,
    "norm_stability": 0.20,
    "sample_stability": 0.15,
}

# Benchmark reference — v2-corrected, averaged across 7 CellxGene datasets
# (excluding ckd_kidney pending KPMP approval).
# outlier_sensitivity is inverted to outlier_robustness (1 - value).
# Source: validation_codex/v2_corrected_results.txt (7-dataset table)
_BENCH_RAW = pd.DataFrame(
    {
        "method": ["ssGSEA", "GSVA", "zscore", "AUCell", "UCell"],
        "bio_relevance": [0.643, 0.609, 0.698, 0.570, 0.587],
        "calc_method_cor": [0.864, 0.728, 0.820, 0.840, 0.958],
        "outlier_robustness": [1 - 0.213, 1 - 0.177, 1 - 0.220, 1 - 0.175, 1 - 0.164],
        "norm_stability": [0.723, 0.849, 0.765, 0.670, 0.798],
        "sample_stability": [0.855, 0.886, 0.860, 0.857, 0.856],
    }
)

CRITERIA = list(CRITERIA_SHORT.keys())


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def scale_01(s: pd.Series) -> pd.Series:
    """Min-max scale a Series to [0, 1]."""
    rng = s.max() - s.min()
    return (s - s.min()) / rng if rng > 0 else pd.Series(0.5, index=s.index)


def load_example(name: str) -> bytes:
    return (EXAMPLE_DIR / name).read_bytes()


def parse_expression(uploaded) -> pd.DataFrame | None:
    """Read an uploaded CSV into a genes-x-samples DataFrame."""
    try:
        df = pd.read_csv(uploaded, index_col=0)
        # Ensure numeric
        df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
        return df
    except Exception as exc:
        st.error(f"Could not parse expression file: {exc}")
        return None


def parse_gene_set(uploaded=None, pasted: str = "") -> list[str]:
    """Return a deduplicated gene list from file or pasted text."""
    genes: list[str] = []
    if uploaded is not None:
        text = uploaded.read().decode("utf-8", errors="ignore")
        genes = [g.strip() for g in text.replace(",", "\n").splitlines() if g.strip()]
    elif pasted.strip():
        genes = [g.strip() for g in pasted.replace(",", "\n").splitlines() if g.strip()]
    return list(dict.fromkeys(genes))  # dedupe preserving order


# ---------------------------------------------------------------------------
# Data profiling
# ---------------------------------------------------------------------------


def profile_data(expr: pd.DataFrame) -> dict:
    """Compute summary statistics for the uploaded expression matrix."""
    info: dict = {}
    info["n_genes"] = expr.shape[0]
    info["n_samples"] = expr.shape[1]
    info["gene_example"] = ", ".join(expr.index[:5].tolist())

    # Normalization heuristic
    sample = expr.iloc[: min(500, expr.shape[0]), : min(10, expr.shape[1])].values
    has_neg = bool(np.nanmin(sample) < 0)
    max_val = float(np.nanmax(sample))
    if has_neg:
        info["norm_guess"] = "Likely voom / sctransform (negative values present)"
    elif max_val < 30:
        info["norm_guess"] = "Likely log2-CPM or log-normalized"
    elif max_val > 1e4:
        info["norm_guess"] = "Likely raw counts (large values)"
    else:
        info["norm_guess"] = "Normalized (CPM or similar)"

    # Library-size outlier detection (IQR)
    lib = expr.sum(axis=0)
    cv = float(lib.std() / lib.mean()) if lib.mean() > 0 else 0.0
    q1, q3 = float(lib.quantile(0.25)), float(lib.quantile(0.75))
    iqr = q3 - q1
    n_out = int(((lib < q1 - 1.5 * iqr) | (lib > q3 + 1.5 * iqr)).sum())
    info["lib_size_cv"] = round(cv, 3)
    info["n_outliers"] = n_out
    info["outlier_pct"] = round(100 * n_out / info["n_samples"], 1)

    # Try to infer condition from column names
    conditions = []
    for col in expr.columns:
        low = col.lower()
        if "disease" in low:
            conditions.append("disease")
        elif "control" in low:
            conditions.append("control")
        else:
            conditions.append("unknown")
    cond_counts = pd.Series(conditions).value_counts().to_dict()
    if "unknown" in cond_counts and len(cond_counts) == 1:
        cond_counts = {}
    info["conditions"] = cond_counts

    return info


# ---------------------------------------------------------------------------
# Z-score scoring (pure numpy — no R required)
# ---------------------------------------------------------------------------


def score_zscore(expr: pd.DataFrame, gene_set: list[str]) -> pd.Series:
    """
    Compute Z-score pathway activity per sample.
    For each gene in the set: z-score across samples, then average.
    """
    overlap = [g for g in gene_set if g in expr.index]
    if len(overlap) < 5:
        return pd.Series(dtype=float)
    sub = expr.loc[overlap].values.astype(float)  # genes x samples
    mu = sub.mean(axis=1, keepdims=True)
    sd = sub.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    z = (sub - mu) / sd
    scores = z.mean(axis=0)
    return pd.Series(scores, index=expr.columns, name="zscore")


# ---------------------------------------------------------------------------
# Recommendation engine
# ---------------------------------------------------------------------------


def compute_recommendation(
    profile: dict,
    zscore_scores: pd.Series | None,
    benchmark: pd.DataFrame = _BENCH_RAW,
) -> dict:
    """
    Weight benchmark criteria, adjust for user-data context,
    and return the recommendation.
    """
    methods = benchmark["method"].tolist()
    metrics = benchmark.copy()

    # Adjust weights based on user data
    weights = dict(DEFAULT_WEIGHTS)
    reasons: list[str] = []

    if profile["n_samples"] < 15:
        weights["sample_stability"] *= 2
        reasons.append(
            f"dataset is small ({profile['n_samples']} samples) "
            "-- sample-size stability weighted higher"
        )
    if profile["outlier_pct"] > 20:
        weights["outlier_robustness"] *= 3
        reasons.append(
            f"{profile['outlier_pct']}% outlier samples detected "
            "-- outlier robustness weighted higher"
        )

    # Normalize weights
    total = sum(weights.values())
    weights = {k: v / total for k, v in weights.items()}

    # Scale each criterion 0-1 across methods
    scaled = metrics[["method"]].copy()
    for c in CRITERIA:
        scaled[c] = scale_01(metrics[c]).values

    # Composite score
    scaled["composite"] = sum(weights[c] * scaled[c] for c in CRITERIA)

    best_idx = int(scaled["composite"].idxmax())
    best_method = scaled.loc[best_idx, "method"]

    # Build explanation
    row = scaled.loc[best_idx]
    strengths = []
    if row["bio_relevance"] > 0.6:
        strengths.append("strong biological relevance")
    if row["calc_method_cor"] > 0.6:
        strengths.append("stable across aggregation methods")
    if row["outlier_robustness"] > 0.6:
        strengths.append("robust to outlier samples")
    if row["norm_stability"] > 0.6:
        strengths.append("consistent across normalizations")
    if row["sample_stability"] > 0.6:
        strengths.append("stable at reduced sample sizes")

    explanation = (
        f"**{best_method}** scored highest overall "
        f"(composite = {row['composite']:.2f})."
    )
    if strengths:
        explanation += f" Key strengths: {'; '.join(strengths)}."
    if reasons:
        explanation += f" Context adjustments: {'; '.join(reasons)}."

    return {
        "metrics": metrics,
        "scaled": scaled,
        "weights": weights,
        "best_method": best_method,
        "best_score": float(row["composite"]),
        "explanation": explanation,
    }


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Priority-based recommendation (matches Figure 9 decision tree)
# ---------------------------------------------------------------------------


def get_priority_based_recommendation(
    priority: str,
    sub_priority_bio: str | None = None,
    sub_priority_robust: str | None = None,
    sctransform_used: bool = False,
    gene_set_size: int = 50,
    outlier_fraction: float = 0.0,
) -> tuple[str, str, list[str]]:
    """Implement the Figure 9 (Practical Guidelines) decision tree.

    Returns (recommended_method, reasoning_text, warnings_list).
    """
    method: str | None = None
    reason = ""

    if priority.startswith("Biological"):
        if sub_priority_bio and "Direction accuracy" in sub_priority_bio:
            method = "Z-score"
            reason = (
                "Z-score has the highest direction accuracy in PathwayBench "
                "(0.694 across 8 datasets, wins outright on 3 of 8)."
            )
        elif sub_priority_bio and "Effect magnitude" in sub_priority_bio:
            method = "AUCell"
            reason = (
                "AUCell produces the largest effect magnitudes "
                "(mean |Cohen's d| = 0.675), useful for ranking pathways. "
                "Also run a magnitude method (ssGSEA, GSVA, or Z-score) to "
                "detect rank-window inversion in diseases with broad "
                "transcriptional remodeling."
            )
    elif priority.startswith("Balanced"):
        method = "GSVA or UCell (run both)"
        reason = (
            "Both achieve Good ratings on 4/5 PathwayBench criteria. "
            "GSVA: stronger biology and sample-size stability. "
            "UCell: stronger aggregation and outlier robustness. "
            "Comparing the two gives a sensitivity-vs-robustness check."
        )
    elif priority.startswith("Technical"):
        if sub_priority_robust == "\u226515 donors":
            method = "UCell"
            reason = (
                "UCell leads aggregation stability (0.960) and outlier "
                "robustness (0.163 deviation). With \u226515 donors, "
                "sample-size stability is not the bottleneck."
            )
        else:
            method = "GSVA"
            reason = (
                "GSVA leads sample-size stability (0.893) and combines it "
                "with low outlier deviation (0.183), making it the safer "
                "choice with <15 donors."
            )

    warnings: list[str] = []
    if sctransform_used and method == "AUCell":
        warnings.append(
            "sctransform normalization detected: AUCell has the lowest "
            "normalization independence (0.669). Prefer UCell or GSVA."
        )
    if gene_set_size < 20 and method == "Z-score":
        warnings.append(
            "Small gene set (<20 genes): Z-score is sensitive to individual "
            "high-variance genes. Prefer GSVA or ssGSEA."
        )
    if outlier_fraction > 20.0 and method not in ("UCell", "GSVA or UCell (run both)"):
        warnings.append(
            "High outlier fraction (>20%): Consider UCell (lowest outlier "
            "deviation, 0.163) or GSVA (0.183)."
        )

    return method or "GSVA or UCell (run both)", reason, warnings


# ---------------------------------------------------------------------------
# Plotting helpers (Plotly)
# ---------------------------------------------------------------------------


def plot_criteria_bars(scaled: pd.DataFrame, best_method: str, weights: dict) -> go.Figure:
    """Grouped bar chart of criteria performance."""
    methods = scaled["method"].tolist()

    fig = go.Figure()

    for method in methods:
        row = scaled[scaled["method"] == method].iloc[0]
        vals = [float(row[c]) for c in CRITERIA]
        labels = [CRITERIA_SHORT[c] for c in CRITERIA]
        border_w = 3 if method == best_method else 0
        border_c = "black" if method == best_method else "rgba(0,0,0,0)"

        fig.add_trace(
            go.Bar(
                name=method,
                x=labels,
                y=vals,
                marker_color=METHOD_COLORS[method],
                marker_line_width=border_w,
                marker_line_color=border_c,
            )
        )

    wt_text = "  ".join(
        f"{s} {round(weights[c]*100)}%"
        for c, s in zip(CRITERIA, ["Bio", "Agg", "Out", "Norm", "Size"])
    )

    fig.update_layout(
        barmode="group",
        title=dict(text="Criteria performance (scaled 0-1)", font_size=16),
        annotations=[
            dict(
                text=f"Recommended method outlined | Weights: {wt_text}",
                xref="paper", yref="paper", x=0, y=1.06,
                showarrow=False, font=dict(size=12, color="grey"),
            )
        ],
        yaxis=dict(title="Scaled score", range=[0, 1.08]),
        xaxis=dict(title=""),
        legend=dict(title="Method"),
        height=420,
        margin=dict(t=80),
    )
    return fig


def plot_benchmark_ref(benchmark: pd.DataFrame) -> go.Figure:
    """Bar chart of the benchmark reference database."""
    # Scale 0-1 for display
    scaled = benchmark[["method"]].copy()
    for c in CRITERIA:
        scaled[c] = scale_01(benchmark[c]).values

    fig = go.Figure()
    for _, row in scaled.iterrows():
        method = row["method"]
        fig.add_trace(
            go.Bar(
                name=method,
                x=[CRITERIA_SHORT[c] for c in CRITERIA],
                y=[float(row[c]) for c in CRITERIA],
                marker_color=METHOD_COLORS[method],
            )
        )

    fig.update_layout(
        barmode="group",
        title=dict(
            text="PathwayBench reference (10 disease datasets)",
            font_size=16,
        ),
        annotations=[
            dict(
                text="Averaged across 10 CellxGene datasets | Scaled per criterion",
                xref="paper", yref="paper", x=0, y=1.06,
                showarrow=False, font=dict(size=12, color="grey"),
            )
        ],
        yaxis=dict(title="Scaled score", range=[0, 1.08]),
        legend=dict(title="Method"),
        height=400,
        margin=dict(t=80),
    )
    return fig


def plot_heatmap(scores_df: pd.DataFrame, title: str = "Pathway scores") -> go.Figure:
    """Heatmap of pathway scores across samples."""
    fig = go.Figure(
        go.Heatmap(
            z=scores_df.values,
            x=scores_df.columns.tolist(),
            y=scores_df.index.tolist(),
            colorscale="RdBu_r",
            zmid=0,
            colorbar=dict(title="Score"),
        )
    )
    fig.update_layout(
        title=dict(text=title, font_size=16),
        xaxis=dict(title="Samples", tickangle=45),
        yaxis=dict(title=""),
        height=max(300, 50 * len(scores_df)),
        margin=dict(l=120, b=100),
    )
    return fig


# ---------------------------------------------------------------------------
# Streamlit App
# ---------------------------------------------------------------------------


def main() -> None:
    st.set_page_config(
        page_title="PathwayBench Advisor",
        page_icon="🧬",
        layout="wide",
    )

    # ---- Header ----
    st.markdown(
        f"""
        <div style="background:linear-gradient(135deg,#1a1a2e 0%,#16213e 50%,#0f3460 100%);
                    color:white; padding:28px 32px 22px; border-radius:0 0 8px 8px;
                    margin:-1rem -1rem 1.5rem -1rem;">
            <h1 style="margin:0 0 4px; font-size:26px; font-weight:700;
                       letter-spacing:-0.3px;">PathwayBench Advisor</h1>
            <p style="margin:0; font-size:13px; color:#a0b4d0;">
                Upload pseudobulk data and tell us your priorities.
                PathwayBench will recommend the best scoring method
                <b>based on what you value most</b> and run all five methods
                so you can compare. v{APP_VERSION}</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # ---- Sidebar ----
    with st.sidebar:
        st.markdown(
            f"""
            **📄 Paper**

            [![DOI](https://zenodo.org/badge/DOI/{ZENODO_DOI}.svg)](https://doi.org/{ZENODO_DOI})

            Alakwaa et al., *Nature Methods* (submitted, 2026).

            [Zenodo](https://doi.org/{ZENODO_DOI}) · [GitHub](https://github.com/fadhlyemen/PathwayBench)
            """
        )
        st.divider()
        # Example data
        st.markdown("### Example Data")
        st.caption(
            "New to PathwayBench? Download example files to see the "
            "required format, then re-upload them to test the tool."
        )
        col_a, col_b = st.columns(2)
        with col_a:
            st.download_button(
                "Expression matrix",
                data=load_example("example_pseudobulk.csv"),
                file_name="example_pseudobulk.csv",
                mime="text/csv",
            )
        with col_b:
            st.download_button(
                "Gene set",
                data=load_example("example_geneset.txt"),
                file_name="example_geneset.txt",
                mime="text/plain",
            )

        st.divider()

        # Upload
        st.markdown("### 1. Upload Data")
        expr_file = st.file_uploader(
            "Expression matrix (CSV)",
            type=["csv"],
            help="Rows = genes, columns = samples. Values = log2-CPM or normalized expression.",
        )
        gs_file = st.file_uploader(
            "Gene set (TXT, one gene per line)",
            type=["txt", "csv", "gmt"],
            help="One gene symbol per line, or comma-separated.",
        )

        st.markdown("**Or paste gene list**")
        gene_paste = st.text_area(
            "Genes",
            placeholder="One gene per line, or comma-separated\ne.g. TP53, BRCA1, EGFR ...",
            height=100,
            label_visibility="collapsed",
        )

        st.divider()

        # ---- STEP 1: User priorities elicitation ----
        st.markdown("### 2. Your priorities")
        st.caption(
            "_PathwayBench shows there is no universally best method "
            "-- the right choice depends on what you value._"
        )

        priority = st.radio(
            "What matters most for your analysis?",
            options=[
                "Biological sensitivity (detect the right direction & effect)",
                "Balanced (good on both sensitivity and robustness)",
                "Technical robustness (stable across pipelines & subsets)",
            ],
            index=1,
            key="priority",
            help=(
                "Sensitivity = does the method find true biological signal? "
                "Robustness = does the method give consistent answers under "
                "different aggregations, normalizations, and donor subsets?"
            ),
        )

        sub_priority_bio = None
        if priority.startswith("Biological"):
            sub_priority_bio = st.radio(
                "Within biological sensitivity, what matters more?",
                options=[
                    "Direction accuracy (is the sign correct?)",
                    "Effect magnitude (how large is the difference?)",
                ],
                index=0,
                key="sub_priority_bio",
                help=(
                    "Direction accuracy: best for hypothesis testing. "
                    "Effect magnitude: best for ranking pathways."
                ),
            )

        sub_priority_robust = None
        if priority.startswith("Technical"):
            n_donors = st.number_input(
                "How many donors are in your dataset?",
                min_value=2, max_value=10000, value=20, step=1,
                key="n_donors_input",
                help="Sample-size stability differs by donor count.",
            )
            sub_priority_robust = (
                "\u226515 donors" if n_donors >= 15 else "<15 donors"
            )
            st.session_state["sub_priority_robust"] = sub_priority_robust

        st.divider()
        run_clicked = st.button("Run Analysis", type="primary", use_container_width=True)

    # ---- Parse inputs ----
    expr_df: pd.DataFrame | None = None
    gene_set: list[str] = []
    profile: dict | None = None
    rec: dict | None = None
    zscore_scores: pd.Series | None = None

    if expr_file is not None:
        expr_df = parse_expression(expr_file)

    gene_set = parse_gene_set(gs_file, gene_paste)

    if expr_df is not None:
        profile = profile_data(expr_df)

    # ---- Run analysis ----
    if run_clicked:
        if expr_df is None:
            st.warning("Please upload an expression matrix first.")
            st.stop()
        if not gene_set:
            st.warning("Please provide at least one gene set.")
            st.stop()

        with st.spinner("Running PathwayBench analysis..."):
            # Z-score scoring (pure Python)
            zscore_scores = score_zscore(expr_df, gene_set)
            # Recommendation from benchmark reference
            rec = compute_recommendation(profile, zscore_scores)

        # Persist in session state so tabs work after re-render
        st.session_state["rec"] = rec
        st.session_state["zscore_scores"] = zscore_scores
        st.session_state["profile"] = profile
        st.session_state["gene_set"] = gene_set
        st.session_state["expr_df"] = expr_df

    # Restore from session state
    rec = st.session_state.get("rec")
    zscore_scores = st.session_state.get("zscore_scores")
    profile = st.session_state.get("profile", profile)
    gene_set = st.session_state.get("gene_set", gene_set)
    expr_df_state = st.session_state.get("expr_df", expr_df)

    # ---- Profile cards ----
    if profile is not None:
        st.markdown("### Data Profile")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Genes", f"{profile['n_genes']:,}")
        c2.metric("Samples", profile["n_samples"])
        c3.metric("Outliers", f"{profile['n_outliers']} ({profile['outlier_pct']}%)")
        c4.metric("Lib-size CV", profile["lib_size_cv"])

        detail_col1, detail_col2 = st.columns(2)
        with detail_col1:
            st.caption(f"**Normalization:** {profile['norm_guess']}")
        with detail_col2:
            if profile["conditions"]:
                cond_str = ", ".join(f"{k}={v}" for k, v in profile["conditions"].items())
                st.caption(f"**Conditions:** {cond_str}")
        st.caption(f"**Example genes:** {profile['gene_example']}")

    # ---- Recommendation banner ----
    if rec is not None:
        best = rec["best_method"]
        color = METHOD_COLORS[best]
        st.markdown(
            f"""
            <div style="background:linear-gradient(135deg,#e8f5e9,#f1f8e9);
                        border:2px solid #66bb6a; border-radius:8px;
                        padding:20px 24px; margin:12px 0 18px;">
                <h2 style="color:#2e7d32; margin:0 0 6px; font-size:20px;">
                    Recommended:
                    <span style="background:{color}; color:white; padding:3px 12px;
                                 border-radius:12px; font-size:16px;">{best}</span>
                </h2>
                <p style="color:#424242; margin:0; font-size:13px; line-height:1.55;">
                    {rec['explanation']}</p>
            </div>
            """,
            unsafe_allow_html=True,
        )

    # ---- Priority-based recommendation (STEP 3) ----
    if rec is not None:
        # Detect sctransform from profiler
        sct_used = (
            profile is not None
            and "sctransform" in profile.get("norm_guess", "").lower()
        )
        gs_size = len(gene_set) if gene_set else 50
        outlier_frac = profile["outlier_pct"] if profile else 0.0

        pri_method, pri_reason, pri_warnings = get_priority_based_recommendation(
            priority=st.session_state.get("priority", "Balanced"),
            sub_priority_bio=st.session_state.get("sub_priority_bio"),
            sub_priority_robust=st.session_state.get("sub_priority_robust"),
            sctransform_used=sct_used,
            gene_set_size=gs_size,
            outlier_fraction=outlier_frac,
        )
        st.success(f"### For your priorities: **{pri_method}**")
        st.markdown(f"**Why:** {pri_reason}")
        if pri_warnings:
            st.warning("**Additional considerations for your data:**")
            for w in pri_warnings:
                st.markdown(f"- {w}")
        st.markdown("---")
        st.markdown("### Full results across all five methods")

    # ---- Tabs ----
    if rec is None:
        st.info("Upload data, provide a gene set, and click **Run Analysis**.")
        return

    tab_crit, tab_scores, tab_compare, tab_bench, tab_rankwin = st.tabs(
        ["Criteria Performance", "Pathway Scores", "Method Comparison",
         "Benchmark Reference", "Rank-Window Competition"]
    )

    # -- Tab 1: Criteria performance --
    with tab_crit:
        fig = plot_criteria_bars(rec["scaled"], rec["best_method"], rec["weights"])
        st.plotly_chart(fig, use_container_width=True)

    # -- Tab 2: Pathway scores --
    with tab_scores:
        if zscore_scores is not None and len(zscore_scores) > 0:
            scores_df = pd.DataFrame({"zscore": zscore_scores}).T
            scores_df.columns = [c.split("_")[0][:8] + ".." if len(c) > 10 else c for c in scores_df.columns]
            fig = plot_heatmap(scores_df, title="Z-score pathway activity (computed on your data)")
            st.plotly_chart(fig, use_container_width=True)

            # Show disease vs control means if detectable
            if expr_df_state is not None:
                disease_cols = [c for c in expr_df_state.columns if "disease" in c.lower()]
                control_cols = [c for c in expr_df_state.columns if "control" in c.lower()]
                if disease_cols and control_cols:
                    d_mean = float(zscore_scores[disease_cols].mean()) if set(disease_cols) <= set(zscore_scores.index) else None
                    c_mean = float(zscore_scores[control_cols].mean()) if set(control_cols) <= set(zscore_scores.index) else None
                    if d_mean is not None and c_mean is not None:
                        mc1, mc2, mc3 = st.columns(3)
                        mc1.metric("Disease mean", f"{d_mean:.3f}")
                        mc2.metric("Control mean", f"{c_mean:.3f}")
                        mc3.metric("Difference", f"{d_mean - c_mean:+.3f}")
        else:
            st.warning("Z-score could not be computed (< 5 genes overlap with gene set).")

        st.info(
            "**Note:** Full multi-method scoring (ssGSEA, GSVA, AUCell, UCell) requires "
            "the R backend. Use the Shiny app (`tool/app.R`) or run the Snakemake pipeline "
            "for complete scoring. The Z-score method is computed here in pure Python."
        )

    # -- Tab 3: Method comparison --
    with tab_compare:
        st.markdown("#### Method performance across 6 criteria")
        st.caption("Values from the PathwayBench benchmark (10 disease datasets).")

        display_df = rec["metrics"].copy()
        display_df = display_df.set_index("method")
        display_df.columns = [CRITERIA_SHORT.get(c, c) for c in display_df.columns]

        # Composite from scaled
        composite = rec["scaled"].set_index("method")["composite"]
        display_df["Composite score"] = composite.values

        # Highlight best
        st.dataframe(
            display_df.style.format("{:.3f}").highlight_max(axis=0, color="#c8e6c9"),
            use_container_width=True,
        )

        # Weight explanation
        st.markdown("##### Applied weights")
        wt_df = pd.DataFrame(
            {
                "Criterion": [CRITERIA_SHORT[c] for c in CRITERIA],
                "Default weight": [DEFAULT_WEIGHTS[c] for c in CRITERIA],
                "Applied weight": [round(rec["weights"][c], 3) for c in CRITERIA],
            }
        )
        st.dataframe(wt_df, use_container_width=True, hide_index=True)

    # -- Tab 4: Benchmark reference --
    with tab_bench:
        fig = plot_benchmark_ref(_BENCH_RAW)
        st.plotly_chart(fig, use_container_width=True)

        st.markdown("##### Raw benchmark values (averaged across 10 datasets)")
        bench_display = _BENCH_RAW.copy().set_index("method")
        bench_display.columns = [CRITERIA_SHORT.get(c, c) for c in bench_display.columns]
        st.dataframe(
            bench_display.style.format("{:.3f}").highlight_max(axis=0, color="#c8e6c9"),
            use_container_width=True,
        )

    # -- Tab 5: Rank-window competition --
    with tab_rankwin:
        st.markdown("### Rank-window competition: a novel failure mode")
        st.markdown(
            "We discovered that rank-based methods (AUCell, UCell) can **invert** "
            "the disease-vs-control direction when non-pathway genes are upregulated "
            "more strongly, displacing pathway genes from the fixed top-N rank window."
        )

        sim_data = pd.DataFrame({
            "Scenario": ["a: Top 10%, clean"] * 5 + ["b: Top 50%, clean"] * 5 +
                        ["c: Top 10%, 200 comp"] * 5 + ["d: Top 10%, 500 comp"] * 5,
            "Method": ["ssGSEA","GSVA","zscore","AUCell","UCell"] * 4,
            "Cohen_d": [9.16, 13.49, 9.28, 7.14, 9.15,
                        8.93, 13.12, 9.07, 0.00, 0.22,
                        6.18, 12.57, 9.17, 4.00, 6.19,
                        2.26, 11.91, 9.40, 0.05, 2.28],
        })
        sim_data["color"] = sim_data.apply(
            lambda r: "#D32F2F" if r["Cohen_d"] < 0.3 and r["Method"] in ("AUCell", "UCell")
            else METHOD_COLORS.get(r["Method"], "#999"), axis=1
        )
        fig_sim = go.Figure()
        for sc in sim_data["Scenario"].unique():
            sub = sim_data[sim_data["Scenario"] == sc]
            fig_sim.add_trace(go.Bar(
                name=sc, x=sub["Method"], y=sub["Cohen_d"],
                marker_color=sub["color"].tolist(),
                showlegend=False,
            ))
        fig_sim.update_layout(
            barmode="group", title="Simulation: Cohen's d under rank-window competition",
            yaxis_title="Cohen's d", height=400,
            annotations=[dict(text="Red = missed or wrong sign", x=0.5, y=1.05,
                              xref="paper", yref="paper", showarrow=False,
                              font=dict(size=11, color="grey"))],
        )
        st.plotly_chart(fig_sim, use_container_width=True)

        st.markdown(
            "**Real-data validation (CKD kidney, ECM remodeling):** "
            "ssGSEA d=+0.39, GSVA d=+0.40, zscore d=+0.44 (detect UP). "
            "AUCell d=+0.03, UCell d=-0.05 (miss or invert)."
        )
        st.info(
            "**Recommendation:** When AUCell/UCell and magnitude methods disagree, "
            "trust the magnitude method. Avoid AUCell/UCell as sole method in diseases "
            "with broad transcriptional remodeling (CKD, fibrotic conditions)."
        )

    # ---- Footer ----
    st.divider()
    st.caption(
        f"PathwayBench Advisor v{APP_VERSION} | Benchmarking 5 pathway scoring methods "
        f"across 5 robustness criteria | Reference: 7 disease datasets from "
        f"CellxGene Census 2025-11-08 | "
        f"[DOI: {ZENODO_DOI}](https://doi.org/{ZENODO_DOI})"
    )


if __name__ == "__main__":
    main()
