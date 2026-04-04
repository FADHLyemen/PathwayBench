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

# Benchmark reference — averaged across 10 CellxGene disease datasets.
# outlier_sensitivity is inverted to outlier_robustness (1 - value).
_BENCH_RAW = pd.DataFrame(
    {
        "method": ["ssGSEA", "GSVA", "zscore", "AUCell", "UCell"],
        "bio_relevance": [0.354, 0.398, 0.340, 0.380, 0.358],
        "calc_method_cor": [0.866, 0.668, 0.796, 0.816, 0.925],
        "outlier_robustness": [1 - 0.251, 1 - 0.255, 1 - 0.256, 1 - 0.195, 1 - 0.195],
        "norm_stability": [0.681, 0.696, 0.705, 0.648, 0.757],
        "sample_stability": [0.888, 0.866, 0.885, 0.856, 0.827],
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
        """
        <div style="background:linear-gradient(135deg,#1a1a2e 0%,#16213e 50%,#0f3460 100%);
                    color:white; padding:28px 32px 22px; border-radius:0 0 8px 8px;
                    margin:-1rem -1rem 1.5rem -1rem;">
            <h1 style="margin:0 0 4px; font-size:26px; font-weight:700;
                       letter-spacing:-0.3px;">PathwayBench Advisor</h1>
            <p style="margin:0; font-size:13px; color:#a0b4d0;">
                Upload your pseudobulk data. Get a benchmarked recommendation
                for the best pathway scoring method.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # ---- Sidebar ----
    with st.sidebar:
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

    # ---- Tabs ----
    if rec is None:
        st.info("Upload data, provide a gene set, and click **Run Analysis**.")
        return

    tab_crit, tab_scores, tab_compare, tab_bench = st.tabs(
        ["Criteria Performance", "Pathway Scores", "Method Comparison", "Benchmark Reference"]
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

    # ---- Footer ----
    st.divider()
    st.caption(
        "PathwayBench Advisor | Benchmarking 5 pathway scoring methods across "
        "6 robustness criteria | Reference database: 10 disease datasets from CellxGene Census"
    )


if __name__ == "__main__":
    main()
