#!/usr/bin/env python3
"""Integrate expression, metadata, and somatic variants into reproducible research outputs."""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from jinja2 import Template

LOG = logging.getLogger("oncoflow")
REQUIRED_META = {"sample_id", "subtype", "purity"}
REQUIRED_VARIANTS = {"sample_id", "gene", "variant_class", "vaf", "impact"}
PALETTE = {"Classical": "#43D9C4", "Basal-like": "#8E7CFF"}
PLOT_BG = "#0B1420"
GRID = "#1D3446"
TEXT = "#EAF4F8"
MUTED = "#8FA8B5"

HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>OncoFlow Atlas · {{ cohort }}</title>
<style>
:root{--bg:#071019;--surface:#0b1622;--surface2:#0e1c2a;--line:#1a3346;--ink:#edf7fa;--muted:#8fa8b5;--teal:#43d9c4;--cyan:#5bb8e8;--violet:#8e7cff;--amber:#f4bf64;--good:#4fd1a5;--warn:#f1b75b}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:radial-gradient(circle at 88% 4%,rgba(26,109,119,.34) 0,transparent 24%),radial-gradient(circle at 14% 16%,rgba(74,71,145,.22) 0,transparent 25%),linear-gradient(180deg,#061019 0%,#08131d 56%,#071019 100%);color:var(--ink);font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.shell{max-width:1360px;margin:auto;padding:30px}.hero{position:relative;overflow:hidden;background:linear-gradient(128deg,rgba(13,36,49,.96) 0%,rgba(13,47,54,.92) 44%,rgba(24,41,74,.94) 100%);border:1px solid #2a5268;border-radius:30px;padding:42px 44px;box-shadow:0 28px 90px rgba(0,0,0,.34),inset 0 1px 0 rgba(255,255,255,.035)}.hero:after{content:"";position:absolute;width:430px;height:430px;border-radius:50%;right:-110px;top:-210px;background:radial-gradient(circle,rgba(67,217,196,.22),transparent 67%)}.eyebrow{font-size:12px;text-transform:uppercase;letter-spacing:.18em;color:var(--teal);font-weight:800}.hero h1{margin:10px 0 9px;font-size:46px;letter-spacing:-.035em}.hero p{max-width:760px;color:#a9bdc7;margin:0;font-size:16px;line-height:1.65}.chips{display:flex;flex-wrap:wrap;gap:9px;margin-top:22px}.chip{font-size:12px;border:1px solid #2b5164;background:rgba(7,16,25,.34);padding:8px 11px;border-radius:999px;color:#c6e5e5}.nav{position:sticky;top:10px;z-index:20;display:flex;gap:8px;overflow:auto;margin:16px 0;padding:7px;border:1px solid rgba(39,78,98,.72);border-radius:15px;background:rgba(6,16,25,.82);backdrop-filter:blur(14px);box-shadow:0 10px 30px rgba(0,0,0,.18)}.nav a{text-decoration:none;color:#bcd0d8;border:1px solid transparent;padding:8px 12px;border-radius:10px;font-size:12px;transition:.2s ease}.nav a:hover{color:#ecffff;border-color:#275069;background:#0f2130}.grid{display:grid;grid-template-columns:repeat(5,1fr);gap:12px}.metric{position:relative;overflow:hidden;background:linear-gradient(165deg,#0d1d2b,#0a1520 70%);border:1px solid #1c3a4d;border-radius:19px;padding:18px;box-shadow:inset 0 1px 0 rgba(255,255,255,.025)}.metric:after{content:"";position:absolute;width:100px;height:100px;border-radius:50%;right:-45px;bottom:-55px;background:radial-gradient(circle,rgba(67,217,196,.18),transparent 68%)}.metric b{display:block;font-size:29px;letter-spacing:-.03em;margin-top:4px}.metric span{font-size:12px;color:var(--muted)}.section{margin-top:18px;background:linear-gradient(180deg,rgba(13,25,37,.98),rgba(8,18,28,.98));border:1px solid #19374a;border-radius:24px;padding:22px;box-shadow:0 18px 46px rgba(0,0,0,.17),inset 0 1px 0 rgba(255,255,255,.02)}.section-head{display:flex;justify-content:space-between;gap:18px;align-items:end;margin-bottom:8px}.section h2{font-size:19px;margin:0}.section p{font-size:12px;color:var(--muted);margin:5px 0 0;line-height:1.55}.twocol{display:grid;grid-template-columns:1fr 1fr;gap:16px}.table{width:100%;border-collapse:collapse;font-size:12px}.table th,.table td{text-align:left;padding:10px 11px;border-bottom:1px solid #183142}.table th{color:#a8c1cb;font-weight:650}.status-PASS{color:var(--good);font-weight:700}.status-REVIEW{color:var(--warn);font-weight:700}.footer{display:flex;justify-content:space-between;gap:20px;color:#6f8995;font-size:11px;padding:20px 2px}.note{border-left:3px solid var(--teal);background:#0d1c26;padding:12px 14px;border-radius:0 12px 12px 0;color:#a9c0ca;font-size:12px;line-height:1.6}@media(max-width:950px){.grid{grid-template-columns:repeat(3,1fr)}.twocol{grid-template-columns:1fr}}@media(max-width:650px){.shell{padding:14px}.hero{padding:28px 24px}.hero h1{font-size:35px}.grid{grid-template-columns:repeat(2,1fr)}}
</style></head><body><main class="shell">
<section class="hero"><div class="eyebrow">Reproducible cancer genomics workflow · v1.0</div><h1>OncoFlow Atlas</h1><p>Integrated multi-omic cohort exploration with transparent quality control, sample-level metrics, somatic-variant summaries, portable execution and verifiable run provenance.</p><div class="chips"><span class="chip">{{ cohort }}</span><span class="chip">{{ n_samples }} samples</span><span class="chip">{{ n_genes }} expression features</span><span class="chip">{{ n_variants }} retained variants</span><span class="chip">VAF ≥ {{ min_vaf }}</span></div></section>
<nav class="nav"><a href="#overview">Overview</a><a href="#expression">Expression</a><a href="#variants">Variants</a><a href="#integration">Integration</a><a href="#qc">QC</a><a href="#provenance">Provenance</a></nav>
<section class="grid" id="overview"><div class="metric"><span>Samples</span><b>{{ n_samples }}</b></div><div class="metric"><span>Retained variants</span><b>{{ n_variants }}</b></div><div class="metric"><span>Genes profiled</span><b>{{ n_genes }}</b></div><div class="metric"><span>Median purity</span><b>{{ purity }}</b></div><div class="metric"><span>QC pass rate</span><b>{{ qc_pass_rate }}%</b></div></section>
<section class="section" id="expression"><div class="section-head"><div><h2>Expression structure</h2><p>PCA on log-transformed expression. PC1 {{ pc1_var }}% · PC2 {{ pc2_var }}% variance explained.</p></div></div>{{ pca|safe }}</section>
<section class="twocol"><div class="section"><div class="section-head"><div><h2>Subtype expression signature</h2><p>Top variable genes, standardized across samples.</p></div></div>{{ heatmap|safe }}</div><div class="section"><div class="section-head"><div><h2>Sample correlation</h2><p>Expression correlation across the cohort.</p></div></div>{{ correlation|safe }}</div></section>
<section class="twocol" id="variants"><div class="section"><div class="section-head"><div><h2>Mutation burden by gene</h2><p>Variants remaining after the configurable VAF threshold.</p></div></div>{{ mutation_bar|safe }}</div><div class="section"><div class="section-head"><div><h2>Variant impact composition</h2><p>Distribution of annotated impact classes.</p></div></div>{{ impact_donut|safe }}</div></section>
<section class="section"><div class="section-head"><div><h2>Mutation landscape</h2><p>Binary sample × gene matrix for the most frequently altered genes.</p></div></div>{{ oncoprint|safe }}</section>
<section class="section" id="integration"><div class="section-head"><div><h2>Sample-level integration</h2><p>Expression signal, variant burden, purity and detected features in one view.</p></div></div>{{ burden|safe }}</section>
<section class="twocol"><div class="section"><div class="section-head"><div><h2>Subtype summary</h2><p>Compact cohort-level comparison for exploratory review.</p></div></div>{{ subtype_table|safe }}</div><div class="section" id="qc"><div class="section-head"><div><h2>Quality-control summary</h2><p>Pre-analysis checks and interpretable review flags.</p></div></div>{{ qc_table|safe }}</div></section>
<section class="section" id="provenance"><div class="section-head"><div><h2>Reproducibility & provenance</h2><p>Execution details are written to <code>run_manifest.json</code> with SHA-256 checksums for all primary inputs.</p></div></div><div class="note">Generated {{ generated }} · Python {{ python_version }} · {{ platform }}. This report is a research workflow demonstration and is not a clinical diagnostic report.</div></section>
<footer class="footer"><span>OncoFlow Atlas · workflow version 1.0</span><span>Portable · inspectable · reproducible</span></footer>
</main></body></html>'''


def configure_plot(fig: go.Figure, height: int = 430) -> go.Figure:
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor=PLOT_BG,
        plot_bgcolor=PLOT_BG,
        font=dict(color=TEXT, family="Inter, Arial, sans-serif"),
        margin=dict(l=40, r=20, t=30, b=40),
        height=height,
        legend_title_text="",
        hoverlabel=dict(bgcolor="#132636", font_color=TEXT),
    )
    fig.update_xaxes(gridcolor=GRID, zerolinecolor=GRID)
    fig.update_yaxes(gridcolor=GRID, zerolinecolor=GRID)
    return fig


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_inputs(expr: pd.DataFrame, var: pd.DataFrame, meta: pd.DataFrame) -> list[str]:
    problems: list[str] = []
    if "gene" not in expr.columns:
        problems.append("expression table must contain a 'gene' column")
    missing_meta = REQUIRED_META - set(meta.columns)
    if missing_meta:
        problems.append(f"metadata is missing columns: {sorted(missing_meta)}")
    missing_var = REQUIRED_VARIANTS - set(var.columns)
    if missing_var:
        problems.append(f"variant table is missing columns: {sorted(missing_var)}")
    if problems:
        return problems

    meta_samples = set(meta["sample_id"].astype(str))
    expr_samples = set(expr.columns) - {"gene"}
    variant_samples = set(var["sample_id"].astype(str)) if len(var) else set()
    absent_expr = meta_samples - expr_samples
    unknown_variants = variant_samples - meta_samples
    if absent_expr:
        problems.append(f"metadata samples absent from expression matrix: {sorted(absent_expr)}")
    if unknown_variants:
        problems.append(f"variant rows reference unknown samples: {sorted(unknown_variants)}")
    if meta["sample_id"].duplicated().any():
        problems.append("metadata contains duplicated sample_id values")
    if expr["gene"].duplicated().any():
        problems.append("expression table contains duplicated gene values")
    if not pd.to_numeric(meta["purity"], errors="coerce").between(0, 1).all():
        problems.append("purity values must be numeric and between 0 and 1")
    if len(var) and not pd.to_numeric(var["vaf"], errors="coerce").between(0, 1).all():
        problems.append("VAF values must be numeric and between 0 and 1")
    return problems


def pca_numpy(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    centered = frame - frame.mean(axis=0)
    u, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
    coords = u[:, :2] * singular_values[:2]
    denom = np.square(singular_values).sum()
    variance = np.square(singular_values[:2]) / denom * 100 if denom else np.array([0.0, 0.0])
    return coords, variance


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--variants", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--cohort-name", default="demo-cohort")
    parser.add_argument("--min-vaf", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=12, help="Number of features shown in top-ranked visualisations")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s · %(message)s")
    if not 0 <= args.min_vaf <= 1:
        raise SystemExit("--min-vaf must be between 0 and 1")
    if args.top_n < 3:
        raise SystemExit("--top-n must be at least 3")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    LOG.info("Loading input tables")
    expr_raw = pd.read_csv(args.expression)
    variants_raw = pd.read_csv(args.variants)
    meta = pd.read_csv(args.metadata)
    errors = validate_inputs(expr_raw, variants_raw, meta)
    if errors:
        raise SystemExit("Input validation failed:\n- " + "\n- ".join(errors))

    expr = expr_raw.set_index("gene")
    sample_cols = list(meta["sample_id"].astype(str))
    numeric_expr = expr[sample_cols].apply(pd.to_numeric, errors="coerce")
    missing_expr = int(numeric_expr.isna().sum().sum())
    numeric_expr = numeric_expr.fillna(0)
    log_expr = np.log1p(numeric_expr.T)

    variants = variants_raw.copy()
    variants["vaf"] = pd.to_numeric(variants["vaf"], errors="coerce")
    variants = variants.loc[variants["vaf"] >= args.min_vaf].copy()
    LOG.info("Retained %s/%s variants at VAF >= %.3f", len(variants), len(variants_raw), args.min_vaf)

    pcs, explained = pca_numpy(log_expr)
    pca = pd.DataFrame({"sample_id": sample_cols, "PC1": pcs[:, 0], "PC2": pcs[:, 1]}).merge(meta, on="sample_id")
    fig_pca = px.scatter(
        pca, x="PC1", y="PC2", color="subtype", symbol="batch" if "batch" in meta.columns else None,
        hover_data=["sample_id", "purity"], color_discrete_map=PALETTE,
    )
    configure_plot(fig_pca, 470)
    fig_pca.update_traces(marker=dict(size=11, line=dict(width=1, color="#DDF7F5")))

    variance = log_expr.var(axis=0).sort_values(ascending=False)
    top_genes = variance.head(min(args.top_n, len(variance))).index.tolist()
    z = log_expr[top_genes].copy()
    z = (z - z.mean(axis=0)) / z.std(axis=0).replace(0, 1)
    order = meta.set_index("sample_id").loc[z.index, "subtype"].sort_values().index
    z = z.loc[order]
    fig_heat = px.imshow(
        z.T, x=z.index, y=z.columns, color_continuous_scale=[[0, "#283B62"], [0.5, "#101C2A"], [1, "#43D9C4"]],
        aspect="auto", labels=dict(color="z-score"), zmin=-2.5, zmax=2.5,
    )
    configure_plot(fig_heat, 470)
    fig_heat.update_xaxes(showticklabels=False, title="Samples")

    corr = log_expr.T.corr()
    fig_corr = px.imshow(corr, color_continuous_scale=[[0, "#111B2C"], [0.5, "#2D3C68"], [1, "#43D9C4"]], zmin=0, zmax=1)
    configure_plot(fig_corr, 470)
    fig_corr.update_xaxes(showticklabels=False, title="Samples")
    fig_corr.update_yaxes(showticklabels=False, title="Samples")

    variant_summary = variants.groupby("gene").size().sort_values(ascending=False).rename("variants").reset_index()
    fig_mut = px.bar(variant_summary.head(args.top_n), x="gene", y="variants", color="variants", color_continuous_scale=["#40507D", "#43D9C4"])
    configure_plot(fig_mut, 410)
    fig_mut.update_layout(coloraxis_showscale=False)

    impact_counts = variants["impact"].value_counts().rename_axis("impact").reset_index(name="variants")
    fig_impact = px.pie(impact_counts, names="impact", values="variants", hole=0.62, color_discrete_sequence=["#43D9C4", "#8E7CFF", "#5BB8E8", "#F4BF64"])
    configure_plot(fig_impact, 410)
    fig_impact.update_traces(textposition="inside", textinfo="percent+label")

    top_mutated = variant_summary.head(min(args.top_n, len(variant_summary)))["gene"].tolist()
    matrix = pd.crosstab(variants["gene"], variants["sample_id"]).reindex(top_mutated).fillna(0).clip(upper=1)
    matrix = matrix.reindex(columns=sample_cols, fill_value=0)
    fig_onco = px.imshow(matrix, color_continuous_scale=[[0, "#101C2A"], [1, "#F0648A"]], zmin=0, zmax=1, aspect="auto")
    configure_plot(fig_onco, 380)
    fig_onco.update_layout(coloraxis_showscale=False)
    fig_onco.update_xaxes(showticklabels=False, title="Samples")

    burden = variants.groupby("sample_id").size().rename("variant_count").reset_index()
    metrics = meta.merge(burden, on="sample_id", how="left").fillna({"variant_count": 0})
    metrics["variant_count"] = metrics["variant_count"].astype(int)
    metrics["expression_signal"] = log_expr.mean(axis=1).reindex(metrics.sample_id).values
    metrics["detected_genes"] = (numeric_expr[sample_cols] > 0).sum(axis=0).reindex(metrics.sample_id).values
    fig_integrated = px.scatter(
        metrics, x="expression_signal", y="variant_count", color="subtype", size="purity",
        hover_name="sample_id", hover_data=["detected_genes", "purity"], color_discrete_map=PALETTE,
    )
    configure_plot(fig_integrated, 450)

    subtype_summary = (
        metrics.groupby("subtype", as_index=False)
        .agg(samples=("sample_id", "count"), median_purity=("purity", "median"), median_variants=("variant_count", "median"), mean_expression_signal=("expression_signal", "mean"))
    )
    subtype_summary["median_purity"] = subtype_summary["median_purity"].round(3)
    subtype_summary["median_variants"] = subtype_summary["median_variants"].round(1)
    subtype_summary["mean_expression_signal"] = subtype_summary["mean_expression_signal"].round(3)

    qc_rows = [
        {"metric": "Expression missing cells", "value": missing_expr, "status": "PASS" if missing_expr == 0 else "REVIEW"},
        {"metric": "Duplicated metadata samples", "value": int(meta.sample_id.duplicated().sum()), "status": "PASS"},
        {"metric": "Samples with purity < 0.40", "value": int((meta.purity < 0.40).sum()), "status": "PASS" if int((meta.purity < 0.40).sum()) == 0 else "REVIEW"},
        {"metric": "Variants below configured VAF", "value": int(len(variants_raw) - len(variants)), "status": "PASS"},
        {"metric": "Median sample purity", "value": round(float(meta.purity.median()), 3), "status": "PASS"},
        {"metric": "Median variants per sample", "value": round(float(metrics.variant_count.median()), 2), "status": "PASS"},
    ]
    qc = pd.DataFrame(qc_rows)
    qc_pass_rate = round(float((qc.status == "PASS").mean() * 100))

    variance.rename("variance").reset_index().rename(columns={"index": "gene"}).head(args.top_n).to_csv(out / "top_genes.csv", index=False)
    variant_summary.to_csv(out / "variant_summary.csv", index=False)
    metrics.to_csv(out / "sample_metrics.csv", index=False)
    pca.to_csv(out / "pca_coordinates.csv", index=False)
    subtype_summary.to_csv(out / "subtype_summary.csv", index=False)
    qc.to_csv(out / "qc_summary.csv", index=False)

    summary = {
        "cohort": args.cohort_name,
        "samples": len(meta),
        "variants_input": len(variants_raw),
        "variants_retained": len(variants),
        "genes": len(expr),
        "median_purity": round(float(meta.purity.median()), 3),
        "subtypes": meta.subtype.value_counts().to_dict(),
        "parameters": {"min_vaf": args.min_vaf, "top_n": args.top_n},
        "pca_variance_explained_pct": {"PC1": round(float(explained[0]), 2), "PC2": round(float(explained[1]), 2)},
        "qc_pass_rate_pct": qc_pass_rate,
    }
    (out / "cohort_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    generated = datetime.now(timezone.utc)
    inputs = {"expression": args.expression, "variants": args.variants, "metadata": args.metadata}
    manifest = {
        "workflow": {"name": "OncoFlow Atlas", "version": "1.0.0"},
        "generated_utc": generated.isoformat(),
        "cohort": args.cohort_name,
        "parameters": {"min_vaf": args.min_vaf, "top_n": args.top_n},
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": version("numpy"),
            "pandas": version("pandas"),
            "plotly": version("plotly"),
        },
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in inputs.items()},
    }
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    qc_html = qc.copy()
    qc_html["status"] = qc_html["status"].map(lambda s: f'<span class="status-{s}">{s}</span>')
    html = Template(HTML).render(
        generated=generated.strftime("%Y-%m-%d %H:%M UTC"),
        cohort=args.cohort_name,
        n_samples=len(meta),
        n_variants=len(variants),
        n_genes=len(expr),
        purity=f"{meta.purity.median():.2f}",
        qc_pass_rate=qc_pass_rate,
        min_vaf=f"{args.min_vaf:.2f}",
        pc1_var=f"{explained[0]:.1f}",
        pc2_var=f"{explained[1]:.1f}",
        pca=fig_pca.to_html(full_html=False, include_plotlyjs=True, config={"displaylogo": False}),
        heatmap=fig_heat.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False}),
        correlation=fig_corr.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False}),
        mutation_bar=fig_mut.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False}),
        impact_donut=fig_impact.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False}),
        oncoprint=fig_onco.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False}),
        burden=fig_integrated.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False}),
        subtype_table=subtype_summary.to_html(index=False, classes="table", border=0),
        qc_table=qc_html.to_html(index=False, classes="table", border=0, escape=False),
        python_version=sys.version.split()[0],
        platform=platform.system(),
    )
    (out / "report.html").write_text(html, encoding="utf-8")

    output_files = [p for p in out.iterdir() if p.is_file() and p.name != "run_manifest.json"]
    manifest["outputs"] = {p.name: {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in sorted(output_files)}
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    LOG.info("Analysis complete: %s", out)


if __name__ == "__main__":
    main()
