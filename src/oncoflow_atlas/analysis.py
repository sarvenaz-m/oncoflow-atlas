"""Numerical analysis and the browser/Python view contract."""

from __future__ import annotations

import hashlib
import json
from collections import Counter

import numpy as np
import pandas as pd

from . import __version__
from .validation import Cohort, InputError, validation_summary

DEFAULTS = {"min_vaf": 0.05, "top_n": 12, "min_purity": 0.4, "subtype": "All", "search": ""}


def parameters(overrides: dict | None = None) -> dict:
    p = DEFAULTS | (overrides or {})
    if set(p) != set(DEFAULTS):
        raise InputError("Unknown analysis parameter")
    for key in ("min_vaf", "min_purity"):
        if (
            isinstance(p[key], bool)
            or not isinstance(p[key], (int, float))
            or not np.isfinite(p[key])
            or not 0 <= p[key] <= 1
        ):
            raise InputError(f"{key} must be a finite number in [0,1]")
    if isinstance(p["top_n"], bool) or not isinstance(p["top_n"], int) or not 3 <= p["top_n"] <= 50:
        raise InputError("top_n must be an integer in [3,50]")
    if not isinstance(p["subtype"], str) or not isinstance(p["search"], str):
        raise InputError("subtype and search must be strings")
    p["search"] = p["search"].strip().lower()
    return p


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean_json(v) for v in value]
    if isinstance(value, (np.integer, np.floating)):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def fingerprint(value) -> str:
    return hashlib.sha256(
        json.dumps(
            clean_json(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
        ).encode()
    ).hexdigest()


def pca_numpy(log_expression: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not np.any(np.ptp(log_expression, axis=0) > 0):
        raise InputError("PCA requires at least one variable expression feature")
    centered = np.asarray(log_expression, dtype=float) - np.mean(log_expression, axis=0)
    u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    total = float(np.square(singular).sum())
    if total <= 0:
        raise InputError("PCA requires at least one variable expression feature")
    coords = u[:, :2] * singular[:2]
    loadings = vt[:2].copy()
    for component in range(2):
        pivot = int(np.argmax(np.abs(loadings[component])))
        if loadings[component, pivot] < 0:
            coords[:, component] *= -1
            loadings[component] *= -1
    return coords, np.square(singular[:2]) / total * 100, loadings


def build_payload(
    cohort: Cohort, name: str, options: dict | None = None, source: dict | None = None
) -> tuple[dict, pd.DataFrame]:
    p = parameters(options)
    values = cohort.expression.to_numpy(dtype=float)
    log_values = np.log2(values + 1)
    coords, explained, loadings = pca_numpy(log_values.T)
    # A constant sample profile has undefined Pearson correlation; keep null, not zero.
    with np.errstate(divide="ignore", invalid="ignore"):
        correlation = np.corrcoef(log_values.T)
    matrix = {
        "genes": cohort.expression.index.tolist(),
        "samples": cohort.expression.columns.tolist(),
        "values": values.tolist(),
    }
    cna = (
        None
        if cohort.cna is None
        else {
            "genes": cohort.cna.index.tolist(),
            "samples": cohort.cna.columns.tolist(),
            "values": cohort.cna.to_numpy(dtype=int).tolist(),
        }
    )
    samples = cohort.metadata.sample_id.tolist()
    payload = clean_json(
        {
            "schema_version": "2.0",
            "software": {"name": "OncoFlow Atlas", "version": __version__},
            "cohort": name,
            "source": source
            or {"kind": "user-supplied", "name": name, "expression_units": "linear-scale expression; user supplied"},
            "parameters": p,
            "metadata": cohort.metadata.to_dict("records"),
            "expression": matrix,
            "variants": cohort.variants.to_dict("records"),
            "cna": cna,
            "inputs": cohort.inputs,
            "validation": validation_summary(cohort),
            "analysis": {
                "transform": "log2(value + 1)",
                "variance_ddof": 1,
                "pca_scope": "Fixed full-cohort coordinates; filtering selects points without refitting PCA",
                "pca": [
                    {"sample_id": sid, "PC1": float(coords[i, 0]), "PC2": float(coords[i, 1])}
                    for i, sid in enumerate(samples)
                ],
                "pca_variance_explained_pct": {"PC1": float(explained[0]), "PC2": float(explained[1])},
                "correlation": correlation,
                "constant_expression_genes": int(np.count_nonzero(np.ptp(log_values, axis=1) == 0)),
            },
        }
    )
    payload["data_fingerprint"] = fingerprint({"inputs": cohort.inputs, "source": payload["source"]})
    payload["default_view"] = compute_view(payload, p)
    payload["sensitivity"] = [
        {"min_vaf": threshold, **variant_totals(payload["variants"], threshold)}
        for threshold in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5]
    ]
    loading_frame = pd.DataFrame({"gene": matrix["genes"], "PC1": loadings[0], "PC2": loadings[1]})
    return payload, loading_frame


def variant_totals(variants: list[dict], threshold: float) -> dict:
    kept = [v for v in variants if v["vaf"] is not None and v["vaf"] >= threshold]
    return {
        "retained_records": len(kept),
        "genes": len({v["gene"] for v in kept}),
        "samples_with_record": len({v["sample_id"] for v in kept}),
        "unknown_vaf": sum(v["vaf"] is None for v in variants),
    }


def compute_view(payload: dict, overrides: dict | None = None) -> dict:
    p = parameters(payload["parameters"] | (overrides or {}))
    metadata = [
        m
        for m in payload["metadata"]
        if (p["subtype"] == "All" or m["subtype"] == p["subtype"])
        and (not p["search"] or p["search"] in m["sample_id"].lower())
    ]
    ids = [m["sample_id"] for m in metadata]
    chosen = set(ids)
    indices = [payload["expression"]["samples"].index(s) for s in ids]
    records = [v for v in payload["variants"] if v["sample_id"] in chosen]
    kept = [v for v in records if v["vaf"] is not None and v["vaf"] >= p["min_vaf"]]
    genes = payload["expression"]["genes"]
    values = np.asarray(payload["expression"]["values"], dtype=float)[:, indices]
    logs = np.log2(values + 1)
    if ids:
        mean = logs.mean(axis=1)
        variance = logs.var(axis=1, ddof=1) if len(ids) > 1 else np.zeros(len(genes))
        variance[np.ptp(logs, axis=1) == 0] = 0
        ranked = sorted(range(len(genes)), key=lambda i: (-variance[i], genes[i]))[: p["top_n"]]
    else:
        mean, variance, ranked = np.zeros(len(genes)), np.zeros(len(genes)), []
    top = [{"gene": genes[i], "variance": float(variance[i])} for i in ranked]
    heat = []
    for i in ranked:
        sd = float(np.sqrt(variance[i]))
        heat.append({"gene": genes[i], "z": ((logs[i] - mean[i]) / sd if sd > 0 else np.zeros(len(ids))).tolist()})
    gene_counts = Counter(v["gene"] for v in kept)
    gene_samples = {}
    for v in kept:
        gene_samples.setdefault(v["gene"], set()).add(v["sample_id"])
    variant_summary = [
        {"gene": g, "records": count, "samples": len(gene_samples[g])}
        for g, count in sorted(gene_counts.items(), key=lambda pair: (-pair[1], pair[0]))
    ]
    counts = Counter(v["sample_id"] for v in kept)
    missing_vaf = Counter(v["sample_id"] for v in records if v["vaf"] is None)
    cna = payload.get("cna")
    cna_values = np.asarray(cna["values"], dtype=int)[:, indices] if cna is not None else None
    cna_summary = []
    if cna is not None and ids:
        amps = (cna_values == 2).sum(axis=1)
        dels = (cna_values == -2).sum(axis=1)
        order = sorted(range(len(cna["genes"])), key=lambda i: (-int(amps[i] + dels[i]), cna["genes"][i]))
        cna_summary = [
            {"gene": cna["genes"][i], "amplified_samples": int(amps[i]), "deep_deleted_samples": int(dels[i])}
            for i in order
            if amps[i] + dels[i] > 0
        ]
    metrics = []
    for i, m in enumerate(metadata):
        reasons = []
        if m.get("purity") is None:
            reasons.append("Purity not reported")
        elif m["purity"] < p["min_purity"]:
            reasons.append("Purity below threshold")
        if missing_vaf[m["sample_id"]]:
            reasons.append("Variant VAF unavailable")
        metrics.append(
            m
            | {
                "variant_count": counts[m["sample_id"]],
                "variants_missing_vaf": missing_vaf[m["sample_id"]],
                "expression_signal": float(logs[:, i].mean()),
                "detected_genes": int((values[:, i] > 0).sum()),
                "amplified_genes": int((cna_values[:, i] == 2).sum()) if cna is not None else None,
                "deep_deleted_genes": int((cna_values[:, i] == -2).sum()) if cna is not None else None,
                "qc_status": "REVIEW" if reasons else "PASS",
                "qc_reasons": reasons,
            }
        )
    known_purity = [m["purity"] for m in metadata if m.get("purity") is not None]
    groups = []
    for subtype in sorted({m["subtype"] for m in metrics}):
        group = [m for m in metrics if m["subtype"] == subtype]
        purities = [m["purity"] for m in group if m["purity"] is not None]
        groups.append(
            {
                "subtype": subtype,
                "samples": len(group),
                "median_purity": float(np.median(purities)) if purities else None,
                "median_variants": float(np.median([m["variant_count"] for m in group])),
                "mean_expression_signal": float(np.mean([m["expression_signal"] for m in group])),
            }
        )
    qc = [
        {"metric": "Input identifiers and sample alignment", "value": "Validated", "status": "PASS"},
        {"metric": "Missing expression cells", "value": 0, "status": "PASS"},
        {
            "metric": "Samples without reported purity",
            "value": len(metadata) - len(known_purity),
            "status": "NOT_REPORTED" if len(metadata) > len(known_purity) else "PASS",
        },
        {
            "metric": "Known purity below threshold",
            "value": sum(x < p["min_purity"] for x in known_purity),
            "status": "REVIEW" if any(x < p["min_purity"] for x in known_purity) else "PASS",
        },
        {
            "metric": "Records with unavailable VAF",
            "value": sum(v["vaf"] is None for v in records),
            "status": "REVIEW" if any(v["vaf"] is None for v in records) else "PASS",
        },
        {
            "metric": "Known VAF below threshold",
            "value": sum(v["vaf"] is not None and v["vaf"] < p["min_vaf"] for v in records),
            "status": "FILTERED",
        },
        {
            "metric": "Copy-number profile",
            "value": len(cna["genes"]) if cna is not None else "Not supplied",
            "status": "PASS" if cna is not None else "NOT_REPORTED",
        },
    ]
    checks = [q for q in qc if q["status"] in ("PASS", "REVIEW")]
    summary = {
        "samples": len(ids),
        "patients": len(ids),
        "genes": len(genes),
        "variants_input": len(records),
        "variants_retained": len(kept),
        "variant_genes": len(gene_counts),
        "median_purity": float(np.median(known_purity)) if known_purity else None,
        "purity_reported": len(known_purity),
        "samples_flagged": sum(m["qc_status"] == "REVIEW" for m in metrics),
        "subtypes": dict(sorted(Counter(m["subtype"] for m in metadata).items())),
        "qc_pass_rate_pct": 100 * sum(q["status"] == "PASS" for q in checks) / len(checks),
        "qc_rate_scope": "Assessable cohort checks only; unavailable fields are listed separately",
        "parameters": p,
        "pca_variance_explained_pct": payload["analysis"]["pca_variance_explained_pct"],
    }
    return clean_json(
        {
            "parameters": p,
            "sample_ids": ids,
            "summary": summary,
            "sample_metrics": metrics,
            "top_genes": top,
            "heatmap": heat,
            "variant_summary": variant_summary,
            "cna_summary": cna_summary,
            "subtype_summary": groups,
            "qc": qc,
            "retained_variant_ids": [v["variant_id"] for v in kept],
        }
    )
