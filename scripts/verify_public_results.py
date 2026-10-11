#!/usr/bin/env python3
"""Independent numeric checks against sklearn/scipy and the original pinned tables."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, zscore
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from oncoflow_atlas.analysis import compute_view  # noqa: E402
from oncoflow_atlas.public_data import prepare  # noqa: E402
from oncoflow_atlas.reporting import write_json  # noqa: E402
from oncoflow_atlas.validation import sha256  # noqa: E402


def verify(results: Path) -> dict:
    report = json.loads((results / "analysis.json").read_text())
    prepared = ROOT / "data/public/prepared"
    source = ROOT / "data/public/source"
    x = pd.read_csv(prepared / "expression.csv", index_col=0).loc[
        report["expression"]["genes"], report["expression"]["samples"]
    ]
    logs = np.log2(x.to_numpy() + 1)
    model = PCA(n_components=2, svd_solver="full")
    reference = model.fit_transform(logs.T)
    actual = np.array([[r["PC1"], r["PC2"]] for r in report["analysis"]["pca"]])
    for i in range(2):
        if np.dot(reference[:, i], actual[:, i]) < 0:
            reference[:, i] *= -1
    np.testing.assert_allclose(actual, reference, rtol=1e-9, atol=1e-9)
    explained = np.array(list(report["analysis"]["pca_variance_explained_pct"].values()))
    np.testing.assert_allclose(explained, model.explained_variance_ratio_ * 100, rtol=1e-10, atol=1e-10)
    top = report["default_view"]["heatmap"]
    indices = [x.index.get_loc(r["gene"]) for r in top]
    heat_reference = zscore(logs[indices], axis=1, ddof=1)
    heat_actual = np.array([r["z"] for r in top])
    np.testing.assert_allclose(heat_actual, heat_reference, rtol=1e-10, atol=1e-10)
    correlation = np.array([[pearsonr(a, b).statistic for b in logs.T] for a in logs.T])
    np.testing.assert_allclose(report["analysis"]["correlation"], correlation, rtol=1e-10, atol=1e-10)
    raw = pd.read_csv(source / "data_mutations.txt", sep="\t", comment="#", low_memory=False)
    raw_vaf = raw.t_alt_count / (raw.t_ref_count + raw.t_alt_count)
    actual_vaf = np.array([r["vaf"] for r in report["variants"]])
    np.testing.assert_allclose(np.sort(actual_vaf), np.sort(raw_vaf), rtol=1e-12, atol=1e-12)
    for row in report["sensitivity"]:
        assert row["retained_records"] == int((raw_vaf >= row["min_vaf"]).sum())
    with tempfile.TemporaryDirectory() as directory:
        tmp = Path(directory)
        rebuilt = prepare(source, ROOT / "config/public-source-lock.json", tmp / "prepared")
        for name, expected in rebuilt["prepared_files"].items():
            assert sha256(prepared / name) == expected, f"Non-reproducible prepared {name}"
        assert rebuilt == json.loads((prepared / "source_manifest.json").read_text())
        scenarios = []
        for path in [results, ROOT / "results/example"]:
            payload = json.loads((path / "analysis.json").read_text())
            options = [
                {},
                *[{"min_vaf": t} for t in [0, 0.1, 0.2, 0.3, 0.5, 1]],
                *[{"subtype": t} for t in sorted({m["subtype"] for m in payload["metadata"]})],
                {"search": payload["metadata"][0]["sample_id"]},
                {"search": "no-sample-matches"},
                {"top_n": 50, "min_purity": 0.8},
                {"min_vaf": 0.2, "top_n": 5},
            ]
            # Keep each payload once and have Node expand references; avoids a large repeated fixture.
            views = [{"options": opt, "expected": compute_view(payload, opt)} for opt in options]
            scenarios.append({"report": payload, "views": views})
        fixture = tmp / "parity.json"
        write_json(fixture, scenarios)
        completed = subprocess.run(
            ["node", str(ROOT / "tests/parity.cjs"), str(fixture)], check=True, text=True, capture_output=True
        )
        parity = json.loads(completed.stdout)
    output = {
        "status": "passed",
        "scope": "Independent numerical and byte-reproduction checks; not external biological validation",
        "libraries": {n: version(n) for n in ["numpy", "pandas", "scipy", "scikit-learn"]},
        "public_analysis_sha256": sha256(results / "analysis.json"),
        "source_revision": report["source"]["revision"],
        "samples": len(x.columns),
        "expression_genes": len(x),
        "variants": len(raw),
        "pca_coordinates_max_abs_error": float(np.max(np.abs(actual - reference))),
        "pca_explained_pct_max_abs_error": float(np.max(np.abs(explained - model.explained_variance_ratio_ * 100))),
        "heatmap_max_abs_error": float(np.max(np.abs(heat_actual - heat_reference))),
        "correlation_max_abs_error": float(np.max(np.abs(np.asarray(report["analysis"]["correlation"]) - correlation))),
        "vaf_max_abs_error": float(np.max(np.abs(np.sort(actual_vaf) - np.sort(raw_vaf)))),
        "prepared_files_reproduced": list(rebuilt["prepared_files"]),
        "threshold_checks": len(report["sensitivity"]),
        "python_javascript_parity": parity,
    }
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results/public")
    parser.add_argument("--output", type=Path, default=ROOT / "validation/independent-checks.json")
    args = parser.parse_args()
    result = verify(args.results)
    write_json(args.output, result)
    print(json.dumps(result, indent=2))
