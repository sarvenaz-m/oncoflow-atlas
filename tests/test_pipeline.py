import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def run_demo(tmp_path: Path, samples: int = 12):
    data = tmp_path / "data"
    results = tmp_path / "results"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "bin/generate_demo.py"),
            "--outdir",
            str(data),
            "--samples",
            str(samples),
            "--seed",
            "19",
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "bin/analyse_multiomics.py"),
            "--expression",
            str(data / "expression.csv"),
            "--variants",
            str(data / "variants.csv"),
            "--metadata",
            str(data / "metadata.csv"),
            "--outdir",
            str(results),
            "--cohort-name",
            "ci-test",
            "--min-vaf",
            "0.08",
            "--top-n",
            "10",
        ],
        check=True,
    )
    return data, results


def test_demo_pipeline_contract(tmp_path):
    _, results = run_demo(tmp_path)
    summary = json.loads((results / "cohort_summary.json").read_text())
    manifest = json.loads((results / "run_manifest.json").read_text())

    assert summary["samples"] == 12
    assert summary["genes"] == 40
    assert summary["variants_retained"] <= summary["variants_input"]
    assert summary["parameters"]["min_vaf"] == 0.08
    assert summary["qc_pass_rate_pct"] >= 80

    expected = {
        "report.html",
        "cohort_summary.json",
        "run_manifest.json",
        "sample_metrics.csv",
        "pca_coordinates.csv",
        "subtype_summary.csv",
        "top_genes.csv",
        "variant_summary.csv",
        "qc_summary.csv",
    }
    assert expected.issubset({p.name for p in results.iterdir()})
    assert (results / "report.html").stat().st_size > 20_000
    assert "outputs" in manifest and "report.html" in manifest["outputs"]
    assert len(manifest["inputs"]["expression"]["sha256"]) == 64


def test_demo_is_deterministic(tmp_path):
    data_a = tmp_path / "a"
    data_b = tmp_path / "b"
    for out in [data_a, data_b]:
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "bin/generate_demo.py"),
                "--outdir",
                str(out),
                "--samples",
                "10",
                "--seed",
                "29",
            ],
            check=True,
        )
    assert (data_a / "expression.csv").read_bytes() == (data_b / "expression.csv").read_bytes()
    assert (data_a / "variants.csv").read_bytes() == (data_b / "variants.csv").read_bytes()


def test_invalid_purity_is_rejected(tmp_path):
    data = tmp_path / "data"
    results = tmp_path / "results"
    subprocess.run(
        [sys.executable, str(ROOT / "bin/generate_demo.py"), "--outdir", str(data), "--samples", "8"], check=True
    )
    metadata = pd.read_csv(data / "metadata.csv")
    metadata.loc[0, "purity"] = 1.4
    metadata.to_csv(data / "metadata.csv", index=False)

    proc = subprocess.run(
        [
            sys.executable,
            str(ROOT / "bin/analyse_multiomics.py"),
            "--expression",
            str(data / "expression.csv"),
            "--variants",
            str(data / "variants.csv"),
            "--metadata",
            str(data / "metadata.csv"),
            "--outdir",
            str(results),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "purity" in (proc.stdout + proc.stderr).lower()
