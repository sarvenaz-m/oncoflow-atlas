#!/usr/bin/env python3
"""Rebuild both report datasets and the offline GitHub Pages workbench."""

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from oncoflow_atlas.reporting import render_report  # noqa: E402


def main():
    for name, folder, output in [
        ("tcga-chol-pancan", "public/prepared", "public"),
        ("synthetic-demo", "demo", "example"),
    ]:
        data = ROOT / "data" / folder
        command = [
            sys.executable,
            str(ROOT / "bin/analyse_multiomics.py"),
            "--expression",
            str(data / "expression.csv"),
            "--variants",
            str(data / "variants.csv"),
            "--metadata",
            str(data / "metadata.csv"),
            "--source-manifest",
            str(data / "source_manifest.json"),
            "--cohort-name",
            name,
            "--outdir",
            str(ROOT / "results" / output),
        ]
        if (data / "cna.csv").exists():
            command += ["--cna", str(data / "cna.csv")]
        subprocess.run(command, check=True, env=os.environ | {"OPENBLAS_NUM_THREADS": "1"})
    reports = [json.loads((ROOT / "results" / name / "analysis.json").read_text()) for name in ["public", "example"]]
    render_report(reports, ROOT / "docs/index.html")
    print("Built docs/index.html and both standalone reports.")


if __name__ == "__main__":
    main()
