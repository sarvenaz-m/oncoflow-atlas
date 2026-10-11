"""Validate inputs, calculate results, and render a reproducible report."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from . import __version__
from .analysis import build_payload
from .reporting import render_report, write_json
from .validation import InputError, load_cohort, sha256, validation_summary

TABLES = {
    "sample_metrics.csv": "sample_metrics",
    "top_genes.csv": "top_genes",
    "variant_summary.csv": "variant_summary",
    "cna_summary.csv": "cna_summary",
    "subtype_summary.csv": "subtype_summary",
    "qc_summary.csv": "qc",
}


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", nargs="?", choices=["validate", "analyse", "render", "run"], default="run")
    p.add_argument("--expression")
    p.add_argument("--variants")
    p.add_argument("--metadata")
    p.add_argument("--cna")
    p.add_argument("--source-manifest")
    p.add_argument("--validation", help="Verify input hashes against this separately generated validation artifact")
    p.add_argument("--analysis", help="Existing analysis.json for the independent render stage")
    p.add_argument("--outdir", required=True)
    p.add_argument("--cohort-name", default="cohort")
    p.add_argument("--min-vaf", type=float, default=0.05)
    p.add_argument("--min-purity", type=float, default=0.4)
    p.add_argument("--top-n", type=int, default=12)
    return p


def csv_frame(path: Path, frame: pd.DataFrame) -> None:
    # Quote spreadsheet formula prefixes in textual display fields, never numeric cells.
    frame = frame.copy()
    for col in frame.select_dtypes(include=["object", "string"]):
        frame[col] = frame[col].map(
            lambda value: (
                "'" + value if isinstance(value, str) and value.startswith(("=", "+", "-", "@", "\t", "\r")) else value
            )
        )
    frame.to_csv(path, index=False, float_format="%.15g")


def finalize(out: Path, payload: dict) -> None:
    names = [
        "analysis.json",
        "validation.json",
        "cohort_summary.json",
        "pca_coordinates.csv",
        "pca_loadings.csv",
        "sample_correlation.csv",
        "sensitivity.csv",
        *TABLES.keys(),
        "report.html",
    ]
    manifest = {
        "workflow": {"name": "OncoFlow Atlas", "version": __version__},
        "cohort": payload["cohort"],
        "generated_utc": datetime.now(UTC).isoformat(),
        "data_fingerprint": payload["data_fingerprint"],
        "parameters": payload["parameters"],
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": version("numpy"),
            "pandas": version("pandas"),
        },
        "inputs": payload["inputs"],
        "software_sha256": {
            str(p.relative_to(Path(__file__).parent)): sha256(p)
            for p in sorted(Path(__file__).parent.rglob("*"))
            if p.is_file() and p.suffix in {".py", ".js", ".html", ".css"}
        },
        "outputs": {
            n: {"bytes": (out / n).stat().st_size, "sha256": sha256(out / n)} for n in names if (out / n).is_file()
        },
    }
    write_json(out / "run_manifest.json", manifest)


def main(argv: list[str] | None = None) -> None:
    args = parser().parse_args(argv)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    try:
        if args.stage == "render":
            if not args.analysis:
                raise InputError("render requires --analysis")
            payload = json.loads(Path(args.analysis).read_text(encoding="utf-8"))
            if payload.get("schema_version") != "2.0":
                raise InputError("Expected a schema-2.0 analysis report")
            render_report([payload], out / "report.html")
            finalize(out, payload)
            print(json.dumps({"stage": "render", "report": str(out / "report.html")}))
            return
        if not all([args.expression, args.variants, args.metadata]):
            raise InputError("Provide --expression, --variants and --metadata")
        cohort = load_cohort(args.expression, args.variants, args.metadata, args.cna)
        validation = validation_summary(cohort)
        if args.validation:
            prior = json.loads(Path(args.validation).read_text(encoding="utf-8"))
            if prior.get("status") != "passed" or prior.get("inputs") != validation["inputs"]:
                raise InputError("Inputs changed after validation; refuse stale validation artifact")
        write_json(out / "validation.json", validation)
        if args.stage == "validate":
            print(json.dumps({"stage": "validate", "status": "passed", "samples": len(cohort.metadata)}))
            return
        source = json.loads(Path(args.source_manifest).read_text(encoding="utf-8")) if args.source_manifest else None
        if source and "prepared_files" in source:
            for key, item in cohort.inputs.items():
                expected = source["prepared_files"].get(item["name"])
                if expected != item["sha256"]:
                    raise InputError(f"Source manifest does not match prepared {key} input")
        payload, loadings = build_payload(
            cohort,
            args.cohort_name,
            {"min_vaf": args.min_vaf, "top_n": args.top_n, "min_purity": args.min_purity},
            source,
        )
        write_json(out / "analysis.json", payload)
        view = payload["default_view"]
        write_json(out / "cohort_summary.json", {"cohort": args.cohort_name, **view["summary"]})
        for filename, key in TABLES.items():
            frame = pd.DataFrame(view[key])
            if frame.empty:
                columns = {
                    "variant_summary": ["gene", "records", "samples"],
                    "cna_summary": ["gene", "amplified_samples", "deep_deleted_samples"],
                }
                frame = pd.DataFrame(columns=columns.get(key, []))
            if key == "sample_metrics":
                frame["qc_reasons"] = frame.qc_reasons.map(lambda items: "; ".join(items))
            csv_frame(out / filename, frame)
        csv_frame(out / "pca_coordinates.csv", pd.DataFrame(payload["analysis"]["pca"]))
        csv_frame(out / "pca_loadings.csv", loadings)
        samples = payload["expression"]["samples"]
        correlation = pd.DataFrame(payload["analysis"]["correlation"], columns=samples)
        correlation.insert(0, "sample_id", samples)
        csv_frame(out / "sample_correlation.csv", correlation)
        csv_frame(out / "sensitivity.csv", pd.DataFrame(payload["sensitivity"]))
        if args.stage == "run":
            render_report([payload], out / "report.html")
        finalize(out, payload)
        print(
            json.dumps(
                {
                    "stage": args.stage,
                    "samples": len(samples),
                    "genes": len(cohort.expression),
                    "retained_variants": view["summary"]["variants_retained"],
                    "outdir": str(out),
                }
            )
        )
    except (InputError, OSError, ValueError) as error:
        raise SystemExit(f"Input/analysis error: {error}") from error


if __name__ == "__main__":
    main()
