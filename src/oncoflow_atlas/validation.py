"""Strict structural/numeric validation, kept separate from analysis."""

from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


class InputError(ValueError):
    """The input cannot be interpreted under the declared data contract."""


def sha256(path: str | Path) -> str:
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.stat().st_size > 150 * 1024 * 1024:
        raise InputError(f"Input exceeds 150 MiB: {path.name}")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader, [])
        if not header or len(header) != len(set(header)) or any(not h.strip() for h in header):
            raise InputError(f"Empty or duplicate column names in {path.name}")
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise InputError(f"Ragged CSV row {line} in {path.name}")
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def identifiers(values: pd.Series, label: str, unique: bool = True) -> None:
    if values.isna().any() or values.astype(str).str.strip().eq("").any():
        raise InputError(f"{label} contains a missing identifier")
    if (values.astype(str) != values.astype(str).str.strip()).any():
        raise InputError(f"{label} has leading/trailing identifier whitespace")
    if unique and values.duplicated().any():
        raise InputError(f"{label} contains duplicated identifiers")


def numeric(values: pd.Series, label: str, nullable: bool = False) -> pd.Series:
    missing = values.astype(str).isin(["", "NA", "NaN", "nan", "null"])
    # Python/NumPy's decimal parser preserves a round-tripped IEEE float. The
    # pandas fast parser can place an exact boundary (e.g. 0.3) one ULP below it.
    try:
        parsed = pd.Series(np.asarray(values.where(~missing, np.nan), dtype=float), index=values.index)
    except (ValueError, TypeError) as error:
        raise InputError(f"{label} must contain finite numeric values") from error
    invalid = (~missing & parsed.isna()) | (~parsed.isna() & ~np.isfinite(parsed))
    if invalid.any() or (not nullable and parsed.isna().any()):
        raise InputError(
            f"{label} must contain finite numeric values" + (" or explicit missing values" if nullable else "")
        )
    return parsed


@dataclass
class Cohort:
    expression: pd.DataFrame
    variants: pd.DataFrame
    metadata: pd.DataFrame
    cna: pd.DataFrame | None
    inputs: dict


def load_cohort(
    expression: str | Path, variants: str | Path, metadata: str | Path, cna: str | Path | None = None
) -> Cohort:
    paths = {"expression": expression, "variants": variants, "metadata": metadata}
    if cna:
        paths["cna"] = cna
    tables = {name: read_table(path) for name, path in paths.items()}
    expr, var, meta = (tables[k] for k in ("expression", "variants", "metadata"))
    for label, frame, required in [
        ("expression", expr, {"gene"}),
        ("metadata", meta, {"sample_id", "subtype"}),
        ("variants", var, {"sample_id", "gene", "variant_class", "vaf", "impact"}),
    ]:
        absent = required - set(frame.columns)
        if absent:
            raise InputError(f"{label} is missing columns: {sorted(absent)}")
    if not 3 <= len(meta) <= 500 or not 2 <= len(expr) <= 60000:
        raise InputError("Expected 3–500 samples and 2–60,000 expression features")
    identifiers(meta.sample_id, "metadata sample_id")
    identifiers(meta.subtype, "metadata subtype", unique=False)
    if "patient_id" not in meta:
        meta["patient_id"] = meta.sample_id
    identifiers(meta.patient_id, "metadata patient_id")
    for name in ("purity",):
        if name not in meta:
            meta[name] = ""
        meta[name] = numeric(meta[name], name, nullable=True)
        if not meta[name].dropna().between(0, 1).all():
            raise InputError(f"{name} must be between 0 and 1")
    if "batch" not in meta:
        meta["batch"] = "Not reported"
    samples = meta.sample_id.tolist()
    sample_set = set(samples)
    for column in ("sample_id", "gene", "variant_class", "impact"):
        identifiers(var[column], f"variants {column}", unique=False)
    if not set(var.sample_id).issubset(sample_set):
        raise InputError("Variant rows reference unknown samples")
    var["vaf"] = numeric(var.vaf, "VAF", nullable=True)
    if not var.vaf.dropna().between(0, 1).all():
        raise InputError("VAF must be between 0 and 1")
    if "variant_id" not in var:
        var["variant_id"] = [f"record-{i + 1:06d}" for i in range(len(var))]
    identifiers(var.variant_id, "variant_id")
    if var.drop(columns="variant_id").duplicated().any():
        raise InputError("Duplicated variant records must be resolved before analysis")
    matrices = {}
    for name in ("expression", "cna"):
        if name not in tables:
            continue
        frame = tables[name]
        if name == "cna" and not 1 <= len(frame) <= 60000:
            raise InputError("Expected 1–60,000 CNA features")
        if "gene" not in frame or set(frame.columns) - {"gene"} != sample_set:
            raise InputError(f"{name} sample columns must exactly match metadata samples")
        identifiers(frame.gene, f"{name} gene")
        values = frame.set_index("gene")[samples].apply(lambda s, label=name: numeric(s, f"{label} values"))
        if name == "expression" and (values.to_numpy() < 0).any():
            raise InputError("Expression values must be nonnegative on the declared linear scale")
        if name == "cna" and not np.isin(values.to_numpy(), [-2, -1, 0, 1, 2]).all():
            raise InputError("CNA values must be discrete -2, -1, 0, 1, 2 calls")
        matrices[name] = values.astype(float).sort_index()
    if (matrices["expression"].sum(axis=0) == 0).any():
        raise InputError("An expression sample has zero total signal")
    inputs = {k: {"name": Path(p).name, "sha256": sha256(p), "bytes": Path(p).stat().st_size} for k, p in paths.items()}
    return Cohort(matrices["expression"], var, meta, matrices.get("cna"), inputs)


def validation_summary(cohort: Cohort) -> dict:
    return {
        "status": "passed",
        "samples": len(cohort.metadata),
        "patients": len(cohort.metadata),
        "expression_genes": len(cohort.expression),
        "variant_records": len(cohort.variants),
        "cna_genes": len(cohort.cna) if cohort.cna is not None else None,
        "missing_expression_cells": 0,
        "purity_unavailable": int(cohort.metadata.purity.isna().sum()),
        "vaf_unavailable": int(cohort.variants.vaf.isna().sum()),
        "checks": [
            "exact_sample_alignment",
            "unique_sample_patient_gene_ids",
            "rectangular_tables",
            "finite_nonnegative_expression",
            "nullable_bounded_purity_vaf",
            "unique_variant_records",
            "discrete_cna_calls",
            "nonzero_sample_signal",
        ],
        "inputs": cohort.inputs,
        "scope": "Tabular contract validation; assay-level sequencing quality is not inferred from processed tables",
    }
