"""Prepare one pinned, attributed cBioPortal/TCGA cohort without fabricating fields."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .analysis import fingerprint
from .reporting import write_json
from .validation import InputError, load_cohort, sha256


def cases(path: Path) -> set[str]:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("case_list_ids:"):
            return set(line.split(":", 1)[1].split())
    raise InputError(f"No case_list_ids in {path.name}")


def clean_profile(frame: pd.DataFrame, samples: list[str], name: str) -> tuple[pd.DataFrame, list[dict]]:
    symbols = frame.Hugo_Symbol.fillna("").astype(str)
    entrez = pd.to_numeric(frame.Entrez_Gene_Id, errors="coerce")
    valid_symbols = symbols.str.strip().ne("") & ~symbols.isin(["?", "NA", "NaN", "nan"])
    valid_ids = entrez.notna() & (entrez > 0) & (entrez % 1 == 0)
    eligible = valid_symbols & valid_ids
    duplicated_symbols = symbols[eligible].duplicated(keep=False).reindex(frame.index, fill_value=False)
    duplicated_ids = entrez[eligible].duplicated(keep=False).reindex(frame.index, fill_value=False)
    keep = eligible & ~duplicated_symbols & ~duplicated_ids
    rejected = []
    for i in frame.index[~keep]:
        reasons = []
        if not valid_symbols.loc[i]:
            reasons.append("missing_gene_symbol")
        if not valid_ids.loc[i]:
            reasons.append("missing_or_invalid_entrez_id")
        if duplicated_symbols.loc[i]:
            reasons.append("ambiguous_gene_symbol")
        if duplicated_ids.loc[i]:
            reasons.append("ambiguous_entrez_id")
        rejected.append(
            {
                "assay": name,
                "source_data_row": int(i + 1),
                "symbol": symbols.loc[i],
                "entrez": int(entrez.loc[i]) if valid_ids.loc[i] else None,
                "reasons": ";".join(reasons),
            }
        )
    selected = frame.loc[keep, samples].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(selected.to_numpy()).all():
        raise InputError(f"Nonfinite {name} values require an explicit preprocessing policy")
    selected.insert(0, "gene", symbols[keep].values)
    return selected.sort_values("gene", kind="stable").reset_index(drop=True), rejected


def prepare(source: Path, lock_path: Path, out: Path) -> dict:
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    for item in lock["files"]:
        file = source / item["path"]
        if not file.is_file() or file.stat().st_size != item["bytes"] or sha256(file) != item["sha256"]:
            raise InputError(f"Pinned source check failed: {item['path']}")
    read = lambda filename: pd.read_csv(source / filename, sep="\t", comment="#", low_memory=False)  # noqa: E731
    raw_expr = read("data_mrna_seq_v2_rsem.txt")
    raw_cna = read("data_cna.txt")
    clinical = read("data_clinical_sample.txt")
    panel = read("data_gene_panel_matrix.txt")
    raw_variants = read("data_mutations.txt")
    measured = {
        "clinical": set(clinical.SAMPLE_ID),
        "rna_matrix": set(raw_expr.columns[2:]),
        "rna_cases": cases(source / "case_lists/cases_RNA_Seq_v2_mRNA.txt"),
        "mutation_cases": cases(source / "case_lists/cases_sequenced.txt"),
        "mutation_coverage": set(panel.loc[panel.mutations == "WXS", "SAMPLE_ID"]),
        "cna_matrix": set(raw_cna.columns[2:]),
        "cna_cases": cases(source / "case_lists/cases_cna.txt"),
    }
    sample_ids = sorted(set.intersection(*measured.values()))
    if not sample_ids:
        raise InputError("No jointly profiled samples")
    if clinical.SAMPLE_ID.duplicated().any():
        raise InputError("Duplicate clinical sample identifiers")
    clinical = clinical.set_index("SAMPLE_ID").loc[sample_ids]
    if clinical.PATIENT_ID.duplicated().any():
        raise InputError("This protocol requires one sample per patient")
    expr, expr_excluded = clean_profile(raw_expr, sample_ids, "expression")
    cna, cna_excluded = clean_profile(raw_cna, sample_ids, "cna")
    if not np.isin(cna.iloc[:, 1:].to_numpy(), [-2, -1, 0, 1, 2]).all():
        raise InputError("Source CNA must contain discrete GISTIC states")
    cna.iloc[:, 1:] = cna.iloc[:, 1:].astype(int)
    metadata = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "patient_id": clinical.PATIENT_ID.values,
            "subtype": clinical.CANCER_TYPE_DETAILED.values,
            "grade": clinical.GRADE.fillna("Not reported").values,
            "purity": [None] * len(sample_ids),
            "batch": ["Not reported"] * len(sample_ids),
        }
    )
    selected = raw_variants[raw_variants.Tumor_Sample_Barcode.isin(sample_ids)].copy()
    keys = [
        "Tumor_Sample_Barcode",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
    ]
    if selected.duplicated(keys).any():
        raise InputError("Repeated source mutation loci require explicit adjudication")
    ref = pd.to_numeric(selected.t_ref_count, errors="coerce")
    alt = pd.to_numeric(selected.t_alt_count, errors="coerce")
    if (
        (ref.dropna() < 0).any()
        or (alt.dropna() < 0).any()
        or (ref.dropna() % 1 != 0).any()
        or (alt.dropna() % 1 != 0).any()
        or not np.isfinite(ref.dropna()).all()
        or not np.isfinite(alt.dropna()).all()
    ):
        raise InputError("Invalid source read counts")
    depth = ref + alt
    vaf = alt / depth.where(depth > 0)
    variants = pd.DataFrame(
        {
            "sample_id": selected.Tumor_Sample_Barcode,
            "gene": selected.Hugo_Symbol,
            "variant_class": selected.Variant_Classification,
            "vaf": vaf,
            "impact": selected.IMPACT.fillna("UNKNOWN"),
            "variant_id": ["tcga-" + fingerprint(row)[:20] for row in selected[keys].astype(str).values.tolist()],
            "chromosome": selected.Chromosome.astype(str),
            "position": selected.Start_Position,
            "ref": selected.Reference_Allele,
            "alt": selected.Tumor_Seq_Allele2,
            "ref_count": ref,
            "alt_count": alt,
        }
    ).sort_values(["sample_id", "variant_id"], kind="stable")
    out.mkdir(parents=True, exist_ok=True)
    for filename, table in [
        ("expression.csv", expr),
        ("cna.csv", cna),
        ("metadata.csv", metadata),
        ("variants.csv", variants),
    ]:
        table.to_csv(out / filename, index=False, float_format="%.17g")
    exclusions = expr_excluded + cna_excluded
    pd.DataFrame(exclusions).to_csv(out / "excluded_features.csv", index=False)
    source_manifest = {
        "kind": "public",
        "name": "TCGA · Cholangiocarcinoma",
        "study_id": lock["study_id"],
        "repository": lock["repository"],
        "revision": lock["revision"],
        "url": "https://www.cbioportal.org/study/summary?id=" + lock["study_id"],
        "source_url": f"https://github.com/{lock['repository']}/tree/{lock['revision']}/public/{lock['study_id']}",
        "citation": "The Cancer Genome Atlas, PanCancer Atlas (Cell, 2018); cBioPortal Datahub.",
        "data_terms": "ODbL-1.0 for the adapted database; retain the upstream TCGA study notice and attribution.",
        "expression_units": "RSEM normalized abundance, upstream batch normalized, linear scale",
        "cna_units": "Discrete GISTIC calls -2,-1,0,1,2",
        "genome_build": sorted(str(v) for v in selected.NCBI_Build.dropna().unique()),
        "samples": len(sample_ids),
        "patients": len(sample_ids),
        "assay_sample_counts": {key: len(value) for key, value in measured.items()},
        "excluded_samples": sorted(set.union(*measured.values()) - set(sample_ids)),
        "feature_counts": {
            "expression_source": len(raw_expr),
            "expression_retained": len(expr),
            "cna_source": len(raw_cna),
            "cna_retained": len(cna),
        },
        "excluded_feature_counts": {"expression": len(expr_excluded), "cna": len(cna_excluded)},
        "variant_records_source": len(raw_variants),
        "variant_records_retained": len(variants),
        "variants_outside_matched_cohort": len(raw_variants) - len(variants),
        "vaf_unavailable": int(vaf.isna().sum()),
        "purity_available": 0,
        "transformations": [
            "Exact sample intersection across all assay coverage sources",
            "Exclude all missing/ambiguous symbol or Entrez rows; no arbitrary duplicate selection",
            "VAF=t_alt_count/(t_ref_count+t_alt_count); retain unknown as missing",
            "Preserve source consequence classes including synonymous/noncoding records",
            "Purity and batch not reported in selected source fields; retain missing explicitly",
        ],
        "source_lock_sha256": sha256(lock_path),
        "source_files": [
            {"path": item["path"], "sha256": item["sha256"], "bytes": item["bytes"]} for item in lock["files"]
        ],
        "prepared_files": {
            name: sha256(out / name)
            for name in ["expression.csv", "cna.csv", "metadata.csv", "variants.csv", "excluded_features.csv"]
        },
        "limitations": [
            "Processed tables only; sequencing read QC and variant calling were not rerun",
            "36 tumors from one study; no normal comparator or external validation cohort",
            "PCA and mutation/CNA summaries are descriptive; no diagnostic classifier or survival model is fitted",
            "Variant counts include all source consequence classes and are not TMB per megabase",
        ],
    }
    write_json(out / "source_manifest.json", source_manifest)
    load_cohort(out / "expression.csv", out / "variants.csv", out / "metadata.csv", out / "cna.csv")
    return source_manifest
