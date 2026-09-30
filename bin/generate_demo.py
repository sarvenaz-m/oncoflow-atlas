#!/usr/bin/env python3
"""Generate a deterministic synthetic multi-omic cohort for local smoke testing."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

GENES = [
    "KRAS", "TP53", "SMAD4", "CDKN2A", "MYC", "GATA6", "KRT19", "MUC1",
    "EPCAM", "VIM", "COL1A1", "CXCL8", "MKI67", "BRCA2", "PALB2", "KRT8",
    "KRT18", "KRT7", "S100A2", "LAMC2", "KRT17", "MMP7", "CEACAM6", "TFF1",
    "TFF2", "AGR2", "MUC5AC", "SPP1", "COL3A1", "FN1", "ITGA2", "TGFB1",
    "CXCL1", "CXCL2", "STAT3", "ERBB2", "EGFR", "PIK3CA", "AKT1", "MTOR",
]

CLASSICAL_MARKERS = {"GATA6": 3.6, "TFF1": 2.4, "TFF2": 2.1, "AGR2": 2.0, "MUC5AC": 1.9}
BASAL_MARKERS = {"VIM": 3.2, "CXCL8": 2.8, "KRT17": 2.5, "S100A2": 2.3, "LAMC2": 2.1}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="data/demo")
    parser.add_argument("--samples", type=int, default=36)
    parser.add_argument("--seed", type=int, default=19)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.samples < 8:
        raise SystemExit("--samples must be at least 8 for a meaningful demonstration cohort")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    samples = [f"OF-{i:03d}" for i in range(1, args.samples + 1)]
    subtypes = np.array(["Classical" if i % 2 == 0 else "Basal-like" for i in range(args.samples)])
    batches = np.array([f"B{(i % 3) + 1}" for i in range(args.samples)])
    purity = np.clip(rng.normal(0.67, 0.11, args.samples), 0.36, 0.93)

    metadata = pd.DataFrame(
        {
            "sample_id": samples,
            "subtype": subtypes,
            "purity": np.round(purity, 3),
            "batch": batches,
        }
    )

    expression = pd.DataFrame({"gene": GENES})
    for j, sample in enumerate(samples):
        values = rng.lognormal(mean=3.1, sigma=0.65, size=len(GENES))
        markers = CLASSICAL_MARKERS if subtypes[j] == "Classical" else BASAL_MARKERS
        for gene, multiplier in markers.items():
            values[GENES.index(gene)] *= multiplier
        values *= rng.normal(1.0, 0.06)
        expression[sample] = np.round(values, 3)

    rows: list[dict[str, object]] = []
    hotspots = {"KRAS": 0.84, "TP53": 0.67, "SMAD4": 0.35, "CDKN2A": 0.30, "BRCA2": 0.12, "PALB2": 0.08}
    impacts = ["MODERATE", "HIGH"]
    for sample in samples:
        for gene, probability in hotspots.items():
            if rng.random() < probability:
                rows.append(
                    {
                        "sample_id": sample,
                        "gene": gene,
                        "variant_class": rng.choice(["SNV", "INDEL"], p=[0.82, 0.18]),
                        "vaf": round(float(rng.uniform(0.06, 0.68)), 3),
                        "impact": rng.choice(impacts, p=[0.62, 0.38]),
                    }
                )
        if rng.random() < 0.18:
            rows.append(
                {
                    "sample_id": sample,
                    "gene": rng.choice(["PIK3CA", "ERBB2", "EGFR", "AKT1"]),
                    "variant_class": "SNV",
                    "vaf": round(float(rng.uniform(0.05, 0.35)), 3),
                    "impact": "MODERATE",
                }
            )

    pd.DataFrame(rows, columns=["sample_id", "gene", "variant_class", "vaf", "impact"]).to_csv(
        out / "variants.csv", index=False
    )
    expression.to_csv(out / "expression.csv", index=False)
    metadata.to_csv(out / "metadata.csv", index=False)
    print(f"Wrote deterministic demo cohort with {args.samples} samples to {out}")


if __name__ == "__main__":
    main()
