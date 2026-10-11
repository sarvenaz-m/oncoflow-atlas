"""Normalize the pinned public cohort and write its full selection ledger."""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from oncoflow_atlas.public_data import prepare  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "data/public/source")
    parser.add_argument("--output", type=Path, default=ROOT / "data/public/prepared")
    args = parser.parse_args()
    result = prepare(args.source, ROOT / "config/public-source-lock.json", args.output)
    print(
        json.dumps(
            {
                k: result[k]
                for k in ["samples", "patients", "feature_counts", "variant_records_retained", "vaf_unavailable"]
            }
        )
    )


if __name__ == "__main__":
    main()
