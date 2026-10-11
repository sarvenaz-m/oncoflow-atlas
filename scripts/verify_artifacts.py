#!/usr/bin/env python3
"""Check that delivered results, browser evidence and current source agree."""
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    checked = 0
    for name in ["public", "example"]:
        folder = ROOT / "results" / name
        manifest = json.loads((folder / "run_manifest.json").read_text())
        for filename, item in manifest["outputs"].items():
            assert sha(folder / filename) == item["sha256"], f"Stale result: {name}/{filename}"
            checked += 1
        for filename, digest in manifest["software_sha256"].items():
            assert sha(ROOT / "src/oncoflow_atlas" / filename) == digest, f"Rebuild results after changing {filename}"
            checked += 1
    lock = json.loads((ROOT / "config/public-source-lock.json").read_text())
    for item in lock["files"]:
        assert sha(ROOT / "data/public/source" / item["path"]) == item["sha256"], item["path"]
        checked += 1
    evidence = json.loads((ROOT / "validation/independent-checks.json").read_text())
    assert evidence["status"] == "passed"
    assert sha(ROOT / "results/public/analysis.json") == evidence["public_analysis_sha256"]
    browser = json.loads((ROOT / "validation/browser-checks.json").read_text())
    assert browser["status"] == "passed"
    assert sha(ROOT / "docs/index.html") == browser["docs_sha256"]
    html = (ROOT / "docs/index.html").read_text()
    embedded = json.loads(re.search(r'<script id="oncoflow-data" type="application/json">(.*?)</script>', html, re.S)[1])
    for report, name in zip(embedded, ["public", "example"], strict=True):
        assert report == json.loads((ROOT / "results" / name / "analysis.json").read_text()), name
    print(json.dumps({"status": "passed", "hash_checks": checked + 2, "embedded_datasets": len(embedded)}))


if __name__ == "__main__":
    main()
