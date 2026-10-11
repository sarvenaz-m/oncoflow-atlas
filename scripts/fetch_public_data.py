"""Fetch pinned source bytes. Git LFS uses the upstream Datahub's configured endpoint."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def download(item: dict, root: Path, endpoint: str) -> str:
    destination = root / item["path"]
    if destination.is_file() and hashlib.sha256(destination.read_bytes()).hexdigest() == item["sha256"]:
        return item["path"]
    url = item["url"]
    headers = {}
    if item.get("lfs_oid_sha256"):
        body = {
            "operation": "download",
            "transfers": ["basic"],
            "objects": [{"oid": item["lfs_oid_sha256"], "size": item["bytes"]}],
        }
        req = urllib.request.Request(
            endpoint,
            data=json.dumps(body).encode(),
            headers={"Content-Type": "application/vnd.git-lfs+json", "Accept": "application/vnd.git-lfs+json"},
        )
        with urllib.request.urlopen(req, timeout=60) as response:
            obj = json.load(response)["objects"][0]
        if "actions" not in obj:
            raise ValueError("Upstream object unavailable: " + item["path"])
        action = obj["actions"]["download"]
        url, headers = action["href"], action.get("header", {})
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part")
    digest = hashlib.sha256()
    with (
        urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=90) as response,
        temporary.open("wb") as output,
    ):
        while chunk := response.read(1024 * 1024):
            digest.update(chunk)
            output.write(chunk)
    if digest.hexdigest() != item["sha256"] or temporary.stat().st_size != item["bytes"]:
        raise ValueError("Source hash/size mismatch: " + item["path"])
    temporary.replace(destination)
    return item["path"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", "--outdir", dest="output", type=Path, required=True)
    args = parser.parse_args()
    lock = json.loads((ROOT / "config/public-source-lock.json").read_text())
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for name in pool.map(lambda item: download(item, args.output, lock["lfs_batch_endpoint"]), lock["files"]):
            print("Verified " + name, flush=True)


if __name__ == "__main__":
    main()
