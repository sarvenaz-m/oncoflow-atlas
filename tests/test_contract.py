import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncoflow_atlas.analysis import build_payload, compute_view, pca_numpy
from oncoflow_atlas.cli import main
from oncoflow_atlas.reporting import render_report, write_json
from oncoflow_atlas.validation import InputError, load_cohort, numeric

ROOT = Path(__file__).resolve().parents[1]


def test_decimal_roundtrip_preserves_closed_vaf_boundary():
    value = 3 / 10
    parsed = numeric(pd.Series([format(value, ".17g")]), "VAF")
    assert parsed.iloc[0] == value
    assert parsed.iloc[0] >= 0.3


@pytest.mark.parametrize(
    "file,old,new,match",
    [
        ("expression", "gene,A,B,C", "gene,A,A,C", "duplicate"),
        ("expression", "gene,A,B,C", "gene,A,B,OTHER", "exactly match"),
        ("expression", "G1,1,3,7", "G1,1,3", "Ragged"),
        ("expression", "G1,1,3,7", "G1,1,,7", "finite numeric"),
        ("expression", "G1,1,3,7", "G1,1,inf,7", "finite numeric"),
        ("expression", "G1,1,3,7", "G1,1,-1,7", "nonnegative"),
        ("expression", "G2,7,3,1", "G1,7,3,1", "duplicated"),
        ("metadata", "B,P2", "B,P1", "duplicated"),
        ("metadata", "0.7", "1.1", "purity"),
        ("variants", "A,G1", "UNKNOWN,G1", "unknown samples"),
        ("variants", "0.5", "oops", "VAF"),
        ("variants", "0.5", "1.2", "VAF"),
        ("variants", "v3", "v1", "duplicated"),
        ("cna", "G1,2,0,-2", "G1,2,0.1,-2", "discrete"),
    ],
)
def test_invalid_inputs_fail_without_imputation(cohort_paths, file, old, new, match):
    path = cohort_paths[file]
    path.write_text(path.read_text().replace(old, new))
    with pytest.raises(InputError, match=match):
        load_cohort(**cohort_paths)


def test_view_math_and_missingness(cohort_paths):
    p, _ = build_payload(load_cohort(**cohort_paths), "manual")
    v = p["default_view"]
    assert v["top_genes"] == [
        {"gene": "G1", "variance": 1.0},
        {"gene": "G2", "variance": 1.0},
        {"gene": "G3", "variance": 0.0},
    ]
    assert v["heatmap"][0]["z"] == [-1, 0, 1]
    assert v["heatmap"][2]["z"] == [0, 0, 0]
    assert v["sample_metrics"][1]["purity"] is None
    assert v["sample_metrics"][1]["variants_missing_vaf"] == 1
    assert v["summary"]["variants_retained"] == 2
    assert compute_view(p, {"min_vaf": 0})["summary"]["variants_retained"] == 2
    boundary = compute_view(p, {"min_vaf": 0.5})
    assert boundary["retained_variant_ids"] == ["v1"]
    single = compute_view(p, {"search": " A "})
    assert single["sample_ids"] == ["A"]
    assert all(row["z"] == [0] for row in single["heatmap"])
    empty = compute_view(p, {"search": "absent"})
    assert empty["sample_ids"] == empty["top_genes"] == empty["cna_summary"] == []
    # Filtering must not claim to refit PCA.
    assert single["summary"]["pca_variance_explained_pct"] == p["analysis"]["pca_variance_explained_pct"]


def test_matrix_column_order_is_aligned(cohort_paths):
    initial, _ = build_payload(load_cohort(**cohort_paths), "cohort")
    for name in ["expression", "cna"]:
        frame = pd.read_csv(cohort_paths[name])
        frame[["gene", "C", "A", "B"]].to_csv(cohort_paths[name], index=False)
    aligned, _ = build_payload(load_cohort(**cohort_paths), "cohort")
    assert initial["default_view"] == aligned["default_view"]
    assert initial["analysis"] == aligned["analysis"]


def test_constant_cohort_has_no_pca():
    with pytest.raises(InputError, match="variable"):
        pca_numpy(np.full((36, 4), np.log2(6)))


def test_empty_summary_has_a_readable_csv_header(cohort_paths, tmp_path):
    flags = [str(v) for key, value in cohort_paths.items() if key != "cna" for v in ("--" + key, value)]
    out = tmp_path / "empty-results"
    main([*flags, "--min-vaf", "1", "--outdir", str(out)])
    assert pd.read_csv(out / "variant_summary.csv").columns.tolist() == ["gene", "records", "samples"]
    assert pd.read_csv(out / "cna_summary.csv").columns.tolist() == [
        "gene",
        "amplified_samples",
        "deep_deleted_samples",
    ]


def test_single_cna_gene_is_supported(cohort_paths):
    cohort_paths["cna"].write_text("gene,A,B,C\nG1,2,0,-2\n")
    assert len(load_cohort(**cohort_paths).cna) == 1


def test_stale_validation_refused(cohort_paths, tmp_path):
    out = tmp_path / "output"
    flags = [str(v) for key, value in cohort_paths.items() for v in ("--" + key, value)]
    main(["validate", *flags, "--outdir", str(out)])
    p = cohort_paths["expression"]
    p.write_text(p.read_text().replace("G1,1,3,7", "G1,1,3,8"))
    with pytest.raises(SystemExit, match="changed after validation"):
        main(["analyse", *flags, "--validation", str(out / "validation.json"), "--outdir", str(out)])
    assert not (out / "analysis.json").exists()


def test_separate_stages_match_direct_run(cohort_paths, tmp_path):
    flags = [str(v) for key, value in cohort_paths.items() for v in ("--" + key, value)]
    direct, staged = tmp_path / "direct", tmp_path / "staged"
    main([*flags, "--outdir", str(direct)])
    main(["validate", *flags, "--outdir", str(staged)])
    main(["analyse", *flags, "--validation", str(staged / "validation.json"), "--outdir", str(staged)])
    assert not (staged / "report.html").exists()
    main(["render", "--analysis", str(staged / "analysis.json"), "--outdir", str(staged)])
    for path in direct.iterdir():
        if path.name != "run_manifest.json":
            assert path.read_bytes() == (staged / path.name).read_bytes(), path.name


def test_render_escapes_embedded_script_and_is_offline(cohort_paths, tmp_path):
    payload, _ = build_payload(load_cohort(**cohort_paths), '</script><script>alert("x")</script>')
    out = tmp_path / "report.html"
    render_report([payload], out)
    html = out.read_text()
    assert payload["cohort"] not in html
    assert "\\u003c/script\\u003e" in html
    assert "connect-src 'none'" in html
    assert "<script src=" not in html


def test_javascript_contract_on_boundary_fixture(cohort_paths, tmp_path):
    payload, _ = build_payload(load_cohort(**cohort_paths), "manual")
    inputs = []
    for options in [{}, {"min_vaf": 0}, {"min_vaf": 0.5}, {"search": "B"}, {"search": "absent"}, {"min_purity": 1}]:
        inputs.append({"report": payload, "options": options, "expected": compute_view(payload, options)})
    fixture = tmp_path / "parity.json"
    write_json(fixture, inputs)
    run = subprocess.run(["node", str(ROOT / "tests/parity.cjs"), str(fixture)], text=True, capture_output=True)
    assert run.returncode == 0, run.stderr
    result = json.loads(run.stdout)
    assert result["scenarios"] == 6
