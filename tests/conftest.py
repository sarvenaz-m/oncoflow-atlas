from pathlib import Path

import pytest


@pytest.fixture
def cohort_paths(tmp_path):
    files = {
        "expression": "gene,A,B,C\nG1,1,3,7\nG2,7,3,1\nG3,2,2,2\n",
        "metadata": "sample_id,patient_id,subtype,purity\nA,P1,X,0.7\nB,P2,X,\nC,P3,Y,0.2\n",
        "variants": "sample_id,gene,variant_class,vaf,impact,variant_id\nA,G1,Missense,0.5,MODERATE,v1\nB,G2,Silent,,LOW,v2\nC,G1,Missense,0.1,MODERATE,v3\n",
        "cna": "gene,A,B,C\nG1,2,0,-2\nG2,1,-1,0\n",
    }
    paths = {}
    for name, content in files.items():
        path: Path = tmp_path / (name + ".csv")
        path.write_text(content)
        paths[name] = path
    return paths
