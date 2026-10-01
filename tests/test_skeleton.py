"""Basic skeleton tests for Gate 0."""

from pathlib import Path

import yaml


def test_repo_structure():
    """Verify essential project files exist."""
    required_files = [
        "README.md",
        "PROTOCOL.md",
        "DEVIATIONS.md",
        "CLAIMS.md",
        "RELEASE_CHECKLIST.md",
        "LICENSE",
        "CITATION.cff",
        ".gitignore",
        "requirements.txt",
        "configs/protocol.yaml",
        "checklists/TRIPOD-AI.md",
        ".github/workflows/ci.yml",
    ]
    for rel_path in required_files:
        assert Path(rel_path).exists(), f"Missing required file: {rel_path}"


def test_protocol_config():
    """Verify protocol.yaml is valid YAML and has key fields."""
    with open("configs/protocol.yaml", "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    assert cfg["study"]["cohort"] == "TCGA-KIRC"
    assert "split_seeds" in cfg["random_seeds"]
    assert len(cfg["random_seeds"]["split_seeds"]) == 5
    assert cfg["models"]["m0"]["name"] == "clinical_cox"
    assert cfg["models"]["m1"]["name"] == "clinical_plus_rna_elasticnet"
