"""Tests for config loading and validation.

These also serve as the guard on the project's central design rule: the
shipped template must contain no real target. If someone fills in
``configs/template.yaml`` with their actual project instead of copying it,
``test_template_has_no_target`` fails.
"""

from __future__ import annotations

import textwrap

import pytest

from phylodonor import config as config_mod

REPO_ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]
TEMPLATE = REPO_ROOT / "configs" / "template.yaml"


def write(tmp_path, body: str):
    path = tmp_path / "cfg.yaml"
    path.write_text(textwrap.dedent(body))
    return path


def test_template_loads():
    cfg = config_mod.load(TEMPLATE, repo_root=REPO_ROOT)
    assert cfg.name == "my_target"


def test_template_has_no_target():
    """The shipped template must stay a template."""
    cfg = config_mod.load(TEMPLATE, repo_root=REPO_ROOT)
    assert cfg.get("target.accession") == "P00000"
    assert cfg.get("domains.definitions") == []
    assert cfg.get("search.database") is None
    assert cfg.get("grafting.candidates.exclude_positions") == []


@pytest.mark.parametrize("cfg_name", ["template.yaml", "demo_ecoli_phrb.yaml"])
def test_shipped_configs_do_not_weight_unvalidated_signals(cfg_name):
    """The likelihood path is implemented but unvalidated in this pipeline.

    Shipped configs must not give it weight, and must not hard-filter on it:
    a filter on an unvalidated signal discards designs silently, which is the
    one failure you cannot see in the output. If this test is failing because
    the PLL path has since been validated, update docs/grafting.md's status
    table in the same change -- the defaults and that table must agree.
    """
    cfg = config_mod.load(REPO_ROOT / "configs" / cfg_name, repo_root=REPO_ROOT)
    weights = cfg.get("grafting.score.weights")
    for term in ("delta_pll", "junction_pll"):
        assert weights.get(term, 0) == 0, f"{cfg_name} weights unvalidated term {term!r}"
    assert cfg.get("grafting.score.filters.max_pll_drop") is None


def test_rejects_bad_project_name(tmp_path):
    path = write(
        tmp_path,
        """
        project: {name: "My Target!"}
    """,
    )
    with pytest.raises(config_mod.ConfigError, match="lowercase slug"):
        config_mod.load(path)


def test_rejects_weights_that_do_not_sum_to_one(tmp_path):
    path = write(
        tmp_path,
        """
        project: {name: t}
        donors:
          weights: {sequence_identity: 0.5, domain_tmscore: 0.9}
    """,
    )
    with pytest.raises(config_mod.ConfigError, match="must sum to 1.0"):
        config_mod.load(path)


def test_accepts_weights_that_sum_to_one(tmp_path):
    path = write(
        tmp_path,
        """
        project: {name: t}
        donors:
          weights: {sequence_identity: 0.4, domain_tmscore: 0.6}
    """,
    )
    assert config_mod.load(path).get("donors.weights.domain_tmscore") == 0.6


@pytest.mark.parametrize(
    "start,end,match",
    [
        (0, 10, "1-based"),
        (50, 10, "end before start"),
    ],
)
def test_rejects_bad_ranges(tmp_path, start, end, match):
    path = write(
        tmp_path,
        f"""
        project: {{name: t}}
        domains:
          definitions:
            - {{name: d, start: {start}, end: {end}}}
    """,
    )
    with pytest.raises(config_mod.ConfigError, match=match):
        config_mod.load(path)


def test_missing_key_raises_with_path_in_message(tmp_path):
    cfg = config_mod.load(write(tmp_path, "project: {name: t}\n"))
    with pytest.raises(config_mod.ConfigError, match="donors.min_ogt_delta"):
        cfg.get("donors.min_ogt_delta")


def test_get_default_allows_none(tmp_path):
    cfg = config_mod.load(write(tmp_path, "project: {name: t}\n"))
    assert cfg.get("nothing.here", None) is None


def test_require_rejects_empty_collections(tmp_path):
    cfg = config_mod.load(
        write(
            tmp_path,
            """
        project: {name: t}
        holo: {cofactors: []}
    """,
        )
    )
    with pytest.raises(config_mod.ConfigError, match="must be set"):
        cfg.require("holo.cofactors")


def test_check_stage_lists_all_missing_keys(tmp_path):
    cfg = config_mod.load(write(tmp_path, "project: {name: t}\n"))
    with pytest.raises(config_mod.ConfigError) as exc:
        cfg.check_stage("homologs")
    assert "target.fasta" in str(exc.value)
    assert "search.database" in str(exc.value)
