"""gecko3-yaml-import tests — YAML GEM vs ecModel auto-detection.

Design Ref: gecko3-yaml-import.plan.md FR-2 / FR-3. Ensures ``load_ecmodel``
sets ``is_ecmodel`` based on the presence of an ``ec-rxns`` section and that
plain GEM YAMLs load cleanly into a fresh project.
"""

from __future__ import annotations

import cobra
import pytest

from cnapy.ecmodel.ecmodel_builder import build_ecmodel
from cnapy.ecmodel.ecmodel_data import ECModelData
from cnapy.ecmodel.exceptions import EcYamlError
from cnapy.ecmodel.yaml_io import load_ecmodel, save_ecmodel


# ── fixtures ─────────────────────────────────────────────────────────────────

def _plain_gem_yaml(path) -> None:
    """Write a minimal GECKO-style YAML *without* any ec-rxns section."""
    path.write_text(
        """---
!!omap
- metaData: !!omap
  - id: plain_gem
  - name: Plain GEM for import tests
- metabolites:
  - !!omap
    - id: A
    - name: A
    - compartment: c
  - !!omap
    - id: B
    - name: B
    - compartment: c
- reactions:
  - !!omap
    - id: r_in
    - name: inflow
    - metabolites: !!omap
      - A: 1
    - lower_bound: 0
    - upper_bound: 1000
    - gene_reaction_rule: ''
  - !!omap
    - id: r_conv
    - name: A to B
    - metabolites: !!omap
      - A: -1
      - B: 1
    - lower_bound: 0
    - upper_bound: 1000
    - gene_reaction_rule: 'gA or gB'
  - !!omap
    - id: r_out
    - name: outflow
    - metabolites: !!omap
      - B: -1
    - lower_bound: 0
    - upper_bound: 1000
    - gene_reaction_rule: ''
- genes:
  - !!omap
    - id: gA
  - !!omap
    - id: gB
- compartments: !!omap
  - c: cytoplasm
""",
        encoding="utf-8",
    )


def _build_tiny_ecmodel():
    m = cobra.Model("ec_src")
    A = cobra.Metabolite("A", compartment="c")
    B = cobra.Metabolite("B", compartment="c")
    m.add_metabolites([A, B])
    r = cobra.Reaction("r1")
    r.add_metabolites({A: -1, B: 1})
    r.bounds = (0, 1000)
    r.gene_reaction_rule = "gA"
    ex_a = cobra.Reaction("EX_A")
    ex_a.add_metabolites({A: -1})
    ex_a.bounds = (-10, 1000)
    ex_b = cobra.Reaction("EX_B")
    ex_b.add_metabolites({B: -1})
    ex_b.bounds = (0, 1000)
    m.add_reactions([r, ex_a, ex_b])
    m.objective = "EX_B"

    ec = ECModelData()
    ec.kcat_entries = [{
        "proteins": ["P001"], "genes": ["gA"], "gene_name": "gA",
        "kcat": 100.0, "rxns": ["r1"], "notes": "", "stoicho": [1],
    }]
    ec.uniprot_data = {
        "P001": {"gene": "gA", "mw_da": 50000.0, "ec": "", "sequence": "M"},
    }
    built, _, _ = build_ecmodel(m, ec)
    return built, ec


# ── plain GEM YAML ────────────────────────────────────────────────────────────

class TestPlainGemYamlLoad:
    """gecko3-yaml-import FR-2 / FR-3 — ``ec-rxns`` absent ⇒ ``is_ecmodel=False``."""

    def test_is_ecmodel_false_when_ec_rxns_absent(self, tmp_path):
        path = tmp_path / "plain.yml"
        _plain_gem_yaml(path)

        model, ec_data = load_ecmodel(path)

        assert ec_data.is_ecmodel is False

    def test_ec_structure_is_empty(self, tmp_path):
        path = tmp_path / "plain.yml"
        _plain_gem_yaml(path)

        _, ec_data = load_ecmodel(path)

        assert ec_data.ec.n_rxns() == 0
        assert ec_data.ec.n_enzymes() == 0

    def test_cobra_model_is_populated(self, tmp_path):
        path = tmp_path / "plain.yml"
        _plain_gem_yaml(path)

        model, _ = load_ecmodel(path)

        assert sorted(m.id for m in model.metabolites) == ["A", "B"]
        assert sorted(r.id for r in model.reactions) == ["r_conv", "r_in", "r_out"]
        assert model.reactions.r_conv.gene_reaction_rule == "gA or gB"

    def test_loaded_gem_supports_subsequent_build(self, tmp_path):
        """Plain GEM YAML → Build ecModel must succeed (FR-3 key benefit)."""
        path = tmp_path / "plain.yml"
        _plain_gem_yaml(path)

        model, ec_data = load_ecmodel(path)
        ec_data.kcat_entries = [{
            "proteins": ["P001"], "genes": ["gA"], "gene_name": "gA",
            "kcat": 50.0, "rxns": ["r_conv"], "notes": "", "stoicho": [1],
        }]
        ec_data.uniprot_data = {
            "P001": {"gene": "gA", "mw_da": 50000.0, "ec": "", "sequence": "M"},
        }

        built, _, _ = build_ecmodel(model, ec_data)

        assert ec_data.is_ecmodel is True
        assert "prot_pool_exchange" in built.reactions
        assert any(r.id.startswith("usage_prot_") for r in built.reactions)


# ── ecModel YAML (regression) ────────────────────────────────────────────────

class TestEcmodelYamlRegression:
    """Round-trip through load_ecmodel still yields ``is_ecmodel=True``."""

    def test_round_trip_keeps_is_ecmodel_true(self, tmp_path):
        built, ec = _build_tiny_ecmodel()
        path = tmp_path / "ec.yml"
        save_ecmodel(built, ec, path)

        _, reloaded_ec = load_ecmodel(path)

        assert reloaded_ec.is_ecmodel is True
        assert reloaded_ec.ec.n_rxns() > 0


# ── None-safety for GUI consumers ────────────────────────────────────────────

class TestGuiNoneSafety:
    """Regression: GUI tooltips concatenate ``rxn.name`` directly, so loaders
    must never hand cobra a ``None`` name / compartment even if the YAML
    omitted the field."""

    def test_reaction_name_never_none_for_missing_field(self, tmp_path):
        path = tmp_path / "no_names.yml"
        path.write_text(
            """---
!!omap
- metaData: !!omap
  - id: nn
- metabolites:
  - !!omap
    - id: A
  - !!omap
    - id: B
- reactions:
  - !!omap
    - id: r1
    - metabolites: !!omap
      - A: -1
      - B: 1
    - lower_bound: 0
    - upper_bound: 1000
""",
            encoding="utf-8",
        )
        model, _ = load_ecmodel(path)
        r = model.reactions.get_by_id("r1")
        assert r.name == "" and r.name is not None
        # GUI pattern from reactions_list.update_tooltips
        assert ("Id: " + r.id + "\nName: " + r.name).startswith("Id: r1")

    def test_metabolite_name_and_compartment_never_none(self, tmp_path):
        path = tmp_path / "no_mnames.yml"
        path.write_text(
            """---
!!omap
- metaData: !!omap
  - id: mn
- metabolites:
  - !!omap
    - id: A
  - !!omap
    - id: B
- reactions:
  - !!omap
    - id: r1
    - metabolites: !!omap
      - A: -1
      - B: 1
    - lower_bound: 0
    - upper_bound: 1000
""",
            encoding="utf-8",
        )
        model, _ = load_ecmodel(path)
        for met in model.metabolites:
            assert met.name is not None, f"{met.id} has None name"
            assert met.compartment is not None, f"{met.id} has None compartment"


# ── malformed inputs (unchanged contract) ────────────────────────────────────

class TestMalformedInputs:
    def test_missing_metabolites_still_errors(self, tmp_path):
        path = tmp_path / "bad.yml"
        path.write_text(
            """---
!!omap
- metaData: !!omap
  - id: empty
""",
            encoding="utf-8",
        )
        with pytest.raises(EcYamlError, match="metabolites|reactions"):
            load_ecmodel(path)
