"""FR-01 tests — pseudoreaction GPR clearance.

Design Ref: §8.2 scenarios 1-2.
"""

from __future__ import annotations

import math

import cobra
import pytest

from cnapy.ecmodel.ecmodel_builder import (
    _stoich_coeff,
    build_ecmodel,
    parse_customkcats_file,
    parse_uniprot_file,
)
from cnapy.ecmodel.ecmodel_data import ECModelData
from cnapy.ecmodel.expansion import (
    PSEUDOREACTION_KEYWORD,
    clear_pseudoreaction_gpr,
)


# ── NaN-safety for DLKcat / UniProt inputs ────────────────────────────────────

class TestNaNGuards:
    """Regression: DLKcat.tsv contains rows whose predicted kcat is NaN
    (substrate out of training distribution). These must be skipped by the
    parser and neutralised by ``_stoich_coeff`` — otherwise they leak into
    the solver as NaN stoichiometric coefficients."""

    def test_stoich_coeff_rejects_nan_kcat(self):
        assert _stoich_coeff(10000, float("nan")) == 0.0

    def test_stoich_coeff_rejects_nan_mw(self):
        assert _stoich_coeff(float("nan"), 10) == 0.0

    def test_stoich_coeff_rejects_inf(self):
        assert _stoich_coeff(float("inf"), 10) == 0.0
        assert _stoich_coeff(10000, float("inf")) == 0.0

    def test_parser_skips_nan_kcat_rows(self, tmp_path):
        tsv = tmp_path / "kcat.tsv"
        tsv.write_text(
            "proteins\tkcat\trxns\tstoicho\n"
            "P001\t100.0\tr1\t1\n"
            "P002\t\tr2\t1\n"          # empty → pandas NaN → must be skipped
            "P003\tnan\tr3\t1\n"         # literal NaN → must be skipped
            "P004\t50.0\tr4\t1\n",
        )
        entries = parse_customkcats_file(str(tsv))
        kcats = [e["kcat"] for e in entries]
        assert all(math.isfinite(k) and k > 0 for k in kcats)
        rxns = {e["rxns"][0] for e in entries}
        assert rxns == {"r1", "r4"}

    def test_parser_zeroes_nan_mw(self, tmp_path):
        tsv = tmp_path / "uni.tsv"
        tsv.write_text(
            "Entry\tMass\tGene Names\n"
            "P001\t50000\tgA\n"
            "P002\t\tgB\n"               # empty mass → pandas NaN
            "P003\tnot_a_number\tgC\n",   # unparseable
        )
        uni = parse_uniprot_file(str(tsv))
        assert all(math.isfinite(v["mw_da"]) for v in uni.values())
        assert uni["P002"]["mw_da"] == 0.0
        assert uni["P003"]["mw_da"] == 0.0


# ── existing FR-01 cases ──────────────────────────────────────────────────────


class TestClearPseudoreactionGpr:
    """Unit tests for clear_pseudoreaction_gpr()."""

    def test_clears_name_matched_pseudoreaction(self, toy_gem: cobra.Model):
        """BIOMASS (name contains 'pseudoreaction') → GPR cleared."""
        assert toy_gem.reactions.BIOMASS.gene_reaction_rule == "g_bio"

        cleared = clear_pseudoreaction_gpr(toy_gem)

        assert "BIOMASS" in cleared
        assert toy_gem.reactions.BIOMASS.gene_reaction_rule == ""

    def test_preserves_non_pseudoreaction_gprs(self, toy_gem: cobra.Model):
        """Reactions with GPR but not matching keyword stay untouched."""
        clear_pseudoreaction_gpr(toy_gem)

        assert toy_gem.reactions.r_iso.gene_reaction_rule == "g1 or g2"
        assert toy_gem.reactions.r_simple.gene_reaction_rule == "g6"

    def test_returns_empty_when_no_pseudoreactions(self, simple_model: cobra.Model):
        """simple_model (from shared conftest) has no pseudoreaction name."""
        cleared = clear_pseudoreaction_gpr(simple_model)
        assert cleared == []

    def test_extra_tsv_flags_by_id(self, toy_gem: cobra.Model, tmp_path):
        """Reactions listed in extra_tsv get GPR cleared even without name match."""
        # r_simple has gene_reaction_rule="g6" and no "pseudoreaction" in its name.
        tsv = tmp_path / "pseudoRxns.tsv"
        tsv.write_text("# comment line\nr_simple\nNON_EXISTENT\n")

        cleared = clear_pseudoreaction_gpr(toy_gem, extra_tsv=str(tsv))

        assert "r_simple" in cleared
        assert toy_gem.reactions.r_simple.gene_reaction_rule == ""

    def test_extra_tsv_missing_file_is_ok(self, toy_gem: cobra.Model, tmp_path):
        """Non-existent TSV should silently fall back to name-only detection."""
        missing = tmp_path / "does_not_exist.tsv"
        cleared = clear_pseudoreaction_gpr(toy_gem, extra_tsv=str(missing))
        # Name-based detection still works.
        assert "BIOMASS" in cleared

    def test_case_insensitive_keyword_match(self):
        """The keyword match must be case-insensitive."""
        model = cobra.Model("t")
        a = cobra.Metabolite("a", compartment="c")
        b = cobra.Metabolite("b", compartment="c")
        model.add_metabolites([a, b])

        r = cobra.Reaction("r1")
        r.name = "Biomass PSEUDOREACTION (all caps)"   # uppercase keyword
        r.add_metabolites({a: -1, b: 1})
        r.gene_reaction_rule = "gX"
        model.add_reactions([r])

        cleared = clear_pseudoreaction_gpr(model)
        assert "r1" in cleared
        assert r.gene_reaction_rule == ""

    def test_pseudoreaction_keyword_constant(self):
        """Sanity check on the keyword constant we're matching against."""
        assert PSEUDOREACTION_KEYWORD == "pseudoreaction"


class TestPseudoreactionClearViaBuilder:
    """Integration — FR-01 runs automatically inside build_ecmodel."""

    def test_build_ecmodel_clears_biomass_gpr(self, toy_gem: cobra.Model):
        """After build_ecmodel, BIOMASS carries no enzyme constraint."""
        ec_data = ECModelData()
        ec_data.kcat_entries = []     # no kcats — only structural checks
        ec_data.uniprot_data = {}

        ecmodel, _, _ = build_ecmodel(toy_gem, ec_data)

        # BIOMASS still exists (pseudoreaction is not deleted, only GPR cleared).
        biomass = ecmodel.reactions.BIOMASS
        assert biomass.gene_reaction_rule == ""
        assert len(biomass.genes) == 0

        # No prot_* metabolite should have been added to BIOMASS.
        prot_coeffs = [c for m, c in biomass.metabolites.items()
                       if m.id.startswith("prot_")]
        assert prot_coeffs == [], "biomass must not draw from protein pool"

    def test_original_source_model_is_not_mutated(self, toy_gem: cobra.Model):
        """FR-01 modifications happen on the internal copy, not the source."""
        original_rule = toy_gem.reactions.BIOMASS.gene_reaction_rule
        ec_data = ECModelData()

        build_ecmodel(toy_gem, ec_data)

        assert toy_gem.reactions.BIOMASS.gene_reaction_rule == original_rule
