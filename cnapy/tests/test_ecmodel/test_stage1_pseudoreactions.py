"""FR-01 tests — pseudoreaction GPR clearance.

Design Ref: §8.2 scenarios 1-2.
"""

from __future__ import annotations

import cobra
import pytest

from cnapy.ecmodel.ecmodel_builder import build_ecmodel
from cnapy.ecmodel.ecmodel_data import ECModelData
from cnapy.ecmodel.expansion import (
    PSEUDOREACTION_KEYWORD,
    clear_pseudoreaction_gpr,
)


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
