"""FR-04 tests — isozyme split (critical correctness fix).

Design Ref: §8.2 scenarios 3-7, §8.3 scenario 2 (isozyme alternative).
"""

from __future__ import annotations

import cobra
import pytest

from cnapy.ecmodel.ecmodel_builder import build_ecmodel, revert_to_gem
from cnapy.ecmodel.ecmodel_data import ECModelData
from cnapy.ecmodel.exceptions import GprParseError
from cnapy.ecmodel.expansion import (
    expand_isozymes,
    parse_gpr_to_isozyme_sets,
)


# ── parse_gpr_to_isozyme_sets ─────────────────────────────────────────────────

class TestParseGpr:
    def test_empty(self):
        assert parse_gpr_to_isozyme_sets("") == []
        assert parse_gpr_to_isozyme_sets("   ") == []

    def test_single_gene(self):
        assert parse_gpr_to_isozyme_sets("g1") == [["g1"]]

    def test_two_isozymes(self):
        assert parse_gpr_to_isozyme_sets("g1 or g2") == [["g1"], ["g2"]]

    def test_single_complex(self):
        assert parse_gpr_to_isozyme_sets("g1 and g2") == [["g1", "g2"]]

    def test_mixed_complex_and_isozyme(self):
        result = parse_gpr_to_isozyme_sets("(g1 and g2) or g3 or (g4 and g5 and g6)")
        assert result == [["g1", "g2"], ["g3"], ["g4", "g5", "g6"]]

    def test_extra_whitespace_tolerated(self):
        assert parse_gpr_to_isozyme_sets("  g1   or    g2  ") == [["g1"], ["g2"]]

    def test_case_insensitive_operators(self):
        # cobra always lowercases but third-party models sometimes use OR / AND.
        assert parse_gpr_to_isozyme_sets("g1 OR g2") == [["g1"], ["g2"]]
        assert parse_gpr_to_isozyme_sets("g1 AND g2") == [["g1", "g2"]]

    def test_nested_parentheses(self):
        # Parens around a single gene should be stripped.
        assert parse_gpr_to_isozyme_sets("(g1)") == [["g1"]]
        assert parse_gpr_to_isozyme_sets("((g1 and g2))") == [["g1", "g2"]]

    def test_unbalanced_parens_raises(self):
        with pytest.raises(GprParseError):
            parse_gpr_to_isozyme_sets("(g1 and g2")
        with pytest.raises(GprParseError):
            parse_gpr_to_isozyme_sets("g1) or g2")

    def test_gene_names_with_underscores_and_digits(self):
        assert parse_gpr_to_isozyme_sets(
            "YFL001W or (Q0045 and P12345)"
        ) == [["YFL001W"], ["Q0045", "P12345"]]


# ── expand_isozymes ───────────────────────────────────────────────────────────

class TestExpandIsozymes:
    def test_simple_or_splits_into_two(self, toy_gem: cobra.Model):
        ec_data = ECModelData()
        split_map = expand_isozymes(toy_gem, ec_data)

        # r_iso had "g1 or g2" → should split
        assert "r_iso" in split_map
        assert set(split_map["r_iso"]) == {"r_iso_EXP_1", "r_iso_EXP_2"}

        # Original r_iso is gone, replaced by two variants.
        assert "r_iso" not in [r.id for r in toy_gem.reactions]
        exp1 = toy_gem.reactions.r_iso_EXP_1
        exp2 = toy_gem.reactions.r_iso_EXP_2
        assert exp1.gene_reaction_rule == "g1"
        assert exp2.gene_reaction_rule == "g2"

    def test_mixed_complex_and_isozyme(self, toy_gem: cobra.Model):
        """(g3 and g4) or g5 → variant 1 has the complex, variant 2 has g5."""
        ec_data = ECModelData()
        expand_isozymes(toy_gem, ec_data)

        rules = [toy_gem.reactions.r_cplx_EXP_1.gene_reaction_rule,
                 toy_gem.reactions.r_cplx_EXP_2.gene_reaction_rule]
        assert "g3 and g4" in rules
        assert "g5" in rules

    def test_single_gene_not_split(self, toy_gem: cobra.Model):
        """r_simple has a single gene g6 — must stay as r_simple."""
        ec_data = ECModelData()
        expand_isozymes(toy_gem, ec_data)

        assert "r_simple" in [r.id for r in toy_gem.reactions]
        assert toy_gem.reactions.r_simple.gene_reaction_rule == "g6"

    def test_no_split_for_reactions_without_gpr(self, toy_gem: cobra.Model):
        """Reactions with empty GPR must be left alone."""
        ec_data = ECModelData()
        expand_isozymes(toy_gem, ec_data)

        assert "r_in" in [r.id for r in toy_gem.reactions]
        assert "r_out" in [r.id for r in toy_gem.reactions]

    def test_split_preserves_stoichiometry_and_bounds(self, toy_gem: cobra.Model):
        """Each _EXP_N variant has identical bounds and metabolite stoich."""
        ec_data = ECModelData()
        orig_lb = toy_gem.reactions.r_iso.lower_bound
        orig_ub = toy_gem.reactions.r_iso.upper_bound
        orig_stoich = dict(toy_gem.reactions.r_iso.metabolites)

        expand_isozymes(toy_gem, ec_data)

        for vid in ["r_iso_EXP_1", "r_iso_EXP_2"]:
            v = toy_gem.reactions.get_by_id(vid)
            assert v.lower_bound == orig_lb
            assert v.upper_bound == orig_ub
            assert {m.id: c for m, c in v.metabolites.items()} == \
                   {m.id: c for m, c in orig_stoich.items()}

    def test_records_split_map_in_ec_data(self, toy_gem: cobra.Model):
        ec_data = ECModelData()
        expand_isozymes(toy_gem, ec_data)

        assert ec_data.isozyme_split_map["r_iso"] == ["r_iso_EXP_1", "r_iso_EXP_2"]

    def test_gecko_light_is_noop(self, toy_gem: cobra.Model):
        ec_data = ECModelData()
        ec_data.gecko_light = True

        split_map = expand_isozymes(toy_gem, ec_data)

        assert split_map == {}
        assert "r_iso" in [r.id for r in toy_gem.reactions]  # not split


# ── integration via build_ecmodel ─────────────────────────────────────────────

class TestBuildWithIsozymes:
    def test_build_ecmodel_splits_isozymes(self, toy_gem: cobra.Model):
        ec_data = ECModelData()
        ecmodel, _, _ = build_ecmodel(toy_gem, ec_data)

        names = {r.id for r in ecmodel.reactions}
        assert "r_iso_EXP_1" in names
        assert "r_iso_EXP_2" in names
        assert "r_iso" not in names

        # Recorded in ec_data for revert.
        assert "r_iso" in ec_data.isozyme_split_map

    def test_revert_merges_isozymes_back(self, toy_gem: cobra.Model):
        """build → revert must give a model equivalent to the original."""
        ec_data = ECModelData()

        ecmodel, _, _ = build_ecmodel(toy_gem, ec_data)
        # Sanity: isozymes were split.
        assert "r_iso_EXP_1" in [r.id for r in ecmodel.reactions]

        gem = revert_to_gem(ecmodel, ec_data)

        # Isozyme variants are gone, original r_iso is back.
        ids = {r.id for r in gem.reactions}
        assert "r_iso" in ids
        assert "r_iso_EXP_1" not in ids
        assert "r_iso_EXP_2" not in ids

        # Restored GPR expresses both isozymes via OR.
        gpr = gem.reactions.r_iso.gene_reaction_rule
        assert "g1" in gpr and "g2" in gpr
        assert " or " in gpr

    def test_isozyme_alternative_is_or_not_and(self, toy_gem: cobra.Model):
        """The critical bug this FR fixes: isozymes must be alternatives (OR).

        If the two isozyme enzymes were required simultaneously (AND) at the
        same reaction node (as the pre-FR-04 code did when customKcats had
        separate entries per isozyme), deleting one variant should block all
        flux. After the fix, the other variant can still carry flux.
        """
        ec_data = ECModelData()
        ecmodel, _, _ = build_ecmodel(toy_gem, ec_data)

        # Baseline: with both variants of the isozymic reaction, FBA is
        # feasible (objective > 0).
        sol_full = ecmodel.optimize()
        assert sol_full.status == "optimal"
        assert sol_full.objective_value > 0

        # Disable EXP_1 — objective should still be reachable via EXP_2.
        with ecmodel as m:
            m.reactions.r_iso_EXP_1.upper_bound = 0
            sol_one = m.optimize()
        assert sol_one.status == "optimal", "EXP_2 alone must sustain flux"
        assert sol_one.objective_value > 0

        # Disable EXP_2 — still feasible via EXP_1.
        with ecmodel as m:
            m.reactions.r_iso_EXP_2.upper_bound = 0
            sol_two = m.optimize()
        assert sol_two.status == "optimal", "EXP_1 alone must sustain flux"
        assert sol_two.objective_value > 0


class TestReversibleIsozymeInteraction:
    """FR-03 + FR-04 composition (reversible split + isozyme split)."""

    def test_reversible_isozyme_produces_four_variants(
        self, toy_reversible_gem: cobra.Model
    ):
        ec_data = ECModelData()
        ecmodel, _, _ = build_ecmodel(toy_reversible_gem, ec_data)

        ids = {r.id for r in ecmodel.reactions}
        # Expected variants: forward EXP_1/EXP_2 and reverse EXP_1/EXP_2.
        assert "r_rev_EXP_1" in ids
        assert "r_rev_EXP_2" in ids
        assert "r_rev_REV_EXP_1" in ids
        assert "r_rev_REV_EXP_2" in ids

    def test_revert_round_trip_restores_reversible(
        self, toy_reversible_gem: cobra.Model
    ):
        ec_data = ECModelData()
        ecmodel, _, _ = build_ecmodel(toy_reversible_gem, ec_data)
        gem = revert_to_gem(ecmodel, ec_data)

        ids = {r.id for r in gem.reactions}
        assert "r_rev" in ids
        assert "r_rev_REV" not in ids
        assert "r_rev_EXP_1" not in ids

        restored = gem.reactions.r_rev
        assert restored.lower_bound < 0 and restored.upper_bound > 0
        # GPR has both isozymes via OR.
        gpr = restored.gene_reaction_rule
        assert "gA" in gpr and "gB" in gpr
        assert " or " in gpr
