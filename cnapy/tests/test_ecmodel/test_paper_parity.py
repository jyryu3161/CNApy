"""gecko3-paper-parity tests — Light mode + FR-02 invert backward-only.

Design Ref: gecko3-paper-parity.plan.md FR-A / FR-B / FR-C. Brings the
GECKO 3.0 implementation to 100% parity with `makeEcModel.m` Step 2,
Step 4 (light), and the light-mode complex cost calculation.
"""

from __future__ import annotations

import cobra
import pytest

from cnapy.ecmodel.ecmodel_builder import build_ecmodel
from cnapy.ecmodel.ecmodel_data import ECModelData
from cnapy.ecmodel.expansion import invert_backward_only


# ── FR-C: invert_backward_only ────────────────────────────────────────────────

class TestInvertBackwardOnly:
    """GECKO 3.0 Box 1 Step 2 — flip backward-only reactions to forward."""

    def _backward_only_model(self) -> cobra.Model:
        m = cobra.Model("bo")
        a = cobra.Metabolite("A", compartment="c")
        b = cobra.Metabolite("B", compartment="c")
        m.add_metabolites([a, b])
        # backward-only: lb<0, ub=0
        r_back = cobra.Reaction("r_back")
        r_back.add_metabolites({a: -1, b: 1})
        r_back.bounds = (-100, 0)
        # forward irreversible
        r_fwd = cobra.Reaction("r_fwd")
        r_fwd.add_metabolites({a: -1, b: 1})
        r_fwd.bounds = (0, 1000)
        # reversible
        r_rev = cobra.Reaction("r_rev")
        r_rev.add_metabolites({a: -1, b: 1})
        r_rev.bounds = (-500, 1000)
        m.add_reactions([r_back, r_fwd, r_rev])
        return m

    def test_backward_only_is_flipped(self):
        m = self._backward_only_model()
        inverted = invert_backward_only(m)

        assert inverted == ["r_back"]
        r = m.reactions.get_by_id("r_back")
        assert r.bounds == (0.0, 100.0)
        # Stoich was {A:-1, B:+1}; after flip {A:+1, B:-1}
        a = m.metabolites.get_by_id("A")
        b = m.metabolites.get_by_id("B")
        assert r.metabolites[a] == 1.0
        assert r.metabolites[b] == -1.0

    def test_forward_and_reversible_untouched(self):
        m = self._backward_only_model()
        invert_backward_only(m)

        # forward irreversible: unchanged
        r_fwd = m.reactions.get_by_id("r_fwd")
        assert r_fwd.bounds == (0.0, 1000.0)
        a, b = m.metabolites.get_by_id("A"), m.metabolites.get_by_id("B")
        assert r_fwd.metabolites[a] == -1.0
        assert r_fwd.metabolites[b] == 1.0

        # reversible: unchanged
        r_rev = m.reactions.get_by_id("r_rev")
        assert r_rev.bounds == (-500.0, 1000.0)

    def test_idempotent(self):
        """After flipping, a second call must be a no-op (now forward)."""
        m = self._backward_only_model()
        first = invert_backward_only(m)
        second = invert_backward_only(m)

        assert first == ["r_back"]
        assert second == []   # already forward; nothing to flip again
        r = m.reactions.get_by_id("r_back")
        assert r.bounds == (0.0, 100.0)

    def test_returns_changed_ids_only(self):
        m = self._backward_only_model()
        inverted = invert_backward_only(m)
        all_ids = {r.id for r in m.reactions}
        assert set(inverted) <= all_ids
        assert set(inverted) == {"r_back"}

    def test_no_backward_reactions_returns_empty(self):
        m = cobra.Model("clean")
        a = cobra.Metabolite("A", compartment="c")
        b = cobra.Metabolite("B", compartment="c")
        m.add_metabolites([a, b])
        r = cobra.Reaction("r1")
        r.add_metabolites({a: -1, b: 1})
        r.bounds = (0, 1000)
        m.add_reactions([r])

        assert invert_backward_only(m) == []


# ── FR-C2: Step 2 ordering inside build_ecmodel ──────────────────────────────

class TestBuildOrdering:
    """invert must run BEFORE reversible split so backward-only is normalised."""

    def test_backward_only_normalised_before_split(self):
        m = cobra.Model("ord")
        a = cobra.Metabolite("A", compartment="c")
        b = cobra.Metabolite("B", compartment="c")
        m.add_metabolites([a, b])
        r = cobra.Reaction("r_back")
        r.add_metabolites({a: -1, b: 1})
        r.bounds = (-100, 0)   # backward-only
        r.gene_reaction_rule = "gA"
        ex_a = cobra.Reaction("EX_A")
        ex_a.add_metabolites({a: -1})
        ex_a.bounds = (-10, 1000)
        ex_b = cobra.Reaction("EX_B")
        ex_b.add_metabolites({b: -1})
        ex_b.bounds = (0, 1000)
        m.add_reactions([r, ex_a, ex_b])

        ec = ECModelData()
        ec.kcat_entries = []
        ec.uniprot_data = {}
        built, _, _ = build_ecmodel(m, ec)

        # After invert + split: r_back should be forward, no _REV variant
        # (because the original was lb<0,ub=0, which becomes lb=0,ub=100
        # after invert — irreversible, so no split needed).
        r_after = built.reactions.get_by_id("r_back")
        assert r_after.bounds == (0.0, 100.0)
        assert "r_back_REV" not in built.reactions


# ── FR-A: Light mode reversible split ────────────────────────────────────────

class TestLightReversibleSplit:
    """gecko3-paper-parity FR-A — Light mode must split reversible reactions
    so reverse flux is enzyme-costed correctly (matches makeEcModel.m Step 4)."""

    def test_light_creates_rev_variants(self):
        m = cobra.Model("light_rev")
        a = cobra.Metabolite("A", compartment="c")
        b = cobra.Metabolite("B", compartment="c")
        m.add_metabolites([a, b])
        r = cobra.Reaction("r1")
        r.add_metabolites({a: -1, b: 1})
        r.bounds = (-1000, 1000)   # reversible
        r.gene_reaction_rule = "gA"
        ex_a = cobra.Reaction("EX_A")
        ex_a.add_metabolites({a: -1})
        ex_a.bounds = (-10, 1000)
        ex_b = cobra.Reaction("EX_B")
        ex_b.add_metabolites({b: -1})
        ex_b.bounds = (0, 1000)
        m.add_reactions([r, ex_a, ex_b])

        ec = ECModelData()
        ec.gecko_light = True
        ec.kcat_entries = []
        ec.uniprot_data = {}
        built, _, _ = build_ecmodel(m, ec)

        # r1 should be split into forward (r1, lb=0) + reverse (r1_REV)
        assert "r1" in built.reactions
        assert "r1_REV" in built.reactions
        assert built.reactions.get_by_id("r1").lower_bound == 0.0


# ── FR-B: Light complex cost SUM ──────────────────────────────────────────────

class TestLightComplexSum:
    """gecko3-paper-parity FR-B — within one entry, complex protein costs sum."""

    def _toy_model(self) -> cobra.Model:
        m = cobra.Model("light_complex")
        a = cobra.Metabolite("A", compartment="c")
        b = cobra.Metabolite("B", compartment="c")
        m.add_metabolites([a, b])
        r = cobra.Reaction("r1")
        r.add_metabolites({a: -1, b: 1})
        r.bounds = (0, 1000)
        r.gene_reaction_rule = "gA and gB"
        ex_a = cobra.Reaction("EX_A")
        ex_a.add_metabolites({a: -1})
        ex_a.bounds = (-10, 1000)
        ex_b = cobra.Reaction("EX_B")
        ex_b.add_metabolites({b: -1})
        ex_b.bounds = (0, 1000)
        m.add_reactions([r, ex_a, ex_b])
        return m

    def test_complex_cost_is_sum_not_min(self):
        """Complex P1+P2 should cost (MW1 + MW2)/(kcat*3600), not min()."""
        m = self._toy_model()
        ec = ECModelData()
        ec.gecko_light = True
        ec.kcat_entries = [{
            "proteins": ["P1", "P2"],
            "genes": ["gA", "gB"],
            "gene_name": "complex",
            "kcat": 100.0,
            "rxns": ["r1"],
            "notes": "",
            "stoicho": [1, 1],
        }]
        ec.uniprot_data = {
            "P1": {"gene": "gA", "mw_da": 50000.0, "ec": "", "sequence": "M"},
            "P2": {"gene": "gB", "mw_da": 30000.0, "ec": "", "sequence": "M"},
        }
        built, _, _ = build_ecmodel(m, ec)

        # Expected coefficient = -(50000 + 30000) / (100 * 3600) ≈ -0.2222
        expected_sum = -(50000 + 30000) / (100 * 3600)
        # Old buggy MIN would give the smallest term: -30000/(100*3600) ≈ -0.0833
        wrong_min = -30000 / (100 * 3600)

        pool_met = built.metabolites.get_by_id("prot_pool")
        # r1 might be split into reversible variants — pick the forward one.
        candidates = [r for r in built.reactions if r.id.startswith("r1")
                      and not r.id.startswith("usage_")]
        assert candidates, "no r1 variant survived build"

        coeff = candidates[0].metabolites.get(pool_met, 0.0)
        assert coeff == pytest.approx(expected_sum, rel=1e-6)
        assert coeff != pytest.approx(wrong_min, rel=1e-6), \
            "regression: still using MIN over complex members"

    def test_dlkcat_single_protein_unchanged(self):
        """For a single-protein entry, SUM == MIN == that value (no regression)."""
        m = self._toy_model()
        m.reactions.get_by_id("r1").gene_reaction_rule = "gA"

        ec = ECModelData()
        ec.gecko_light = True
        ec.kcat_entries = [{
            "proteins": ["P1"],
            "genes": ["gA"],
            "gene_name": "gA",
            "kcat": 100.0,
            "rxns": ["r1"],
            "notes": "",
            "stoicho": [1],
        }]
        ec.uniprot_data = {
            "P1": {"gene": "gA", "mw_da": 50000.0, "ec": "", "sequence": "M"},
        }
        built, _, _ = build_ecmodel(m, ec)

        expected = -50000 / (100 * 3600)
        pool_met = built.metabolites.get_by_id("prot_pool")
        candidates = [r for r in built.reactions if r.id.startswith("r1")
                      and not r.id.startswith("usage_")]
        coeff = candidates[0].metabolites.get(pool_met, 0.0)
        assert coeff == pytest.approx(expected, rel=1e-6)
