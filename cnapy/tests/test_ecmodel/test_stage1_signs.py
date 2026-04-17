"""FR-07 + FR-08 tests — GECKO 3.0 sign convention for pool & usage reactions.

Design Ref: §8.2 scenarios 11 (unit), §8.3 scenario 4 (integration).

The v2 (GECKO 3.0) convention:
  * ``prot_pool_exchange`` : stoich ``{prot_pool: -1}``,
                             bounds ``(-pool_bound, 0)``.
  * ``usage_prot_{id}``    : stoich ``{prot_i: -1, prot_pool: +1}``,
                             bounds ``(-1000, 0)``.

The v1 (legacy CNApy) convention is the opposite sign / bound polarity.
``flip_reaction_signs_v1_to_v2`` is the migration helper used by
``ECModelData.upgrade`` (see ``test_cna_migration.py``).
"""

from __future__ import annotations

import cobra

from cnapy.ecmodel.ecmodel_builder import build_ecmodel
from cnapy.ecmodel.ecmodel_data import ECModelData
from cnapy.ecmodel.expansion import (
    add_enzyme_metabolite_and_usage,
    add_protein_pool,
    flip_reaction_signs_v1_to_v2,
    set_usage_capacity,
    usage_capacity,
)


# ── helpers ────────────────────────────────────────────────────────────────────

def _bare_model() -> cobra.Model:
    """Minimal cobra model with a single A→B reaction — no prot_pool yet."""
    m = cobra.Model("bare")
    a = cobra.Metabolite("A", compartment="c")
    b = cobra.Metabolite("B", compartment="c")
    m.add_metabolites([a, b])
    r = cobra.Reaction("r1")
    r.add_metabolites({a: -1, b: 1})
    r.bounds = (0, 1000)
    m.add_reactions([r])
    return m


def _make_v1_pool_and_usage(model: cobra.Model,
                             uniprot_id: str = "P00001",
                             pool_ub: float = 42.0,
                             usage_ub: float = 1000.0) -> None:
    """Manually build v1 pool + usage reactions on ``model`` (for migration tests).

    v1 convention (legacy CNApy):
      * ``prot_pool_exchange`` stoich ``{pool: +1}``, bounds ``(0, pool_ub)``.
      * ``usage_prot_{id}``    stoich ``{pool: -1, prot: +1}``, bounds ``(0, usage_ub)``.
    """
    pool_met = cobra.Metabolite("prot_pool", name="Protein pool", compartment="c")
    enz_met = cobra.Metabolite(
        f"prot_{uniprot_id}", name=f"Enzyme {uniprot_id}", compartment="c",
    )
    model.add_metabolites([pool_met, enz_met])

    pool_rxn = cobra.Reaction(
        "prot_pool_exchange", name="Protein pool exchange",
        lower_bound=0.0, upper_bound=pool_ub,
    )
    pool_rxn.add_metabolites({pool_met: 1.0})

    usage_rxn = cobra.Reaction(
        f"usage_prot_{uniprot_id}", name=f"Enzyme usage {uniprot_id}",
        lower_bound=0.0, upper_bound=usage_ub,
    )
    usage_rxn.add_metabolites({pool_met: -1.0, enz_met: 1.0})

    model.add_reactions([pool_rxn, usage_rxn])


# ── add_protein_pool ──────────────────────────────────────────────────────────

class TestAddProteinPool:
    """FR-08: protein pool exchange uses negative-flux (v2) convention."""

    def test_adds_pool_metabolite(self):
        m = _bare_model()
        add_protein_pool(m, pool_bound_mg_gDCW=100.0)
        assert "prot_pool" in m.metabolites

    def test_adds_exchange_reaction(self):
        m = _bare_model()
        add_protein_pool(m, pool_bound_mg_gDCW=100.0)
        assert "prot_pool_exchange" in m.reactions

    def test_exchange_stoichiometry_is_negative_one(self):
        m = _bare_model()
        add_protein_pool(m, pool_bound_mg_gDCW=100.0)
        rxn = m.reactions.get_by_id("prot_pool_exchange")
        pool_met = m.metabolites.get_by_id("prot_pool")
        assert rxn.metabolites == {pool_met: -1.0}

    def test_exchange_bounds_are_negative_flux(self):
        m = _bare_model()
        add_protein_pool(m, pool_bound_mg_gDCW=100.0)
        rxn = m.reactions.get_by_id("prot_pool_exchange")
        assert rxn.lower_bound == -100.0
        assert rxn.upper_bound == 0.0

    def test_negative_pool_bound_is_taken_as_magnitude(self):
        """Callers should be able to pass pool_bound either sign."""
        m = _bare_model()
        add_protein_pool(m, pool_bound_mg_gDCW=-100.0)
        rxn = m.reactions.get_by_id("prot_pool_exchange")
        assert rxn.lower_bound == -100.0
        assert rxn.upper_bound == 0.0

    def test_idempotent_and_refreshes_bound(self):
        m = _bare_model()
        add_protein_pool(m, pool_bound_mg_gDCW=50.0)
        add_protein_pool(m, pool_bound_mg_gDCW=200.0)

        # Metabolite + reaction exist exactly once.
        assert sum(1 for met in m.metabolites if met.id == "prot_pool") == 1
        assert sum(1 for rxn in m.reactions if rxn.id == "prot_pool_exchange") == 1

        # Bound is refreshed to the most recent call.
        rxn = m.reactions.get_by_id("prot_pool_exchange")
        assert rxn.lower_bound == -200.0
        assert rxn.upper_bound == 0.0


# ── add_enzyme_metabolite_and_usage ───────────────────────────────────────────

class TestAddEnzymeMetaboliteAndUsage:
    """FR-07: usage reaction uses v2 convention."""

    def test_returns_expected_ids(self):
        m = _bare_model()
        add_protein_pool(m)
        met_id, rxn_id = add_enzyme_metabolite_and_usage(m, "P12345")
        assert met_id == "prot_P12345"
        assert rxn_id == "usage_prot_P12345"

    def test_usage_reaction_stoichiometry_is_v2(self):
        m = _bare_model()
        add_protein_pool(m)
        add_enzyme_metabolite_and_usage(m, "P12345")

        rxn = m.reactions.get_by_id("usage_prot_P12345")
        pool_met = m.metabolites.get_by_id("prot_pool")
        prot_met = m.metabolites.get_by_id("prot_P12345")

        # {prot_i: -1, prot_pool: +1} — usage consumes enzyme, releases pool.
        assert rxn.metabolites[prot_met] == -1.0
        assert rxn.metabolites[pool_met] == 1.0

    def test_usage_reaction_bounds_are_negative_flux(self):
        m = _bare_model()
        add_protein_pool(m)
        add_enzyme_metabolite_and_usage(m, "P12345")
        rxn = m.reactions.get_by_id("usage_prot_P12345")
        assert rxn.lower_bound == -1000.0
        assert rxn.upper_bound == 0.0

    def test_idempotent(self):
        m = _bare_model()
        add_protein_pool(m)
        add_enzyme_metabolite_and_usage(m, "P12345")
        add_enzyme_metabolite_and_usage(m, "P12345")

        assert sum(1 for met in m.metabolites if met.id == "prot_P12345") == 1
        assert sum(1 for rxn in m.reactions if rxn.id == "usage_prot_P12345") == 1


# ── usage_capacity / set_usage_capacity helpers ───────────────────────────────

class TestUsageCapacityHelpers:
    """Sign-agnostic readout + writer for legacy/modern callers."""

    def test_usage_capacity_reads_v2(self):
        """v2 convention: capacity is |lower_bound|."""
        m = _bare_model()
        add_protein_pool(m)
        add_enzyme_metabolite_and_usage(m, "P1")
        rxn = m.reactions.get_by_id("usage_prot_P1")
        rxn.lower_bound = -250.0
        rxn.upper_bound = 0.0
        assert usage_capacity(rxn) == 250.0

    def test_usage_capacity_reads_v1(self):
        """v1 convention (positive flux): capacity is upper_bound."""
        m = _bare_model()
        _make_v1_pool_and_usage(m, "P1", usage_ub=500.0)
        rxn = m.reactions.get_by_id("usage_prot_P1")
        assert usage_capacity(rxn) == 500.0

    def test_set_usage_capacity_writes_v2_bounds(self):
        m = _bare_model()
        add_protein_pool(m)
        add_enzyme_metabolite_and_usage(m, "P1")
        rxn = m.reactions.get_by_id("usage_prot_P1")

        set_usage_capacity(rxn, 123.0)
        assert rxn.lower_bound == -123.0
        assert rxn.upper_bound == 0.0

    def test_set_usage_capacity_ignores_sign(self):
        m = _bare_model()
        add_protein_pool(m)
        add_enzyme_metabolite_and_usage(m, "P1")
        rxn = m.reactions.get_by_id("usage_prot_P1")

        set_usage_capacity(rxn, -77.0)
        assert rxn.lower_bound == -77.0
        assert rxn.upper_bound == 0.0


# ── flip_reaction_signs_v1_to_v2 ──────────────────────────────────────────────

class TestFlipReactionSignsV1ToV2:
    """FR-20 migration helper — flips sign convention in place."""

    def test_flip_pool_and_usage_stoichiometry(self):
        m = _bare_model()
        _make_v1_pool_and_usage(m, "P1", pool_ub=42.0, usage_ub=1000.0)

        flip_reaction_signs_v1_to_v2(m)

        pool_rxn = m.reactions.get_by_id("prot_pool_exchange")
        usage_rxn = m.reactions.get_by_id("usage_prot_P1")
        pool_met = m.metabolites.get_by_id("prot_pool")
        prot_met = m.metabolites.get_by_id("prot_P1")

        # v2 stoich: pool exchange {pool: -1}, usage {prot: -1, pool: +1}
        assert pool_rxn.metabolites == {pool_met: -1.0}
        assert usage_rxn.metabolites[prot_met] == -1.0
        assert usage_rxn.metabolites[pool_met] == 1.0

    def test_flip_pool_and_usage_bounds(self):
        m = _bare_model()
        _make_v1_pool_and_usage(m, "P1", pool_ub=42.0, usage_ub=1000.0)

        flip_reaction_signs_v1_to_v2(m)

        pool_rxn = m.reactions.get_by_id("prot_pool_exchange")
        usage_rxn = m.reactions.get_by_id("usage_prot_P1")

        # (0, pool_ub) → (-pool_ub, 0) ; (0, 1000) → (-1000, 0)
        assert pool_rxn.lower_bound == -42.0
        assert pool_rxn.upper_bound == 0.0
        assert usage_rxn.lower_bound == -1000.0
        assert usage_rxn.upper_bound == 0.0

    def test_flip_returns_list_of_flipped_ids(self):
        m = _bare_model()
        _make_v1_pool_and_usage(m, "P1")
        _make_v1_pool_and_usage.__wrapped__ = None  # noqa (just marker)

        flipped = flip_reaction_signs_v1_to_v2(m)

        assert set(flipped) == {"prot_pool_exchange", "usage_prot_P1"}

    def test_flip_preserves_untargeted_reactions(self):
        """Non-usage / non-pool reactions must be left alone."""
        m = _bare_model()
        _make_v1_pool_and_usage(m, "P1")
        original_r1 = {
            met.id: coeff for met, coeff in m.reactions.get_by_id("r1").metabolites.items()
        }
        original_r1_bounds = m.reactions.get_by_id("r1").bounds

        flip_reaction_signs_v1_to_v2(m)

        r1 = m.reactions.get_by_id("r1")
        assert {met.id: coeff for met, coeff in r1.metabolites.items()} == original_r1
        assert r1.bounds == original_r1_bounds

    def test_flip_on_proteomics_applied_bounds(self):
        """v1 bounds (0, level) from applied proteomics → (-level, 0)."""
        m = _bare_model()
        # Simulate proteomics: usage upper bound set to a specific measured level.
        _make_v1_pool_and_usage(m, "P1", usage_ub=17.5)

        flip_reaction_signs_v1_to_v2(m)

        usage_rxn = m.reactions.get_by_id("usage_prot_P1")
        assert usage_rxn.lower_bound == -17.5
        assert usage_rxn.upper_bound == 0.0

    def test_flip_empty_list_when_no_targets(self):
        m = _bare_model()  # no pool / usage reactions exist
        flipped = flip_reaction_signs_v1_to_v2(m)
        assert flipped == []


# ── integration: build_ecmodel produces v2 sign convention (L2 #4) ────────────

class TestBuildEcmodelV2Signs:
    """L2 scenario #4 — sign convention is enforced end-to-end."""

    def test_pool_exchange_uses_negative_lower_bound(self, toy_gem):
        ec_data = ECModelData()
        # Trigger pool bound computation (ptot*f*sigma*1000) to be non-trivial.
        ec_data.ptot, ec_data.f, ec_data.sigma = 0.5, 0.5, 0.5
        ec_data.kcat_entries = []
        ec_data.uniprot_data = {}

        ecmodel, _, _ = build_ecmodel(toy_gem, ec_data)

        pool_rxn = ecmodel.reactions.get_by_id(ec_data.protein_pool_rxn_id)
        assert pool_rxn.lower_bound < 0
        assert pool_rxn.upper_bound == 0.0

    def test_pool_exchange_stoich_is_minus_one(self, toy_gem):
        ec_data = ECModelData()
        ec_data.kcat_entries = []
        ec_data.uniprot_data = {}

        ecmodel, _, _ = build_ecmodel(toy_gem, ec_data)

        pool_rxn = ecmodel.reactions.get_by_id(ec_data.protein_pool_rxn_id)
        pool_met = ecmodel.metabolites.get_by_id(ec_data.protein_pool_met_id)
        assert pool_rxn.metabolites == {pool_met: -1.0}

    def test_usage_reactions_use_v2_convention(self, toy_gem):
        """Every usage_prot_* reaction obeys v2 sign + bounds."""
        ec_data = ECModelData()
        # Provide one kcat entry so build_ecmodel instantiates at least one
        # usage reaction via _ensure_enzyme_met_and_usage.
        ec_data.kcat_entries = [{
            "proteins": ["P1"],
            "genes": ["g1"],
            "gene_name": "g1",
            "kcat": 10.0,
            "rxns": [],
            "notes": "test",
            "stoicho": [1],
        }]
        ec_data.uniprot_data = {
            "P1": {"gene": "g1", "mw_da": 50000.0, "ec": "", "sequence": ""},
        }

        ecmodel, _, _ = build_ecmodel(toy_gem, ec_data)

        usage_rxns = [r for r in ecmodel.reactions if r.id.startswith("usage_prot_")]
        assert usage_rxns, "build_ecmodel should have created at least one usage rxn"

        for rxn in usage_rxns:
            # v2 bounds.
            assert rxn.lower_bound == -1000.0
            assert rxn.upper_bound == 0.0

            # v2 stoich: one enzyme (-1), one pool (+1), nothing else.
            pool_met = ecmodel.metabolites.get_by_id(ec_data.protein_pool_met_id)
            prot_mets = [m for m in rxn.metabolites if m.id.startswith("prot_")
                         and m.id != "prot_pool"]
            assert len(prot_mets) == 1
            assert rxn.metabolites[prot_mets[0]] == -1.0
            assert rxn.metabolites[pool_met] == 1.0
