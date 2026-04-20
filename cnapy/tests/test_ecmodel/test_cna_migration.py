"""FR-20 tests — v1 .cna ECModelData migration to v2 sign convention.

Design Ref: §8.2 scenario 11, §8.4 scenario 2.

``ECModelData.from_dict_with_migration`` loads a legacy v1 dict and, if a
cobra model is supplied alongside, flips the sign convention on its
``prot_pool_exchange`` and ``usage_prot_*`` reactions via
``flip_reaction_signs_v1_to_v2``.
"""

from __future__ import annotations

import cobra

from cnapy.ecmodel.ecmodel_data import (
    CURRENT_SCHEMA_VERSION,
    ECModelData,
)


# ── helpers ────────────────────────────────────────────────────────────────────

def _build_v1_ecmodel() -> cobra.Model:
    """Cobra model that mimics a CNApy v1 ecModel on disk.

    v1 convention:
      * ``prot_pool_exchange`` stoich ``{pool: +1}``, bounds ``(0, pool_ub)``.
      * ``usage_prot_{id}``    stoich ``{pool: -1, prot: +1}``, bounds ``(0, usage_ub)``.
    """
    m = cobra.Model("v1_ecmodel")
    a = cobra.Metabolite("A", compartment="c")
    b = cobra.Metabolite("B", compartment="c")
    pool = cobra.Metabolite("prot_pool", name="Protein pool", compartment="c")
    prot = cobra.Metabolite("prot_P1", name="Enzyme P1", compartment="c")
    m.add_metabolites([a, b, pool, prot])

    r = cobra.Reaction("r1")
    r.add_metabolites({a: -1, b: 1})
    r.bounds = (0, 1000)

    pool_rxn = cobra.Reaction(
        "prot_pool_exchange", name="Protein pool exchange",
        lower_bound=0.0, upper_bound=100.0,
    )
    pool_rxn.add_metabolites({pool: 1.0})

    usage_rxn = cobra.Reaction(
        "usage_prot_P1", name="Enzyme usage P1",
        lower_bound=0.0, upper_bound=1000.0,
    )
    usage_rxn.add_metabolites({pool: -1.0, prot: 1.0})

    m.add_reactions([r, pool_rxn, usage_rxn])
    return m


def _v1_dict(**overrides) -> dict:
    """Produce a dict shaped like ecmodel_data.json from a v1 .cna file.

    v1 files were written before FR-20 → ``schema_version`` key is absent.
    """
    base = {
        "is_ecmodel": True,
        "gecko_light": False,
        "sigma": 0.5,
        "f": 0.5,
        "ptot": 0.5,
        "kcat_entries": [],
        "uniprot_data": {
            "P1": {"gene": "g1", "mw_da": 50000.0, "ec": "", "sequence": ""},
        },
        "enzyme_met_ids": {"P1": "prot_P1"},
        "enzyme_rxn_ids": {"P1": "usage_prot_P1"},
        "original_reaction_ids": ["r1"],
        "split_rxn_map": {},
        # FR-20 note: v1 predates isozyme_split_map entirely.
    }
    base.update(overrides)
    return base


# ── from_dict (baseline, schema_version default) ──────────────────────────────

class TestFromDictSchemaVersion:
    def test_missing_schema_version_defaults_to_v1(self):
        """v1 .cna files have no ``schema_version`` key."""
        data = _v1_dict()
        obj = ECModelData.from_dict(data)
        assert obj.schema_version == 1

    def test_explicit_v2_is_kept(self):
        data = _v1_dict(schema_version=2)
        obj = ECModelData.from_dict(data)
        assert obj.schema_version == 2

    def test_future_schema_is_kept(self):
        """A future schema version must not be silently downgraded."""
        data = _v1_dict(schema_version=99)
        obj = ECModelData.from_dict(data)
        assert obj.schema_version == 99


# ── from_dict_with_migration ──────────────────────────────────────────────────

class TestFromDictWithMigration:
    """Load legacy v1 data + cobra model → upgraded in place."""

    def test_was_migrated_flag_is_true_for_v1(self):
        data = _v1_dict()
        model = _build_v1_ecmodel()

        obj, was_migrated = ECModelData.from_dict_with_migration(data, model)

        assert was_migrated is True
        assert obj.schema_version == CURRENT_SCHEMA_VERSION

    def test_was_migrated_flag_is_false_for_v2(self):
        data = _v1_dict(schema_version=2)
        model = _build_v1_ecmodel()   # model state doesn't matter; upgrade is skipped

        _, was_migrated = ECModelData.from_dict_with_migration(data, model)

        assert was_migrated is False

    def test_migration_flips_pool_signs(self):
        """Post-migration pool_rxn obeys v2 convention (Design §8.2 #11)."""
        data = _v1_dict()
        model = _build_v1_ecmodel()

        ECModelData.from_dict_with_migration(data, model)

        pool_rxn = model.reactions.get_by_id("prot_pool_exchange")
        assert pool_rxn.lower_bound < 0
        assert pool_rxn.upper_bound == 0.0

        pool_met = model.metabolites.get_by_id("prot_pool")
        assert pool_rxn.metabolites == {pool_met: -1.0}

    def test_migration_flips_usage_signs(self):
        data = _v1_dict()
        model = _build_v1_ecmodel()

        ECModelData.from_dict_with_migration(data, model)

        usage_rxn = model.reactions.get_by_id("usage_prot_P1")
        pool_met = model.metabolites.get_by_id("prot_pool")
        prot_met = model.metabolites.get_by_id("prot_P1")

        # v2 bounds.
        assert usage_rxn.lower_bound == -1000.0
        assert usage_rxn.upper_bound == 0.0

        # v2 stoich: {prot: -1, pool: +1}.
        assert usage_rxn.metabolites[prot_met] == -1.0
        assert usage_rxn.metabolites[pool_met] == 1.0

    def test_migration_preserves_non_ec_reactions(self):
        """r1 (no prot_*) must not be touched by migration."""
        data = _v1_dict()
        model = _build_v1_ecmodel()
        r1 = model.reactions.get_by_id("r1")
        original_bounds = r1.bounds
        original_stoich = {m.id: c for m, c in r1.metabolites.items()}

        ECModelData.from_dict_with_migration(data, model)

        r1_after = model.reactions.get_by_id("r1")
        assert r1_after.bounds == original_bounds
        assert {m.id: c for m, c in r1_after.metabolites.items()} == original_stoich

    def test_migration_skipped_when_not_ecmodel(self):
        """Non-ecModel projects (is_ecmodel=False) skip the cobra flip.

        The schema still bumps to v2 (metadata-only), but the cobra model
        must not be mutated because it has no prot_pool / usage reactions.
        """
        data = _v1_dict(is_ecmodel=False)
        # Plain GEM without any pool / usage reactions.
        plain = cobra.Model("plain")
        a = cobra.Metabolite("A", compartment="c")
        b = cobra.Metabolite("B", compartment="c")
        plain.add_metabolites([a, b])
        r = cobra.Reaction("r1")
        r.add_metabolites({a: -1, b: 1})
        r.bounds = (0, 1000)
        plain.add_reactions([r])

        obj, was_migrated = ECModelData.from_dict_with_migration(data, plain)

        assert was_migrated is True  # schema bump did run
        assert obj.schema_version == CURRENT_SCHEMA_VERSION
        # Plain cobra model untouched.
        assert plain.reactions.r1.bounds == (0, 1000)

    def test_migration_without_model_only_bumps_schema(self):
        """Callers may pass cobra_model=None for metadata-only upgrade."""
        data = _v1_dict()

        obj, was_migrated = ECModelData.from_dict_with_migration(data, None)

        assert was_migrated is True
        assert obj.schema_version == CURRENT_SCHEMA_VERSION


# ── ECModelData.upgrade ───────────────────────────────────────────────────────

class TestECModelDataUpgrade:
    def test_upgrade_from_v1_flips_signs(self):
        """Direct call on an ECModelData instance with v1 data."""
        obj = ECModelData()
        obj.schema_version = 1
        obj.is_ecmodel = True
        model = _build_v1_ecmodel()

        obj.upgrade(model)

        assert obj.schema_version == CURRENT_SCHEMA_VERSION
        pool_rxn = model.reactions.get_by_id("prot_pool_exchange")
        assert pool_rxn.lower_bound < 0
        assert pool_rxn.upper_bound == 0.0

    def test_upgrade_is_noop_for_current_version(self):
        """upgrade() on a v2 object must not touch the cobra model."""
        obj = ECModelData()   # schema_version == CURRENT_SCHEMA_VERSION
        model = _build_v1_ecmodel()

        before_pool_bounds = model.reactions.get_by_id(
            "prot_pool_exchange").bounds
        before_usage_bounds = model.reactions.get_by_id(
            "usage_prot_P1").bounds

        obj.upgrade(model)

        # Not flipped (we would have flipped if upgrade ran).
        assert model.reactions.get_by_id(
            "prot_pool_exchange").bounds == before_pool_bounds
        assert model.reactions.get_by_id(
            "usage_prot_P1").bounds == before_usage_bounds

    def test_upgrade_is_idempotent(self):
        """Calling upgrade() twice is safe and has no extra effect."""
        obj = ECModelData()
        obj.schema_version = 1
        obj.is_ecmodel = True
        model = _build_v1_ecmodel()

        obj.upgrade(model)
        first_pool_bounds = model.reactions.get_by_id(
            "prot_pool_exchange").bounds
        first_usage_bounds = model.reactions.get_by_id(
            "usage_prot_P1").bounds

        # Second call — schema_version is now 2, so this is a no-op.
        obj.upgrade(model)

        assert obj.schema_version == CURRENT_SCHEMA_VERSION
        assert model.reactions.get_by_id(
            "prot_pool_exchange").bounds == first_pool_bounds
        assert model.reactions.get_by_id(
            "usage_prot_P1").bounds == first_usage_bounds

    def test_upgrade_without_model_still_bumps_schema(self):
        obj = ECModelData()
        obj.schema_version = 1
        obj.is_ecmodel = True

        obj.upgrade(None)

        assert obj.schema_version == CURRENT_SCHEMA_VERSION
