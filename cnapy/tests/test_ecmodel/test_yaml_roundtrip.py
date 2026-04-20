"""FR-19 tests — GECKO 3.0 YAML save / load round-trip.

Design Ref: §8.2 #8-10 (EcStructure omap + yaml_io round-trip),
§8.4 #1 (GECKO-compatible yeast-GEM round-trip — scoped to a toy model
since the full yeast-GEM fixture is ~40 MB).
"""

from __future__ import annotations

import cobra
import numpy as np
import pytest
from scipy.sparse import csr_matrix

from cnapy.ecmodel.ec_structure import EcStructure
from cnapy.ecmodel.ecmodel_builder import build_ecmodel
from cnapy.ecmodel.ecmodel_data import ECModelData
from cnapy.ecmodel.exceptions import EcYamlError
from cnapy.ecmodel.yaml_io import load_ecmodel, save_ecmodel


# ── EcStructure !!omap ─────────────────────────────────────────────────────────

class TestEcStructureOmap:
    """Design §8.2 #8-9 — EcStructure to/from !!omap structure."""

    def _sample(self) -> EcStructure:
        return EcStructure(
            rxns=["r1_EXP_1", "r1_EXP_2"],
            kcat=[100.0, 200.0],
            source=["brenda", "dlkcat"],
            notes=["", "auto"],
            eccodes=["1.2.3.4", "1.2.3.4;1.2.3.5"],
            genes=["gA", "gB"],
            enzymes=["P001", "P002"],
            mw=[50000.0, 60000.0],
            sequence=["MSEQA", "MSEQB"],
            concs=[0.0, 12.5],
            rxnEnzMat=csr_matrix(np.array([[1, 0], [0, 2]], dtype=float)),
        )

    def test_to_omap_produces_expected_shape(self):
        ec = self._sample()
        ec_rxns, ec_enzymes = ec.to_yaml_omap()
        assert len(ec_rxns) == 2
        assert len(ec_enzymes) == 2

    def test_ec_rxns_entry_keys_are_ordered(self):
        """GECKO canonical order: id, kcat, source, eccodes[, notes], enzymes."""
        ec = self._sample()
        ec_rxns, _ = ec.to_yaml_omap()
        keys0 = [k for k, _ in ec_rxns[0]]
        assert keys0[:4] == ["id", "kcat", "source", "eccodes"]
        assert keys0[-1] == "enzymes"

    def test_eccodes_single_emits_as_string(self):
        ec = self._sample()
        ec_rxns, _ = ec.to_yaml_omap()
        entry0 = dict(ec_rxns[0])
        assert entry0["eccodes"] == "1.2.3.4"

    def test_eccodes_multi_emits_as_list(self):
        ec = self._sample()
        ec_rxns, _ = ec.to_yaml_omap()
        entry1 = dict(ec_rxns[1])
        assert entry1["eccodes"] == ["1.2.3.4", "1.2.3.5"]

    def test_enzymes_subunit_preserved(self):
        ec = self._sample()
        ec_rxns, _ = ec.to_yaml_omap()
        enz1 = dict(ec_rxns[1])["enzymes"]
        # [(uid, subunit_count), ...] — r1_EXP_2 has P002 with 2 subunits.
        assert enz1 == [("P002", 2)]

    def test_from_omap_inverse(self):
        """Round-trip: to → from preserves all fields."""
        ec = self._sample()
        ec_rxns, ec_enzymes = ec.to_yaml_omap()
        ec2 = EcStructure.from_yaml_omap(ec_rxns, ec_enzymes)

        assert ec2.rxns == ec.rxns
        assert ec2.kcat == ec.kcat
        assert ec2.source == ec.source
        assert ec2.enzymes == ec.enzymes
        assert ec2.mw == ec.mw
        assert ec2.sequence == ec.sequence

        # rxnEnzMat shape + non-zero pattern preserved.
        assert ec2.rxnEnzMat is not None
        np.testing.assert_array_equal(
            ec2.rxnEnzMat.toarray(), ec.rxnEnzMat.toarray()
        )

    def test_from_omap_tolerates_missing_ec_enzymes(self):
        """Enzymes mentioned in ec-rxns but absent from ec-enzymes still
        produce placeholder enzyme rows (accepted permissively)."""
        ec_rxns = [
            [("id", "rxX"), ("kcat", 10.0), ("source", "custom"),
             ("eccodes", ""), ("enzymes", [("P999", 1)])],
        ]
        ec2 = EcStructure.from_yaml_omap(ec_rxns, [])
        assert "P999" in ec2.enzymes
        assert ec2.mw[ec2.enzymes.index("P999")] == 0.0

    def test_validate_catches_misaligned_columns(self):
        ec = EcStructure(
            rxns=["r1", "r2"], kcat=[1.0], source=["", ""],
            notes=["", ""], eccodes=["", ""],
        )
        with pytest.raises(ValueError, match="kcat"):
            ec.validate()

    def test_plain_dict_round_trip(self):
        """to_dict → from_dict preserves CSR entries for .cna storage."""
        ec = self._sample()
        data = ec.to_dict()
        restored = EcStructure.from_dict(data)
        assert restored.rxns == ec.rxns
        assert restored.kcat == ec.kcat
        assert restored.rxnEnzMat is not None
        np.testing.assert_array_equal(
            restored.rxnEnzMat.toarray(), ec.rxnEnzMat.toarray()
        )


# ── yaml_io save/load round-trip ──────────────────────────────────────────────

def _toy_gem_with_isozyme() -> cobra.Model:
    m = cobra.Model("toy_rt")
    A = cobra.Metabolite("A", compartment="c")
    B = cobra.Metabolite("B", compartment="c")
    m.add_metabolites([A, B])

    r1 = cobra.Reaction("r1")
    r1.add_metabolites({A: -1, B: 1})
    r1.bounds = (0, 1000)
    r1.gene_reaction_rule = "gA or gB"

    r_in = cobra.Reaction("EX_A")
    r_in.add_metabolites({A: -1})
    r_in.bounds = (-10, 1000)

    r_out = cobra.Reaction("EX_B")
    r_out.add_metabolites({B: -1})
    r_out.bounds = (0, 1000)

    m.add_reactions([r1, r_in, r_out])
    m.objective = "EX_B"
    return m


def _toy_ec_data() -> ECModelData:
    ec = ECModelData()
    ec.sigma = 0.4
    ec.f = 0.3
    ec.ptot = 0.5
    ec.kcat_entries = [
        {"proteins": ["P001"], "genes": ["gA"], "gene_name": "gA",
         "kcat": 100.0, "rxns": ["r1"], "notes": "brenda", "stoicho": [1]},
        {"proteins": ["P002"], "genes": ["gB"], "gene_name": "gB",
         "kcat": 200.0, "rxns": ["r1"], "notes": "custom", "stoicho": [1]},
    ]
    ec.uniprot_data = {
        "P001": {"gene": "gA", "mw_da": 50000.0, "ec": "", "sequence": "MSEQ"},
        "P002": {"gene": "gB", "mw_da": 60000.0, "ec": "", "sequence": "MSEQ2"},
    }
    return ec


class TestSaveLoadRoundTrip:
    """Design §8.4 #1 — save → load → bit-for-bit equivalence (toy scope)."""

    def test_save_creates_yaml_file(self, tmp_path):
        model = _toy_gem_with_isozyme()
        ec = _toy_ec_data()
        built, _, _ = build_ecmodel(model, ec)

        path = tmp_path / "out.yml"
        save_ecmodel(built, ec, path)

        assert path.exists()
        text = path.read_text(encoding="utf-8")
        # Accept either of PyYAML's valid emissions for an !!omap-tagged root:
        #   "--- !!omap\n..."   (explicit_start with tag inline)
        #   "---\n!!omap\n..."  (tag on its own line)
        assert text.startswith("---") and "!!omap" in text.split("\n", 3)[1] + text.split("\n", 3)[0]
        # Required top-level sections present.
        for section in ("metaData", "metabolites", "reactions",
                         "genes", "compartments", "ec-rxns", "ec-enzymes"):
            assert section in text

    def test_save_refuses_non_ecmodel(self, tmp_path):
        """Saving a non-ecmodel ECModelData must raise EcYamlError."""
        model = _toy_gem_with_isozyme()
        ec = ECModelData()   # is_ecmodel=False

        with pytest.raises(EcYamlError, match="ecModel"):
            save_ecmodel(model, ec, tmp_path / "no.yml")

    def test_save_refuses_missing_directory(self, tmp_path):
        model = _toy_gem_with_isozyme()
        ec = _toy_ec_data()
        built, _, _ = build_ecmodel(model, ec)

        with pytest.raises(EcYamlError, match="directory"):
            save_ecmodel(built, ec, tmp_path / "missing_dir" / "out.yml")

    def test_load_missing_file(self, tmp_path):
        with pytest.raises(EcYamlError, match="not found"):
            load_ecmodel(tmp_path / "nope.yml")

    def test_load_invalid_yaml(self, tmp_path):
        bad = tmp_path / "bad.yml"
        bad.write_text("this is not: - valid: yaml: at all:\n")
        with pytest.raises(EcYamlError):
            load_ecmodel(bad)

    def test_load_missing_required_section(self, tmp_path):
        """Missing 'metabolites' / 'reactions' must surface a clear error."""
        path = tmp_path / "no_rxns.yml"
        path.write_text("""---
!!omap
- metaData: !!omap
  - id: dummy
""")
        with pytest.raises(EcYamlError, match="metabolites|reactions"):
            load_ecmodel(path)

    def test_round_trip_preserves_reactions_and_bounds(self, tmp_path):
        model = _toy_gem_with_isozyme()
        ec = _toy_ec_data()
        built, _, _ = build_ecmodel(model, ec)

        path = tmp_path / "rt.yml"
        save_ecmodel(built, ec, path)
        reloaded, _ec_data = load_ecmodel(path)

        # Reaction set parity.
        assert sorted(r.id for r in reloaded.reactions) == sorted(
            r.id for r in built.reactions
        )
        # Metabolite set parity.
        assert sorted(m.id for m in reloaded.metabolites) == sorted(
            m.id for m in built.metabolites
        )
        # Pool exchange bounds preserved (v2 sign convention).
        pool_before = built.reactions.get_by_id("prot_pool_exchange")
        pool_after = reloaded.reactions.get_by_id("prot_pool_exchange")
        assert pool_after.lower_bound == pytest.approx(pool_before.lower_bound)
        assert pool_after.upper_bound == pytest.approx(pool_before.upper_bound)
        # Pool stoichiometry preserved.
        assert pool_after.metabolites[
            pool_after.metabolites.__iter__().__next__()
        ] == -1.0

    def test_round_trip_preserves_reaction_metabolites(self, tmp_path):
        model = _toy_gem_with_isozyme()
        ec = _toy_ec_data()
        built, _, _ = build_ecmodel(model, ec)

        path = tmp_path / "rt.yml"
        save_ecmodel(built, ec, path)
        reloaded, _ = load_ecmodel(path)

        # Stoichiometry for each reaction should round-trip.
        for rxn in built.reactions:
            reloaded_rxn = reloaded.reactions.get_by_id(rxn.id)
            before = {m.id: c for m, c in rxn.metabolites.items()}
            after = {m.id: c for m, c in reloaded_rxn.metabolites.items()}
            assert after == before, f"stoichiometry mismatch on {rxn.id}"

    def test_round_trip_preserves_gpr(self, tmp_path):
        model = _toy_gem_with_isozyme()
        ec = _toy_ec_data()
        built, _, _ = build_ecmodel(model, ec)

        path = tmp_path / "rt.yml"
        save_ecmodel(built, ec, path)
        reloaded, _ = load_ecmodel(path)

        for rxn in built.reactions:
            reloaded_rxn = reloaded.reactions.get_by_id(rxn.id)
            assert reloaded_rxn.gene_reaction_rule == rxn.gene_reaction_rule

    def test_round_trip_preserves_ec_metadata(self, tmp_path):
        model = _toy_gem_with_isozyme()
        ec = _toy_ec_data()
        built, _, _ = build_ecmodel(model, ec)

        path = tmp_path / "rt.yml"
        save_ecmodel(built, ec, path)
        _, reloaded_ec = load_ecmodel(path)

        assert reloaded_ec.ec.rxns == ec.ec.rxns
        assert reloaded_ec.ec.enzymes == ec.ec.enzymes
        # kcats survive as floats.
        assert reloaded_ec.ec.kcat == pytest.approx(ec.ec.kcat)
        # MW survives.
        assert reloaded_ec.ec.mw == pytest.approx(ec.ec.mw)
        # schema_version stamped correctly.
        assert reloaded_ec.schema_version == 2
        assert reloaded_ec.is_ecmodel is True

    def test_round_trip_preserves_pool_parameters(self, tmp_path):
        model = _toy_gem_with_isozyme()
        ec = _toy_ec_data()
        built, _, _ = build_ecmodel(model, ec)

        path = tmp_path / "rt.yml"
        save_ecmodel(built, ec, path)
        _, reloaded_ec = load_ecmodel(path)

        assert reloaded_ec.sigma == pytest.approx(ec.sigma)
        assert reloaded_ec.f == pytest.approx(ec.f)
        assert reloaded_ec.ptot == pytest.approx(ec.ptot)

    def test_round_trip_twice_is_stable(self, tmp_path):
        """save → load → save → load — second load should equal first load."""
        model = _toy_gem_with_isozyme()
        ec = _toy_ec_data()
        built, _, _ = build_ecmodel(model, ec)

        p1 = tmp_path / "rt1.yml"
        p2 = tmp_path / "rt2.yml"

        save_ecmodel(built, ec, p1)
        m1, ec1 = load_ecmodel(p1)
        save_ecmodel(m1, ec1, p2)
        m2, ec2 = load_ecmodel(p2)

        # ec metadata stable.
        assert ec2.ec.rxns == ec1.ec.rxns
        assert ec2.ec.kcat == pytest.approx(ec1.ec.kcat)
        assert ec2.ec.enzymes == ec1.ec.enzymes
        assert ec2.ec.mw == pytest.approx(ec1.ec.mw)
        # Reaction stoich stable.
        for rxn in m1.reactions:
            before = {m.id: c for m, c in rxn.metabolites.items()}
            after = {m.id: c
                     for m, c in m2.reactions.get_by_id(rxn.id).metabolites.items()}
            assert after == before
