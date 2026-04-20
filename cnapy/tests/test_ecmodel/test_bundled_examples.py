"""gecko-data-bundle tests — bundled GECKO example datasets.

Design Ref: gecko-data-bundle.plan.md §3.3 — verifies that the Human
and Yeast bundles ship a valid manifest, their referenced files exist,
and each file can be consumed by CNApy's loaders.
"""

from __future__ import annotations

import json

import pytest

from cnapy.data.examples.gecko import (
    SUPPORTED_SPECIES,
    gecko_bundle_dir,
    gecko_bundle_file,
    gecko_bundle_manifest,
)
from cnapy.ecmodel.ecmodel_builder import (
    parse_customkcats_file,
    parse_uniprot_file,
)
from cnapy.ecmodel.yaml_io import load_ecmodel


# ── Manifest + file presence ──────────────────────────────────────────────────

class TestBundleManifest:
    @pytest.mark.parametrize("species", SUPPORTED_SPECIES)
    def test_manifest_json_parses(self, species):
        m = gecko_bundle_manifest(species)
        assert m["species"] == species
        for key in ("display_name", "files", "version", "source",
                    "license_summary"):
            assert key in m, f"manifest missing '{key}'"
        for role in ("model", "kcat", "uniprot"):
            assert role in m["files"], f"manifest.files missing '{role}'"

    @pytest.mark.parametrize("species", SUPPORTED_SPECIES)
    def test_bundled_files_exist(self, species):
        """Every file referenced by manifest.json is physically present."""
        for role in ("model", "kcat", "uniprot"):
            path = gecko_bundle_file(species, role)
            assert path.exists(), f"{species}/{role} → missing {path}"
            assert path.stat().st_size > 0, f"{path} is empty"

    def test_bundle_dir_rejects_unknown_species(self):
        with pytest.raises(ValueError, match="Unknown"):
            gecko_bundle_dir("mouse")


# ── Loader smoke tests ───────────────────────────────────────────────────────

class TestBundleLoaders:
    """Each bundled file must be consumable by the CNApy loaders."""

    @pytest.mark.parametrize("species", SUPPORTED_SPECIES)
    def test_model_yaml_loads(self, species):
        yml_path = gecko_bundle_file(species, "model")
        model, ec = load_ecmodel(str(yml_path))
        assert len(model.reactions) > 0
        # Plain GEM YAML — no ec-rxns section, so is_ecmodel=False.
        assert ec.is_ecmodel is False

    @pytest.mark.parametrize("species", SUPPORTED_SPECIES)
    def test_kcat_tsv_parses(self, species):
        path = gecko_bundle_file(species, "kcat")
        entries = parse_customkcats_file(str(path))
        assert len(entries) > 0
        # Every parsed entry has a finite positive kcat (NaN guard).
        for e in entries:
            assert e["kcat"] > 0

    @pytest.mark.parametrize("species", SUPPORTED_SPECIES)
    def test_uniprot_tsv_parses(self, species):
        path = gecko_bundle_file(species, "uniprot")
        data = parse_uniprot_file(str(path))
        assert len(data) > 0
        # Every loaded protein has a non-negative MW.
        for uid, info in data.items():
            assert info["mw_da"] >= 0, f"{uid}: negative MW"


# ── Manifest coverage fields line up with actual data ────────────────────────

class TestManifestCoverage:
    """If manifest lists expected coverage counts, they should be consistent
    with the bundled files (order-of-magnitude check, not exact)."""

    @pytest.mark.parametrize("species", SUPPORTED_SPECIES)
    def test_kcat_entry_count_plausible(self, species):
        manifest = gecko_bundle_manifest(species)
        expected = manifest.get("coverage", {}).get("kcat_entries")
        if expected is None:
            pytest.skip("manifest has no kcat_entries coverage hint")
        actual = len(parse_customkcats_file(str(gecko_bundle_file(species, "kcat"))))
        # Allow ±10% drift before tests fail — catches regressions without
        # being brittle to minor dataset refreshes.
        assert 0.9 * expected <= actual <= 1.1 * expected, (
            f"{species}: manifest says {expected} kcat entries, found {actual}"
        )
