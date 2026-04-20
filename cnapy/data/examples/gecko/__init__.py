"""Bundled GECKO 3.0 example data (Human-GEM, yeast-GEM).

Each species subdirectory carries a ``manifest.json`` describing file
roles and metadata. Access the bundle via :func:`gecko_bundle_dir`.
"""

from __future__ import annotations

import json
from importlib import resources
from pathlib import Path


SUPPORTED_SPECIES = ("human", "yeast")


def gecko_bundle_dir(species: str) -> Path:
    """Return the on-disk path to the bundled GECKO example for ``species``.

    Works both for editable/dev clones and pip-installed wheels because
    ``importlib.resources.as_file`` transparently extracts the package
    data when needed.
    """
    if species not in SUPPORTED_SPECIES:
        raise ValueError(
            f"Unknown GECKO example species {species!r}. "
            f"Supported: {SUPPORTED_SPECIES}"
        )
    with resources.as_file(
        resources.files("cnapy.data.examples.gecko") / species
    ) as path:
        return Path(path)


def gecko_bundle_manifest(species: str) -> dict:
    """Load and return the ``manifest.json`` for ``species``."""
    return json.loads((gecko_bundle_dir(species) / "manifest.json").read_text())


def gecko_bundle_file(species: str, role: str) -> Path:
    """Return the bundled file path for a given role (``model``/``kcat``/``uniprot``)."""
    manifest = gecko_bundle_manifest(species)
    filename = manifest["files"][role]
    return gecko_bundle_dir(species) / filename
