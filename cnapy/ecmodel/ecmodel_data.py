"""ECModelData: data class for GECKO enzyme-constrained model metadata."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from cnapy.ecmodel.ec_structure import EcStructure

if TYPE_CHECKING:
    import cobra

# Design Ref: §3.2 — current schema version. v1 files get auto-migrated on load.
CURRENT_SCHEMA_VERSION = 2


@dataclass
class KcatEntry:
    """One row from customKcats.tsv: a user-supplied kcat constraint."""
    proteins: list   # UniProt IDs; multiple = enzyme complex ('+' separated in file)
    genes: list      # gene IDs (informational)
    gene_name: str   # short gene name (informational)
    kcat: float      # turnover number (1/s)
    rxns: list       # reaction IDs to constrain; empty = use protein matching
    notes: str       # source / comment
    stoicho: list    # subunit copies per protein, same length as proteins


@dataclass
class ECModelData:
    """
    Stores all metadata needed to build, track, and revert a GECKO ecModel.

    Populated by gecko_dialog when the user loads CSV files and converts the
    model.  Serialised as ecmodel_data.json inside the .cna project ZIP so
    the ecModel state survives save/load.
    """

    # ── schema version ──────────────────────────────────────────────────────
    # Design Ref: §3.2 FR-20 — absent or 1 in old .cna files triggers automatic
    # v1→v2 sign-convention flip on load.
    schema_version: int = CURRENT_SCHEMA_VERSION

    # ── status ──────────────────────────────────────────────────────────────
    is_ecmodel: bool = False
    gecko_light: bool = False   # True → light formulation (no per-enzyme reactions)

    # ── protein-pool parameters ──────────────────────────────────────────────
    sigma: float = 0.5   # average in-vivo enzyme saturation
    f: float = 0.5       # fraction of metabolic enzymes in total proteome
    ptot: float = 0.5    # total protein content (g/gDCW)

    # ── user-supplied enzyme data ─────────────────────────────────────────────
    # Parsed from customKcats.tsv – list of KcatEntry dicts (serialised as plain dicts)
    kcat_entries: list = field(default_factory=list)

    # Parsed from uniprot.tsv – { uniprot_id: {gene, mw_da, ec, sequence} }
    uniprot_data: dict = field(default_factory=dict)

    # ── internal mappings set by ecmodel_builder ──────────────────────────────
    # { uniprot_id -> prot_<id> metabolite id }
    enzyme_met_ids: dict = field(default_factory=dict)
    # { uniprot_id -> usage_prot_<id> reaction id }
    enzyme_rxn_ids: dict = field(default_factory=dict)
    # original reaction IDs before any splitting
    original_reaction_ids: list = field(default_factory=list)
    # mapping { original_rxn_id -> [split_rxn_ids] }  — reversible fwd/rev split
    split_rxn_map: dict = field(default_factory=dict)
    # mapping { original_rxn_id -> [_EXP_N ids] }     — isozyme split (FR-04)
    # Design Ref: §2.1 — stored separately from split_rxn_map so reversible
    # and isozyme splits don't collide on the same key.
    isozyme_split_map: dict = field(default_factory=dict)

    # ── fixed identifiers ─────────────────────────────────────────────────────
    protein_pool_met_id: str = "prot_pool"
    protein_pool_rxn_id: str = "prot_pool_exchange"

    # ── FR-05 GECKO 3.0 model.ec metadata ─────────────────────────────────────
    # Design Ref: §3.1 — populated by build_ecmodel and persisted to YAML/.cna.
    ec: EcStructure = field(default_factory=EcStructure)

    # ──────────────────────────────────────────────────────────────────────────
    def pool_bound(self) -> float:
        """Upper bound for protein pool exchange (mg/gDCW)."""
        return self.ptot * self.f * self.sigma * 1000

    def to_dict(self) -> dict:
        """Serialise to a plain dict for JSON storage."""
        return {
            "schema_version": self.schema_version,
            "is_ecmodel": self.is_ecmodel,
            "gecko_light": self.gecko_light,
            "sigma": self.sigma,
            "f": self.f,
            "ptot": self.ptot,
            "kcat_entries": [
                {
                    "proteins": e["proteins"],
                    "genes": e["genes"],
                    "gene_name": e["gene_name"],
                    "kcat": e["kcat"],
                    "rxns": e["rxns"],
                    "notes": e["notes"],
                    "stoicho": e["stoicho"],
                }
                for e in self.kcat_entries
            ],
            "uniprot_data": self.uniprot_data,
            "enzyme_met_ids": self.enzyme_met_ids,
            "enzyme_rxn_ids": self.enzyme_rxn_ids,
            "original_reaction_ids": self.original_reaction_ids,
            "split_rxn_map": self.split_rxn_map,
            "isozyme_split_map": self.isozyme_split_map,
            "ec": self.ec.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "ECModelData":
        """Deserialise from JSON storage. Keeps the on-disk schema_version."""
        obj = cls()
        # schema_version defaults to 1 for files written before FR-20.
        obj.schema_version = int(data.get("schema_version", 1))
        obj.is_ecmodel = data.get("is_ecmodel", False)
        obj.gecko_light = data.get("gecko_light", False)
        obj.sigma = data.get("sigma", 0.5)
        obj.f = data.get("f", 0.5)
        obj.ptot = data.get("ptot", 0.5)
        obj.kcat_entries = data.get("kcat_entries", [])
        obj.uniprot_data = data.get("uniprot_data", {})
        obj.enzyme_met_ids = data.get("enzyme_met_ids", {})
        obj.enzyme_rxn_ids = data.get("enzyme_rxn_ids", {})
        obj.original_reaction_ids = data.get("original_reaction_ids", [])
        obj.split_rxn_map = data.get("split_rxn_map", {})
        obj.isozyme_split_map = data.get("isozyme_split_map", {})
        obj.ec = EcStructure.from_dict(data.get("ec"))
        return obj

    @classmethod
    def from_dict_with_migration(
        cls, data: dict, cobra_model: "cobra.Model | None" = None,
    ) -> "tuple[ECModelData, bool]":
        """Load from dict and apply any schema migrations.

        Design Ref: §3.2 FR-20 — v1 files (missing or schema_version == 1)
        get their sign convention flipped via
        :func:`cnapy.ecmodel.expansion.flip_reaction_signs_v1_to_v2`. The
        cobra model (if given) is modified in place to reflect the new
        sign convention.

        Returns
        -------
        (ec_data, was_migrated) : the loaded/upgraded ECModelData and a bool
            indicating whether migration ran. Callers can use this to surface
            a user-facing notice and make a ``.bak`` backup of the source.
        """
        obj = cls.from_dict(data)
        was_migrated = False
        if obj.schema_version < CURRENT_SCHEMA_VERSION:
            obj.upgrade(cobra_model)
            was_migrated = True
        return obj, was_migrated

    def upgrade(self, cobra_model: "cobra.Model | None" = None) -> None:
        """Upgrade this ECModelData (and its cobra_model) to the current schema.

        v1 → v2: flip sign convention on ``usage_prot_*`` and
        ``prot_pool_exchange`` reactions (FR-07 + FR-08).
        """
        if self.schema_version >= CURRENT_SCHEMA_VERSION:
            return

        if self.schema_version == 1:
            if cobra_model is not None and self.is_ecmodel:
                from cnapy.ecmodel.expansion import flip_reaction_signs_v1_to_v2
                flip_reaction_signs_v1_to_v2(cobra_model)
            self.schema_version = 2

    def reset(self):
        """Reset to default (non-ecModel) state.

        Note: ``gecko_light`` is a user-selected build mode (not derived
        state), so it is preserved across resets. Clearing it here would
        make the GUI's Light-radio selection silently ignored whenever
        build_ecmodel is re-run.
        """
        self.is_ecmodel = False
        self.enzyme_met_ids.clear()
        self.enzyme_rxn_ids.clear()
        self.original_reaction_ids.clear()
        self.split_rxn_map.clear()
        self.isozyme_split_map.clear()
        self.ec = EcStructure()
