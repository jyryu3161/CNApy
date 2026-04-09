"""ECModelData: data class for GECKO enzyme-constrained model metadata."""

from dataclasses import dataclass, field


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
    # mapping { original_rxn_id -> [split_rxn_ids] }  (only populated when splitting occurred)
    split_rxn_map: dict = field(default_factory=dict)

    # ── fixed identifiers ─────────────────────────────────────────────────────
    protein_pool_met_id: str = "prot_pool"
    protein_pool_rxn_id: str = "prot_pool_exchange"

    # ──────────────────────────────────────────────────────────────────────────
    def pool_bound(self) -> float:
        """Upper bound for protein pool exchange (mg/gDCW)."""
        return self.ptot * self.f * self.sigma * 1000

    def to_dict(self) -> dict:
        """Serialise to a plain dict for JSON storage."""
        return {
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
        }

    @classmethod
    def from_dict(cls, data: dict) -> "ECModelData":
        """Deserialise from JSON storage."""
        obj = cls()
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
        return obj

    def reset(self):
        """Reset to default (non-ecModel) state."""
        self.is_ecmodel = False
        self.gecko_light = False
        self.enzyme_met_ids.clear()
        self.enzyme_rxn_ids.clear()
        self.original_reaction_ids.clear()
        self.split_rxn_map.clear()
