"""FR-05: ``EcStructure`` — GECKO 3.0 ``model.ec`` metadata container.

Design Ref: §3.1 — mirrors the MATLAB GECKO 3.0 ``model.ec`` substructure
(Nat Protoc 2024, Table 3, p.642) so CNApy ecModels are round-trip
compatible with the ``ecYeastGEM.yml`` / ``human-GEM`` distributions.

Layout
------
Per-reaction arrays (length m = len(rxns)):
    rxns, kcat, source, notes, eccodes
Per-enzyme arrays (length n = len(enzymes)):
    genes, enzymes, mw, sequence, concs
Sparse matrix:
    rxnEnzMat — (m × n) CSR with subunit copies of enzyme j in rxn i.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from scipy.sparse import coo_matrix, csr_matrix


# Design Ref: §10.1 — source labels that round-trip to the GECKO YAML.
SOURCE_LABELS = {"brenda", "dlkcat", "custom", "standard", ""}


@dataclass
class EcStructure:
    """GECKO 3.0 ``model.ec`` structure.

    See module docstring for the overall layout. All list fields are
    aligned column-wise; ``len(rxns) == len(kcat) == len(source) == ...``
    and ``len(genes) == len(enzymes) == len(mw) == ...``.
    """

    # ── per-reaction (m entries) ───────────────────────────────────────────
    rxns:     list[str] = field(default_factory=list)
    kcat:     list[float] = field(default_factory=list)
    source:   list[str] = field(default_factory=list)
    notes:    list[str] = field(default_factory=list)
    eccodes:  list[str] = field(default_factory=list)

    # ── per-enzyme (n entries) ─────────────────────────────────────────────
    genes:    list[str] = field(default_factory=list)
    enzymes:  list[str] = field(default_factory=list)
    mw:       list[float] = field(default_factory=list)
    sequence: list[str] = field(default_factory=list)
    concs:    list[float] = field(default_factory=list)

    # ── sparse matrix ──────────────────────────────────────────────────────
    # Design Ref: §3.1 — entry[i, j] is the subunit-copy count of enzyme j
    # in reaction i. ``None`` means "not populated yet"; use an empty CSR
    # to mean "definitively zero enzymes".
    rxnEnzMat: csr_matrix | None = None

    # ── light/full flag ────────────────────────────────────────────────────
    # Design Ref: §3.1 — parity with MATLAB's ``geckoLight`` boolean.
    gecko_light: bool = False

    # ── size helpers ───────────────────────────────────────────────────────

    def n_rxns(self) -> int:
        return len(self.rxns)

    def n_enzymes(self) -> int:
        return len(self.enzymes)

    # ── validation ────────────────────────────────────────────────────────

    def validate(self) -> None:
        """Sanity-check column alignment.

        Raises
        ------
        ValueError
            if per-reaction or per-enzyme lists are not equal length,
            or the sparse matrix shape disagrees with the array sizes.
        """
        m = self.n_rxns()
        n = self.n_enzymes()

        per_rxn = {
            "kcat": len(self.kcat),
            "source": len(self.source),
            "notes": len(self.notes),
            "eccodes": len(self.eccodes),
        }
        for name, length in per_rxn.items():
            if length != m:
                raise ValueError(
                    f"per-reaction field '{name}' has length {length}, "
                    f"expected {m} (== len(rxns))"
                )

        per_enz = {
            "genes": len(self.genes),
            "mw": len(self.mw),
            "sequence": len(self.sequence),
            "concs": len(self.concs),
        }
        for name, length in per_enz.items():
            if length != n:
                raise ValueError(
                    f"per-enzyme field '{name}' has length {length}, "
                    f"expected {n} (== len(enzymes))"
                )

        if self.rxnEnzMat is not None and self.rxnEnzMat.shape != (m, n):
            raise ValueError(
                f"rxnEnzMat shape {self.rxnEnzMat.shape} does not match "
                f"(n_rxns, n_enzymes) = ({m}, {n})"
            )

    # ── YAML !!omap serialization ─────────────────────────────────────────
    # Design Ref: §3.1 + §8.2 #8-9. The GECKO YAML uses ``!!omap`` (ordered
    # mapping) tags. Both sections are *lists of single-key dicts* at the
    # outer level, and each enzyme entry inside an ec-rxn is also an omap.
    # For CNApy's I/O layer we represent omaps as lists of ``(key, value)``
    # pairs; the yaml_io layer attaches the ``!!omap`` tag on dump.

    def to_yaml_omap(self) -> tuple[list, list]:
        """Return ``(ec_rxns_list, ec_enzymes_list)`` for YAML emission.

        Structure of each ``ec_rxns`` entry::

            [("id", rxn_id),
             ("kcat", kcat_value),
             ("source", source_label),
             ("eccodes", eccodes_value),      # str or list[str]
             ("notes", notes_str),
             ("enzymes", [(uniprot_id, subunit_count), ...])]

        Structure of each ``ec_enzymes`` entry::

            [("genes", gene),
             ("enzymes", uniprot_id),
             ("mw", mw_da),
             ("sequence", seq),
             ("concs", measured_conc_mg_gDCW)]
        """
        self.validate()

        # --- ec_rxns ------------------------------------------------------
        ec_rxns: list[list[tuple[str, Any]]] = []
        for i, rxn_id in enumerate(self.rxns):
            entry: list[tuple[str, Any]] = [
                ("id", rxn_id),
                ("kcat", self.kcat[i]),
                ("source", self.source[i]),
            ]
            # eccodes: YAML convention — single code as str, multiple as list
            # matching GECKO 3.0 output.
            codes = _split_eccodes(self.eccodes[i])
            if len(codes) == 1:
                entry.append(("eccodes", codes[0]))
            elif codes:
                entry.append(("eccodes", codes))
            else:
                entry.append(("eccodes", ""))

            if self.notes[i]:
                entry.append(("notes", self.notes[i]))

            # enzymes omap: { uniprot_id: subunit_count }.
            enz_pairs = self._rxn_enzymes_pairs(i)
            entry.append(("enzymes", enz_pairs))

            ec_rxns.append(entry)

        # --- ec_enzymes ---------------------------------------------------
        ec_enzymes: list[list[tuple[str, Any]]] = []
        for j in range(self.n_enzymes()):
            ec_enzymes.append([
                ("genes", self.genes[j]),
                ("enzymes", self.enzymes[j]),
                ("mw", self.mw[j]),
                ("sequence", self.sequence[j]),
                ("concs", self.concs[j]),
            ])

        return ec_rxns, ec_enzymes

    def _rxn_enzymes_pairs(self, rxn_index: int) -> list[tuple[str, float]]:
        """Extract the enzyme list for reaction ``rxn_index`` as omap pairs.

        Walks the CSR row's ``indices`` / ``data`` arrays directly so we
        stay O(nnz) per row instead of materialising the dense row.
        """
        if self.rxnEnzMat is None:
            return []
        row = self.rxnEnzMat.getrow(rxn_index)
        pairs: list[tuple[str, float]] = []
        for j, stoich in zip(row.indices, row.data):
            val = float(stoich)
            if val == 0:
                continue
            val_out: Any = int(val) if val.is_integer() else val
            pairs.append((self.enzymes[j], val_out))
        return pairs

    @classmethod
    def from_yaml_omap(cls,
                        ec_rxns: list | None,
                        ec_enzymes: list | None,
                        gecko_light: bool = False) -> "EcStructure":
        """Inverse of :meth:`to_yaml_omap`.

        Accepts the relaxed shapes produced by PyYAML's ``safe_load`` with
        either pairs (list of ``(key, value)``) or plain dicts.
        """
        ec_rxns = ec_rxns or []
        ec_enzymes = ec_enzymes or []

        # --- enzymes first so we can index them by UniProt when filling
        # rxnEnzMat below.
        genes: list[str] = []
        enzymes: list[str] = []
        mw: list[float] = []
        sequence: list[str] = []
        concs: list[float] = []

        for entry in ec_enzymes:
            d = _to_dict(entry)
            genes.append(_as_str(d.get("genes", "")))
            enzymes.append(_as_str(d.get("enzymes", "")))
            mw.append(_as_float(d.get("mw", 0.0)))
            sequence.append(_as_str(d.get("sequence", "")))
            concs.append(_as_float(d.get("concs", 0.0)))

        enz_index = {eid: j for j, eid in enumerate(enzymes)}

        # --- reactions ----------------------------------------------------
        rxns: list[str] = []
        kcat: list[float] = []
        source: list[str] = []
        notes: list[str] = []
        eccodes: list[str] = []
        rxn_enz_rows: list[tuple[int, int, float]] = []

        for i, entry in enumerate(ec_rxns):
            d = _to_dict(entry)
            rxns.append(_as_str(d.get("id", "")))
            kcat.append(_as_float(d.get("kcat", 0.0)))
            source.append(_as_str(d.get("source", "")))
            notes.append(_as_str(d.get("notes", "")))
            eccodes.append(_join_eccodes(d.get("eccodes", "")))

            enz_entry = d.get("enzymes", [])
            enz_pairs = _to_pairs(enz_entry)
            for uid, stoich in enz_pairs:
                # Tolerate enzymes mentioned in ec-rxns but missing from
                # ec-enzymes (seen occasionally in partial exports).
                if uid not in enz_index:
                    enz_index[uid] = len(enzymes)
                    enzymes.append(uid)
                    genes.append("")
                    mw.append(0.0)
                    sequence.append("")
                    concs.append(0.0)
                rxn_enz_rows.append((i, enz_index[uid], _as_float(stoich)))

        rxn_enz_mat = _triples_to_csr(rxn_enz_rows, len(rxns), len(enzymes))

        return cls(
            rxns=rxns, kcat=kcat, source=source, notes=notes, eccodes=eccodes,
            genes=genes, enzymes=enzymes, mw=mw, sequence=sequence, concs=concs,
            rxnEnzMat=rxn_enz_mat,
            gecko_light=gecko_light,
        )

    # ── plain-dict round-trip (for .cna JSON storage) ─────────────────────

    def to_dict(self) -> dict:
        """Serialise to a JSON-compatible dict for the .cna archive."""
        self.validate()
        if self.rxnEnzMat is None or self.rxnEnzMat.nnz == 0:
            mat_entries: list[list[Any]] = []
        else:
            coo = self.rxnEnzMat.tocoo()
            mat_entries = [
                [int(r), int(c), float(v)]
                for r, c, v in zip(coo.row, coo.col, coo.data)
            ]
        return {
            "rxns": list(self.rxns),
            "kcat": list(self.kcat),
            "source": list(self.source),
            "notes": list(self.notes),
            "eccodes": list(self.eccodes),
            "genes": list(self.genes),
            "enzymes": list(self.enzymes),
            "mw": list(self.mw),
            "sequence": list(self.sequence),
            "concs": list(self.concs),
            "rxnEnzMat_coo": mat_entries,
            "gecko_light": bool(self.gecko_light),
        }

    @classmethod
    def from_dict(cls, data: dict | None) -> "EcStructure":
        """Inverse of :meth:`to_dict`. Accepts ``None`` → empty EcStructure."""
        if not data:
            return cls()
        obj = cls(
            rxns=list(data.get("rxns", [])),
            kcat=[float(v) for v in data.get("kcat", [])],
            source=[_as_str(v) for v in data.get("source", [])],
            notes=[_as_str(v) for v in data.get("notes", [])],
            eccodes=[_as_str(v) for v in data.get("eccodes", [])],
            genes=[_as_str(v) for v in data.get("genes", [])],
            enzymes=[_as_str(v) for v in data.get("enzymes", [])],
            mw=[float(v) for v in data.get("mw", [])],
            sequence=[_as_str(v) for v in data.get("sequence", [])],
            concs=[float(v) for v in data.get("concs", [])],
            gecko_light=bool(data.get("gecko_light", False)),
        )
        coo_entries = data.get("rxnEnzMat_coo", [])
        triples = [(int(r), int(c), float(v)) for r, c, v in coo_entries]
        obj.rxnEnzMat = _triples_to_csr(triples, obj.n_rxns(), obj.n_enzymes())
        return obj


# ── helpers ────────────────────────────────────────────────────────────────────

def _triples_to_csr(triples: list[tuple[int, int, float]],
                     m: int, n: int) -> csr_matrix | None:
    """Build a CSR matrix from ``(row, col, value)`` triples.

    Returns ``None`` for empty dimensions so callers can keep the "not
    populated" semantics. Goes via COO construction in one shot, avoiding
    the ``lil_matrix`` → incremental assignment → ``tocsr`` round-trip.
    """
    if not triples or not m or not n:
        return None
    rows, cols, data = zip(*triples)
    return coo_matrix((data, (rows, cols)), shape=(m, n), dtype=float).tocsr()


def _to_dict(entry: Any) -> dict:
    """Normalise a YAML !!omap (list of 1-key dicts, or list of pairs, or
    already-built dict) into a plain dict."""
    if isinstance(entry, dict):
        return dict(entry)
    if isinstance(entry, list):
        out: dict = {}
        for item in entry:
            if isinstance(item, dict):
                # single-key dict as produced by yaml.safe_load on !!omap
                out.update(item)
            elif isinstance(item, (list, tuple)) and len(item) == 2:
                out[item[0]] = item[1]
        return out
    return {}


def _to_pairs(entry: Any) -> list[tuple[str, float]]:
    """Normalise an enzymes !!omap into ``[(uniprot_id, stoich), ...]``."""
    if isinstance(entry, dict):
        return [(str(k), _as_float(v)) for k, v in entry.items()]
    if isinstance(entry, list):
        pairs: list[tuple[str, float]] = []
        for item in entry:
            if isinstance(item, dict):
                for k, v in item.items():
                    pairs.append((str(k), _as_float(v)))
            elif isinstance(item, (list, tuple)) and len(item) == 2:
                pairs.append((str(item[0]), _as_float(item[1])))
        return pairs
    return []


def _as_str(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def _as_float(value: Any) -> float:
    if value is None or value == "":
        return 0.0
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _split_eccodes(value: str) -> list[str]:
    """Split our internal ';'-joined eccodes into a YAML-friendly list."""
    if not value:
        return []
    return [code.strip() for code in value.split(";") if code.strip()]


def _join_eccodes(value: Any) -> str:
    """Inverse of :func:`_split_eccodes`. Accepts str or list input."""
    if value is None:
        return ""
    if isinstance(value, list):
        return ";".join(str(code).strip() for code in value if str(code).strip())
    return str(value)
