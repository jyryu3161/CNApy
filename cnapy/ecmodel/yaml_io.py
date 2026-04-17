"""FR-19: GECKO 3.0 compatible YAML save/load for CNApy ecModels.

Design Ref: §2.2 + §3.1 + §5.1 — serialise a ``cobra.Model`` together with
its ``ECModelData`` to the ``!!omap`` layout used by SysBioChalmers
(``ecYeastGEM.yml``, ``human-GEM.yml`` and friends).

Top-level YAML sections (in order)::

    - metaData
    - metabolites
    - reactions
    - genes
    - compartments
    - ec-rxns      (from ``ec_data.ec``)
    - ec-enzymes   (from ``ec_data.ec``)

CNApy-specific notes
--------------------
*   ``usage_prot_*`` and ``prot_pool_exchange`` reactions are serialised
    as regular cobra reactions with their v2 sign convention preserved.
*   ``yaml.safe_load`` is used for parsing (Design §7) — no arbitrary
    Python object instantiation is possible from a malicious YAML.
*   ``!!omap`` is represented internally with ``yaml.add_representer`` so
    ordered dicts round-trip cleanly; on load the reverse resolver maps
    ``!!omap`` back to plain dicts for ergonomic consumption.
"""

from __future__ import annotations

import logging
from collections import OrderedDict
from pathlib import Path
from typing import Any

import cobra
import yaml

from cnapy.ecmodel.ec_structure import EcStructure, _as_float, _to_dict
from cnapy.ecmodel.ecmodel_data import CURRENT_SCHEMA_VERSION, ECModelData
from cnapy.ecmodel.exceptions import EcYamlError

logger = logging.getLogger("cnapy.ecmodel.yaml_io")


# ── !!omap representer / constructor ──────────────────────────────────────────
# PyYAML emits plain ``dict``s as unordered mappings; GECKO requires ordered
# maps tagged ``!!omap``. We register a custom representer for OrderedDict
# and a constructor that turns ``!!omap`` back into a plain ``dict`` (Python
# 3.7+ dicts preserve insertion order, matching GECKO's ordering guarantees).


def _ordered_dict_representer(dumper: yaml.SafeDumper, data: OrderedDict):
    # YAML `!!omap` is spec'd as a *sequence of single-key mappings*. We
    # build that sequence explicitly so the emitted file matches GECKO's
    # ``ecYeastGEM.yml`` layout line-for-line.
    return dumper.represent_sequence(
        "tag:yaml.org,2002:omap",
        [{k: v} for k, v in data.items()],
    )


def _omap_constructor(loader: yaml.SafeLoader, node):
    """Constructor for `!!omap` — tolerant of both sequence and mapping nodes.

    Sequence form (GECKO canonical): list of single-key mappings.
    Mapping form: some producers emit plain mappings with the ``!!omap`` tag.
    """
    if isinstance(node, yaml.SequenceNode):
        return _to_dict(loader.construct_sequence(node, deep=True))
    if isinstance(node, yaml.MappingNode):
        return dict(loader.construct_mapping(node, deep=True))
    raise yaml.constructor.ConstructorError(
        None, None,
        f"unexpected node type {type(node).__name__} for !!omap", node.start_mark,
    )


class _EcDumper(yaml.SafeDumper):
    """SafeDumper with OrderedDict → ``!!omap`` support."""


class _EcLoader(yaml.SafeLoader):
    """SafeLoader that resolves ``!!omap`` to plain dict."""


_EcDumper.add_representer(OrderedDict, _ordered_dict_representer)
_EcLoader.add_constructor("tag:yaml.org,2002:omap", _omap_constructor)


# ── metabolite / reaction / gene → omap serialisers ───────────────────────────

def _omap(pairs: list[tuple[str, Any]]) -> OrderedDict:
    return OrderedDict(pairs)


def _metabolite_to_omap(met: cobra.Metabolite) -> OrderedDict:
    pairs: list[tuple[str, Any]] = [
        ("id", met.id),
        ("name", met.name or ""),
        ("compartment", met.compartment or ""),
    ]
    if met.formula:
        pairs.append(("formula", met.formula))
    if met.charge is not None:
        pairs.append(("charge", met.charge))
    if met.annotation:
        pairs.append(("annotation", _omap(sorted(met.annotation.items()))))
    return _omap(pairs)


def _reaction_to_omap(rxn: cobra.Reaction) -> OrderedDict:
    met_pairs = [(m.id, _as_number(coeff))
                 for m, coeff in rxn.metabolites.items()]
    pairs: list[tuple[str, Any]] = [
        ("id", rxn.id),
        ("name", rxn.name or ""),
        ("metabolites", _omap(met_pairs)),
        ("lower_bound", _as_number(rxn.lower_bound)),
        ("upper_bound", _as_number(rxn.upper_bound)),
        ("gene_reaction_rule", rxn.gene_reaction_rule or ""),
    ]
    if rxn.subsystem:
        pairs.append(("subsystem", rxn.subsystem))
    if rxn.annotation:
        pairs.append(("annotation", _omap(sorted(rxn.annotation.items()))))
    return _omap(pairs)


def _gene_to_omap(gene: cobra.Gene) -> OrderedDict:
    pairs: list[tuple[str, Any]] = [("id", gene.id)]
    if gene.name and gene.name != gene.id:
        pairs.append(("name", gene.name))
    if gene.annotation:
        pairs.append(("annotation", _omap(sorted(gene.annotation.items()))))
    return _omap(pairs)


def _ec_rxns_to_omap(ec: EcStructure) -> list[OrderedDict]:
    """Convert EcStructure to list-of-omap for the ``ec-rxns`` section."""
    raw_rxns, _ = ec.to_yaml_omap()
    out: list[OrderedDict] = []
    for entry in raw_rxns:
        # entry is a list of (k, v) pairs; nested 'enzymes' is also pair-list.
        pairs: list[tuple[str, Any]] = []
        for key, value in entry:
            if key == "enzymes" and isinstance(value, list):
                pairs.append(("enzymes", _omap(value)))
            else:
                pairs.append((key, value))
        out.append(_omap(pairs))
    return out


def _ec_enzymes_to_omap(ec: EcStructure) -> list[OrderedDict]:
    _, raw_enz = ec.to_yaml_omap()
    return [_omap(entry) for entry in raw_enz]


def _as_number(value: Any) -> Any:
    """Emit ints as ints and everything else as float (GECKO convention)."""
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value
    try:
        f = float(value)
    except (TypeError, ValueError):
        return value
    if f.is_integer() and abs(f) < 1e15:
        return int(f)
    return f


# ── save ───────────────────────────────────────────────────────────────────────

def save_ecmodel(model: cobra.Model,
                  ec_data: ECModelData,
                  path: str | Path) -> None:
    """Serialise (``model``, ``ec_data``) to a GECKO-compatible YAML file.

    Parameters
    ----------
    model    : the cobra.Model containing metabolites / reactions / genes.
    ec_data  : ECModelData. Must have ``is_ecmodel=True`` and a populated
               ``ec`` structure (done by :func:`build_ecmodel`).
    path     : output file path. Parent directory must exist.

    Raises
    ------
    EcYamlError
        if ``ec_data`` is not an active ecModel, or if writing fails.
    """
    if not ec_data.is_ecmodel:
        raise EcYamlError(
            "Cannot export YAML: ECModelData is not an active ecModel. "
            "Run build_ecmodel() first."
        )

    out_path = Path(path)
    if not out_path.parent.exists():
        raise EcYamlError(
            f"Output directory does not exist: {out_path.parent}"
        )

    # Build the nested !!omap document.
    meta = _omap([
        ("id", model.id or ""),
        ("name", model.name or model.id or ""),
        ("schema_version", CURRENT_SCHEMA_VERSION),
        ("geckoLight", "true" if ec_data.gecko_light else "false"),
        ("sigma", _as_number(ec_data.sigma)),
        ("f", _as_number(ec_data.f)),
        ("ptot", _as_number(ec_data.ptot)),
    ])

    compartments = _omap(sorted(model.compartments.items()))

    document: list[tuple[str, Any]] = [
        ("metaData", meta),
        ("metabolites", [_metabolite_to_omap(m) for m in model.metabolites]),
        ("reactions", [_reaction_to_omap(r) for r in model.reactions]),
        ("genes", [_gene_to_omap(g) for g in model.genes]),
        ("compartments", compartments),
        ("ec-rxns", _ec_rxns_to_omap(ec_data.ec)),
        ("ec-enzymes", _ec_enzymes_to_omap(ec_data.ec)),
    ]
    # Top-level document is itself an ``!!omap`` — an OrderedDict routes
    # through ``_ordered_dict_representer`` so the output matches GECKO's
    # ``ecYeastGEM.yml`` layout exactly.
    root = OrderedDict(document)

    with out_path.open("w", encoding="utf-8") as fp:
        yaml.dump(
            root, fp, Dumper=_EcDumper,
            default_flow_style=False, sort_keys=False, allow_unicode=True,
            explicit_start=True,
        )

    logger.info("Saved ecModel YAML to %s (%d reactions, %d ec-rxns)",
                out_path, len(model.reactions), ec_data.ec.n_rxns())


# ── load ───────────────────────────────────────────────────────────────────────

def load_ecmodel(path: str | Path) -> tuple[cobra.Model, ECModelData]:
    """Load a GECKO-compatible YAML file into ``(cobra.Model, ECModelData)``.

    Accepts both CNApy-written YAML (with a ``metaData`` block carrying
    ``schema_version``) and foreign GECKO YAML (in which case the loaded
    ``ECModelData`` is initialised with sensible defaults — callers should
    rebuild it if they need CNApy-specific state).

    Raises
    ------
    EcYamlError
        on missing sections, malformed structure, or I/O errors.
    """
    in_path = Path(path)
    if not in_path.exists():
        raise EcYamlError(f"YAML file not found: {in_path}")

    try:
        with in_path.open("r", encoding="utf-8") as fp:
            raw = yaml.load(fp, Loader=_EcLoader)
    except yaml.YAMLError as exc:
        raise EcYamlError(f"Failed to parse YAML: {exc}") from exc

    if not isinstance(raw, (dict, list)):
        raise EcYamlError(
            f"Unexpected YAML root type: {type(raw).__name__}. "
            "Expected top-level !!omap."
        )
    document = _to_dict(raw)

    if "metabolites" not in document or "reactions" not in document:
        raise EcYamlError(
            "Not a valid GECKO YAML: missing 'metabolites' or 'reactions'."
        )

    meta = document.get("metaData", {}) or {}

    # ── rebuild cobra.Model ────────────────────────────────────────────────
    model = cobra.Model(id_or_model=str(meta.get("id", "")) or "ecModel")
    if meta.get("name"):
        model.name = str(meta["name"])

    mets_by_id: dict[str, cobra.Metabolite] = {}
    for m in document.get("metabolites", []) or []:
        met = cobra.Metabolite(
            id=str(m.get("id", "")),
            name=str(m.get("name", "")) or None,
            compartment=str(m.get("compartment", "")) or None,
        )
        formula = m.get("formula")
        if formula:
            met.formula = str(formula)
        charge = m.get("charge")
        if charge is not None:
            try:
                met.charge = int(charge)
            except (TypeError, ValueError):
                pass
        annotation = m.get("annotation")
        if isinstance(annotation, dict):
            met.annotation.update(annotation)
        mets_by_id[met.id] = met
    if mets_by_id:
        model.add_metabolites(list(mets_by_id.values()))

    # Compartments: cobra derives from metabolites but allow explicit overrides.
    compartments = document.get("compartments")
    if isinstance(compartments, dict) and compartments:
        model.compartments = dict(compartments)

    # Reactions
    for r in document.get("reactions", []) or []:
        rxn = cobra.Reaction(
            id=str(r.get("id", "")),
            name=str(r.get("name", "")) or None,
            lower_bound=_as_float(r.get("lower_bound", 0.0)),
            upper_bound=_as_float(r.get("upper_bound", 1000.0)),
        )
        model.add_reactions([rxn])

        met_dict: dict[cobra.Metabolite, float] = {}
        mets_payload = r.get("metabolites") or {}
        if isinstance(mets_payload, dict):
            for met_id, coeff in mets_payload.items():
                met = mets_by_id.get(str(met_id))
                if met is not None:
                    met_dict[met] = _as_float(coeff)
        if met_dict:
            rxn.add_metabolites(met_dict)

        gpr = r.get("gene_reaction_rule") or ""
        if gpr:
            rxn.gene_reaction_rule = str(gpr)

        subsystem = r.get("subsystem")
        if subsystem:
            rxn.subsystem = str(subsystem) if not isinstance(
                subsystem, list) else "; ".join(str(s) for s in subsystem)

        annotation = r.get("annotation")
        if isinstance(annotation, dict):
            rxn.annotation.update(annotation)

    # ── rebuild ECModelData ────────────────────────────────────────────────
    ec_data = ECModelData()
    ec_data.is_ecmodel = True
    ec_data.schema_version = int(meta.get("schema_version",
                                           CURRENT_SCHEMA_VERSION))
    ec_data.gecko_light = str(meta.get("geckoLight", "false")).lower() == "true"
    ec_data.sigma = _as_float(meta.get("sigma", 0.5))
    ec_data.f = _as_float(meta.get("f", 0.5))
    ec_data.ptot = _as_float(meta.get("ptot", 0.5))

    ec_data.ec = EcStructure.from_yaml_omap(
        document.get("ec-rxns"),
        document.get("ec-enzymes"),
        gecko_light=ec_data.gecko_light,
    )

    # Rebuild uniprot_data from the ec-enzymes section so downstream CNApy
    # code (apply_proteomics, Usage page, etc.) has MW / gene available.
    ec_data.uniprot_data = {
        uid: {
            "gene": ec_data.ec.genes[j] if j < len(ec_data.ec.genes) else "",
            "mw_da": ec_data.ec.mw[j] if j < len(ec_data.ec.mw) else 0.0,
            "ec": "",
            "sequence": ec_data.ec.sequence[j]
                if j < len(ec_data.ec.sequence) else "",
        }
        for j, uid in enumerate(ec_data.ec.enzymes)
    }

    # Map enzyme metabolites / reactions if present in the cobra model —
    # gives the rest of CNApy something to look up via ec_data.enzyme_met_ids.
    ec_data.enzyme_met_ids = {
        uid: f"prot_{uid}"
        for uid in ec_data.ec.enzymes
        if f"prot_{uid}" in mets_by_id
    }
    ec_data.enzyme_rxn_ids = {
        uid: f"usage_prot_{uid}"
        for uid in ec_data.ec.enzymes
        if any(r.id == f"usage_prot_{uid}" for r in model.reactions)
    }

    logger.info("Loaded ecModel YAML from %s (%d reactions, %d ec-rxns)",
                in_path, len(model.reactions), ec_data.ec.n_rxns())

    return model, ec_data


