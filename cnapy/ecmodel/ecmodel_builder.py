"""
ecmodel_builder.py
==================
Translates GECKO 3.0 algorithm (makeEcModel / applyKcatConstraints /
setProtPoolSize) from MATLAB into COBRApy-compatible Python.

Public API
----------
parse_customkcats_file(path)  → list[dict]
parse_uniprot_file(path)      → dict
compute_coverage(model, ec_data) → dict
build_ecmodel(model, ec_data) → (cobra.Model, list, list)
revert_to_gem(ecmodel, ec_data) → cobra.Model
get_enzyme_usage(ecmodel, solution) → dict
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import cobra
import pandas as pd

if TYPE_CHECKING:
    from cnapy.ecmodel.ecmodel_data import ECModelData

# ── file parsers ──────────────────────────────────────────────────────────────

def parse_customkcats_file(path: str) -> list[dict]:
    """
    Parse a customKcats TSV/CSV/Excel file into a list of kcat entry dicts.

    Column names are case-insensitive. Missing optional columns default to
    empty values and do not cause errors.

    REQUIRED:
        kcat        — turnover number (1/s), must be > 0
        proteins    — UniProt ID(s) of the catalysing enzyme; must match
                      'Entry' in the UniProt file (for MW lookup).
                      Use '+' for complexes (e.g. P0A9P0+P12345).
        rxns        — reaction ID(s) in the model, comma-separated.
                      Must match the loaded model's reaction identifiers
                      (e.g. MAR03905 for Human-GEM, PFK for BiGG models).

    NOTE: if 'rxns' is omitted, reactions are found via proteins → gene
    name (from UniProt file) → model GPR. This only works when gene names
    in the UniProt file match the model's gene ID format. When they differ
    (e.g. gene symbols vs Ensembl IDs), provide 'rxns' explicitly.

    OPTIONAL:
        stoicho     — subunit copy number per protein ('+'-separated;
                      default 1)
        genes       — gene ID(s) matching the model's gene identifiers;
                      used as additional targets for reaction matching
        gene_name   — informational only, not used in calculations
        notes       — free-text comment; not used in calculations

    Returns list of dicts with keys:
        proteins (list[str]), genes (list[str]), gene_name (str),
        kcat (float), rxns (list[str]), notes (str), stoicho (list[int])
    """
    sep = "\t" if path.endswith(".tsv") else ","
    if path.endswith((".xlsx", ".xls")):
        df = pd.read_excel(path, header=0)
    else:
        df = pd.read_csv(path, sep=sep, header=0)

    df.columns = [c.strip().lower() for c in df.columns]

    entries = []
    for _, row in df.iterrows():
        proteins_raw = str(row.get("proteins", "")).strip()
        genes_raw    = str(row.get("genes", "")).strip()
        gene_name    = str(row.get("gene_name", "")).strip()
        kcat_raw     = row.get("kcat", None)
        rxns_raw     = str(row.get("rxns", "")).strip()
        notes        = str(row.get("notes", "")).strip()
        stoicho_raw  = str(row.get("stoicho", "1")).strip()

        # Skip blank or NaN kcat
        try:
            kcat = float(kcat_raw)
        except (TypeError, ValueError):
            continue
        if kcat <= 0:
            continue

        proteins = [p.strip() for p in proteins_raw.split("+") if p.strip() and p.strip() != "nan"]
        genes    = [g.strip() for g in re.split(r"\band\b|\+", genes_raw) if g.strip() and g.strip() != "nan"]
        rxns     = [r.strip() for r in rxns_raw.split(",") if r.strip() and r.strip() != "nan"]

        # stoicho – one integer per protein
        stoicho_parts = [s.strip() for s in stoicho_raw.split("+") if s.strip() and s.strip() != "nan"]
        stoicho = []
        for part in stoicho_parts:
            try:
                stoicho.append(int(float(part)))
            except ValueError:
                stoicho.append(1)
        # pad / truncate to match proteins length
        if proteins:
            while len(stoicho) < len(proteins):
                stoicho.append(1)
            stoicho = stoicho[: len(proteins)]
        else:
            stoicho = [1]

        entries.append({
            "proteins": proteins,
            "genes":    genes,
            "gene_name": gene_name,
            "kcat":     kcat,
            "rxns":     rxns,
            "notes":    notes,
            "stoicho":  stoicho,
        })

    return entries


def parse_uniprot_file(path: str) -> dict:
    """
    Parse a UniProt TSV/CSV/Excel file into a protein data dict.

    Column names are matched flexibly (case-insensitive, partial match).
    Missing optional columns default to empty values.

    REQUIRED:
        Entry                — UniProt accession ID (e.g. P0A9P0);
                               rows without a valid entry are skipped
        Mass                 — molecular weight in Daltons (e.g. 86376);
                               also recognised as "Molecular weight", "MW", etc.
                               Enzymes with MW = 0 cannot be constrained.

    OPTIONAL:
        Gene Names           — gene name for display and GPR-based reaction
                               matching (e.g. "thrA"); any column whose name
                               contains both "gene" and "name" is accepted.
                               NOTE: for reaction matching via kcat 'proteins',
                               these must use the same ID format as the model's
                               genes (e.g. gene symbols or Ensembl IDs).
        Sequence             — amino acid sequence; only needed for UniKP
                               kcat prediction; not used in GECKO calculations
        EC number            — enzyme classification number; stored but not
                               used in GECKO calculations

    Returns dict: { uniprot_id: { gene, mw_da, ec, sequence } }
    """
    sep = "\t" if path.endswith(".tsv") else ","
    if path.endswith((".xlsx", ".xls")):
        df = pd.read_excel(path, header=0)
    else:
        df = pd.read_csv(path, sep=sep, header=0)

    df.columns = [c.strip() for c in df.columns]

    # Flexible column name matching
    col_map = {}
    for c in df.columns:
        lc = c.lower()
        if lc == "entry":
            col_map["entry"] = c
        elif "gene" in lc and "name" in lc:
            col_map["gene"] = c
        elif "ec" in lc and "number" in lc:
            col_map["ec"] = c
        elif "mass" in lc or "molecular weight" in lc or lc in ("mw", "mw (da)", "mw(da)"):
            if "mass" not in col_map:   # first match wins
                col_map["mass"] = c
        elif "sequence" in lc:
            col_map["sequence"] = c

    result = {}
    for _, row in df.iterrows():
        uid = str(row.get(col_map.get("entry", "Entry"), "")).strip()
        if not uid or uid == "nan":
            continue
        gene = str(row.get(col_map.get("gene", ""), "")).strip()
        ec   = str(row.get(col_map.get("ec", ""), "")).strip()
        seq  = str(row.get(col_map.get("sequence", ""), "")).strip()
        try:
            mw = float(str(row.get(col_map.get("mass", ""), 0)).replace(",", ""))
        except (ValueError, TypeError):
            mw = 0.0

        result[uid] = {
            "gene":     gene if gene != "nan" else "",
            "mw_da":   mw,
            "ec":      ec if ec != "nan" else "",
            "sequence": seq if seq != "nan" else "",
        }

    return result


# ── coverage analysis ─────────────────────────────────────────────────────────

def _build_gene_to_uniprot(ec_data: ECModelData) -> dict:
    """Build { gene_id: uniprot_id } from uniprot_data and kcat_entries."""
    g2u = {}
    for uid, info in ec_data.uniprot_data.items():
        gene = info.get("gene", "")
        if gene:
            g2u[gene] = uid
    # also capture from kcat_entries (proteins vs genes fields)
    for entry in ec_data.kcat_entries:
        for prot, gene in zip(entry["proteins"], entry["genes"]):
            if gene and prot:
                g2u[gene] = prot
    return g2u


def _reactions_for_entry(model: cobra.Model, entry: dict,
                          gene_to_uniprot: dict) -> list[str]:
    """Return the model reaction IDs matched by a KcatEntry."""
    rxn_ids = []

    # Priority 1: explicit reaction IDs in entry["rxns"]
    if entry["rxns"]:
        for rid in entry["rxns"]:
            for rxn in model.reactions:
                if rxn.id == rid or rxn.id == rid + "_REV":
                    rxn_ids.append(rxn.id)
        return list(dict.fromkeys(rxn_ids))  # deduplicate, preserve order

    # Priority 2: match by proteins → genes → reactions
    if entry["proteins"]:
        # For each protein, find the gene, then find reactions with that gene
        target_genes = set()
        for prot in entry["proteins"]:
            # direct gene from kcat entry
            for g, u in gene_to_uniprot.items():
                if u == prot:
                    target_genes.add(g)
            # also try genes field of the same entry
        for g in entry["genes"]:
            target_genes.add(g)

        for rxn in model.reactions:
            rxn_gene_ids = {g.id for g in rxn.genes}
            if rxn_gene_ids & target_genes:
                rxn_ids.append(rxn.id)
        return rxn_ids

    return []


def compute_coverage(model: cobra.Model, ec_data: ECModelData) -> dict:
    """
    Compute kcat/MW coverage statistics.

    Returns dict:
        total_enzyme_rxns   – reactions with at least one gene
        covered_rxns        – reactions that have a matching kcat entry
        covered_rxn_ids     – list of those reaction IDs
        missing_rxn_ids     – reactions with genes but no kcat
        missing_mw_proteins – UniProt IDs referenced in kcat_entries but not in uniprot_data
    """
    gene_to_uniprot = _build_gene_to_uniprot(ec_data)

    enzyme_rxn_ids = {r.id for r in model.reactions if r.genes}
    covered = set()
    missing_mw_proteins = set()

    for entry in ec_data.kcat_entries:
        matched = _reactions_for_entry(model, entry, gene_to_uniprot)
        covered.update(matched)
        # check MW availability
        for prot in entry["proteins"]:
            if prot and prot not in ec_data.uniprot_data:
                missing_mw_proteins.add(prot)
            elif prot and ec_data.uniprot_data.get(prot, {}).get("mw_da", 0) <= 0:
                missing_mw_proteins.add(prot)

    missing_rxn_ids = sorted(enzyme_rxn_ids - covered)
    covered_rxn_ids = sorted(covered & enzyme_rxn_ids)

    return {
        "total_enzyme_rxns":   len(enzyme_rxn_ids),
        "covered_rxns":        len(covered_rxn_ids),
        "covered_rxn_ids":     covered_rxn_ids,
        "missing_rxn_ids":     missing_rxn_ids,
        "missing_mw_proteins": sorted(missing_mw_proteins),
    }


# ── ecModel construction ──────────────────────────────────────────────────────

def _split_reversible_reactions(ecmodel: cobra.Model, ec_data: ECModelData):
    """
    Split every reversible (lb < 0, ub > 0) enzyme-catalyzed reaction into
    a forward copy (original id, lb=0) and a reverse copy (id_REV, lb=0,
    negated stoichiometry).

    Populates ec_data.split_rxn_map.
    """
    to_split = [
        r for r in list(ecmodel.reactions)
        if r.lower_bound < 0 and r.upper_bound > 0 and r.genes
    ]

    for rxn in to_split:
        rev_id = rxn.id + "_REV"
        rev_rxn = cobra.Reaction(
            id=rev_id,
            name=rxn.name + " (reverse)",
            lower_bound=0.0,
            upper_bound=-rxn.lower_bound,
        )
        rev_rxn.gene_reaction_rule = rxn.gene_reaction_rule
        # negate stoichiometry for reverse direction
        rev_stoich = {met: -coeff for met, coeff in rxn.metabolites.items()}
        rev_rxn.add_metabolites(rev_stoich)
        ecmodel.add_reactions([rev_rxn])

        # restrict original to forward only
        rxn.lower_bound = 0.0

        ec_data.split_rxn_map[rxn.id] = [rxn.id, rev_id]


def _add_protein_pool(ecmodel: cobra.Model, ec_data: ECModelData):
    """Add prot_pool pseudo-metabolite and exchange reaction."""
    prot_pool_met = cobra.Metabolite(
        id=ec_data.protein_pool_met_id,
        name="Protein pool",
        compartment="c",
    )
    ecmodel.add_metabolites([prot_pool_met])

    pool_rxn = cobra.Reaction(
        id=ec_data.protein_pool_rxn_id,
        name="Protein pool exchange",
        lower_bound=0.0,
        upper_bound=ec_data.pool_bound(),
    )
    pool_rxn.add_metabolites({prot_pool_met: 1.0})
    ecmodel.add_reactions([pool_rxn])


def _ensure_enzyme_met_and_usage(ecmodel: cobra.Model, ec_data: ECModelData,
                                  uniprot_id: str):
    """
    Add prot_{uniprot_id} pseudo-metabolite and usage_prot_{uniprot_id}
    reaction (if not already present).  Only used for the full ecModel.
    """
    met_id = f"prot_{uniprot_id}"
    rxn_id = f"usage_prot_{uniprot_id}"

    if met_id in ec_data.enzyme_met_ids:
        return  # already added

    prot_pool_met = ecmodel.metabolites.get_by_id(ec_data.protein_pool_met_id)

    enz_met = cobra.Metabolite(
        id=met_id,
        name=f"Enzyme {uniprot_id}",
        compartment="c",
    )
    ecmodel.add_metabolites([enz_met])

    usage_rxn = cobra.Reaction(
        id=rxn_id,
        name=f"Enzyme usage {uniprot_id}",
        lower_bound=0.0,
        upper_bound=1000.0,
    )
    usage_rxn.add_metabolites({
        prot_pool_met: -1.0,
        enz_met:       1.0,
    })
    ecmodel.add_reactions([usage_rxn])

    ec_data.enzyme_met_ids[uniprot_id] = met_id
    ec_data.enzyme_rxn_ids[uniprot_id] = rxn_id


def _stoich_coeff(mw_da: float, kcat_per_s: float, subunits: int = 1) -> float:
    """
    Stoichiometric coefficient for the enzyme pseudo-metabolite.
    Follows GECKO formula: -(MW_Da / (kcat_1/s * 3600)) * subunit_copies
    Units are internally consistent with GECKO's protein pool in mg/gDCW.
    """
    if kcat_per_s <= 0 or mw_da <= 0:
        return 0.0
    return -(mw_da / (kcat_per_s * 3600.0)) * subunits


def _apply_entry_full(ecmodel: cobra.Model, ec_data: ECModelData,
                       entry: dict, gene_to_uniprot: dict) -> list[str]:
    """
    Apply one kcat entry to matching reactions (full ecModel).
    Returns list of reaction IDs that could not be constrained.
    """
    proteins = entry["proteins"]
    kcat     = entry["kcat"]
    stoicho  = entry["stoicho"]
    failed   = []

    matched_rxn_ids = _reactions_for_entry(ecmodel, entry, gene_to_uniprot)
    if not matched_rxn_ids:
        return failed  # nothing to apply

    # Build per-protein metabolites and compute total stoich coeff
    # For a complex: coeff = sum_i(stoicho_i * MW_i) / kcat (one enzyme unit)
    if proteins:
        total_coeff = 0.0
        valid = True
        for prot, units in zip(proteins, stoicho):
            info = ec_data.uniprot_data.get(prot, {})
            mw = info.get("mw_da", 0.0)
            if mw <= 0:
                valid = False
                break
            _ensure_enzyme_met_and_usage(ecmodel, ec_data, prot)
            coeff = _stoich_coeff(mw, kcat, units)
            # We'll store per-protein coefficients on individual reactions below
            total_coeff += abs(coeff)

        if not valid:
            return matched_rxn_ids  # can't constrain without MW

        for rxn_id in matched_rxn_ids:
            try:
                rxn = ecmodel.reactions.get_by_id(rxn_id)
            except KeyError:
                failed.append(rxn_id)
                continue

            for prot, units in zip(proteins, stoicho):
                info = ec_data.uniprot_data[prot]
                mw = info["mw_da"]
                coeff = _stoich_coeff(mw, kcat, units)
                enz_met = ecmodel.metabolites.get_by_id(ec_data.enzyme_met_ids[prot])
                # Use minimum-cost (= highest kcat) entry if this enzyme-reaction
                # pair appears in multiple kcat rows.
                existing = rxn.metabolites.get(enz_met, 0.0)
                if existing == 0.0 or abs(coeff) < abs(existing):
                    rxn.add_metabolites({enz_met: coeff}, combine=False)

    else:
        # No proteins specified – only rxns.  We need the gene → protein mapping.
        for rxn_id in matched_rxn_ids:
            try:
                rxn = ecmodel.reactions.get_by_id(rxn_id)
            except KeyError:
                failed.append(rxn_id)
                continue

            # Find enzymes catalysing this reaction via GPR
            for gene in rxn.genes:
                prot = gene_to_uniprot.get(gene.id)
                if prot is None:
                    continue
                info = ec_data.uniprot_data.get(prot, {})
                mw = info.get("mw_da", 0.0)
                if mw <= 0:
                    continue
                _ensure_enzyme_met_and_usage(ecmodel, ec_data, prot)
                coeff = _stoich_coeff(mw, kcat, 1)
                enz_met = ecmodel.metabolites.get_by_id(ec_data.enzyme_met_ids[prot])
                rxn.add_metabolites({enz_met: coeff}, combine=False)

    return failed


def _apply_entry_light(ecmodel: cobra.Model, ec_data: ECModelData,
                        entry: dict, gene_to_uniprot: dict) -> list[str]:
    """
    Apply one kcat entry to matching reactions (GECKO-light ecModel).
    Protein cost goes directly to prot_pool with the minimum MW/kcat coefficient.
    """
    proteins = entry["proteins"]
    kcat     = entry["kcat"]
    stoicho  = entry["stoicho"]
    failed   = []

    matched_rxn_ids = _reactions_for_entry(ecmodel, entry, gene_to_uniprot)
    if not matched_rxn_ids:
        return failed

    prot_pool_met = ecmodel.metabolites.get_by_id(ec_data.protein_pool_met_id)

    if proteins:
        # Use minimum protein cost across all proteins in the entry
        min_coeff = None
        for prot, units in zip(proteins, stoicho):
            info = ec_data.uniprot_data.get(prot, {})
            mw = info.get("mw_da", 0.0)
            if mw <= 0:
                continue
            c = _stoich_coeff(mw, kcat, units)
            if min_coeff is None or abs(c) < abs(min_coeff):
                min_coeff = c

        if min_coeff is None:
            return matched_rxn_ids

        for rxn_id in matched_rxn_ids:
            try:
                rxn = ecmodel.reactions.get_by_id(rxn_id)
            except KeyError:
                failed.append(rxn_id)
                continue
            # Use the minimum-cost (= highest kcat) isozyme for each reaction.
            # If a previous entry already set a coefficient, only overwrite if
            # this entry's cost is strictly lower (i.e. this enzyme is faster).
            existing = rxn.metabolites.get(prot_pool_met, 0.0)
            if existing == 0.0 or abs(min_coeff) < abs(existing):
                rxn.add_metabolites({prot_pool_met: min_coeff}, combine=False)

    return failed


def build_ecmodel(model: cobra.Model,
                  ec_data: ECModelData) -> tuple[cobra.Model, list, list]:
    """
    Convert a conventional GEM into an enzyme-constrained model.

    Implements the GECKO 3.0 algorithm (makeEcModel + applyKcatConstraints +
    setProtPoolSize) in COBRApy.

    Parameters
    ----------
    model    : original cobra.Model (not modified)
    ec_data  : ECModelData populated with kcat_entries, uniprot_data, sigma, f, ptot

    Returns
    -------
    ecmodel          : new cobra.Model with enzyme constraints
    unconstrained    : list of reaction IDs with GPR but no kcat applied
    missing_mw_prots : list of UniProt IDs whose MW is missing
    """
    ec_data.reset()
    ec_data.original_reaction_ids = [r.id for r in model.reactions]

    ecmodel = model.copy()

    gene_to_uniprot = _build_gene_to_uniprot(ec_data)

    # ── Stage 1: Split reversible enzyme-catalysed reactions ──────────────────
    if not ec_data.gecko_light:
        _split_reversible_reactions(ecmodel, ec_data)

    # ── Stage 1 cont.: Add protein pool ──────────────────────────────────────
    _add_protein_pool(ecmodel, ec_data)

    # ── Stage 2: Apply kcat constraints ──────────────────────────────────────
    unconstrained_set: set[str] = set()
    missing_mw_prots: set[str] = set()

    for entry in ec_data.kcat_entries:
        # check MW availability
        for prot in entry["proteins"]:
            if prot and (prot not in ec_data.uniprot_data
                         or ec_data.uniprot_data[prot].get("mw_da", 0) <= 0):
                missing_mw_prots.add(prot)

        if ec_data.gecko_light:
            failed = _apply_entry_light(ecmodel, ec_data, entry, gene_to_uniprot)
        else:
            failed = _apply_entry_full(ecmodel, ec_data, entry, gene_to_uniprot)

        unconstrained_set.update(failed)

    # ── Stage 3: Set protein pool bound ──────────────────────────────────────
    pool_rxn = ecmodel.reactions.get_by_id(ec_data.protein_pool_rxn_id)
    pool_rxn.upper_bound = ec_data.pool_bound()

    ec_data.is_ecmodel = True

    # Reactions with GPR that have no kcat constraint at all
    enzyme_rxn_ids = {r.id for r in ecmodel.reactions if r.genes
                      and not r.id.startswith(("usage_prot_", "prot_pool"))}
    # Check which have no enzyme metabolite added
    constrained = set()
    for rxn in ecmodel.reactions:
        if any(m.id.startswith("prot_") and m.id != "prot_pool"
               for m in rxn.metabolites):
            constrained.add(rxn.id)
    if ec_data.gecko_light:
        ppm = ec_data.protein_pool_met_id
        for rxn in ecmodel.reactions:
            if ppm in {m.id for m in rxn.metabolites}:
                constrained.add(rxn.id)

    unconstrained = sorted(enzyme_rxn_ids - constrained)

    return ecmodel, unconstrained, sorted(missing_mw_prots)


# ── revert ────────────────────────────────────────────────────────────────────

def revert_to_gem(ecmodel: cobra.Model, ec_data: ECModelData) -> cobra.Model:
    """
    Remove all GECKO additions and return a model equivalent to the original GEM.

    Removes:
    - prot_pool metabolite + exchange reaction
    - prot_{id} metabolites + usage_prot_{id} reactions (full mode)
    - enzyme stoichiometric coefficients from metabolic reactions
    - _REV split reactions (restoring original reversibility)
    """
    gem = ecmodel.copy()

    # ── remove enzyme pseudo-metabolite stoich from metabolic reactions ───────
    prot_mets = [m for m in gem.metabolites if m.id.startswith("prot_")]
    for met in prot_mets:
        rxns_to_clean = list(met.reactions)
        for rxn in rxns_to_clean:
            rxn.add_metabolites({met: 0}, combine=False)

    # ── remove usage reactions ────────────────────────────────────────────────
    usage_rxns = [r for r in gem.reactions if r.id.startswith("usage_prot_")]
    gem.remove_reactions(usage_rxns, remove_orphans=True)

    # ── remove pool exchange ──────────────────────────────────────────────────
    try:
        pool_rxn = gem.reactions.get_by_id(ec_data.protein_pool_rxn_id)
        gem.remove_reactions([pool_rxn], remove_orphans=True)
    except KeyError:
        pass

    # ── remove _REV split reactions and restore original bounds ──────────────
    rev_rxns = [r for r in gem.reactions if r.id.endswith("_REV")]
    for rev_rxn in rev_rxns:
        base_id = rev_rxn.id[:-4]  # strip _REV
        try:
            base_rxn = gem.reactions.get_by_id(base_id)
            # restore original reversibility
            base_rxn.lower_bound = -rev_rxn.upper_bound
        except KeyError:
            pass
    gem.remove_reactions(rev_rxns, remove_orphans=True)

    ec_data.reset()
    return gem


# ── enzyme usage analysis ─────────────────────────────────────────────────────

def get_enzyme_usage(ecmodel: cobra.Model,
                     solution: cobra.Solution) -> dict:
    """
    Calculate enzyme usage from an FBA/pFBA solution on a full ecModel.

    Returns dict: { uniprot_id: { abs_usage, cap_usage, ub } }
    abs_usage – protein used (model units, proportional to mg/gDCW)
    cap_usage – fraction of available capacity used (0–1)
    ub        – upper bound of the usage reaction
    """
    result = {}
    for rxn in ecmodel.reactions:
        if not rxn.id.startswith("usage_prot_"):
            continue
        uid = rxn.id[len("usage_prot_"):]
        flux = solution.fluxes.get(rxn.id, 0.0)
        ub = rxn.upper_bound
        abs_usage = abs(flux)
        cap_usage = (abs_usage / ub) if ub > 0 else 0.0
        result[uid] = {
            "abs_usage": abs_usage,
            "cap_usage": cap_usage,
            "ub": ub,
        }
    return result


# ── proteomics integration ────────────────────────────────────────────────────

def parse_proteomics_file(path: str) -> dict:
    """
    Parse a proteomics TSV/CSV file.

    Expected format (tab or comma separated, with or without header):
        uniprot_id  |  level (mg/gDCW)

    Returns dict: { uniprot_id: level_mg_per_gDCW }
    """
    sep = "\t" if path.endswith(".tsv") else ","
    if path.endswith((".xlsx", ".xls")):
        df = pd.read_excel(path, header=0)
    else:
        # Try with header first, fall back to no-header
        df = pd.read_csv(path, sep=sep, header=0)

    # If first column looks like a header word (not a UniProt ID), keep header
    # UniProt IDs are typically 6-10 alphanumeric chars
    cols = list(df.columns)
    if len(cols) >= 2:
        try:
            float(cols[1])
            # Second column name is a number → no header, re-read
            df = pd.read_csv(path, sep=sep, header=None)
        except (ValueError, TypeError):
            pass

    result = {}
    for _, row in df.iterrows():
        uid = str(row.iloc[0]).strip()
        if not uid or uid.lower() in ("nan", "uniprot", "entry", "id"):
            continue
        try:
            level = float(str(row.iloc[1]).replace(",", ""))
            if level > 0:
                result[uid] = level
        except (ValueError, TypeError, IndexError):
            continue
    return result


def apply_proteomics(ecmodel: cobra.Model, ec_data: "ECModelData",
                     prot_data: dict) -> tuple[int, list]:
    """
    Apply measured protein abundances to constrain individual enzyme usage
    reactions (GECKO constrainEnzConcs equivalent).

    For each enzyme with a measured level (mg/gDCW):
      1. Disconnect from prot_pool (stoich coefficient → 0)
      2. Set usage_prot_{id}.upper_bound = measured_level [mg/gDCW]
      3. Reduce prot_pool UB by the sum of applied abundances, so the pool
         budget only covers unmeasured enzymes (GECKO paper §2.3).

    For enzymes without measured data: unchanged (still draw from prot_pool).

    Parameters
    ----------
    ecmodel   : ecModel (modified in-place)
    ec_data   : ECModelData
    prot_data : { uniprot_id: level_mg_per_gDCW }

    Returns
    -------
    n_applied : number of enzymes successfully constrained
    missing   : UniProt IDs in prot_data that have no usage reaction
    """
    prot_pool_met_id = ec_data.protein_pool_met_id
    n_applied = 0
    missing = []
    total_applied_mg = 0.0

    try:
        prot_pool_met = ecmodel.metabolites.get_by_id(prot_pool_met_id)
    except KeyError:
        raise ValueError("prot_pool metabolite not found – is this an ecModel?")

    try:
        pool_rxn = ecmodel.reactions.get_by_id(ec_data.protein_pool_rxn_id)
    except KeyError:
        raise ValueError("prot_pool_exchange reaction not found – is this an ecModel?")

    for uid, level in prot_data.items():
        usage_rxn_id = f"usage_prot_{uid}"
        try:
            usage_rxn = ecmodel.reactions.get_by_id(usage_rxn_id)
        except KeyError:
            missing.append(uid)
            continue

        # 1. Disconnect from protein pool
        usage_rxn.add_metabolites({prot_pool_met: 0}, combine=False)
        # 2. Constrain by measured abundance [mg/gDCW]
        usage_rxn.upper_bound = level
        total_applied_mg += level
        n_applied += 1

    # 3. Reduce pool UB so it covers only the remaining (unmeasured) enzymes
    new_pool_ub = max(0.0, pool_rxn.upper_bound - total_applied_mg)
    pool_rxn.upper_bound = new_pool_ub

    return n_applied, missing


def remove_proteomics(ecmodel: cobra.Model, ec_data: "ECModelData") -> int:
    """
    Remove proteomics constraints: reconnect all usage reactions to prot_pool,
    restore individual upper bounds to 1000, and restore the pool UB to its
    original value (ec_data.pool_bound()).

    Returns number of reactions restored.
    """
    prot_pool_met_id = ec_data.protein_pool_met_id
    try:
        prot_pool_met = ecmodel.metabolites.get_by_id(prot_pool_met_id)
    except KeyError:
        return 0

    n = 0
    for rxn in ecmodel.reactions:
        if not rxn.id.startswith("usage_prot_"):
            continue
        # Only restore if this reaction was disconnected from pool
        if prot_pool_met not in rxn.metabolites or rxn.metabolites[prot_pool_met] == 0:
            rxn.add_metabolites({prot_pool_met: -1.0}, combine=False)
            rxn.upper_bound = 1000.0
            n += 1

    # Restore pool UB to the value set at ecModel build time
    try:
        pool_rxn = ecmodel.reactions.get_by_id(ec_data.protein_pool_rxn_id)
        pool_rxn.upper_bound = ec_data.pool_bound()
    except KeyError:
        pass

    return n


def flexibilize_enz_concs(ecmodel: cobra.Model, prot_data: dict) -> list[str]:
    """
    Relax the minimum set of proteomics constraints needed to restore FBA
    feasibility (GECKO flexibilizeEnzConcs equivalent).

    Strategy
    --------
    1. Relax ALL proteomics constraints to 1000 (unconstrained).
    2. Verify FBA is feasible in this state (baseline check — if not, the
       model itself is broken and we cannot help).
    3. Re-tighten each constraint one by one (highest abundance first, as
       highly abundant enzymes are less likely to be the bottleneck).
       - If tightening keeps FBA feasible → keep it tight (enzyme stays
         constrained).
       - If tightening causes infeasibility → leave it relaxed.

    This correctly handles cases where *multiple* enzymes simultaneously
    contribute to infeasibility (the previous one-at-a-time approach failed
    in those cases).

    Parameters
    ----------
    ecmodel   : ecModel (modified in-place)
    prot_data : { uniprot_id: level_mg_per_gDCW } — the applied proteomics data

    Returns
    -------
    relaxed : list of UniProt IDs whose constraints were relaxed
    """
    # Collect enzymes that are currently proteomics-constrained
    constrained: list[tuple[str, float, cobra.Reaction]] = []
    for uid, level in prot_data.items():
        try:
            usage_rxn = ecmodel.reactions.get_by_id(f"usage_prot_{uid}")
        except KeyError:
            continue
        constrained.append((uid, level, usage_rxn))

    if not constrained:
        return []

    # Step 1: relax all proteomics constraints
    for uid, level, usage_rxn in constrained:
        usage_rxn.upper_bound = 1000.0

    # Step 2: baseline feasibility check (model without any proteomics UBs)
    sol = run_fba_on_ecmodel(ecmodel)
    if sol.status != "optimal":
        # Model is infeasible even without proteomics — restore and give up
        for uid, level, usage_rxn in constrained:
            usage_rxn.upper_bound = level
        return []

    # Step 3: re-tighten one by one, highest abundance first
    # (most abundant enzymes are least likely to be the bottleneck)
    relaxed_uids = {uid for uid, _, _ in constrained}  # start: all relaxed
    for uid, level, usage_rxn in sorted(constrained, key=lambda x: x[1], reverse=True):
        usage_rxn.upper_bound = level          # tighten
        sol = run_fba_on_ecmodel(ecmodel)
        if sol.status == "optimal":
            relaxed_uids.discard(uid)          # tightening is OK → keep tight
        else:
            usage_rxn.upper_bound = 1000.0     # tightening broke FBA → stay relaxed

    return sorted(relaxed_uids)


# ── kcat what-if analysis ─────────────────────────────────────────────────────

def get_kcat_entries_for_enzyme(ec_data: "ECModelData",
                                uniprot_id: str) -> list[int]:
    """Return indices into ec_data.kcat_entries that involve uniprot_id."""
    indices = []
    for i, entry in enumerate(ec_data.kcat_entries):
        if uniprot_id in entry["proteins"]:
            indices.append(i)
    return indices


def apply_kcat_multiplier(ecmodel: cobra.Model, ec_data: "ECModelData",
                          uniprot_id: str, multiplier: float) -> int:
    """
    Scale the kcat of all reactions catalysed by uniprot_id by `multiplier`.

    Modifies stoichiometric coefficients of prot_{uniprot_id} in all
    metabolic reactions in-place.

    New coeff = old_coeff / multiplier  (less protein needed per flux unit)

    Returns number of reactions updated.
    """
    if multiplier <= 0:
        raise ValueError("multiplier must be positive")

    met_id = f"prot_{uniprot_id}"
    try:
        prot_met = ecmodel.metabolites.get_by_id(met_id)
    except KeyError:
        raise ValueError(f"No enzyme metabolite found for {uniprot_id}. "
                         "Check that this enzyme is in the ecModel.")

    n = 0
    for rxn in list(prot_met.reactions):
        # Skip usage reaction itself
        if rxn.id.startswith("usage_prot_"):
            continue
        old_coeff = rxn.metabolites[prot_met]
        new_coeff = old_coeff / multiplier
        rxn.add_metabolites({prot_met: new_coeff}, combine=False)
        n += 1
    return n


def run_fba_on_ecmodel(ecmodel: cobra.Model) -> cobra.Solution:
    """Run FBA and return the solution."""
    with ecmodel:
        return ecmodel.optimize()
