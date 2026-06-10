"""Stage 1 expansion utilities for GECKO 3.0 ecModel construction.

Design Ref: §2.1 expansion.py — houses the transformations that are new or
changed vs. the legacy CNApy implementation:

* FR-01 clear_pseudoreaction_gpr    — Box 1 Step 1
* FR-04 expand_isozymes              — Box 1 Step 5
* FR-07/08 add_protein_pool, add_usage_reactions (STUBS — implemented in p1-signs)

Stage 1 Box 1 reference: Chen et al. Nat Protoc 19:629–667 (2024), p.642.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING

import cobra

from cnapy.ecmodel.exceptions import GprParseError

if TYPE_CHECKING:
    from cnapy.ecmodel.ecmodel_data import ECModelData

logger = logging.getLogger("cnapy.ecmodel.expansion")

# Design Ref: §10.1 — GECKO uses "pseudoreaction" keyword in rxnNames
PSEUDOREACTION_KEYWORD = "pseudoreaction"


# ── FR-01: pseudoreaction GPR clearance ──────────────────────────────────────

def clear_pseudoreaction_gpr(model: cobra.Model,
                              extra_tsv: str | Path | None = None) -> list[str]:
    """Clear gene associations from pseudoreactions.

    GECKO 3.0 Box 1 Step 1 (Nat Protoc 2024, p.642): pseudoreactions such as
    biomass assembly or ATP maintenance should not be catalysed by enzymes and
    therefore must not draw from the protein pool. This function sets their
    gene_reaction_rule to an empty string.

    A reaction is treated as pseudoreaction when *any* of the following holds:
      * ``rxn.name`` contains the substring "pseudoreaction" (case-insensitive)
      * ``rxn.id`` is listed in ``extra_tsv`` first column (optional user file)

    Parameters
    ----------
    model : cobra.Model
        The metabolic model. Modified in place.
    extra_tsv : str or Path, optional
        TSV file whose first column contains reaction IDs to treat as
        pseudoreactions, regardless of the name heuristic.

    Returns
    -------
    cleared : list of reaction IDs whose GPR was cleared.
    """
    extra_ids: set[str] = set()
    if extra_tsv is not None:
        p = Path(extra_tsv)
        if p.exists():
            with p.open() as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    rxn_id = line.split("\t", 1)[0].strip()
                    if rxn_id:
                        extra_ids.add(rxn_id)

    cleared: list[str] = []
    for rxn in model.reactions:
        is_pseudo = PSEUDOREACTION_KEYWORD in (rxn.name or "").lower()
        is_extra = rxn.id in extra_ids
        if (is_pseudo or is_extra) and rxn.gene_reaction_rule:
            rxn.gene_reaction_rule = ""
            cleared.append(rxn.id)

    if cleared:
        logger.info("Cleared GPR from %d pseudoreaction(s): %s",
                    len(cleared), ", ".join(cleared[:10])
                    + (" …" if len(cleared) > 10 else ""))
    return cleared


# ── FR-02: invert backward-only irreversible reactions ───────────────────────

def invert_backward_only(model: cobra.Model) -> list[str]:
    """Flip stoichiometry + bounds of reactions that only run backwards.

    GECKO 3.0 Box 1 Step 2 (Nat Protoc 2024, p.642). A reaction with
    ``lb < 0 and ub <= 0`` is "backward-only" — only the reverse
    direction is allowed. GECKO normalises these so every irreversible
    reaction runs in the forward direction:

    * Multiply every stoichiometric coefficient by ``-1``.
    * Set new bounds to ``(0, -old_lb)``.

    Reactions that are already forward (``lb >= 0``) or reversible
    (``lb < 0 < ub``) are untouched.

    Parameters
    ----------
    model : cobra.Model (mutated in place)

    Returns
    -------
    inverted : list of reaction IDs that were flipped.
    """
    inverted: list[str] = []
    for rxn in model.reactions:
        if rxn.lower_bound < 0 and rxn.upper_bound <= 0:
            # Negate every coefficient in one batched call (cobra emits a
            # solver-side update per add_metabolites, so batching matters
            # on large models).
            rxn.add_metabolites(
                {met: -2 * coeff for met, coeff in rxn.metabolites.items()},
                combine=True,
            )
            old_lb = rxn.lower_bound
            rxn.lower_bound = 0.0
            rxn.upper_bound = -old_lb
            inverted.append(rxn.id)

    if inverted:
        logger.info("Inverted %d backward-only reaction(s)", len(inverted))
    return inverted


# ── FR-04: isozyme split (expandModel equivalent) ────────────────────────────

_AND_PATTERN = re.compile(r"\s+and\s+", re.IGNORECASE)
_OR_PATTERN = re.compile(r"\s+or\s+", re.IGNORECASE)


def parse_gpr_to_isozyme_sets(gpr: str) -> list[list[str]]:
    """Parse a gene_reaction_rule string into isozyme gene sets.

    GECKO 3.0 Box 1 Step 5 helper — an ``OR`` at the top level separates
    alternative isozymes (or complexes), while ``AND`` inside each isozyme
    lists the proteins that form a complex.

    Examples
    --------
    >>> parse_gpr_to_isozyme_sets("(g1 and g2) or g3 or (g4 and g5 and g6)")
    [['g1', 'g2'], ['g3'], ['g4', 'g5', 'g6']]
    >>> parse_gpr_to_isozyme_sets("g1")
    [['g1']]
    >>> parse_gpr_to_isozyme_sets("")
    []

    Parameters
    ----------
    gpr : str
        cobra-style gene_reaction_rule. May be empty.

    Returns
    -------
    list[list[str]]
        Outer list = isozymes (OR). Inner list = genes (AND). Empty if gpr blank.

    Raises
    ------
    GprParseError
        If the gpr contains unbalanced parentheses.
    """
    gpr = (gpr or "").strip()
    if not gpr:
        return []
    gpr = re.sub(r"\s+", " ", gpr)

    # Validate parenthesis balance up-front.
    depth = 0
    for ch in gpr:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                raise GprParseError(gpr, "unbalanced parenthesis")
    if depth != 0:
        raise GprParseError(gpr, "unbalanced parenthesis")

    # Recursively distribute to disjunctive normal form (DNF) so that nested
    # OR groups inside an AND complex — e.g. "(g1 or g2) and g3" — are fully
    # expanded into separate isozyme sets [['g1', 'g3'], ['g2', 'g3']] rather
    # than left as a single set containing the literal token "(g1 or g2)".
    isozymes = _expand_to_dnf(gpr)

    result: list[list[str]] = []
    for gene_set in isozymes:
        genes = [g for g in gene_set if g]
        if genes:
            result.append(genes)
    return result


def _expand_to_dnf(text: str) -> list[list[str]]:
    """Recursively expand a GPR fragment to disjunctive normal form.

    Returns a list of isozyme sets (each an AND-list of gene tokens). An
    ``OR`` produces the union of its operands' alternatives; an ``AND``
    produces the Cartesian product of its operands' alternatives so that any
    nested OR group is distributed over the surrounding complex.
    """
    text = _strip_outer_parens(text.strip())
    if not text:
        return []

    or_parts = _split_top_level(text, _OR_PATTERN)
    if len(or_parts) > 1:
        alternatives: list[list[str]] = []
        for part in or_parts:
            alternatives.extend(_expand_to_dnf(part))
        return alternatives

    and_parts = _split_top_level(text, _AND_PATTERN)
    if len(and_parts) > 1:
        # Cartesian product of each AND operand's alternatives.
        product: list[list[str]] = [[]]
        for part in and_parts:
            sub_alts = _expand_to_dnf(part)
            if not sub_alts:
                continue
            product = [combo + alt for combo in product for alt in sub_alts]
        return [combo for combo in product if combo]

    # Single leaf token (a gene id), possibly wrapped in redundant parens.
    leaf = _strip_outer_parens(text.strip()).strip()
    return [[leaf]] if leaf else []


def _split_top_level(text: str, sep_pattern: re.Pattern) -> list[str]:
    """Split ``text`` by ``sep_pattern`` only at paren-depth 0."""
    parts: list[str] = []
    depth = 0
    last = 0
    i = 0
    while i < len(text):
        ch = text[i]
        if ch == "(":
            depth += 1
            i += 1
            continue
        if ch == ")":
            depth -= 1
            i += 1
            continue
        if depth == 0:
            m = sep_pattern.match(text, i)
            if m:
                parts.append(text[last:i])
                last = m.end()
                i = m.end()
                continue
        i += 1
    parts.append(text[last:])
    return parts


def _strip_outer_parens(part: str) -> str:
    """Remove a matched pair of outermost parentheses, if any."""
    while part.startswith("(") and part.endswith(")"):
        depth = 0
        matched = True
        for i, ch in enumerate(part):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0 and i != len(part) - 1:
                    matched = False
                    break
        if matched:
            part = part[1:-1].strip()
        else:
            break
    return part


def expand_isozymes(model: cobra.Model,
                    ec_data: "ECModelData | None" = None
                    ) -> dict[str, list[str]]:
    """Split isozyme-catalysed reactions into ``_EXP_N`` variants.

    GECKO 3.0 Box 1 Step 5 (Nat Protoc 2024, p.642): reactions whose GPR
    contains ``OR`` are duplicated so that each resulting reaction is catalysed
    by exactly one isozyme (or complex). This ensures that later kcat
    constraints treat isozymes as alternatives (OR), not as co-requirements
    (AND).

    Skipped for ``gecko_light`` models — the light formulation collapses
    isozymes to the minimum-cost one at kcat-application time.

    Parameters
    ----------
    model : cobra.Model
        ecModel being built. Modified in place.
    ec_data : ECModelData, optional
        If given, populates ``ec_data.isozyme_split_map``.

    Returns
    -------
    split_map : ``{original_rxn_id: [exp_rxn_ids]}``
        Mapping used by ``revert_to_gem`` to reverse the operation.
    """
    if ec_data is not None and getattr(ec_data, "gecko_light", False):
        logger.debug("expand_isozymes skipped for GECKO-light model")
        return {}

    split_map: dict[str, list[str]] = {}

    # Snapshot the reaction list: we will add new and remove existing ones.
    candidates = [r for r in list(model.reactions) if r.gene_reaction_rule]

    for rxn in candidates:
        try:
            isozyme_sets = parse_gpr_to_isozyme_sets(rxn.gene_reaction_rule)
        except GprParseError as exc:
            logger.warning("%s: GPR parse failed, skipping isozyme split (%s)",
                           rxn.id, exc.reason)
            continue

        if len(isozyme_sets) <= 1:
            continue  # Single isozyme or complex — no split needed.

        exp_ids: list[str] = []
        for idx, gene_set in enumerate(isozyme_sets, start=1):
            new_id = f"{rxn.id}_EXP_{idx}"
            if new_id in model.reactions:
                logger.warning("%s: %s already exists, skipping entire split",
                               rxn.id, new_id)
                exp_ids = []
                break

            new_rxn = cobra.Reaction(
                id=new_id,
                name=rxn.name,
                subsystem=rxn.subsystem or "",
                lower_bound=rxn.lower_bound,
                upper_bound=rxn.upper_bound,
            )
            new_rxn.add_metabolites({m: c for m, c in rxn.metabolites.items()})
            # Copy annotation and notes so downstream tools (SBML export etc.)
            # still see the same metadata on each variant.
            try:
                new_rxn.annotation = dict(rxn.annotation)
            except Exception:  # noqa: BLE001
                pass
            try:
                new_rxn.notes = dict(rxn.notes)
            except Exception:  # noqa: BLE001
                pass
            model.add_reactions([new_rxn])
            # gene_reaction_rule must be set *after* add_reactions so that
            # model.genes gets the new Gene objects registered.
            new_rxn.gene_reaction_rule = " and ".join(gene_set)
            exp_ids.append(new_id)

        if not exp_ids:
            continue

        split_map[rxn.id] = exp_ids
        # Remove the original — its GPR has been distributed to _EXP variants.
        # Leave orphan genes alone; other reactions may still reference them.
        model.remove_reactions([rxn], remove_orphans=False)

    if ec_data is not None:
        # Store on ec_data so revert_to_gem can reverse.
        existing = getattr(ec_data, "isozyme_split_map", None)
        if existing is None:
            # Field doesn't exist yet (pre-schema-v2) — attach as attribute.
            ec_data.isozyme_split_map = {}  # type: ignore[attr-defined]
            existing = ec_data.isozyme_split_map  # type: ignore[attr-defined]
        existing.update(split_map)

    if split_map:
        logger.info("Expanded %d reaction(s) into isozyme variants",
                    len(split_map))
    return split_map


# ── FR-07/08: GECKO 3.0 sign convention for usage / pool_exchange ───────────

# Design Ref: §3.3 — GECKO 3.0 convention
#   usage_prot_{id}     : stoich {prot_i: -1, prot_pool: +1}, lb=-1000, ub=0
#   prot_pool_exchange  : stoich {prot_pool: -1},             lb=-pool_bound, ub=0
#
# Both reactions carry flux in the NEGATIVE direction; ``abs(flux)`` gives the
# protein usage in mg/gDCW.


def add_protein_pool(model: cobra.Model,
                     pool_met_id: str = "prot_pool",
                     pool_rxn_id: str = "prot_pool_exchange",
                     pool_bound_mg_gDCW: float = 0.0,
                     compartment: str = "c") -> None:
    """Add the prot_pool pseudo-metabolite and its exchange reaction.

    GECKO 3.0 Box 1 Step 10 (metabolite) + Step 12 (exchange reaction).

    The exchange reaction is defined as ``prot_pool ->`` with stoichiometric
    coefficient -1, so a negative flux produces protein pool mass (replenishes
    enzyme usage). The ``lb`` sets the maximum pool size in mg/gDCW.

    Idempotent: if the metabolite or reaction already exist they are left
    untouched (except for the ``lb`` being refreshed on the pool exchange).
    """
    if pool_met_id not in model.metabolites:
        met = cobra.Metabolite(
            id=pool_met_id, name="Protein pool", compartment=compartment,
        )
        model.add_metabolites([met])
    met = model.metabolites.get_by_id(pool_met_id)

    if pool_rxn_id not in model.reactions:
        rxn = cobra.Reaction(
            id=pool_rxn_id, name="Protein pool exchange",
            lower_bound=-abs(pool_bound_mg_gDCW), upper_bound=0.0,
        )
        rxn.add_metabolites({met: -1.0})
        model.add_reactions([rxn])
    else:
        rxn = model.reactions.get_by_id(pool_rxn_id)
        rxn.lower_bound = -abs(pool_bound_mg_gDCW)
        rxn.upper_bound = 0.0


def add_enzyme_metabolite_and_usage(model: cobra.Model,
                                     uniprot_id: str,
                                     pool_met_id: str = "prot_pool",
                                     compartment: str = "c"
                                     ) -> tuple[str, str]:
    """Add one ``prot_{uniprot_id}`` metabolite and its ``usage_prot_{id}`` rxn.

    GECKO 3.0 Box 1 Steps 9 & 11.

    Returns
    -------
    (met_id, rxn_id) : IDs of the added (or existing) metabolite and usage rxn.
    """
    met_id = f"prot_{uniprot_id}"
    rxn_id = f"usage_prot_{uniprot_id}"

    if met_id not in model.metabolites:
        met = cobra.Metabolite(
            id=met_id, name=f"Enzyme {uniprot_id}", compartment=compartment,
        )
        model.add_metabolites([met])
    enz_met = model.metabolites.get_by_id(met_id)

    if rxn_id not in model.reactions:
        pool_met = model.metabolites.get_by_id(pool_met_id)
        usage_rxn = cobra.Reaction(
            id=rxn_id, name=f"Enzyme usage {uniprot_id}",
            lower_bound=-1000.0, upper_bound=0.0,
        )
        # GECKO 3.0: stoich {prot_i: -1, prot_pool: +1}
        usage_rxn.add_metabolites({enz_met: -1.0, pool_met: 1.0})
        model.add_reactions([usage_rxn])

    return met_id, rxn_id


# ── FR-20: v1 → v2 schema migration helpers ──────────────────────────────────

def flip_reaction_signs_v1_to_v2(model: cobra.Model) -> list[str]:
    """Flip the sign convention of usage / pool reactions for schema v1 → v2.

    v1 (legacy CNApy) used stoich {prot_pool: -1, prot_i: +1} with bounds
    ``(0, 1000)`` on usage reactions and ``(0, pool_bound)`` on the pool
    exchange. v2 (GECKO 3.0 standard) uses the opposite sign convention so
    both reactions carry flux in the negative direction.

    Algorithm
    ---------
    For every reaction with id starting with ``usage_prot_`` or equal to
    ``prot_pool_exchange``:

    1. Multiply every stoichiometric coefficient by ``-1``.
    2. Swap bounds: ``new_lb = -old_ub``, ``new_ub = -old_lb``.

    This handles all v1 states — default ``(0, 1000)`` becomes ``(-1000, 0)``,
    proteomics-applied ``(0, level)`` becomes ``(-level, 0)``, flex-relaxed
    ``(0, 1000)`` becomes ``(-1000, 0)``, and ``(0, pool_bound)`` becomes
    ``(-pool_bound, 0)``.

    Parameters
    ----------
    model : cobra.Model (modified in place)

    Returns
    -------
    flipped : list of reaction IDs that were modified.
    """
    flipped: list[str] = []
    for rxn in model.reactions:
        is_usage = rxn.id.startswith("usage_prot_")
        is_pool = rxn.id == "prot_pool_exchange"
        if not (is_usage or is_pool):
            continue

        # Negate all stoichiometric coefficients in a single call — cobra
        # emits one solver-side update per add_metabolites, so batching
        # avoids O(n) notifications on large (yeast-GEM scale) models.
        rxn.add_metabolites(
            {met: -2 * coeff for met, coeff in rxn.metabolites.items()},
            combine=True,
        )

        # Swap bounds: new_lb = -old_ub, new_ub = -old_lb.
        old_lb, old_ub = rxn.lower_bound, rxn.upper_bound
        rxn.lower_bound = -old_ub
        rxn.upper_bound = -old_lb

        flipped.append(rxn.id)

    if flipped:
        logger.info("Flipped signs on %d reaction(s) during v1→v2 migration",
                    len(flipped))
    return flipped


# ── Sign-convention helpers for code that needs "capacity" regardless ───────

def usage_capacity(rxn: cobra.Reaction) -> float:
    """Return the upper bound on enzyme usage magnitude, in mg/gDCW.

    Works for both v1 (positive flux: ``upper_bound``) and v2 (negative flux:
    ``-lower_bound``) sign conventions. Call sites that need "how much can
    this enzyme consume of the pool" should use this instead of reading
    ``rxn.upper_bound`` directly.
    """
    if rxn.lower_bound < 0 and rxn.upper_bound <= 0:
        return -rxn.lower_bound  # v2
    return rxn.upper_bound        # v1 (or unconstrained)


def set_usage_capacity(rxn: cobra.Reaction, cap_mg_gDCW: float) -> None:
    """Set the enzyme usage capacity (v2 convention: lb=-cap, ub=0)."""
    rxn.lower_bound = -abs(cap_mg_gDCW)
    rxn.upper_bound = 0.0
