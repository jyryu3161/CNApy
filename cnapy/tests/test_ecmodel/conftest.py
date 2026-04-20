"""Fixtures for ecmodel tests.

Design Ref: §8.5 — in-memory toy models avoid shipping XML fixtures while
still exercising all Stage 1 transformations.
"""

from __future__ import annotations

import cobra
import pytest


@pytest.fixture
def toy_gem() -> cobra.Model:
    """A small GEM with enough structure to test pseudoreaction + isozyme logic.

    Design:
      * ``r_in``    : exchange-like import of A (no GPR)
      * ``r_iso``   : A → B, catalysed by two isozymes ``g1`` and ``g2``
                      (classic OR case → should split into ``_EXP_1/2``).
      * ``r_cplx``  : B → C, catalysed by a complex ``(g3 and g4) or g5``
                      (mixed case → ``_EXP_1`` is complex, ``_EXP_2`` is g5).
      * ``r_simple``: C → D, single gene g6 (must NOT be split).
      * ``r_out``   : D → (no GPR).
      * ``BIOMASS`` : D → biomass, with a fake GPR ``g_bio``, marked as a
                      pseudoreaction via its name so FR-01 clears it.
    """
    model = cobra.Model("toy_gem")

    mets = {
        m: cobra.Metabolite(m, compartment="c") for m in ["A", "B", "C", "D", "BIO"]
    }
    model.add_metabolites(list(mets.values()))

    r_in = cobra.Reaction("r_in")
    r_in.add_metabolites({mets["A"]: 1})
    r_in.bounds = (0, 1000)

    r_iso = cobra.Reaction("r_iso")
    r_iso.name = "Isozymic conversion A to B"
    r_iso.add_metabolites({mets["A"]: -1, mets["B"]: 1})
    r_iso.bounds = (0, 1000)
    r_iso.gene_reaction_rule = "g1 or g2"

    r_cplx = cobra.Reaction("r_cplx")
    r_cplx.name = "Complex or single B to C"
    r_cplx.add_metabolites({mets["B"]: -1, mets["C"]: 1})
    r_cplx.bounds = (0, 1000)
    r_cplx.gene_reaction_rule = "(g3 and g4) or g5"

    r_simple = cobra.Reaction("r_simple")
    r_simple.name = "Single gene C to D"
    r_simple.add_metabolites({mets["C"]: -1, mets["D"]: 1})
    r_simple.bounds = (0, 1000)
    r_simple.gene_reaction_rule = "g6"

    r_out = cobra.Reaction("r_out")
    r_out.add_metabolites({mets["D"]: -1})
    r_out.bounds = (0, 1000)

    biomass = cobra.Reaction("BIOMASS")
    biomass.name = "biomass pseudoreaction"   # triggers FR-01 via name match
    biomass.add_metabolites({mets["D"]: -1, mets["BIO"]: 1})
    biomass.bounds = (0, 1000)
    # Intentionally a fake GPR to verify FR-01 clears it.
    biomass.gene_reaction_rule = "g_bio"

    ex_bio = cobra.Reaction("EX_BIO")
    ex_bio.add_metabolites({mets["BIO"]: -1})
    ex_bio.bounds = (0, 1000)

    model.add_reactions([r_in, r_iso, r_cplx, r_simple, r_out, biomass, ex_bio])
    model.objective = "EX_BIO"
    return model


@pytest.fixture
def toy_reversible_gem() -> cobra.Model:
    """A tiny GEM with a reversible isozymic reaction.

    Used to verify that reversible-split + isozyme-split compose correctly
    (FR-03 + FR-04 interaction).
    """
    model = cobra.Model("toy_reversible")
    a = cobra.Metabolite("A", compartment="c")
    b = cobra.Metabolite("B", compartment="c")
    model.add_metabolites([a, b])

    r = cobra.Reaction("r_rev")
    r.add_metabolites({a: -1, b: 1})
    r.bounds = (-1000, 1000)   # reversible
    r.gene_reaction_rule = "gA or gB"

    ex_a = cobra.Reaction("EX_A")
    ex_a.add_metabolites({a: -1})
    ex_a.bounds = (-10, 1000)

    ex_b = cobra.Reaction("EX_B")
    ex_b.add_metabolites({b: -1})
    ex_b.bounds = (0, 1000)

    model.add_reactions([r, ex_a, ex_b])
    model.objective = "EX_B"
    return model
