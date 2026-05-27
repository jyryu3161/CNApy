"""Tests for omics integration: GPR evaluation and E-Flux2."""

import ast

import pytest

from cnapy.gui_elements.omics_integration_dialog import (
    evaluate_gpr_expression,
    gene_expression_to_reaction_weights,
    run_eflux2,
)


def _gpr(expr: str) -> ast.AST:
    """Parse a Python boolean expression into a GPR-style AST node."""
    return ast.parse(expr, mode="eval").body


class TestEvaluateGPRExpression:
    def test_single_gene(self):
        assert evaluate_gpr_expression(_gpr("g1"), {"g1": 5.0}) == 5.0

    def test_missing_gene_returns_none(self):
        assert evaluate_gpr_expression(_gpr("g1"), {}) is None

    def test_and_takes_min(self):
        assert evaluate_gpr_expression(_gpr("g1 and g2"), {"g1": 3.0, "g2": 7.0}) == 3.0

    def test_and_with_missing_subunit_is_none(self):
        # Enzyme complex cannot form if any subunit is missing.
        assert evaluate_gpr_expression(_gpr("g1 and g2"), {"g1": 3.0}) is None

    def test_or_takes_max(self):
        assert evaluate_gpr_expression(_gpr("g1 or g2"), {"g1": 3.0, "g2": 7.0}) == 7.0

    def test_or_with_partial_missing_uses_available(self):
        # Missing isozyme is tolerated as long as at least one alternative exists.
        assert evaluate_gpr_expression(_gpr("g1 or g2"), {"g1": 3.0}) == 3.0

    def test_or_all_missing_is_none(self):
        assert evaluate_gpr_expression(_gpr("g1 or g2"), {}) is None

    def test_nested_and_of_or(self):
        # (g1 or g2) and g3 → min(max(g1, g2), g3)
        node = _gpr("(g1 or g2) and g3")
        assert evaluate_gpr_expression(node, {"g1": 1, "g2": 5, "g3": 3}) == 3.0


class TestGeneExpressionToReactionWeights:
    def test_uniform_expression(self, ecoli_core_model):
        with ecoli_core_model as model:
            expr = {g.id: 1.0 for g in model.genes}
            weights = gene_expression_to_reaction_weights(model, expr)
        assert len(weights) > 0
        assert all(abs(w - 1.0) < 1e-9 for w in weights.values())

    def test_reaction_without_genes_is_skipped(self, simple_model):
        # simple_model has no genes, so no weights are produced.
        weights = gene_expression_to_reaction_weights(simple_model, {"gA": 1.0})
        assert weights == {}


class TestRunEflux2:
    def test_empty_weights(self, ecoli_core_model):
        with ecoli_core_model as model:
            status, obj, fluxes, method = run_eflux2(model, {})
        assert status == "no_targets"
        assert obj is None
        assert fluxes is None
        assert method == ""

    def test_unknown_reactions_only(self, ecoli_core_model):
        with ecoli_core_model as model:
            status, *_ = run_eflux2(model, {"NOT_A_REAL_RXN": 1.0})
        assert status == "no_targets"

    def test_uniform_expression_preserves_growth(self, ecoli_core_model):
        # With uniform expression on every internal reaction, normalized = 1.0
        # everywhere → no effective bound tightening → objective matches FBA.
        with ecoli_core_model as model:
            fba_obj = model.optimize().objective_value
            weights = {r.id: 1.0 for r in model.reactions if not r.id.startswith("EX_")}
            status, obj, fluxes, method = run_eflux2(model, weights)
        assert status == "optimal"
        assert obj == pytest.approx(fba_obj, rel=1e-3)
        assert method in {"qp", "pfba", "fba"}
        assert fluxes is not None

    def test_downregulation_reduces_flux(self, ecoli_core_model):
        # Pick a reaction with substantial FBA flux and a GPR, downregulate
        # it severely, and confirm its flux drops below baseline.
        with ecoli_core_model as model:
            target_id = "PFK"
            if target_id not in {r.id for r in model.reactions}:
                pytest.skip("PFK not in model")
            baseline = abs(model.optimize().fluxes[target_id])
            if baseline < 1e-6:
                pytest.skip("PFK baseline flux is essentially zero")

            weights = {r.id: 1.0 for r in model.reactions}
            weights[target_id] = 0.001  # 0.1% of max → tight bound
            status, _obj, fluxes, _method = run_eflux2(model, weights, min_scale=0.001)

        assert status == "optimal"
        assert fluxes is not None
        assert abs(fluxes[target_id]) < baseline

    def test_input_model_not_mutated(self, ecoli_core_model):
        # Run E-Flux2 with non-trivial bound changes and confirm the input
        # model's reaction bounds and objective are restored on return.
        with ecoli_core_model as model:
            # Pick any reaction with a GPR.
            target = next(r for r in model.reactions if r.gpr and r.gpr.body)
            original_bounds = target.bounds
            original_obj = str(model.objective.expression)

            weights = {r.id: 1.0 for r in model.reactions}
            weights[target.id] = 0.01
            run_eflux2(model, weights)

            assert target.bounds == original_bounds
            assert str(model.objective.expression) == original_obj

    def test_percentile_normalization_invalid_value(self, ecoli_core_model):
        with ecoli_core_model as model:
            weights = {r.id: 1.0 for r in model.reactions}
            with pytest.raises(ValueError):
                run_eflux2(model, weights, normalization_percentile=0.0)
            with pytest.raises(ValueError):
                run_eflux2(model, weights, normalization_percentile=150.0)

    def test_percentile_parameter_changes_result_under_outlier(self, ecoli_core_model):
        # Inject a single outlier weight; the denominator at p100 (max) and at
        # p99 must differ → bounds differ → either the objective or the flux
        # distribution must observably change.  Whether biomass *increases*
        # at p99 depends on which bounds bind in the L2 problem; the test
        # below only verifies the mechanism is wired through, leaving the
        # biological direction to the JS66 integration sanity-check.
        with ecoli_core_model as model:
            weights = {r.id: 1.0 for r in model.reactions}
            target_id = next(iter(weights))
            weights[target_id] = 100.0  # strong outlier

            _, obj_max, f_max, _ = run_eflux2(
                model, weights, weight_threshold=0.0, normalization_percentile=100.0
            )
            _, obj_p99, f_p99, _ = run_eflux2(
                model, weights, weight_threshold=0.0, normalization_percentile=99.0
            )
        assert obj_max is not None and obj_p99 is not None
        fluxes_differ = (
            f_max is not None
            and f_p99 is not None
            and any(abs(f_max[r] - f_p99[r]) > 1e-6 for r in f_max)
        )
        obj_differs = abs(obj_max - obj_p99) > 1e-6
        assert fluxes_differ or obj_differs, (
            "normalization_percentile must change either fluxes or objective when "
            "outliers are present (mechanism appears un-wired)"
        )

    def test_percentile_paper_default_matches_max_normalization(self, ecoli_core_model):
        # Sanity check: the default percentile=100 must reproduce the
        # paper-faithful behavior identical to explicitly calling with
        # normalization_percentile=100.0.
        with ecoli_core_model as model:
            weights = {r.id: float(i + 1) for i, r in enumerate(model.reactions)}
            _, obj_default, _, _ = run_eflux2(model, weights, weight_threshold=0.0)
            _, obj_100, _, _ = run_eflux2(
                model, weights, weight_threshold=0.0, normalization_percentile=100.0
            )
        assert obj_default == pytest.approx(obj_100, rel=1e-9)

    def test_scenario_constraints_take_precedence(self, ecoli_core_model):
        # A reaction listed in flux_constraints should keep those bounds even
        # if its expression weight would otherwise scale it.
        with ecoli_core_model as model:
            target_id = "PFK"
            if target_id not in {r.id for r in model.reactions}:
                pytest.skip("PFK not in model")
            forced = (2.0, 2.0)
            weights = {r.id: 1.0 for r in model.reactions}
            weights[target_id] = 0.0001  # would normally tighten bound

            status, _obj, fluxes, _method = run_eflux2(
                model,
                weights,
                flux_constraints={target_id: forced},
            )

        assert status == "optimal"
        assert fluxes is not None
        assert fluxes[target_id] == pytest.approx(2.0, abs=1e-6)
