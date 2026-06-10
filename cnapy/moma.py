#!/usr/bin/env python3
#
# Copyright 2022 CNApy organization
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -*- coding: utf-8 -*-

"""
MOMA (Minimization of Metabolic Adjustment) and ROOM (Regulatory On/Off Minimization) implementations.

This module includes:
- linear_moma: Linear MOMA analysis
- room: ROOM analysis (requires MILP solver)
- has_milp_solver: Check for MILP solver availability

This module was enhanced as part of CNApy improvements.
"""

import cobra
from optlang.symbolics import Add


def linear_moma(model: cobra.Model, reference_fluxes: dict[str, float]) -> cobra.Solution:
    """
    Linear MOMA (Minimization of Metabolic Adjustment) analysis.

    Finds a flux distribution that minimizes the Manhattan distance (L1 norm)
    to a reference flux distribution.

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model to analyze
    reference_fluxes : Dict[str, float]
        Reference flux distribution (typically wild-type)

    Returns:
    --------
    cobra.Solution
        The optimized solution
    """
    with model as model:
        # We need to minimize sum(|v - v_ref|)
        # v - v_ref = delta_pos - delta_neg
        # minimize sum(delta_pos + delta_neg)
        # delta_pos, delta_neg >= 0

        objective_terms = []
        vars_to_add = []
        cons_to_add = []

        for rid, v_ref in reference_fluxes.items():
            if rid not in model.reactions:
                continue

            rxn = model.reactions.get_by_id(rid)

            # Skip if reaction is fixed to 0 in both reference and current model (optimization)
            if v_ref == 0 and rxn.lower_bound == 0 and rxn.upper_bound == 0:
                continue

            pos = model.problem.Variable(f"{rid}_pos_slack", lb=0)
            neg = model.problem.Variable(f"{rid}_neg_slack", lb=0)

            vars_to_add.extend([pos, neg])
            objective_terms.extend([pos, neg])

            # Constraint: v - pos + neg = v_ref
            # => v - pos + neg - v_ref = 0
            # In optlang: Constraint(expression, lb=v_ref, ub=v_ref) where expression is v - pos + neg

            expression = rxn.flux_expression - pos + neg
            constraint = model.problem.Constraint(expression, lb=v_ref, ub=v_ref)
            cons_to_add.append(constraint)

        model.add_cons_vars(vars_to_add)
        model.add_cons_vars(cons_to_add)

        model.objective = model.problem.Objective(Add(*objective_terms), direction="min")

        solution = model.optimize()
        return solution


def has_milp_solver() -> tuple[bool, str]:
    """
    Check if a MILP (Mixed-Integer Linear Programming) solver is available.

    Returns:
    --------
    Tuple[bool, str]
        (is_available, solver_name or error_message)
    """

    # Check current solver
    try:
        test_model = cobra.Model()
        solver_interface = test_model.solver.interface.__name__

        # GLPK, CPLEX, and Gurobi support MILP
        milp_solvers = ["glpk", "cplex", "gurobi"]
        for milp in milp_solvers:
            if milp in solver_interface.lower():
                return True, milp

        return False, f"Current solver ({solver_interface}) does not support MILP"
    except Exception as e:
        return False, str(e)


def room(
    model: cobra.Model, reference_fluxes: dict[str, float], delta: float = 0.03, epsilon: float = 0.001
) -> cobra.Solution:
    """
    Regulatory On/Off Minimization (ROOM) analysis.

    ROOM finds a flux distribution that minimizes the number of significant flux changes
    compared to a reference (wild-type) flux distribution. It uses binary variables to
    count reactions that deviate significantly from the reference state.

    Note: Requires a MILP-capable solver (CPLEX, Gurobi, or GLPK).

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model to analyze
    reference_fluxes : Dict[str, float]
        Reference flux distribution (typically wild-type)
    delta : float, optional (default=0.03)
        Relative tolerance for flux changes. A reaction is considered "changed"
        if |v - w| > delta * |w| where v is mutant flux and w is wild-type flux.
    epsilon : float, optional (default=0.001)
        Absolute tolerance for near-zero fluxes.

    Returns:
    --------
    cobra.Solution
        The optimized solution

    Raises:
    -------
    RuntimeError
        If no MILP solver is available
    """
    # Check for MILP solver
    has_milp, msg = has_milp_solver()
    if not has_milp:
        raise RuntimeError(f"ROOM requires a MILP solver. {msg}")

    with model as model:
        binary_vars = []
        vars_to_add = []
        cons_to_add = []

        for rid, w in reference_fluxes.items():
            if rid not in model.reactions:
                continue

            rxn = model.reactions.get_by_id(rid)
            lb, ub = rxn.bounds

            # Create binary indicator variable
            y = model.problem.Variable(f"{rid}_y", type="binary")
            binary_vars.append(y)
            vars_to_add.append(y)

            # Calculate allowable flux range
            if abs(w) < epsilon:
                w_upper = epsilon
                w_lower = -epsilon
            else:
                w_upper = w + delta * abs(w)
                w_lower = w - delta * abs(w)

            # Big-M values
            # Cap the Big-M at a finite value when bounds are large or infinite
            # (e.g. exchange reactions with +/-inf bounds), mirroring the
            # reference implementation in test_data/Simulator.py. An infinite
            # coefficient would otherwise be rejected/mis-scaled by the backend.
            if ub is not None and ub < 1000:
                M_upper = max(abs(ub - w_upper), 1000)
            else:
                M_upper = 2000
            if lb is not None and lb > -1000:
                M_lower = max(abs(w_lower - lb), 1000)
            else:
                M_lower = 2000

            # Constraints using Big-M formulation:
            # v <= w_upper + M_upper * y  (when y=0: v <= w_upper)
            # v >= w_lower - M_lower * y  (when y=0: v >= w_lower)

            upper_constraint = model.problem.Constraint(
                rxn.flux_expression - M_upper * y, ub=w_upper, name=f"room_upper_{rid}"
            )

            lower_constraint = model.problem.Constraint(
                rxn.flux_expression + M_lower * y, lb=w_lower, name=f"room_lower_{rid}"
            )

            cons_to_add.extend([upper_constraint, lower_constraint])

        model.add_cons_vars(vars_to_add)
        model.add_cons_vars(cons_to_add)

        # Objective: minimize number of significantly changed reactions
        model.objective = model.problem.Objective(Add(*binary_vars), direction="min")

        solution = model.optimize()
        return solution
