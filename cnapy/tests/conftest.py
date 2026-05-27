"""Shared pytest fixtures for CNApy tests."""

import cobra
import numpy as np
import pandas as pd
import pytest

# Pin tests to GLPK so the suite is reproducible regardless of the developer's
# licensed solver state (e.g. an expired Gurobi license would otherwise break
# bare ``cobra.Model(...)`` calls).
try:
    cobra.Configuration().solver = "glpk"
except Exception:
    pass


@pytest.fixture(scope="session")
def ecoli_core_model():
    """Load E. coli core model for testing.

    This fixture is session-scoped for performance as model loading is expensive.
    Tests should use the model within a context manager to avoid modifying it.
    """
    try:
        model = cobra.io.load_model("textbook")
        return model
    except Exception:
        pytest.skip("E. coli core model not available")


@pytest.fixture
def empty_model():
    """Create an empty cobra.Model for testing."""
    return cobra.Model("empty_test_model")


@pytest.fixture
def simple_model():
    """Create a simple model with a few reactions for testing."""
    model = cobra.Model("simple_model")

    # Create metabolites
    A = cobra.Metabolite("A", compartment="c")
    B = cobra.Metabolite("B", compartment="c")
    C = cobra.Metabolite("C", compartment="c")

    # Create reactions
    r1 = cobra.Reaction("R1")
    r1.add_metabolites({A: -1, B: 1})
    r1.bounds = (-1000, 1000)

    r2 = cobra.Reaction("R2")
    r2.add_metabolites({B: -1, C: 1})
    r2.bounds = (-1000, 1000)

    r3 = cobra.Reaction("EX_A")
    r3.add_metabolites({A: -1})
    r3.bounds = (-10, 1000)

    r4 = cobra.Reaction("EX_C")
    r4.add_metabolites({C: -1})
    r4.bounds = (0, 1000)

    model.add_reactions([r1, r2, r3, r4])
    model.objective = "R2"

    return model


@pytest.fixture
def wildtype_fluxes(ecoli_core_model):
    """Get FBA solution fluxes for the E. coli core model."""
    with ecoli_core_model as model:
        solution = model.optimize()
        if solution.status == "optimal":
            return solution.fluxes.to_dict()
    return {}


@pytest.fixture
def sample_dataframe():
    """Create a sample DataFrame for testing sampling statistics."""
    np.random.seed(42)
    n_samples = 100
    n_reactions = 5

    data = {f"R{i}": np.random.normal(loc=i, scale=0.5, size=n_samples) for i in range(n_reactions)}
    return pd.DataFrame(data)


@pytest.fixture
def flux_vector_data():
    """Create sample flux vector data for FluxVectorContainer tests."""
    np.random.seed(42)
    n_vectors = 10
    n_reactions = 5

    fv_mat = np.random.randn(n_vectors, n_reactions)
    reac_id = [f"R{i}" for i in range(n_reactions)]

    return fv_mat, reac_id


@pytest.fixture
def simple_scenario_dict():
    """Create a simple scenario dictionary for testing."""
    return {
        "R1": (5.0, 5.0),
        "R2": (0.0, 10.0),
        "R3": (-5.0, 5.0),
    }
