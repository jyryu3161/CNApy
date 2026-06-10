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

"""Flux Sampling module for CNApy

This module provides two sampling methods:

1. Random Sampling:
   - Standard flux sampling using hit-and-run algorithm
   - Samples uniformly from the feasible flux space

2. Predicted Flux-Based Sampling:
   - Uses a reference flux distribution (from FBA, MOMA, etc.) as the mean
   - Samples with Gaussian noise around the predicted fluxes
   - Useful for uncertainty quantification around a predicted solution
"""

import cobra
import numpy as np
import pandas as pd
from cobra.sampling import sample


def perform_sampling(model: cobra.Model, n: int, thinning: int, processes: int) -> pd.DataFrame:
    """
    Perform standard random flux sampling.

    Uses the OptGP (Optim-GP) sampler by default for efficient sampling.

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    n : int
        Number of samples
    thinning : int
        Thinning factor (keep every nth sample)
    processes : int
        Number of parallel processes

    Returns:
    --------
    pd.DataFrame
        DataFrame with samples (rows) and reactions (columns)
    """
    with model as model:
        s = sample(model, n, thinning=thinning, processes=processes)
    return s


def perform_predicted_flux_sampling(
    model: cobra.Model,
    reference_fluxes: dict[str, float],
    n: int,
    std_fraction: float = 0.1,
    constraint_mode: str = "bounds",
    min_fraction: float = 0.8,
    max_fraction: float = 1.2,
    thinning: int = 100,
    processes: int = 4,
) -> tuple[pd.DataFrame, dict[str, tuple[float, float]]]:
    """
    Perform flux sampling centered around predicted flux values.

    This method constrains the flux space around a reference solution
    (e.g., from FBA, MOMA, LAD) and samples within that constrained space.

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    reference_fluxes : Dict[str, float]
        Reference flux values (e.g., from FBA solution)
    n : int
        Number of samples
    std_fraction : float
        Standard deviation as fraction of flux value for Gaussian noise
        (used in post-processing)
    constraint_mode : str
        How to constrain the space:
        - "bounds": Constrain bounds to [min_fraction * flux, max_fraction * flux]
        - "fixed": Fix fluxes to reference values (no sampling)
        - "free": Don't constrain (same as random sampling)
    min_fraction : float
        Minimum fraction of reference flux (for "bounds" mode)
    max_fraction : float
        Maximum fraction of reference flux (for "bounds" mode)
    thinning : int
        Thinning factor
    processes : int
        Number of parallel processes

    Returns:
    --------
    Tuple[pd.DataFrame, Dict[str, Tuple[float, float]]]
        (samples DataFrame, applied bounds for each reaction)
    """
    applied_bounds = {}

    with model as m:
        if constraint_mode == "bounds":
            # Constrain each reaction to a range around the reference flux
            for rid, flux in reference_fluxes.items():
                if rid in m.reactions:
                    rxn = m.reactions.get_by_id(rid)

                    # Keep the reaction's original feasible bounds so the
                    # constrained window can never be relaxed outside them.
                    orig_lb = rxn.lower_bound
                    orig_ub = rxn.upper_bound

                    if abs(flux) < 1e-6:
                        # Near-zero flux: constrain to small range around zero
                        new_lb = max(rxn.lower_bound, -0.1)
                        new_ub = min(rxn.upper_bound, 0.1)
                    else:
                        # Non-zero flux: constrain to fraction range
                        if flux >= 0:
                            new_lb = max(rxn.lower_bound, flux * min_fraction)
                            new_ub = min(rxn.upper_bound, flux * max_fraction)
                        else:
                            new_lb = max(rxn.lower_bound, flux * max_fraction)  # Note: reversed for negative
                            new_ub = min(rxn.upper_bound, flux * min_fraction)

                    # An empty intersection (new_lb > new_ub) means the
                    # requested window does not overlap the reaction's
                    # existing bounds. Swapping the values would advertise
                    # (and sample over) a range outside the model's own
                    # feasible bounds, silently relaxing a hard constraint.
                    # Instead, clamp the window into [orig_lb, orig_ub] so
                    # the applied bounds always stay within the model's
                    # original feasible range.
                    if new_lb > new_ub:
                        new_lb = min(max(new_lb, orig_lb), orig_ub)
                        new_ub = max(min(new_ub, orig_ub), orig_lb)
                        if new_lb > new_ub:
                            new_lb, new_ub = new_ub, new_lb

                    rxn.lower_bound = new_lb
                    rxn.upper_bound = new_ub
                    applied_bounds[rid] = (new_lb, new_ub)

        elif constraint_mode == "fixed":
            # Fix all reference fluxes
            for rid, flux in reference_fluxes.items():
                if rid in m.reactions:
                    rxn = m.reactions.get_by_id(rid)
                    rxn.lower_bound = flux
                    rxn.upper_bound = flux
                    applied_bounds[rid] = (flux, flux)

        # Perform sampling
        try:
            samples = sample(m, n, thinning=thinning, processes=processes)
        except Exception as e:
            # If sampling fails, try with fewer samples
            try:
                samples = sample(m, min(n, 100), thinning=thinning, processes=1)
            except Exception:
                raise RuntimeError(f"Sampling failed: {str(e)}")

    return samples, applied_bounds


def compute_sampling_statistics(samples: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.DataFrame]:
    """
    Compute statistics from sampling results.

    Parameters:
    -----------
    samples : pd.DataFrame
        Samples from flux sampling

    Returns:
    --------
    Tuple[pd.Series, pd.Series, pd.DataFrame]
        (mean, std, correlation matrix)
    """
    mean = samples.mean()
    std = samples.std()

    # Compute correlation for reactions with variance
    high_var_cols = std[std > 1e-6].index
    if len(high_var_cols) > 0:
        corr = samples[high_var_cols].corr()
    else:
        corr = pd.DataFrame()

    return mean, std, corr


def add_gaussian_noise_to_samples(
    samples: pd.DataFrame, reference_fluxes: dict[str, float], std_fraction: float = 0.1
) -> pd.DataFrame:
    """
    Add Gaussian noise to samples centered around reference fluxes.

    This post-processing step can be used to add uncertainty
    to a deterministic solution.

    Parameters:
    -----------
    samples : pd.DataFrame
        Original samples
    reference_fluxes : Dict[str, float]
        Reference flux values
    std_fraction : float
        Standard deviation as fraction of absolute flux value

    Returns:
    --------
    pd.DataFrame
        Samples with added Gaussian noise
    """
    noisy_samples = samples.copy()

    for rid, flux in reference_fluxes.items():
        if rid in noisy_samples.columns:
            std = abs(flux) * std_fraction if abs(flux) > 1e-6 else std_fraction
            noise = np.random.normal(0, std, len(noisy_samples))
            noisy_samples[rid] = noisy_samples[rid] + noise

    return noisy_samples
