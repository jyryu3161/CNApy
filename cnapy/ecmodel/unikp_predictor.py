"""
unikp_predictor.py
==================
Optional UniKP-based kcat prediction for reactions that lack experimental data.

UniKP uses:
  - ProtT5-XL-UniRef50  (protein language model, ~3 GB download on first use)
  - SMILES Transformer  (substrate language model)
  - Extra Trees ensemble (scikit-learn)

Reference: Han Yu et al., Nature Communications 14, 8211 (2023)
GitHub:    https://github.com/Luo-SynBioLab/UniKP

This module wraps the UniKP Python API.  If UniKP or its dependencies are not
installed the functions return gracefully with an informative message.
"""

from __future__ import annotations

from typing import Callable

# ── dependency check ──────────────────────────────────────────────────────────

def check_unikp_available() -> tuple[bool, str]:
    """
    Verify that all UniKP runtime dependencies are importable.

    Returns (available: bool, message: str).
    """
    missing = []
    for pkg, import_name in [
        ("torch",        "torch"),
        ("transformers", "transformers"),
        ("scikit-learn", "sklearn"),
        ("UniKP",        "UniKP"),
    ]:
        try:
            __import__(import_name)
        except ImportError:
            missing.append(pkg)

    if missing:
        return False, (
            "The following packages are required for UniKP but are not installed:\n"
            + "  " + ", ".join(missing) + "\n\n"
            "Install with:\n"
            "  pip install torch transformers scikit-learn UniKP"
        )
    return True, "UniKP is available."


# ── prediction ────────────────────────────────────────────────────────────────

def predict_kcat_unikp(
    sequences: dict[str, str],
    smiles: dict[str, str],
    progress_callback: Callable[[int, int], None] | None = None,
) -> tuple[dict[str, float], list[str]]:
    """
    Predict kcat values using UniKP.

    Parameters
    ----------
    sequences : { identifier -> amino_acid_sequence }
        Enzyme amino acid sequences (one-letter code).
    smiles    : { identifier -> SMILES_string }
        Substrate SMILES strings, keyed by the same identifiers.
    progress_callback : callable(current, total) or None
        Called after each prediction so a progress bar can be updated.

    Returns
    -------
    predictions : { identifier -> predicted_kcat_1/s }
    failed      : list of identifiers for which prediction failed
    """
    available, msg = check_unikp_available()
    if not available:
        raise ImportError(msg)

    try:
        from UniKP import UniKP as _UniKP  # noqa: N813
    except ImportError as exc:
        raise ImportError(
            "UniKP package could not be imported.  "
            "Install it with: pip install UniKP"
        ) from exc

    common_ids = set(sequences) & set(smiles)
    predictions: dict[str, float] = {}
    failed: list[str] = []
    total = len(common_ids)

    for idx, uid in enumerate(sorted(common_ids)):
        if progress_callback:
            progress_callback(idx, total)
        try:
            kcat_log10 = _UniKP.predict_kcat(
                protein_sequence=sequences[uid],
                smiles=smiles[uid],
            )
            # UniKP returns log10(kcat); convert to kcat in 1/s
            import math
            predictions[uid] = 10 ** kcat_log10
        except Exception:  # noqa: BLE001
            failed.append(uid)

    if progress_callback:
        progress_callback(total, total)

    return predictions, failed
