"""Nearest-epoch IXPE IRF calibration resolver.

ixpeobssim carries no concept of IRF validity dates -- its own CALDB helper
(``ixpeobssim/irf/caldb.py``) states plainly "we don't really have a concept of
validity date" -- so the reduction caller must choose which ``obssimYYYYMMDD``
calibration epoch to use. Hand-filling that date is error prone: the wrong epoch
silently biases the recovered polarization (a Cyg X-1 May-2022 reproduction run
that used a mid-2023 epoch recovered a polarization degree ~2 sigma below the
published value). This module removes the choice from the caller. Given an
observation's date it selects the installed calibration epoch NEAREST in time,
so omitting ``irf_date`` yields the temporally appropriate response and an
explicit ``irf_date`` remains an override.

The epoch strings are the ``obssimYYYYMMDD`` tokens embedded in the installed
CALDB MRF filenames (each the "first date of validity" of that calibration, per
the HEASARC naming convention ixpeobssim approximately follows). The active
ixpeobssim installation is located from ``ixpeobssim.__file__`` -- never a
hardcoded path -- mirroring ``red_ixpe._installed_du_rmf``.
"""

from __future__ import annotations

import glob
import os
import re
from datetime import date

# MJD 0 == 1858-11-17; ``date.toordinal`` counts from 0001-01-01, so subtracting
# that day's ordinal converts a calendar date to an (integer-day) MJD with only
# the standard library -- keeping this module importable without the science stack.
_MJD_ORIGIN_ORDINAL = date(1858, 11, 17).toordinal()

_EPOCH_RE = re.compile(r"obssim(\d{8})_")


def _epoch_to_mjd(epoch: str) -> float:
    """Convert an ``'YYYYMMDD'`` epoch token to its MJD (integer day)."""
    d = date(int(epoch[0:4]), int(epoch[4:6]), int(epoch[6:8]))
    return float(d.toordinal() - _MJD_ORIGIN_ORDINAL)


def _caldb_mrf_dir() -> str:
    """MRF directory of the currently active ixpeobssim installation.

    Imported lazily so this module stays importable (and unit-testable) in an
    environment without the science stack; only auto-resolution needs it.
    """
    import ixpeobssim

    return os.path.join(
        os.path.dirname(os.path.abspath(ixpeobssim.__file__)),
        "caldb", "ixpe", "gpd", "cpf", "mrf",
    )


def installed_irf_epochs(caldb_mrf_dir: str | None = None) -> list[str]:
    """Sorted, de-duplicated ``'YYYYMMDD'`` epochs present in the CALDB MRF dir.

    Every calibration flavour (``alpha075``/``gray``/...) and DU repeats the same
    epoch token, so the set collapses to one entry per epoch. ``caldb_mrf_dir``
    is injectable for testing; it defaults to the active ixpeobssim installation.
    """
    mrf_dir = caldb_mrf_dir if caldb_mrf_dir is not None else _caldb_mrf_dir()
    epochs: set[str] = set()
    for name in glob.glob(os.path.join(mrf_dir, "ixpe_d*_obssim*_v0*.mrf")):
        match = _EPOCH_RE.search(os.path.basename(name))
        if match:
            epochs.add(match.group(1))
    return sorted(epochs)


def resolve_irf_date(
    obs_mjd: float,
    epochs: list[str] | None = None,
    caldb_mrf_dir: str | None = None,
) -> str:
    """Return the installed IRF epoch (``'YYYYMMDD'``) nearest the observation.

    Nearest by absolute day-distance between ``obs_mjd`` and each epoch's first
    date of validity; the nearest epoch may fall slightly *after* the observation
    (IXPE calibration drifts smoothly, so the closest-in-time calibration best
    represents the detector state). A tie resolves to the earlier epoch, whose
    calibration window has already opened. ``epochs`` is injectable for testing;
    it defaults to :func:`installed_irf_epochs`.

    Raises
    ------
    RuntimeError
        No IRF epoch is installed, so nothing can be selected.
    """
    epochs = epochs if epochs is not None else installed_irf_epochs(caldb_mrf_dir)
    if not epochs:
        raise RuntimeError(
            "No installed IXPE IRF epochs found (looked for obssimYYYYMMDD MRFs in "
            "the active ixpeobssim CALDB); cannot auto-resolve irf_date -- pass one "
            "explicitly."
        )
    # ``sorted`` ascending + ``min`` keeps the first (earliest) entry on a tie.
    return min(sorted(epochs), key=lambda e: abs(obs_mjd - _epoch_to_mjd(e)))
