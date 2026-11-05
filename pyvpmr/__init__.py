#  Copyright (C) 2024-2026 Theodore Chang
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Python interface to the VPMR algorithm.

VPMR approximates a smooth kernel function :math:`g(t)` by a finite sum of
exponentials,

.. math::

    g(t) \\approx \\sum_{j=1}^{N} m_j \\exp(-s_j t).

The main entry point is :func:`vpmr`, which returns the coefficients ``m`` and
the exponents ``s`` as complex NumPy arrays.  :func:`evaluate`, :func:`plot`,
:func:`split`, and the ``to_*_damping`` helpers build on that result.
"""

from __future__ import annotations

import operator
import re
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Callable, NamedTuple

import numpy as np

try:  # Normal wheel layout: the extension is installed as a top-level module.
    import _pyvpmr
except ImportError:  # pragma: no cover - supports alternative build layouts.
    try:
        from . import _pyvpmr  # type: ignore[no-redef]
    except ImportError as exc:
        raise ImportError(
            "The compiled extension '_pyvpmr' is not available. "
            "Install pyvpmr from a wheel or build it from source."
        ) from exc

__all__ = [
    "DEFAULT_KERNEL",
    "VPMRResult",
    "evaluate",
    "plot",
    "split",
    "to_elemental_damping",
    "to_global_damping",
    "vpmr",
]

try:
    __version__ = version("pyvpmr")
except PackageNotFoundError:  # pragma: no cover - running from a source tree.
    __version__ = "0.0.0"

#: Kernel expression used by :func:`vpmr` when no kernel is supplied.
DEFAULT_KERNEL = "exp(-t*t/4)"

_COMPLEX_LINE = re.compile(
    r"^\s*"
    r"([+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?)"
    r"([+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?)"
    r"[jJ]\s*$"
)


class VPMRResult(NamedTuple):
    """Result of :func:`vpmr`.

    The result behaves like the historical ``(m, s)`` tuple, so both
    ``m, s = vpmr(...)`` and ``result = vpmr(...); result.m`` are supported.
    """

    m: np.ndarray
    """Coefficients :math:`m_j` as a one-dimensional complex array."""

    s: np.ndarray
    """Exponents :math:`s_j` as a one-dimensional complex array."""


def _as_int(name: str, value: object, minimum: int) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an integer, not bool")
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer, got {type(value).__name__}") from exc
    if result < minimum:
        raise ValueError(f"{name} must be greater than or equal to {minimum}, got {result}")
    return result


def _as_float(name: str, value: object, minimum: float, *, inclusive: bool = True) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a real number, not bool")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a real number, got {type(value).__name__}") from exc
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite, got {result}")
    if result < minimum or (not inclusive and result == minimum):
        relation = ">=" if inclusive else ">"
        raise ValueError(f"{name} must be {relation} {minimum}, got {result}")
    return result


def vpmr(
    *,
    n: int = 20,
    c: int = 4,
    d: int = 0,
    q: int = 500,
    m: float = 1.05,
    e: float = 1e-8,
    k: str = "",
    omit: bool = True,
    verbose: bool = False,
) -> VPMRResult:
    """Run the VPMR algorithm.

    Args:
        n: Maximum number of exponential terms (at least ``1``).
        c: Maximum exponent/order used by the algorithm (at least ``1``).
        d: Precision in bits.  Use ``0`` for automatic precision selection.
        q: Gauss-Legendre quadrature order (at least ``1``).
        m: Precision multiplier (at least ``1.05``).
        e: Convergence tolerance, strictly positive.
        k: Kernel expression in terms of ``t``.  An empty string uses
            :data:`DEFAULT_KERNEL`.
        omit: Whether to omit numerically trivial terms from the result.
        verbose: Whether to print progress information and the effective
            parameters.  The Python API is quiet by default.

    Returns:
        A :class:`VPMRResult` containing the coefficient array ``m`` and the
        exponent array ``s``.

    Raises:
        TypeError: If an argument has an incompatible type.
        ValueError: If an argument is outside its supported range.
        RuntimeError: If the kernel expression cannot be compiled or the
            algorithm fails.
    """
    n = _as_int("n", n, 1)
    c = _as_int("c", c, 1)
    d = _as_int("d", d, 0)
    q = _as_int("q", q, 1)
    m = _as_float("m", m, 1.05)
    e = _as_float("e", e, 0.0, inclusive=False)
    if not isinstance(k, str):
        raise TypeError(f"k must be a string, got {type(k).__name__}")

    result = _pyvpmr.vpmr(
        n=n,
        c=c,
        d=d,
        q=q,
        m=m,
        e=e,
        k=k,
        omit=bool(omit),
        verbose=bool(verbose),
    )

    try:
        raw_m, raw_s = result
    except (TypeError, ValueError) as exc:  # pragma: no cover - defensive.
        raise RuntimeError("The VPMR extension returned an invalid result.") from exc

    return VPMRResult(
        np.asarray(raw_m, dtype=np.complex128).reshape(-1),
        np.asarray(raw_s, dtype=np.complex128).reshape(-1),
    )


def evaluate(
    m: np.ndarray | list[complex] | VPMRResult,
    s: np.ndarray | list[complex],
    t: float | np.ndarray,
) -> complex | np.ndarray:
    """Evaluate the exponential approximation :math:`\\sum_j m_j e^{-s_j t}`.

    Args:
        m: Coefficients of the approximation.
        s: Exponents of the approximation.
        t: Scalar or array of time values.

    Returns:
        A complex scalar when ``t`` is scalar, otherwise a complex array with
        the same shape as ``t``.
    """
    m_array, s_array = _process_args((m, s))
    time = np.asarray(t, dtype=float)
    result = m_array @ np.exp(-np.multiply.outer(s_array, time))
    return result.item() if result.ndim == 0 else result


def split(result: str) -> VPMRResult:
    """Parse the textual output of the standalone ``vpmr`` program.

    The standalone program prints the coefficients after an ``M = `` line and
    the exponents after an ``S = `` line.  This helper extracts those two
    blocks without depending on whitespace or fixed-width formatting.

    Args:
        result: Raw standard output produced by the ``vpmr`` executable.

    Returns:
        A :class:`VPMRResult` parsed from the output.

    Raises:
        TypeError: If ``result`` is not a string.
        ValueError: If no values, or an odd number of values, can be parsed.
    """
    if not isinstance(result, str):
        raise TypeError(f"result must be a string, got {type(result).__name__}")

    values: list[complex] = []
    for line in result.splitlines():
        match = _COMPLEX_LINE.match(line)
        if match is None:
            continue
        real, imag = match.groups()
        values.append(complex(float(real), float(imag)))

    if not values:
        raise ValueError("No complex values found in the VPMR output.")
    if len(values) % 2 != 0:
        raise ValueError(
            f"Expected an even number of values (M and S blocks), got {len(values)}."
        )

    half = len(values) // 2
    return VPMRResult(
        np.asarray(values[:half], dtype=np.complex128),
        np.asarray(values[half:], dtype=np.complex128),
    )


def _evaluate_kernel(kernel: Callable, x: np.ndarray) -> np.ndarray:
    """Evaluate a user kernel, supporting both vectorized and scalar callables."""
    try:
        values = kernel(x)
    except (TypeError, ValueError):
        values = np.fromiter((kernel(item) for item in x), dtype=np.complex128, count=x.size)

    values = np.asarray(values)
    if values.ndim == 0:
        values = np.full(x.shape, values)
    else:
        values = np.broadcast_to(values, x.shape)
    return values


def plot(
    m: np.ndarray | list[complex] | VPMRResult,
    s: np.ndarray | list[complex],
    kernel: Callable,
    *,
    size: tuple[float, float] = (6, 4),
    xlim: tuple[float, float] = (0, 10),
    show: bool = True,
    save_to: str | Path | None = None,
    title: str | None = None,
) -> tuple:
    """Plot a kernel and its exponential approximation.

    Args:
        m: Coefficients of the approximation.
        s: Exponents of the approximation.
        kernel: Reference kernel callable.  It may accept a NumPy array or a
            scalar ``t``.
        size: Figure size in inches.
        xlim: Lower and upper limits of the time axis.
        show: Whether to call :func:`matplotlib.pyplot.show`.
        save_to: Optional path used to save the figure.
        title: Optional figure title.

    Returns:
        A ``(figure, axes, error_axes)`` tuple.  Matplotlib is imported lazily
        so that importing :mod:`pyvpmr` does not require a display backend.
    """
    import matplotlib.pyplot as plt

    m_array, s_array = _process_args((m, s))
    x = np.linspace(*xlim, 401)

    y_ref = _evaluate_kernel(kernel, x)
    y = evaluate(m_array, s_array, x)

    fig, ax = plt.subplots(figsize=size)
    ax.plot(x, np.real_if_close(y_ref), "b-", label="kernel", linewidth=2)
    ax.plot(x, np.real(y), "r", linestyle="dashdot", label="approximation", linewidth=3)
    ax.set_xlabel("time $t$ [s]")
    ax.set_ylabel("kernel function $g(t)$")
    ax.set_xlim(np.min(x), np.max(x))
    ax.legend(loc="upper right", handlelength=4)
    if title is not None:
        ax.set_title(title)

    error_ax = ax.twinx()
    # Log axes cannot display exact zeros; clip them to the smallest positive
    # normal number while keeping the visual representation meaningful.
    error = np.maximum(np.abs(y_ref - y), np.finfo(float).tiny)
    error_ax.plot(x, error, "g--", label="absolute error", linewidth=1)
    error_ax.set_yscale("log")
    error_ax.set_ylabel("absolute error")
    error_ax.legend(loc="center right")

    fig.tight_layout(pad=0.05)
    if save_to is not None:
        fig.savefig(Path(save_to))
    if show:
        plt.show()

    return fig, ax, error_ax


def _process_args(*args) -> tuple[np.ndarray, np.ndarray]:
    """Normalize the supported ``(m, s)`` calling conventions."""
    if len(args) == 1:
        value = args[0]
        if isinstance(value, VPMRResult):
            m, s = value.m, value.s
        else:
            try:
                m, s = value
            except (TypeError, ValueError) as exc:
                raise TypeError(
                    "Expected a VPMRResult or a pair of (m, s) arrays."
                ) from exc
    elif len(args) == 2:
        m, s = args
    else:
        raise TypeError(f"Expected 1 or 2 arguments, got {len(args)}.")

    m_array = np.asarray(m, dtype=np.complex128).reshape(-1)
    s_array = np.asarray(s, dtype=np.complex128).reshape(-1)
    if m_array.size != s_array.size:
        raise ValueError(
            f"m and s must have the same length, got {m_array.size} and {s_array.size}."
        )
    return m_array, s_array


def _format_terms(m: np.ndarray, s: np.ndarray) -> list[str]:
    return [
        f"{mi.real:+.15e} {mi.imag:+.15e} {si.real:+.15e} {si.imag:+.15e}"
        for mi, si in zip(m, s)
    ]


def to_global_damping(*args) -> str:
    """Generate a suanPan command for global nonviscous damping.

    Args:
        *args: Either ``(m, s)``, ``((m, s),)``, or a single
            :class:`VPMRResult`.

    Returns:
        A suanPan ``integrator NonviscousNewmark`` command.  The generated
        comments explain which values may need to be adjusted for a particular
        model.
    """
    m, s = _process_args(*args)
    command = (
        "# The following can be used as a global nonviscous damping with the Newmark time integration.\n"
        "# You may need to modify the first line to change tag and integration parameters.\n"
        "integrator NonviscousNewmark 1 .25 .5"
    )
    terms = _format_terms(m, s)
    if terms:
        command += " \\\n" + " \\\n".join(terms)
    return command + "\n"


def to_elemental_damping(*args) -> str:
    """Generate a suanPan command for per-element nonviscous damping.

    Args:
        *args: Either ``(m, s)``, ``((m, s),)``, or a single
            :class:`VPMRResult`.

    Returns:
        A suanPan ``modifier ElementalNonviscous`` command.  The generated
        comments explain how to apply the command to a multiplier element
        group instead.
    """
    m, s = _process_args(*args)
    command = (
        "# The following can be used as a per-element based nonviscous damping.\n"
        "# You may need to modify the first line to change tags.\n"
        "# Use the alternative form to apply to multiplier elements.\n"
        "# modifier ElementalNonviscousGroup {unique_modifier_tag} {associated_element_group_tag}\n"
        "modifier ElementalNonviscous {unique_modifier_tag} {associated_element_tag}"
    )
    terms = _format_terms(m, s)
    if terms:
        command += " \\\n" + " \\\n".join(terms)
    return command + "\n"
