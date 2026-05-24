from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Callable

import numpy as np

try:
    import _pyvpmr
except ImportError as exc:  # pragma: no cover - only exercised when extension is unavailable.
    _pyvpmr = None
    _PYVPMR_IMPORT_ERROR = exc
else:
    _PYVPMR_IMPORT_ERROR = None

__all__ = [
    "VPMROptions",
    "VPMRResult",
    "plot",
    "split",
    "to_elemental_damping",
    "to_global_damping",
    "vpmr",
]


def _require_backend() -> None:
    if _pyvpmr is None:
        raise ImportError(
            "The compiled _pyvpmr extension is unavailable. Install pyvpmr from a wheel "
            "or build the extension before calling vpmr()."
        ) from _PYVPMR_IMPORT_ERROR


def _ensure_positive_int(name: str, value: int) -> int:
    if not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer.")
    return value


def _ensure_non_negative_int(name: str, value: int) -> int:
    if not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a non-negative integer.")
    return value


def _ensure_positive_float(name: str, value: float) -> float:
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a positive finite number.")
    return value


@dataclass(frozen=True)
class VPMROptions:
    """User-facing options for the VPMR algorithm."""

    terms: int = 20
    max_exponent: int = 4
    precision_bits: int = 0
    quadrature_order: int = 500
    precision_multiplier: float = 1.05
    tolerance: float = 1e-8
    kernel: str = ""
    omit_trivial_terms: bool = True

    def __post_init__(self) -> None:
        object.__setattr__(self, "terms", _ensure_positive_int("terms", self.terms))
        object.__setattr__(
            self, "max_exponent", _ensure_positive_int("max_exponent", self.max_exponent)
        )
        object.__setattr__(
            self,
            "precision_bits",
            _ensure_non_negative_int("precision_bits", self.precision_bits),
        )
        object.__setattr__(
            self,
            "quadrature_order",
            _ensure_positive_int("quadrature_order", self.quadrature_order),
        )
        object.__setattr__(
            self,
            "precision_multiplier",
            _ensure_positive_float("precision_multiplier", self.precision_multiplier),
        )
        object.__setattr__(self, "tolerance", _ensure_positive_float("tolerance", self.tolerance))
        if not isinstance(self.kernel, str):
            raise ValueError("kernel must be a string expression.")
        if not isinstance(self.omit_trivial_terms, bool):
            raise ValueError("omit_trivial_terms must be a boolean.")


@dataclass(frozen=True)
class VPMRResult:
    """Result of the VPMR algorithm."""

    weights: np.ndarray
    poles: np.ndarray

    def __post_init__(self) -> None:
        weights = np.atleast_1d(np.asarray(self.weights, dtype=complex))
        poles = np.atleast_1d(np.asarray(self.poles, dtype=complex))
        if weights.shape != poles.shape:
            raise ValueError("weights and poles must have the same shape.")
        object.__setattr__(self, "weights", weights)
        object.__setattr__(self, "poles", poles)

    def evaluate(self, x: float | np.ndarray) -> complex | np.ndarray:
        """Evaluate the approximation at one or more time points."""
        sample = np.asarray(x, dtype=float)
        if sample.ndim == 0:
            return np.sum(self.weights * np.exp(-self.poles * sample)).item()
        weights = self.weights.reshape((-1,) + (1,) * sample.ndim)
        poles = self.poles.reshape((-1,) + (1,) * sample.ndim)
        return np.sum(weights * np.exp(-poles * sample), axis=0)

    def plot(self, kernel: Callable, **kwargs):
        """Plot the kernel and approximation."""
        return plot(self, kernel, **kwargs)

    def to_global_damping(self) -> str:
        """Format as a suanPan global nonviscous damping command."""
        return to_global_damping(self)

    def to_elemental_damping(self) -> str:
        """Format as a suanPan elemental nonviscous damping command."""
        return to_elemental_damping(self)


def vpmr(options: VPMROptions | None = None, /, **kwargs) -> VPMRResult:
    """
    Run the VPMR algorithm.

    Provide either an options object or keyword arguments matching VPMROptions.
    """
    _require_backend()
    if options is not None and kwargs:
        raise ValueError("Provide either 'options' or keyword arguments, not both.")
    if options is None:
        options = VPMROptions(**kwargs)
    elif not isinstance(options, VPMROptions):
        raise TypeError("options must be a VPMROptions instance.")

    raw_weights, raw_poles = _pyvpmr.vpmr(
        n=options.terms,
        c=options.max_exponent,
        d=options.precision_bits,
        q=options.quadrature_order,
        m=options.precision_multiplier,
        e=options.tolerance,
        k=options.kernel,
        omit=options.omit_trivial_terms,
    )
    return VPMRResult(raw_weights, raw_poles)


def split(result: str) -> VPMRResult:
    """Parse standalone VPMR output text into a VPMRResult."""
    split_r = result.strip().split("\n")
    regex = re.compile(r"([+\-]\d+\.\d+e[+\-]\d+){2}j")
    items = [i for i in split_r if regex.match(i)]
    if len(items) % 2 != 0:
        raise ValueError("Cannot parse standalone output into paired weights and poles.")

    half = len(items) // 2
    weights = [complex(i) for i in items[:half]]
    poles = [complex(i) for i in items[half:]]
    return VPMRResult(weights, poles)


def _evaluate_kernel(kernel: Callable, x: np.ndarray) -> np.ndarray:
    try:
        y_ref = np.asarray(kernel(x), dtype=complex)
        if y_ref.shape == x.shape:
            return y_ref
    except (TypeError, ValueError):
        pass

    try:
        return np.asarray(np.vectorize(kernel, otypes=[complex])(x), dtype=complex)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "The kernel function must accept either a NumPy array or scalar float values."
        ) from exc


def plot(
    result: VPMRResult,
    kernel: Callable,
    *,
    size: tuple[float, float] = (6, 4),
    xlim: tuple[float, float] = (0, 10),
    show: bool = True,
    save_to: str | None = None,
):
    """
    Plot kernel reference and approximation for a VPMRResult.

    The rendered lines use the real parts of the evaluated values.
    """
    if not isinstance(result, VPMRResult):
        raise TypeError("plot expects a VPMRResult instance as the first argument.")

    from matplotlib import pyplot as plt

    x = np.linspace(*xlim, 401)
    y_ref = _evaluate_kernel(kernel, x)
    y = result.evaluate(x)

    fig = plt.figure(figsize=size)
    ax1 = plt.gca()
    ax1.plot(x, y_ref.real, "b-", label="kernel", linewidth=2)
    ax1.plot(x, y.real, "r", linestyle="dashdot", label="approximation", linewidth=3)
    ax1.set_xlabel("time $t$ [s]")
    ax1.set_ylabel("kernel function $g(t)$")
    ax1.legend(loc="upper right", handlelength=4)

    ax2 = plt.twinx()
    ax2.plot(x, np.abs(y_ref - y), "g--", label="absolute error", linewidth=1)
    ax2.set_yscale("log")
    ax2.set_ylabel("absolute error")
    ax2.legend(loc="center right")

    plt.xlim(np.min(x), np.max(x))
    plt.tight_layout(pad=0.05)
    if show:
        plt.show()
    if save_to:
        fig.savefig(save_to)
    return fig, (ax1, ax2)


def to_global_damping(result: VPMRResult) -> str:
    """Generate a suanPan global nonviscous damping command."""
    if not isinstance(result, VPMRResult):
        raise TypeError("to_global_damping expects a VPMRResult instance.")

    command = "# The following can be used as a global nonviscous damping with the Newmark time integration.\n"
    command += "# You may need to modify the first line to change tag and integration parameters.\n"
    command += "integrator NonviscousNewmark 1 .25 .5"

    for weight, pole in zip(result.weights, result.poles):
        command += (
            f" \\\n{weight.real:+.15e} {weight.imag:+.15e} "
            f"{pole.real:+.15e} {pole.imag:+.15e}"
        )
    command += "\n"
    return command


def to_elemental_damping(result: VPMRResult) -> str:
    """Generate a suanPan elemental nonviscous damping command."""
    if not isinstance(result, VPMRResult):
        raise TypeError("to_elemental_damping expects a VPMRResult instance.")

    command = "# The following can be used as a per-element based nonviscous damping.\n"
    command += "# You may need to modify the first line to change tags.\n"
    command += "# Use the alternative form to apply to multiplier elements.\n"
    command += "# modifier ElementalNonviscousGroup {unique_modifier_tag} {associated_element_group_tag}\n"
    command += "modifier ElementalNonviscous {unique_modifier_tag} {associated_element_tag}"

    for weight, pole in zip(result.weights, result.poles):
        command += (
            f" \\\n{weight.real:+.15e} {weight.imag:+.15e} "
            f"{pole.real:+.15e} {pole.imag:+.15e}"
        )
    command += "\n"
    return command
