from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Callable

import numpy as np

try:
    import _pyvpmr
except ImportError as exc:  # pragma: no cover - exercised in environments without the extension.
    _pyvpmr = None
    _PYVPMR_IMPORT_ERROR = exc
else:
    _PYVPMR_IMPORT_ERROR = None

__all__ = [
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


def _coalesce_argument(
    canonical,
    alias,
    *,
    canonical_name: str,
    alias_name: str,
    default,
):
    if canonical is not None and alias is not None and canonical != alias:
        raise ValueError(
            f"{canonical_name!r} and {alias_name!r} refer to the same option and must match."
        )
    if canonical is not None:
        return canonical
    if alias is not None:
        return alias
    return default


@dataclass(frozen=True)
class VPMRResult:
    """
    Container returned by :func:`vpmr`.

    The result can still be unpacked as ``m, s = vpmr(...)`` for backwards
    compatibility, while also exposing named attributes and convenience methods.
    """

    m: np.ndarray
    s: np.ndarray

    def __post_init__(self) -> None:
        m = np.atleast_1d(np.asarray(self.m, dtype=complex))
        s = np.atleast_1d(np.asarray(self.s, dtype=complex))
        if m.shape != s.shape:
            raise ValueError("The length of m and s must be the same.")
        object.__setattr__(self, "m", m)
        object.__setattr__(self, "s", s)

    @property
    def weights(self) -> np.ndarray:
        """Alias for ``m``."""
        return self.m

    @property
    def poles(self) -> np.ndarray:
        """Alias for ``s``."""
        return self.s

    def as_tuple(self) -> tuple[np.ndarray, np.ndarray]:
        """Return the legacy ``(m, s)`` tuple."""
        return self.m, self.s

    def evaluate(self, x: float | np.ndarray) -> complex | np.ndarray:
        """Evaluate the exponential approximation at one or more points."""
        sample = np.asarray(x, dtype=float)
        weights = self.m.reshape((-1,) + (1,) * sample.ndim)
        poles = self.s.reshape((-1,) + (1,) * sample.ndim)
        values = np.sum(weights * np.exp(-poles * sample), axis=0)
        return values.item() if sample.ndim == 0 else values

    def plot(self, kernel: Callable, **kwargs):
        """Plot the reference kernel and the approximation."""
        return plot(self, kernel, **kwargs)

    def to_global_damping(self) -> str:
        """Format the approximation as a suanPan global damping command."""
        return to_global_damping(self)

    def to_elemental_damping(self) -> str:
        """Format the approximation as a suanPan elemental damping command."""
        return to_elemental_damping(self)

    def __iter__(self):
        yield self.m
        yield self.s

    def __len__(self) -> int:
        return 2

    def __getitem__(self, index: int) -> np.ndarray:
        return self.as_tuple()[index]


def vpmr(
    *,
    n: int | None = None,
    terms: int | None = None,
    c: int | None = None,
    max_exponent: int | None = None,
    d: int | None = None,
    precision_bits: int | None = None,
    q: int | None = None,
    quadrature_order: int | None = None,
    m: float | None = None,
    precision_multiplier: float | None = None,
    e: float | None = None,
    tolerance: float | None = None,
    k: str | None = None,
    kernel: str | None = None,
    omit: bool | None = None,
    omit_trivial_terms: bool | None = None,
) -> VPMRResult:
    """
    Run the VPMR algorithm.

    Short option names match the C++ binding. Descriptive aliases are also
    accepted for better readability:

    ``terms`` → ``n``, ``max_exponent`` → ``c``, ``precision_bits`` → ``d``,
    ``quadrature_order`` → ``q``, ``precision_multiplier`` → ``m``,
    ``tolerance`` → ``e``, ``kernel`` → ``k`` and ``omit_trivial_terms`` → ``omit``.
    """
    _require_backend()

    options = {
        "n": _coalesce_argument(
            n, terms, canonical_name="n", alias_name="terms", default=20
        ),
        "c": _coalesce_argument(
            c,
            max_exponent,
            canonical_name="c",
            alias_name="max_exponent",
            default=4,
        ),
        "d": _coalesce_argument(
            d,
            precision_bits,
            canonical_name="d",
            alias_name="precision_bits",
            default=0,
        ),
        "q": _coalesce_argument(
            q,
            quadrature_order,
            canonical_name="q",
            alias_name="quadrature_order",
            default=500,
        ),
        "m": _coalesce_argument(
            m,
            precision_multiplier,
            canonical_name="m",
            alias_name="precision_multiplier",
            default=1.05,
        ),
        "e": _coalesce_argument(
            e,
            tolerance,
            canonical_name="e",
            alias_name="tolerance",
            default=1e-8,
        ),
        "k": _coalesce_argument(
            k,
            kernel,
            canonical_name="k",
            alias_name="kernel",
            default="",
        ),
        "omit": _coalesce_argument(
            omit,
            omit_trivial_terms,
            canonical_name="omit",
            alias_name="omit_trivial_terms",
            default=True,
        ),
    }

    return VPMRResult(*_pyvpmr.vpmr(**options))


def split(result: str) -> VPMRResult | None:
    """
    Split the standalone VPMR output into a :class:`VPMRResult`.

    :param result: The raw output of the vpmr program.
    :return: Parsed VPMR result.
    """
    split_r = result.strip().split("\n")
    regex = re.compile(r"([+\-]\d+\.\d+e[+\-]\d+){2}j")
    items = [i for i in split_r if regex.match(i)]
    if len(items) % 2 != 0:
        print("something wrong with the output")
        return None

    m_complex = [complex(i) for i in items[: len(items) // 2]]
    s_complex = [complex(i) for i in items[len(items) // 2 :]]
    return VPMRResult(m_complex, s_complex)


def _evaluate_kernel(kernel: Callable, x: np.ndarray) -> np.ndarray:
    try:
        y_ref = np.asarray(kernel(x), dtype=complex)
        if y_ref.shape == x.shape:
            return y_ref
    except (TypeError, ValueError):
        pass

    try:
        return np.asarray([kernel(float(item)) for item in x], dtype=complex)
    except Exception as exc:
        raise ValueError(
            "The kernel function must accept either a NumPy array or scalar float values."
        ) from exc


def plot(
    m: VPMRResult | list | np.ndarray,
    s: list | np.ndarray | Callable | None = None,
    kernel: Callable | None = None,
    *,
    size: tuple[float, float] = (6, 4),
    xlim: tuple[float, float] = (0, 10),
    show: bool = True,
    save_to: str | None = None,
):
    """
    Plot the kernel function and the approximation.

    Accepts either ``plot(result, kernel, ...)`` or ``plot(m, s, kernel, ...)``.
    Returns the created figure and axes tuple.
    """
    if isinstance(m, VPMRResult):
        result = m
        if kernel is None and callable(s):
            kernel = s
            s = None
        elif s is not None:
            raise ValueError(
                "When the first argument is a VPMRResult, provide the kernel as the "
                "second positional argument or as the 'kernel' keyword argument."
            )
    else:
        if s is None or callable(s):
            raise ValueError("Both m and s must be provided unless a VPMRResult is used.")
        result = VPMRResult(*_process_args(m, s))

    if kernel is None:
        raise ValueError("A kernel function is required.")

    from matplotlib import pyplot as plt

    x = np.linspace(*xlim, 401)
    y_ref = _evaluate_kernel(kernel, x)
    y = result.evaluate(x)

    fig = plt.figure(figsize=size)

    ax1 = plt.gca()
    ax1.plot(x, y_ref.real, "b-", label="kernel", linewidth=2)
    ax1.plot(
        x,
        y.real,
        "r",
        linestyle="dashdot",
        label="approximation",
        linewidth=3,
    )
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


def _process_args(*args) -> tuple[np.ndarray, np.ndarray]:
    if len(args) == 1:
        if isinstance(args[0], VPMRResult):
            return args[0].as_tuple()
        if 2 != len(args[0]):
            raise ValueError("Expected a two-item result tuple.")
        m, s = args[0]
    elif len(args) == 2:
        m, s = args
    else:
        raise ValueError("Wrong number of arguments.")

    result = VPMRResult(m, s)
    return result.as_tuple()


def to_global_damping(*args):
    """
    Generate a command to use the kernel as a global nonviscous damping model in suanPan.

    :param args: A result object, a ``(m, s)`` tuple, or separate ``m`` and ``s`` arrays.
    :return: The command.
    """
    command = "# The following can be used as a global nonviscous damping with the Newmark time integration.\n"
    command += "# You may need to modify the first line to change tag and integration parameters.\n"
    command += "integrator NonviscousNewmark 1 .25 .5"

    for m, s in zip(*_process_args(*args)):
        command += f" \\\n{m.real:+.15e} {m.imag:+.15e} {s.real:+.15e} {s.imag:+.15e}"

    command += "\n"

    return command


def to_elemental_damping(*args):
    """
    Generate a command to use the kernel as a per-element nonviscous damping model in suanPan.

    :param args: A result object, a ``(m, s)`` tuple, or separate ``m`` and ``s`` arrays.
    :return: The command.
    """
    command = "# The following can be used as a per-element based nonviscous damping.\n"
    command += "# You may need to modify the first line to change tags.\n"
    command += "# Use the alternative form to apply to multiplier elements.\n"
    command += "# modifier ElementalNonviscousGroup {unique_modifier_tag} {associated_element_group_tag}"
    command += (
        "modifier ElementalNonviscous {unique_modifier_tag} {associated_element_tag}"
    )

    for m, s in zip(*_process_args(*args)):
        command += f" \\\n{m.real:+.15e} {m.imag:+.15e} {s.real:+.15e} {s.imag:+.15e}"

    command += "\n"

    return command
