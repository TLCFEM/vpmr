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

from __future__ import annotations

import math

import numpy as np
import pytest

import _pyvpmr
from pyvpmr import (
    DEFAULT_KERNEL,
    VPMRResult,
    evaluate,
    plot,
    split,
    to_elemental_damping,
    to_global_damping,
    vpmr,
)

HAS_EXTENSION = getattr(_pyvpmr, "__file__", None) is not None
requires_extension = pytest.mark.skipif(
    not HAS_EXTENSION, reason="the compiled _pyvpmr extension is not available"
)


def test_default_kernel_is_available():
    assert DEFAULT_KERNEL == "exp(-t*t/4)"


def test_split_parses_standalone_output():
    output = """
Using the following parameters:
       terms = 2.

M =
+2.3132817597168739e+01-2.8586221856566439e-105j
-1.0000000000000000e+00+2.0000000000000000e+00j
S =
+0.0000000000000000e+00-0.0000000000000000e+00j
+3.1400000000000000e+00+0.0000000000000000e+00j
"""
    result = split(output)

    assert isinstance(result, VPMRResult)
    np.testing.assert_allclose(result.m, [23.132817597168739 - 2.8586221856566439e-105j, -1 + 2j])
    np.testing.assert_allclose(result.s, [0 + 0j, 3.14 + 0j])


def test_split_accepts_compact_numbers():
    result = split("M =\n-1+2j\nS =\n0-3.5j\n")
    np.testing.assert_array_equal(result.m, [-1 + 2j])
    np.testing.assert_array_equal(result.s, [-3.5j])


@pytest.mark.parametrize("output", ["", "no numbers here", "M =\n+1.0+2.0j\n"])
def test_split_rejects_incomplete_output(output):
    with pytest.raises(ValueError):
        split(output)


def test_split_rejects_non_string():
    with pytest.raises(TypeError):
        split(None)  # type: ignore[arg-type]


def test_evaluate_scalar_and_array():
    m = np.array([1.0, 2.0])
    s = np.array([0.0, 1.0])

    assert evaluate(m, s, 0.0) == pytest.approx(3.0)
    np.testing.assert_allclose(
        evaluate(m, s, np.array([0.0, 1.0])),
        [3.0, 1.0 + 2.0 / math.e],
    )


def test_evaluate_accepts_vpmr_result():
    result = VPMRResult(np.array([1 + 1j]), np.array([0 + 0j]))
    assert evaluate(result.m, result.s, 2.0) == pytest.approx(1 + 1j)


def test_chained_helpers_accept_result_objects():
    result = VPMRResult(np.array([2.0]), np.array([3.0]))
    command = to_global_damping(result)
    assert "2.000000000000000e+00" in command
    assert "3.000000000000000e+00" in command


def test_to_global_damping_formats_all_terms():
    m = np.array([1 + 2j, -3 - 4j])
    s = np.array([0.5 + 0j, 1.5 - 0.25j])

    command = to_global_damping(m, s)
    lines = command.splitlines()

    assert lines[0].startswith("#")
    assert lines[2].startswith("integrator NonviscousNewmark")
    assert command.count(" \\\n") == len(m)
    assert any(
        line.startswith("+1.000000000000000e+00 +2.000000000000000e+00") for line in lines
    )
    assert command.endswith("\n")


def test_to_elemental_damping_keeps_commands_uncommented():
    command = to_elemental_damping([1 + 0j], [2 + 0j])
    lines = command.splitlines()

    assert lines[3] == (
        "# modifier ElementalNonviscousGroup {unique_modifier_tag} "
        "{associated_element_group_tag}"
    )
    assert lines[4].startswith("modifier ElementalNonviscous")
    assert not lines[4].startswith("#")


def test_damping_helpers_validate_shapes():
    with pytest.raises(ValueError):
        to_global_damping([1 + 0j], [1 + 0j, 2 + 0j])
    with pytest.raises(TypeError):
        to_global_damping(1, 2, 3)  # type: ignore[call-arg]


def test_plot_writes_a_figure(tmp_path):
    import matplotlib

    matplotlib.use("Agg")

    m = np.array([1.0])
    s = np.array([1.0])
    file = tmp_path / "approximation.png"

    fig, ax, error_ax = plot(
        m,
        s,
        lambda t: np.exp(-t),
        xlim=(0, 2),
        show=False,
        save_to=file,
        title="test",
    )

    assert file.exists()
    assert ax.get_title() == "test"
    assert error_ax.get_yscale() == "log"

    import matplotlib.pyplot as plt

    plt.close(fig)


def test_plot_supports_scalar_only_kernels(tmp_path):
    fig, _, _ = plot(
        [1.0],
        [1.0],
        lambda t: math.exp(-t),
        xlim=(0, 1),
        show=False,
        save_to=tmp_path / "scalar-kernel.png",
    )

    import matplotlib.pyplot as plt

    plt.close(fig)


@pytest.mark.parametrize(
    ("argument", "value"),
    [
        ("n", 0),
        ("c", 0),
        ("d", -1),
        ("q", 0),
        ("m", 1.0),
        ("e", 0.0),
    ],
)
def test_vpmr_validates_numeric_ranges(argument, value):
    with pytest.raises(ValueError):
        vpmr(**{argument: value})


def test_vpmr_validates_types():
    with pytest.raises(TypeError):
        vpmr(n=2.5)
    with pytest.raises(TypeError):
        vpmr(k=1)  # type: ignore[arg-type]


def test_vpmr_returns_arrays():
    result = vpmr(n=8, q=100, verbose=False)

    assert isinstance(result, VPMRResult)
    assert isinstance(result.m, np.ndarray)
    assert isinstance(result.s, np.ndarray)
    assert result.m.dtype == np.complex128
    assert result.s.dtype == np.complex128
    assert result.m.shape == result.s.shape
    assert result.m.size > 0


@requires_extension
def test_vpmr_repeated_calls_do_not_reuse_kernel_cache():
    # Regression test: the first call used to leave both the kernel expression
    # and the sampled-kernel cache in place for later calls.
    first = vpmr(n=10, q=120, e=1e-6, k="exp(-t*t/4)", verbose=False)
    second = vpmr(n=10, q=120, e=1e-6, k="exp(-t*t/10)", verbose=False)
    third = vpmr(n=10, q=120, e=1e-6, verbose=False)

    assert first.m.size > 0 and second.m.size > 0 and third.m.size > 0
    assert first.m.shape == third.m.shape
    np.testing.assert_allclose(first.m, third.m)
    assert first.m.shape != second.m.shape or not np.allclose(first.m, second.m)


@requires_extension
def test_vpmr_omission_keeps_terms_above_tolerance():
    # Regression test: trivial-term detection used the squared magnitude, so
    # terms around 3e-5 were discarded against a 1e-8 tolerance and the n=20
    # approximation was only accurate to about 5e-5.
    result = vpmr(n=20, e=1e-8, omit=True, verbose=False)
    t = np.linspace(0.0, 10.0, 401)
    error = np.abs(np.exp(-(t**2) / 4.0) - evaluate(result.m, result.s, t))

    assert error.max() < 1e-6


@requires_extension
def test_vpmr_invalid_kernel_raises_runtime_error():
    with pytest.raises(RuntimeError):
        vpmr(n=4, q=40, k="this is not an expression", verbose=False)
