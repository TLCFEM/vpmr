#  Copyright (C) 2023Theodore Chang
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

import re
from typing import Callable

import numpy as np
from _pyvpmr import vpmr as vpmr_impl
from matplotlib import pyplot as plt

vpmr = vpmr_impl


def split(result: str) -> tuple | None:
    """
    Split the output of the vpmr program into two arrays of complex numbers.

    :param result: The raw output of the vpmr program.
    :return: a tuple of two arrays of complex numbers.
    """
    split_r = result.strip().split('\n')
    regex = re.compile(r'([+\-]\d+\.\d+e[+\-]\d+){2}j')
    items = [i for i in split_r if regex.match(i)]
    if len(items) % 2 != 0:
        print('something wrong with the output')
        return None

    m_complex = [complex(i) for i in items[:len(items) // 2]]
    s_complex = [complex(i) for i in items[len(items) // 2:]]
    return np.array(m_complex), np.array(s_complex)


def plot(
        m: list | np.ndarray, s: list | np.ndarray, kernel: Callable, *,
        size: tuple = (6, 4),
        xlim: tuple = (0, 10),
        show: bool = True,
        save_to: str = None
):
    """
    Plot the kernel function and the approximation.

    :param m: The list of m values.
    :param s: The list of s values.
    :param kernel: The kernel function.
    :param size: The size of the figure.
    :param xlim: The x-axis limits.
    :param show: If to show the figure.
    :param save_to: Where to save the figure.
    """
    x = np.linspace(*xlim, 401)

    y_ref = np.zeros(len(x))
    for i in range(len(x)):
        y_ref[i] = kernel(x[i])

    y = np.zeros(len(x), dtype=complex)
    for ml, sl in zip(m, s):
        y += ml * np.exp(-sl * x)

    fig = plt.figure(figsize=size)

    ax1 = plt.gca()
    ax1.plot(x, y_ref, 'b-', label='kernel', linewidth=2)
    ax1.plot(x, y.real, 'r', linestyle='dashdot', label='approximation', linewidth=3)
    ax1.set_xlabel('time $t$ [s]')
    ax1.set_ylabel('kernel function $g(t)$')
    ax1.legend(loc='upper right', handlelength=4)

    ax2 = plt.twinx()
    ax2.plot(x, np.abs(y_ref - y), 'g--', label='absolute error', linewidth=1)
    ax2.set_yscale('log')
    ax2.set_ylabel('absolute error')
    ax2.legend(loc='center right')

    plt.xlim(np.min(x), np.max(x))
    plt.tight_layout(pad=.05)
    if show:
        plt.show()
    if save_to:
        fig.savefig(save_to)
