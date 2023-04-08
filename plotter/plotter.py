################################################################################
#  Copyright (C) 2023 Theodore Chang
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
################################################################################
import os
import re
import subprocess
from os.path import exists
from shutil import which

import click
import numpy as np
from matplotlib import pyplot as plt


# change this kernel before plotting
def kernel(x):
    return np.exp(-x ** 2 / 4)


def split(r: str):
    split_r = r.strip().split('\n')
    regex = re.compile(r'([+\-]\d+\.\d+e[+\-]\d+){2}j')
    items = [i for i in split_r if regex.match(i)]
    if len(items) % 2 != 0:
        print('something wrong with the output')
        return None

    m_complex = [complex(i) for i in items[:len(items) // 2]]
    s_complex = [complex(i) for i in items[len(items) // 2:]]
    return np.array(m_complex), np.array(s_complex)


def plotter(output: str, save: str):
    if (result := split(output)) is None:
        return
    x = np.linspace(0, 10, 401)
    yy = np.zeros(len(x))
    for i in range(len(x)):
        yy[i] = kernel(x[i])
    y = np.zeros(len(x), dtype=complex)
    for ml, sl in zip(*result):
        y += ml * np.exp(-sl * x)

    fig = plt.figure(figsize=(6, 4))

    ax1 = plt.gca()
    ax1.plot(x, yy, 'b-', label='kernel', linewidth=2)
    ax1.plot(x, y.real, 'r', linestyle='dashdot', label='approximation', linewidth=3)
    ax1.set_xlabel('time $t$')
    ax1.set_ylabel('kernel function $g(t)$')
    ax1.legend(loc='upper right', handlelength=4)

    ax2 = plt.twinx()
    ax2.plot(x, np.abs(yy - y), 'g--', label='absolute error', linewidth=1)
    ax2.set_yscale('log')
    ax2.set_ylabel('absolute error')
    ax2.legend(loc='center right')

    plt.xlim(np.min(x), np.max(x))
    plt.tight_layout(pad=.05)
    plt.show()
    if save:
        fig.savefig(save)


@click.command()
@click.option('-n', default=30, help='number of terms')
@click.option('-q', default=500, help='quadrature order')
@click.option('-d', default=100, help='number of precision digits')
@click.option('-nc', default=4, help='controls the maximum exponents')
@click.option('-e', default=1e-8, help='tolerance')
@click.option('-o', default='', help='save plot to file')
@click.option('-k', default='', help='file name of kernel function')
@click.option('--exe', default='', help='path to vpmr executable')
def execute(n, nc, d, q, e, k, exe, o):
    exe = exe if exe else 'vpmr'

    if exists(exe):
        exe = os.path.abspath(exe)
    elif which(exe) is None:
        print('cannot find the vpmr executable')
        return

    command = [exe, '-n', str(n), '-nc', str(nc), '-d', str(d), '-q', str(q), '-e', str(e)]
    if k:
        command.extend(['-k', k])
    result = subprocess.check_output(command).decode('utf-8')
    print(result)
    plotter(result, o)


if __name__ == '__main__':
    execute()
