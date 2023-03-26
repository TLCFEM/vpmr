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
    return np.exp(-x * x / 4)


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


def plotter(output: str):
    if (result := split(output)) is None:
        return
    x = np.linspace(0, 5, 401)
    yy = kernel(x)
    y = np.zeros(len(x), dtype=complex)
    for ml, sl in zip(*result):
        y += ml * np.exp(-sl * x)

    fig = plt.figure(figsize=(6, 4))
    plt.plot(x, yy, 'b-', label='kernel', linewidth=2)
    plt.plot(x, y.real, 'r--', label='approximation', linewidth=3)

    plt.xlabel('t')
    plt.ylabel('f(t)')
    plt.legend(loc='lower left')

    ax2 = plt.twinx()
    ax2.plot(x, np.abs(yy - y), 'g-', label='Error', linewidth=1)
    ax2.set_yscale('log')
    ax2.set_ylabel('Absolute Error')

    plt.title('Approximation of the kernel via VPMR')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.show()


@click.command()
@click.option('-n', default=30, help='number of terms')
@click.option('-nc', default=4, help='number of exponents')
@click.option('-d', default=512, help='number of precision digits')
@click.option('-q', default=500, help='quadrature order')
@click.option('-e', default=1e-8, help='tolerance')
@click.option('-k', default=None, help='file name of kernel function')
@click.option('-exe', default=None, help='path of vpmr executable')
def execute(n, nc, d, q, e, k, exe):
    if not exe:
        exe = 'vpmr'

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
    plotter(result)


if __name__ == '__main__':
    execute()
    # make a command line tool using click
