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

import numpy as np
from matplotlib import pyplot as plt

weights = '''
M = 
+1.0589202279681926e-11+0.0000000000000000e+00j
-5.4905134221689723e-03+2.2104939243740062e-03j
-5.4905134221689723e-03-2.2104939243740062e-03j
+1.1745193571738945e+01+2.3424418474362436e-153j
-5.5143304351134406e+00-5.7204056791636839e+00j
-5.5143304351134406e+00+5.7204056791636839e+00j
-1.6161617424833765e-02+2.3459542440459513e+00j
-1.6161617424833765e-02-2.3459542440459513e+00j
+1.6338578576177490e-01-1.9308431539218421e-01j
+1.6338578576177490e-01+1.9308431539218421e-01j
S = 
+0.0000000000000000e+00+0.0000000000000000e+00j
+1.7655956664692953e+00-2.7555720406099038e+00j
+1.7655956664692953e+00+2.7555720406099038e+00j
+1.8757961592204051e+00-0.0000000000000000e+00j
+1.8700580506914817e+00-6.2013413918954552e-01j
+1.8700580506914817e+00+6.2013413918954552e-01j
+1.8521958553280000e+00-1.2601975249082220e+00j
+1.8521958553280000e+00+1.2601975249082220e+00j
+1.8197653300065935e+00-1.9494562062795735e+00j
+1.8197653300065935e+00+1.9494562062795735e+00j
'''


def kernel(x):
    return np.exp(-x * x / 4)


def split(r):
    split_r = r.strip().split('\n')
    m = split_r[:len(split_r) // 2]
    s = split_r[len(split_r) // 2:]
    m_complex = []
    s_complex = []
    for i in m:
        try:
            m_complex.append(complex(i))
        except ValueError:
            pass
    for i in s:
        try:
            s_complex.append(complex(i))
        except ValueError:
            pass
    return np.array(m_complex), np.array(s_complex)


def plotter():
    x = np.linspace(0, 5, 401)
    yy = kernel(x)
    y = np.zeros(len(x), dtype=complex)
    for ml, sl in zip(*split(weights)):
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


if __name__ == '__main__':
    plotter()
