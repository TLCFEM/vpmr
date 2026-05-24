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

import numpy as np

from pyvpmr import vpmr


def kernel(x):
    return np.exp(-(x**2) / 4) * x**3


if __name__ == "__main__":
    result = vpmr(
        terms=60,
        kernel="exp(-t^2/4)*t^3",
        tolerance=1e-12,
        max_exponent=10,
    )
    print(result.to_global_damping())
    result.plot(kernel)
