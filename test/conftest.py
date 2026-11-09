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

"""Allow the pure-Python tests to run when the compiled extension is absent."""

from __future__ import annotations

import sys
import types

try:
    import _pyvpmr  # noqa: F401
except ImportError:
    module = types.ModuleType("_pyvpmr")

    def _fake_vpmr(**_: object):
        # A deterministic stand-in so API-level tests remain independent of
        # the optional build step.  Tests that need real numerics use the
        # HAS_EXTENSION flag in test_pyvpmr.py.
        return [complex(index + 1.0) for index in range(3)], [
            complex(index + 1.0) for index in range(3)
        ]

    module.vpmr = _fake_vpmr  # type: ignore[attr-defined]
    sys.modules["_pyvpmr"] = module
