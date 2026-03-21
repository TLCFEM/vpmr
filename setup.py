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

import os
import platform
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup


def initialize():
    include_dirs = [
        "eigen",
        "eigen/unsupported",
        "mpreal",
        "exprtk",
        "exprtk-custom-types",
    ]

    library_dirs = []

    def add_to_include(target: str):
        if not os.path.exists(target):
            return

        for root, dirs, _ in os.walk(target):
            if "include" in dirs:
                include_dirs.append(root + "/include")
                break

    def add_to_library(target: str, header: str):
        if not os.path.exists(target):
            return

        for root, _, files in os.walk(target):
            if any(header in file for file in files):
                library_dirs.append(root)
                break

    if "Darwin" == platform.system():
        os.environ["CC"] = "gcc-13"
        os.environ["CXX"] = "g++-13"

    for path in ["/opt/homebrew/Cellar", "/usr/local/Cellar"]:
        add_to_include(f"{path}/mpfr")
        add_to_include(f"{path}/gmp")
        add_to_include(f"{path}/tbb")

        add_to_library(f"{path}/mpfr", "libmpfr")
        add_to_library(f"{path}/gmp", "libgmp")
        add_to_library(f"{path}/tbb", "libtbb")

    libraries = ["mpfr", "gmp", "tbb"]
    extra = []

    # for path in [
    #     "/usr/lib/libmpfr.a",
    #     "/usr/lib/x86_64-linux-gnu/libmpfr.a",
    # ]:
    #     if Path(path).exists():
    #         extra.append(path)
    #         break
    # else:
    #     libraries.append("mpfr")
    #
    # for path in [
    #     "/usr/lib/x86_64-linux-gnu/libgmp.a",
    #     "/usr/lib64/libgmp.a",
    #     "/usr/lib/libgmp.a",
    # ]:
    #     if Path(path).exists():
    #         extra.append(path)
    #         break
    # else:
    #     libraries.append("gmp")

    return Pybind11Extension(
        "_pyvpmr",
        ["src/VPMR.cpp"],
        include_dirs=include_dirs,
        library_dirs=library_dirs,
        libraries=libraries,
        extra_objects=extra,
        define_macros=[("PYVPMR", 1)],
        extra_compile_args=["-fPIC"],
    )


# noinspection PyTypeChecker
setup(
    long_description=(Path(__file__).parent / "README.md")
    .read_text(encoding="utf8")
    .replace(
        "resource/", "https://raw.githubusercontent.com/TLCFEM/vpmr/master/resource/"
    ),
    long_description_content_type="text/markdown",
    ext_modules=[
        initialize(),
    ],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    packages=["pyvpmr"],
)
