import os
import platform
from datetime import datetime
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

if "Darwin" == platform.system():
    os.environ["CC"] = "gcc-13"
    os.environ["CXX"] = "g++-13"

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


for root in ["/opt/homebrew/Cellar", "/usr/local/Cellar"]:
    add_to_include(f"{root}/mpfr")
    add_to_include(f"{root}/gmp")
    add_to_include(f"{root}/tbb")

    add_to_library(f"{root}/mpfr", "libmpfr")
    add_to_library(f"{root}/gmp", "libgmp")
    add_to_library(f"{root}/tbb", "libtbb")

# noinspection PyTypeChecker
setup(
    version=datetime.now().strftime("%y%m%d"),
    long_description=(Path(__file__).parent / "README.md")
    .read_text(encoding="utf8")
    .replace(
        "resource/", "https://raw.githubusercontent.com/TLCFEM/vpmr/master/resource/"
    ),
    long_description_content_type="text/markdown",
    ext_modules=[
        Pybind11Extension(
            "_pyvpmr",
            ["src/VPMR.cpp"],
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            libraries=["mpfr", "gmp", "tbb"],
            define_macros=[("PYVPMR", 1)],
        ),
    ],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    packages=["pyvpmr"],
)
