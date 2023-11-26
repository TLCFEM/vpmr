import platform
from datetime import datetime
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# noinspection PyTypeChecker
setup(
    name="pyvpmr",
    version=datetime.now().strftime("%y%m%d"),
    author="Theodore Chang",
    author_email="tlcfem@gmail.com",
    url="https://github.com/TLCFEM/vpmr",
    description="The VPMR Algorithm",
    long_description=(Path(__file__).parent / "README.md").read_text().replace(
        "resource/", "https://raw.githubusercontent.com/TLCFEM/vpmr/master/resource/"),
    long_description_content_type='text/markdown',
    install_requires=["numpy", "matplotlib"],
    ext_modules=[
        Pybind11Extension(
            "_pyvpmr",
            ["src/VPMR.cpp"],
            include_dirs=["eigen", "eigen/unsupported", "mpreal", "exprtk", "exprtk-custom-types"],
            libraries=["mpfr", "gmp", "tbb"],
            define_macros=[("PYVPMR", 1)],
            extra_compile_args=['-fno-aligned-allocation'] if "Darwin" == platform.system() else [],
        ),
    ],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
    packages=["pyvpmr"],
)
