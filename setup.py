from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
from datetime import datetime
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# noinspection PyTypeChecker
setup(
    name="pyvpmr",
    version=datetime.now().strftime("%y%m%d"),
    author="Theodore Chang",
    author_email="tlcfem@gmail.com",
    url="https://github.com/TLCFEM/vpmr",
    description="The VPMR Algorithm",
    long_description=long_description,
    long_description_content_type='text/markdown',
    ext_modules=[
        Pybind11Extension(
            "pyvpmr",
            ["src/VPMR.cpp"],
            include_dirs=["eigen", "eigen/unsupported", "mpreal", "exprtk", "exprtk-custom-types"],
            libraries=["mpfr", "gmp", "tbb"],
            define_macros=[("PYVPMR", 1)],
        ),
    ],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)
