from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
from datetime import datetime

# noinspection PyTypeChecker
setup(
    name="pyvpmr",
    version=datetime.now().strftime("%y%m%d"),
    author="Theodore Chang",
    author_email="tlcfem@gmail.com",
    url="https://github.com/TLCFEM/vpmr",
    description="The VPMR Algorithm",
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
