# setup.py (example)
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import pybind11
ext_modules = [
    Pybind11Extension(
        "triangulation",
        ["bindings.cpp"],   # include both files
        # include_dirs=[pybind11.get_include(),'.'],
        cxx_std=17,
    ),
]

setup(
    name="triangulation",
    version="0.0.1",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
