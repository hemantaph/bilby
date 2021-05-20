""" Setup for Coeff.pyx and WaVel.pyx """

from setuptools import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize(["Coeff.pyx", "WaVel.pyx", "function.pyx"]))

#python setup.py build_ext --inplace
