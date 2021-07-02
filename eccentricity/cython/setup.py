""" Setup for hphc.pyx """

from setuptools import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize(["hphc.pyx"]))