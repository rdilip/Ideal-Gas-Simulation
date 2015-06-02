from distutils.core import setup
from Cython.Build import cythonize
import sys

file_name = "idealGas"
inp_file = file_name + ".pyx"

setup(
		name = file_name,
		ext_modules = cythonize(inp_file),
	)
