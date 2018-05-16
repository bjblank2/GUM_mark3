__author__ = 'brian'
from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize('clustersCY.pyx'),include_dirs=[numpy.get_include()])
setup(ext_modules = cythonize('clustermag_rulesCY.pyx'),include_dirs=[numpy.get_include()])
setup(ext_modules = cythonize('jsCY.pyx'),include_dirs=[numpy.get_include()])
setup(ext_modules = cythonize('mc_neighborObjCY.pyx'),include_dirs=[numpy.get_include()])
setup(ext_modules = cythonize('mc_siteObjCY.pyx'),include_dirs=[numpy.get_include()])
#setup(ext_modules = cythonize('mc_supercellCY.pyx'),include_dirs=[numpy.get_include()])
#setup(ext_modules = cythonize('mc_functions_2CY.pyx'),include_dirs=[numpy.get_include()])