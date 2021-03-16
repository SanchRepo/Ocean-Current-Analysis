# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 13:15:35 2017

@author: qzh14
"""

try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Distutils import build_ext
import numpy as np

ext_modules = [Extension("SoccerHycomHF",["SoccerHycomHF.pyx"])]
               

setup(
    name= 'Generic model class',
    cmdclass = {'Soccer_ext': build_ext},
    include_dirs = [np.get_include()],
    ext_modules = ext_modules)