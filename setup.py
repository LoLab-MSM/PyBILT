#!/usr/bin/env python

from distutils.core import setup

setup(name='VORBILT',
      version='1.0.0',
      description='Lipid bilayer molecular simulation trajectory analysis toolkit',
      author='Blake Wilson',
      author_email='blake.a.wilson@vanderbilt.edu',
      url='https://github.com/blakeaw/VORBILT',
      packages=['MDAnalysis', 'numpy', 'scipy', 'matplotlib', 'seaborn'],
     )
