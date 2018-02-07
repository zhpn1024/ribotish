#!/usr/bin/python
#Last-modified: 15 Jan 2018

#         Module/Scripts Description
# 
# Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
# 
# This code is free software; you can redistribute it and/or modify it
# under the terms of the BSD License (see the file COPYING included with
# the distribution).
# 
# @status:  Alpha
# @version: 0.2.0
# @author:  Peng Zhang
# @contact: zhpn1024@163.com

# ------------------------------------
# python modules
# ------------------------------------

import os,sys
from setuptools import setup, find_packages, Extension

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__ == '__main__':
    if float(sys.version[:3])<2.7 : #  or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.7 or higher!\n")
        sys.exit(1)

    # includepy = "%s/include/python%s" % (sys.prefix, sys.version[:3])
    with open("README.rst",'r') as fh:
        long_description = fh.read()
    # version
    from src import __version__ as VERSION
    PROG = 'ribotish'
    #with open('RELEASE','r') as fh:
        #PROG, VERSION = fh.next().split()[:2]

    if 'clean' in sys.argv:
        print >>sys.stderr, "Clean dist and egg info ..."
        os.system('if [ -d dist ]; then rm -rf dist; fi')
        os.system('if [ -f ngslib.egg-info ]; then rm ngslib.egg-info; fi')
        os.system('if [ -d ngslib.egg-info ]; then rm -rf ngslib.egg-info; fi')
        sys.exit(0)
    
    # install requirement
    install_requires = [ ["pysam >= 0.8.3"],
                         ["matplotlib >= 1.4.3"],
                         ["scipy >= 0.15.1"],
                         ]

    setup(name=PROG,
          version=VERSION,
          author='Peng Zhang',
          author_email='zhpn1024@163.com',
          url='https://github.com/zhpn1024/ribotish',
          license="GNU General Public License (GPL)",
          keywords = "Python, Riboseq, translation",
          description = ("Python Modules for Riboseq data analysis."),
          long_description = long_description,
          package_dir={PROG:'src'},
          packages = [PROG,'ribotish.run','ribotish.zbio'],
          scripts=['bin/ribotish',
                   #'bin/DumpTargetRNAs.py',
                   #'bin/EVDPermutation.py',
                    ],
          classifiers=['Environment :: Console',
                       'Development Status :: 3 - Alpha',
                       'Intended Audience :: Developers',
                       'License :: OSI Approved :: GNU General Public License (GPL)',
                       'License :: Free for non-commercial use',
                       'Operating System :: Unix',
                       'Programming Language :: Python :: 2.7',
                       'Programming Language :: Python :: 3',
                       'Topic :: Scientific/Engineering :: Bio-Informatics'],
          install_requires=install_requires)

