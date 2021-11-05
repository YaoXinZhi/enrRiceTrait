# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 28/03/2021 22:58
@Author: XINZHI YAO
"""

# from distutils.core import setup
import setuptools
from setuptools import find_packages, setup

with open("README.md", "r") as fh:
  long_description = fh.read()

setup(
  name='enrRiceTrait',
  version='1.1.6',
  author='xinzhi_yao',
  author_email='xinzhi_bioinfo@163.com',
  description="A python package for Rice Trait Enrichment.",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/YaoXinZhi/enrRiceTrait",

  packages=setuptools.find_packages(),

  # include_package_data is bdist, and MANIFEST.in is sdist.
  # include_package_data=True,
  package_data={
    'enrRiceTrait': ['data/to.wto.ro.funRiceGene.obo',
                     'data/Total_Association.txt'],
    'data': ['*']
  },

  classifiers=[
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
  ],
  python_requires='>3.7'
)


