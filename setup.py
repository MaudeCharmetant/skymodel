#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from setuptools import setup

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

setup(
    name="skymodel",
    version="1.0",
    author="Maude Charmetant & Jens Erler",
    author_email="mcharmetant@astro.uni-bonn.de, jens@astro.uni-bonn.de",
    packages=["skymodel"],
    url="https://github.com/MaudeCharmetant/skymodel",
    license="MIT License",
    description=("A template-based model of the microwave sky for CCAT-p science forecasts"),
    long_description=open("README.rst").read(),
    package_data={"skymodel": ["LICENSE", "masks/*.fits"]},
    include_package_data=True,
    install_requires=["numpy", "healpy", "astropy", "pysm", "tqdm"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    zip_safe=False,
)
