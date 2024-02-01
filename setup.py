#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import io
import os.path
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext
from typing import Any

from setuptools import find_packages
from setuptools import setup


def read(*names: Any, **kwargs: dict[str, str]):
    with io.open(join(dirname(__file__), *names), encoding=kwargs.get("encoding", "utf8")) as fh:  # type: ignore
        return fh.read()


setup(
    name="nsp13",
    version="0.1.0",
    license="MIT",
    description="Interdiciplinary Group Project on finding SARS-CoV-2 NSP13 inhibitors",
    long_description="{}".format(
        re.compile("^.. start-badges.*^.. end-badges", re.M | re.S).sub(
            "", read("README.md")
        ),
    ),
    long_description_content_type="text/markdown",
    author="Anastazja Avdonina, Julia Byrska, Jakub Guzek, Paulina Kucharewicz, Michalina Wysocka",
    url="file://" + os.path.abspath(dirname(__file__)),
    packages=find_packages("src"),
    package_dir={"": "src"},
    py_modules=[splitext(basename(path))[0] for path in glob("src/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.9",
        "Topic :: Utilities",
        "Private :: Do Not Upload",
    ],
    keywords=["Machine Learning", "Drug Discovery", "PCM", "nsp13", "SARS-CoV-2"],
    python_requires=">=3.9",
    install_requires=[
        "biopython>=1.83",
        "matplotlib>=3.8.2",
        "numpy>=1.26.3",
        "pandas>=2.2.0",
        "rdkit>=2023.9.4",
        "requests>=2.31.0",
        "scikit-learn==1.2.2",
        "scipy>=1.12.0",
        "seaborn>=0.13.2",
        "umap-learn>=0.5.5"
    ],
)
