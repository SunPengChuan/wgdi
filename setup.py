#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from setuptools import setup, find_packages

with open("README.md", "r",encoding='utf-8') as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name="WGDI",
    version="0.2.1",
    author="Pengchuan Sun",
    author_email="sunpengchuan@gmail.com",
    description="Whole Genome Duplication Identification",
    license="BSD License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SunPengChuan/wgdi",
    packages=find_packages(),
    package_data={'':['*.conf'],'wgdi': ['*.ini']},
    classifiers=(
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ),
    entry_points={
        'console_scripts':[
            'wgdi = wgdi.run:main',
        ]
      },
    zip_safe=True,
    install_requires=required
)