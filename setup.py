#!/usr/bin/env python3
"""
setup.py: installer for exfi
"""

from setuptools import setup

from exfi import __version__

setup(
    name='exfi',
    version=__version__,
    packages=[
        'exfi',
        'exfi.io'
    ],
    license='MIT',
    install_requires=[
        'numpy',
        'Biopython',
        'networkx',
        'pandas',
        'natsort',
    ],
    long_description=open('README.md').read(),
    test_suite='nose.collector',
    tests_require=['nose'],
    scripts=[
        "bin/build_baited_bloom_filter",
        "bin/build_splice_graph",
        "bin/gfa1_to_exons",
        "bin/gfa1_to_gapped_transcripts"
    ],
    include_package_data=True,
    zip_safe=False
)
