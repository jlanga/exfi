from setuptools import setup

from exfi import __version__

setup(
    name='exfi',
    version=__version__,
    packages=['exfi'],
    license='MIT',
    install_requires=[
        'numpy',
        'Biopython',
        'networkx',
        'pandas'
    ],
    long_description=open('README.md').read(),
    test_suite='nose.collector',
    tests_require=['nose'],
    scripts=[
        "bin/build_baited_bloom_filter",
        "bin/build_splicegraph",
        "bin/exons_to_gapped_transcript"
    ],
    include_package_data=True,
    zip_safe=False
)
