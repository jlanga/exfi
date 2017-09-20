from setuptools import setup


setup(
    name='exfi',
    version='1.1.0',
    packages=['exfi'],
    license='MIT',
    install_requires=[
        'numpy',
        'Biopython'
    ],
    long_description=open('README.md').read(),
    test_suite='nose.collector',
    tests_require=['nose'],
    scripts=[
        "bin/build_baited_bloom_filter",
        "bin/find_exons",
        "bin/reduce_exons",
        "bin/exons_to_gapped_transcript"
    ],
    include_package_data=True,
    zip_safe=False
)
