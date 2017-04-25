from setuptools import setup


setup(
    name='exfi',
    version='0..1',
    packages=['exfi'],
    license='MIT',
    install_requires=[
        'numpy',
        'Biopython'
    ],
    long_description=open('README.md').read(),
    test_suite='nose.collector',
    tests_require=['nose'],
    scripts = [
        "bin/find_exons",
        "bin/filter_exons_by_length",
        "bin/filter_exons_by_extensibility",
        "bin/paste_exons",
        "bin/fill_microindels"
    ],
    include_package_data=True,
    zip_safe=False
)
