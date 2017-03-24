from setuptools import setup


setup(
    name='exfi',
    version='0.1dev',
    packages=['exfi'],
    license='MIT',
    install_requires=[
        'Biopython'
    ],
    long_description=open('README.rst').read(),
    test_suite='nose.collector',
    tests_require=['nose'],
    include_package_data=True,
    zip_safe=False
)
