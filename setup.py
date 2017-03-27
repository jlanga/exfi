from setuptools import setup


setup(
    name='exfi',
    version='0.1dev0',
    packages=['exfi'],
    license='MIT',
    install_requires=[
        'Biopython'
    ],
    long_description=open('README.rst').read(),
    test_suite='nose.collector',
    tests_require=['nose'],
    entry_points={
          'console_scripts': ['funniest-joke=funniest.command_line:main'],
      },
    scripts = [
        "bin/find_exons",
    ],
    include_package_data=True,
    zip_safe=False
)
