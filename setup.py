from setuptools import setup


setup(
    name='exfi',
    version='0.1dev0',
    packages=['exfi'],
    license='MIT',
    install_requires=[
        'Biopython'
    ],
    long_description=open('README.md').read(),
    test_suite='nose.collector',
    tests_require=['nose'],
    entry_points={
          'console_scripts': ['funniest-joke=funniest.command_line:main'],
      },
    scripts = [
        "bin/find_exons",
        "bin/filter_exons",
        "bin/paste_exons",
        "bin/fill_microindels"
    ],
    include_package_data=True,
    zip_safe=False
)
