#!/usr/bin/env bash
# Conda packages
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda install --yes abyss bedtools biopython

# Manually installed repos
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    brew tap homebrew/science
    brew install biobloomtools
else
    git clone https://github.com/bcgsc/biobloom src/biobloom
    pushd src/biobloom
    ./autogen.sh && ./configure && make -j 2 && sudo make install -j 2
    popd
fi
