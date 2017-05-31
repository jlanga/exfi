#!/usr/bin/env bash

# Install conda packages
export PATH="$HOME/miniconda3/bin:$PATH"
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda install --yes abyss bedtools biopython

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    brew tap homebrew/science
    brew install biobloomtools gcc
else
    git clone https://github.com/bcgsc/biobloom ~/download/biobloom
    pushd ~/download/biobloom
    ./autogen.sh && ./configure && make -j 2 && sudo make install
    popd
fi
