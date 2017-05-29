#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    brew tap homebrew/science
    brew install biobloom curl
    curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b
else
    sudo apt-get -qq update
    sudo apt-get install -y build-essential curl git libboost-dev
    curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b
    git clone https://github.com/bcgsc/biobloom src/biobloom
    pushd src/biobloom
    ./autogen.sh && ./configure && make -j 2 && sudo make install-j 2
    popd

fi

export PATH="/home/travis/miniconda3/bin:$PATH"
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda install --yes abyss bedtools biopython
