#!/usr/bin/env bash

# Install conda packages
export PATH="$HOME/miniconda3_$TRAVIS_OS_NAME/bin:$PATH"
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda install --yes abyss biopython bedtools pandas networkx

git clone https://github.com/bcgsc/biobloom ~/download/biobloom
pushd ~/download/biobloom
./autogen.sh
./configure
make -j
sudo make install
popd
