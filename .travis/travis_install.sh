#!/usr/bin/env bash

# Install conda packages
export PATH="$HOME/miniconda3_$TRAVIS_OS_NAME/bin:$PATH"
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda install --yes abyss biopython bedtools pandas networkx


# SDSL-lite
# https://hub.docker.com/r/adamnovak/sequence-graphs/~/dockerfile/
git clone https://github.com/simongog/sdsl-lite.git && \
pushd sdsl-lite && \
./install.sh /usr && \
popd && \
rm -Rf sdsl-lite

# biobloomtools
# as in https://github.com/bcgsc/biobloom
git clone https://github.com/bcgsc/biobloom.git && \
pushd biobloom/ && \
git submodule update --init && \
./autogen.sh && \
./configure && \
make -j 2 && \
make install && \
popd && \
rm -rf biobloom/ 
