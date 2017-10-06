#!/usr/bin/env bash

# Install conda packages
export PATH="$HOME/miniconda3_$TRAVIS_OS_NAME/bin:$PATH"
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda install --yes \
    abyss \
    biopython \
    bedtools \
    pandas \
    networkx
conda clean --all --yes

# SDSL-lite
# https://hub.docker.com/r/adamnovak/sequence-graphs/~/dockerfile/
if [[ ! -d $HOME/sdsl-lite ]]; then
    git clone https://github.com/simongog/sdsl-lite.git
fi
pushd $HOME/sdsl-lite && \
sudo ./install.sh /usr/local && \
popd

# biobloomtools
# as in https://github.com/bcgsc/biobloom
if [[ ! -d $HOME/biobloom ]]; then
    git clone https://github.com/bcgsc/biobloom.git
fi
pushd $HOME/biobloom/ && \
git submodule update --init && \
./autogen.sh && \
./configure --prefix=/usr/local && \
make -j 2 && \
sudo make install && \
popd
