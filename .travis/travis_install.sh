#!/usr/bin/env bash

# Install conda packages
<<<<<<< HEAD
export PATH="$HOME/miniconda3/bin:$PATH"
=======
export PATH="$HOME/miniconda3_$TRAVIS_OS_NAME/bin:$PATH"
>>>>>>> devel
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
<<<<<<< HEAD
conda install --yes abyss bedtools biopython

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    brew tap homebrew/science
    brew install biobloomtools
else
    git clone https://github.com/bcgsc/biobloom ~/download/biobloom
    pushd ~/download/biobloom
    ./autogen.sh && ./configure && make -j 2 && sudo make install
    popd
fi
=======
conda install --yes abyss biopython bedtools

git clone https://github.com/bcgsc/biobloom ~/download/biobloom
pushd ~/download/biobloom
./autogen.sh
./configure
make -j
sudo make install
popd
>>>>>>> devel
