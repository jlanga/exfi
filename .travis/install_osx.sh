brew tap homebrew/science
brew install biobloom curl
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
export PATH="/home/travis/miniconda3/bin:$PATH"
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda install --yes abyss bedtools biopython
