#!/usr/bin/env bash

# Install miniconda

if test -e $HOME/miniconda3_$TRAVIS_OS_NAME/bin ; then
    echo "miniconda already installed."
else
    echo "Installing miniconda."
    rm -rf $HOME/miniconda3
    mkdir -p $HOME/download
    if [[ -d $HOME/download/miniconda_$TRAVIS_OS_NAME.sh ]] ; then
        rm -rf $HOME/download/miniconda_$TRAVIS_OS_NAME.sh
    fi
    if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
        url="https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    else
        url="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    fi
    wget \
        --continue \
        --output-document $HOME/download/miniconda_$TRAVIS_OS_NAME.sh \
        $url
    chmod +x $HOME/download/miniconda_$TRAVIS_OS_NAME.sh
    $HOME/download/miniconda_$TRAVIS_OS_NAME.sh \
        -b \
        -p $HOME/miniconda3_$TRAVIS_OS_NAME
fi


# Install brew
## https://gist.github.com/njlr/ac278ee448dc18002763092724d4e4d4
if test -e $HOME/.linuxbrew ; then
    echo "linubrew already installed"
else
    echo "Installing linuxbrew from GitHub"
    rm -rf $HOME/.linuxbrewbefore_install:
    git clone https://github.com/Linuxbrew/brew.git $HOME/.linuxbrew
    export PATH="$HOME/.linuxbrew/bin:$PATH"
    echo 'export PATH="$HOME/.linuxbrew/bin:$PATH"' >>~/.bash_profile
    export MANPATH="$(brew --prefix)/share/man:$MANPATH"
    export INFOPATH="$(brew --prefix)/share/info:$INFOPATH"
fi
brew --version
brew tap homebrew/science
brew install biobloomtools
