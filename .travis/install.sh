#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    bash ./travis/install_osx.sh
else
    bash ./travis/install_linux.sh
fi
