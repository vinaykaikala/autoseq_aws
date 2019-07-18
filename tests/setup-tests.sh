#!/usr/bin/env bash

# Get and install anaconda for custom Python installation
echo -n "Downloading conda... "
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh 2> /dev/null > /dev/null
echo "Done!"
bash Miniconda-latest-Linux-x86_64.sh -b -p $HOME/miniconda2 -f

conda config --add channels r
conda config --add channels bioconda
conda config --add channels dakl

conda config --set always_yes yes --set changeps1 no

# Useful for debugging any issues with conda
conda info -a

conda install --file conda-list-tests.txt --quiet
