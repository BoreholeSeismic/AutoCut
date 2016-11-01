#!/bin/bash

cd ~/Download
wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
bash Miniconda2-latest-MacOSX-x86_64.sh
conda update conda
conda install -y scipy pandas
conda install -y -c conda-forge colorlog=2.7.0
conda install -y -c conda-forge tqdm=4.8.4
conda install -y -c siruix autocut