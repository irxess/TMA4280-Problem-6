#! /bin/bash
export LANG=en_US.utf8
export LC_ALL=en_US.utf8

module load intelcomp
module load openmpi/1.4.3-intel

rm ex6
make clean && make && chmod +x ex6
