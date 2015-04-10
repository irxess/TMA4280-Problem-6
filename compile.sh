#! /bin/bash
export LANG=en_US.utf8
export LC_ALL=en_US.utf8

module load intelcomp
module load openmpi/1.4.3-intel

rm ex6
make clean && make

if [ -f ex6 ]; then
	chmod +x ex6
	mpirun -np 4 ./ex6 3 4
	mpirun -np 4 ./ex6 3 6
	mpirun -np 4 ./ex6 3 16
fi
