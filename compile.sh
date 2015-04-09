export LANG=en_US.utf8
export LC_ALL=en_US.utf8

#Create CMAKE sub-folders, use -p to force success if folder already exists
#mkdir debug_mpi release_mpi -p

#Delete all that potentially already was there
#rm -rf debug_mpi/* release_mpi/*

module load intelcomp
module load openmpi/1.4.3-intel
#module load cmake

#Create makefiles for each sub-folder
#cd debug_mpi/
#CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Debug
#CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Debug
make clean && make && chmod +x ex6

#cd ../release_mpi/
#CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
#make clean && make && chmod +x ex6
