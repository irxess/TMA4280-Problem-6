CC=mpicc
CXX=icpc
FC=ifort

poisson: poisson.o fst.o
	$(CC) -openmp  $^ -lm -o ex6

fst.o:  fst.f
	$(FC) -c  $< 

%.o : %.c 
	$(CC) -c -g -openmp $<

clean:	
	rm -f *.o 

.PHONY:	clean
