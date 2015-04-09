CC=mpicc
FC=ifort

poisson: poisson.o fst.o
	$(CC) -openmp  $^ -lm -o ex6

fst.o:  fst.f
	$(FC) -c  $< 

%.o : %.c 
	$(CC) -openmp -c $< 

.PHONY:	
clean:	
	rm -f *.o 
