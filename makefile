CC=clang
FC=gfortran

poisson: poisson.o fst.o
	$(CC) -fopenmp  $^ $(shell pkg-config ompi-c --libs) -lm

fst.o:  fst.f
	$(FC) -c  $< 

%.o : %.c 
	$(CC) -c $< $(shell pkg-config ompi-c --cflags) 

.PHONY:	
clean:	
	rm -f *.o 
