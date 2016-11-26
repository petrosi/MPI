MPI=mpicc
CFLAGS=-std=c11 -O3 -lm -Wall
RES=-DPRINT_RESULTS
CONV=-DTEST_CONV

all: Jacobi

Jacobi: Jacobi.c utils.c
	$(MPI) $(CFLAGS) $(RES) $(CONV) Jacobi.c utils.c -o Jacobi

clean:
	rm mpi_skeleton

