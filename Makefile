MPI=mpicc
CFLAGS=-std=c11 -O3 -lm -Wall
RES=-DPRINT_RESULTS
CONV=-DTEST_CONV

all: mpi

mpi: mpi_skeleton.c utils.c
	$(MPI) $(CFLAGS) $(RES) $(CONV) mpi_skeleton.c utils.c -o mpi_skeleton

clean:
	rm mpi_skeleton

