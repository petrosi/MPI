MPI=mpicc
CFLAGS=-std=c11 -O3 -Wall
MATH=-lm
RES=-DPRINT_RESULTS
CONV=-DTEST_CONV

all: jacobi seidelsor redblacksor jacobi_mpi seidelsor_mpi redblacksor_mpi

jacobi: Jacobi.c utils.c
	$(MPI) $(CFLAGS) $(RES) $(CONV) Jacobi.c utils.c -o jacobi

seidelsor: GaussSeidelSOR.c utils.c
	$(MPI) $(CFLAGS) $(RES) $(CONV) GaussSeidelSOR.c utils.c -o seidelsor $(MATH)

redblacksor: RedBlackSOR.c utils.c
	$(MPI) $(CFLAGS) $(RES) $(CONV) RedBlackSOR.c utils.c -o redblacksor $(MATH)

jacobi_mpi: Jacobi.c utils.c
	$(MPI) $(CFLAGS) Jacobi.c utils.c -o jacobi_mpi

seidelsor_mpi: GaussSeidelSOR.c utils.c
	$(MPI) $(CFLAGS) GaussSeidelSOR.c utils.c -o seidelsor_mpi $(MATH)

redblacksor_mpi: RedBlackSOR.c utils.c
	$(MPI) $(CFLAGS) RedBlackSOR.c utils.c -o redblacksor_mpi $(MATH)


clean:
	rm jacobi seidelsor redblacksor jacobi_mpi seidelsor_mpi redblacksor_mpi

