MPI=mpicc
CFLAGS=-std=c11 -O3 -Wall
MATH=-lm
RES=-DPRINT_RESULTS
CONV=-DTEST_CONV

all: Jacobi GaussSeidelSOR RedBlackSOR

Jacobi: Jacobi.c utils.c
	$(MPI) $(CFLAGS) $(RES) $(CONV) Jacobi.c utils.c -o Jacobi

GaussSeidelSOR: GaussSeidelSOR.c utils.c
	$(MPI) $(CFLAGS) $(RES) $(CONV) GaussSeidelSOR.c utils.c -o GaussSeidelSOR $(MATH)

RedBlackSOR: RedBlackSOR.c utils.c
	$(MPI) $(CFLAGS) $(RES) $(CONV) RedBlackSOR.c utils.c -o RedBlackSOR $(MATH)

clean:
	rm Jacobi GaussSeidelSOR RedBlackSOR

