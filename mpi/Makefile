CC = mpicc
CFLAGS = -lm -fopenmp -O2

all: life3d

life3d: life3d-mpi.c
	$(CC) -o life3d-mpi life3d-mpi.c $(CFLAGS)

clean:
	rm -f life3d-mpi
