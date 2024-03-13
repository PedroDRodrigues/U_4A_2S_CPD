CC = gcc
CFLAGS = -lm -fopenmp -O2

all: life3d

life3d: life3d-omp.c
	$(CC) -o life3d life3d-omp.c $(CFLAGS)

clean:
	rm -f life3d
