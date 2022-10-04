CFLAGS=-O3 -g
CC=gcc
CILKCC=/usr/local/OpenCilk-2.0.0-x86_64-Linux-Ubuntu-22.04/bin/clang
MPICC=mpicc
default: all

vptree_build:
	@mkdir -p build
	$(CC) $(CFLAGS) -c ./src/sequential/sequential_vptree.c -o ./build/sequential_vptree.o
	$(CC) $(CFLAGS) -c ./src/quick_select.c -o ./build/quick_select.o
	$(CC) $(CFLAGS) -o ./build/sequential_vptree.out ./build/quick_select.o ./build/sequential_vptree.o

build_quick_select:
	@mkdir -p build
	$(CC) $(CFLAGS) -c ./src/quick_select.c -o ./build/quick_select.o
	$(CC) $(CFLAGS) -o ./build/quick_select.out ./build/quick_select.o

cilk_vptree:
	@mkdir -p build
	$(CILKCC) $(CFLAGS) -c ./src/cilk_vptree.c -o ./build/cilk_vptree.o -fopencilk
	$(CILKCC) $(CFLAGS) -c ./src/quick_select.c -o ./build/quick_select.o -fopencilk
	$(CILKCC) $(CFLAGS) -o ./build/cilk_vptree.out ./build/quick_select.o ./build/cilk_vptree.o -fopencilk

build_knn:
	@mkdir -p build
	$(CC) $(CFLAGS) -c ./src/KNNSearch.c -o ./build/KNNSearch.o
	$(CC) $(CFLAGS) -c ./src/quick_select.c -o ./build/quick_select.o
	$(CC) $(CFLAGS) -o ./build/KNNSearch.out ./build/quick_select.o ./build/KNNSearch.o

build_mpi:
	@mkdir -p build
	$(MPICC) $(CFLAGS) -c ./src/mpi_knn.c -o ./build/mpi_knn.o
	$(MPICC) $(CFLAGS) -c ./src/quick_select.c -o ./build/quick_select.o
	$(MPICC) $(CFLAGS) -o ./build/mpi_knn.out ./build/mpi_knn.o ./build/quick_select.o -lm


.PHONY: clean

all: vptree_build

run_quick_select: build_quick_select
	./build/quick_select.out

run_vptree: vptree_build
	valgrind --leak-check=yes --track-origins=yes --log-file=check.rpt ./build/sequential_vptree.out ./src/data

run_cilk_vptree: cilk_vptree
	./build/cilk_vptree.out ./src/data

run_knn: build_knn
	valgrind --leak-check=yes --track-origins=yes --log-file=check.rpt ./build/KNNSearch.out ./src/data 5

run_mpi_knn: build_mpi
	#mpirun -hostfile hosts  ./build/mpi_knn.out ./src/data
	mpirun -np 4 ./build/mpi_knn.out ./src/data


clean:
	rm -rf ./build/*.out
	rm -rf ./build/*.o
