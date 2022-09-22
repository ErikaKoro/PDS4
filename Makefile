CFLAGS=-O3 -g
CC=gcc
default: all

vptree_build:
	@mkdir -p build
	$(CC) $(CFLAGS) -c ./src/vptree.c -o ./build/vptree.o
	$(CC) $(CFLAGS) -c ./src/quick_select.c -o ./build/quick_select.o
	$(CC) $(CFLAGS) -o ./build/vptree.out ./build/quick_select.o ./build/vptree.o


build_quick_select:
	@mkdir -p build
	$(CC) $(CFLAGS) -c ./src/quick_select.c -o ./build/quick_select.o
	$(CC) $(CFLAGS) -o ./build/quick_select.out ./build/quick_select.o

.PHONY: clean

all: vptree_build

run_quick_select: build_quick_select
	./build/quick_select.out

run_vptree: vptree_build
	valgrind --leak-check=yes --track-origins=yes --log-file=check.rpt ./build/vptree.out ./src/data

clean:
	rm -rf ./build/*.out
	rm -rf ./build/*.o
