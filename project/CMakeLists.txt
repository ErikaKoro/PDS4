cmake_minimum_required(VERSION 3.21)
project(PDS4 C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS -O3)

include_directories(src)

add_executable(PDS4
        src/quick_select.c
        src/quick_select.h
        src/sequential/sequential_vptree.c
        src/sequential/sequential_vptree.h
        src/timer.c)

add_executable(data
        src/Data/main.c)

add_executable(CILK
        src/quick_select.c
        src/quick_select.h
        src/cilk/cilk_vptree.c
        src/timer.c src/cilk/cilk_vptree.h)

target_link_libraries(CILK -fopencilk)

add_executable(KNN2
        src/quick_select.c
        src/quick_select.h
        src/timer.c
        src/knn/KNNSearch.c)

add_executable(MPI
        src/mpi/mpi_knn.c
        src/mpi/mpi_quick_select.h
        src/mpi/mpi_quick_select.c
        src/mpi/mpi_tree.c src/mpi/mpi_master.c src/mpi/mpi_master.h src/mpi/mpi_slave.c src/mpi/mpi_slave.h src/mpi/mpi_exchanges.c src/mpi/mpi_exchanges.h src/mpi/distribute_by_median.c src/mpi/distribute_by_median.h)

target_link_libraries(MPI -lm)

