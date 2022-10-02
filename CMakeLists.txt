cmake_minimum_required(VERSION 3.21)
project(PDS4 C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS -O3)

include_directories(src)

add_executable(PDS4
        src/quick_select.c
        src/quick_select.h
        src/sequential_vptree.c
        src/sequential_vptree.h)

add_executable(data
        src/main.c)

add_executable(CILK
        src/quick_select.c
        src/quick_select.h
        src/cilk_vptree.c
        src/cilk_vptree.h)

target_link_libraries(CILK -fopencilk)

add_executable(KNN2
        src/quick_select.c
        src/quick_select.h
        src/KNNSearch.c)

add_executable(MPI
        src/mpi_knn.h
        src/mpi_knn.c
        src/mpi_quick_select.h
        src/mpi_quick_select.c)

target_link_libraries(MPI -lm)

add_executable(Cilk_test src/cilk_test.c)

target_link_libraries(Cilk_test -fopencilk)
