cmake_minimum_required(VERSION 3.21)
project(PDS4 C)

set(CMAKE_C_STANDARD 11)

include_directories(src)

add_executable(PDS4
        src/quick_select.c
        src/quick_select.h
        src/vptree.c
        src/vptree.h)

add_executable(data
        src/main.c)
