<div id="top"></div>

<br />
<div align="center">
  <h1 align="center">Parallel and Distributed Systems Assignment 4</h1>
  <h3 align="center">Aristotle University of Thessaloniki</h3>
  <h4 align="center">School of Electrical & Computer Engineering</h4>
  <p align="center">
    Koro Erika
    <br />
    Winter Semester 2021 - 2022
    <br />
    <br />
    <br />
    <br />
  </p>
</div>

- [1. About this project](#1-about-this-project)
- [2. Getting started](#2-getting-started)
- [3. Dependencies](#3-dependencies)
    - [3.0.1. Make or Cmake](#301-make-or-cmake)
    - [3.0.2. OpenMPI](#302-openmpi)
    - [3.0.3. OpenCilk-9.0.1](#303-opencilk-901)
- [4. Project directory layout](#4-project-directory-layout)
  - [4.1. PDS4](#41-pds4)
- [5. Compile and run](#5-compile-and-run)
  - [5.1. Linux](#51-linux)
  - [Sequential](#sequential)
  - [OpenCilk](#opencilk)
  - [KNN_Search](#knn_search)
  - [mpi](#mpi)
  - [Clean](#clean)
  - [5.3. Command line argumnets](#53-command-line-argumnets)
<br/>

## 1. About this project

<p align="justify">
  The objective of this assignment is to parallelize the KNN algorithm
  using Cilk and MPI C libraries.
<br/>
<br/>
  The first part of the project is to build sequentially the vantage point tree and then calculate its k-nearest neighbours.Secondly, this algorithm is used for calculating all points' k-nearest neighbours.Finally, cilk is used to parallelize the construction of the vantage point tree and the distances' calculations and MPI is implemented in order to parallelize the KNN algorithm.
<br/>
<br/>
</p>

## 2. Getting started

To setup this repository on your local machine run the following commands on the terminal:

```console
git clone git@github.com:ErikaKoro/PDS4.git
```

Or alternatively [*download*](https://github.com/ErikaKoro/PDS4.git) and extract the zip file of the repository
<br/>
<br/>

## 3. Dependencies
#### 3.0.1. Make or Cmake

This project uses make utilities to build and run the excecutables. Alternatively you can use Cmake but openCilk and MPI are not supported in the build script.

#### 3.0.2. OpenMPI

Make sure you have an mpicc capable compiler.

#### 3.0.3. OpenCilk-9.0.1

You can install OpenCilk by following the instructions of the official website. The official support is for Ubuntu but the binaries have been tested on Manjaro Linux as well and they are functional.

IMPORTANT! The path to the openCilk clang binary should be updated in the Makefile


<br/>

## 4. Project directory layout

### 4.1. PDS4

```
.PDS4
├── src                     main code
|   ├── Data                # dataset creation related source files
|   ├── CILK                # OpenCilk related source files
|   ├── knn                 # sequential knn related source files
|   ├── mpi                 # OpenMPI related source files
|   ├── sequential          # sequential related source files
|   ├── quick_select.c      # quick_select algorithm
|   ├── quick_select.h      # quick_select header file 
|   ├── timer.c             # calculate execution time source code  
|   └── timer.h             # timer.c header file
├── .gitignore          
├── CMakeLists.txt          # Cmake build script for Clion
├── Makefile                # Makefile for linux base distros
├── README.md
└── hosts                   # file that defines number of slots for MPI

```
<br/>
<br/>

## 5. Compile and run

### 5.1. Linux

### Sequential
Simply run `make run_sequential` in the [*PDS4*](PDS4) directory in order to run the executable for sequential code that creates the vantage point tree. The executable files will be created in the `BUILD_DIR` directory along with the `object` files.

### OpenCilk
Run `make run_cilk` in order to run the executable for the OpenCilk parallelized code.

### KNN_Search
Run `make run_knn` in order to run the executable for the KNN algorithm for all points of the dataset.

### mpi
Run `make run_mpi` in order to run the executable for the distributed solution for the KNN algorithm.

### Clean
Run `make clean` to clear all the object and executable files. For more information on the command line arguments read bellow.


<br/>

### 5.3. Command line argumnets

* Data executable
    ```C
    argv[1] = dimension of points
    argv[2] = number of points
    argv[3] = 'name of the file'
    ```
* Sequential executable
    ```C
    argv[1] = 'path to the dataset file'
    ```
* OpenCilk executable
    ```C
    argv[1] = 'Path to the dataset file'
    ```
* KNN executable
    ```C
    argv[1] = 'Path to the dataset file'
    argv[2] = number of searching neighbours for all points
    ```
* OpenMPI executable
    ```C
    argv[1] = 'Path to the dataset file'
    argv[2] = number of searching neighbours for all points
    ```
<br/>


