#include "mpi_master.h"
#include "mpi_tree.h"
#include <stdlib.h>
#include <mpi.h>

/**
 * This function finds which processes has the closest points to the pivot, builds its VPTree according to the median distance from the pivot
 * and send its points to the master. Master process calls the KNNSearch function and announces the k-nearest neighbours
 * @param pivot The initial chosen pivot
 * @param rank The process id
 * @param size The world size of the communicator
 * @param pointsPerProc the number of points per process
 * @param holdPoints The array that holds each process's points
 * @param dimension The dimension of the points
 * @param k_neighbours The number of neighbours needed to be found
 */
void calculateMasterKNN(KNN perPivot, const double *pivot, int rank, int size, int64_t pointsPerProc, double **holdPoints, int64_t dimension, int k_neighbours){


    // Find which processes will construct their VPTree in order to find the k-nearest neighbours
    int64_t limit = 0;
    if(k_neighbours % pointsPerProc != 0){
        limit = k_neighbours / pointsPerProc + 1;
    }
    else{
        limit = k_neighbours / pointsPerProc;
    }


    vptree initial;
    initial.inner = NULL;
    initial.outer = NULL;
    initial.start = 0;
    initial.stop = pointsPerProc - 1;

    // Allocate memory for the vantage point and fill it with the last point's coordinates(random choice)
    initial.vpPoint = (double *) malloc(dimension * sizeof(double));     // FREEMEM (check)
    for (int i = 0; i < dimension; ++i) {
        initial.vpPoint[i] = pivot[i];
    }

    // Update the distances
    double *distancesPerProc = (double *) calloc(pointsPerProc, sizeof(double));        // FREEMEM(CHECK)

    findDistance(distancesPerProc, holdPoints, dimension, initial.vpPoint, pointsPerProc);

    buildVPTree(&initial, holdPoints, distancesPerProc, dimension, pointsPerProc, k_neighbours);      // FREEMEM(CHECK)

    freeMpiMemory(initial.inner);
    if(initial.outer != NULL) {
        freeMpiMemory(initial.outer);
    }
    free(initial.vpPoint);
    free(distancesPerProc);

    // This condition checks whether the master's points are enough for the needed neighbours
    if(pointsPerProc > k_neighbours) {
        knn_mpi_search(k_neighbours, perPivot.nearest, holdPoints, pointsPerProc, dimension);

    }
    else{
        double *recvbuffer;

        // The number of elements that have been requested to be sent
        int *elements = (int *) malloc((limit - 1) * sizeof(int));  // FREEMEM(CHECK)
        MPI_Request *requests = (MPI_Request *) malloc(limit * sizeof(MPI_Request));    // FREEMEM(CHECK)

        unsigned long int reallocSize = 0;

        for (int i = 0; i < limit - 1; i++) {

            MPI_Status status;
            MPI_Probe(i + 1, 12, MPI_COMM_WORLD, &status);  // Ask MPI to give you the size of the message

            MPI_Get_count(&status, MPI_DOUBLE, &elements[i]);

            if(i == 0){
                // malloc the first time
                recvbuffer = (double *)malloc(elements[i] * dimension * sizeof(double));    // FREEMEM (check)
                reallocSize += elements[i] * dimension * sizeof(double);
//                printf("receive buffer initial address: %p\n", recvbuffer);
            }
            else {
                // realloc
                recvbuffer = realloc(recvbuffer, reallocSize + elements[i] * dimension * sizeof(double));   // FREEMEM(CHECK)
                reallocSize += elements[i] * dimension * sizeof(double);

//                printf("receive buffer new address: %p\n", recvbuffer);
            }

            MPI_Irecv(recvbuffer + i * pointsPerProc * dimension,
                      elements[i] * (int) dimension,
                      MPI_DOUBLE,
                      i + 1,
                      12,
                      MPI_COMM_WORLD,
                      &requests[i]
            );

        }

        for (int i = 0; i < limit - 1; ++i) {
            MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
        }

        // Copy the points of the master rank to the nearest array
        for (int i = 0; i < pointsPerProc; ++i) {
            for (int j = 0; j < dimension; ++j) {
                perPivot.nearest[i][j] = holdPoints[i][j];
            }
        }

        // Copy the received points to the nearest array
        int index = 0;
        for (int i = 1; i < limit; ++i) {
            for (int j = 0; j < elements[i - 1] / dimension; ++j) {
                for (int k = 0; k < dimension; ++k) {
                    perPivot.nearest[pointsPerProc * i + j][k] = recvbuffer[index];
                    index++;
                }
            }
        }

//        free(recvbuffer);
        free(elements);
    }

}