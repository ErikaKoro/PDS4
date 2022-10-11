#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>
#include <math.h>

#include "mpi_master.h"
#include "mpi_slave.h"
#include "mpi_quick_select.h"
#include "distribute_by_median.h"


/**
 * check if the processes with ranks between 0 and worldSize(exclusive) have distances bigger than the biggest distance they received from the previous process(sorted processes)
 * so as to know if the processes in the end have the right points
 *
 * @param points
 * @param pointsPerProc
 * @param dimension
 * @param pivot
 */
void mpi_test_function(double **points, int pointsPerProc, int dimension, double *pivot) {

    int rank;
    int worldSize;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    double *dist = (double *) malloc(pointsPerProc * sizeof (double));

    findMPIDistances(rank, dist, points, dimension, pivot, pointsPerProc);

    double maxDist = dist[0];
    double minDist = dist[0];
    for (int i = 0; i < pointsPerProc; ++i) {
        if(dist[i] > maxDist){
            maxDist = dist[i];
        }

        if(dist[i] < minDist) {
            minDist = dist[i];
        }
    }

    if (rank == 0){
        MPI_Send(&maxDist, 1, MPI_DOUBLE, rank + 1, 50, MPI_COMM_WORLD);
    }
    else if (rank == worldSize - 1){
        double prevMax;

        MPI_Recv(&prevMax, 1, MPI_DOUBLE, rank - 1, 50, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (prevMax < maxDist) {
            printf("Rank: %d SUCCESS\n", rank - 1);
            printf("Rank: %d SUCCESS\n", rank);
        } else {
            printf("Rank: %d FAILURE\n", rank - 1);

        }


    } else {
        double prevMax;
        MPI_Request request_1;
        MPI_Request request_2;

        MPI_Irecv(&prevMax, 1, MPI_DOUBLE, rank - 1, 50, MPI_COMM_WORLD, &request_1);
        MPI_Isend(&maxDist, 1, MPI_DOUBLE, rank + 1, 50, MPI_COMM_WORLD, &request_2);

        MPI_Wait(&request_1, NULL);
        MPI_Wait(&request_2, NULL);

        if (prevMax < maxDist) {
            printf("Rank: %d SUCCESS\n", rank - 1);
        } else {
            printf("Rank: %d FAILURE\n", rank - 1);

        }
    }
    free(dist);
}

void printArray(double **points, int64_t elements, int64_t dimension){
    for (int i = 0; i < elements; ++i) {
        for (int j = 0; j < dimension; ++j) {
            printf("%.1f ", points[i][j]);
        }
        printf("\n");
    }
}

bool isPowerOfTwo(int64_t number){
    return (ceil(log2((double)number)) == floor(log2((double)number)));
}


/**
 * This function copies one array to another
 * @param points The array that need to be copied
 * @param copied The array where the other arrays is copied to
 * @param numberOfPoints The number of array's rows
 * @param dimension The number of array's cols
 */
void deepCopyArray(double **points, double **copied, int64_t numberOfPoints, int64_t dimension){

    for (int i = 0; i < numberOfPoints; ++i) {
        for (int j = 0; j < dimension; ++j) {
            copied[i][j] = points[i][j];
        }
    }
}


int main(int argc, char **argv) {
    int size;
    int rank;


    MPI_Init(&argc, &argv);
    double start, end;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (argc != 3)
        printf("\nGimme args!\n");

    /* Read the binary file with our data back in */
    // The first two 64-bit numbers from the file give us the number of our points and the space's dimension
    int64_t numberOfPoints;
    int64_t dimension;
    int64_t pointsPerProc;
    double **holdThePoints;

    FILE *fh = fopen(argv[1], "rb");
    if (fh != NULL) {
        fread(&dimension, sizeof(int64_t), 1, fh);
        fread(&numberOfPoints, sizeof(int64_t), 1, fh);

        // hold number of points that is a power of 2
        if (!isPowerOfTwo(numberOfPoints)) {
            int64_t power = 1;
            while (power < numberOfPoints)
                power *= 2;

            numberOfPoints = power / 2;
            if (rank == 0)
                printf("The number of points is: %ld\n", numberOfPoints);
        } else {
            if (rank == 0)
                printf("The number of points is: %ld\n", numberOfPoints);
        }


        pointsPerProc = numberOfPoints / size;
//        printf("Rank: %d, Points per proc: %ld\n", rank, pointsPerProc);

        // hold our points with their coordinates read from the file
        holdThePoints = (double **) malloc(pointsPerProc * sizeof(double *));   // FREEMEM
        for (int i = 0; i < pointsPerProc; i++) {
            holdThePoints[i] = (double *) malloc(dimension * sizeof(double));
        }

        // Seek where is the pointer in our file
        fseek(fh, (int)(rank * dimension * pointsPerProc * sizeof(double)), SEEK_CUR);

        for (int i = 0; i < pointsPerProc; i++) {
            double *x2 = (double *) malloc(sizeof(double) * dimension); // represents only one point

            fread(x2, sizeof(double), dimension, fh);
            holdThePoints[i] = x2;
        }
        fclose(fh);


        double *pivot = (double *) calloc(dimension, sizeof(double));

        /*for (int l = 0; l < size; l++) {
            if (l == rank) {
                // Select pivot
                for (int i = 0; i < pointsPerProc; ++i) {


                    printf("The pivot is: ");
                    for (int j = 0; j < dimension; j++) {
                        pivot[j] = holdThePoints[i][j];
                        printf("%.1f ", pivot[j]);
                    }


                }
            }

            int k_neighbours = atoi(argv[2]);
            if (k_neighbours > numberOfPoints)
                k_neighbours = (int) numberOfPoints;

            if (rank == l) {
                start = MPI_Wtime();
            }
            //Broadcast the pivot to the processes
            MPI_Bcast(pivot, (int) dimension, MPI_DOUBLE, l, MPI_COMM_WORLD);

            distributeByMedian(pivot,
                               0,
                               rank,
                               (int) dimension,
                               holdThePoints,
                               (int) pointsPerProc,
                               size,
                               MPI_COMM_WORLD
            );

            if (rank == l) {
                end = MPI_Wtime();
                printf("The time is: %.4f\n", end - start);
            }

            mpi_test_function(holdThePoints, (int) pointsPerProc, (int) dimension, pivot);

            KNN *totalNearest;
            if (l == rank) {
                totalNearest = (KNN *) malloc(numberOfPoints * sizeof(double));  // FREEMEM

                for (int i = 0; i < numberOfPoints; ++i) {

                    // for each point, each totalNearest array has k points
                    totalNearest[i].nearest = (double **) malloc(sizeof(double *) * k_neighbours);

                    for (int j = 0; j < k_neighbours; ++j) {
                        totalNearest[i].nearest[j] = (double *) malloc(sizeof(double) * dimension);
                    }
                }
                calculateMasterKNN(totalNearest[0],
                                   pivot,
                                   rank,
                                   size,
                                   pointsPerProc,
                                   holdThePoints,
                                   dimension,
                                   k_neighbours
                );
            }

            calculateSlaveKNN(pivot, rank, size, pointsPerProc, holdThePoints, dimension, k_neighbours);

            if (rank == l) {
                printf("\nThe nearest points are:\n");
                for (int i = 0; i < k_neighbours; ++i) {
                    for (int j = 0; j < dimension; ++j) {
                        printf("%.1f ", totalNearest->nearest[i][j]);
                    }
                    printf("\n");
                }

                for (int m = 0; m < numberOfPoints; ++m) {
                    for (int n = 0; n < k_neighbours; ++n) {
                        free(totalNearest[m].nearest[n]);
                    }
                }

                for (int m = 0; m < numberOfPoints; ++m) {
                    free(totalNearest[m].nearest);
                }
                free(totalNearest);

            }


            for (int i = 0; i < pointsPerProc; ++i) {
                free(holdThePoints[i]);
            }
            free(holdThePoints);

            MPI_Barrier(MPI_COMM_WORLD);

            free(pivot);
        }*/

        int k_neighbours = atoi(argv[2]);

        if(k_neighbours > numberOfPoints)
            k_neighbours = (int)numberOfPoints;

        double **copied = (double **)malloc(sizeof (double *) * pointsPerProc);        // FREEMEM
        for (int i = 0; i < pointsPerProc; ++i) {
            copied[i] = (double *) malloc(sizeof (double ) * dimension);        // FREEMEM
        }

        deepCopyArray(holdThePoints, copied, pointsPerProc, dimension);

        for (int l = 0; l < numberOfPoints; ++l) {

            if(rank == l / pointsPerProc){
                for (int j = 0; j < dimension; ++j) {
                    pivot[j] = copied[l % pointsPerProc][j];
                }
//                printf("\nThe rank is %d\n", l / (int)pointsPerProc);
            }

            //Broadcast the pivot to the processes
            MPI_Bcast(pivot, (int) dimension, MPI_DOUBLE, l / (int)pointsPerProc, MPI_COMM_WORLD);

            if (rank == 0) {
                start = MPI_Wtime();
            }

            distributeByMedian(pivot,
                               0,
                               rank,
                               (int) dimension,
                               holdThePoints,
                               (int) pointsPerProc,
                               size,
                               MPI_COMM_WORLD
            );

            if (rank == 0) {
                end = MPI_Wtime();
                printf("The time is: %.4f\n", end - start);
            }


            mpi_test_function(holdThePoints, (int) pointsPerProc, (int) dimension, pivot);

            KNN *totalNearest;
            if(rank == 0) {
                totalNearest = (KNN *)malloc(numberOfPoints * sizeof(double));  // FREEMEM

                for (int i = 0; i < numberOfPoints; ++i) {

                    // for each point, each totalNearest array has k points
                    totalNearest[i].nearest = (double **) malloc(sizeof(double *) * k_neighbours);

                    for (int j = 0; j < k_neighbours; ++j) {
                        totalNearest[i].nearest[j] = (double *) malloc(sizeof(double) * dimension);
                    }
                }
                calculateMasterKNN(totalNearest[0], pivot, rank, size, pointsPerProc, holdThePoints, dimension, k_neighbours);
            }

            calculateSlaveKNN(pivot, rank, size, pointsPerProc, holdThePoints, dimension, k_neighbours);

            if(rank == 0){
                for (int k = 0; k < numberOfPoints; ++k) {
                    for (int j = 0; j < k_neighbours; ++j) {
                        free(totalNearest[k].nearest[j]);
                    }
                }

                for (int i = 0; i < numberOfPoints; ++i) {
                    free(totalNearest[i].nearest);
                }
                free(totalNearest);
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }

        for (int i = 0; i < pointsPerProc; ++i) {
            free(holdThePoints[i]);
        }

        free(holdThePoints);
        free(pivot);
        MPI_Finalize();

        return 0;
    }
}
