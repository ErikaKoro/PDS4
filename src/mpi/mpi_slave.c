#include "mpi_slave.h"
#include "mpi_tree.h"
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void calculateSlaveKNN(const double *pivot, int rank, int size, int64_t pointsPerProc, double **holdPoints, int64_t dimension, int k_neighbours){
    if (rank * pointsPerProc < k_neighbours) {
        double *sendPoints;

        vptree initial;
        initial.inner = NULL;
        initial.outer = NULL;
        initial.start = 0;
        initial.stop = pointsPerProc - 1;

        // Allocate memory for the vantage point and fill it with the last point's coordinates(random choice)
        initial.vpPoint = (double *) malloc(dimension * sizeof(double));     // FREEMEM(CHECK)
        for (int i = 0; i < dimension; ++i) {
            initial.vpPoint[i] = pivot[i];
        }

        // Update the distances
        double *distancesPerProc = (double *) calloc(pointsPerProc, sizeof(double));        // FREEMEM(CHECK)

        findDistance(distancesPerProc, holdPoints, dimension, initial.vpPoint, pointsPerProc);

        buildVPTree(&initial, holdPoints, distancesPerProc, dimension, pointsPerProc);      // FREEMEM(CHECK)

        freeMpiMemory(initial.inner);
        freeMpiMemory(initial.outer);
        free(distancesPerProc);


        if ((rank + 1) * pointsPerProc <= k_neighbours) {

            // The buffer with the points that will be sent to the master  // FREEMEM(CHECK)
            sendPoints = (double *) malloc(pointsPerProc * dimension * sizeof(double));

            for (int i = 0; i < pointsPerProc; i++) {
                memcpy(sendPoints + i * dimension, holdPoints[i],
                       dimension * sizeof(double));  // copy the exchange points from holdPoints to the sendPoints
                //printf("\nThe sendPoints[i] is: %.1f\n", sendPoints[i]);
            }

            MPI_Send(sendPoints, (int) (pointsPerProc * dimension), MPI_DOUBLE, 0, 12, MPI_COMM_WORLD);

        } else {

            int64_t pointsToGive = pointsPerProc - ((rank + 1) * pointsPerProc) % k_neighbours;

            sendPoints = (double *) malloc(pointsToGive * dimension * sizeof(double));

            for (int i = 0; i < pointsToGive; i++) {
                memcpy(sendPoints + i * dimension, holdPoints[i], dimension * sizeof(double));
            }

            MPI_Send(sendPoints, (int) (pointsToGive * dimension), MPI_DOUBLE, 0, 12, MPI_COMM_WORLD);
        }

        free(initial.vpPoint);
        free(sendPoints);
    }
    else{
        return;
    }
}