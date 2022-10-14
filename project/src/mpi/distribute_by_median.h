#ifndef PDS4_DISTRIBUTE_BY_MEDIAN_H
#define PDS4_DISTRIBUTE_BY_MEDIAN_H
#include <mpi.h>


/**
 * Finds the distance of each element, per process, from the pivot
 *
 * @param dist array to hold the results
 * @param points points per process with their coordinates
 * @param dimension the dimension of our space (e.g 3D)
 * @param pivot array with the coordinates of the pivot
 * @param pointsPerProc
 */
void findMPIDistances(int rank, double *dist, double **points, int dimension, const double *pivot, int pointsPerProc);

/**
 * Calculates the distances, gathers them to the master, which finds their median, then broadcasts it again to the
 * other processes in order to calculate how many points do they want to exchange and sort their array with their points based on their distances from the median. After that,
 * the master gathers from each process the number of points each one of them wants to exchange. The master based on the array with the counters finds which process should communicate with
 * with which process in order to the processes with rank < worldSize/2 have distances(from the pivot) < median and the processes with rank >= worldSize/2 have
 * distances(from the pivot) > median. Continuously, a buffer info is sent to each process with which one should communicate aiming at the completion of the exchanges.
 * Finally, the processes communicate and exchange their points.
 *
 * @param master master rank
 * @param rank process id
 * @param dimension
 * @param holdPoints the array with the points and their coordinates
 * @param pointsPerProc
 * @param worldSize the number of processes
 * @param communicator "holds" the processes with which the function will be called with, recursively
 */
void distributeByMedian(double *pivot,int master, int rank, int dimension, double **holdPoints, int pointsPerProc, int worldSize, MPI_Comm communicator);


#endif //PDS4_DISTRIBUTE_BY_MEDIAN_H
