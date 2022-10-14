#ifndef PDS4_MPI_EXCHANGES_H
#define PDS4_MPI_EXCHANGES_H
#include <mpi.h>

/**
 *
 * This function finds the exchanges between the processes according to the following algorithm:
 * counterReceiver: The array in which the master has gathered the number of points each process wants to exchange
 * The main idea is to compare one by one the processes with ranks smaller than worldSize/2 to the processes with ranks bigger than worldSize/2,
 * based on how many points each one wants to exchange. For example, let's consider the array with processes' ids
 * and the counterReceiver array to be:
 * ranks = [0 1 2 3 ]
 * counterReceiver = [2 3 5 0]
 * The aim is to give the points of the first half of the array to the second half. So, in the beginning we compare the counterReceiver[0] = 2 to
 * counterReceiver[0 + worldSize/2] = counterReceiver[0 + 4/2] = counterReceiver[2] = 5.
 * We see that 2 < 5 which means process with rank 0 should give to process with rank 2 2 points and vice versa.
 * The process with rank 0 has completed its exchanges in contrast to the process with rank 2. The last one wants to give 5 - 2 = 3 more points to another process.
 * So, now we have to compare the (i + 1)th element of the counterReceiver to the left points of the incomplete . We observe that counterReceiver[1] = 3
 * is equal to the number of points that has remained in the process with rank 2. So, the process with rank 1 should give to process with rank 2
 * 3 points and take from it 3 points as well.
 * Then, both processes have completed their exchanges.This algorithm continues till all processes have exchanged as many elements as they wanted.
 *
 * @param rank
 * @param counterReceiver The array with the number of points, each process wants to give, that has bees sent to the master
 * @param worldSize
 * @param master
 * @param communicator the cluster of processes that calls the function
 * @param numberOfExchanges points to the length of the array with the infos for each rank
 * @return the array with the infos per process
 *
 */
int *findExchanges(int rank, int *counterReceiver, int worldSize, int master, MPI_Comm communicator, int *numberOfExchanges);

#endif //PDS4_MPI_EXCHANGES_H
