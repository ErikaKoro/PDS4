#include "mpi_exchanges.h"
#include <stdlib.h>

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
int *findExchanges(int rank, int *counterReceiver, int worldSize, int master, MPI_Comm communicator, int *numberOfExchanges){
    int *helperIndex;

    /// IN THE EVEN INDEXES WILL BE STORED THE NUMBER OF POINTS THAT WILL SEND IN THE RANK THAT LIES IN THE ODD INDEXES
    int **info;

    if(rank == master) {
        // Allocate memory for the array that will contain with which process should the current process communicate and how many points will they exchange
        info = (int **) calloc(worldSize, sizeof(int *));      // FREEMEM (CHECK)
        for (int i = 0; i < worldSize; i++) {
            info[i] = (int *) calloc(2 * worldSize,  sizeof(int));
        }

        // This array holds the index for its line in the info array in which we should write the next element
        helperIndex = (int *) calloc(worldSize, sizeof(int));       // FREEMEM (CHECK)

        int nextIndex = worldSize / 2;  // Indicates the index of the second half of the array
        for (int i = 0; i < worldSize / 2; i++) {

            while (counterReceiver[i] > 0) {  // Repeat until each process gets the info for its exchanges. When counterReceiver[i] = 0 the rank doesn't want to exchange anymore or not at all

                if (counterReceiver[i] == counterReceiver[nextIndex]) {  // Compare elements of the first half of array to the second half's elements
//                    printf("Case 1 equals Rank %d will send %d points to the Rank %d\n", i, counterReceiver[i], nextIndex);
                    info[i][helperIndex[i]] = nextIndex;  // Store for the rank i that it should give points to the nextIndex(it lies in the second half of the array as it is initialized with the value worldSize/2

                    helperIndex[i]++;  // Increase the index which indicates the last index we wrote sth about the i rank

                    info[i][helperIndex[i]] = counterReceiver[i];  // In the new index store how many points will the i rank send to the newIndex rank

                    helperIndex[i]++;

                    info[nextIndex][helperIndex[nextIndex]] = i;  // Store the info exchanges for the nextIndex rank as well
                    // Send to the rank i

                    helperIndex[nextIndex]++;  // Increase the index which indicates the last index we wrote sth about the nextIndex rank

                    info[nextIndex][helperIndex[nextIndex]] = counterReceiver[i];  // Send to rank i as many elements as it sent to rank nextIndex

                    helperIndex[nextIndex]++;

                    counterReceiver[i] = 0; // The i rank has gotten the info for its exchanges

                    counterReceiver[nextIndex] = 0;  // The nextIndex has gotten the info for its exchanges too

                    nextIndex++;  // Point to the new index we' re going to write
                } else if (counterReceiver[i] < counterReceiver[nextIndex]) {  // If the current rank wants to give less elements than the nextIndex requires
//                    printf("Case 2 sender has less Rank %d will send %d points to the Rank %d\n", i, counterReceiver[i], nextIndex);
                    info[i][helperIndex[i]] = nextIndex;

                    helperIndex[i]++;

                    info[i][helperIndex[i]] = counterReceiver[i];

                    helperIndex[i]++;

                    info[nextIndex][helperIndex[nextIndex]] = i;

                    helperIndex[nextIndex]++;

                    info[nextIndex][helperIndex[nextIndex]] = counterReceiver[i];

                    helperIndex[nextIndex]++;

                    counterReceiver[nextIndex] -= counterReceiver[i]; // the points remained that the nextIndex rank should exchange to finish its exchanges

                    counterReceiver[i] = 0; // The i rank has gotten all the info
                } else if (counterReceiver[i] > counterReceiver[nextIndex] && counterReceiver[nextIndex] != 0) {
//                    printf("Case 3 sender has more Rank %d will send %d points to the Rank %d\n", i, counterReceiver[nextIndex], nextIndex);
                    info[i][helperIndex[i]] = nextIndex;

                    helperIndex[i]++;

                    info[i][helperIndex[i]] = counterReceiver[nextIndex];

                    helperIndex[i]++;

                    info[nextIndex][helperIndex[nextIndex]] = i;

                    helperIndex[nextIndex]++;


                    info[nextIndex][helperIndex[nextIndex]] = counterReceiver[nextIndex];

                    helperIndex[nextIndex]++;

                    counterReceiver[i] -= counterReceiver[nextIndex];  // The points left in i rank for which we should find with which next process should they be exchanged

                    counterReceiver[nextIndex] = 0;  // The nextIndex has gotten all its info

                    nextIndex++;
                } else {
                    nextIndex++;  // In case we have counterReceiver[i] > counterReceiver[nextIndex] && counterReceiver[nextIndex] == 0
                }
            }
        }

        // Allocate memory for the request objects
        MPI_Request *requests = (MPI_Request *) malloc((worldSize - 1) * sizeof (MPI_Request));     // FREEMEM (CHECK)

        // For each rank in the worldSize send a request from the master to the other processes
        for (int i = 1; i < worldSize; i++) {
            MPI_Isend(info[i], helperIndex[i], MPI_INT, i, 10, communicator, &requests[i - 1]);
        }

//        printf("\nThe rank %d has receive buffer: ", rank);
//        for(int i = 0; i < helperIndex[0]; i++){
//            printf("%d ", info[0][i]);
//        }
//        printf("\n");

        for(int i = 0; i < worldSize - 1; i++){
            MPI_Wait(&requests[i], MPI_STATUS_IGNORE);  // Wait for the requests to finish
        }
        *numberOfExchanges = helperIndex[0];  // The length of the array for the infos of the rank master

        free(helperIndex);

        // TODO possible snake
        for (int i = 1; i < worldSize; i++) {
            free(info[i]);
        }

        return info[0];

    } else {
        // Allocate memory for each process for the info it will receive
        int *infoPerProc;

        int elements;  // The number of elements that have been requested to be sent
        MPI_Status status;
        MPI_Probe(master, 10, communicator, &status);  // Ask MPI to give you the size of the message

        MPI_Get_count(&status, MPI_INT, &elements);

        infoPerProc = (int *) calloc(elements, sizeof (int));  // Allocate memory for the receiver buffer of the processes that will get infos from the master     // FREEMEM (CHECK)
        MPI_Recv(infoPerProc, elements, MPI_INT, master, 10, communicator, MPI_STATUS_IGNORE);  // Receive the array per process!

//        printf("\nThe rank %d has receive buffer: ", rank);
//        for(int i = 0; i < elements; i++){
//            printf("%d ", infoPerProc[i]);
//        }
//        printf("\n");
        *numberOfExchanges = elements;

        return infoPerProc;
    }
}