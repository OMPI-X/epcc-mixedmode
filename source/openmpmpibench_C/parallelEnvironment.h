/*****************************************************************************
 *                                                                           *
 *             Mixed-mode OpenMP/MPI MicroBenchmark Suite - Version 1.0      *
 *                                                                           *
 *                            produced by                                    *
 *                                                                           *
 *                Mark Bull, Jim Enright and Fiona Reid                      *
 *                                                                           *
 *                                at                                         *
 *                                                                           *
 *                Edinburgh Parallel Computing Centre                        *
 *                                                                           *
 *   email: markb@epcc.ed.ac.uk, fiona@epcc.ed.ac.uk                         *
 *                                                                           *
 *                                                                           *
 *              Copyright 2012, The University of Edinburgh                  *
 *                                                                           *
 *                                                                           *
 *  Licensed under the Apache License, Version 2.0 (the "License");          *
 *  you may not use this file except in compliance with the License.         *
 *  You may obtain a copy of the License at                                  *
 *                                                                           *
 *      http://www.apache.org/licenses/LICENSE-2.0                           *
 *                                                                           *
 *  Unless required by applicable law or agreed to in writing, software      *
 *  distributed under the License is distributed on an "AS IS" BASIS,        *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 *  See the License for the specific language governing permissions and      *
 *  limitations under the License.                                           *
 *                                                                           *
 ****************************************************************************/

/*-----------------------------------------------------------*/
/* Header file for parallelEnvironment.c.                    */
/* Contains variables declarations and function prototypes   */
/* used to setup the MPI and OpenMP programming environment. */
/*-----------------------------------------------------------*/

#ifndef PARALLELENVIRONMENT_H_
#define PARALLELENVIRONMENT_H_

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "benchmarkSetup.h"

/* function prototypes */
int initParallelEnv();
int finaliseParallelEnv();
int findRank(int rankIn);
int findNeighbours();
int compareProcNames(int rankA, int rankB);
int setupCommunicators();
int procNameToHash();
int exchangeWorldRanks(int nodeA, int nodeB, int *otherWorldRank);
int sendProcName(int destNode, int srcNode, char *destProcName);
int crossCommBalance(int nodeA, int nodeB);

/* variable declaration */
/*MPI variables */
#define TAG 1 /* set tag to match messages */
int myMPIRank, numMPIprocs;
MPI_Comm comm, commCart;
MPI_Comm crossComm, localComm;
int localCommRank, localCommSize, crossCommRank;
MPI_Status status;
char myProcName[MPI_MAX_PROCESSOR_NAME];
int procNameLen;
MPI_Request requestID;
MPI_Status statusArray[4]; /* for haloexchange */
MPI_Request requestArray[4]; /* for haloexchange */
int leftNeighbour, rightNeighbour;
int sizeInteger;

int PPRanks[2]; /* ranks for pingpong or pingping */

/* OpenMP variables */
int myThreadID, numThreads;
/* make myThreadID a thread private variable */
#pragma omp threadprivate(myThreadID)
/*Array to hold the global ID for each thread */
int *globalIDarray;

int threadSupport;
#endif /* PARALLELENVIRONMENT_H_ */
