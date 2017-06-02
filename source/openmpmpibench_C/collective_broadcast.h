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

#ifndef COLLECTIVE_BROADCAST_H_
#define COLLECTIVE_BROADCAST_H_

#include "benchmarkSetup.h"
#include "output.h"
#include "parallelEnvironment.h"
#include <stdio.h>

#define BROADCASTNUM 100
#define BROADCASTROOT 0 /* Source of broadcast */


/* function prototypes */
int broadcast();
int broadcastKernel(int totalReps, int dataSize);
int allocateBroadcastData(int bufferSize);
int freeBroadcastData();
int testBroadcast(int bufferSize);

/* variable declaration */
int *broadcastBuf;
int *finalBroadcastBuf;

#endif /* COLLECTIVE_BROADCAST_H_ */
