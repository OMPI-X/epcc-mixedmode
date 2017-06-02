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

#ifndef PT_TO_PT_PINGPING_H_
#define PT_TO_PT_PINGPING_H_

#include "parallelEnvironment.h"
#include "benchmarkSetup.h"
#include "output.h"
#include <stdio.h>

/* function prototypes */
int pingPing(int benchmarkType);
int allocatePingpingData(int sizeofBuffer);
int freePingpingData();
int masteronlyPingping(int totalReps, int dataSize);
int funnelledPingping(int totalReps, int dataSize);
int multiplePingping(int totalReps, int dataSize);
int testPingping(int sizeofBuffer, int dataSize);

/* variable declaration */
int pingRankA, pingRankB;
int sizeofBuffer;
int *pingSendBuf, *pingRecvBuf;
int *finalRecvBuf;
int *testBuf;

#endif /* PT_TO_PT_PINGPING_H_ */
