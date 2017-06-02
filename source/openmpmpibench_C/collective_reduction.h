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

#ifndef COLLECTION_REDUCTION_H_
#define COLLECTION_REDUCTION_H_

#include "benchmarkSetup.h"
#include "output.h"
#include "parallelEnvironment.h"
#include <stdio.h>

/* function prototypes */
int reduction(int benchmarkType);
int reduceKernel(int totalReps, int dataSize);
int allReduceKernel(int totalReps, int dataSize);
int allocateReduceData(int bufferSize);
int freeReduceData();
int testReduce(int bufferSize, int benchmarkType);

/* variable declaration */
int *localReduceBuf, *globalReduceBuf;
int *tempBuf, *finalReduceBuf;

#endif /* COLLECTION_REDUCTION_H_ */
