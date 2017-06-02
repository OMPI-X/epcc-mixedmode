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
/* Header file for benchmarkSetup.c.                         */
/*-----------------------------------------------------------*/

#ifndef BENCHMARKSETUP_H_
#define BENCHMARKSETUP_H_
#include "parallelEnvironment.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FALSE 0
#define TRUE 1

#define NUM_BENCHMARKS 22
#define MAXSTRING 30

#define FINISHED 999 
#define ERROR 100

#define LASTPPID 5 /* id of last pingpong/pingping bench */
#define LAST_PT_PT_ID 8 /* id of last point-to-point bench */
#define LASTMULTIPPID 14 /* id of last multi pt-to-pt bench */

/*Benchmark types */
#define MASTERONLY 1
#define FUNNELLED  2
#define MULTIPLE 3
#define REDUCE 4
#define ALLREDUCE 5
#define SCATTER 6
#define GATHER 7

/* function prototypes */
int openFile(char *fileName);
int closeFile();
int setupBenchmarkList();
int readBenchmarkParams();
int findBenchmarkNumber();
int convertToLowercase(char *convertString);
int repTimeCheck(double time, int numReps);
int max(int a, int b);

/* variable declaration */
FILE *inputFile;
char benchmarkList[NUM_BENCHMARKS][MAXSTRING];

int warmUpIters; /* number of iterations of warmup loop */
int defaultReps; /* the default number of repetitions */
int repsToDo; /* reps to do for a benchmark */
int minDataSize;
int maxDataSize;
double targetTime; /* threshold time for benchmark */

int benchComplete;
int benchmarkNumber;

/* variables for  timing */
double startTime,finishTime, totalTime;



#endif /* BENCHMARKSETUP_H_ */
