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

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <stdio.h>
#include <string.h>
#include "benchmarkSetup.h"

/* function prototypes */
int printHeader();
int setParallelInfo(int numMPIprocs, int threadSupport, int numThreads);
int setBenchName(char *name, int number, int support);
int printNodeReport(int sameNode, int rankA, int rankB);
int printBenchHeader();
int setTestOutcome(int outcome);
int setReportParams(int size, int reps, double time);
int printBenchName();
int printReport();
int threadSupportToString(int threadSupport, char *string);
int printMultiProcInfo(int printNode, int pairWorldRank, char *pairProcName);
int printBalanceError();

/* define report data type */
struct report
     {
		char benchName[30];
		int benchNumber;
		int supported;
		int dataSize;
		int bytes;
		char testOutcome[5];
		int numReps;
		double benchTime;
		double timePerRep;
		/* parallel environment info */
		int numMPIprocs;
		int numThreads;
		int supportLevel;
     };

/* variable declaration */
struct report benchReport;

#endif /* OUTPUT_H_ */
