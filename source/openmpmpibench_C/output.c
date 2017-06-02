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
/* Routines for output of benchmark timings.                 */
/*-----------------------------------------------------------*/
#include "output.h"

/*-----------------------------------------------------------*/
/* printHeader                                               */
/*                                                           */
/* Prints a header in the output.                            */
/*-----------------------------------------------------------*/
int printHeader(){
	char string[MAXSTRING];

    /* Convert threadSupport to a string for output */
    threadSupportToString(benchReport.supportLevel, string);

    printf("----------------------------------------------\n");
    printf("Mixed mode MPI/OpenMP benchmark suite v1.0\n");
    printf("----------------------------------------------\n");
    printf("Number of MPI processes = %d\n", benchReport.numMPIprocs);
    printf("Number of OpenMP threads = %d\n", benchReport.numThreads);
    printf("Thread support = %s\n", string);

    return 0;
}

/*-----------------------------------------------------------*/
/* setParallelInfo                                           */
/*                                                           */
/* Sets the number of MPI processes, the OpenMP thread count */
/* and the level of thread support provided by the MPI       */
/* implementation.                                           */
/*-----------------------------------------------------------*/
int setParallelInfo(int numMPIprocs, int threadSupport, int numThreads){
	/* set the number of MPI processes in report data type*/
	benchReport.numMPIprocs = numMPIprocs;
	/* set thread support level in report data type */
	benchReport.supportLevel = threadSupport;
	/* set number of OpenMP threads in report data type */
	benchReport.numThreads = numThreads;

	return 0;
}

/*-----------------------------------------------------------*/
/* setBenchName                                              */
/* Sets the benchName, benchNumber and if the benchmark is   */
/* supported.                                                */
/*-----------------------------------------------------------*/
int setBenchName(char *name, int number, int support){
	strcpy(benchReport.benchName,name);
	benchReport.benchNumber = number;
	benchReport.supported = support;

	printBenchName();

	return 0;
}

/*-----------------------------------------------------------*/
/* printBenchName                                            */
/*                                                           */
/* Print header for benchmark - name of benchmark and        */
/* list of names of each column.                             */
/*-----------------------------------------------------------*/
int printBenchName(){

	printf("--------------------------------------------\n");
	printf("# %s\n", benchReport.benchName);
	printf("--------------------------------------------\n");

	if (benchReport.supported == FALSE){
		printf("WARNING: Implementation does not support benchmark.\n");
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* printNodeReport                                           */
/* 															 */
/* For pingpong and pingping benchmarks prints out if the    */
/* two MPI processes are on the same node or not.            */
/*-----------------------------------------------------------*/
int printNodeReport(int sameNode, int rankA, int rankB){
	if (sameNode == TRUE){
		printf("Intra node benchmark between process %d and process %d\n",rankA,rankB);
	}
	else if (sameNode == FALSE){
		printf("Inter node benchmark between process %d and process %d\n",rankA,rankB);
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* printBenchHeader                                          */
/* 															 */
/* Prints the column headings for the benchmark report.      */
/*-----------------------------------------------------------*/
int printBenchHeader(){

	printf(" Data Size     Msg Size (bytes)     No. Reps     ");
	printf("Time (sec)     Time/Rep (s)     Test\n");

	printf("-----------   ------------------   ----------   ");
	printf("------------   --------------   ------\n");

	return 0;
}

/*-----------------------------------------------------------*/
/* setTestOutcome                                            */
/*															 */
/* Sets benchReport's testOutcome element.                   */
/* Called in test routine of each benchmark.                 */
/*-----------------------------------------------------------*/
int setTestOutcome(int outcome){

	if (outcome == TRUE){
		strcpy(benchReport.testOutcome,"Pass");
	}
	else if (outcome == FALSE){
		strcpy(benchReport.testOutcome,"Fail");
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* setReportParams                                           */
/*															 */
/* Sets the numReps and benchTime for a certain datasize     */
/*-----------------------------------------------------------*/
int setReportParams(int size, int reps, double time){
	benchReport.dataSize = size;
	benchReport.numReps = reps;
	benchReport.benchTime = time;

	/* Calculate and set time for 1 rep */
	benchReport.timePerRep = time / reps;
	/* Calculate the size of message in bytes */
	if (benchReport.benchNumber <= LAST_PT_PT_ID){
		/* if point to point benchmark size of msg is
		 * dataSize x numThreads x sizeof(int)
		 */
		benchReport.bytes = size * benchReport.numThreads * sizeInteger;
	}
	else if (benchReport.benchNumber <= LASTMULTIPPID){
		benchReport.bytes = size * benchReport.numThreads * sizeInteger * localCommSize;
	}
	else {
		benchReport.bytes = size * sizeInteger;
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* printMultiProcInfo                                        */
/*                                                           */
/* This prints the comm world ranks and processor names for  */
/* each pair of processes in the multi-pingpong or           */
/* multi-pingping benchmarks.                                */
/*-----------------------------------------------------------*/
int printMultiProcInfo(int printNode, int pairWorldRank, char *pairProcName){

	/* MPI Processes under printNode of crossComm perform the output */
	if (crossCommRank == printNode){
		printf("MPI process %d on %s ", myMPIRank,myProcName);
		printf("communicating with MPI process %d on %s\n",\
				pairWorldRank,pairProcName);
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* printReport                                               */
/*															 */
/* Prints out the a column of information after each         */
/* data size iteration.                                      */
/*-----------------------------------------------------------*/
int printReport(){

	/* Write to output */
	printf("d %d\t\t%d\t\t   %d\t\t%lf\t%lf\t%s\n", \
			benchReport.dataSize, benchReport.bytes, benchReport.numReps, \
	        benchReport.benchTime, benchReport.timePerRep,benchReport.testOutcome);

	return 0;
}

/*-----------------------------------------------------------*/
/* printBalanceError                                         */
/*                                                           */
/* Prints an error if there isn't the same number of MPI     */
/* processes in the nodes selected for the multi-pingpong    */
/* or multi-pingping benchmarks.                             */
/*-----------------------------------------------------------*/
int printBalanceError(){

	printf("\nERROR: Nodes selected for this benchmark do not ");
	printf("have same number of MPI processes per node.");
	printf("Skipping benchmark...\n");

	return 0;
}

/*-----------------------------------------------------------*/
/* threadSupportToString                                     */
/*															 */
/* Converts the threadSupport integer variable to a          */
/* string for output.                                        */
/*-----------------------------------------------------------*/
int threadSupportToString(int threadSupport, char *string){

	if (threadSupport == MPI_THREAD_SINGLE){
		strcpy(string,"MPI_THREAD_SINGLE");
	}
	else if (threadSupport == MPI_THREAD_FUNNELED){
		strcpy(string,"MPI_THREAD_FUNNELED");
	}
	else if (threadSupport == MPI_THREAD_SERIALIZED){
			strcpy(string,"MPI_THREAD_SERIALIZED");
	}
	else if (threadSupport == MPI_THREAD_MULTIPLE){
				strcpy(string,"MPI_THREAD_MULTIPLE");
	}

	return 0;
}


