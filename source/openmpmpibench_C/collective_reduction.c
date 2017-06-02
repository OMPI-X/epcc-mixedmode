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
/* Implements the collective reduce and allreduce mixed      */
/* mode OpenMP/MPI benchmarks.                               */
/*-----------------------------------------------------------*/
#include "collective_reduction.h"

/*-----------------------------------------------------------*/
/* reduction                                                 */
/*                                                           */
/* Driver subroutine for the reduce and allReduce            */
/* benchmarks.                                               */
/*-----------------------------------------------------------*/
int reduction(int benchmarkType){
	int dataSizeIter, sizeofBuf;

	/* Initialise repsToDo to defaultReps */
	repsToDo = defaultReps;

	/* Start loop over data sizes */
	dataSizeIter = minDataSize; /* initialise dataSizeIter */
	while (dataSizeIter <= maxDataSize){
		/* allocate space for the main data arrays.. */
		allocateReduceData(dataSizeIter);

		/* Perform benchmark warm-up */
		if (benchmarkType == REDUCE){
			reduceKernel(warmUpIters, dataSizeIter);
			/* Master process tests if reduce was a success */
			if (myMPIRank == 0){
				testReduce(dataSizeIter, benchmarkType);
			}
		}
		else if (benchmarkType == ALLREDUCE){
			/* calculate sizeofBuf for test */
			sizeofBuf = dataSizeIter * numThreads;
			allReduceKernel(warmUpIters, dataSizeIter);
			/* all processes need to perform unit test */
			testReduce(sizeofBuf, benchmarkType);
		}

		/* Initialise the benchmark */
		benchComplete = FALSE;
		/* Execute benchmark until target time is reached */
		while (benchComplete != TRUE){
			/* Start timer */
			MPI_Barrier(comm);
			startTime = MPI_Wtime();

			/* Execute reduce for repsToDo repetitions */
			if (benchmarkType == REDUCE){
				reduceKernel(repsToDo, dataSizeIter);
			}
			else if (benchmarkType == ALLREDUCE){
				allReduceKernel(repsToDo, dataSizeIter);
			}

			/* Stop timer */
			MPI_Barrier(comm);
			finishTime = MPI_Wtime();
			totalTime = finishTime - startTime;

			/* Test if target time was reached with the number of reps */
			if (myMPIRank==0){
			  benchComplete = repTimeCheck(totalTime, repsToDo);
			}
			/* Ensure all procs have the same value of benchComplete */
			/* and repsToDo */
			MPI_Bcast(&benchComplete, 1, MPI_INT, 0, comm);
			MPI_Bcast(&repsToDo, 1, MPI_INT, 0, comm);
		}

		/* Master process sets benchmark result for reporting */
		if (myMPIRank == 0){
			setReportParams(dataSizeIter, repsToDo, totalTime);
			printReport();
		}

		/* Free allocated data */
		freeReduceData();

		/* Double dataSize and loop again */
		dataSizeIter = dataSizeIter * 2;

	}

	return 0;
}

/*-----------------------------------------------------------*/
/* reduceKernel                                              */
/*                                                           */
/* Implements the reduce mixed mode benchmark.               */
/* Each thread under every MPI process combines its local    */
/* buffer. This is then sent to the master MPI process to    */
/* get the overall reduce value.                             */
/*-----------------------------------------------------------*/
int reduceKernel(int totalReps, int dataSize){
	int repIter, i, j;

	for (repIter=1; repIter<totalReps; repIter++){

		/* Manually perform the reduction between OpenMP threads */
#pragma omp parallel default(none) \
	private(i,j) \
	shared(tempBuf,globalIDarray,dataSize,numThreads) \
	shared(localReduceBuf)
		{
		/* 1) Intialise the tempBuf array */
#pragma omp for schedule(static,dataSize)
		for(i=0; i<(numThreads * dataSize); i++){
			tempBuf[i] = globalIDarray[myThreadID];
		}
		/* 2) Reduce tempBuf into localReduceBuf */
#pragma omp for
		for(i=0; i<dataSize; i++){
			localReduceBuf[i] = 0;
			for (j=0; j<numThreads; j++){
				localReduceBuf[i] += tempBuf[(j*dataSize)+i];
			}
		}
		}
		/* Perform a reduce of localReduceBuf across the
		 * MPI processes.
		 */
		MPI_Reduce(localReduceBuf, globalReduceBuf, dataSize,\
				MPI_INT, MPI_SUM, 0, comm);

		/* Copy globalReduceBuf into master Threads portion
		 * of finalReduceBuf.
		 */
		// FR this should only happen on rank == 0 thus added if to ensure this
		if (myMPIRank==0) {
		        for (i=0; i<dataSize; i++){
			        finalReduceBuf[i] = globalReduceBuf[i];
			}
		}

	}
	return 0;
}

/*-----------------------------------------------------------*/
/* allReduce                                                 */
/*                                                           */
/* Implements the allreduce mixed mode benchmark.            */
/* Each thread under every MPI process combines its local    */
/* buffer. All MPI processes then combine this value to      */
/* the overall reduction value at each process.              */
/*-----------------------------------------------------------*/
int allReduceKernel(int totalReps, int dataSize){
	int repIter, i, j;
	int startPos;


	for (repIter=0; repIter<totalReps; repIter++){

		/* Manually perform the reduction between OpenMP threads */
#pragma omp parallel default(none) \
	private(i,j) \
	shared(tempBuf,globalIDarray,dataSize,numThreads) \
	shared(localReduceBuf)
		{
		/* 1) Intialise the tempBuf array */
#pragma omp for schedule(static,dataSize)
		for(i=0; i<(numThreads * dataSize); i++){
			tempBuf[i] = globalIDarray[myThreadID];
		}
		/* 2) Reduce tempBuf into localReduceBuf */
#pragma omp for
		for(i=0; i<dataSize; i++){
			localReduceBuf[i] = 0;
			for (j=0; j<numThreads; j++){
				localReduceBuf[i] += tempBuf[(j*dataSize)+i];
			}
		}
		}

		/* Perform an all reduce of localReduceBuf across
		 * the MPI processes.
		 */
		 MPI_Allreduce(localReduceBuf, globalReduceBuf, \
				 dataSize, MPI_INTEGER, MPI_SUM, comm);

		 /* Each thread copies globalReduceBuf into its portion
		  * of finalReduceBuf.
		  */
#pragma omp parallel default(none) \
	private(i,startPos) \
	shared(dataSize,finalReduceBuf,globalReduceBuf)
		 {
		 /* Calculate the start of each threads portion
		  * of finalReduceBuf.
		  */
		 startPos = (myThreadID * dataSize);
		 for (i=0; i<dataSize; i++){
			 finalReduceBuf[startPos + i] = globalReduceBuf[i];
		 }
		 }

	}
	
	return 0;
}

/*-----------------------------------------------------------*/
/* allocateReduceData                                        */
/*                                                           */
/* Allocate memory for the main data arrays in the           */
/* reduction operation.                                      */
/*-----------------------------------------------------------*/
int allocateReduceData(int bufferSize){

	localReduceBuf = (int *) malloc(bufferSize * sizeof(int));
	globalReduceBuf = (int *) malloc(bufferSize * sizeof(int));
	/* tempBuf and Final reduce is of size dataSize*numThreads */
	tempBuf = (int *) malloc((bufferSize * numThreads) * sizeof(int));
	finalReduceBuf = (int *) malloc((bufferSize * numThreads) * sizeof(int));

	return 0;
}

/*-----------------------------------------------------------*/
/* freeReduceData                                            */
/*                                                           */
/* Free allocated memory for main data arrays.               */
/*-----------------------------------------------------------*/
int freeReduceData(){

	free(localReduceBuf);
	free(globalReduceBuf);
	free(tempBuf);
	free(finalReduceBuf);
	return 0;
}

/*-----------------------------------------------------------*/
/* testReduce                                                */
/*                                                           */
/* Verifies that the reduction benchmarks worked correctly.  */
/*-----------------------------------------------------------*/
int testReduce(int bufferSize, int benchmarkType){
	int i, testFlag, reduceFlag;
	int correctReduce, lastGlobalID;

	/* Initialise correctReduce to 0.. */
	correctReduce = 0;
	/* ..and testFlag to true */
	testFlag = TRUE;

	/* set lastGlobalID */
	lastGlobalID = (numMPIprocs * numThreads);

	/* Now find correctReduce value by summing to lastGlobalID */
	for (i=0; i<lastGlobalID; i++){
		correctReduce = correctReduce + i;
	}

	/* Compare each element of finalRecvBuf to correctRedcue */
	for (i=0; i<bufferSize; i++){
		if (finalReduceBuf[i] != correctReduce){
			testFlag = FALSE;
		}
	}

	/* For allReduce, combine testFlag into master
	 * with logical AND.
	 */
	if (benchmarkType == ALLREDUCE){
		MPI_Reduce(&testFlag, &reduceFlag, 1, MPI_INT, MPI_LAND, 0, comm);
		/* then master sets testOutcome using reduceFlag */
		if (myMPIRank == 0){
			setTestOutcome(reduceFlag);
		}
	}
	else{
		/* for reduce master process just sets testOurcome using testFlag */
		setTestOutcome(testFlag);
	}

	return 0;
}
