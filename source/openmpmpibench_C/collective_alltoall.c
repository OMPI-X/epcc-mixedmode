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
/* Implements the alltoall mixed mode OpenMP/MPI benchmark.  */
/*-----------------------------------------------------------*/
#include "collective_alltoall.h"

/*-----------------------------------------------------------*/
/* alltoall                                                  */
/*                                                           */
/* Driver routine for the alltoall benchmark.                */
/*-----------------------------------------------------------*/
int alltoall(){
	int dataSizeIter;
	int bufferSize;

	/* Initialise repsToDo to defaultReps */
	repsToDo = defaultReps;

	/* Start loop over data sizes */
	dataSizeIter = minDataSize; /* initialise dataSizeIter */
	while (dataSizeIter <= maxDataSize){
		/* Calculate bufferSize and allocate space for
		 * the data arrays.
		 */
		bufferSize = dataSizeIter * numThreads * \
			numMPIprocs * numThreads;

		allocateAlltoallData(bufferSize);

		/* Perform warm-up of benchmark */
		alltoallKernel(warmUpIters, dataSizeIter);

		/* Test if alltoall was successful */
		testAlltoall(dataSizeIter);

		/* Initialise the benchmark */
		benchComplete = FALSE;

		/* Execute benchmark until target time is reached */
		while (benchComplete != TRUE){
			/* Start timer */
			MPI_Barrier(comm);
			startTime = MPI_Wtime();

			/* Execute alltoall for repsToDo repetitions */
			alltoallKernel(repsToDo, dataSizeIter);

			/* Stop timer */
			MPI_Barrier(comm);
			finishTime = MPI_Wtime();
			totalTime = finishTime - startTime;

			/* Test if target time was reached */
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
		freeAlltoallData();

		/* Double data size and loop again */
		dataSizeIter = dataSizeIter * 2;
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* alltoallKernel                                            */
/*                                                           */
/* Implements the all to all benchmark.                      */
/* Each thread sends/receives dataSize items to/from         */
/* every other process.                                      */
/*-----------------------------------------------------------*/
int alltoallKernel(int totalReps, int dataSize){
	int repIter, i, j;
	int dataForEachProc, numsToWrite;
	int blockNum, startOffset;

	/* Calculate how much data each thread sends to each process */
	numsToWrite = numThreads * dataSize;
	/* Calculate total amount of data each process receives
	 * from any other process....
	 * ...each thread gets dataSize items from every other thread.
	 */
	dataForEachProc = numThreads * numThreads * dataSize;

	for (repIter=0; repIter<totalReps; repIter++){

		/* Each thread writes numsToWrite items for each
		 * MPI process to alltoallSendBuf.
		 */
#pragma omp parallel default(none) \
	private(blockNum,i,j) \
	shared(numsToWrite,dataForEachProc,globalIDarray) \
	shared(alltoallSendBuf,numMPIprocs)
		{
			/* Calculate the blockNum of each thread.
			 * This is used to find which portion of the
			 * dataForEachProc elements a thread will
			 * be responsible for.
			 */
			blockNum = (myThreadID)* numsToWrite;

			/* Write threadID to correct location in
			 * alltoallSendBuf.
			 */
			for (i=0; i<numMPIprocs; i++){ /* loop over MPI processes */
				for (j=0; j<numsToWrite; j++){ /* loop over data to write */
					alltoallSendBuf[blockNum +(i * dataForEachProc) + j] = \
						globalIDarray[myThreadID];
				}
			}
		}

		/* Call MPI_AlltoAll */
		MPI_Alltoall(alltoallSendBuf, dataForEachProc, MPI_INT, \
		            alltoallRecvBuf, dataForEachProc, MPI_INT, \
		            comm);

		/* Each thread now reads the receive buffer so that it gets
		 * dataSize values from every other thread in its portion
		 * of alltoallFinalBuf.
		 */
#pragma omp parallel default(none) \
	private(blockNum,startOffset,i,j) \
	shared(alltoallRecvBuf,alltoallFinalBuf,numMPIprocs) \
	shared(dataForEachProc,numsToWrite,dataSize,globalIDarray) \
	shared(numThreads)
		{
			/* Calculate the blockNum.
			 * This identifies which portion of the data from each
			 * process a thread is responsible for.
			 */
			blockNum = myThreadID * dataSize;

			/* Calculate the offset into each MPI processes finalBuf
			 * where each thread will start to read its data.
			 */
			startOffset = (numsToWrite * numMPIprocs) * myThreadID;

			/* Loop over all processors (threads & proceeses). */
			for (i=0; i<(numThreads * numMPIprocs); i++){
				for (j=0; j<dataSize; j++){
					alltoallFinalBuf[startOffset + (i * dataSize) + j] = \
						alltoallRecvBuf[blockNum + (i * numsToWrite) + j];
				}
			}

		}

	}
	return 0;
}

/*-----------------------------------------------------------*/
/* allocateAlltoallData                                      */
/*                                                           */
/* Allocates memory for the main data arrays used in the     */
/* alltoall benchmark.                                       */
/*-----------------------------------------------------------*/
int allocateAlltoallData(int bufferSize){

	alltoallSendBuf = (int *) malloc(bufferSize * sizeof(int));
	alltoallRecvBuf = (int *) malloc(bufferSize * sizeof(int));
	alltoallFinalBuf = (int *) malloc(bufferSize * sizeof(int));

	return 0;
}

/*-----------------------------------------------------------*/
/* freeAlltoallData                                          */
/*                                                           */
/* Free memory of the main data arrays.                      */
/*-----------------------------------------------------------*/
int freeAlltoallData(){

	free(alltoallSendBuf);
	free(alltoallRecvBuf);
	free(alltoallFinalBuf);

	return 0;
}

/*-----------------------------------------------------------*/
/* testAlltoall                                              */
/*                                                           */
/* Verifies that the all to all completed successfully.      */
/*-----------------------------------------------------------*/
int testAlltoall(int dataSize){
	int sizeofBuffer, i, j;
	int dataForEachThread, startElem;
	int testFlag, reduceFlag;
	int *testBuf;

	/* Set testFlag to true */
	testFlag = TRUE;

	/* calculate the size of the buffer on each process and allocate */
	sizeofBuffer = dataSize * numThreads * numMPIprocs * numThreads;
	testBuf = (int *) malloc(sizeofBuffer * sizeof(int));

	/* Calculate how many elements each thread will work with */
	dataForEachThread = dataSize * numThreads * numMPIprocs;

	/* Fill buffer with expected values. */
#pragma omp parallel default(none) \
	private(i,j,startElem) \
	shared(testBuf,globalIDarray,sizeofBuffer,dataSize) \
	shared(numThreads,numMPIprocs,dataForEachThread)
	{
	/* Calculate start element for each thread */
	startElem = (myThreadID) * dataForEachThread;

	for (i=0; i<(numThreads * numMPIprocs); i++){
		for (j=0; j<dataSize; j++){
			testBuf[startElem + (i * dataSize) + j] = i;
		}
	}
	}

	/* Compare */
	for (i=0; i<sizeofBuffer; i++){
		if (alltoallFinalBuf[i] != testBuf[i]){
			testFlag = FALSE;
		}
	}

	/* Reduce testFlag with logical AND operator to
	 * get overall test result.
	 */
	MPI_Reduce(&testFlag, &reduceFlag, 1, MPI_INT, MPI_LAND, 0, comm);

	/* Master then sets testOutcome flag */
	if (myMPIRank == 0){
		setTestOutcome(reduceFlag);
	}

	/* Free space for testBuf */
	free(testBuf);

	return 0;
}

