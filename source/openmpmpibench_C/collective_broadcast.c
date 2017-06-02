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
/* Implements the mixed mode OpenMP/MPI collective           */
/* broadcast benchmark.                                      */
/*-----------------------------------------------------------*/
#include "collective_broadcast.h"

/*-----------------------------------------------------------*/
/* broadcast                                   	             */
/* 															 */
/* Driver subroutine for the broadcast benchmark.            */
/*-----------------------------------------------------------*/
int broadcast(){
	int dataSizeIter, sizeofFinalBuf;

	/* initialise repsToDo to defaultReps */
	repsToDo = defaultReps;

	/* Start loop over data sizes */
	dataSizeIter = minDataSize;

	while (dataSizeIter <= maxDataSize){
		/* allocate space for main data arrays */
		allocateBroadcastData(dataSizeIter);

		/* perform benchmark warm-up */
		broadcastKernel(warmUpIters,dataSizeIter);

		/* set sizeofFinalBuf and test if broadcast was a success */
		sizeofFinalBuf = dataSizeIter * numThreads;
		testBroadcast(sizeofFinalBuf);

		/* Initialise the benchmark */
		benchComplete = FALSE;

		/* Execute benchmark until target time is reached */
		while (benchComplete != TRUE){
			/* Start timer */
			MPI_Barrier(comm);
			startTime = MPI_Wtime();

			/* Execute broadcast for repsToDo repetitions */
			broadcastKernel(repsToDo, dataSizeIter);

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

		/* master process sets benchmark result for reporting */
		if (myMPIRank == 0){
			setReportParams(dataSizeIter, repsToDo, totalTime);
			printReport();

		}

		/* Free allocated data */
		freeBroadcastData();

		/* Double dataSize and loop again */
		dataSizeIter = dataSizeIter * 2;

	}

	return 0;
}

/*-----------------------------------------------------------*/
/* broadcastKernel                                           */
/*                                                           */
/* The broadcast benchmark.                                  */
/* At the start one process owns the data. After, all        */
/* processes and threads have a copy of the data.            */
/*-----------------------------------------------------------*/
int broadcastKernel(int totalReps, int dataSize){
	int repIter, i;
	int startPos; /* Start position in finalBroadcastBuf for each thread */

	for (repIter=0; repIter<totalReps; repIter++){

		/* Master MPI process writes to broadcastBuf */
		if (myMPIRank == BROADCASTROOT){
			for (i=0; i<dataSize; i++){
				broadcastBuf[i] = BROADCASTNUM;
			}
		}
		/* Broadcast array to all other processes */
		MPI_Bcast(broadcastBuf, dataSize, MPI_INT, BROADCASTROOT, comm);

		/* Each thread copies broadcastBuf to its portion of finalBroadcastBuf */
#pragma omp parallel default(none) \
	private(i,startPos) \
	shared(dataSize,finalBroadcastBuf,broadcastBuf)
		{
		startPos = ((myThreadID) * dataSize);
		for (i=0; i<dataSize; i++){
			finalBroadcastBuf[startPos + i] = broadcastBuf[i];
		}
		}
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* allocateBroadcastData                                     */
/*                                                           */
/* Allocate memory for the main data arrays in the           */
/* broadcast operation.                                      */
/*-----------------------------------------------------------*/
int allocateBroadcastData(int bufferSize){

	broadcastBuf = (int *)malloc(bufferSize * sizeof(int));

	/* finalBroadcastBuf is of size dataSize*numThreads */
	finalBroadcastBuf = (int *)malloc((bufferSize*numThreads)*sizeof(int));

	return 0;
}

/*-----------------------------------------------------------*/
/* freeBroadcastData                                         */
/*                                                           */
/* Free memory of main data arrays.                          */
/*-----------------------------------------------------------*/
int freeBroadcastData(){

	free(broadcastBuf);
	free(finalBroadcastBuf);

	return 0;
}

/*-----------------------------------------------------------*/
/* testBroadcast                                             */
/*                                                           */
/* Verifies that the broadcast benchmark worked correctly.   */
/*-----------------------------------------------------------*/
int testBroadcast(int bufferSize){
	int i, testFlag, reduceFlag;

	/* Initialise testFlag to true */
	testFlag = TRUE;

	/* Compare each element of finalBroadcast with BROADCASTNUM */
	for (i=0; i<bufferSize; i++){
		if (finalBroadcastBuf[i] != BROADCASTNUM){
			testFlag = FALSE;
		}
	}

	/* Reduce testFlag to master with logical AND operation */
	MPI_Reduce(&testFlag, &reduceFlag, 1, MPI_INT, MPI_LAND, 0, comm);

	/* Master then sets testOutcome using reduceFlag */
	if (myMPIRank == 0){
		setTestOutcome(testFlag);
	}

	return 0;
}
