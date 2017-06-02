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
/* Contains the point-to-point pingping mixed mode           */
/* OpenMP/MPI benchmarks.                                    */
/* This includes: -masteronly pingping                       */
/*                -funnelled pingping                        */
/*                -multiple pingping                         */
/*-----------------------------------------------------------*/
#include "pt_to_pt_pingping.h"


/*-----------------------------------------------------------*/
/* pingPing                                    				 */
/*                                                           */
/* Driver subroutine for the pingping benchmark.             */
/*-----------------------------------------------------------*/
int pingPing(int benchmarkType){
	int dataSizeIter;
	int sameNode;

	pingRankA = PPRanks[0];
	pingRankB = PPRanks[1];

	/* Check if pingRankA and pingRankB are on the same node */
	sameNode = compareProcNames(pingRankA, pingRankB);


	if (myMPIRank == 0){
		/* print message saying if benchmark is inter or intra node */
		printNodeReport(sameNode,pingRankA,pingRankB);
		/* then print report column headings. */
		printBenchHeader();
	}

	/* initialise repsToDo to defaultReps at start of benchmark */
	repsToDo = defaultReps;

	/* Loop over data sizes */
	dataSizeIter = minDataSize; /* initialise dataSizeIter to minDataSize */
	while (dataSizeIter <= maxDataSize){
		/* set sizeofBuffer */
		sizeofBuffer = dataSizeIter * numThreads;

		/* Allocate space for main data arrays */
		allocatePingpingData(sizeofBuffer);

		/* warm-up for benchmarkType */
		if (benchmarkType == MASTERONLY){
			/* Masteronly warmp sweep */
			masteronlyPingping(warmUpIters, dataSizeIter);
		}
		else if (benchmarkType == FUNNELLED){
			/* perform funnelled warm-up sweep */
			funnelledPingping(warmUpIters, dataSizeIter);
		}
		else if (benchmarkType == MULTIPLE){
			multiplePingping(warmUpIters, dataSizeIter);
		}

		/* perform verification test for the pingping */
		testPingping(sizeofBuffer, dataSizeIter);

		/* Initialise benchmark */
		benchComplete = FALSE;

		/* keep executing benchmark until target time is reached */
		while (benchComplete != TRUE){
			/* Start the timer...MPI_Barrier to synchronise */
			MPI_Barrier(comm);
			startTime = MPI_Wtime();

			if (benchmarkType == MASTERONLY){
				/* execute for repsToDo repetitions */
				masteronlyPingping(repsToDo, dataSizeIter);
			}
			else if (benchmarkType == FUNNELLED){
				funnelledPingping(repsToDo, dataSizeIter);
			}
			else if (benchmarkType == MULTIPLE){
				multiplePingping(repsToDo, dataSizeIter);
			}

			/* Stop the timer...MPI_Barrier to synchronise processes */
			MPI_Barrier(comm);
			finishTime = MPI_Wtime();
			totalTime = finishTime - startTime;

			/* Call repTimeCheck function to test if target time is reached */
			if (myMPIRank==0){
			  benchComplete = repTimeCheck(totalTime, repsToDo);
			}
			/* Ensure all procs have the same value of benchComplete */
			/* and repsToDo */
			MPI_Bcast(&benchComplete, 1, MPI_INT, 0, comm);
			MPI_Bcast(&repsToDo, 1, MPI_INT, 0, comm);
		}

		/* Master process sets benchmark results */
		if (myMPIRank == 0){
			setReportParams(dataSizeIter, repsToDo, totalTime);
			printReport();
		}

		/* Free the allocated space for the main data arrays */
		freePingpingData();

		/* Update dataSize before the next iteration */
		dataSizeIter = dataSizeIter * 2; /* double data size */
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* masteronlyPingping										 */
/* 															 */
/* Two processes send a message to each other using the      */
/* MPI_Isend, MPI_Recv and MPI_Wait routines.				 */
/* Inter-process communication takes place outside of the    */
/* parallel region.											 */
/*-----------------------------------------------------------*/
int masteronlyPingping(int totalReps, int dataSize){
	int repIter, i;
	int destRank;

	/* set destRank to ID of other process */
	if (myMPIRank == pingRankA){
		destRank = pingRankB;
	}
	else if (myMPIRank == pingRankB){
		destRank = pingRankA;
	}

	for (repIter = 0; repIter < totalReps; repIter++){

		if (myMPIRank == pingRankA || myMPIRank == pingRankB){

			/* Each thread writes its globalID to pingSendBuf
			 * using a PARALLEL DO directive.
			 */
#pragma omp parallel for default(none) \
	private(i) \
	shared(pingSendBuf,dataSize,sizeofBuffer,globalIDarray) \
	schedule(static,dataSize)
			for (i=0; i<sizeofBuffer; i++){
				pingSendBuf[i] = globalIDarray[myThreadID];
			}

			/* Process calls non-bloacking send to start transfer of pingSendBuf
			 * to other process.
			 */
			MPI_Isend(pingSendBuf, sizeofBuffer, MPI_INT, destRank, \
					TAG, comm, &requestID);

			/* Process then waits for message from other process. */
			MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INT, destRank, \
					TAG, comm, &status);

			/* Finish the Send operation with an MPI_Wait */
			MPI_Wait(&requestID, &status);

			/* Each thread under the MPI process now reads its part of the
			 * received buffer.
			 */
#pragma omp parallel for default(none) \
	private(i) \
	shared(finalRecvBuf,dataSize,sizeofBuffer,pingRecvBuf) \
	schedule(static,dataSize)
			for (i=0; i<sizeofBuffer; i++){
				finalRecvBuf[i] = pingRecvBuf[i];
			}

		}
	}
		return 0;
}

/*-----------------------------------------------------------*/
/* funnelledPingPing                            		     */
/* 															 */
/* Two processes send a message to each other using the      */
/* MPI_Isend, MPI_Recv and MPI_Wait routines.                */
/* Inter-process communication takes place inside the        */
/* OpenMP parallel region.                                   */
/*-----------------------------------------------------------*/
int funnelledPingping(int totalReps, int dataSize){
	int repIter, i;
	int destRank;

    /* set destRank to ID of other process */
    if (myMPIRank == pingRankA){
    	destRank = pingRankB;
    }
    else if (myMPIRank == pingRankB){
    	destRank = pingRankA;
    }

	/* Open the parallel region */
#pragma omp parallel default(none) \
	private(i, repIter) \
	shared(dataSize,sizeofBuffer,pingSendBuf,globalIDarray) \
	shared(pingRecvBuf,finalRecvBuf,status,requestID) \
	shared(destRank,comm,myMPIRank,pingRankA,pingRankB,totalReps)

	for (repIter = 0; repIter < totalReps; repIter++){


		if (myMPIRank == pingRankA || myMPIRank == pingRankB){

			/* Each thread writes its globalID to its part of
			 * pingSendBuf.
			 */
#pragma omp for schedule(static,dataSize)
			for (i=0; i<sizeofBuffer; i++){
				pingSendBuf[i] = globalIDarray[myThreadID];
			}
/* Implicit barrier here takes care of necessary synchronisation */

#pragma omp master
			{
			/* Master thread starts send of buffer */
			MPI_Isend(pingSendBuf, sizeofBuffer, MPI_INT, destRank, \
					TAG, comm, &requestID);

			/* then waits for message from other process */
			MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INT, destRank, \
					TAG, comm, &status);

			/* Master thread then completes send using an MPI_Wait */
			MPI_Wait(&requestID, &status);
			}

/* Barrier needed to ensure master thread has completed transfer */
#pragma omp barrier

			/* Each thread reads its part of the received buffer */
#pragma omp for schedule(static,dataSize)
			for (i=0; i<sizeofBuffer; i++){
				finalRecvBuf[i] = pingRecvBuf[i];
			}
		}
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* multiplePingping                            				 */
/* 															 */
/* With this algorithm multiple threads take place in the    */
/* communication and computation.                            */
/* Each thread sends its portion of the pingSendBuf to the   */
/* other process using MPI_Isend/ MPI_Recv/ MPI_Wait         */
/* routines.                                                 */
/*-----------------------------------------------------------*/
int multiplePingping(int totalReps, int dataSize){
	int repIter, i;
	int destRank;
	int lBound;

    /* set destRank to ID of other process */
    if (myMPIRank == pingRankA){
    	destRank = pingRankB;
    }
    else if (myMPIRank == pingRankB){
    	destRank = pingRankA;
    }

    /* Open parallel region */
#pragma omp parallel default(none) \
	private(i,lBound,requestID,status,repIter) \
	shared(pingSendBuf,pingRecvBuf,finalRecvBuf,sizeofBuffer) \
	shared(destRank,myMPIRank,pingRankA,pingRankB,totalReps) \
	shared(dataSize,globalIDarray,comm)
    {
    for (repIter = 0; repIter < totalReps; repIter++){

    	if (myMPIRank == pingRankA || myMPIRank == pingRankB){

    		/* Calculate the lower bound of each threads
    		 * portion of the data arrays.
    		 */
    		lBound = (myThreadID * dataSize);

    		/* Each thread writes to its part of pingSendBuf */
#pragma omp for nowait schedule(static,dataSize)
    		for (i=0; i<sizeofBuffer; i++){
    			pingSendBuf[i] = globalIDarray[myThreadID];
    		}

    		/* Each thread starts send of dataSize items of
    		 * pingSendBuf to process with rank = destRank.
    		 */
    		MPI_Isend(&pingSendBuf[lBound], dataSize, MPI_INT, destRank, \
    				myThreadID, comm, &requestID);

    		/* Thread then waits for message from destRank with
    		 * tag equal to it thread id.
    		 */
    		MPI_Recv(&pingRecvBuf[lBound], dataSize, MPI_INT, destRank, \
    				myThreadID, comm, &status);

    		/* Thread completes send using MPI_Wait */
    		MPI_Wait(&requestID, &status);

    		/* Each thread reads its part of received buffer. */
#pragma omp for nowait schedule(static,dataSize)
    		for (i=0; i<sizeofBuffer; i++){
    			finalRecvBuf[i] = pingRecvBuf[i];
    		}
    	}
    }
    }

    return 0;
}

/*-----------------------------------------------------------*/
/* allocatePingpingData                                      */
/*															 */
/* Allocates space for the main data arrays.                 */
/* Size of each array is specified by subroutine argument.   */
/*-----------------------------------------------------------*/
int allocatePingpingData(int sizeofBuffer){

	pingSendBuf = (int *)malloc(sizeofBuffer * sizeof(int));
	pingRecvBuf = (int *)malloc(sizeofBuffer * sizeof(int));
	finalRecvBuf = (int *)malloc(sizeofBuffer * sizeof(int));

	return 0;
}

/*-----------------------------------------------------------*/
/* freePingpingData                                          */
/*															 */
/* Deallocates the storage space for the main data arrays.   */
/*-----------------------------------------------------------*/
int freePingpingData(){

	free(pingSendBuf);
	free(pingRecvBuf);
	free(finalRecvBuf);

	return 0;
}

/*-----------------------------------------------------------*/
/* testPingping												 */
/*															 */
/* Verifies that the PingPing benchmark worked correctly.    */
/*-----------------------------------------------------------*/
int testPingping(int sizeofBuffer,int dataSize){
	int otherPingRank, i, testFlag, reduceFlag;
	int *testBuf;

	/* initialise testFlag to true (test passed) */
	testFlag = TRUE;

	/* Testing only needs to be done by pingRankA & pingRankB */
	if (myMPIRank == pingRankA || myMPIRank == pingRankB){
		/* allocate space for testBuf */
		testBuf = (int *)malloc(sizeofBuffer * sizeof(int));

		/* set the ID of other pingRank */
		if (myMPIRank == pingRankA){
			otherPingRank = pingRankB;
		}
		else if (myMPIRank == pingRankB){
			otherPingRank = pingRankA;
		}

		/* construct testBuf array with correct values.
		 * These are the values that should be in finalRecvBuf.
		 */
#pragma omp parallel for default(none) \
	private(i) \
	shared(otherPingRank,numThreads,testBuf,dataSize,sizeofBuffer) \
	schedule(static,dataSize)
		for (i=0; i<sizeofBuffer; i++){
			/* calculate globalID of thread expected in finalRecvBuf
			 * This is done by using otherPingRank
			 */
			testBuf[i] = (otherPingRank * numThreads) + myThreadID;
		}

		/* compare each element of testBuf and finalRecvBuf */
		for (i=0; i<sizeofBuffer; i++){
			if (testBuf[i] != finalRecvBuf[i]){
				testFlag = FALSE;
			}
		}

		/* free space for testBuf */
		free(testBuf);
	}


	MPI_Reduce(&testFlag, &reduceFlag, 1, MPI_INT, MPI_LAND, 0, comm);

	/* Master process sets the testOutcome using testFlag. */
	 if (myMPIRank == 0){
		 setTestOutcome(reduceFlag);
	 }

	 return 0;
}
