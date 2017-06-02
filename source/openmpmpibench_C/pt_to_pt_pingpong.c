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
/* Contains the point-to-point pingpong mixed mode           */
/* OpenMP/MPI benchmarks.                                    */
/* This includes: -masteronly pingpong                       */
/*                -funnelled pingpong                        */
/*                -multiple pingpong                         */
/*-----------------------------------------------------------*/
#include "pt_to_pt_pingpong.h"


/*-----------------------------------------------------------*/
/* pingPong                                    				 */
/*                                                           */
/* Driver subroutine for the pingpong benchmark.             */
/*-----------------------------------------------------------*/
int pingPong(int benchmarkType){
	int dataSizeIter;
	int sameNode;

	pingRank = PPRanks[0];
	pongRank = PPRanks[1];

	/* Check if pingRank and pongRank are on the same node */
	sameNode = compareProcNames(pingRank,pongRank);

	/* Master process then does some reporting */
	if (myMPIRank == 0){
		/* print message saying if benchmark is inter or intra node */
		printNodeReport(sameNode,pingRank,pongRank);
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

		/* allocate space for the main data arrays */
		allocatePingpongData(sizeofBuffer);

		/* warm-up for either masteronly, funnelled or multiple */
		if (benchmarkType == MASTERONLY){
			/* perform masteronly warm-up sweep */
			masteronlyPingpong(warmUpIters, dataSizeIter);
		}
		else if (benchmarkType == FUNNELLED){
			/* perform funnelled warm-up sweep */
			funnelledPingpong(warmUpIters, dataSizeIter);
		}
		else if (benchmarkType == MULTIPLE){
			multiplePingpong(warmUpIters, dataSizeIter);
		}

		/* perform verification test for the pingpong */
		testPingpong(sizeofBuffer, dataSizeIter);

		/* Initialise benchmark */
		benchComplete = FALSE;

		/* keep executing benchmark until target time is reached */
		while (benchComplete != TRUE){
			/* Start the timer...MPI_Barrier to synchronise */
			MPI_Barrier(comm);
			startTime = MPI_Wtime();

			if (benchmarkType == MASTERONLY){
				/* execute for repsToDo repetitions */
				masteronlyPingpong(repsToDo, dataSizeIter);
			}
			else if (benchmarkType == FUNNELLED){
				funnelledPingpong(repsToDo, dataSizeIter);
			}
			else if (benchmarkType == MULTIPLE){
				multiplePingpong(repsToDo, dataSizeIter);
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
		freePingpongData();

		/* Update dataSize before the next iteration */
		dataSizeIter = dataSizeIter * 2; /* double data size */
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* masteronlyPingpong										 */
/* 															 */
/* One MPI process sends single fixed length message to      */
/* another MPI process.                                      */
/* This other process then sends it back to the first        */
/* process.                                                  */
/*-----------------------------------------------------------*/
int masteronlyPingpong(int totalReps, int dataSize){
	int repIter, i;

	for (repIter = 0; repIter < totalReps; repIter++){
		/* All threads under MPI process with rank = pingRank
		 * write to their part of the pingBuf array using a
		 * parallel for directive.
		 */
		if (myMPIRank == pingRank){
#pragma omp parallel for default(none) \
	private(i) \
	shared(pingSendBuf,dataSize,sizeofBuffer,globalIDarray) \
	schedule(static,dataSize)

			for(i=0; i<sizeofBuffer; i++){
				pingSendBuf[i] = globalIDarray[myThreadID];
			}

			/* Ping process sends buffer to MPI process with rank equal to
			 * pongRank.
			 */
			MPI_Send(pingSendBuf, sizeofBuffer, MPI_INT, pongRank, TAG, comm);

			/* Process then waits for a message from pong process and
			 * each thread reads its part of received buffer.
			 */
			MPI_Recv(pongRecvBuf, sizeofBuffer, MPI_INT, pongRank, \
					TAG, comm, &status);

#pragma omp parallel for default(none) \
	private(i) \
	shared(pongRecvBuf,finalRecvBuf,dataSize,sizeofBuffer) \
	schedule(static,dataSize)
			for(i=0; i<sizeofBuffer; i++){
				finalRecvBuf[i] = pongRecvBuf[i];
			}
		}
		else if (myMPIRank == pongRank){
			/* pongRank receives the message from the ping process */
			MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INT, pingRank, \
					TAG, comm, &status);

			/* each thread under the pongRank MPI process now copies
			 * its part of the received buffer to pongSendBuf.
			 */
#pragma omp parallel for default(none) \
	private(i) \
        shared(pongSendBuf,pingRecvBuf,dataSize,sizeofBuffer) \
	schedule(static,dataSize)
			for(i=0; i< sizeofBuffer; i++){
				pongSendBuf[i] = pingRecvBuf[i];
			}

			/* pongRank process now sends pongSendBuf to ping process. */
			MPI_Send(pongSendBuf, sizeofBuffer, MPI_INTEGER, pingRank, \
					TAG, comm);
		}
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* funnelledPingpong										 */
/*															 */
/* One MPI process sends single fixed length message to      */
/* another MPI process.                                      */
/* This other process then sends it back to the first        */
/* process.                                                  */
/* All communication takes place within the OpenMP           */
/* region for this benchmark.                                */
/*-----------------------------------------------------------*/
int funnelledPingpong(int totalReps, int dataSize){
	int repIter, i;

	/* Open parallel region for threads */
#pragma omp parallel default(none) \
	private(i,repIter) \
	shared(pingRank,pongRank,pingSendBuf,pingRecvBuf) \
	shared(pongSendBuf,pongRecvBuf,finalRecvBuf,sizeofBuffer) \
	shared(dataSize,globalIDarray,comm,status,totalReps,myMPIRank)
	{
	for (repIter=0; repIter< totalReps; repIter++){
		/* All threads under MPI process with rank = pingRank
		 * write to its part of the pingBuf array using a
		 * parallel do directive.
		 */
		if (myMPIRank == pingRank){
#pragma omp for schedule(static,dataSize)
			for (i=0; i<sizeofBuffer; i++){
				pingSendBuf[i] = globalIDarray[myThreadID];
			}
/* Implicit barrier at end of for takes care of synchronisation */

			/*Master thread under ping process sends buffer to
			 * MPI process with rank equal to pongRank.
			 */
#pragma omp master
			{
			MPI_Send(pingSendBuf, sizeofBuffer, MPI_INT, pongRank, TAG, comm);

			/* Master thread then waits for a message from pong process */
			MPI_Recv(pongRecvBuf, sizeofBuffer, MPI_INT, pongRank, TAG, \
					comm, &status);
			}
/* Barrier needed to wait for master thread to complete MPI_Recv */
#pragma omp barrier

		/*Each thread reads its part of received buffer */
#pragma omp for schedule(static,dataSize)
		for (i=0; i<sizeofBuffer; i++){
			finalRecvBuf[i] = pongRecvBuf[i];
		}
	}
	else if (myMPIRank == pongRank){
		/* Master thread under pongRank receives the message
		 * from the ping process.
		 */
#pragma omp master
		{
		MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INT, pingRank, TAG, comm, &status);
		}

/* Barrier needed to wait on master thread */
#pragma omp barrier

		/* Each thread under the pongRank MPI process now copies its part
		 * of the received buffer to pongSendBuf.
		 */
#pragma omp for schedule(static,dataSize)
		for (i=0; i<sizeofBuffer; i++){
			pongSendBuf[i] = pingRecvBuf[i];
		}
/* Implicit barrier at end of DO */

		/* Master thread of pongRank process now sends pongSendBuf
		 * to ping process.
		 */
#pragma omp master
		{
		MPI_Send(pongSendBuf, sizeofBuffer, MPI_INT, pingRank, TAG, comm);
		}
	}
	} /* end of repetitions */
	} /* end of parallel region */

	return 0;
}

/*-----------------------------------------------------------*/
/* multiplePingpong											 */
/*      													 */
/* With this algorithm multiple threads take place in the    */
/* communication and computation.                            */
/* Each thread under the MPI ping process sends a portion    */
/* of the message to the other MPI process.                  */
/* Each thread of the other process then sends it back to    */
/* the first process.                                        */
/*-----------------------------------------------------------*/
int multiplePingpong(int totalReps, int dataSize){
	int repIter, i, lBound, uBound;

	/* Open parallel region for threads under pingRank */
#pragma omp parallel default(none) \
	private(i,repIter,lBound) \
	shared(pingRank,pongRank,pingSendBuf,pingRecvBuf) \
	shared(pongSendBuf,pongRecvBuf,finalRecvBuf,sizeofBuffer) \
	shared(dataSize,globalIDarray,comm,status,totalReps,myMPIRank)
	{
	for (repIter=0; repIter < totalReps; repIter++){

		if (myMPIRank == pingRank){
			/* Calculate lower bound of data array for the thread */
			lBound = (myThreadID * dataSize);

			/* All threads under MPI process with rank = pingRank
			 * write to their part of the pingBuf array using
			 * a parallel for directive.
			 */
#pragma omp for nowait schedule(static,dataSize)
			for (i=0; i<sizeofBuffer; i++){
				pingSendBuf[i] = globalIDarray[myThreadID];
			}
/* Implicit barrier not needed for multiple*/

			/*Each thread under ping process sends dataSize items
			 * to MPI process with rank equal to pongRank.
			 * myThreadID is used as tag to ensure data goes to correct
			 * place in recv buffer.
			 */
			MPI_Send(&pingSendBuf[lBound], dataSize, MPI_INT, pongRank, \
					myThreadID, comm);

			/* Thread then waits for a message from pong process. */
			MPI_Recv(&pongRecvBuf[lBound], dataSize, MPI_INT, pongRank, \
					myThreadID, comm, &status);

			/* Each thread reads its part of the received buffer */
#pragma omp for nowait schedule(static,dataSize)
			for (i=0; i<sizeofBuffer; i++){
				finalRecvBuf[i] = pongRecvBuf[i];
			}
		}
		else if (myMPIRank == pongRank){
			/* Calculate lower bound of the data array */
			lBound = (myThreadID * dataSize);

			/* Each thread under pongRank receives a message
			 * from the ping process.
			 */
			MPI_Recv(&pingRecvBuf[lBound], dataSize, MPI_INT, pingRank, \
					myThreadID, comm, &status);

			/* Each thread now copies its part of the received buffer
			 * to pongSendBuf.
			 */
#pragma omp for nowait schedule(static,dataSize)
			for (i=0; i<sizeofBuffer; i++)
			{
				pongSendBuf[i] = pingRecvBuf[i];
			}

			/* Each thread now sends pongSendBuf to ping process. */
			MPI_Send(&pongSendBuf[lBound], dataSize, MPI_INT, pingRank, \
					myThreadID, comm);
		}

	}/* end of repetitions */
	} /* end of parallel region */
	return 0;
}



/*-----------------------------------------------------------*/
/* allocateData                                              */
/*															 */
/* Allocates space for the main data arrays.                 */
/* Size of each array is specified by subroutine argument.   */
/*-----------------------------------------------------------*/
int allocatePingpongData(int sizeofBuffer){

	pingSendBuf = (int *)malloc(sizeofBuffer * sizeof(int));
	pingRecvBuf = (int *)malloc(sizeofBuffer * sizeof(int));
	pongSendBuf = (int *)malloc(sizeofBuffer * sizeof(int));
	pongRecvBuf = (int *)malloc(sizeofBuffer * sizeof(int));
	finalRecvBuf = (int *)malloc(sizeofBuffer * sizeof(int));

	return 0;
}

/*-----------------------------------------------------------*/
/* freeData                                                  */
/*															 */
/* Deallocates the storage space for the main data arrays.   */
/*-----------------------------------------------------------*/
int freePingpongData(){

	free(pingSendBuf);
	free(pingRecvBuf);
	free(pongSendBuf);
	free(pongRecvBuf);
	free(finalRecvBuf);

	return 0;
}

/*-----------------------------------------------------------*/
/* testPingpong												 */
/*															 */
/* Verifies that the Ping Pong benchmark worked correctly.   */
/*-----------------------------------------------------------*/
int testPingpong(int sizeofBuffer,int dataSize){
	int i, testFlag;
	int *testBuf;

	/* PingRank process checks if pingpong worked ok. */
	if (myMPIRank == pingRank){
		/* initialise testFlag to true (test passed) */
		testFlag = TRUE;

		/* allocate space for the testBuf */
		testBuf = (int *)malloc(sizeofBuffer * sizeof(int));

		/* construct testBuf array with correct values.
		 * These are the values that should be in finalRecvBuf.
		 */
#pragma omp parallel for default(none) \
	private(i) \
	shared(testBuf,dataSize,sizeofBuffer,globalIDarray) \
	schedule(static,dataSize)
		for (i=0; i<sizeofBuffer; i++){
			testBuf[i] = globalIDarray[myThreadID];
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
	/* pingRank broadcasts testFlag to the other processes */
    MPI_Bcast(&testFlag, 1, MPI_INT, pingRank, comm);

    /* Master process sets the testOutcome using testFlag. */
    if (myMPIRank == 0){
    	setTestOutcome(testFlag);
    }

    return 0;
}






