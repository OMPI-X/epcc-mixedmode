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
/* Contains the point-to-point multi-pingpong mixed mode     */
/* OpenMP/MPI benchmarks.                                    */
/* This includes: -masteronly multiPingpong                  */
/*                -funnelled multiPingpong                   */
/*                -multiple multiPingpong                    */
/*-----------------------------------------------------------*/
#include "pt_to_pt_multiPingpong.h"

/*-----------------------------------------------------------*/
/* multiPingPong                                             */
/*                                                           */
/* Driver subroutine for the multi-pingpong benchmark.       */
/*-----------------------------------------------------------*/
int multiPingPong(int benchmarkType){
	int dataSizeIter;
	int pongWorldRank;
	char pongProcName[MPI_MAX_PROCESSOR_NAME];
	int balance;

	pingNode = 0;
	pongNode = 1;

	/* Check if there's a balance in num of MPI processes
	  on pingNode and pongNode. */
	balance = crossCommBalance(pingNode, pongNode);
	/* If not balanced.. */
	if (balance == FALSE){
		/* ..master prints error */
		if (myMPIRank == 0){
			printBalanceError();
		}
		/* ..and all process exit function. */
		return 1;
	}

	/* Exchange MPI_COMM_WORLD ranks for processes in same crossComm */
	exchangeWorldRanks(pingNode, pongNode, &pongWorldRank);

	/* Processes on pongNode send processor name to pingNode procs. */
	sendProcName(pingNode, pongNode, pongProcName);

	/* Print comm world ranks & processor name of processes
	 * taking part in multi-pingpong benchmark.
	 */
	printMultiProcInfo(pingNode, pongWorldRank, pongProcName);

	/* Barrier to ensure that all procs have completed
	 * printMultiProcInfo before prinring column headings.
	 */
	MPI_Barrier(comm);

	/* Master process then prints report column headings */
	if (myMPIRank == 0){
		printBenchHeader();
	}

	/* Initialise repsToDo to defaultReps at start of benchmark */
	repsToDo = defaultReps;
	dataSizeIter = minDataSize; /* initialise dataSizeIter to minDataSize */

	/* Loop over data sizes */
	while (dataSizeIter <= maxDataSize){
		/* set sizeofBuffer */
		sizeofBuffer = dataSizeIter * numThreads;

		/* Allocate space for the main data arrays */
		allocateMultiPingpongData(sizeofBuffer);

		/* warm-up */
		if (benchmarkType == MASTERONLY){
			/* Masteronly warm-up */
			masteronlyMultiPingpong(warmUpIters, dataSizeIter);
		}
		else if (benchmarkType == FUNNELLED){
			/* Funnelled warm-up sweep */
			funnelledMultiPingpong(warmUpIters, dataSizeIter);
		}
		else if (benchmarkType == MULTIPLE){
			/* Multiple pingpong warm-up */
			multipleMultiPingpong(warmUpIters, dataSizeIter);
		}

		/* Verification test for multi-pingpong */
		testMultiPingpong(sizeofBuffer, dataSizeIter);

		/* Initialise benchmark */
		benchComplete = FALSE;

		/* Keep executing benchmark until target time is reached */
		while (benchComplete != TRUE){

			/* MPI_Barrier to synchronise processes.
			   Then start the timer. */
			MPI_Barrier(comm);
			startTime = MPI_Wtime();

			if (benchmarkType == MASTERONLY){
				/* Execute masteronly multipingpong repsToDo times */
				masteronlyMultiPingpong(repsToDo, dataSizeIter);
			}
			else if (benchmarkType == FUNNELLED){
				/* Execute funnelled multipingpong */
				funnelledMultiPingpong(repsToDo, dataSizeIter);
			}
			else if (benchmarkType == MULTIPLE){
				multipleMultiPingpong(repsToDo, dataSizeIter);
			}

			/* Stop the timer..MPI_Barrier to synchronise processes
			 * for more accurate timing.
			 */
			MPI_Barrier(comm);
			finishTime = MPI_Wtime();
			totalTime = finishTime - startTime;

			/* Call repTimeCheck to check if target time is reached. */
			if (myMPIRank==0){
			  benchComplete = repTimeCheck(totalTime, repsToDo);
			}
			/* Ensure all procs have the same value of benchComplete */
			/* and repsToDo */
			MPI_Bcast(&benchComplete, 1, MPI_INT, 0, comm);
			MPI_Bcast(&repsToDo, 1, MPI_INT, 0, comm);
		} /* End of loop to check if benchComplete is true */

		/* Master process sets benchmark results */
		if (myMPIRank == 0){
			setReportParams(dataSizeIter, repsToDo, totalTime);
			printReport();
		}

		/* Free the allocated space for the main data arrays */
		freeMultiPingpongData();

		/* Update dataSize before next iteration */
		dataSizeIter = dataSizeIter * 2;

	} /* end loop over data sizes */

	return 0;
}

/*-----------------------------------------------------------*/
/* masteronlyMultiPingpong                                   */
/*                                                           */
/* All MPI processes in crossComm = pingNode sends a single  */
/* fixed length message to the neighbouring process in       */
/* crossComm = pongNode.                                     */
/* The neighbouring processes then sends the message back    */
/* to the first process.                                     */
/*-----------------------------------------------------------*/
int masteronlyMultiPingpong(int totalReps, int dataSize){
	int repIter, i;

	for (repIter = 1; repIter <= totalReps; repIter++){

		/* Threads under each MPI process with
		 * crossCommRank = pingNode write to pingSendBuf
		 * array with a PARALLEL FOR directive.
		 */
		if (crossCommRank == pingNode){

#pragma omp parallel for default(none) \
	private(i) \
	shared(pingSendBuf,dataSize,sizeofBuffer,globalIDarray) \
	schedule(static,dataSize)

			for (i=0; i<sizeofBuffer; i++){
				pingSendBuf[i] = globalIDarray[myThreadID];
			}

			/* Each process with crossCommRank = pingNode sends
			 * buffer to MPI process with rank = pongNode in crossComm.
			 */
			MPI_Send(pingSendBuf, sizeofBuffer, MPI_INT, pongNode, TAG, crossComm);

			/* The processes then wait for a message from pong process
			 * and each thread reads its part of the received buffer.
			 */
			MPI_Recv(pongRecvBuf, sizeofBuffer, MPI_INT, pongNode, \
					TAG, crossComm, &status);

#pragma omp parallel for default(none) \
	private(i) \
	shared(pongRecvBuf,finalRecvBuf,dataSize,sizeofBuffer) \
	schedule(static,dataSize)

			for (i=0; i<sizeofBuffer; i++){
				finalRecvBuf[i] = pongRecvBuf[i];
			}
		}
		else if (crossCommRank == pongNode){

			/* Each process with crossCommRank = pongNode receives
			 * the message from the pingNode processes.
			 */
			MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INT, pingNode,\
					TAG, crossComm, &status);

			/* Each thread copies its part of the received buffer
			 * to pongSendBuf.
			 */
#pragma omp parallel for default(none) \
	private(i) \
	shared(pongSendBuf,pingRecvBuf,dataSize,sizeofBuffer) \
	schedule(static,dataSize)

			for (i=0; i<sizeofBuffer; i++){
				pongSendBuf[i] = pingRecvBuf[i];
			}

			/* The processes now send pongSendBuf to processes
			 * with crossCommRank = pingNode.
			 */
			MPI_Send(pongSendBuf, sizeofBuffer, MPI_INT, pingNode, \
					TAG, crossComm);
		}
	} /* End repetitions loop */

	return 0;
}

/*-----------------------------------------------------------*/
/* funnelledMultiPingpong                                    */
/*                                                           */
/* All MPI processes in crossComm = pingNode sends a single  */
/* fixed length message to the neighbouring process in       */
/* crossComm = pongNode.                                     */
/* The neighbouring processes then sends the message back    */
/* to the first process.                                     */
/* All communication takes place within the OpenMP parallel  */
/* region for this benchmark.                                */
/*-----------------------------------------------------------*/
int funnelledMultiPingpong(int totalReps, int dataSize){
	int repIter, i;

	/* Open the parallel region for threads */
#pragma omp parallel default(none) \
	private(i,repIter) \
	shared(pingNode,pongNode,pingSendBuf,pingRecvBuf) \
	shared(pongSendBuf,pongRecvBuf,finalRecvBuf,sizeofBuffer) \
	shared(dataSize,globalIDarray,crossComm,status) \
	shared(totalReps,myMPIRank,crossCommRank)
	{

		/* loop totalRep times */
		for (repIter = 1; repIter <= totalReps; repIter++){

			/* All threads under each MPI process with
			 * crossCommRank = pingNode write to pingSendBuf
			 * array using a parallel for directive.
			 */
			if (crossCommRank == pingNode){
#pragma omp for schedule(static,dataSize)

				for (i=0; i<sizeofBuffer; i++){
					pingSendBuf[i] = globalIDarray[myThreadID];
				}
/* Implicit barrier at end of omp for takes care of synchronisation */

				/* Master thread under each pingNode process sends
				 * buffer to corresponding MPI process in pongNode
				 * using crossComm.
				 */
#pragma omp master
			{
				MPI_Send(pingSendBuf, sizeofBuffer, MPI_INT, pongNode, TAG, crossComm);

				/* Master thread then waits for a message from the pong process. */
				MPI_Recv(pongRecvBuf, sizeofBuffer, MPI_INT, pongNode, TAG, \
					crossComm, &status);
			}
/* Barrier needed to wait for master thread to complete MPI_Recv */
#pragma omp barrier

			/* Each thread then reads its part of the received buffer. */
#pragma omp for schedule(static,dataSize)

				for (i=0; i<sizeofBuffer; i++){
					finalRecvBuf[i] = pongRecvBuf[i];
				}

			}
			else if (crossCommRank == pongNode){

				/* Master thread under each pongNode process receives
				 * the message from the pingNode processes.
				 */
#pragma omp master
			{
					MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INT, pingNode,\
							TAG, crossComm, &status);
			}
/* Barrier needed to wait on master thread */
#pragma omp barrier

			/* Each thread reads its part of the received buffer. */
#pragma omp for schedule(static,dataSize)

				for (i=0; i<sizeofBuffer; i++){
					pongSendBuf[i] = pingRecvBuf[i];
				}
/* Implicit barrier at end of omp for */

			/* Master threads send their pongSendBuf to processes
			 * with crossCommRank = pingNode.
			 */
#pragma omp master
			{

				MPI_Send(pongSendBuf, sizeofBuffer, MPI_INT, pingNode, TAG, crossComm);

			}

			}
		} /* End of repetitions loop. */
	} /* End of parallel region */

	return 0;
}

/*-----------------------------------------------------------*/
/* multipleMultiPingpong                                     */
/*                                                           */
/* Multiple threads take place in the communication and      */
/* computation.                                              */
/* Each thread of all MPI processes in crossComm = pingNode  */
/* sends a portion of the message to the neighbouring        */
/* process in crossComm = pongNode.                          */
/* Each thread of the neighbouring processes then sends      */
/* the message back to the first process.                    */
/*-----------------------------------------------------------*/
int multipleMultiPingpong(int totalReps, int dataSize){
	int repIter, i;
	int lBound;

	/* Open parallel region for threads */
#pragma omp parallel default(none) \
	private(i,repIter,status,lBound) \
	shared(pingNode,pongNode,pingSendBuf,pingRecvBuf) \
	shared(pongSendBuf,pongRecvBuf,finalRecvBuf,sizeofBuffer) \
	shared(dataSize,globalIDarray,crossComm) \
	shared(totalReps,myMPIRank,crossCommRank)
	{
		for (repIter=1; repIter<=totalReps; repIter++){ /* loop totalRep times */

			if (crossCommRank == pingNode){

				/* Calculate lower bound of data array for the thread */
				lBound = (myThreadID * dataSize);

				/* All threads write to its part of the pingBuf
				 * array using a parallel for directive.
				 */
#pragma omp for nowait schedule(static,dataSize)

				for (i=0; i<sizeofBuffer; i++){
					pingSendBuf[i] = globalIDarray[myThreadID];
				}
/* Implicit barrier at end of for not needed for multiple */

				/* Each thread under ping process sends dataSize items
				 * to pongNode process in crossComm.
				 * myThreadID is used as tag to ensure data goes to
				 * correct place in buffer.
				 */
				MPI_Send(&pingSendBuf[lBound], dataSize, MPI_INT, pongNode, \
						myThreadID, crossComm);

				/* Thread then waits for a message from pongNode. */
				MPI_Recv(&pongRecvBuf[lBound], dataSize, MPI_INT, pongNode, \
						myThreadID, crossComm, &status);

				/* Each thread reads its part of the received buffer. */
#pragma omp for nowait schedule(static,dataSize)

				for (i=0; i<sizeofBuffer; i++){
					finalRecvBuf[i] = pongRecvBuf[i];
				}

			}
			else if (crossCommRank == pongNode){
				/* Calculate lower and upper bound of data array */
				lBound = (myThreadID * dataSize);

				/* Each thread under pongRank receives a message from
				 * the ping process.
				 */
				MPI_Recv(&pingRecvBuf[lBound], dataSize, MPI_INT, pingNode, \
						myThreadID, crossComm, &status);

				/* Each thread now copies its part of the received buffer
				 * to pongSendBuf.
				 */
#pragma omp for nowait schedule(static,dataSize)
				for (i=0; i<sizeofBuffer; i++){
					pongSendBuf[i] = pingRecvBuf[i];
				}

				/* Each thread now sends pongSendBuf to ping process. */
				MPI_Send(&pongSendBuf[lBound], dataSize, MPI_INT, pingNode, \
						myThreadID, crossComm);

			}
		} /* End repetitions loop */
	} /* End parallel region */

	return 0;
}

/*-----------------------------------------------------------*/
/* allocateMultiPingpongData                                 */
/*                                                           */
/* Allocates space for the main data arrays.                 */
/* Size of each array is specified by subroutine argument.   */
/*-----------------------------------------------------------*/
int allocateMultiPingpongData(int sizeofBuffer){

	if (crossCommRank == pingNode){
		/* allocate space for arrays that MPI processes
		 * with crossCommRank = pingRank will use.
		 */
		pingSendBuf = (int *)malloc(sizeof(int) * sizeofBuffer);
		pongRecvBuf = (int *)malloc(sizeof(int) * sizeofBuffer);
		finalRecvBuf = (int *)malloc(sizeof(int) * sizeofBuffer);
	}
	else if (crossCommRank == pongNode){
		/* allocate space for arrays that MPI processes
		 * with crossCommRank = pongNode will use.
		 */
		pingRecvBuf = (int *)malloc(sizeof(int) * sizeofBuffer);
		pongSendBuf = (int *)malloc(sizeof(int) * sizeofBuffer);
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* freeMultiPingpongData                                     */
/*                                                           */
/* Deallocates the storage space for the main data arrays.   */
/*-----------------------------------------------------------*/
int freeMultiPingpongData(){

	if (crossCommRank == pingNode){
		free(pingSendBuf);
		free(pongRecvBuf);
		free(finalRecvBuf);
	}
	else if (crossCommRank == pongNode){
		free(pingRecvBuf);
		free(pongSendBuf);
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* testMultiPingpong                                         */
/*                                                           */
/* Verifies the the multi pingpong benchmark worked          */
/* correctly.                                                */
/*-----------------------------------------------------------*/
int testMultiPingpong(int sizeofBuffer, int dataSize){
	int i;
	int testFlag, localTestFlag;

	/* Initialise localTestFlag to true */
	localTestFlag = TRUE;

	/* All processes with crossCommRank = pingNode check
	 * if multi-pingpong worked ok.
	 */
	if (crossCommRank == pingNode){

		/* allocate space for testBuf */
		testBuf = (int *)malloc(sizeof(int) * sizeofBuffer);

		/* Construct testBuf array with correct values.
		 * These are the values that should be in finalRecvBuf.
		 */
#pragma omp parallel for default(none) \
	private(i) \
	shared(testBuf,dataSize,sizeofBuffer,globalIDarray)\
	schedule(static,dataSize)

		for (i=0; i<sizeofBuffer; i++){
			testBuf[i] = globalIDarray[myThreadID];
		}


		/* Compare each element of testBuf and finalRecvBuf */
		for (i=0; i<sizeofBuffer; i++){
			if (testBuf[i] != finalRecvBuf[i]){
				localTestFlag = FALSE;
			}
		}

		/* Free space for testBuf */
		free(testBuf);
	}

	/* Reduce localTestFlag to master */
	MPI_Reduce(&localTestFlag, &testFlag, 1, MPI_INT,MPI_LAND, 0, comm);

	/* Master then sets testOutcome using reduceFlag */
	if (myMPIRank == 0){
		setTestOutcome(testFlag);
	}

	return 0;
}

