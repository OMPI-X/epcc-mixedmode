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
/* Contains the point-to-point multi-pingping mixed mode     */
/* OpenMP/MPI benchmarks.                                    */
/* This includes: -masteronly multiPingping                  */
/*                -funnelled multiPingping                   */
/*                -multiple multiPingping                    */
/*-----------------------------------------------------------*/
#include "pt_to_pt_multiPingping.h"

/*-----------------------------------------------------------*/
/* multiPingPing                                             */
/*                                                           */
/* Driver subroutine for the multi-pingping benchmark.       */
/*-----------------------------------------------------------*/
int multiPingping(int benchmarkType){
  int dataSizeIter;
  char otherProcName[MPI_MAX_PROCESSOR_NAME];
  int balance;

  pingNodeA = 0;
  pingNodeB = 1;

  /* Check if there's a balance in num of MPI processes
     on pingNodeA and pingNodeB. */
  balance = crossCommBalance(pingNodeA, pingNodeB);
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
  exchangeWorldRanks(pingNodeA, pingNodeB, &otherPingRank);

  /* Processes on pongNode send processor name to pingNode procs. */
  sendProcName(pingNodeA, pingNodeB, otherProcName);

  /* Print comm world ranks & processor name of processes
   * taking part in multi-pingpong benchmark.
   */
  printMultiProcInfo(pingNodeA, otherPingRank, otherProcName);

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
  /* Initialise dataSizeIter */
  dataSizeIter = minDataSize;

  /* Start loop over data sizes */
  while (dataSizeIter <= maxDataSize){
    /* set size of buffer */
    sizeofBuffer = dataSizeIter * numThreads;

    /* Allocate space for the main data arrays */
    allocateMultiPingpingData(sizeofBuffer);

    /* warm-up */
    if (benchmarkType == MASTERONLY){
      /* Masteronly warm-up */
      masteronlyMultiPingping(warmUpIters, dataSizeIter);
    }
    else if (benchmarkType == FUNNELLED){
      /* Funnelled warm-up sweep */
      funnelledMultiPingping(warmUpIters, dataSizeIter);
    }
    else if (benchmarkType == MULTIPLE){
      /* Multiple pingpong warm-up */
      multipleMultiPingping(warmUpIters, dataSizeIter);
    }

    /* Verification test for multi-pingpong */
    testMultiPingping(sizeofBuffer, dataSizeIter);

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
	masteronlyMultiPingping(repsToDo, dataSizeIter);
      }
      else if (benchmarkType == FUNNELLED){
	/* Execute funnelled multipingpong */
	funnelledMultiPingping(repsToDo, dataSizeIter);
      }
      else if (benchmarkType == MULTIPLE){
	multipleMultiPingping(repsToDo, dataSizeIter);
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
    freeMultiPingpingData();

    /* Update dataSize before next iteration */
    dataSizeIter = dataSizeIter * 2;
  }

  return 0;
}

/*-----------------------------------------------------------*/
/* masteronlyMultiPingping                                   */
/*                                                           */
/* All Processes with rank of pingNodeA or pingNodeB in      */
/* crossComm send a message to each other.                   */
/* MPI communication takes place outside of the parallel     */
/* region.                                                   */
/*-----------------------------------------------------------*/
int masteronlyMultiPingping(int totalReps, int dataSize){
  int repIter, i;
  int destRank;

  /* set destRank to ID of other process */
  if (crossCommRank == pingNodeA){
    destRank = pingNodeB;
  }
  else if (crossCommRank == pingNodeB){
    destRank = pingNodeA;
  }

  /* loop totalRep times */
  for (repIter=1; repIter<=totalReps; repIter++){

    if ((crossCommRank == pingNodeA) || (crossCommRank == pingNodeB) ){

      /* Each thread writes its globalID to pingSendBuf
       * with a parallel for directive.
       */
#pragma omp parallel for default(none)				\
  private(i)							\
  shared(pingSendBuf,dataSize,sizeofBuffer,globalIDarray)	\
  schedule(static,dataSize)

      for (i=0; i<sizeofBuffer; i++){
	pingSendBuf[i] = globalIDarray[myThreadID];
      }

      /* Process calls non-blocking send to start transfer of
       * pingSendBuf to other process.
       */
      MPI_Isend(pingSendBuf, sizeofBuffer, MPI_INT, destRank, TAG,\
		crossComm, &requestID);

      /* Processes then wait for message from other process. */
      MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INT, destRank, TAG, \
	       crossComm, &status);

      /* Finish the send operation with an MPI_Wait */
      MPI_Wait(&requestID, &status);

      /* Threads under the MPI processes read their part of the
       * received buffer.
       */
#pragma omp parallel for default(none)				\
  private(i)							\
  shared(finalRecvBuf,dataSize,sizeofBuffer,pingRecvBuf)	\
  schedule(static,dataSize)

      for (i=0; i<sizeofBuffer; i++){
	finalRecvBuf[i] = pingRecvBuf[i];
      }

    }
  } /* End repetitions loop */

  return 0;
}

/*-----------------------------------------------------------*/
/* funnelledMultiPingping                                    */
/*                                                           */
/* All processes with rank of pingNodeA or pingNodeB in      */
/* crossComm send a message to each other.                   */
/* Inter-process communication takes place inside the        */
/* OpenMP parallel region by the master thread.              */
/*-----------------------------------------------------------*/
int funnelledMultiPingping(int totalReps, int dataSize){
  int repIter, i;
  int destRank;

  /* Set destRank to id of other process */
  if (crossCommRank == pingNodeA){
    destRank = pingNodeB;
  }
  else if (crossCommRank == pingNodeB){
    destRank = pingNodeA;
  }

  /* Open the parallel region */
#pragma omp parallel default(none)				\
  private(i,repIter)						\
  shared(dataSize,sizeofBuffer,pingSendBuf,globalIDarray)	\
  shared(pingRecvBuf,finalRecvBuf,status,requestID,destRank)	\
  shared(crossComm,crossCommRank,pingNodeA,pingNodeB,totalReps)
  {
    /* loop totalRep times */
    for (repIter = 1; repIter <= totalReps; repIter++){

      if (crossCommRank == pingNodeA || crossCommRank == pingNodeB){
	/* Each thread writes its globalID to its part of
	 * pingSendBuf with an omp for.
	 */
#pragma omp for schedule(static,dataSize)

	for (i=0; i<sizeofBuffer; i++){
	  pingSendBuf[i] = globalIDarray[myThreadID];
	}
	/* Implicit barrier here takes care of necessary synchronisation. */

#pragma omp master
	{
	  /* Master thread of each process starts send. */
	  MPI_Isend(pingSendBuf, sizeofBuffer, MPI_INT, \
		    destRank, TAG, crossComm, &requestID);

	  /* Processes then wait for message. */
	  MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INT, \
		   destRank, TAG, crossComm, &status);

	  /* Finish the send operation with an MPI_Wait */
	  MPI_Wait(&requestID, &status);
	}
	/* Barrier to ensure master thread has completed transfer. */
#pragma omp barrier

	/* Each thread reads its part of the received buffer */
#pragma omp for schedule(static,dataSize)
	for (i=0; i<sizeofBuffer; i++){
	  finalRecvBuf[i] = pingRecvBuf[i];
	}

      }
    } /* End repetitions loop */

  } /* End parallel region */

  return 0;

}

/*-----------------------------------------------------------*/
/* multipleMultiPingping                                     */
/* 															 */
/* All processes with crossCommRank of pingNodeA and         */
/* pingNodeB in crossComm send a message to each other.      */
/* Multiple threads take part in the communication.          */
/*-----------------------------------------------------------*/
int multipleMultiPingping(int totalReps, int dataSize){
  int repIter, i;
  int destRank;
  int lBound;

  /* set destRank to be id of other process */
  if (crossCommRank == pingNodeA){
    destRank = pingNodeB;
  }
  else if (crossCommRank == pingNodeB){
    destRank = pingNodeA;
  }

  /* Open parallel region */
#pragma omp parallel default(none)				\
  private(i,repIter,lBound,requestID,status)			\
  shared(dataSize,sizeofBuffer,pingSendBuf,globalIDarray)	\
  shared(pingRecvBuf,finalRecvBuf,destRank,crossComm)		\
  shared(crossCommRank,pingNodeA,pingNodeB,totalReps)
  {

    /* loop totalRep times */
    for (repIter = 1; repIter <= totalReps; repIter++){

      if (crossCommRank == pingNodeA || crossCommRank == pingNodeB){

	/* Calculate lower bound of each threads portion
	 * of the data array.
	 */
	lBound = (myThreadID * dataSize);

	/* Each thread writes to its part of pingSendBuf */
#pragma omp for nowait schedule(static,dataSize)

	for (i=0; i<sizeofBuffer; i++){
	  pingSendBuf[i] = globalIDarray[myThreadID];
	}

	/* Each thread starts send of dataSize items from
	 * pingSendBuf.
	 */
	MPI_Isend(&pingSendBuf[lBound], dataSize, MPI_INT, \
		  destRank, myThreadID, crossComm, &requestID);

	/* Thread then waits for message from destRank
	 * with tag equal to its threadID.
	 */
	MPI_Recv(&pingRecvBuf[lBound], dataSize, MPI_INT, destRank, \
		 myThreadID, crossComm, &status);


	/* Thread completes send using MPI_Wait */
	MPI_Wait(&requestID, &status);

	/* Each thread reads its part of received buffer. */
#pragma omp for nowait schedule(static,dataSize)

	for (i=0; i<sizeofBuffer; i++){
	  finalRecvBuf[i] = pingRecvBuf[i];
	}


      }

    } /* End repetitions loop */
  }

  return 0;
}

/*-----------------------------------------------------------*/
/* allocateMultiPingpingData                                 */
/*                                                           */
/* Allocates space for the main data arrays.                 */
/* Size of each array is specified by subroutine argument.   */
/*-----------------------------------------------------------*/
int allocateMultiPingpingData(int sizeofBuffer){

  if (crossCommRank == pingNodeA || crossCommRank == pingNodeB){

    pingSendBuf = (int *)malloc(sizeof(int) * sizeofBuffer);
    pingRecvBuf = (int *)malloc(sizeof(int) * sizeofBuffer);
    finalRecvBuf = (int *)malloc(sizeof(int) * sizeofBuffer);
  }

  return 0;
}

/*-----------------------------------------------------------*/
/* freeMultiPingpingData                                     */
/*                                                           */
/* Free allocated memory for main data arrays.               */
/*-----------------------------------------------------------*/
int freeMultiPingpingData(){

  if (crossCommRank == pingNodeA || crossCommRank == pingNodeB){

    free(pingSendBuf);
    free(pingRecvBuf);
    free(finalRecvBuf);

  }

  return 0;
}

/*-----------------------------------------------------------*/
/* testMultiPingping                                         */
/*                                                           */
/* Verifies the the multi-pingping benchmark worked          */
/* correctly.                                                */
/*-----------------------------------------------------------*/
int testMultiPingping(int sizeofBuffer, int dataSize){
  int i;
  int testFlag, localTestFlag;

  /* set localTestFlag to true */
  localTestFlag = TRUE;

  /* Testing done for processes on pingNodeA & pingNodeB */
  if (crossCommRank == pingNodeA || crossCommRank == pingNodeB) {

    /* allocate space for testBuf */
    testBuf = (int *)malloc(sizeof(int) * sizeofBuffer);

    /* Construct testBuf with correct values */
#pragma omp parallel for default(none)					\
  private(i)								\
  shared(otherPingRank,numThreads,dataSize,sizeofBuffer,testBuf)	\
  schedule(static,dataSize)

    for (i=0; i<sizeofBuffer; i++){
      /* calculate globalID of thread expected in finalRecvBuf.
       * This is done using otherPingRank.
       */
      testBuf[i] = (otherPingRank * numThreads) + myThreadID;
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

  /* Reduce testFlag into master with logical AND */
  MPI_Reduce(&localTestFlag, &testFlag, 1, MPI_INT, MPI_LAND, 0, comm);


  /* master sets testOutcome flag */
  if (myMPIRank == 0){
    setTestOutcome(testFlag);
  }

  return 0;
}
