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
/* Implements the scatter and gather mixed mode OpenMP/MPI   */
/* benchmarks.                                               */
/*-----------------------------------------------------------*/
#include "collective_scatterGather.h"

/*-----------------------------------------------------------*/
/* scatterGather                                             */
/*                                                           */
/* Driver routine for the scatter benchmark.                 */
/*-----------------------------------------------------------*/
int scatterGather(int benchmarkType){
  int dataSizeIter, bufferSize;

  /* Initialise repsToDo to defaultReps */
  repsToDo = defaultReps;

  /* Start loop over data sizes */
  dataSizeIter = minDataSize; /* initialise dataSizeIter */
  while (dataSizeIter <= maxDataSize){

    /* Calculate buffer size and allocate spaces for
     * data arrays.
     */
    bufferSize = dataSizeIter * numThreads;

    if (benchmarkType == SCATTER){
      allocateScatterGatherData(bufferSize, benchmarkType);
      /* perform benchmark warm-up */
      scatterKernel(warmUpIters, dataSizeIter);
      /* Test if scatter was successful */
      testScatterGather(bufferSize, benchmarkType);
    }
    else if (benchmarkType == GATHER){
      allocateScatterGatherData(bufferSize, benchmarkType);
      /* Perform benchmark warm-up */
      gatherKernel(warmUpIters, dataSizeIter);
      /* Test if gather was successful */
      if (myMPIRank == GATHERROOT){
	testScatterGather(bufferSize*numMPIprocs, benchmarkType);
      }
    }

    /* Initialise the benchmark */
    benchComplete = FALSE;
    /* Execute benchmark until target time is reached */
    while (benchComplete != TRUE){

      /* Start timer */
      MPI_Barrier(comm);
      startTime = MPI_Wtime();

      if (benchmarkType == SCATTER){
	/* Execute scatter for repsToDo repetitions */
	scatterKernel(repsToDo, dataSizeIter);
      }
      else if (benchmarkType == GATHER){
	/* Execute gather for repsToDo repetitions */
	gatherKernel(repsToDo, dataSizeIter);
      }

      /* Stop timer */
      MPI_Barrier(comm);
      finishTime = MPI_Wtime();
      totalTime = finishTime - startTime;

      /* Test if target time was reached */
      if (myMPIRank==0) {
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
    freeScatterGatherData(benchmarkType);

    /* Dobule data size and loop again */
    dataSizeIter = dataSizeIter * 2;
  }

  return 0;
}

/*-----------------------------------------------------------*/
/* scatterKernel                                             */
/*                                                           */
/* Implement the scatter benchmark.                          */
/* Root process first scatters send buffer to other          */
/* processes.                                                */
/* Each thread under a MPI process then reads its portion    */
/* of scatterRecvBuf.                                        */
/*-----------------------------------------------------------*/
int scatterKernel(int totalReps, int dataSize){
  int repIter, i;
  int totalSendBufElems, sendCount, recvCount;

  /* Calculate totalSendBufElems */
  totalSendBufElems = numMPIprocs * numThreads * dataSize;

  /* Calculate sendCount */
  sendCount = dataSize * numThreads;
  recvCount = sendCount;

  for (repIter=0; repIter<totalReps; repIter++){
    /* Master process writes to scatterSendBuf */

    if (myMPIRank == SCATTERROOT){
      for (i=0; i<totalSendBufElems; i++){
	scatterSendBuf[i] = SCATTERSTARTVAL + i;
      }
    }

    /* Scatter the data to other processes */
    MPI_Scatter(scatterSendBuf, sendCount, MPI_INT, \
		scatterRecvBuf, recvCount, MPI_INT, \
		SCATTERROOT, comm);

    /* Each thread now reads its portion of scatterRecvBuf */
#pragma omp parallel for default(none)			\
    private(i)						\
    shared(dataSize,recvCount,finalBuf,scatterRecvBuf)	\
    schedule(static,dataSize)
    for (i=0; i<recvCount; i++){ /* loop over all data in recv buffer */
	finalBuf[i] = scatterRecvBuf[i];
    }
  } /* End of loop over reps */

  return 0;
}

/*-----------------------------------------------------------*/
/* gatherKernel                                              */
/*														     */
/* Implements the gather benchmark.                          */
/* Each thread writes part of its buffer then all data       */
/* is gathered to the master process.                        */
/*-----------------------------------------------------------*/
int gatherKernel(int totalReps, int dataSize){
  int repIter,i;
  int totalRecvBufElems, sendCount, recvCount;
  int startVal;

  /* Calculate totalRecvBufElems */
  totalRecvBufElems = dataSize * numThreads * numMPIprocs;
  /* Each process calculates its send and recv count */
  sendCount = dataSize * numThreads;
  recvCount = sendCount;

  /* Calculate startVal for each process.
   * This is used to find the values of gatherSendBuf.
   */
  startVal = (myMPIRank * sendCount) + GATHERSTARTVAL;

  for (repIter=0; repIter<totalReps; repIter++){
    /* Each thread writes to its portion of gatherSendBuf */
#pragma omp parallel for default(none)			\
  private(i)						\
  shared(gatherSendBuf,startVal,dataSize,sendCount)	\
  schedule(static,dataSize)
    for (i=0; i<sendCount; i++){
      gatherSendBuf[i] = startVal + i;
    }

    /* Gather the data to GATHERROOT */
    MPI_Gather(gatherSendBuf, sendCount, MPI_INT,\
	       gatherRecvBuf, recvCount, MPI_INT,\
	       GATHERROOT, comm);

    /* GATHERROOT process then copies its received data
     * to finalBuf.
     */
    if (myMPIRank == GATHERROOT){
      for (i=0; i<totalRecvBufElems; i++){
	finalBuf[i] = gatherRecvBuf[i];
      }
    }
  }

  return 0;
}

/*-----------------------------------------------------------*/
/* allocateScatterGatherData                                 */
/*                                                           */
/* Allocate memory for main data arrays                      */
/*-----------------------------------------------------------*/
int allocateScatterGatherData(int bufferSize, int benchmarkType){

  if (benchmarkType == SCATTER){
    /* scatterSendBuf is size (bufferSize * numMPIprocs) */
    if (myMPIRank == SCATTERROOT){
      scatterSendBuf = (int *) malloc((bufferSize * numMPIprocs) * sizeof(int));
    }
    scatterRecvBuf = (int *) malloc(bufferSize * sizeof(int));
    finalBuf = (int *)malloc(bufferSize * sizeof(int));
  }
  else if (benchmarkType == GATHER){
    gatherSendBuf = (int *) malloc(bufferSize * sizeof(int));
    if (myMPIRank == GATHERROOT){
      gatherRecvBuf = (int *) malloc((bufferSize * numMPIprocs) * sizeof(int));
      finalBuf = (int *) malloc((bufferSize * numMPIprocs) * sizeof(int));
    }
  }

  return 0;
}

/*-----------------------------------------------------------*/
/* freeScatterGatherData                                     */
/*                                                           */
/* Free memory of main data arrays.                          */
/*-----------------------------------------------------------*/
int freeScatterGatherData(int benchmarkType){

  if (benchmarkType == SCATTER){
    if (myMPIRank == SCATTERROOT){
      free(scatterSendBuf);
    }
    free(scatterRecvBuf);
    free(finalBuf);
  }
  else if (benchmarkType == GATHER){
    free(gatherSendBuf);
    if (myMPIRank == GATHERROOT){
      free(gatherRecvBuf);
      free(finalBuf);
    }
  }


  return 0;
}

/*-----------------------------------------------------------*/
/* testScatterGather                                         */
/*                                                           */
/* Verifies that the scatter and gahter benchmarks worked    */
/* correctly.                                                */
/*-----------------------------------------------------------*/
int testScatterGather(int sizeofBuffer, int benchmarkType){
  int i, startVal;
  int testFlag, reduceFlag;
  int *testBuf;

  /* Initialise testFlag to true */
  testFlag = TRUE;

  /* Allocate space for testBuf */
  testBuf = (int *) malloc (sizeofBuffer * sizeof(int));

  if (benchmarkType == SCATTER){
    /* Find the start scatter value for each MPI process */
    startVal = (myMPIRank * sizeofBuffer) + SCATTERSTARTVAL;
  }
  else if (benchmarkType == GATHER){
    /* startVal is GATHERSTARTVAL */
    startVal = GATHERSTARTVAL;
  }

  /* Fill testBuf with correct values */
  for (i=0; i<sizeofBuffer; i++){
    testBuf[i] = startVal + i;
  }

  /* Compare each element of finalBuf with testBuf */
  for (i=0; i<sizeofBuffer; i++){
    if (finalBuf[i] != testBuf[i]){
      testFlag = FALSE;
    }
  }

  /* For scatter: reduce testFlag into master with
   * logical AND operator.
   */
  if (benchmarkType == SCATTER){
    MPI_Reduce(&testFlag, &reduceFlag, 1, MPI_INT, MPI_LAND, 0, comm);
    /* Master then sets testOutcome using reduceFlag */
    if (myMPIRank == 0){
      setTestOutcome(reduceFlag);
    }
  }
  else if (benchmarkType == GATHER){
    setTestOutcome(testFlag);
  }

  free(testBuf);

  return 0;
}
