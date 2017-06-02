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
/* Implements the collective barrier mixed mode OpenMP/MPI   */
/* benchmark.                                                */
/*-----------------------------------------------------------*/
#include "collective_barrier.h"

/*-----------------------------------------------------------*/
/* barrierDriver                                             */
/*                                                           */
/* Driver subroutine for the barrier benchmark.              */
/*-----------------------------------------------------------*/
int barrierDriver(){
	/* initialise repsToDo to defaultReps */
	repsToDo = defaultReps;

	/* perform warm-up for benchmark */
	barrierKernel(warmUpIters);

	/* Initialise the benchmark */
	benchComplete = FALSE;

	/* Execute benchmark until target time is reached */
	while (benchComplete != TRUE){
		/* Start timer */
		MPI_Barrier(comm);
		startTime = MPI_Wtime();

		/* Execute benchmark for repsToDo repetitions */
		barrierKernel(repsToDo);

		/* Stop timer */
		MPI_Barrier(comm);
		finishTime = MPI_Wtime();
		totalTime = finishTime - startTime;

		/* Test if target time was reached with number of
		 * repetitions.
		 */
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
		/* no unit test, hardwire test result to pass */
		setTestOutcome(TRUE);
		setReportParams(1,repsToDo,totalTime);
		printReport();
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* barrierKernel                                             */
/* 															 */
/* Main kernel for barrier benchmark.                        */
/* First threads under each process synchronise with an      */
/* OMP BARRIER. Then a MPI barrier synchronises each MPI     */
/* process. MPI barrier is called within a OpenMP master     */
/* directive.                                                */
/*-----------------------------------------------------------*/
int barrierKernel(int totalReps){
	int repIter;

	/* Open the parallel region */
#pragma omp parallel default(none) \
	private(repIter) \
	shared(totalReps,comm)
	{
	for (repIter=0; repIter<totalReps; repIter++){

		/* Threads synchronise with an OpenMP barrier */
#pragma omp barrier

		/* Master threads on each process now synchronise */
#pragma omp master
		{
		MPI_Barrier(comm);
		}
	}
	}

	return 0;
}




