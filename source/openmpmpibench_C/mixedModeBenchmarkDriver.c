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
/* Main driver for mixed mode benchmark program.             */
/* Reads benchmark input file.                               */
/* Initialises the parallel environment.                     */
/* Calls each benchmark.                                     */
/*-----------------------------------------------------------*/
#include "parallelEnvironment.h"
#include "benchmarkSetup.h"
#include "pt_to_pt_pingpong.h"
#include "pt_to_pt_pingping.h"
#include "pt_to_pt_multiPingpong.h"
#include "pt_to_pt_multiPingping.h"
#include "pt_to_pt_haloexchange.h"
#include "collective_barrier.h"
#include "collective_broadcast.h"
#include "collective_reduction.h"
#include "collective_alltoall.h"
#include "collective_scatterGather.h"

int main(int argc, char *argv[]){

	/* Flag to check if benchmark is supported */
	int supportFlag;
	/* String for setting benchmark name for output */
	char name[MAXSTRING];

	/* Initialise the parallel execution environment */
	initParallelEnv();

	/* Master MPI process..... */
	if (myMPIRank == 0){
		/* Check command line argument */
		if (argc != 2){
			printf("ERROR Reading input file from command line.\n");
			printf("Usage: %s <filename>", argv[0] );
			/* Finalise programming environment and exit */
			finaliseParallelEnv();
			exit(-1);
		}
		else{
			/* Print header and parallel environment info. */
			printHeader(numMPIprocs,numThreads,threadSupport);
			/* Open the input file */
			openFile(argv[1]);
			/* Setup the list of all possible benchmarks */
			setupBenchmarkList();
		}
	}

	/* Master reads parameters from input file and
	broadcasts them to the other processes. */
	readBenchmarkParams();

	/* Execute benchmarks by reading list from
	input file. */
	 while (findBenchmarkNumber() != FINISHED){
		 switch (benchmarkNumber){
			/* Masteronly Pingpong */
			case 0:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
				    strcpy(name,"Masteronly Pingpong");
				    setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				pingPong(MASTERONLY);
				break;
			/* Funnelled Pingpong */
			case 1:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Funnelled Pingpong");
					setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				pingPong(FUNNELLED);
				break;
			/* Multiple Pingpong */
			case 2:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_MULTIPLE);
					strcpy(name,"Multiple Pingpong");
					setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				pingPong(MULTIPLE);
				break;

			/* Masteronly Pingping */
			case 3:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Masteronly Pingping");
					setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				pingPing(MASTERONLY);
				break;
			/* Funnelled Pingping */
			case 4:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Funnelled Pingping");
					setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				pingPing(FUNNELLED);
				break;
			/* Multiple Pingping */
			case 5:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_MULTIPLE);
					strcpy(name,"Multiple Pingping");
					setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				pingPing(MULTIPLE);
				break;

			/* Masteronly Haloexchange */
			case 6:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Masteronly Haloexchange");
					setBenchName(name, benchmarkNumber, supportFlag);
					printBenchHeader();
				}
				/* Execute benchmark */
				haloExchange(MASTERONLY);
				break;
			/* Funnelled Haloexchange */
			case 7:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Funnelled Haloexchange");
					setBenchName(name, benchmarkNumber, supportFlag);
					printBenchHeader();
				}
				/* Execute benchmark */
				haloExchange(FUNNELLED);
				break;
			/* Multiple Haloexchange */
			case 8:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_MULTIPLE);
					strcpy(name,"Multiple Haloexchange");
					setBenchName(name, benchmarkNumber, supportFlag);
					printBenchHeader();
				}
				/* Execute benchmark */
				haloExchange(MULTIPLE);
				break;
			/* Masteronly Multi-Pingpong */
			case 9:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Masteronly MultiPingpong");
					setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				multiPingPong(MASTERONLY);
				break;
			/* Funnelled Multi-Pingpong */
			case 10:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Funnelled MultiPingpong");
					setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				multiPingPong(FUNNELLED);
				break;
			/* Multiple Multi-Pingpong */
			case 11:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_MULTIPLE);
					strcpy(name,"Multiple MultiPingpong");
					setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				multiPingPong(MULTIPLE);
				break;
			/* Masteronly Multi-Pingping */
			case 12:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Masteronly MultiPingping");
					setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				multiPingping(MASTERONLY);
				break;
			/* Funnelled Multi-Pingping */
			case 13:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Funnelled MultiPingping");
					setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				multiPingping(FUNNELLED);
				break;
			/* Multiple Multi-Pingping */
			case 14:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_MULTIPLE);
					strcpy(name, "Multiple MultiPingping");
					setBenchName(name, benchmarkNumber, supportFlag);
				}
				/* Execute benchmark */
				multiPingping(MULTIPLE);
				break;
			/* Barrier */
			case 15:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Barrier");
					setBenchName(name, benchmarkNumber, supportFlag);
					printBenchHeader();
				}
				/* Execute benchmark */
				barrierDriver();
				break;
			/* Reduce */
			case 16:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Reduce");
					setBenchName(name, benchmarkNumber, supportFlag);
					printBenchHeader();
				}
				/* Execute benchmark */
				reduction(REDUCE);
				break;
			/* All-Reduce */
			case 17:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"All Reduce");
					setBenchName(name, benchmarkNumber, supportFlag);
					printBenchHeader();
				}
				/* Execute benchmark */
				reduction(ALLREDUCE);
				break;
			/* Broadcast */
			case 18:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Broadcast");
					setBenchName(name, benchmarkNumber, supportFlag);
					printBenchHeader();
				}
				/* Execute benchmark */
				broadcast();
				break;
			/* Scatter */
			case 19:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Scatter");
					setBenchName(name, benchmarkNumber, supportFlag);
					printBenchHeader();
				}
				/* Execute benchmark */
				scatterGather(SCATTER);
				break;
			/* Gather */
			case 20:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"Gather");
					setBenchName(name, benchmarkNumber, supportFlag);
					printBenchHeader();
				}
				/* Execute benchmark */
				scatterGather(GATHER);
				break;
			/* All-to-all */
			case 21:
				/* Set name */
				if (myMPIRank == 0){
					supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED);
					strcpy(name,"All to All");
					setBenchName(name, benchmarkNumber, supportFlag);
					printBenchHeader();
				}
				/* Execute benchmark */
				alltoall();
				break;
		}

	}

	/* Finalise programming environment */
	finaliseParallelEnv();
	if (myMPIRank == 0){
	     closeFile();
	}

}
