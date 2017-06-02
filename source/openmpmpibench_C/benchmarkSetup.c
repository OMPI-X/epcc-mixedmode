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
/* Setup functions for Mixed mode benchmark program.         */
/* Routines to setup the benchmark.                          */
/*-----------------------------------------------------------*/
#include "benchmarkSetup.h"
#include "parallelEnvironment.h"
/*-----------------------------------------------------------*/
/* openFile                                                  */
/*                                                           */
/* Attempts to open the file passed as an argument.          */
/* --------------------------------------------------------- */
int openFile(char *fileName){

	printf("Attempting to open %s ....",fileName);

	inputFile = fopen(fileName, "r");

	/* Check that file opened successfully */
	if (inputFile == NULL){
		printf("ERROR.\n");
	}
	else{
		printf("Success.\n");
	}

return 0;
}


/*-----------------------------------------------------------*/
/* closeFile                                                 */
/*                                                           */
/* Closes the input file.                                    */
/*-----------------------------------------------------------*/
int closeFile(){

	/* close the input file */
	fclose(inputFile);

return 0;
}

/*-----------------------------------------------------------*/
/* setupBenchmarkList                          				 */
/*                                                           */
/* Subroutine to setup the benchmarkList array with the      */
/* list of all possible benchmarks.                          */
/*-----------------------------------------------------------*/
 int setupBenchmarkList(){

	 /* Pingpong benchmarks */
	 strcpy (benchmarkList[0], "masteronlypingpong");
	 strcpy (benchmarkList[1], "funnelledpingpong");
	 strcpy (benchmarkList[2], "multiplepingpong");
	 /* Pingping benchmarks */
	 strcpy (benchmarkList[3], "masteronlypingping");
	 strcpy (benchmarkList[4], "funnelledpingping");
	 strcpy (benchmarkList[5], "multiplepingping");
	 /* Haloexchange benchmarks */
	 strcpy (benchmarkList[6], "masteronlyhaloexchange");
	 strcpy (benchmarkList[7], "funnelledhaloexchange");
	 strcpy (benchmarkList[8], "multiplehaloexchange");
	 /* Multi-pingpong benchmarks */
	 strcpy (benchmarkList[9], "masteronlymultipingpong");
	 strcpy (benchmarkList[10], "funnelledmultipingpong");
	 strcpy (benchmarkList[11], "multiplemultipingpong");
	 /* Multi-pingpong benchmarks */
	 strcpy (benchmarkList[12], "masteronlymultipingping");
	 strcpy (benchmarkList[13], "funnelledmultipingping");
	 strcpy (benchmarkList[14], "multiplemultipingping");
	 /* Collective benchmarks */
	 strcpy (benchmarkList[15], "barrier");
	 strcpy (benchmarkList[16], "reduce");
	 strcpy (benchmarkList[17], "allreduce");
	 strcpy (benchmarkList[18], "broadcast");
	 strcpy (benchmarkList[19], "scatter");
	 strcpy (benchmarkList[20], "gather");
	 strcpy (benchmarkList[21], "alltoall");

return 0;
 }

/*-----------------------------------------------------------*/
/* readBenchmarkParams                                       */
/*															 */
/* Initialises the benchmark parameters.                     */
/* Reads the minimum and maximum data for the benchmarks     */
/* from the input file (unit 10).                            */
/*-----------------------------------------------------------*/
int readBenchmarkParams(){

	/* Rank 0 reads parameters from input file */
    if (myMPIRank == 0){
    	printf ("Reading parameters from input file....\n");
        /* read minimum data size from input file */
        fscanf(inputFile, "%d", &minDataSize);
        /* read maximum data size from input file */
        fscanf(inputFile, "%d", &maxDataSize);
        /* read target time from input file */
        fscanf(inputFile, "%lf", &targetTime);

        /* set other benchmark parameters */
        warmUpIters = 2;
        defaultReps = 1000;

        /* Report benchmark parameters */
        printf("------------------------------------------\n");
        printf("           Benchmark parameters           \n");
        printf("------------------------------------------\n");
        printf("Minimum data size %d\n", minDataSize);
        printf("Maximum data size %d\n", maxDataSize);
        printf("Target time (sec) %lf\n", targetTime);
        printf("Default Repetitions %d\n", defaultReps);
        printf("No. Warmup iterations %d\n", warmUpIters);

    }
	/*Initialise benchmarkNumber to 0 so that the  WHILE loop in
	the driver is entered the first time */
    benchmarkNumber = 0;

    /* Broadcast benchmark parameters from master to all
    other MPI processes. */
    MPI_Bcast(&minDataSize, 1, MPI_INT, 0, comm);
    MPI_Bcast(&maxDataSize, 1, MPI_INT, 0, comm);
    MPI_Bcast(&targetTime, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&defaultReps, 1, MPI_INT, 0, comm);
    MPI_Bcast(&warmUpIters, 1, MPI_INT, 0, comm);

    return 0;
}

/*-----------------------------------------------------------*/
/* findBenchmarkNumber                                       */
/*                                                           */
/* Finds the ID of the next benchmark which will be          */
/* executed. Master MPI process reads next line from input   */
/* file. It then compares it to the benchmark list to find   */
/* the ID and broadcasts this to the other MPI processes.    */
/* 															 */
/* The function sets the benchmarkNumber variable and also   */
/* returns the benchmarkNumber.  							 */
/*-----------------------------------------------------------*/
int findBenchmarkNumber(){
	char benchmarkName[MAXSTRING];
	int rankInA, rankInB;
	int i;

	/* Master MPI process reads next line from file */
	if (myMPIRank == 0){
		/* set benchmarkNumber to ERROR before read to allow error
		check */
		benchmarkNumber = ERROR;

		/* read next benchmark from file */
		if (fscanf(inputFile, "%s", benchmarkName) == EOF){
			benchmarkNumber = FINISHED;
		}
		else {
			/* convert benchmarkName to lowercase characters */
			convertToLowercase(benchmarkName);
			/* ..and check if benchmark name matches. */
			for (i = 0; i< NUM_BENCHMARKS; i++){
				if (strcmp(benchmarkName,benchmarkList[i]) == 0){
					benchmarkNumber = i;
				}
			}
		}

		/* Check if benchmark Name does not match */
		if (benchmarkNumber == ERROR){
		   printf("ERROR: %s does not match any possible benchmarks\n",benchmarkName);
		}

		/* Check if pingpong or pingping benchmark */
		if (benchmarkNumber <= LASTPPID){
			/* Read ranks from input file */
			if (fscanf(inputFile, "%d %d",&rankInA, &rankInB) != 2){
				printf("ERROR: expecting ranks after %s\n",benchmarkName);
			}
			else {
				PPRanks[0] = findRank(rankInA);
				PPRanks[1] = findRank(rankInB);
			}
			/* Check if PPRanks are the same */
			if (PPRanks[0] == PPRanks[1]){
				printf("Warning: Ranks are the same; benchmark will not work.\n");
			}

		}
	}

	/* Broadcast benchmarkNumber to other MPI processes */
	MPI_Bcast(&benchmarkNumber, 1, MPI_INT, 0, comm);

	/* If pingpong or pingping benchmark then broadcast ranks of participating processes */
	if (benchmarkNumber <= LASTPPID) {
		MPI_Bcast(PPRanks, 2, MPI_INT, 0, comm);
	}

	return benchmarkNumber;
}

/*-----------------------------------------------------------*/
/* convertToLowerCase             							 */
/*                                                           */
/* Takes a string as an agrument and converts all            */
/* uppercase characters to lowercase using its ASCII value.  */
/*-----------------------------------------------------------*/
int convertToLowercase(char *convertString){
	int i;
	int len;

	len = strlen(convertString);

	for (i=0; i<len; i++){
		convertString[i] = tolower(convertString[i]);
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* repTimeCheck                                              */
/*                                                           */
/* Checks if the time for the benchmark reached the target   */
/* time. Changes the number of repetitions for the next      */
/* data size based on the difference between the time        */
/* taken and the target time.                                */
/*-----------------------------------------------------------*/
int repTimeCheck(double time, int numReps){
	int repCheck;

	if (time < targetTime){
		/* double totalReps and repeat benchmark */
		repsToDo = 2 * numReps;
		repCheck = FALSE;
	}
	else if (time > (2 * targetTime)){
		/* finish benchmark and half number of reps for next dataSize */
		repsToDo = max(numReps/2,1);
		repCheck = TRUE;
	}
	else { /* time is >= targetTime */
		/* finish benchmark and keep reps for next data size */
		repCheck = TRUE;
	}

	return repCheck;
}

/*-----------------------------------------------------------*/
/* max                                                       */
/*                                                           */
/* Finds the maximum value of two integers passed as         */
/* arguments.												 */
/*-----------------------------------------------------------*/
int max(int a, int b){

	return (a > b) ? a : b;
}


