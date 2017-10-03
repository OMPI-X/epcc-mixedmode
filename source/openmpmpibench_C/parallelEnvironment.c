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
/* Contains variables and routines for the parallel          */
/* environment.                                              */
/* Routines to setup the MPI and OpenMP programming          */
/* environment.                                              */
/* Header file = parallelEnvironment.h                       */
/*-----------------------------------------------------------*/
#define _GNU_SOURCE
#include "parallelEnvironment.h"
#include "output.h"
#include <pmix.h>
#include <unistd.h>
#include <stdio.h>

#define PMIX_PROGRAMMING_MODEL      "pmix.pgm.model"        // (char*) programming model being initialized (e.g., "MPI" or "OpenMP")
#define PMIX_MODEL_LIBRARY_NAME     "pmix.mdl.name"         // (char*) programming model implementation ID (e.g., "OpenMPI" or "MPICH")
#define PMIX_MODEL_LIBRARY_VERSION  "pmix.mld.vrs"          // (char*) programming model version string (e.g., "2.1.1")
#define PMIX_THREADING_MODEL        "pmix.threads"          // (char*) threading model used (e.g., "pthreads")
#define PMIX_LOCAL_SIZE             "pmix.local.size"

/*-----------------------------------------------------------*/
/* initParallelEnv			                                 */
/*                                                           */
/* Initialises the MPI and OpenMP environments.              */
/* Finds the total number of MPI processes and               */
/* OpenMP threads.                                           */
/* Also finds the ID of each MPI process and OpenMP thread.  */
/*-----------------------------------------------------------*/

typedef struct {
    pmix_info_t *info;
    size_t      ninfo;
    int         *active;
    char        *nspace;
} mydata_t;

int active  = -1;
int ready   = -1;
int all_set = -1;
pmix_proc_t myproc;

uint32_t    n_node_local_ranks = 0;

static void release_fn(size_t evhdlr_registration_id,
                       pmix_status_t status,
                       const pmix_proc_t *source,
                       pmix_info_t info[], size_t ninfo,
                       pmix_info_t results[], size_t nresults,
                       pmix_event_notification_cbfunc_fn_t cbfunc,
                       void *cbdata)
{
    /* tell the event handler state machine that we are the last step */
    if (NULL != cbfunc) {
        cbfunc(PMIX_EVENT_ACTION_COMPLETE, NULL, 0, NULL, NULL, cbdata);
    }

    /* do whatever we want/need to do to coordinate */
}

static void
info_event_cb (pmix_status_t status, void *cbdata)
{
    mydata_t        *cd     = (mydata_t*)cbdata;
    size_t          n_keys;
    int             i;
    pmix_proc_t     proc;
    pmix_value_t    *val    = NULL;

    //fprintf (stdout, "%s - Start\n", __func__);

    n_keys = cd->ninfo;
    //fprintf (stdout, "%s - Nb of info keys: %d\n", __func__, (int)n_keys);
    for (i = 0; i < n_keys; i++)
    {
        //fprintf (stdout, "%s - %s/%s\n", __func__, cd->info[i].key, cd->info[i].value.data.string);
        if (strcmp (cd->info[i].key, "pmix.pgm.model") == 0 &&
            strcmp (cd->info[i].value.data.string, "OpenMP") == 0)
        {
            fprintf (stdout, "%s - Just detected an OpenMP runtime!!\n", __func__);
        }

        if (strcmp (cd->info[i].key, "pmix.local.size") == 0)
        {
            n_node_local_ranks = cd->info[i].value.data.uint32;
            fprintf (stdout, "%s - Just detected the number of node local MPI ranks: %d\n", __func__, (int)n_node_local_ranks);
        }
    }    
    PMIX_INFO_FREE (cd->info, cd->ninfo);
    free (cd);

    //fprintf (stdout, "%s - End\n", __func__);
}

static void
mpiinfo_event_cb (pmix_status_t status, void *cbdata)
{
    mydata_t    *cd = (mydata_t*)cbdata;
    
}

static void
evhandler_reg_callbk_2 (pmix_status_t status, size_t evhandler_ref, void *cbdata)
{
   ready = 1; 
}

static void
evhandler_reg_callbk (pmix_status_t status, size_t evhandler_ref, void *cbdata)
{
    char            *prog_model = NULL;
    pmix_proc_t     proc;
    pmix_value_t    *val        = NULL;
    pmix_rank_t     npeers;

    active = status;
    fprintf (stdout, "%s - Start\n", __func__);

#if 0
    PMIX_PROC_CONSTRUCT (&proc);
    (void)strncpy (proc.nspace, myproc.nspace, PMIX_MAX_NSLEN);
    proc.rank = PMIX_RANK_WILDCARD;

    fprintf (stdout, "%s - Getting PMIX_PROGRAMMING_MODEL\n", __func__);
    PMIx_Get (&proc, PMIX_PROGRAMMING_MODEL, NULL, 0, &val);
    if (val == NULL)
    {
        fprintf (stderr, "PMIX_PROGRAMMING_MODEL undefined\n");
    } else {
        prog_model = val->data.string;
        fprintf (stdout, "PMIX_PROGRAMMING_MODEL: %s\n", prog_model);
    }

    fprintf (stdout, "%s - Getting PMIX_LOCAL_SIZE\n", __func__);
    PMIx_Get (&proc, PMIX_LOCAL_SIZE, NULL, 0, &val);
    if (val == NULL)
    {
        fprintf (stderr, "PMIX_LOCAL_SIZE undefined\n");
    } else {
        npeers = val->data.uint32;
        fprintf (stdout, "PMIX_LOCAL_SIZE: %d\n", (int)npeers);
    }
#endif

    fprintf (stdout, "%s - End\n", __func__);
}

int initParallelEnv(){
    int             i;
    int             numCPU;
    mydata_t        *mpi_data;
    mydata_t        *proc_data;
    pmix_status_t   rc;
    uint32_t        node_local_procs    = 0;
    pmix_info_t     *mpi_info;
    char            *mpi_model          = "MPI";
    char            *mpi_modelname      = "openMPI";
    char            *mpi_version        = "master";
    pmix_status_t   code                = PMIX_MODEL_DECLARED;
    char            *s_procNames        = NULL;
    char            *r_procNames        = NULL;
    char            *n_threads_str      = NULL;

    /* Setup MPI programming environment */
    PMIx_Init (&myproc, NULL, 0);

    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &threadSupport);
    fprintf (stderr, "[%s:%d:%s] Check\n", __FILE__, __LINE__, __func__);

    /* Find the processor name of each MPI process */
    MPI_Get_processor_name(myProcName, &procNameLen);

    mpi_data = (mydata_t*) malloc (sizeof (mydata_t));

    /* Register a handler specifically for notifying us if other programming models declare themselves */
    active = -1;
    mpi_data->active = &active;
    //mpi_data->nspace = 
    PMIx_Register_event_handler (&code, 1, NULL, 0, release_fn, evhandler_reg_callbk, NULL);
    fprintf (stderr, "[%s:%d:%s] Check\n", __FILE__, __LINE__, __func__);

#if 0
    /* Setup to declare our programming model */
    mpi_data = (mydata_t*) malloc (sizeof (mydata_t));
    mpi_data->ninfo = 3;
    PMIX_INFO_CREATE (mpi_data->info, mpi_data->ninfo);
    PMIX_INFO_LOAD (&mpi_data->info[0], PMIX_PROGRAMMING_MODEL, mpi_model, PMIX_STRING);
    PMIX_INFO_LOAD (&mpi_data->info[1], PMIX_MODEL_LIBRARY_NAME, mpi_modelname, PMIX_STRING);
    PMIX_INFO_LOAD (&mpi_data->info[2], PMIX_MODEL_LIBRARY_VERSION, mpi_version, PMIX_STRING);
    fprintf (stderr, "[%s:%d:%s] Check\n", __FILE__, __LINE__, __func__);

    PMIx_Notify_event (PMIX_MODEL_DECLARED, NULL, PMIX_RANGE_PROC_LOCAL, mpi_data->info, mpi_data->ninfo, info_event_cb, (void*)mpi_data);
#endif
    fprintf (stdout, "Client ns %s rank %d: Running\n", myproc.nspace, myproc.rank);

#if 0
    while (-1 == active)
        usleep (10);

    fprintf (stderr, "[%s:%d:%s] Check\n", __FILE__, __LINE__, __func__);
    if (0 != active)
        exit (active);
#endif

	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &numMPIprocs);
	MPI_Comm_rank(comm, &myMPIRank);

	/*Find the number of bytes for an int */
	sizeInteger = sizeof(int);

	/* Find the processor name of each MPI process */
	MPI_Get_processor_name(myProcName, &procNameLen);

    /* Use processor name to create a communicator
     * across node boundaries.
     */
    setupCommunicators();

    /* Figure out the nunber of ranks on the node */
    /* Not portable but good enough for now */
    s_procNames = (char*) malloc (procNameLen * numMPIprocs);
    for (i = 0; i < numMPIprocs; i++)
    {
        strncpy (&s_procNames[i*procNameLen], myProcName, procNameLen);
    }
    r_procNames = (char*) malloc (procNameLen * numMPIprocs);
    MPI_Alltoall (s_procNames, procNameLen, MPI_CHAR,
                  r_procNames, procNameLen, MPI_CHAR,
                  comm);
    for (i = 0; i < numMPIprocs; i++)
    {
        if (strncmp (&r_procNames[i*procNameLen], myProcName, procNameLen) == 0)
        {
            node_local_procs++;
        }
    }

    /* Figure out the number of processors */
    numCPU = sysconf (_SC_NPROCESSORS_ONLN);
    printf ("numCPU: %d\n", numCPU);

    asprintf (&n_threads_str, "%d", (int)(numCPU - node_local_procs));
    setenv ("OMP_NUM_THREADS", n_threads_str, 1); 
    free (n_threads_str);
    n_threads_str = NULL;

    /* Setup to declare our programming model */
    mpi_data = (mydata_t*) malloc (sizeof (mydata_t));
    mpi_data->ninfo = 4;
    PMIX_INFO_CREATE (mpi_data->info, mpi_data->ninfo);
    PMIX_INFO_LOAD (&mpi_data->info[0], PMIX_PROGRAMMING_MODEL, mpi_model, PMIX_STRING);
    PMIX_INFO_LOAD (&mpi_data->info[1], PMIX_MODEL_LIBRARY_NAME, mpi_modelname, PMIX_STRING);
    PMIX_INFO_LOAD (&mpi_data->info[2], PMIX_MODEL_LIBRARY_VERSION, mpi_version, PMIX_STRING);
    PMIX_INFO_LOAD (&mpi_data->info[3], PMIX_LOCAL_SIZE, &node_local_procs, PMIX_UINT32);
    fprintf (stderr, "[%s:%d:%s] Check\n", __FILE__, __LINE__, __func__);

    PMIx_Notify_event (PMIX_MODEL_DECLARED, NULL, PMIX_RANGE_PROC_LOCAL, mpi_data->info, mpi_data->ninfo, info_event_cb, (void*)mpi_data);

    while (-1 == active)
        usleep (10);
    if (0 != active)
        exit (active);
#if 0
    PMIx_Register_event_handler (&code, 1, NULL, 0, release_fn, evhandler_reg_callbk_2, NULL);

    proc_data = (mydata_t*) malloc (sizeof (mydata_t));
    proc_data->ninfo = 1;
    PMIX_INFO_CREATE (proc_data->info, proc_data->ninfo);
    PMIX_INFO_LOAD (&proc_data->info[0], PMIX_LOCAL_SIZE, &node_local_procs, PMIX_UINT32);

    PMIx_Notify_event (PMIX_MODEL_DECLARED, NULL, PMIX_RANGE_PROC_LOCAL, proc_data->info, proc_data->ninfo, mpiinfo_event_cb, (void*)proc_data);
#endif
    

    free (r_procNames);
    free (s_procNames);
    r_procNames = NULL;
    s_procNames = NULL;

/*
    while (-1 == ready)
        usleep (10);
*/

	/* setup OpenMP programming environment */
{
        mydata_t    *omp_data;
        char        *omp_model      = "OpenMP";
        char        *omp_modelname  = "StdOpenMP";
        char        *omp_version    = "TBD";

        omp_data = (mydata_t*) malloc (sizeof (mydata_t));

        omp_data->ninfo = 4;
        PMIX_INFO_CREATE (omp_data->info, omp_data->ninfo);
        PMIX_INFO_LOAD (&omp_data->info[0], PMIX_PROGRAMMING_MODEL, omp_model, PMIX_STRING);
        PMIX_INFO_LOAD (&omp_data->info[1], PMIX_MODEL_LIBRARY_NAME, omp_modelname, PMIX_STRING);
        PMIX_INFO_LOAD (&omp_data->info[2], PMIX_MODEL_LIBRARY_VERSION, omp_version, PMIX_STRING);
        PMIX_INFO_LOAD (&omp_data->info[3], PMIX_THREADING_MODEL, "openmp", PMIX_STRING);

        PMIx_Notify_event (PMIX_MODEL_DECLARED, NULL, PMIX_RANGE_PROC_LOCAL, omp_data->info, omp_data->ninfo, info_event_cb, (void*)omp_data);

        if (n_node_local_ranks > 0)
        {
            numThreads = numCPU - n_node_local_ranks;
            omp_set_dynamic (0);
            omp_set_num_threads (numThreads);
        }
}
#pragma omp parallel default(none) \
shared(numThreads,globalIDarray,myMPIRank,n_node_local_ranks,numCPU)
   {
        mydata_t    *omp_data;
        char        *omp_model      = "OpenMP";
        char        *omp_modelname  = "StdOpenMP";
        char        *omp_version    = "TBD";

        omp_data = (mydata_t*) malloc (sizeof (mydata_t));

        omp_data->ninfo = 4;
        PMIX_INFO_CREATE (omp_data->info, omp_data->ninfo);
        PMIX_INFO_LOAD (&omp_data->info[0], PMIX_PROGRAMMING_MODEL, omp_model, PMIX_STRING);
        PMIX_INFO_LOAD (&omp_data->info[1], PMIX_MODEL_LIBRARY_NAME, omp_modelname, PMIX_STRING);
        PMIX_INFO_LOAD (&omp_data->info[2], PMIX_MODEL_LIBRARY_VERSION, omp_version, PMIX_STRING);
        PMIX_INFO_LOAD (&omp_data->info[3], PMIX_THREADING_MODEL, "openmp", PMIX_STRING);

        PMIx_Notify_event (PMIX_MODEL_DECLARED, NULL, PMIX_RANGE_PROC_LOCAL, omp_data->info, omp_data->ninfo, info_event_cb, (void*)omp_data);

	    numThreads = omp_get_num_threads();
	    myThreadID = omp_get_thread_num();

	   /* Allocate space for globalIDarray */
#pragma omp single
   {
	   globalIDarray = (int *)malloc(numThreads * sizeof(int));
   }

	   /*calculate the globalID for each thread */
	   globalIDarray[myThreadID] = (myMPIRank * numThreads) + myThreadID;
   }

   /* set parallel info in benchmark report type */
   setParallelInfo(numMPIprocs,threadSupport,numThreads);

return 0;
}

/*-----------------------------------------------------------*/
/* finaliseParallelEnv                                       */
/*                                                           */
/* Closes the MPI programming environment.                   */
/*                                                         	 */
/*-----------------------------------------------------------*/
int finaliseParallelEnv(){

	/* finalise the MPI programming environment */
    MPI_Finalize();
    /*free the space created for globalIDarray...*/
    free(globalIDarray);

return 0;
}

/*-----------------------------------------------------------*/
/* findRank                                                  */
/*                                                           */
/* Finds the MPI ranks which will take part in the pingping  */
/* or pingpong benchmarks based on the numbers read from the */
/* input file.                                               */
/*-----------------------------------------------------------*/
int findRank(int rankIn){
	int CalcRank;

    /* Figure out actual MPI rank */
    if (rankIn < 0){
    	CalcRank = numMPIprocs + rankIn;
    }
    else{
    	CalcRank = rankIn;
    }

    /* Check if findRank is too big or still -ve */
    if (CalcRank > (numMPIprocs-1)){
       printf("Warning: Rank input greater than total process count.\n");
       printf("Using Rank = %d ", numMPIprocs-1);
       CalcRank = numMPIprocs - 1;
    }
    else if(CalcRank < 0){
    	printf("Warning: MPI process offset greater than total process count.\n");
    	printf("Using Rank = 0 ");
    	CalcRank = 0;
    }

    return CalcRank;
}

/*-----------------------------------------------------------*/
/* findNeighbourRanks                          				 */
/*															 */
/* This creates a cartesian topology and finds the left      */
/* and right neighbours of each process.                     */
/*-----------------------------------------------------------*/
int findNeighbours(){
	int dims[1];
	int periods[1];
	int reorder;

	/* find a good process distribution */
	dims[0] = 0; /* zero so that dims_create tries to rearrange */
	MPI_Dims_create(numMPIprocs, 1, dims);

	/* set periods equal to TURE for periodic boundary conditions ... */
	periods[0] = TRUE;
	/* ...and reorder = FALSE */
	reorder = FALSE;

	/* Create the cartesian topology */
	MPI_Cart_create(comm, 1, dims, periods, reorder, &commCart);

	/* Find the ranks of the left and right neighbour */
	MPI_Cart_shift(commCart, 0, 1, &leftNeighbour, &rightNeighbour);

	return 0;
}

/*-----------------------------------------------------------*/
/* benchmarkSupport                                          */
/*                                                           */
/* This function compares the level of thread support        */
/* needed by a particular benchmark with the level provided  */
/* by the implementation.                                    */
/*-----------------------------------------------------------*/
int benchmarkSupport(int required){
    int benchSupport;

    if (required <= threadSupport){
       benchSupport = TRUE;
    }
    else {
       benchSupport = FALSE;
    }

    return benchSupport;
}

/*-----------------------------------------------------------*/
/* compareProcNames                                          */
/*                                                           */
/* Compares the names of 2 processes to check if they are on */
/* the same node or not.                                     */
/*-----------------------------------------------------------*/
int compareProcNames(int rankA, int rankB){
	int sameNode;
	char recvProcName[MPI_MAX_PROCESSOR_NAME];

	/* Rank B sends procName to Rank A */
	if (myMPIRank == rankB){
		MPI_Send(myProcName, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, rankA, TAG, comm);
	}
	else if (myMPIRank == rankA){
		MPI_Recv(recvProcName, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, rankB, TAG, comm, &status);
		/* Rank B compares the two processor names */
		if (strcmp(myProcName,recvProcName) == 0){
			sameNode = TRUE;
		}
		else{
			sameNode = FALSE;
		}
	}

	/* Rank A then broadcasts its sameNode value to the other processes */
	MPI_Bcast(&sameNode, 1, MPI_INT, rankA, comm);

	return sameNode;
}

/*-----------------------------------------------------------*/
/* setupCommunicators                                        */
/*                                                           */
/* This creates two new communicators.                       */
/* The first gives a local communicator for processes on     */
/* the same node.                                            */
/* The second uses the local rank to give a communicator     */
/* across node boundaries.                                   */
/*                                                           */
/* e.g. for 16 nodes each with 2 processors, this routine    */
/*      will give 16 local communicators of size 2 and       */
/*      2 communicators of size 2 across nodes.              */
/*-----------------------------------------------------------*/
int setupCommunicators(){
	int procHash;

	/* Get hash from processor name */
	procHash = procNameToHash();

	/* Comm_split using procHash as colour to get
	 * local communicator.
	 */
	MPI_Comm_split(comm, procHash, 0, &localComm);

	/* Find ranks of processes in localComm */
	MPI_Comm_rank(localComm, &localCommRank);

	/* Find the size of localComm (for use in calculating multi datasize) */
	MPI_Comm_size(localComm, &localCommSize);

	/* Use localRank as colour to get communicator across nodes. */
	MPI_Comm_split(comm, localCommRank, 0, &crossComm);

	/* Find ranks of processes in crossComm */
	MPI_Comm_rank(crossComm, &crossCommRank);

	return 0;
}

/*-----------------------------------------------------------*/
/* procNameToHash                                            */
/*                                                           */
/* Creates an integer hash for each process.                 */
/* Each process on the same node will have the same hash     */
/* value.                                                    */
/*-----------------------------------------------------------*/
int procNameToHash(){
	int procHash,i;

	/* Initialise hash to 0 */
	procHash = 0;

	for (i=0; i<procNameLen; i++){
		procHash = (7 * procHash) + (int)(myProcName[i]);
	}

	return procHash;
}


/*-----------------------------------------------------------*/
/* exchangeWorldRanks                                        */
/*                                                           */
/* Finds the MPI_COMM_WORLD ranks of the processes           */
/* participating in the multi-pingpong and multi-pingping    */
/* benchmarks.                                               */
/*-----------------------------------------------------------*/
int exchangeWorldRanks(int nodeA, int nodeB, int *otherWorldRank){
	int destRank;

	if (crossCommRank == nodeA){
		destRank = nodeB;
	}
	else if (crossCommRank == nodeB){
		destRank = nodeA;
	}

	if (crossCommRank == nodeA || crossCommRank == nodeB){
		/* Start send of comm_world rank to destRank in crossComm */
		MPI_Isend(&myMPIRank, 1, MPI_INT, destRank, TAG, crossComm, &requestID);

		/* Then wait for message from destRank and store in otherWorldRank. */
		MPI_Recv(otherWorldRank, 1, MPI_INT, destRank, TAG, crossComm, &status);

		MPI_Wait(&requestID, &status);
	}

	return 0;
}

/*-----------------------------------------------------------*/
/* sendProcName                                              */
/*                                                           */
/* Sends the processor name from processes in destNode       */
/* of crossComm to srcNode.                                  */
/*-----------------------------------------------------------*/
int sendProcName(int destNode, int srcNode, char *destProcName){

	/* MPI processes under srcNode of crossComm send their
	 * processor name to destNode.
	 */
	if (crossCommRank == srcNode){
		MPI_Send(myProcName, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, \
				destNode, TAG, crossComm);
	}
	else if (crossCommRank == destNode){
		MPI_Recv(destProcName, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, \
				srcNode, TAG, crossComm, &status);
	}
}

/*-----------------------------------------------------------*/
/* checkCrossCommBalance                                     */
/*                                                           */
/* Checks if there's a balance in the number of processes    */
/* in crossComm nodes.                                       */
/*-----------------------------------------------------------*/
int crossCommBalance(int nodeA, int nodeB){
	int localCommSize, otherLocalCommSize;
	int crossCommBalance;

	/* Find the size of localComm */
	MPI_Comm_size(localComm, &localCommSize);

	/* Master process on nodeB sends localCommSize */
	if ((crossCommRank == nodeB) && (localCommRank == 0)){
		MPI_Send(&localCommSize, 1, MPI_INT, nodeA, TAG, crossComm);
	}
	/* Master process on nodeA... */
	else if ((crossCommRank == nodeA) && (localCommRank == 0)){
		/* 1) receives nodeB's localCommSize */
		MPI_Recv(&otherLocalCommSize, 1, MPI_INT, nodeB, TAG, \
				crossComm, &status);

		/* 2) Test for balance by comparing otherLocalCommSize
		 * to localCommSize.
		 */
		if (localCommSize == otherLocalCommSize){
			/* Set balance to TRUE */
			crossCommBalance = TRUE;
		}
		else{
			crossCommBalance = FALSE;
		}

		/* 3) Send balance to master commWorld process.
		 * Only need explicit send if commWorld is not same
		 * process as master process on nodeA.
		 */
		if (myMPIRank != 0){
			MPI_Send(&crossCommBalance, 1, MPI_INT, 0, TAG, comm);
		}
	}

	/* Master commWorld process.. */
	if (myMPIRank == 0){
		/* Receives balance variable if not same process as
		 * master process on nodeA.
		 */
		if ((crossCommRank != nodeA) && (localCommRank != 0)){
			MPI_Recv(&crossCommRank, 1, MPI_INT, MPI_ANY_SOURCE, \
					TAG, comm, &status);
		}
	}

	/* Broadcast balance to all processes */
	MPI_Bcast(&crossCommBalance, 1, MPI_INT, 0, comm);

	return crossCommBalance;
}
