# EPCC Summary/Notes

There is both a C and Fortran implementation of the benchmarks.  I have
only played with the C version.

There are three different forms used:
  1. "master-only" - master MPI communicates outside omp-parallel
  2. "funnelled" - master MPI communicates inside omp-parallel
  3. "multiple" - all (threads) MPI communicate inside omp-parallel

There are 22 tests, many are just the same operation using the three
different forms (above). For example, MasteronlyHaloExchange,
FunnelledHaloExchange, MultipleHaloExchange.

Also, the PingPong and PingPing tests have a standard two rank and a
multiple (all) ranks version, e.g., MasteronlyPingPong,
MasteronlyMultiPingPong.

PingPong vs. PingPing: The PingPing is just that messages are exchanged in
both directions between the two processes concurrently.

The general groupings & descriptions for the tests (excerpts from [2]):

 - {masteronly,funnelled,multiple} Ping Pong

     Pingpong between two MPI processes with reads and writes by all OpenMP
     threads under each process.  
     With variant ...
        masteronly* - MPI communications outside OpenMP parallel region
         funnelled* - MPI communications inside  OpenMP parallel region,
                      by master thread
          multiple* - all threads take part in MPI communication

  - {masteronly,funnelled,multiple} Ping Ping

     Two MPI processes send a message to each other with parallel reads and
     writes by all OpenMP threads under each process. 
     With variant ...
        masteronly* - MPI communications outside OpenMP parallel region
         funnelled* - MPI communications inside  OpenMP parallel region,
                      by master thread
          multiple* - all threads take part in MPI communication

  - {masteronly,funnelled,multiple} Multi Ping Pong

     Every core on a node performs a pingpong with the corresponding core
     on another node. 
     With variant ...
        masteronly* - MPI communications outside OpenMP parallel region
         funnelled* - MPI communications inside  OpenMP parallel region,
                      by master thread
          multiple* - all threads take part in MPI communication

  - {masteronly,funnelled,multiple} Multi Ping Ping

     Every core on a node performs a pingping with the corresponding core
     on another node. 
     With variant ...
        masteronly* - MPI communications outside OpenMP parallel region
         funnelled* - MPI communications inside  OpenMP parallel region,
                      by master thread
          multiple* - all threads take part in MPI communication

  - {masteronly,funnelled,multiple} Halo Exchange

     All MPI processes are arranged in a ring and each process exchanges
     messages with its two neighbouring processes with parallel reads and
     writes the OpenMP threads under each process. 
     With variant ...
        masteronly* - MPI communications outside OpenMP parallel region
         funnelled* - MPI communications inside  OpenMP parallel region,
                      by master thread
          multiple* - all threads take part in MPI communication

 - Barrier

    Mixed-mode barrier where threads under each process first synchronise
    with an OMP BARRIER followed by an MPI_Barrier to synchronise each MPI
    process.

 - Reduce

    Mixed mode reduce benchmark. All threads under every MPI process
    combines its local buffer. All MPI processes then combined their
    values to get the overall reduction value which is stored on the
    master MPI process.

 - AllReduce

    Mixed mode allreduce benchmark. All threads under every MPI process
    combines its local buffer. All MPI processes then perform and
    MPI_Allreduce on these values to give the overall reduction value at
    each process.

 - Broadcast

    Data is broadcast from the master MPI process to all other MPI
    processes and then copied by each OpenMP thread. 

 - Scatter

    Mixed mode scatter benchmark. Master MPI process scatters its buffer
    to all other MPI processes. Each OpenMP thread then reads its portion
    of the buffer received.

 - Gather

    Mixed mode gather benchmark. All threads under an MPI process writes
    to a specific portion of a buffer. The master MPI process then gathers
    all the data using an MPI_Gather.

 - AlltoAll

    Mixed mode all to all benchmark. Each OpenMP thread sends/receives a
    portion of data to/from every other thread.


## Config File Notes

 - PingPong and PingPing tests include args indicating two MPI ranks to
   participate in exchange.  The intent is to allow for intra-node and
   inter-node rank selection (i.e., on same node or different nodes).

   Note, they use MPI_GET_PROCESSOR_NAME to differentiate ranks on 
   same/different nodes inside the benchmark code itself.

 - HaloExchange - all MPI procs participate (arranged in a ring, exchanging
   messages with two neighboring procs. 

 - Specify Min (integer), Max (integer), Time (seconds) at top of input file,
   first three lines.
        o Minimum data size
        o Maximum data size
        o Target time to roughly bound execution time for each test

    "The top three lines, which specify the minimum data size, maximum data 
     size and target time for each benchmark are required (and in that order).   
     The data size starts at the minimum data size and is doubled until the 
     maximum data size is reached. The target time parameter is used to keep 
     the execution time for each test approximately constant. It is used as 
     follows: for a given data size the benchmark is run for a certain number 
     of iterations. If the execution time is less than the target time the 
     test is re-run with twice the number of iterations until the target 
     time is met."

## Output 

 - Written to stdout
 - The Amount of data involved in the benchmark (labelled `Msg Size`)
 - The Total number of repetitions needed to execute for at least the target
   time (`No. Reps`) 
 - The Exact time for this number of repetitions (`Time (sec)`); 
 - The calculated time per repetition (` Time/Rep (s) `); 
 - The outcome of a correctness test carried on the microbenchmark (`Test`) 

## Misc

 - There are parse & plot python scripts included to process benchmark
   output.  Note, I had problem with plotting script so may just re-write
   as bash script (would avoid need for python-gnuplot to talk to gnuplot).

 - I had to hack a few bits in C code to get things to run with OpenMPI
   because the OMP parallel regions are marked as `shared(none)`.
   I had to add things like `shared(ompi_mpi_comm_world,ompi_mpi_int)`
   to parallel section specs.
   
 - My changes are currently as separate patch files in `scripts/patches/`
   [https://code-int.ornl.gov/3t4/t3ompix/tree/master/apps/epcc/scripts/patches]

 - I wrote a helper script to run the benchmark, `RUN-EPCC.sh` that takes
   three args: `-t NUM_THREADS`, `-r NUM_RANKS`, `-f INPUT_FILE`.
   It defaults to using localhost, but you can edit setting at top of
   script to add a hostfile to run over multiple nodes, set MPI binding,
   etc.

 - I had problems running the *MultiPingPong and *MultiPingPing tests.
   (See `scripts/inputs/FIXME_1` and `scripts/inputs/TODO.md`)
   Worked with `-t 2 -r 2 -f FIXME_1`.
   Failed with `-t 4 -r 2 -f FIXME_1`.


## References
 [1] "A Microbenchmark Suite for Mixed-Mode OpenMP/MPI",
     J. M. Bull, J. P. Enright, N. Ameer,  IWOMP 2009,
     Springer. 
     [http://rdcu.be/tbW6].

 [2] "Intro: OpenMPI/MPI Mixed-Mode Microbenchmarks",
     [http://www2.epcc.ed.ac.uk/~markb/mpiopenmpbench/intro.html]


## Related Links
 - EPCC OpenMP/MPI micro-benchmark webpage
   [https://www.epcc.ed.ac.uk/research/computing/performance-characterisation-and-benchmarking/epcc-openmpmpi-micro-benchmark]

 - My local git repo with source & scripts
   [https://code-int.ornl.gov/3t4/t3ompix/tree/master/apps/epcc]
 
