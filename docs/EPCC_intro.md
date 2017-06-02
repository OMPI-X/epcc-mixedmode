
# OpenMP/MPI Mixed-Mode Microbenchmarks

[Source](http://www2.epcc.ed.ac.uk/~markb/mpiopenmpbench/intro.html "Permalink to OpenMP/MPI Mixed-Mode Microbenchmarks")

## About

With the current prevalence of multi-core processors in HPC architectures, mixed-mode programming, using both MPI and OpenMP in the same application is becoming increasingly important. We have designed and implemented a set of low-level microbenchmarks to test the performance of this programming model. 

The mixed-mode microbenchmarks provide analogues for the typical operations found in MPI microbenchmark suites, for both point-to-point and collective communication patterns.   
Two main considerations are captured in the microbenchmarks: 

* The effects of intra-thread communication and synchronisation by measuring both the times of MPI library calls and the reading and writing of buffers. 
* The ability to easily compare the effects of changing the OpenMP thread to MPI process ratio while using the same total number of cores. This is done by the appropriate choice of buffer sizes and message lengths. 

## Download

A Fortran and a C version of the mixed-mode microbenchmarks are available for download. 

  - Fortran version: [http://www2.epcc.ed.ac.uk/mixedMode_Fortran.tgz]
  - C version: [http://www2.epcc.ed.ac.uk/mixedMode_C.tgz]

## Compiling

#### Fortran version

1. Unpack the tar file with the command ` tar -xzf mixedMode_Fortran.tgz `
2. Edit the ` Makefile ` as follows: 
    - Set ` FC ` to the Fortran compiler you want to use. 
    - Set ` FFLAGS ` to any required compiler flags.   
    - These should include flags to process OpenMP directives and MPI calls and also some 
      standard optimisation flags (e.g. ` -O3 `) 
3. Type ` make ` to build the code. 

#### C version

1. Unpack the tar file with the command ` tar -xzf mixedMode_C.tgz `
2. Edit the ` Makefile ` as follows: 
    - Set ` CC ` to the C compiler you want to use. 
    - Set ` CFLAGS ` to any required compiler flags.   
    - These should include flags to process OpenMP directives and MPI calls and also some 
      standard optimisation flags (e.g. ` -O3 `) 
3. Type ` make ` to build the code. 

### Running the benchmark suite

An executable called ` mixedModeBenchmark ` is produced after successful compilation. 

To run the benchmark suite set the number of OpenMP threads via the ` OMP_NUM_THREADS ` 
environment variable and specify the number of MPI processes via the MPI job launcher(?).   
e.g. To run with 4 OpenMP threads and 2 MPI processes:   

    ```
     export OMP_NUM_THREADS=4
     mpirun -np 2 ./mixedModeBenchmark inputFile.txt
    ```
    
The executable takes an input file as an argument. This contains benchmark
setup information and a list of the benchmarks to be run. The layout of the
input file is described below. 

#### Input File Format

A sample input file is shown below: 

    ```
    1 # min data size 
    4194304 # max data size 
    1.00 # target time 
    masteronlyPingpong 
    0 1 
    funnelledPingping 
    0 -1 
    multipleMultiPingpong 
    funnelledHaloExchange 
    alltoall
    gather
    allReduce
    barrier
    scatter
    ```
    

The top three lines, which specify the minimum data size, maximum data size
and target time for each benchmark are required (and in that order).   The
data size starts at the minimum data size and is doubled until the maximum
data size is reached. The target time parameter is used to keep the
execution time for each test approximately constant. It is used as follows:
for a given data size the benchmark is run for a certain number of
iterations. If the execution time is less than the target time the test is
re-run with twice the number of iterations until the target time is met. 

Following this the list of benchmarks to run are specified.   For all
pingpong and pingping benchmarks the ranks of two MPI processes to
participate in the benchmark must be specified in the line following the
benchmark name. The numbers may be negative, in which case it is added to
the total number of MPI processes to get a valid MPI rank.   

A full list of all benchmarks available is given below: 

 1. `MasteronlyPingpong`  
    - Pingpong between two MPI processes with reads and writes by all OpenMP
      threads under each process. MPI communication takes place outside
      OpenMP parallel regions. 

 2. `FunnelledPingpong`  
    - Same as ` MasteronlyPingpong ` but MPI communication takes inside
      OpenMP parallel region by master thread. 

 3. `MultiplePingpong` 
    - Same as ` MasteronlyPingpong ` except all threads take part in MPI
      communication.

 4. `MasteronlyPingping`  
    - Two MPI processes send a message to each other with parallel reads and
      writes by all OpenMP threads under each process. Inter-process
      communication takes place outside of the parallel regions.  

 5. `FunnelledPingping`  
    - Same as ` MasteronlyPingping ` but inter-process communication takes
      place inside of the parallel regions by the master thread. 

 6. `MultiplePingping` 
    - Same as ` MultiplePingping ` but all threads take part in
      inter-process communication.

 7. `MasteronlyHaloexchange` 
    - All MPI processes are arranged in a ring and each process exchanges
      messages with its two neighbouring processes with parallel reads and
      writes the OpenMP threads under each process. Communication takes
      place outside of the parallel region. 

 8. `FunnelledHaloexchange` 
    - Same as above but communication is carried out inside the the OpenMP
      parallel regions by the master thread.

 9. `MultipleHaloexchange` 
    - Same as ` MasteronlyHaloexchange ` except inter-process communication
      is carried out inside the the OpenMP by all threads. 

 10. `MasteronlyMultipingpong` 
    - Every core on a node performs a pingpong with the corresponding core
      on another node. MPI messages are sent/received outside of the
      parallel regions.

 11. `FunnelledMultipingpong`  
    - Same as ` MasteronlyMultipingpong ` but MPI messages are carried out
      inside the parallel region by the master thread. 

 12. `MultipleMultipingpong` 
    - Same as ` MasteronlyMultipingpong ` with MPI messages carried out
      inside the parallel region by all threads. 

 13. `MasteronlyMultipingping` 
    - Every core on a node performs a pingping with the corresponding core
      on another node. MPI messages are sent/received outside of the
      parallel regions. 

 14. `FunnelledMultipingping`
    - Same as ` MasteronlyMultipingping ` but MPI messages are carried out
      inside the parallel region by the master thread.

 15. `MultipleMultipingping`  
    - Same as ` MasteronlyMultipingping ` with MPI messages carried out
      inside the parallel region by all threads. 

 16.  `Barrier` 
    - Mixed-mode barrier where threads under each process first synchronise
      with an ` OMP BARRIER ` followed by an ` MPI_Barrier ` to synchronise
      each MPI process.

 17. `Reduce`  
    - Mixed mode reduce benchmark. All threads under every MPI process
      combines its local buffer. All MPI processes then combined their
      values to get the overall reduction value which is stored on the
      master MPI process.

 18. `AllReduce`  
    - Mixed mode allreduce benchmark. All threads under every MPI process
      combines its local buffer. All MPI processes then perform and
      MPI_Allreduce on these values to give the overall reduction value at
      each process.

 19. `Broadcast`  
    - Data is broadcast from the master MPI process to all other MPI
      processes and then copied by each OpenMP thread. 

 20. `Scatter` 
    - Mixed mode scatter benchmark. Master MPI process scatters its buffer
      to all other MPI processes. Each OpenMP thread then reads its portion
      of the buffer received.

 21. `Gather`
    - Mixed mode gather benchmark. All threads under an MPI process writes
      to a specific portion of a buffer. The master MPI process then gathers
      all the data using an MPI_Gather.

 22. `AlltoAll`  
    - Mixed mode all to all benchmark. Each OpenMP thread sends/receives a
      portion of data to/from every other thread.

#### Output Format

The output from a run of the mixed-mode microbenchmarks is written to
stdout.   First general information about the run is given, the number of
MPI processes and OpenMP threads being used; the level of threading claimed
by the MPI implementation; the min and max data sizes; target time etc..
The output from each microbenchmark kernel follows. This includes: 

    - the amount of data involved in the benchmark (labelled ` Msg Size `); 
    - the total number of repetitions needed to execute for at least the target time (` No. Reps `) ; 
    - the exact time for this number of repetitions (` Time (sec) `); 
    - the calculated time per repetition (` Time/Rep (s) `); 
    - the outcome of a correctness test carried on the microbenchmark (` Test `) 
  
A sample output is shown below: 

`
    
    ---------------------------------------------- 
       Mixed mode MPI/OpenMP benchmark suite v1.0 
    ---------------------------------------------- 
      Number of MPI processes =           2 
      Number of OpenMP threads =           4 
      Thread support = MPI_THREAD_MULTIPLE 
     Reading parameters from input file.... 
    ------------------------------------------ 
               Benchmark parameters 
    ------------------------------------------ 
    Minimum data size                1
    Maximum data size          4194304
    Target time (sec)             1.00
    Default Repetitions           1000
    No. Warmup iterations            2
    --------------------------------------------
     # Masteronly Pingpong
    --------------------------------------------
     Inter node benchmark between process           0 and process           1
      Data Size     Msg Size (bytes)     No. Reps     Time (sec)     Time/Rep (s)     Test
     -----------   ------------------   ----------   ------------   --------------   ------
    d         1                   16         8000       2.083506       0.000260438     Pass
    d         2                   32         4000       1.142429       0.000285607     Pass
    d         4                   64         8000       2.124921       0.000265615     Pass
    d         8                  128         4000       1.161517       0.000290379     Pass
    d        16                  256         8000       2.144729       0.000268091     Pass
    --------------------------------------------
     # Funnelled Haloexchange
    --------------------------------------------
      Data Size     Msg Size (bytes)     No. Reps     Time (sec)     Time/Rep (s)     Test
     -----------   ------------------   ----------   ------------   --------------   ------
    d         1                   16         8000       2.116878       0.000264610     Pass
    d         2                   32         4000       1.183819       0.000295955     Pass
    d         4                   64         8000       2.158311       0.000269789     Pass
    d         8                  128         4000       1.219967       0.000304992     Pass
    d        16                  256         8000       2.188064       0.000273508     Pass
    --------------------------------------------
     # Scatter
    --------------------------------------------
      Data Size     Msg Size (bytes)     No. Reps     Time (sec)     Time/Rep (s)     Test
     -----------   ------------------   ----------   ------------   --------------   ------
    d         1                    4        16000       1.589212       0.000099326     Pass
    d         2                    8        16000       1.384335       0.000086521     Pass
    d         4                   16        16000       1.634958       0.000102185     Pass
    d         8                   32        16000       1.644211       0.000102763     Pass
    d        16                   64        16000       1.439673       0.000089980     Pass
    

`

Two Python programs to parse the output file and to generate plots using gnuplot are also available: 

  - `mixedModeFileParser.py` [http://www2.epcc.ed.ac.uk/mixedModeFileParser.py]
      takes the output from a run of the benchmark suite, extracts the
      message size and execution time fields and writes them to individual
      files for each benchmark. 

  - `mixedModePlotter.py` [http://www2.epcc.ed.ac.uk/mixedModePlotter.py]
      uses these files to generate plots for each set of benchmarks. For the
      barrier benchmark it plots execution time versus number of threads.
      For all other benchmarks the ratio of execution time to the execution
      time for 1 thread per process is plotted against message size. To use
      the plotter script you'll need to install the gnuplot.py package which
      can be found at [http://gnuplot-py.sourceforge.net]


Files
------
 - [http://www2.epcc.ed.ac.uk/mixedMode_Fortran.tgz]
 - [http://www2.epcc.ed.ac.uk/mixedMode_C.tgz]
 - [http://www2.epcc.ed.ac.uk/mixedModeFileParser.py]
 - [http://www2.epcc.ed.ac.uk/mixedModePlotter.py]
 - [http://gnuplot-py.sourceforge.net/]

  
