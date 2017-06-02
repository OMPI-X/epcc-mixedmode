EPCC Notes
----------


Download the EPCC hybrid benchmark (OpenMP + MPI)

    https://www.epcc.ed.ac.uk/research/computing/performance-characterisation-and-benchmarking/epcc-openmpmpi-micro-benchmark


If using OpenMPI, must hack the source files to add below to `#omp parallel` regions

    ```
       shared(ompi_mpi_comm_world,ompi_mpi_int)
    ```

The `RUN-EPCC.sh` wrapper script can be use to run the benchmark.  
To run EPCC by hand:

    ```
      export OMP_NUM_THREADS=4
      mpirun -np 2 ./mixedModeBenchmark inputfile
    ```

The `mixedModeFileParser.py` and `mixedModePlotter.py` Python scripts can be
used to generate plots with gnuplot for the results of the benchmark.  The
`mixedModePlotter.py` script requires the gnuplot-py package.  On
Ubuntu-16.04 this package can be installed like this:

    ```
      sudo apt-get install python-gnuplot gnuplot
    ```


References
----------
 - Main webpage
   https://www.epcc.ed.ac.uk/research/computing/performance-characterisation-and-benchmarking/epcc-openmpmpi-micro-benchmark

 - Intro to Mixed-Mode Microbenchmark webpage (saved locally as EPCC_intro.md)
   http://www2.epcc.ed.ac.uk/~markb/mpiopenmpbench/intro.html

 - Benchmark Summary paper
   https://doi.org/10.1007/978-3-642-02303-3_10
