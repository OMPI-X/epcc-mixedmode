!*****************************************************************************
!*                                                                           *
!*             Mixed-mode OpenMP/MPI MicroBenchmark Suite - Version 1.0      *
!*                                                                           *
!*                            produced by                                    *
!*                                                                           *
!*                Mark Bull, Jim Enright and Fiona Reid                      *
!*                                                                           *
!*                                at                                         *
!*                                                                           *
!*                Edinburgh Parallel Computing Centre                        *
!*                                                                           *
!*   email: markb@epcc.ed.ac.uk, fiona@epcc.ed.ac.uk                         *
!*                                                                           *
!*                                                                           *
!*              Copyright 2012, The University of Edinburgh                  *
!*                                                                           *
!*                                                                           *
!*  Licensed under the Apache License, Version 2.0 (the "License");          *
!*  you may not use this file except in compliance with the License.         *
!*  You may obtain a copy of the License at                                  *
!*                                                                           *
!*      http://www.apache.org/licenses/LICENSE-2.0                           *
!*                                                                           *
!*  Unless required by applicable law or agreed to in writing, software      *
!*  distributed under the License is distributed on an "AS IS" BASIS,        *
!*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
!*  See the License for the specific language governing permissions and      *
!*  limitations under the License.                                           *
!*                                                                           *
!*****************************************************************************

!-----------------------------------------------------------!
! Main driver for mixed mode benchmark program.             !
! Reads benchmark input file.                               !
! Initialises the parallel environment.                     !
! Calls each benchmark.                                     !
!-----------------------------------------------------------!
PROGRAM mixedModeBenchmark
  use pt_to_pt_pingpong
  use pt_to_pt_pingping
  use pt_to_pt_multiPingPong
  use pt_to_pt_multiPingPing
  use pt_to_pt_haloExchange
  use collective_barrier
  use collective_reduction
  use collective_broadcast
  use collective_scatterGather
  use collective_alltoall
  use parallelEnvironment
  use benchmarkSetup
  use output

  implicit none
  
  !String for setting benchmark name for output
  character (len = MAXSTRING) :: name
  
  !Flag to check if benchmark is supported 
  logical :: supportFlag 

  !Initialise the parallel execution environment
  CALL initParallelEnv()
  
  !Master MPI process.....
  IF (myMPIRank == 0) THEN
     !1) Ptints header and parallel environment info.
     CALL printHeader(numMPIprocs,numThreads,threadSupport)
     !2) Opens the input file
     CALL openFile()
     !3) Setup the list of all possible benchmarks
     CALL setupBenchmarkList()
  END IF
  
  !Master reads parameters from input file and
  !broadcasts them to the other processes.
  CALL readBenchmarkParams()
  
  !Execute bencmarks by reading list from
  !input file.
  CALL findBenchmarkNumber()
  DO WHILE(benchmarkNumber /= FINISHED)
  
     benchmarks : SELECT CASE (benchmarkNumber)
     !Masteronly Pingpong
     CASE(1) 
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Masteronly Pingpong"
           CALL setBenchName(name, benchmarkNumber, supportFlag) 
        END IF
        !Execute benchmark
        CALL pingPong(MASTERONLY)
     
     !Funnelled Pingpong   
     CASE(2)
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Funnelled Pingpong"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
        END IF
        !Execute benchmark
        CALL pingPong(FUNNELLED)

     !Multiple Pingpong
     CASE(3) 
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_MULTIPLE)
           name = "Multiple Pingpong"
           CALL setBenchName(name, benchmarkNumber, supportFlag) 
        END IF
        !Execute benchmark 
        CALL pingPong(MULTIPLE)

     !Masteronly Pingping
     CASE(4)
        !Set name
         IF (myMPIRank == 0) THEN
            supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
            name = "Masteronly Pingping"
           CALL setBenchName(name, benchmarkNumber, supportFlag) 
        END IF
        !Execute benchmark
        CALL pingPing(MASTERONLY)

     !Funnelled Pingping
     CASE(5)
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Funnelled Pingping"
           CALL setBenchName(name, benchmarkNumber, supportFlag) 
        END IF
        !Execute benchmark
        CALL pingPing(FUNNELLED)

     !Multiple Pingping
     CASE(6)
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_MULTIPLE)
           name = "Multiple Pingping"
           CALL setBenchName(name, benchmarkNumber, supportFlag) 
        END IF
        !Execute benchmark
        CALL pingPing(MULTIPLE)

          !Masteronly Haloexchange
     CASE(7)
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Masteronly Haloexchange"
           CALL setBenchName(name, benchmarkNumber, supportFlag) 
           CALL printBenchHeader()
        END IF
        !Execute benchmark
        CALL haloExchange(MASTERONLY)

     !Funnelled Haloexchange    
     CASE(8) 
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Funnelled Haloexchange"
           CALL setBenchName(name, benchmarkNumber, supportFlag) 
           CALL printBenchHeader()
        END IF
        !Execute benchmark
        CALL haloExchange(FUNNELLED)

     !Multiple Haloexchange
     CASE(9) 
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_MULTIPLE)
           name = "Multiple Haloexchange"
           CALL setBenchName(name, benchmarkNumber, supportFlag) 
           CALL printBenchHeader()
        END IF
        !Execute benchmark
        CALL haloExchange(MULTIPLE)

     !Masteronly Multipingpong
     CASE(10)
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Masteronly MultiPingpong"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
        END IF
        !Execute benchmark
        CALL multiPingPong(MASTERONLY)

     !Funnelled Multipingpong
     CASE(11)
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Funnelled MultiPingpong"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
        END IF
        !Execute benchmark
        CALL multiPingPong(FUNNELLED)           

     !Multiple Multipingpong
     CASE(12)
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_MULTIPLE)
           name = "Multiple MultiPingpong"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
        END IF
        !Execute benchmark
        CALL multiPingPong(MULTIPLE)

     !Masteronly Multipingping
     CASE(13)
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Masteronly MultiPingping"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
        END IF
        !Execute benchmark
        CALL multiPingPing(MASTERONLY)

     !Funnelled Multipingping
     CASE(14)
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Funnelled MultiPingping"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
        END IF
        !Execute benchmark
        CALL multiPingPing(FUNNELLED)

     !Multiple Multipingping
     CASE(15)
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_MULTIPLE)
           name = "Multiple MultiPingping"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
        END IF
        !Execute benchmark
        CALL multiPingPing(MULTIPLE)

     !Barrier   
     CASE(16) 
        !Set name
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Barrier"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
           CALL printBenchHeader()
        END IF
        !Execute benchmark
        CALL barrierDriver()

     !Reduce   
     CASE(17)
        !Set name
        IF (myMPIRank == 0) THEN
            supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Reduce"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
           CALL printBenchHeader()
        END IF
        !Execute benchmark
        CALL reduction(REDUCE)

     !All-reduce   
     CASE(18)
        !Set name
        IF (myMPIRank == 0) THEN
            supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "All Reduce"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
           CALL printBenchHeader()
        END IF
        !Execute benchmark
        CALL reduction(ALLREDUCE)

     !Broadcast   
     CASE(19) 
        !Set name
        IF (myMPIRank == 0) THEN
            supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Broadcast"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
           CALL printBenchHeader()
        END IF
        !Execute benchmark
        CALL broadcast()
     
     !Scatter   
     CASE(20)
        IF (myMPIRank == 0) THEN
           supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Scatter"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
           CALL printBenchHeader()
        END IF
        !Execute banehmark
        CALL scatterGather(SCATTER)

     !Gather    
     CASE(21)
        IF (myMPIRank == 0) THEN
            supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "Gather"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
           CALL printBenchHeader()
        END IF
        !Execute benchmark
        CALL scatterGather(GATHER)

     !All to all
     CASE(22)
        !Set name
        IF (myMPIRank == 0) THEN
            supportFlag = benchmarkSupport(MPI_THREAD_FUNNELED)
           name = "All to all"
           CALL setBenchName(name, benchmarkNumber, supportFlag)
           CALL printBenchHeader()
        END IF
        !Execute benchmark
        CALL alltoall()
     
     !Default..file read error
     CASE default
        !..error message will already be printed out.

     END SELECT benchmarks
     
     !Read next benchmark from file
     CALL findBenchmarkNumber()
  END DO

  !Finalise programming environment
  CALL finaliseParallelEnv()
  !Master process closes file
  IF (myMPIRank == 0) THEN
     CALL closeFile()
  END IF
END PROGRAM mixedModeBenchmark

 
