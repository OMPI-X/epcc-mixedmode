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
! Implements the collective reduce and allreduce mixed      !
! mode OpenMP/MPI benchmarks.                               !
!-----------------------------------------------------------!
MODULE collective_reduction

  use parallelEnvironment
  use omp_lib
  use benchmarkSetup
  use output

  implicit none 
  !Data array declarations
  integer, dimension(:), allocatable :: globalReduceBuf
  integer, dimension(:), allocatable :: finalReduceBuf

CONTAINS
  !---------------------------------------------------------!
  ! Subroutine: reduction                                   !
  !                                                         !
  ! Driver subroutine for the reduce and allReduce          !
  ! benchmarks.                                             !
  !---------------------------------------------------------!
  SUBROUTINE reduction(benchmarkType)
    integer, intent(in) :: benchmarkType
    integer :: dataSizeIter
    integer :: sizeofBuf !for allReduce operation

    !initialise repsToDo to defaultReps
    repsToDo = defaultReps

    !Start loop over data sizes
    dataSizeIter = minDataSize !initialise dataSizeIter
    DO WHILE (dataSizeIter <= maxDataSize)
       !allocate space for the main data arrays..
       CALL allocateData(dataSizeIter)
       
       !Perform benchmark warm-up
       IF (benchmarkType == REDUCE) THEN
          CALL reduceKernel(warmUpIters,dataSizeIter)
          !Master process tests if reduce was success
          IF (myMPIRank == 0) THEN
             CALL testReduce(dataSizeIter,benchmarkType)
          END IF
       ELSE IF (benchmarkType == ALLREDUCE) THEN
          !calculate sizeofBuf for test
          sizeofBuf = dataSizeIter * numThreads
          CALL allReduceKernel(warmUpIters, dataSizeIter)
          !All processes need to perform unit test
          CALL testReduce(sizeofBuf,benchmarkType)
       END IF

       !Initialise the benchmark
       benchComplete = .false.
       !Execute benchmark until target time is reached
       DO WHILE (benchComplete .NEQV. .true.)
          !Start timer
          CALL MPI_Barrier(comm, ierr)
          startTime = MPI_WTime()
          
          !Execute reduce for repsToDo repetitions
          IF (benchmarkType == REDUCE) THEN
             CALL reduceKernel(repsToDo, dataSizeIter)
          ELSE IF (benchmarkType == ALLREDUCE) THEN
             CALL allReduceKernel(repsToDo, dataSizeIter)
          END IF
       
          !Stop timer
          CALL MPI_Barrier(comm, ierr)
          finishTime = MPI_WTime()
          totalTime = finishTime - startTime

          !Test if target time was reached with the number of reps
          if (myMPIRank==0) then 
             benchComplete = repTimeCheck(totalTime, repsToDo)
          end if
          !Ensure all procs have the same value of benchComplete
          !and repsToDo
          call MPI_Bcast(benchComplete, 1, MPI_INTEGER, 0, comm, ierr)
          call MPI_Bcast(repsToDo, 1, MPI_INTEGER, 0, comm, ierr)

       END DO
       
       !Master process sets benchmark result for reporting
       IF (myMPIRank == 0) THEN
          CALL setReportParams(dataSizeIter,repsToDo,totalTime)
          CALL printReport()
       END IF

       !Free allocated data
       CALL freeData()
       
       !Double dataSize and loop again
       dataSizeIter = dataSizeIter * 2
       
    END DO !End loop over data sizes

  END SUBROUTINE reduction

  !---------------------------------------------------------!
  ! Subroutine: reduce                                      !
  !                                                         !
  ! Implements the reduce mixed mode benchmark.             !
  ! Each thread under every MPI process combines its local  !
  ! buffer. This is then sent to the master MPI process to  !
  ! get the overall reduce value.                           !
  !---------------------------------------------------------!
  SUBROUTINE reduceKernel(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter,i
    !Decalre array which each thread reduces into
    integer, dimension(dataSize) :: localReduceBuf

    DO repIter = 1, totalReps !loop for totalReps

       !initialise all reduce arrays to ensure correct results
       localReduceBuf = 0
       globalReduceBuf = 0
       finalReduceBuf = 0

       !Open the parallel region and declare localReduceBuf
       !as a reduction variable.
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i),&
!$OMP SHARED(dataSize,globalIDarray),&
!$OMP REDUCTION(+:localReduceBuf)
       DO i = 1, dataSize
          localReduceBuf(i) = localReduceBuf(i) + globalIDarray(myThreadID)
       END DO
!$OMP END PARALLEL

       !Perform a reduce of localReduceBuf across the
       !MPI processes.
       CALL MPI_Reduce(localReduceBuf, globalReduceBuf, &
            dataSize, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
          
       !Copy globalReduceBuf into master Threads portion
       !of finalReduceBuf.
       ! FR this should only happen on rank==0 
       if (myMPIRank==0) then 
          finalReduceBuf(1:dataSize) = globalReduceBuf
       end if

    END DO !End repetitions loop

  END SUBROUTINE reduceKernel

  !---------------------------------------------------------!
  ! Subroutine: allReduce                                   !
  !                                                         !
  ! Implements the allreduce mixed mode benchmark.          !
  ! Each thread under every MPI process combines its local  !
  ! buffer. All MPI processes then combine this value to    !
  ! the overall reduction value at each process.            !
  !---------------------------------------------------------!
  SUBROUTINE allReduceKernel(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter,i
    integer :: startPos
    !Decalre array which each thread reduces into
    integer, dimension(dataSize) :: localReduceBuf

     DO repIter = 1, totalReps !loop for totalReps

       !initialise all reduce arrays to ensure correct results
       localReduceBuf = 0
       globalReduceBuf = 0
       finalReduceBuf = 0

       !Open the parallel region and declare localReduceBuf
       !as a reduction variable.
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i),&
!$OMP SHARED(dataSize,globalIDarray),&
!$OMP REDUCTION(+:localReduceBuf)
       DO i = 1, dataSize
          localReduceBuf(i) = localReduceBuf(i) + globalIDarray(myThreadID)
       END DO
!$OMP END PARALLEL

       !Perform an all reduce of localReduceBuf across 
       !the MPI processes.
       CALL MPI_Allreduce(localReduceBuf, globalReduceBuf, &
            dataSize, MPI_INTEGER, MPI_SUM, comm, ierr)
       
       !Each thread copies globalReduceBuf into its portion 
       !of finalReduceBuf
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,startPos), &
!$OMP SHARED(dataSize,finalReduceBuf,globalReduceBuf)

       !Calculate the start of each threads portion of finalReduceBuf
       startPos = ((myThreadID-1) * dataSize)
       DO i = 1, dataSize
          finalReduceBuf(startPos + i) = globalReduceBuf(i)
       END DO
!$OMP END PARALLEL

    END DO !End repetitions loop

  END SUBROUTINE allReduceKernel

  !---------------------------------------------------------!
  ! Subroutine: allocateData                                !
  !                                                         !
  ! Allocate memory for the main data arrays in the         !
  ! reduction operation.                                    !
  !---------------------------------------------------------!
  SUBROUTINE allocateData(bufferSize)
    integer, intent(in) :: bufferSize
    
    allocate(globalReduceBuf(bufferSize))
    !Final reduce is of size dataSize*numThreads
    allocate(finalReduceBuf(bufferSize*numThreads))

  END SUBROUTINE allocateData

  !---------------------------------------------------------!
  ! Subroutine: freeData                                    !
  !                                                         !
  ! Free allocated memory for main data arrays.             !
  !---------------------------------------------------------!
  SUBROUTINE freeData()
    
    deallocate(globalReduceBuf)
    deallocate(finalReduceBuf)
    
  END SUBROUTINE freeData

  !---------------------------------------------------------!
  ! Subroutine: testReduce                                  !
  !                                                         !
  ! Verifies that the reduction benchmarks worked correctly.!
  !---------------------------------------------------------!
  SUBROUTINE testReduce(bufferSize,benchmarkType)
    integer, intent(in) :: bufferSize, benchmarkType
    integer :: i
    integer :: correctReduce, lastGlobalID
    logical :: testFlag, reduceFlag
    
    !Initialise correctReduce to 0..
    correctReduce = 0
    !..and testFlag to true
    testFlag = .true.

    !set lastGlobalID
    lastGlobalID = (numMPIprocs * numThreads)

    !Now find correctReduce value by summing to lastGlobalID
    DO i = 1, lastGlobalID
       correctReduce = correctReduce + i
    END DO

    !Compare each element of finalRecvBuf to correctReduce 
    DO i = 1, bufferSize
       IF (finalReduceBuf(i) /= correctReduce) THEN
          testFlag = .false.
       END IF
    END DO

    !For allReduce, combine testFlag into master with logical AND
    IF (benchmarkType == ALLREDUCE) THEN
       CALL MPI_Reduce(testFlag, reduceFlag, 1, MPI_LOGICAL, &
            MPI_LAND, 0, comm, ierr)
       !then master sets testOutcome using reduceFlag
       IF (myMPIRank == 0) THEN
          CALL setTestOutcome(reduceFlag)
       END IF
    ELSE
       !For reduce master process just sets testOutcome using testFlag
       CALL setTestOutcome(testFlag)
    END IF

  END SUBROUTINE testReduce

END MODULE collective_reduction
