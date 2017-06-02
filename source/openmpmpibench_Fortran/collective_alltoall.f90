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
! Implements the alltoall mixed mode OpenMP/MPI benchmark.  !   
!-----------------------------------------------------------!
MODULE collective_alltoall

  use parallelEnvironment
  use omp_lib
  use benchmarkSetup
  use output

  implicit none

  !Data arrays
  integer, dimension(:), allocatable :: alltoallSendBuf
  integer, dimension(:), allocatable :: alltoallRecvBuf
  integer, dimension(:), allocatable :: alltoallFinalBuf
  integer, dimension(:), allocatable :: testBuf

CONTAINS
  !---------------------------------------------------------!
  ! Subroutine: alltoall                                    !
  !                                                         !
  ! Driver routine for the alltoall benchmark.              !
  !---------------------------------------------------------!
  SUBROUTINE alltoall()
    integer :: dataSizeIter
    integer :: bufferSize

    !Initialise repsToDo to defaultReps
    repsToDo = defaultReps

    !Start loop over data sizes
    dataSizeIter = minDataSize !initialise dataSizeIter
    DO WHILE (dataSizeIter <= maxDataSize) 
       !Calculate buffer size and allocate space for 
       !the data arrays.
       bufferSize = dataSizeIter * (numThreads * numMPIprocs) &
            * numThreads
       
       CALL allocateData(bufferSize)
       
       !Perform warm-up of benchmark
       CALL alltoallKernel(warmUpIters,dataSizeIter)
       
       !Test if alltoall was successful
       CALL testAlltoall(dataSizeIter)

       !Initialise the benchmark
       benchComplete = .false.

       !Execute benchmark until target time is reached
       DO WHILE (benchComplete .NEQV. .true.) 
          !Start timer
          CALL MPI_Barrier(comm, ierr)
          startTime = MPI_Wtime()

          !Execute alltoall for repsToDo repetitions
          CALL alltoallKernel(repsToDo, dataSizeIter)

          !Stop timer
          CALL MPI_Barrier(comm, ierr)
          finishTime = MPI_Wtime()
          totalTime = finishTime - startTime

          !Test if target time was reached
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
       
       !Double data size and loop again
       dataSizeIter = dataSizeIter * 2

    END DO !End loop over data sizes

  END SUBROUTINE alltoall
  !---------------------------------------------------------!
  ! Subroutine: alltoallKernel                              !
  !                                                         !
  ! Implements the all to all benchmark.                    !
  ! Each thread sends/receives dataSize items to/from       !
  ! every other process.                                    !
  !---------------------------------------------------------!
  SUBROUTINE alltoallKernel(totalReps,dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i, j
    integer :: dataForEachProc, numsToWrite
    integer :: blockNum, startOffset

    !Calculate how much data each thread sends to each process
    numsToWrite = numThreads * dataSize
    !Calculate total amount of data each process gets 
    !from any other process...
    !..each thread gets dataSize items from every other thread.
    dataForEachProc = numThreads * numThreads * dataSize
    
    DO repIter = 1, totalReps
       
       !Each thread writes to numsToWrite items for each 
       !MPI process to alltoallSendBuf.
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(blockNum,i,j), &
!$OMP SHARED(numsToWrite,dataForEachProc,globalIDarray), &
!$OMP SHARED(alltoallSendBuf,numMPIprocs)
       !Calculate the blockNum of  a thread.
       !This is used to find which portion of the
       !dataForEachProc elements a thread will be responsible.
       blockNum = (myThreadID - 1)*numsToWrite

       !Write threadID to correct location in alltoallSendBuf
       DO i = 1, numMPIprocs !loop over MPI processes
          DO j = 1, numsToWrite !loop over data to write
             alltoallSendBuf(blockNum + ((i-1)*dataForEachProc + j)) = &
                  globalIDarray(myThreadID) 
          END DO
       END DO
!$OMP END PARALLEL
       
       !Call MPI_AlltoAll
       CALL MPI_Alltoall(alltoallSendBuf, dataForEachProc, MPI_INTEGER, &
            alltoallRecvBuf, dataForEachProc, MPI_INTEGER, &
            comm, ierr)

       !Each thread now reads the receive buffer so that 
       !it gets dataSize values from every other thread 
       !in its portion of alltoallFinalBuf
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(blockNum,startOffset,i,j),&
!$OMP SHARED(alltoallRecvBuf,alltoallFinalBuf,numMPIprocs),&
!$OMP SHARED(dataForEachProc,numsToWrite,dataSize,globalIDarray),&
!$OMP SHARED(numThreads)

       !Calculate the blockNum.
       !This which portion of the data from each process
       !a thread is responsible for.
       blockNum = (myThreadID-1)*dataSize
       
       !Calculate offset into the each MPI processes finalBuf where
       !each thread will start to write its data...
       !1) Calculate amount of data for each thread..
       startOffset = (numsToWrite * numMPIprocs)
       !2) Find which block in finalBuf for each thread
       startOffset = startOffset * (myThreadID -1)

       !Loop over all processors (threads & processes)
       DO i = 1, (numThreads*numMPIprocs)
          DO j = 1, dataSize
             alltoallFinalBuf(startOffset + ((i-1)*dataSize) + j ) = &
                  alltoallRecvBuf(blockNum + ((i-1)*numsToWrite) + j)
          END DO
       END DO

!$OMP END PARALLEL


    END DO !End loop over repetitions

  END SUBROUTINE alltoallKernel
  !---------------------------------------------------------!
  ! Subroutine: allocateData                                !
  !                                                         !
  ! Allocates memory for the main data arrays used in the   !
  ! alltoall benchmark.                                     !
  !---------------------------------------------------------!
  SUBROUTINE allocateData(bufferSize)
    integer, intent(in) :: bufferSize

    allocate(alltoallSendBuf(bufferSize))
    allocate(alltoallRecvBuf(bufferSize))
    allocate(alltoallFinalBuf(bufferSize))

  END SUBROUTINE allocateData

  !---------------------------------------------------------!
  ! Subroutine: freeData                                    !
  !                                                         !
  ! Free memory of the main data arrays.                    !
  !---------------------------------------------------------!
  SUBROUTINE freeData()
    
    deallocate(alltoallSendBuf)
    deallocate(alltoallRecvBuf)
    deallocate(alltoallFinalBuf)

  END SUBROUTINE freeData
  
  !---------------------------------------------------------!
  ! Subroutine: testAlltoall                                !
  !                                                         !
  ! Verifies that the all to all completed successfully.    !
  !---------------------------------------------------------!
  SUBROUTINE testAlltoall(dataSize)
    integer, intent(in) :: dataSize
    integer :: sizeofBuffer, i,j
    integer :: dataForEachThread, startElem
    logical :: testFlag, reduceFlag
    
    !Set testFlag to true
    testFlag = .true.

    !calculate the size of buffer on each process and allocte
    sizeofBuffer = dataSize * numThreads * numMPIprocs * numThreads
    allocate(testBuf(sizeofBuffer))

    !Calculate how many elements each thread will work with
    dataForEachThread = dataSize * numThreads * numMPIProcs

    !Fill buffer with expected values.
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,j,startElem), &
!$OMP SHARED(testBuf,globalIDarray,sizeofBuffer,dataSize),&
!$OMP SHARED(numThreads,numMPIprocs,dataForEachThread)

    !Calculate start element for each thread
    startElem = (myThreadID - 1)* dataForEachThread

    DO i = 1, (numThreads * numMPIprocs)
       DO j = 1, dataSize
          testBuf(startElem + (i-1)*dataSize + j) = i
       END DO
    END DO
    
!$OMP END PARALLEL
    
    !Compare
    DO i = 1, sizeofBuffer
       IF(alltoallFinalBuf(i) /= testBuf(i)) THEN
          testFlag = .false.
       END IF
    END DO

     !print *, myMPIRank,"Recvbuf =", alltoallFinalBuf
     
     !print *, myMPIRank,"testbuf =", testBuf

    !Reduce testFlag with logical AND operator to 
    !get overall test result.
    CALL MPI_Reduce(testFlag, reduceFlag, 1, MPI_LOGICAL, &
         MPI_LAND, 0, comm, ierr)

    !Master then sets testOutcome flag
    IF (myMPIRank == 0) THEN
       CALL setTestOutcome(reduceFlag)
    END IF
    
    !free space for testBuf
    deallocate(testBuf)

  END SUBROUTINE testAlltoall

END MODULE collective_alltoall
          
  
