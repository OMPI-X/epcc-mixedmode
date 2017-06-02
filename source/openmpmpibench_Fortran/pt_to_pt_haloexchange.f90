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
! Contains the point-to-point halo exchange mixed mode      !
! OpenMP/MPI benchmarks.                                    !
! This includes: -masteronly haloexchange                   !
!                -funnelled haloexchange                    !
!                -multiple haloexchange                     !
!-----------------------------------------------------------!
MODULE pt_to_pt_haloExchange

  use parallelEnvironment
  use omp_lib
  use benchmarkSetup
  use output

  implicit none

!Module variables
  integer :: sizeofBuffer
  integer, dimension(:), allocatable :: leftSendBuf, leftRecvBuf
  integer, dimension(:), allocatable :: rightSendBuf, rightRecvBuf
  integer, dimension(:), allocatable :: finalLeftBuf, finalRightBuf
  integer, dimension(:), allocatable :: testLeftBuf, testRightBuf

CONTAINS
  !---------------------------------------------------------!
  ! Subroutine: haloExchange                                !
  !                                                         !
  ! Driver subroutine for the haloExchange benchmark.       !
  !---------------------------------------------------------!
  SUBROUTINE haloExchange(benchmarkType)
    integer, intent(in) :: benchmarkType
    integer :: dataSizeIter

    !find the ranks of the left and right neighbour
    CALL findNeighbourRanks()
    
    !initialise repsToDo to defaultReps
    repsToDo = defaultReps

    !Start loop over data sizes
    dataSizeIter = minDataSize !initialise dataSizeIter
    DO WHILE (dataSizeIter <= maxDataSize)
       !set sizeofBuffer
       sizeofBuffer = dataSizeIter * numThreads
       
       !Allocate space for the main data arrays
       CALL allocateData(sizeofBuffer)

       !perform benchmark warm-up 
       IF (benchmarkType == MASTERONLY) THEN
          CALL masteronlyHaloexchange(warmUpIters, dataSizeIter)
       ELSE IF (benchmarkType == FUNNELLED) THEN
          CALL funnelledHaloexchange(warmUpIters, dataSizeIter)
       ELSE IF (benchmarkType == MULTIPLE) THEN
          CALL multipleHaloexchange(warmUpIters, dataSizeIter)
       END IF

       !Each process performs a verification test 
       CALL testHaloexchange(sizeofBuffer, dataSizeIter)
       
       !Initialise the benchmark
       benchComplete = .false.
       !Execute benchmark until target time is reached
       DO WHILE (benchComplete .NEQV. .true.) 
          !Start timer
          CALL MPI_Barrier(comm, ierr)
          startTime = MPI_Wtime()
          
          !Execute benchmarkType for repsToDo repetitions
          IF (benchmarkType == MASTERONLY) THEN
             CALL masteronlyHaloexchange(repsToDo, dataSizeIter)
          ELSE IF (benchmarkType == FUNNELLED) THEN
             CALL funnelledHaloexchange(repsToDo, dataSizeIter)
          ELSE IF (benchmarkType == MULTIPLE) THEN
             CALL multipleHaloexchange(repsToDo, dataSizeIter)
          END IF

          !Stop timer.
          CALL MPI_Barrier(comm, ierr)
          finishTime = MPI_Wtime()
          totalTime = finishTime - startTime

          !Test if target time reached with the number of reps
          if (myMPIRank==0) then 
             benchComplete = repTimeCheck(totalTime, repsToDo)
          end if
          !Ensure all procs have the same value of benchComplete
          !and repsToDo
          call MPI_Bcast(benchComplete, 1, MPI_INTEGER, 0, comm, ierr)
          call MPI_Bcast(repsToDo, 1, MPI_INTEGER, 0, comm, ierr)
          
       END DO !end of benchcomplete loop

       !Master process sets benchmark results
       IF (myMPIRank == 0) THEN
          CALL setReportParams(dataSizeIter,repsToDo,totalTime)
          CALL printReport()
       END IF
       
       !Free allocated data
       CALL freeData()

       !Double dataSize and loop again
       dataSizeIter = dataSizeIter * 2

    END DO !End loop over data sizes

  END SUBROUTINE haloExchange

  !---------------------------------------------------------!
  ! Subroutine: masteronlyHaloexchange                      !
  !                                                         !
  ! Each process exchanges a message with its left and      !
  ! right neighbour.                                        !
  ! Communication takes place outside of the parallel       !
  ! region.                                                 !
  !---------------------------------------------------------!
  SUBROUTINE masteronlyHaloexchange(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    
    DO repIter = 1, totalReps !loop totalReps times
       
       !Each thread writes its globalID to rightSendBuf
       !and leftSendBuf with a parallel do directive
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(leftSendBuf,rightSendBuf,dataSize,sizeofBuffer), &
!$OMP SHARED(globalIDarray),&
!$OMP SCHEDULE(STATIC,dataSize)
       DO i = 1, sizeofBuffer
          leftSendBuf(i) = globalIDarray(myThreadID)
          rightSendBuf(i) = globalIDarray(myThreadID)
       END DO
!$OMP END PARALLEL DO

       !Process starts send of data to leftNeighbour and
       !rightNeighbour using non-blocking send...
       !..to leftNeighbour
       CALL MPI_ISend(leftSendBuf, sizeofBuffer, MPI_INTEGER, &
            leftNeighbour, tag, comm, requestArray(1), ierr)
       !..to rightNeighbour
       CALL MPI_ISend(rightSendBuf, sizeofBuffer, MPI_INTEGER, &
            rightNeighbour, tag, comm, requestArray(2), ierr)

       !Process then waits for messages from leftNeighbour and
       !rightNeighbour.
       !Receive leftRecvBuf from leftNeighbour
       CALL MPI_IRecv(leftRecvBuf, sizeofBuffer, MPI_INTEGER, &
            leftNeighbour, tag, comm, requestArray(3), ierr)
       !Receive rightRecvBuf from rightNeighbour
       CALL MPI_IRecv(rightRecvBuf, sizeofBuffer, MPI_INTEGER, &
            rightNeighbour, tag, comm, requestArray(4), ierr)

       !Finish the sends with an MPI_Waitall on the requests
       CALL MPI_Waitall(4, requestArray, statusArray, ierr) 

       !Each thread now reads its part of the left and right
       !received buffers.
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(leftRecvBuf,rightRecvBuf,dataSize,sizeofBuffer),&
!$OMP SHARED(finalLeftBuf,finalRightBuf), &
!$OMP SCHEDULE(STATIC,dataSize)
       DO i = 1, sizeofBuffer
          finalLeftBuf(i) = leftRecvBuf(i)
          finalRightBuf(i) = rightRecvBuf(i)
       END DO
!$OMP END PARALLEL DO

    END DO !End repetitions loop

  END SUBROUTINE masteronlyHaloexchange

  !---------------------------------------------------------!
  ! Subroutine: funnelledHaloexchange                       !
  !                                                         ! 
  ! Each process exchanges a message with its left and      !
  ! right neighbour.                                        !
  ! Communication takes place by one thread inside of the   !
  ! parallel region.                                        !
  !---------------------------------------------------------!
  SUBROUTINE funnelledHaloexchange(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    
     !Open the parallel region
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,repIter), &
!$OMP SHARED(dataSize,sizeofBuffer,leftSendBuf,rightSendBuf), &
!$OMP SHARED(rightRecvBuf,leftRecvBuf,finalLeftBuf), &
!$OMP SHARED(finalRightBuf,leftNeighbour,rightNeighbour), &
!$OMP SHARED(globalIDarray,ierr,comm,status), &
!$OMP SHARED(totalReps,requestArray,statusArray)

    DO repIter = 1, totalReps !loop for totalReps 

       !Each thread writes its globalID to rightSendBuf and
       !leftSendBuf.
!$OMP DO SCHEDULE(STATIC,dataSize)
       DO i = 1, sizeofBuffer
          leftSendBuf(i) = globalIDarray(myThreadID)
          rightSendBuf(i) = globalIDarray(myThreadID)
       END DO
!$OMP END DO
!Implicit barrier here takes care of necessary synchronisation
    
!$OMP MASTER
       !Master thread starts send of data to left and right
       !neighbours with a non-blocking send
       !..to leftNeighbour
       CALL MPI_ISend(leftSendBuf, sizeofBuffer, MPI_INTEGER, &
            leftNeighbour, tag, comm, requestArray(1), ierr)
       !..to rightNeighbour
       CALL MPI_ISend(rightSendBuf, sizeofBuffer, MPI_INTEGER, &
            rightNeighbour, tag, comm, requestArray(2), ierr)

       !Thread then starts receive of messages from leftNeighbour 
       !and rightNeighbour.
       !Receive leftRecvBuf from leftNeighbour
       CALL MPI_IRecv(leftRecvBuf, sizeofBuffer, MPI_INTEGER, &
            leftNeighbour, tag, comm, requestArray(3), ierr)
       !Receive rightRecvBuf from rightNeighbour
       CALL MPI_IRecv(rightRecvBuf, sizeofBuffer, MPI_INTEGER, &
            rightNeighbour, tag, comm, requestArray(4), ierr)
       
       !Finish the sends and receives with an MPI_Waitall 
       !on the requests.
       CALL MPI_Waitall(4, requestArray, statusArray, ierr) 
!$OMP END MASTER

!Barrier to ensure master thread has completed transfer.
!$OMP BARRIER
       
       !Each thread now reads its part of the left and right
       !received buffers.
!$OMP DO SCHEDULE(STATIC,dataSize)
       DO i = 1, sizeofBuffer
          finalLeftBuf(i) = leftRecvBuf(i)
          finalRightBuf(i) = rightRecvBuf(i)
       END DO
!$OMP END DO

    END DO !end of repetitions loop
!$OMP END PARALLEL
  END SUBROUTINE funnelledHaloexchange
  
  !---------------------------------------------------------!
  ! Subroutine: multipleHaloexchange                        !
  !                                                         ! 
  ! Each process exchanges a message with its left and      !
  ! right neighbour.                                        !
  ! All threads take part in the inter-porcess              !
  ! communication.                                          !
  !---------------------------------------------------------!
  SUBROUTINE multipleHaloexchange(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    integer :: lBound, uBound
    
     !Open the parallel region
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,requestArray,statusArray,status,ierr), &
!$OMP PRIVATE(lBound,uBound,repIter),&
!$OMP SHARED(dataSize,sizeofBuffer,leftSendBuf,rightSendBuf), &
!$OMP SHARED(rightRecvBuf,leftRecvBuf,finalLeftBuf), &
!$OMP SHARED(finalRightBuf,leftNeighbour,rightNeighbour), &
!$OMP SHARED(globalIDarray,comm,totalReps)

    DO repIter = 1, totalReps !repetition loop      

       !Calculate lower and upper bound of each threads
       !portion of the data arrays.
       lBound = ((myThreadID-1) * dataSize) + 1
       uBound = (myThreadID * dataSize)
       
       !Each thread writes its globalID to rightSendBuf and
       !leftSendBuf.
!$OMP DO SCHEDULE(STATIC,dataSize)
       DO i = 1, sizeofBuffer
          leftSendBuf(i) = globalIDarray(myThreadID)
          rightSendBuf(i) = globalIDarray(myThreadID)
       END DO
!$OMP END DO NOWAIT
!Implicit barrier not needed for multiple version

       !Each thread starts send of dataSize items to leftNeighbour
       !and to rightNeighbour
       CALL MPI_Isend(leftSendBuf(lBound:uBound), dataSize, &
            MPI_INTEGER, leftNeighbour, myThreadID, comm, &
            requestArray(1), ierr)
       CALL MPI_Isend(rightSendBuf(lBound:uBound), dataSize, &
            MPI_INTEGER, rightNeighbour, myThreadID, comm, &
             requestArray(2), ierr)
       
       !Each Thread then starts receive of messages from 
       !leftNeighbour and rightNeighbour.
       !Receive leftRecvBuf from leftNeighbour...
       CALL MPI_IRecv(leftRecvBuf(lBound:uBound), dataSize, &
            MPI_INTEGER, leftNeighbour, myThreadID, comm, &
            requestArray(3), ierr)
       !Receive rightRecvBuf from rightNeighbour
       CALL MPI_IRecv(rightRecvBuf(lBound:uBound), dataSize, &
            MPI_INTEGER, rightNeighbour, myThreadID, comm, &
            requestArray(4), ierr)
       
       !Finish the sends with an MPI_Waitall on the requests
       CALL MPI_Waitall(4, requestArray, statusArray, ierr)

       !Each thread now reads its part of the left and right
       !received buffers.
!$OMP DO SCHEDULE(STATIC,dataSize)
       DO i = 1, sizeofBuffer
          finalLeftBuf(i) = leftRecvBuf(i)
          finalRightBuf(i) = rightRecvBuf(i)
       END DO
!$OMP END DO NOWAIT

    END DO !End repetitions loop
!$OMP END PARALLEL

  END SUBROUTINE multipleHaloexchange

  !---------------------------------------------------------!
  ! Subroutine: allocateData                                !
  !                                                         !
  ! Allocate memory for the main data arrays in the         !
  ! haloexchange.                                           !
  !---------------------------------------------------------!
  SUBROUTINE allocateData(bufferSize)
    integer, intent(in) :: bufferSize

    allocate(leftSendBuf(bufferSize), leftRecvBuf(bufferSize))
    allocate(rightSendBuf(bufferSize), rightRecvBuf(bufferSize))
    allocate(finalLeftBuf(bufferSize), finalRightBuf(bufferSize))

  END SUBROUTINE allocateData

  !---------------------------------------------------------!
  ! Subroutine: freeData                                    !
  !                                                         !
  ! Free allocated memory for main data arrays.             !
  !---------------------------------------------------------!
  SUBROUTINE freeData()
    
    deallocate(leftSendBuf, leftRecvBuf)
    deallocate(rightSendBuf, rightRecvBuf)
    deallocate(finalLeftBuf, finalRightBuf)

  END SUBROUTINE freeData
    
  !---------------------------------------------------------!
  ! Subroutine: testHaloexchange                            !
  !                                                         !
  ! Verifies that the halo exchange benchmark worked        !
  ! correctly.                                              !
  !---------------------------------------------------------!
  SUBROUTINE testHaloexchange(sizeofBuffer, dataSize)
    integer, intent(in) :: sizeofBuffer, dataSize
    integer :: i
    logical :: testFlag, reduceFlag

    !set testFlag to true
    testFlag = .true.

    !allocate space for testLeftBuf and testRightBuf
    allocate(testLeftBuf(sizeofBuffer),testRightBuf(sizeofBuffer))

    !Construct testLeftBuf and testRightBuf with correct values
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(leftNeighbour,rightNeighbour,numThreads), &
!$OMP SHARED(dataSize,sizeofBuffer,testLeftBuf,testRightBuf),&
!$OMP SCHEDULE(STATIC,dataSize)
    DO i = 1, sizeofBuffer
       !Calculate globalID of thread expected in finalLeftBuf..
       testLeftBuf(i) = (leftNeighbour * numThreads) + myThreadID
       !...and in finalRightBuf.
       testRightBuf(i) = (rightNeighbour * numThreads) + myThreadID
    END DO
!OMP END PARALLEL DO

    !Compare..
    DO i = 1, sizeofBuffer
       !1) values from left neighbour
       IF (testLeftBuf(i) /= finalLeftBuf(i)) THEN
          testFlag = .false.
       END IF
       !2) values from right neighbour
       IF (testRightBuf(i) /= finalRightBuf(i)) THEN
          testFlag = .false.
       END IF
    END DO

    !Reduce testFlag into master with logical AND operator
    CALL MPI_Reduce(testFlag, reduceFlag, 1, MPI_LOGICAL, &
         MPI_LAND, 0, comm, ierr)

    !Master then sets testOutcome flag
    IF (myMPIRank == 0) THEN
       CALL setTestOutcome(reduceFlag)
    END IF

    !free space for testLeftBuf and testRightBuf
    deallocate(testLeftBuf, testRightBuf)

  END SUBROUTINE testHaloexchange

END MODULE pt_to_pt_haloExchange
