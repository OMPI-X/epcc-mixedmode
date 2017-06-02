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
! Contains the point-to-point multi-pingping mixed mode     !
! OpenMP/MPI benchmarks.                                    !
! This includes: -masteronly multiPingping                  !
!                -funnelled multiPingping                   !
!                -multiple multiPingping                    !
!-----------------------------------------------------------!
MODULE pt_to_pt_multiPingPing
  
  use parallelEnvironment
  use omp_lib
  use benchmarkSetup
  use output

  implicit none

!Module variables
  integer :: pingNodeA, pingNodeB
  integer :: otherPingRank
  integer :: sizeofBuffer
  integer, dimension(:), allocatable :: pingSendBuf, pingRecvBuf
  integer, dimension(:), allocatable :: finalRecvBuf, testBuf

CONTAINS
!Module procedure specification

  !---------------------------------------------------------!
  ! Subroutine: multiPingPing                               !
  !                                                         !
  ! Driver subroutine for the multi-pingping benchmark.     !
  !---------------------------------------------------------!
  SUBROUTINE multiPingPing(benchmarkType)
    integer, intent(in) :: benchmarkType
    integer :: dataSizeIter
    character (len = MPI_MAX_PROCESSOR_NAME) :: otherProcName
    logical :: balance

    pingNodeA = 0
    pingNodeB = 1

    !Check if there's a balance in num of MPI processes 
    !in pingNodeA and pingNodeB.
    balance = crossCommBalance(pingNodeA, pingNodeB)
    !If not balanced..
    IF (balance .EQV. .false.) THEN
       !..Master prints error
       IF (myMPIRank == 0) THEN
          CALL printBalanceError()
       END IF
       !..and all processes return from subroutine.
       RETURN
    END IF

    !Exchange MPI_COMM_WORLD ranks for processes in same crossComm
    CALL exchangeWorldRanks(pingNodeA, pingNodeB, otherPingRank)

    !Processes on pingNodeB send processor name to pingNodeA procs
    CALL sendProcName(pingNodeA, pingNodeB, otherProcName)

    !Print comm world ranks & processor names of 
    !processes taking part in multi-pingping benchmark.
    CALL printMultiProcInfo(pingNodeA, otherPingRank, otherProcName)

    !Barrier to ensure that all procs have completed
    !printMultiProcInfo before printing column headings
    CALL MPI_Barrier(comm, ierr)
    !Master process then prints report column headings
    IF(myMPIRank == 0) THEN
       CALL printBenchHeader()
    END IF

    !initialise repsToDo to defaultReps at start of benchmark
    repsToDo = defaultReps

    !Start loop over data sizes
    dataSizeIter = minDataSize !initialise dataSizeIter
    DO WHILE (dataSizeIter <= maxDataSize)

       !set size of buffer
       sizeofBuffer = dataSizeIter * numThreads

       !Allocate space for the main data arrays
       CALL allocateData(sizeofBuffer)

       !Warm-up for benchmarkType
       IF (benchmarkType == MASTERONLY) THEN
          !Masteronly warm-up sweep
          CALL masteronlyMultiPingping(warmUpIters, dataSizeIter)
       ELSE IF (benchmarkType == FUNNELLED) THEN
          !Funnelled warm-up sweep
          CALL funnelledMultiPingping(warmUpIters, dataSizeIter)
       ELSE IF (benchmarkType == MULTIPLE) THEN
          !Multiple warm-up
          CALL multipleMultiPingping(warmUpIters, dataSizeIter)
       END IF
  
       !Verification test for multi-pingping
       CALL testMultiPingping(sizeofBuffer, dataSizeIter)

       !Initialise benchmark
       benchComplete = .false.
       !Keep executing benchmark until target time is reached
       DO WHILE (benchComplete .NEQV. .true.)

          !Start the timer..MPI_Barrier to synchronise
          !processes for more accurate timing.
          CALL MPI_Barrier(comm, ierr)
          startTime = MPI_Wtime()

          IF (benchmarkType == MASTERONLY) THEN
             !Execute masteronly multipingping repsToDo times
             CALL masteronlyMultiPingping(repsToDo, dataSizeIter)
          ELSE IF (benchmarkType == FUNNELLED) THEN
             !Execute funnelled multipingping
             CALL funnelledMultiPingping(repsToDo, dataSizeIter)
          ELSE IF (benchmarkType == MULTIPLE) THEN
             !Execute multiple multipingping
             CALL multipleMultiPingping(repsToDo, dataSizeIter)
          END IF

          !Stop the timer..MPI_Barrier to synchronise processes
          !for more accurate timing.
          CALL MPI_Barrier(comm, ierr)
          finishTime = MPI_Wtime()
          totalTime = finishTime - startTime

          !Call repTimeCheck function to test if target time
          !is reached.
          if (myMPIRank==0) then 
             benchComplete = repTimeCheck(totalTime, repsToDo)
          end if
          !Ensure all procs have the same value of benchComplete
          !and repsToDo
          call MPI_Bcast(benchComplete, 1, MPI_INTEGER, 0, comm, ierr)
          call MPI_Bcast(repsToDo, 1, MPI_INTEGER, 0, comm, ierr)

       END DO !end of loop to check if benchComplete is true

       !Master process sets benchmark results
       IF (myMPIRank == 0) THEN
          CALL setReportParams(dataSizeIter, repsToDo, totalTime)
          CALL printReport()
       END IF

       !Free the allocated space for the main data arrays
       CALL freeData()

       !Update dataSize before next iteration
       dataSizeIter = dataSizeIter * 2 !double data size

    END DO !end of loop over data sizes.

  END SUBROUTINE multiPingPing

  !---------------------------------------------------------!
  ! Subroutine: masteronlyMultiPingping                     !
  !                                                         !
  ! All Processes with rank of pingNodeA or pingNodeB in    !
  ! crossComm send a message to each other.                 !
  ! MPI communication takes place outside of the parallel   !
  ! region.                                                 !
  !---------------------------------------------------------!
  SUBROUTINE masteronlyMultiPingping(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    integer :: destRank

    !set destRank to ID of other process
    IF (crossCommRank == pingNodeA) THEN
       destRank = pingNodeB
    ELSE IF (crossCommRank == pingNodeB) THEN
       destRank = pingNodeA
    END IF

    DO repIter = 1, totalReps !loop totalRep times

       IF (crossCommRank == pingNodeA .or. &
            crossCommRank == pingNodeB) THEN
          
          !Each thread writes its globalID to pingSendBuf
          !with a PARALLEL DO directive.
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(pingSendBuf,dataSize,sizeofBuffer,globalIDarray), &
!$OMP SCHEDULE(STATIC,dataSize)

          DO i = 1, sizeofBuffer
             pingSendBuf(i) = globalIDarray(myThreadID)
          END DO

!$OMP END PARALLEL DO

          !Process calls non-blocking send to start transfer of
          !pingSendBuf to other process.
          CALL MPI_Isend(pingSendBuf, sizeofBuffer, MPI_INTEGER, &
               destRank, tag, crossComm, requestID, ierr)

          !Processes then waits for message from other process.
          CALL MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INTEGER, &
               destRank, tag, crossComm, status, ierr)

          !Finish the send operation with an MPI_Wait
          CALL MPI_Wait(requestID, status, ierr)

          !The threads under the MPI processes read their 
          !part of the received buffer.
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(finalRecvBuf,dataSize,sizeofBuffer,pingRecvBuf), &
!$OMP SCHEDULE(STATIC,dataSize)
          DO i = 1, sizeofBuffer
             finalRecvBuf(i) = pingRecvBuf(i)
          END DO
!$OMP END PARALLEL DO
          
       END IF
       
    END DO !End repetitions loop

  END SUBROUTINE masteronlyMultiPingping

  !---------------------------------------------------------!
  ! Subroutine: funnelledMultiPingping                      !
  !                                                         !
  ! All processes with rank of pingNodeA or pingNodeB in    !
  ! crossComm send a message to each other.                 !
  ! Inter-process communication takes place inside the      !
  ! OpenMP parallel region by the master thread.            !
  !---------------------------------------------------------!
  SUBROUTINE funnelledMultiPingping(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    integer :: destRank

    !set destRank to ID of other process
    IF (crossCommRank == pingNodeA) THEN
       destRank = pingNodeB
    ELSE IF (crossCommRank == pingNodeB) THEN
       destRank = pingNodeA
    END IF

    !Open the parallel region
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,repIter), &
!$OMP SHARED(dataSize,sizeofBuffer,pingSendBuf,globalIDarray), &
!$OMP SHARED(pingRecvBuf,finalRecvBuf,status,requestID,ierr), &
!$OMP SHARED(destRank,crossComm,crossCommRank,pingNodeA), &
!$OMP SHARED(pingNodeB,totalReps)

    DO repIter = 1, totalReps !loop totalRep times

       IF (crossCommRank == pingNodeA .or. &
            crossCommRank == pingNodeB) THEN
          
          !Each thread writes its globalID to its part 
          !of pingSendBuf with a OMP DO.
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1, sizeofBuffer
             pingSendBuf(i) = globalIDarray(myThreadID)
          END DO
!$OMP END DO
!Implicit barrier here takes care of necessary synchronisation.

!$OMP MASTER
          !Master thread of each process starts send.
          CALL MPI_Isend(pingSendBuf, sizeofBuffer, MPI_INTEGER, &
               destRank, tag, crossComm, requestID, ierr)

          !Processes then wait for message.
          CALL MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INTEGER, &
               destRank, tag, crossComm, status, ierr)

          !Finish the send operation with an MPI_Wait
          CALL MPI_Wait(requestID, status, ierr)
!$OMP END MASTER

!Barrier to ensure master thread has completed transfer.
!$OMP BARRIER

          !Each thread reads its part of the received buffer
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1, sizeofBuffer
             finalRecvBuf(i) = pingRecvBuf(i)
          END DO
!$OMP END DO

       END IF

    END DO !End repetitions loop
!$OMP END PARALLEL

  END SUBROUTINE funnelledMultiPingping

  !---------------------------------------------------------!
  ! Subroutine: multipleMultiPingping                       !
  !                                                         !
  ! All processes with crossCommRank of pingNodeA and       !
  ! pingNodeB in crossComm send a message to each other.    !
  ! Multiple threads take part in the communication.        !
  !---------------------------------------------------------!
  SUBROUTINE multipleMultiPingping(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    integer :: destRank
    integer :: lBound, uBound

    !set destRank to be ID of other process
    IF (crossCommRank == pingNodeA) THEN
       destRank = pingNodeB
    ELSE IF (crossCommRank == pingNodeB) THEN
       destRank = pingNodeA
    END IF

    !Open parallel region
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,repIter,lBound,uBound,requestID,status,ierr), &
!$OMP SHARED(dataSize,sizeofBuffer,pingSendBuf,globalIDarray), &
!$OMP SHARED(pingRecvBuf,finalRecvBuf,destRank,crossComm), &
!$OMP SHARED(crossCommRank,pingNodeA,pingNodeB,totalReps)

    DO repIter = 1, totalReps !loop totalRep times

       IF (crossCommRank == pingNodeA .or. &
            crossCommRank == pingNodeB) THEN
          
          !Calculate lower and upper bound of each threads
          !portion of the data arrays
          lBound = ((myThreadID-1) * dataSize) + 1
          uBound = (myThreadID * dataSize)

          !Each thread writes to its part of pingSendBuf
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1, sizeofBuffer
             pingSendBuf(i) = globalIDarray(myThreadID)
          END DO
!$OMP END DO NOWAIT

          !Each thread starts send of dataSize items from
          !pingSendBuf.
          CALL MPI_Isend(pingSendBuf(lBound:uBound), dataSize, &
               MPI_INTEGER, destRank, myThreadID, crossComm, &
               requestID, ierr)

          !Thread then waits for message from destRank
          !with tag equal to its threadID.
          CALL MPI_Recv(pingRecvBuf(lBound:uBound), dataSize, &
               MPI_INTEGER, destRank, myThreadID, crossComm, &
               status, ierr)

          !Thread completes send using MPI_Wait
          CALL MPI_Wait(requestID, status, ierr)

          !Each thread reads its part of received buffer.
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1, sizeofBuffer
             finalRecvBuf(i) = pingRecvBuf(i)
          END DO
!$OMP END DO NOWAIT

       END IF

    END DO !End reps
!$OMP END PARALLEL

  END SUBROUTINE multipleMultiPingping
          
  !---------------------------------------------------------!
  ! Subroutine: allocateData                                !
  !                                                         !
  ! Allocates space for the main data arrays.               !
  ! Size of each array is specified by subroutine argument. ! 
  !---------------------------------------------------------!
  SUBROUTINE allocateData(sizeofBuffer)
    integer, intent(in) :: sizeofBuffer

    IF (crossCommRank == pingNodeA .or. &
         crossCommRank == pingNodeB) THEN

       allocate(pingSendBuf(sizeofBuffer))
       allocate(pingRecvBuf(sizeofBuffer))
       allocate(finalRecvBuf(sizeofBuffer))

    END IF

  END SUBROUTINE allocateData

  !---------------------------------------------------------!
  ! Subroutine: freeData                                    !
  !                                                         !
  ! Free allocated memory for main data arrays.             !
  !---------------------------------------------------------!
  SUBROUTINE freeData()
    
    IF (crossCommRank == pingNodeA .or. &
         crossCommRank == pingNodeB) THEN
       
       deallocate(pingSendBuf, pingRecvBuf)
       deallocate(finalRecvBuf)

    END IF
 
  END SUBROUTINE freeData

  !---------------------------------------------------------!
  ! Subroutine: testMultiPingping                           !
  !                                                         !
  ! Verifies the the multi pingping benchmark worked        !
  ! correctly.                                              !
  !---------------------------------------------------------!
  SUBROUTINE testMultiPingping(sizeofBuffer, dataSize)
    integer, intent(in) :: sizeofBuffer, dataSize
    integer :: i
    logical :: testFlag, localTestFlag

    !set localTestFlag to true
    localTestFlag = .true.

    !Testing done for processes on pingNodeA & pingNodeB
    IF (crossCommRank == pingNodeA .or. &
         crossCommRank == pingNodeB) THEN

       !allocate space for testBuf
       allocate(testBuf(sizeofBuffer))

       !Construct testBuf with correct values
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(otherPingRank,numThreads,dataSize), &
!$OMP SHARED(sizeofBuffer,testBuf), &
!$OMP SCHEDULE(STATIC,dataSize)
       DO i = 1, sizeofBuffer
          !calculate globalID of thread expected in finalRecvBuf.
          !This is done using otherPingRank
          testBuf(i) = (otherPingRank * numThreads) + myThreadID
       END DO
!$OMP END PARALLEL DO

       !Compare each element of testBuf and finalRecvBuf
       DO i = 1, sizeofBuffer
          IF (testBuf(i) /= finalRecvBuf(i)) THEN
             localTestFlag = .false.
          END IF
       END DO

       !free space for testBuf
       deallocate(testBuf)

    END IF

    !Reduce testFlag into master with logical AND 
    CALL MPI_Reduce(localTestFlag, testFlag, 1, MPI_LOGICAL, &
         MPI_LAND, 0, comm, ierr)

    !master sets testOutcome flag
    IF (myMPIRank == 0) THEN
       CALL setTestOutcome(testFlag)
    END IF

  END SUBROUTINE testMultiPingping


END MODULE pt_to_pt_multiPingPing
