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
! Contains the point-to-point pingping mixed mode           !
! OpenMP/MPI benchmarks.                                    !
! This includes: -masteronly pingping                       !
!                -funnelled pingping                        !
!                -multiple pingping                         ! 
!-----------------------------------------------------------!
MODULE pt_to_pt_pingPing

  use parallelEnvironment
  use omp_lib
  use benchmarkSetup
  use output

  implicit none

!Module variables
  integer :: pingRankA, pingRankB !ranks of 2 pingping processes
  integer :: sizeofBuffer
  integer, dimension(:), allocatable :: pingSendBuf, pingRecvBuf
  integer, dimension(:), allocatable :: finalRecvBuf, testBuf
  
CONTAINS

  !---------------------------------------------------------!
  ! Subroutine: pingPing                                    !
  !                                                         !
  ! Driver subroutine for the pingping benchmark.           !
  !---------------------------------------------------------!
  SUBROUTINE pingPing(benchmarkType)
    integer, intent(in) :: benchmarkType
    integer :: dataSizeIter
    logical :: sameNode

    pingRankA = PPRanks(1)
    pingRankB = PPRanks(2)

    !Check if pingRankA and pingRankB are on the same node
    sameNode = compareProcNames(pingRankA, pingRankB)
    
    IF (myMPIRank == 0) THEN
       !print message saying if benchmark is inter or intra node
       CALL printNodeReport(sameNode,pingRankA,pingRankB)
       !..then print benchmark report column headings.
       CALL printBenchHeader()
    END IF

    !initialise repsToDo to defaultReps
    repsToDo = defaultReps

    !Start loop over data sizes
    dataSizeIter = minDataSize !initialise dataSizeIter
    DO WHILE (dataSizeIter <= maxDataSize)
 
       !set size of buffer
       sizeofBuffer = dataSizeIter * numThreads

       !Allocate space for main data arrays
       CALL allocateData(sizeofBuffer)

       !warm-up for benchmarkType
       IF (benchmarkType == MASTERONLY) THEN
          !Masteronly warmup sweep
          CALL masteronlyPingping(warmUpIters, dataSizeIter)
       ELSE IF (benchmarkType == FUNNELLED) THEN
          !funnelled warmup
          CALL funnelledPingping(warmUpIters, dataSizeIter)
       ELSE IF (benchmarkType == MULTIPLE) THEN
          !perform multiple pinping warmup
          CALL multiplePingping(warmUpIters, dataSizeIter)
       END IF

       !Perform verification test for pingping
       !this is only done by pingRankA and pingRankB
       CALL testPingping(sizeofBuffer, dataSizeIter)
       
       !Initialise the benchmark
       benchComplete = .false.
       !Execute benchmark until target time is reached
       DO WHILE (benchComplete .NEQV. .true.)
          !Start timer.
          CALL MPI_Barrier(comm, ierr)
          startTime = MPI_Wtime()
          
          !Execute benchmarkType for repsToDo repetitions
          IF (benchmarkType == MASTERONLY) THEN
             CALL masteronlyPingping(repsToDo, dataSizeIter)
          ELSE IF (benchmarkType == FUNNELLED) THEN
             CALL funnelledPingping(repsToDo, dataSizeIter)
          ELSE IF (benchmarkType == MULTIPLE) THEN
             CALL multiplePingping(repsToDo, dataSizeIter)
          END IF
          
          !Stop timer.
          CALL MPI_Barrier(comm, ierr)
          finishTime = MPI_Wtime()
          totalTime = finishTime - startTime

          !Test if target time was reached with number of
          !repetitions.
          if (myMPIRank==0) then 
             benchComplete = repTimeCheck(totalTime, repsToDo)
          end if
          !Ensure all procs have the same value of benchComplete
          !and repsToDo
          call MPI_Bcast(benchComplete, 1, MPI_INTEGER, 0, comm, ierr)
          call MPI_Bcast(repsToDo, 1, MPI_INTEGER, 0, comm, ierr)
          
       END DO !End of benchComplete loop

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

  END SUBROUTINE pingPing

  !---------------------------------------------------------!
  !Subroutine: masteronlyPingPing                           !
  !                                                         !
  ! Two processes send a message to each other using the    !
  ! MPI_Isend, MPI_Recv and MPI_Wait routines.              !
  ! Inter-process communication takes place outside of the  !
  ! parallel region.                                        !
  !---------------------------------------------------------!
  SUBROUTINE masteronlyPingping(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer ::  repIter, i
    integer :: destRank

    !set destRank to ID of other process
    IF (myMPIRank == pingRankA) THEN
       destRank = pingRankB
    ELSE IF (myMPIRank == pingRankB) THEN
       destRank = pingRankA
    END IF

    DO repIter = 1, totalReps !loop totalReps times

       IF(myMPIRank == pingRankA .or. myMPIRank == pingRankB) THEN

          !Each thread writes its globalID to pingSendBuf
          !using a PARALLEL DO directive.
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
          CALL MPI_ISend(pingSendBuf, sizeofBuffer, MPI_INTEGER, &
            destRank, tag, comm, requestID, ierr)

          !Process then waits for message from other process.
          CALL MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INTEGER, &
               destRank, tag, comm, status, ierr)

          !Finish the Send operation with an MPI_Wait
          CALL MPI_Wait(requestID, status, ierr)

          !Each thread under the MPI process now reads its part 
          !of the received buffer.
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
  END SUBROUTINE masteronlyPingping

  !---------------------------------------------------------!
  !Subroutine: funnelledPingPing                            !
  !                                                         !
  ! Two processes send a message to each other using the    !
  ! MPI_Isend, MPI_Recv and MPI_Wait routines.              !
  ! Inter-process communication takes place inside the      !
  ! OpenMP parallel region.                                 !
  !---------------------------------------------------------!
  SUBROUTINE funnelledPingping(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    integer :: destRank

    !set destRank to ID of other process
    IF (myMPIRank == pingRankA) THEN
       destRank = pingRankB
    ELSE IF (myMPIRank == pingRankB) THEN
       destRank = pingRankA
    END IF

    !Open the parallel region
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,repIter), &
!$OMP SHARED(dataSize,sizeofBuffer,pingSendBuf,globalIDarray), &
!$OMP SHARED(pingRecvBuf,finalRecvBuf,status,requestID,ierr), &
!$OMP SHARED(destRank,comm,myMPIRank,pingRankA,pingRankB,totalReps)

    DO repIter = 1, totalReps !loop totalRep times
       
       IF(myMPIRank == pingRankA .or. myMPIRank == pingRankB) THEN
          
          !Each threads write its globalID to its part of
          !pingSendBuf.
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1,sizeofBuffer
             pingSendBuf(i) = globalIDarray(myThreadID)
          END DO
!$OMP END DO
!Implicit barrier here takes care of necessary synchronisation.
          
!$OMP MASTER
          !Master thread starts send of buffer
          CALL MPI_Isend(pingSendBuf, sizeofBuffer, MPI_INTEGER,&
               destRank, tag, comm, requestID, ierr)
          !then waits for message from other process.
          CALL MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INTEGER,&
               destRank, tag, comm, status, ierr)
          !Master thread then completes send using MPI_Wait.
          CALL MPI_Wait(requestID, status, ierr)
!$OMP END MASTER

!Barrier needed to ensure master thread has completed transfer.
!$OMP BARRIER
          
          !Each thread reads its part of the received buffer
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1,sizeofBuffer
             finalRecvBuf(i) = pingRecvBuf(i)
          END DO
!$OMP END DO

       END IF

    END DO !End repetitions loop
!$OMP END PARALLEL

  END SUBROUTINE funnelledPingping
               
  !---------------------------------------------------------!
  ! Subroutine: multiplePingping                            !
  !                                                         !
  ! With this algorithmn multiple threads take part in the  !
  ! communication and computation.                          !
  ! Each thread sends its portion of the pingSendBuf to the !
  ! other process using MPI_Isend/ MPI_Recv/ MPI_Wait       !
  ! routines.                                               !
  !---------------------------------------------------------!
  SUBROUTINE multiplePingping(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    integer :: destRank 
    integer :: lBound, uBound

     !set destRank to be ID of other process
    IF (myMPIRank == pingRankA) THEN
       destRank = pingRankB
    ELSE IF (myMPIRank == pingRankB) THEN
       destRank = pingRankA
    END IF
       
     !Open parallel region
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,lBound,uBound,requestID,status,ierr,repIter), &
!$OMP SHARED(pingSendBuf,pingRecvBuf,finalRecvBuf,sizeofBuffer),&
!$OMP SHARED(destRank,myMPIRank,pingRankA,pingRankB,totalReps),&
!$OMP SHARED(dataSize,globalIDarray,comm)

    DO repIter = 1,totalReps !Loop totalReps times
       IF(myMPIRank == pingRankA .or. myMPIRank == pingRankB) THEN
          
           !Calculate lower and upper bound of each threads 
           !portion of the data arrays
           lBound = ((myThreadID-1)* dataSize) + 1
           uBound = (myThreadID * dataSize)
           
           !Each thread writes to its part of pingSendBuf
!$OMP DO SCHEDULE(STATIC,dataSize)
           DO i = 1,sizeofBuffer
              pingSendBuf(i) = globalIDarray(myThreadID)
           END DO
!$OMP END DO NOWAIT

           !Each thread starts send of dataSize items of 
           !pingSendBuf to process with rank = destRank
           CALL MPI_Isend(pingSendBuf(lBound:uBound), dataSize,&
                MPI_INTEGER, destRank, myThreadID, comm,&
                requestID, ierr)
           !Thread then waits for message from destRank with 
           !tag equal to its thread id.
           CALL MPI_Recv(pingRecvBuf(lBound:uBound), dataSize,&
                MPI_INTEGER, destRank, myThreadID, comm,&
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

     END DO !End repetitions loop
!$OMP END PARALLEL

   END SUBROUTINE multiplePingping
    
  !---------------------------------------------------------!
  ! Subroutine: allocateData                                !
  !                                                         !
  ! Allocate memory for the main data arrays in pingping.   !
  !---------------------------------------------------------!
  SUBROUTINE allocateData(bufferSize)
    integer, intent(in) :: bufferSize

    allocate(pingSendBuf(bufferSize), pingRecvBuf(bufferSize))
    allocate(finalRecvBuf(bufferSize))

  END SUBROUTINE allocateData
  
  !---------------------------------------------------------!
  ! Subroutine: freeData                                    !
  !                                                         !
  ! Free allocated memory for main data arrays.             !
  !---------------------------------------------------------!
  SUBROUTINE freeData()
    
    deallocate(pingSendBuf, pingRecvBuf)
    deallocate(finalRecvBuf)

  END SUBROUTINE freeData

  !---------------------------------------------------------!
  ! Subroutine: testPingPing                                !
  !                                                         !
  ! Verifies that the Ping Ping benchmark worked correctly. !
  !---------------------------------------------------------!
  SUBROUTINE testPingPing(sizeofBuffer, dataSize)
    integer, intent(in) :: sizeofBuffer, dataSize
    integer :: otherPingRank, i
    logical :: testFlag, reduceFlag
    
    !set testFlag to true
    testFlag = .true.
    
    !Testing only needs to be done by pingRankA & pingRankB
    IF (myMPIRank == pingRankA .or. myMPIRank == pingRankB) THEN
       !allocate space for testBuf
       allocate(testBuf(sizeofBuffer))
    
       !set the ID of other pingRank
       IF (myMPIRank == pingRankA) THEN
          otherPingRank = pingRankB
       ELSE IF (myMPIRank == pingRankB) THEN
          otherPingRank = pingRankA
       END IF
    
       !Construct testBuf with correct values
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(otherPingRank,numThreads,dataSize), &
!$OMP SHARED(sizeofBuffer,testBuf), &
!$OMP SCHEDULE(STATIC,dataSize)
       DO i = 1, sizeofBuffer
          !calculate globalID of thread expected in finalRecvBuf
          !this is done by using otherPingRank
          testBuf(i) = (otherPingRank * numThreads) + myThreadID
       END DO
!$OMP END PARALLEL DO

       !Compare each element of testBuf and finalRecvBuf
       DO i = 1, sizeofBuffer
          IF (testBuf(i) /= finalRecvBuf(i)) THEN
             testFlag = .false.
          END IF
       END DO

       !free space for testBuf
       deallocate(testBuf)
       
    END IF !End test loop
    
    !Reduce testFlag into master with logical AND
    CALL MPI_Reduce(testFlag, reduceFlag, 1, MPI_LOGICAL, &
         MPI_LAND, 0, comm, ierr)

    !master sets testOutcome flag
    IF (myMPIRank == 0) THEN
       CALL setTestOutcome(reduceFlag)
    END IF

  END SUBROUTINE testPingPing

END MODULE pt_to_pt_pingPing
