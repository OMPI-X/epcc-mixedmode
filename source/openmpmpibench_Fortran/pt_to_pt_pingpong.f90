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
! Contains the point-to-point pingpong mixed mode           !
! OpenMP/MPI benchmarks.                                    !
! This includes: -masteronly pingpong                       !
!                -funnelled pingpong                        !
!                -multiple pingpong                         !
!-----------------------------------------------------------!
MODULE pt_to_pt_pingPong
  
  use parallelEnvironment
  use omp_lib
  use benchmarkSetup
  use output

  implicit none

!Module variable specification 
  integer :: pingRank,pongRank
  integer :: sizeofBuffer
  integer, dimension(:), allocatable :: pingSendBuf, pingRecvBuf
  integer, dimension(:), allocatable :: pongSendBuf, pongRecvBuf
  integer, dimension(:), allocatable :: finalRecvBuf

  integer, dimension(:), allocatable :: testBuf
   
CONTAINS
!Module procedure specification

  !---------------------------------------------------------!
  ! Subroutine: pingPong                                    !
  !                                                         !
  ! Driver subroutine for the pingpong benchmark.           !
  !---------------------------------------------------------!
  SUBROUTINE pingPong(benchmarkType)
    integer, intent(in) :: benchmarkType
    integer :: dataSizeIter
    logical :: sameNode

    pingRank = PPRanks(1) 
    pongRank = PPRanks(2)
    
    !Check if pingRank and pongRank are on the same node
    sameNode = compareProcNames(pingRank,pongRank)
  
    !Master process then does reporting...
    IF (myMPIRank == 0) THEN
       !print message saying if benchmark is inter or intra node
       CALL printNodeReport(sameNode,pingRank,pongRank)
       !..then print report column headings.
       CALL printBenchHeader()
    END IF

    !initialise repsToDo to defaultReps at start of benchmark
    repsToDo = defaultReps

    !Loop over data sizes
    dataSizeIter = minDataSize !initialise dataSizeIter to minDataSize
    DO WHILE (dataSizeIter <= maxDataSize)
       
       !set sizeofBuffer 
       sizeofBuffer = dataSizeIter*numThreads

       !Allocate space for the main data arrays
       CALL allocateData(sizeofBuffer)

       !warm-up for either masteronly, funnelled or multiple
       IF (benchmarkType == MASTERONLY) THEN
          !Perform masteronly warm-up sweep
          CALL masteronlyPingpong(warmUpIters, dataSizeIter)
       ELSE IF (benchmarkType == FUNNELLED) THEN
          !Perform funnelled warm-up sweep
          CALL funnelledPingpong(warmUpIters, dataSizeIter)
       ELSE IF (benchmarkType == MULTIPLE) THEN
          !Perform multiple pingpong warmup
          CALL multiplePingpong(warmUpIters, dataSizeIter)
       END IF

       !Perform verification test for the pingpong
       CALL testPingpong(sizeofBuffer, dataSizeIter)
       
       !Initialise benchmark
       benchComplete = .false.
       !Keep executing benchmark until target time is reached
       DO WHILE (benchComplete .NEQV. .true.)

          !Start the timer...MPI_Barrier to synchronise 
          !processes for more accurate timing.
          CALL MPI_Barrier(comm, ierr)
          startTime = MPI_Wtime()
       
          IF (benchmarkType == MASTERONLY) THEN
             !Execute for repsToDo repetitions
             CALL masteronlyPingpong(repsToDo, dataSizeIter)
          ELSE IF (benchmarkType == FUNNELLED) THEN
             CALL funnelledPingpong(repsToDo, dataSizeIter)
          ELSE IF (benchmarkType == MULTIPLE) THEN
             CALL multiplePingpong(repsToDo, dataSizeIter)
          END IF

          !Stop the timer...MPI_Barrier to synchronise processes 
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
          CALL setReportParams(dataSizeIter,repsToDo,totalTime)
          CALL printReport()
       END IF

       !Free the allocated space for the main data arrays
       CALL freeData()

       !Update dataSize before next iteration
       dataSizeIter = dataSizeIter * 2 !double data size
       
    END DO !end of loop over data sizes
  
  END SUBROUTINE pingPong
  
  !---------------------------------------------------------!
  ! Subroutine: masteronlyPingpong                          !
  !                                                         !
  ! One MPI process sends single fixed length message to    !
  ! another MPI process.                                    !
  ! This other process then sends it back to the first      !
  ! process.                                                !
  !---------------------------------------------------------!
  SUBROUTINE masteronlyPingpong(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    DO repIter = 1, totalReps !loop totalRep times

       !All threads under MPI process with rank = pingRank
       !write to its part of the pingBuf array using a 
       !PARALLEL DO directive.
       IF (myMPIRank == pingRank) THEN
          
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(pingSendBuf,dataSize,sizeofBuffer,globalIDarray), &
!$OMP SCHEDULE(STATIC,dataSize)

          DO i = 1,sizeofBuffer
             pingSendBuf(i) = globalIDarray(myThreadID)
          ENDDO
!$OMP END PARALLEL DO
             
          !Ping process sends buffer to MPI process with rank equal
          !to pongRank.
          CALL MPI_Send(pingSendBuf, sizeofBuffer, MPI_INTEGER,&
               pongRank, tag, comm, ierr)
       
          !Process then waits for a message from pong process and
          !each thread reads its part of received buffer.
          !This completes the pingpong operation.
          CALL MPI_Recv(pongRecvBuf, sizeofBuffer, MPI_INTEGER,&
               pongRank, tag, comm, status, ierr)

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(pongRecvBuf,finalRecvBuf,dataSize,sizeofBuffer), &
!$OMP SCHEDULE(STATIC,dataSize)
          DO i = 1,sizeofBuffer
             finalRecvBuf(i) = pongRecvBuf(i)
          ENDDO
!$OMP END PARALLEL DO
          
       ELSEIF (myMPIRank == pongRank) THEN
          
          !pongRank receives the message from the ping process
          CALL MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INTEGER, &
               pingRank, tag, comm, status, ierr)

          !Each thread under the pongRank MPI process now copies 
          !its part of the received buffer to pongSendBuf

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(pongSendBuf,pingRecvBuf,dataSize,sizeofBuffer), &
!$OMP SCHEDULE(STATIC,dataSize)
          DO i = 1,sizeofBuffer
             pongSendBuf(i) = pingRecvBuf(i)
          END DO
!$OMP END PARALLEL DO
             
          !pongRank process now sends pongSendBuf to ping process.
          CALL MPI_Send(pongSendBuf, sizeofBuffer, MPI_INTEGER, &
               pingRank, tag, comm, ierr)
       END IF

    END DO !End repetitions loop
  
  END SUBROUTINE masteronlyPingpong
     
  !---------------------------------------------------------!
  ! Subroutine: funnelledPingpong                           !
  !                                                         !
  ! One MPI process sends single fixed length message to    !
  ! another MPI process.                                    !
  ! This other process then sends it back to the first      !
  ! process.                                                !
  ! All communication takes place within the OpenMP         !
  ! region for this benchmark.                              !
  !---------------------------------------------------------!
  SUBROUTINE funnelledPingpong(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i

    !Open parallel region for threads
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,repIter), &
!$OMP SHARED(pingRank,pongRank,pingSendBuf,pingRecvBuf), &
!$OMP SHARED(pongSendBuf,pongRecvBuf,finalRecvBuf,sizeofBuffer), &
!$OMP SHARED(dataSize,globalIDarray,comm,ierr,status), &
!$OMP SHARED(totalReps,myMPIRank)

    DO repIter = 1, totalReps !loop totalRep times

       !All threads under MPI process with rank = pingRank
       !write to its part of the pingBuf array using a 
       !PARALLEL DO directive.
       IF (myMPIRank == pingRank) THEN
       
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1,sizeofBuffer
             pingSendBuf(i) = globalIDarray(myThreadID)
          ENDDO
!$OMP END DO
!Implicit barrier at end of DO takes care of synchronisation here

          !Master thread under Ping process sends buffer 
          !to MPI process with rank equal to pongRank.
!$OMP MASTER
          CALL MPI_Send(pingSendBuf, sizeofBuffer, MPI_INTEGER,&
               pongRank, tag, comm, ierr)
       
          !Thread then waits for a message from pong process.
          CALL MPI_Recv(pongRecvBuf, sizeofBuffer, MPI_INTEGER,&
               pongRank, tag, comm, status, ierr)
!$OMP END MASTER

!Barrier needed to wait for master thread to complete MPI_Recv
!$OMP BARRIER 

          !Each thread reads its part of received buffer.
          !This completes the pingpong operation.
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1,sizeofBuffer
             finalRecvBuf(i) = pongRecvBuf(i)
          ENDDO
!$OMP END DO

       ELSEIF (myMPIRank == pongRank) THEN
          
          !Master thread under pongRank receives the message 
          !from the ping process
!$OMP MASTER
          CALL MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INTEGER, &
               pingRank, tag, comm, status, ierr)
!$OMP END MASTER

!Barrier needed to wait on master thread
!$OMP BARRIER

          !Each thread under the pongRank MPI process now copies 
          !its part of the received buffer to pongSendBuf

!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1,sizeofBuffer
             pongSendBuf(i) = pingRecvBuf(i)
          END DO
!$OMP END DO
!Implicit barrier at end of DO
   
          !Master thread pongRank process now sends pongSendBuf 
          !to ping process.
!$OMP MASTER
          CALL MPI_Send(pongSendBuf, sizeofBuffer, MPI_INTEGER, &
               pingRank, tag, comm, ierr)
!$OMP END MASTER

       END IF
    END DO !End repetitions loop
!$OMP END PARALLEL 
 
  END SUBROUTINE funnelledPingpong

  !---------------------------------------------------------!
  ! Subroutine: multiplePingpong                            !
  !                                                         !
  ! With this algorithmn multiple threads take place in the !
  ! communication and computation.                          !
  ! Each thread under the MPI ping process sends a portion  !
  ! of the message to the other MPI process.                !
  ! Each thread of the other process then sends it back to  !
  ! the first process.                                      !
  !---------------------------------------------------------!
  SUBROUTINE multiplePingpong(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    integer :: lBound, uBound

    !Open parallel region for threads under pingRank
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,lBound,uBound,repIter,ierr,status), &
!$OMP SHARED(myMPIRank,pingRank,pongRank,pingSendBuf,pingRecvBuf), &
!$OMP SHARED(pongSendBuf,pongRecvBuf,finalRecvBuf,sizeofBuffer),&
!$OMP SHARED(totalReps,dataSize,globalIDarray,comm)

    DO repIter = 1, totalReps !loop totalRep times

       IF (myMPIRank == pingRank) THEN
       
          !Calculate lower and upper bound of data array
          lBound = ((myThreadID-1)* dataSize) + 1
          uBound = (myThreadID * dataSize)
          
          !All threads under MPI process with rank = pingRank
          !write to its part of the pingBuf array using a 
          !PARALLEL DO directive
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1,sizeofBuffer
             pingSendBuf(i) = globalIDarray(myThreadID)
          ENDDO
!$OMP END DO NOWAIT
!Implicit barrier at end of DO not needed for multiple

          !Each thread under Ping process sends dataSize items
          !to MPI process with rank equal to pongRank.
          !myThreadID is used as tag to ensure data goes to correct 
          !place in buffer.
          CALL MPI_Send(pingSendBuf(lBound:uBound), dataSize,&
               MPI_INTEGER, pongRank, myThreadID, comm, ierr)
       
          !Thread then waits for a message from pong process.
          CALL MPI_Recv(pongRecvBuf(lBound:uBound), dataSize, &
               MPI_INTEGER, pongRank, myThreadID, comm, &
               status, ierr)

          !Each thread reads its part of received buffer.
          !This completes the pingpong operation.
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1,sizeofBuffer
             finalRecvBuf(i) = pongRecvBuf(i)
          ENDDO
!$OMP END DO NOWAIT


       ELSEIF (myMPIRank == pongRank) THEN

          !Calculate lower and upper bound of data array
          lBound = ((myThreadID-1)* dataSize)+1
          uBound = (myThreadID * dataSize)


          !Each thread under pongRank receives a message 
          !from the ping process
          CALL MPI_Recv(pingRecvBuf(lBound:uBound), dataSize, &
               MPI_INTEGER, pingRank, myThreadID, comm, &
               status, ierr)

          !Each thread under the pongRank MPI process now copies 
          !its part of the received buffer to pongSendBuf

!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1,sizeofBuffer
             pongSendBuf(i) = pingRecvBuf(i)
          END DO
!$OMP END DO NOWAIT
   
          !Each thread now sends pongSendBuf to ping process.
          CALL MPI_Send(pongSendBuf(lBound:uBound), dataSize, &
               MPI_INTEGER, pingRank, myThreadID, comm, ierr)

       END IF
    END DO !End repetitions loop
!$OMP END PARALLEL 
 
  END SUBROUTINE multiplePingpong

  !---------------------------------------------------------!
  ! Subroutine: allocateData                                !
  !                                                         !
  ! Allocates space for the main data arrays.               !
  ! Size of each array is specified by subroutine argument. ! 
  !---------------------------------------------------------!
  SUBROUTINE allocateData(sizeofBuffer)
    integer, intent(in) :: sizeofBuffer
    
    allocate(pingSendBuf(sizeofBuffer), pingRecvBuf(sizeofBuffer))
    allocate(pongSendBuf(sizeofBuffer), pongRecvBuf(sizeofBuffer))
    allocate(finalRecvBuf(sizeofBuffer))

  END SUBROUTINE allocateData

  !---------------------------------------------------------!
  ! Subroutine: freeData                                    !
  !                                                         !
  ! Deallocates the storage space for the main data arrays. ! 
  !---------------------------------------------------------!
  SUBROUTINE freeData()
    
    deallocate(pingSendBuf, pingRecvBuf)
    deallocate(pongSendBuf, pongRecvBuf)
    deallocate(finalRecvBuf)

  END SUBROUTINE freeData

  !---------------------------------------------------------!
  ! Subroutine: testPingPong                                !
  !                                                         !
  ! Verifies that the Ping Pong benchmark worked correctly. !
  !---------------------------------------------------------!
  SUBROUTINE testPingPong(sizeofBuffer, dataSize)
    integer, intent(in) :: sizeofBuffer, dataSize
    integer :: i
    logical :: testFlag

    !PingRank process checks if pingpong worked ok
    IF (myMPIRank == pingRank) THEN
       !initialise testFlag to true (test passed)
       testFlag = .true.

       !allocate space for the testBuf
       allocate(testBuf(sizeofBuffer))
    
       !Construct testBuf array with correct values.
       !These are the values that should be in finalRecvBuf.
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(testBuf,dataSize,sizeofBuffer,globalIDarray), &
!$OMP SCHEDULE(STATIC,dataSize)          

       DO i = 1,sizeofBuffer
          testBuf(i) = globalIDarray(myThreadID)
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
    
    END IF

    !pingRank broadcasts testFlag to the other processes
    CALL MPI_Bcast(testFlag, 1, MPI_LOGICAL, pingRank, comm, ierr)

    !Master process receives sets the testOutcome using testFlag.
    IF (myMPIRank == 0) THEN
       CALL setTestOutcome(testFlag)
    END IF
    
  END SUBROUTINE testPingPong

END MODULE pt_to_pt_pingPong
