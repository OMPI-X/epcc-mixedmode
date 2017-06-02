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
! Contains the point-to-point multi-pingpong mixed mode     !
! OpenMP/MPI benchmarks.                                    !
! This includes: -masteronly multiPingpong                  !
!                -funnelled multiPingpong                   !
!                -multiple multiPingpong                    !
!-----------------------------------------------------------!
MODULE pt_to_pt_multiPingPong
  
  use parallelEnvironment
  use omp_lib
  use benchmarkSetup
  use output

  implicit none

!Module variable specification
  integer :: pingNode, pongNode
  integer :: sizeofBuffer
  integer, dimension(:), allocatable :: pingSendBuf, pingRecvBuf
  integer, dimension(:), allocatable :: pongSendBuf, pongRecvBuf
  integer, dimension(:), allocatable :: finalRecvBuf
  integer, dimension(:), allocatable :: testBuf

CONTAINS
!Module procedure specification

  !---------------------------------------------------------!
  ! Subroutine: multiPingPong                               !
  !                                                         !
  ! Driver subroutine for the multi-pingpong benchmark.     !
  !---------------------------------------------------------!
  SUBROUTINE multiPingPong(benchmarkType)
    integer, intent(in) :: benchmarkType
    integer :: dataSizeIter
    integer :: pongWorldRank
    character (len = MPI_MAX_PROCESSOR_NAME) :: pongProcName
    logical :: balance

    pingNode = 0
    pongNode = 1

    !Check if there's a balance in num of MPI processes on
    !pingNode and pongNode.
    balance = crossCommBalance(pingNode, pongNode)
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
    CALL exchangeWorldRanks(pingNode, pongNode, pongWorldRank)

    !Processes on pongNode send processor name to pingNode procs
    CALL sendProcName(pingNode, pongNode, pongProcName)

    !Print comm world ranks & processor names of
    !processes taking part in multi-pingpong benchmark.
    CALL printMultiProcInfo(pingNode, pongWorldRank, pongProcName)

    !Barrier to ensure that all procs have completed
    !printMultiProcInfo before printing column headings
    CALL MPI_Barrier(comm, ierr)
    !Master process then prints report column headings
    IF (myMPIRank == 0) THEN
       CALL printBenchHeader()
    END IF

    !initialise repsToDo to defaultReps at start of benchmark
    repsToDo = defaultReps

    !Loop over data sizes
    dataSizeIter = minDataSize !initialise dataSizeIter to minDataSize
    DO WHILE (dataSizeIter <= maxDataSize)
       
       !set sizeofBuffer
       sizeofBuffer = dataSizeIter * numThreads

       !Allocate space for the main data arrays
       CALL allocateData(sizeofBuffer)

       !Warm-up for either masteronly, funnelled or multiple
       IF (benchmarkType == MASTERONLY) THEN
          !Perform masteronly warm-up sweep
          CALL masteronlyMultiPingpong(warmUpIters, dataSizeIter)
       ELSE IF (benchmarkType == FUNNELLED) THEN
          !Funnelled warm-up sweep
          CALL funnelledMultiPingpong(warmUpIters, dataSizeIter)
       ELSE IF (benchmarkType == MULTIPLE) THEN
          !Mutiple pingpong warmup
          CALL multipleMultiPingpong(warmUpIters, dataSizeIter)
       END IF

       !Verification test for multi-pingpong
       CALL testMultiPingpong(sizeofBuffer, dataSizeIter)

       !Initialise benchmark
       benchComplete = .false.
       !Keep executing benchmark until target time is reached
       DO WHILE (benchComplete .NEQV. .true.)
          
          !Start the timer...MPI_Barrier to synchronise
          !processes for more accurate timing.
          CALL MPI_Barrier(comm, ierr)
          startTime = MPI_Wtime()
          
          IF (benchmarkType == MASTERONLY) THEN
             !Execute masteronly multipingpong repsToDo times
             CALL masteronlyMultiPingpong(repsToDo, dataSizeIter)
          ELSE IF (benchmarkType == FUNNELLED) THEN
             !Execute funnelled multipingpong
             CALL funnelledMultiPingpong(repsToDo, dataSizeIter)
          ELSE IF (benchmarkType == MULTIPLE) THEN
             !Exexute multiple multipingpong
             CALL multipleMultiPingpong(repsToDo, dataSizeIter)
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
          CALL setReportParams(dataSizeIter, repsToDo, totalTime)
          CALL printReport()
       END IF

       !Free the allocated space for the main data arrays
       CALL freeData()

       !Update dataSize before next iteration
       dataSizeIter = dataSizeIter * 2 !double data size

    END DO !end of loop over data sizes.

  END SUBROUTINE multiPingPong
    
  !---------------------------------------------------------!
  ! Subroutine: masteronlyMultiPingpong                     !
  !                                                         !
  ! All MPI processes in crossComm = pingNode sends a single!
  ! fixed length message to the neighbouring process in     !
  ! crossComm = pongNode.                                   !
  ! The neighbouring processes then sends the message back  !
  ! to the first process.                                   !
  !---------------------------------------------------------!
  SUBROUTINE masteronlyMultiPingpong(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i

    DO repIter = 1, totalReps !loop totalRep times

       !All threads under each MPI process with 
       !crossCommRank = pingNode write to pingSendBuf array
       !using a PARALLEL DO directive.
       IF (crossCommRank == pingNode) THEN
          
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(pingSendBuf,dataSize,sizeofBuffer,globalIDarray), &
!$OMP SCHEDULE(STATIC,dataSize)

          DO i = 1, sizeofBuffer
             pingSendBuf(i) = globalIDarray(myThreadID)
          END DO
!$OMP END PARALLEL DO

          !Each process with crossCommRank = pingNode sends
          !buffer to MPI process with rank = pongNode in crossComm.
          CALL MPI_Send(pingSendBuf, sizeofBuffer, MPI_INTEGER,&
               pongNode, tag, crossComm, ierr)

          !The processes then wait for a message from pong process
          !and each thread reads it part of the recieved buffer.
          CALL MPI_Recv(pongRecvBuf, sizeofBuffer, MPI_INTEGER,&
               pongNode, tag, crossComm, status, ierr)

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(pongRecvBuf,finalRecvBuf,dataSize,sizeofBuffer), &
!$OMP SCHEDULE(STATIC,dataSize)
          DO i = 1, sizeofBuffer
             finalRecvBuf(i) = pongRecvBuf(i)
          END DO
!$OMP END PARALLEL DO

       ELSEIF (crossCommRank == pongNode) THEN

          !Each process with crossCommRank = pongNode receives
          !the message from the pingNode processes.
          CALL MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INTEGER,&
               pingNode, tag, crossComm, status, ierr)

          !Each thread copies its part of the received buffer
          !to pongSendBuf.
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(pongSendBuf,pingRecvBuf,dataSize,sizeofBuffer),&
!$OMP SCHEDULE(STATIC,dataSize)
          DO i = 1,sizeofBuffer
             pongSendBuf(i) = pingRecvBuf(i)
          END DO
!$OMP END PARALLEL DO

          !The processes now send pongSendBuf to processes
          !with crossCommRank = pingNode.
          CALL MPI_Send(pongSendBuf, sizeofBuffer, MPI_INTEGER, &
               pingNode, tag, crossComm, ierr)
          
       END IF

    END DO !End repetitions loop

  END SUBROUTINE masteronlyMultiPingpong

  !---------------------------------------------------------!
  ! Subroutine: funnelledMultiPingpong                      !
  !                                                         !
  ! All MPI processes in crossComm = pingNode sends a single!
  ! fixed length message to the neighbouring process in     !
  ! crossComm = pongNode.                                   !
  ! The neighbouring processes then sends the message back  !
  ! to the first process.                                   !
  ! All communication takes place within the OpenMP parallel!
  ! region for this benchmark.                              !
  !---------------------------------------------------------!
  SUBROUTINE funnelledMultiPingpong(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i

    !Open parallel region for threads
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,repIter), &
!$OMP SHARED(pingNode,pongNode,pingSendBuf,pingRecvBuf),&
!$OMP SHARED(pongSendBuf,pongRecvBuf,finalRecvBuf,sizeofBuffer),&
!$OMP SHARED(dataSize,globalIDarray,crossComm,ierr,status), &
!$OMP SHARED(totalReps,myMPIRank,crossCommRank)

    DO repIter = 1, totalReps !loop totalRep times
       
       !All threads under each MPI process with
       !crossCommRank = pingNode write to pingSendBuf array
       !using a PARALLEL DO directive.
       IF (crossCommRank == pingNode) THEN

!$OMP DO SCHEDULE(STATIC,dataSize)

          DO i = 1, sizeofBuffer
             pingSendBuf(i) = globalIDarray(myThreadID)
          END DO

!$OMP END DO
!Implicit barrier at end of DO takes care of synchronisation here
          
          !Master thread under each pingNode process sends 
          !buffer to corresponding MPI process in pongNode
          !using the crossComm.
!$OMP MASTER
          CALL MPI_Send(pingSendBuf, sizeofBuffer, MPI_INTEGER,&
               pongNode, tag, crossComm, ierr)

          !The Master thread then waits for a message from the
          !pong process.
          CALL MPI_Recv(pongRecvBuf, sizeofBuffer, MPI_INTEGER,&
               pongNode, tag, crossComm, status, ierr)
!$OMP END MASTER

!Barrier needed to wait for master thread to complete MPI_Recv
!$OMP BARRIER

          !Each thread then reads its part of the received buffer.
!$OMP DO SCHEDULE(STATIC,dataSize)

          DO i = 1, sizeofBuffer
             finalRecvBuf(i) = pongRecvBuf(i)
          END DO

!$OMP END DO
       
       ELSE IF (crossCommRank == pongNode) THEN

          !Master thread under each pongNode process receives
          !the message from the pingNode processes.
!$OMP MASTER          
          CALL MPI_Recv(pingRecvBuf, sizeofBuffer, MPI_INTEGER,&
               pingNode, tag, crossComm, status, ierr)
!$OMP END MASTER

!Barrier needed to wait on master thread
!$OMP BARRIER

          !Each thread reads its part of the received buffer.
!$OMP DO SCHEDULE(STATIC,dataSize)

          DO i = 1, sizeofBuffer
             pongSendBuf(i) = pingRecvBuf(i)
          END DO

!$OMP END DO
!Implicit barrier at end of DO

          !Master threads sends pongSendBuf to processes
          !with crossCommRank = pingNode.
!$OMP MASTER
          CALL MPI_Send(pongSendBuf, sizeofBuffer, MPI_INTEGER,&
               pingNode, tag, crossComm, ierr)
!$OMP END MASTER

       END IF
    END DO !End repetitions loop

!$OMP END PARALLEL

  END SUBROUTINE funnelledMultiPingpong

  !---------------------------------------------------------!
  ! Subroutine: multipleMultiPingpong                       !
  !                                                         !
  ! Multiple threads take place in the communication and    !
  ! computation.                                            !
  ! Each thread of all MPI processes in crossComm = pingNode!
  ! sends a portion of the message to the neighbouring      !
  ! process in crossComm = pongNode.                        !
  ! Each thread of the neighbouring processes then sends    !
  ! the message back to the first process.                  !
  !---------------------------------------------------------!
  SUBROUTINE multipleMultiPingpong(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    integer :: lBound, uBound

    !Open parallel region for threads
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,repIter,ierr,status,lBound,uBound), &
!$OMP SHARED(pingNode,pongNode,pingSendBuf,pingRecvBuf),&
!$OMP SHARED(pongSendBuf,pongRecvBuf,finalRecvBuf,sizeofBuffer),&
!$OMP SHARED(dataSize,globalIDarray,crossComm), &
!$OMP SHARED(totalReps,myMPIRank,crossCommRank)

    DO repIter = 1, totalReps !loop totalReps time

       IF (crossCommRank == pingNode) THEN
          
          !Calculate lower and upper bound of data array
          lBound = ((myThreadID-1)* dataSize) + 1
          uBound = (myThreadID * dataSize)

          !All threads write to its part of the pingBuf
          !array using a PARALLEL DO directive
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1, sizeofBuffer
             pingSendBuf(i) = globalIDarray(myThreadID)
          END DO
!$OMP END DO NOWAIT
!Implicit barrier at end of DO not needed for multiple

          !Each thread under Ping process sends dataSize items
          !to pongNode process in crossComm.
          !myThreadID is used as tag to ensure data goes to 
          !correct place in buffer.
          CALL MPI_Send(pingSendBuf(lBound:uBound), dataSize,&
               MPI_INTEGER, pongNode, myThreadID, crossComm, ierr)

          !Thread then waits for a message from pongNode.
          CALL MPI_Recv(pongRecvBuf(lBound:uBound), dataSize,&
               MPI_INTEGER, pongNode, myThreadID, crossComm, &
               status, ierr)

          !Each thread reads its part of the received buffer,
!$OMP DO SCHEDULE(STATIC,dataSize)
          DO i = 1, sizeofBuffer
             finalRecvBuf(i) = pongRecvBuf(i)
          END DO
!$OMP END DO NOWAIT

          ELSEIF (crossCommRank == pongNode) THEN

             !Calculate lower and upper bound of data array
             lBound = ((myThreadID - 1) * dataSize) + 1
             uBound = (myThreadID * dataSize)

             !Each thread under pongRank receives a message
             !from the ping process
             CALL MPI_Recv(pingRecvBuf(lBound:uBound), dataSize, &
                  MPI_INTEGER, pingNode, myThreadID, crossComm, &
                  status, ierr)

             !Each thread now copies its part of the received
             !buffer to pongSendBuf
!$OMP DO SCHEDULE(STATIC,dataSize)
             DO i = 1, sizeofBuffer
                pongSendBuf(i) = pingRecvBuf(i)
             END DO
!$OMP END DO NOWAIT

                !Each thread now sends pongSendBuf to ping process.
                CALL MPI_Send(pongSendBuf(lBound:uBound), dataSize,&
                     MPI_INTEGER, pingNode, myThreadID, &
                     crossComm, ierr)

             END IF
          END DO !End repetitions loop
!$OMP END PARALLEL

        END SUBROUTINE multipleMultiPingpong

  !---------------------------------------------------------!
  ! Subroutine: allocateData                                !
  !                                                         !
  ! Allocates space for the main data arrays.               !
  ! Size of each array is specified by subroutine argument. ! 
  !---------------------------------------------------------
  SUBROUTINE allocateData(sizeofBuffer)
    integer, intent(in) :: sizeofBuffer
    
    IF (crossCommRank == pingNode) THEN
       !allocate space for arrays that MPI processes
       !with crossCommRank = pingNode will use
       allocate(pingSendBuf(sizeofBuffer))
       allocate(pongRecvBuf(sizeofBuffer))
       allocate(finalRecvBuf(sizeofBuffer))
    ELSE IF (crossCommRank == pongNode) THEN
       !allocate space for arrays that MPI processes
       !with crossCommRank = pongNode will use.
       allocate(pingRecvBuf(sizeofBuffer))
       allocate(pongSendBuf(sizeofBuffer))
    END IF
      
  END SUBROUTINE allocateData
          
  !---------------------------------------------------------!
  ! Subroutine: freeData                                    !
  !                                                         !
  ! Deallocates the storage space for the main data arrays. ! 
  !---------------------------------------------------------!
  SUBROUTINE freeData()
    
    IF (crossCommRank == pingNode) THEN
       deallocate(pingSendBuf)
       deallocate(pongRecvBuf)
       deallocate(finalRecvBuf)
    ELSE IF (crossCommRank == pongNode) THEN
       deallocate(pingRecvBuf)
       deallocate(pongSendBuf)
    END IF

  END SUBROUTINE freeData
        
  !---------------------------------------------------------!
  ! Subroutine: testMultiPingpong                           !
  !                                                         !
  ! Verifies the the multi pingpong benchmark worked        !
  ! correctly.                                              !
  !---------------------------------------------------------!
  SUBROUTINE testMultiPingpong(sizeofBuffer, dataSize)
    integer, intent(in) :: sizeofBuffer, dataSize
    integer :: i
    logical :: testFlag, localTestFlag

    !Initialise localtestFlag to true
    localTestFlag = .true.

    !All processes with crossCommRank = pingNode check
    !if multi-pingpong worked ok.
    IF (crossCommRank == pingNode) THEN
       
       !allocate space for testBuf
       allocate(testBuf(sizeofBuffer))
       
       !Construct testBuf array with correct values.
       !These are the values that should be in finalRecvBuf.
!$OMP PARALLEL DO DEFAULT(NONE),&
!$OMP PRIVATE(i), &
!$OMP SHARED(testBuf,dataSize,sizeofBuffer,globalIDarray),&
!$OMP SCHEDULE(STATIC,dataSize)

       DO i = 1,sizeofBuffer
          testBuf(i) = globalIDarray(myThreadID)
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

    !Reduce localTestFlag to master with logical AND operator
    CALL MPI_Reduce(localTestFlag, testFlag, 1, MPI_LOGICAL, &
         MPI_LAND, 0, comm, ierr)
    !Master then sets testOutcome using reduceFlag
    IF (myMPIRank == 0) THEN
       CALL setTestOutcome(testFlag)
    END IF

  END SUBROUTINE testMultiPingpong

END MODULE pt_to_pt_multiPingpong
          
          
