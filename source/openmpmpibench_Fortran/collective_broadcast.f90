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
! Implements the mixed mode OpenMP/MPI collective           !
!broadcast benchmark.                                       !
!-----------------------------------------------------------!
MODULE collective_broadcast

  use parallelEnvironment
  use omp_lib
  use benchmarkSetup
  use output

  implicit none 
  !Data array declarations
  integer, dimension(:), allocatable :: broadcastBuf
  integer, dimension(:), allocatable :: finalBroadcastBuf
  
  integer, parameter :: BROADCASTNUM = 100 !Number to broadcast
    
CONTAINS
  !---------------------------------------------------------!
  ! Subroutine: broadcast                                   !
  !                                                         !
  ! Driver subroutine for the broadcast benchmark.          !
  !---------------------------------------------------------!
  SUBROUTINE broadcast()
    integer :: dataSizeIter
    integer :: sizeofFinalBuf !needed for test

    !initialise repsToDo to defaultReps
    repsToDo = defaultReps

    !Start loop over data sizes
    dataSizeIter = minDataSize 
    DO WHILE (dataSizeIter <= maxDataSize) 
       !allocate space for main data arrays
       CALL allocateData(dataSizeIter)

       !Perform benchmark warm-up
       CALL broadcastKernel(warmUpIters,dataSizeIter)
       
       !Set sizeofFinalBuf and test if broadcast was a success
       sizeofFinalBuf = dataSizeIter * numThreads
       CALL testBroadcast(sizeofFinalBuf)

       !Initialise the benchmark 
       benchComplete = .false.
       !Execute benchmark until target time is reached
       DO WHILE (benchComplete .NEQV. .true.)
          !Start timer
          CALL MPI_Barrier(comm, ierr)
          startTime = MPI_WTime()

          !Execute broadcast for repsToDo repetitions
          CALL broadcastKernel(repsToDo, dataSizeIter)

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

       !Double dataSize and loop again
       dataSizeIter = dataSizeIter * 2
      
    END DO !End loop over data sizes

  END SUBROUTINE broadcast
  !---------------------------------------------------------!
  ! Subroutine: broadcastKernel                             !
  !                                                         !
  ! The broadcast benchmark.                                !
  ! At the start one process owns the data. After, all      !
  ! processes and threads have a copy of the data.          !
  !---------------------------------------------------------!
  SUBROUTINE broadcastKernel(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    !Set source of broadcast
    integer, parameter :: BROADCASTROOT = 0
    !Start position in finalBroadcastBuf of each thread.
    integer :: startPos 

    DO repIter = 1, totalReps

       !Master MPI process writes to broadcastBuf
       IF (myMPIRank == BROADCASTROOT) THEN
          DO i = 1, dataSize
             broadcastBuf(i) = BROADCASTNUM
          END DO
       END IF
    
       !Broadcast array to all other processes.
       CALL MPI_Bcast(broadcastBuf, dataSize, MPI_INTEGER, &
            BROADCASTROOT, comm, ierr)

       !Each thread copies broadcastBuf to its portion of
       !finalBroadcastBuf
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,startPos), &
!$OMP SHARED(dataSize,finalBroadcastBuf,broadcastBuf)
    
       !Calculate start of each threads portions of finalBroadcastBuf
       startPos = ((myThreadID-1) * dataSize)
       DO i = 1, dataSize
          finalBroadcastBuf(startPos + i) = broadcastBuf(i)
       END DO

!$OMP END PARALLEL
   
    END DO !End of repetitions loop
  
  END SUBROUTINE broadcastKernel
    
  !---------------------------------------------------------!
  ! Subroutine: allocateData                                !
  !                                                         !
  ! Allocate memory for the main data arrays in the         !
  ! broadcast operation.                                    !
  !---------------------------------------------------------!
  SUBROUTINE allocateData(bufferSize)
    integer, intent(in) :: bufferSize

    allocate(broadcastBuf(bufferSize))
    !finalBroadcastBuf is of size dataSize*numThreads
    allocate(finalBroadcastBuf(bufferSize*numThreads))
    
  END SUBROUTINE allocateData

  !---------------------------------------------------------!
  ! Subroutine: freeData                                    !
  !                                                         !
  ! Free memory of main data arrays.                        !
  !---------------------------------------------------------!
  SUBROUTINE freeData()
    
    deallocate(broadcastBuf)
    deallocate(finalBroadcastBuf)

  END SUBROUTINE freeData
  
  !---------------------------------------------------------!
  ! Subroutine: testBroadcast                               !
  !                                                         !
  ! Verifies that the broadcast benchmark worked correctly. !
  !---------------------------------------------------------!
  SUBROUTINE testBroadcast(bufferSize)
    integer, intent(in) :: bufferSize
    integer :: i
    logical :: testFlag, reduceFlag

    !Initialise testFlag to true
    testFlag = .true.

    !Compare each element of finalBroadcastBuf with BROADCASTNUM
    DO i = 1, bufferSize
       IF (finalBroadcastBuf(i) /= BROADCASTNUM) THEN
          testFlag = .false.
       END IF
    END DO

    !Reduce testFlag to master with logical AND operation
    CALL MPI_Reduce(testFlag, reduceFlag, 1, MPI_LOGICAL, &
         MPI_LAND, 0, comm, ierr)
    !Master then sets testOutcome using reduceFlag
    IF (myMPIRank == 0) THEN
       CALL setTestOutcome(reduceFlag)
    END IF

  END SUBROUTINE testBroadcast

END MODULE collective_broadcast
