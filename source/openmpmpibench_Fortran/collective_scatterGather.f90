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
! Implements the scatter and gather mixed mode OpenMP/MPI   !
! benchmarks.                                               !   
!-----------------------------------------------------------!
MODULE collective_scatterGather

  use parallelEnvironment
  use omp_lib
  use benchmarkSetup
  use output

  implicit none 
  !Data array declarations for scatter
  integer, dimension(:), allocatable :: scatterSendBuf
  integer, dimension(:), allocatable :: scatterRecvBuf

  !Data array declarations for gather
  integer, dimension(:), allocatable :: gatherSendBuf
  integer, dimension(:), allocatable :: gatherRecvBuf

  !Data arrays for final read and testing
  integer, dimension(:), allocatable :: finalBuf
  integer, dimension(:), allocatable :: testBuf

  integer, parameter :: SCATTERROOT = 0
  integer, parameter :: SCATTERSTARTVAL = 100

  integer, parameter :: GATHERROOT = 0
  integer, parameter :: GATHERSTARTVAL = 100


CONTAINS
  !---------------------------------------------------------!
  ! Subroutine: scatterGather                               !
  !                                                         !
  ! Driver routine for the scatter benchmark.               !
  !---------------------------------------------------------!
  SUBROUTINE scatterGather(benchmarkType)
    integer, intent(in) :: benchmarkType
    integer :: dataSizeIter
    integer :: bufferSize

    !Initialise repsToDo to defaultReps
    repsToDo = defaultReps
    
    !Start loop over data sizes
    dataSizeIter = minDataSize !initialise dataSizeIter
    DO WHILE (dataSizeIter <= maxDataSize)
       !Calculate buffer size and allocate space for 
       !scatter data arrays.
       bufferSize = dataSizeIter * numThreads

       
       IF (benchmarkType == SCATTER) THEN !Scatter
          CALL allocateData(bufferSize,benchmarkType)
          !Perform benchmark warm-up
          CALL scatterKernel(warmUpIters, dataSizeIter)
          !Test if scatter was successful
          CALL testScatterGather(bufferSize, benchmarkType)

       ELSE IF (benchmarkType == GATHER) THEN !Gather
          CALL allocateData(bufferSize,benchmarkType)
          !Perform benchmark warm-up
          CALL gatherKernel(warmUpIters, dataSizeIter)
          !Test if gather was successful
          IF (myMPIRank == GATHERROOT) THEN
             CALL testScatterGather(bufferSize*numMPIprocs, benchmarkType)
          END IF
       END IF
       
       !Initialise the benchmark
       benchComplete = .false.
       !Execute benchmark until target time is reached
       DO WHILE (benchComplete .NEQV. .true.)
          !Start timer
          CALL MPI_Barrier(comm, ierr)
          startTime = MPI_Wtime()

          IF (benchmarkType == SCATTER) THEN !Scatter
             !Execute scatter for repsToDo repetitions
             CALL scatterKernel(repsToDo, dataSizeIter)

          ELSE IF (benchmarkType == GATHER) THEN !Gather
             !Execute gather for repsToDo repetitions
             CALL gatherKernel(repsToDo, dataSizeIter)
          END IF
          
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
       CALL freeData(benchmarkType)

       !Double data size and loop again
       dataSizeIter = dataSizeIter * 2
       
    END DO !End loop over data sizes

  END SUBROUTINE scatterGather
  !---------------------------------------------------------!
  ! Subroutine: scatterKernel                               !
  !                                                         !
  ! Implement the scatter benchmark.                        !
  ! Root process first scatters send buffer to other        !
  ! processes.                                              !
  ! Each thread under a MPI process then reads its portion  !
  ! of scatterRecvBuf.                                      !
  !---------------------------------------------------------!
  SUBROUTINE scatterKernel(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    integer :: totalSendBufElems, sendCount, recvCount

    !Calculate totalSendBufElems
    totalSendBufElems = numMPIprocs * numThreads * dataSize
    !Calculate sendCount
    sendCount = dataSize * numThreads
    recvCount = sendCount

    DO repIter = 1, totalReps
       
       !Master process writes to scatterSendBuf
       IF (myMPIRank == SCATTERROOT) THEN
          DO i = 1, totalSendBufElems
             scatterSendBuf(i) = SCATTERSTARTVAL + i
          END DO
       END IF
       
       !Scatter the data to other processes
       CALL MPI_Scatter(scatterSendBuf, sendCount, MPI_INTEGER, &
            scatterRecvBuf, recvCount, MPI_INTEGER, &
            SCATTERROOT, comm, ierr)

       !Each thread now reads its portion of scatterRecvBuf
!$OMP PARALLEL DO DEFAULT(NONE),&
!$OMP PRIVATE(i),&
!$OMP SHARED(dataSize,recvCount,finalBuf,scatterRecvBuf),&
!$OMP SCHEDULE(STATIC,dataSize)

       DO i = 1, recvCount !looping over all data in recv buffer
          finalBuf(i) = scatterRecvBuf(i)
       END DO

!$OMP END PARALLEL DO

    END DO !End repetitions loop

  END SUBROUTINE scatterKernel

  !---------------------------------------------------------!
  ! Subroutine: gatherKernel                                !
  !                                                         !
  ! Implements the gather benchmark.                        !
  ! Each thread writes part of its buffer then all data     !
  ! is gathered to the master process.                      !
  !---------------------------------------------------------!
  SUBROUTINE gatherKernel(totalReps, dataSize)
    integer, intent(in) :: totalReps, dataSize
    integer :: repIter, i
    integer :: totalRecvBufElems
    integer :: sendCount, recvCount
    integer :: startVal

    !Calculate totalRecvBufElems
    totalRecvBufElems = dataSize * numThreads * numMPIprocs
    !Each process calculates its send and recv count
    sendCount = dataSize * numThreads
    recvCount = sendCount

    !Calculate startVal for each process
    !This is used to find the values to put in gatherSendBuf
    startVal = (myMPIRank * sendCount) + GATHERSTARTVAL

    DO repIter = 1, totalReps
       
       !Each thread writes to its portion of gatherSendBuf
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP PRIVATE(i), &
!$OMP SHARED(gatherSendBuf,startVal,dataSize,sendCount),&
!$OMP SCHEDULE(STATIC,dataSize)
       DO i = 1, sendCount
          gatherSendBuf(i) = startVal + i
       END DO
!$OMP END PARALLEL DO
       
       !Gather the data to GATHERROOT
       CALL MPI_Gather(gatherSendBuf, sendCount, MPI_INTEGER, &
            gatherRecvBuf, recvCount, MPI_INTEGER, &
            GATHERROOT, comm, ierr)
      
       !GATHERROOT process then copies its received data to
       !finalBuf
       IF (myMPIRank == GATHERROOT) THEN
          DO i = 1, totalRecvBufElems
             finalBuf(i) = gatherRecvBuf(i)
          END DO
       END IF

    END DO !End of repetitions loop

  END SUBROUTINE gatherKernel

  !---------------------------------------------------------!
  ! Subroutine: allocateData                                !
  !                                                         !
  ! Allocate memory for main data arrays                    !
  !---------------------------------------------------------!
  SUBROUTINE allocateData(bufferSize, benchmarkType)
    integer, intent(in) :: bufferSize, benchmarkType

    IF (benchmarkType == SCATTER) THEN !Allocate for scatter 

       !scatterSendBuf is size (bufferSize * numMPIprocs)
       IF (myMPIRank == SCATTERROOT) THEN
          allocate(scatterSendBuf(bufferSize*numMPIprocs))
       END IF
       allocate(scatterRecvBuf(bufferSize))
       allocate(finalBuf(bufferSize))

    ELSE IF (benchmarkType == GATHER) THEN !Allocate for gather

       allocate(gatherSendBuf(bufferSize))
       IF (myMPIRank == GATHERROOT) THEN
          allocate(gatherRecvBuf(bufferSize*numMPIprocs))
          allocate(finalBuf(bufferSize*numMPIprocs))
       END IF

    END IF
       
  END SUBROUTINE allocateData
          
  !---------------------------------------------------------!
  ! Subroutine: freeData                                    !
  !                                                         !
  ! Free memory of main data arrays.                        !
  !---------------------------------------------------------!
  SUBROUTINE freeData(benchmarkType)
    integer, intent(in) :: benchmarkType
    
    IF (benchmarkType == SCATTER) THEN
       
       IF (myMPIRank == SCATTERROOT) THEN
          deallocate(scatterSendBuf)
       END IF
       deallocate(scatterRecvBuf)
       deallocate(finalBuf)

    ELSE IF (benchmarkType == GATHER) THEN

       deallocate(gatherSendBuf)
       IF (myMPIRank == GATHERROOT) THEN
          deallocate(gatherRecvBuf) 
          deallocate(finalBuf)
       END IF

    END IF

  END SUBROUTINE freeData

  !---------------------------------------------------------!
  ! Subroutine: testScatterGather                           !
  !                                                         !
  ! Verifies that the scatter and gahter benchmarks worked  !
  ! correctly.                                              !
  !---------------------------------------------------------!
  SUBROUTINE testScatterGather(sizeofBuffer, benchmarkType)
    integer, intent(in) :: sizeofBuffer, benchmarkType
    integer :: i
    integer :: startVal
    logical :: testFlag, reduceFlag

    !initialise testFlag to true
    testFlag = .true.

    !Allocate space for testBuf
    allocate(testBuf(sizeofBuffer))
    
    IF (benchmarkType == SCATTER) THEN
       !Find the start scatter value for each MPI process
       startVal = (myMPIRank*sizeofBuffer) + SCATTERSTARTVAL
    ELSE IF (benchmarkType == GATHER) THEN
       !startVal is GATHERSTARTVAL
       startVal = GATHERSTARTVAL
    END IF

    !Fill testBuf with correct values
    DO i = 1, sizeofBuffer
       testBuf(i) = startVal + i
    END DO

    !Compare each element of finalBuf with testBuf 
    DO i = 1, sizeofBuffer
       IF (finalBuf(i) /= testBuf(i)) THEN
          testFlag = .false.
       END IF
    END DO

    !For scatter reduce testFlag into master with 
    !logical AND operator
    IF (benchmarkType == SCATTER) THEN
       CALL MPI_Reduce(testFlag, reduceFlag, 1, MPI_LOGICAL, &
            MPI_LAND, 0, comm, ierr)
       !Master then sets testOutcome using reduceFlag
       IF (myMPIRank == 0) THEN
          CALL setTestOutcome(reduceFlag)
       END IF
    ELSE IF (benchmarkType == GATHER) THEN
       CALL setTestOutcome(testFlag)
    END IF

    deallocate(testBuf)

  END SUBROUTINE testScatterGather

END MODULE collective_scatterGather
