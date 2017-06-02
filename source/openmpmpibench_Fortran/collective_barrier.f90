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
! Implements the collective barrier mixed mode OpenMP/MPI   !   
! benchmark.                                                !
!-----------------------------------------------------------!
MODULE collective_barrier

  use parallelEnvironment
  use omp_lib
  use benchmarkSetup
  use output

  implicit none 
  
CONTAINS
  !---------------------------------------------------------!
  ! Subroutine: barrierDriver                               !
  !                                                         !
  ! Driver subroutine for the barrier benchmark.            !
  !---------------------------------------------------------!
  SUBROUTINE barrierDriver()
    
    !initialise repsToDo to defaultReps
    repsToDo = defaultReps
    
    !Perform warm-up for benchmark
    CALL barrierKernel(warmUpIters)
    
    !Initialise the benchmark 
    benchComplete = .false.
    !Execute benchmark until target time is reached
    DO WHILE (benchComplete .NEQV. .true.)
       !Start timer
       CALL MPI_Barrier(comm,ierr)
       startTime = MPI_Wtime()
       
       !Execute benchmark for repsToDo repetitions
       CALL barrierKernel(repsToDo)

       !Stop timer 
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
       !No unit test, hardwire test result to pass
       CALL setTestOutcome(.true.)
       CALL setReportParams(1,repsToDo,totalTime)
       CALL printReport()
    END IF
  
  END SUBROUTINE barrierDriver

  !---------------------------------------------------------!
  ! Subroutine: barrierKernel                               !
  !                                                         !
  ! Main kernel for barrier benchmark.                      !
  ! First threads under each process synchronise with an    !
  ! OMP BARRIER. Then a MPI barrier synchronises each MPI   !
  ! process. MPI barrier is called within a OpenMP master   !
  ! directive.                                              !
  !---------------------------------------------------------!
  SUBROUTINE barrierKernel(totalReps)
    integer, intent(in) :: totalReps
    integer :: repIter

    !Open the parallel region
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(repIter), &
!$OMP SHARED(totalReps,comm,ierr) 

    DO repIter = 1, totalReps

       !Threads synchronise with an OpenMP barrier
!$OMP BARRIER

       !Master threads on each process now synhronise
!$OMP MASTER
       CALL MPI_Barrier(comm, ierr)
!$OMP END MASTER
       
    END DO !End repetitions loop

!$OMP END PARALLEL

  END SUBROUTINE barrierKernel

END MODULE collective_barrier
