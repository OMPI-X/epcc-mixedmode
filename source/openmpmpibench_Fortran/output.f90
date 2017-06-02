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
! Module with data and routines for benchmarks timing       !
! reports.                                                  !
!-----------------------------------------------------------!
MODULE output

  use benchmarkSetup
  use parallelEnvironment

  implicit none
 
  !Define data type for reporting
  type report
     character (len = MAXSTRING) :: benchName
     integer :: benchNumber
     logical :: supported
     integer :: dataSize
     integer :: bytes 
     logical :: testOutcome
     integer :: numReps
     DOUBLE PRECISION :: benchTime
     DOUBLE PRECISION :: timePerRep
  end type report

  !Declare benchReport, of type report
  type(report) :: benchReport

CONTAINS
  !---------------------------------------------------------!
  ! Subroutine: setBenchName                                !
  !                                                         !
  ! Sets the benchName, benchNumber and if the benchmark is !
  ! supported.                                              !
  !                                                         !
  !---------------------------------------------------------!
  SUBROUTINE setBenchName(name,number,support)
    character (len = MAXSTRING), intent (in) :: name
    integer, intent(in) :: number
    logical, intent(in) :: support

    benchReport%benchName = name
    benchReport%benchNumber = number
    benchReport%supported = support

    CALL printBenchName() 

  END SUBROUTINE setBenchName

  !---------------------------------------------------------!
  ! Subroutine: setTestOutcome                              !
  !                                                         !
  ! Sets benchReport's testOutcome element.                 !
  ! Called in test routine of each benchmark.               !
  !---------------------------------------------------------!
  SUBROUTINE setTestOutcome(outcome)
    logical, intent(in) :: outcome
    
    benchReport%testOutcome = outcome

  END SUBROUTINE setTestOutcome

  !---------------------------------------------------------!
  ! Subroutine: setReportParams                             !
  !                                                         !
  ! Sets the numReps and benchTime for a certain datasize   !  
  !---------------------------------------------------------!
  SUBROUTINE setReportParams(size,reps,time)
    integer, intent(in) :: size, reps
    DOUBLE PRECISION, intent(in) :: time
    
    benchReport%dataSize = size
    benchReport%numReps = reps
    benchReport%benchTime = time
    !Calculate and set time for 1 rep
    benchReport%timePerRep = time/reps
    !Calculate the size of message in bytes
    IF (benchReport%benchNumber <= LAST_PT_PT_ID) THEN
       !If point to point benchmark size of msg is 
       !dataSize * numThreads * sizeof(integer)
       benchReport%bytes = size * numThreads * sizeInteger
    ELSE IF (benchReport%benchNumber <= LASTMULTIPPID) THEN
       !If multi point to point benchmark size of msg is 
       !the size of the message leaving each node.
       benchReport%bytes = size * numThreads * sizeInteger * localCommSize
    ELSE
       benchReport%bytes = size * sizeInteger
    END IF

  END SUBROUTINE setReportParams
    
  !---------------------------------------------------------!
  ! Subroutine: printHeader                                 !
  !                                                         !
  ! Prints a header in the output.                          !
  !---------------------------------------------------------!
  SUBROUTINE printHeader(numProcs, numThreads, threadSupport)
    integer, intent(in) :: numProcs, numThreads, threadSupport
    character (len = MAXSTRING) :: string

    !Convert threadSupport to a string for output
    CALL threadSupportToString(threadSupport, string)

    write(*,*) "----------------------------------------------"
    write(*,*) "  Mixed mode MPI/OpenMP benchmark suite v1.0  "
    write(*,*) "----------------------------------------------"
    write(*,*) " Number of MPI processes =", numProcs
    write(*,*) " Number of OpenMP threads =", numThreads
    write(*,*) " Thread support = ", trim(string)
    
  END SUBROUTINE printHeader
  
  !---------------------------------------------------------!
  ! Subroutine: printBenchName                              !
  !                                                         !
  ! Print header for benchmark - name of benchmark and      !
  ! list of names of each column.                           !
  !---------------------------------------------------------!
  SUBROUTINE printBenchName()
    write(*,*) "--------------------------------------------"
    write(*,*) "# ", benchReport%benchName
    write(*,*) "--------------------------------------------"

    !print warning if benchmark not supported
    IF (benchReport%supported .EQV. .false.) THEN
       write(*,*) "WARNING: Implementation does not ",&
            "support benchmark"
    END IF
    
    !Flush output buffer
    !CALL flush(6)

  END SUBROUTINE printBenchName

  !----------------------------------------------------------!
  ! Subroutine: printBenchHeader                             !
  !                                                          !
  ! Prints the column headings for the benchmark report.     !
  !----------------------------------------------------------!
  SUBROUTINE printBenchHeader()


    write(*,fmt="(2x,a9,5x,a16,5x,a8,5x,a10,5x,a12,5x,a4)")&
         "Data Size","Msg Size (bytes)","No. Reps",&
         "Time (sec)","Time/Rep (s)","Test"
    write(*,fmt="(1x,a11,3x,a18,3x,a10,3x,a12,3x,a14,3x,a6)")&
         "-----------","------------------","----------",&
         "------------","--------------","------"


  END SUBROUTINE printBenchHeader

  !---------------------------------------------------------!
  ! Subroutine: printNodeReport                             !
  !                                                         !
  ! For pingpong and pingping benchmarks prints out if the  !
  ! two MPI processes are on the same node or not.          !
  !---------------------------------------------------------!
  SUBROUTINE printNodeReport(sameNode,rankA,rankB)
    integer, intent(in) :: rankA, rankB
    logical, intent(in) :: sameNode
    
    IF (sameNode .EQV. .true.) THEN
       write(*,*) "Intra node benchmark between process",&
            rankA, "and process", rankB 
    ELSE IF (sameNode .EQV. .false.) THEN
       write(*,*) "Inter node benchmark between process",&
            rankA, "and process", rankB
    END IF
    
  END SUBROUTINE printNodeReport

  !---------------------------------------------------------!
  ! Subroutine: printMultiProcInfo                          !
  !                                                         !
  ! This prints the comm world ranks and processor names for!
  ! each pair of processes in the multi-pingpong or         !
  ! multi-pingping benchmarks.                              !
  !---------------------------------------------------------!
  SUBROUTINE printMultiProcInfo(printNode, pairWorldRank, pairProcName)
    integer, intent(in) :: printNode, pairWorldRank
    character (len = MPI_MAX_PROCESSOR_NAME) :: pairProcName
    
    !MPI processes under printNode of crossComm perform the output
    IF (crossCommRank == printNode) THEN
       print *, "MPI process ", myMPIRank, "on ", trim(myProcName), &
            " commumicating with MPI process ", pairWorldRank, &
            "on ", trim(pairProcName)
    END IF
    
  END SUBROUTINE printMultiProcInfo

  !---------------------------------------------------------!
  ! Subroutine: printReport                                 !
  !                                                         !
  ! Prints out the a column of information after each       !
  ! data size iteration.                                    !
  !---------------------------------------------------------!
  SUBROUTINE printReport()
    character (len =4) testString
    
    !Set testString
    IF(benchReport%testOutcome .EQV. .true.) THEN
       testString = "Pass"
    ELSE
       testString = "Fail"
    END IF
    
    !Write to output
    write(*,fmt="('d',i10,5x,i16,5x,i8,5x,f10.6,4x,f14.9,5x,a4)")&
         benchReport%dataSize, benchReport%bytes,&
         benchReport%numReps,benchReport%benchTime,&
         benchReport%timePerRep,testString
    
  END SUBROUTINE printReport

  !--------------------------------------------------------!
  ! Subroutine: printBalanceError                          !
  !                                                        !
  ! Prints an error if there isn't the same number of MPI  !
  ! processes in the nodes selected for the multi-pingpong !
  ! or multi-pingping benchmarks.                          !
  !--------------------------------------------------------!
  SUBROUTINE printBalanceError()
    
    print *, ""
    print *, "ERROR: Nodes selected for this benchmark do not",&
         "have same number of MPI processes per node.", &
         "Skipping benchmark..."
    print *, ""

  END SUBROUTINE printBalanceError
  !--------------------------------------------------------!
  ! Subroutine: threadSupportToString                      !
  !                                                        !
  ! Converts the threadSupport integer variable to a       !
  ! string for output.                                     !
  !--------------------------------------------------------!
  SUBROUTINE threadSupportToString(threadSupport, string)
    integer, intent(in) :: threadSupport
    character (len = MAXSTRING), intent(out) :: string

    IF (threadSupport == MPI_THREAD_SINGLE) THEN
       string = "MPI_THREAD_SINGLE"
    ELSE IF (threadSupport == MPI_THREAD_FUNNELED) THEN
       string = "MPI_THREAD_FUNNELED"
    ELSE IF (threadSupport == MPI_THREAD_SERIALIZED) THEN
       string = "MPI_THREAD_SERIALIZED"
    ELSE IF (threadSupport == MPI_THREAD_MULTIPLE) THEN
       string = "MPI_THREAD_MULTIPLE"
    END IF
    
  END SUBROUTINE threadSupportToString

END MODULE output
