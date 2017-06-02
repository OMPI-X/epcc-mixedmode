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
! Setup module for Mixed mode benchmark program.            !
! Contains variables for the settings of the benchmark and  !
! routines to setup the benchmark.                          !
!-----------------------------------------------------------!
MODULE benchmarkSetup
  use parallelEnvironment
  implicit none

!Module variables
  integer :: warmUpIters !number of iterations of warmup loop
  integer :: defaultReps !the default number of repetitions
  integer :: repsToDo !reps to do for a benchmark
  integer :: minDataSize 
  integer :: maxDataSize
  
  logical :: benchComplete
  DOUBLE PRECISION :: targetTime !threshold time for benchmark
  !variables for timing
  DOUBLE PRECISION :: startTime,finishTime, totalTime
 
  !benchmark types (constants, don't need to change)
  integer, parameter :: MASTERONLY = 1
  integer, parameter :: FUNNELLED = 2
  integer, parameter :: MULTIPLE = 3
  integer, parameter :: REDUCE = 4
  integer, parameter :: ALLREDUCE = 5
  integer, parameter :: SCATTER = 6
  integer, parameter :: GATHER = 7

  !Longest string constant..used in string decalrations
  integer, parameter :: MAXSTRING = 30

  !for parsing list of benchmarks in the input file
  integer, parameter :: NUM_BENCHMARKS = 22
  character (len = MAXSTRING), dimension(NUM_BENCHMARKS) :: benchmarkList
  integer, parameter :: LASTPPID = 6 !id of last pingpong/pingping bench
  integer, parameter :: LAST_PT_PT_ID = 9 !id of last pt-to-pt bench
  integer, parameter :: LASTMULTIPPID = 15 !id of last multi pt-to-pt bench
  integer :: benchmarkNumber
  integer, parameter :: FINISHED = 999, ERROR = 100
  
  !Input file variables
  integer :: ioStatus
  
CONTAINS
  
  !---------------------------------------------------------!
  ! Subroutine: openFile                                    !
  !                                                         !
  ! Looks for filename as a program argument and attempts   !
  ! to open this file.                                      !
  !---------------------------------------------------------!
  SUBROUTINE openFile()
    integer :: argStatus
    character (len = MAXSTRING) :: fileName

    !Read input file name from command line
    !CALL getarg(1,fileName)
    CALL GET_COMMAND_ARGUMENT(1,fileName,STATUS=argStatus)
    IF(argStatus > 0) THEN
       print *, "ERROR Reading input file from command line."
       print *, "Usage: mpiexec -np mixedModeBenchmark <fileName>"
    END IF

    write(*,fmt="(2x,A,A,A)",advance="no") "Attempting to open '"&
         ,trim(fileName),"'...."

    !Open fileName; assign it to unit 10
    OPEN(UNIT=10, FILE=fileName, STATUS="OLD", &
         ACTION="READ", FORM="FORMATTED", IOSTAT=ioStatus)

    !Check that file opened successfully
    IF(ioStatus > 0) THEN
       write(*,*) "ERROR. IOSTAT =", ioStatus
    ELSE
       write(*,*) "Success."
    END IF

  END SUBROUTINE openFile
  
  !---------------------------------------------------------!
  ! Subroutine: closeFile                                   !
  !                                                         !
  ! Closes the input file.                                  !
  !---------------------------------------------------------!
  SUBROUTINE closeFile()

    !Close file
    CLOSE(10)

  END SUBROUTINE closeFile

  !---------------------------------------------------------!
  ! Subroutine: readBenchmarkParams                         !
  !                                                         !  
  ! Initialises the benchmark parameters.                   !
  ! Reads the minimum and maximum data for the benchmarks   !
  ! from the input file (unit 10).                          !       
  !---------------------------------------------------------!
  SUBROUTINE readBenchmarkParams()

    !Rank 0 reads parameters from input file
    IF (myMPIRank == 0) THEN
       write (*,*) "Reading parameters from input file...."
       !read minimum data size from input file
       read(10,*) minDataSize 
       !read maximum data size from input file
       read(10,*) maxDataSize
       !read target time from input file
       read(10,*) targetTime

       !set other benchmark parameters
       warmUpIters = 2
       defaultReps = 1000
    
       !Report benchmark parameters
       write(*,fmt="(A)") "------------------------------------------"
       write(*,fmt="(A)") "           Benchmark parameters           "
       write(*,fmt="(A)") "------------------------------------------"
       write(*,fmt="(A,t25,i10)") "Minimum data size", minDataSize
       write(*,fmt="(A,t25,i10)") "Maximum data size", maxDataSize
       write(*,fmt="(A,t25,f10.2)") "Target time (sec)", targetTime
       write(*,fmt="(A,t25,i10)") "Default Repetitions", defaultReps
       write(*,fmt="(A,t25,i10)") "No. Warmup iterations", warmUpIters

    END IF
    !Initialise benchmarkNumber to 0 so that the
    !DO WHILE loop in the driver is entered the first time
    benchmarkNumber = 0

    !Broadcast benchmark parameters from master to all 
    !other MPI processes.
    CALL MPI_Bcast(minDataSize, 1, MPI_INTEGER, 0, comm, ierr)
    CALL MPI_Bcast(maxDataSize, 1, MPI_INTEGER, 0, comm, ierr)
    CALL MPI_Bcast(targetTime, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    CALL MPI_Bcast(defaultReps, 1, MPI_INTEGER, 0, comm, ierr)
    CALL MPI_Bcast(warmUpIters, 1, MPI_INTEGER, 0, comm, ierr)

  END SUBROUTINE readBenchmarkParams

  !---------------------------------------------------------!
  ! Subroutine: setupBenchmarkList                          !
  !                                                         !  
  ! Subroutine to setup the benchmarkList array with the    !
  ! list of all possible benchmarks.                        ! 
  !---------------------------------------------------------!
  SUBROUTINE setupBenchmarkList()
    
    !Pingpong benchmarks
    benchmarkList(1) = "masteronlypingpong"
    benchmarkList(2) = "funnelledpingpong"
    benchmarkList(3) = "multiplepingpong"
    !Pingping benchmarks
    benchmarkList(4) = "masteronlypingping"
    benchmarkList(5) = "funnelledpingping"
    benchmarkList(6) = "multiplepingping"
    !Haloexchange benchmarks
    benchmarkList(7) = "masteronlyhaloexchange"
    benchmarkList(8) = "funnelledhaloexchange"
    benchmarkList(9) = "multiplehaloexchange"
    !Multi-Pingpong benchmarks
    benchmarkList(10) = "masteronlymultipingpong"
    benchmarkList(11) = "funnelledmultipingpong"
    benchmarkList(12) = "multiplemultipingpong"
    !Multi-Pingping benchmarks
    benchmarkList(13) = "masteronlymultipingping" 
    benchmarkList(14) = "funnelledmultipingping"
    benchmarkList(15) = "multiplemultipingping"
    !Collective benchmarks
    benchmarkList(16) = "barrier"
    benchmarkList(17) = "reduce"
    benchmarkList(18) = "allreduce"
    benchmarkList(19) = "broadcast"
    benchmarkList(20) = "scatter"
    benchmarkList(21) = "gather"
    benchmarkList(22) = "alltoall"

  END SUBROUTINE setupBenchmarkList
  
  !---------------------------------------------------------!
  ! Subroutine: findBenchmarkNumber                         !
  !                                                         !
  ! Finds the ID of the next benchmark which will be        !
  ! executed. Master MPI process reads next line from input !
  ! file. It then compares it to the benchmark list to find !
  ! the ID and broadcasts this to the other MPI processes.  !
  !---------------------------------------------------------!
  SUBROUTINE findBenchmarkNumber()
    character (len = MAXSTRING) :: benchmarkName
    integer :: rankInA, rankInB
    integer :: i

    !Master MPI process reads next line from file
    IF (myMPIRank == 0) THEN
       !set benchmarkNumber to ERROR before read 
       !to allow error check
       benchmarkNumber = ERROR
       
       !Read next benchmark from file
       READ(10,*,IOSTAT=ioStatus) benchmarkName
       
       !Check if EOF is reached
       IF (ioStatus < 0) THEN
          benchmarkNumber = FINISHED
       ELSE
          !Convert benchmarkName to lowercase characters
          CALL ConvertTolowercase(benchmarkName)
          !..and check if benchmark name matches.
          DO i = 1,NUM_BENCHMARKS
             IF (benchmarkName == benchmarkList(i)) THEN
                benchmarkNumber = i
             END IF
          END DO
       END IF

       !Check if benchmark Name does not match
       IF (benchmarkNumber == ERROR) THEN
          write(*,*) "ERROR: ", trim(benchmarkName) , &
               " does not match any possible benchmarks"
       END IF

       !Check if pingpong or pingping benchmark
       IF (benchmarkNumber <= LASTPPID) THEN 
          !Read ranks from input file
          READ(10,*,IOSTAT=ioStatus) rankInA, rankInB
          
          !Check if error in read
          IF (ioStatus < 0) THEN
             write(*,*) "ERROR: expecting ranks after ",&
                  trim(benchmarkName)
          ELSE !if no error find actual MPI ranks
             PPRanks(1) = findRank(rankInA)
             PPRanks(2) = findRank(rankInB)
          END IF
          !Check if PPRanks are the same
          IF (PPRanks(1) == PPRanks(2)) THEN
             write(*,*) "Warning: Ranks are the same; benchmark will",&
                  " not work."
          END IF
       END IF
       
    END IF

    !Broadcast benchmarkNumber to other MPI processes
    CALL MPI_Bcast(benchmarkNumber, 1, MPI_INTEGER, 0, &
         comm, ierr)

    !If pingpong or pingping benchmark broadcast ranks of
    !participating processes.
    IF (benchmarkNumber <= LASTPPID) THEN
       CALL MPI_Bcast(PPRanks, 2, MPI_INTEGER, 0, comm, ierr)
    END IF

  END SUBROUTINE findBenchmarkNumber

  !---------------------------------------------------------!
  ! Subroutine: convertToLowerCase                          !
  !                                                         !
  ! Takes a string as an agrument and converts all          !
  ! uppercase characters to lowercase using its ASCII value.!
  !                                                         !
  ! Adopted from: Fortran 90/95 for Scientists & Engineers  !
  !               by S.J. Chapman.(2nd Ed. Chp 10)          !
  !---------------------------------------------------------!
  SUBROUTINE convertToLowerCase(string)
    character (len = *), intent(inout) :: string
    integer :: i, length

    !Find length of string.
    length = LEN(string)
    
    !Loop through each character of string.
    DO i = 1, length
       !If character between A and Z...
       IF((string(i:i) >= 'A') .AND. (string(i:i)) <= 'Z') THEN
          !Make character lowercase
          string(i:i) = ACHAR(IACHAR(string(i:i)) + 32)
       END IF
    END DO

  END SUBROUTINE convertToLowerCase

  !---------------------------------------------------------!
  ! Function: repTimeCheck                                  !
  !                                                         !
  ! Checks if the time for the benchmark reached the target !
  ! time. Changes the number of repetitions for the next    !
  ! data size based on the difference between the time      !
  ! taken and the target time.                              !
  !---------------------------------------------------------!
  FUNCTION repTimeCheck(time, numReps)
    DOUBLE PRECISION, intent(in) :: time
    integer, intent(in) :: numReps
    logical :: repTimeCheck

    IF (time < targetTime) THEN
       !double repsToDo and repeat benchmark
       repsToDo = 2*numReps
       repTimeCheck = .false.
    ELSE IF (time > (2*targetTime)) THEN
       !finish benchmark and half number of reps for next dataSize
       repsToDo = MAX(numReps/2,1) !repsToDo is at least 1 
       repTimeCheck = .true.
    ELSE !time is >= targetTime
       !finish benchmark and keep reps for next data size
       repTimeCheck = .true.
    END IF

   
  END FUNCTION repTimeCheck

 
END MODULE benchmarkSetup
