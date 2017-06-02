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
! Contains variables and routines for the parallel          !
! environment.                                              !
! Routines to setup the MPI and OpenMP programming          ! 
! environment.                                              !
!-----------------------------------------------------------!
MODULE parallelEnvironment
  
  use omp_lib
  use mpi
  !INCLUDE 'mpif.h'
  implicit none
  
!Module variable specification
  !MPI variables
  integer :: myMPIRank, numMPIprocs 
  integer :: comm, commCart
  integer :: localComm, crossComm, localCommSize
  integer :: localCommRank, crossCommRank
  integer :: ierr
  integer :: status(MPI_STATUS_SIZE)
  character (len = MPI_MAX_PROCESSOR_NAME) :: myProcName
  integer :: procNameLen
  integer :: requestID
  integer :: requestArray(4) !for haloexchange
  integer :: statusArray(MPI_STATUS_SIZE,4) !for haloexchange
  integer, parameter :: tag = 1 !set tag to match messages
  integer :: leftNeighbour, rightNeighbour
  
  integer :: threadSupport !Level of thread support by implementation
  integer, dimension(2) :: PPRanks !ranks for pingpong or pingping

  !OpenMP variables
  integer :: myThreadID, numThreads
  !make myThreadID a thread private variable
!$OMP THREADPRIVATE(myThreadID)

  !Array to hold the global ID for each thread
  integer, dimension(:), allocatable :: globalIDarray

  integer :: sizeInteger !holds number of bytes for an int
CONTAINS
!Module procedure specification

  !---------------------------------------------------------!
  ! Subroutine: initParallelEnv                             !
  !                                                         !
  ! Initialises the MPI and OpenMP environments.            !
  ! Finds the total number of MPI processes and             !
  ! OpenMP threads.                                         !
  ! Also finds the ID of each MPI process and OpenMP thread.! 
  !---------------------------------------------------------!
  SUBROUTINE initParallelEnv()

    !setup MPI programming environment 
    CALL MPI_Init_thread(MPI_THREAD_MULTIPLE,threadSupport,ierr)

    comm = MPI_COMM_WORLD
    CALL MPI_Comm_size(comm, numMPIprocs, ierr)
    CALL MPI_Comm_rank(comm, myMPIRank, ierr)

    !Find the number of bytes for an int (numMPIprocs)
    CALL MPI_Type_size(MPI_INTEGER, sizeInteger, ierr)

    !Find the processor name of each MPI process.
    CALL  MPI_Get_processor_name(myProcName, procNameLen, ierr)

    !Use processor name to create a communicator across 
    !node boundaries.
    CALL setupCommunicators()

    !setup OpenMP programming environment
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP SHARED(numThreads,globalIDarray,myMPIRank)

    numThreads = omp_get_num_threads()
    myThreadID = omp_get_thread_num() + 1 !threadID from 1 to totalThreads

    !Allocate space for globalIDarray
!$OMP SINGLE
    allocate(globalIDarray(numThreads))
!$OMP END SINGLE

    !Calculate the globalID for each thread
    globalIDarray(myThreadID) = (myMPIRank * numThreads) + myThreadID

!$OMP END PARALLEL
  
  END SUBROUTINE initParallelEnv

  !---------------------------------------------------------!
  ! Subroutine: finaliseParallelEnv                         !
  !                                                         !
  ! Closes the MPI programming environment.                 !
  !                                                         ! 
  !---------------------------------------------------------!
  SUBROUTINE finaliseParallelEnv()

    !finalise the MPI programming environment
    CALL MPI_Finalize(ierr)
    !free the space created for globalIDarray...
    deallocate(globalIDarray)
    
  END SUBROUTINE finaliseParallelEnv
  
  !---------------------------------------------------------!
  ! Subroutine: findNeighbourRanks                          !
  !                                                         !
  ! This creates a cartesian topology and finds the left    !
  ! and right neighbours of each process.                   !
  !---------------------------------------------------------!
  SUBROUTINE findNeighbourRanks()
    integer :: dims(1) !dims array for MPI_Dims_Create
    logical, parameter :: PERIODS(1) = (/.true./), REORDER = .false.

    !find a good process distribution
    dims = 0 !zero so that dims_create tries to rearrange
    CALL MPI_Dims_Create(numMPIProcs,1,dims,ierr)

    !Create the cartesian topology
    CALL MPI_Cart_Create(comm,1,dims,PERIODS,REORDER, &
         commCart, ierr)

    !Find the ranks of left and right neighbour
    CALL MPI_Cart_Shift(commCart, 0, 1, leftNeighbour, &
         rightNeighbour, ierr)

  END SUBROUTINE findNeighbourRanks
  
  !----------------------------------------------------------!
  ! Function: benchmarkSupport                               !
  !                                                          !
  ! This function compares the level of thread support       !
  ! needed by a particular benchmark with the level provided !
  ! by the implementation.                                   !
  !----------------------------------------------------------!
  FUNCTION benchmarkSupport(required)
    integer, intent(in) :: required
    logical :: benchmarkSupport

    IF (required <= threadSupport) THEN
       benchmarkSupport = .true.
    ELSE
       benchmarkSupport = .false.
    END IF

  END FUNCTION benchmarkSupport

  !-----------------------------------------------------------!
  ! Function: findRank                                        !
  !                                                           !
  ! Finds the MPI ranks which will take part in the pingping  !
  ! or pingpong benchmarks based on the numbers read from the !
  ! input file.                                               !
  !-----------------------------------------------------------!
  FUNCTION findRank(rankIn)
    integer, intent(in) ::rankIn
    integer :: findRank

    !Figure out actual MPI rank
    IF (rankin < 0) THEN
       findRank = numMPIprocs + rankin
    ELSE
       findRank = rankIn
    END IF

    !Check if findRank is too big or still -ve
    IF (findRank > (numMPIprocs-1)) THEN
       !write(*,*) "Warning: Rank input greater than total",&
       !     "process count. Using Rank = ", numMPIprocs-1
       findRank = numMPIprocs - 1
    ELSE IF (findRank < 0) THEN
       !write(*,*) "Warning: MPI process offset greater than",&
       !     "total process count. Using Rank = 0"
       findRank = 0
    END IF

  END FUNCTION findRank
  
  !-----------------------------------------------------------!
  ! Function: compareProcNames                                !
  !                                                           !
  ! Compares the names of 2 processes to check if they are on !
  ! the same node or not.                                     !
  !-----------------------------------------------------------!
  FUNCTION compareProcNames(rankA, rankB)
    integer, intent(in) :: rankA, rankB
    logical compareProcNames
    character (len = MPI_MAX_PROCESSOR_NAME) :: recvProcName

    !Rank B sends procName to Rank A
    IF (myMPIRank == rankB) THEN
       CALL MPI_Send(myProcName, MPI_MAX_PROCESSOR_NAME, &
            MPI_CHARACTER, rankA, tag, comm, ierr)
    ELSE IF (myMPIRank == rankA) THEN
       CALL MPI_Recv(recvProcName, MPI_MAX_PROCESSOR_NAME, &
            MPI_CHARACTER, rankB, tag, comm, status, ierr)
       !rankB compares the two processor names
       IF (myProcName == recvProcName) THEN
          compareProcNames = .true.
       ELSE
          compareProcNames = .false.
       END IF
    END IF

    !Rank A then broadcasts its compareProcNames value to 
    !the other processes
    CALL MPI_Bcast(compareProcNames, 1, MPI_LOGICAL, rankA, &
         comm, ierr)

  END FUNCTION compareProcNames

  !---------------------------------------------------------!
  ! Subroutine: setupCommunicators                          !
  !                                                         !
  ! This creates two new communicators.                     !
  ! The first gives a local communicator for processes on   !
  ! the same node.                                          !
  ! The second uses the local rank to give a communicator   !
  ! across node boundaries.                                 !
  !                                                         !
  ! e.g. for 16 nodes each with 2 processors, this routine  !
  !      will give 16 local communicators of size 2 and     !
  !      2 communicators of size 2 across nodes.            !
  !---------------------------------------------------------!
  SUBROUTINE setupCommunicators()
    integer :: procHash
    
    !Get hash based on processor name
    procHash = procNameToHash()

    !Comm_split using procHash as colour to get 
    !local cmmunicator.
    CALL MPI_Comm_split(comm, procHash, 0, localComm, ierr)

    !Find ranks of processes in localComm
    CALL MPI_Comm_rank(localComm, localCommRank, ierr)

    !Find the size of localComm (for use in calculating multi datasize)
    CALL MPI_Comm_size(localComm, localCommSize, ierr)
    
    !Use localRank as colour to get communicator across nodes.
    CALL MPI_Comm_split(comm, localCommRank, 0, crossComm, ierr)

    !Find ranks of processes in crossComm
    CALL MPI_Comm_rank(crossComm, crossCommRank, ierr)
    
  END SUBROUTINE setupCommunicators

  !---------------------------------------------------------!
  ! Function: procNameToHash                                !
  !                                                         !
  ! Creates an integer hash for each process.               !
  ! Each process on the same node will have the same hash   !
  ! value.                                                  !
  !---------------------------------------------------------!
  FUNCTION procNameToHash()
    integer :: procNameToHash
    integer :: i

    !Initialise hash to 0
    procNameToHash = 0

    DO i = 1, procNameLen
       
       procNameToHash = 7 * procNameToHash + &
            ICHAR(myProcName(i:i))
    END DO
    
  END FUNCTION procNameToHash
  
  !---------------------------------------------------------!
  ! Subroutine: exchangeWorldRanks                          !
  !                                                         !
  ! Finds the MPI_COMM_WORLD ranks of the processes         !
  ! participating in the multi-pingpong and multi-pingping  !
  ! benchmarks.                                             !
  !---------------------------------------------------------!
  SUBROUTINE exchangeWorldRanks(nodeA, nodeB, otherWorldRank)
    integer, intent(in) :: nodeA, nodeB
    integer, intent(out) :: otherWorldRank
    integer :: destRank

    IF (crossCommRank == nodeA) THEN
       destRank = nodeB
    ELSE IF (crossCommRank == nodeB) THEN
       destRank = nodeA
    END IF

    IF (crossCommRank == nodeA .or. crossCommRank == nodeB) THEN
       
       !Start send of COMM_WORLD rank to destRank in crossComm
       CALL MPI_Isend(myMPIRank, 1, MPI_INTEGER, destRank, &
            tag, crossComm, requestID, ierr)

       !Then wait for message from destRank and store in 
       !otherWorldRank.
       CALL MPI_Recv(otherWorldRank, 1, MPI_INTEGER, destRank, &
            tag, crossComm, status, ierr)

       !Finish the send with an MPI_Wait
       CALL MPI_Wait(requestID, status, ierr)

    END IF
    
  END SUBROUTINE exchangeWorldRanks
  
  !---------------------------------------------------------!
  ! Subroutine: sendProcName                                !
  !                                                         !
  ! Sends the processor name from processes in destNode     !
  ! of crossComm to srcNode.                                !
  !---------------------------------------------------------!
  SUBROUTINE sendProcName(destNode, srcNode, destProcName)
    integer, intent(in) :: srcNode, destNode
    character (len = MPI_MAX_PROCESSOR_NAME), intent(out) :: destProcName
    
    !MPI Processes under srcNode of crossComm send their 
    !processor name to destNode.
    IF (crossCommRank == srcNode) THEN
       CALL MPI_Send(myProcName, MPI_MAX_PROCESSOR_NAME, &
            MPI_CHARACTER, destNode, tag, crossComm, ierr)
    ELSE IF (crossCommRank == destNode) THEN
       CALL MPI_Recv(destProcName, MPI_MAX_PROCESSOR_NAME, &
            MPI_CHARACTER, srcNode, tag, crossComm, status, ierr)
    END IF
    
  END SUBROUTINE sendProcName

  !---------------------------------------------------------!
  ! Function: checkCrossCommBalance                         !
  !                                                         !
  ! Checks if there's a balance in the number of processes  !
  ! in crossComm nodes.                                     !
  !---------------------------------------------------------!
  FUNCTION crossCommBalance(nodeA, nodeB)
    integer, intent(in) :: nodeA, nodeB
    integer :: localCommSize, otherLocalCommSize
    logical :: crossCommBalance

    !Find the size of localComm
    CALL MPI_Comm_size(localComm, localCommSize, ierr)

    !master process on nodeB sends localCommSize
    IF (crossCommRank == nodeB .and. localCommRank == 0) THEN
       CALL MPI_Send(localCommSize, 1, MPI_INTEGER, nodeA, &
            tag, crossComm, ierr)
    !master process on nodeA.....
    ELSEIF (crossCommRank == nodeA .and. localCommRank == 0) THEN
       !1) receives nodeB's localCommSize
       CALL MPI_Recv(otherLocalCommSize, 1, MPI_INTEGER, nodeB, &
            tag, crossComm, status, ierr)
       
       !2) Test for balance by comparing otherLocalCommSize 
       !   to localCommSize.
       IF (localCommSize == otherLocalCommSize) THEN
          !Set balance to .true.
          crossCommBalance = .true.
       ELSE 
          crossCommBalance = .false.
       END IF

       !3) Send balance to master commWorld process.
       !Only need to send if not master commWorld process
       IF (myMPIRank /= 0) THEN
          CALL MPI_Send(crossCommBalance, 1, MPI_LOGICAL, &
               0, tag, comm, ierr)
       END IF

    END IF
    
    !Master commWorld process....
    IF (myMPIRank == 0) THEN
       !Receives balance variable if not same process as 
       !master process on nodeA.
       IF (crossCommRank /= nodeA .and. localCommRank /= 0) THEN

          CALL MPI_Recv(crossCommBalance, 1, MPI_LOGICAL, &
               MPI_ANY_SOURCE, tag, comm, status, ierr)
       END IF

    END IF
       
    !Broadcast balance to all processes
    CALL MPI_Bcast(crossCommBalance, 1, MPI_LOGICAL, 0, comm, ierr)

  END FUNCTION crossCommBalance
  
END MODULE parallelEnvironment
