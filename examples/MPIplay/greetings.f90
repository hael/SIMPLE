program greetings
use mpi_f08
implicit none
integer            :: myid, ntasks, source, dest, tag, ierr, islave
integer            :: size
TYPE(MPI_Status)   :: status
character(len=100) :: message
character(len=10)  :: digit_string
integer, parameter :: master = 0

! always start with initialising the MPI env
call MPI_Init(ierr) 
call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
! first argument: MPI_COMM_WORLD is a communicator (collection of processes that can send messages to each other)
! myid is now process rank; N processes means myid 0..N-1; In SIMPLE terminology: part = myid + 1
call MPI_Comm_size(MPI_COMM_WORLD, ntasks, ierr)
! takes the communicator and returns the number of processors p = nparts

if( myid == 0 )then
    write(digit_string,'(i3)') myid
    message ='Greetings from process '//trim(digit_string)//' !'
    tag     = 0
    do islave=1,ntasks-1
        call MPI_Send(message, 100, MPI_CHARACTER, islave, tag, MPI_COMM_WORLD, ierr)
    end do
    ! message sent
    ! message        = the message sent
    ! 100 = count    = # characters
    ! MPI_CHARACTER  = datatype
    ! islave         = id of recieving process (cannot be wildcard)
    ! tag            = integer tag (cannot be wildcard)
    ! MPI_COMM_WORLD = communicator
else    
    tag = 0
    call MPI_Recv(message, 100, MPI_CHARACTER, master, tag, MPI_COMM_WORLD, status, ierr)
    ! message recieved
    ! message        = the message recieved
    ! 100 = count    = # characters
    ! MPI_CHARACTER  = datatype
    ! master         = id of sending process (can be wildcard: MPI_ANY_SOURCE)
    ! tag            = integer tag (can be wildcard: MPI_ANY_TAG)
    ! MPI_COMM_WORLD = communicator
    ! status         = indicates whether data was recieved
    write(*,'(a)') trim(message)
endif
! in order for process A (any rank /= 0) to send a message to process B (rank = 0), the
! argument comm that A uses in MPI_Send must be identical to the argument that B uses in
! MPI_Recv

! always end with Finalize
call MPI_Finalize(ierr) 

end program greetings

! MPI DATATYPES
! MPI_INTEGER          = integer
! MPI_REAL             = real
! MPI_DOUBLE_PRECISION = double precision
! MPI_COMPLEX          = complex
! MPI_LOGICAL          = logical
! MPI_CHARACTER        = character
! MPI_BYTE
! MPI_PACKED