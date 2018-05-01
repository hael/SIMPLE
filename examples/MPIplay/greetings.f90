program greetings
use mpi_f08
implicit none
integer            :: my_rank, p, source, dest, tag, ierr
integer            :: size
TYPE(MPI_Status)   :: status
character(len=100) :: message
character(len=10)  :: digit_string
call MPI_Init(ierr) ! always start with Init
call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
! first argument: MPI_COMM_WORLD is a communicator (collection of processes that can send messages to each other)
! my_rank is now process rank; N processes means my_rank 0..N-1; In SIMPLE terminology: part = my_rank + 1
call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)
! takes the communicator and returns the number of processors p = nparts
if( my_rank /= 0 )then
    write(digit_string,'(i3)') my_rank
    message ='Greetings from process '//trim(digit_string)//' !'
    dest    = 0
    tag     = 0
    call MPI_Send(message, 100, MPI_CHARACTER, dest, tag, MPI_COMM_WORLD, ierr)
    ! message sent
    ! message        = the message sent
    ! 100 = count    = # characters
    ! MPI_CHARACTER  = datatype
    ! dest           = rank of recieving process (cannot be wildcard)
    ! tag            = integer tag (cannot be wildcard)
    ! MPI_COMM_WORLD = communicator
else
    do source=1,p-1
        tag = 0
        call MPI_Recv(message, 100, MPI_CHARACTER, source, tag, MPI_COMM_WORLD, status, ierr)
        ! message recieved
        ! message        = the message recieved
        ! 100 = count    = # characters
        ! MPI_CHARACTER  = datatype
        ! source         = rank of sending process (can be wildcard: MPI_ANY_SOURCE)
        ! tag            = integer tag (can be wildcard: MPI_ANY_TAG)
        ! MPI_COMM_WORLD = communicator
        ! status         = indicates whether data was recieved
        write(*,'(a)') trim(message)
    end do
endif
! in order for process A (any rank /= 0) to send a message to process B (rank = 0), the
! argument comm that A uses in MPI_Send must be identical to the argument that B uses in
! MPI_Recv
call MPI_Finalize(ierr) ! always end with Finalize

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