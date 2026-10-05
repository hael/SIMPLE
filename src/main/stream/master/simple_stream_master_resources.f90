!@descr: the computing resources the stream master gives each stage: named defaults, overridden per stage by environment variables
!==============================================================================
! MODULE: simple_stream_master_resources
!
! PURPOSE:
!   One table of the threads and parts of every stage and of the jobs the
!   stages run (stream fix plan, decision 14). Each value has a named default
!   below; a stage's environment variables override them, read once by the
!   master (stream_resources_from_env), which passes the values on the stage
!   command lines and logs the table at start. The stages hold no resource
!   literals of their own.
!
!   Environment variables, one set per stage (SIMPLE_STREAM_<STAGE>_...):
!     PREPROC  NTHR, NPARTS  preprocessing jobs (p01)
!     REFGEN   NTHR, NPARTS  the initial analysis' jobs (p03): its 2D jobs and
!                            its sieve's chunks, its 3D job and reprojection
!                            (threads), and the parts of its 2D and 3D jobs
!     PICK     NTHR, NPARTS  reference-picking jobs (p04)
!     CHUNK    NTHR, NPARTS  particle sieving's chunks (p05): threads per chunk
!                            job, and the number of chunks run at once
!     POOL     NTHR, NPARTS  the 2D pool's iterations (p06)
!     SOLVE3D  NTHR, NPARTS  the multistate 3D jobs (p07)
!   Each stage reads its own PARTITION variable for its queue environment.
!   A value the user gives on the master's command line wins over both
!   (the picking references' 3D threads, nthr3D_pickrefs).
!==============================================================================
module simple_stream_master_resources
use simple_defs,             only: logfhandle
use simple_string_utils,     only: str2int
use simple_defs_environment, only: SIMPLE_STREAM_PREPROC_NTHR, SIMPLE_STREAM_PREPROC_NPARTS, SIMPLE_STREAM_REFGEN_NTHR,&
                                  &SIMPLE_STREAM_REFGEN_NPARTS, SIMPLE_STREAM_PICK_NTHR, SIMPLE_STREAM_PICK_NPARTS,&
                                  &SIMPLE_STREAM_CHUNK_NTHR, SIMPLE_STREAM_CHUNK_NPARTS, SIMPLE_STREAM_POOL_NTHR,&
                                  &SIMPLE_STREAM_POOL_NPARTS, SIMPLE_STREAM_SOLVE3D_NTHR, SIMPLE_STREAM_SOLVE3D_NPARTS
implicit none

public :: stream_resources, stream_resources_from_env
private

! the defaults
integer, parameter :: PREPROCESS_NTHR          = 4   ! per preprocessing job
integer, parameter :: PREPROCESS_NPARTS        = 16  ! preprocessing jobs at once
integer, parameter :: OPTICS_NTHR              = 1
integer, parameter :: INITIAL_ANALYSIS_NTHR    = 32  ! the initial analysis' own process (in-process picking)
integer, parameter :: INITIAL_ANALYSIS_NTHR2D  = 16  ! its 2D jobs and its sieve's chunks
integer, parameter :: INITIAL_ANALYSIS_NTHR3D  = 16  ! its 3D job and reprojection
integer, parameter :: INITIAL_ANALYSIS_NPARTS  = 1   ! parts of its 2D and 3D jobs
integer, parameter :: INITIAL_ANALYSIS_NCHUNKS = 4   ! its sieve's chunks at once
integer, parameter :: REFPICK_NTHR             = 8
integer, parameter :: REFPICK_NPARTS           = 8
integer, parameter :: SIEVE_NTHR               = 16  ! per chunk job
integer, parameter :: SIEVE_NCHUNKS            = 4   ! chunks at once
integer, parameter :: POOL2D_NTHR              = 8
integer, parameter :: POOL2D_NPARTS            = 6
integer, parameter :: SOLVE3D_NTHR             = 8   ! the 3D stage's own process
integer, parameter :: SOLVE3D_NTHR3D           = 8   ! per 3D job
integer, parameter :: SOLVE3D_NPARTS3D         = 8   ! parts of each 3D job

!> The resources of every stage and of its jobs
type :: stream_resources
    integer :: preprocess_nthr          = PREPROCESS_NTHR
    integer :: preprocess_nparts        = PREPROCESS_NPARTS
    integer :: optics_nthr              = OPTICS_NTHR
    integer :: initial_analysis_nthr    = INITIAL_ANALYSIS_NTHR
    integer :: initial_analysis_nthr2D  = INITIAL_ANALYSIS_NTHR2D
    integer :: initial_analysis_nthr3D  = INITIAL_ANALYSIS_NTHR3D
    integer :: initial_analysis_nparts  = INITIAL_ANALYSIS_NPARTS
    integer :: initial_analysis_nchunks = INITIAL_ANALYSIS_NCHUNKS
    integer :: refpick_nthr             = REFPICK_NTHR
    integer :: refpick_nparts           = REFPICK_NPARTS
    integer :: sieve_nthr               = SIEVE_NTHR
    integer :: sieve_nchunks            = SIEVE_NCHUNKS
    integer :: pool2D_nthr              = POOL2D_NTHR
    integer :: pool2D_nparts            = POOL2D_NPARTS
    integer :: solve3D_nthr             = SOLVE3D_NTHR
    integer :: solve3D_nthr3D           = SOLVE3D_NTHR3D
    integer :: solve3D_nparts3D         = SOLVE3D_NPARTS3D
contains
    procedure :: log => log_resources
end type stream_resources

contains

    !> The defaults, each overridden by its stage's environment variable when that is set to a
    !! positive integer; a variable that is not one is ignored with a warning.
    function stream_resources_from_env() result( res )
        type(stream_resources) :: res
        call env_override(SIMPLE_STREAM_PREPROC_NTHR,   res%preprocess_nthr)
        call env_override(SIMPLE_STREAM_PREPROC_NPARTS, res%preprocess_nparts)
        call env_override(SIMPLE_STREAM_REFGEN_NTHR,    res%initial_analysis_nthr2D)
        call env_override(SIMPLE_STREAM_REFGEN_NTHR,    res%initial_analysis_nthr3D)
        call env_override(SIMPLE_STREAM_REFGEN_NPARTS,  res%initial_analysis_nparts)
        call env_override(SIMPLE_STREAM_PICK_NTHR,      res%refpick_nthr)
        call env_override(SIMPLE_STREAM_PICK_NPARTS,    res%refpick_nparts)
        call env_override(SIMPLE_STREAM_CHUNK_NTHR,     res%sieve_nthr)
        call env_override(SIMPLE_STREAM_CHUNK_NPARTS,   res%sieve_nchunks)
        call env_override(SIMPLE_STREAM_POOL_NTHR,      res%pool2D_nthr)
        call env_override(SIMPLE_STREAM_POOL_NPARTS,    res%pool2D_nparts)
        call env_override(SIMPLE_STREAM_SOLVE3D_NTHR,   res%solve3D_nthr3D)
        call env_override(SIMPLE_STREAM_SOLVE3D_NPARTS, res%solve3D_nparts3D)

    contains

        subroutine env_override( name, val )
            character(len=*), intent(in)    :: name
            integer,          intent(inout) :: val
            character(len=64) :: env_val
            integer           :: envlen, ival, ios
            call get_environment_variable(name, env_val, envlen)
            if( envlen <= 0 ) return
            ival = str2int(trim(env_val), ios)
            if( ios /= 0 .or. ival <= 0 )then
                write(logfhandle,'(A,A,A,A)') '>>> WARNING: IGNORING ', name, '=', trim(env_val)
                return
            endif
            val = ival
        end subroutine env_override

    end function stream_resources_from_env

    !> The table, to the log.
    subroutine log_resources( self )
        class(stream_resources), intent(in) :: self
        write(logfhandle,'(A)')       '>>> STREAM RESOURCES (THREADS / PARTS)'
        write(logfhandle,'(A,I4,A,I4)') '>>>   preprocessing jobs            ', self%preprocess_nthr, ' / ', self%preprocess_nparts
        write(logfhandle,'(A,I4)')      '>>>   optics assignment             ', self%optics_nthr
        write(logfhandle,'(A,I4)')      '>>>   initial analysis              ', self%initial_analysis_nthr
        write(logfhandle,'(A,I4,A,I4)') '>>>     2D jobs and sieve chunks    ', self%initial_analysis_nthr2D, ' / ',&
            &self%initial_analysis_nparts
        write(logfhandle,'(A,I4,A,I4)') '>>>     3D job and reprojection     ', self%initial_analysis_nthr3D, ' / ',&
            &self%initial_analysis_nparts
        write(logfhandle,'(A,I4)')      '>>>     sieve chunks at once        ', self%initial_analysis_nchunks
        write(logfhandle,'(A,I4,A,I4)') '>>>   reference-picking jobs        ', self%refpick_nthr, ' / ', self%refpick_nparts
        write(logfhandle,'(A,I4,A,I4)') '>>>   sieve chunks (each / at once) ', self%sieve_nthr, ' / ', self%sieve_nchunks
        write(logfhandle,'(A,I4,A,I4)') '>>>   pool 2D iterations            ', self%pool2D_nthr, ' / ', self%pool2D_nparts
        write(logfhandle,'(A,I4)')      '>>>   multistate 3D                 ', self%solve3D_nthr
        write(logfhandle,'(A,I4,A,I4)') '>>>     3D jobs                     ', self%solve3D_nthr3D, ' / ', self%solve3D_nparts3D
    end subroutine log_resources

end module simple_stream_master_resources
