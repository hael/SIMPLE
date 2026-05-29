!@descr: coarray-backend matrix-sum parallelization strategy test
module simple_coarray_matrix_sum_strategy
use simple_core_module_api
use simple_cmdline,      only: cmdline
use simple_parameters,   only: parameters
use simple_qsys_env,     only: qsys_env
use simple_qsys_funs,    only: qsys_job_finished
implicit none

public :: coarray_matrix_sum_strategy
public :: coarray_matrix_sum_master_strategy
public :: coarray_matrix_sum_worker_strategy
public :: create_coarray_matrix_sum_strategy
private
#include "simple_local_flags.inc"

integer, parameter :: DEFAULT_MATRIX_N = 4
integer, parameter :: DEFAULT_NPARTS   = 4
integer, parameter :: JOB_DESCR_SIZE   = 16

type, abstract :: coarray_matrix_sum_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
    procedure(endmsg_interface),   deferred :: end_message
end type coarray_matrix_sum_strategy

type, extends(coarray_matrix_sum_strategy) :: coarray_matrix_sum_master_strategy
    integer :: n       = DEFAULT_MATRIX_N
    integer :: nparts  = DEFAULT_NPARTS
    integer :: ncunits = DEFAULT_NPARTS
    integer :: element_sum = 0
    integer :: total_sum   = 0
contains
    procedure :: initialize   => master_initialize
    procedure :: execute      => master_execute
    procedure :: finalize_run => master_finalize_run
    procedure :: cleanup      => master_cleanup
    procedure :: end_message  => master_end_message
end type coarray_matrix_sum_master_strategy

type, extends(coarray_matrix_sum_strategy) :: coarray_matrix_sum_worker_strategy
    integer :: n       = DEFAULT_MATRIX_N
    integer :: part    = 1
    integer :: nparts  = 1
    integer :: image   = 0
    integer :: nimages = 0
    integer :: part_sum = 0
contains
    procedure :: initialize   => worker_initialize
    procedure :: execute      => worker_execute
    procedure :: finalize_run => worker_finalize_run
    procedure :: cleanup      => worker_cleanup
    procedure :: end_message  => worker_end_message
end type coarray_matrix_sum_worker_strategy

abstract interface
    subroutine init_interface( self, cline )
        import :: coarray_matrix_sum_strategy, cmdline
        class(coarray_matrix_sum_strategy), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface( self )
        import :: coarray_matrix_sum_strategy
        class(coarray_matrix_sum_strategy), intent(inout) :: self
    end subroutine exec_interface

    subroutine finalize_interface( self )
        import :: coarray_matrix_sum_strategy
        class(coarray_matrix_sum_strategy), intent(inout) :: self
    end subroutine finalize_interface

    subroutine cleanup_interface( self )
        import :: coarray_matrix_sum_strategy
        class(coarray_matrix_sum_strategy), intent(inout) :: self
    end subroutine cleanup_interface

    function endmsg_interface( self ) result( msg )
        import :: coarray_matrix_sum_strategy
        class(coarray_matrix_sum_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
    end function endmsg_interface
end interface

contains

    function create_coarray_matrix_sum_strategy( cline ) result( strategy )
        class(cmdline), intent(in) :: cline
        class(coarray_matrix_sum_strategy), allocatable :: strategy
        if( cline%defined('part') )then
            allocate(coarray_matrix_sum_worker_strategy :: strategy)
        else
            allocate(coarray_matrix_sum_master_strategy :: strategy)
        endif
    end function create_coarray_matrix_sum_strategy

    subroutine master_initialize( self, cline )
        class(coarray_matrix_sum_master_strategy), intent(inout) :: self
        class(cmdline),                            intent(inout) :: cline
        self%n       = DEFAULT_MATRIX_N
        self%nparts  = DEFAULT_NPARTS
        if( cline%defined('box')    ) self%n      = cline%get_iarg('box')
        if( cline%defined('nparts') ) self%nparts = cline%get_iarg('nparts')
        if( self%n < 1 )      THROW_HARD('box must be > 0 for coarray matrix sum test')
        if( self%nparts < 1 ) THROW_HARD('nparts must be > 0 for coarray matrix sum test')
        self%ncunits = self%nparts
        if( cline%defined('ncunits') ) self%ncunits = cline%get_iarg('ncunits')
        if( self%ncunits < 1 ) THROW_HARD('ncunits must be > 0 for coarray matrix sum test')
        self%element_sum = 0
        self%total_sum   = 0
    end subroutine master_initialize

    subroutine master_execute( self )
        class(coarray_matrix_sum_master_strategy), intent(inout) :: self
        type(parameters) :: params
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        integer          :: ipart

        call init_qsys_params(params, self%nparts, self%ncunits)
        call cleanup_result_files(self%nparts)
        call build_job_descr(job_descr, self%n, self%nparts, self%ncunits)
        call qenv%new(params, self%nparts, nptcls=self%nparts, exec_bin=string('simple_test_exec'), &
            &qsys_name=string('coarray'), qsys_nthr=1)
        write(logfhandle,'(a,i0)') '>>> COARRAY MATRIX SIZE N: ', self%n
        write(logfhandle,'(a,i0)') '>>> COARRAY PARTITIONS:    ', self%nparts
        write(logfhandle,'(a,i0)') '>>> COARRAY NCUNITS:       ', self%ncunits
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=.false., extra_params=params)

        do ipart = 1,self%nparts
            call accumulate_result(self, ipart)
        end do
        self%total_sum = self%n * self%n * self%element_sum
        write(logfhandle,'(a,i0)') '>>> MATRIX ELEMENT SUM:    ', self%element_sum
        write(logfhandle,'(a,i0)') '>>> TOTAL MATRIX SUM:      ', self%total_sum
        if( self%element_sum < self%nparts ) THROW_HARD('coarray matrix sum contains unexpected values')
        call qenv%kill
        call job_descr%kill
    end subroutine master_execute

    subroutine master_finalize_run( self )
        class(coarray_matrix_sum_master_strategy), intent(inout) :: self
    end subroutine master_finalize_run

    subroutine master_cleanup( self )
        class(coarray_matrix_sum_master_strategy), intent(inout) :: self
    end subroutine master_cleanup

    function master_end_message( self ) result( msg )
        class(coarray_matrix_sum_master_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_TEST_COARRAY_MATRIX_SUM NORMAL STOP ****'
    end function master_end_message

    subroutine worker_initialize( self, cline )
        class(coarray_matrix_sum_worker_strategy), intent(inout) :: self
        class(cmdline),                            intent(inout) :: cline
        self%n      = DEFAULT_MATRIX_N
        self%part   = 1
        self%nparts = 1
        if( cline%defined('box')    ) self%n      = cline%get_iarg('box')
        if( cline%defined('part')   ) self%part   = cline%get_iarg('part')
        if( cline%defined('nparts') ) self%nparts = cline%get_iarg('nparts')
        if( self%n < 1 )      THROW_HARD('box must be > 0 for coarray matrix sum worker')
        if( self%part < 1 )   THROW_HARD('part must be > 0 for coarray matrix sum worker')
        if( self%nparts < 1 ) THROW_HARD('nparts must be > 0 for coarray matrix sum worker')
        self%image   = read_env_int('SIMPLE_COARRAY_IMAGE')
        self%nimages = read_env_int('SIMPLE_COARRAY_IMAGES')
        if( self%image < 1 )   THROW_HARD('SIMPLE_COARRAY_IMAGE is not defined for coarray matrix sum worker')
        if( self%nimages < 1 ) THROW_HARD('SIMPLE_COARRAY_IMAGES is not defined for coarray matrix sum worker')
    end subroutine worker_initialize

    subroutine worker_execute( self )
        class(coarray_matrix_sum_worker_strategy), intent(inout) :: self
        integer, allocatable :: matrix(:,:)
        integer              :: ios, funit
        character(len=:), allocatable :: fname
        allocate(matrix(self%n,self%n), source=self%image)
        self%part_sum = sum(matrix)
        fname = result_fname(self%part)
        open(newunit=funit, file=fname, status='replace', action='write', iostat=ios)
        if( ios /= 0 ) THROW_HARD('could not write '//fname)
        write(funit,'(i0,1x,i0,1x,i0,1x,i0,1x,i0)') self%part, self%image, self%nimages, self%n, self%part_sum
        close(funit)
        call mark_worker_finished(self)
        deallocate(matrix)
    end subroutine worker_execute

    subroutine worker_finalize_run( self )
        class(coarray_matrix_sum_worker_strategy), intent(inout) :: self
    end subroutine worker_finalize_run

    subroutine worker_cleanup( self )
        class(coarray_matrix_sum_worker_strategy), intent(inout) :: self
    end subroutine worker_cleanup

    function worker_end_message( self ) result( msg )
        class(coarray_matrix_sum_worker_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_TEST_COARRAY_MATRIX_SUM_WORKER NORMAL STOP ****'
    end function worker_end_message

    subroutine init_qsys_params( params, nparts, ncunits )
        type(parameters), intent(inout) :: params
        integer,          intent(in)    :: nparts, ncunits
        call simple_getcwd(params%cwd)
        if( .not. allocated(CWD_GLOB_ORIG) ) allocate(CWD_GLOB_ORIG, source=params%cwd%to_char())
        CWD_GLOB = params%cwd%to_char()
        params%prg        = 'coarray_matrix_sum'
        params%projfile   = ''
        params%qsys_name  = 'coarray'
        params%split_mode = 'even'
        params%nparts     = nparts
        params%nptcls     = nparts
        params%ncunits    = ncunits
        params%nthr       = 1
        params%numlen     = len(int2str(nparts))
    end subroutine init_qsys_params

    subroutine build_job_descr( job_descr, n, nparts, ncunits )
        type(chash), intent(inout) :: job_descr
        integer,     intent(in)    :: n, nparts, ncunits
        call job_descr%new(JOB_DESCR_SIZE)
        call job_descr%set('test',      'coarray_matrix_sum')
        call job_descr%set('box',       int2str(n))
        call job_descr%set('nparts',    int2str(nparts))
        call job_descr%set('ncunits',   int2str(ncunits))
        call job_descr%set('qsys_name', 'coarray')
        call job_descr%move_key_to_front('test')
    end subroutine build_job_descr

    subroutine cleanup_result_files( nparts )
        integer, intent(in) :: nparts
        integer :: ipart
        do ipart = 1,nparts
            call del_file(result_fname(ipart))
        end do
    end subroutine cleanup_result_files

    subroutine accumulate_result( self, ipart )
        class(coarray_matrix_sum_master_strategy), intent(inout) :: self
        integer,                                   intent(in)    :: ipart
        character(len=:), allocatable :: fname
        integer :: funit, ios, part_read, image_read, nimages_read, n_read, part_sum_read
        fname = result_fname(ipart)
        if( .not. file_exists(fname) ) THROW_HARD('missing coarray matrix sum result file: '//fname)
        open(newunit=funit, file=fname, status='old', action='read', iostat=ios)
        if( ios /= 0 ) THROW_HARD('could not read '//fname)
        read(funit,*,iostat=ios) part_read, image_read, nimages_read, n_read, part_sum_read
        close(funit)
        if( ios /= 0 ) THROW_HARD('could not parse '//fname)
        if( part_read /= ipart ) THROW_HARD('coarray matrix sum part mismatch in '//fname)
        if( n_read /= self%n ) THROW_HARD('coarray matrix sum size mismatch in '//fname)
        if( image_read < 1 .or. nimages_read < 1 ) THROW_HARD('coarray image metadata mismatch in '//fname)
        if( part_sum_read /= self%n * self%n * image_read ) THROW_HARD('coarray matrix part sum mismatch in '//fname)
        self%element_sum = self%element_sum + image_read
    end subroutine accumulate_result

    subroutine mark_worker_finished( self )
        class(coarray_matrix_sum_worker_strategy), intent(in) :: self
        type(parameters) :: params
        params%nparts         = self%nparts
        params%part           = self%part
        params%numlen         = len(int2str(self%nparts))
        params%l_distr_worker = .true.
        call qsys_job_finished(params, string('simple_coarray_matrix_sum_strategy :: worker_execute'))
    end subroutine mark_worker_finished

    function read_env_int( key ) result( val )
        character(len=*), intent(in) :: key
        integer :: val
        character(len=STDLEN) :: env
        integer :: envlen, envstat, ios
        val = 0
        call get_environment_variable(key, value=env, length=envlen, status=envstat)
        if( envstat /= 0 .or. envlen < 1 ) return
        read(env(:envlen),*,iostat=ios) val
        if( ios /= 0 ) val = 0
    end function read_env_int

    function result_fname( ipart ) result( fname )
        integer, intent(in) :: ipart
        character(len=:), allocatable :: fname
        fname = 'coarray_matrix_sum_part_'//int2str(ipart)//'.dat'
    end function result_fname

end module simple_coarray_matrix_sum_strategy
