! particles watcher for stream processing

module simple_extractwatcher
#include "simple_lib.f08"
use simple_defs
use simple_syslib
use simple_fileio
use simple_params, only: params
use simple_timer
implicit none

public :: extractwatcher
private

type extractwatcher
    private
    character(len=STDLEN), allocatable :: history(:)             !< history of micrographs detected
    character(len=STDLEN), allocatable :: new_stks(:)            !< new extracted stacks
    character(len=STDLEN)              :: cwd            = ''    !< CWD
    character(len=STDLEN)              :: watch_dir      = ''    !< movies directory to watch
    character(len=4)                   :: ext            = ''    !< 
    character(len=STDLEN)              :: mic_ext        = ''    !< 
    integer                            :: n_history      = 0     !< history of movies detected
    integer                            :: report_time    = 10    !< time ellapsed prior to processing
    integer                            :: starttime      = 0     !< time of first watch
    integer                            :: ellapsedtime   = 0     !< time ellapsed between last and first watch
    integer                            :: lastreporttime = 0     !< time ellapsed between last and first watch
    integer                            :: n_watch        = 0     !< number of times the folder has been watched
    logical                            :: ctf            = .false.
    logical                            :: doprint        = .false.
contains

    ! DOERS
    procedure          :: watch
    procedure, private :: is_past
    procedure, private :: add2history
    procedure, private :: to_process
    procedure, private :: deftab_from_stk
    ! GETTERS
    procedure          :: get_new_stacks
    procedure          :: get_new_deftabs
    ! DESTRUCTOR
    procedure          :: kill
end type

interface extractwatcher
    module procedure constructor
end interface extractwatcher

character(len=STDLEN), parameter   :: FILETABNAME = 'extractwatcher.txt'
integer,               parameter   :: FAIL_THRESH = 50
integer,               parameter   :: FAIL_TIME   = 7200 ! 2 hours

contains
    
    !>  \brief  is a constructor
    function constructor( p, report_time, print )result( self )
        class(params),     intent(in) :: p
        integer,           intent(in) :: report_time  ! in seconds
        logical, optional, intent(in) :: print
        type(extractwatcher)      :: self
        character(len=STDLEN) :: cwd
        call self%kill
        call simple_getcwd(cwd)
        self%cwd       = trim(cwd)
        self%watch_dir = trim(p%dir_ptcls)//'/'
        if( .not. file_exists(trim(adjustl(self%watch_dir))) )then
            print *, 'Directory does not exist: ', trim(adjustl(self%watch_dir))
            stop
        endif
        self%report_time = report_time
        self%ext         = trim(p%ext)
        if( present(print) )then
            self%doprint = print
        else
            self%doprint = .false.
        endif
    end function constructor

    !>  \brief  is the watching procedure
    subroutine watch( self, n_stks )
#ifdef PGI
        include 'lib3f.h'
        procedure :: cast_time_char => ctime
#endif
        class(extractwatcher), intent(inout) :: self
        integer,               intent(out)   :: n_stks
        character(len=STDLEN),   allocatable :: farray(:)
        integer,                 allocatable :: fileinfo(:)
        logical,                 allocatable :: is_new_stk(:)
        integer               :: tnow, last_accessed, last_modified, last_status_change ! in seconds
        integer               :: i, io_stat, alloc_stat, n_lsfiles, cnt, fail_cnt
        character(len=STDLEN) :: fname, abs_fname, fbody, ext
        logical               :: is_closed
        ! init
        self%n_watch = self%n_watch + 1
        tnow = simple_gettime()
        if( self%n_watch .eq. 1 )then
            ! first call
            self%starttime  = tnow
        endif
        self%ellapsedtime = tnow - self%starttime
        n_stks   = 0
        fail_cnt = 0
        ! builds mics array
        fbody = trim(self%watch_dir) // trim('/ptcls_from_')
        call ls_filetab(trim(fbody), trim(self%ext), trim(FILETABNAME))
        call read_filetable( trim(FILETABNAME), farray )
        n_lsfiles = size(farray)
        if( n_lsfiles .eq. 0 )then
            ! nothing to report
            if(allocated(self%history))deallocate(self%history)
            return
        endif
        ! identifies closed and untouched files
        allocate(is_new_stk(n_lsfiles), source=.false.)
        do i = 1, n_lsfiles
            fname = trim(adjustl(farray(i)))
            if( self%is_past(fname) )then
                is_new_stk(i) = .false.
            else
                abs_fname = trim(self%cwd)//'/'//trim(adjustl(fname))
                call simple_file_stat(abs_fname, io_stat, fileinfo, doprint=.false.)
                is_closed = .not. is_file_open(abs_fname)
                if( io_stat.eq.0 )then
                    ! new movie
                    last_accessed      = tnow - fileinfo( 9)
                    last_modified      = tnow - fileinfo(10)
                    last_status_change = tnow - fileinfo(11)
                    if(    (last_accessed      > self%report_time)&
                    &.and. (last_modified      > self%report_time)&
                    &.and. (last_status_change > self%report_time)&
                    &.and. is_closed ) is_new_stk(i) = .true.
                else
                    ! some error occured
                    fail_cnt = fail_cnt + 1
                    print *,'Error watching file: ', trim(fname), ' with code: ',io_stat
                endif
                if(allocated(fileinfo))deallocate(fileinfo)
            endif
        enddo
        ! identifies already processed files for restart
        do i = 1, n_lsfiles
            if( is_new_stk(i) )then
                is_new_stk(i) = self%to_process( farray(i) )
                if( is_new_stk(i) )write(*,'(A,A,A,A)')'>>> NEW EXTRACTED STACK: ',&
                &trim(adjustl(farray(i))), '; ', cast_time_char(tnow)
            endif
        enddo
        ! report
        n_stks = count(is_new_stk)
        if( n_stks > 0 )then
            if(allocated(self%new_stks))deallocate(self%new_stks)
            allocate(self%new_stks(n_stks))
            cnt = 0
            do i = 1, n_lsfiles
                if( .not.is_new_stk(i) )cycle
                cnt   = cnt + 1
                fname = trim(adjustl(farray(i)))
                call self%add2history( fname )
                self%new_stks(cnt) = trim(fname)
            enddo
        endif
    end subroutine watch

    !>  \brief whether the movie should be processed or not
    !!         if one file is missing it is re-processed  
    logical function to_process( self, fname )
        class(extractwatcher), intent(inout) :: self
        character(len=*),    intent(in)  :: fname
        character(len=STDLEN) :: ctf_name
        integer :: nl, ldim(3), nptcls
        logical :: stack_done, ctf_done
        to_process  = .false.
        ctf_name    = self%deftab_from_stk(fname)
        if( .not.file_exists(trim(ctf_name)) )return
        nl = nlines(trim(ctf_name))
        call find_ldim_nptcls(trim(fname), ldim, nptcls)
        if( nptcls.ne.nl )return
        to_process = .true.
    end function to_process

    !>  \brief  is for adding to the history of already reported files
    subroutine add2history( self, fname )
        class(extractwatcher),   intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        character(len=STDLEN), allocatable :: tmp_farr(:)
        integer :: n, alloc_stat
        if( .not.allocated(self%history) )then
            n = 0
            allocate(self%history(1), stat=alloc_stat)
            allocchk("In: simple_moviewatcher%add2history 1")
        else
            n = size(self%history)
            allocate(tmp_farr(n), stat=alloc_stat)
            allocchk("In: simple_moviewatcher%add2history 2")
            tmp_farr(:) = self%history
            deallocate(self%history)
            allocate(self%history(n+1), stat=alloc_stat)
            allocchk("In: simple_moviewatcher%add2history 3")
            self%history(:n) = tmp_farr
            deallocate(tmp_farr)
        endif
        self%history(n+1) = trim(adjustl(fname))
        self%n_history    = self%n_history + 1
    end subroutine add2history

    !>  \brief  is for checking a file has already been reported
    logical function is_past( self, fname )
        class(extractwatcher), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        integer :: i
        is_past = .false.
        if( allocated(self%history) )then
            do i = 1, size(self%history)
                if( trim(adjustl(fname)) .eq. trim(adjustl(self%history(i))) )then
                    is_past = .true.
                    exit
                endif
            enddo
        endif
    end function is_past

    character(len=STDLEN) function deftab_from_stk( self, fname )
        class(extractwatcher), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        character(len=STDLEN) :: fname_here, fbody
        fname_here = remove_abspath(trim(adjustl(fname)))
        fbody      = get_fbody(trim(fname_here), trim(self%ext), separator=.false.)
        fbody      = fbody(12:) ! removes the 'ptcls_from_'
        deftab_from_stk = trim(self%watch_dir) // trim('extract_params_') // trim(fbody) // trim('.txt')
    end function deftab_from_stk

    ! GETTERS

    subroutine get_new_deftabs( self, deftabs )
        class(extractwatcher),             intent(inout) :: self
        character(len=*), allocatable, intent(out)   :: deftabs(:)
        integer :: i, n
        if( allocated(deftabs) )deallocate(deftabs)
        if( allocated(self%new_stks) )then
            n = size(self%new_stks)
            allocate(deftabs(n))
            do i = 1, n
                deftabs(i) = self%deftab_from_stk( self%new_stks(i) )
            enddo
        endif
    end subroutine get_new_deftabs

    subroutine get_new_stacks( self, stacks )
        class(extractwatcher),             intent(inout) :: self
        character(len=*), allocatable, intent(out)   :: stacks(:)
        integer :: i, n
        if( allocated(stacks) )deallocate(stacks)
        if( allocated(self%new_stks) )then
            stacks = self%new_stks
        endif
    end subroutine get_new_stacks

    !>  \brief  is a destructor
    subroutine kill( self )
        class(extractwatcher), intent(inout) :: self
        self%cwd        = ''
        self%watch_dir  = ''
        if(allocated(self%history))deallocate(self%history)
        if(allocated(self%new_stks))deallocate(self%new_stks)
        self%report_time    = 0
        self%starttime      = 0
        self%ellapsedtime   = 0
        self%lastreporttime = 0
        self%n_watch        = 0
        self%n_history      = 0
        self%doprint        = .false.
    end subroutine kill

end module simple_extractwatcher
