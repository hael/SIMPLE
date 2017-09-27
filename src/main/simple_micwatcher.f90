! particles watcher for stream processing

module simple_micwatcher
#include "simple_lib.f08"
use simple_defs
use simple_syslib
use simple_fileio
use simple_params, only: params
use simple_timer
implicit none

public :: micwatcher
private

type micwatcher
    private
    character(len=STDLEN), allocatable :: history(:)             !< history of micrographs detected
    character(len=STDLEN), allocatable :: new_mics(:)             !< new micrographs
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
    procedure, private :: stack_from_mic
    procedure, private :: deftab_from_mic
    ! GETTERS
    procedure          :: get_new_mics
    procedure          :: get_new_stacks
    procedure          :: get_new_deftabs
    ! DESTRUCTOR
    procedure          :: kill
end type

interface micwatcher
    module procedure constructor
end interface micwatcher

character(len=STDLEN), parameter   :: FILETABNAME = 'mics_stream.txt'
integer,               parameter   :: FAIL_THRESH = 50
integer,               parameter   :: FAIL_TIME   = 7200 ! 2 hours

contains
    
    !>  \brief  is a constructor
    function constructor( p, report_time, print )result( self )
        class(params),     intent(in) :: p
        integer,           intent(in) :: report_time  ! in seconds
        logical, optional, intent(in) :: print
        type(micwatcher)            :: self
        character(len=STDLEN)         :: cwd
        call self%kill
        if( .not. file_exists(trim(adjustl(p%dir_mics))) )then
            print *, 'Directory does not exist: ', trim(adjustl(p%dir_mics))
            stop
        endif
        call simple_getcwd(cwd)
        self%cwd         = trim(cwd)
        self%watch_dir   = trim(adjustl(p%dir_mics))
        self%report_time = report_time
        self%ext         = p%ext
        self%mic_ext     = trim('_intg')//trim(self%ext)
        if( present(print) )then
            self%doprint = print
        else
            self%doprint = .false.
        endif
    end function constructor

    !>  \brief  is the watching procedure
    subroutine watch( self, n_mics )
#ifdef PGI
        include 'lib3f.h'
        procedure :: cast_time_char => ctime
#endif
        class(micwatcher),   intent(inout) :: self
        integer,             intent(out)   :: n_mics
        character(len=STDLEN), allocatable :: farray(:)
        integer,               allocatable :: fileinfo(:)
        logical,               allocatable :: is_new_mic(:)
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
        n_mics   = 0
        fail_cnt = 0
        ! builds mics array
        fbody = trim(self%watch_dir) // trim('/')
        call ls_filetab(trim(fbody), trim(self%mic_ext), trim(FILETABNAME))
        call read_filetable( trim(FILETABNAME), farray )
        n_lsfiles = size(farray)
        if( n_lsfiles .eq. 0 )then
            ! nothing to report
            if(allocated(self%history))deallocate(self%history)
            return
        endif
        ! identifies closed and untouched files
        allocate(is_new_mic(n_lsfiles), source=.false.)
        do i = 1, n_lsfiles
            fname = trim(adjustl(farray(i)))
            if( self%is_past(fname) )then
                is_new_mic(i) = .false.
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
                    &.and. is_closed ) is_new_mic(i) = .true.
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
            if( is_new_mic(i) )then
                is_new_mic(i) = self%to_process( farray(i) )
                if( is_new_mic(i) )write(*,'(A,A,A,A)')'>>> NEW MICROGRAPH: ',&
                &trim(adjustl(farray(i))), '; ', cast_time_char(tnow)
            endif
        enddo
        ! report
        n_mics = count(is_new_mic)
        if( n_mics > 0 )then
            if(allocated(self%new_mics))deallocate(self%new_mics)
            allocate(self%new_mics(n_mics))
            cnt = 0
            do i = 1, n_lsfiles
                if( .not.is_new_mic(i) )cycle
                cnt   = cnt + 1
                fname = trim(adjustl(farray(i)))
                call self%add2history( fname )
                self%new_mics(cnt) = trim(fname)
            enddo
        endif
    end subroutine watch

    !>  \brief whether the movie should be processed or not
    !!         if one file is missing it is re-processed  
    logical function to_process( self, fname )
        class(micwatcher), intent(inout) :: self
        character(len=*),    intent(in)  :: fname
        character(len=STDLEN) :: ctf_name, picker_name
        integer :: n
        logical :: stack_done, ctf_done
        to_process  = .false.
        ctf_name    = self%deftab_from_mic(fname)
        picker_name = self%stack_from_mic(fname)
        print *, trim(ctf_name), file_exists(ctf_name)
        print *, trim(picker_name), file_exists(picker_name)
        if( .not.file_exists(trim(ctf_name)) )return
        if( .not.file_exists(trim(picker_name)) )return
        n = nlines(trim(ctf_name))
        to_process = .true.
    end function to_process

    !>  \brief  is for adding to the history of already reported files
    subroutine add2history( self, fname )
        class(micwatcher),   intent(inout) :: self
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
        class(micwatcher), intent(inout) :: self
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

    character(len=STDLEN) function stack_from_mic( self, fname )
        class(micwatcher), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        character(len=STDLEN) :: fname_here, fbody
        fname_here = remove_abspath(trim(adjustl(fname)))
        fbody      = get_fbody(trim(fname_here), trim(self%mic_ext), separator=.false.)
        stack_from_mic = trim(self%watch_dir) // trim('/particles/ptcls_from_') // trim(fbody) // trim(self%ext)
    end function stack_from_mic

    character(len=STDLEN) function deftab_from_mic( self, fname )
        class(micwatcher), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        character(len=STDLEN) :: fname_here, fbody
        fname_here = remove_abspath(trim(adjustl(fname)))
        fbody      = get_fbody(trim(fname_here), trim(self%mic_ext), separator=.false.)
        deftab_from_mic = trim(self%watch_dir) // trim('/particles/extract_params_') // trim(fbody) // trim('.txt')
    end function deftab_from_mic

    ! GETTERS

    subroutine get_new_mics( self, mics )
        class(micwatcher),             intent(inout) :: self
        character(len=*), allocatable, intent(out)   :: mics(:)
        if( allocated(mics) )deallocate(mics)
        if( allocated(self%new_mics) )then
            mics = self%new_mics
        endif
    end subroutine get_new_mics

    subroutine get_new_deftabs( self, deftabs )
        class(micwatcher),             intent(inout) :: self
        character(len=*), allocatable, intent(out)   :: deftabs(:)
        integer :: i, n
        if( allocated(deftabs) )deallocate(deftabs)
        if( allocated(self%new_mics) )then
            n = size(self%new_mics)
            allocate(deftabs(n))
            do i = 1, n
                deftabs(i) = self%deftab_from_mic( self%new_mics(i) )
            enddo
        endif
    end subroutine get_new_deftabs

    subroutine get_new_stacks( self, stacks )
        class(micwatcher),             intent(inout) :: self
        character(len=*), allocatable, intent(out)   :: stacks(:)
        integer :: i, n
        if( allocated(stacks) )deallocate(stacks)
        if( allocated(self%new_mics) )then
            n = size(self%new_mics)
            allocate(stacks(n))
            do i = 1, n
                stacks(i) = self%stack_from_mic( self%new_mics(i) )
            enddo
        endif
    end subroutine get_new_stacks

    !>  \brief  is a destructor
    subroutine kill( self )
        class(micwatcher), intent(inout) :: self
        self%cwd        = ''
        self%watch_dir  = ''
        if(allocated(self%history))deallocate(self%history)
        self%report_time    = 0
        self%starttime      = 0
        self%ellapsedtime   = 0
        self%lastreporttime = 0
        self%n_watch        = 0
        self%n_history      = 0
        self%doprint        = .false.
    end subroutine kill

end module simple_micwatcher
