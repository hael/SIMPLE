
module simple_moviewatcher
use simple_defs
use simple_syscalls
use simple_filehandling
use simple_jiffys
implicit none

public :: moviewatcher
private

type moviewatcher
    private
        character(len=STDLEN)              :: cwd       = ''    !< CWD
        character(len=STDLEN)              :: watch_dir = ''    !< movies directory to watch
        character(len=STDLEN), allocatable :: history(:)        !< history of movies detected
        integer                            :: n_history = 0     !< history of movies detected
        integer                            :: report_time       !< time ellapsed prior to processing
        integer                            :: starttime         !< time of first watch
        integer                            :: ellapsedtime      !< time ellapsed between last and first watch
        integer                            :: lastreporttime    !< time ellapsed between last and first watch
        integer                            :: n_watch   = 0     !< number of times the folder has been watched
        integer                            :: n_fails   = 0     !< number of times a watch has resulted in some error
        logical                            :: doprint   = .false.

  contains

    ! doers
    procedure          :: watch
    procedure, private :: is_past
    procedure, private :: add2history
    ! destructor
    procedure          :: kill

end type

interface moviewatcher
    module procedure constructor
end interface moviewatcher

character(len=STDLEN), parameter   :: FILETABNAME='movieftab_preproc_stream.txt'
integer,               parameter   :: FAIL_THRESH = 50

contains
    
    !>  \brief  is a constructor
    function constructor( workdir, report_time, print )result( self )
        character(len=*),  intent(in) :: workdir
        integer,           intent(in) :: report_time  ! in seconds
        logical, optional, intent(in) :: print
        type(moviewatcher)            :: self
        character(len=STDLEN)         :: cwd
        call self%kill
        if( .not. file_exists(trim(adjustl(workdir))) )then
            print *, 'Directory does not exist: ', trim(adjustl(workdir))
            stop
        endif
        call getcwd(cwd)
        self%cwd         = trim(cwd)
        self%watch_dir   = trim(adjustl(workdir))
        self%report_time = report_time
        if( present(print) )then
            self%doprint = print
        else
            self%doprint = .false.
        endif
    end function constructor

    !>  \brief  is the watching procedure
    subroutine watch( self, n_files, files )
        class(moviewatcher),           intent(inout) :: self
        integer,                       intent(out)   :: n_files
        character(len=*), allocatable, intent(out)   :: files(:)
        character(len=STDLEN), allocatable :: farray(:)
        integer,               allocatable :: stats(:)
        integer               :: tnow, last_accessed, last_modified, last_status_change ! in seconds
        integer               :: i, fstat, fail_cnt, alloc_stat, n_lsfiles
        character(len=STDLEN) :: fname, abs_fname
        logical               :: is_closed, has_new
        has_new      = .false.
        fail_cnt     = 0
        self%n_watch = self%n_watch + 1
        ! updates timing
        if( self%n_watch .eq. 1 )then
            ! first call
            self%starttime = time()
            tnow = self%starttime
        else
            tnow = time()
        endif
        self%ellapsedtime = tnow - self%starttime
        ! builds files array
        call sys_gen_mrcfiletab(self%watch_dir, trim(FILETABNAME))
        call read_filetable( trim(FILETABNAME), farray )
        n_lsfiles  = size(farray)
        if( n_lsfiles .eq. 0 )then
            if(allocated(self%history))deallocate(self%history)
            return
        endif
        ! identifies files to report
        do i = 1, n_lsfiles
            fname = trim(adjustl(farray(i)))
            if( self%is_past(fname) )then
                ! already been reported, nothing to do
            else
                abs_fname = trim(self%cwd)//'/'//trim(adjustl(fname))
                call file_stats(trim(abs_fname), fstat, stats)
                is_closed = .not. is_file_open(abs_fname)
                if( fstat.eq.0 )then
                    ! new file identified
                    last_accessed      = tnow - stats( 9)
                    last_modified      = tnow - stats(10)
                    last_status_change = tnow - stats(11)
                    if(    (last_accessed      > self%report_time)&
                    &.and. (last_modified      > self%report_time)&
                    &.and. (last_status_change > self%report_time)&
                    &.and. is_closed )then
                        ! new file to report
                        call self%add2history( fname )
                        self%lastreporttime = tnow
                        print *,'New_movie: ', trim(abs_fname), ', ',ctime(maxval(stats(9:11)))
                    endif
                else
                    ! some error occured
                    fail_cnt = fail_cnt + 1
                    print *,'Error watching file: ', trim(fname), ' with code: ',fstat
                endif
                deallocate(stats)
            endif
        enddo
        ! Failures handling
        if( fail_cnt > 0 )self%n_fails = self%n_fails + 1
        if( self%n_fails > FAIL_THRESH )then
            ! DO SOMETHING
            print *, 'Warning: unknown failure occured, simple_moviewatcher%watch'
        endif
        ! report
        if(allocated(files))deallocate(files)
        n_files = self%n_history
        if( n_files > 0 )then
            allocate(files(n_files), stat=alloc_stat)
            call alloc_err("In: simple_moviewatcher%watch", alloc_stat)
            files = self%history
        endif
        ! cleanup
        deallocate(farray)
    end subroutine watch

    !>  \brief  is for adding to the history of already reported files
    subroutine add2history( self, fname )
        class(moviewatcher), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        character(len=STDLEN), allocatable :: tmp_farr(:)
        integer :: n, alloc_stat
        if( .not.allocated(self%history) )then
            n = 0
            allocate(self%history(1), stat=alloc_stat)
            call alloc_err("In: simple_moviewatcher%add2history 2", alloc_stat)
        else
            n = size(self%history)
            allocate(tmp_farr(n), stat=alloc_stat)
            call alloc_err("In: simple_moviewatcher%add2history 2", alloc_stat)
            tmp_farr(:) = self%history
            deallocate(self%history)
            allocate(self%history(n+1), stat=alloc_stat)
            call alloc_err("In: simple_moviewatcher%add2history 3", alloc_stat)
            self%history(:n) = tmp_farr
            deallocate(tmp_farr)
        endif
        self%history(n+1) = trim(adjustl(fname))
        self%n_history    = self%n_history + 1
    end subroutine add2history

    !>  \brief  is for checking a file has already been reported
    logical function is_past( self, fname )
        class(moviewatcher), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        integer :: i
        is_past = .false.
        if( allocated(self%history) )then
            do i = 1, size(self%history)
                if( trim(adjustl(fname)) .eq. trim(adjustl(self%history(i))) )then
                    is_past = .true.
                    return
                endif
            enddo
        else
            return
        endif
    end function is_past

    !>  \brief  is a destructor
    subroutine kill( self )
        class(moviewatcher), intent(inout) :: self
        self%cwd       = ''
        self%watch_dir = ''
        if(allocated(self%history))deallocate(self%history)
        self%report_time    = 0
        self%starttime      = 0
        self%ellapsedtime   = 0
        self%lastreporttime = 0
        self%n_watch        = 0
        self%n_fails        = 0
        self%n_history      = 0
        self%doprint        = .false.
    end subroutine kill

end module simple_moviewatcher
