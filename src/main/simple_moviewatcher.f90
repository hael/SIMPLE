! movie watcher for stream processing
module simple_moviewatcher
#include "simple_lib.f08"
use simple_params,        only: params
use simple_timer
implicit none

public :: moviewatcher
private

type moviewatcher
    private
    character(len=STDLEN), allocatable :: history(:)             !< history of movies detected
    character(len=STDLEN)              :: cwd            = ''    !< CWD
    character(len=STDLEN)              :: watch_dir      = ''    !< movies directory to watch
    character(len=STDLEN)              :: ext            = ''    !< target directory
    character(len=STDLEN)              :: fbody          = ''    !< template name
    integer                            :: n_history      = 0     !< history of movies detected
    integer                            :: report_time    = 600   !< time ellapsed prior to processing
    integer                            :: starttime      = 0     !< time of first watch
    integer                            :: ellapsedtime   = 0     !< time ellapsed between last and first watch
    integer                            :: lastreporttime = 0     !< time ellapsed between last and first watch
    integer                            :: n_watch        = 0     !< number of times the folder has been watched
    logical                            :: dopick         = .false.
    logical                            :: doprint        = .false.
contains

    ! doers
    procedure          :: watch
    procedure, private :: is_past
    procedure, private :: add2history
    procedure, private :: to_process
    ! destructor
    procedure          :: kill
end type

interface moviewatcher
    module procedure constructor
end interface moviewatcher

character(len=STDLEN), parameter   :: FILETABNAME = 'movieftab_preprocess_stream.txt'
integer,               parameter   :: FAIL_THRESH = 50
integer,               parameter   :: FAIL_TIME   = 7200 ! 2 hours

contains

    !>  \brief  is a constructor
    function constructor( p, report_time, print )result( self )
        class(params),     intent(in) :: p
        integer,           intent(in) :: report_time  ! in seconds
        logical, optional, intent(in) :: print
        type(moviewatcher)            :: self
        character(len=STDLEN)         :: cwd
        call self%kill
        if( .not. file_exists(trim(adjustl(p%dir_movies))) )then
            print *, 'Directory does not exist: ', trim(adjustl(p%dir_movies))
            stop
        endif
        call simple_getcwd(cwd)
        self%cwd         = trim(cwd)
        self%watch_dir   = trim(adjustl(p%dir_movies))
        self%report_time = report_time
        self%dopick      = p%l_pick
        self%ext         = trim(adjustl(p%ext))
        self%fbody       = trim(adjustl(p%fbody))
        if( present(print) )then
            self%doprint = print
        else
            self%doprint = .false.
        endif
    end function constructor

    !>  \brief  is the watching procedure
    subroutine watch( self, n_movies, movies )
         use simple_timer
        class(moviewatcher),           intent(inout) :: self
        integer,                       intent(out)   :: n_movies
        character(len=*), allocatable, intent(out)   :: movies(:)
        character(len=STDLEN), allocatable :: farray(:)
        integer,               allocatable :: fileinfo(:)
        logical,               allocatable :: is_new_movie(:)
        integer               :: tnow, last_accessed, last_modified, last_status_change ! in seconds
        integer               :: i, io_stat, n_lsfiles, cnt, fail_cnt
        character(len=STDLEN) :: fname, abs_fname
        logical               :: is_closed
        integer(timer_int_kind) :: twatch
        ! init
        self%n_watch = self%n_watch + 1
        tnow = simple_gettime()
        if( self%n_watch .eq. 1 )then
            ! first call
            self%starttime  = tnow
        endif
        self%ellapsedtime = tnow - self%starttime
        if(allocated(movies))deallocate(movies)
        n_movies = 0
        fail_cnt = 0
        ! builds files array
        call ls_mrcfiletab(self%watch_dir, trim(FILETABNAME))
        call read_filetable( trim(FILETABNAME), farray )
        n_lsfiles = size(farray)
        if( n_lsfiles .eq. 0 )then
            ! nothing to report
            if(allocated(self%history))deallocate(self%history)
            return
        endif
        ! identifies closed and untouched files
        allocate(is_new_movie(n_lsfiles), source=.false.)
        do i = 1, n_lsfiles
            fname = trim(adjustl(farray(i)))
            if( self%is_past(fname) )then
                is_new_movie(i) = .false.
            else
                if( fname(1:1).eq."/" )then
                    ! already has absolute path
                    abs_fname = trim(adjustl(fname))
                else
                    ! generate absolute path
                    abs_fname = trim(self%cwd)//'/'//trim(adjustl(fname))
                endif
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
                    &.and. is_closed ) is_new_movie(i) = .true.
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
            if( is_new_movie(i) )then
                is_new_movie(i) = self%to_process( farray(i) )
                if( is_new_movie(i) )write(*,'(A,A,A,A)')'>>> NEW MOVIE: ',&
                &trim(adjustl(farray(i))), '; ', cast_time_char(tnow)
            endif
        enddo
        ! report
        n_movies = count(is_new_movie)
        if( n_movies > 0 )then
            allocate(movies(n_movies))
            cnt = 0
            do i = 1, n_lsfiles
                if( .not.is_new_movie(i) )cycle
                cnt   = cnt + 1
                fname = trim(adjustl(farray(i)))
                call self%add2history( fname )
                movies(cnt) = trim(fname)
            enddo
        endif
    end subroutine watch

    !>  \brief whether the movie should be processed or not
    !!         if one file is missing it is re-processed
    logical function to_process( self, fname )
        class(moviewatcher), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        character(len=STDLEN) :: fname_here, fbody, ext, motion_correct_name, unidoc_name, picker_name
        logical :: motion_correct_done, ctffind_done, picker_done
        fname_here = remove_abspath(trim(adjustl(fname)))
        ext        = trim(fname2ext(fname_here))
        fbody      = get_fbody(trim(fname_here), trim(ext))
        if( trim(self%fbody) .ne. '' )fbody = trim(self%fbody)//trim(adjustl(fbody))
        ! motion_correct
        motion_correct_name = motion_correct_STREAM_DIR//trim(adjustl(fbody))//THUMBNAIL_SUFFIX//trim(self%ext)
        motion_correct_done = file_exists(trim(motion_correct_name))
        if( motion_correct_done )then
            motion_correct_name = motion_correct_STREAM_DIR//trim(adjustl(fbody))//INTGMOV_SUFFIX//trim(self%ext)
            motion_correct_done = file_exists(trim(motion_correct_name))
        endif
        ! ctffind
        unidoc_name  = UNIDOC_STREAM_DIR//UNIDOC_OUTPUT//trim(adjustl(fbody))//'.txt'
        ctffind_done = file_exists(trim(unidoc_name))
        to_process = .not. (motion_correct_done .and. ctffind_done)
    end function to_process

    !>  \brief  is for adding to the history of already reported files
    subroutine add2history( self, fname )
        class(moviewatcher), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        character(len=STDLEN), allocatable :: tmp_farr(:)
        integer :: n
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
        class(moviewatcher), intent(inout) :: self
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

    !>  \brief  is a destructor
    subroutine kill( self )
        class(moviewatcher), intent(inout) :: self
        self%cwd        = ''
        self%watch_dir  = ''
        self%ext        = ''
        if(allocated(self%history))deallocate(self%history)
        self%report_time    = 0
        self%starttime      = 0
        self%ellapsedtime   = 0
        self%lastreporttime = 0
        self%n_watch        = 0
        self%n_history      = 0
        self%doprint        = .false.
    end subroutine kill

end module simple_moviewatcher
