!@descr: the numbered sets a stream stage runs as streaming jobs: naming, writing, queueing, collecting, completing and restoring them
!==============================================================================
! MODULE: simple_stream_job_sets
!
! PURPOSE:
!   A stream stage that works in batches (preprocessing: five movies; reference
!   picking: the accepted micrographs of one upstream project) writes each
!   batch as a small project, a "set", in its job folder, runs one streaming
!   job on it, and moves it to the completed folder once the job is done. This
!   type keeps the bookkeeping the stages used to repeat:
!     - sets are numbered 1, 2, ... and named <number padded to numlen>.simple
!     - a set can record where it came from (an upstream project), so that a
!       restarted stage knows which inputs it has already handled
!     - on restart the numbering continues after the highest set of the
!       previous run, completed or left unfinished, whether or not anything in
!       it was accepted; the job folder with the unfinished sets is set aside
!       (a job left running reads its set by path), and their inputs are
!       submitted again
!     - a stopping stage cancels the jobs in flight (cancel)
!   What goes into a set, what a stage does with a finished one, and whether
!   a finished set is moved to the completed folder stay in the stage.
!
! LIFECYCLE:
!   new(job_dir, completed_dir, numlen) -> [restore] ->
!   { write_set -> submit ; schedule ; collect -> complete } -> kill
!
!   It could move beside qsys_async_job in src/utils/qsys, which would then
!   depend on sp_project.
!==============================================================================
module simple_stream_job_sets
use simple_defs_fname,   only: METADATA_EXT, STDERROUT_DIR
use simple_string,       only: string
use simple_string_utils, only: int2str_pad, str2int
use simple_fileio,       only: basename, get_fbody, simple_abspath, simple_list_files_regexp, simple_rename, simple_rmdir
use simple_syslib,       only: simple_mkdir
use simple_cmdline,      only: cmdline
use simple_sp_project,   only: sp_project
use simple_qsys_env,     only: qsys_env
use simple_qsys_job_record, only: set_aside_dir
implicit none

public :: stream_job_sets
private

character(len=*), parameter :: SOURCE_KEY = 'stream_set_source' ! projinfo key of a set's origin

type :: stream_job_sets
    private
    type(string) :: job_dir             ! absolute; sets are written and their jobs run here
    type(string) :: completed_dir       ! absolute; finished sets are moved here
    integer      :: numlen   = 5
    integer      :: counter  = 0        ! number of the last set
    logical      :: l_exists = .false.
contains
    procedure :: new
    procedure :: write_set
    procedure :: submit
    procedure :: schedule
    procedure :: cancel
    procedure :: collect
    procedure :: complete
    procedure :: restore
    procedure :: get_counter
    procedure :: get_job_dir
    procedure :: get_completed_dir
    procedure :: kill
end type stream_job_sets

contains

    !> Makes the job folder (with the folder the job scripts write their output to) and the
    !! completed folder; set names are padded to @p numlen digits.
    subroutine new( self, job_dir, completed_dir, numlen )
        class(stream_job_sets), intent(inout) :: self
        class(string),          intent(in)    :: job_dir, completed_dir
        integer,                intent(in)    :: numlen
        call self%kill()
        call make_job_dir(job_dir)
        call simple_mkdir(completed_dir)
        self%job_dir       = simple_abspath(job_dir)
        self%completed_dir = simple_abspath(completed_dir)
        self%numlen        = numlen
        self%counter       = 0
        self%l_exists      = .true.
    end subroutine new

    !> Writes @p set_proj as the next set in the job folder and points @p cline_worker at it:
    !! projname, projfile, and items 1 to @p nitems. @p source, when given, is recorded in the
    !! set as its origin and returned by restore().
    subroutine write_set( self, set_proj, cline_worker, nitems, source )
        class(stream_job_sets),  intent(inout) :: self
        class(sp_project),       intent(inout) :: set_proj
        class(cmdline),          intent(inout) :: cline_worker
        integer,                 intent(in)    :: nitems
        class(string), optional, intent(in)    :: source
        type(string) :: projname, projfile
        self%counter = self%counter + 1
        projname     = int2str_pad(self%counter, self%numlen)
        projfile     = projname//METADATA_EXT
        call set_proj%projinfo%new(1, is_ptcl=.false.)
        call set_proj%projinfo%set(1, 'projname', projname)
        call set_proj%projinfo%set(1, 'projfile', projfile)
        call set_proj%projinfo%set(1, 'cwd',      self%job_dir)
        if( present(source) ) call set_proj%projinfo%set(1, SOURCE_KEY, source)
        call set_proj%write(self%job_dir//'/'//projfile)
        call cline_worker%set('projname', projname)
        call cline_worker%set('projfile', projfile)
        call cline_worker%set('fromp',    1)
        call cline_worker%set('top',      nitems)
    end subroutine write_set

    !> Queues one streaming job with @p cline_worker (as write_set left it).
    subroutine submit( self, qenv, cline_worker )
        class(stream_job_sets), intent(inout) :: self
        class(qsys_env),        intent(inout) :: qenv
        class(cmdline),         intent(in)    :: cline_worker
        call qenv%qscripts%add_to_streaming(cline_worker)
    end subroutine submit

    !> Cancels the jobs in flight and drops the queued ones (a stopping stage).
    subroutine cancel( self, qenv )
        class(stream_job_sets), intent(inout) :: self
        class(qsys_env),        intent(inout) :: qenv
        call qenv%qscripts%cancel_streaming(self%job_dir)
    end subroutine cancel

    !> Starts queued jobs on the computing units that are free; the jobs run in the job folder.
    subroutine schedule( self, qenv )
        class(stream_job_sets), intent(inout) :: self
        class(qsys_env),        intent(inout) :: qenv
        call qenv%qscripts%schedule_streaming(qenv%qdescr, path=self%job_dir)
    end subroutine schedule

    !> The sets whose jobs have finished since the last call (absolute paths, in the job folder),
    !! and the items of the jobs that failed (the items each set's job was given, fromp to top:
    !! the movies of a preprocessing set, partial or not).
    subroutine collect( self, qenv, done, nfailed_items )
        class(stream_job_sets),    intent(inout) :: self
        class(qsys_env),           intent(inout) :: qenv
        type(string), allocatable, intent(inout) :: done(:)
        integer,                   intent(out)   :: nfailed_items
        class(cmdline), allocatable :: done_clines(:), failed_clines(:)
        integer :: i, n
        if( allocated(done) ) deallocate(done)
        allocate(done(0))
        nfailed_items = 0
        if( qenv%qscripts%get_done_stacksz() > 0 )then
            call qenv%qscripts%get_stream_done_stack(done_clines)
            n = size(done_clines)
            deallocate(done)
            allocate(done(n))
            do i = 1,n
                done(i) = self%job_dir//'/'//done_clines(i)%get_carg('projfile')
                call done_clines(i)%kill
            enddo
            deallocate(done_clines)
        endif
        if( qenv%qscripts%get_failed_stacksz() > 0 )then
            call qenv%qscripts%get_stream_fail_stack(failed_clines, n)
            if( n > 0 )then
                do i = 1,n
                    nfailed_items = nfailed_items + max(1, failed_clines(i)%get_iarg('top') - failed_clines(i)%get_iarg('fromp') + 1)
                    call failed_clines(i)%kill
                enddo
                deallocate(failed_clines)
            endif
        endif
    end subroutine collect

    !> Moves the finished set @p fname to the completed folder; @p completed_fname is its new path.
    subroutine complete( self, fname, completed_fname )
        class(stream_job_sets), intent(in)    :: self
        class(string),          intent(in)    :: fname
        type(string),           intent(inout) :: completed_fname
        completed_fname = self%completed_dir//'/'//basename(fname)
        call simple_rename(fname, completed_fname)
    end subroutine complete

    !> Restart: the completed sets of the previous run (absolute paths) and, with @p sources,
    !! each one's recorded origin ('' when none). Numbering continues after the highest set,
    !! completed or left unfinished; a job folder holding unfinished sets is set aside, since a
    !! job left running reads its set by path, and the job folder is made afresh.
    subroutine restore( self, completed, sources )
        class(stream_job_sets),              intent(inout) :: self
        type(string), allocatable,           intent(inout) :: completed(:)
        type(string), allocatable, optional, intent(inout) :: sources(:)
        type(string),     allocatable :: unfinished(:)
        type(sp_project) :: proj
        integer          :: i
        logical          :: l_unfinished
        call simple_list_files_regexp(self%completed_dir, '\.simple$', completed)
        if( .not. allocated(completed) ) allocate(completed(0))
        call simple_list_files_regexp(self%job_dir, '\.simple$', unfinished)
        if( .not. allocated(unfinished) ) allocate(unfinished(0))
        self%counter = 0
        do i = 1,size(completed)
            self%counter = max(self%counter, set_number(completed(i)))
        enddo
        do i = 1,size(unfinished)
            self%counter = max(self%counter, set_number(unfinished(i)))
        enddo
        l_unfinished = size(unfinished) > 0
        if( present(sources) )then
            if( allocated(sources) ) deallocate(sources)
            allocate(sources(size(completed)))
            do i = 1,size(completed)
                sources(i) = ''
                call proj%read_segment('projinfo', completed(i))
                if( proj%projinfo%get_noris() > 0 )then
                    if( proj%projinfo%isthere(1, SOURCE_KEY) ) sources(i) = proj%projinfo%get_str(1, SOURCE_KEY)
                endif
                call proj%kill
            enddo
        endif
        if( l_unfinished )then
            call set_aside_dir(self%job_dir)
        else
            call simple_rmdir(self%job_dir)
        endif
        call make_job_dir(self%job_dir)
    end subroutine restore

    ! the number of the set file @p fname (0 when its name is not a number)
    integer function set_number( fname )
        class(string), intent(in) :: fname
        type(string) :: fbody
        integer      :: id, iostat
        set_number = 0
        fbody = basename(fname)
        fbody = get_fbody(fbody, METADATA_EXT, separator=.false.)
        id    = str2int(fbody, iostat)
        if( iostat == 0 ) set_number = id
    end function set_number

    integer function get_counter( self )
        class(stream_job_sets), intent(in) :: self
        get_counter = self%counter
    end function get_counter

    function get_job_dir( self ) result( dir )
        class(stream_job_sets), intent(in) :: self
        type(string) :: dir
        dir = self%job_dir
    end function get_job_dir

    function get_completed_dir( self ) result( dir )
        class(stream_job_sets), intent(in) :: self
        type(string) :: dir
        dir = self%completed_dir
    end function get_completed_dir

    subroutine kill( self )
        class(stream_job_sets), intent(inout) :: self
        if( .not. self%l_exists ) return
        call self%job_dir%kill
        call self%completed_dir%kill
        self%numlen   = 5
        self%counter  = 0
        self%l_exists = .false.
    end subroutine kill

    ! the job folder and the folder its job scripts write their output to
    subroutine make_job_dir( dir )
        class(string), intent(in) :: dir
        call simple_mkdir(dir)
        call simple_mkdir(dir//'/'//STDERROUT_DIR)
    end subroutine make_job_dir

end module simple_stream_job_sets
