module simple_commander_stream_wflows
use simple_defs
use simple_cmdline,           only: cmdline
use simple_chash,             only: chash
use simple_qsys_base,         only: qsys_base
use simple_qsys_factory,      only: qsys_factory
use simple_qsys_ctrl,         only: qsys_ctrl
use simple_params,            only: params
use simple_commander_base,    only: commander_base
use simple_strings,           only: real2str
use simple_commander_preproc, only: preproc_commander
use simple_filehandling       ! use all in there
use simple_qsys_funs          ! use all in there
use simple_syscalls           ! use all in there
implicit none

public :: preproc_stream_commander
private

type, extends(commander_base) :: preproc_stream_commander
  contains
    procedure :: execute      => exec_preproc_stream
end type preproc_stream_commander

contains

    ! PRE-PROCESS SINGLE-PARTICLE DDDs IN STREAMING MODE

    subroutine exec_preproc_stream( self, cline )
        use simple_commander_preproc, only: preproc_commander
        class(preproc_stream_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! constants
        logical,               parameter   :: DEBUG = .true.
        character(len=STDLEN), parameter   :: FILETABNAME='movieftab_preproc_stream.txt'
        integer,               parameter   :: SHORTTIME = 5
        ! commanders
        type(preproc_commander)            :: xpreproc
        ! other vars
        character(len=STDLEN), allocatable :: movienames(:)
        type(params)                       :: p_master
        type(qsys_ctrl)                    :: qscripts
        type(chash)                        :: myq_descr, job_descr
        type(chash), allocatable           :: part_params(:)
        integer, allocatable               :: parts(:,:)
        logical, allocatable               :: lmask_stream(:), ltmp(:)
        type(qsys_factory)                 :: qsys_fac
        class(qsys_base), pointer          :: myqsys
        integer                            :: nmovies, nmovies_prev, i, nmovies_in_queue, imovies
        ! make master parameters
        p_master = params(cline, checkdistr=.false.)
        ! set defaults
        p_master%split_mode = 'singles'
        p_master%numlen     = 5
        if( cline%defined('fbody') )then
            call cline%set('fbody', trim(p_master%dir_target)//'/'//trim(p_master%fbody))
        else
            call cline%set('fbody', trim(p_master%dir_target)//'/')
        endif
        ! make target directory
        call exec_cmdline('mkdir -p '//trim(adjustl(p_master%dir_target)))

        

        nmovies = 0
        do
            call create_stream_filetable
            if( nmovies > nmovies_prev )then
                ! there are new movies to process
                if( DEBUG ) print *, 'nmovies updated to: ', nmovies
                ! prepare stream mask
                if( allocated(lmask_stream) )then
                    nmovies_in_queue = size(lmask_stream)
                else
                    nmovies_in_queue = 0
                endif
                if( nmovies > nmovies_in_queue )then
                    allocate(ltmp(nmovies))
                    ltmp = .true.
                    ltmp(1:nmovies_in_queue) = lmask_stream(1:nmovies_in_queue)
                    if( allocated(lmask_stream) ) deallocate(lmask_stream)
                    allocate(lmask_stream(1:nmovies), source=ltmp)
                    deallocate(ltmp)
                    if( DEBUG )then
                        do i=1,nmovies
                            print *, i, 'lmask: ', lmask_stream(i)
                        end do
                    endif
                    call create_individual_filetables
                    call setup_distr_env
                    call qscripts%schedule_jobs
                else
                    write(*,'(a)') 'WEIRD! nmovies <= nmovies_in_queue, despite that&
                    &there are new movies in filetab'
                endif
            endif
            call sleep(SHORTTIME)
        end do

        contains

            subroutine create_stream_filetable
                ! generate the filetable via system call
                call sys_gen_mrcfiletab(p_master%dir_movies, FILETABNAME)
                nmovies_prev = nmovies
                nmovies      = nlines(FILETABNAME)
                ! read the movienames back in
                call read_filetable(FILETABNAME, movienames)
                ! sanity check
                if( size(movienames) /= nmovies ) print *, 'WEIRD! nlines of streaming filetable .ne. size of arr'
                ! need to set p_master%nptcls = nmovies
                p_master%nptcls = nmovies
            end subroutine create_stream_filetable

            subroutine create_individual_filetables
                integer :: imovie, fnr, file_stat
                character(len=STDLEN), allocatable :: individual_filetabs(:)
                individual_filetabs = make_filenames('movie_to_process', nmovies, '.txt', p_master%numlen)
                if( allocated(part_params) ) deallocate(part_params)
                allocate(part_params(size(individual_filetabs)))
                do imovie=1,nmovies
                    call part_params(imovie)%new(1)
                    call part_params(imovie)%set('filetab', trim(individual_filetabs(imovie)))
                    if( lmask_stream(imovie) )then
                        fnr = get_fileunit()
                        open(unit=fnr, status='replace', action='write', file=trim(individual_filetabs(imovie)), iostat=file_stat)
                        call fopen_err('exec_preproc_stream :: create_individual_filetables', file_stat)
                        write(fnr,'(a)') trim(movienames(imovie))
                        close(unit=fnr)
                        call flush(unit=fnr)
                    else
                        call del_file(individual_filetabs(imovie))
                    endif
                end do
            end subroutine create_individual_filetables

            subroutine setup_distr_env
                ! setup the environment for distributed execution
                call setup_qsys_env(p_master, qsys_fac, myqsys, parts, qscripts, myq_descr, lmask_stream)
                ! prepare job description
                call cline%gen_job_descr(job_descr)
                ! prepare scripts
                call qsys_cleanup(p_master)
                call qscripts%generate_scripts(job_descr, p_master%ext, myq_descr, part_params=part_params)
            end subroutine setup_distr_env

    end subroutine exec_preproc_stream

end module simple_commander_stream_wflows