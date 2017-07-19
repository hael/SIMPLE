module simple_commander_stream_wflows
use simple_defs
use simple_cmdline,           only: cmdline
use simple_chash,             only: chash
use simple_params,            only: params
use simple_commander_base,    only: commander_base
use simple_strings,           only: real2str
use simple_commander_preproc, only: preproc_commander
use simple_qsys_env,          only: qsys_env
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

    ! subroutine exec_preproc_stream( self, cline )
    !     use simple_commander_preproc, only: preproc_commander
    !     class(preproc_stream_commander), intent(inout) :: self
    !     class(cmdline),                  intent(inout) :: cline
    !     logical,               parameter   :: DEBUG = .true.
    !     character(len=STDLEN), parameter   :: FILETABNAME='movieftab_preproc_stream.txt'
    !     integer,               parameter   :: SHORTTIME = 15
    !     character(len=STDLEN), allocatable :: movienames(:)
    !     type(qsys_env)           :: qenv
    !     type(params)             :: p_master
    !     type(chash)              :: job_descr
    !     type(chash), allocatable :: part_params(:)
    !     logical,     allocatable :: jobs_done(:), jobs_submitted(:)
    !     integer                  :: nmovies, nmovies_prev
    !     integer, parameter       :: TRAILING=5
    !     ! make master parameters
    !     p_master = params(cline, checkdistr=.false.)
    !     ! set defaults
    !     p_master%split_mode = 'singles'
    !     p_master%numlen     = 5
    !     call cline%set('numlen', real(p_master%numlen))
    !     if( cline%defined('fbody') )then
    !         call cline%set('fbody', trim(p_master%dir_target)//'/'//trim(p_master%fbody))
    !     else
    !         call cline%set('fbody', trim(p_master%dir_target)//'/')
    !     endif
    !     ! make target directory
    !     call exec_cmdline('mkdir -p '//trim(adjustl(p_master%dir_target)))
    !     nmovies = 0
    !     nmovies_prev = 0
    !     do
    !         ! generate the filetable via system call
    !         call sys_gen_mrcfiletab(p_master%dir_movies, FILETABNAME)
    !         ! read the movienames back in
    !         call read_filetable(FILETABNAME, movienames)
    !         nmovies_prev = nmovies
    !         nmovies      = size(movienames)
    !         call cline%set('nparts', real(nmovies)) ! need dyn update of nparts 4 stream
    !         p_master%nparts = nmovies               ! need dyn update of nparts 4 stream
    !         p_master%nptcls = nmovies               ! need dyn update of nparts 4 stream
    !         call create_individual_filetables       ! 1-of-1 ftab but index comes from part
    !         call setup_distr_env                    ! just to reduce complexity
    !         ! manage job scheduling
    !         if( qenv%exists() )then
    !             call qenv%qscripts%update_queue
    !             call qenv%qscripts%submit_scripts
    !         endif
    !         call simple_sleep(SHORTTIME)
    !     end do

    !     contains

    !         subroutine create_individual_filetables
    !             integer :: imovie, fnr, file_stat
    !             character(len=STDLEN), allocatable :: individual_filetabs(:)
    !             individual_filetabs = make_filenames('movie_to_process', nmovies, '.txt', p_master%numlen)
    !             if( allocated(part_params) ) deallocate(part_params)
    !             allocate(part_params(size(individual_filetabs)))
    !             do imovie=1,nmovies
    !                 call part_params(imovie)%new(1)
    !                 call part_params(imovie)%set('filetab', trim(individual_filetabs(imovie)))
    !                 fnr = get_fileunit()
    !                 open(unit=fnr, status='replace', action='write', file=trim(individual_filetabs(imovie)), iostat=file_stat)
    !                 call fopen_err('exec_preproc_stream :: create_individual_filetables', file_stat)
    !                 write(fnr,'(a)') trim(movienames(imovie))
    !                 close(unit=fnr)
    !             end do
    !         end subroutine create_individual_filetables

    !         subroutine setup_distr_env
    !             ! get the old jobs status
    !             if( qenv%exists() ) call qenv%qscripts%get_jobs_status( jobs_done, jobs_submitted )
    !             ! setup the environment for distributed execution
    !             call qenv%new(p_master, stream=.true.)
    !             ! put back the old jobs status
    !             if( allocated(jobs_done) ) call qenv%qscripts%set_jobs_status( jobs_done, jobs_submitted )
    !             ! prepare job description
    !             call cline%gen_job_descr(job_descr)
    !             ! prepare scripts
    !             if( nmovies > nmovies_prev )then
    !                 call qenv%qscripts%generate_scripts(job_descr, p_master%ext, qenv%qdescr, part_params=part_params)
    !             endif
    !         end subroutine setup_distr_env

    ! end subroutine exec_preproc_stream

    subroutine exec_preproc_stream( self, cline )
        use simple_commander_preproc, only: preproc_commander
        use simple_moviewatcher,      only: moviewatcher
        class(preproc_stream_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        logical,               parameter   :: DEBUG = .true.
        character(len=STDLEN), parameter   :: FILETABNAME='movieftab_preproc_stream.txt'
        integer,               parameter   :: SHORTTIME = 15
        character(len=STDLEN), allocatable :: movienames(:)
        type(qsys_env)           :: qenv
        type(params)             :: p_master
        type(chash)              :: job_descr
        type(moviewatcher)       :: movie_buff
        type(chash), allocatable :: part_params(:)
        logical,     allocatable :: jobs_done(:), jobs_submitted(:)
        integer                  :: nmovies, nmovies_prev
        integer, parameter       :: TRAILING=5
        ! make master parameters
        p_master = params(cline, checkdistr=.false.)
        ! set defaults
        p_master%split_mode = 'singles'
        p_master%numlen     = 5
        call cline%set('numlen', real(p_master%numlen))
        if( cline%defined('fbody') )then
            call cline%set('fbody', trim(p_master%dir_target)//'/'//trim(p_master%fbody))
        else
            call cline%set('fbody', trim(p_master%dir_target)//'/')
        endif
        ! make target directory
        call exec_cmdline('mkdir -p '//trim(adjustl(p_master%dir_target)))
        ! movie watcher init
        movie_buff   = moviewatcher(p_master%dir_movies, SHORTTIME, print=.true.)
        nmovies      = 0
        do
            nmovies_prev = nmovies
            call movie_buff%watch( nmovies, movienames )
            call cline%set('nparts', real(nmovies)) ! need dyn update of nparts 4 stream
            p_master%nparts = nmovies               ! need dyn update of nparts 4 stream
            p_master%nptcls = nmovies               ! need dyn update of nparts 4 stream
            call create_individual_filetables       ! 1-of-1 ftab but index comes from part
            call setup_distr_env                    ! just to reduce complexity
            ! manage job scheduling
            if( qenv%exists() )then
                call qenv%qscripts%update_queue
                call qenv%qscripts%submit_scripts
            endif
            ! wait...
            call simple_sleep(SHORTTIME)
        end do

        contains

            subroutine create_individual_filetables
                character(len=STDLEN), allocatable :: individual_filetabs(:)
                integer               :: imovie, fnr, file_stat
                character(len=STDLEN) :: fname
                individual_filetabs = make_filenames('movie_to_process', nmovies, '.txt', p_master%numlen)
                if( allocated(part_params) ) deallocate(part_params)
                allocate(part_params(size(individual_filetabs)))
                do imovie=1,nmovies
                    fname = trim(individual_filetabs(imovie))
                    call part_params(imovie)%new(1)
                    call part_params(imovie)%set('filetab', fname)
                    fnr = get_fileunit()
                    open(unit=fnr, status='replace', action='write', file=fname, iostat=file_stat)
                    call fopen_err('exec_preproc_stream :: create_individual_filetables', file_stat)
                    write(fnr,'(a)') trim(movienames(imovie))
                    close(unit=fnr)
                end do
            end subroutine create_individual_filetables

            subroutine setup_distr_env
                ! get the old jobs status
                if( qenv%exists() ) call qenv%qscripts%get_jobs_status( jobs_done, jobs_submitted )
                ! setup the environment for distributed execution
                call qenv%new(p_master, stream=.true.)
                ! put back the old jobs status
                if( allocated(jobs_done) ) call qenv%qscripts%set_jobs_status( jobs_done, jobs_submitted )
                ! prepare job description
                call cline%gen_job_descr(job_descr)
                ! prepare scripts
                if( nmovies > nmovies_prev )then
                    call qenv%qscripts%generate_scripts(job_descr, p_master%ext, qenv%qdescr, part_params=part_params)
                endif
            end subroutine setup_distr_env

    end subroutine exec_preproc_stream

end module simple_commander_stream_wflows