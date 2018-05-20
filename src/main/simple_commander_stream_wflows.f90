! concrete commander: stream processing routines
module simple_commander_stream_wflows
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_sp_project,     only: sp_project
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs,      only: qsys_cleanup
use simple_parameters,     only: parameters
implicit none

public :: preprocess_stream_commander
public :: cluster2D_stream_distr_commander
private

type, extends(commander_base) :: preprocess_stream_commander
  contains
    procedure :: execute      => exec_preprocess_stream
end type preprocess_stream_commander
type, extends(commander_base) :: cluster2D_stream_distr_commander
  contains
    procedure :: execute      => exec_cluster2D_stream_distr
end type cluster2D_stream_distr_commander

contains

    subroutine exec_preprocess_stream( self, cline )
        use simple_moviewatcher,         only: moviewatcher
        use simple_qsys_funs,            only: qsys_cleanup
        use simple_commander_preprocess, only: preprocess_commander
        class(preprocess_stream_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)                       :: params
        integer,                   parameter   :: SHORTTIME = 30   ! folder watched every minute
        integer,                   parameter   :: LONGTIME  = 10
        class(cmdline),            allocatable :: completed_jobs_clines(:)
        type(qsys_env)                         :: qenv
        type(moviewatcher)                     :: movie_buff
        type(sp_project)                       :: spproj, stream_spproj
        character(len=LONGSTRLEN), allocatable :: movies(:), prev_movies(:)
        character(len=:),          allocatable :: output_dir, output_dir_ctf_estimate, output_dir_picker
        character(len=:),          allocatable :: output_dir_motion_correct, output_dir_extract, stream_spprojfile
        character(len=LONGSTRLEN)              :: movie
        integer                                :: nmovies, imovie, stacksz, prev_stacksz, iter, icline
        integer                                :: nptcls, nptcls_prev, nmovs, nmovs_prev
        logical                                :: l_pick
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
        params = parameters(cline)
        params%stream     = 'yes'
        params%split_mode = 'stream'
        params%numlen     = 5
        params%ncunits    = params%nparts
        call cline%set('numlen', real(params%numlen))
        call cline%set('mkdir', 'no')
        call cline%set('prg',   'preprocess')
        l_pick = cline%defined('refs')
        ! read in movies
        call spproj%read( params%projfile )
        ! check for previously processed movies
        call spproj%get_movies_table(prev_movies)
        ! output directories
        output_dir = './'
        output_dir_ctf_estimate   = trim(output_dir)//trim(DIR_CTF_ESTIMATE)
        output_dir_motion_correct = trim(output_dir)//trim(DIR_MOTION_CORRECT)
        call simple_mkdir(output_dir)
        call simple_mkdir(output_dir_ctf_estimate)
        call simple_mkdir(output_dir_motion_correct)
        if( l_pick )then
            output_dir_picker  = trim(output_dir)//trim(DIR_PICKER)
            output_dir_extract = trim(output_dir)//trim(DIR_EXTRACT)
            call simple_mkdir(output_dir_picker)
            call simple_mkdir(output_dir_extract)
        endif
        ! setup the environment for distributed execution
        call qenv%new( stream=.true. )
        ! movie watcher init
        movie_buff = moviewatcher(longtime, prev_movies)
        ! start watching
        prev_stacksz = 0
        nmovies      = 0
        iter         = 0
        do
            if( file_exists(trim(TERM_STREAM)) )then
                write(*,'(A)')'>>> TERMINATING PREPROCESSING STREAM'
                exit
            endif
            iter = iter + 1
            call movie_buff%watch( nmovies, movies )
            ! append movies to processing stack
            if( nmovies > 0 )then
                do imovie = 1, nmovies
                    movie = trim(adjustl(movies(imovie)))
                    call create_individual_project
                    call qenv%qscripts%add_to_streaming( cline )
                enddo
            endif
            ! stream scheduling
            call qenv%qscripts%schedule_streaming( qenv%qdescr )
            stacksz = qenv%qscripts%get_stacksz()
            if( stacksz .ne. prev_stacksz )then
                prev_stacksz = stacksz
                write(*,'(A,I5)')'>>> MOVIES TO PROCESS: ', stacksz
            endif
            ! completed jobs update the current project
            if( qenv%qscripts%get_done_stacksz() > 0 )then
                ! append new processed movies to project
                call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
                nptcls_prev = spproj%get_nptcls()
                nmovs_prev  = spproj%get_nmovies()
                do icline=1,size(completed_jobs_clines)
                    stream_spprojfile = completed_jobs_clines(icline)%get_carg('projfile')
                    call stream_spproj%read( stream_spprojfile )
                    call spproj%append_project(stream_spproj, 'mic')
                    if( l_pick )call spproj%append_project(stream_spproj, 'stk')
                    call stream_spproj%kill()
                    deallocate(stream_spprojfile)
                enddo
                nptcls = spproj%get_nptcls()
                nmovs  = spproj%get_nmovies()
                ! write
                if( nmovs == 0 )then
                    ! first write
                    call spproj%write
                else
                    ! write inside
                    call spproj%write_segment_inside('mic', fromto=[nmovs_prev+1, nmovs])
                    if( l_pick )then
                        call spproj%write_segment_inside('stk',    fromto=[nmovs_prev+1, nmovs])
                        call spproj%write_segment_inside('ptcl2D', fromto=[nptcls_prev+1, nptcls])
                        call spproj%write_segment_inside('ptcl3D', fromto=[nptcls_prev+1, nptcls])
                    endif
                endif
                call update_projects_list
                deallocate(completed_jobs_clines)
            endif
            ! wait
            call simple_sleep(SHORTTIME)
        end do
        ! termination
        call spproj%write
        call spproj%kill
        ! cleanup
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_PREPROCESS_STREAM NORMAL STOP ****')
        contains

            subroutine update_projects_list
                type(cmdline) :: cline_mov
                character(len=:), allocatable :: fname
                character(len=LONGSTRLEN), allocatable :: old_fnames(:), fnames(:)
                integer :: i, n_spprojs, n_old
                n_spprojs = size(completed_jobs_clines)
                if( n_spprojs == 0 )return
                if( file_exists(STREAM_SPPROJFILES) )then
                    ! append
                    call read_filetable(STREAM_SPPROJFILES, old_fnames)
                    n_old = size(old_fnames)
                    allocate(fnames(n_spprojs+n_old))
                    fnames(1:n_old) = old_fnames(:)
                    do i=1,n_spprojs
                        cline_mov       = completed_jobs_clines(i)
                        fname           = trim(cline_mov%get_carg('projfile'))
                        fnames(n_old+i) = trim(fname)
                    enddo
                else
                    ! first write
                    allocate(fnames(n_spprojs))
                    do i=1,n_spprojs
                        cline_mov = completed_jobs_clines(i)
                        fname     = trim(cline_mov%get_carg('projfile'))
                        fnames(i) = trim(fname)
                    enddo
                endif
                call write_filetable(STREAM_SPPROJFILES, fnames)
            end subroutine update_projects_list

            subroutine create_individual_project
                type(sp_project)              :: spproj_here
                type(cmdline)                 :: cline_here
                type(ctfparams)               :: ctfvars
                character(len=STDLEN)         :: ext, movie_here
                character(len=LONGSTRLEN)     :: projname, projfile
                movie_here = basename(trim(movie))
                ext        = fname2ext(trim(movie_here))
                projname   = 'preprocess_'//trim(get_fbody(trim(movie_here), trim(ext)))
                projfile   = trim(projname)//trim(METADATA_EXT)
                call cline_here%set('projname', trim(projname))
                call cline_here%set('projfile', trim(projfile))
                call spproj_here%update_projinfo(cline_here)
                spproj_here%compenv  = spproj%compenv
                spproj_here%jobproc  = spproj%jobproc
                ctfvars%ctfflag      = CTFFLAG_YES
                ctfvars%smpd         = params%smpd
                ctfvars%cs           = params%cs
                ctfvars%kv           = params%kv
                ctfvars%fraca        = params%fraca
                ctfvars%l_phaseplate = params%phaseplate.eq.'yes'
                call spproj_here%add_single_movie(trim(movie), ctfvars)
                call spproj_here%write()
                call spproj_here%kill()
                call cline%set('projname', trim(projname))
                call cline%set('projfile', trim(projfile))
            end subroutine create_individual_project

    end subroutine exec_preprocess_stream

    subroutine exec_cluster2D_stream_distr( self, cline )
        use simple_commander_distr_wflows, only: cluster2D_distr_commander, make_cavgs_distr_commander,scale_project_distr_commander
        use simple_ori,                    only: ori
        use simple_oris,                   only: oris
        use simple_image,                  only: image
        class(cluster2D_stream_distr_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        character(len=*), parameter   :: STK_FILETAB        = 'stkstreamtab.txt'
        character(len=:), allocatable :: WORK_PROJFILE
        integer,          parameter   :: SHIFTSRCH_PTCLSLIM = 2000  ! # of ptcls required to turn on shift search
        integer,          parameter   :: SHIFTSRCH_ITERLIM  = 5     ! # of iterations prior to turn on shift search
        integer,          parameter   :: CCRES_NPTCLS_LIM   = 10000 ! # of ptcls required to turn on objfun=ccres
        integer,          parameter   :: WAIT_WATCHER       = 60    ! seconds prior to new stack detection
        integer,          parameter   :: MAXNCLS            = 1000  ! maximum # of classes
        type(parameters)                    :: params
        type(cluster2D_distr_commander)     :: xcluster2D_distr
        type(make_cavgs_distr_commander)    :: xmake_cavgs
        type(scale_project_distr_commander) :: xscale_distr
        type(cmdline)                       :: cline_cluster2D, cline_scale, cline_make_cavgs
        type(sp_project)                    :: orig_proj, work_proj
        character(len=:),       allocatable :: orig_projfile
        character(len=STDLEN)               :: str_iter, refs_glob
        real    :: orig_smpd, msk, scale_factor, orig_msk
        integer :: iter, n_newstks, orig_box, box, nptcls_glob, ncls_glob, tnow
        integer :: nptcls_glob_prev, ncls_glob_prev, last_injection, n_stk, n_stk_prev
        logical :: do_autoscale
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        call cline%set('stream','yes') ! only for parameters determination
        params = parameters(cline)
        if( .not.file_exists(params%projfile) )stop 'Project does not exist!'
        call cline%set('stream','no') ! was only for parameters determination
        orig_projfile = params%projfile
        allocate(WORK_PROJFILE, source='cluster2D_stream_tmproj.simple')
        ! init command-lines
        cline_cluster2D  = cline
        call cline_cluster2D%set('prg',       'cluster2D')
        call cline_cluster2D%set('autoscale', 'no')
        call cline_cluster2D%set('extr_iter', 100.) ! no extremal randomization
        call cline_cluster2D%delete('frac_update')
        call cline_cluster2D%delete('projfile_target')
        call cline_cluster2D%delete('projname')
        call cline_cluster2D%set('projfile', trim(WORK_PROJFILE))
        call cline_cluster2D%set('projname', trim(get_fbody(trim(WORK_PROJFILE),trim('simple'))))
        cline_make_cavgs = cline
        call cline_make_cavgs%set('prg',    'make_cavgs')
        call cline_make_cavgs%set('refs',   'cavgs_final'//params%ext)
        call cline_make_cavgs%delete('autoscale')
        call cline_make_cavgs%set('projfile', trim(WORK_PROJFILE)) ! TO WORK OUT
        ! wait for the first stacks to trickle in...
        do
            if( .not.is_file_open(params%projfile_target) )then
                call orig_proj%read(params%projfile_target)
                n_stk = orig_proj%os_stk%get_noris()
                if( n_stk > 0)then
                    nptcls_glob = orig_proj%get_nptcls()
                    if( nptcls_glob > params%ncls_start * params%nptcls_per_cls )then
                        ! Enough particles to start 2D clustering...
                        exit
                    endif
                endif
            endif
            call simple_sleep(WAIT_WATCHER) ! parameter instead
        enddo
        ! transfer of general info
        call work_proj%read(params%projfile)
        orig_proj%projinfo = work_proj%projinfo
        orig_proj%compenv  = work_proj%compenv
        if( orig_proj%jobproc%get_noris()>0 ) orig_proj%jobproc = work_proj%jobproc
        call work_proj%kill()
        call orig_proj%projinfo%delete_entry('projname')
        call orig_proj%projinfo%delete_entry('projfile')
        call orig_proj%update_projinfo(cline_cluster2D) ! name change
        call orig_proj%write()                          ! & write
        call orig_proj%kill()
        ! scaling
        call work_proj%read(trim(WORK_PROJFILE))
        orig_box  = work_proj%get_box()
        orig_smpd = work_proj%get_smpd()
        orig_msk  = params%msk
        params%smpd_targets2D(1) = max(orig_smpd, params%lp*LP2SMPDFAC2D)
        do_autoscale = params%autoscale.eq.'yes' .and. params%smpd_targets2D(1) > orig_smpd
        if( do_autoscale )then
            deallocate(WORK_PROJFILE)
            call work_proj%scale_projfile(params%smpd_targets2D(1), WORK_PROJFILE, cline_cluster2D, cline_scale)
            scale_factor = cline_scale%get_rarg('scale')
            box          = nint(cline_scale%get_rarg('newbox'))
            msk          = cline_cluster2D%get_rarg('msk')
            call xscale_distr%execute( cline_scale )
        else
            scale_factor = 1.
            box          = orig_box
            msk          = orig_msk
        endif
        call cline_cluster2D%set('projfile', trim(WORK_PROJFILE))
        call cline_cluster2D%set('box',      real(box))
        call cline_cluster2D%set('msk',      real(msk))
        call work_proj%kill()
        ! Main loop
        ncls_glob      = nint(real(nptcls_glob) / real(params%nptcls_per_cls))
        last_injection = simple_gettime()
        do iter = 1, 999
            str_iter = int2str_pad(iter,3)
            tnow     = simple_gettime()
            if(tnow-last_injection > params%time_inactive)then
                write(*,*)'>>> TIME LIMIT WITHOUT NEW PARTICLES ADDITION REACHED'
                exit
            endif
            n_stk_prev       = n_stk
            nptcls_glob_prev = nptcls_glob
            ! CLUSTER2D
            call cline_cluster2D%delete('endit')
            call cline_cluster2D%set('ncls',    real(ncls_glob))
            call cline_cluster2D%set('startit', real(iter))
            call cline_cluster2D%set('maxits',  real(iter))
            call cline_cluster2D%set('ncls',    real(ncls_glob))
            call cline_cluster2D%set('nparts',  real(min(ncls_glob,params%nparts)))
            if( nptcls_glob > SHIFTSRCH_PTCLSLIM .and. iter > SHIFTSRCH_ITERLIM )then
                call cline_cluster2D%set('trs', MINSHIFT)
            endif
            if( nptcls_glob > CCRES_NPTCLS_LIM )then
                call cline_cluster2D%set('objfun', 'ccres')
            endif
            if(nptcls_glob > SHIFTSRCH_PTCLSLIM .and. iter > SHIFTSRCH_ITERLIM)then
                call cline_cluster2D%set('trs', MINSHIFT)
            endif
            ! execute
            call xcluster2D_distr%execute(cline_cluster2D)
            ! update
            call work_proj%kill
            call work_proj%read(trim(WORK_PROJFILE))
            ! current refernce file name
            refs_glob = trim(CAVGS_ITER_FBODY)//trim(str_iter)//trim(params%ext)
            ! remap zero-population classes
            call remap_empty_classes
            ! user termination
            if( file_exists(trim(TERM_STREAM)) )then
                    write(*,'(A)')'>>> TERMINATING CLUSTER2D STREAM'
                    exit
                endif
            ! detects whether new images have been added to the original project
            if( .not.is_file_open(params%projfile_target) )then
                call orig_proj%read_segment('stk', params%projfile_target)
                n_stk     = orig_proj%os_stk%get_noris()
                n_newstks = n_stk - n_stk_prev
            endif
            if(n_newstks > 0)then
                call work_proj%kill
                call work_proj%read(params%projfile)
                call append_newstacks( orig_proj%os_stk )
                ! updates # of classes
                ncls_glob_prev = ncls_glob
                ncls_glob      = nint(real(nptcls_glob) / params%nptcls_per_cls)
                ncls_glob      = min(ncls_glob, MAXNCLS)
                ! allocate new particles to previous references
                call map_new_ptcls
                last_injection = simple_gettime()
            endif

            ! wait
            call simple_sleep(WAIT_WATCHER)
        enddo
        ! un-scaling
        if( do_autoscale )then
            call work_proj%os_ptcl2D%mul_shifts( 1./scale_factor )
            call work_proj%os_stk%set_all2single('box',  real(orig_box))
            call work_proj%os_stk%set_all2single('smpd', real(orig_smpd))
        endif
        call work_proj%write()
        ! class averages at original sampling
        call cline_make_cavgs%set('ncls', real(ncls_glob))
        call xmake_cavgs%execute(cline_make_cavgs) ! should be distributed
        ! cleanup
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_CLUSTER2D_STREAM NORMAL STOP ****')

        contains

            subroutine append_newstacks( os_stk_in )
                use simple_commander_imgproc, only: scale_commander
                class(oris),  intent(inout) :: os_stk_in
                character(len=*), parameter :: SCALE_FILETAB = 'stkscale.txt'
                type(oris)                         :: os_ptcls
                type(ori)                          :: o_stk
                type(ctfparams)                    :: ctfvars
                type(qsys_env)                     :: qenv
                character(len=LONGSTRLEN), allocatable :: new_stks(:)
                character(len=:),      allocatable     :: stkname, ctfstr, phaseplate
                integer :: istk, nptcls, iptcl
                allocate(new_stks(n_newstks))
                do istk=1,n_newstks
                    call os_stk_in%getter(n_stk+istk, 'stk', stkname)
                    new_stks(istk) = trim(stkname)
                enddo
                ! scale stacks to append
                if( params%autoscale.eq.'yes' )then
                    call qenv%new()
                    call cline_scale%delete('stk')
                    call cline_scale%delete('stktab')
                    call cline_scale%set('dir_target', trim(SCSTK_DIR))
                    call cline_scale%set('filetab', trim(STK_FILETAB))
                    call cline_scale%set('stream', 'yes')
                    call cline_scale%set('numlen', 1.)
                    call write_filetable(SCALE_FILETAB, new_stks)
                    call qenv%exec_simple_prg_in_queue(cline_scale, 'OUT1','JOB_FINISHED_1')
                    call qsys_cleanup
                    call del_file(SCALE_FILETAB)
                    do istk=1,n_newstks
                        new_stks(istk) = SCSTK_DIR//add2fbody(new_stks(istk), params%ext, SCALE_SUFFIX)
                    enddo
                endif
                ! append to working project
                do istk=1,n_newstks
                    o_stk  = os_stk_in%get_ori(n_stk+istk)
                    nptcls = nint(o_stk%get('nptcls'))
                    call os_ptcls%new(nptcls)
                    do iptcl=1,nptcls
                        call os_ptcls%set_ori(iptcl, o_stk)
                    enddo
                    call o_stk%getter('ctf', ctfstr)
                    ctfvars%ctfflag = 1
                    select case( trim(ctfstr) )
                        case('no')
                            ctfvars%ctfflag = 0
                        case('yes')
                            ctfvars%ctfflag = 1
                        case('flip')
                            ctfvars%ctfflag = 2
                    end select
                    ctfvars%kv    = o_stk%get('kv')
                    ctfvars%cs    = o_stk%get('cs')
                    ctfvars%fraca = o_stk%get('fraca')
                    if( o_stk%isthere('phaseplate'))then
                        call o_stk%getter('phaseplate', phaseplate)
                        ctfvars%l_phaseplate = trim(phaseplate) .eq. 'yes'
                    else
                        ctfvars%l_phaseplate = .false.
                    endif
                    call work_proj%add_stk(new_stks(istk), ctfvars, os_ptcls)
                enddo
                ! update global parameters
                n_stk       = work_proj%os_stk%get_noris()
                nptcls_glob = work_proj%get_nptcls()
                write(*,'(A,I6)')'>>> NEW PARTICLES COUNT: ', nptcls_glob
            end subroutine append_newstacks

            subroutine remap_empty_classes
                integer, allocatable  :: fromtocls(:,:)
                type(image)           :: img_cavg
                character(len=STDLEN) :: stk
                real                  :: smpd
                integer               :: icls
                call work_proj%os_ptcl2D%fill_empty_classes(ncls_glob, fromtocls)
                if( allocated(fromtocls) )then
                    smpd = work_proj%get_smpd()
                    ! updates document & classes
                    call work_proj%write()
                    call img_cavg%new([box, box,1], smpd)
                    do icls = 1, size(fromtocls, dim=1)
                        ! cavg
                        call img_cavg%read(trim(refs_glob), fromtocls(icls, 1))
                        call img_cavg%write(trim(refs_glob), fromtocls(icls, 2))
                        ! even & odd
                        stk = add2fbody(trim(refs_glob),params%ext,'_even')
                        call img_cavg%read(stk, fromtocls(icls, 1))
                        call img_cavg%write(stk, fromtocls(icls, 2))
                        stk = add2fbody(trim(refs_glob),params%ext,'_odd')
                        call img_cavg%read(stk, fromtocls(icls, 1))
                        call img_cavg%write(stk, fromtocls(icls, 2))
                    enddo
                    ! stack size preservation
                    call img_cavg%read(trim(refs_glob), ncls_glob)
                    call img_cavg%write(trim(refs_glob), ncls_glob)
                    stk = add2fbody(trim(refs_glob),params%ext,'_even')
                    call img_cavg%read(stk, ncls_glob)
                    call img_cavg%write(stk, ncls_glob)
                    stk = add2fbody(trim(refs_glob),params%ext,'_odd')
                    call img_cavg%read(stk, ncls_glob)
                    call img_cavg%write(stk, ncls_glob)
                endif
            end subroutine remap_empty_classes

            subroutine map_new_ptcls
                use simple_projection_frcs, only: projection_frcs
                type(ran_tabu)        :: rt
                type(projection_frcs) :: frcs_prev, frcs
                type(image)           :: img_cavg
                integer, allocatable  :: fromtocls(:,:), cls(:)
                real                  :: smpd
                integer               :: icls, nptcls, iptcl, i, state
                character(len=STDLEN) :: stk
                nptcls = nptcls_glob - nptcls_glob_prev
                state  = 1
                ! randomize new ptcls to previous references
                allocate(cls(nptcls))
                rt = ran_tabu(nptcls)
                call rt%balanced(ncls_glob_prev, cls)
                do i = 1, nptcls
                    iptcl = nptcls_glob_prev + i
                    call work_proj%os_ptcl2D%set(iptcl, 'class', real(cls(i)))
                    call work_proj%os_ptcl2D%set(iptcl, 'w',     1.)
                    call work_proj%os_ptcl2D%set(iptcl, 'corr',  0.)
                enddo
                deallocate(cls)
                ! updates references &FRCs
                if( ncls_glob.eq.ncls_glob_prev )then
                    ! nothing to do
                else
                    write(*,'(A,I6)')'>>> NEW CLASSES COUNT: ', ncls_glob
                    call work_proj%os_ptcl2D%fill_empty_classes(ncls_glob, fromtocls)
                    if( allocated(fromtocls) )then
                        ! references
                        smpd = work_proj%get_smpd()
                        call img_cavg%new([box,box,1], smpd)
                        do icls = 1, size(fromtocls, dim=1)
                            ! cavg
                            call img_cavg%read(trim(refs_glob), fromtocls(icls, 1))
                            call img_cavg%write(trim(refs_glob), fromtocls(icls, 2))
                            ! even & odd
                            stk = add2fbody(trim(refs_glob),params%ext,'_even')
                            call img_cavg%read(stk, fromtocls(icls, 1))
                            call img_cavg%write(stk, fromtocls(icls, 2))
                            stk = add2fbody(trim(refs_glob),params%ext,'_odd')
                            call img_cavg%read(stk, fromtocls(icls, 1))
                            call img_cavg%write(stk, fromtocls(icls, 2))
                        enddo
                        ! stack size preservation
                        call img_cavg%read(trim(refs_glob), ncls_glob)
                        call img_cavg%write(trim(refs_glob), ncls_glob)
                        stk = add2fbody(trim(refs_glob),params%ext,'_even')
                        call img_cavg%read(stk, ncls_glob)
                        call img_cavg%write(stk, ncls_glob)
                        stk = add2fbody(trim(refs_glob),params%ext,'_odd')
                        call img_cavg%read(stk, ncls_glob)
                        call img_cavg%write(stk, ncls_glob)
                        ! FRCs
                        call frcs_prev%new(ncls_glob_prev, box, smpd, state)
                        call frcs%new(ncls_glob, box, smpd, state)
                        call frcs_prev%read(FRCS_FILE)
                        do icls = 1, ncls_glob_prev
                            call frcs%set_frc(icls,&
                            &frcs_prev%get_frc(icls, box, state), state)
                        enddo
                        do icls = 1, size(fromtocls, dim=1)
                            call frcs%set_frc( fromtocls(icls,2),&
                            &frcs%get_frc(fromtocls(icls,1), box, state), state)
                        enddo
                        call frcs%write(FRCS_FILE)
                    endif
                endif
            end subroutine map_new_ptcls

    end subroutine exec_cluster2D_stream_distr

end module simple_commander_stream_wflows
