! concrete commander: stream processing routines
module simple_commander_stream_wflows
include 'simple_lib.f08'
use simple_cmdline,              only: cmdline
use simple_params,               only: params
use simple_commander_base,       only: commander_base
use simple_oris,                 only: oris
use simple_image,                only: image
use simple_sp_project,           only: sp_project
use simple_binoris_io              ! use all in there
use simple_qsys_env,             only: qsys_env
use simple_qsys_funs,            only: qsys_cleanup
use simple_binoris_io              ! use all in there

!use simple_commander_distr,      only: merge_algndocs_commander! use all in there
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
        use simple_commander_preprocess, only: preprocess_commander
        use simple_moviewatcher,         only: moviewatcher
        class(preprocess_stream_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        integer,               parameter   :: SHORTTIME = 60   ! folder watched every minute
        integer,               parameter   :: LONGTIME  = 900  ! 15 mins before processing a new movie
        class(cmdline),       allocatable  :: completed_jobs(:)
        type(qsys_env)                     :: qenv
        type(params)                       :: p_master
        type(moviewatcher)                 :: movie_buff
        type(sp_project)                   :: spproj, stream_spproj
        character(len=STDLEN), allocatable :: movies(:)
        character(len=:),      allocatable :: output_dir, output_dir_ctf_estimate, output_dir_picker
        character(len=:),      allocatable :: output_dir_motion_correct, output_dir_extract, stream_spprojfile
        character(len=LONGSTRLEN)          :: movie
        integer                            :: nmovies, imovie, stacksz, prev_stacksz, iter, icline
        logical                            :: l_pick
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
        call cline%set('prg', 'preprocess')
        ! make master parameters
        p_master        = params(cline)
        p_master%stream = 'yes'
        l_pick          = cline%defined('refs')
        ! output directories
        if( cline%defined('dir') )then
            output_dir = trim(p_master%dir)//'/'
        else
            output_dir = trim(DIR_PREPROC)
        endif
        output_dir_ctf_estimate   = trim(output_dir)//trim(DIR_CTF_ESTIMATE)
        output_dir_motion_correct = trim(output_dir)//trim(DIR_MOTION_CORRECT)
        call mkdir(output_dir)
        call mkdir(output_dir_ctf_estimate)
        call mkdir(output_dir_motion_correct)
        if( l_pick )then
            output_dir_picker  = trim(output_dir)//trim(DIR_PICKER)
            output_dir_extract = trim(output_dir)//trim(DIR_EXTRACT)
            call mkdir(output_dir_picker)
            call mkdir(output_dir_extract)
        endif
        ! read in movies
        call spproj%read( p_master%projfile )
        ! set defaults
        p_master%numlen     = 5
        call cline%set('numlen', real(p_master%numlen))
        p_master%split_mode = 'stream'
        p_master%ncunits    = p_master%nparts
        ! setup the environment for distributed execution
        call qenv%new(p_master, stream=.true.)
        ! movie watcher init
        movie_buff = moviewatcher(p_master, longtime, print=.true.)
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
                    call create_individual_project()
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
                call qenv%qscripts%get_stream_done_stack( completed_jobs )
                do icline=1,size(completed_jobs)
                    stream_spprojfile = completed_jobs(icline)%get_carg('projfile')
                    call stream_spproj%read( stream_spprojfile )
                    call spproj%append_project( stream_spproj, 'mic' )
                    if( l_pick )call spproj%append_project( stream_spproj, 'stk' )
                    call stream_spproj%kill()
                    call del_file(stream_spprojfile)
                    deallocate(stream_spprojfile)
                enddo
                call spproj%write()
                deallocate(completed_jobs)
            endif
            ! wait
            call simple_sleep(SHORTTIME)
        end do
        ! cleanup
        call qsys_cleanup(p_master)
        ! end gracefully
        call simple_end('**** SIMPLE_PREPROCESS_STREAM NORMAL STOP ****')
        contains

            subroutine create_individual_project
                type(sp_project)              :: spproj_here
                type(cmdline)                 :: cline_here
                type(ctfparams)               :: ctfvars
                character(len=STDLEN)         :: ext, movie_here
                character(len=LONGSTRLEN)     :: projname
                movie_here = remove_abspath(trim(movie))
                ext        = fname2ext(trim(movie_here))
                projname   = 'preprocess_'//trim(get_fbody(trim(movie_here), trim(ext)))
                call cline_here%set('projname', trim(projname)) ! necessary?
                call spproj_here%update_projinfo(cline_here)
                spproj_here%compenv  = spproj%compenv
                spproj_here%jobproc  = spproj%jobproc
                ctfvars%ctfflag      = CTFFLAG_YES
                ctfvars%smpd         = p_master%smpd
                ctfvars%cs           = p_master%cs
                ctfvars%kv           = p_master%kv
                ctfvars%fraca        = p_master%fraca
                ctfvars%l_phaseplate = p_master%phaseplate.eq.'yes'
                call spproj_here%add_single_movie(trim(movie), ctfvars)
                call spproj_here%write()
                call spproj_here%kill()
                call cline%set('projfile', trim(projname)//'.simple')
            end subroutine create_individual_project

    end subroutine exec_preprocess_stream

    subroutine exec_cluster2D_stream_distr( self, cline )
        use simple_commander_distr_wflows, only: cluster2D_distr_commander, make_cavgs_distr_commander,scale_project_distr_commander
        use simple_oris,                   only: oris
        use simple_image,                  only: image
        use simple_commander_distr         ! use all in there
        class(cluster2D_stream_distr_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        character(len=*), parameter   :: STK_FILETAB        = 'stkstreamtab.txt'
        character(len=:), allocatable :: WORK_PROJFILE
        integer,          parameter   :: SHIFTSRCH_PTCLSLIM = 2000  ! # of ptcls required to turn on shift search
        integer,          parameter   :: SHIFTSRCH_ITERLIM  = 5     ! # of iterations prior to turn on shift search
        integer,          parameter   :: CCRES_NPTCLS_LIM   = 10000 ! # of ptcls required to turn on objfun=ccres
        integer,          parameter   :: WAIT_WATCHER       = 60    ! seconds prior to new stack detection
        integer,          parameter   :: MAXNCLS            = 1000  ! maximum # of classes
        type(cluster2D_distr_commander)     :: xcluster2D_distr
        type(make_cavgs_distr_commander)    :: xmake_cavgs
        type(scale_project_distr_commander) :: xscale_distr
        type(cmdline)                       :: cline_cluster2D, cline_scale, cline_make_cavgs
        type(params)                        :: p_master
        type(sp_project)                    :: orig_proj, work_proj
        character(len=:),       allocatable :: orig_projfile
        character(len=STDLEN)               :: str_iter, refs_glob
        real    :: orig_smpd, msk, scale_factor, orig_msk
        integer :: iter, n_newstks, orig_box, box, nptcls_glob, ncls_glob, tnow
        integer :: nptcls_glob_prev, ncls_glob_prev, last_injection, n_stk, n_stk_prev
        logical :: do_autoscale

        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        ! make master parameters
        call cline%set('stream','yes') ! only for parameters determination
        p_master      = params(cline)
        if( .not.file_exists(p_master%projfile) )stop 'Project does not exist!'
        call cline%set('stream','no') ! was only for parameters determination
        orig_projfile = p_master%projfile
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
        call cline_make_cavgs%set('refs',   'cavgs_final'//p_master%ext)
        call cline_make_cavgs%delete('autoscale')
        call cline_make_cavgs%set('projfile', trim(WORK_PROJFILE)) ! TO WORK OUT
        ! wait for the first stacks to trickle in...
        do
            if( .not.is_file_open(p_master%projfile_target) )then
                call orig_proj%read(p_master%projfile_target)
                n_stk = orig_proj%os_stk%get_noris()
                if( n_stk > 0)then
                    nptcls_glob = orig_proj%get_nptcls()
                    if( nptcls_glob > p_master%ncls_start * p_master%nptcls_per_cls )then
                        ! Enough particles to start 2D clustering...
                        exit
                    endif
                endif
            endif
            call simple_sleep(WAIT_WATCHER) ! parameter instead
        enddo
        ! transfer of general info
        call work_proj%read(p_master%projfile)
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
        orig_msk  = p_master%msk
        p_master%smpd_targets2D(1) = max(orig_smpd, p_master%lp*LP2SMPDFAC2D)
        do_autoscale = p_master%autoscale.eq.'yes' .and. p_master%smpd_targets2D(1) > orig_smpd
        if( do_autoscale )then
            deallocate(WORK_PROJFILE)
            call work_proj%scale_projfile(p_master%smpd_targets2D(1), WORK_PROJFILE, cline_cluster2D, cline_scale)
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
        ncls_glob      = nint(real(nptcls_glob) / real(p_master%nptcls_per_cls))
        last_injection = simple_gettime()
        do iter = 1, 999
            str_iter = int2str_pad(iter,3)
            tnow     = simple_gettime()
            if(tnow-last_injection > p_master%time_inactive)then
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
            call cline_cluster2D%set('nparts',  real(min(ncls_glob,p_master%nparts)))
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
            refs_glob = trim(CAVGS_ITER_FBODY)//trim(str_iter)//trim(p_master%ext)
            ! remap zero-population classes
            call remap_empty_classes
            ! user termination
            if( file_exists(trim(TERM_STREAM)) )then
                    write(*,'(A)')'>>> TERMINATING CLUSTER2D STREAM'
                    exit
                endif
            ! detects whether new images have been added to the original project
            if( .not.is_file_open(p_master%projfile_target) )then
                call orig_proj%read_segment('stk', p_master%projfile_target)
                n_stk     = orig_proj%os_stk%get_noris()
                n_newstks = n_stk - n_stk_prev
            endif
            if(n_newstks > 0)then
                call work_proj%kill
                call work_proj%read(p_master%projfile)
                call append_newstacks( orig_proj%os_stk )
                ! updates # of classes
                ncls_glob_prev = ncls_glob
                ncls_glob      = nint(real(nptcls_glob) / p_master%nptcls_per_cls)
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
        call qsys_cleanup(p_master)
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
                character(len=STDLEN), allocatable :: new_stks(:)
                character(len=:),      allocatable :: stkname, val, ctfstr, phaseplate
                integer :: istk, nptcls, iptcl
                allocate(new_stks(n_newstks))
                do istk=1,n_newstks
                    call os_stk_in%getter(n_stk+istk, 'stk', stkname)
                    new_stks(istk) = trim(stkname)
                enddo
                ! scale stacks to append
                if( p_master%autoscale.eq.'yes' )then
                    call qenv%new(p_master)
                    call cline_scale%delete('stk')
                    call cline_scale%delete('stktab')
                    call cline_scale%set('dir_target', trim(SCSTK_DIR))
                    call cline_scale%set('filetab', trim(STK_FILETAB))
                    call cline_scale%set('stream', 'yes')
                    call cline_scale%set('numlen', 1.)
                    call write_filetable(SCALE_FILETAB, new_stks)
                    call qenv%exec_simple_prg_in_queue(cline_scale, 'OUT1','JOB_FINISHED_1')
                    call qsys_cleanup(p_master)
                    call del_file(SCALE_FILETAB)
                    do istk=1,n_newstks
                        new_stks(istk) = SCSTK_DIR//add2fbody(new_stks(istk), p_master%ext, SCALE_SUFFIX)
                    enddo
                endif
                ! append to working project
                do istk=1,n_newstks
                    o_stk  = os_stk_in%get_ori(n_stk+istk)
                    nptcls = nint(o_stk%get('nptcls'))
                    call os_ptcls%new_clean(nptcls)
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
                        stk = add2fbody(trim(refs_glob),p_master%ext,'_even')
                        call img_cavg%read(stk, fromtocls(icls, 1))
                        call img_cavg%write(stk, fromtocls(icls, 2))
                        stk = add2fbody(trim(refs_glob),p_master%ext,'_odd')
                        call img_cavg%read(stk, fromtocls(icls, 1))
                        call img_cavg%write(stk, fromtocls(icls, 2))
                    enddo
                    ! stack size preservation
                    call img_cavg%read(trim(refs_glob), ncls_glob)
                    call img_cavg%write(trim(refs_glob), ncls_glob)
                    stk = add2fbody(trim(refs_glob),p_master%ext,'_even')
                    call img_cavg%read(stk, ncls_glob)
                    call img_cavg%write(stk, ncls_glob)
                    stk = add2fbody(trim(refs_glob),p_master%ext,'_odd')
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
                            stk = add2fbody(trim(refs_glob),p_master%ext,'_even')
                            call img_cavg%read(stk, fromtocls(icls, 1))
                            call img_cavg%write(stk, fromtocls(icls, 2))
                            stk = add2fbody(trim(refs_glob),p_master%ext,'_odd')
                            call img_cavg%read(stk, fromtocls(icls, 1))
                            call img_cavg%write(stk, fromtocls(icls, 2))
                        enddo
                        ! stack size preservation
                        call img_cavg%read(trim(refs_glob), ncls_glob)
                        call img_cavg%write(trim(refs_glob), ncls_glob)
                        stk = add2fbody(trim(refs_glob),p_master%ext,'_even')
                        call img_cavg%read(stk, ncls_glob)
                        call img_cavg%write(stk, ncls_glob)
                        stk = add2fbody(trim(refs_glob),p_master%ext,'_odd')
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
