! concrete commander: high-level workflows
module simple_commander_cluster3D
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters
use simple_sp_project,     only: sp_project
use simple_qsys_funs
implicit none

public :: cluster3D_commander
public :: cluster3D_refine_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: cluster3D_commander
  contains
    procedure :: execute      => exec_cluster3D
end type cluster3D_commander
type, extends(commander_base) :: cluster3D_refine_commander
  contains
    procedure :: execute      => exec_cluster3D_refine
end type cluster3D_refine_commander

contains

    subroutine exec_cluster3D( self, cline )
        use simple_oris,               only: oris
        use simple_sym,                only: sym
        use simple_cluster_seed,       only: gen_labelling
        use simple_commander_refine3D, only: refine3D_commander_distr
        use simple_commander_rec,      only: reconstruct3D_commander_distr
        class(cluster3D_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        ! constants
        integer,           parameter :: MAXITS1        = 50
        integer,           parameter :: MAXITS2        = 40
        character(len=*),  parameter :: one            = '01'
        character(len=12), parameter :: cls3D_projfile = 'cls3D.simple'
        ! distributed commanders
        type(refine3D_commander_distr)         :: xrefine3D_distr
        type(reconstruct3D_commander_distr)    :: xreconstruct3D_distr
        ! command lines
        type(cmdline)                          :: cline_refine3D1, cline_refine3D2
        type(cmdline)                          :: cline_reconstruct3D_mixed_distr
        type(cmdline)                          :: cline_reconstruct3D_multi_distr
        ! other variables
        type(parameters)                       :: params
        type(sym)                              :: symop
        type(sp_project)                       :: spproj, work_proj
        type(oris)                             :: os, opeaks
        type(ctfparams)                        :: ctfparms
        character(len=:),          allocatable :: cavg_stk, orig_projfile, prev_vol, target_name
        character(len=LONGSTRLEN), allocatable :: list(:)
        real,                      allocatable :: corrs(:), x(:), z(:), res(:), tmp_rarr(:)
        integer,                   allocatable :: labels(:), states(:), tmp_iarr(:)
        real     :: trs, extr_init, lp_cls3D, smpdfoo
        integer  :: i, iter, startit, rename_stat, ncls, boxfoo, iptcl, ipart
        integer  :: nptcls_part, istate, n_nozero
        logical  :: fall_over, cavgs_import
        if( nint(cline%get_rarg('nstates')) <= 1 ) THROW_HARD('Non-sensical NSTATES argument for heterogeneity analysis!')
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('refine') ) call cline%set('refine', 'cluster')
        if( .not. cline%defined('oritype')) call cline%set('oritype', 'ptcl3D')
        ! make master parameters
        call params%new(cline)
        orig_projfile   = trim(params%projfile)
        params%projfile = trim(params%cwd)//'/'//trim(params%projname)//trim(METADATA_EXT)
        call cline%set('projfile',params%projfile)
        ! set mkdir to no
        call cline%set('mkdir', 'no')
        ! prep project
        cavgs_import = .false.
        fall_over    = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                call work_proj%read(params%projfile)
                fall_over = work_proj%get_nptcls() == 0
            case('cls3D')
                call spproj%read(params%projfile)
                fall_over = spproj%os_out%get_noris() == 0
        case DEFAULT
            write(logfhandle,*)'Unsupported ORITYPE; simple_commander_cluster3D::exec_cluster3D'
        end select
        if( fall_over ) THROW_HARD('no particles found! exec_cluster3D')
        if( params%oritype.eq.'ptcl3D' )then
            ! just splitting
            call work_proj%split_stk(params%nparts, dir=PATH_PARENT)
        else
            ! class-averages
            params%projfile = trim(cls3d_projfile)
            call cline%set('oritype', 'ptcl3D')
            call spproj%get_cavgs_stk(cavg_stk, ncls, ctfparms%smpd)
            cavgs_import = spproj%os_ptcl2D%get_noris() == 0
            if( cavgs_import )then
                ! start from import
                if(.not.params%l_lpset ) THROW_HARD('need LP=XXX for imported class-averages; cluster3D')
                lp_cls3D = params%lp
                allocate(states(ncls), source=1)
            else
                ! start from previous 2D
                states = nint(spproj%os_cls2D%get_all('state'))
                ! determines resolution limit
                if( params%l_lpset )then
                    lp_cls3D = params%lp
                else
                    tmp_rarr  = spproj%os_cls2D%get_all('res')
                    tmp_iarr  = nint(spproj%os_cls2D%get_all('state'))
                    res       = pack(tmp_rarr, mask=(tmp_iarr>0))
                    lp_cls3D  = median_nocopy(res)
                    deallocate(res, tmp_iarr, tmp_rarr)
                endif
                if(cline%defined('lpstop')) lp_cls3D = max(lp_cls3D, params%lpstop)
            endif
            if( count(states==0) .eq. ncls )then
                THROW_HARD('no class averages detected in project file: '//trim(params%projfile)//'; cluster3D')
            endif
            work_proj%projinfo = spproj%projinfo
            work_proj%compenv  = spproj%compenv
            if(spproj%jobproc%get_noris()  > 0) work_proj%jobproc = spproj%jobproc
            call work_proj%add_single_stk(trim(cavg_stk), ctfparms, spproj%os_cls3D)
            ! takes care of states
            call work_proj%os_ptcl3D%set_all('state', real(states))
            ! name change
            call work_proj%projinfo%delete_entry('projname')
            call work_proj%projinfo%delete_entry('projfile')
            call cline%set('projfile', trim(params%projfile))
            call cline%set('projname', trim(get_fbody(trim(params%projfile),trim('simple'))))
            call work_proj%update_projinfo(cline)
            ! splitting in CURRENT directory
            call work_proj%split_stk(params%nparts, dir=PATH_HERE)
            ! write
            call work_proj%write
        endif
        ! fetch project oris
        call work_proj%get_sp_oris('ptcl3D', os)
        ! wipe previous states
        labels = nint(os%get_all('state'))
        if( any(labels > 1) )then
            where(labels > 0) labels = 1
            call os%set_all('state', real(labels))
        endif
        deallocate(labels)

        ! e/o partition
        if( .not.params%l_lpset )then
            if( os%get_nevenodd() == 0 ) call os%partition_eo
        else
            call os%set_all2single('eo', -1.)
        endif
        if( trim(params%refine) .eq. 'sym' )then
            ! randomize projection directions with respect to symmetry
            symop = sym(params%pgrp)
            call symop%symrandomize(os)
            call symop%kill
        endif

        ! prepare command lines from prototype
        call cline%delete('refine')
        ! resolution limits
        if( trim(params%oritype).eq.'cls3D' )then
            params%l_lpset = .true.
            call cline%set('lp',lp_cls3D)
        else
            if(.not.cline%defined('lplim_crit'))call cline%set('lplim_crit', 0.5)
        endif
        cline_refine3D1                 = cline ! first stage, extremal optimization
        cline_refine3D2                 = cline ! second stage, stochastic refinement
        cline_reconstruct3D_mixed_distr = cline
        cline_reconstruct3D_multi_distr = cline
        ! first stage
        call cline_refine3D1%set('prg',       'refine3D')
        call cline_refine3D1%set('match_filt','no')
        call cline_refine3D1%set('maxits',     real(MAXITS1))
        call cline_refine3D1%set('neigh',     'yes') ! always consider neighbours
        if( .not.cline_refine3D1%defined('nnn') )then
            call cline_refine3D1%set('nnn', 0.05*real(params%nspace))
        endif
        call cline_refine3D1%delete('update_frac')  ! no update frac for extremal optimization
        ! second stage
        call cline_refine3D2%set('prg', 'refine3D')
        call cline_refine3D2%set('match_filt','no')
        call cline_refine3D2%set('refine', 'multi')
        if( .not.cline%defined('update_frac') )call cline_refine3D2%set('update_frac', 0.5)
        ! reconstructions
        call cline_reconstruct3D_mixed_distr%set('prg',    'reconstruct3D')
        call cline_reconstruct3D_mixed_distr%set('nstates', 1.)
        call cline_reconstruct3D_mixed_distr%delete('lp')
        call cline_reconstruct3D_multi_distr%set('prg', 'reconstruct3D')
        call cline_reconstruct3D_multi_distr%delete('lp')
        if( trim(params%refine) .eq. 'sym' )then
            call cline_reconstruct3D_multi_distr%set('pgrp','c1')
            call cline_reconstruct3D_mixed_distr%set('pgrp','c1')
        endif
        if( cline%defined('trs') )then
            ! all good
        else
            ! works out shift limits for in-plane search
            trs = MSK_FRAC*real(params%msk)
            trs = min(MAXSHIFT, max(MINSHIFT, trs))
            call cline_refine3D1%set('trs',trs)
            call cline_refine3D2%set('trs',trs)
        endif
        ! refinement specific section
        select case(trim(params%refine))
            case('sym')
                call cline_refine3D1%set('refine','clustersym')
                call cline_refine3D2%set('pgrp','c1')
                call cline_refine3D2%delete('neigh') ! no neighbour mode for symmetry
                call cline_refine3D2%delete('nnn')
            case DEFAULT
                call cline_refine3D1%set('refine', 'cluster')
        end select

        ! MIXED MODEL RECONSTRUCTION
        ! retrieve mixed model Fourier components, normalization matrix, FSC & anisotropic filter
        if( .not.params%l_lpset )then
            work_proj%os_ptcl3D = os
            call work_proj%write
            call xreconstruct3D_distr%execute(cline_reconstruct3D_mixed_distr)
            rename_stat = simple_rename(trim(VOL_FBODY)//one//params%ext, trim(CLUSTER3D_VOL)//params%ext)
            rename_stat = simple_rename(trim(VOL_FBODY)//one//'_even'//params%ext, trim(CLUSTER3D_VOL)//'_even'//params%ext)
            rename_stat = simple_rename(trim(VOL_FBODY)//one//'_odd'//params%ext,  trim(CLUSTER3D_VOL)//'_odd'//params%ext)
            rename_stat = simple_rename(trim(FSC_FBODY)//one//BIN_EXT, trim(CLUSTER3D_FSC))
        endif

        ! calculate extremal initial ratio
        if( os%isthere('corr') )then
            labels    = nint(os%get_all('state'))
            corrs     = os%get_all('corr')
            x         = pack(corrs, mask=(labels>0))
            z         = robust_z_scores(x)
            extr_init = 2.*real(count(z<-1.)) / real(count(labels>0))
            extr_init = max(0.1,extr_init)
            extr_init = min(extr_init,EXTRINITHRESH)
            deallocate(x,z,corrs,labels)
        else
            extr_init = EXTRINITHRESH
        endif
        call cline_refine3D1%set('extr_init', extr_init)
        write(logfhandle,'(A,F5.2)') '>>> INITIAL EXTREMAL RATIO: ',extr_init

        ! randomize state labels
        write(logfhandle,'(A)') '>>>'
        call gen_labelling(os, params%nstates, 'squared_uniform')
        work_proj%os_ptcl3D = os
        ! writes for reconstruct3D,refine3D
        call work_proj%write
        call work_proj%kill
        call os%kill

        ! STAGE1: extremal optimization, frozen orientation parameters
        write(logfhandle,'(A)')    '>>>'
        write(logfhandle,'(A,I3)') '>>> 3D CLUSTERING - STAGE 1'
        write(logfhandle,'(A)')    '>>>'
        call xrefine3D_distr%execute(cline_refine3D1)
        iter = nint(cline_refine3D1%get_rarg('endit'))
        ! for analysis purpose only
        call work_proj%read_segment('ptcl3D', params%projfile)
        call work_proj%kill

        ! STAGE2: multi-states refinement
        startit = iter + 1
        call cline_refine3D2%set('startit', real(startit))
        call cline_refine3D2%set('maxits',  real(min(params%maxits,startit+MAXITS2)))
        write(logfhandle,'(A)')    '>>>'
        write(logfhandle,'(A,I3)') '>>> 3D CLUSTERING - STAGE 2'
        write(logfhandle,'(A)')    '>>>'
        call xrefine3D_distr%execute(cline_refine3D2)

        ! class-averages mapping
        if( params%oritype.eq.'cls3D' )then
            call work_proj%read(params%projfile)
            spproj%os_cls3D = work_proj%os_ptcl3D
            if( cavgs_import )then
                ! no mapping
            else
                ! map to ptcl3D
                call spproj%map2ptcls
            endif
            call spproj%write
            call spproj%kill
            call del_file(cls3d_projfile)
        endif

        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER3D NORMAL STOP ****')
    end subroutine exec_cluster3D

    subroutine exec_cluster3D_refine( self, cline )
        use simple_oris,               only: oris
        use simple_ori,                only: ori
        use simple_parameters,         only: params_glob
        use simple_commander_refine3D, only: refine3D_commander_distr
        class(cluster3D_refine_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! constants
        integer,                     parameter :: MAXITS = 40
        character(len=12),           parameter :: cls3D_projfile = 'cls3D.simple'
        ! distributed commanders
        type(refine3D_commander_distr)         :: xrefine3D_distr
        ! command lines
        type(cmdline),             allocatable :: cline_refine3D(:)
        ! other variables
        integer,                   allocatable :: state_pops(:), states(:), master_states(:)
        character(len=STDLEN),     allocatable :: dirs(:), projfiles(:)
        character(len=LONGSTRLEN), allocatable :: rel_stks(:), stks(:)
        character(len=:),          allocatable :: projname, cavg_stk, frcs_fname, orig_projfile, stk
        type(parameters)         :: params
        type(ctfparams)          :: ctfparms
        type(sp_project)         :: spproj, spproj_master
        class(oris),     pointer :: pos => null()
        type(ori)                :: o_tmp
        integer                  :: state, iptcl, nstates, single_state, ncls, istk, nstks
        logical                  :: l_singlestate, cavgs_import, fall_over
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call params%new(cline)
        ! set mkdir to no
        call cline%set('mkdir', 'no')
        ! sanity checks
        if( .not.cline%defined('maxits') )call cline%set('maxits',real(MAXITS))
        l_singlestate = cline%defined('state')
        if( l_singlestate )then
            single_state = nint(cline%get_rarg('state'))
        else
            single_state = 0
        endif
        orig_projfile = trim(params%projfile)
        cavgs_import  = .false.
        fall_over     = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                call spproj_master%read(params%projfile)
                fall_over = spproj_master%get_nptcls() == 0
            case('cls3D')
                call spproj%read(params%projfile)
                fall_over = spproj%os_out%get_noris() == 0
        case DEFAULT
            write(logfhandle,*)'Unsupported ORITYPE; simple_commander_cluster3D::exec_cluster3D_refine'
        end select
        if( fall_over ) THROW_HARD('no particles found! exec_cluster3D_refine')
        ! stash states
        if(params%oritype.eq.'cls3D')then
            master_states  = nint(spproj%os_cls3D%get_all('state'))
            call spproj%os_cls3D%get_pops(state_pops, 'state', consider_w=.false.)
        else
            master_states  = nint(spproj_master%os_ptcl3D%get_all('state'))
            call spproj_master%os_ptcl3D%get_pops(state_pops, 'state', consider_w=.false.)
        endif
        nstates        = maxval(master_states)
        params%nstates = nstates
        if( params%nstates==1 )then
            THROW_HARD('non-sensical # states: '//int2str(params%nstates)//' for multi-particle refinement')
        endif
        if( state_pops(params%state) == 0 )then
            THROW_HARD('state: '//int2str(params%state)//' is empty')
        endif
        ! state dependent variables
        allocate(projfiles(params%nstates), dirs(params%nstates), cline_refine3D(params%nstates))
        do state = 1, params%nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. single_state.ne.state )cycle
            ! name & directory
            projname         = 'state_'//trim(int2str_pad(state,2))
            projfiles(state) = trim(projname)//trim(METADATA_EXT)
            dirs(state)      = trim(int2str(state))//'_refine3D'
            ! command line
            cline_refine3D(state) = cline
            call cline_refine3D(state)%set('prg',     'refine3D')
            call cline_refine3D(state)%set('projname',trim(projname))
            call cline_refine3D(state)%set('projfile',trim(projfiles(state)))
            call cline_refine3D(state)%set('mkdir',   'yes')
            call cline_refine3D(state)%set('refine',  'single')
            call cline_refine3D(state)%delete('state')
            call cline_refine3D(state)%delete('nstates')
            if(params%oritype.eq.'cls3D') call cline_refine3D(state)%set('oritype', 'ptcl3D')
        enddo

        ! transfer cavgs to ptcl3D
        if(params%oritype.eq.'cls3D')then
            call spproj%get_cavgs_stk(cavg_stk, ncls, ctfparms%smpd)
            states       = nint(spproj%os_cls3D%get_all('state'))
            cavgs_import = spproj%os_ptcl2D%get_noris() == 0
            if( cavgs_import )then
                ! start from import
                if(.not.cline%defined('lp')) THROW_HARD('need LP=XXX for imported class-averages; cluster3D_refine')
            else
                call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
                if( .not.file_exists(frcs_fname) )then
                    THROW_HARD('the project file does not contain enough information for e/o alignment, use a low-pass instead')
                endif
            endif
            if( count(states==0) .eq. ncls )then
                THROW_HARD('no class averages detected in project file: '//trim(params%projfile)//'; cluster3D_refine')
            endif
            spproj_master%projinfo = spproj%projinfo
            spproj_master%compenv  = spproj%compenv
            if(spproj%jobproc%get_noris()  > 0) spproj_master%jobproc = spproj%jobproc
            if( cavgs_import )then
                call spproj_master%add_single_stk(trim(cavg_stk), ctfparms, spproj%os_cls3D)
                call spproj_master%os_ptcl3D%set_all('state', real(states))
            else
                call prep_eo_stks
                params_glob%nptcls = spproj_master%get_nptcls()
            endif
            ! name & oritype change
            call spproj_master%projinfo%delete_entry('projname')
            call spproj_master%projinfo%delete_entry('projfile')
            call cline%set('projfile', cls3D_projfile)
            call cline%set('projname', trim(get_fbody(trim(cls3D_projfile),trim('simple'))))
            call spproj_master%update_projinfo(cline)
            ! splitting in CURRENT directory
            call spproj_master%split_stk(params%nparts, dir=PATH_HERE)
            ! write
            call spproj_master%write
        endif

        ! states are lost from the project after this loop and stored in master_states
        nstks = spproj_master%os_stk%get_noris()
        allocate(stks(nstks), rel_stks(nstks))
        do istk = 1,nstks
            stk = spproj_master%get_stkname(istk)
            stks(istk) = NIL
            if( file_exists(stk) )then
                rel_stks(istk) = trim(stk)
                stks(istk)     = simple_abspath(stk)
            endif
            ! turns to absolute paths
            call spproj_master%os_stk%set(istk,'stk',stks(istk))
        enddo
        do state = 1, params%nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. single_state.ne.state )cycle
            ! states
            states = master_states
            where(states /= state) states = 0
            where(states /= 0)     states = 1
            call spproj_master%os_ptcl3D%set_all('state', real(states))
            ! write
            call spproj_master%update_projinfo(cline_refine3D(state))
            call spproj_master%write(projfiles(state))
            deallocate(states)
        enddo
        do istk = 1,nstks
            stk = spproj_master%get_stkname(istk)
            if( trim(stks(istk)) /= NIL )then
                ! restores path
                call spproj_master%os_stk%set(istk,'stk',rel_stks(istk))
            endif
        enddo
        deallocate(stks,rel_stks)
        ! restores name
        call spproj_master%update_projinfo(cline)

        ! Execute individual refine3D jobs
        do state = 1, nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. state.ne.single_state )cycle
            write(logfhandle,'(A)')   '>>>'
            write(logfhandle,'(A,I2,A,A)')'>>> REFINING STATE: ', state
            write(logfhandle,'(A)')   '>>>'
            params_glob%projname = 'state_'//trim(int2str_pad(state,2))
            params_glob%projfile = projfiles(state)
            params_glob%nstates = 1
            params_glob%state   = 1
            call xrefine3D_distr%execute(cline_refine3D(state))
            call simple_chdir(PATH_PARENT,errmsg="commander_hlev_wflows :: exec_cluster3D_refine;")
        enddo
        ! restores original values
        params_glob%projname = trim(get_fbody(trim(orig_projfile),trim('simple')))
        params_glob%projfile = trim(orig_projfile)
        params_glob%nstates  = nstates

        ! consolidates new orientations parameters & files
        ! gets original project back
        if(params%oritype.eq.'cls3D')then
            params_glob%nptcls = ncls
            call spproj_master%kill
            call spproj_master%read(params%projfile)
        endif
        do state=1,params%nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. state.ne.single_state ) cycle
            ! renames volumes and updates in os_out
            call stash_state(state)
            ! updates orientations
            call spproj%read_segment('ptcl3D',filepath(dirs(state),projfiles(state)))
            call spproj_master%ptr2oritype(params%oritype, pos)
            do iptcl=1,params%nptcls
                if( master_states(iptcl)==state )then
                    call spproj%os_ptcl3D%get_ori(iptcl, o_tmp)
                    call pos%set_ori(iptcl, o_tmp)
                    ! reset original states
                    call pos%set(iptcl,'state',real(state))
                endif
            enddo
            call spproj%kill
        enddo
        ! map to ptcls for non-imported class-averages
        if(params%oritype.eq.'cls3D' .and. .not.cavgs_import) call spproj_master%map2ptcls
        ! final write
        call spproj_master%write
        ! cleanup
        call spproj%kill
        call spproj_master%kill
        do state=1,params%nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. state.ne.single_state )cycle
            call del_file(projfiles(state))
        enddo
        if(params%oritype.eq.'cls3D') call del_file(cls3D_projfile)
        deallocate(master_states, dirs, projfiles)
        call o_tmp%kill
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER3D_REFINE NORMAL STOP ****')

        contains

            ! stash docs, volumes , etc.
            subroutine stash_state(s)
                integer, intent(in) :: s
                character(len=2),            parameter :: one = '01'
                character(len=LONGSTRLEN), allocatable :: files(:)
                character(len=LONGSTRLEN) :: src, dest, vol, fsc!, volfilt
                character(len=2)          :: str_state
                character(len=8)          :: str_iter
                integer                   :: i, final_it, stat, l, pos
                final_it  = nint(cline_refine3D(s)%get_rarg('endit'))
                str_state = int2str_pad(s,2)
                str_iter  = '_ITER'//int2str_pad(final_it,3)
                if( s == 1 )then
                    vol     = filepath(dirs(s), trim(VOL_FBODY)//one//trim(params%ext))
                    fsc     = filepath(dirs(s), trim(FSC_FBODY)//one//BIN_EXT)
                else
                    ! renames all *state01* files
                     call simple_list_files(trim(dirs(s))//'/*state01*', files)
                     do i=1,size(files)
                         src  = files(i)
                         dest = basename(files(i))
                         l    = len_trim(dest)
                         pos  = index(dest(1:l),'_state01',back=.true.)
                         dest(pos:pos+7) = '_state' // str_state
                         stat = simple_rename(src, filepath(trim(dirs(s)),dest))
                     enddo
                     deallocate(files)
                     call simple_list_files(trim(dirs(s))//'/*STATE01*', files)
                     do i=1,size(files)
                         src  = files(i)
                         dest = basename(files(i))
                         l    = len_trim(dest)
                         pos  = index(dest(1:l),'_STATE01',back=.true.)
                         dest(pos:pos+7) = '_STATE' // str_state
                         stat = simple_rename(src, filepath(trim(dirs(s)),dest))
                     enddo
                     deallocate(files)
                     vol     = filepath(dirs(s), trim(VOL_FBODY)//str_state//trim(params%ext))
                     fsc     = filepath(dirs(s), trim(FSC_FBODY)//str_state//BIN_EXT)
                endif
                ! updates os_out
                if(params%oritype.eq.'cls3D')then
                    call spproj%add_vol2os_out(vol, params%smpd, s, 'vol_cavg')
                else
                    call spproj_master%add_vol2os_out(vol, params%smpd, s, 'vol')
                endif
                call spproj_master%add_fsc2os_out(fsc, s, params%box)
                call spproj_master%add_frcs2os_out(filepath(dirs(s),FRCS_FILE),'frc3D')
            end subroutine stash_state

            subroutine prep_eo_stks
                use simple_ori, only: ori
                type(ori)                     :: o, o_even, o_odd
                character(len=:), allocatable :: eostk, ext
                integer :: even_ind, odd_ind, state, icls
                do state=1,params%nstates
                    call cline_refine3D(state)%delete('lp')
                    call cline_refine3D(state)%set('frcs',trim(frcs_fname))
                    call cline_refine3D(state)%set('lplim_crit', 0.5)
                    call cline_refine3D(state)%set('clsfrcs',   'yes')
                enddo
                ! add stks
                ext   = '.'//fname2ext( cavg_stk )
                eostk = add2fbody(trim(cavg_stk), trim(ext), '_even')
                call spproj_master%add_stk(eostk, ctfparms)
                eostk = add2fbody(trim(cavg_stk), trim(ext), '_odd')
                call spproj_master%add_stk(eostk, ctfparms)
                ! update orientations parameters
                if(allocated(master_states))deallocate(master_states)
                allocate(master_states(2*ncls), source=0)
                do icls=1,ncls
                    even_ind = icls
                    odd_ind  = ncls+icls
                    call spproj%os_cls3D%get_ori(icls, o)
                    state    = spproj%os_cls3D%get_state(icls)
                    call o%set('class', real(icls)) ! for mapping frcs in 3D
                    call o%set('state', real(state))
                    ! even
                    o_even = o
                    call o_even%set('eo', 0.)
                    call o_even%set('stkind', spproj_master%os_ptcl3D%get(even_ind,'stkind'))
                    call spproj_master%os_ptcl3D%set_ori(even_ind, o_even)
                    master_states(even_ind) = state
                    ! odd
                    o_odd = o
                    call o_odd%set('eo', 1.)
                    call o_odd%set('stkind', spproj_master%os_ptcl3D%get(odd_ind,'stkind'))
                    call spproj_master%os_ptcl3D%set_ori(odd_ind, o_odd)
                    master_states(odd_ind) = state
                enddo
                ! cleanup
                deallocate(eostk, ext)
                call o%kill
                call o_even%kill
                call o_odd%kill
            end subroutine prep_eo_stks

    end subroutine exec_cluster3D_refine

end module simple_commander_cluster3D
