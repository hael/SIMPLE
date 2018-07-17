! concrete commander: cluster2D for simultanous 2D alignment and clustering of single-particle images
module simple_commander_cluster2D
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_qsys_funs,      only: qsys_job_finished
use simple_parameters,     only: parameters
use simple_builder,        only: builder
implicit none

public :: make_cavgs_commander
public :: cluster2D_commander
public :: cavgassemble_commander
public :: check_2Dconv_commander
public :: rank_cavgs_commander
public :: cluster_cavgs_commander
private

type, extends(commander_base) :: make_cavgs_commander
  contains
    procedure :: execute      => exec_make_cavgs
end type make_cavgs_commander
type, extends(commander_base) :: cluster2D_commander
  contains
    procedure :: execute      => exec_cluster2D
end type cluster2D_commander
type, extends(commander_base) :: cavgassemble_commander
  contains
    procedure :: execute      => exec_cavgassemble
end type cavgassemble_commander
type, extends(commander_base) :: check_2Dconv_commander
  contains
    procedure :: execute      => exec_check_2Dconv
end type check_2Dconv_commander
type, extends(commander_base) :: rank_cavgs_commander
  contains
    procedure :: execute      => exec_rank_cavgs
end type rank_cavgs_commander
type, extends(commander_base) :: cluster_cavgs_commander
  contains
    procedure :: execute      => exec_cluster_cavgs
end type cluster_cavgs_commander

#include "simple_local_flags.inc"

contains

    subroutine exec_make_cavgs( self, cline )
        use simple_classaverager
        class(make_cavgs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: ncls_here
        integer(timer_int_kind)::t1, t2
        t1=tic()
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params)
        DebugPrint ' exec_make_cavgs init and build                              ', toc(t1)
        t2=tic()
        write(*,'(a)') '>>> GENERATING CLUSTER CENTERS'
        ! deal with the orientations
        ncls_here = build%spproj_field%get_n('class')
        if( .not. cline%defined('ncls') ) params%ncls = build%spproj_field%get_n('class')
         t2=tic()
        if( params%l_remap_cls )then
            call build%spproj_field%remap_cls()
            DebugPrint ' exec_make_cavgs remap_cls                               ', toc(t1)
            t2=tic()
            if( cline%defined('ncls') )then
                if( params%ncls < ncls_here ) stop 'ERROR, inputted ncls < ncls_in_oritab; not allowed!'
                if( params%ncls > ncls_here )then
                    call build%spproj_field%expand_classes(params%ncls)
                endif
            endif
            DebugPrint ' exec_make_cavgs expand_classes                          ', toc(t1)
        else if( params%tseries .eq. 'yes' )then
            if( .not. cline%defined('ncls') )then
                stop '# class averages (ncls) need to be part of command line when tseries=yes'
            endif
            call build%spproj_field%ini_tseries(params%ncls, 'class')
            DebugPrint ' exec_make_cavgs init tseries                            ', toc(t1)
            t2=tic()
            call build%spproj_field%partition_eo(tseries=.true.)
            DebugPrint ' exec_make_cavgs partition_eo tseries                    ', toc(t1)
        endif
        ! shift multiplication
        if( params%mul > 1. )then
            t2=tic()
            call build%spproj_field%mul_shifts(params%mul)
             DebugPrint ' exec_make_cavgs mul_shifts                             ', toc(t2)
        endif
        ! setup weights in case the 2D was run without them (specscore will still be there)
        if( params%weights2D.eq.'yes' )then
            call build%spproj_field%calc_spectral_weights
        else
            t2=tic()
            call build%spproj_field%set_all2single('w', 1.0)
            if( params%shellw.eq.'yes' ) call build%spproj_field%calc_bfac_rec
             DebugPrint ' exec_make_cavgs spproj_field%calc_bfac_rec             ', toc(t2)
        endif
        ! even/odd partitioning
        t2=tic()
        if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
        DebugPrint ' exec_make_cavgs spproj_field%partition_eo                  ', toc(t2)
        ! write
        if( nint(cline%get_rarg('part')) .eq. 1 )then
            call build%spproj%write_segment_inside(params%oritype)
        endif
        t2=tic()
        ! create class averager
        call cavger_new('class')
        DebugPrint ' exec_make_cavgs cavger                                     ', toc(t2)
        ! transfer ori data to object
        call cavger_transf_oridat(build%spproj)
        DebugPrint ' exec_make_cavgs cavger transf_oridat                       ', toc(t2)
        ! standard cavg assembly
        call cavger_assemble_sums( .false. )
        DebugPrint ' exec_make_cavgs cavger assemble_sums                       ', toc(t2)
        ! write sums
        call cavger_readwrite_partial_sums('write')
        DebugPrint ' exec_make_cavgs cavger readwrite_partial_sums              ', toc(t2)
        call qsys_job_finished(  'simple_commander_cluster2D :: exec_make_cavgs' )
        call cavger_kill
        ! end gracefully
        call simple_end('**** SIMPLE_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
        DebugPrint ' exec_make_cavgs total time                                 ', toc(t1)
    end subroutine exec_make_cavgs

     subroutine exec_cluster2D( self, cline )
        use simple_strategy2D_matcher, only: cluster2D_exec
        class(cluster2D_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: startit, ncls_from_refs, lfoo(3)
        integer(timer_int_kind)::t1, t2
        t1=tic()
        call cline%set('oritype', 'ptcl2D')
        t2=tic()
        call build%init_params_and_build_strategy2D_tbox(cline, params)
        DebugPrint ' exec_cluster2D init and build                              ', toc(t2)
        t2=tic()
        if( cline%defined('refs') )then
            call find_ldim_nptcls(params%refs, lfoo, ncls_from_refs)
            ! consistency check
            if( params%ncls /=  ncls_from_refs ) stop 'nrefs /= inputted ncls'
        endif
        startit = 1
        DebugPrint ' exec_cluster2D refs                                        ', toc(t2)
        t2=tic()
        if( cline%defined('startit') )startit = params%startit
        if( startit == 1 )call build%spproj_field%clean_updatecnt
        DebugPrint ' exec_cluster2D clean updatecnt                             ', toc(t2)
        ! execute
        if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
        t2=tic()
        call cluster2D_exec( cline, startit) ! partition or not, depending on 'part'
        DebugPrint ' exec_cluster2D exec time                                   ', toc(t2)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
        call qsys_job_finished('simple_commander_cluster2D :: exec_cluster2D')
        DebugPrint ' exec_cluster2D total time                                  ', toc(t1)
     end subroutine exec_cluster2D

    subroutine exec_cavgassemble( self, cline )
        use simple_classaverager
        use simple_strategy2D3D_common, only: gen2Dclassdoc
        class(cavgassemble_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer(timer_int_kind)::t1, t2
        t1=tic()
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params)
        DebugPrint ' exec_cavgassemble init and build                            ', toc(t1)
        call cavger_new( 'class')
        t2=tic()
        call cavger_assemble_sums_from_parts()
        DebugPrint ' exec_cavgassemble assemble_sums_from_parts                 ', toc(t2)
        if( cline%defined('which_iter') )then
            params%refs      = 'cavgs_iter'//int2str_pad(params%which_iter,3)//params%ext
            params%refs_even = 'cavgs_iter'//int2str_pad(params%which_iter,3)//'_even'//params%ext
            params%refs_odd  = 'cavgs_iter'//int2str_pad(params%which_iter,3)//'_odd'//params%ext
        else if( .not. cline%defined('refs') )then
            params%refs      = 'start2Drefs'//params%ext
            params%refs_even = 'start2Drefs_even'//params%ext
            params%refs_odd  = 'start2Drefs_odd'//params%ext
        endif
        t2=tic()
        call cavger_calc_and_write_frcs_and_eoavg(params%frcs)
        DebugPrint ' exec_cavgassemble assemble_sums_from_parts                  ', toc(t2)
        ! classdoc gen needs to be after calc of FRCs
        t2=tic()
        call gen2Dclassdoc
        DebugPrint ' exec_cavgassemble gen2Dclassdoc                            ', toc(t2)
        ! write references
        t2=tic()
        call cavger_write(trim(params%refs),      'merged')
        call cavger_write(trim(params%refs_even), 'even'  )
        call cavger_write(trim(params%refs_odd),  'odd'   )
        DebugPrint ' exec_cavgassemble write                                     ', toc(t2)
        call cavger_kill()
        ! write project
        t2=tic()
        call build%spproj%write_segment_inside('cls2D', params%projfile)
        DebugPrint ' exec_cavgassemble write_segment_inside                     ', toc(t2)
        ! end gracefully
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        call simple_touch('CAVGASSEMBLE_FINISHED', errmsg='In: commander_rec :: eo_cavgassemble ')
         DebugPrint ' exec_cavgassemble Completed in                             ', toc(t1), ' secs'
    end subroutine exec_cavgassemble

    subroutine exec_check_2Dconv( self, cline )
        use simple_convergence, only: convergence
        use simple_parameters,  only: params_glob
        class(check_2Dconv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        type(convergence) :: conv
        logical :: converged
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        ! convergence check
        converged = conv%check_conv2D(cline, build%spproj_field%get_n('class'), params%msk)
        call cline%set('frac', conv%get('frac'))
        if( params_glob%l_doshift )then
            ! activates shift search
            call cline%set('trs', params_glob%trs)
        endif
        if( converged )then
            call cline%set('converged', 'yes')
        else
            call cline%set('converged', 'no')
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_2DCONV NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_2Dconv

    subroutine exec_rank_cavgs( self, cline )
        use simple_oris, only: oris
        class(rank_cavgs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        integer, allocatable :: order(:)
        real,    allocatable :: res(:)
        integer    :: ldim(3), ncls, iclass
        type(oris) :: clsdoc_ranked
        call cline%set('oritype', 'cls2D')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call find_ldim_nptcls(params%stk, ldim, ncls)
        params%ncls = ncls
        if( build%spproj_field%get_noris() == params%ncls )then
            ! all we need to do is fetch from classdoc in projfile &
            ! order according to resolution
            call clsdoc_ranked%new(params%ncls)
            res = build%spproj%os_cls2D%get_all('res')
            allocate(order(params%ncls))
            order = (/(iclass,iclass=1,params%ncls)/)
            call hpsort(res, order)
            do iclass=1,params%ncls
                call clsdoc_ranked%set(iclass, 'class', real(order(iclass)))
                call clsdoc_ranked%set(iclass, 'rank',  real(iclass))
                call clsdoc_ranked%set(iclass, 'pop',   build%spproj_field%get(order(iclass),  'pop'))
                call clsdoc_ranked%set(iclass, 'res',   build%spproj_field%get(order(iclass),  'res'))
                call clsdoc_ranked%set(iclass, 'corr',  build%spproj_field%get(order(iclass), 'corr'))
                call clsdoc_ranked%set(iclass, 'w',     build%spproj_field%get(order(iclass),    'w'))
                write(*,'(a,1x,i5,1x,a,1x,i5,1x,a,i5,1x,a,1x,f6.2)') 'CLASS:', order(iclass),&
                    &'RANK:', iclass ,'POP:', nint(build%spproj_field%get(order(iclass), 'pop')),&
                    &'RES:', build%spproj_field%get(order(iclass), 'res')
                call build%img%read(params%stk, order(iclass))
                call build%img%write(params%outstk, iclass)
            end do
            call clsdoc_ranked%write('classdoc_ranked.txt')
        else
            ! nothing to do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_RANK_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_rank_cavgs

    subroutine exec_cluster_cavgs( self, cline )
        use simple_cluster_cavgs, only: cluster_cavgs_exec
        class(cluster_cavgs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call cluster_cavgs_exec( )
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_cluster_cavgs

end module simple_commander_cluster2D
