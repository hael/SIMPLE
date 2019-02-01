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
#include "simple_local_flags.inc"

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

contains

    subroutine exec_make_cavgs( self, cline )
        use simple_classaverager
        class(make_cavgs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: ncls_here
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params)
        write(logfhandle,'(a)') '>>> GENERATING CLUSTER CENTERS'
        ! deal with the orientations
        ncls_here = build%spproj_field%get_n('class')
        if( .not. cline%defined('ncls') ) params%ncls = build%spproj_field%get_n('class')
        if( params%l_remap_cls )then
            call build%spproj_field%remap_cls()
            if( cline%defined('ncls') )then
                if( params%ncls < ncls_here ) THROW_HARD('inputted ncls < ncls_in_oritab not allowed!')
                if( params%ncls > ncls_here )then
                    call build%spproj_field%expand_classes(params%ncls)
                endif
            endif
        else if( params%tseries .eq. 'yes' )then
            if( .not. cline%defined('ncls') )then
                THROW_HARD('# class averages (ncls) need to be part of command line when tseries=yes')
            endif
            call build%spproj_field%ini_tseries(params%ncls, 'class')
            call build%spproj_field%partition_eo(tseries=.true.)
        endif
        ! shift multiplication
        if( params%mul > 1. )then
            call build%spproj_field%mul_shifts(params%mul)
        endif
        ! setup weights
        if( params%softpw2D.eq.'yes' )then
            call build%spproj_field%calc_soft_weights_specscore
        else
            call build%spproj_field%calc_hard_weights2D(params%frac, params%ncls)
        endif
        ! shell weighted classes
        if( params%shellw.eq.'yes' )then
            call build%spproj_field%calc_bfac_rec_specscore(params%bfac_sdev)
        else
            call build%spproj_field%set_all2single('bfac_rec', 0.)
        endif
        ! even/odd partitioning
        if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
        ! write
        if( nint(cline%get_rarg('part')) .eq. 1 )then
            call build%spproj%write_segment_inside(params%oritype)
        endif
        ! create class averager
        call cavger_new('class')
        ! transfer ori data to object
        call cavger_transf_oridat(build%spproj)
        ! standard cavg assembly
        call cavger_assemble_sums( .false. )
        ! write sums
        call cavger_readwrite_partial_sums('write')
        call qsys_job_finished(  'simple_commander_cluster2D :: exec_make_cavgs' )
        call cavger_kill
        ! end gracefully
        call simple_end('**** SIMPLE_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_make_cavgs

     subroutine exec_cluster2D( self, cline )
        use simple_strategy2D_matcher, only: cluster2D_exec
        class(cluster2D_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: startit, ncls_from_refs, lfoo(3)
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params)
        if( cline%defined('refs') )then
            call find_ldim_nptcls(params%refs, lfoo, ncls_from_refs)
            ! consistency check
            if( params%ncls /=  ncls_from_refs ) THROW_HARD('nrefs /= inputted ncls')
        endif
        startit = 1
        if( cline%defined('startit') )startit = params%startit
        if( startit == 1 )call build%spproj_field%clean_updatecnt
        ! execute
        if( .not. cline%defined('outfile') ) THROW_HARD('need unique output file for parallel jobs')
        call cluster2D_exec( cline, startit ) ! partition or not, depending on 'part'
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
        call qsys_job_finished('simple_commander_cluster2D :: exec_cluster2D')
     end subroutine exec_cluster2D

    subroutine exec_cavgassemble( self, cline )
        use simple_classaverager
        class(cavgassemble_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params)
        call cavger_new( 'class')
        call cavger_transf_oridat( build%spproj )
        call cavger_assemble_sums_from_parts()
        if( cline%defined('which_iter') )then
            params%refs      = 'cavgs_iter'//int2str_pad(params%which_iter,3)//params%ext
            params%refs_even = 'cavgs_iter'//int2str_pad(params%which_iter,3)//'_even'//params%ext
            params%refs_odd  = 'cavgs_iter'//int2str_pad(params%which_iter,3)//'_odd'//params%ext
        else if( .not. cline%defined('refs') )then
            params%refs      = 'start2Drefs'//params%ext
            params%refs_even = 'start2Drefs_even'//params%ext
            params%refs_odd  = 'start2Drefs_odd'//params%ext
        endif
        call cavger_calc_and_write_frcs_and_eoavg(params%frcs, params%l_locres)
        ! classdoc gen needs to be after calc of FRCs
        call cavger_gen2Dclassdoc(build%spproj)
        ! write references
        call cavger_write(trim(params%refs),      'merged')
        call cavger_write(trim(params%refs_even), 'even'  )
        call cavger_write(trim(params%refs_odd),  'odd'   )
        call cavger_kill()
        ! write project
        call build%spproj%write_segment_inside('cls2D', params%projfile)
        ! end gracefully
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        call simple_touch('CAVGASSEMBLE_FINISHED', errmsg='In: commander_rec :: eo_cavgassemble ')
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
                call clsdoc_ranked%set(iclass, 'class',     real(order(iclass)))
                call clsdoc_ranked%set(iclass, 'rank',      real(iclass))
                call clsdoc_ranked%set(iclass, 'pop',       build%spproj_field%get(order(iclass),  'pop'))
                call clsdoc_ranked%set(iclass, 'res',       build%spproj_field%get(order(iclass),  'res'))
                call clsdoc_ranked%set(iclass, 'corr',      build%spproj_field%get(order(iclass), 'corr'))
                call clsdoc_ranked%set(iclass, 'w',         build%spproj_field%get(order(iclass),    'w'))
                call clsdoc_ranked%set(iclass, 'specscore', build%spproj_field%get(order(iclass), 'specscore'))
                write(logfhandle,'(a,1x,i5,1x,a,1x,i5,1x,a,i5,1x,a,1x,f6.2)') 'CLASS:', order(iclass),&
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
