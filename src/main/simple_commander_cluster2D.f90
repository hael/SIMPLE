! concrete commander: cluster2D for simultanous 2D alignment and clustering of single-particle images
module simple_commander_cluster2D
#include "simple_lib.f08"
use simple_cmdline,             only: cmdline
use simple_params,              only: params
use simple_build,               only: build
use simple_commander_base,      only: commander_base
use simple_imghead,             only: find_ldim_nptcls
use simple_strategy2D3D_common, only: gen2Dclassdoc
use simple_qsys_funs,           only: qsys_job_finished
use simple_projection_frcs,     only: projection_frcs
use simple_binoris_io           ! use all in there
implicit none

public :: make_cavgs_commander
public :: cluster2D_commander
public :: cavgassemble_commander
public :: check2D_conv_commander
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
type, extends(commander_base) :: check2D_conv_commander
  contains
    procedure :: execute      => exec_check2D_conv
end type check2D_conv_commander
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
        class(cmdline),             intent(inout) :: cline
        type(params)  :: p
        type(build)   :: b
        integer :: ncls_here, icls, fnr, file_stat, j
        p = params(cline, spproj_a_seg=PTCL2D_SEG)        ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        call b%build_strategy2D_tbox(p)                   ! 2D Hadamard matcher built
        write(*,'(a)') '>>> GENERATING CLUSTER CENTERS'
        ! deal with the orientations
        ncls_here = b%a%get_n('class')
        if( .not. cline%defined('ncls') ) p%ncls = b%a%get_n('class')
        if( p%l_remap_cls )then
            call b%a%remap_cls()
            if( cline%defined('ncls') )then
                if( p%ncls < ncls_here ) stop 'ERROR, inputted ncls < ncls_in_oritab; not allowed!'
                if( p%ncls > ncls_here )then
                    call b%a%expand_classes(p%ncls)
                endif
            endif
        else if( p%tseries .eq. 'yes' )then
            if( .not. cline%defined('ncls') )then
                stop '# class averages (ncls) need to be part of command line when tseries=yes'
            endif
            call b%a%ini_tseries(p%ncls, 'class')
            call b%a%partition_eo(tseries=.true.)
        endif
        ! shift multiplication
        if( p%mul > 1. )then
            call b%a%mul_shifts(p%mul)
        endif
        ! setup weights in case the 2D was run without them (specscore will still be there)
        if( p%weights2D.eq.'yes' )then
            if( p%nptcls <= SPECWMINPOP )then
                call b%a%set_all2single('w', 1.0)
            else
                ! frac is one by default in cluster2D (no option to set frac)
                ! so spectral weighting is done over all images
                call b%a%calc_spectral_weights(1.0)
            endif
        else
            call b%a%set_all2single('w', 1.0)
        endif
        ! even/odd partitioning
        if( b%a%get_nevenodd() == 0 ) call b%a%partition_eo
        ! write
        if( nint(cline%get_rarg('part')) .eq. 1 ) call b%spproj%write()
        ! create class averager
        call cavger_new(b, p, 'class')
        ! transfer ori data to object
        call cavger_transf_oridat(b%a)
        ! standard cavg assembly
        call cavger_assemble_sums( .false. )
        ! write sums
        call cavger_readwrite_partial_sums('write')
        call qsys_job_finished( p, 'simple_commander_cluster2D :: exec_make_cavgs' )
        call cavger_kill
        ! end gracefully
        call simple_end('**** SIMPLE_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_make_cavgs

    subroutine exec_cluster2D( self, cline )
        use simple_strategy2D_matcher, only: cluster2D_exec
        class(cluster2D_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: i, startit, ncls_from_refs, lfoo(3)
        p = params(cline, spproj_a_seg=PTCL2D_SEG)        ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        call b%build_strategy2D_tbox(p)             ! 2D Hadamard matcher built
        if( cline%defined('refs') )then
            call find_ldim_nptcls(p%refs, lfoo, ncls_from_refs)
            ! consistency check
            if( p%ncls /=  ncls_from_refs ) stop 'nrefs /= inputted ncls'
        endif
        startit = 1
        if( cline%defined('startit') )startit = p%startit
        if( startit == 1 )call b%a%clean_updatecnt
        ! execute
        if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
        call cluster2D_exec(b, p, cline, startit) ! partition or not, depending on 'part'
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
        call qsys_job_finished(p, 'simple_commander_cluster2D :: exec_cluster2D')
    end subroutine exec_cluster2D

    subroutine exec_cavgassemble( self, cline )
        use simple_classaverager
        class(cavgassemble_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: fnr, file_stat, j
        p = params(cline, spproj_a_seg=PTCL2D_SEG)        ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        call b%build_strategy2D_tbox(p)
        call cavger_new(b, p, 'class')
        call cavger_assemble_sums_from_parts()
        p%frcs  = 'frcs.bin'
        if( cline%defined('which_iter') )then
            p%refs      = 'cavgs_iter'//int2str_pad(p%which_iter,3)//p%ext
            p%refs_even = 'cavgs_iter'//int2str_pad(p%which_iter,3)//'_even'//p%ext
            p%refs_odd  = 'cavgs_iter'//int2str_pad(p%which_iter,3)//'_odd'//p%ext
            p%frcs  = 'frcs_iter'//int2str_pad(p%which_iter,3)//'.bin'
        else if( .not. cline%defined('refs') )then
            p%refs      = 'start2Drefs'//p%ext
            p%refs_even = 'start2Drefs_even'//p%ext
            p%refs_odd  = 'start2Drefs_odd'//p%ext
        endif
        call cavger_calc_and_write_frcs_and_eoavg(p%frcs)
        ! classdoc gen needs to be after calc of FRCs
        call gen2Dclassdoc( b, p )
        ! write references
        call cavger_write(trim(p%refs),      'merged')
        call cavger_write(trim(p%refs_even), 'even'  )
        call cavger_write(trim(p%refs_odd),  'odd'   )
        call cavger_kill()
        ! write project
        call b%spproj%write(p%projfile)
        ! end gracefully
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        call simple_touch('CAVGASSEMBLE_FINISHED', errmsg='In: commander_rec :: eo_cavgassemble ')
    end subroutine exec_cavgassemble

    subroutine exec_check2D_conv( self, cline )
        class(check2D_conv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        logical      :: converged
        p = params(cline, spproj_a_seg=PTCL2D_SEG)        ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        p%ncls    = b%a%get_n('class')
        converged = b%conv%check_conv2D() ! convergence check
        call cline%set('frac', b%conv%get('frac'))
        if( p%l_doshift )then
            ! activates shift search
            call cline%set('trs', p%trs)
        endif
        if( converged )then
            call cline%set('converged', 'yes')
        else
            call cline%set('converged', 'no')
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK2D_CONV NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check2D_conv

    subroutine exec_rank_cavgs( self, cline )
        use simple_oris, only: oris
        use simple_math, only: hpsort
        class(rank_cavgs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: iclass, pop
        type(oris)   :: clsdoc_ranked
        integer, allocatable :: order(:)
        real,    allocatable :: res(:)
        integer :: ldim(3), ncls
        p = params(cline, spproj_a_seg=CLS2D_SEG)         ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        call find_ldim_nptcls(p%stk, ldim, ncls)
        p%ncls = ncls
        if( b%a%get_noris() == p%ncls )then
            ! all we need to do is fetch from classdoc in projfile &
            ! order according to resolution
            call clsdoc_ranked%new_clean(p%ncls)
            res = b%spproj%os_cls2D%get_all('res')
            allocate(order(p%ncls))
            order = (/(iclass,iclass=1,p%ncls)/)
            call hpsort(res, order)
            do iclass=1,p%ncls
                call clsdoc_ranked%set(iclass, 'class', real(order(iclass)))
                call clsdoc_ranked%set(iclass, 'rank',  real(iclass))
                call clsdoc_ranked%set(iclass, 'pop',   b%a%get(order(iclass),  'pop'))
                call clsdoc_ranked%set(iclass, 'res',   b%a%get(order(iclass),  'res'))
                call clsdoc_ranked%set(iclass, 'corr',  b%a%get(order(iclass), 'corr'))
                call clsdoc_ranked%set(iclass, 'w',     b%a%get(order(iclass),    'w'))
                write(*,'(a,1x,i5,1x,a,1x,i5,1x,a,i5,1x,a,1x,f6.2)') 'CLASS:', order(iclass),&
                &'RANK:', iclass ,'POP:', nint(b%a%get(order(iclass), 'pop')), 'RES:', b%a%get(order(iclass), 'res')
                call b%img%read(p%stk, order(iclass))
                call b%img%write(p%outstk, iclass)
            end do
            call clsdoc_ranked%write('classdoc_ranked.txt')
        else
            ! nothing to do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_RANK_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_rank_cavgs

    subroutine exec_cluster_cavgs( self, cline )
        use simple_cluster_cavgs
        class(cluster_cavgs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        p = params(cline, spproj_a_seg=PTCL2D_SEG)        ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        call cluster_cavgs_exec( b, p )
         ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_cluster_cavgs

end module simple_commander_cluster2D
