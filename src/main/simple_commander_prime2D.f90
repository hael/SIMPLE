! concrete commander: prime2D for simultanous 2D alignment and clustering of single-particle images
module simple_commander_prime2D
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_fileio          ! use all in there
use simple_jiffys          ! use all in there
implicit none

public :: makecavgs_commander
public :: prime2D_commander
public :: cavgassemble_commander
public :: check2D_conv_commander
public :: rank_cavgs_commander
private

type, extends(commander_base) :: makecavgs_commander 
 contains
    procedure :: execute      => exec_makecavgs
end type makecavgs_commander
type, extends(commander_base) :: prime2D_commander 
  contains
    procedure :: execute      => exec_prime2D
end type prime2D_commander
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

contains

    subroutine exec_makecavgs( self, cline )
        use simple_classaverager, only: classaverager
        use simple_binoris_io,    only: binwrite_oritab
        use simple_qsys_funs,     only: qsys_job_finished
        class(makecavgs_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(params)        :: p
        type(build)         :: b
        type(classaverager) :: cavger
        integer :: ncls_in_oritab, icls, fnr, file_stat
        p = params(cline)                                 ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        call b%build_hadamard_prime2D_tbox(p)             ! 2D Hadamard matcher built
        write(*,'(a)') '>>> GENERATING CLUSTER CENTERS'
        ! deal with the orientations
        if( cline%defined('oritab') .and. p%l_remap_classes )then
            call b%a%remap_classes
            ncls_in_oritab = b%a%get_n('class')
            if( cline%defined('ncls') )then
                if( p%ncls < ncls_in_oritab ) stop 'ERROR, inputted ncls < ncls_in_oritab; not allowed!'
                if( p%ncls > ncls_in_oritab )then
                    call b%a%expand_classes(p%ncls)
                endif
            else
                p%ncls = ncls_in_oritab
            endif
        else if( cline%defined('oritab') )then
            if( .not. cline%defined('ncls') ) p%ncls = b%a%get_n('class')
        else if( p%tseries .eq. 'yes' )then
            if( .not. cline%defined('ncls') )then
                stop '# class averages (ncls) need to be part of command line when tseries=yes'
            endif
            call b%a%ini_tseries(p%ncls, 'class')
        else
            if( .not. cline%defined('ncls') )then
                stop 'If no oritab is provided ncls (# class averages) need to be part of command line'
            endif
            call b%a%rnd_cls(p%ncls)
        endif
        if( cline%defined('outfile') )then
            p%oritab = p%outfile
        else
            p%oritab = 'prime2D_startdoc.txt'
        endif
        ! shift multiplication
        if( p%mul > 1. ) call b%a%mul_shifts(p%mul)
        ! setup weights in case the 2D was run without them (specscore will still be there)
        if( p%weights2D.eq.'yes' )then
            if( p%nptcls <= SPECWMINPOP )then
                call b%a%set_all2single('w', 1.0)
            else
                ! frac is one by default in prime2D (no option to set frac)
                ! so spectral weighting is done over all images
                call b%a%calc_spectral_weights(1.0)
            endif
        else
            call b%a%set_all2single('w', 1.0)
        endif
        ! even/odd partitioning
        if( b%a%get_nevenodd() == 0 ) call b%a%partition_eo('class', [p%fromp,p%top])    
        ! write
        if( p%l_distr_exec .and. nint(cline%get_rarg('part')) .eq. 1 )then
            call binwrite_oritab(p%oritab, b%a, [1,p%nptcls])
        else
            call binwrite_oritab(p%oritab, b%a, [1,p%nptcls])
        endif
        ! create class averager
        call cavger%new(b, p, 'class')
        ! transfer ori data to object
        call cavger%transf_oridat(b%a)
        if( cline%defined('filwidth') )then
            ! filament option
            if( p%l_distr_exec)then
                stop 'filwidth mode not implemented for distributed mode; simple_commander_prime2D.f90; exec_makecavgs'
            endif
            call b%img%bin_filament(p%filwidth)
            do icls=1,p%ncls
                call cavger%set_cavg(icls, 'merged', b%img)
            end do
        else
            ! standard cavg assembly
            call cavger%assemble_sums()
        endif
        ! write sums
        if( p%l_distr_exec)then
            call cavger%write_partial_sums()
            call qsys_job_finished( p, 'simple_commander_prime2D :: exec_makecavgs' )
        else
            if( cline%defined('refs') )then
                call cavger%write(p%refs, 'merged')
            else
                call cavger%write('startcavgs'//p%ext, 'merged')
            endif
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_MAKECAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_makecavgs
    
    subroutine exec_prime2D( self, cline )
        use simple_hadamard2D_matcher, only: prime2D_exec
        use simple_qsys_funs,          only: qsys_job_finished
        use simple_imgfile,            only: find_ldim_nptcls
        class(prime2D_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer      :: i, startit, ncls_from_refs, lfoo(3)
        logical      :: converged, l_distr_exec
        converged    = .false.
        l_distr_exec = .false.
        p = params(cline)                                 ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        call b%build_hadamard_prime2D_tbox(p)             ! 2D Hadamard matcher built
        if( cline%defined('refs') )then
            call find_ldim_nptcls(p%refs, lfoo, ncls_from_refs)
            ! consistency check
            if( p%ncls /=  ncls_from_refs ) stop 'nrefs /= inputted ncls'
        endif
        ! execute
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
            call prime2D_exec(b, p, cline, p%startit, converged) ! partition or not, depending on 'part'       
        else
            startit = 1
            if( cline%defined('startit') ) startit = p%startit
            ! extremal dynamics
            if( cline%defined('extr_iter') )then
                ! all is well
            else
                p%extr_iter = startit
            endif
            do i=startit,p%maxits
                call prime2D_exec(b, p, cline, i, converged)
                ! updates extremal dynamics
                p%extr_iter = p%extr_iter + 1
                if(converged) exit
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME2D NORMAL STOP ****')
        ! this is needed for chunk-based prime2D parallellisation
        call qsys_job_finished(p, 'simple_commander_prime2D :: exec_prime2D')
    end subroutine exec_prime2D

    subroutine exec_cavgassemble( self, cline )
        use simple_classaverager, only: classaverager
        class(cavgassemble_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        type(classaverager)  :: cavger
        integer              :: fnr, file_stat
        p = params(cline)                                 ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        call b%build_hadamard_prime2D_tbox(p)
        call cavger%new(b, p, 'class')
        call cavger%assemble_sums_from_parts()
        if( cline%defined('which_iter') )then
            p%refs = 'cavgs_iter'//int2str_pad(p%which_iter,3)//p%ext
        else if( .not. cline%defined('refs') )then
            p%refs = 'startcavgs'//p%ext
        endif
        call cavger%write(trim(p%refs), 'merged')
        ! end gracefully
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        if(.not.fopen(fnr, FILE='CAVGASSEMBLE_FINISHED', STATUS='REPLACE', action='WRITE', iostat=file_stat))&
        call fileio_errmsg('In: commander_rec :: eo_volassemble', file_stat )
        if(.not.fclose( fnr , iostat=file_stat))&
        call fileio_errmsg('In: commander_rec :: eo_volassemble fclose', file_stat )
    end subroutine exec_cavgassemble
    
    subroutine exec_check2D_conv( self, cline )
        class(check2D_conv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        logical      :: converged
        p = params(cline) ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        p%ncls    = b%a%get_n('class')
        converged = b%conv%check_conv2D() ! convergence check
        call cline%set('frac', b%conv%get('frac'))
        if( p%doshift )then
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
        use simple_oris,       only: oris
        use simple_binoris_io, only: binread_oritab, binread_nlines
        class(rank_cavgs_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(params)         :: p
        type(build)          :: b
        integer              :: iclass
        integer, allocatable :: order(:)
        p = params(cline) ! parameters generated
        call b%build_general_tbox(p, cline, do3d=.false.) ! general objects built
        p%ncls   = p%nptcls
        p%nptcls = binread_nlines(p%oritab) 
        call b%a%new(p%nptcls)
        call binread_oritab(p%oritab, b%a, [1,p%nptcls])
        order = b%a%order_cls(p%ncls)
        do iclass=1,p%ncls
            write(*,'(a,1x,i5,1x,a,1x,i5,1x,a,i5)') 'CLASS:', order(iclass),&
            &'CLASS_RANK:', iclass ,'POP:', b%a%get_pop(order(iclass), 'class') 
            call b%img%read(p%stk, order(iclass))
            call b%img%write(p%outstk, iclass)
        end do
        call simple_end('**** SIMPLE_RANK_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_rank_cavgs

end module simple_commander_prime2D
