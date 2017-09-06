! projection-matching based on Hadamard products, high-level search routines for PRIME2D
module simple_hadamard2D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_prime2D_srch,     only: prime2D_srch
use simple_classaverager,    only: classaverager
use simple_ori,              only: ori
use simple_build,            only: build
use simple_params,           only: params
use simple_cmdline,          only: cmdline
use simple_strings,          only: int2str_pad
use simple_jiffys            ! use all in there
use simple_fileio            ! use all in there
use simple_hadamard_common   ! use all in there
use simple_filterer          ! use all in there
use simple_defs              ! use all in there
use simple_syslib            ! use all in there
implicit none

public :: prime2D_exec, preppftcc4align, pftcc
private
#include "simple_local_flags.inc"

type(polarft_corrcalc)          :: pftcc
type(prime2D_srch), allocatable :: primesrch2D(:)
type(classaverager)             :: cavger

contains

    !>  \brief  is the prime2D algorithm
    subroutine prime2D_exec( b, p, cline, which_iter, converged )
        use simple_qsys_funs,   only: qsys_job_finished
        use simple_strings,     only: str_has_substr
        use simple_procimgfile, only: random_selection_from_imgfile
        use simple_binoris_io,  only: binwrite_oritab
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: converged
        logical, allocatable :: ptcl_mask(:)
        real,    allocatable :: frc(:), res(:)
        integer :: iptcl, icls, j
        real    :: corr_thresh, frac_srch_space, skewness, extr_thresh, frc05, frc0143
        logical :: l_do_read

        ! PREP REFERENCES
        call cavger%new(b, p, 'class')
        l_do_read = .true.
        if( p%l_distr_exec )then
            if( .not. cline%defined('refs') )&
            &stop 'need refs to be part of command line for distributed prime2D execution'
        else
            if( which_iter == p%startit )then
                if( .not. cline%defined('refs') .and. cline%defined('oritab') )then
                    ! we make references
                    call cavger%transf_oridat(b%a)
                    call cavger%assemble_sums()
                    l_do_read = .false.
                else
                    ! we randomly select particle images as initial references
                    p%refs = 'start2Drefs'//p%ext
                    if( p%chunktag .ne. '' ) p%refs = trim(p%chunktag)//trim(p%refs)
                    ptcl_mask = b%a%included()
                    call random_selection_from_imgfile(p%stk, p%refs, p%ncls, p%box, p%smpd, ptcl_mask)
                    deallocate(ptcl_mask)
                endif
            endif
        endif
        if( l_do_read )then
            if( .not. file_exists(p%refs) ) stop 'input references (refs) does not exist in cwd'
            call cavger%read(p%refs, 'merged')
        endif

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = b%a%get_avg('frac')

        ! SETUP WEIGHTS
        ! this needs to be done prior to search such that each part
        ! sees the same information in distributed execution
        if( p%weights2D .eq. 'yes' .and. frac_srch_space >= FRAC_INTERPOL )then
            if( p%nptcls <= SPECWMINPOP )then
                call b%a%set_all2single('w', 1.0)
            else
                ! frac is one by default in prime2D (no option to set frac)
                ! so spectral weighting is done over all images
                call b%a%calc_spectral_weights(1.0)
            endif
        else
            ! defaults to unitary weights
            call b%a%set_all2single('w', 1.0)
        endif

        ! READ FOURIER RING CORRELATIONS
        if( which_iter >= LPLIM1ITERBOUND .and. frac_srch_space >= FRAC_INTERPOL )then
            if( file_exists(p%fsc) )then
                call b%projfrcs%read(p%fsc) ! spectral whitening + optiml low-pass filter activated on read
            else
                write(*,*) 'the FRC file does not exist in cwd: ', trim(p%fsc)
                stop 'simple_hadamard2D_matcher :: prime2D_exec'
            endif
        endif

        ! POPULATION BALANCING LOGICS
        ! this needs to be done prior to search such that each part
        ! sees the same information in distributed execution
        if( p%balance > 0 )then
            call b%a%balance(p%balance, skewness)
            write(*,'(A,F8.2)') '>>> CLASS DISTRIBUTION SKEWNESS(%):', 100. * skewness
        else
            call b%a%set_all2single('state_balance', 1.0)
        endif

        ! EXTREMAL LOGICS
        if( frac_srch_space < 98. .and. p%extr_iter <= 15 )then
            extr_thresh = EXTRINITHRESH * (1.-EXTRTHRESH_CONST)**real(p%extr_iter-1)  ! factorial decay
            extr_thresh = min(EXTRINITHRESH, max(0., extr_thresh))
            corr_thresh = b%a%extremal_bound(extr_thresh)
        else
            corr_thresh = -huge(corr_thresh)
        endif

        ! SET FOURIER INDEX RANGE
        call set_bp_range2D( b, p, cline, which_iter, frac_srch_space )

        ! GENERATE REFERENCE & PARTICLE POLAR FTs
        call preppftcc4align( b, p )

        ! INITIALIZE
        write(*,'(A,1X,I3)') '>>> PRIME2D DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter
        if( .not. p%l_distr_exec )then
            p%outfile = 'prime2Ddoc_'//int2str_pad(which_iter,3)//'.txt'
            if( p%chunktag .ne. '' ) p%outfile= trim(p%chunktag)//trim(p%outfile)
        endif

        ! STOCHASTIC IMAGE ALIGNMENT
        allocate( primesrch2D(p%fromp:p%top) )
        do iptcl=p%fromp,p%top
            call primesrch2D(iptcl)%new(p, pftcc)
        end do
        ! calculate CTF matrices
        if( p%ctf .ne. 'no' ) call pftcc%create_polar_ctfmats(b%a)
        ! execute the search
        call del_file(p%outfile)
        select case(trim(p%refine))
            case('neigh')
                !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                do iptcl=p%fromp,p%top
                    call primesrch2D(iptcl)%nn_srch(pftcc, iptcl, b%a, b%nnmat)
                end do
                !$omp end parallel do
            case DEFAULT
                if( p%oritab .eq. '' )then
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        call primesrch2D(iptcl)%exec_prime2D_srch(pftcc, iptcl, b%a, greedy=.true.)
                    end do
                    !$omp end parallel do
                else
                    if( corr_thresh > 0. )then
                        write(*,'(A,F8.2)') '>>> PARTICLE RANDOMIZATION(%):', 100.*extr_thresh
                        write(*,'(A,F8.2)') '>>> CORRELATION THRESHOLD:    ', corr_thresh
                    endif
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        call primesrch2D(iptcl)%exec_prime2D_srch(pftcc, iptcl, b%a, extr_bound=corr_thresh)
                    end do
                    !$omp end parallel do
                endif
        end select
        DebugPrint ' hadamard2D_matcher; completed alignment'

        ! SETUP (OVERWRITES) EVEN/ODD PARTITION
        ! needs to be here since parameters just updated
        ! need to be in the [p%fromp, p%top] range or it will be based on previous params
        call b%a%partition_eo('class', [p%fromp, p%top])

        ! REMAPPING OF HIGHEST POPULATED CLASSES
        ! needs to be here since parameters just updated
        if( p%l_distr_exec )then
            ! this is done in the distributed workflow
        else
            call b%a%fill_empty_classes(p%ncls)
        endif

        ! OUTPUT ORIENTATIONS
        call binwrite_oritab(p%outfile, b%a, [p%fromp,p%top])
        p%oritab = p%outfile

        ! WIENER RESTORATION OF CLASS AVERAGES
        call cavger%transf_oridat(b%a)
        call cavger%set_grid_flag(frac_srch_space >= FRAC_INTERPOL)
        call cavger%assemble_sums()
        if( p%l_distr_exec)then
            call cavger%write_partial_sums()
        else
            p%refs = 'cavgs_iter'//int2str_pad(which_iter,3)//p%ext
            if( p%chunktag .ne. '' ) p%refs = trim(p%chunktag)//trim(p%refs)
            call cavger%write(p%refs, 'merged')
        endif

        ! CALCULATE & WRITE TO DISK FOURIER RING CORRELATION
        if( p%l_distr_exec )then
            ! this is done in prime2D commander; cavgassemble
        else
            p%fsc = 'frcs_iter'//int2str_pad(which_iter,3)//'.bin'
            if( p%chunktag .ne. '' ) p%fsc = trim(p%chunktag)//trim(p%fsc)
            call cavger%calc_and_write_frcs(p%fsc)
            call b%projfrcs%estimate_res(frc, res, frc05, frc0143)
            do j=1,size(res)
                write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', frc(j)
            end do
            write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FRC=0.500 DETERMINED TO:', frc05
            write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FRC=0.143 DETERMINED TO:', frc0143 
            deallocate(frc, res)
        endif

        ! DESTRUCT
        do iptcl=p%fromp,p%top
            call primesrch2D(iptcl)%kill
        end do
        deallocate( primesrch2D )
        call pftcc%kill
        call cavger%kill

        ! REPORT CONVERGENCE
        if( p%l_distr_exec )then
            call qsys_job_finished(p, 'simple_hadamard2D_matcher :: prime2D_exec')
        else
            converged = b%conv%check_conv2D()
        endif
    end subroutine prime2D_exec

    !>  \brief  prepares the polarft corrcalc object for search
    subroutine preppftcc4align( b, p )
        use simple_syslib,       only: alloc_errchk
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        type(ori) :: o
        integer   :: cnt, iptcl, icls, pop, istate
        integer   :: filtsz, alloc_stat, filnum, io_stat
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PRIME2D SEARCH ENGINE'
        ! must be done here since constants in p are dynamically set
        call pftcc%new(p%ncls, [p%fromp,p%top], [p%boxmatch,p%boxmatch,1], p%smpd, p%kfromto, p%ring2, p%ctf)
        ! prepare the polarizers
        call b%img_match%init_polarizer(pftcc)
        ! prepare the automasker
        if( p%l_envmsk .and. p%automsk .eq. 'cavg' ) call b%mskimg%init2D(p, p%ncls)

        ! PREPARATION OF REFERENCES IN PFTCC
        ! read references and transform into polar coordinates
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING REFERENCES'
        do icls=1,p%ncls
            call progress(icls, p%ncls)
            pop = 1
            if( p%oritab /= '' ) pop = b%a%get_pop(icls, 'class')
            if( pop > 0 )then
                ! prepare the reference
                call cavger%get_cavg(icls, 'merged', b%img)
                if( p%oritab /= '' )then
                    call prep2Dref(b, p, icls, center=(pop > MINCLSPOPLIM))
                else
                    call prep2Dref(b, p, icls)
                endif
                ! transfer to polar coordinates
                call b%img_match%polarize(pftcc, icls, isptcl=.false.)
            endif
        end do

        ! PREPARATION OF PARTICLES IN PFTCC
        ! read particle images and create polar projections
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PARTICLES'
        cnt = 0
        do iptcl=p%fromp,p%top
            cnt = cnt+1
            call progress(cnt, p%top-p%fromp+1)
            call read_img_from_stk( b, p, iptcl )
            o = b%a%get_ori(iptcl)
            call prepimg4align(b, p, o)
            ! transfer to polar coordinates
            call b%img_match%polarize(pftcc, iptcl)
        end do
        DebugPrint '*** hadamard2D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align

end module simple_hadamard2D_matcher
