! projection-matching based on Hadamard products, high-level search routines for PRIME2D
module simple_hadamard2D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
#include "simple_lib.f08"
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_ori,              only: ori
use simple_build,            only: build
use simple_params,           only: params
use simple_cmdline,          only: cmdline
use simple_hadamard_common   ! use all in there
use simple_filterer          ! use all in there
use simple_timer             ! use all in there
use simple_classaverager,    ! use all in there
use simple_prime2D_srch,     ! use all in there
implicit none

public :: prime2D_exec, preppftcc4align, pftcc
private
#include "simple_local_flags.inc"

logical, parameter              :: L_BENCH         = .true.
logical, parameter              :: L_BENCH_PRIME2D = .true.
type(polarft_corrcalc)          :: pftcc
type(prime2D_srch), allocatable :: primesrch2D(:)
integer(timer_int_kind)         :: t_init, t_prep_pftcc, t_align, t_cavg, t_tot
real(timer_int_kind)            :: rt_init, rt_prep_pftcc, rt_align, rt_cavg
real(timer_int_kind)            :: rt_tot, rt_refloop, rt_inpl, rt_tot_sum, rt_refloop_sum
real(timer_int_kind)            :: rt_inpl_sum
character(len=STDLEN)           :: benchfname

contains

    !>  \brief  is the prime2D algorithm
    subroutine prime2D_exec( b, p, cline, which_iter, converged )
        use simple_qsys_funs,   only: qsys_job_finished
        use simple_procimgfile, only: random_selection_from_imgfile, copy_imgfile
        use simple_binoris_io,  only: binwrite_oritab
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: converged
        integer, allocatable :: prev_pops(:), pinds(:)
        logical, allocatable :: ptcl_mask(:)
        integer :: iptcl, icls, j, fnr
        real    :: corr_thresh, frac_srch_space, skewness, extr_thresh
        logical :: l_do_read, doprint

        ! PREP REFERENCES
        if( L_BENCH )then
            t_init = tic()
            t_tot  = t_init
        endif
        call cavger_new(b, p, 'class')
        l_do_read = .true.
        if( p%l_distr_exec )then
            if( b%a%get_nevenodd() == 0 )then
                stop 'ERROR! no eo partitioning available; hadamard3D_matcher :: prime2D_exec'
            endif
            if( .not. cline%defined('refs') )&
            &stop 'need refs to be part of command line for distributed prime2D execution'
        else
            if( b%a%get_nevenodd() == 0 )then
                call b%a%partition_eo
            endif
            if( which_iter == p%startit )then
                if( .not. cline%defined('refs') .and. cline%defined('oritab') )then
                    ! we make references
                    call cavger_transf_oridat(b%a)
                    call cavger_assemble_sums()
                    l_do_read = .false.
                else if( cline%defined('refs') )then
                    ! do nothing
                else
                    ! we randomly select particle images as initial references
                    p%refs      = 'start2Drefs'//p%ext
                    p%refs_even = 'start2Drefs_even'//p%ext
                    p%refs_odd  = 'start2Drefs_odd'//p%ext
                    ptcl_mask   = b%a%included()
                    if( cline%defined('stktab') )then
                        call random_selection_from_imgfile(p%stktab, p%refs, p%ncls, p%box, p%smpd, ptcl_mask)
                    else
                        call random_selection_from_imgfile(p%stk, p%refs, p%ncls, p%box, p%smpd, ptcl_mask)
                    endif
                    deallocate(ptcl_mask)
                    call copy_imgfile(trim(p%refs), trim(p%refs_even), p%smpd, [1,p%ncls])
                    call copy_imgfile(trim(p%refs), trim(p%refs_odd),  p%smpd, [1,p%ncls])
                endif
            endif
        endif
        if( l_do_read )then
            if( .not. file_exists(p%refs) ) stop 'input references (refs) does not exist in cwd'
            call cavger_read(p%refs,      'merged')
            call cavger_read(p%refs_even, 'even'  )
            call cavger_read(p%refs_odd,  'odd'   )
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
                if( p%weights2D .eq. 'yes' .and. which_iter > 3 )then
                    call b%a%get_pops(prev_pops, 'class', consider_w=.true., maxn=p%ncls)
                else
                    call b%a%get_pops(prev_pops, 'class', consider_w=.false., maxn=p%ncls)
                endif
                ! frac is one by default in prime2D (no option to set frac)
                ! so spectral weighting is done over all images
                call b%a%calc_spectral_weights(1.0)
                if( any(prev_pops == 0) )then
                    ! now ensuring the spectral re-ranking does not re-populates 
                    ! zero-populated classes, for congruence with empty cavgs
                    do icls = 1, p%ncls
                        if( prev_pops(icls) > 0 )cycle
                        if( b%a%get_pop(icls, 'class', consider_w=.false.) == 0 )cycle
                        pinds = b%a%get_pinds(icls, 'class', consider_w=.false.)
                        do iptcl = 1, size(pinds)
                            call b%a%set(pinds(iptcl), 'w', 0.)
                        enddo
                        deallocate(pinds)
                    enddo
                endif
                deallocate(prev_pops)
            endif
        else
            ! defaults to unitary weights
            call b%a%set_all2single('w', 1.0)
        endif

        ! READ FOURIER RING CORRELATIONS
        if( file_exists(p%frcs) ) call b%projfrcs%read(p%frcs)

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
        if( L_BENCH ) rt_init = toc(t_init)

        ! GENERATE REFERENCE & PARTICLE POLAR FTs
        if( L_BENCH ) t_prep_pftcc = tic()
        call preppftcc4align( b, p, which_iter )
        if( L_BENCH ) rt_prep_pftcc = toc(t_prep_pftcc)

        ! INITIALIZE
        write(*,'(A,1X,I3)') '>>> PRIME2D DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter
        if( .not. p%l_distr_exec )then
            p%outfile = 'prime2Ddoc_'//int2str_pad(which_iter,3)//METADATEXT
        endif

        ! STOCHASTIC IMAGE ALIGNMENT
        ! create the search objects, need to re-create every round because parameters are changing
        call prep4prime2D_srch( b, p, which_iter )
        allocate( primesrch2D(p%fromp:p%top), stat=alloc_stat)
        allocchk("In hadamard2D_matcher::prime2D_exec primesrch2D objects ")
        do iptcl=p%fromp,p%top
            call primesrch2D(iptcl)%new(iptcl, pftcc, b%a, p)
        end do
        ! generate CTF matrices
        if( p%ctf .ne. 'no' ) call pftcc%create_polar_ctfmats(b%a)
        ! memoize FFTs for improved performance
        call pftcc%memoize_ffts
        ! execute the search
        call del_file(p%outfile)
        if( L_BENCH ) t_align = tic()
        select case(trim(p%refine))
            case('neigh')
                !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                do iptcl=p%fromp,p%top
                    call primesrch2D(iptcl)%nn_srch(b%nnmat)
                end do
                !$omp end parallel do
            case DEFAULT
                if( p%oritab .eq. '' )then
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        call primesrch2D(iptcl)%exec_prime2D_srch(greedy=.true.)
                    end do
                    !$omp end parallel do
                else
                    if( corr_thresh > 0. )then
                        write(*,'(A,F8.2)') '>>> PARTICLE RANDOMIZATION(%):', 100.*extr_thresh
                        write(*,'(A,F8.2)') '>>> CORRELATION THRESHOLD:    ', corr_thresh
                    endif
                    if( L_BENCH_PRIME2D )then
                        rt_refloop_sum = 0.
                        rt_tot_sum     = 0.
                    endif
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        call primesrch2D(iptcl)%exec_prime2D_srch(extr_bound=corr_thresh)
                        if( L_BENCH_PRIME2D )then
                            call primesrch2D(iptcl)%get_times(rt_refloop, rt_inpl, rt_tot)
                            rt_refloop_sum      = rt_refloop_sum + rt_refloop
                            rt_inpl_sum         = rt_inpl_sum    + rt_inpl
                            rt_tot_sum          = rt_tot_sum     + rt_tot
                        endif
                    end do
                    !$omp end parallel do
                endif
        end select
        ! pftcc & primesrch2D not needed anymore
        call pftcc%kill
        deallocate( primesrch2D )
        if( L_BENCH ) rt_align = toc(t_align)
        DebugPrint ' hadamard2D_matcher; completed alignment'

        ! REMAPPING OF HIGHEST POPULATED CLASSES
        ! needs to be here since parameters just updated
        if( p%dyncls.eq.'yes' )then
            if( p%l_distr_exec )then
                ! this is done in the distributed workflow
            else
                call b%a%fill_empty_classes(p%ncls)
            endif
        endif

        ! OUTPUT ORIENTATIONS
        call binwrite_oritab(p%outfile, b%a, [p%fromp,p%top])
        p%oritab = p%outfile

        ! WIENER RESTORATION OF CLASS AVERAGES
        if( L_BENCH ) t_cavg = tic()
        call cavger_transf_oridat(b%a)
        call cavger_assemble_sums()
        ! write results to disk
        if( p%l_distr_exec )then
            call cavger_write_partial_sums()
        else
            p%frcs = 'frcs_iter'//int2str_pad(which_iter,3)//'.bin'
            call cavger_calc_and_write_frcs_and_eoavg(p%frcs)
            call gen2Dclassdoc( b, p, 'classdoc_'//int2str_pad(which_iter,3)//'.txt')
            p%refs      = 'cavgs_iter'//int2str_pad(which_iter,3)//p%ext
            p%refs_even = 'cavgs_iter'//int2str_pad(which_iter,3)//'_even'//p%ext
            p%refs_odd  = 'cavgs_iter'//int2str_pad(which_iter,3)//'_odd'//p%ext
            call cavger_write(p%refs,      'merged')
            call cavger_write(p%refs_even, 'even'  )
            call cavger_write(p%refs_odd,  'odd'   )
        endif
        call cavger_kill
        if( L_BENCH ) rt_cavg = toc(t_cavg)

        ! REPORT CONVERGENCE
        if( p%l_distr_exec )then
            call qsys_job_finished(p, 'simple_hadamard2D_matcher :: prime2D_exec')
        else
            converged = b%conv%check_conv2D()
        endif
        if( L_BENCH )then
            rt_tot  = toc(t_tot)
            doprint = .true.
            if( p%l_distr_exec .and. p%part /= 1 ) doprint = .false.
            if( doprint )then
                benchfname = 'HADAMARD2D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
                call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation       : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'pftcc preparation    : ', rt_prep_pftcc
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', rt_align
                write(fnr,'(a,1x,f9.2)') 'class averaging      : ', rt_cavg
                write(fnr,'(a,1x,f9.2)') 'total time           : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation       : ', (rt_init/rt_tot)        * 100.
                write(fnr,'(a,1x,f9.2)') 'pftcc preparation    : ', (rt_prep_pftcc/rt_tot)  * 100.
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', (rt_align/rt_tot)       * 100.
                write(fnr,'(a,1x,f9.2)') 'class averaging      : ', (rt_cavg/rt_tot)        * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for      : ',&
                    &((rt_init+rt_prep_pftcc+rt_align+rt_cavg)/rt_tot) * 100.
                call fclose(fnr)
            endif
        endif
        if( L_BENCH_PRIME2D )then
            doprint = .true.
            if( p%l_distr_exec .and. p%part /= 1 ) doprint = .false.
            if( doprint )then
                benchfname = 'PRIME2D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
                call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'refloop : ', rt_refloop_sum
                write(fnr,'(a,1x,f9.2)') 'in-plane: ', rt_inpl_sum
                write(fnr,'(a,1x,f9.2)') 'tot     : ', rt_tot_sum
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'refloop : ', (rt_refloop_sum/rt_tot_sum)      * 100.
                write(fnr,'(a,1x,f9.2)') 'in-plane: ', (rt_inpl_sum/rt_tot_sum)         * 100.
            endif
        endif
    end subroutine prime2D_exec

    !>  \brief  prepares the polarft corrcalc object for search
    subroutine preppftcc4align( b, p, which_iter )
        use simple_polarizer, only: polarizer
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        integer,       intent(in)    :: which_iter
        type(polarizer), allocatable :: match_imgs(:)
        type(ori) :: o
        integer   :: iptcl, icls, pop, pop_even, pop_odd
        integer   :: batchlims(2), batchsz, imatch, batchsz_max, iptcl_batch
        logical   :: do_center
        real      :: xyz(3)
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PRIME2D SEARCH ENGINE'

        ! create the polarft_corrcalc object
        call pftcc%new(p%ncls, p, eoarr=nint(b%a%get_all('eo', [p%fromp,p%top])))
        ! prepare the polarizer images
        call b%img_match%init_polarizer(pftcc, p%alpha)
        ! this is tro avoid excessive allocation, allocate what is the upper bound on the 
        ! # matchimgs needed for both parallel loops
        batchsz_max = max(MAXIMGBATCHSZ,p%ncls)
        allocate(match_imgs(batchsz_max))
        do imatch=1,batchsz_max
            call match_imgs(imatch)%new([p%boxmatch, p%boxmatch, 1], p%smpd)
            call match_imgs(imatch)%copy_polarizer(b%img_match)
        end do

        ! PREPARATION OF REFERENCES IN PFTCC
        ! read references and transform into polar coordinates
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING REFERENCES'
        !$omp parallel do default(shared) private(icls,pop,pop_even,pop_odd,do_center,xyz)&
        !$omp schedule(static) proc_bind(close)
        do icls=1,p%ncls
            pop      = 1
            pop_even = 0
            pop_odd  = 0
            if( p%oritab /= '' )then
                pop      = b%a%get_pop(icls, 'class'      )
                pop_even = b%a%get_pop(icls, 'class', eo=0)
                pop_odd  = b%a%get_pop(icls, 'class', eo=1)
            endif
            if( pop > 0 )then
                ! prepare the references
                do_center = (p%oritab /= '' .and. (pop > MINCLSPOPLIM) .and. (which_iter > 2))
                ! here we are determining the shifts and map them back to classes
                call prep2Dref(b, p, cavgs_merged(1,icls), match_imgs(icls), icls, center=do_center, xyz_out=xyz)
                if( pop_even >= MINCLSPOPLIM .and. pop_odd >= MINCLSPOPLIM )then
                    ! here we are passing in the shifts and do NOT map them back to classes
                    call prep2Dref(b, p, cavgs_even(1,icls), match_imgs(icls), icls, center=do_center, xyz_in=xyz) 
                    call match_imgs(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.true.)  ! 2 polar coords
                    ! here we are passing in the shifts and do NOT map them back to classes
                    call prep2Dref(b, p, cavgs_odd(1,icls), match_imgs(icls), icls, center=do_center, xyz_in=xyz)
                    call match_imgs(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.false.) ! 2 polar coords
                else
                    ! put the merged class average in both even and odd positions
                    call match_imgs(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.true. ) ! 2 polar coords
                    call pftcc%cp_even2odd_ref(icls)
                endif   
            endif
        end do
        !$omp end parallel do
        
        ! PREPARATION OF PARTICLES IN PFTCC
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PARTICLES'
        call prepimgbatch(b, p, batchsz_max)
        do iptcl_batch=p%fromp,p%top,batchsz_max
            batchlims = [iptcl_batch,min(p%top,iptcl_batch + batchsz_max - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call read_imgbatch( b, p, batchlims)
            !$omp parallel do default(shared) private(iptcl,imatch)&
            !$omp schedule(static) proc_bind(close)
            do iptcl=batchlims(1),batchlims(2)
                imatch = iptcl - batchlims(1) + 1
                call prepimg4align(b, p, iptcl, b%imgbatch(imatch), match_imgs(imatch), is3D=.false.)
                ! transfer to polar coordinates
                call match_imgs(imatch)%polarize(pftcc, iptcl, .true., .true.)
            end do
            !$omp end parallel do 
        end do
        ! DESTRUCT
        do imatch=1,batchsz_max
            call match_imgs(imatch)%kill_polarizer
            call match_imgs(imatch)%kill
            call b%imgbatch(imatch)%kill
        end do
        deallocate(match_imgs, b%imgbatch)

        DebugPrint '*** hadamard2D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align

end module simple_hadamard2D_matcher
