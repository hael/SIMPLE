! projection-matching based on Hadamard products, high-level search routines for CLUSTER2D
module simple_strategy2D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
#include "simple_lib.f08"
use simple_polarft_corrcalc,      only: polarft_corrcalc
use simple_ori,                   only: ori
use simple_build,                 only: build
use simple_params,                only: params
use simple_cmdline,               only: cmdline
use simple_strategy2D,            only: strategy2D
use simple_strategy2D_srch,       only: strategy2D_spec
use simple_strategy2D_greedy,     only: strategy2D_greedy
use simple_strategy2D_neigh,      only: strategy2D_neigh
use simple_strategy2D_stochastic, only: strategy2D_stochastic
use simple_strategy2D_alloc       ! use all in there
use simple_strategy2D3D_common    ! use all in there
use simple_filterer               ! use all in there
use simple_timer                  ! use all in there
use simple_classaverager          ! use all in there

implicit none

public :: cluster2D_exec, preppftcc4align, pftcc
private
#include "simple_local_flags.inc"

logical, parameter             :: L_BENCH         = .false.
logical, parameter             :: L_BENCH_CLUSTER2D = .false.
type(polarft_corrcalc), target :: pftcc
integer                        :: nptcls2update
integer(timer_int_kind)        :: t_init, t_prep_pftcc, t_align, t_cavg, t_tot
real(timer_int_kind)           :: rt_init, rt_prep_pftcc, rt_align, rt_cavg
real(timer_int_kind)           :: rt_tot, rt_refloop, rt_inpl, rt_tot_sum, rt_refloop_sum
real(timer_int_kind)           :: rt_inpl_sum
character(len=STDLEN)          :: benchfname

contains

    !>  \brief  is the prime2D algorithm
    subroutine cluster2D_exec( b, p, cline, which_iter )
        use simple_qsys_funs,    only: qsys_job_finished
        use simple_binoris_io,   only: binwrite_oritab
        class(build),  target, intent(inout) :: b
        class(params), target, intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        integer, allocatable       :: prev_pops(:), pinds(:)
        logical, allocatable       :: ptcl_mask(:)
        class(strategy2D), pointer :: strategy2Dsrch(:)
        type(strategy2D_spec)      :: strategy2Dspec
        integer :: iptcl, icls, i, j, fnr, cnt
        real    :: corr_bound, frac_srch_space, skewness, extr_thresh
        logical :: doprint, l_partial_sums, l_extr, l_frac_update

        if( L_BENCH )then
            t_init = tic()
            t_tot  = t_init
        endif

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = b%a%get_avg('frac')

        ! SWITCHES
        if( p%extr_iter == 1 )then
            ! greedy start
            l_partial_sums = .false.
            l_extr         = .false.
            l_frac_update  = .false.
        else if( p%extr_iter <= MAX_EXTRLIM2D )then
            ! extremal opt without fractional update
            l_partial_sums = .false.
            l_extr         = .true.
            l_frac_update  = .false.
        else
            ! optional fractional update, no extremal opt
            l_partial_sums = p%l_frac_update
            l_extr         = .false.
            l_frac_update  = p%l_frac_update
        endif

        ! EXTREMAL LOGICS
        if( l_extr  )then
            ! factorial decay, -2 because first step is always greedy
            extr_thresh = EXTRINITHRESH * (1.-EXTRTHRESH_CONST)**real(p%extr_iter-2)
            extr_thresh = min(EXTRINITHRESH, max(0., extr_thresh))
            corr_bound  = b%a%extremal_bound(extr_thresh)
            write(*,'(A,F8.2)') '>>> PARTICLE RANDOMIZATION(%):', 100.*extr_thresh
            write(*,'(A,F8.2)') '>>> CORRELATION THRESHOLD:    ', corr_bound
        else
            extr_thresh = 0.
            corr_bound  = -huge(corr_bound)
        endif

        ! PARTICLE INDEX SAMPLING FOR FRACTIONAL UPDATE (OR NOT)
        if( allocated(pinds) )     deallocate(pinds)
        if( allocated(ptcl_mask) ) deallocate(ptcl_mask)
        if( l_frac_update )then
            allocate(ptcl_mask(p%fromp:p%top))
            call b%a%sample4update_and_incrcnt2D(p%ncls, [p%fromp,p%top], p%update_frac, nptcls2update, pinds, ptcl_mask)
            ! correct convergence stats
            do iptcl=p%fromp,p%top
                if( .not. ptcl_mask(iptcl) )then
                    ! these are not updated
                    call b%a%set(iptcl, 'mi_class',    1.0)
                    call b%a%set(iptcl, 'mi_inpl',     1.0)
                    call b%a%set(iptcl, 'mi_joint',    1.0)
                    call b%a%set(iptcl, 'dist_inpl',   0.0)
                    call b%a%set(iptcl, 'frac',      100.0)
                endif
            end do
        else
            nptcls2update = p%top - p%fromp + 1
            allocate(pinds(nptcls2update), ptcl_mask(p%fromp:p%top))
            pinds = (/(i,i=p%fromp,p%top)/)
            ptcl_mask = .true.
        endif

        ! PREP REFERENCES
        call cavger_new(b, p, 'class', ptcl_mask)
        if( b%a%get_nevenodd() == 0 )then
            stop 'ERROR! no eo partitioning available; strategy2D_matcher :: cluster2D_exec'
        endif
        if( .not. cline%defined('refs') ) stop 'need refs to be part of command line for cluster2D execution'
        if( .not. file_exists(p%refs) )   stop 'input references (refs) does not exist in cwd'
        call cavger_read(p%refs, 'merged')
        if( file_exists(p%refs_even) )then
            call cavger_read(p%refs_even, 'even')
        else
            call cavger_read(p%refs, 'even')
        endif
        if( file_exists(p%refs_odd) )then
            call cavger_read(p%refs_odd, 'odd')
        else
            call cavger_read(p%refs, 'odd')
        endif

        ! SETUP WEIGHTS
        ! this needs to be done prior to search such that each part
        ! sees the same information in distributed execution
        if( p%weights2D .eq. 'yes' .and. frac_srch_space >= FRAC_INTERPOL )then
            if( p%nptcls <= SPECWMINPOP )then
                call b%a%set_all2single('w', 1.0)
                ! call b%a%calc_hard_weights(p%frac)
            else
                if( p%weights2D .eq. 'yes' .and. which_iter > 3 )then
                    call b%a%get_pops(prev_pops, 'class', consider_w=.true., maxn=p%ncls)
                else
                    call b%a%get_pops(prev_pops, 'class', consider_w=.false., maxn=p%ncls)
                endif
                ! frac is one by default in cluster2D (no option to set frac)
                ! so spectral weighting is done over all images
                call b%a%calc_spectral_weights(1.0)
                ! call b%a%calc_spectral_weights(p%frac)
                if( any(prev_pops == 0) )then
                    ! now ensuring the spectral re-ranking does not re-populates
                    ! zero-populated classes, for congruence with empty cavgs
                    do icls = 1, p%ncls
                        if( prev_pops(icls) > 0 ) cycle
                        call b%a%get_pinds(icls, 'class', pinds, consider_w=.false.)
                        if( .not.allocated(pinds) )cycle
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
            ! call b%a%set_all2single('w', p%frac) ! should be done by class
        endif

        ! READ FOURIER RING CORRELATIONS
        if( file_exists(p%frcs) ) call b%projfrcs%read(p%frcs)

        ! SET FOURIER INDEX RANGE
        call set_bp_range2D( b, p, cline, which_iter, frac_srch_space )
        if( L_BENCH ) rt_init = toc(t_init)

        ! GENERATE REFERENCE & PARTICLE POLAR FTs
        if( L_BENCH ) t_prep_pftcc = tic()
        call preppftcc4align( b, p, which_iter )
        if( L_BENCH ) rt_prep_pftcc = toc(t_prep_pftcc)

        ! INITIALIZE
        write(*,'(A,1X,I3)') '>>> CLUSTER2D DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter

        ! STOCHASTIC IMAGE ALIGNMENT
        ! array allocation for strategy2D
        call prep_strategy2D( b, p, ptcl_mask, which_iter )
        ! switch for polymorphic strategy2D construction
        select case(trim(p%neigh))
            case('yes')
                allocate(strategy2D_neigh :: strategy2Dsrch(p%fromp:p%top))
            case DEFAULT
                if( b%spproj%is_virgin_field('ptcl2D') )then
                    allocate(strategy2D_greedy :: strategy2Dsrch(p%fromp:p%top))
                else
                    allocate(strategy2D_stochastic :: strategy2Dsrch(p%fromp:p%top))
                endif
        end select
        ! actual construction
        cnt = 0
        do iptcl = p%fromp, p%top
            if( ptcl_mask(iptcl) )then
                cnt = cnt + 1
                ! search spec
                strategy2Dspec%iptcl      =  iptcl
                strategy2Dspec%iptcl_map  =  cnt
                strategy2Dspec%corr_bound =  corr_bound
                strategy2Dspec%pp         => p
                strategy2Dspec%ppftcc     => pftcc
                strategy2Dspec%pa         => b%a
                if( allocated(b%nnmat) ) strategy2Dspec%nnmat => b%nnmat
                ! search object
                call strategy2Dsrch(iptcl)%new(strategy2Dspec)
            endif
        end do
        ! memoize CTF matrices
        if( b%spproj%get_ctfflag('ptcl2D').ne.'no' ) call pftcc%create_polar_ctfmats(b%spproj, 'ptcl2D')
        ! memoize FFTs for improved performance
        call pftcc%memoize_ffts
        
        ! SEARCH
        if( L_BENCH ) t_align = tic()
        if( L_BENCH_CLUSTER2D )then
            rt_refloop_sum = 0.
            rt_inpl_sum    = 0.
            rt_tot_sum     = 0.
        endif
        !$omp parallel do default(shared) schedule(guided) private(i,iptcl) proc_bind(close)
        do i=1,nptcls2update
           iptcl = pinds(i)
           call strategy2Dsrch(iptcl)%srch
        end do
        !$omp end parallel do

        ! CLEAN-UP
        call clean_strategy2D()
        call pftcc%kill
        do i=1,nptcls2update
           iptcl = pinds(i)
           call strategy2Dsrch(iptcl)%kill
        end do
        deallocate( strategy2Dsrch )
        if( L_BENCH ) rt_align = toc(t_align)
        DebugPrint ' strategy2D_matcher; completed alignment'

        ! OUTPUT ORIENTATIONS
        call binwrite_oritab(p%outfile, b%spproj, b%a, [p%fromp,p%top], isegment=PTCL2D_SEG)
        p%oritab = p%outfile

        ! WIENER RESTORATION OF CLASS AVERAGES
        if( L_BENCH ) t_cavg = tic()
        call cavger_transf_oridat( b%a )
        call cavger_assemble_sums( l_partial_sums )
        ! write results to disk
        call cavger_readwrite_partial_sums('write')
        call cavger_kill
        if( L_BENCH ) rt_cavg = toc(t_cavg)

        call qsys_job_finished(p, 'simple_strategy2D_matcher :: cluster2D_exec')
        if( L_BENCH )then
            rt_tot  = toc(t_tot)
            doprint = .true.
            if( p%part /= 1 ) doprint = .false.
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
        if( L_BENCH_CLUSTER2D )then
            doprint = .true.
            if( p%part /= 1 ) doprint = .false.
            if( doprint )then
                benchfname = 'CLUSTER2D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
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
    end subroutine cluster2D_exec

    !>  \brief  prepares the polarft corrcalc object for search
    subroutine preppftcc4align( b, p, which_iter )
        use simple_polarizer, only: polarizer
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        integer,       intent(in)    :: which_iter
        type(polarizer), allocatable :: match_imgs(:)
        integer   :: iptcl, icls, pop, pop_even, pop_odd
        integer   :: batchlims(2), imatch, batchsz_max, iptcl_batch
        logical   :: do_center
        real      :: xyz(3)
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
                do_center = (p%oritab /= '' .and. (pop > MINCLSPOPLIM) .and.&
                    &(which_iter > 2) .and. .not.p%l_frac_update)
                ! here we are determining the shifts and map them back to classes
                call prep2Dref(b, p, cavgs_merged(icls), match_imgs(icls), icls, center=do_center, xyz_out=xyz)
                if( pop_even >= MINCLSPOPLIM .and. pop_odd >= MINCLSPOPLIM )then
                    ! here we are passing in the shifts and do NOT map them back to classes
                    call prep2Dref(b, p, cavgs_even(icls), match_imgs(icls), icls, center=do_center, xyz_in=xyz)
                    call match_imgs(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.true.)  ! 2 polar coords
                    ! here we are passing in the shifts and do NOT map them back to classes
                    call prep2Dref(b, p, cavgs_odd(icls), match_imgs(icls), icls, center=do_center, xyz_in=xyz)
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
        call build_pftcc_particles( b, p, pftcc, batchsz_max, match_imgs, .false.)
        ! DESTRUCT
        do imatch=1,batchsz_max
            call match_imgs(imatch)%kill_polarizer
            call match_imgs(imatch)%kill
            call b%imgbatch(imatch)%kill
        end do
        deallocate(match_imgs, b%imgbatch)
        DebugPrint '*** strategy2D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align

end module simple_strategy2D_matcher
