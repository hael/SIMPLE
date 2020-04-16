! projection-matching based on Hadamard products, high-level search routines for CLUSTER2D
module simple_strategy2D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_cmdline,          only: cmdline
use simple_builder,          only: build_glob
use simple_parameters,       only: params_glob
use simple_classaverager
implicit none

public :: cluster2D_exec
private
#include "simple_local_flags.inc"

logical, parameter      :: L_BENCH = .false.

type(polarft_corrcalc)  :: pftcc
logical,    allocatable :: ptcl_mask(:)
integer                 :: nptcls2update
integer(timer_int_kind) :: t_init, t_prep_pftcc, t_align, t_cavg, t_tot
real(timer_int_kind)    :: rt_init, rt_prep_pftcc, rt_align, rt_cavg
real(timer_int_kind)    :: rt_tot, rt_tot_sum, rt_refloop_sum
real(timer_int_kind)    :: rt_inpl_sum
character(len=STDLEN)   :: benchfname

contains

    !>  \brief  is the prime2D algorithm
    subroutine cluster2D_exec( cline, which_iter )
        use simple_qsys_funs,    only: qsys_job_finished
        use simple_binoris_io,   only: binwrite_oritab
        use simple_strategy2D3D_common,   only: set_bp_range2d
        use simple_strategy2D,            only: strategy2D, strategy2D_per_ptcl
        use simple_strategy2D_srch,       only: strategy2D_spec
        use simple_strategy2D_alloc,      only: prep_strategy2d,clean_strategy2d
        use simple_strategy2D_greedy,     only: strategy2D_greedy
        use simple_strategy2D_tseries,    only: strategy2D_tseries
        use simple_strategy2D_neigh,      only: strategy2D_neigh
        use simple_strategy2D_snhc,       only: strategy2D_snhc
        use simple_strategy2D_inpl,       only: strategy2D_inpl
        use simple_strategy2D_eval,       only: strategy2D_eval
        class(cmdline),          intent(inout) :: cline
        integer,                 intent(in)    :: which_iter
        type(strategy2D_per_ptcl), allocatable :: strategy2Dsrch(:)
        type(strategy2D_spec) :: strategy2Dspec
        integer, allocatable  :: pinds(:)
        real                  :: snhc_sz(params_glob%fromp:params_glob%top)
        integer               :: chunk_id(params_glob%fromp:params_glob%top)
        real                  :: frac_srch_space, rrnd
        integer               :: iptcl, i, fnr, cnt, updatecnt
        logical               :: doprint, l_partial_sums, l_frac_update
        logical               :: l_snhc, l_greedy, l_stream, l_np_cls_defined
        if( L_BENCH )then
            t_init = tic()
            t_tot  = t_init
        endif

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = build_glob%spproj_field%get_avg('frac')

        ! SWITCHES
        l_partial_sums = .false.
        l_snhc         = .false.
        l_greedy       = .false.
        l_frac_update  = .false.
        l_stream       = trim(params_glob%stream).eq.'yes'
        if( params_glob%extr_iter == 1 )then
            ! greedy start
            l_partial_sums = .false.
            l_snhc         = .false.
            l_greedy       = .true.
            l_frac_update  = .false.
        else if( params_glob%extr_iter <= MAX_EXTRLIM2D )then
            ! extremal opt without fractional update
            l_partial_sums = .false.
            l_frac_update  = .false.
            l_greedy       = .false.
            l_snhc         = .true.
            if(params_glob%refine.eq.'greedy')then
                l_greedy = .true.
                l_snhc   = .false.
            endif
        else
            ! optional fractional update, no snhc opt
            l_partial_sums = params_glob%l_frac_update
            l_frac_update  = params_glob%l_frac_update
            l_snhc         = .false.
            l_greedy       = (params_glob%refine.eq.'greedy') .or.(params_glob%cc_objfun.eq.OBJFUN_EUCLID)
        endif

        ! PARTICLE INDEX SAMPLING FOR FRACTIONAL UPDATE (OR NOT)
        if( allocated(pinds) )     deallocate(pinds)
        if( allocated(ptcl_mask) ) deallocate(ptcl_mask)
        if( l_frac_update )then
            allocate(ptcl_mask(params_glob%fromp:params_glob%top))
            call build_glob%spproj_field%sample4update_and_incrcnt2D(params_glob%ncls, &
                [params_glob%fromp,params_glob%top], params_glob%update_frac, nptcls2update, pinds, ptcl_mask)
        else
            call build_glob%spproj_field%mask_from_state(1, ptcl_mask, pinds, fromto=[params_glob%fromp,params_glob%top])
            call build_glob%spproj_field%incr_updatecnt([params_glob%fromp,params_glob%top], mask=ptcl_mask)
            nptcls2update = count(ptcl_mask)
        endif

        ! SNHC LOGICS
        if( l_snhc )then
            ! factorial decay, -2 because first step is always greedy
            snhc_sz = min(SNHC2D_INITFRAC,&
                &max(0.,SNHC2D_INITFRAC*(1.-SNHC2D_DECAY)**real(params_glob%extr_iter-2)))
            write(logfhandle,'(A,F8.2)') '>>> STOCHASTIC NEIGHBOURHOOD SIZE(%):',&
                &100.*(1.-snhc_sz(params_glob%fromp))
        else
            snhc_sz = 0. ! full neighbourhood
        endif

        ! CHUNKS
        chunk_id = 1

        ! ARRAY ALLOCATION FOR STRATEGY2D prior to weights
        call prep_strategy2D( ptcl_mask, which_iter )
        write(logfhandle,'(A)') '>>> STRATEGY2D OBJECTS ALLOCATED'

        ! SETUP WEIGHTS
        ! this needs to be done prior to search such that each part
        ! sees the same information in distributed execution
        ! has to be done prior to classaverager initialization
        if( which_iter > 3 )then
            if( trim(params_glob%ptclw) .eq. 'yes' )then
                call build_glob%spproj_field%calc_soft_weights2D
            else
                call build_glob%spproj_field%calc_hard_weights2D(params_glob%frac, params_glob%ncls)
            endif
        else
            call build_glob%spproj_field%set_all2single('w', 1.0)
        endif
        write(logfhandle,'(A)') '>>> WEIGHTS ASSIGNED'

        ! PREP REFERENCES
        call cavger_new( 'class', ptcl_mask)
        if( build_glob%spproj_field%get_nevenodd() == 0 )then
            THROW_HARD('no eo partitioning available; cluster2D_exec')
        endif
        if( .not. cline%defined('refs') )         THROW_HARD('need refs to be part of command line for cluster2D execution')
        if( .not. file_exists(params_glob%refs) ) THROW_HARD('input references (refs) does not exist in cwd')
        call cavger_read(params_glob%refs, 'merged')
        if( file_exists(params_glob%refs_even) )then
            call cavger_read(params_glob%refs_even, 'even')
        else
            call cavger_read(params_glob%refs, 'even')
        endif
        if( file_exists(params_glob%refs_odd) )then
            call cavger_read(params_glob%refs_odd, 'odd')
        else
            call cavger_read(params_glob%refs, 'odd')
        endif

        ! READ FOURIER RING CORRELATIONS
        if( file_exists(params_glob%frcs) ) call build_glob%projfrcs%read(params_glob%frcs)
        if( params_glob%l_pssnr )then
            if( file_exists(trim(PSSNR_FBODY)//int2str_pad(1,2)//BIN_EXT) )then
                call build_glob%projpssnrs%read(trim(PSSNR_FBODY)//int2str_pad(1,2)//BIN_EXT)
            endif
        endif
        ! SET FOURIER INDEX RANGE
        call set_bp_range2D(cline, which_iter, frac_srch_space )
        if( L_BENCH ) rt_init = toc(t_init)
        ! GENERATE REFERENCE & PARTICLE POLAR FTs
        if( L_BENCH ) t_prep_pftcc = tic()
        call preppftcc4align( which_iter )
        if( L_BENCH ) rt_prep_pftcc = toc(t_prep_pftcc)

        ! INITIALIZE STOCHASTIC IMAGE ALIGNMENT
        write(logfhandle,'(A,1X,I3)') '>>> CLUSTER2D DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter
        ! switch for polymorphic strategy2D construction
        l_np_cls_defined = cline%defined('nptcls_per_cls')
        allocate(strategy2Dsrch(params_glob%fromp:params_glob%top))
        select case(trim(params_glob%neigh))
        case('yes')
            do iptcl=params_glob%fromp,params_glob%top
                if( ptcl_mask(iptcl) ) allocate(strategy2D_neigh :: strategy2Dsrch(iptcl)%ptr)
            end do
        case DEFAULT
            do iptcl=params_glob%fromp,params_glob%top
                if( ptcl_mask(iptcl) )then
                    updatecnt = nint(build_glob%spproj_field%get(iptcl,'updatecnt'))
                    if( l_stream )then
                        ! online mode, based on particle history
                        if( updatecnt <= STREAM_SRCHLIM )then
                            ! reproduces offline mode
                            if( l_greedy .or. (.not.build_glob%spproj_field%has_been_searched(iptcl) .or. updatecnt==1) )then
                                allocate(strategy2D_greedy :: strategy2Dsrch(iptcl)%ptr, stat=alloc_stat)
                            else
                                allocate(strategy2D_snhc   :: strategy2Dsrch(iptcl)%ptr, stat=alloc_stat)
                            endif
                        else
                            ! stochastic per particle decision: stochastic/greedy or in-plane searches or evaluation
                            rrnd = ran3()
                            if( rrnd < STREAM_SRCHFRAC )then
                                if( l_greedy)then
                                    allocate(strategy2D_greedy :: strategy2Dsrch(iptcl)%ptr, stat=alloc_stat)
                                else
                                    allocate(strategy2D_snhc   :: strategy2Dsrch(iptcl)%ptr, stat=alloc_stat)
                                endif
                            else if( rrnd < STREAM_SRCHFRAC+STREAM_INPLFRAC )then
                                allocate(strategy2D_inpl :: strategy2Dsrch(iptcl)%ptr, stat=alloc_stat)
                            else
                                allocate(strategy2D_eval :: strategy2Dsrch(iptcl)%ptr, stat=alloc_stat)
                                call build_glob%spproj_field%set(iptcl,'updatecnt',real(updatecnt-1))
                            endif
                        endif
                    else
                        ! offline mode, based on iteration
                        if( l_greedy .or. (.not.build_glob%spproj_field%has_been_searched(iptcl) .or. updatecnt==1) )then
                            if( trim(params_glob%tseries).eq.'yes' .and. l_np_cls_defined )then
                                allocate(strategy2D_tseries :: strategy2Dsrch(iptcl)%ptr, stat=alloc_stat)
                            else
                                allocate(strategy2D_greedy  :: strategy2Dsrch(iptcl)%ptr, stat=alloc_stat)
                            endif
                        else
                            allocate(strategy2D_snhc :: strategy2Dsrch(iptcl)%ptr, stat=alloc_stat)
                        endif
                    endif
                    if(alloc_stat/=0)call allocchk("In strategy2D_matcher:: cluster2D_exec strategy2Dsrch objects ")
                endif
            enddo
        end select
        ! actual construction
        cnt = 0
        do iptcl = params_glob%fromp, params_glob%top
            if( ptcl_mask(iptcl) )then
                cnt = cnt + 1
                ! search spec
                strategy2Dspec%iptcl       = iptcl
                strategy2Dspec%iptcl_map   = cnt
                strategy2Dspec%chunk_id    = chunk_id(iptcl)
                strategy2Dspec%stoch_bound = snhc_sz(iptcl)
                ! search object
                call strategy2Dsrch(iptcl)%ptr%new(strategy2Dspec)
            endif
        end do
        ! memoize CTF matrices
        if( build_glob%spproj%get_ctfflag('ptcl2D').ne.'no' )then
            call pftcc%create_polar_absctfmats(build_glob%spproj, 'ptcl2D')
        endif
        ! memoize FFTs for improved performance
        call pftcc%memoize_ffts
        ! SEARCH
        if( L_BENCH ) t_align = tic()
        !$omp parallel do default(shared) schedule(guided) private(i,iptcl) proc_bind(close)
        do i=1,nptcls2update
            iptcl = pinds(i)
            call strategy2Dsrch(iptcl)%ptr%srch
        end do
        !$omp end parallel do
        ! cleanup
        call clean_strategy2D
        do i=1,nptcls2update
            iptcl = pinds(i)
            call strategy2Dsrch(iptcl)%ptr%kill
            nullify(strategy2Dsrch(iptcl)%ptr)
        end do
        deallocate(strategy2Dsrch, pinds)

        ! CLEAN-UP
        call pftcc%kill
        deallocate(ptcl_mask)
        if( L_BENCH ) rt_align = toc(t_align)

        ! OUTPUT ORIENTATIONS
        call binwrite_oritab(params_glob%outfile, build_glob%spproj, build_glob%spproj_field, &
            &[params_glob%fromp,params_glob%top], isegment=PTCL2D_SEG)
        params_glob%oritab = params_glob%outfile

        ! WIENER RESTORATION OF CLASS AVERAGES
        if( L_BENCH ) t_cavg = tic()
        call cavger_transf_oridat( build_glob%spproj )
        call cavger_assemble_sums( l_partial_sums )
        ! write results to disk
        call cavger_readwrite_partial_sums('write')
        call cavger_kill
        if( L_BENCH ) rt_cavg = toc(t_cavg)
        call qsys_job_finished('simple_strategy2D_matcher :: cluster2D_exec')
        if( L_BENCH )then
            rt_tot  = toc(t_tot)
            doprint = .true.
            if( params_glob%part /= 1 ) doprint = .false.
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
    end subroutine cluster2D_exec

    !>  \brief  prepares the polarft corrcalc object for search
    subroutine preppftcc4align( which_iter )
        use simple_polarizer,           only: polarizer
        use simple_strategy2D3D_common, only: prep2dref,build_pftcc_particles
        integer,       intent(in)    :: which_iter
        type(polarizer), allocatable :: match_imgs(:)
        real      :: xyz(3)
        integer   :: icls, pop, pop_even, pop_odd, imatch, batchsz_max
        logical   :: do_center, has_been_searched
        ! create the polarft_corrcalc object
        call pftcc%new(params_glob%ncls, [params_glob%fromp,params_glob%top], ptcl_mask,&
            &eoarr=nint(build_glob%spproj_field%get_all('eo',[params_glob%fromp,params_glob%top])))
        ! prepare the polarizer images
        call build_glob%img_match%init_polarizer(pftcc, params_glob%alpha)
        ! this is tro avoid excessive allocation, allocate what is the upper bound on the
        ! # matchimgs needed for both parallel loops
        batchsz_max = max(MAXIMGBATCHSZ,params_glob%ncls)
        allocate(match_imgs(batchsz_max))
        do imatch=1,batchsz_max
            call match_imgs(imatch)%new([params_glob%boxmatch, params_glob%boxmatch, 1], params_glob%smpd)
            call match_imgs(imatch)%copy_polarizer(build_glob%img_match)
        end do
        ! PREPARATION OF REFERENCES IN PFTCC
        has_been_searched = .not.build_glob%spproj%is_virgin_field(params_glob%oritype)
        ! read references and transform into polar coordinates
        !$omp parallel do default(shared) private(icls,pop,pop_even,pop_odd,do_center,xyz)&
        !$omp schedule(static) proc_bind(close)
        do icls=1,params_glob%ncls
            pop      = 1
            pop_even = 0
            pop_odd  = 0
            if( has_been_searched )then
                pop      = build_glob%spproj_field%get_pop(icls, 'class'      )
                pop_even = build_glob%spproj_field%get_pop(icls, 'class', eo=0)
                pop_odd  = build_glob%spproj_field%get_pop(icls, 'class', eo=1)
            endif
            if( pop > 0 )then
                ! prepare the references
                ! here we are determining the shifts and map them back to classes
                do_center = (has_been_searched .and. (pop > MINCLSPOPLIM) .and. (which_iter > 2)&
                    &.and. .not.params_glob%l_frac_update)
                call prep2Dref(pftcc, cavgs_merged(icls), match_imgs(icls), icls, center=do_center, xyz_out=xyz)
                if( .not.params_glob%l_lpset )then
                    if( pop_even >= MINCLSPOPLIM .and. pop_odd >= MINCLSPOPLIM )then
                        ! here we are passing in the shifts and do NOT map them back to classes
                        call prep2Dref(pftcc, cavgs_even(icls), match_imgs(icls), icls, center=do_center, xyz_in=xyz)
                        call match_imgs(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.true., mask=build_glob%l_resmsk)  ! 2 polar coords
                        ! here we are passing in the shifts and do NOT map them back to classes
                        call prep2Dref(pftcc, cavgs_odd(icls), match_imgs(icls), icls, center=do_center, xyz_in=xyz)
                        call match_imgs(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.false., mask=build_glob%l_resmsk) ! 2 polar coords
                    else
                        ! put the merged class average in both even and odd positions
                        call match_imgs(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.true., mask=build_glob%l_resmsk) ! 2 polar coords
                        call pftcc%cp_even2odd_ref(icls)
                    endif
                else
                    call prep2Dref(pftcc, cavgs_merged(icls), match_imgs(icls), icls, center=do_center, xyz_in=xyz)
                    call match_imgs(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.true., mask=build_glob%l_resmsk) ! 2 polar coords
                    call pftcc%cp_even2odd_ref(icls)
                endif
            endif
        end do
        !$omp end parallel do
        ! PREPARATION OF PARTICLES IN PFTCC
        call build_pftcc_particles( pftcc, batchsz_max, match_imgs, ptcl_mask)
        ! DESTRUCT
        do imatch=1,batchsz_max
            call match_imgs(imatch)%kill_polarizer
            call match_imgs(imatch)%kill
            call build_glob%imgbatch(imatch)%kill
        end do
        deallocate(match_imgs, build_glob%imgbatch)
    end subroutine preppftcc4align

end module simple_strategy2D_matcher
