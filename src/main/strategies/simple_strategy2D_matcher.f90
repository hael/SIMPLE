! projection-matching based on Hadamard products, high-level search routines for CLUSTER2D
module simple_strategy2D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_binoris_io
use simple_polarft_corrcalc,    only: polarft_corrcalc
use simple_cmdline,             only: cmdline
use simple_builder,             only: build_glob
use simple_parameters,          only: params_glob
use simple_polarizer,           only: polarizer
use simple_qsys_funs,           only: qsys_job_finished
use simple_convergence,         only: convergence
use simple_strategy2D3D_common, only: set_bp_range2d, prepimgbatch
use simple_strategy2D,          only: strategy2D, strategy2D_per_ptcl
use simple_strategy2D_srch,     only: strategy2D_spec
use simple_strategy2D_alloc,    only: prep_strategy2d_batch, clean_strategy2d, prep_strategy2D_glob
use simple_strategy2D_greedy,   only: strategy2D_greedy
use simple_strategy2D_tseries,  only: strategy2D_tseries
use simple_strategy2D_snhc,     only: strategy2D_snhc
use simple_strategy2D_eval,     only: strategy2D_eval
use simple_opt_filter,          only: uni_filt2D_sub
use simple_masker,              only: automask2D
use simple_euclid_sigma2,       only: euclid_sigma2
use simple_classaverager
use simple_progress
implicit none

public :: cluster2D_exec
private
#include "simple_local_flags.inc"

type(polarft_corrcalc)       :: pftcc
type(euclid_sigma2)          :: eucl_sigma
logical,         allocatable :: ptcl_mask(:)
integer                      :: batchsz_max
real(timer_int_kind)         :: rt_init, rt_prep_pftcc, rt_align, rt_cavg, rt_projio, rt_tot
integer(timer_int_kind)      ::  t_init,  t_prep_pftcc,  t_align,  t_cavg,  t_projio,  t_tot
character(len=STDLEN)        :: benchfname

contains

    !>  \brief  is the prime2D algorithm
    subroutine cluster2D_exec( cline, which_iter, converged )
        class(cmdline),          intent(inout) :: cline
        integer,                 intent(in)    :: which_iter
        logical,                 intent(inout) :: converged
        type(strategy2D_per_ptcl), allocatable :: strategy2Dsrch(:)
        type(strategy2D_spec),     allocatable :: strategy2Dspecs(:)
        integer,                   allocatable :: pinds(:), batches(:,:)
        real,                      allocatable :: states(:)
        type(convergence) :: conv
        type(cmdline)     :: cline_opt_filt
        type(ori)         :: orientation
        real    :: frac_srch_space, snhc_sz, frac
        integer :: iptcl, fnr, updatecnt, iptcl_map, nptcls2update, min_nsamples, icls
        integer :: batchsz, nbatches, batch_start, batch_end, iptcl_batch, ibatch
        logical :: doprint, l_partial_sums, l_frac_update, l_ctf, have_frcs
        logical :: l_snhc, l_greedy, l_stream, l_np_cls_defined
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = build_glob%spproj_field%get_avg('frac')

        ! SWITCHES
        l_partial_sums     = .false.
        l_snhc             = .false.
        l_greedy           = .false.
        l_frac_update      = .false.
        l_stream           = trim(params_glob%stream).eq.'yes'
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
            if( (params_glob%refine.eq.'greedy') )then
                l_greedy   = .true.
                l_snhc     = .false.
            endif
        else
            ! optional fractional update, no snhc opt
            l_partial_sums = params_glob%l_frac_update
            l_frac_update  = params_glob%l_frac_update
            l_snhc         = .false.
            l_greedy       = (params_glob%refine.eq.'greedy') .or.(params_glob%cc_objfun.eq.OBJFUN_EUCLID)
        endif
        if( l_stream )then
            l_frac_update             = .false.
            l_partial_sums            = .false.
            params_glob%l_frac_update = .false.
            if( which_iter > 1 )then
                if( params_glob%update_frac < 0.99 )then
                    l_partial_sums            = .true.
                    params_glob%l_frac_update = .true.
                else
                    params_glob%update_frac = 1.
                endif
            else
                params_glob%update_frac = 1.
            endif
        endif

        ! PARTICLE INDEX SAMPLING FOR FRACTIONAL UPDATE (OR NOT)
        if( allocated(pinds) )     deallocate(pinds)
        if( allocated(ptcl_mask) ) deallocate(ptcl_mask)
        allocate(ptcl_mask(params_glob%fromp:params_glob%top))
        if( l_frac_update )then
            call build_glob%spproj_field%sample4update_and_incrcnt([params_glob%fromp,params_glob%top],&
            &params_glob%update_frac, nptcls2update, pinds, ptcl_mask)
        else
            call build_glob%spproj_field%sample4update_and_incrcnt([params_glob%fromp,params_glob%top],&
            &1.0, nptcls2update, pinds, ptcl_mask)
        endif

        ! SNHC LOGICS
        if( l_snhc )then
            ! factorial decay, -2 because first step is always greedy
            snhc_sz = min(SNHC2D_INITFRAC,&
                &max(0.,SNHC2D_INITFRAC*(1.-SNHC2D_DECAY)**real(params_glob%extr_iter-2)))
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A,F8.2)') '>>> STOCHASTIC NEIGHBOURHOOD SIZE(%):', 100.*(1.-snhc_sz)
        else
            snhc_sz = 0. ! full neighbourhood
        endif

        ! ARRAY ALLOCATION FOR STRATEGY2D prior to weights
        call prep_strategy2D_glob
        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> STRATEGY2D OBJECTS ALLOCATED'

        ! SETUP WEIGHTS
        call build_glob%spproj_field%set_all2single('w', 1.0)

        ! READ FOURIER RING CORRELATIONS
        have_frcs = .false.
        if( file_exists(params_glob%frcs) )then
            call build_glob%clsfrcs%read(params_glob%frcs)
            have_frcs = .true.
        endif

        ! PREP REFERENCES
        call cavger_new(ptcl_mask)
        if( build_glob%spproj_field%get_nevenodd() == 0 )then
            if( l_distr_exec_glob ) THROW_HARD('no eo partitioning available; cluster2D_exec')
            call build_glob%spproj_field%partition_eo
            call build_glob%spproj%write_segment_inside(params_glob%oritype)
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
        if( params_glob%l_automsk .and. which_iter > AMSK2D_ITERLIM )then
            do icls = 1,params_glob%ncls
                call build_glob%env_masks(icls)%copy(cavgs_even(icls))
                call build_glob%env_masks(icls)%add(cavgs_odd(icls))
                call build_glob%env_masks(icls)%mul(0.5)
            enddo
            call automask2D(build_glob%env_masks, params_glob%ngrow, nint(params_glob%winsz), params_glob%edge, build_glob%diams)
            call uni_filt2D_sub(cavgs_even, cavgs_odd, build_glob%env_masks, build_glob%cutoff_finds_eo)
            if( params_glob%part.eq.1 )then
                do icls = 1,params_glob%ncls
                    call cavgs_even(icls)%write('filtered_refs_iter'//int2str_pad(which_iter,2)//'.mrc', icls)
                enddo
            endif
            ! read back the unfiltered class averages
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
            ! apply the envelope masks
            do icls = 1,params_glob%ncls
                call cavgs_even(icls)%mul(build_glob%env_masks(icls))
                call cavgs_odd(icls)%mul(build_glob%env_masks(icls))
            enddo
        endif

        ! SET FOURIER INDEX RANGE
        call set_bp_range2D(cline, which_iter, frac_srch_space)

        ! PREP BATCH ALIGNEMENT
        batchsz_max = min(nptcls2update,params_glob%nthr*BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls2update)/real(batchsz_max))
        batches     = split_nobjs_even(nptcls2update, nbatches)
        batchsz_max = maxval(batches(:,2)-batches(:,1)+1)

        ! GENERATE REFERENCES
        if( L_BENCH_GLOB )then
            rt_init = toc(t_init)
            t_prep_pftcc = tic()
        endif
        call preppftcc4align( which_iter )

        ! GENERATE PARTICLES IMAGE OBJECTS
        allocate(strategy2Dspecs(batchsz_max),strategy2Dsrch(batchsz_max))
        if( .not. params_glob%l_cartesian ) call build_glob%img_match%init_polarizer(pftcc, params_glob%alpha)
        call prepimgbatch(batchsz_max)
        if( L_BENCH_GLOB ) rt_prep_pftcc = toc(t_prep_pftcc)

        ! STOCHASTIC IMAGE ALIGNMENT
        rt_align         = 0.
        l_ctf            = build_glob%spproj%get_ctfflag('ptcl2D',iptcl=params_glob%fromp).ne.'no'
        l_np_cls_defined = cline%defined('nptcls_per_cls')
        write(logfhandle,'(A,1X,I3)') '>>> CLUSTER2D DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter
        ! Batch loop
        do ibatch=1,nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            ! Prep particles in pftcc
            if( L_BENCH_GLOB ) t_prep_pftcc = tic()
            call build_pftcc_batch_particles(batchsz, pinds(batch_start:batch_end))
            if( l_ctf ) call pftcc%create_polar_absctfmats(build_glob%spproj, 'ptcl2D')
            if( L_BENCH_GLOB ) rt_prep_pftcc = rt_prep_pftcc + toc(t_prep_pftcc)
            ! batch strategy2D objects
            if( L_BENCH_GLOB ) t_init = tic()
            call prep_strategy2D_batch( pftcc, which_iter, batchsz, pinds(batch_start:batch_end))
            if( L_BENCH_GLOB ) rt_init = rt_init + toc(t_init)
            ! Particles threaded loop
            if( L_BENCH_GLOB ) t_align = tic()
            !$omp parallel do default(shared) private(iptcl,iptcl_batch,iptcl_map,updatecnt)&
            !$omp schedule(static) proc_bind(close)
            do iptcl_batch = 1,batchsz                     ! particle batch index
                iptcl_map  = batch_start + iptcl_batch - 1 ! masked global index (cumulative batch index)
                iptcl      = pinds(iptcl_map)              ! global index
                ! Search strategy (polymorphic strategy2D construction)
                updatecnt = nint(build_glob%spproj_field%get(iptcl,'updatecnt'))
                if( l_stream )then
                    ! online mode, based on history
                    if( updatecnt==1 .or. (.not.build_glob%spproj_field%has_been_searched(iptcl)) )then
                        ! brand new particles
                        allocate(strategy2D_greedy          :: strategy2Dsrch(iptcl_batch)%ptr)
                    else
                        ! other particles
                        if( l_greedy )then
                            allocate(strategy2D_greedy      :: strategy2Dsrch(iptcl_batch)%ptr)
                        else
                            allocate(strategy2D_snhc        :: strategy2Dsrch(iptcl_batch)%ptr)
                        endif
                    endif
                else
                    ! offline mode, based on iteration
                    if( l_greedy .or. (.not.build_glob%spproj_field%has_been_searched(iptcl) .or. updatecnt==1) )then
                        if( trim(params_glob%tseries).eq.'yes' )then
                            if( l_np_cls_defined )then
                                allocate(strategy2D_tseries :: strategy2Dsrch(iptcl_batch)%ptr)
                            else
                                allocate(strategy2D_greedy  :: strategy2Dsrch(iptcl_batch)%ptr)
                            endif
                        else
                            allocate(strategy2D_greedy      :: strategy2Dsrch(iptcl_batch)%ptr)
                        endif
                    else
                        allocate(strategy2D_snhc            :: strategy2Dsrch(iptcl_batch)%ptr)
                    endif
                endif
                ! Search specification & object
                strategy2Dspecs(iptcl_batch)%iptcl       = iptcl
                strategy2Dspecs(iptcl_batch)%iptcl_map   = iptcl_batch
                strategy2Dspecs(iptcl_batch)%stoch_bound = snhc_sz
                call strategy2Dsrch(iptcl_batch)%ptr%new(strategy2Dspecs(iptcl_batch))
                call strategy2Dsrch(iptcl_batch)%ptr%srch
                ! calculate sigma2 for ML-based refinement
                if ( params_glob%l_needs_sigma ) then
                    call build_glob%spproj_field%get_ori(iptcl, orientation)
                    if( orientation%isstatezero() ) cycle
                    call eucl_sigma%calc_sigma2(build_glob%spproj_field, iptcl, orientation)
                end if
                ! cleanup
                call strategy2Dsrch(iptcl_batch)%ptr%kill
            enddo ! Particles threaded loop
            !$omp end parallel do
            if( L_BENCH_GLOB ) rt_align = rt_align + toc(t_align)
        enddo ! Batch loop

        ! CLEAN-UP
        call clean_strategy2D
        do iptcl_batch = 1,batchsz_max
            nullify(strategy2Dsrch(iptcl_batch)%ptr)
        end do
        deallocate(strategy2Dsrch,pinds,strategy2Dspecs,batches,ptcl_mask)

        ! WRITE SIGMAS FOR ML-BASED REFINEMENT
        if( params_glob%l_needs_sigma ) call eucl_sigma%write_sigma2

        ! OUTPUT ORIENTATIONS
        if( L_BENCH_GLOB ) t_projio = tic()
        call binwrite_oritab(params_glob%outfile, build_glob%spproj, build_glob%spproj_field, &
            &[params_glob%fromp,params_glob%top], isegment=PTCL2D_SEG)
        params_glob%oritab = params_glob%outfile
        if( L_BENCH_GLOB ) rt_projio = toc(t_projio)

        ! WIENER RESTORATION OF CLASS AVERAGES
        if( L_BENCH_GLOB ) t_cavg = tic()
        call cavger_transf_oridat( build_glob%spproj )
        call cavger_assemble_sums( l_partial_sums )
        if( l_distr_exec_glob )then
            ! write results to disk
            call cavger_readwrite_partial_sums('write')
        else
            ! check convergence
            converged = conv%check_conv2D(cline, build_glob%spproj_field, build_glob%spproj_field%get_n('class'), params_glob%msk)
            ! Update progress file if not stream
            if(.not. l_stream) call progressfile_update(conv%get('progress'))
            ! finalize classes
            if( converged .or. which_iter == params_glob%maxits)then
                call cavger_readwrite_partial_sums('write')
            endif
            call cavger_merge_eos_and_norm
            if( cline%defined('which_iter') )then
                params_glob%refs      = trim(CAVGS_ITER_FBODY)//int2str_pad(params_glob%which_iter,3)//params_glob%ext
                params_glob%refs_even = trim(CAVGS_ITER_FBODY)//int2str_pad(params_glob%which_iter,3)//'_even'//params_glob%ext
                params_glob%refs_odd  = trim(CAVGS_ITER_FBODY)//int2str_pad(params_glob%which_iter,3)//'_odd'//params_glob%ext
            else
                THROW_HARD('which_iter expected to be part of command line in shared-memory execution')
            endif
            call cavger_calc_and_write_frcs_and_eoavg(params_glob%frcs, params_glob%which_iter)
            ! classdoc gen needs to be after calc of FRCs
            call cavger_gen2Dclassdoc(build_glob%spproj)
            ! write references
            call cavger_write(trim(params_glob%refs),      'merged')
            call cavger_write(trim(params_glob%refs_even), 'even'  )
            call cavger_write(trim(params_glob%refs_odd),  'odd'   )
            ! update command line
            call cline%set('refs', trim(params_glob%refs))
            ! write project: cls2D and state congruent cls3D
            call build_glob%spproj%os_cls3D%new(params_glob%ncls, is_ptcl=.false.)
            states = build_glob%spproj%os_cls2D%get_all('state')
            call build_glob%spproj%os_cls3D%set_all('state',states)
            call build_glob%spproj%write_segment_inside('cls2D', params_glob%projfile)
            call build_glob%spproj%write_segment_inside('cls3D', params_glob%projfile)
            deallocate(states)
        endif
        call cavger_kill
        call eucl_sigma%kill
        call pftcc%kill ! necessary for shared mem implementation, which otherwise bugs out when the bp-range changes
        if( L_BENCH_GLOB ) rt_cavg = toc(t_cavg)
        call qsys_job_finished('simple_strategy2D_matcher :: cluster2D_exec')
        if( L_BENCH_GLOB )then
            rt_tot  = toc(t_tot)
            doprint = .true.
            if( params_glob%part /= 1 ) doprint = .false.
            if( doprint )then
                benchfname = 'CLUSTER2D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
                call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation       : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'pftcc preparation    : ', rt_prep_pftcc
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', rt_align
                write(fnr,'(a,1x,f9.2)') 'class averaging      : ', rt_cavg
                write(fnr,'(a,1x,f9.2)') 'project file I/O     : ', rt_projio
                write(fnr,'(a,1x,f9.2)') 'total time           : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation       : ', (rt_init/rt_tot)       * 100.
                write(fnr,'(a,1x,f9.2)') 'pftcc preparation    : ', (rt_prep_pftcc/rt_tot) * 100.
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', (rt_align/rt_tot)      * 100.
                write(fnr,'(a,1x,f9.2)') 'class averaging      : ', (rt_cavg/rt_tot)       * 100.
                write(fnr,'(a,1x,f9.2)') 'project file I/O     : ', (rt_projio/rt_tot)     * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for      : ',&
                    &((rt_init+rt_prep_pftcc+rt_align+rt_cavg+rt_projio)/rt_tot) * 100.
                call fclose(fnr)
            endif
        endif
    end subroutine cluster2D_exec

    !>  \brief  prepares batch particle images for alignment
    subroutine build_pftcc_batch_particles( nptcls_here, pinds )
        use simple_strategy2D3D_common, only: read_imgbatch, prepimg4align
        integer, intent(in) :: nptcls_here
        integer, intent(in) :: pinds(nptcls_here)
        integer :: iptcl_batch, iptcl
        call read_imgbatch( nptcls_here, pinds, [1,nptcls_here] )
        ! reassign particles indices & associated variables
        call pftcc%reallocate_ptcls(nptcls_here, pinds)
        !$omp parallel do default(shared) private(iptcl,iptcl_batch)&
        !$omp schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            iptcl = pinds(iptcl_batch)
            call prepimg4align(iptcl, build_glob%imgbatch(iptcl_batch))
            ! transfer to polar coordinates
            call build_glob%imgbatch(iptcl_batch)%polarize(pftcc, iptcl, .true., .true., mask=build_glob%l_resmsk)
            ! e/o flag
            call pftcc%set_eo(iptcl, nint(build_glob%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
        ! Memoize particles FFT parameters
        call pftcc%memoize_ffts
    end subroutine build_pftcc_batch_particles

    !>  \brief  prepares the polarft corrcalc object for search and imports the references
    subroutine preppftcc4align( which_iter )
        use simple_strategy2D3D_common, only: prep2dref
        integer,           intent(in) :: which_iter
        character(len=:), allocatable :: fname
        real      :: xyz(3)
        integer   :: icls, pop, pop_even, pop_odd
        logical   :: do_center, has_been_searched
        has_been_searched = .not.build_glob%spproj%is_virgin_field(params_glob%oritype)
        ! create the polarft_corrcalc object
        if( params_glob%l_bfac )then
            call pftcc%new(params_glob%ncls, [1,batchsz_max], params_glob%l_match_filt, bfac=params_glob%bfac)
        else
            call pftcc%new(params_glob%ncls, [1,batchsz_max], params_glob%l_match_filt)
        endif
        if( params_glob%l_needs_sigma )then
            fname = SIGMA2_FBODY//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'
            call eucl_sigma%new(fname, params_glob%box)
            call eucl_sigma%read_part(  build_glob%spproj_field, ptcl_mask)
            call eucl_sigma%read_groups(build_glob%spproj_field, ptcl_mask)
        end if
        ! prepare the polarizer images
        call build_glob%img_match%init_polarizer(pftcc, params_glob%alpha)
        ! PREPARATION OF REFERENCES IN PFTCC
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
                call cavgs_even(icls)%copy_polarizer(build_glob%img_match)
                call cavgs_odd(icls)%copy_polarizer(build_glob%img_match)
                call cavgs_merged(icls)%copy_polarizer(build_glob%img_match)
                ! prepare the references
                ! here we are determining the shifts and map them back to classes
                do_center = (has_been_searched .and. (pop > MINCLSPOPLIM) .and. (which_iter > 2)&
                    &.and. .not.params_glob%l_frac_update)
                call prep2Dref(cavgs_merged(icls), icls, .false., center=do_center, xyz_out=xyz)
                if( .not.params_glob%l_lpset )then
                    if( pop_even >= MINCLSPOPLIM .and. pop_odd >= MINCLSPOPLIM )then
                        ! here we are passing in the shifts and do NOT map them back to classes
                        call prep2Dref(cavgs_even(icls), icls, .true., center=do_center, xyz_in=xyz)
                        call cavgs_even(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.true., mask=build_glob%l_resmsk)  ! 2 polar coords
                        ! here we are passing in the shifts and do NOT map them back to classes
                        call prep2Dref(cavgs_odd(icls), icls, .false., center=do_center, xyz_in=xyz)
                        call cavgs_odd(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.false., mask=build_glob%l_resmsk) ! 2 polar coords
                    else
                        ! put the merged class average in both even and odd positions
                        call cavgs_merged(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.true., mask=build_glob%l_resmsk)  ! 2 polar coords
                        call pftcc%cp_even2odd_ref(icls)
                    endif
                else
                    call prep2Dref(cavgs_merged(icls), icls, .false., center=do_center, xyz_in=xyz)
                    call cavgs_merged(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.true., mask=build_glob%l_resmsk)     ! 2 polar coords
                    call pftcc%cp_even2odd_ref(icls)
                endif
            endif
        end do
        !$omp end parallel do
    end subroutine preppftcc4align

end module simple_strategy2D_matcher
