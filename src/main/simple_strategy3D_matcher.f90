! projection-matching by stochastic hill-climbing, high-level search routines for REFINE3D
module simple_strategy3D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_strategy3D_alloc ! singleton s3D
use simple_timer
use simple_binoris_io
use simple_oris,                    only: oris
use simple_qsys_funs,               only: qsys_job_finished
use simple_kbinterpol,              only: kbinterpol
use simple_ori,                     only: ori
use simple_sym,                     only: sym
use simple_image,                   only: image
use simple_cmdline,                 only: cmdline
use simple_parameters,              only: params_glob
use simple_builder,                 only: build_glob
use simple_polarizer,               only: polarizer
use simple_polarft_corrcalc,        only: polarft_corrcalc
use simple_cartft_corrcalc,         only: cartft_corrcalc
use simple_strategy3D_cluster,      only: strategy3D_cluster
use simple_strategy3D_shc,          only: strategy3D_shc
use simple_strategy3D_shcc,         only: strategy3D_shcc
use simple_strategy3D_snhc,         only: strategy3D_snhc
use simple_strategy3D_snhcc,        only: strategy3D_snhcc
use simple_strategy3D_greedy,       only: strategy3D_greedy
use simple_strategy3D_greedy_neigh, only: strategy3D_greedy_neigh
use simple_strategy3D_neigh,        only: strategy3D_neigh
use simple_strategy3D_neighc,       only: strategy3D_neighc
use simple_strategy3D_cont,         only: strategy3D_cont
use simple_strategy3D,              only: strategy3D
use simple_strategy3D_srch,         only: strategy3D_spec, set_ptcl_stats, eval_ptcl
use simple_convergence,             only: convergence
use simple_euclid_sigma2,           only: euclid_sigma2
use simple_strategy2D3D_common
implicit none

public :: refine3D_exec, preppftcc4align, pftcc, calc_3Drec
private
#include "simple_local_flags.inc"

logical, parameter             :: DEBUG_HERE = .false.
logical                        :: has_been_searched
type(polarft_corrcalc), target :: pftcc
type(cartft_corrcalc),  target :: cftcc
type(polarizer),   allocatable :: match_imgs(:)
integer,           allocatable :: prev_states(:), pinds(:)
logical,           allocatable :: ptcl_mask(:)
type(sym)                      :: c1_symop
integer                        :: nptcls2update
type(euclid_sigma2)            :: eucl_sigma
! benchmarking
integer(timer_int_kind)        :: t_init,   t_prep_pftcc,  t_prep_orisrch,  t_align,  t_rec,  t_tot,  t_projio
real(timer_int_kind)           :: rt_init, rt_prep_pftcc, rt_prep_orisrch, rt_align, rt_rec, rt_tot, rt_projio
character(len=STDLEN)          :: benchfname

contains

    subroutine refine3D_exec( cline, which_iter, converged )
        class(cmdline),        intent(inout) :: cline
        integer,               intent(in)    :: which_iter
        logical,               intent(inout) :: converged
        integer, target, allocatable :: symmat(:,:)
        logical,         allocatable :: het_mask(:)
        !---> The below is to allow particle-dependent decision about which 3D strategy to use
        type :: strategy3D_per_ptcl
            class(strategy3D), pointer :: ptr  => null()
        end type strategy3D_per_ptcl
        type(strategy3D_per_ptcl), allocatable :: strategy3Dsrch(:)
        !<---- hybrid or combined search strategies can then be implemented as extensions of the
        !      relevant strategy3D base class
        type(strategy3D_spec),     allocatable :: strategy3Dspecs(:)
        real,                      allocatable :: resarr(:)
        integer,                   allocatable :: batches(:,:)
        character(len=:),          allocatable :: maps_dir, iter_dir
        type(convergence) :: conv
        type(ori)         :: orientation
        real    :: frac_srch_space, extr_thresh, extr_score_thresh, mi_proj, anneal_ratio
        integer :: nbatches, batchsz_max, batch_start, batch_end, batchsz, imatch
        integer :: iptcl, fnr, ithr, state, n_nozero, iptcl_batch, iptcl_map
        integer :: ibatch, iextr_lim, lpind_anneal, lpind_start, ncavgs
        logical :: doprint, do_extr, l_ctf, l_cartesian
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif

        ! CARTESIAN REFINEMENT FLAG
        select case(trim(params_glob%refine))
            case('shcc','neighc','snhcc')
                l_cartesian = .true.
            case DEFAULT
                l_cartesian = .false.
        end select
        
        ! CHECK THAT WE HAVE AN EVEN/ODD PARTITIONING
        if( build_glob%spproj_field%get_nevenodd() == 0 )then
            if( l_distr_exec_glob ) THROW_HARD('no eo partitioning available; refine3D_exec')
            call build_glob%spproj_field%partition_eo
            call build_glob%spproj%write_segment_inside(params_glob%oritype)
        endif

        ! CHECK WHETHER WE HAVE PREVIOUS 3D ORIENTATIONS
        has_been_searched = .not.build_glob%spproj%is_virgin_field(params_glob%oritype)

        ! SET FOURIER INDEX RANGE
        call set_bp_range(cline)

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = build_glob%spproj_field%get_avg('frac')

        ! PARTICLE INDEX SAMPLING FOR FRACTIONAL UPDATE (OR NOT)
        if( allocated(pinds) )     deallocate(pinds)
        if( allocated(ptcl_mask) ) deallocate(ptcl_mask)
        allocate(ptcl_mask(params_glob%fromp:params_glob%top))
        if( params_glob%l_frac_update )then
            call build_glob%spproj_field%sample4update_and_incrcnt([params_glob%fromp,params_glob%top],&
            &params_glob%update_frac, nptcls2update, pinds, ptcl_mask)
        else
            call build_glob%spproj_field%sample4update_and_incrcnt([params_glob%fromp,params_glob%top],&
            &1.0, nptcls2update, pinds, ptcl_mask)
        endif

        ! EXTREMAL LOGICS
        do_extr           = .false.
        extr_score_thresh = -huge(extr_score_thresh)
        select case(trim(params_glob%refine))
            case('cluster','clustersym')
                ! general logics
                if(allocated(het_mask))deallocate(het_mask)
                allocate(het_mask(params_glob%fromp:params_glob%top), source=ptcl_mask)
                call build_glob%spproj_field%set_extremal_vars(params_glob%extr_init, params_glob%extr_iter,&
                    &which_iter, frac_srch_space, do_extr, iextr_lim, update_frac=params_glob%update_frac)
                if( do_extr )then
                    anneal_ratio      = max(0., cos(PI/2.*real(params_glob%extr_iter-1)/real(iextr_lim)))
                    extr_thresh       = params_glob%extr_init * anneal_ratio
                    extr_score_thresh = build_glob%spproj_field%extremal_bound(extr_thresh, 'corr')
                    if( cline%defined('lpstart') )then
                        ! resolution limit update
                        lpind_start       = calc_fourier_index(params_glob%lpstart,params_glob%box,params_glob%smpd)
                        lpind_anneal      = nint(real(lpind_start) + (1.-anneal_ratio)*real(params_glob%kstop-lpind_start))
                        params_glob%kstop = min(lpind_anneal, params_glob%kstop)
                        resarr            = build_glob%img%get_res()
                        params_glob%lp    = resarr(params_glob%kstop)
                        if( params_glob%cc_objfun .ne. OBJFUN_EUCLID ) params_glob%kfromto(2) = params_glob%kstop
                        call build_glob%spproj_field%set_all2single('lp',params_glob%lp)
                        deallocate(resarr)
                    endif
                else
                    het_mask = .false.
                endif
                ! refinement mode specifics
                select case(trim(params_glob%refine))
                    case('clustersym')
                       ! symmetry pairing matrix
                        c1_symop = sym('c1')
                        params_glob%nspace = min(params_glob%nspace*build_glob%pgrpsyms%get_nsym(), 2500)
                        call build_glob%eulspace%new(params_glob%nspace, is_ptcl=.false.)
                        call build_glob%eulspace%spiral
                        call build_glob%pgrpsyms%nearest_sym_neighbors(build_glob%eulspace, symmat)
                end select
        end select
        ! PREP BATCH ALIGNEMENT
        batchsz_max = min(nptcls2update,params_glob%nthr*BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls2update)/real(batchsz_max))
        batches     = split_nobjs_even(nptcls2update, nbatches)
        batchsz_max = maxval(batches(:,2)-batches(:,1)+1)

        ! PREPARE THE POLARFT_CORRCALC DATA STRUCTURE
        if( L_BENCH_GLOB )then
            rt_init = toc(t_init)
            t_prep_pftcc = tic()
        endif

        if( l_cartesian )then
            call prepcftcc4align(cline, batchsz_max)
            ! pftcc still needs to be initiated due to build_glob%img_match%init_polarizer
            call pftcc%new(params_glob%nrefs, [1,batchsz_max], params_glob%l_match_filt)
        else
            call preppftcc4align(cline, batchsz_max)
        endif
        if( L_BENCH_GLOB ) rt_prep_pftcc = toc(t_prep_pftcc)
        if( L_BENCH_GLOB ) t_prep_orisrch = tic()
        ! clean big objects before starting to allocate new big memory chunks
        ! cannot kill build_glob%vol since used in continuous search
        call build_glob%vol2%kill
        ! call build_glob%vol_odd%kill ! cannot kill when we support cftcc refinement (Cartesian)
        ! array allocation for strategy3D
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: array allocation for strategy3D'
        call prep_strategy3D(ptcl_mask) ! allocate s3D singleton
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: array allocation for strategy3D, DONE'
        ! generate particles image objects
        call build_glob%img_match%init_polarizer(pftcc, params_glob%alpha)
        allocate(match_imgs(batchsz_max),strategy3Dspecs(batchsz_max),strategy3Dsrch(batchsz_max))
        call prepimgbatch(batchsz_max)
        !$omp parallel do default(shared) private(imatch) schedule(static) proc_bind(close)
        do imatch=1,batchsz_max
            call match_imgs(imatch)%new([params_glob%box, params_glob%box, 1], params_glob%smpd)
            call match_imgs(imatch)%copy_polarizer(build_glob%img_match)
        end do
        !$omp end parallel do

        ! STOCHASTIC IMAGE ALIGNMENT
        if( trim(params_glob%oritype) .eq. 'ptcl3D' )then
            l_ctf = build_glob%spproj%get_ctfflag('ptcl3D').ne.'no'
        else
            ! class averages have no CTF
            l_ctf = .false.
        endif
        write(logfhandle,'(A,1X,I3)') '>>> REFINE3D SEARCH, ITERATION:', which_iter
        if( L_BENCH_GLOB )then
            rt_prep_orisrch = toc(t_prep_orisrch)
            rt_align        = 0.
        endif
        ! Batch loop
        do ibatch=1,nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            ! Prep particles in pftcc
            if( L_BENCH_GLOB ) t_prep_pftcc = tic()
            if( l_cartesian )then
                call build_cftcc_batch_particles(batchsz, pinds(batch_start:batch_end))
                if( l_ctf ) call cftcc%create_absctfmats(build_glob%spproj, 'ptcl3D')
            else
                call build_pftcc_batch_particles(batchsz, pinds(batch_start:batch_end))
                if( l_ctf ) call pftcc%create_polar_absctfmats(build_glob%spproj, 'ptcl3D')
            endif
            if( L_BENCH_GLOB ) rt_prep_pftcc = rt_prep_pftcc + toc(t_prep_pftcc)
            ! Particles loop
            if( L_BENCH_GLOB ) t_align = tic()
            !$omp parallel do default(shared) private(iptcl,iptcl_batch,iptcl_map,ithr,orientation)&
            !$omp schedule(static) proc_bind(close)
            do iptcl_batch = 1,batchsz                     ! particle batch index
                iptcl_map  = batch_start + iptcl_batch - 1 ! masked global index (cumulative)
                iptcl      = pinds(iptcl_map)              ! global index
                ithr       = omp_get_thread_num() + 1
                ! switch for per-particle polymorphic strategy3D construction
                select case(trim(params_glob%refine))
                    case('snhc')
                        allocate(strategy3D_snhc                 :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('snhcc')
                        allocate(strategy3D_snhcc                :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('shc')
                        if( .not. has_been_searched )then
                            allocate(strategy3D_greedy           :: strategy3Dsrch(iptcl_batch)%ptr)
                        else
                            if( ran3() < GREEDY_FREQ )then
                                allocate(strategy3D_greedy       :: strategy3Dsrch(iptcl_batch)%ptr)
                            else
                                allocate(strategy3D_shc          :: strategy3Dsrch(iptcl_batch)%ptr)
                            endif
                        endif
                    case('shcc')
                        allocate(strategy3D_shcc                 :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('neigh')
                        if( ran3() < GLOB_FREQ )then
                            allocate(strategy3D_shc              :: strategy3Dsrch(iptcl_batch)%ptr)
                        else
                            if( ran3() < GREEDY_FREQ )then
                                allocate(strategy3D_greedy_neigh :: strategy3Dsrch(iptcl_batch)%ptr)
                            else
                                allocate(strategy3D_neigh        :: strategy3Dsrch(iptcl_batch)%ptr)
                            endif
                        endif
                    case('neighc')
                        allocate(strategy3D_neighc               :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('greedy')
                        allocate(strategy3D_greedy               :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('greedy_neigh')
                        if( ran3() < GLOB_FREQ )then
                            allocate(strategy3D_greedy           :: strategy3Dsrch(iptcl_batch)%ptr)
                        else
                            allocate(strategy3D_greedy_neigh     :: strategy3Dsrch(iptcl_batch)%ptr)
                        endif
                    case('cont')
                        THROW_HARD('refine=cont mode (continuous refinement in polar coordinates) is currently nonfunctional')
                        ! allocate(strategy3D_cont                 :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('cluster','clustersym')
                        allocate(strategy3D_cluster              :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('eval')
                        call eval_ptcl(pftcc, iptcl)
                        cycle
                    case DEFAULT
                        THROW_HARD('refinement mode: '//trim(params_glob%refine)//' unsupported')
                end select
                strategy3Dspecs(iptcl_batch)%iptcl =  iptcl
                strategy3Dspecs(iptcl_batch)%szsn  =  params_glob%szsn
                strategy3Dspecs(iptcl_batch)%extr_score_thresh = extr_score_thresh
                if( allocated(het_mask) ) strategy3Dspecs(iptcl_batch)%do_extr =  het_mask(iptcl)
                if( allocated(symmat)   ) strategy3Dspecs(iptcl_batch)%symmat  => symmat
                ! search object(s) & search
                call strategy3Dsrch(iptcl_batch)%ptr%new(strategy3Dspecs(iptcl_batch))
                call strategy3Dsrch(iptcl_batch)%ptr%srch(ithr)
                call strategy3Dsrch(iptcl_batch)%ptr%kill
                ! calculate sigma2 for ML-based refinement
                if ( params_glob%l_needs_sigma ) then
                    call build_glob%spproj_field%get_ori(iptcl, orientation)
                    if( orientation%isstatezero() ) cycle
                    call eucl_sigma%calc_sigma2(build_glob%spproj_field, iptcl, orientation)
                end if
            enddo ! Particles loop
            !$omp end parallel do
            if( L_BENCH_GLOB ) rt_align = rt_align + toc(t_align)
        enddo
        ! cleanup
        do iptcl_batch = 1,batchsz_max
            nullify(strategy3Dsrch(iptcl_batch)%ptr)
        end do
        deallocate(strategy3Dsrch,strategy3Dspecs,batches)

        ! WRITE SIGMAS FOR ML-BASED REFINEMENT
        if( params_glob%l_needs_sigma ) call eucl_sigma%write_sigma2

        ! UPDATE PARTICLE STATS
        if( .not. l_cartesian ) call calc_ptcl_stats( batchsz_max, l_ctf )

        ! CALCULATE PARTICLE WEIGHTS
        select case(trim(params_glob%ptclw))
            case('yes')
                ! weights are set at search time, so nothing to do here.
            case DEFAULT
                call build_glob%spproj_field%calc_hard_weights(params_glob%frac)
        end select

        ! CLEAN
        call clean_strategy3D ! deallocate s3D singleton
        call pftcc%kill
        call cftcc%kill
        call build_glob%vol%kill
        call build_glob%vol_odd%kill
        call orientation%kill
        do ibatch=1,batchsz_max
            call match_imgs(ibatch)%kill_polarizer
            call match_imgs(ibatch)%kill
        end do
        deallocate(match_imgs)
        if( allocated(symmat)   ) deallocate(symmat)
        if( allocated(het_mask) ) deallocate(het_mask)

        ! OUTPUT ORIENTATIONS
        if( L_BENCH_GLOB ) t_projio = tic()
        select case(trim(params_glob%oritype))
            case('ptcl3D')
                call binwrite_oritab(params_glob%outfile, build_glob%spproj, &
                    &build_glob%spproj_field, [params_glob%fromp,params_glob%top], isegment=PTCL3D_SEG)
            case('cls3D')
                call binwrite_oritab(params_glob%outfile, build_glob%spproj, &
                    &build_glob%spproj_field, [params_glob%fromp,params_glob%top], isegment=CLS3D_SEG)
            case DEFAULT
                THROW_HARD('unsupported oritype: '//trim(params_glob%oritype)//'; refine3D_exec')
        end select
        params_glob%oritab = params_glob%outfile
        if( L_BENCH_GLOB ) rt_projio = toc(t_projio)

        ! VOLUMETRIC 3D RECONSTRUCTION
        if( L_BENCH_GLOB ) t_rec = tic()
        call calc_3Drec( cline, which_iter )
        call eucl_sigma%kill
        if( L_BENCH_GLOB ) rt_rec = toc(t_rec)

        ! REPORT CONVERGENCE
        call qsys_job_finished(  'simple_strategy3D_matcher :: refine3D_exec')
        if( .not. params_glob%l_distr_exec ) converged = conv%check_conv3D(cline, params_glob%msk)
        if( L_BENCH_GLOB )then
            rt_tot  = toc(t_tot)
            doprint = .true.
            if( params_glob%part /= 1 ) doprint = .false.
            if( doprint )then
                benchfname = 'REFINE3D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
                call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation        : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'pftcc preparation     : ', rt_prep_pftcc
                write(fnr,'(a,1x,f9.2)') 'orisrch3D preparation : ', rt_prep_orisrch
                write(fnr,'(a,1x,f9.2)') '3D alignment          : ', rt_align
                write(fnr,'(a,1x,f9.2)') 'project file I/O      : ', rt_projio
                write(fnr,'(a,1x,f9.2)') 'reconstruction        : ', rt_rec
                write(fnr,'(a,1x,f9.2)') 'total time            : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation        : ', (rt_init/rt_tot)         * 100.
                write(fnr,'(a,1x,f9.2)') 'pftcc preparation     : ', (rt_prep_pftcc/rt_tot)   * 100.
                write(fnr,'(a,1x,f9.2)') 'orisrch3D preparation : ', (rt_prep_orisrch/rt_tot) * 100.
                write(fnr,'(a,1x,f9.2)') '3D alignment          : ', (rt_align/rt_tot)        * 100.
                write(fnr,'(a,1x,f9.2)') 'project file I/O      : ', (rt_projio/rt_tot)       * 100.
                write(fnr,'(a,1x,f9.2)') 'reconstruction        : ', (rt_rec/rt_tot)          * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for       : ',&
                    &((rt_init+rt_prep_pftcc+rt_prep_orisrch+rt_align+rt_projio+rt_rec)/rt_tot) * 100.
                call fclose(fnr)
            endif
        endif
    end subroutine refine3D_exec

    !> Prepare discrete search using polar projection Fourier cross correlation
    subroutine preppftcc4align( cline, batchsz_max )
        class(cmdline), intent(inout) :: cline !< command line
        integer,        intent(in)    :: batchsz_max
        type(ori) :: o_tmp
        real      :: xyz(3)
        integer   :: cnt, s, ind, iref, nrefs
        logical   :: do_center
        character(len=:), allocatable :: fname
        nrefs = params_glob%nspace * params_glob%nstates
        ! must be done here since params_glob%kfromto is dynamically set
        call pftcc%new(nrefs, [1,batchsz_max], params_glob%l_match_filt)
        if( params_glob%l_needs_sigma )then
            fname = SIGMA2_FBODY//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'
            call eucl_sigma%new(fname, params_glob%box)
            call eucl_sigma%read_part(  build_glob%spproj_field, ptcl_mask)
            call eucl_sigma%read_groups(build_glob%spproj_field, ptcl_mask)
        end if
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read reference volumes and create polar projections
        cnt = 0
        do s=1,params_glob%nstates
            if( str_has_substr(params_glob%refine,'greedy') )then
                if( .not.file_exists(params_glob%vols(s)) )then
                    cnt = cnt + params_glob%nspace
                    call progress(cnt, nrefs)
                    cycle
                endif
            else
                if( has_been_searched )then
                    if( build_glob%spproj_field%get_pop(s, 'state') == 0 )then
                        ! empty state
                        cnt = cnt + params_glob%nspace
                        call progress(cnt, nrefs)
                        cycle
                    endif
                endif
            endif
            call calcrefvolshift_and_mapshifts2ptcls( cline, s, params_glob%vols(s), do_center, xyz)
            call read_and_filter_refvols( cline, params_glob%vols_even(s), params_glob%vols_odd(s) )
            ! PREPARE ODD REFERENCES
            call preprefvol(pftcc, cline, s, do_center, xyz, .false.)
            !$omp parallel do default(shared) private(iref, o_tmp) schedule(static) proc_bind(close)
            do iref=1,params_glob%nspace
                call build_glob%eulspace%get_ori(iref, o_tmp)
                call build_glob%vol_odd%fproject_polar((s - 1) * params_glob%nspace + iref, &
                    &o_tmp, pftcc, iseven=.false., mask=build_glob%l_resmsk)
                call o_tmp%kill
            end do
            !$omp end parallel do
            ! PREPARE EVEN REFERENCES
            call preprefvol(pftcc,  cline, s, do_center, xyz, .true.)
            !$omp parallel do default(shared) private(iref, o_tmp) schedule(static) proc_bind(close)
            do iref=1,params_glob%nspace
                call build_glob%eulspace%get_ori(iref, o_tmp)
                call build_glob%vol%fproject_polar((s - 1) * params_glob%nspace + iref, &
                    &o_tmp, pftcc, iseven=.true., mask=build_glob%l_resmsk)
                call o_tmp%kill
            end do
            !$omp end parallel do
        end do
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align

    !> Prepare continuous search using Cartesian projection Fourier cross correlation
    subroutine prepcftcc4align( cline, batchsz_max )
        class(cmdline), intent(inout) :: cline !< command line
        integer,        intent(in)    :: batchsz_max
        real      :: xyz(3)
        integer   :: s
        logical   :: do_center
        character(len=:), allocatable :: fname
        ! must be done here since params_glob%kfromto is dynamically set
        call cftcc%new(build_glob%vol, build_glob%vol_odd, [1,batchsz_max], params_glob%l_match_filt)
        do s = 1,params_glob%nstates
            call calcrefvolshift_and_mapshifts2ptcls( cline, s, params_glob%vols(s), do_center, xyz)
            call read_and_filter_refvols(cline, params_glob%vols_even(s), params_glob%vols_odd(s))
            ! odd refvol
            call preprefvol(cftcc, cline, s, do_center, xyz, .false.)
            ! even refvol
            call preprefvol(cftcc, cline, s, do_center, xyz, .true.)
        end do
    end subroutine prepcftcc4align

    !>  \brief  prepares batch particle images for alignment
    subroutine build_pftcc_batch_particles( nptcls_here, pinds_here )
        use simple_strategy2D3D_common, only: read_imgbatch, prepimg4align
        integer, intent(in) :: nptcls_here
        integer, intent(in) :: pinds_here(nptcls_here)
        integer :: iptcl_batch, iptcl
        call read_imgbatch( nptcls_here, pinds_here, [1,nptcls_here] )
        ! reassign particles indices & associated variables
        call pftcc%reallocate_ptcls(nptcls_here, pinds_here)
        !$omp parallel do default(shared) private(iptcl,iptcl_batch) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            iptcl = pinds_here(iptcl_batch)
            ! prep
            call match_imgs(iptcl_batch)%zero_and_unflag_ft
            call prepimg4align(iptcl, build_glob%imgbatch(iptcl_batch), match_imgs(iptcl_batch))
            ! transfer to polar coordinates
            call match_imgs(iptcl_batch)%polarize(pftcc, iptcl, .true., .true., mask=build_glob%l_resmsk)
            ! e/o flag
            call pftcc%set_eo(iptcl, nint(build_glob%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
        ! Memoize particles FFT parameters
        call pftcc%memoize_ffts
    end subroutine build_pftcc_batch_particles

    !>  \brief  prepares batch particle images for alignment
    subroutine build_cftcc_batch_particles( nptcls_here, pinds_here )
        use simple_strategy2D3D_common, only: read_imgbatch, prepimg4align
        integer, intent(in) :: nptcls_here
        integer, intent(in) :: pinds_here(nptcls_here)
        integer :: iptcl_batch, iptcl
        call read_imgbatch( nptcls_here, pinds_here, [1,nptcls_here] )
        ! reassign particles indices & associated variables
        call cftcc%reallocate_ptcls(nptcls_here, pinds_here)
        !$omp parallel do default(shared) private(iptcl,iptcl_batch) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            iptcl = pinds_here(iptcl_batch)
            ! prep
            call match_imgs(iptcl_batch)%zero_and_unflag_ft
            call prepimg4align(iptcl, build_glob%imgbatch(iptcl_batch), match_imgs(iptcl_batch))
            call cftcc%set_ptcl(iptcl, match_imgs(iptcl_batch))
            ! e/o flag
            call cftcc%set_eo(iptcl, nint(build_glob%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
    end subroutine build_cftcc_batch_particles

    !> Prepare alignment search using polar projection Fourier cross correlation
    subroutine calc_ptcl_stats( batchsz_max, l_ctf )
        use simple_strategy2D3D_common, only: prepimg4align
        integer,   intent(in) :: batchsz_max
        logical,   intent(in) :: l_ctf
        integer, allocatable  :: pinds_here(:), batches(:,:)
        integer :: nptcls, iptcl_batch, iptcl, nbatches, ibatch, batch_start, batch_end, iptcl_map, batchsz
        if( .not.params_glob%l_frac_update ) return
        select case(params_glob%refine)
            case('cluster','clustersym','eval')
                return
            case DEFAULT
                ! all good
        end select
        ! build local particles index map
        nptcls = 0
        do iptcl = params_glob%fromp,params_glob%top
            if( .not.ptcl_mask(iptcl) )then
                if( build_glob%spproj_field%get_state(iptcl) > 0 ) nptcls = nptcls + 1
            endif
        enddo
        if( nptcls == 0 ) return
        allocate(pinds_here(nptcls),source=0)
        nptcls = 0
        do iptcl = params_glob%fromp,params_glob%top
            if( .not.ptcl_mask(iptcl) )then
                if( build_glob%spproj_field%get_state(iptcl) > 0 )then
                    nptcls = nptcls + 1
                    pinds_here(nptcls) = iptcl
                endif
            endif
        enddo
        ! Batch loop
        nbatches = ceiling(real(nptcls)/real(batchsz_max))
        batches  = split_nobjs_even(nptcls, nbatches)
        do ibatch=1,nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            call build_pftcc_batch_particles(batchsz, pinds_here(batch_start:batch_end))
            if( l_ctf ) call pftcc%create_polar_absctfmats(build_glob%spproj, 'ptcl3D')
            !$omp parallel do default(shared) private(iptcl,iptcl_batch,iptcl_map)&
            !$omp schedule(static) proc_bind(close)
            do iptcl_batch = 1,batchsz                     ! particle batch index
                iptcl_map  = batch_start + iptcl_batch - 1 ! masked global index (cumulative batch index)
                iptcl      = pinds_here(iptcl_map)         ! global index
                call set_ptcl_stats(pftcc, iptcl)
            enddo
            !$omp end parallel do
        enddo
    end subroutine calc_ptcl_stats

    !> volumetric 3d reconstruction
    subroutine calc_3Drec( cline, which_iter )
        use simple_fplane, only: fplane
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        type(fplane),    allocatable  :: fpls(:)
        type(ctfparams), allocatable  :: ctfparms(:)
        type(ori)        :: orientation
        type(kbinterpol) :: kbwin
        real    :: sdev_noise
        integer :: batchlims(2), iptcl, i, i_batch, ibatch
        if( trim(params_glob%dorec) .eq. 'no' ) return
        select case(trim(params_glob%refine))
            case('eval')
                ! no reconstruction
            case DEFAULT
                c1_symop = sym('c1')
                ! make the gridding prepper
                kbwin = build_glob%eorecvols(1)%get_kbwin()
                ! init volumes
                call preprecvols
                ! prep batch imgs
                call prepimgbatch(MAXIMGBATCHSZ)
                ! allocate array
                allocate(fpls(MAXIMGBATCHSZ),ctfparms(MAXIMGBATCHSZ))
                ! prep batch imgs
                call prepimgbatch(MAXIMGBATCHSZ)
                ! gridding batch loop
                do i_batch=1,nptcls2update,MAXIMGBATCHSZ
                    batchlims = [i_batch,min(nptcls2update,i_batch + MAXIMGBATCHSZ - 1)]
                    call read_imgbatch( nptcls2update, pinds, batchlims)
                    !$omp parallel do default(shared) private(i,iptcl,ibatch) schedule(static) proc_bind(close)
                    do i=batchlims(1),batchlims(2)
                        iptcl  = pinds(i)
                        ibatch = i - batchlims(1) + 1
                        if( .not.fpls(ibatch)%does_exist() ) call fpls(ibatch)%new(build_glob%imgbatch(1))
                        call build_glob%imgbatch(ibatch)%noise_norm(build_glob%lmsk, sdev_noise)
                        call build_glob%imgbatch(ibatch)%fft
                        ctfparms(ibatch) = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
                        call fpls(ibatch)%gen_planes(build_glob%imgbatch(ibatch), ctfparms(ibatch), iptcl=iptcl, serial=.true.)
                    end do
                    !$omp end parallel do
                    ! gridding
                    do i=batchlims(1),batchlims(2)
                        iptcl       = pinds(i)
                        ibatch      = i - batchlims(1) + 1
                        call build_glob%spproj_field%get_ori(iptcl, orientation)
                        if( orientation%isstatezero() ) cycle
                        select case(trim(params_glob%refine))
                            case('clustersym')
                                ! always C1 reconstruction
                                call grid_ptcl(fpls(ibatch), c1_symop, orientation)
                            case DEFAULT
                                call grid_ptcl(fpls(ibatch), build_glob%pgrpsyms, orientation)
                        end select
                    end do
                end do
                ! normalise structure factors
                call norm_struct_facts( cline, which_iter)
                ! destruct
                call killrecvols()
                do ibatch=1,MAXIMGBATCHSZ
                    call fpls(ibatch)%kill
                end do
                deallocate(fpls,ctfparms)
       end select
       call orientation%kill
    end subroutine calc_3Drec

end module simple_strategy3D_matcher
