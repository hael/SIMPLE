module simple_strategy3D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_strategy3D_alloc ! singleton s3D
use simple_timer
use simple_binoris_io
use simple_qsys_funs,               only: qsys_job_finished
use simple_image,                   only: image
use simple_cmdline,                 only: cmdline
use simple_parameters,              only: params_glob
use simple_builder,                 only: build_glob
use simple_polarft_corrcalc,        only: polarft_corrcalc
use simple_cartft_corrcalc,         only: cartft_corrcalc
use simple_strategy3D_cluster,      only: strategy3D_cluster
use simple_strategy3D_shc,          only: strategy3D_shc
use simple_strategy3D_shcc,         only: strategy3D_shcc
use simple_strategy3D_snhc,         only: strategy3D_snhc
use simple_strategy3D_greedy,       only: strategy3D_greedy
use simple_strategy3D_greedyc,      only: strategy3D_greedyc
use simple_strategy3D_greedy_neigh, only: strategy3D_greedy_neigh
use simple_strategy3D_greedy_sub,   only: strategy3D_greedy_sub
use simple_strategy3D_shc_sub,      only: strategy3D_shc_sub
use simple_strategy3D_neigh,        only: strategy3D_neigh
use simple_strategy3D_neighc,       only: strategy3D_neighc
use simple_strategy3D,              only: strategy3D
use simple_strategy3D_srch,         only: strategy3D_spec
use simple_convergence,             only: convergence
use simple_euclid_sigma2,           only: euclid_sigma2
use simple_strategy2D3D_common
implicit none

public :: refine3D_exec, prep_ccobjs4align, pftcc
private
#include "simple_local_flags.inc"

logical, parameter             :: DEBUG_HERE = .false.
logical                        :: has_been_searched
type(polarft_corrcalc), target :: pftcc
type(cartft_corrcalc),  target :: cftcc
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
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: converged
        integer, target, allocatable  :: symmat(:,:)
        logical,         allocatable  :: het_mask(:)
        !---> The below is to allow particle-dependent decision about which 3D strategy to use
        type :: strategy3D_per_ptcl
            class(strategy3D), pointer :: ptr  => null() ! 1st polar greedy discrete search
            class(strategy3D), pointer :: ptr2 => null() ! 2nd Cartesian stochastic continuous search
        end type strategy3D_per_ptcl
        type(strategy3D_per_ptcl), allocatable :: strategy3Dsrch(:)
        !<---- hybrid or combined search strategies can then be implemented as extensions of the
        !      relevant strategy3D base class
        type(strategy3D_spec), allocatable :: strategy3Dspecs(:)
        real,                  allocatable :: resarr(:)
        integer,               allocatable :: batches(:,:)
        character(len=:),      allocatable :: maps_dir, iter_dir
        type(convergence) :: conv
        type(ori)         :: orientation
        real    :: frac_srch_space, extr_thresh, extr_score_thresh, mi_proj, anneal_ratio
        integer :: nbatches, batchsz_max, batch_start, batch_end, batchsz, imatch
        integer :: iptcl, fnr, ithr, state, n_nozero, iptcl_batch, iptcl_map
        integer :: ibatch, iextr_lim, lpind_anneal, lpind_start, ncavgs
        logical :: doprint, do_extr, l_ctf
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif

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
                        lpind_anneal      = nint(real(lpind_start) + (1.-anneal_ratio)*real(params_glob%kfromto(2)-lpind_start))
                        params_glob%kfromto(2) = min(lpind_anneal, params_glob%kfromto(2))
                        resarr            = build_glob%img%get_res()
                        params_glob%lp    = resarr(params_glob%kfromto(2))
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

        ! PREPARE THE FT_CORRCALC DATA STRUCTURES
        if( L_BENCH_GLOB )then
            rt_init = toc(t_init)
            t_prep_pftcc = tic()
        endif

        call prep_ccobjs4align(cline, batchsz_max)
        if( L_BENCH_GLOB ) rt_prep_pftcc = toc(t_prep_pftcc)
        if( L_BENCH_GLOB ) t_prep_orisrch = tic()
        ! clean big objects before starting to allocate new big memory chunks
        ! cannot kill build_glob%vol since used in continuous search
        ! cannot kill build_glob%vol_odd when we support cftcc refinement (Cartesian)
        call build_glob%vol2%kill
        ! array allocation for strategy3D
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: array allocation for strategy3D'
        call prep_strategy3D(ptcl_mask) ! allocate s3D singleton
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: array allocation for strategy3D, DONE'
        ! generate particles image objects
        allocate(strategy3Dspecs(batchsz_max),strategy3Dsrch(batchsz_max))
        if( pftcc%exists() )call build_glob%img_match%init_polarizer(pftcc, params_glob%alpha)
        call prepimgbatch(batchsz_max)

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
            call build_batch_particles(batchsz, pinds(batch_start:batch_end))
            call pftcc%create_polar_absctfmats(build_glob%spproj, 'ptcl3D')
            call cftcc%create_absctfmats(build_glob%spproj, 'ptcl3D')
            if( params_glob%l_obj_reg ) call pftcc%memoize_ptcl_prob(pinds(batch_start:batch_end))
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
                        if( .not. has_been_searched )then
                            allocate(strategy3D_greedyc          :: strategy3Dsrch(iptcl_batch)%ptr)
                        else
                            allocate(strategy3D_shcc             :: strategy3Dsrch(iptcl_batch)%ptr)
                        endif
                    case('neigh')
                        allocate(strategy3D_greedy_sub           :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('shc_neigh')
                        allocate(strategy3D_shc_sub              :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('neigh_test')
                        ! only do shifting in the ptr2
                        params_glob%l_doshift = .false.
                        allocate(strategy3D_greedy_sub           :: strategy3Dsrch(iptcl_batch)%ptr)
                        allocate(strategy3D_neighc               :: strategy3Dsrch(iptcl_batch)%ptr2)
                    case('neighc')
                        if( ran3() < GLOB_FREQ )then
                            allocate(strategy3D_shcc             :: strategy3Dsrch(iptcl_batch)%ptr)
                        else
                            allocate(strategy3D_neighc           :: strategy3Dsrch(iptcl_batch)%ptr)
                        endif
                    case('greedy')
                        allocate(strategy3D_greedy               :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('greedyc')
                        allocate(strategy3D_greedyc              :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('greedy_neigh')
                        if( ran3() < GLOB_FREQ )then
                            allocate(strategy3D_greedy           :: strategy3Dsrch(iptcl_batch)%ptr)
                        else
                            allocate(strategy3D_greedy_neigh     :: strategy3Dsrch(iptcl_batch)%ptr)
                        endif
                    case('cluster','clustersym')
                        allocate(strategy3D_cluster              :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('sigma')
                        ! first sigma estimation (done below)
                    case DEFAULT
                        THROW_HARD('refinement mode: '//trim(params_glob%refine)//' unsupported')
                end select
                strategy3Dspecs(iptcl_batch)%iptcl =  iptcl
                strategy3Dspecs(iptcl_batch)%szsn  =  params_glob%szsn
                strategy3Dspecs(iptcl_batch)%extr_score_thresh = extr_score_thresh
                if( allocated(het_mask) ) strategy3Dspecs(iptcl_batch)%do_extr =  het_mask(iptcl)
                if( allocated(symmat)   ) strategy3Dspecs(iptcl_batch)%symmat  => symmat
                ! search object(s) & search
                if( associated(strategy3Dsrch(iptcl_batch)%ptr) )then
                    call strategy3Dsrch(iptcl_batch)%ptr%new(strategy3Dspecs(iptcl_batch))
                    call strategy3Dsrch(iptcl_batch)%ptr%srch(ithr)
                    call strategy3Dsrch(iptcl_batch)%ptr%kill
                endif
                if( associated(strategy3Dsrch(iptcl_batch)%ptr2) )then
                    ! do shifting in the ptr2 now
                    params_glob%l_doshift = .true.
                    ! updating sigma with new orientation (same reference though), when ptr2 is set
                    call build_glob%spproj_field%get_ori(iptcl, orientation)
                    if( params_glob%l_cartesian )then
                        call eucl_sigma%update_sigma2(cftcc, iptcl, orientation, 'proj')
                    else
                        call eucl_sigma%update_sigma2(pftcc, iptcl, orientation, 'proj')
                    endif
                    call strategy3Dsrch(iptcl_batch)%ptr2%new(strategy3Dspecs(iptcl_batch))
                    call strategy3Dsrch(iptcl_batch)%ptr2%srch(ithr)
                    call strategy3Dsrch(iptcl_batch)%ptr2%kill
                endif
                ! calculate sigma2 for ML-based refinement
                if ( params_glob%l_needs_sigma ) then
                    call build_glob%spproj_field%get_ori(iptcl, orientation)
                    if( params_glob%l_cartesian )then
                        call eucl_sigma%calc_sigma2(cftcc, iptcl, orientation, 'proj')
                    else
                        call eucl_sigma%calc_sigma2(pftcc, iptcl, orientation, 'proj')
                    endif
                end if
            enddo ! Particles loop
            !$omp end parallel do
            if( L_BENCH_GLOB ) rt_align = rt_align + toc(t_align)
        enddo
        ! cleanup
        do iptcl_batch = 1,batchsz_max
            nullify(strategy3Dsrch(iptcl_batch)%ptr)
            nullify(strategy3Dsrch(iptcl_batch)%ptr2)
        end do
        deallocate(strategy3Dsrch,strategy3Dspecs,batches)

        ! WRITE SIGMAS FOR ML-BASED REFINEMENT
        if( params_glob%l_needs_sigma ) call eucl_sigma%write_sigma2

        ! CALCULATE PARTICLE WEIGHTS
        select case(trim(params_glob%ptclw))
            case('yes')
                ! weights are set at search time, so nothing to do here.
            case DEFAULT
                call build_glob%spproj_field%calc_hard_weights(params_glob%frac)
        end select

        ! CLEAN
        call clean_strategy3D ! deallocate s3D singleton
        if( params_glob%l_cartesian )then
            call cftcc%kill
        else
            call pftcc%kill
        endif
        call build_glob%vol%kill
        call orientation%kill
        if( allocated(symmat)   ) deallocate(symmat)
        if( allocated(het_mask) ) deallocate(het_mask)

        ! OUTPUT ORIENTATIONS
        select case(trim(params_glob%refine))
            case('sigma')
                ! nothing to do
            case DEFAULT
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
        end select

        ! VOLUMETRIC 3D RECONSTRUCTION
        select case(trim(params_glob%refine))
            case('eval','sigma')
                ! no reconstruction
            case DEFAULT
                if( L_BENCH_GLOB ) t_rec = tic()
                call calc_3Drec( cline, nptcls2update, pinds, which_iter )
                call eucl_sigma%kill
                call killimgbatch
                if( L_BENCH_GLOB ) rt_rec = toc(t_rec)
        end select

        ! REPORT CONVERGENCE
        call qsys_job_finished(  'simple_strategy3D_matcher :: refine3D_exec')
        if( .not. params_glob%l_distr_exec )then
            if( params_glob%l_cartesian )then
                converged = conv%check_conv3Dc(cline, params_glob%msk)
            else
                converged = conv%check_conv3D(cline, params_glob%msk)
            endif
        endif
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

    subroutine prep_ccobjs4align( cline, batchsz_max )
        class(cmdline), intent(inout) :: cline !< command line
        integer,        intent(in)    :: batchsz_max
        character(len=:), allocatable :: fname
        type(ori) :: o_tmp
        real      :: xyz(3)
        integer   :: cnt, s, ind, iref, nrefs
        logical   :: do_center
        ! first the polar
        nrefs = params_glob%nspace * params_glob%nstates
        ! must be done here since params_glob%kfromto is dynamically set
        call pftcc%new(nrefs, [1,batchsz_max], params_glob%kfromto)
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
            ! PREPARE E/O VOLUMES
            call preprefvol(cline, s, do_center, xyz, .false.)
            call preprefvol(cline, s, do_center, xyz, .true.)
            ! PREPARE REFERENCES
            !$omp parallel do default(shared) private(iref, o_tmp) schedule(static) proc_bind(close)
            do iref=1,params_glob%nspace
                call build_glob%eulspace%get_ori(iref, o_tmp)
                call build_glob%vol_odd%fproject_polar((s - 1) * params_glob%nspace + iref,&
                    &o_tmp, pftcc, iseven=.false., mask=build_glob%l_resmsk)
                call build_glob%vol%fproject_polar(    (s - 1) * params_glob%nspace + iref,&
                    &o_tmp, pftcc, iseven=.true.,  mask=build_glob%l_resmsk)
                call o_tmp%kill
            end do
        end do
        ! then the Cartesian
        ! must be done here since params_glob%kfromto is dynamically set
        call cftcc%new(build_glob%vol, build_glob%vol_odd, [1,batchsz_max])
        if( params_glob%l_needs_sigma )then
            fname = SIGMA2_FBODY//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'
            call eucl_sigma%new(fname, params_glob%box)
            call eucl_sigma%read_part(  build_glob%spproj_field, ptcl_mask)
            call eucl_sigma%read_groups(build_glob%spproj_field, ptcl_mask)
        end if
        if( params_glob%l_obj_reg )then
            call pftcc%calc_ptcl_reg(.true.)
            call pftcc%calc_ptcl_reg(.false.)
        endif
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: finished prep_ccobjs4align'
    end subroutine prep_ccobjs4align

    subroutine build_batch_particles( nptcls_here, pinds_here )
        use simple_strategy2D3D_common, only: read_imgbatch, prepimg4align
        integer, intent(in) :: nptcls_here
        integer, intent(in) :: pinds_here(nptcls_here)
        integer :: iptcl_batch, iptcl
        call read_imgbatch( nptcls_here, pinds_here, [1,nptcls_here] )
        ! reassign particles indices & associated variables
        call pftcc%reallocate_ptcls(nptcls_here, pinds_here)
        call cftcc%reallocate_ptcls(nptcls_here, pinds_here)
        !$omp parallel do default(shared) private(iptcl,iptcl_batch) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            iptcl = pinds_here(iptcl_batch)
            ! prep
            call prepimg4align(iptcl, build_glob%imgbatch(iptcl_batch))
            ! transfer to polar coordinates
            call build_glob%img_match%polarize(pftcc, build_glob%imgbatch(iptcl_batch), iptcl, .true., .true., mask=build_glob%l_resmsk)
            ! set Cartesian
            call cftcc%set_ptcl(iptcl, build_glob%imgbatch(iptcl_batch))
            ! e/o flags
            call pftcc%set_eo(iptcl, nint(build_glob%spproj_field%get(iptcl,'eo'))<=0 )
            call cftcc%set_eo(iptcl, nint(build_glob%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
        ! Memoize particles FFT parameters
        call pftcc%memoize_ffts
    end subroutine build_batch_particles

end module simple_strategy3D_matcher
