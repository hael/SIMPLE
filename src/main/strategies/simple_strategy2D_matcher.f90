! projection-matching based on Hadamard products, high-level search routines for CLUSTER2D
module simple_strategy2D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_builder,                only: build_glob
use simple_cmdline,                only: cmdline
use simple_euclid_sigma2,          only: euclid_sigma2
use simple_image,                  only: image
use simple_parameters,             only: params_glob
use simple_polarft_calc,       only: polarft_calc
use simple_qsys_funs,              only: qsys_job_finished
use simple_strategy2D,             only: strategy2D, strategy2D_per_ptcl
use simple_strategy2D3D_common,    only: set_bp_range2d, prepimgbatch, killimgbatch
use simple_strategy2D_greedy,      only: strategy2D_greedy
use simple_strategy2D_greedy_smpl, only: strategy2D_greedy_smpl
use simple_strategy2D_inpl,        only: strategy2D_inpl
use simple_strategy2D_inpl_smpl,   only: strategy2D_inpl_smpl
use simple_strategy2D_prob,        only: strategy2D_prob
use simple_strategy2D_snhc,        only: strategy2D_snhc
use simple_strategy2D_snhc_smpl,   only: strategy2D_snhc_smpl
use simple_strategy2D_srch,        only: strategy2D_spec
use simple_strategy2D_tseries,     only: strategy2D_tseries
use simple_binoris_io
use simple_classaverager
use simple_polarops
use simple_progress
use simple_strategy2D_alloc
implicit none

public :: cluster2D_exec
public :: sample_ptcls4update2D, preppftc4align2D, prep_batch_particles2D
public :: build_batch_particles2D, clean_batch_particles2D
private
#include "simple_local_flags.inc"

type(polarft_calc)   :: pftc
type(euclid_sigma2)      :: eucl_sigma
type(image), allocatable :: ptcl_match_imgs(:)
real(timer_int_kind)     :: rt_init, rt_prep_pftc, rt_align, rt_cavg, rt_projio, rt_tot
integer(timer_int_kind)  ::  t_init,  t_prep_pftc,  t_align,  t_cavg,  t_projio,  t_tot
type(string)             :: benchfname

contains

    !>  \brief  is the prime2D algorithm
    subroutine cluster2D_exec( cline, which_iter, converged )
        use simple_convergence,    only: convergence
        use simple_eul_prob_tab2D, only: eul_prob_tab2D
        use simple_decay_funs,     only: inv_cos_decay, extremal_decay2D
        class(cmdline),          intent(inout) :: cline
        integer,                 intent(in)    :: which_iter
        logical,                 intent(inout) :: converged
        type(strategy2D_per_ptcl), allocatable :: strategy2Dsrch(:)
        character(len=STDLEN)                  :: refine_flag
        real,                      allocatable :: states(:), incr_shifts(:,:)
        integer,                   allocatable :: pinds(:), batches(:,:)
        type(eul_prob_tab2D),           target :: probtab
        type(ori)             :: orientation
        type(convergence)     :: conv
        type(strategy2D_spec) :: strategy2Dspec
        real    :: frac_srch_space, neigh_frac, clinw
        integer :: iptcl, fnr, updatecnt, iptcl_map, iptcl_batch, ibatch, nptcls2update
        integer :: batchsz_max, batchsz, nbatches, batch_start, batch_end
        logical :: l_partial_sums, l_update_frac, l_ctf, l_prob, l_snhc, l_polar
        logical :: l_stream, l_greedy, l_np_cls_defined, l_alloc_read_cavgs, l_clin
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = build_glob%spproj_field%get_avg('frac')

        ! SWITCHES
        refine_flag    = trim(params_glob%refine)
        l_snhc         = str_has_substr(refine_flag, 'snhc')
        l_greedy       = str_has_substr(refine_flag, 'greedy')
        l_prob         = str_has_substr(refine_flag, 'prob')
        l_stream       = trim(params_glob%stream).eq.'yes'
        l_update_frac  = params_glob%l_update_frac  ! refers to particles sampling
        l_partial_sums = l_update_frac              ! to apply fractional momentum to class averages
        if( params_glob%extr_iter == 1 )then
            l_greedy       = .true.     ! greedy start
            l_snhc         = .false.
            l_partial_sums = .false.
        else if( params_glob%extr_iter > params_glob%extr_lim )then
            ! snhc_smpl turns to snhc after extremal phase
            if( trim(refine_flag)=='snhc_smpl' ) refine_flag = 'snhc'
        endif
        if( l_stream )then
            l_update_frac = .false.
            if( (which_iter>1) .and. (params_glob%update_frac<0.99) )then
                params_glob%l_update_frac = .true.
                l_partial_sums            = .true.
            else
                params_glob%update_frac   = 1.
                params_glob%l_update_frac = .false.
                l_partial_sums            = .false.
            endif
        endif
        if( l_prob )then
            ! all search decisions are made beforehand in prob_tab2D
            l_snhc         = .false.
            l_greedy       = .false.
            l_update_frac  = params_glob%l_update_frac
            l_partial_sums = l_update_frac .and. (params_glob%extr_iter>1)
        endif
        l_polar = trim(params_glob%polar).eq.'yes'
        l_clin  = .false.
        if( l_polar .and. trim(params_glob%ref_type)=='comlin_hybrid' )then
            if( l_snhc .or. (params_glob%extr_iter==1 .and.l_greedy))then
                l_clin =.true.
            else
                THROW_HARD('REF_TYPE=COMLIN_HYBRID only supported with refine=snhc|snhc_smpl')
            endif
        endif

        ! PARTICLE SAMPLING
        if( allocated(pinds) ) deallocate(pinds)
        if( l_prob )then
            ! generation of random sample and incr of updatecnts delegated to prob_tab2D_distr
            call build_glob%spproj_field%sample4update_reprod([params_glob%fromp,params_glob%top],&
            &nptcls2update, pinds )
        else
            call sample_ptcls4update2D([params_glob%fromp,params_glob%top], l_update_frac, nptcls2update, pinds)
        endif

        ! SNHC LOGICS
        neigh_frac = 0.
        if( params_glob%extr_iter > params_glob%extr_lim )then
            ! done
        else
            if( l_snhc )then
                neigh_frac = extremal_decay2D( params_glob%extr_iter, params_glob%extr_lim )
                if( L_VERBOSE_GLOB ) write(logfhandle,'(A,F8.2)') '>>> STOCHASTIC NEIGHBOURHOOD SIZE(%):', 100.*(1.-neigh_frac)
            endif
        endif

        ! READ FOURIER RING CORRELATIONS
        if( file_exists(params_glob%frcs) ) call build_glob%clsfrcs%read(params_glob%frcs)

        ! PREP REFERENCES
        if( build_glob%spproj_field%get_nevenodd() == 0 )then
            if( l_distr_exec_glob ) THROW_HARD('no eo partitioning available; cluster2D_exec')
            call build_glob%spproj_field%partition_eo
            call build_glob%spproj%write_segment_inside(params_glob%oritype)
        endif
        if( l_polar .and. which_iter>1 )then
            ! references are read in prep_polar_pftc4align2D below
            ! On first iteration the references are taken from the input images
        else
            l_alloc_read_cavgs = .true.
            if( .not.l_distr_exec_glob )then
                l_alloc_read_cavgs = which_iter==1
            endif
            call cavger_new(pinds, alloccavgs=l_alloc_read_cavgs)
            if( l_alloc_read_cavgs )then
                if( .not. cline%defined('refs') )then
                    THROW_HARD('need refs to be part of command line for cluster2D execution')
                endif
                call cavger_read_all
            endif
        endif

        ! SET FOURIER INDEX RANGE
        call set_bp_range2D(cline, which_iter, frac_srch_space)

        ! PREP BATCH ALIGNEMENT
        batchsz_max = min(nptcls2update,params_glob%nthr*BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls2update)/real(batchsz_max))
        batches     = split_nobjs_even(nptcls2update, nbatches)
        batchsz_max = maxval(batches(:,2)-batches(:,1)+1)
        allocate(incr_shifts(2,batchsz_max),source=0.)

        ! GENERATE POLAR REFERENCES
        if( L_BENCH_GLOB )then
            rt_init = toc(t_init)
            t_prep_pftc = tic()
        endif
        if( l_polar .and. which_iter>1)then
            ! Polar references, on first iteration the references are taken from the input images
            call prep_polar_pftc4align2D( pftc, batchsz_max, which_iter, l_stream )
        else
            ! Cartesian references
            call preppftc4align2D( pftc, batchsz_max, which_iter, l_stream )
        endif
        if( l_polar )then
            ! for restoration
            if( which_iter == 1 ) call polar_cavger_new(pftc, l_clin)
            call polar_cavger_zero_pft_refs
        endif

        ! ARRAY ALLOCATION FOR STRATEGY2D after pftc initialization
        call prep_strategy2D_glob( neigh_frac )
        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> STRATEGY2D OBJECTS ALLOCATED'

        ! SETUP WEIGHTS
        call build_glob%spproj_field%set_all2single('w', 1.0)

        ! GENERATE PARTICLES IMAGE OBJECTS
        allocate(strategy2Dsrch(batchsz_max))
        call prep_batch_particles2D(batchsz_max)
        if( L_BENCH_GLOB ) rt_prep_pftc = toc(t_prep_pftc)

        ! READ THE ASSIGNMENT FOR PROB MODE
        if( l_prob )then
            call probtab%new(pinds)
            call probtab%read_assignment(string(ASSIGNMENT_FBODY)//'.dat')
            s2D%probtab => probtab ! table accessible to strategies
        endif

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
            ! Prep particles in pftc
            if( L_BENCH_GLOB ) t_prep_pftc = tic()
            call build_batch_particles2D(pftc, batchsz, pinds(batch_start:batch_end), l_ctf)
            if( L_BENCH_GLOB ) rt_prep_pftc = rt_prep_pftc + toc(t_prep_pftc)
            ! batch strategy2D objects
            if( L_BENCH_GLOB ) t_init = tic()
            call prep_strategy2D_batch( pftc, which_iter, batchsz, pinds(batch_start:batch_end))
            if( L_BENCH_GLOB ) rt_init = rt_init + toc(t_init)
            ! Particles threaded loop
            if( L_BENCH_GLOB ) t_align = tic()
            !$omp parallel do private(iptcl,iptcl_batch,iptcl_map,updatecnt,orientation,strategy2Dspec)&
            !$omp default(shared) schedule(static) proc_bind(close)
            do iptcl_batch = 1,batchsz                     ! particle batch index
                iptcl_map  = batch_start + iptcl_batch - 1 ! masked global index (cumulative batch index)
                iptcl      = pinds(iptcl_map)              ! global index
                ! Search strategy (polymorphic strategy2D construction)
                updatecnt = build_glob%spproj_field%get_updatecnt(iptcl)
                if( l_stream )then
                    ! online mode, based on history
                    if( updatecnt==1 .or. (.not.build_glob%spproj_field%has_been_searched(iptcl)) )then
                        ! brand new particles
                        allocate(strategy2D_greedy                :: strategy2Dsrch(iptcl_batch)%ptr)
                    else
                        select case(trim(refine_flag))
                            case('greedy')
                                allocate(strategy2D_greedy        :: strategy2Dsrch(iptcl_batch)%ptr)
                            case('greedy_smpl')
                                allocate(strategy2D_greedy_smpl   :: strategy2Dsrch(iptcl_batch)%ptr)
                            case('snhc_smpl')
                                allocate(strategy2D_snhc_smpl     :: strategy2Dsrch(iptcl_batch)%ptr)
                            case DEFAULT ! is refine=snhc
                                allocate(strategy2D_snhc          :: strategy2Dsrch(iptcl_batch)%ptr)
                        end select
                    endif
                else
                    ! offline mode, based on iteration
                    if( l_prob )then
                        allocate(strategy2D_prob                    :: strategy2Dsrch(iptcl_batch)%ptr)
                    else
                        if( str_has_substr(refine_flag,'inpl') )then
                            if( refine_flag.eq.'inpl' )then
                                allocate(strategy2D_inpl            :: strategy2Dsrch(iptcl_batch)%ptr)
                            else if( refine_flag.eq.'inpl_smpl' )then
                                allocate(strategy2D_inpl_smpl       :: strategy2Dsrch(iptcl_batch)%ptr)
                            endif
                        else if( l_greedy .or. (updatecnt==1 .or. (.not.build_glob%spproj_field%has_been_searched(iptcl))) )then
                            ! first iteration | refine=*greedy*
                            if( trim(params_glob%tseries).eq.'yes' )then
                                if( l_np_cls_defined )then
                                    allocate(strategy2D_tseries     :: strategy2Dsrch(iptcl_batch)%ptr)
                                else
                                    allocate(strategy2D_greedy      :: strategy2Dsrch(iptcl_batch)%ptr)
                                endif
                            else
                                select case(trim(refine_flag))
                                case('greedy_smpl')
                                    allocate(strategy2D_greedy_smpl :: strategy2Dsrch(iptcl_batch)%ptr)
                                case DEFAULT ! is refine=greedy
                                    allocate(strategy2D_greedy      :: strategy2Dsrch(iptcl_batch)%ptr)
                                end select
                            endif
                        else
                            ! iteration>1 & refine/=*greedy*
                            select case(trim(refine_flag))
                                case('snhc_smpl')
                                    allocate(strategy2D_snhc_smpl   :: strategy2Dsrch(iptcl_batch)%ptr)
                                case DEFAULT ! is refine=snhc
                                    allocate(strategy2D_snhc        :: strategy2Dsrch(iptcl_batch)%ptr)
                            end select
                        endif
                    endif
                endif
                ! Search specification & object
                strategy2Dspec%iptcl       = iptcl
                strategy2Dspec%iptcl_batch = iptcl_batch
                strategy2Dspec%iptcl_map   = iptcl_map
                strategy2Dspec%stoch_bound = neigh_frac
                call strategy2Dsrch(iptcl_batch)%ptr%new(strategy2Dspec)
                call strategy2Dsrch(iptcl_batch)%ptr%srch
                ! keep track of incremental shift
                incr_shifts(:,iptcl_batch) = strategy2Dsrch(iptcl_batch)%ptr%s%best_shvec
                ! calculate sigma2 for ML-based refinement
                if ( params_glob%l_needs_sigma ) then
                    call build_glob%spproj_field%get_ori(iptcl, orientation)
                    call orientation%set_shift(incr_shifts(:,iptcl_batch)) ! incremental shift
                    call eucl_sigma%calc_sigma2(pftc, iptcl, orientation, 'class')
                end if
                ! cleanup
                call strategy2Dsrch(iptcl_batch)%ptr%kill
            enddo ! Particles threaded loop
            !$omp end parallel do
            if( L_BENCH_GLOB ) rt_align = rt_align + toc(t_align)
            ! restore polar cavgs
            if( l_polar )then
                call polar_cavger_update_sums(batchsz, pinds(batch_start:batch_end),&
                    &build_glob%spproj, pftc, incr_shifts(:,1:batchsz))
            endif
        enddo ! Batch loop

        ! BALANCING OF NUMBER OF PARTICLES PER CLASS
        if( l_prob )then
            ! done before, when assigning class from table
        else
            if( l_stream )then
                if( params_glob%l_update_frac .and. params_glob%maxpop>0 )then
                    call build_glob%spproj_field%balance_ptcls_within_cls(nptcls2update, pinds,&
                        &params_glob%maxpop, params_glob%nparts)
                endif
            else
                if( params_glob%maxpop>0 )then
                    call build_glob%spproj_field%balance_ptcls_within_cls(nptcls2update, pinds,&
                        &params_glob%maxpop, params_glob%nparts)
                endif
            endif
        endif

        ! CLEAN-UP
        call clean_strategy2D
        call orientation%kill
        call probtab%kill
        do iptcl_batch = 1,batchsz_max
            nullify(strategy2Dsrch(iptcl_batch)%ptr)
        end do
        call clean_batch_particles2D
        deallocate(strategy2Dsrch,pinds,batches)

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
        if( l_distr_exec_glob )then
            if( trim(params_glob%restore_cavgs).eq.'yes' )then
                if( l_polar )then
                    call polar_cavger_readwrite_partial_sums('write')
                else
                    call cavger_transf_oridat( build_glob%spproj )
                    call cavger_assemble_sums( l_partial_sums )
                    call cavger_readwrite_partial_sums('write')
                endif
            endif
            call cavger_kill
            call polar_cavger_kill
        else
            ! check convergence
            converged = conv%check_conv2D(cline, build_glob%spproj_field, build_glob%spproj_field%get_n('class'), params_glob%msk)
            converged = converged .and. (params_glob%which_iter >= params_glob%minits)
            converged = converged .or.  (params_glob%which_iter >= params_glob%maxits)
            ! Update progress file if not stream
            if(.not. l_stream) call progressfile_update(conv%get('progress'))
            if( trim(params_glob%restore_cavgs).eq.'yes' )then
                if( cline%defined('which_iter') )then
                    params_glob%refs      = CAVGS_ITER_FBODY//int2str_pad(params_glob%which_iter,3)//params_glob%ext%to_char()
                    params_glob%refs_even = CAVGS_ITER_FBODY//int2str_pad(params_glob%which_iter,3)//'_even'//params_glob%ext%to_char()
                    params_glob%refs_odd  = CAVGS_ITER_FBODY//int2str_pad(params_glob%which_iter,3)//'_odd'//params_glob%ext%to_char()
                else
                    THROW_HARD('which_iter expected to be part of command line in shared-memory execution')
                endif
                if( l_polar )then
                    if( which_iter == 1) call cavger_kill
                    ! polar restoration
                    if( l_clin )then
                        clinw = min(1.0, max(0.0, 1.0-max(0.0, real(params_glob%extr_iter-4)/real(params_glob%extr_lim-3))))
                        call polar_cavger_merge_eos_and_norm(build_glob%eulspace, clinw)
                    else
                        call polar_cavger_merge_eos_and_norm2D
                    endif
                    call polar_cavger_calc_and_write_frcs_and_eoavg(string(FRCS_FILE), cline)
                    call polar_cavger_writeall(string(POLAR_REFS_FBODY))
                    call polar_cavger_write_cartrefs(pftc, get_fbody(params_glob%refs,params_glob%ext,separator=.false.), 'merged')
                    call polar_cavger_gen2Dclassdoc(build_glob%spproj)
                    call polar_cavger_kill
                else
                    ! cartesian restoration
                    call cavger_transf_oridat( build_glob%spproj )
                    call cavger_assemble_sums( l_partial_sums )
                    call cavger_merge_eos_and_norm
                    call cavger_calc_and_write_frcs_and_eoavg(params_glob%frcs, params_glob%which_iter)
                    ! classdoc gen needs to be after calc of FRCs
                    call cavger_gen2Dclassdoc(build_glob%spproj)
                    ! write references
                    call cavger_write(params_glob%refs,'merged')
                    if( l_stream )then
                        call cavger_write(params_glob%refs_even,'even')
                        call cavger_write(params_glob%refs_odd, 'odd')
                        call cavger_readwrite_partial_sums('write')
                    endif
                    call cavger_kill(dealloccavgs=.false.)
                endif
                ! update command line
                call cline%set('refs', params_glob%refs)
                ! write project: cls2D and state congruent cls3D
                call build_glob%spproj%os_cls3D%new(params_glob%ncls, is_ptcl=.false.)
                states = build_glob%spproj%os_cls2D%get_all('state')
                call build_glob%spproj%os_cls3D%set_all('state',states)
                call build_glob%spproj%write_segment_inside('cls2D', params_glob%projfile)
                call build_glob%spproj%write_segment_inside('cls3D', params_glob%projfile)
                deallocate(states)
            endif
        endif
        call eucl_sigma%kill
        ! necessary for shared mem implementation, which otherwise bugs out when the bp-range changes
        call pftc%kill
        if( L_BENCH_GLOB ) rt_cavg = toc(t_cavg)
        call qsys_job_finished(string('simple_strategy2D_matcher :: cluster2D_exec'))
        if( L_BENCH_GLOB )then
            if( params_glob%part == 1 )then
                rt_tot  = toc(t_tot)
                benchfname = 'CLUSTER2D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
                call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation       : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'pftc preparation    : ', rt_prep_pftc
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', rt_align
                write(fnr,'(a,1x,f9.2)') 'class averaging      : ', rt_cavg
                write(fnr,'(a,1x,f9.2)') 'project file I/O     : ', rt_projio
                write(fnr,'(a,1x,f9.2)') 'total time           : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation       : ', (rt_init/rt_tot)       * 100.
                write(fnr,'(a,1x,f9.2)') 'pftc preparation    : ', (rt_prep_pftc/rt_tot) * 100.
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', (rt_align/rt_tot)      * 100.
                write(fnr,'(a,1x,f9.2)') 'class averaging      : ', (rt_cavg/rt_tot)       * 100.
                write(fnr,'(a,1x,f9.2)') 'project file I/O     : ', (rt_projio/rt_tot)     * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for      : ',&
                    &((rt_init+rt_prep_pftc+rt_align+rt_cavg+rt_projio)/rt_tot) * 100.
                call fclose(fnr)
            endif
        endif
    end subroutine cluster2D_exec

    subroutine sample_ptcls4update2D( pfromto, l_updatefrac, nptcls, pinds )
        logical,              intent(in)    :: l_updatefrac
        integer,              intent(in)    :: pfromto(2)
        integer,              intent(inout) :: nptcls
        integer, allocatable, intent(inout) :: pinds(:)
        if( l_updatefrac )then
            ! fractional sampling
            call build_glob%spproj_field%sample4update_rnd(pfromto, params_glob%update_frac,&
                &nptcls, pinds, .true.)
        else
            ! we sample all state > 0
            call build_glob%spproj_field%sample4update_all(pfromto, nptcls, pinds, .true.)
        endif
    end subroutine sample_ptcls4update2D

    !>  \brief  initializes convenience objects for particles polar alignment
    subroutine prep_batch_particles2D( batchsz_max )
        integer, intent(in) :: batchsz_max
        integer :: ithr
        call prepimgbatch(batchsz_max)
        allocate(ptcl_match_imgs(params_glob%nthr))
        !$omp parallel do private(ithr) default(shared) proc_bind(close) schedule(static)
        do ithr = 1,params_glob%nthr
            call ptcl_match_imgs(ithr)%new([params_glob%box_crop, params_glob%box_crop, 1],&
                &params_glob%smpd_crop, wthreads=.false.)
        enddo
        !$omp end parallel do
    end subroutine prep_batch_particles2D

    subroutine clean_batch_particles2D
        integer :: ithr
        call killimgbatch
        do ithr = 1,params_glob%nthr
            call ptcl_match_imgs(ithr)%kill
        enddo
        deallocate(ptcl_match_imgs)
    end subroutine clean_batch_particles2D

    !>  \brief  fills batch particle images for polar alignment
    subroutine build_batch_particles2D( pftc, nptcls_here, pinds, l_ctf_here )
        use simple_strategy2D3D_common, only: discrete_read_imgbatch, prepimg4align
        class(polarft_calc), intent(inout) :: pftc
        integer, intent(in) :: nptcls_here
        integer, intent(in) :: pinds(nptcls_here)
        logical, intent(in) :: l_ctf_here
        integer :: iptcl_batch, iptcl, ithr
        call discrete_read_imgbatch( nptcls_here, pinds, [1,nptcls_here])
        ! reassign particles indices & associated variables
        call pftc%reallocate_ptcls(nptcls_here, pinds)
        if( .not.build_glob%img_crop_polarizer%polarizer_initialized() )then
            call build_glob%img_crop_polarizer%init_polarizer(pftc, params_glob%alpha)
        endif
        ! call build_glob%spproj_field%write('ptcl2Dfield_going_in.txt')
        !$omp parallel do default(shared) private(iptcl,iptcl_batch,ithr)&
        !$omp schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds(iptcl_batch)
            call prepimg4align(iptcl, build_glob%imgbatch(iptcl_batch), ptcl_match_imgs(ithr))
            ! transfer to polar coordinates
            call build_glob%img_crop_polarizer%polarize(pftc, ptcl_match_imgs(ithr), iptcl, .true., .true., mask=build_glob%l_resmsk)
            ! e/o flag
            call pftc%set_eo(iptcl, nint(build_glob%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
        ! Memoize particles FFT parameters
        ! always create this one, CTF logic internal
        call pftc%create_polar_absctfmats(build_glob%spproj, 'ptcl2D')
        call pftc%memoize_ptcls
    end subroutine build_batch_particles2D

    !>  \brief  prepares the polarft corrcalc object for search and imports the references
    subroutine preppftc4align2D( pftc, batchsz_max, which_iter, l_stream )
        use simple_strategy2D3D_common, only: prep2dref
        class(polarft_calc), intent(inout) :: pftc
        integer,                 intent(in)    :: batchsz_max, which_iter
        logical,                 intent(in)    :: l_stream
        type(image),      allocatable :: match_imgs(:), tmp_imgs(:)
        type(string) :: fname
        real         :: xyz(3)
        integer      :: icls, pop, pop_even, pop_odd
        logical      :: do_center, has_been_searched
        has_been_searched = .not.build_glob%spproj%is_virgin_field(params_glob%oritype)
        ! create the polarft_calc object
        call pftc%new(params_glob%ncls, [1,batchsz_max], params_glob%kfromto)
        ! objective functions & sigma
        if( params_glob%l_needs_sigma )then
            fname = SIGMA2_FBODY//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'
            call eucl_sigma%new(fname, params_glob%box)
            if( l_stream )then
                call eucl_sigma%read_groups(build_glob%spproj_field)
                call eucl_sigma%allocate_ptcls
            else
                call eucl_sigma%read_part(  build_glob%spproj_field)
                if( params_glob%cc_objfun == OBJFUN_EUCLID ) call eucl_sigma%read_groups(build_glob%spproj_field)
            endif
        endif
        ! prepare the polarizer images
        call build_glob%img_crop_polarizer%init_polarizer(pftc, params_glob%alpha)
        allocate(match_imgs(params_glob%ncls),tmp_imgs(params_glob%ncls))
        call cavgs_merged(1)%construct_thread_safe_tmp_imgs(nthr_glob)
        ! PREPARATION OF REFERENCES IN pftc
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
                call match_imgs(icls)%new([params_glob%box_crop, params_glob%box_crop, 1], params_glob%smpd_crop, wthreads=.false.)
                call tmp_imgs(icls)%new([params_glob%box_crop, params_glob%box_crop, 1], params_glob%smpd_crop, wthreads=.false.)
                ! here we are determining the shifts and map them back to classes
                do_center = (has_been_searched .and. (pop > MINCLSPOPLIM) .and. (which_iter > 2)&
                    &.and. .not.params_glob%l_update_frac)
                call tmp_imgs(icls)%copy_fast(cavgs_merged(icls))
                call prep2Dref(tmp_imgs(icls), match_imgs(icls), icls, center=do_center, xyz_out=xyz)
                if( .not.params_glob%l_lpset )then
                    if( pop_even >= MINCLSPOPLIM .and. pop_odd >= MINCLSPOPLIM )then
                        ! here we are passing in the shifts and do NOT map them back to classes
                        call tmp_imgs(icls)%copy_fast(cavgs_even(icls))
                        call prep2Dref(tmp_imgs(icls), match_imgs(icls), icls, center=do_center, xyz_in=xyz)
                        call build_glob%img_crop_polarizer%polarize(pftc, match_imgs(icls), icls, isptcl=.false., iseven=.true., mask=build_glob%l_resmsk)  ! 2 polar coords
                        ! here we are passing in the shifts and do NOT map them back to classes
                        call tmp_imgs(icls)%copy_fast(cavgs_odd(icls))
                        call prep2Dref(tmp_imgs(icls), match_imgs(icls), icls, center=do_center, xyz_in=xyz)
                        call build_glob%img_crop_polarizer%polarize(pftc, match_imgs(icls), icls, isptcl=.false., iseven=.false., mask=build_glob%l_resmsk)  ! 2 polar coords
                    else
                        ! put the merged class average in both even and odd positions
                        call build_glob%img_crop_polarizer%polarize(pftc, match_imgs(icls), icls, isptcl=.false., iseven=.true., mask=build_glob%l_resmsk)  ! 2 polar coords
                        call pftc%cp_even2odd_ref(icls)
                    endif
                else
                    call tmp_imgs(icls)%copy_fast(cavgs_merged(icls))
                    call prep2Dref(cavgs_merged(icls), match_imgs(icls), icls, center=do_center, xyz_in=xyz)
                    call build_glob%img_crop_polarizer%polarize(pftc, match_imgs(icls), icls, isptcl=.false., iseven=.true., mask=build_glob%l_resmsk)  ! 2 polar coords
                    call pftc%cp_even2odd_ref(icls)
                endif
                call match_imgs(icls)%kill
                call tmp_imgs(icls)%kill
            endif
        end do
        !$omp end parallel do
        call pftc%memoize_refs
        ! CLEANUP
        deallocate(match_imgs,tmp_imgs)
        call cavgs_merged(1)%kill_thread_safe_tmp_imgs
    end subroutine preppftc4align2D

    !>  \brief  prepares the polarft corrcalc object for search and imports polar references
    subroutine prep_polar_pftc4align2D( pftc, batchsz_max, which_iter, l_stream )
        use simple_strategy2D3D_common, only: prep2dref
        class(polarft_calc), intent(inout) :: pftc
        integer,                 intent(in)    :: batchsz_max, which_iter
        logical,                 intent(in)    :: l_stream
        type(image),      allocatable :: tmp_imgs(:)
        type(string) :: fname
        real         :: xyz(3)
        integer      :: icls, pop, pop_even, pop_odd
        logical      :: has_been_searched, do_center, l_center
        ! pftc instantiation
        call pftc%new(params_glob%ncls, [1,batchsz_max], params_glob%kfromto)
        ! Sigma2
        if( params_glob%l_needs_sigma )then
            fname = SIGMA2_FBODY//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'
            call eucl_sigma%new(fname, params_glob%box)
            if( l_stream )then
                call eucl_sigma%read_groups(build_glob%spproj_field)
                call eucl_sigma%allocate_ptcls
            else
                call eucl_sigma%read_part(  build_glob%spproj_field)
                if( params_glob%cc_objfun == OBJFUN_EUCLID )then
                    call eucl_sigma%read_groups(build_glob%spproj_field)
                endif
            endif
        endif
        ! Read polar references
        call polar_cavger_new(pftc, trim(params_glob%ref_type)=='comlin_hybrid')
        call polar_cavger_read_all(string(POLAR_REFS_FBODY)//BIN_EXT)
        has_been_searched = .not.build_glob%spproj%is_virgin_field(params_glob%oritype)
        ! Centering-related objects
        do_center = (params_glob%center .eq. 'yes') .and. has_been_searched&
             &.and. (which_iter > 2) .and. (.not.params_glob%l_update_frac)
        if( do_center )then
            allocate(tmp_imgs(params_glob%ncls))
            call polar_cavger_refs2cartesian(pftc, tmp_imgs, 'merged')
            call tmp_imgs(1)%construct_thread_safe_tmp_imgs(nthr_glob)
        endif
        ! PREPARATION OF REFERENCES IN pftc
        !$omp parallel do default(shared) private(icls,pop,pop_even,pop_odd,xyz,l_center)&
        !$omp schedule(static) proc_bind(close)
        do icls=1,params_glob%ncls
            ! populations
            pop      = 1
            pop_even = 0
            pop_odd  = 0
            if( has_been_searched )then
                pop      = build_glob%spproj_field%get_pop(icls, 'class'      )
                pop_even = build_glob%spproj_field%get_pop(icls, 'class', eo=0)
                pop_odd  = build_glob%spproj_field%get_pop(icls, 'class', eo=1)
            endif
            if( pop > 0 )then
                ! centering
                l_center = do_center .and. (pop > MINCLSPOPLIM)
                if( l_center )then
                    call polar_prep2Dref( icls, cavg=tmp_imgs(icls), center=.true., xyz=xyz )
                else
                    call polar_prep2Dref( icls )
                    xyz = 0.0
                endif
                ! transfer to pftc
                if( .not.params_glob%l_lpset )then
                    if( pop_even >= MINCLSPOPLIM .and. pop_odd >= MINCLSPOPLIM )then
                        ! transfer e/o refs to pftc
                        call polar_cavger_set_ref_pftc(icls, 'even', pftc)
                        call polar_cavger_set_ref_pftc(icls, 'odd',  pftc)
                    else
                        ! put the merged class average in both even and odd positions
                        call polar_cavger_set_ref_pftc(icls, 'merged', pftc)
                        call pftc%cp_even2odd_ref(icls)
                    endif
                else
                    ! put the merged class average in both even and odd positions
                    call polar_cavger_set_ref_pftc(icls, 'merged', pftc)
                    call pftc%cp_even2odd_ref(icls)
                endif
                ! centering within the pftc
                if( l_center .and. (arg(xyz) > CENTHRESH) )then
                    call build_glob%spproj_field%add_shift2class(icls, -xyz(1:2))
                    call pftc%shift_ref(icls, xyz(1:2))
                endif
            endif
            if( do_center ) call tmp_imgs(icls)%kill
        end do
        !$omp end parallel do
        call pftc%memoize_refs
        ! cleanup
        if( do_center )then
            call tmp_imgs(1)%kill_thread_safe_tmp_imgs
            deallocate(tmp_imgs)
        endif
    end subroutine prep_polar_pftc4align2D

end module simple_strategy2D_matcher