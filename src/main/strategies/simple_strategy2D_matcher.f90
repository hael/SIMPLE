!@descr: high-level search routines for the cluster2D and abinitio2D applications
module simple_strategy2D_matcher
use simple_pftc_srch_api
use simple_binoris_io
use simple_classaverager
use simple_new_classaverager
use simple_progress
use simple_strategy2D_alloc
use simple_builder,                only: build_glob
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
implicit none

public :: cluster2D_exec
public :: sample_ptcls4update2D, preppftc4align2D, prep_batch_particles2D
public :: build_batch_particles2D, clean_batch_particles2D
private
#include "simple_local_flags.inc"

type(polarft_calc)       :: pftc
type(euclid_sigma2)      :: eucl_sigma
type(image), allocatable :: ptcl_match_imgs(:), ptcl_match_imgs_pad(:)
real(timer_int_kind)     :: rt_init, rt_prep_pftc, rt_align, rt_cavg, rt_projio, rt_tot
integer(timer_int_kind)  :: t, t_init,  t_prep_pftc,  t_align,  t_cavg,  t_projio,  t_tot
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
            if( L_NEW_CAVGER )then
                call cavger_new_new(build_glob, pinds, alloccavgs=l_alloc_read_cavgs)
                if( l_alloc_read_cavgs )then
                    if( .not. cline%defined('refs') )then
                        THROW_HARD('need refs to be part of command line for cluster2D execution')
                    endif
                    call cavger_new_read_all
                endif
            else
                call cavger_new(build_glob, pinds, alloccavgs=l_alloc_read_cavgs)
                if( l_alloc_read_cavgs )then
                    if( .not. cline%defined('refs') )then
                        THROW_HARD('need refs to be part of command line for cluster2D execution')
                    endif
                    call cavger_read_all
                endif
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
        ! generate particles/references image objects
        call prep_batch_particles2D(batchsz_max)
        if( l_polar .and. which_iter>1)then
            ! Polar references, on first iteration the references are taken from the input images
            call prep_polar_pftc4align2D( pftc, batchsz_max, which_iter, l_stream )
        else
            ! Cartesian references
            call preppftc4align2D( pftc, batchsz_max, which_iter, l_stream )
        endif
        if( l_polar )then
            ! for restoration
            if( which_iter == 1 ) call pftc%polar_cavger_new(l_clin)
            call pftc%polar_cavger_zero_pft_refs
        endif

        ! ARRAY ALLOCATION FOR STRATEGY2D after pftc initialization
        call prep_strategy2D_glob( build_glob%spproj, neigh_frac )
        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> STRATEGY2D OBJECTS ALLOCATED'

        ! SETUP WEIGHTS
        call build_glob%spproj_field%set_all2single('w', 1.0)

        ! GENERATE PARTICLES SEARCH OBJECTS
        allocate(strategy2Dsrch(batchsz_max))
        if( L_BENCH_GLOB ) rt_prep_pftc = toc(t_prep_pftc)

        ! READ THE ASSIGNMENT FOR PROB MODE
        if( l_prob )then
            call probtab%new(build_glob, pinds)
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
            call build_batch_particles2D(pftc, batchsz, pinds(batch_start:batch_end))
            if( L_BENCH_GLOB ) rt_prep_pftc = rt_prep_pftc + toc(t_prep_pftc)
            ! batch strategy2D objects
            if( L_BENCH_GLOB ) t_init = tic()
            call prep_strategy2D_batch( build_glob%spproj, pftc, which_iter, batchsz, pinds(batch_start:batch_end))
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
                call pftc%polar_cavger_update_sums(batchsz, pinds(batch_start:batch_end),&
                    &build_glob%spproj, incr_shifts(:,1:batchsz))
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
                    call pftc%polar_cavger_readwrite_partial_sums('write')
                else
                    if( L_NEW_CAVGER )then
                        call cavger_new_transf_oridat( build_glob%spproj )
                        call cavger_new_assemble_sums( l_partial_sums )
                        call cavger_new_readwrite_partial_sums('write')
                    else
                        call cavger_transf_oridat( build_glob%spproj )
                        call cavger_assemble_sums( l_partial_sums )
                        call cavger_readwrite_partial_sums('write')
                    endif
                endif
            endif
            if( L_NEW_CAVGER )then
                call cavger_new_kill
            else
                call cavger_kill
            endif
            call pftc%polar_cavger_kill
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
                    if( L_NEW_CAVGER )then
                        if( which_iter == 1) call cavger_new_kill
                    else
                        if( which_iter == 1) call cavger_kill
                    endif
                    ! polar restoration
                    if( l_clin )then
                        clinw = min(1.0, max(0.0, 1.0-max(0.0, real(params_glob%extr_iter-4)/real(params_glob%extr_lim-3))))
                        call pftc%polar_cavger_merge_eos_and_norm(build_glob%eulspace, build_glob%pgrpsyms, clinw)
                    else
                        call pftc%polar_cavger_merge_eos_and_norm2D
                    endif
                    call pftc%polar_cavger_calc_and_write_frcs_and_eoavg(build_glob%clsfrcs, build_glob%spproj_field%get_update_frac(), string(FRCS_FILE), cline)
                    call pftc%polar_cavger_writeall(string(POLAR_REFS_FBODY))
                    call pftc%polar_cavger_gen2Dclassdoc(build_glob%spproj, build_glob%clsfrcs)
                    call pftc%polar_cavger_kill
                else
                    ! cartesian restoration
                    if( L_NEW_CAVGER )then
                        call cavger_new_transf_oridat( build_glob%spproj )
                        call cavger_new_assemble_sums( l_partial_sums )
                        call cavger_new_restore_cavgs( params_glob%frcs )
                        ! classdoc gen needs to be after calc of FRCs
                        call cavger_new_gen2Dclassdoc( build_glob%spproj )
                        ! write references
                        call cavger_new_write_merged( params_glob%refs )
                        if( l_stream )then
                            call cavger_new_write_eo( params_glob%refs_even, params_glob%refs_odd )
                            call cavger_new_readwrite_partial_sums( 'write' )
                        endif
                        call cavger_new_kill(dealloccavgs=.false.)
                    else
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
                write(fnr,'(a,1x,f9.2)') 'pftc preparation     : ', rt_prep_pftc
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', rt_align
                write(fnr,'(a,1x,f9.2)') 'class averaging      : ', rt_cavg
                write(fnr,'(a,1x,f9.2)') 'project file I/O     : ', rt_projio
                write(fnr,'(a,1x,f9.2)') 'total time           : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation       : ', (rt_init/rt_tot)      * 100.
                write(fnr,'(a,1x,f9.2)') 'pftc preparation     : ', (rt_prep_pftc/rt_tot) * 100.
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', (rt_align/rt_tot)     * 100.
                write(fnr,'(a,1x,f9.2)') 'class averaging      : ', (rt_cavg/rt_tot)      * 100.
                write(fnr,'(a,1x,f9.2)') 'project file I/O     : ', (rt_projio/rt_tot)    * 100.
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
        allocate(ptcl_match_imgs(params_glob%nthr), ptcl_match_imgs_pad(params_glob%nthr))
        !$omp parallel do private(ithr) default(shared) proc_bind(close) schedule(static)
        do ithr = 1,params_glob%nthr
            call ptcl_match_imgs(ithr)%new(    [params_glob%box_crop,   params_glob%box_crop,   1],&
            &params_glob%smpd_crop, wthreads=.false.)
            call ptcl_match_imgs_pad(ithr)%new([params_glob%box_croppd, params_glob%box_croppd, 1],&
            &params_glob%smpd_crop, wthreads=.false.)
        enddo
        !$omp end parallel do
    end subroutine prep_batch_particles2D

    subroutine clean_batch_particles2D
        use simple_imgarr_utils, only: dealloc_imgarr
        call killimgbatch
        call dealloc_imgarr(ptcl_match_imgs)
        call dealloc_imgarr(ptcl_match_imgs_pad)
    end subroutine clean_batch_particles2D

    !>  \brief  fills batch particle images for polar alignment
    subroutine build_batch_particles2D( pftc, nptcls_here, pinds )
        use simple_strategy2D3D_common, only: discrete_read_imgbatch, prepimg4align!, prepimg4align_bench
        class(polarft_calc), intent(inout) :: pftc
        integer,             intent(in)    :: nptcls_here
        integer,             intent(in)    :: pinds(nptcls_here)
        complex, allocatable :: pft(:,:)
        ! real(timer_int_kind)    :: rt_prep1, rt_prep2, rt_prep, rt_polarize, rt_sum, rt_loop
        ! integer(timer_int_kind) :: t_polarize, t_loop
        integer     :: iptcl_batch, iptcl, ithr
        call discrete_read_imgbatch( nptcls_here, pinds, [1,nptcls_here])
        ! reassign particles indices & associated variables
        call pftc%reallocate_ptcls(nptcls_here, pinds)
        ! memoization for polarize_oversamp
        call ptcl_match_imgs_pad(1)%memoize4polarize_oversamp(pftc%get_pdim())
        ! mask memoization for prepimg4align
        call ptcl_match_imgs(1)%memoize_mask_coords
        ! memoize FT mapping stuff
        call memoize_ft_maps(ptcl_match_imgs(1)%get_ldim(), ptcl_match_imgs(1)%get_smpd())
        ! allocate pft
        ! rt_prep1    = 0.
        ! rt_prep2    = 0.
        ! rt_prep     = 0.
        ! rt_polarize = 0.
        ! rt_sum      = 0.
        !$omp parallel do default(shared) private(iptcl,iptcl_batch,ithr,pft) schedule(static) proc_bind(close)
        ! t_loop = tic()
        do iptcl_batch = 1,nptcls_here
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds(iptcl_batch)
            call prepimg4align(iptcl, build_glob%imgbatch(iptcl_batch), ptcl_match_imgs(ithr), ptcl_match_imgs_pad(ithr))
            ! t_polarize = tic()
            ! call prepimg4align_bench(iptcl, build_glob%imgbatch(iptcl_batch), ptcl_match_imgs(ithr), ptcl_match_imgs_pad(ithr),&
            ! &rt_prep1, rt_prep2, rt_prep)
            ! t_polarize = tic()
            pft = pftc%allocate_pft()
            call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft, mask=build_glob%l_resmsk)
            call pftc%set_ptcl_pft(iptcl, pft)
            deallocate(pft)
            ! rt_polarize = rt_polarize + toc(t_polarize)
            ! e/o flag
            call pftc%set_eo(iptcl, nint(build_glob%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
        ! rt_loop = toc(t_loop)
        ! rt_sum = rt_prep + rt_polarize

        ! print *, 'rt_prep1    =', rt_prep1,    ' % ', 100.*(rt_prep1/rt_sum)
        ! print *, 'rt_prep2    =', rt_prep2,    ' % ', 100.*(rt_prep2/rt_sum)
        ! print *, 'rt_prep     =', rt_prep,     ' % ', 100.*(rt_prep/rt_sum)
        ! print *, 'rt_polarize =', rt_polarize, ' % ', 100.*(rt_polarize/rt_sum)
        ! print *, 'accounted for % ', 100.*(rt_sum/rt_loop)
        ! print *, ''

        ! rt_prep1    =   1.1228021550000002       %    47.656603255671037         
        ! rt_prep2    =  0.78998149100000004       %    33.530247807469209     
        ! rt_prep     =   2.2056463639999984       %    93.617217622334607     
        ! rt_polarize =  0.15038003799999988       %    6.3827823776653929     
        ! accounted for %    99.978044541257489  

        ! always create this one, CTF logic internal
        call pftc%create_polar_absctfmats(build_glob%spproj, 'ptcl2D')
        call pftc%memoize_ptcls
        ! destruct
        call forget_ft_maps
    end subroutine build_batch_particles2D

    !>  \brief  prepares the polarft corrcalc object for search and imports the references
    subroutine preppftc4align2D( pftc, batchsz_max, which_iter, l_stream )
        use simple_strategy2D3D_common, only: prep2dref, calc_2Dref_offset
        use simple_image,               only: unmemoize_mask_coords
        class(polarft_calc), intent(inout) :: pftc
        integer,             intent(in)    :: batchsz_max, which_iter
        logical,             intent(in)    :: l_stream
        class(image),    pointer :: cavgs_m(:), cavgs_e(:), cavgs_o(:)
        type(image), allocatable :: match_imgs(:)
        complex,     allocatable :: pft(:,:)
        type(string) :: fname
        real         :: xyz(3)
        integer      :: icls, pop, pop_even, pop_odd, centype, ithr
        logical      :: do_center, has_been_searched, input_center
        has_been_searched = .not.build_glob%spproj%is_virgin_field(params_glob%oritype)
        input_center      = trim(params_glob%center) .eq. 'yes'
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
        call ptcl_match_imgs_pad(1)%memoize4polarize_oversamp(pftc%get_pdim())
        allocate(match_imgs(params_glob%ncls))
        if( L_NEW_CAVGER )then
            cavgs_m => cavgs_merged_new
            cavgs_e => cavgs_even_new
            cavgs_o => cavgs_odd_new
        else
            cavgs_m => cavgs_merged
            cavgs_e => cavgs_even
            cavgs_o => cavgs_odd
        endif
        call cavgs_m(1)%construct_thread_safe_tmp_imgs(nthr_glob)
        ! mask memoization
        call cavgs_m(1)%memoize_mask_coords
        ! mode of cavg centering
        centype = get_centype(params_glob%center_type)
        ! PREPARATION OF REFERENCES IN pftc
        ! read references and transform into polar coordinates
        !$omp parallel do default(shared) private(icls,ithr,pop,pop_even,pop_odd,do_center,xyz,pft)&
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
            ithr = omp_get_thread_num() + 1
            if( pop > 0 )then
                call match_imgs(icls)%new([params_glob%box_crop, params_glob%box_crop, 1], params_glob%smpd_crop, wthreads=.false.)
                ! Calculate center
                do_center = (has_been_searched .and. (pop > MINCLSPOPLIM) .and. (which_iter > 2)&
                    &.and. .not.params_glob%l_update_frac)
                do_center = input_center .and. do_center
                if( do_center )then
                    call match_imgs(icls)%copy_fast(cavgs_m(icls))
                    call calc_2Dref_offset(match_imgs(icls), icls, centype, xyz)
                else
                    xyz = 0.0
                endif
                ! Prepare the references
                ! allocte pft
                pft = pftc%allocate_pft()
                if( params_glob%l_lpset )then
                    ! merged class average in both even and odd positions
                    call match_imgs(icls)%copy_fast(cavgs_m(icls))
                    call prep2Dref(match_imgs(icls), icls, xyz, ptcl_match_imgs_pad(ithr))
                    call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft, mask=build_glob%l_resmsk)
                    call pftc%set_ref_pft(icls, pft, iseven=.true.)
                    call pftc%cp_even2odd_ref(icls)
                else
                    if( pop_even >= MINCLSPOPLIM .and. pop_odd >= MINCLSPOPLIM )then
                        ! even & odd
                        call match_imgs(icls)%copy_fast(cavgs_e(icls))
                        call prep2Dref(match_imgs(icls), icls, xyz, ptcl_match_imgs_pad(ithr))
                        call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft, mask=build_glob%l_resmsk)
                        call pftc%set_ref_pft(icls, pft, iseven=.true.)
                        call match_imgs(icls)%copy_fast(cavgs_o(icls))
                        call prep2Dref(match_imgs(icls), icls, xyz, ptcl_match_imgs_pad(ithr))
                        call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft, mask=build_glob%l_resmsk)
                        call pftc%set_ref_pft(icls, pft, iseven=.false.)
                    else
                        ! merged class average in both even and odd positions
                        call match_imgs(icls)%copy_fast(cavgs_m(icls))
                        call prep2Dref(match_imgs(icls), icls, xyz, ptcl_match_imgs_pad(ithr))
                        call ptcl_match_imgs_pad(ithr)%polarize_oversamp(pft, mask=build_glob%l_resmsk)
                        call pftc%set_ref_pft(icls, pft, iseven=.true.)
                        call pftc%cp_even2odd_ref(icls)
                    endif
                endif
                deallocate(pft)
                call match_imgs(icls)%kill
            endif
        end do
        !$omp end parallel do
        call pftc%memoize_refs
        ! CLEANUP
        deallocate(match_imgs)
        call cavgs_m(1)%kill_thread_safe_tmp_imgs
        nullify(cavgs_m,cavgs_e,cavgs_o)
    end subroutine preppftc4align2D

    !>  \brief  prepares the polarft corrcalc object for search and imports polar references
    subroutine prep_polar_pftc4align2D( pftc, batchsz_max, which_iter, l_stream )
        use simple_strategy2D3D_common, only: calc_2Dref_offset
        class(polarft_calc), intent(inout) :: pftc
        integer,             intent(in)    :: batchsz_max, which_iter
        logical,             intent(in)    :: l_stream
        type(image),      allocatable :: tmp_imgs(:)
        type(string) :: fname
        real         :: xyz(3)
        integer      :: icls, pop, pop_even, pop_odd, centype
        logical      :: has_been_searched, do_center, l_center, l_gaufilt
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
        call pftc%polar_cavger_new(trim(params_glob%ref_type)=='comlin_hybrid')
        call pftc%polar_cavger_read_all(string(POLAR_REFS_FBODY)//BIN_EXT)
        has_been_searched = .not.build_glob%spproj%is_virgin_field(params_glob%oritype)
        ! Centering-related objects
        do_center = (params_glob%center .eq. 'yes') .and. has_been_searched&
             &.and. (which_iter > 2) .and. (.not.params_glob%l_update_frac)
        if( do_center )then
            allocate(tmp_imgs(params_glob%ncls))
            call pftc%polar_cavger_refs2cartesian(tmp_imgs, 'merged')
            call tmp_imgs(1)%construct_thread_safe_tmp_imgs(nthr_glob)
            ! mask memoization
            call tmp_imgs(1)%memoize_mask_coords
        endif
        ! Filtering
        l_gaufilt = trim(params_glob%gauref)=='yes'
        ! Mode of cavg centering
        centype = get_centype(params_glob%center_type)
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
                xyz      = 0.
                if( l_center ) call calc_2Dref_offset(tmp_imgs(icls), icls, centype, xyz)
                ! Prep for alignment
                call pftc%polar_prep2Dref(build_glob%clsfrcs, icls, l_gaufilt)
                ! transfer to pftc
                if( params_glob%l_lpset )then
                    ! merged class average in both even and odd positions
                    call pftc%polar_cavger_set_ref_pft(icls, 'merged')
                    call pftc%cp_even2odd_ref(icls)
                else
                    if( pop_even >= MINCLSPOPLIM .and. pop_odd >= MINCLSPOPLIM )then
                        ! transfer e/o refs to pftc
                        call pftc%polar_cavger_set_ref_pft(icls, 'even')
                        call pftc%polar_cavger_set_ref_pft(icls, 'odd')
                    else
                        ! merged class average in both even and odd positions
                        call pftc%polar_cavger_set_ref_pft(icls, 'merged')
                        call pftc%cp_even2odd_ref(icls)
                    endif
                endif
                ! centering cavg & particles within the pftc
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
