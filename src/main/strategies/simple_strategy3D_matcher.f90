!@descr: high-level search routines for the refine3D application
module simple_strategy3D_matcher
use simple_pftc_srch_api
use simple_strategy3D_alloc
use simple_strategy2D3D_common
use simple_binoris_io
use simple_builder,                 only: builder
use simple_convergence,             only: convergence
use simple_euclid_sigma2,           only: euclid_sigma2
use simple_eul_prob_tab,            only: eul_prob_tab
use simple_qsys_funs,               only: qsys_job_finished
use simple_strategy3D,              only: strategy3D
use simple_strategy3D_eval,         only: strategy3D_eval
use simple_strategy3D_greedy,       only: strategy3D_greedy
use simple_strategy3D_greedy_smpl,  only: strategy3D_greedy_smpl
use simple_strategy3D_greedy_sub,   only: strategy3D_greedy_sub
use simple_strategy3D_prob,         only: strategy3D_prob
use simple_strategy3D_shc,          only: strategy3D_shc
use simple_strategy3D_shc_smpl,     only: strategy3D_shc_smpl
use simple_strategy3D_snhc_smpl,    only: strategy3D_snhc_smpl
use simple_strategy3D_srch,         only: strategy3D_spec
implicit none

public :: refine3D_exec
private
#include "simple_local_flags.inc"

logical                    :: has_been_searched
type(eul_prob_tab), target :: eulprob_obj_part
type(image),   allocatable :: ptcl_match_imgs(:), ptcl_match_imgs_pad(:)
integer,       allocatable :: pinds(:)
type(string)               :: fname
integer                    :: nptcls2update
type(polarft_calc)         :: pftc
class(parameters), pointer :: p_ptr => null()
class(builder),    pointer :: b_ptr => null()
! benchmarking
integer(timer_int_kind)    :: t_init, t_build_batch_particles, t_prep_orisrch, t_align, t_rec, t_tot, t_projio
integer(timer_int_kind)    :: t_prepare_refs_sigmas_ptcls, t_prepare_polar_references
real(timer_int_kind)       :: rt_init, rt_build_batch_particles, rt_prep_orisrch, rt_align, rt_rec, rt_tot, rt_projio
real(timer_int_kind)       :: rt_prepare_refs_sigmas_ptcls, rt_prepare_polar_references
type(string)               :: benchfname

contains

    subroutine refine3D_exec( params, build, cline, which_iter, converged )
        class(parameters), target, intent(inout) :: params
        class(builder),    target, intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
        integer,                   intent(in)    :: which_iter
        logical,                   intent(inout) :: converged
        !---> The below is to allow particle-dependent decision about which 3D strategy to use
        type :: strategy3D_per_ptcl
            class(strategy3D), pointer :: ptr  => null()
        end type strategy3D_per_ptcl
        type(strategy3D_per_ptcl), allocatable :: strategy3Dsrch(:)
        !<---- hybrid or combined search strategies can then be implemented as extensions of the
        !      relevant strategy3D base class
        type(strategy3D_spec),     allocatable :: strategy3Dspecs(:)
        type(class_sample),        allocatable :: clssmp(:)
        integer,                   allocatable :: batches(:,:), cnt_greedy(:), cnt_all(:)
        real,                      allocatable :: incr_shifts(:,:)
        type(convergence) :: conv
        type(ori)         :: orientation
        real    :: frac_greedy
        integer :: nbatches, batchsz_max, batch_start, batch_end, batchsz
        integer :: iptcl, fnr, ithr, iptcl_batch, iptcl_map
        integer :: ibatch
        logical :: doprint, l_polar, l_restore
        
        ! assign parameters pointer
        p_ptr => params
        ! assign builder pointer
        b_ptr => build
        
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif
        l_polar = trim(p_ptr%polar).eq.'yes'
        select case(trim(p_ptr%refine))
            case('eval','sigma')
                l_restore = .false.
            case DEFAULT
                l_restore = .true.
        end select

        ! CHECK THAT WE HAVE AN EVEN/ODD PARTITIONING
        if( b_ptr%spproj_field%get_nevenodd() == 0 )then
            if( l_distr_exec_glob ) THROW_HARD('no eo partitioning available; refine3D_exec')
            call b_ptr%spproj_field%partition_eo
            call b_ptr%spproj%write_segment_inside(p_ptr%oritype)
        endif

        ! CHECK WHETHER WE HAVE PREVIOUS 3D ORIENTATIONS
        has_been_searched = .not.b_ptr%spproj%is_virgin_field(p_ptr%oritype)

        ! SET FOURIER INDEX RANGE
        call set_bp_range(p_ptr, b_ptr, cline)

        ! PARTICLE INDEX SAMPLING FOR FRACTIONAL UPDATE (OR NOT)
        if( allocated(pinds) ) deallocate(pinds)
        if( str_has_substr(p_ptr%refine, 'prob') )then
            ! generation of random sample and incr of updatecnts delegated to prob_align
            call b_ptr%spproj_field%sample4update_reprod([p_ptr%fromp,p_ptr%top],&
            &nptcls2update, pinds )
        else
            ! sampled incremented
            if( p_ptr%l_fillin .and. mod(which_iter,5) == 0 )then
                call sample_ptcls4fillin(b_ptr, [p_ptr%fromp,p_ptr%top], .true., nptcls2update, pinds)
            else
                call sample_ptcls4update(p_ptr, b_ptr, [p_ptr%fromp,p_ptr%top], .true., nptcls2update, pinds)
            endif
        endif

        ! PREP BATCH ALIGNMENT
        batchsz_max = min(nptcls2update,p_ptr%nthr*BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls2update)/real(batchsz_max))
        batches     = split_nobjs_even(nptcls2update, nbatches)
        batchsz_max = maxval(batches(:,2)-batches(:,1)+1)

        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, POLARIZER, PTCLS
        if( L_BENCH_GLOB )then
            rt_init = toc(t_init)
            t_prepare_refs_sigmas_ptcls = tic()
        endif
        call prepare_refs_sigmas_ptcls( p_ptr, b_ptr, pftc, cline, ptcl_match_imgs, ptcl_match_imgs_pad,&
                                        &batchsz_max, which_iter, do_polar=(l_polar .and. .not.cline%defined('vol1')) )
        if( L_BENCH_GLOB )then
            rt_prepare_refs_sigmas_ptcls = toc(t_prepare_refs_sigmas_ptcls)
            t_prepare_polar_references   = tic()
        endif
        if( l_polar .and. l_restore )then
            ! for restoration
            if( cline%defined('vol1') )then
                call pftc%polar_cavger(.true.)
                if( p_ptr%l_trail_rec )then
                    ! In the first iteration the polarized cartesian references are written down
                    call pftc%polar_cavger_writeall_pftcrefs(string(POLAR_REFS_FBODY))
                endif
            endif
            call pftc%polar_cavger_zero_pft_refs
            if( file_exists(p_ptr%frcs) )then
                call b_ptr%clsfrcs%read(p_ptr%frcs)
            else
                call b_ptr%clsfrcs%new(p_ptr%nspace, p_ptr%box_crop,&
                    &p_ptr%smpd_crop, p_ptr%nstates)
            endif
        endif
        if( L_BENCH_GLOB )then
            rt_prepare_polar_references = toc(t_prepare_polar_references)
            t_prep_orisrch              = tic()
        endif

        ! PREPARE STRATEGY3D
        call prep_strategy3D(p_ptr, b_ptr) ! allocate s3D singleton
        allocate(strategy3Dspecs(batchsz_max),strategy3Dsrch(batchsz_max))

        ! READING THE ASSIGNMENT FOR PROB MODE
        if( str_has_substr(p_ptr%refine, 'prob') .and. .not.(trim(p_ptr%refine) .eq. 'sigma') )then
            call eulprob_obj_part%new(p_ptr, b_ptr, pinds)
            call eulprob_obj_part%read_assignment(string(ASSIGNMENT_FBODY)//'.dat')
        endif

        if( L_BENCH_GLOB )then
            rt_prep_orisrch          = toc(t_prep_orisrch)
            rt_build_batch_particles = 0.
            rt_align                 = 0.
            rt_rec                   = 0.
        endif

        ! BATCH LOOP
        write(logfhandle,'(A,1X,I3)') '>>> REFINE3D SEARCH, ITERATION:', which_iter
        allocate(cnt_greedy(p_ptr%nthr), cnt_all(p_ptr%nthr), source=0)
        allocate(incr_shifts(2,batchsz_max),source=0.)
        do ibatch=1,nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            ! Prep particles in pftc
            if( L_BENCH_GLOB ) t_build_batch_particles = tic()
            call build_batch_particles(p_ptr, b_ptr, pftc, batchsz, pinds(batch_start:batch_end), ptcl_match_imgs, ptcl_match_imgs_pad)
            if( L_BENCH_GLOB ) rt_build_batch_particles = rt_build_batch_particles + toc(t_build_batch_particles)
            ! Particles loop
            if( L_BENCH_GLOB ) t_align = tic()
            !$omp parallel do default(shared) private(iptcl,iptcl_batch,iptcl_map,ithr,orientation)&
            !$omp schedule(static) proc_bind(close)
            do iptcl_batch    = 1,batchsz                     ! particle batch index
                iptcl_map     = batch_start + iptcl_batch - 1 ! masked global index (cumulative)
                iptcl         = pinds(iptcl_map)              ! global index
                ithr          = omp_get_thread_num() + 1
                cnt_all(ithr) = cnt_all(ithr) + 1
                ! switch for per-particle polymorphic strategy3D construction
                select case(trim(p_ptr%refine))
                    case('shc')
                        if( .not. has_been_searched )then
                            allocate(strategy3D_greedy           :: strategy3Dsrch(iptcl_batch)%ptr)
                            cnt_greedy(ithr) = cnt_greedy(ithr) + 1
                        else
                            if( ran3() < GREEDY_FREQ )then
                                allocate(strategy3D_greedy       :: strategy3Dsrch(iptcl_batch)%ptr)
                                cnt_greedy(ithr) = cnt_greedy(ithr) + 1
                            else
                                allocate(strategy3D_shc          :: strategy3Dsrch(iptcl_batch)%ptr)
                            endif
                        endif
                    case('shc_smpl')
                        if( b_ptr%spproj_field%is_first_update(which_iter, iptcl) )then
                            allocate(strategy3D_greedy_smpl      :: strategy3Dsrch(iptcl_batch)%ptr)
                            cnt_greedy(ithr) = cnt_greedy(ithr) + 1
                        else
                            allocate(strategy3D_shc_smpl         :: strategy3Dsrch(iptcl_batch)%ptr)
                        endif
                    case('snhc_smpl')
                        if( b_ptr%spproj_field%is_first_update(which_iter, iptcl) )then
                            allocate(strategy3D_greedy_smpl      :: strategy3Dsrch(iptcl_batch)%ptr)
                            cnt_greedy(ithr) = cnt_greedy(ithr) + 1
                        else
                            allocate(strategy3D_snhc_smpl        :: strategy3Dsrch(iptcl_batch)%ptr)
                        endif
                    case('eval')
                        allocate(strategy3D_eval                 :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('neigh')
                        allocate(strategy3D_greedy_sub           :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('greedy')
                        allocate(strategy3D_greedy               :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('prob','prob_state')
                        allocate(strategy3D_prob                 :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('sigma')
                        ! first sigma estimation (done below)
                        call b_ptr%spproj_field%get_ori(iptcl, orientation)
                        call b_ptr%spproj_field%set(iptcl, 'proj', b_ptr%eulspace%find_closest_proj(orientation))
                    case DEFAULT
                        THROW_HARD('refinement mode: '//trim(p_ptr%refine)//' unsupported')
                end select
                strategy3Dspecs(iptcl_batch)%iptcl     = iptcl
                strategy3Dspecs(iptcl_batch)%iptcl_map = iptcl_map
                if( str_has_substr(p_ptr%refine, 'prob') ) strategy3Dspecs(iptcl_batch)%eulprob_obj_part => eulprob_obj_part
                ! search
                if( associated(strategy3Dsrch(iptcl_batch)%ptr) )then
                    ! instance & search
                    call strategy3Dsrch(iptcl_batch)%ptr%new(p_ptr, strategy3Dspecs(iptcl_batch), b_ptr)
                    call strategy3Dsrch(iptcl_batch)%ptr%srch(b_ptr%spproj_field, ithr)
                    ! keep track of incremental shift
                    incr_shifts(:,iptcl_batch) = b_ptr%spproj_field%get_2Dshift(iptcl) - strategy3Dsrch(iptcl_batch)%ptr%s%prev_shvec
                    ! cleanup
                    call strategy3Dsrch(iptcl_batch)%ptr%kill
                endif
                ! calculate sigma2 for ML-based refinement
                if ( p_ptr%l_needs_sigma ) then
                    call b_ptr%spproj_field%get_ori(iptcl, orientation)
                    call orientation%set_shift(incr_shifts(:,iptcl_batch))
                    call b_ptr%esig%calc_sigma2(pftc, iptcl, orientation, 'proj')
                end if
            enddo ! Particles loop
            !$omp end parallel do
            if( L_BENCH_GLOB ) rt_align = rt_align + toc(t_align)
            ! restore polar cavgs
            if( l_polar .and. l_restore )then
                if( L_BENCH_GLOB ) t_rec = tic()
                call pftc%polar_cavger_update_sums(batchsz, pinds(batch_start:batch_end),&
                    &b_ptr%spproj, b_ptr%esig%sigma2_noise, incr_shifts(:,1:batchsz), is3D=.true.)
                if( L_BENCH_GLOB ) rt_rec = rt_rec + toc(t_rec)
            endif
        enddo
        ! report fraction of greedy searches
        frac_greedy = 0.
        if( any(cnt_greedy > 0) .and. any(cnt_all > 0) )then
            frac_greedy = real(sum(cnt_greedy)) / real(sum(cnt_all))
        endif
        call b_ptr%spproj_field%set_all2single('frac_greedy', frac_greedy)
        ! cleanup
        do iptcl_batch = 1,batchsz_max
            nullify(strategy3Dsrch(iptcl_batch)%ptr)
        end do
        deallocate(strategy3Dsrch,strategy3Dspecs,batches)
        call eulprob_obj_part%kill
        call deallocate_class_samples(clssmp)

        ! WRITE SIGMAS FOR ML-BASED REFINEMENT
        if( p_ptr%l_needs_sigma ) call b_ptr%esig%write_sigma2

        ! CALCULATE PARTICLE WEIGHTS
        if( trim(p_ptr%ptclw).eq.'yes' )then
            ! not supported
        else
            if( trim(p_ptr%cavgw).eq.'yes' )then
                ! class averages
                call b_ptr%spproj_field%calc_cavg_soft_weights(p_ptr%frac)
            else
                ! particles
                call b_ptr%spproj_field%calc_hard_weights(p_ptr%frac)
            endif
        endif

        ! CLEAN
        call clean_strategy3D ! deallocate s3D singleton
        call b_ptr%vol%kill
        call orientation%kill
        call dealloc_imgarr(ptcl_match_imgs)
        call dealloc_imgarr(ptcl_match_imgs_pad)

        ! OUTPUT ORIENTATIONS
        select case(trim(p_ptr%refine))
            case('sigma')
                ! nothing to do
            case DEFAULT
                if( L_BENCH_GLOB ) t_projio = tic()
                select case(trim(p_ptr%oritype))
                    case('ptcl3D')
                        call binwrite_oritab(p_ptr%outfile, b_ptr%spproj, &
                            &b_ptr%spproj_field, [p_ptr%fromp,p_ptr%top], isegment=PTCL3D_SEG)
                    case('cls3D')
                        call binwrite_oritab(p_ptr%outfile, b_ptr%spproj, &
                            &b_ptr%spproj_field, [p_ptr%fromp,p_ptr%top], isegment=CLS3D_SEG)
                    case DEFAULT
                        THROW_HARD('unsupported oritype: '//trim(p_ptr%oritype)//'; refine3D_exec')
                end select
                p_ptr%oritab = p_ptr%outfile
                if( L_BENCH_GLOB ) rt_projio = toc(t_projio)
        end select

        ! VOLUMETRIC 3D RECONSTRUCTION
        if( l_restore )then
            if( L_BENCH_GLOB ) t_rec = tic()
            if( l_polar )then
                ! Polar representation
                call killimgbatch(b_ptr)
                call b_ptr%esig%kill
                if( l_distr_exec_glob )then
                    call pftc%polar_cavger_readwrite_partial_sums('write')
                else
                    call polar_restoration
                endif
                call pftc%kill
            else
                ! Cartesian volume
                call pftc%kill
                if( trim(p_ptr%volrec).eq.'yes' )then
                     call calc_3Drec( p_ptr, b_ptr, cline, nptcls2update, pinds )
                endif
                call b_ptr%esig%kill
                call killimgbatch(b_ptr)
            endif
            if( L_BENCH_GLOB ) rt_rec = rt_rec + toc(t_rec)
        else
            ! cleanup
            call killimgbatch(b_ptr)
            call b_ptr%esig%kill
            call pftc%kill
        endif

        ! REPORT CONVERGENCE
        call qsys_job_finished(p_ptr, string('simple_strategy3D_matcher :: refine3D_exec'))
        if( .not. p_ptr%l_distr_exec .and. trim(p_ptr%refine).ne.'sigma' )then
            converged = conv%check_conv3D(p_ptr, cline, b_ptr%spproj_field, p_ptr%msk)
        endif
        if( L_BENCH_GLOB )then
            rt_tot  = toc(t_tot)
            doprint = .true.
            if( p_ptr%part /= 1 ) doprint = .false.
            if( doprint )then
                benchfname = 'REFINE3D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
                call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation            : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'build_batch_particles     : ', rt_build_batch_particles
                write(fnr,'(a,1x,f9.2)') 'prepare_refs_sigmas_ptcls : ', rt_prepare_refs_sigmas_ptcls
                write(fnr,'(a,1x,f9.2)') 'prepare_polar_references  : ', rt_prepare_polar_references
                write(fnr,'(a,1x,f9.2)') 'orisrch3D preparation     : ', rt_prep_orisrch
                write(fnr,'(a,1x,f9.2)') '3D alignment              : ', rt_align
                write(fnr,'(a,1x,f9.2)') 'project file I/O          : ', rt_projio
                write(fnr,'(a,1x,f9.2)') 'reconstruction            : ', rt_rec
                write(fnr,'(a,1x,f9.2)') 'total time                : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation            : ', (rt_init/rt_tot)                      * 100.
                write(fnr,'(a,1x,f9.2)') 'build_batch_particles     : ', (rt_build_batch_particles/rt_tot)     * 100.
                write(fnr,'(a,1x,f9.2)') 'prepare_refs_sigmas_ptcls : ', (rt_prepare_refs_sigmas_ptcls/rt_tot) * 100.
                write(fnr,'(a,1x,f9.2)') 'prepare_polar_references  : ', (rt_prepare_polar_references/rt_tot)  * 100.
                write(fnr,'(a,1x,f9.2)') 'orisrch3D preparation     : ', (rt_prep_orisrch/rt_tot)              * 100.
                write(fnr,'(a,1x,f9.2)') '3D alignment              : ', (rt_align/rt_tot)                     * 100.
                write(fnr,'(a,1x,f9.2)') 'project file I/O          : ', (rt_projio/rt_tot)                    * 100.
                write(fnr,'(a,1x,f9.2)') 'reconstruction            : ', (rt_rec/rt_tot)                       * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for           : ',&
                    &((rt_init+rt_build_batch_particles+rt_prepare_refs_sigmas_ptcls+rt_prepare_polar_references+rt_prep_orisrch+rt_align+rt_projio+rt_rec)/rt_tot) * 100.
                call fclose(fnr)
            endif
        endif

      contains

        subroutine polar_restoration()
            p_ptr%refs = CAVGS_ITER_FBODY//int2str_pad(p_ptr%which_iter,3)//MRC_EXT
            call pftc%polar_cavger_merge_eos_and_norm(reforis=b_ptr%eulspace, symop=b_ptr%pgrpsyms)
            call pftc%polar_cavger_calc_and_write_frcs_and_eoavg(b_ptr%clsfrcs, b_ptr%spproj_field%get_update_frac(), string(FRCS_FILE), cline)
            call pftc%polar_cavger_writeall(string(POLAR_REFS_FBODY))
            call pftc%polar_cavger_kill
        end subroutine polar_restoration

    end subroutine refine3D_exec
    
end module simple_strategy3D_matcher
