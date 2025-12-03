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
use simple_eul_prob_tab,            only: eul_prob_tab
use simple_polarft_calc,        only: polarft_calc
use simple_strategy3D_shc,          only: strategy3D_shc
use simple_strategy3D_shc_smpl,     only: strategy3D_shc_smpl
use simple_strategy3D_snhc_smpl,    only: strategy3D_snhc_smpl
use simple_strategy3D_greedy,       only: strategy3D_greedy
use simple_strategy3D_greedy_smpl,  only: strategy3D_greedy_smpl
use simple_strategy3D_greedy_sub,   only: strategy3D_greedy_sub
use simple_strategy3D_prob,         only: strategy3D_prob
use simple_strategy3D_eval,         only: strategy3D_eval
use simple_strategy3D,              only: strategy3D
use simple_strategy3D_srch,         only: strategy3D_spec
use simple_convergence,             only: convergence
use simple_euclid_sigma2,           only: euclid_sigma2
use simple_strategy2D3D_common
use simple_polarops
implicit none

public :: refine3D_exec
private
#include "simple_local_flags.inc"

logical                    :: has_been_searched
type(eul_prob_tab), target :: eulprob_obj_part
type(image),   allocatable :: ptcl_match_imgs(:)
integer,       allocatable :: pinds(:)
type(string)               :: fname
integer                    :: nptcls2update
type(euclid_sigma2)        :: eucl_sigma
type(polarft_calc)     :: pftc
! benchmarking
integer(timer_int_kind)    :: t_init, t_build_batch_particles, t_prep_orisrch, t_align, t_rec, t_tot, t_projio
integer(timer_int_kind)    :: t_prepare_polar_references
real(timer_int_kind)       :: rt_init, rt_build_batch_particles, rt_prep_orisrch, rt_align, rt_rec, rt_tot, rt_projio
real(timer_int_kind)       :: rt_prepare_polar_references
type(string)               :: benchfname

contains

    subroutine refine3D_exec( cline, which_iter, converged )
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: converged
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
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif
        l_polar = trim(params_glob%polar).eq.'yes'
        select case(trim(params_glob%refine))
            case('eval','sigma')
                l_restore = .false.
            case DEFAULT
                l_restore = .true.
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

        ! PARTICLE INDEX SAMPLING FOR FRACTIONAL UPDATE (OR NOT)
        if( allocated(pinds) ) deallocate(pinds)
        if( str_has_substr(params_glob%refine, 'prob') )then
            ! generation of random sample and incr of updatecnts delegated to prob_align
            call build_glob%spproj_field%sample4update_reprod([params_glob%fromp,params_glob%top],&
            &nptcls2update, pinds )
        else
            ! sampled incremented
            if( params_glob%l_fillin .and. mod(which_iter,5) == 0 )then
                call sample_ptcls4fillin([params_glob%fromp,params_glob%top], .true., nptcls2update, pinds)
            else
                call sample_ptcls4update([params_glob%fromp,params_glob%top], .true., nptcls2update, pinds)
            endif
        endif

        ! PREP BATCH ALIGNMENT
        batchsz_max = min(nptcls2update,params_glob%nthr*BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls2update)/real(batchsz_max))
        batches     = split_nobjs_even(nptcls2update, nbatches)
        batchsz_max = maxval(batches(:,2)-batches(:,1)+1)

        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, POLARIZER, PTCLS
        if( L_BENCH_GLOB )then
            rt_init                    = toc(t_init)
            t_prepare_polar_references = tic()
        endif
        call prepare_refs_sigmas_ptcls( pftc, cline, eucl_sigma, ptcl_match_imgs, batchsz_max, which_iter,&
                                        &do_polar=(l_polar .and. .not.cline%defined('vol1')) )
        if( l_polar .and. l_restore )then
            ! for restoration
            if( cline%defined('vol1') )then
                call polar_cavger_new(pftc, .true.)
                if( params_glob%l_trail_rec )then
                    ! In the first iteration the polarized cartesian references are written down
                    call polar_cavger_writeall_pftcrefs(string(POLAR_REFS_FBODY))
                endif
            endif
            call polar_cavger_zero_pft_refs
            if( file_exists(params_glob%frcs) )then
                call build_glob%clsfrcs%read(params_glob%frcs)
            else
                call build_glob%clsfrcs%new(params_glob%nspace, params_glob%box_crop,&
                    &params_glob%smpd_crop, params_glob%nstates)
            endif
        endif
        if( L_BENCH_GLOB )then
            rt_prepare_polar_references = toc(t_prepare_polar_references)
            t_prep_orisrch              = tic()
        endif

        ! PREPARE STRATEGY3D
        call prep_strategy3D ! allocate s3D singleton
        allocate(strategy3Dspecs(batchsz_max),strategy3Dsrch(batchsz_max))

        ! READING THE ASSIGNMENT FOR PROB MODE
        if( str_has_substr(params_glob%refine, 'prob') .and. .not.(trim(params_glob%refine) .eq. 'sigma') )then
            call eulprob_obj_part%new(pinds)
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
        allocate(cnt_greedy(params_glob%nthr), cnt_all(params_glob%nthr), source=0)
        allocate(incr_shifts(2,batchsz_max),source=0.)
        do ibatch=1,nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            ! Prep particles in pftc
            if( L_BENCH_GLOB ) t_build_batch_particles = tic()
            call build_batch_particles(pftc, batchsz, pinds(batch_start:batch_end), ptcl_match_imgs)
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
                select case(trim(params_glob%refine))
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
                        if( build_glob%spproj_field%is_first_update(which_iter, iptcl) )then
                            allocate(strategy3D_greedy_smpl      :: strategy3Dsrch(iptcl_batch)%ptr)
                            cnt_greedy(ithr) = cnt_greedy(ithr) + 1
                        else
                            allocate(strategy3D_shc_smpl         :: strategy3Dsrch(iptcl_batch)%ptr)
                        endif
                    case('snhc_smpl')
                        if( build_glob%spproj_field%is_first_update(which_iter, iptcl) )then
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
                        call build_glob%spproj_field%get_ori(iptcl, orientation)
                        call build_glob%spproj_field%set(iptcl, 'proj', build_glob%eulspace%find_closest_proj(orientation))
                    case DEFAULT
                        THROW_HARD('refinement mode: '//trim(params_glob%refine)//' unsupported')
                end select
                strategy3Dspecs(iptcl_batch)%iptcl     = iptcl
                strategy3Dspecs(iptcl_batch)%iptcl_map = iptcl_map
                if( str_has_substr(params_glob%refine, 'prob') ) strategy3Dspecs(iptcl_batch)%eulprob_obj_part => eulprob_obj_part
                ! search
                if( associated(strategy3Dsrch(iptcl_batch)%ptr) )then
                    ! instance & search
                    call strategy3Dsrch(iptcl_batch)%ptr%new(strategy3Dspecs(iptcl_batch))
                    call strategy3Dsrch(iptcl_batch)%ptr%srch(ithr)
                    ! keep track of incremental shift
                    incr_shifts(:,iptcl_batch) = build_glob%spproj_field%get_2Dshift(iptcl) - strategy3Dsrch(iptcl_batch)%ptr%s%prev_shvec
                    ! cleanup
                    call strategy3Dsrch(iptcl_batch)%ptr%kill
                endif
                ! calculate sigma2 for ML-based refinement
                if ( params_glob%l_needs_sigma ) then
                    call build_glob%spproj_field%get_ori(iptcl, orientation)
                    call orientation%set_shift(incr_shifts(:,iptcl_batch))
                    call eucl_sigma%calc_sigma2(pftc, iptcl, orientation, 'proj')
                end if
            enddo ! Particles loop
            !$omp end parallel do
            if( L_BENCH_GLOB ) rt_align = rt_align + toc(t_align)
            ! restore polar cavgs
            if( l_polar .and. l_restore )then
                if( L_BENCH_GLOB ) t_rec = tic()
                call polar_cavger_update_sums(batchsz, pinds(batch_start:batch_end),&
                    &build_glob%spproj, pftc, incr_shifts(:,1:batchsz), is3D=.true.)
                if( L_BENCH_GLOB ) rt_rec = rt_rec + toc(t_rec)
            endif
        enddo
        ! report fraction of greedy searches
        frac_greedy = 0.
        if( any(cnt_greedy > 0) .and. any(cnt_all > 0) )then
            frac_greedy = real(sum(cnt_greedy)) / real(sum(cnt_all))
        endif
        call build_glob%spproj_field%set_all2single('frac_greedy', frac_greedy)
        ! cleanup
        do iptcl_batch = 1,batchsz_max
            nullify(strategy3Dsrch(iptcl_batch)%ptr)
        end do
        deallocate(strategy3Dsrch,strategy3Dspecs,batches)
        call eulprob_obj_part%kill
        call deallocate_class_samples(clssmp)

        ! WRITE SIGMAS FOR ML-BASED REFINEMENT
        if( params_glob%l_needs_sigma ) call eucl_sigma%write_sigma2

        ! CALCULATE PARTICLE WEIGHTS
        if( trim(params_glob%ptclw).eq.'yes' )then
            ! not supported
        else
            if( trim(params_glob%cavgw).eq.'yes' )then
                ! class averages
                call build_glob%spproj_field%calc_cavg_soft_weights(params_glob%frac)
            else
                ! particles
                call build_glob%spproj_field%calc_hard_weights(params_glob%frac)
            endif
        endif

        ! CLEAN
        call clean_strategy3D ! deallocate s3D singleton
        call build_glob%vol%kill
        call orientation%kill
        do ithr = 1,params_glob%nthr
            call ptcl_match_imgs(ithr)%kill
        enddo
        deallocate(ptcl_match_imgs)

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
        if( l_restore )then
            if( L_BENCH_GLOB ) t_rec = tic()
            if( l_polar )then
                ! Polar representation
                call killimgbatch
                call eucl_sigma%kill
                if( l_distr_exec_glob )then
                    call polar_cavger_readwrite_partial_sums('write')
                else
                    call polar_restoration
                endif
                call pftc%kill
            else
                ! Cartesian volume
                call pftc%kill
                if( trim(params_glob%volrec).eq.'yes' )then
                    if( trim(params_glob%projrec).eq.'yes' )then
                        call calc_projdir3Drec( cline, nptcls2update, pinds )
                    else
                        call calc_3Drec( cline, nptcls2update, pinds )
                    endif
                endif
                call eucl_sigma%kill
                call killimgbatch
            endif
            if( L_BENCH_GLOB ) rt_rec = rt_rec + toc(t_rec)
        else
            ! cleanup
            call killimgbatch
            call eucl_sigma%kill
            call pftc%kill
        endif

        ! REPORT CONVERGENCE
        call qsys_job_finished(string('simple_strategy3D_matcher :: refine3D_exec'))
        if( .not. params_glob%l_distr_exec .and. trim(params_glob%refine).ne.'sigma' )then
            converged = conv%check_conv3D(cline, params_glob%msk)
        endif
        if( L_BENCH_GLOB )then
            rt_tot  = toc(t_tot)
            doprint = .true.
            if( params_glob%part /= 1 ) doprint = .false.
            if( doprint )then
                benchfname = 'REFINE3D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
                call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation           : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'build_batch_particles    : ', rt_build_batch_particles
                write(fnr,'(a,1x,f9.2)') 'prepare_polar_references : ', rt_prepare_polar_references
                write(fnr,'(a,1x,f9.2)') 'orisrch3D preparation    : ', rt_prep_orisrch
                write(fnr,'(a,1x,f9.2)') '3D alignment             : ', rt_align
                write(fnr,'(a,1x,f9.2)') 'project file I/O         : ', rt_projio
                write(fnr,'(a,1x,f9.2)') 'reconstruction           : ', rt_rec
                write(fnr,'(a,1x,f9.2)') 'total time               : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation           : ', (rt_init/rt_tot)                     * 100.
                write(fnr,'(a,1x,f9.2)') 'build_batch_particles    : ', (rt_build_batch_particles/rt_tot)    * 100.
                write(fnr,'(a,1x,f9.2)') 'prepare_polar_references : ', (rt_prepare_polar_references/rt_tot) * 100.
                write(fnr,'(a,1x,f9.2)') 'orisrch3D preparation    : ', (rt_prep_orisrch/rt_tot)             * 100.
                write(fnr,'(a,1x,f9.2)') '3D alignment             : ', (rt_align/rt_tot)                    * 100.
                write(fnr,'(a,1x,f9.2)') 'project file I/O         : ', (rt_projio/rt_tot)                   * 100.
                write(fnr,'(a,1x,f9.2)') 'reconstruction           : ', (rt_rec/rt_tot)                      * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for          : ',&
                    &((rt_init+rt_build_batch_particles+rt_prepare_polar_references+rt_prep_orisrch+rt_align+rt_projio+rt_rec)/rt_tot) * 100.
                call fclose(fnr)
            endif
        endif

      contains

        subroutine polar_restoration()
            params_glob%refs = CAVGS_ITER_FBODY//int2str_pad(params_glob%which_iter,3)//params_glob%ext%to_char()
            call polar_cavger_merge_eos_and_norm(reforis=build_glob%eulspace)
            call polar_cavger_calc_and_write_frcs_and_eoavg(string(FRCS_FILE), cline)
            call polar_cavger_writeall(string(POLAR_REFS_FBODY))
            call polar_cavger_write_cartrefs(pftc, get_fbody(params_glob%refs,params_glob%ext,separator=.false.), 'merged')
            call polar_cavger_kill
        end subroutine polar_restoration

    end subroutine refine3D_exec
    
end module simple_strategy3D_matcher
