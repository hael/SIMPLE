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
use simple_polarft_corrcalc,        only: polarft_corrcalc
use simple_strategy3D_shc,          only: strategy3D_shc
use simple_strategy3D_shc_smpl,     only: strategy3D_shc_smpl
use simple_strategy3D_smpl,         only: strategy3D_smpl
use simple_strategy3D_smpl_sub,     only: strategy3D_smpl_sub
use simple_strategy3D_greedy,       only: strategy3D_greedy
use simple_strategy3D_greedy_smpl,  only: strategy3D_greedy_smpl
use simple_strategy3D_greedy_sub,   only: strategy3D_greedy_sub
use simple_strategy3D_prob,         only: strategy3D_prob
use simple_strategy3D,              only: strategy3D
use simple_strategy3D_srch,         only: strategy3D_spec
use simple_convergence,             only: convergence
use simple_euclid_sigma2,           only: euclid_sigma2
use simple_strategy2D3D_common
implicit none

public :: refine3D_exec, prepare_polar_references, pftcc
private
#include "simple_local_flags.inc"

logical, parameter             :: DEBUG_HERE = .false.
logical                        :: has_been_searched
type(polarft_corrcalc), target :: pftcc
type(eul_prob_tab),     target :: eulprob_obj_part
type(image),       allocatable :: ptcl_match_imgs(:)
integer,           allocatable :: prev_states(:), pinds(:)
logical,           allocatable :: ptcl_mask(:)
integer                        :: nptcls2update
type(euclid_sigma2)            :: eucl_sigma
! benchmarking
integer(timer_int_kind)        :: t_init,   t_build_batch_particles,  t_prep_orisrch,  t_align,  t_rec,  t_tot,  t_projio
integer(timer_int_kind)        :: t_prepare_polar_references, t_read_and_filter_refvols
real(timer_int_kind)           :: rt_init, rt_build_batch_particles, rt_prep_orisrch, rt_align, rt_rec, rt_tot, rt_projio
real(timer_int_kind)           :: rt_prepare_polar_references, rt_read_and_filter_refvols
character(len=STDLEN)          :: benchfname

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
        type(strategy3D_spec), allocatable :: strategy3Dspecs(:)
        real,                  allocatable :: resarr(:)
        integer,               allocatable :: batches(:,:)
        type(convergence) :: conv
        type(oris)        :: prev_oris
        type(ori)         :: orientation
        real    :: frac_srch_space, extr_thresh, extr_score_thresh, anneal_ratio
        integer :: nbatches, batchsz_max, batch_start, batch_end, batchsz
        integer :: iptcl, fnr, ithr, iptcl_batch, iptcl_map
        integer :: ibatch, iextr_lim, lpind_anneal, lpind_start
        logical :: doprint, do_extr
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
        if( params_glob%l_batchfrac )then
            if( trim(params_glob%refine).eq.'prob' )then
                call build_glob%spproj_field%sample4batchupdate_reprod([params_glob%fromp,params_glob%top],&
                &nptcls2update, pinds, ptcl_mask )
            else
                call build_glob%spproj_field%sample4batchupdate([params_glob%fromp,params_glob%top],&
                &params_glob%batchfrac, nptcls2update, pinds, ptcl_mask)
                call build_glob%spproj_field%incr_updatecnt([params_glob%fromp,params_glob%top], ptcl_mask)
            endif
            call prev_oris%copy(build_glob%spproj_field, build_glob%spproj_field%is_particle())
        else
            if( trim(params_glob%refine).eq.'prob' )then
                ! generation of random sample and incr of updatecnts delegated to prob_align
                call build_glob%spproj_field%sample4update_reprod([params_glob%fromp,params_glob%top],&
                &nptcls2update, pinds, ptcl_mask )
            else
                if( params_glob%l_frac_update )then
                    if( params_glob%l_stoch_update )then
                        call build_glob%spproj_field%sample4update_rnd([params_glob%fromp,params_glob%top],&
                            &params_glob%update_frac, nptcls2update, pinds, ptcl_mask, .true.) ! sampled incremented
                    else
                        if( build_glob%spproj_field%has_been_sampled() )then ! we have a random subset
                            call build_glob%spproj_field%sample4update_reprod([params_glob%fromp,params_glob%top],&
                                                    &nptcls2update, pinds, ptcl_mask)
                        else                                                 ! we generate a random subset
                            call build_glob%spproj_field%sample4update_rnd2([params_glob%fromp,params_glob%top],&
                            &params_glob%update_frac, nptcls2update, pinds, ptcl_mask, .true.) ! sampled incremented
                        endif
                    endif
                else                                                     ! we sample all state > 0
                    call build_glob%spproj_field%sample4update_all([params_glob%fromp,params_glob%top],&
                                                &nptcls2update, pinds, ptcl_mask, .true.) ! sampled incremented
                endif
                ! increment update counter
                call build_glob%spproj_field%incr_updatecnt([params_glob%fromp,params_glob%top], ptcl_mask)
            endif
        endif

        ! PREP BATCH ALIGNEMENT
        batchsz_max = min(nptcls2update,params_glob%nthr*BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls2update)/real(batchsz_max))
        batches     = split_nobjs_even(nptcls2update, nbatches)
        batchsz_max = maxval(batches(:,2)-batches(:,1)+1)

        ! PREPARE THE POLARFT DATA STRUCTURES
        if( L_BENCH_GLOB )then
            rt_init                    = toc(t_init)
            t_prepare_polar_references = tic()
        endif
        call prepare_polar_references(cline, batchsz_max)
        if( L_BENCH_GLOB )then
            rt_prepare_polar_references = toc(t_prepare_polar_references)
            t_prep_orisrch              = tic()
        endif
        call build_glob%img_crop_polarizer%init_polarizer(pftcc, params_glob%alpha)
        call build_glob%vol%kill
        call build_glob%vol_odd%kill
        call build_glob%vol2%kill
        ! array allocation for strategy3D
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: array allocation for strategy3D'
        call prep_strategy3D ! allocate s3D singleton
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: array allocation for strategy3D, DONE'
        ! generate particles image objects
        allocate(strategy3Dspecs(batchsz_max),strategy3Dsrch(batchsz_max))
        call prepimgbatch(batchsz_max)
        allocate(ptcl_match_imgs(params_glob%nthr))
        do ithr = 1,params_glob%nthr
            call ptcl_match_imgs(ithr)%new([params_glob%box_crop, params_glob%box_crop, 1], params_glob%smpd_crop, wthreads=.false.)
        enddo
        write(logfhandle,'(A,1X,I3)') '>>> REFINE3D SEARCH, ITERATION:', which_iter
        if( str_has_substr(params_glob%refine, 'prob') .and. .not.(trim(params_glob%refine) .eq. 'sigma') )then
            call eulprob_obj_part%new(pinds)
            call eulprob_obj_part%read_assignment(trim(ASSIGNMENT_FBODY)//'.dat')
        endif
        if( L_BENCH_GLOB )then
            rt_prep_orisrch          = toc(t_prep_orisrch)
            rt_build_batch_particles = 0.
            rt_align                 = 0.
        endif
        ! Batch loop
        do ibatch=1,nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            ! Prep particles in pftcc
            if( L_BENCH_GLOB ) t_build_batch_particles = tic()
            call build_batch_particles(batchsz, pinds(batch_start:batch_end))
            if( L_BENCH_GLOB ) rt_build_batch_particles = rt_build_batch_particles + toc(t_build_batch_particles)
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
                    case('shc_smpl')
                        allocate(strategy3D_shc_smpl             :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('neigh')
                        allocate(strategy3D_greedy_sub           :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('greedy')
                        allocate(strategy3D_greedy               :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('prob')
                        allocate(strategy3D_prob                 :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('smpl')
                        allocate(strategy3D_smpl                 :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('smpl_neigh')
                        allocate(strategy3D_smpl_sub             :: strategy3Dsrch(iptcl_batch)%ptr)
                    case('sigma')
                        ! first sigma estimation (done below)
                        call build_glob%spproj_field%get_ori(iptcl, orientation)
                        call build_glob%spproj_field%set(iptcl, 'proj', real(build_glob%eulspace%find_closest_proj(orientation)))
                    case DEFAULT
                        THROW_HARD('refinement mode: '//trim(params_glob%refine)//' unsupported')
                end select
                strategy3Dspecs(iptcl_batch)%iptcl     = iptcl
                strategy3Dspecs(iptcl_batch)%iptcl_map = iptcl_map
                if( str_has_substr(params_glob%refine, 'prob') ) strategy3Dspecs(iptcl_batch)%eulprob_obj_part => eulprob_obj_part
                ! search object(s) & search
                if( associated(strategy3Dsrch(iptcl_batch)%ptr) )then
                    call strategy3Dsrch(iptcl_batch)%ptr%new(strategy3Dspecs(iptcl_batch))
                    call strategy3Dsrch(iptcl_batch)%ptr%srch(ithr)
                    call strategy3Dsrch(iptcl_batch)%ptr%kill
                endif
                ! calculate sigma2 for ML-based refinement
                if ( params_glob%l_needs_sigma ) then
                    call build_glob%spproj_field%get_ori(iptcl, orientation)
                    call eucl_sigma%calc_sigma2(pftcc, iptcl, orientation, 'proj')
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
        call eulprob_obj_part%kill

        ! WRITE SIGMAS FOR ML-BASED REFINEMENT
        if( params_glob%l_needs_sigma ) call eucl_sigma%write_sigma2

        ! CALCULATE PARTICLE WEIGHTS
        if( .not.(trim(params_glob%ptclw).eq.'yes') )then
            call build_glob%spproj_field%calc_hard_weights(params_glob%frac)
        endif

        ! CLEAN
        call clean_strategy3D ! deallocate s3D singleton
        call pftcc%kill
        if( str_has_substr(params_glob%refine, 'prob') ) call eulprob_obj_part%kill
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
        select case(trim(params_glob%refine))
            case('eval','sigma')
                ! no reconstruction
            case DEFAULT
                if( L_BENCH_GLOB ) t_rec = tic()
                if( params_glob%l_batchfrac )then
                    call calc_3Dbatchrec( cline, nptcls2update, pinds, prev_oris, which_iter )
                    call prev_oris%kill
                else
                    call calc_3Drec( cline, nptcls2update, pinds, which_iter )
                endif
                call eucl_sigma%kill
                call killimgbatch
                if( L_BENCH_GLOB ) rt_rec = toc(t_rec)
        end select

        ! REPORT CONVERGENCE
        call qsys_job_finished('simple_strategy3D_matcher :: refine3D_exec')
        if( .not. params_glob%l_distr_exec .and. trim(params_glob%refine).ne.'sigma' )then
            converged = conv%check_conv3D(cline, params_glob%msk)
        endif
        if( L_BENCH_GLOB )then
            rt_tot  = toc(t_tot)
            doprint = .true.
            if( params_glob%part /= 1 ) doprint = .false.
            if( doprint )then
                benchfname = 'REFINE3D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
                call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation           : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'build_batch_particles    : ', rt_build_batch_particles
                write(fnr,'(a,1x,f9.2)') 'prepare_polar_references : ', rt_prepare_polar_references
                write(fnr,'(a,1x,f9.2)') 'read_and_filter_refvols  : ', rt_read_and_filter_refvols
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
                write(fnr,'(a,1x,f9.2)') 'read_and_filter_refvols  : ', (rt_read_and_filter_refvols/rt_tot)  * 100.
                write(fnr,'(a,1x,f9.2)') 'orisrch3D preparation    : ', (rt_prep_orisrch/rt_tot)             * 100.
                write(fnr,'(a,1x,f9.2)') '3D alignment             : ', (rt_align/rt_tot)                    * 100.
                write(fnr,'(a,1x,f9.2)') 'project file I/O         : ', (rt_projio/rt_tot)                   * 100.
                write(fnr,'(a,1x,f9.2)') 'reconstruction           : ', (rt_rec/rt_tot)                      * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for          : ',&
                    &((rt_init+rt_build_batch_particles+rt_prepare_polar_references+rt_prep_orisrch+rt_align+rt_projio+rt_rec)/rt_tot) * 100.
                call fclose(fnr)
            endif
        endif
    end subroutine refine3D_exec

    subroutine prepare_polar_references( cline, batchsz_max )
        class(cmdline), intent(in)    :: cline !< command line
        integer,        intent(in)    :: batchsz_max
        character(len=:), allocatable :: fname
        type(ori) :: o_tmp
        real      :: xyz(3)
        integer   :: cnt, s, iref, nrefs
        logical   :: do_center
        if( L_BENCH_GLOB ) rt_read_and_filter_refvols = 0.
        ! first the polar
        nrefs = params_glob%nspace * params_glob%nstates
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
            if( trim(params_glob%refine).eq.'prob' )then
                ! already mapping shifts in prob_tab with shared-memory execution
                call calcrefvolshift_and_mapshifts2ptcls( cline, s, params_glob%vols(s), do_center, xyz, map_shift=l_distr_exec_glob)
            else
                call calcrefvolshift_and_mapshifts2ptcls( cline, s, params_glob%vols(s), do_center, xyz, map_shift=.true.)
            endif
            if( L_BENCH_GLOB ) t_read_and_filter_refvols = tic()
            call read_and_filter_refvols(s)
            if( L_BENCH_GLOB ) rt_read_and_filter_refvols = rt_read_and_filter_refvols + toc(t_read_and_filter_refvols)
            ! must be done here since params_glob%kfromto is dynamically set when lp_auto=yes
            if( s == 1 )then
                call pftcc%new(nrefs, [1,batchsz_max], params_glob%kfromto)
                if( params_glob%l_needs_sigma )then
                    fname = SIGMA2_FBODY//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'
                    call eucl_sigma%new(fname, params_glob%box)
                    call eucl_sigma%read_part(  build_glob%spproj_field, ptcl_mask)
                    call eucl_sigma%read_groups(build_glob%spproj_field, ptcl_mask)
                end if
            endif
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
            !$omp end parallel do
        end do
        call pftcc%memoize_refs
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: finished prepare_polar_references'
    end subroutine prepare_polar_references

    subroutine build_batch_particles( nptcls_here, pinds_here )
        use simple_strategy2D3D_common, only: read_imgbatch, prepimg4align
        integer, intent(in) :: nptcls_here
        integer, intent(in) :: pinds_here(nptcls_here)
        integer :: iptcl_batch, iptcl, ithr
        call discrete_read_imgbatch( nptcls_here, pinds_here, [1,nptcls_here], params_glob%l_use_denoised )
        ! reassign particles indices & associated variables
        call pftcc%reallocate_ptcls(nptcls_here, pinds_here)
        !$omp parallel do default(shared) private(iptcl,iptcl_batch,ithr) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds_here(iptcl_batch)
            ! prep
            call prepimg4align(iptcl, build_glob%imgbatch(iptcl_batch), ptcl_match_imgs(ithr))
            ! transfer to polar coordinates
            call build_glob%img_crop_polarizer%polarize(pftcc, ptcl_match_imgs(ithr), iptcl, .true., .true., mask=build_glob%l_resmsk)
            ! e/o flags
            call pftcc%set_eo(iptcl, nint(build_glob%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
        call pftcc%create_polar_absctfmats(build_glob%spproj, 'ptcl3D')
        ! Memoize particles FFT parameters
        call pftcc%memoize_ptcls
    end subroutine build_batch_particles
    
end module simple_strategy3D_matcher
