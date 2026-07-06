!@descr: SIMPLE-native linear 3D variability analysis from fixed particle poses
module simple_flex_eigenvol_strategy
use simple_core_module_api
use simple_builder,          only: builder
use simple_cmdline,          only: cmdline
use simple_image,            only: image
use simple_imgarr_utils,     only: dealloc_imgarr
use simple_matcher_3Drec,    only: init_rec
use simple_matcher_ptcl_io,  only: discrete_read_imgbatch, discrete_read_imgbatch_source, prepimgbatch, killimgbatch
use simple_memoize_ft_maps,  only: memoize_ft_maps, forget_ft_maps
use simple_parameters,       only: parameters
use simple_reconstructor,    only: reconstructor, insert_plane_oversamp_multi_scaled, project_fplanes_mean_basis
implicit none

public :: run_flex_eigenvol_linear
private
#include "simple_local_flags.inc"

real(dp), parameter :: LATENT_RIDGE = 1.0d-3
real(dp), parameter :: MODE_VAR_FLOOR = 1.0d-3
real,     parameter :: TRAJ_SIGMAS(7) = [-3., -2., -1., 0., 1., 2., 3.]

contains

    subroutine run_flex_eigenvol_linear( params, build, cline )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        type(reconstructor)              :: mean_rec
        type(reconstructor), allocatable :: basis_recs(:)
        type(fplane_type),  allocatable  :: fpls(:)
        type(string)                     :: sigma2_fname, zfname
        integer, allocatable             :: pinds(:)
        real(dp), allocatable            :: z(:,:), z_postvar(:,:), resid_energy(:), resid_mean_energy(:), mode_vars(:)
        integer                          :: nptcls, ncomp, niters, iter, q, kfromto_eff(2)
        integer(timer_int_kind)          :: t_total, t_step
        t_total = tic()
        call validate_inputs(params, cline, ncomp, niters)
        call select_particles(params, build, pinds, nptcls)
        allocate(z(nptcls,ncomp), z_postvar(nptcls,ncomp), resid_energy(nptcls), resid_mean_energy(nptcls), &
            &mode_vars(ncomp), basis_recs(ncomp))
        call initialize_latents(z, nptcls, ncomp)
        call orthonormalize_latents(z, nptcls, ncomp)
        z_postvar = 0.d0
        mode_vars = 1.d0
        kfromto_eff = flex_kfromto(params)
        write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL PARTICLES  : ', nptcls
        write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL COMPONENTS : ', ncomp
        write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL PPCA ITERS : ', niters
        write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL THREADS    : ', nthr_glob
        write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL BOX        : ', params%box_crop
        write(logfhandle,'(A,F8.3)') '>>> FLEX_EIGENVOL SMPD       : ', params%smpd_crop
        write(logfhandle,'(A,F8.3)') '>>> FLEX_EIGENVOL MASK PX    : ', params%msk_crop
        write(logfhandle,'(A,F8.3)') '>>> FLEX_EIGENVOL BASIS LP A : ', params%lp
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_EIGENVOL KFROMTO    : ', kfromto_eff(1), ' / ', kfromto_eff(2)
        write(logfhandle,'(A,A)')  '>>> FLEX_EIGENVOL PTCL SRC   : ', trim(params%ptcl_src)
        write(logfhandle,'(A,A)')  '>>> FLEX_EIGENVOL MEAN MAP  : ', params%vols(1)%to_char()
        call flush(logfhandle)
        t_step = tic()
        sigma2_fname = string('flex_eigenvol_sigma2.tmp')
        call build%esig%new(params, build%pftc, sigma2_fname, params%box_crop)
        call init_mean_reconstructor(params, build, mean_rec)
        do q = 1, ncomp
            call init_output_reconstructor(params, build, basis_recs(q))
        end do
        call log_seconds('>>> FLEX_EIGENVOL SETUP SECONDS', toc(t_step))
        do iter = 1, niters
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_EIGENVOL LINEAR ITERATION ', iter, ' / ', niters
            call flush(logfhandle)
            t_step = tic()
            call update_basis_from_latents(params, build, mean_rec, basis_recs, z, z_postvar, pinds, nptcls, ncomp, fpls)
            call log_iter_seconds('>>> FLEX_EIGENVOL ITER M-STEP SECONDS', iter, toc(t_step))
            t_step = tic()
            call infer_latents_from_basis(params, build, mean_rec, basis_recs, z, mode_vars, &
                &z_postvar, resid_energy, resid_mean_energy, pinds, nptcls, ncomp, fpls)
            call log_iter_seconds('>>> FLEX_EIGENVOL ITER E-STEP SECONDS', iter, toc(t_step))
            call log_residual_stats('>>> FLEX_EIGENVOL ITER MEAN-ONLY RESIDUAL ENERGY', resid_mean_energy, nptcls)
            call log_residual_stats('>>> FLEX_EIGENVOL ITER MODE RESIDUAL ENERGY', resid_energy, nptcls)
            call log_residual_reduction('>>> FLEX_EIGENVOL ITER RESIDUAL REDUCTION', resid_mean_energy, resid_energy, nptcls)
            call log_latent_sdevs('>>> FLEX_EIGENVOL ITER RAW LATENT SD', z, nptcls, ncomp)
            call log_mode_vars('>>> FLEX_EIGENVOL ITER PPCA MODE VAR', mode_vars, ncomp)
            call orthonormalize_latents(z, nptcls, ncomp)
            call log_latent_sdevs('>>> FLEX_EIGENVOL ITER NORM LATENT SD', z, nptcls, ncomp)
        end do
        write(logfhandle,'(A)') '>>> FLEX_EIGENVOL FINAL REFIT'
        call flush(logfhandle)
        t_step = tic()
        call update_basis_from_latents(params, build, mean_rec, basis_recs, z, z_postvar, pinds, nptcls, ncomp, fpls)
        call log_seconds('>>> FLEX_EIGENVOL FINAL M-STEP SECONDS', toc(t_step))
        t_step = tic()
        call infer_latents_from_basis(params, build, mean_rec, basis_recs, z, mode_vars, &
            &z_postvar, resid_energy, resid_mean_energy, pinds, nptcls, ncomp, fpls)
        call log_seconds('>>> FLEX_EIGENVOL FINAL E-STEP SECONDS', toc(t_step))
        call log_residual_stats('>>> FLEX_EIGENVOL FINAL MEAN-ONLY RESIDUAL ENERGY', resid_mean_energy, nptcls)
        call log_residual_stats('>>> FLEX_EIGENVOL FINAL MODE RESIDUAL ENERGY', resid_energy, nptcls)
        call log_residual_reduction('>>> FLEX_EIGENVOL FINAL RESIDUAL REDUCTION', resid_mean_energy, resid_energy, nptcls)
        call log_latent_sdevs('>>> FLEX_EIGENVOL RAW FINAL LATENT SD', z, nptcls, ncomp)
        call log_latent_means('>>> FLEX_EIGENVOL RAW FINAL LATENT MEAN', z, nptcls, ncomp)
        call log_mode_vars('>>> FLEX_EIGENVOL FINAL PPCA MODE VAR', mode_vars, ncomp)
        call log_latent_covariance_eigs('>>> FLEX_EIGENVOL FINAL LATENT COV EIGVAL', z, nptcls, ncomp)
        call orthonormalize_latents(z, nptcls, ncomp)
        call log_latent_sdevs('>>> FLEX_EIGENVOL OUTPUT LATENT SD', z, nptcls, ncomp)
        call log_latent_means('>>> FLEX_EIGENVOL OUTPUT LATENT MEAN', z, nptcls, ncomp)
        call log_basis_fourier_norms('>>> FLEX_EIGENVOL OUTPUT BASIS FOURIER NORM', basis_recs, ncomp)
        zfname = output_prefix(params)//'_zcoords.txt'
        t_step = tic()
        call write_basis_outputs(params, basis_recs, ncomp)
        call write_particle_coordinates(zfname, pinds, z, resid_energy, nptcls, ncomp)
        call write_trajectory_volumes(params, basis_recs, z, nptcls, ncomp)
        call log_seconds('>>> FLEX_EIGENVOL OUTPUT SECONDS', toc(t_step))
        write(logfhandle,'(A,A)') '>>> FLEX_EIGENVOL WROTE Z TABLE : ', zfname%to_char()
        call log_seconds('>>> FLEX_EIGENVOL TOTAL SECONDS', toc(t_total))
        call flush(logfhandle)
        call cleanup_planes(fpls)
        if( allocated(pinds) ) deallocate(pinds)
        if( allocated(z) ) deallocate(z)
        if( allocated(z_postvar) ) deallocate(z_postvar)
        if( allocated(resid_energy) ) deallocate(resid_energy)
        if( allocated(resid_mean_energy) ) deallocate(resid_mean_energy)
        if( allocated(mode_vars) ) deallocate(mode_vars)
        call build%esig%kill
        call mean_rec%dealloc_rho
        call mean_rec%kill
        if( allocated(basis_recs) )then
            do q = 1, size(basis_recs)
                call basis_recs(q)%dealloc_rho
                call basis_recs(q)%kill
            end do
            deallocate(basis_recs)
        endif
        call sigma2_fname%kill
        call zfname%kill
    end subroutine run_flex_eigenvol_linear

    subroutine validate_inputs( params, cline, ncomp, niters )
        class(parameters), intent(in)    :: params
        class(cmdline),    intent(inout) :: cline
        integer,           intent(out)   :: ncomp, niters
        if( trim(params%oritype) /= 'ptcl3D' )then
            THROW_HARD('flex_eigenvol currently requires oritype=ptcl3D')
        endif
        if( .not. cline%defined('vol1') )then
            THROW_HARD('flex_eigenvol requires vol1=<mean map> from abinitio3D/refine3D')
        endif
        if( params%nstates /= 1 )then
            THROW_HARD('flex_eigenvol currently supports a single ptcl3D state')
        endif
        ncomp  = max(1, params%neigs)
        niters = max(1, params%maxits)
        if( ncomp > 16 )then
            THROW_WARN('flex_eigenvol capping neigs at 16 for this first PPCA implementation')
            ncomp = 16
        endif
    end subroutine validate_inputs

    subroutine select_particles( params, build, pinds, nptcls )
        class(parameters), intent(in)       :: params
        class(builder),    intent(inout)    :: build
        integer, allocatable, intent(inout) :: pinds(:)
        integer, intent(out)                :: nptcls
        call build%spproj_field%sample4rec([params%fromp, params%top], nptcls, pinds)
        if( nptcls <= 0 ) THROW_HARD('flex_eigenvol found no active particles')
    end subroutine select_particles

    subroutine init_mean_reconstructor( params, build, mean_rec )
        class(parameters),   intent(inout) :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        call mean_rec%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
        call mean_rec%fft
        call mean_rec%alloc_rho(params, build%spproj, expand=.true.)
        call mean_rec%expand_exp
    end subroutine init_mean_reconstructor

    subroutine init_output_reconstructor( params, build, basis_rec )
        class(parameters),   intent(inout) :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: basis_rec
        call basis_rec%new([params%box_crop, params%box_crop, params%box_crop], params%smpd_crop)
        call basis_rec%alloc_rho(params, build%spproj, expand=.true.)
        call basis_rec%reset
        call basis_rec%reset_exp
    end subroutine init_output_reconstructor

    subroutine update_basis_from_latents( params, build, mean_rec, basis_recs, z, z_postvar, pinds, nptcls, ncomp, fpls )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: nptcls, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(in)    :: z(nptcls,ncomp), z_postvar(nptcls,ncomp)
        integer,             intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        type(fplane_type) :: mean_fpl
        type(ori)         :: orientation
        integer           :: batchlims(2), batchsz, ibatch, i, iptcl, q, row, progress_stride
        integer(timer_int_kind) :: t_total, t_phase, t_comp
        t_total = tic()
        write(logfhandle,'(A)') '>>> FLEX_EIGENVOL M-STEP: UPDATING EIGENVOLUMES'
        call flush(logfhandle)
        progress_stride = max(1, 5 * MAXIMGBATCHSZ)
        do q = 1, ncomp
            call basis_recs(q)%reset
            call basis_recs(q)%reset_exp
        end do
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, init_volumes=.false.)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        t_phase = tic()
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call read_particles(params, build, nptcls, pinds, batchlims, batchsz)
            call prep_imgs4flex(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            do i = 1, batchsz
                row   = batchlims(1) + i - 1
                iptcl = pinds(row)
                call build%spproj_field%get_ori(iptcl, orientation)
                if( orientation%isstatezero() ) cycle
                call mean_rec%project_fplane(orientation, fpls(i), mean_fpl, apply_ctf_amp=.true.)
                call subtract_plane(fpls(i), mean_fpl)
                call insert_plane_oversamp_multi_scaled(basis_recs, build%pgrpsyms, orientation, fpls(i), &
                    &z(row,:), z(row,:) * z(row,:) + z_postvar(row,:))
            end do
            if( batchlims(2) == nptcls .or. mod(batchlims(2), progress_stride) == 0 )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_EIGENVOL M-STEP PARTICLES: ', batchlims(2), ' / ', nptcls
                call flush(logfhandle)
            endif
        end do
        call log_seconds('>>> FLEX_EIGENVOL M-STEP INSERT SECONDS', toc(t_phase))
        call orientation%kill
        call cleanup_runtime_batch(build, fpls)
        call cleanup_plane(mean_fpl)
        t_phase = tic()
        do q = 1, ncomp
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_EIGENVOL M-STEP FINALIZE COMPONENT ', q, ' / ', ncomp
            call flush(logfhandle)
            t_comp = tic()
            call finalize_basis_for_projection(params, basis_recs(q))
            call log_comp_seconds('>>> FLEX_EIGENVOL M-STEP FINALIZE SECONDS', q, toc(t_comp))
        end do
        call log_seconds('>>> FLEX_EIGENVOL M-STEP FINALIZE TOTAL SECONDS', toc(t_phase))
        call log_seconds('>>> FLEX_EIGENVOL M-STEP TOTAL SECONDS', toc(t_total))
    end subroutine update_basis_from_latents

    subroutine infer_latents_from_basis( params, build, mean_rec, basis_recs, z, mode_vars, &
        &z_postvar, resid_energy, resid_mean_energy, pinds, nptcls, ncomp, fpls )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: nptcls, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(inout) :: z(nptcls,ncomp)
        real(dp),            intent(inout) :: mode_vars(ncomp)
        real(dp),            intent(out)   :: z_postvar(nptcls,ncomp)
        real(dp),            intent(out)   :: resid_energy(nptcls), resid_mean_energy(nptcls)
        integer,             intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpls(:)
        type(ori),         allocatable :: orientations(:)
        real(dp), allocatable :: gram(:,:,:), rhs(:,:), zrow(:,:), post_diag(:,:), mode_second(:,:)
        integer           :: batchlims(2), batchsz, ibatch, i, iptcl, q, r, row, ithr, progress_stride
        integer(timer_int_kind) :: t_total, t_phase
        t_total = tic()
        write(logfhandle,'(A)') '>>> FLEX_EIGENVOL E-STEP: INFERRING LATENTS'
        call flush(logfhandle)
        progress_stride = max(1, 5 * MAXIMGBATCHSZ)
        allocate(basis_fpls(ncomp,nthr_glob), mean_fpls(nthr_glob), orientations(nthr_glob), &
            &gram(ncomp,ncomp,nthr_glob), rhs(ncomp,nthr_glob), zrow(ncomp,nthr_glob), &
            &post_diag(ncomp,nthr_glob), mode_second(ncomp,nthr_glob))
        resid_energy = 0.d0
        resid_mean_energy = 0.d0
        z_postvar = 0.d0
        mode_second = 0.d0
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, init_volumes=.false.)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        t_phase = tic()
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call read_particles(params, build, nptcls, pinds, batchlims, batchsz)
            call prep_imgs4flex(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            !$omp parallel do default(shared) private(i,row,iptcl,q,r,ithr) schedule(static) proc_bind(close)
            do i = 1, batchsz
                row   = batchlims(1) + i - 1
                iptcl = pinds(row)
                ithr  = omp_get_thread_num() + 1
                call build%spproj_field%get_ori(iptcl, orientations(ithr))
                if( orientations(ithr)%isstatezero() ) cycle
                call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(ithr), fpls(i), &
                    &mean_fpls(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
                call subtract_plane(fpls(i), mean_fpls(ithr))
                resid_mean_energy(row) = plane_energy(fpls(i))
                gram(:,:,ithr) = 0.d0
                rhs(:,ithr)    = 0.d0
                do q = 1, ncomp
                    rhs(q,ithr) = plane_inner_product(fpls(i), basis_fpls(q,ithr))
                    do r = q, ncomp
                        gram(q,r,ithr) = plane_inner_product(basis_fpls(q,ithr), basis_fpls(r,ithr))
                        gram(r,q,ithr) = gram(q,r,ithr)
                    end do
                    gram(q,q,ithr) = gram(q,q,ithr) + ppca_prior_precision(mode_vars(q))
                end do
                call solve_ppca_posterior(gram(:,:,ithr), rhs(:,ithr), zrow(:,ithr), post_diag(:,ithr))
                z(row,:) = zrow(:,ithr)
                do q = 1, ncomp
                    z_postvar(row,q) = post_diag(q,ithr)
                    mode_second(q,ithr) = mode_second(q,ithr) + zrow(q,ithr) * zrow(q,ithr) + post_diag(q,ithr)
                end do
                do q = 1, ncomp
                    call subtract_scaled_plane(fpls(i), basis_fpls(q,ithr), zrow(q,ithr))
                end do
                resid_energy(row) = plane_energy(fpls(i))
            end do
            !$omp end parallel do
            if( batchlims(2) == nptcls .or. mod(batchlims(2), progress_stride) == 0 )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_EIGENVOL E-STEP PARTICLES: ', batchlims(2), ' / ', nptcls
                call flush(logfhandle)
            endif
        end do
        call log_seconds('>>> FLEX_EIGENVOL E-STEP INFERENCE SECONDS', toc(t_phase))
        do q = 1, ncomp
            mode_vars(q) = max(MODE_VAR_FLOOR, sum(mode_second(q,:)) / real(max(1,nptcls), dp))
        end do
        do ithr = 1, nthr_glob
            call orientations(ithr)%kill
            call cleanup_plane(mean_fpls(ithr))
        end do
        call cleanup_runtime_batch(build, fpls)
        do ithr = 1, nthr_glob
            do q = 1, ncomp
                call cleanup_plane(basis_fpls(q,ithr))
            end do
        end do
        deallocate(basis_fpls, mean_fpls, orientations, gram, rhs, zrow, post_diag, mode_second)
        call log_seconds('>>> FLEX_EIGENVOL E-STEP TOTAL SECONDS', toc(t_total))
    end subroutine infer_latents_from_basis

    subroutine finalize_basis_for_projection( params, basis_rec )
        class(parameters),   intent(in)    :: params
        type(reconstructor), intent(inout) :: basis_rec
        call basis_rec%compress_exp
        call basis_rec%sampl_dens_correct
        call basis_rec%ifft
        call basis_rec%div(real(params%box))
        call regularize_basis_volume(params, basis_rec)
        call basis_rec%fft
        call basis_rec%expand_exp
    end subroutine finalize_basis_for_projection

    subroutine regularize_basis_volume( params, basis_rec )
        class(parameters),   intent(in)    :: params
        type(reconstructor), intent(inout) :: basis_rec
        if( basis_rec%is_ft() ) call basis_rec%ifft
        if( params%msk_crop > TINY )then
            call basis_rec%mask3D_soft(params%msk_crop, backgr=0.)
        endif
        if( params%lp > 2.0 * params%smpd_crop + TINY )then
            call basis_rec%bp(0., params%lp)
        endif
        if( params%msk_crop > TINY )then
            call basis_rec%mask3D_soft(params%msk_crop, backgr=0.)
        endif
    end subroutine regularize_basis_volume

    subroutine write_basis_outputs( params, basis_recs, ncomp )
        class(parameters),   intent(in)    :: params
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        type(string) :: fname
        integer      :: q
        do q = 1, ncomp
            fname = basis_volume_name(params, q)
            write(logfhandle,'(A,I0,A,I0,A,A)') '>>> FLEX_EIGENVOL WRITING COMPONENT ', q, ' / ', ncomp, ': ', fname%to_char()
            call flush(logfhandle)
            call write_projection_ready_volume(basis_recs(q), fname)
            write(logfhandle,'(A,I0,A,A)') '>>> FLEX_EIGENVOL WROTE COMPONENT ', q, ': ', fname%to_char()
            call flush(logfhandle)
            call fname%kill
        end do
    end subroutine write_basis_outputs

    subroutine write_projection_ready_volume( basis_rec, fname )
        type(reconstructor), intent(inout) :: basis_rec
        class(string),       intent(in)    :: fname
        type(image) :: img
        call img%copy(basis_rec)
        if( img%is_ft() ) call img%ifft
        call log_image_stats('>>> FLEX_EIGENVOL COMPONENT VOLUME STATS', img)
        call img%write(fname, del_if_exists=.true.)
        call img%kill
    end subroutine write_projection_ready_volume

    subroutine write_particle_coordinates( zfname, pinds, z, resid_energy, nptcls, ncomp )
        class(string), intent(in) :: zfname
        integer,      intent(in) :: nptcls, ncomp
        integer,      intent(in) :: pinds(nptcls)
        real(dp),     intent(in) :: z(nptcls,ncomp), resid_energy(nptcls)
        integer :: funit, ierr, i, q
        call del_file(zfname)
        call fopen(funit, file=zfname, status='NEW', action='WRITE', iostat=ierr)
        call fileiochk('flex_eigenvol opening '//zfname%to_char(), ierr)
        write(funit,'(A)', advance='no') '# particle'
        do q = 1, ncomp
            write(funit,'(A,I0)', advance='no') ' z', q
        end do
        write(funit,'(A)') ' residual_energy'
        do i = 1, nptcls
            write(funit,'(I10)', advance='no') pinds(i)
            do q = 1, ncomp
                write(funit,'(1X,ES14.6)', advance='no') z(i,q)
            end do
            write(funit,'(1X,ES14.6)') resid_energy(i)
        end do
        call fclose(funit)
    end subroutine write_particle_coordinates

    subroutine write_trajectory_volumes( params, basis_recs, z, nptcls, ncomp )
        class(parameters),   intent(inout) :: params
        integer,             intent(in)    :: nptcls, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(in)    :: z(nptcls,ncomp)
        type(image)  :: mean_img, basis_img, traj_img, diff_img
        type(string) :: fname
        real         :: sig
        integer      :: q, isig
        call mean_img%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
        call log_image_stats('>>> FLEX_EIGENVOL MEAN VOLUME STATS', mean_img)
        do q = 1, ncomp
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_EIGENVOL WRITING TRAJECTORY COMPONENT ', q, ' / ', ncomp
            call flush(logfhandle)
            call basis_img%copy(basis_recs(q))
            if( basis_img%is_ft() ) call basis_img%ifft
            call log_image_stats('>>> FLEX_EIGENVOL TRAJECTORY BASIS STATS', basis_img)
            sig = latent_sdev(z(:,q), nptcls)
            if( sig <= TINY ) sig = 1.
            write(logfhandle,'(A,I0,A,ES12.4)') '>>> FLEX_EIGENVOL TRAJECTORY LATENT SD ', q, ': ', sig
            call flush(logfhandle)
            do isig = 1, size(TRAJ_SIGMAS)
                call traj_img%copy(mean_img)
                call traj_img%add(basis_img, TRAJ_SIGMAS(isig) * sig)
                call log_image_stats('>>> FLEX_EIGENVOL TRAJECTORY VOLUME STATS', traj_img)
                call diff_img%copy(traj_img)
                call diff_img%subtr(mean_img)
                call log_image_stats('>>> FLEX_EIGENVOL TRAJECTORY MINUS MEAN STATS', diff_img)
                fname = trajectory_diff_volume_name(params, q, TRAJ_SIGMAS(isig))
                call diff_img%write(fname, del_if_exists=.true.)
                call fname%kill
                call diff_img%kill
                fname = trajectory_volume_name(params, q, TRAJ_SIGMAS(isig))
                call traj_img%write(fname, del_if_exists=.true.)
                call traj_img%kill
                call fname%kill
            end do
            call basis_img%kill
            write(logfhandle,'(A,I0)') '>>> FLEX_EIGENVOL WROTE TRAJECTORY COMPONENT ', q
            call flush(logfhandle)
        end do
        call mean_img%kill
    end subroutine write_trajectory_volumes

    subroutine read_particles( params, build, nptcls, pinds, batchlims, batchsz )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer, intent(in) :: nptcls, batchlims(2), batchsz
        integer, intent(in) :: pinds(nptcls)
        if( trim(params%ptcl_src) == 'raw' )then
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
        else
            call discrete_read_imgbatch_source(params, build, trim(params%ptcl_src), &
                &nptcls, pinds, batchlims, build%imgbatch(:batchsz))
        endif
    end subroutine read_particles

    subroutine prep_imgs4flex( params, build, nptcls, ptcl_imgs, pinds, fplanes )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        class(image),      intent(inout) :: ptcl_imgs(nptcls)
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        type(ctfparams) :: ctfparms(nthr_glob)
        real    :: shift(2)
        integer :: iptcl, i, ithr, kfromto(2)
        call memoize_ft_maps([params%boxpd, params%boxpd, 1], params%smpd)
        kfromto = flex_kfromto(params)
        !$omp parallel do default(shared) private(i,ithr,iptcl,shift) schedule(static) proc_bind(close)
        do i = 1, nptcls
            ithr   = omp_get_thread_num() + 1
            iptcl  = pinds(i)
            call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft(build%lmsk, build%img_pad_heap(ithr))
            ctfparms(ithr) = build%spproj%get_ctfparams(params%oritype, iptcl)
            shift = build%spproj_field%get_2Dshift(iptcl)
            if( params%l_ml_reg .and. allocated(build%esig%sigma2_noise) )then
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr), &
                    &shift, fplanes(i), build%esig%sigma2_noise(kfromto(1):kfromto(2),iptcl), store_transfer=.true.)
            else
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr), &
                    &shift, fplanes(i), store_transfer=.true.)
            endif
            call cap_fplane_for_flex(fplanes(i), kfromto)
        end do
        !$omp end parallel do
    end subroutine prep_imgs4flex

    subroutine cap_fplane_for_flex( fpl, kfromto )
        type(fplane_type), intent(inout) :: fpl
        integer,           intent(in)    :: kfromto(2)
        integer :: nyq_eff
        nyq_eff = max(OSMPL_PAD_FAC, OSMPL_PAD_FAC * kfromto(2))
        if( fpl%nyq > 0 ) fpl%nyq = min(fpl%nyq, nyq_eff)
    end subroutine cap_fplane_for_flex

    function flex_kfromto( params ) result( kfromto )
        class(parameters), intent(in) :: params
        integer :: kfromto(2), kto_full
        real    :: dstep_crop
        kto_full = max(1, fdim(params%box_crop) - 1)
        kfromto(1) = 1
        kfromto(2) = kto_full
        if( params%lp > 2.0 * params%smpd_crop + TINY )then
            dstep_crop = real(max(1, params%box_crop - 1)) * params%smpd_crop
            kfromto(2) = max(1, min(kto_full, int(dstep_crop / params%lp)))
        endif
    end function flex_kfromto

    subroutine initialize_latents( z, nptcls, ncomp )
        integer,  intent(in) :: nptcls, ncomp
        real(dp), intent(out) :: z(nptcls,ncomp)
        integer :: i, q
        real(dp) :: phase
        do q = 1, ncomp
            do i = 1, nptcls
                phase  = DTWOPI * real(q * (i - 1), dp) / real(max(1,nptcls), dp)
                z(i,q) = sin(phase)
            end do
        end do
    end subroutine initialize_latents

    subroutine orthonormalize_latents( z, nptcls, ncomp )
        integer,  intent(in)    :: nptcls, ncomp
        real(dp), intent(inout) :: z(nptcls,ncomp)
        integer  :: q, r
        real(dp) :: avg, denom, coeff, sd
        do q = 1, ncomp
            avg = sum(z(:,q)) / real(max(1,nptcls), dp)
            z(:,q) = z(:,q) - avg
            do r = 1, q - 1
                denom = sum(z(:,r) * z(:,r))
                if( denom > DTINY )then
                    coeff = sum(z(:,q) * z(:,r)) / denom
                    z(:,q) = z(:,q) - coeff * z(:,r)
                endif
            end do
            sd = sqrt(sum(z(:,q) * z(:,q)) / real(max(1,nptcls - 1), dp))
            if( sd <= DTINY )then
                call reseed_latent_column(z(:,q), nptcls, q)
                sd = sqrt(sum(z(:,q) * z(:,q)) / real(max(1,nptcls - 1), dp))
            endif
            if( sd > DTINY ) z(:,q) = z(:,q) / sd
        end do
    end subroutine orthonormalize_latents

    subroutine reseed_latent_column( zcol, nptcls, q )
        integer,  intent(in) :: nptcls, q
        real(dp), intent(out) :: zcol(nptcls)
        integer :: i
        real(dp) :: phase, avg
        do i = 1, nptcls
            phase   = DTWOPI * real((q + 1) * (i - 1), dp) / real(max(1,nptcls), dp)
            zcol(i) = cos(phase)
        end do
        avg = sum(zcol) / real(max(1,nptcls), dp)
        zcol = zcol - avg
    end subroutine reseed_latent_column

    function latent_sdev( zcol, nptcls ) result( sd )
        integer,  intent(in) :: nptcls
        real(dp), intent(in) :: zcol(nptcls)
        real :: sd
        real(dp) :: avg
        avg = sum(zcol) / real(max(1,nptcls), dp)
        sd  = real(sqrt(sum((zcol - avg) * (zcol - avg)) / real(max(1,nptcls - 1), dp)))
    end function latent_sdev

    subroutine log_latent_covariance_eigs( label, z, nptcls, ncomp )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: nptcls, ncomp
        real(dp),         intent(in) :: z(nptcls,ncomp)
        real(dp), allocatable :: cov(:,:), rot(:,:), eigvals(:)
        integer :: nrot, q
        if( ncomp <= 0 ) return
        allocate(cov(ncomp,ncomp), rot(ncomp,ncomp), eigvals(ncomp), source=0.d0)
        call latent_covariance(z, nptcls, ncomp, cov)
        nrot = 0
        call jacobi(cov, ncomp, ncomp, eigvals, rot, nrot)
        call eigsrt(eigvals, rot, ncomp, ncomp)
        eigvals = max(eigvals, 0.d0)
        do q = 1, ncomp
            write(logfhandle,'(A,1X,I0,A,ES12.4)') trim(label), q, ': ', eigvals(q)
        end do
        call flush(logfhandle)
        deallocate(cov, rot, eigvals)
    end subroutine log_latent_covariance_eigs

    subroutine latent_covariance( z, nptcls, ncomp, cov )
        integer,  intent(in)  :: nptcls, ncomp
        real(dp), intent(in)  :: z(nptcls,ncomp)
        real(dp), intent(out) :: cov(ncomp,ncomp)
        real(dp) :: avg(ncomp), denom
        integer  :: q, r
        do q = 1, ncomp
            avg(q) = sum(z(:,q)) / real(max(1,nptcls), dp)
        end do
        denom = real(max(1,nptcls - 1), dp)
        cov = 0.d0
        do q = 1, ncomp
            do r = q, ncomp
                cov(q,r) = sum((z(:,q) - avg(q)) * (z(:,r) - avg(r))) / denom
                cov(r,q) = cov(q,r)
            end do
        end do
    end subroutine latent_covariance

    subroutine log_seconds( label, seconds )
        character(len=*),      intent(in) :: label
        real(timer_int_kind),  intent(in) :: seconds
        write(logfhandle,'(A,A,F10.3)') trim(label), ': ', seconds
        call flush(logfhandle)
    end subroutine log_seconds

    subroutine log_iter_seconds( label, iter, seconds )
        character(len=*),     intent(in) :: label
        integer,              intent(in) :: iter
        real(timer_int_kind), intent(in) :: seconds
        write(logfhandle,'(A,1X,I0,A,F10.3)') trim(label), iter, ': ', seconds
        call flush(logfhandle)
    end subroutine log_iter_seconds

    subroutine log_comp_seconds( label, icomp, seconds )
        character(len=*),     intent(in) :: label
        integer,              intent(in) :: icomp
        real(timer_int_kind), intent(in) :: seconds
        write(logfhandle,'(A,1X,I0,A,F10.3)') trim(label), icomp, ': ', seconds
        call flush(logfhandle)
    end subroutine log_comp_seconds

    subroutine log_residual_stats( label, resid_energy, nptcls )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: nptcls
        real(dp),         intent(in) :: resid_energy(nptcls)
        real(dp) :: avg, sd
        avg = sum(resid_energy) / real(max(1,nptcls), dp)
        sd  = sqrt(sum((resid_energy - avg) * (resid_energy - avg)) / real(max(1,nptcls - 1), dp))
        write(logfhandle,'(A,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') trim(label), &
            &' mean=', avg, ' sd=', sd, ' min=', minval(resid_energy), ' max=', maxval(resid_energy)
        call flush(logfhandle)
    end subroutine log_residual_stats

    subroutine log_residual_reduction( label, resid_mean_energy, resid_energy, nptcls )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: nptcls
        real(dp),         intent(in) :: resid_mean_energy(nptcls), resid_energy(nptcls)
        real(dp) :: mean_only, with_modes, reduction
        mean_only  = sum(resid_mean_energy) / real(max(1,nptcls), dp)
        with_modes = sum(resid_energy)      / real(max(1,nptcls), dp)
        reduction  = 0.d0
        if( mean_only > DTINY ) reduction = (mean_only - with_modes) / mean_only
        write(logfhandle,'(A,A,ES12.4,A,ES12.4,A,ES12.4)') trim(label), &
            &' mean_only=', mean_only, ' with_modes=', with_modes, ' frac=', reduction
        call flush(logfhandle)
    end subroutine log_residual_reduction

    subroutine log_image_stats( label, img )
        character(len=*), intent(in)    :: label
        type(image),      intent(inout) :: img
        real :: ave, sdev, maxv, minv
        logical :: err
        if( img%is_ft() ) call img%ifft
        call img%stats(ave, sdev, maxv, minv, errout=err)
        write(logfhandle,'(A,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') trim(label), &
            &' mean=', ave, ' sd=', sdev, ' min=', minv, ' max=', maxv, ' range=', maxv - minv
        call flush(logfhandle)
    end subroutine log_image_stats

    subroutine log_latent_sdevs( label, z, nptcls, ncomp )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: nptcls, ncomp
        real(dp),         intent(in) :: z(nptcls,ncomp)
        integer :: q
        do q = 1, ncomp
            write(logfhandle,'(A,1X,I0,A,ES12.4)') trim(label), q, ': ', latent_sdev(z(:,q), nptcls)
        end do
        call flush(logfhandle)
    end subroutine log_latent_sdevs

    subroutine log_latent_means( label, z, nptcls, ncomp )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: nptcls, ncomp
        real(dp),         intent(in) :: z(nptcls,ncomp)
        integer :: q
        do q = 1, ncomp
            write(logfhandle,'(A,1X,I0,A,ES12.4)') trim(label), q, ': ', &
                &sum(z(:,q)) / real(max(1,nptcls), dp)
        end do
        call flush(logfhandle)
    end subroutine log_latent_means

    subroutine log_mode_vars( label, mode_vars, ncomp )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: ncomp
        real(dp),         intent(in) :: mode_vars(ncomp)
        integer :: q
        do q = 1, ncomp
            write(logfhandle,'(A,1X,I0,A,ES12.4)') trim(label), q, ': ', mode_vars(q)
        end do
        call flush(logfhandle)
    end subroutine log_mode_vars

    subroutine log_basis_fourier_norms( label, basis_recs, ncomp )
        character(len=*),   intent(in) :: label
        integer,            intent(in) :: ncomp
        type(reconstructor), intent(in) :: basis_recs(ncomp)
        real(dp) :: energy
        integer  :: q
        do q = 1, ncomp
            energy = basis_fourier_energy(basis_recs(q))
            write(logfhandle,'(A,1X,I0,A,ES12.4,A,ES12.4)') trim(label), q, &
                &': norm=', sqrt(max(0.d0, energy)), ' energy=', energy
        end do
        call flush(logfhandle)
    end subroutine log_basis_fourier_norms

    function basis_fourier_energy( basis_rec ) result( val )
        type(reconstructor), intent(in) :: basis_rec
        real(dp) :: val
        integer :: h, k, m
        val = 0.d0
        if( .not. allocated(basis_rec%cmat_exp) ) return
        do m = lbound(basis_rec%cmat_exp,3), ubound(basis_rec%cmat_exp,3)
            do k = lbound(basis_rec%cmat_exp,2), ubound(basis_rec%cmat_exp,2)
                do h = lbound(basis_rec%cmat_exp,1), ubound(basis_rec%cmat_exp,1)
                    val = val + real(basis_rec%cmat_exp(h,k,m) * conjg(basis_rec%cmat_exp(h,k,m)), dp)
                end do
            end do
        end do
    end function basis_fourier_energy

    subroutine subtract_plane( data_fpl, model_fpl )
        type(fplane_type), intent(inout) :: data_fpl
        type(fplane_type), intent(in)    :: model_fpl
        if( .not. allocated(data_fpl%cmplx_plane) .or. .not. allocated(model_fpl%cmplx_plane) )then
            THROW_HARD('subtract_plane received unallocated Fourier plane')
        endif
        data_fpl%cmplx_plane = data_fpl%cmplx_plane - model_fpl%cmplx_plane
    end subroutine subtract_plane

    subroutine subtract_scaled_plane( data_fpl, model_fpl, scale )
        type(fplane_type), intent(inout) :: data_fpl
        type(fplane_type), intent(in)    :: model_fpl
        real(dp),          intent(in)    :: scale
        if( .not. allocated(data_fpl%cmplx_plane) .or. .not. allocated(model_fpl%cmplx_plane) )then
            THROW_HARD('subtract_scaled_plane received unallocated Fourier plane')
        endif
        data_fpl%cmplx_plane = data_fpl%cmplx_plane - real(scale) * model_fpl%cmplx_plane
    end subroutine subtract_scaled_plane

    subroutine scale_plane_for_latent_insert( fpl, data_scale, density_scale )
        type(fplane_type), intent(inout) :: fpl
        real(dp),          intent(in)    :: data_scale, density_scale
        real :: data_scale_sp, density_scale_sp
        data_scale_sp    = real(data_scale)
        density_scale_sp = real(max(0.d0, density_scale))
        if( allocated(fpl%cmplx_plane) ) fpl%cmplx_plane = data_scale_sp * fpl%cmplx_plane
        if( allocated(fpl%ctfsq_plane) ) fpl%ctfsq_plane = density_scale_sp * fpl%ctfsq_plane
        if( allocated(fpl%transfer_plane) ) fpl%transfer_plane = data_scale_sp * fpl%transfer_plane
    end subroutine scale_plane_for_latent_insert

    function plane_inner_product( lhs_fpl, rhs_fpl ) result( val )
        use simple_math, only: ceil_div, floor_div
        type(fplane_type), intent(in) :: lhs_fpl, rhs_fpl
        real(dp) :: val
        complex(dp) :: acc
        integer :: h, k, hmin, hmax, kmin, kmax, nyq_eff, sample_stride
        acc = DCMPLX_ZERO
        sample_stride = OSMPL_PAD_FAC
        nyq_eff = lhs_fpl%nyq
        if( rhs_fpl%nyq > 0 ) nyq_eff = min(nyq_eff, rhs_fpl%nyq)
        if( nyq_eff <= 0 ) nyq_eff = ubound(lhs_fpl%cmplx_plane,1)
        hmin = max(sample_stride * ceil_div(lbound(lhs_fpl%cmplx_plane,1), sample_stride), &
            &sample_stride * ceil_div(-nyq_eff, sample_stride))
        hmax = min(sample_stride * floor_div(ubound(lhs_fpl%cmplx_plane,1), sample_stride), &
            &sample_stride * floor_div(nyq_eff, sample_stride))
        kmin = max(sample_stride * ceil_div(lbound(lhs_fpl%cmplx_plane,2), sample_stride), &
            &sample_stride * ceil_div(-nyq_eff, sample_stride))
        kmax = min(0, sample_stride * floor_div(nyq_eff, sample_stride))
        do k = kmin, kmax, sample_stride
            do h = hmin, hmax, sample_stride
                if( nint(sqrt(real(h*h + k*k))) > nyq_eff ) cycle
                acc = acc + cmplx(lhs_fpl%cmplx_plane(h,k), kind=dp) * conjg(cmplx(rhs_fpl%cmplx_plane(h,k), kind=dp))
            end do
        end do
        val = real(acc, dp)
    end function plane_inner_product

    function plane_energy( fpl ) result( val )
        use simple_math, only: ceil_div, floor_div
        type(fplane_type), intent(in) :: fpl
        real(dp) :: val
        integer :: h, k, hmin, hmax, kmin, kmax, nyq_eff, sample_stride
        val = 0.d0
        sample_stride = OSMPL_PAD_FAC
        nyq_eff = fpl%nyq
        if( nyq_eff <= 0 ) nyq_eff = ubound(fpl%cmplx_plane,1)
        hmin = max(sample_stride * ceil_div(lbound(fpl%cmplx_plane,1), sample_stride), &
            &sample_stride * ceil_div(-nyq_eff, sample_stride))
        hmax = min(sample_stride * floor_div(ubound(fpl%cmplx_plane,1), sample_stride), &
            &sample_stride * floor_div(nyq_eff, sample_stride))
        kmin = max(sample_stride * ceil_div(lbound(fpl%cmplx_plane,2), sample_stride), &
            &sample_stride * ceil_div(-nyq_eff, sample_stride))
        kmax = min(0, sample_stride * floor_div(nyq_eff, sample_stride))
        do k = kmin, kmax, sample_stride
            do h = hmin, hmax, sample_stride
                if( nint(sqrt(real(h*h + k*k))) > nyq_eff ) cycle
                val = val + real(fpl%cmplx_plane(h,k) * conjg(fpl%cmplx_plane(h,k)), dp)
            end do
        end do
    end function plane_energy

    real(dp) function ppca_prior_precision( mode_var ) result( prec )
        real(dp), intent(in) :: mode_var
        prec = 1.d0 / max(MODE_VAR_FLOOR, mode_var)
        prec = max(prec, LATENT_RIDGE)
    end function ppca_prior_precision

    subroutine solve_ppca_posterior( gram, rhs, x, post_diag )
        real(dp), intent(in)  :: gram(:,:), rhs(:)
        real(dp), intent(out) :: x(:), post_diag(:)
        real(dp) :: unit(size(rhs)), sol(size(rhs))
        integer  :: q, n
        n = size(rhs)
        call solve_latent_system(gram, rhs, x)
        post_diag = 0.d0
        do q = 1, n
            unit = 0.d0
            unit(q) = 1.d0
            call solve_latent_system(gram, unit, sol)
            post_diag(q) = max(0.d0, sol(q))
        end do
    end subroutine solve_ppca_posterior

    subroutine solve_latent_system( gram, rhs, x )
        real(dp), intent(in)  :: gram(:,:), rhs(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: a(size(rhs),size(rhs)), b(size(rhs)), rowtmp(size(rhs))
        real(dp) :: pivot, factor, tmp
        integer  :: n, i, j, k, piv
        n = size(rhs)
        a = gram
        b = rhs
        x = 0.d0
        do k = 1, n - 1
            piv   = k
            pivot = abs(a(k,k))
            do i = k + 1, n
                if( abs(a(i,k)) > pivot )then
                    pivot = abs(a(i,k))
                    piv   = i
                endif
            end do
            if( pivot <= DTINY ) return
            if( piv /= k )then
                rowtmp(:) = a(k,:)
                a(k,:)    = a(piv,:)
                a(piv,:)  = rowtmp(:)
                tmp       = b(k)
                b(k)      = b(piv)
                b(piv)    = tmp
            endif
            do i = k + 1, n
                factor = a(i,k) / a(k,k)
                a(i,k:n) = a(i,k:n) - factor * a(k,k:n)
                b(i)     = b(i)     - factor * b(k)
            end do
        end do
        if( abs(a(n,n)) <= DTINY ) return
        do i = n, 1, -1
            tmp = b(i)
            do j = i + 1, n
                tmp = tmp - a(i,j) * x(j)
            end do
            if( abs(a(i,i)) <= DTINY )then
                x = 0.d0
                return
            endif
            x(i) = tmp / a(i,i)
        end do
    end subroutine solve_latent_system

    subroutine copy_plane( src, dst )
        type(fplane_type), intent(in)    :: src
        type(fplane_type), intent(inout) :: dst
        call cleanup_plane(dst)
        dst%frlims  = src%frlims
        dst%shconst = src%shconst
        dst%nyq     = src%nyq
        allocate(dst%cmplx_plane(lbound(src%cmplx_plane,1):ubound(src%cmplx_plane,1), &
            &lbound(src%cmplx_plane,2):ubound(src%cmplx_plane,2)))
        allocate(dst%ctfsq_plane(lbound(src%ctfsq_plane,1):ubound(src%ctfsq_plane,1), &
            &lbound(src%ctfsq_plane,2):ubound(src%ctfsq_plane,2)))
        dst%cmplx_plane = src%cmplx_plane
        dst%ctfsq_plane = src%ctfsq_plane
        if( allocated(src%transfer_plane) )then
            allocate(dst%transfer_plane(lbound(src%transfer_plane,1):ubound(src%transfer_plane,1), &
                &lbound(src%transfer_plane,2):ubound(src%transfer_plane,2)))
            dst%transfer_plane = src%transfer_plane
        endif
    end subroutine copy_plane

    subroutine cleanup_runtime_batch( build, fpls )
        class(builder), intent(inout) :: build
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        call cleanup_planes(fpls)
        call dealloc_imgarr(build%img_pad_heap)
        call killimgbatch(build)
        call forget_ft_maps
    end subroutine cleanup_runtime_batch

    subroutine cleanup_planes( fpls )
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        integer :: i
        if( allocated(fpls) )then
            do i = 1, size(fpls)
                call cleanup_plane(fpls(i))
            end do
            deallocate(fpls)
        endif
    end subroutine cleanup_planes

    subroutine cleanup_plane( fpl )
        type(fplane_type), intent(inout) :: fpl
        if( allocated(fpl%cmplx_plane) ) deallocate(fpl%cmplx_plane)
        if( allocated(fpl%ctfsq_plane) ) deallocate(fpl%ctfsq_plane)
        if( allocated(fpl%transfer_plane) ) deallocate(fpl%transfer_plane)
        fpl%frlims  = 0
        fpl%shconst = 0.
        fpl%nyq     = 0
    end subroutine cleanup_plane

    function output_prefix( params ) result( prefix )
        class(parameters), intent(in) :: params
        type(string) :: prefix
        type(string) :: ext
        character(len=:), allocatable :: stem
        ext    = fname2ext(params%outvol)
        prefix = get_fbody(params%outvol, ext)
        stem   = prefix%to_char()
        if( len_trim(stem) >= 4 )then
            if( stem(len_trim(stem)-3:len_trim(stem)) == '_001' ) stem = stem(:len_trim(stem)-4)
        endif
        prefix = string(trim(stem))
        call ext%kill
    end function output_prefix

    function basis_volume_name( params, icomp ) result( fname )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: icomp
        type(string) :: fname, prefix
        character(len=3) :: tag
        if( icomp == 1 )then
            fname = params%outvol
            return
        endif
        write(tag,'(I3.3)') icomp
        prefix = output_prefix(params)
        fname  = prefix//'_'//tag//'.mrc'
        call prefix%kill
    end function basis_volume_name

    function trajectory_volume_name( params, icomp, sigma ) result( fname )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: icomp
        real,              intent(in) :: sigma
        type(string) :: fname, prefix
        character(len=3) :: ctag
        character(len=2) :: stag
        write(ctag,'(I3.3)') icomp
        if( sigma < 0. )then
            write(stag,'(A,I1)') 'm', nint(abs(sigma))
        else if( sigma > 0. )then
            write(stag,'(A,I1)') 'p', nint(abs(sigma))
        else
            stag = 'z0'
        endif
        prefix = output_prefix(params)
        fname  = prefix//'_mode'//ctag//'_'//stag//'.mrc'
        call prefix%kill
    end function trajectory_volume_name

    function trajectory_diff_volume_name( params, icomp, sigma ) result( fname )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: icomp
        real,              intent(in) :: sigma
        type(string) :: fname, prefix
        character(len=3) :: ctag
        character(len=2) :: stag
        write(ctag,'(I3.3)') icomp
        if( sigma < 0. )then
            write(stag,'(A,I1)') 'm', nint(abs(sigma))
        else if( sigma > 0. )then
            write(stag,'(A,I1)') 'p', nint(abs(sigma))
        else
            stag = 'z0'
        endif
        prefix = output_prefix(params)
        fname  = prefix//'_mode'//ctag//'_'//stag//'_diff.mrc'
        call prefix%kill
    end function trajectory_diff_volume_name

end module simple_flex_eigenvol_strategy
