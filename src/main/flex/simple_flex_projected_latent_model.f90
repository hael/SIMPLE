!@descr: Projection-aware latent volume model kernels for flex_analysis
module simple_flex_projected_latent_model
use simple_core_module_api
use simple_builder,          only: builder
use simple_image,            only: image
use simple_imgarr_utils,     only: dealloc_imgarr
use simple_gridding,         only: prep3D_inv_instrfun4mul
use simple_linalg,           only: eigsrt, hermitian_solve, jacobi
use simple_matcher_3Drec,    only: init_rec, prep_imgs4rec, cleanup_rec_buffers
use simple_matcher_ptcl_io,  only: discrete_read_imgbatch, discrete_read_imgbatch_source, prepimgbatch, killimgbatch
use simple_memoize_ft_maps,  only: memoize_ft_maps, forget_ft_maps
use simple_parameters,       only: parameters
use simple_reconstructor,    only: reconstructor
use simple_flex_reconstructor_latent_ops, only: insert_planes_oversamp_coupled_batch_scaled, &
    &accumulate_planes_oversamp_coupled_stats_batch, project_fplane_mean, project_fplanes_mean_basis
use simple_map_reduce,       only: split_nobjs_even
implicit none

public :: update_basis_from_latents, infer_latents_from_basis
public :: write_mstep_stats_part_file, update_basis_from_mstep_stats_part_files
public :: write_estep_latent_part_file, reduce_estep_latent_part_files
public :: initialize_latents, orthonormalize_latents, latent_sdev, latent_covariance
public :: basis_fourier_energy, cleanup_planes, projected_model_kfromto, prep_imgs4projected_model
public :: solve_coupled_basis_exp
public :: canonicalize_projected_latent_basis
public :: test_projected_latent_mstep_stats_io, test_projected_latent_canonicalization
public :: test_projected_model_plane_preparation
private
#include "simple_local_flags.inc"

real(dp), parameter :: COUPLED_MSTEP_RIDGE_REL = 1.0d-8
! Keep the same Fourier-cell observability convention as
! simple_image_arith::div_cmat_at_1, which zeroes a reconstruction cell when
! its sampling/CTF density is at or below 1.e-6.  The coupled solve must not
! turn those unobserved cells into high-amplitude noise by dividing by DTINY.
real(dp), parameter :: COUPLED_DENSITY_FLOOR = 1.0d-6
real(dp), parameter :: CANON_METRIC_REL_TOL = 1.0d-8
real(dp), parameter :: CANON_CHECK_TOL      = 5.0d-8
integer,  parameter :: MSTEP_STATS_MAGIC = 1180053581
integer,  parameter :: MSTEP_STATS_VERSION = 1
integer,  parameter :: ESTEP_PART_MAGIC = 1180053580
integer,  parameter :: ESTEP_PART_VERSION = 3
integer(longer), parameter :: MSTEP_STATS_IO_TARGET_BYTES = 64_longer * 1024_longer * 1024_longer

type :: projected_latent_mstep_2d_block
    integer :: nrecords = 0
    integer :: ncomp    = 0
    integer,  allocatable :: rows(:), pinds(:)
    logical,  allocatable :: valid(:)
    type(ori), allocatable :: orientations(:)
    real(dp), allocatable :: zrows(:,:)
    real(dp), allocatable :: latent_second(:,:,:)
end type projected_latent_mstep_2d_block

type :: projected_latent_mstep_stats
    integer :: nrecords  = 0
    integer :: ncomp     = 0
    integer :: npairs    = 0
    integer :: exp_lb(3) = 0
    integer :: exp_ub(3) = -1
    integer :: model_nyq = 0
    complex, allocatable :: basis_rhs(:,:,:,:)
    real,    allocatable :: rho_cross(:,:,:,:)
end type projected_latent_mstep_stats

type :: projected_latent_estep_part
    integer :: nrecords = 0
    integer :: ncomp    = 0
    integer :: nmetric  = 0
    integer, allocatable :: rows(:), pinds(:)
    logical, allocatable :: valid(:)
    real(dp), allocatable :: zrows(:,:)
    real(dp), allocatable :: resid_energy(:), resid_mean_energy(:)
    real(dp), allocatable :: basis_metric(:,:)
end type projected_latent_estep_part

contains

    !> Verify that the observation/operator representation used by the
    !! projected latent model is algebraically identical to the standard
    !! refine3D reconstruction plane representation.  In particular,
    !! observation * transfer must equal the standard CTF-weighted numerator
    !! and both paths must have the same CTF-squared denominator.
    subroutine test_projected_model_plane_preparation( params, build, pinds )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:)
        type(parameters) :: test_params
        type(fplane_type), allocatable :: standard_fpls(:), projected_fpls(:)
        integer :: batchlims(2), batchsz, ibatch, i
        integer :: funit, nplanes
        real(dp) :: numerator_err2, numerator_ref2, density_err2, density_ref2
        real(dp) :: numerator_rel, density_rel, numerator_max, density_max
        complex, allocatable :: numerator_diff(:,:)
        real,    allocatable :: density_diff(:,:)
        real(dp), parameter :: REL_TOL = 1.0d-5
        character(len=*), parameter :: METRICS_FILE = 'flex_preimage_plane_preparation_metrics_test.txt'
        if( size(pinds)<1 ) THROW_HARD('flex plane-preparation test requires at least one particle')
        test_params = params
        test_params%lp       = 0.
        test_params%ml_reg   = 'no'
        test_params%l_ml_reg = .false.
        numerator_err2 = 0.d0
        numerator_ref2 = 0.d0
        density_err2   = 0.d0
        density_ref2   = 0.d0
        numerator_max  = 0.d0
        density_max    = 0.d0
        nplanes        = 0
        write(logfhandle,'(A,I0)') '>>> FLEX PRE-IMAGE PLANE CONTRACT TEST PARTICLES: ',size(pinds)
        call init_rec(test_params, build, MAXIMGBATCHSZ, standard_fpls)
        allocate(projected_fpls(MAXIMGBATCHSZ))
        call prepimgbatch(test_params, build, MAXIMGBATCHSZ)
        do ibatch = 1,size(pinds),MAXIMGBATCHSZ
            batchlims = [ibatch,min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
            batchsz   = batchlims(2)-batchlims(1)+1
            call read_particles(test_params, build, size(pinds), pinds, batchlims, batchsz)
            call prep_imgs4rec(test_params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), standard_fpls(:batchsz))
            ! norm_noise_taper_edge_pad_fft intentionally modifies its image
            ! input, so re-read the batch before preparing the flex planes.
            call read_particles(test_params, build, size(pinds), pinds, batchlims, batchsz)
            call prep_imgs4projected_model(test_params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), projected_fpls(:batchsz))
            do i = 1,batchsz
                if( .not.allocated(standard_fpls(i)%cmplx_plane) .or. &
                    &.not.allocated(projected_fpls(i)%cmplx_plane) .or. &
                    &.not.allocated(projected_fpls(i)%transfer_plane) )then
                    THROW_HARD('flex plane-preparation test found an incomplete Fourier plane')
                endif
                if( any(shape(standard_fpls(i)%cmplx_plane)/=shape(projected_fpls(i)%cmplx_plane)) .or. &
                    &any(shape(standard_fpls(i)%ctfsq_plane)/=shape(projected_fpls(i)%ctfsq_plane)) )then
                    THROW_HARD('flex plane-preparation test found incompatible Fourier-plane shapes')
                endif
                numerator_diff = projected_fpls(i)%cmplx_plane * projected_fpls(i)%transfer_plane - &
                    &standard_fpls(i)%cmplx_plane
                density_diff = projected_fpls(i)%ctfsq_plane - standard_fpls(i)%ctfsq_plane
                numerator_err2 = numerator_err2 + sum(real(numerator_diff*conjg(numerator_diff),dp))
                numerator_ref2 = numerator_ref2 + sum(real(standard_fpls(i)%cmplx_plane * &
                    &conjg(standard_fpls(i)%cmplx_plane),dp))
                density_err2 = density_err2 + sum(real(density_diff,dp)*real(density_diff,dp))
                density_ref2 = density_ref2 + sum(real(standard_fpls(i)%ctfsq_plane,dp) * &
                    &real(standard_fpls(i)%ctfsq_plane,dp))
                numerator_max = max(numerator_max,maxval(abs(numerator_diff)))
                density_max   = max(density_max,real(maxval(abs(density_diff)),dp))
                nplanes = nplanes + 1
                deallocate(numerator_diff,density_diff)
            end do
        end do
        numerator_rel = sqrt(numerator_err2/max(numerator_ref2,DTINY))
        density_rel   = sqrt(density_err2/max(density_ref2,DTINY))
        open(newunit=funit,file=METRICS_FILE,status='replace',action='write')
        write(funit,'(A,I0)') 'planes=',nplanes
        write(funit,'(A,ES16.8)') 'numerator_relative_l2=',numerator_rel
        write(funit,'(A,ES16.8)') 'numerator_max_abs=',numerator_max
        write(funit,'(A,ES16.8)') 'ctfsq_relative_l2=',density_rel
        write(funit,'(A,ES16.8)') 'ctfsq_max_abs=',density_max
        write(funit,'(A,ES16.8)') 'relative_tolerance=',REL_TOL
        close(funit)
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
            &'>>> FLEX PRE-IMAGE PLANE CONTRACT numerator_rel=',numerator_rel, &
            &' numerator_max=',numerator_max,' ctfsq_rel=',density_rel,' ctfsq_max=',density_max
        call cleanup_planes(projected_fpls)
        call cleanup_rec_buffers(build,standard_fpls)
        if( numerator_rel>REL_TOL .or. density_rel>REL_TOL )then
            THROW_HARD('flex observation/operator planes differ from standard refine3D preparation; inspect '//METRICS_FILE)
        endif
        write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE PLANE CONTRACT TEST PASSED'
    end subroutine test_projected_model_plane_preparation

    subroutine update_basis_from_latents( params, build, mean_rec, basis_recs, z, pinds, nptcls, ncomp, fpls, log_label )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: nptcls, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(in)    :: z(nptcls,ncomp)
        integer,             intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        character(len=*), optional, intent(in) :: log_label
        type(fplane_type), allocatable :: mean_fpls(:)
        type(image) :: gridcorr_img
        type(projected_latent_mstep_2d_block) :: mstep_block
        real,    allocatable :: rho_cross_exp(:,:,:,:)
        integer, allocatable :: parts(:,:)
        character(len=:), allocatable :: log_prefix
        integer              :: exp_shape(3), npairs
        integer           :: batchlims(2), batchsz, ibatch, ipart, ithr, nparts_eff, partlims(2), q
        integer(timer_int_kind) :: t_total, t_phase, t_comp, t_batch
        real(dp) :: read_prep_seconds, mean_project_seconds, insert_seconds
        t_total = tic()
        if( present(log_label) )then
            log_prefix = projected_model_log_prefix(log_label)
        else
            log_prefix = projected_model_log_prefix()
        endif
        write(logfhandle,'(A)') log_prefix//' M-STEP: UPDATING EIGENVOLUMES WITH COUPLED BLOCK SOLVE'
        call flush(logfhandle)
        read_prep_seconds    = 0.
        mean_project_seconds = 0.
        insert_seconds       = 0.
        nparts_eff      = max(1, min(max(1, params%nparts), nptcls))
        parts           = split_nobjs_even(nptcls, nparts_eff)
        if( nparts_eff > 1 )then
            write(logfhandle,'(A,I0)') log_prefix//' M-STEP LOCAL PARTITIONS: ', nparts_eff
            call flush(logfhandle)
        endif
        do q = 1, ncomp
            call basis_recs(q)%reset
            call basis_recs(q)%reset_exp
        end do
        npairs    = (ncomp * (ncomp + 1)) / 2
        exp_shape = shape(basis_recs(1)%cmat_exp)
        allocate(rho_cross_exp(npairs, exp_shape(1), exp_shape(2), exp_shape(3)), source=0.)
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call init_projected_latent_mstep_2d_block(mstep_block, MAXIMGBATCHSZ, ncomp)
        allocate(mean_fpls(nthr_glob))
        t_phase = tic()
        do ipart = 1, nparts_eff
            partlims = parts(ipart,:)
            do ibatch = partlims(1), partlims(2), MAXIMGBATCHSZ
                batchlims = [ibatch, min(partlims(2), ibatch + MAXIMGBATCHSZ - 1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                t_batch = tic()
                call read_particles(params, build, nptcls, pinds, batchlims, batchsz)
                call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                    &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
                read_prep_seconds = read_prep_seconds + toc(t_batch)
                t_batch = tic()
                call prepare_projected_latent_mstep_2d_block(params, build, mean_rec, fpls(:batchsz), z, &
                    &pinds, batchlims, batchsz, ncomp, mstep_block, mean_fpls)
                mean_project_seconds = mean_project_seconds + toc(t_batch)
                t_batch = tic()
                call insert_projected_latent_mstep_2d_block(build, basis_recs, rho_cross_exp, ncomp, &
                    &mstep_block, fpls(:batchsz))
                insert_seconds = insert_seconds + toc(t_batch)
                write(logfhandle,'(A,I0,A,I0,A,F9.3,A,F9.3,A,F9.3)') log_prefix//' M-STEP BATCH: ', &
                    &batchlims(2), ' / ', nptcls, ' read_prep_total_s=', read_prep_seconds, &
                    &' mean_project_total_s=', mean_project_seconds, ' insert_total_s=', insert_seconds
                call flush(logfhandle)
            end do
        end do
        call log_seconds(log_prefix//' M-STEP INSERT SECONDS', toc(t_phase))
        call log_seconds(log_prefix//' M-STEP READ/PREP SECONDS', read_prep_seconds)
        call log_seconds(log_prefix//' M-STEP MEAN PROJECT SECONDS', mean_project_seconds)
        call log_seconds(log_prefix//' M-STEP COUPLED ACCUMULATE SECONDS', insert_seconds)
        call kill_projected_latent_mstep_2d_block(mstep_block)
        if( allocated(parts) ) deallocate(parts)
        call cleanup_runtime_batch(build, fpls)
        do ithr = 1, size(mean_fpls)
            call cleanup_plane(mean_fpls(ithr))
        end do
        deallocate(mean_fpls)
        t_phase = tic()
        call solve_coupled_basis_exp(basis_recs, rho_cross_exp, ncomp)
        call log_seconds(log_prefix//' M-STEP COUPLED SOLVE SECONDS', toc(t_phase))
        deallocate(rho_cross_exp)
        t_phase = tic()
        gridcorr_img = prep3D_inv_instrfun4mul([params%box_crop,params%box_crop,params%box_crop], &
            &OSMPL_PAD_FAC*[params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        do q = 1, ncomp
            write(logfhandle,'(A,I0,A,I0)') log_prefix//' M-STEP FINALIZE COMPONENT ', q, ' / ', ncomp
            call flush(logfhandle)
            t_comp = tic()
            call finalize_basis_for_projection(params, basis_recs(q), gridcorr_img, density_corrected=.true.)
            call log_comp_seconds(log_prefix//' M-STEP FINALIZE SECONDS', q, toc(t_comp))
        end do
        call gridcorr_img%kill
        call log_seconds(log_prefix//' M-STEP FINALIZE TOTAL SECONDS', toc(t_phase))
        call log_seconds(log_prefix//' M-STEP TOTAL SECONDS', toc(t_total))
    end subroutine update_basis_from_latents

    subroutine init_projected_latent_mstep_2d_block( block, nrecords_max, ncomp )
        type(projected_latent_mstep_2d_block), intent(inout) :: block
        integer, intent(in) :: nrecords_max, ncomp
        call kill_projected_latent_mstep_2d_block(block)
        block%nrecords = 0
        block%ncomp    = ncomp
        allocate(block%rows(nrecords_max), block%pinds(nrecords_max), block%valid(nrecords_max), &
            &block%orientations(nrecords_max), block%zrows(ncomp,nrecords_max), &
            &block%latent_second(ncomp,ncomp,nrecords_max))
        block%rows          = 0
        block%pinds         = 0
        block%valid         = .false.
        block%zrows         = 0.d0
        block%latent_second = 0.d0
    end subroutine init_projected_latent_mstep_2d_block

    subroutine kill_projected_latent_mstep_2d_block( block )
        type(projected_latent_mstep_2d_block), intent(inout) :: block
        integer :: i
        if( allocated(block%orientations) )then
            do i = 1, size(block%orientations)
                call block%orientations(i)%kill
            end do
            deallocate(block%orientations)
        endif
        if( allocated(block%rows) ) deallocate(block%rows)
        if( allocated(block%pinds) ) deallocate(block%pinds)
        if( allocated(block%valid) ) deallocate(block%valid)
        if( allocated(block%zrows) ) deallocate(block%zrows)
        if( allocated(block%latent_second) ) deallocate(block%latent_second)
        block%nrecords = 0
        block%ncomp    = 0
    end subroutine kill_projected_latent_mstep_2d_block

    subroutine reset_projected_latent_mstep_2d_block( block, nrecords )
        type(projected_latent_mstep_2d_block), intent(inout) :: block
        integer, intent(in) :: nrecords
        integer :: i
        if( .not. allocated(block%valid) ) THROW_HARD('unallocated M-step 2D block')
        if( nrecords > size(block%valid) ) THROW_HARD('M-step 2D block capacity exceeded')
        do i = 1, size(block%orientations)
            call block%orientations(i)%kill
        end do
        block%nrecords = nrecords
        block%rows(:nrecords)          = 0
        block%pinds(:nrecords)         = 0
        block%valid(:nrecords)         = .false.
        block%zrows(:,:nrecords)       = 0.d0
        block%latent_second(:,:,:nrecords) = 0.d0
    end subroutine reset_projected_latent_mstep_2d_block

    subroutine prepare_projected_latent_mstep_2d_block( params, build, mean_rec, fpls_batch, z, &
        &pinds, batchlims, batchsz, ncomp, block, mean_fpls )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: batchlims(2), batchsz, ncomp
        type(fplane_type),   intent(inout) :: fpls_batch(batchsz)
        real(dp),            intent(in)    :: z(:,:)
        integer,             intent(in)    :: pinds(:)
        type(projected_latent_mstep_2d_block), intent(inout) :: block
        type(fplane_type), intent(inout) :: mean_fpls(nthr_glob)
        integer :: i, iptcl, q, r, row, ithr
        call reset_projected_latent_mstep_2d_block(block, batchsz)
        !$omp parallel do default(shared) private(i,row,iptcl,q,r,ithr) schedule(static) proc_bind(close)
        do i = 1, batchsz
            row   = batchlims(1) + i - 1
            iptcl = pinds(row)
            ithr  = omp_get_thread_num() + 1
            block%rows(i)  = row
            block%pinds(i) = iptcl
            call build%spproj_field%get_ori(iptcl, block%orientations(i))
            if( block%orientations(i)%isstatezero() ) cycle
            call project_fplane_mean(mean_rec, block%orientations(i), fpls_batch(i), mean_fpls(ithr), apply_ctf_amp=.true.)
            call subtract_plane(fpls_batch(i), mean_fpls(ithr))
            block%zrows(:,i) = z(row,:)
            block%latent_second(:,:,i) = 0.d0
            do q = 1, ncomp
                do r = 1, ncomp
                    block%latent_second(q,r,i) = block%latent_second(q,r,i) + z(row,q) * z(row,r)
                end do
            end do
            block%valid(i) = .true.
        end do
        !$omp end parallel do
    end subroutine prepare_projected_latent_mstep_2d_block

    subroutine insert_projected_latent_mstep_2d_block( build, basis_recs, rho_cross_exp, ncomp, block, fpls_batch )
        class(builder),      intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real,                intent(inout) :: rho_cross_exp(:,:,:,:)
        type(projected_latent_mstep_2d_block), intent(inout) :: block
        type(fplane_type),   intent(in)    :: fpls_batch(:)
        call insert_planes_oversamp_coupled_batch_scaled(basis_recs,rho_cross_exp,build%pgrpsyms, &
            &block%orientations,fpls_batch,block%zrows,block%latent_second,block%valid,block%nrecords)
    end subroutine insert_projected_latent_mstep_2d_block

    subroutine init_projected_latent_mstep_stats( stats, exp_lb, exp_ub, model_nyq, ncomp )
        type(projected_latent_mstep_stats), intent(inout) :: stats
        integer, intent(in) :: exp_lb(3), exp_ub(3), model_nyq, ncomp
        integer :: exp_shape(3)
        call kill_projected_latent_mstep_stats(stats)
        exp_shape = exp_ub - exp_lb + 1
        if( ncomp < 1 .or. model_nyq < 1 .or. any(exp_shape < 1) )then
            THROW_HARD('invalid projected latent M-step statistics dimensions')
        endif
        stats%ncomp     = ncomp
        stats%npairs    = (ncomp * (ncomp + 1)) / 2
        stats%exp_lb    = exp_lb
        stats%exp_ub    = exp_ub
        stats%model_nyq = model_nyq
        allocate(stats%basis_rhs(exp_shape(1),exp_shape(2),exp_shape(3),ncomp), source=CMPLX_ZERO)
        allocate(stats%rho_cross(stats%npairs,exp_shape(1),exp_shape(2),exp_shape(3)), source=0.)
    end subroutine init_projected_latent_mstep_stats

    subroutine kill_projected_latent_mstep_stats( stats )
        type(projected_latent_mstep_stats), intent(inout) :: stats
        if( allocated(stats%basis_rhs) ) deallocate(stats%basis_rhs)
        if( allocated(stats%rho_cross) ) deallocate(stats%rho_cross)
        stats%nrecords  = 0
        stats%ncomp     = 0
        stats%npairs    = 0
        stats%exp_lb    = 0
        stats%exp_ub    = -1
        stats%model_nyq = 0
    end subroutine kill_projected_latent_mstep_stats

    integer(longer) function projected_latent_mstep_stats_nbytes( stats ) result( nbytes )
        type(projected_latent_mstep_stats), intent(in) :: stats
        integer(longer) :: ncells, cbytes, rbytes
        integer :: exp_shape(3)
        nbytes = 0_longer
        if( stats%ncomp < 1 .or. stats%npairs < 1 ) return
        exp_shape = stats%exp_ub - stats%exp_lb + 1
        if( any(exp_shape < 1) ) return
        ncells = product(int(exp_shape, longer))
        cbytes = int(storage_size(CMPLX_ZERO) / 8, longer)
        rbytes = int(storage_size(0.) / 8, longer)
        nbytes = ncells * (cbytes * int(stats%ncomp,longer) + rbytes * int(stats%npairs,longer))
    end function projected_latent_mstep_stats_nbytes

    subroutine accumulate_projected_latent_mstep_2d_block( build, stats, block, fpls_batch )
        class(builder), intent(inout) :: build
        type(projected_latent_mstep_stats), intent(inout) :: stats
        type(projected_latent_mstep_2d_block), intent(inout) :: block
        type(fplane_type), intent(in) :: fpls_batch(:)
        if( block%ncomp /= stats%ncomp ) THROW_HARD('M-step block/statistics component mismatch')
        call accumulate_planes_oversamp_coupled_stats_batch(stats%basis_rhs,stats%rho_cross,stats%exp_lb, &
            &stats%model_nyq,build%pgrpsyms,block%orientations,fpls_batch,block%zrows, &
            &block%latent_second,block%valid,block%nrecords)
        stats%nrecords=stats%nrecords+count(block%valid(:block%nrecords))
    end subroutine accumulate_projected_latent_mstep_2d_block

    subroutine write_projected_latent_mstep_stats( fname, stats )
        class(string), intent(in) :: fname
        type(projected_latent_mstep_stats), intent(in) :: stats
        type(string) :: tmp_fname
        integer :: funit, io_stat, q, header(12)
        if( .not. allocated(stats%basis_rhs) .or. .not. allocated(stats%rho_cross) )then
            THROW_HARD('unallocated projected latent M-step statistics')
        endif
        header = [MSTEP_STATS_MAGIC, MSTEP_STATS_VERSION, stats%nrecords, stats%ncomp, stats%npairs, &
            &stats%exp_lb, stats%exp_ub, stats%model_nyq]
        tmp_fname = fname%to_char()//'.tmp'
        call del_file(tmp_fname)
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_projected_latent_mstep_stats; open '//tmp_fname%to_char(), io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_projected_latent_mstep_stats; header '//tmp_fname%to_char(), io_stat)
        do q = 1, stats%ncomp
            write(funit, iostat=io_stat) stats%basis_rhs(:,:,:,q)
            call fileiochk('write_projected_latent_mstep_stats; RHS '//tmp_fname%to_char(), io_stat)
        end do
        write(funit, iostat=io_stat) stats%rho_cross
        call fileiochk('write_projected_latent_mstep_stats; cross density '//tmp_fname%to_char(), io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname, overwrite=.true.)
        call tmp_fname%kill
    end subroutine write_projected_latent_mstep_stats

    subroutine write_mstep_stats_part_file( params, build, mean_rec, z, pinds, nptcls, ncomp, &
        &partlims, fname, fpls, log_label )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: nptcls, ncomp, partlims(2)
        real(dp),            intent(in)    :: z(nptcls,ncomp)
        integer,             intent(in)    :: pinds(nptcls)
        class(string),       intent(in)    :: fname
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        character(len=*), optional, intent(in) :: log_label
        type(fplane_type), allocatable :: mean_fpls(:)
        type(projected_latent_mstep_2d_block) :: mstep_block
        type(projected_latent_mstep_stats) :: stats
        character(len=:), allocatable :: log_prefix
        integer :: batchlims(2), batchsz, ibatch, ithr
        integer(timer_int_kind) :: t_total, t_batch
        real(dp) :: read_prep_seconds, mean_project_seconds, insert_seconds
        t_total = tic()
        if( present(log_label) )then
            log_prefix = projected_model_log_prefix(log_label)
        else
            log_prefix = projected_model_log_prefix()
        endif
        if( partlims(1) < 1 .or. partlims(2) > nptcls .or. partlims(1) > partlims(2) )then
            THROW_HARD('invalid M-step statistics part limits')
        endif
        write(logfhandle,'(A,I0,A,I0,A,A)') log_prefix//' M-STEP WORKER ROWS: ', partlims(1), ' / ', partlims(2), &
            &' -> ', fname%to_char()
        call flush(logfhandle)
        read_prep_seconds    = 0.
        mean_project_seconds = 0.
        insert_seconds       = 0.
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call init_projected_latent_mstep_2d_block(mstep_block, MAXIMGBATCHSZ, ncomp)
        allocate(mean_fpls(nthr_glob))
        call init_projected_latent_mstep_stats(stats, lbound(mean_rec%cmat_exp), ubound(mean_rec%cmat_exp), &
            &mean_rec%get_lfny(1), ncomp)
        write(logfhandle,'(A,I0)') log_prefix//' M-STEP WORKER STATISTICS BYTES: ', &
            &projected_latent_mstep_stats_nbytes(stats)
        call flush(logfhandle)
        do ibatch = partlims(1), partlims(2), MAXIMGBATCHSZ
            batchlims = [ibatch, min(partlims(2), ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            t_batch = tic()
            call read_particles(params, build, nptcls, pinds, batchlims, batchsz)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            read_prep_seconds = read_prep_seconds + toc(t_batch)
            t_batch = tic()
            call prepare_projected_latent_mstep_2d_block(params, build, mean_rec, fpls(:batchsz), z, &
                &pinds, batchlims, batchsz, ncomp, mstep_block, mean_fpls)
            mean_project_seconds = mean_project_seconds + toc(t_batch)
            t_batch = tic()
            call accumulate_projected_latent_mstep_2d_block(build, stats, mstep_block, fpls(:batchsz))
            insert_seconds = insert_seconds + toc(t_batch)
            write(logfhandle,'(A,I0,A,I0,A,F9.3,A,F9.3,A,F9.3)') log_prefix//' M-STEP WORKER BATCH: ', &
                &batchlims(2) - partlims(1) + 1, ' / ', partlims(2) - partlims(1) + 1, &
                &' read_prep_total_s=', read_prep_seconds, ' mean_project_total_s=', mean_project_seconds, &
                &' insert_total_s=', insert_seconds
            call flush(logfhandle)
        end do
        call write_projected_latent_mstep_stats(fname, stats)
        call kill_projected_latent_mstep_2d_block(mstep_block)
        call kill_projected_latent_mstep_stats(stats)
        call cleanup_runtime_batch(build, fpls)
        do ithr = 1, size(mean_fpls)
            call cleanup_plane(mean_fpls(ithr))
        end do
        deallocate(mean_fpls)
        call log_seconds(log_prefix//' M-STEP WORKER TOTAL SECONDS', toc(t_total))
    end subroutine write_mstep_stats_part_file

    subroutine update_basis_from_mstep_stats_part_files( params, basis_recs, ncomp, part_fnames, nparts, log_label )
        class(parameters),   intent(in)    :: params
        integer,             intent(in)    :: ncomp, nparts
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        class(string),       intent(in)    :: part_fnames(nparts)
        character(len=*), optional, intent(in) :: log_label
        real, allocatable :: rho_cross_exp(:,:,:,:)
        type(image) :: gridcorr_img
        character(len=:), allocatable :: log_prefix
        integer :: exp_shape(3), npairs, ipart, q
        integer(timer_int_kind) :: t_total, t_phase
        t_total = tic()
        if( present(log_label) )then
            log_prefix = projected_model_log_prefix(log_label)
        else
            log_prefix = projected_model_log_prefix()
        endif
        write(logfhandle,'(A,I0)') log_prefix//' M-STEP MASTER REDUCING STATISTICS PARTS: ', nparts
        call flush(logfhandle)
        do q = 1, ncomp
            call basis_recs(q)%reset
            call basis_recs(q)%reset_exp
        end do
        npairs    = (ncomp * (ncomp + 1)) / 2
        exp_shape = shape(basis_recs(1)%cmat_exp)
        allocate(rho_cross_exp(npairs, exp_shape(1), exp_shape(2), exp_shape(3)), source=0.)
        t_phase = tic()
        do ipart = 1, nparts
            if( .not. file_exists(part_fnames(ipart)) )then
                THROW_HARD('missing M-step statistics part: '//part_fnames(ipart)%to_char())
            endif
            call reduce_projected_latent_mstep_stats_file(part_fnames(ipart), basis_recs, rho_cross_exp, ncomp)
            call del_file(part_fnames(ipart))
        end do
        call log_seconds(log_prefix//' M-STEP MASTER REDUCE SECONDS', toc(t_phase))
        t_phase = tic()
        call solve_coupled_basis_exp(basis_recs, rho_cross_exp, ncomp)
        call log_seconds(log_prefix//' M-STEP MASTER COUPLED SOLVE SECONDS', toc(t_phase))
        deallocate(rho_cross_exp)
        t_phase = tic()
        gridcorr_img = prep3D_inv_instrfun4mul([params%box_crop,params%box_crop,params%box_crop], &
            &OSMPL_PAD_FAC*[params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        do q = 1, ncomp
            write(logfhandle,'(A,I0,A,I0)') log_prefix//' M-STEP MASTER FINALIZE COMPONENT ', q, ' / ', ncomp
            call flush(logfhandle)
            call finalize_basis_for_projection(params, basis_recs(q), gridcorr_img, density_corrected=.true.)
        end do
        call gridcorr_img%kill
        call log_seconds(log_prefix//' M-STEP MASTER FINALIZE SECONDS', toc(t_phase))
        call log_seconds(log_prefix//' M-STEP MASTER TOTAL SECONDS', toc(t_total))
    end subroutine update_basis_from_mstep_stats_part_files

    subroutine reduce_projected_latent_mstep_stats_file( fname, basis_recs, rho_cross_exp, ncomp )
        class(string), intent(in) :: fname
        integer, intent(in) :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real, intent(inout) :: rho_cross_exp(:,:,:,:)
        complex, allocatable :: rhs_buf(:,:,:)
        real,    allocatable :: rho_buf(:,:,:,:)
        integer :: funit, io_stat, header(12), exp_lb(3), exp_ub(3), exp_shape(3)
        integer :: npairs, q, im1, im2, nread, mlo, mhi, nz_rhs, nz_rho
        integer(longer) :: rhs_plane_bytes, rho_plane_bytes, expected_bytes, file_bytes
        if( .not. file_exists(fname) ) THROW_HARD('missing M-step statistics file: '//fname%to_char())
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('reduce_projected_latent_mstep_stats_file; open '//fname%to_char(), io_stat)
        read(funit, iostat=io_stat) header
        call fileiochk('reduce_projected_latent_mstep_stats_file; header '//fname%to_char(), io_stat)
        if( header(1) /= MSTEP_STATS_MAGIC ) THROW_HARD('bad M-step statistics magic: '//fname%to_char())
        if( header(2) /= MSTEP_STATS_VERSION ) THROW_HARD('bad M-step statistics version: '//fname%to_char())
        if( header(3) < 0 ) THROW_HARD('negative M-step statistics record count: '//fname%to_char())
        if( header(4) /= ncomp ) THROW_HARD('M-step statistics component mismatch: '//fname%to_char())
        npairs = (ncomp * (ncomp + 1)) / 2
        if( header(5) /= npairs ) THROW_HARD('M-step statistics pair-count mismatch: '//fname%to_char())
        exp_lb    = header(6:8)
        exp_ub    = header(9:11)
        exp_shape = exp_ub - exp_lb + 1
        if( any(exp_shape < 1) ) THROW_HARD('invalid M-step statistics bounds: '//fname%to_char())
        if( any(exp_lb /= lbound(basis_recs(1)%cmat_exp)) .or. &
            &any(exp_ub /= ubound(basis_recs(1)%cmat_exp)) )then
            THROW_HARD('M-step statistics grid mismatch: '//fname%to_char())
        endif
        if( header(12) /= basis_recs(1)%get_lfny(1) )then
            THROW_HARD('M-step statistics Nyquist mismatch: '//fname%to_char())
        endif
        expected_bytes = int(size(header),longer) * int(storage_size(header(1))/8,longer) + &
            &product(int(exp_shape,longer)) * (int(storage_size(CMPLX_ZERO)/8,longer) * int(ncomp,longer) + &
            &int(storage_size(0.)/8,longer) * int(npairs,longer))
        inquire(file=fname%to_char(), size=file_bytes, iostat=io_stat)
        call fileiochk('reduce_projected_latent_mstep_stats_file; size '//fname%to_char(), io_stat)
        if( file_bytes /= expected_bytes ) THROW_HARD('truncated or oversized M-step statistics file: '//fname%to_char())
        rhs_plane_bytes = int(storage_size(CMPLX_ZERO)/8,longer) * int(exp_shape(1),longer) * int(exp_shape(2),longer)
        rho_plane_bytes = int(storage_size(0.)/8,longer) * int(npairs,longer) * &
            &int(exp_shape(1),longer) * int(exp_shape(2),longer)
        nz_rhs = max(1, min(exp_shape(3), int(MSTEP_STATS_IO_TARGET_BYTES / max(1_longer,rhs_plane_bytes))))
        nz_rho = max(1, min(exp_shape(3), int(MSTEP_STATS_IO_TARGET_BYTES / max(1_longer,rho_plane_bytes))))
        allocate(rhs_buf(exp_shape(1),exp_shape(2),nz_rhs))
        do q = 1, ncomp
            do im1 = 1, exp_shape(3), nz_rhs
                im2   = min(exp_shape(3), im1 + nz_rhs - 1)
                nread = im2 - im1 + 1
                read(funit, iostat=io_stat) rhs_buf(:,:,1:nread)
                call fileiochk('reduce_projected_latent_mstep_stats_file; RHS '//fname%to_char(), io_stat)
                mlo = exp_lb(3) + im1 - 1
                mhi = exp_lb(3) + im2 - 1
                basis_recs(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),mlo:mhi) = &
                    &basis_recs(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),mlo:mhi) + rhs_buf(:,:,1:nread)
            end do
        end do
        deallocate(rhs_buf)
        allocate(rho_buf(npairs,exp_shape(1),exp_shape(2),nz_rho))
        do im1 = 1, exp_shape(3), nz_rho
            im2   = min(exp_shape(3), im1 + nz_rho - 1)
            nread = im2 - im1 + 1
            read(funit, iostat=io_stat) rho_buf(:,:,:,1:nread)
            call fileiochk('reduce_projected_latent_mstep_stats_file; cross density '//fname%to_char(), io_stat)
            rho_cross_exp(:,:,:,im1:im2) = rho_cross_exp(:,:,:,im1:im2) + rho_buf(:,:,:,1:nread)
        end do
        deallocate(rho_buf)
        call fclose(funit)
    end subroutine reduce_projected_latent_mstep_stats_file

    subroutine test_projected_latent_mstep_stats_io()
        integer, parameter :: TEST_BOX = 16, TEST_NCOMP = 3
        type(projected_latent_mstep_stats) :: stats
        type(reconstructor) :: basis_recs(TEST_NCOMP)
        type(string) :: fname
        character(len=256) :: padded_fname
        real, allocatable :: rho_cross(:,:,:,:)
        integer :: exp_lb(3), exp_ub(3), exp_shape(3), npairs, q, ipair, i, j, k
        real :: rhs_err, rho_err
        exp_lb    = [-3, -TEST_BOX/2-2, -TEST_BOX/2-2]
        exp_ub    = [ TEST_BOX/2+2, TEST_BOX/2+2, TEST_BOX/2+2]
        exp_shape = exp_ub - exp_lb + 1
        npairs    = (TEST_NCOMP * (TEST_NCOMP + 1)) / 2
        do q = 1, TEST_NCOMP
            call basis_recs(q)%new([TEST_BOX,TEST_BOX,TEST_BOX], 1.)
            allocate(basis_recs(q)%cmat_exp(exp_lb(1):exp_ub(1),exp_lb(2):exp_ub(2),exp_lb(3):exp_ub(3)), &
                &source=CMPLX_ZERO)
        end do
        call init_projected_latent_mstep_stats(stats, exp_lb, exp_ub, basis_recs(1)%get_lfny(1), TEST_NCOMP)
        stats%nrecords = 7
        do q = 1, TEST_NCOMP
            do k = 1, exp_shape(3)
                do j = 1, exp_shape(2)
                    do i = 1, exp_shape(1)
                        stats%basis_rhs(i,j,k,q) = cmplx(0.001*real(i+2*j+3*k+5*q), -0.002*real(2*i-j+k+q))
                    end do
                end do
            end do
        end do
        do ipair = 1, npairs
            stats%rho_cross(ipair,:,:,:) = 0.01 * real(ipair)
        end do
        allocate(rho_cross(npairs,exp_shape(1),exp_shape(2),exp_shape(3)), source=0.)
        padded_fname = 'test_projected_latent_mstep_stats_io.bin'
        fname = padded_fname
        call write_projected_latent_mstep_stats(fname, stats)
        call reduce_projected_latent_mstep_stats_file(fname, basis_recs, rho_cross, TEST_NCOMP)
        rhs_err = 0.
        do q = 1, TEST_NCOMP
            rhs_err = max(rhs_err, maxval(abs(basis_recs(q)%cmat_exp - stats%basis_rhs(:,:,:,q))))
        end do
        rho_err = maxval(abs(rho_cross - stats%rho_cross))
        if( rhs_err > epsilon(1.) ) THROW_HARD('projected latent M-step statistics RHS I/O mismatch')
        if( rho_err > epsilon(1.) ) THROW_HARD('projected latent M-step statistics density I/O mismatch')
        call del_file(fname)
        call fname%kill
        call kill_projected_latent_mstep_stats(stats)
        deallocate(rho_cross)
        do q = 1, TEST_NCOMP
            if( allocated(basis_recs(q)%cmat_exp) ) deallocate(basis_recs(q)%cmat_exp)
            call basis_recs(q)%kill
        end do
    end subroutine test_projected_latent_mstep_stats_io

    subroutine write_estep_latent_part_file( params, build, mean_rec, basis_recs, pinds, nptcls, ncomp, &
        &partlims, fname, fpls, log_label )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: nptcls, ncomp, partlims(2)
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        integer,             intent(in)    :: pinds(nptcls)
        class(string),       intent(in)    :: fname
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        character(len=*), optional, intent(in) :: log_label
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpls(:)
        type(ori),         allocatable :: orientations(:)
        type(projected_latent_estep_part) :: estep_part, estep_out
        complex(dp), allocatable :: gram_h(:,:,:), rhs_h(:,:)
        real(dp),    allocatable :: gram(:,:,:), rhs(:,:), zrow(:,:), basis_metric_thread(:,:,:)
        integer,     allocatable :: metric_count_thread(:)
        character(len=:), allocatable :: log_prefix
        integer :: batchlims(2), batchsz, ibatch, ithr, q, irec, i, progress_stride
        integer(timer_int_kind) :: t_total
        t_total = tic()
        if( present(log_label) )then
            log_prefix = projected_model_log_prefix(log_label)
        else
            log_prefix = projected_model_log_prefix()
        endif
        if( partlims(1) < 1 .or. partlims(2) > nptcls .or. partlims(1) > partlims(2) )then
            THROW_HARD('invalid E-step latent part limits')
        endif
        write(logfhandle,'(A,I0,A,I0,A,A)') log_prefix//' E-STEP WORKER ROWS: ', partlims(1), ' / ', partlims(2), &
            &' -> ', fname%to_char()
        call flush(logfhandle)
        progress_stride = max(1, 5 * MAXIMGBATCHSZ)
        allocate(basis_fpls(ncomp,nthr_glob), mean_fpls(nthr_glob), orientations(nthr_glob), &
            &gram_h(ncomp,ncomp,nthr_glob), rhs_h(ncomp,nthr_glob), gram(ncomp,ncomp,nthr_glob), &
            &rhs(ncomp,nthr_glob), zrow(ncomp,nthr_glob), &
            &basis_metric_thread(ncomp,ncomp,nthr_glob), metric_count_thread(nthr_glob))
        basis_metric_thread = 0.d0
        metric_count_thread = 0
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call init_projected_latent_estep_part(estep_part, MAXIMGBATCHSZ, ncomp)
        call init_projected_latent_estep_part(estep_out, partlims(2) - partlims(1) + 1, ncomp)
        estep_out%nrecords = 0
        do ibatch = partlims(1), partlims(2), MAXIMGBATCHSZ
            batchlims = [ibatch, min(partlims(2), ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call read_particles(params, build, nptcls, pinds, batchlims, batchsz)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            call prepare_projected_latent_estep_part(build, mean_rec, basis_recs, fpls(:batchsz), pinds, &
                &batchlims, batchsz, ncomp, estep_part, basis_fpls, mean_fpls, orientations, &
                &gram_h, rhs_h, gram, rhs, zrow, basis_metric_thread, metric_count_thread)
            do i = 1, estep_part%nrecords
                if( .not. estep_part%valid(i) ) cycle
                if( estep_out%nrecords >= size(estep_out%valid) ) THROW_HARD('E-step output part capacity exceeded')
                irec = estep_out%nrecords + 1
                estep_out%nrecords = irec
                estep_out%rows(irec)              = estep_part%rows(i)
                estep_out%pinds(irec)             = estep_part%pinds(i)
                estep_out%valid(irec)             = .true.
                estep_out%zrows(:,irec)           = estep_part%zrows(:,i)
                estep_out%resid_energy(irec)      = estep_part%resid_energy(i)
                estep_out%resid_mean_energy(irec) = estep_part%resid_mean_energy(i)
            end do
            if( batchlims(2) == partlims(2) .or. mod(batchlims(2) - partlims(1) + 1, progress_stride) == 0 )then
                write(logfhandle,'(A,I0,A,I0)') log_prefix//' E-STEP WORKER PARTICLES: ', &
                    &batchlims(2) - partlims(1) + 1, ' / ', partlims(2) - partlims(1) + 1
                call flush(logfhandle)
            endif
        end do
        estep_out%basis_metric = sum(basis_metric_thread, dim=3)
        estep_out%nmetric      = sum(metric_count_thread)
        call del_file(fname)
        call write_projected_latent_estep_part(fname, estep_out)
        call kill_projected_latent_estep_part(estep_part)
        call kill_projected_latent_estep_part(estep_out)
        do ithr = 1, nthr_glob
            call orientations(ithr)%kill
            call cleanup_plane(mean_fpls(ithr))
            do q = 1, ncomp
                call cleanup_plane(basis_fpls(q,ithr))
            end do
        end do
        call cleanup_runtime_batch(build, fpls)
        deallocate(basis_fpls, mean_fpls, orientations, gram_h, rhs_h, gram, rhs, zrow, &
            &basis_metric_thread, metric_count_thread)
        call log_seconds(log_prefix//' E-STEP WORKER TOTAL SECONDS', toc(t_total))
    end subroutine write_estep_latent_part_file

    subroutine reduce_estep_latent_part_files( part_fnames, nparts, z, resid_energy, resid_mean_energy, &
        &nptcls, ncomp, log_label, basis_metric, metric_valid_count )
        integer,       intent(in)    :: nparts, nptcls, ncomp
        class(string), intent(in)    :: part_fnames(nparts)
        real(dp),      intent(inout) :: z(nptcls,ncomp)
        real(dp),      intent(inout) :: resid_energy(nptcls), resid_mean_energy(nptcls)
        character(len=*), optional, intent(in) :: log_label
        real(dp), optional, intent(out) :: basis_metric(ncomp,ncomp)
        integer,  optional, intent(out) :: metric_valid_count
        type(projected_latent_estep_part) :: part
        real(dp) :: basis_metric_sum(ncomp,ncomp)
        character(len=:), allocatable :: log_prefix
        integer :: ipart, nmetric
        integer(timer_int_kind) :: t_total
        t_total = tic()
        if( present(log_label) )then
            log_prefix = projected_model_log_prefix(log_label)
        else
            log_prefix = projected_model_log_prefix()
        endif
        write(logfhandle,'(A,I0)') log_prefix//' E-STEP MASTER REDUCING LATENT PARTS: ', nparts
        call flush(logfhandle)
        z = 0.d0
        resid_energy = 0.d0
        resid_mean_energy = 0.d0
        basis_metric_sum = 0.d0
        nmetric = 0
        do ipart = 1, nparts
            if( .not. file_exists(part_fnames(ipart)) )then
                THROW_HARD('missing E-step latent part: '//part_fnames(ipart)%to_char())
            endif
            call read_projected_latent_estep_part(part_fnames(ipart), part)
            if( part%ncomp /= ncomp ) THROW_HARD('E-step latent part component count mismatch')
            call reduce_projected_latent_estep_part(part, z, resid_energy, resid_mean_energy)
            basis_metric_sum = basis_metric_sum + part%basis_metric
            nmetric = nmetric + part%nmetric
            call kill_projected_latent_estep_part(part)
        end do
        if( present(basis_metric) )then
            if( nmetric <= 0 ) THROW_HARD('E-step reduction produced an empty basis metric')
            basis_metric = basis_metric_sum / real(nmetric, dp)
        endif
        if( present(metric_valid_count) ) metric_valid_count = nmetric
        call log_seconds(log_prefix//' E-STEP MASTER REDUCE SECONDS', toc(t_total))
    end subroutine reduce_estep_latent_part_files

    subroutine write_projected_latent_estep_part( fname, part )
        class(string), intent(in) :: fname
        type(projected_latent_estep_part), intent(in) :: part
        integer :: funit, io_stat, header(5), nrec
        nrec = part%nrecords
        header = [ESTEP_PART_MAGIC, ESTEP_PART_VERSION, nrec, part%ncomp, part%nmetric]
        call fopen(funit, file=fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_projected_latent_estep_part; open '//fname%to_char(), io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_projected_latent_estep_part; header '//fname%to_char(), io_stat)
        write(funit, iostat=io_stat) part%basis_metric
        call fileiochk('write_projected_latent_estep_part; basis metric '//fname%to_char(), io_stat)
        if( nrec > 0 )then
            write(funit, iostat=io_stat) part%rows(:nrec), part%pinds(:nrec), part%valid(:nrec)
            call fileiochk('write_projected_latent_estep_part; particle fields '//fname%to_char(), io_stat)
            write(funit, iostat=io_stat) part%zrows(:,:nrec)
            call fileiochk('write_projected_latent_estep_part; latent fields '//fname%to_char(), io_stat)
            write(funit, iostat=io_stat) part%resid_energy(:nrec), part%resid_mean_energy(:nrec)
            call fileiochk('write_projected_latent_estep_part; residual fields '//fname%to_char(), io_stat)
        endif
        call fclose(funit)
    end subroutine write_projected_latent_estep_part

    subroutine read_projected_latent_estep_part( fname, part )
        class(string), intent(in) :: fname
        type(projected_latent_estep_part), intent(inout) :: part
        integer :: funit, io_stat, header(5), nrec, ncomp
        call kill_projected_latent_estep_part(part)
        if( .not. file_exists(fname) ) THROW_HARD('missing projected latent E-step part: '//fname%to_char())
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_projected_latent_estep_part; open '//fname%to_char(), io_stat)
        read(funit, iostat=io_stat) header
        call fileiochk('read_projected_latent_estep_part; header '//fname%to_char(), io_stat)
        if( header(1) /= ESTEP_PART_MAGIC ) THROW_HARD('bad projected latent E-step part magic: '//fname%to_char())
        if( header(2) /= ESTEP_PART_VERSION ) THROW_HARD('bad projected latent E-step part version: '//fname%to_char())
        nrec  = header(3)
        ncomp = header(4)
        if( nrec < 0 .or. ncomp < 1 .or. header(5) < 0 )then
            THROW_HARD('invalid projected latent E-step part header: '//fname%to_char())
        endif
        call init_projected_latent_estep_part(part, nrec, ncomp)
        part%nrecords = nrec
        part%nmetric  = header(5)
        read(funit, iostat=io_stat) part%basis_metric
        call fileiochk('read_projected_latent_estep_part; basis metric '//fname%to_char(), io_stat)
        if( nrec > 0 )then
            read(funit, iostat=io_stat) part%rows(:nrec), part%pinds(:nrec), part%valid(:nrec)
            call fileiochk('read_projected_latent_estep_part; particle fields '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) part%zrows(:,:nrec)
            call fileiochk('read_projected_latent_estep_part; latent fields '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) part%resid_energy(:nrec), part%resid_mean_energy(:nrec)
            call fileiochk('read_projected_latent_estep_part; residual fields '//fname%to_char(), io_stat)
        endif
        call fclose(funit)
    end subroutine read_projected_latent_estep_part

    subroutine init_projected_latent_estep_part( part, nrecords_max, ncomp )
        type(projected_latent_estep_part), intent(inout) :: part
        integer, intent(in) :: nrecords_max, ncomp
        call kill_projected_latent_estep_part(part)
        part%nrecords = 0
        part%ncomp    = ncomp
        part%nmetric  = 0
        allocate(part%rows(nrecords_max), part%pinds(nrecords_max), part%valid(nrecords_max), &
            &part%zrows(ncomp,nrecords_max), part%resid_energy(nrecords_max), &
            &part%resid_mean_energy(nrecords_max), part%basis_metric(ncomp,ncomp))
        part%rows              = 0
        part%pinds             = 0
        part%valid             = .false.
        part%zrows             = 0.d0
        part%resid_energy      = 0.d0
        part%resid_mean_energy = 0.d0
        part%basis_metric      = 0.d0
    end subroutine init_projected_latent_estep_part

    subroutine kill_projected_latent_estep_part( part )
        type(projected_latent_estep_part), intent(inout) :: part
        if( allocated(part%rows) ) deallocate(part%rows)
        if( allocated(part%pinds) ) deallocate(part%pinds)
        if( allocated(part%valid) ) deallocate(part%valid)
        if( allocated(part%zrows) ) deallocate(part%zrows)
        if( allocated(part%resid_energy) ) deallocate(part%resid_energy)
        if( allocated(part%resid_mean_energy) ) deallocate(part%resid_mean_energy)
        if( allocated(part%basis_metric) ) deallocate(part%basis_metric)
        part%nrecords = 0
        part%ncomp    = 0
        part%nmetric  = 0
    end subroutine kill_projected_latent_estep_part

    subroutine reset_projected_latent_estep_part( part, nrecords )
        type(projected_latent_estep_part), intent(inout) :: part
        integer, intent(in) :: nrecords
        if( .not. allocated(part%valid) ) THROW_HARD('unallocated E-step part')
        if( nrecords > size(part%valid) ) THROW_HARD('E-step part capacity exceeded')
        part%nrecords = nrecords
        part%nmetric  = 0
        part%rows(:nrecords)              = 0
        part%pinds(:nrecords)             = 0
        part%valid(:nrecords)             = .false.
        part%zrows(:,:nrecords)           = 0.d0
        part%resid_energy(:nrecords)      = 0.d0
        part%resid_mean_energy(:nrecords) = 0.d0
        part%basis_metric                 = 0.d0
    end subroutine reset_projected_latent_estep_part

    subroutine prepare_projected_latent_estep_part( build, mean_rec, basis_recs, fpls_batch, pinds, &
        &batchlims, batchsz, ncomp, part, basis_fpls, mean_fpls, orientations, gram_h, rhs_h, gram, rhs, zrow, &
        &basis_metric_thread, metric_count_thread )
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: batchlims(2), batchsz, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        type(fplane_type),   intent(inout) :: fpls_batch(batchsz)
        integer,             intent(in)    :: pinds(:)
        type(projected_latent_estep_part), intent(inout) :: part
        type(fplane_type), intent(inout) :: basis_fpls(ncomp,nthr_glob), mean_fpls(nthr_glob)
        type(ori),         intent(inout) :: orientations(nthr_glob)
        complex(dp),       intent(inout) :: gram_h(ncomp,ncomp,nthr_glob), rhs_h(ncomp,nthr_glob)
        real(dp),          intent(inout) :: gram(ncomp,ncomp,nthr_glob), rhs(ncomp,nthr_glob)
        real(dp),          intent(inout) :: zrow(ncomp,nthr_glob)
        real(dp),          intent(inout) :: basis_metric_thread(ncomp,ncomp,nthr_glob)
        integer,           intent(inout) :: metric_count_thread(nthr_glob)
        integer :: i, iptcl, q, r, row, ithr
        call reset_projected_latent_estep_part(part, batchsz)
        !$omp parallel do default(shared) private(i,row,iptcl,q,r,ithr) schedule(static) proc_bind(close)
        do i = 1, batchsz
            row   = batchlims(1) + i - 1
            iptcl = pinds(row)
            ithr  = omp_get_thread_num() + 1
            part%rows(i)  = row
            part%pinds(i) = iptcl
            call build%spproj_field%get_ori(iptcl, orientations(ithr))
            if( orientations(ithr)%isstatezero() ) cycle
            call project_fplanes_mean_basis(mean_rec, basis_recs, orientations(ithr), fpls_batch(i), &
                &mean_fpls(ithr), basis_fpls(:,ithr), apply_ctf_amp=.true.)
            call subtract_plane(fpls_batch(i), mean_fpls(ithr))
            part%resid_mean_energy(i) = plane_energy(fpls_batch(i))
            gram_h(:,:,ithr) = DCMPLX_ZERO
            rhs_h(:,ithr)    = DCMPLX_ZERO
            gram(:,:,ithr)   = 0.d0
            rhs(:,ithr)      = 0.d0
            do q = 1, ncomp
                rhs_h(q,ithr) = hermitian_plane_inner_product(basis_fpls(q,ithr), fpls_batch(i))
                rhs(q,ithr)   = real(rhs_h(q,ithr), dp)
                do r = q, ncomp
                    gram_h(q,r,ithr) = hermitian_plane_inner_product(basis_fpls(q,ithr), basis_fpls(r,ithr))
                    gram_h(r,q,ithr) = conjg(gram_h(q,r,ithr))
                    gram(q,r,ithr)   = real(gram_h(q,r,ithr), dp)
                    gram(r,q,ithr)   = gram(q,r,ithr)
                end do
            end do
            basis_metric_thread(:,:,ithr) = basis_metric_thread(:,:,ithr) + gram(:,:,ithr)
            metric_count_thread(ithr) = metric_count_thread(ithr) + 1
            call solve_latent_least_squares(gram(:,:,ithr), rhs(:,ithr), zrow(:,ithr))
            part%zrows(:,i) = zrow(:,ithr)
            do q = 1, ncomp
                call subtract_scaled_plane(fpls_batch(i), basis_fpls(q,ithr), zrow(q,ithr))
            end do
            part%resid_energy(i) = plane_energy(fpls_batch(i))
            part%valid(i) = .true.
        end do
        !$omp end parallel do
    end subroutine prepare_projected_latent_estep_part

    subroutine reduce_projected_latent_estep_part( part, z, resid_energy, resid_mean_energy )
        type(projected_latent_estep_part), intent(in) :: part
        real(dp), intent(inout) :: z(:,:)
        real(dp), intent(inout) :: resid_energy(:), resid_mean_energy(:)
        integer :: i, row
        do i = 1, part%nrecords
            if( .not. part%valid(i) ) cycle
            row = part%rows(i)
            z(row,:)                 = part%zrows(:,i)
            resid_mean_energy(row)   = part%resid_mean_energy(i)
            resid_energy(row)        = part%resid_energy(i)
        end do
    end subroutine reduce_projected_latent_estep_part

    subroutine infer_latents_from_basis( params, build, mean_rec, basis_recs, z, &
        &resid_energy, resid_mean_energy, pinds, nptcls, ncomp, fpls, log_label, &
        &basis_metric, metric_valid_count )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: nptcls, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(inout) :: z(nptcls,ncomp)
        real(dp),            intent(out)   :: resid_energy(nptcls), resid_mean_energy(nptcls)
        integer,             intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        character(len=*), optional, intent(in) :: log_label
        real(dp), optional, intent(out) :: basis_metric(ncomp,ncomp)
        integer,  optional, intent(out) :: metric_valid_count
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpls(:)
        type(ori),         allocatable :: orientations(:)
        type(projected_latent_estep_part) :: estep_part
        complex(dp), allocatable :: gram_h(:,:,:), rhs_h(:,:)
        real(dp),    allocatable :: gram(:,:,:), rhs(:,:), zrow(:,:)
        real(dp),    allocatable :: basis_metric_thread(:,:,:)
        integer,     allocatable :: metric_count_thread(:)
        integer,     allocatable :: parts(:,:)
        character(len=:), allocatable :: log_prefix
        integer           :: batchlims(2), batchsz, ibatch, ipart, ithr, nparts_eff, partlims(2), q, progress_stride
        integer(timer_int_kind) :: t_total, t_phase
        t_total = tic()
        if( present(log_label) )then
            log_prefix = projected_model_log_prefix(log_label)
        else
            log_prefix = projected_model_log_prefix()
        endif
        write(logfhandle,'(A)') log_prefix//' E-STEP: INFERRING LATENTS'
        call flush(logfhandle)
        progress_stride = max(1, 5 * MAXIMGBATCHSZ)
        nparts_eff      = max(1, min(max(1, params%nparts), nptcls))
        parts           = split_nobjs_even(nptcls, nparts_eff)
        if( nparts_eff > 1 )then
            write(logfhandle,'(A,I0)') log_prefix//' E-STEP LOCAL PARTITIONS: ', nparts_eff
            call flush(logfhandle)
        endif
        allocate(basis_fpls(ncomp,nthr_glob), mean_fpls(nthr_glob), orientations(nthr_glob), &
            &gram_h(ncomp,ncomp,nthr_glob), rhs_h(ncomp,nthr_glob), gram(ncomp,ncomp,nthr_glob), &
            &rhs(ncomp,nthr_glob), zrow(ncomp,nthr_glob), &
            &basis_metric_thread(ncomp,ncomp,nthr_glob), metric_count_thread(nthr_glob))
        resid_energy = 0.d0
        resid_mean_energy = 0.d0
        basis_metric_thread = 0.d0
        metric_count_thread = 0
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call init_projected_latent_estep_part(estep_part, MAXIMGBATCHSZ, ncomp)
        t_phase = tic()
        do ipart = 1, nparts_eff
            partlims = parts(ipart,:)
            do ibatch = partlims(1), partlims(2), MAXIMGBATCHSZ
                batchlims = [ibatch, min(partlims(2), ibatch + MAXIMGBATCHSZ - 1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                call read_particles(params, build, nptcls, pinds, batchlims, batchsz)
                call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                    &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
                call prepare_projected_latent_estep_part(build, mean_rec, basis_recs, fpls(:batchsz), pinds, &
                    &batchlims, batchsz, ncomp, estep_part, basis_fpls, mean_fpls, orientations, &
                    &gram_h, rhs_h, gram, rhs, zrow, basis_metric_thread, metric_count_thread)
                call reduce_projected_latent_estep_part(estep_part, z, resid_energy, resid_mean_energy)
                if( batchlims(2) == nptcls .or. mod(batchlims(2), progress_stride) == 0 )then
                    write(logfhandle,'(A,I0,A,I0)') log_prefix//' E-STEP PARTICLES: ', batchlims(2), ' / ', nptcls
                    call flush(logfhandle)
                endif
            end do
        end do
        call log_seconds(log_prefix//' E-STEP INFERENCE SECONDS', toc(t_phase))
        if( present(basis_metric) )then
            if( sum(metric_count_thread) <= 0 ) THROW_HARD('E-step produced an empty basis metric')
            basis_metric = sum(basis_metric_thread, dim=3) / real(sum(metric_count_thread), dp)
        endif
        if( present(metric_valid_count) ) metric_valid_count = sum(metric_count_thread)
        call kill_projected_latent_estep_part(estep_part)
        do ithr = 1, nthr_glob
            call orientations(ithr)%kill
            call cleanup_plane(mean_fpls(ithr))
        end do
        if( allocated(parts) ) deallocate(parts)
        call cleanup_runtime_batch(build, fpls)
        do ithr = 1, nthr_glob
            do q = 1, ncomp
                call cleanup_plane(basis_fpls(q,ithr))
            end do
        end do
        deallocate(basis_fpls, mean_fpls, orientations, gram_h, rhs_h, gram, rhs, zrow, &
            &basis_metric_thread, metric_count_thread)
        call log_seconds(log_prefix//' E-STEP TOTAL SECONDS', toc(t_total))
    end subroutine infer_latents_from_basis

    subroutine canonicalize_projected_latent_basis( basis_recs, z, eigvals, basis_metric, &
        &pinds, nptcls, ncomp, log_label )
        integer,             intent(in)    :: nptcls, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(inout) :: z(nptcls,ncomp)
        real(dp),            intent(out)   :: eigvals(ncomp)
        real(dp),            intent(in)    :: basis_metric(ncomp,ncomp)
        integer,             intent(in)    :: pinds(nptcls)
        character(len=*), optional, intent(in) :: log_label
        real(dp), allocatable :: transform(:,:), inv_transform(:,:)
        real(dp), allocatable :: latent_cov(:,:), check_mat(:,:), identity(:,:), metric_can(:,:)
        real(dp) :: metric_cond, inverse_err, metric_err, covariance_err, explained_sum, anchor_abs, val_abs
        real(dp) :: zrow_work(ncomp)
        real(dp) :: signs(ncomp)
        character(len=:), allocatable :: log_prefix
        integer :: i, q, anchor
        if( nptcls <= 0 .or. ncomp <= 0 ) THROW_HARD('canonicalization requires a nonempty projected latent model')
        if( present(log_label) )then
            log_prefix = projected_model_log_prefix(log_label)
        else
            log_prefix = projected_model_log_prefix()
        endif
        allocate(transform(ncomp,ncomp), inv_transform(ncomp,ncomp), latent_cov(ncomp,ncomp), &
            &check_mat(ncomp,ncomp), identity(ncomp,ncomp), metric_can(ncomp,ncomp), source=0.d0)
        identity  = 0.d0
        do q = 1, ncomp
            identity(q,q)  = 1.d0
        end do
        call latent_covariance(z, nptcls, ncomp, latent_cov)
        call compute_canonical_transform(basis_metric, latent_cov, ncomp, transform, inv_transform, eigvals, metric_cond)
        check_mat = matmul(transform, inv_transform)
        inverse_err = maxval(abs(check_mat - identity))
        metric_can = matmul(transpose(transform), matmul(basis_metric, transform))
        metric_err = maxval(abs(metric_can - identity))
        check_mat = matmul(inv_transform, matmul(latent_cov, transpose(inv_transform)))
        do q = 1, ncomp
            check_mat(q,q) = check_mat(q,q) - eigvals(q)
        end do
        covariance_err = maxval(abs(check_mat)) / max(1.d0, maxval(eigvals))
        if( inverse_err > CANON_CHECK_TOL .or. metric_err > CANON_CHECK_TOL .or. covariance_err > CANON_CHECK_TOL )then
            THROW_HARD('projected latent canonicalization failed its matrix invariants')
        endif

        !$omp parallel do default(shared) schedule(static) private(i,zrow_work) proc_bind(close)
        do i = 1, nptcls
            zrow_work = matmul(inv_transform, z(i,:))
            z(i,:) = zrow_work
        end do
        !$omp end parallel do

        ! A latent eigenvector has arbitrary sign. Anchor each well-defined
        ! component to the largest-magnitude particle coordinate, breaking an
        ! exact tie by the smallest project particle index.
        signs = 1.d0
        do q = 1, ncomp
            anchor = 1
            anchor_abs = abs(z(1,q))
            do i = 2, nptcls
                val_abs = abs(z(i,q))
                if( val_abs > anchor_abs .or. (val_abs == anchor_abs .and. pinds(i) < pinds(anchor)) )then
                    anchor = i
                    anchor_abs = val_abs
                endif
            end do
            if( anchor_abs > DTINY .and. z(anchor,q) < 0.d0 ) signs(q) = -1.d0
            transform(:,q)     = signs(q) * transform(:,q)
            inv_transform(q,:) = signs(q) * inv_transform(q,:)
            z(:,q)             = signs(q) * z(:,q)
        end do
        call mix_projected_latent_basis(basis_recs, transform, ncomp)
        explained_sum = sum(max(0.d0, eigvals))
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') log_prefix//' CANONICAL METRIC CONDITION=', &
            &metric_cond, ' inverse_err=', inverse_err, ' metric_err=', metric_err, ' covariance_err=', covariance_err
        do q = 1, ncomp
            if( explained_sum > DTINY )then
                write(logfhandle,'(A,I0,A,ES12.4,A,F10.6,A,F4.1)') log_prefix//' CANONICAL EIGENVALUE ', q, ': ', &
                    &eigvals(q), ' explained_fraction=', eigvals(q) / explained_sum, ' sign=', signs(q)
            else
                write(logfhandle,'(A,I0,A,ES12.4,A,F4.1)') log_prefix//' CANONICAL EIGENVALUE ', q, ': ', &
                    &eigvals(q), ' sign=', signs(q)
            endif
        end do
        call flush(logfhandle)
        deallocate(transform, inv_transform, latent_cov, check_mat, identity, metric_can)
    end subroutine canonicalize_projected_latent_basis

    subroutine compute_canonical_transform( basis_metric, latent_cov, ncomp, transform, inv_transform, eigvals, metric_cond )
        integer,  intent(in)  :: ncomp
        real(dp), intent(in)  :: basis_metric(ncomp,ncomp), latent_cov(ncomp,ncomp)
        real(dp), intent(out) :: transform(ncomp,ncomp), inv_transform(ncomp,ncomp), eigvals(ncomp), metric_cond
        real(dp) :: metric_work(ncomp,ncomp), metric_vecs(ncomp,ncomp), metric_vals(ncomp)
        real(dp) :: metric_half(ncomp,ncomp), metric_invhalf(ncomp,ncomp)
        real(dp) :: signal_work(ncomp,ncomp), signal_vecs(ncomp,ncomp)
        real(dp) :: metric_max, metric_min, metric_tol, signal_tol, coeff
        integer  :: i, j, q, nrot
        metric_work = 0.5d0 * (basis_metric + transpose(basis_metric))
        nrot = 0
        call jacobi(metric_work, ncomp, ncomp, metric_vals, metric_vecs, nrot)
        call eigsrt(metric_vals, metric_vecs, ncomp, ncomp)
        metric_max = maxval(metric_vals)
        if( metric_max <= DTINY ) THROW_HARD('canonicalization basis metric has no observed support')
        metric_tol = CANON_METRIC_REL_TOL * metric_max
        metric_min = minval(metric_vals)
        if( metric_min <= metric_tol )then
            THROW_HARD('canonicalization basis metric is rank deficient; reduce neigs')
        endif
        metric_cond = metric_max / metric_min
        metric_half    = 0.d0
        metric_invhalf = 0.d0
        do q = 1, ncomp
            do j = 1, ncomp
                do i = 1, ncomp
                    coeff = metric_vecs(i,q) * metric_vecs(j,q)
                    metric_half(i,j)    = metric_half(i,j)    + sqrt(metric_vals(q)) * coeff
                    metric_invhalf(i,j) = metric_invhalf(i,j) + coeff / sqrt(metric_vals(q))
                end do
            end do
        end do
        signal_work = matmul(metric_half, matmul(latent_cov, metric_half))
        signal_work = 0.5d0 * (signal_work + transpose(signal_work))
        nrot = 0
        call jacobi(signal_work, ncomp, ncomp, eigvals, signal_vecs, nrot)
        call eigsrt(eigvals, signal_vecs, ncomp, ncomp)
        signal_tol = CANON_CHECK_TOL * max(1.d0, maxval(abs(eigvals)))
        if( minval(eigvals) < -signal_tol ) THROW_HARD('canonicalization signal covariance is not positive semidefinite')
        eigvals = max(0.d0, eigvals)
        transform     = matmul(metric_invhalf, signal_vecs)
        inv_transform = matmul(transpose(signal_vecs), metric_half)
    end subroutine compute_canonical_transform

    subroutine mix_projected_latent_basis( basis_recs, transform, ncomp )
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(in)    :: transform(ncomp,ncomp)
        complex(dp) :: old_vals(ncomp), new_vals(ncomp)
        integer :: lb(3), ub(3), q, h, k, m
        if( .not. allocated(basis_recs(1)%cmat_exp) ) THROW_HARD('canonicalization received an unallocated basis')
        lb = lbound(basis_recs(1)%cmat_exp)
        ub = ubound(basis_recs(1)%cmat_exp)
        do q = 2, ncomp
            if( .not. allocated(basis_recs(q)%cmat_exp) ) THROW_HARD('canonicalization received an unallocated basis')
            if( any(lbound(basis_recs(q)%cmat_exp) /= lb) .or. any(ubound(basis_recs(q)%cmat_exp) /= ub) )then
                THROW_HARD('canonicalization basis Fourier-grid shape mismatch')
            endif
        end do
        !$omp parallel do collapse(3) default(shared) schedule(static) &
        !$omp private(h,k,m,q,old_vals,new_vals) proc_bind(close)
        do m = lb(3), ub(3)
            do k = lb(2), ub(2)
                do h = lb(1), ub(1)
                    do q = 1, ncomp
                        old_vals(q) = cmplx(basis_recs(q)%cmat_exp(h,k,m), kind=dp)
                    end do
                    new_vals = matmul(old_vals, transform)
                    do q = 1, ncomp
                        basis_recs(q)%cmat_exp(h,k,m) = cmplx(real(new_vals(q),sp), real(aimag(new_vals(q)),sp))
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        do q = 1, ncomp
            call basis_recs(q)%compress_exp
        end do
    end subroutine mix_projected_latent_basis

    subroutine test_projected_latent_canonicalization()
        integer, parameter :: TEST_NCOMP = 3, TEST_NOBS = 4, TEST_NFEAT = 5
        type(projected_latent_estep_part) :: part_in, part_out
        type(string) :: part_fname
        real(dp) :: metric(TEST_NCOMP,TEST_NCOMP), covariance(TEST_NCOMP,TEST_NCOMP)
        real(dp) :: transform(TEST_NCOMP,TEST_NCOMP), inv_transform(TEST_NCOMP,TEST_NCOMP), eigvals(TEST_NCOMP)
        real(dp) :: basis(TEST_NFEAT,TEST_NCOMP), z(TEST_NOBS,TEST_NCOMP), basis_new(TEST_NFEAT,TEST_NCOMP)
        real(dp) :: z_new(TEST_NOBS,TEST_NCOMP), identity(TEST_NCOMP,TEST_NCOMP)
        real(dp) :: check(TEST_NCOMP,TEST_NCOMP), pred(TEST_NOBS,TEST_NFEAT), pred_new(TEST_NOBS,TEST_NFEAT)
        real(dp) :: metric_cond, err
        integer :: q
        metric = reshape([4.d0, 1.d0, 0.5d0, 1.d0, 3.d0, 0.25d0, 0.5d0, 0.25d0, 2.d0], shape(metric))
        covariance = reshape([2.5d0,0.3d0,0.1d0, 0.3d0,1.2d0,0.2d0, 0.1d0,0.2d0,0.4d0], shape(covariance))
        basis = reshape([1.d0,2.d0,3.d0,4.d0,5.d0, 2.d0,-1.d0,0.5d0,1.5d0,-2.d0, &
            &0.25d0,1.25d0,-0.75d0,2.25d0,0.8d0], shape(basis))
        z = reshape([1.d0,0.2d0,-0.4d0,0.5d0, -0.3d0,1.1d0,0.7d0,-0.2d0, &
            &0.8d0,-0.6d0,0.1d0,1.3d0], shape(z))
        call compute_canonical_transform(metric, covariance, TEST_NCOMP, transform, inv_transform, eigvals, metric_cond)
        identity = 0.d0
        do q = 1, TEST_NCOMP
            identity(q,q) = 1.d0
        end do
        err = maxval(abs(matmul(transform, inv_transform) - identity))
        if( err > CANON_CHECK_TOL ) THROW_HARD('canonicalization test inverse mismatch')
        check = matmul(transpose(transform), matmul(metric, transform))
        if( maxval(abs(check - identity)) > CANON_CHECK_TOL ) THROW_HARD('canonicalization test metric mismatch')
        check = matmul(inv_transform, matmul(covariance, transpose(inv_transform)))
        do q = 1, TEST_NCOMP
            check(q,q) = check(q,q) - eigvals(q)
        end do
        if( maxval(abs(check)) > CANON_CHECK_TOL ) THROW_HARD('canonicalization test covariance mismatch')
        if( any(eigvals(2:) > eigvals(:TEST_NCOMP-1)) ) THROW_HARD('canonicalization test eigenvalue ordering mismatch')
        basis_new = matmul(basis, transform)
        z_new     = matmul(z, transpose(inv_transform))
        pred      = matmul(z, transpose(basis))
        pred_new  = matmul(z_new, transpose(basis_new))
        if( maxval(abs(pred - pred_new)) > CANON_CHECK_TOL ) THROW_HARD('canonicalization test model prediction mismatch')
        call init_projected_latent_estep_part(part_in, 2, TEST_NCOMP)
        part_in%nrecords = 2
        part_in%nmetric  = 2
        part_in%rows     = [1, 2]
        part_in%pinds    = [11, 12]
        part_in%valid    = .true.
        part_in%basis_metric = metric
        part_in%zrows(:,1) = z(1,:)
        part_in%zrows(:,2) = z(2,:)
        part_fname = 'test_projected_latent_canonicalization_estep.bin'
        call write_projected_latent_estep_part(part_fname, part_in)
        call read_projected_latent_estep_part(part_fname, part_out)
        if( part_out%nmetric /= part_in%nmetric ) THROW_HARD('canonicalization test E-step metric count mismatch')
        if( maxval(abs(part_out%basis_metric - part_in%basis_metric)) > CANON_CHECK_TOL )then
            THROW_HARD('canonicalization test E-step metric I/O mismatch')
        endif
        if( maxval(abs(part_out%zrows - part_in%zrows)) > CANON_CHECK_TOL )then
            THROW_HARD('canonicalization test E-step latent I/O mismatch')
        endif
        call del_file(part_fname)
        call part_fname%kill
        call kill_projected_latent_estep_part(part_in)
        call kill_projected_latent_estep_part(part_out)
        write(logfhandle,'(A,ES12.4)') '>>> PROJECTED_MODEL CANONICALIZATION TEST PASSED; METRIC CONDITION=', metric_cond
        call flush(logfhandle)
    end subroutine test_projected_latent_canonicalization

    subroutine solve_coupled_basis_exp( basis_recs, rho_cross_exp, ncomp )
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real,                intent(in)    :: rho_cross_exp(:,:,:,:)
        complex(dp) :: rhs(ncomp), sol(ncomp)
        real(dp)    :: amat(ncomp,ncomp)
        real(dp)    :: diag_sum, diag_max, ridge, denom
        integer     :: lb(3), ub(3), h, k, m, ih, ik, im, q, r, flag, shell, nyq
        lb = lbound(basis_recs(1)%cmat_exp)
        ub = ubound(basis_recs(1)%cmat_exp)
        nyq = basis_recs(1)%get_lfny(1)
        !$omp parallel do collapse(3) default(shared) schedule(static) &
        !$omp private(h,k,m,ih,ik,im,q,r,amat,rhs,sol,diag_sum,diag_max,ridge,denom,flag,shell) proc_bind(close)
        do m = lb(3), ub(3)
            do k = lb(2), ub(2)
                do h = lb(1), ub(1)
                    ih = h - lb(1) + 1
                    ik = k - lb(2) + 1
                    im = m - lb(3) + 1
                    ! Match reconstructor%sampl_dens_correct: values outside
                    ! the spherical Nyquist support are not reconstructable,
                    ! even though they lie inside the Cartesian FFT cube.
                    shell = nint(sqrt(real(h*h+k*k+m*m)))
                    if( shell > nyq )then
                        do q = 1, ncomp
                            basis_recs(q)%cmat_exp(h,k,m) = CMPLX_ZERO
                        end do
                        cycle
                    endif
                    rhs  = DCMPLX_ZERO
                    diag_sum = 0.d0
                    diag_max = 0.d0
                    do q = 1, ncomp
                        rhs(q) = cmplx(basis_recs(q)%cmat_exp(h,k,m), kind=dp)
                        denom = max(0.d0,real(rho_cross_exp(pair_index(q,q),ih,ik,im),dp))
                        diag_sum = diag_sum + denom
                        diag_max = max(diag_max,denom)
                    end do
                    if( diag_max <= COUPLED_DENSITY_FLOOR )then
                        do q = 1, ncomp
                            basis_recs(q)%cmat_exp(h,k,m) = CMPLX_ZERO
                        end do
                        cycle
                    endif
                    amat = 0.d0
                    do q = 1, ncomp
                        do r = q, ncomp
                            amat(q,r) = real(rho_cross_exp(pair_index(q,r),ih,ik,im), dp)
                            amat(r,q) = amat(q,r)
                        end do
                    end do
                    ridge = COUPLED_MSTEP_RIDGE_REL * diag_sum / real(max(1,ncomp), dp)
                    do q = 1, ncomp
                        amat(q,q) = amat(q,q) + ridge
                    end do
                    call solve_real_spd_complex(amat, rhs, sol, ncomp, flag)
                    if( flag /= 0 )then
                        do q = 1, ncomp
                            denom = max(abs(amat(q,q)), ridge)
                            if( denom > DTINY )then
                                sol(q) = rhs(q) / denom
                            else
                                sol(q) = DCMPLX_ZERO
                            endif
                        end do
                    endif
                    do q = 1, ncomp
                        basis_recs(q)%cmat_exp(h,k,m) = cmplx(real(sol(q), sp), real(aimag(sol(q)), sp))
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine solve_coupled_basis_exp

    integer pure function pair_index( q, r ) result( ipair )
        integer, intent(in) :: q, r
        ipair = (r * (r - 1)) / 2 + q
    end function pair_index

    subroutine solve_real_spd_complex( amat_in, rhs, sol, n, flag )
        integer,     intent(in)  :: n
        real(dp),    intent(in)  :: amat_in(n,n)
        complex(dp), intent(in)  :: rhs(n)
        complex(dp), intent(out) :: sol(n)
        integer,     intent(out) :: flag
        real(dp) :: chol(n,n), yr(n), yi(n), xr(n), xi(n)
        real(dp) :: sumr, sumi, sumv, tol
        integer  :: i, j, l
        flag = 0
        sol  = DCMPLX_ZERO
        chol = 0.d0
        tol  = max(DTINY, epsilon(1.d0) * max(1.d0, maxval(abs(amat_in))))
        do j = 1, n
            sumv = amat_in(j,j)
            do l = 1, j - 1
                sumv = sumv - chol(j,l) * chol(j,l)
            end do
            if( sumv <= tol )then
                flag = 1
                return
            endif
            chol(j,j) = sqrt(sumv)
            do i = j + 1, n
                sumv = amat_in(i,j)
                do l = 1, j - 1
                    sumv = sumv - chol(i,l) * chol(j,l)
                end do
                chol(i,j) = sumv / chol(j,j)
            end do
        end do
        do i = 1, n
            sumr = real(rhs(i), dp)
            sumi = aimag(rhs(i))
            do l = 1, i - 1
                sumr = sumr - chol(i,l) * yr(l)
                sumi = sumi - chol(i,l) * yi(l)
            end do
            yr(i) = sumr / chol(i,i)
            yi(i) = sumi / chol(i,i)
        end do
        do i = n, 1, -1
            sumr = yr(i)
            sumi = yi(i)
            do l = i + 1, n
                sumr = sumr - chol(l,i) * xr(l)
                sumi = sumi - chol(l,i) * xi(l)
            end do
            xr(i) = sumr / chol(i,i)
            xi(i) = sumi / chol(i,i)
        end do
        do i = 1, n
            sol(i) = cmplx(xr(i), xi(i), kind=dp)
        end do
    end subroutine solve_real_spd_complex

    subroutine finalize_basis_for_projection( params, basis_rec, gridcorr_img, density_corrected )
        class(parameters),   intent(in)    :: params
        type(reconstructor), intent(inout) :: basis_rec
        class(image),        intent(in)    :: gridcorr_img
        logical, optional,   intent(in)    :: density_corrected
        logical :: l_density_corrected
        l_density_corrected = .false.
        if( present(density_corrected) ) l_density_corrected = density_corrected
        call basis_rec%compress_exp
        if( .not. l_density_corrected ) call basis_rec%sampl_dens_correct
        call basis_rec%ifft
        call basis_rec%div(real(params%box))
        call basis_rec%mul(gridcorr_img)
        call regularize_basis_volume(params, basis_rec)
        call basis_rec%fft
        call basis_rec%expand_exp
    end subroutine finalize_basis_for_projection

    subroutine regularize_basis_volume( params, basis_rec )
        class(parameters),   intent(in)    :: params
        type(reconstructor), intent(inout) :: basis_rec
        if( basis_rec%is_ft() ) call basis_rec%ifft
        if( params%lp > 2.0 * params%smpd_crop + TINY )then
            call basis_rec%bp(0., params%lp)
        endif
        if( params%msk_crop > TINY )then
            call basis_rec%mask3D_soft(params%msk_crop, backgr=0.)
        endif
    end subroutine regularize_basis_volume

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

    subroutine prep_imgs4projected_model( params, build, nptcls, ptcl_imgs, pinds, fplanes )
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
        kfromto = projected_model_kfromto(params)
        !$omp parallel do default(shared) private(i,ithr,iptcl,shift) schedule(static) proc_bind(close)
        do i = 1, nptcls
            ithr   = omp_get_thread_num() + 1
            iptcl  = pinds(i)
            call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft(build%lmsk, build%img_pad_heap(ithr))
            ctfparms(ithr) = build%spproj%get_ctfparams(params%oritype, iptcl)
            shift = build%spproj_field%get_2Dshift(iptcl)
            call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr), &
                &shift, fplanes(i), store_transfer=.true., observation_model=.true.)
            call cap_fplane_for_projected_model(fplanes(i), kfromto)
        end do
        !$omp end parallel do
    end subroutine prep_imgs4projected_model

    subroutine cap_fplane_for_projected_model( fpl, kfromto )
        type(fplane_type), intent(inout) :: fpl
        integer,           intent(in)    :: kfromto(2)
        integer :: nyq_eff
        nyq_eff = max(OSMPL_PAD_FAC, OSMPL_PAD_FAC * kfromto(2))
        if( fpl%nyq > 0 ) fpl%nyq = min(fpl%nyq, nyq_eff)
    end subroutine cap_fplane_for_projected_model

    function projected_model_kfromto( params ) result( kfromto )
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
    end function projected_model_kfromto

    subroutine initialize_latents( z, nptcls, ncomp )
        use simple_rnd, only: mnorm_smp
        integer,  intent(in) :: nptcls, ncomp
        real(dp), intent(out) :: z(nptcls,ncomp)
        real :: identity(ncomp,ncomp), means(ncomp)
        integer :: i, q
        identity = 0.
        means    = 0.
        do q = 1, ncomp
            identity(q,q) = 1.
        end do
        do i = 1, nptcls
            z(i,:) = real(mnorm_smp(identity, ncomp, means), dp)
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

    function hermitian_plane_inner_product( lhs_fpl, rhs_fpl ) result( val )
        use simple_math, only: ceil_div, floor_div
        type(fplane_type), intent(in) :: lhs_fpl, rhs_fpl
        complex(dp) :: val
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
                acc = acc + conjg(cmplx(lhs_fpl%cmplx_plane(h,k), kind=dp)) * cmplx(rhs_fpl%cmplx_plane(h,k), kind=dp)
            end do
        end do
        val = acc
    end function hermitian_plane_inner_product

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

    subroutine solve_latent_least_squares( gram, rhs, x )
        real(dp), intent(in)  :: gram(:,:), rhs(:)
        real(dp), intent(out) :: x(:)
        integer :: flag
        call hermitian_solve(gram, rhs, x, flag)
        if( flag /= 0 ) x = 0.d0
    end subroutine solve_latent_least_squares

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

    function projected_model_log_prefix( log_label ) result( prefix )
        character(len=*), optional, intent(in) :: log_label
        character(len=:), allocatable :: prefix
        if( present(log_label) )then
            prefix = '>>> '//trim(log_label)
        else
            prefix = '>>> PROJECTED_LATENT_MODEL'
        endif
    end function projected_model_log_prefix

    subroutine log_seconds( label, seconds )
        character(len=*),      intent(in) :: label
        real(timer_int_kind),  intent(in) :: seconds
        write(logfhandle,'(A,A,F10.3)') trim(label), ': ', seconds
        call flush(logfhandle)
    end subroutine log_seconds

    subroutine log_comp_seconds( label, icomp, seconds )
        character(len=*),     intent(in) :: label
        integer,              intent(in) :: icomp
        real(timer_int_kind), intent(in) :: seconds
        write(logfhandle,'(A,1X,I0,A,F10.3)') trim(label), icomp, ': ', seconds
        call flush(logfhandle)
    end subroutine log_comp_seconds

end module simple_flex_projected_latent_model
