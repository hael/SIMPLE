!@descr: Projection-aware latent volume model kernels
module simple_projected_latent_model
use simple_core_module_api
use simple_builder,          only: builder
use simple_image,            only: image
use simple_imgarr_utils,     only: dealloc_imgarr
use simple_linalg,           only: hermitian_invert, hermitian_solve
use simple_matcher_3Drec,    only: init_rec
use simple_matcher_ptcl_io,  only: discrete_read_imgbatch, discrete_read_imgbatch_source, prepimgbatch, killimgbatch
use simple_memoize_ft_maps,  only: memoize_ft_maps, forget_ft_maps
use simple_parameters,       only: parameters
use simple_reconstructor,    only: reconstructor
use simple_reconstructor_latent_ops, only: insert_plane_oversamp_coupled_scaled, accumulate_plane_oversamp_coupled_stats, &
    &project_fplane_mean, project_fplanes_mean_basis
use simple_map_reduce,       only: split_nobjs_even
implicit none

public :: update_basis_from_latents, infer_latents_from_basis
public :: write_mstep_stats_part_file, update_basis_from_mstep_stats_part_files
public :: write_estep_latent_part_file, reduce_estep_latent_part_files
public :: initialize_latents, orthonormalize_latents, latent_sdev, latent_covariance
public :: basis_fourier_energy, cleanup_planes, projected_model_kfromto
public :: test_projected_latent_mstep_stats_io
private
#include "simple_local_flags.inc"

real(dp), parameter :: LATENT_RIDGE = 1.0d-3
real(dp), parameter :: MODE_VAR_FLOOR = 1.0d-3
real(dp), parameter :: COUPLED_MSTEP_RIDGE_REL = 1.0d-8
real(dp), parameter :: COUPLED_MSTEP_RIDGE_ABS = 1.0d-10
integer,  parameter :: MSTEP_STATS_MAGIC = 1180053581
integer,  parameter :: MSTEP_STATS_VERSION = 1
integer,  parameter :: ESTEP_PART_MAGIC = 1180053580
integer,  parameter :: ESTEP_PART_VERSION = 1
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
    integer, allocatable :: rows(:), pinds(:)
    logical, allocatable :: valid(:)
    real(dp), allocatable :: zrows(:,:)
    real(dp), allocatable :: z_postcov(:,:,:)
    real(dp), allocatable :: resid_energy(:), resid_mean_energy(:)
    real(dp), allocatable :: mode_second(:,:)
end type projected_latent_estep_part

contains

    subroutine update_basis_from_latents( params, build, mean_rec, basis_recs, z, z_postcov, pinds, nptcls, ncomp, fpls, log_label )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: nptcls, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(in)    :: z(nptcls,ncomp), z_postcov(nptcls,ncomp,ncomp)
        integer,             intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        character(len=*), optional, intent(in) :: log_label
        type(fplane_type) :: mean_fpl
        type(projected_latent_mstep_2d_block) :: mstep_block
        real,    allocatable :: rho_cross_exp(:,:,:,:)
        integer, allocatable :: parts(:,:)
        character(len=:), allocatable :: log_prefix
        integer              :: exp_shape(3), npairs
        integer           :: batchlims(2), batchsz, ibatch, ipart, nparts_eff, partlims(2), q, progress_stride
        integer(timer_int_kind) :: t_total, t_phase, t_comp
        t_total = tic()
        if( present(log_label) )then
            log_prefix = projected_model_log_prefix(log_label)
        else
            log_prefix = projected_model_log_prefix()
        endif
        write(logfhandle,'(A)') log_prefix//' M-STEP: UPDATING EIGENVOLUMES WITH COUPLED BLOCK SOLVE'
        call flush(logfhandle)
        progress_stride = max(1, 5 * MAXIMGBATCHSZ)
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
        t_phase = tic()
        do ipart = 1, nparts_eff
            partlims = parts(ipart,:)
            do ibatch = partlims(1), partlims(2), MAXIMGBATCHSZ
                batchlims = [ibatch, min(partlims(2), ibatch + MAXIMGBATCHSZ - 1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                call read_particles(params, build, nptcls, pinds, batchlims, batchsz)
                call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                    &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
                call prepare_projected_latent_mstep_2d_block(params, build, mean_rec, fpls(:batchsz), z, z_postcov, &
                    &pinds, batchlims, batchsz, ncomp, mstep_block, mean_fpl)
                call insert_projected_latent_mstep_2d_block(build, basis_recs, rho_cross_exp, ncomp, &
                    &mstep_block, fpls(:batchsz))
                if( batchlims(2) == nptcls .or. mod(batchlims(2), progress_stride) == 0 )then
                    write(logfhandle,'(A,I0,A,I0)') log_prefix//' M-STEP PARTICLES: ', batchlims(2), ' / ', nptcls
                    call flush(logfhandle)
                endif
            end do
        end do
        call log_seconds(log_prefix//' M-STEP INSERT SECONDS', toc(t_phase))
        call kill_projected_latent_mstep_2d_block(mstep_block)
        if( allocated(parts) ) deallocate(parts)
        call cleanup_runtime_batch(build, fpls)
        call cleanup_plane(mean_fpl)
        t_phase = tic()
        call solve_coupled_basis_exp(basis_recs, rho_cross_exp, ncomp)
        call log_seconds(log_prefix//' M-STEP COUPLED SOLVE SECONDS', toc(t_phase))
        deallocate(rho_cross_exp)
        t_phase = tic()
        do q = 1, ncomp
            write(logfhandle,'(A,I0,A,I0)') log_prefix//' M-STEP FINALIZE COMPONENT ', q, ' / ', ncomp
            call flush(logfhandle)
            t_comp = tic()
            call finalize_basis_for_projection(params, basis_recs(q), density_corrected=.true.)
            call log_comp_seconds(log_prefix//' M-STEP FINALIZE SECONDS', q, toc(t_comp))
        end do
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

    subroutine prepare_projected_latent_mstep_2d_block( params, build, mean_rec, fpls_batch, z, z_postcov, &
        &pinds, batchlims, batchsz, ncomp, block, mean_fpl )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: batchlims(2), batchsz, ncomp
        type(fplane_type),   intent(inout) :: fpls_batch(batchsz)
        real(dp),            intent(in)    :: z(:,:), z_postcov(:,:,:)
        integer,             intent(in)    :: pinds(:)
        type(projected_latent_mstep_2d_block), intent(inout) :: block
        type(fplane_type), intent(inout) :: mean_fpl
        integer :: i, iptcl, q, r, row
        call reset_projected_latent_mstep_2d_block(block, batchsz)
        do i = 1, batchsz
            row   = batchlims(1) + i - 1
            iptcl = pinds(row)
            block%rows(i)  = row
            block%pinds(i) = iptcl
            call build%spproj_field%get_ori(iptcl, block%orientations(i))
            if( block%orientations(i)%isstatezero() ) cycle
            call project_fplane_mean(mean_rec, block%orientations(i), fpls_batch(i), mean_fpl, apply_ctf_amp=.true.)
            call subtract_plane(fpls_batch(i), mean_fpl)
            block%zrows(:,i) = z(row,:)
            block%latent_second(:,:,i) = z_postcov(row,:,:)
            do q = 1, ncomp
                do r = 1, ncomp
                    block%latent_second(q,r,i) = block%latent_second(q,r,i) + z(row,q) * z(row,r)
                end do
            end do
            block%valid(i) = .true.
        end do
    end subroutine prepare_projected_latent_mstep_2d_block

    subroutine insert_projected_latent_mstep_2d_block( build, basis_recs, rho_cross_exp, ncomp, block, fpls_batch )
        class(builder),      intent(inout) :: build
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real,                intent(inout) :: rho_cross_exp(:,:,:,:)
        type(projected_latent_mstep_2d_block), intent(inout) :: block
        type(fplane_type),   intent(in)    :: fpls_batch(:)
        integer :: i
        do i = 1, block%nrecords
            if( .not. block%valid(i) ) cycle
            call insert_plane_oversamp_coupled_scaled(basis_recs, rho_cross_exp, build%pgrpsyms, &
                &block%orientations(i), fpls_batch(i), block%zrows(:,i), block%latent_second(:,:,i))
        end do
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
        integer :: i
        if( block%ncomp /= stats%ncomp ) THROW_HARD('M-step block/statistics component mismatch')
        do i = 1, block%nrecords
            if( .not. block%valid(i) ) cycle
            call accumulate_plane_oversamp_coupled_stats(stats%basis_rhs, stats%rho_cross, stats%exp_lb, &
                &stats%model_nyq, build%pgrpsyms, block%orientations(i), fpls_batch(i), &
                &block%zrows(:,i), block%latent_second(:,:,i))
            stats%nrecords = stats%nrecords + 1
        end do
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

    subroutine write_mstep_stats_part_file( params, build, mean_rec, z, z_postcov, pinds, nptcls, ncomp, &
        &partlims, fname, fpls, log_label )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: nptcls, ncomp, partlims(2)
        real(dp),            intent(in)    :: z(nptcls,ncomp), z_postcov(nptcls,ncomp,ncomp)
        integer,             intent(in)    :: pinds(nptcls)
        class(string),       intent(in)    :: fname
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        character(len=*), optional, intent(in) :: log_label
        type(fplane_type) :: mean_fpl
        type(projected_latent_mstep_2d_block) :: mstep_block
        type(projected_latent_mstep_stats) :: stats
        character(len=:), allocatable :: log_prefix
        integer :: batchlims(2), batchsz, ibatch, progress_stride
        integer(timer_int_kind) :: t_total
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
        progress_stride = max(1, 5 * MAXIMGBATCHSZ)
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call init_projected_latent_mstep_2d_block(mstep_block, MAXIMGBATCHSZ, ncomp)
        call init_projected_latent_mstep_stats(stats, lbound(mean_rec%cmat_exp), ubound(mean_rec%cmat_exp), &
            &mean_rec%get_lfny(1), ncomp)
        write(logfhandle,'(A,I0)') log_prefix//' M-STEP WORKER STATISTICS BYTES: ', &
            &projected_latent_mstep_stats_nbytes(stats)
        call flush(logfhandle)
        do ibatch = partlims(1), partlims(2), MAXIMGBATCHSZ
            batchlims = [ibatch, min(partlims(2), ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call read_particles(params, build, nptcls, pinds, batchlims, batchsz)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            call prepare_projected_latent_mstep_2d_block(params, build, mean_rec, fpls(:batchsz), z, z_postcov, &
                &pinds, batchlims, batchsz, ncomp, mstep_block, mean_fpl)
            call accumulate_projected_latent_mstep_2d_block(build, stats, mstep_block, fpls(:batchsz))
            if( batchlims(2) == partlims(2) .or. mod(batchlims(2) - partlims(1) + 1, progress_stride) == 0 )then
                write(logfhandle,'(A,I0,A,I0)') log_prefix//' M-STEP WORKER PARTICLES: ', &
                    &batchlims(2) - partlims(1) + 1, ' / ', partlims(2) - partlims(1) + 1
                call flush(logfhandle)
            endif
        end do
        call write_projected_latent_mstep_stats(fname, stats)
        call kill_projected_latent_mstep_2d_block(mstep_block)
        call kill_projected_latent_mstep_stats(stats)
        call cleanup_runtime_batch(build, fpls)
        call cleanup_plane(mean_fpl)
        call log_seconds(log_prefix//' M-STEP WORKER TOTAL SECONDS', toc(t_total))
    end subroutine write_mstep_stats_part_file

    subroutine update_basis_from_mstep_stats_part_files( params, basis_recs, ncomp, part_fnames, nparts, log_label )
        class(parameters),   intent(in)    :: params
        integer,             intent(in)    :: ncomp, nparts
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        class(string),       intent(in)    :: part_fnames(nparts)
        character(len=*), optional, intent(in) :: log_label
        real, allocatable :: rho_cross_exp(:,:,:,:)
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
        do q = 1, ncomp
            write(logfhandle,'(A,I0,A,I0)') log_prefix//' M-STEP MASTER FINALIZE COMPONENT ', q, ' / ', ncomp
            call flush(logfhandle)
            call finalize_basis_for_projection(params, basis_recs(q), density_corrected=.true.)
        end do
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

    subroutine write_estep_latent_part_file( params, build, mean_rec, basis_recs, mode_vars, pinds, nptcls, ncomp, &
        &partlims, fname, fpls, log_label )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: nptcls, ncomp, partlims(2)
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(in)    :: mode_vars(ncomp)
        integer,             intent(in)    :: pinds(nptcls)
        class(string),       intent(in)    :: fname
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        character(len=*), optional, intent(in) :: log_label
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpls(:)
        type(ori),         allocatable :: orientations(:)
        type(projected_latent_estep_part) :: estep_part, estep_out
        complex(dp), allocatable :: gram_h(:,:,:), rhs_h(:,:)
        real(dp),    allocatable :: gram(:,:,:), rhs(:,:), zrow(:,:), post_cov(:,:,:)
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
            &rhs(ncomp,nthr_glob), zrow(ncomp,nthr_glob), post_cov(ncomp,ncomp,nthr_glob))
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
            call prepare_projected_latent_estep_part(build, mean_rec, basis_recs, fpls(:batchsz), mode_vars, pinds, &
                &batchlims, batchsz, ncomp, estep_part, basis_fpls, mean_fpls, orientations, &
                &gram_h, rhs_h, gram, rhs, zrow, post_cov)
            do i = 1, estep_part%nrecords
                if( .not. estep_part%valid(i) ) cycle
                if( estep_out%nrecords >= size(estep_out%valid) ) THROW_HARD('E-step output part capacity exceeded')
                irec = estep_out%nrecords + 1
                estep_out%nrecords = irec
                estep_out%rows(irec)              = estep_part%rows(i)
                estep_out%pinds(irec)             = estep_part%pinds(i)
                estep_out%valid(irec)             = .true.
                estep_out%zrows(:,irec)           = estep_part%zrows(:,i)
                estep_out%z_postcov(:,:,irec)     = estep_part%z_postcov(:,:,i)
                estep_out%resid_energy(irec)      = estep_part%resid_energy(i)
                estep_out%resid_mean_energy(irec) = estep_part%resid_mean_energy(i)
                estep_out%mode_second(:,irec)     = estep_part%mode_second(:,i)
            end do
            if( batchlims(2) == partlims(2) .or. mod(batchlims(2) - partlims(1) + 1, progress_stride) == 0 )then
                write(logfhandle,'(A,I0,A,I0)') log_prefix//' E-STEP WORKER PARTICLES: ', &
                    &batchlims(2) - partlims(1) + 1, ' / ', partlims(2) - partlims(1) + 1
                call flush(logfhandle)
            endif
        end do
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
        deallocate(basis_fpls, mean_fpls, orientations, gram_h, rhs_h, gram, rhs, zrow, post_cov)
        call log_seconds(log_prefix//' E-STEP WORKER TOTAL SECONDS', toc(t_total))
    end subroutine write_estep_latent_part_file

    subroutine reduce_estep_latent_part_files( part_fnames, nparts, z, z_postcov, resid_energy, resid_mean_energy, &
        &mode_vars, nptcls, ncomp, log_label )
        integer,       intent(in)    :: nparts, nptcls, ncomp
        class(string), intent(in)    :: part_fnames(nparts)
        real(dp),      intent(inout) :: z(nptcls,ncomp), z_postcov(nptcls,ncomp,ncomp)
        real(dp),      intent(inout) :: resid_energy(nptcls), resid_mean_energy(nptcls), mode_vars(ncomp)
        character(len=*), optional, intent(in) :: log_label
        type(projected_latent_estep_part) :: part
        real(dp) :: mode_second(ncomp)
        character(len=:), allocatable :: log_prefix
        integer :: ipart, q
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
        z_postcov = 0.d0
        resid_energy = 0.d0
        resid_mean_energy = 0.d0
        mode_second = 0.d0
        do ipart = 1, nparts
            if( .not. file_exists(part_fnames(ipart)) )then
                THROW_HARD('missing E-step latent part: '//part_fnames(ipart)%to_char())
            endif
            call read_projected_latent_estep_part(part_fnames(ipart), part)
            call reduce_projected_latent_estep_part(part, z, z_postcov, resid_energy, resid_mean_energy, mode_second)
            call kill_projected_latent_estep_part(part)
        end do
        do q = 1, ncomp
            mode_vars(q) = max(MODE_VAR_FLOOR, mode_second(q) / real(max(1,nptcls), dp))
        end do
        call log_seconds(log_prefix//' E-STEP MASTER REDUCE SECONDS', toc(t_total))
    end subroutine reduce_estep_latent_part_files

    subroutine write_projected_latent_estep_part( fname, part )
        class(string), intent(in) :: fname
        type(projected_latent_estep_part), intent(in) :: part
        integer :: funit, io_stat, header(4), nrec
        nrec = part%nrecords
        header = [ESTEP_PART_MAGIC, ESTEP_PART_VERSION, nrec, part%ncomp]
        call fopen(funit, file=fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_projected_latent_estep_part; open '//fname%to_char(), io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_projected_latent_estep_part; header '//fname%to_char(), io_stat)
        if( nrec > 0 )then
            write(funit, iostat=io_stat) part%rows(:nrec), part%pinds(:nrec), part%valid(:nrec)
            call fileiochk('write_projected_latent_estep_part; particle fields '//fname%to_char(), io_stat)
            write(funit, iostat=io_stat) part%zrows(:,:nrec), part%z_postcov(:,:,:nrec)
            call fileiochk('write_projected_latent_estep_part; latent fields '//fname%to_char(), io_stat)
            write(funit, iostat=io_stat) part%resid_energy(:nrec), part%resid_mean_energy(:nrec), part%mode_second(:,:nrec)
            call fileiochk('write_projected_latent_estep_part; residual fields '//fname%to_char(), io_stat)
        endif
        call fclose(funit)
    end subroutine write_projected_latent_estep_part

    subroutine read_projected_latent_estep_part( fname, part )
        class(string), intent(in) :: fname
        type(projected_latent_estep_part), intent(inout) :: part
        integer :: funit, io_stat, header(4), nrec, ncomp
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
        if( nrec < 0 .or. ncomp < 1 ) THROW_HARD('invalid projected latent E-step part header: '//fname%to_char())
        call init_projected_latent_estep_part(part, nrec, ncomp)
        part%nrecords = nrec
        if( nrec > 0 )then
            read(funit, iostat=io_stat) part%rows(:nrec), part%pinds(:nrec), part%valid(:nrec)
            call fileiochk('read_projected_latent_estep_part; particle fields '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) part%zrows(:,:nrec), part%z_postcov(:,:,:nrec)
            call fileiochk('read_projected_latent_estep_part; latent fields '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) part%resid_energy(:nrec), part%resid_mean_energy(:nrec), part%mode_second(:,:nrec)
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
        allocate(part%rows(nrecords_max), part%pinds(nrecords_max), part%valid(nrecords_max), &
            &part%zrows(ncomp,nrecords_max), part%z_postcov(ncomp,ncomp,nrecords_max), &
            &part%resid_energy(nrecords_max), part%resid_mean_energy(nrecords_max), &
            &part%mode_second(ncomp,nrecords_max))
        part%rows              = 0
        part%pinds             = 0
        part%valid             = .false.
        part%zrows             = 0.d0
        part%z_postcov         = 0.d0
        part%resid_energy      = 0.d0
        part%resid_mean_energy = 0.d0
        part%mode_second       = 0.d0
    end subroutine init_projected_latent_estep_part

    subroutine kill_projected_latent_estep_part( part )
        type(projected_latent_estep_part), intent(inout) :: part
        if( allocated(part%rows) ) deallocate(part%rows)
        if( allocated(part%pinds) ) deallocate(part%pinds)
        if( allocated(part%valid) ) deallocate(part%valid)
        if( allocated(part%zrows) ) deallocate(part%zrows)
        if( allocated(part%z_postcov) ) deallocate(part%z_postcov)
        if( allocated(part%resid_energy) ) deallocate(part%resid_energy)
        if( allocated(part%resid_mean_energy) ) deallocate(part%resid_mean_energy)
        if( allocated(part%mode_second) ) deallocate(part%mode_second)
        part%nrecords = 0
        part%ncomp    = 0
    end subroutine kill_projected_latent_estep_part

    subroutine reset_projected_latent_estep_part( part, nrecords )
        type(projected_latent_estep_part), intent(inout) :: part
        integer, intent(in) :: nrecords
        if( .not. allocated(part%valid) ) THROW_HARD('unallocated E-step part')
        if( nrecords > size(part%valid) ) THROW_HARD('E-step part capacity exceeded')
        part%nrecords = nrecords
        part%rows(:nrecords)              = 0
        part%pinds(:nrecords)             = 0
        part%valid(:nrecords)             = .false.
        part%zrows(:,:nrecords)           = 0.d0
        part%z_postcov(:,:,:nrecords)     = 0.d0
        part%resid_energy(:nrecords)      = 0.d0
        part%resid_mean_energy(:nrecords) = 0.d0
        part%mode_second(:,:nrecords)     = 0.d0
    end subroutine reset_projected_latent_estep_part

    subroutine prepare_projected_latent_estep_part( build, mean_rec, basis_recs, fpls_batch, mode_vars, pinds, &
        &batchlims, batchsz, ncomp, part, basis_fpls, mean_fpls, orientations, gram_h, rhs_h, gram, rhs, zrow, post_cov )
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: batchlims(2), batchsz, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        type(fplane_type),   intent(inout) :: fpls_batch(batchsz)
        real(dp),            intent(in)    :: mode_vars(ncomp)
        integer,             intent(in)    :: pinds(:)
        type(projected_latent_estep_part), intent(inout) :: part
        type(fplane_type), intent(inout) :: basis_fpls(ncomp,nthr_glob), mean_fpls(nthr_glob)
        type(ori),         intent(inout) :: orientations(nthr_glob)
        complex(dp),       intent(inout) :: gram_h(ncomp,ncomp,nthr_glob), rhs_h(ncomp,nthr_glob)
        real(dp),          intent(inout) :: gram(ncomp,ncomp,nthr_glob), rhs(ncomp,nthr_glob)
        real(dp),          intent(inout) :: zrow(ncomp,nthr_glob), post_cov(ncomp,ncomp,nthr_glob)
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
                gram(q,q,ithr) = gram(q,q,ithr) + ppca_prior_precision(mode_vars(q))
            end do
            call solve_ppca_posterior(gram(:,:,ithr), rhs(:,ithr), zrow(:,ithr), post_cov(:,:,ithr))
            part%zrows(:,i)         = zrow(:,ithr)
            part%z_postcov(:,:,i)   = post_cov(:,:,ithr)
            do q = 1, ncomp
                part%mode_second(q,i) = zrow(q,ithr) * zrow(q,ithr) + max(0.d0, post_cov(q,q,ithr))
            end do
            do q = 1, ncomp
                call subtract_scaled_plane(fpls_batch(i), basis_fpls(q,ithr), zrow(q,ithr))
            end do
            part%resid_energy(i) = plane_energy(fpls_batch(i))
            part%valid(i) = .true.
        end do
        !$omp end parallel do
    end subroutine prepare_projected_latent_estep_part

    subroutine reduce_projected_latent_estep_part( part, z, z_postcov, resid_energy, resid_mean_energy, mode_second )
        type(projected_latent_estep_part), intent(in) :: part
        real(dp), intent(inout) :: z(:,:), z_postcov(:,:,:)
        real(dp), intent(inout) :: resid_energy(:), resid_mean_energy(:), mode_second(:)
        integer :: i, row
        do i = 1, part%nrecords
            if( .not. part%valid(i) ) cycle
            row = part%rows(i)
            z(row,:)                 = part%zrows(:,i)
            z_postcov(row,:,:)       = part%z_postcov(:,:,i)
            resid_mean_energy(row)   = part%resid_mean_energy(i)
            resid_energy(row)        = part%resid_energy(i)
            mode_second(:)           = mode_second(:) + part%mode_second(:,i)
        end do
    end subroutine reduce_projected_latent_estep_part

    subroutine infer_latents_from_basis( params, build, mean_rec, basis_recs, z, mode_vars, &
        &z_postcov, resid_energy, resid_mean_energy, pinds, nptcls, ncomp, fpls, log_label )
        class(parameters),   intent(in)    :: params
        class(builder),      intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        integer,             intent(in)    :: nptcls, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real(dp),            intent(inout) :: z(nptcls,ncomp)
        real(dp),            intent(inout) :: mode_vars(ncomp)
        real(dp),            intent(out)   :: z_postcov(nptcls,ncomp,ncomp)
        real(dp),            intent(out)   :: resid_energy(nptcls), resid_mean_energy(nptcls)
        integer,             intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        character(len=*), optional, intent(in) :: log_label
        type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpls(:)
        type(ori),         allocatable :: orientations(:)
        type(projected_latent_estep_part) :: estep_part
        complex(dp), allocatable :: gram_h(:,:,:), rhs_h(:,:)
        real(dp),    allocatable :: gram(:,:,:), rhs(:,:), zrow(:,:), post_cov(:,:,:), mode_second(:)
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
            &rhs(ncomp,nthr_glob), zrow(ncomp,nthr_glob), post_cov(ncomp,ncomp,nthr_glob), mode_second(ncomp))
        resid_energy = 0.d0
        resid_mean_energy = 0.d0
        z_postcov = 0.d0
        mode_second = 0.d0
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
                call prepare_projected_latent_estep_part(build, mean_rec, basis_recs, fpls(:batchsz), mode_vars, pinds, &
                    &batchlims, batchsz, ncomp, estep_part, basis_fpls, mean_fpls, orientations, &
                    &gram_h, rhs_h, gram, rhs, zrow, post_cov)
                call reduce_projected_latent_estep_part(estep_part, z, z_postcov, resid_energy, resid_mean_energy, mode_second)
                if( batchlims(2) == nptcls .or. mod(batchlims(2), progress_stride) == 0 )then
                    write(logfhandle,'(A,I0,A,I0)') log_prefix//' E-STEP PARTICLES: ', batchlims(2), ' / ', nptcls
                    call flush(logfhandle)
                endif
            end do
        end do
        call log_seconds(log_prefix//' E-STEP INFERENCE SECONDS', toc(t_phase))
        do q = 1, ncomp
            mode_vars(q) = max(MODE_VAR_FLOOR, mode_second(q) / real(max(1,nptcls), dp))
        end do
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
        deallocate(basis_fpls, mean_fpls, orientations, gram_h, rhs_h, gram, rhs, zrow, post_cov, mode_second)
        call log_seconds(log_prefix//' E-STEP TOTAL SECONDS', toc(t_total))
    end subroutine infer_latents_from_basis

    subroutine solve_coupled_basis_exp( basis_recs, rho_cross_exp, ncomp )
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real,                intent(in)    :: rho_cross_exp(:,:,:,:)
        complex(dp) :: rhs(ncomp), sol(ncomp)
        real(dp)    :: amat(ncomp,ncomp)
        real(dp)    :: diag_sum, ridge, denom
        integer     :: lb(3), ub(3), h, k, m, ih, ik, im, q, r, flag
        lb = lbound(basis_recs(1)%cmat_exp)
        ub = ubound(basis_recs(1)%cmat_exp)
        !$omp parallel do collapse(3) default(shared) schedule(static) &
        !$omp private(h,k,m,ih,ik,im,q,r,amat,rhs,sol,diag_sum,ridge,denom,flag) proc_bind(close)
        do m = lb(3), ub(3)
            do k = lb(2), ub(2)
                do h = lb(1), ub(1)
                    ih = h - lb(1) + 1
                    ik = k - lb(2) + 1
                    im = m - lb(3) + 1
                    amat = 0.d0
                    rhs  = DCMPLX_ZERO
                    diag_sum = 0.d0
                    do q = 1, ncomp
                        rhs(q) = cmplx(basis_recs(q)%cmat_exp(h,k,m), kind=dp)
                        do r = q, ncomp
                            amat(q,r) = real(rho_cross_exp(pair_index(q,r),ih,ik,im), dp)
                            amat(r,q) = amat(q,r)
                        end do
                        diag_sum = diag_sum + max(0.d0, amat(q,q))
                    end do
                    if( diag_sum <= DTINY .and. sum(abs(rhs)) <= DTINY )then
                        do q = 1, ncomp
                            basis_recs(q)%cmat_exp(h,k,m) = CMPLX_ZERO
                        end do
                        cycle
                    endif
                    ridge = max(COUPLED_MSTEP_RIDGE_ABS, COUPLED_MSTEP_RIDGE_REL * diag_sum / real(max(1,ncomp), dp))
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

    subroutine finalize_basis_for_projection( params, basis_rec, density_corrected )
        class(parameters),   intent(in)    :: params
        type(reconstructor), intent(inout) :: basis_rec
        logical, optional,   intent(in)    :: density_corrected
        logical :: l_density_corrected
        l_density_corrected = .false.
        if( present(density_corrected) ) l_density_corrected = density_corrected
        call basis_rec%compress_exp
        if( .not. l_density_corrected ) call basis_rec%sampl_dens_correct
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
            if( params%l_ml_reg .and. allocated(build%esig%sigma2_noise) )then
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr), &
                    &shift, fplanes(i), build%esig%sigma2_noise(kfromto(1):kfromto(2),iptcl), store_transfer=.true.)
            else
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr), &
                    &shift, fplanes(i), store_transfer=.true.)
            endif
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

    real(dp) function ppca_prior_precision( mode_var ) result( prec )
        real(dp), intent(in) :: mode_var
        prec = 1.d0 / max(MODE_VAR_FLOOR, mode_var)
        prec = max(prec, LATENT_RIDGE)
    end function ppca_prior_precision

    subroutine solve_ppca_posterior( gram, rhs, x, post_cov )
        real(dp), intent(in)  :: gram(:,:), rhs(:)
        real(dp), intent(out) :: x(:), post_cov(:,:)
        integer :: flag
        integer  :: n
        n = size(rhs)
        x = 0.d0
        post_cov = 0.d0
        call hermitian_solve(gram, rhs, x, flag)
        if( flag /= 0 ) return
        call hermitian_invert(gram, post_cov, flag)
        if( flag /= 0 ) post_cov = 0.d0
    end subroutine solve_ppca_posterior

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

end module simple_projected_latent_model
