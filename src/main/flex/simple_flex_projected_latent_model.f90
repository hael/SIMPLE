!@descr: Projection-aware latent volume model kernels for flex_analysis
module simple_flex_projected_latent_model
use simple_core_module_api
use simple_builder,          only: builder
use simple_image,            only: image
use simple_linalg,           only: eigsrt, jacobi
use simple_memoize_ft_maps,  only: memoize_ft_maps
use simple_parameters,       only: parameters
use simple_reconstructor,    only: reconstructor
implicit none

public :: prep_imgs4projected_model, solve_coupled_basis_exp, projected_model_kfromto
public :: test_projected_latent_mstep_stats_io, test_projected_latent_canonicalization
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
! Version 2 stores one ordinary reconstruction density shared by all
! independently reconstructed residual modes.
integer,  parameter :: MSTEP_STATS_VERSION = 2
integer,  parameter :: ESTEP_PART_MAGIC = 1180053580
integer,  parameter :: ESTEP_PART_VERSION = 3
integer(longer), parameter :: MSTEP_STATS_IO_TARGET_BYTES = 64_longer * 1024_longer * 1024_longer

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
        stats%npairs    = 1
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
        npairs = 1
        if( header(5) /= npairs ) THROW_HARD('M-step statistics density-count mismatch: '//fname%to_char())
        if( size(rho_cross_exp,1) /= npairs ) THROW_HARD('M-step reduction requires one shared density volume')
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
        integer :: exp_lb(3), exp_ub(3), exp_shape(3), q, i, j, k
        real :: rhs_err, rho_err
        exp_lb    = [-3, -TEST_BOX/2-2, -TEST_BOX/2-2]
        exp_ub    = [ TEST_BOX/2+2, TEST_BOX/2+2, TEST_BOX/2+2]
        exp_shape = exp_ub - exp_lb + 1
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
        stats%rho_cross(1,:,:,:) = 0.01
        allocate(rho_cross(1,exp_shape(1),exp_shape(2),exp_shape(3)), source=0.)
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
        logical     :: diagonal_density
        ! Same shape convention insert_planes_oversamp_coupled_batch_scaled uses to pick its
        ! accumulation mode: a leading extent of ncomp means only the diagonal of the coupled normal
        ! matrix was accumulated, so the per-voxel system decouples into ncomp scalar divisions.
        diagonal_density = size(rho_cross_exp,1) == ncomp .and. ncomp /= (ncomp*(ncomp+1))/2
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
                        if( diagonal_density )then
                            denom = max(0.d0,real(rho_cross_exp(q,ih,ik,im),dp))
                        else
                            denom = max(0.d0,real(rho_cross_exp(pair_index(q,q),ih,ik,im),dp))
                        endif
                        diag_sum = diag_sum + denom
                        diag_max = max(diag_max,denom)
                    end do
                    if( diag_max <= COUPLED_DENSITY_FLOOR )then
                        do q = 1, ncomp
                            basis_recs(q)%cmat_exp(h,k,m) = CMPLX_ZERO
                        end do
                        cycle
                    endif
                    ridge = COUPLED_MSTEP_RIDGE_REL * diag_sum / real(max(1,ncomp), dp)
                    if( diagonal_density )then
                        ! No off-diagonals were accumulated, so the normal matrix IS its diagonal and
                        ! the ncomp x ncomp Cholesky collapses to ncomp divisions. Same ridge, same
                        ! floor, same fallback -- only the cross terms between components are dropped.
                        do q = 1, ncomp
                            denom = max(0.d0,real(rho_cross_exp(q,ih,ik,im),dp)) + ridge
                            if( denom > DTINY )then
                                sol(q) = rhs(q) / denom
                            else
                                sol(q) = DCMPLX_ZERO
                            endif
                        end do
                        do q = 1, ncomp
                            basis_recs(q)%cmat_exp(h,k,m) = cmplx(real(sol(q), sp), real(aimag(sol(q)), sp))
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

    !!  mskrad (optional, pixels at params%box): when present the particle is soft-masked to
    !!  that radius after noise normalization instead of edge-tapered. This is SIMPLE's
    !!  equivalent of the reference's mask_images_in_H_B/mask_images_in_proj (covariance_estimation
    !!  options, both default True), which masks each image to the projected molecular
    !!  envelope before the covariance accumulation. Solvent outside the particle contributes
    !!  only noise, and at box_crop=64/mskdiam=200 the disc keeps ~21% of the frame, so the
    !!  noise in every per-image inner product drops by roughly the same factor. Left absent
    !!  the behaviour is exactly as before.
    subroutine prep_imgs4projected_model( params, build, nptcls, ptcl_imgs, pinds, fplanes, &
        &mskrad, force_cpu, resident, cached )
        use simple_flex_gpu, only: flex_gpu_prep_ready
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        class(image),      intent(inout) :: ptcl_imgs(nptcls)
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        real, optional,    intent(in)    :: mskrad
        logical, optional, intent(in)    :: force_cpu   !< cross-check reference building
        logical, optional, intent(in)    :: resident    !< device path: leave planes resident only
        logical, optional, intent(in)    :: cached      !< serve reads from the downscaled cache
        type(ctfparams) :: ctfparms(nthr_glob)
        real    :: shift(2), crop_factor
        integer :: iptcl, i, ithr, kfromto(2)
        logical :: l_mask, l_cpu, l_res, l_cached
        l_mask = .false.
        if( present(mskrad) ) l_mask = mskrad > 0.0
        l_cpu = .false.
        if( present(force_cpu) ) l_cpu = force_cpu
        l_res = .false.
        if( present(resident) ) l_res = resident
        l_cached = .false.
        if( present(cached) ) l_cached = cached
        ! A cache entry is the noise-normalised, Fourier-cropped particle at box_crop. That prefix
        ! is equivalent to the full-box path only for the TAPER variant, which is the one
        ! prep_imgs4rec certified: norm_noise_mask_pad_fft has no renorm= switch, so a masked run
        ! would noise-normalise a second time. Refuse rather than change the numerics silently.
        if( l_cached .and. l_mask ) THROW_HARD('particle cache is incompatible with image masking &
            &(SIMPLE_COV_MASK_IMAGES); prep_imgs4projected_model')
        ! device path: when a stage driver has begun the GPU prep lifecycle, the whole
        ! taper->norm->pad->FFT->plane chain runs on device and the planes are fetched packed
        ! (taper variant only; the mask variant stays on the CPU)
        if( flex_gpu_prep_ready() .and. .not. l_mask .and. .not. l_cpu .and. .not. l_cached )then
            call prep_imgs4projected_model_dev(params, build, nptcls, ptcl_imgs, pinds, &
                &fplanes, fetch=.not. l_res)
            return
        endif
        if( l_res ) THROW_HARD('resident prep requested without the device prep lifecycle')
        ! logical/physical address mapping for padded Fourier planes: a cached particle already
        ! lives on the cropped grid, so the pad heap and the map must both be box_croppd
        if( l_cached )then
            call memoize_ft_maps([params%box_croppd, params%box_croppd, 1], params%smpd_crop)
        else
            call memoize_ft_maps([params%boxpd, params%boxpd, 1], params%smpd)
        endif
        kfromto = projected_model_kfromto(params)
        if( l_cached ) kfromto(2) = min(kfromto(2), params%box_crop/2)
        crop_factor = real(params%box_crop) / real(params%box)
        !$omp parallel do default(shared) private(i,ithr,iptcl,shift) schedule(static) proc_bind(close)
        do i = 1, nptcls
            ithr   = omp_get_thread_num() + 1
            iptcl  = pinds(i)
            if( l_mask )then
                call ptcl_imgs(i)%norm_noise_mask_pad_fft(build%lmsk, mskrad, build%img_pad_heap(ithr))
            else if( l_cached )then
                ! already noise-normalised when the entry was written; renorm=.false. stops the
                ! cropped particle being rescaled a second time (see simple_ptcl_cache)
                call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft(build%lmsk_crop, &
                    &build%img_pad_heap(ithr), renorm=.false.)
            else
                call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft(build%lmsk, build%img_pad_heap(ithr))
            endif
            ctfparms(ithr) = build%spproj%get_ctfparams(params%oritype, iptcl)
            shift = build%spproj_field%get_2Dshift(iptcl)
            if( l_cached )then
                ! shconst is in pixels of the padded box the image actually has, and the CTF kernel
                ! reads cycles/pixel of the current grid -- both must move to the cropped grid
                ctfparms(ithr)%smpd = ctfparms(ithr)%smpd / crop_factor   ! = smpd_crop
                shift               = shift * crop_factor
            endif
            if( params%l_ml_reg )then
                if( .not. allocated(build%esig%sigma2_noise) )then
                    THROW_HARD('projected covariance model requested whitening without loaded sigma2 spectra')
                endif
                if( iptcl < lbound(build%esig%sigma2_noise,2) .or. &
                    &iptcl > ubound(build%esig%sigma2_noise,2) )then
                    THROW_HARD('projected covariance particle index is outside the sigma2 table')
                endif
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr), &
                    &shift, fplanes(i), build%esig%sigma2_noise(kfromto(1):kfromto(2),iptcl), &
                    &store_transfer=.true., observation_model=.true.)
            else
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr), &
                    &shift, fplanes(i), store_transfer=.true., observation_model=.true.)
            endif
            call cap_fplane_for_projected_model(fplanes(i), kfromto)
        end do
        !$omp end parallel do
    end subroutine prep_imgs4projected_model

    !> device variant of the stage prep: gathers the per-record host scalars (CTF, shift,
    !! sigma2), runs the taper->norm->pad->FFT->plane chain on the GPU, fetches the packed
    !! separated observation-model planes and unpacks them into fplane_type. Byte-equivalent
    !! to the CPU generator: values live on the OSMPL_PAD_FAC-multiple lattice, zeros
    !! elsewhere; plane reuse across batches skips the re-zero of the never-written gaps.
    subroutine prep_imgs4projected_model_dev( params, build, nptcls, ptcl_imgs, pinds, &
        &fplanes, fetch )
        use simple_flex_gpu, only: flex_gpu_prep_batch_f, flex_gpu_prep_fetch_sep_f
        use simple_ftiter,   only: ftiter
        use simple_math,     only: ceil_div, floor_div
        use simple_math_ft,  only: resample_sigma2
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        class(image),      intent(inout) :: ptcl_imgs(nptcls)
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        logical, optional, intent(in)    :: fetch   !< .false. = planes stay resident only
        type(ctfparams) :: ctfp_arr(nptcls)
        real            :: shf(2,nptcls)
        logical         :: vld(nptcls)
        complex(sp), allocatable :: plcy(:,:,:)
        real(sp),    allocatable :: plt(:,:,:), plct(:,:,:)
        real,        allocatable :: sig2_ups(:,:)
        type(ftiter) :: fit_pd, fit_cr
        real    :: shconst_pd(3)
        integer :: kfromto(2), frlims_pd(3,2), i, h, k, hlo, hhi, klo, nyqpd, signyq
        integer :: hmin, hmax, kmin, iptcl, pf
        logical :: l_fresh, l_fetch
        l_fetch = .true.
        if( present(fetch) ) l_fetch = fetch
        pf      = OSMPL_PAD_FAC
        kfromto = projected_model_kfromto(params)
        call fit_pd%new([params%boxpd, params%boxpd, 1], params%smpd_crop)
        frlims_pd = fit_pd%loop_lims(3)
        ! fill to the FULL padded band like the CPU generator; the working-band cap acts
        ! through fpl%nyq downstream, not through the stored values
        nyqpd     = fit_pd%get_lfny(1)
        if( params%l_ml_reg )then
            if( .not. allocated(build%esig%sigma2_noise) )then
                THROW_HARD('projected covariance model requested whitening without loaded sigma2 spectra')
            endif
            call fit_cr%new([params%box_crop, params%box_crop, 1], params%smpd_crop)
            signyq = fit_cr%get_lfny(1)
            allocate(sig2_ups(0:nyqpd, nptcls), source=1.0)
        endif
        vld = .true.
        !$omp parallel do default(shared) private(i,iptcl) schedule(static) proc_bind(close)
        do i = 1, nptcls
            iptcl       = pinds(i)
            ctfp_arr(i) = build%spproj%get_ctfparams(params%oritype, iptcl)
            shf(:,i)    = build%spproj_field%get_2Dshift(iptcl)
            if( params%l_ml_reg )then
                if( iptcl < lbound(build%esig%sigma2_noise,2) .or. &
                    &iptcl > ubound(build%esig%sigma2_noise,2) )then
                    THROW_HARD('projected covariance particle index is outside the sigma2 table')
                endif
                call resample_sigma2(kfromto(1), signyq, &
                    &build%esig%sigma2_noise(kfromto(1):kfromto(2), iptcl), nyqpd, &
                    &real(signyq)/real(nyqpd), sig2_ups(:,i))
            endif
        end do
        !$omp end parallel do
        if( params%l_ml_reg )then
            call flex_gpu_prep_batch_f(ptcl_imgs, ctfp_arr, shf, vld, nptcls, params%box, &
                &frlims_pd, nyqpd, sig2_ups=sig2_ups)
        else
            call flex_gpu_prep_batch_f(ptcl_imgs, ctfp_arr, shf, vld, nptcls, params%box, &
                &frlims_pd, nyqpd)
        endif
        if( .not. l_fetch )then
            ! consumer is a resident device stage (fused columns); planes stay on device
            if( allocated(sig2_ups) ) deallocate(sig2_ups)
            return
        endif
        hlo = ceil_div (frlims_pd(1,1), pf); hhi = floor_div(frlims_pd(1,2), pf)
        klo = ceil_div (frlims_pd(2,1), pf)
        allocate(plcy(hlo:hhi, klo:0, nptcls), plt(hlo:hhi, klo:0, nptcls), &
            &plct(hlo:hhi, klo:0, nptcls))
        call flex_gpu_prep_fetch_sep_f(plcy, plt, plct, nptcls, hlo, hhi, klo)
        hmin = frlims_pd(1,1); hmax = frlims_pd(1,2); kmin = frlims_pd(2,1)
        shconst_pd      = 0.
        shconst_pd(1:2) = PI/real(params%boxpd/2)
        !$omp parallel do default(shared) private(i,h,k,l_fresh) schedule(static) proc_bind(close)
        do i = 1, nptcls
            l_fresh = .not. allocated(fplanes(i)%cmplx_plane)
            if( .not. l_fresh )then
                if( lbound(fplanes(i)%cmplx_plane,1) /= hmin .or. &
                    &ubound(fplanes(i)%cmplx_plane,1) /= hmax .or. &
                    &lbound(fplanes(i)%cmplx_plane,2) /= kmin .or. &
                    &ubound(fplanes(i)%cmplx_plane,2) /= 0    .or. &
                    &.not. allocated(fplanes(i)%transfer_plane) )then
                    call cleanup_plane(fplanes(i))
                    l_fresh = .true.
                endif
            endif
            if( l_fresh )then
                allocate(fplanes(i)%cmplx_plane(hmin:hmax, kmin:0),    source=cmplx(0.,0.))
                allocate(fplanes(i)%ctfsq_plane(hmin:hmax, kmin:0),    source=0.)
                allocate(fplanes(i)%transfer_plane(hmin:hmax, kmin:0), source=cmplx(0.,0.))
            endif
            do k = klo, 0
                do h = hlo, hhi
                    fplanes(i)%cmplx_plane(pf*h, pf*k)    = plcy(h,k,i)
                    fplanes(i)%transfer_plane(pf*h, pf*k) = cmplx(plt(h,k,i), 0.)
                    fplanes(i)%ctfsq_plane(pf*h, pf*k)    = plct(h,k,i)
                end do
            end do
            fplanes(i)%frlims  = frlims_pd
            fplanes(i)%nyq     = nyqpd
            fplanes(i)%shconst = shconst_pd
            call cap_fplane_for_projected_model(fplanes(i), kfromto)
        end do
        !$omp end parallel do
        deallocate(plcy, plt, plct)
        if( allocated(sig2_ups) ) deallocate(sig2_ups)
    end subroutine prep_imgs4projected_model_dev

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

    subroutine cleanup_plane( fpl )
        type(fplane_type), intent(inout) :: fpl
        if( allocated(fpl%cmplx_plane) ) deallocate(fpl%cmplx_plane)
        if( allocated(fpl%ctfsq_plane) ) deallocate(fpl%ctfsq_plane)
        if( allocated(fpl%transfer_plane) ) deallocate(fpl%transfer_plane)
        fpl%frlims  = 0
        fpl%shconst = 0.
        fpl%nyq     = 0
    end subroutine cleanup_plane

end module simple_flex_projected_latent_model
