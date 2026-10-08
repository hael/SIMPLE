!@descr: flex_pca projection-aware latent model: Fourier projection/backprojection helpers, particle prep, the coupled M-step solve
module simple_flex_reconstructor_latent_ops
!$ use omp_lib, only: omp_get_thread_num
use simple_flex_pca_planes,      only: flex_plane_store
use simple_core_module_api, only: cmplx_zero, dcmplx_zero, dp, dtiny, fdim, fplane_type, kbwinsz, ori, &
    &simple_exception, sp, tic, timer_int_kind, tiny, toc
use simple_reconstructor,        only: reconstructor, exp_samples
use simple_builder,              only: builder
use simple_image,                only: image
use simple_linalg,               only: solve_real_spd_complex
use simple_matcher_ptcl_io,      only: discrete_read_imgbatch
use simple_matcher_3Drec,        only: prep_imgs4rec, gen_rec_plane
use simple_memoize_ft_maps,      only: memoize_ft_maps
use simple_parameters,           only: parameters
implicit none

public :: project_fplanes_mean_basis, LATENT_WDIM
!> the projection-aware latent model (merged from simple_flex_projected_latent_model)
public :: planes_batch_load, prep_imgs4projected_model, solve_coupled_basis_exp, projected_model_kfromto
public :: add_invtausq2rho_coupled, pair_index
private
#include "simple_local_flags.inc"

!> the KB window width of the particle-plane samplers (the polar bank's 2D gather)
integer, parameter :: LATENT_WDIM = 2 * ceiling(KBWINSZ - 0.5) + 1

real(dp), parameter :: COUPLED_MSTEP_RIDGE_REL = 1.0d-8
real(dp), parameter :: COUPLED_DENSITY_FLOOR   = 1.0d-6

contains

    !> The mean's and the basis volumes' central sections at the samples of plane fpl_ref (orientation o),
    !! times the forward transfer with apply_ctf_amp: one sample set (window geometry once per sample),
    !! gathered from every volume. All volumes share the mean's expanded lattice.
    subroutine project_fplanes_mean_basis( mean_rec, basis_recs, o, fpl_ref, mean_fpl, basis_fpls, apply_ctf_amp )
        type(reconstructor), intent(in)    :: mean_rec
        type(reconstructor), intent(in)    :: basis_recs(:)
        class(ori),          intent(inout) :: o
        class(fplane_type),  intent(in)    :: fpl_ref
        type(fplane_type),   intent(inout) :: mean_fpl
        type(fplane_type),   intent(inout) :: basis_fpls(:)
        logical,             intent(in)    :: apply_ctf_amp
        type(exp_samples)    :: samples
        complex, allocatable :: vals(:)
        integer :: q, ncomp
        if( .not. allocated(mean_rec%cmat_exp) )then
            THROW_HARD('expanded mean matrix does not exist; project_fplanes_mean_basis')
        endif
        if( .not. allocated(fpl_ref%cmplx_plane) )then
            THROW_HARD('reference Fourier plane does not exist; project_fplanes_mean_basis')
        endif
        ncomp = size(basis_recs)
        if( size(basis_fpls) < ncomp )then
            THROW_HARD('basis output plane array too small; project_fplanes_mean_basis')
        endif
        do q = 1, ncomp
            if( .not. allocated(basis_recs(q)%cmat_exp) )then
                THROW_HARD('expanded basis matrix does not exist; project_fplanes_mean_basis')
            endif
        end do
        call samples%new_plane(mean_rec, o, fpl_ref)
        allocate(vals(samples%get_n()))
        call samples%gather(mean_rec, vals)
        call samples%put_plane(vals, fpl_ref, mean_fpl, apply_ctf_amp)
        do q = 1, ncomp
            call samples%gather(basis_recs(q), vals)
            call samples%put_plane(vals, fpl_ref, basis_fpls(q), apply_ctf_amp)
        end do
        deallocate(vals)
        call samples%kill
    end subroutine project_fplanes_mean_basis

    subroutine solve_coupled_basis_exp( basis_recs, rho_cross_exp, ncomp )
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real,                intent(in)    :: rho_cross_exp(:,:,:,:)
        complex(dp) :: rhs(ncomp), sol(ncomp)
        real(dp)    :: amat(ncomp,ncomp)
        real(dp)    :: diag_sum, diag_max, ridge, denom
        integer     :: lb(3), ub(3), h, k, m, ih, ik, im, q, r, flag, shell, nyq
        logical     :: diagonal_density
        ! Same shape convention insert_planes_multi (simple_reconstructor) uses to pick its
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

    !> Flex analog of reconstructor::add_invtausq2rho: adds the per-component, per-shell invtau2(q,sh)
    !! (crossfsc_to_invtau2) to the diagonal rows pair_index(q,q) of the packed coupled density.
    !! Call after any distributed reduction and before solve_coupled_basis_exp.
    subroutine add_invtausq2rho_coupled( basis_recs, rho_cross_exp, ncomp, invtau2 )
        integer,             intent(in)    :: ncomp
        type(reconstructor), intent(in)    :: basis_recs(ncomp)  !< lattice geometry reference only
        real,                intent(inout) :: rho_cross_exp(:,:,:,:)
        real,                intent(in)    :: invtau2(:,:)       !< (ncomp, nshells) per component, per shell
        integer :: lb(3), ub(3), h, k, m, ih, ik, im, q, shell, nyq, shmax
        if( size(rho_cross_exp,1) /= (ncomp*(ncomp+1))/2 ) &
            &THROW_HARD('add_invtausq2rho_coupled requires the FULL packed coupled density')
        if( size(invtau2,1) < ncomp ) THROW_HARD('add_invtausq2rho_coupled: invtau2 rank mismatch')
        lb    = lbound(basis_recs(1)%cmat_exp)
        ub    = ubound(basis_recs(1)%cmat_exp)
        nyq   = basis_recs(1)%get_lfny(1)
        shmax = min(nyq, size(invtau2,2))
        !$omp parallel do collapse(3) default(shared) schedule(static) &
        !$omp private(h,k,m,ih,ik,im,q,shell) proc_bind(close)
        do m = lb(3), ub(3)
            do k = lb(2), ub(2)
                do h = lb(1), ub(1)
                    ih = h - lb(1) + 1
                    ik = k - lb(2) + 1
                    im = m - lb(3) + 1
                    ! same shell convention and spherical-Nyquist support rule as the solve
                    shell = nint(sqrt(real(h*h+k*k+m*m)))
                    if( shell < 1 .or. shell > shmax ) cycle
                    do q = 1, ncomp
                        rho_cross_exp(pair_index(q,q),ih,ik,im) = &
                            &rho_cross_exp(pair_index(q,q),ih,ik,im) + invtau2(q,shell)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine add_invtausq2rho_coupled

    !> Serve a batch from the resident store when possible; otherwise read, prepare and retain it.
    subroutine planes_batch_load( plane_store, params, build, n, pinds, batchlims, fpls, cached, sec_read, sec_prep )
        class(flex_plane_store),        intent(inout) :: plane_store
        class(parameters),              intent(in)    :: params
        class(builder),                 intent(inout) :: build
        integer,                        intent(in)    :: n, pinds(n), batchlims(2)
        type(fplane_type),              intent(inout) :: fpls(:)
        logical,                        intent(in)    :: cached
        real(timer_int_kind), optional, intent(inout) :: sec_read, sec_prep
        integer(timer_int_kind) :: t
        integer :: batchsz
        logical :: found
        batchsz = batchlims(2) - batchlims(1) + 1
        call plane_store%fetch(pinds(batchlims(1):batchlims(2)), fpls(:batchsz), cached, found, sec_read)
        if( found ) return
        t = tic()
        if( cached )then
            call plane_store%read_cache_batch(params, n, pinds, batchlims)
        else
            call discrete_read_imgbatch(params, build, n, pinds, batchlims)
        endif
        if( present(sec_read) ) sec_read = sec_read + toc(t)
        t = tic()
        call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
            &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), cached=cached, plane_store=plane_store)
        if( present(sec_prep) ) sec_prep = sec_prep + toc(t)
        call plane_store%store(pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
    end subroutine planes_batch_load

    !> The projected model's particle planes, prepared as the reconstruction prepares them (prep_imgs4rec):
    !! on the projected model's band, as observation-model planes (whitened observation and forward
    !! transfer), whitened when ML regularization is on. cached: the batch comes from the plane cache, which
    !! holds the full-box preparation's padded transform on the box_crop grid, so only the plane generation
    !! (gen_rec_plane) runs, on the cropped grid.
    subroutine prep_imgs4projected_model( params, build, nptcls, ptcl_imgs, pinds, fplanes, cached, plane_store )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        class(image),      intent(inout) :: ptcl_imgs(nptcls)
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        logical, optional, intent(in)    :: cached      !< serve reads from the downscaled cache
        class(flex_plane_store), optional, intent(in) :: plane_store
        integer :: i, ithr, kfromto(2)
        logical :: l_cached
        l_cached = .false.
        if( present(cached) ) l_cached = cached
        kfromto = projected_model_kfromto(params%box_crop, params%smpd_crop, params%lp)
        if( .not. l_cached )then
            call prep_imgs4rec(params, build, nptcls, ptcl_imgs, pinds, fplanes, kfromto=kfromto, &
                &observation_model=.true., whiten=params%l_ml_reg)
            return
        endif
        if( .not. present(plane_store) ) THROW_HARD('cached particle preparation requested without its plane store')
        if( params%l_ml_reg .and. .not. allocated(build%esig%sigma2_noise) )then
            THROW_HARD('projected covariance model requested whitening without loaded sigma2 spectra')
        endif
        ! a cached particle already lives on the cropped grid, so the pad heap and the map are box_croppd
        call memoize_ft_maps([params%box_croppd, params%box_croppd, 1], params%smpd_crop)
        kfromto(2) = min(kfromto(2), params%box_crop/2)
        !$omp parallel do default(shared) private(i,ithr) schedule(static) proc_bind(close)
        do i = 1, nptcls
            ithr = omp_get_thread_num() + 1
            call plane_store%fill_cached_image(i, build%img_pad_heap(ithr))
            call gen_rec_plane(params, build, pinds(i), build%img_pad_heap(ithr), .true., kfromto, &
                &params%l_ml_reg, .true., .true., fplanes(i))
        end do
        !$omp end parallel do
    end subroutine prep_imgs4projected_model

    function projected_model_kfromto( box_crop, smpd_crop, lp ) result( kfromto )
        integer, intent(in) :: box_crop
        real,    intent(in) :: smpd_crop, lp
        integer :: kfromto(2), kto_full
        real    :: dstep_crop
        kto_full = max(1, fdim(box_crop) - 1)
        kfromto(1) = 1
        kfromto(2) = kto_full
        if( lp > 2.0 * smpd_crop + TINY )then
            dstep_crop = real(max(1, box_crop - 1)) * smpd_crop
            kfromto(2) = max(1, min(kto_full, int(dstep_crop / lp)))
        endif
    end function projected_model_kfromto

end module simple_flex_reconstructor_latent_ops
