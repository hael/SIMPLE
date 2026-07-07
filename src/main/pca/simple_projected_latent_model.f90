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
use simple_reconstructor_latent_ops, only: insert_plane_oversamp_coupled_scaled, project_fplanes_mean_basis
implicit none

public :: update_basis_from_latents, infer_latents_from_basis
public :: initialize_latents, orthonormalize_latents, latent_sdev, latent_covariance
public :: basis_fourier_energy, cleanup_planes, projected_model_kfromto
private
#include "simple_local_flags.inc"

real(dp), parameter :: LATENT_RIDGE = 1.0d-3
real(dp), parameter :: MODE_VAR_FLOOR = 1.0d-3
real(dp), parameter :: COUPLED_MSTEP_RIDGE_REL = 1.0d-8
real(dp), parameter :: COUPLED_MSTEP_RIDGE_ABS = 1.0d-10

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
        type(ori)         :: orientation
        real,    allocatable :: rho_cross_exp(:,:,:,:)
        real(dp)             :: latent_second(ncomp,ncomp)
        character(len=:), allocatable :: log_prefix
        integer              :: exp_shape(3), npairs
        integer           :: batchlims(2), batchsz, ibatch, i, iptcl, q, r, row, progress_stride
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
        do q = 1, ncomp
            call basis_recs(q)%reset
            call basis_recs(q)%reset_exp
        end do
        npairs    = (ncomp * (ncomp + 1)) / 2
        exp_shape = shape(basis_recs(1)%cmat_exp)
        allocate(rho_cross_exp(npairs, exp_shape(1), exp_shape(2), exp_shape(3)), source=0.)
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, init_volumes=.false.)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        t_phase = tic()
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call read_particles(params, build, nptcls, pinds, batchlims, batchsz)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            do i = 1, batchsz
                row   = batchlims(1) + i - 1
                iptcl = pinds(row)
                call build%spproj_field%get_ori(iptcl, orientation)
                if( orientation%isstatezero() ) cycle
                call mean_rec%project_fplane(orientation, fpls(i), mean_fpl, apply_ctf_amp=.true.)
                call subtract_plane(fpls(i), mean_fpl)
                latent_second = z_postcov(row,:,:)
                do q = 1, ncomp
                    do r = 1, ncomp
                        latent_second(q,r) = latent_second(q,r) + z(row,q) * z(row,r)
                    end do
                end do
                call insert_plane_oversamp_coupled_scaled(basis_recs, rho_cross_exp, build%pgrpsyms, &
                    &orientation, fpls(i), z(row,:), latent_second)
            end do
            if( batchlims(2) == nptcls .or. mod(batchlims(2), progress_stride) == 0 )then
                write(logfhandle,'(A,I0,A,I0)') log_prefix//' M-STEP PARTICLES: ', batchlims(2), ' / ', nptcls
                call flush(logfhandle)
            endif
        end do
        call log_seconds(log_prefix//' M-STEP INSERT SECONDS', toc(t_phase))
        call orientation%kill
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
        complex(dp), allocatable :: gram_h(:,:,:), rhs_h(:,:)
        real(dp),    allocatable :: gram(:,:,:), rhs(:,:), zrow(:,:), post_cov(:,:,:), mode_second(:,:)
        character(len=:), allocatable :: log_prefix
        integer           :: batchlims(2), batchsz, ibatch, i, iptcl, q, r, row, ithr, progress_stride
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
        allocate(basis_fpls(ncomp,nthr_glob), mean_fpls(nthr_glob), orientations(nthr_glob), &
            &gram_h(ncomp,ncomp,nthr_glob), rhs_h(ncomp,nthr_glob), gram(ncomp,ncomp,nthr_glob), &
            &rhs(ncomp,nthr_glob), zrow(ncomp,nthr_glob), post_cov(ncomp,ncomp,nthr_glob), mode_second(ncomp,nthr_glob))
        resid_energy = 0.d0
        resid_mean_energy = 0.d0
        z_postcov = 0.d0
        mode_second = 0.d0
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, init_volumes=.false.)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        t_phase = tic()
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call read_particles(params, build, nptcls, pinds, batchlims, batchsz)
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
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
                gram_h(:,:,ithr) = DCMPLX_ZERO
                rhs_h(:,ithr)    = DCMPLX_ZERO
                gram(:,:,ithr)   = 0.d0
                rhs(:,ithr)      = 0.d0
                do q = 1, ncomp
                    rhs_h(q,ithr) = hermitian_plane_inner_product(basis_fpls(q,ithr), fpls(i))
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
                z(row,:) = zrow(:,ithr)
                z_postcov(row,:,:) = post_cov(:,:,ithr)
                do q = 1, ncomp
                    mode_second(q,ithr) = mode_second(q,ithr) + zrow(q,ithr) * zrow(q,ithr) + &
                        &max(0.d0, post_cov(q,q,ithr))
                end do
                do q = 1, ncomp
                    call subtract_scaled_plane(fpls(i), basis_fpls(q,ithr), zrow(q,ithr))
                end do
                resid_energy(row) = plane_energy(fpls(i))
            end do
            !$omp end parallel do
            if( batchlims(2) == nptcls .or. mod(batchlims(2), progress_stride) == 0 )then
                write(logfhandle,'(A,I0,A,I0)') log_prefix//' E-STEP PARTICLES: ', batchlims(2), ' / ', nptcls
                call flush(logfhandle)
            endif
        end do
        call log_seconds(log_prefix//' E-STEP INFERENCE SECONDS', toc(t_phase))
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
