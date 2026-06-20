!@descr: CSR graph coefficient-projection denoising helpers
module simple_diff_map_denoise
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_diff_map_graphs, only: diffmap_graph, graph_matvec
use simple_image,           only: image
use simple_imgarr_utils,    only: copy_imgarr, dealloc_imgarr
use simple_linalg,          only: sparse_eigh
use simple_parameters,      only: parameters
use simple_polarft_calc,    only: polarft_calc
implicit none
#include "simple_local_flags.inc"

private
public :: graph_coeffproj_denoise

integer, parameter :: SE2_NTRANS_MODES = 2

type :: graph_mode_context
    type(diffmap_graph), pointer :: graph => null()
    integer :: mode = 0
end type graph_mode_context

contains

    subroutine graph_coeffproj_denoise(params, imgs, avg, graph, coeff_ptcls, rank_keep_override)
        type(parameters),         intent(in)  :: params
        class(image),             intent(in)  :: imgs(:)
        real,                     intent(in)  :: avg(:)
        type(diffmap_graph), target, intent(in) :: graph
        type(image), allocatable, intent(out) :: coeff_ptcls(:)
        integer, optional,        intent(in)  :: rank_keep_override
        select case(lowercase(trim(graph%steering)))
            case('none')
                call so2_coeffproj_denoise(params, imgs, avg, graph, coeff_ptcls, force_zero_theta=.true., &
                                           rank_keep_override=rank_keep_override)
            case('so2')
                call so2_coeffproj_denoise(params, imgs, avg, graph, coeff_ptcls, rank_keep_override=rank_keep_override)
            case('se2')
                call se2_sync_coeffproj_denoise(params, imgs, avg, graph, coeff_ptcls, rank_keep_override=rank_keep_override)
            case DEFAULT
                write(logfhandle,'(A,A)') 'Graph coefficient denoising skipped: unknown steering=', trim(graph%steering)
                call flush(logfhandle)
        end select
    end subroutine graph_coeffproj_denoise

    subroutine so2_coeffproj_denoise(params, imgs, avg, graph, coeff_ptcls, force_zero_theta, rank_keep_override)
        type(parameters),         intent(in)  :: params
        class(image),             intent(in)  :: imgs(:)
        real,                     intent(in)  :: avg(:)
        type(diffmap_graph), target, intent(in) :: graph
        type(image), allocatable, intent(out) :: coeff_ptcls(:)
        logical, optional,        intent(in)  :: force_zero_theta
        integer, optional,        intent(in)  :: rank_keep_override
        type(parameters) :: params_pft
        type(polarft_calc) :: pftc
        type(image), allocatable :: imgs_work(:)
        type(image) :: avg_img, synth_img
        complex(sp), allocatable :: pfts(:,:,:), pfts_hat(:,:,:), pft_tmp(:,:)
        complex(sp), allocatable :: mode_coeff(:,:), mode_hat(:,:)
        integer :: nptcls, ldim(3), kfromto(2), pdim_srch(3), pftsz, nrots, nk
        integer :: i, k, mode, mode_max, rank_keep
        real :: smpd
        logical :: zero_theta
        nptcls = size(imgs)
        if( nptcls < 2 ) return
        if( graph%n /= nptcls ) return
        zero_theta = .false.
        if( present(force_zero_theta) ) zero_theta = force_zero_theta
        if( .not. zero_theta .and. .not. graph%has_theta() ) return
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        if( size(avg) /= product(ldim) ) return
        call prepare_steerable_pft_params(params, ldim, smpd, params_pft)
        kfromto(1) = max(2, calc_fourier_index(params_pft%hp, params_pft%box, params_pft%smpd))
        kfromto(2) =        calc_fourier_index(params_pft%lp, params_pft%box, params_pft%smpd)
        if( kfromto(2) < kfromto(1) ) return
        call avg_img%new(ldim, smpd)
        call avg_img%unserialize(avg)
        call pftc%new(params_pft, 1, [1,nptcls], kfromto)
        pdim_srch = pftc%get_pdim_srch()
        pftsz     = pftc%get_pftsz()
        nrots     = pftc%get_nrots()
        nk        = kfromto(2) - kfromto(1) + 1
        mode_max  = min(max(0, params%steerable_nmodes), max(0, pftsz - 1))
        if( present(rank_keep_override) )then
            rank_keep = min(max(1, rank_keep_override), nptcls)
        else
            rank_keep = min(max(2, params%k_nn), nptcls)
        endif
        allocate(pfts(pftsz,nk,nptcls), pfts_hat(pftsz,nk,nptcls))
        pfts     = cmplx(0.,0.,kind=sp)
        pfts_hat = cmplx(0.,0.,kind=sp)
        allocate(pft_tmp(pftsz,kfromto(1):kfromto(2)))
        imgs_work = copy_imgarr(imgs)
        call imgs_work(1)%memoize4polarize(pdim_srch)
        do i = 1,nptcls
            call imgs_work(i)%subtr(avg_img)
            call imgs_work(i)%fft()
            call pftc%polarize_ptcl_pft(imgs_work(i), i, pdim_srch, oversamp=.false.)
            call pftc%get_ptcl_pft(i, pft_tmp)
            do k = 1,nk
                pfts(:,k,i) = pft_tmp(:,kfromto(1) + k - 1)
            end do
        end do
        call pftc%kill
        call dealloc_imgarr(imgs_work)
        allocate(mode_coeff(nptcls,nk), mode_hat(nptcls,nk))
        do mode = 0,mode_max
            call extract_pft_angular_mode(pfts, kfromto, nrots, mode, mode_coeff)
            call project_pft_angular_mode(graph, mode, mode_coeff, mode_hat, rank_keep, zero_theta)
            call accumulate_pft_angular_mode(pfts_hat, kfromto, nrots, mode, mode_hat)
            if( mode > 0 )then
                call extract_pft_angular_mode(pfts, kfromto, nrots, -mode, mode_coeff)
                call project_pft_angular_mode(graph, -mode, mode_coeff, mode_hat, rank_keep, zero_theta)
                call accumulate_pft_angular_mode(pfts_hat, kfromto, nrots, -mode, mode_hat)
            endif
        end do
        coeff_ptcls = copy_imgarr(imgs)
        do i = 1,nptcls
            call synthesize_pft_residual(pfts_hat(:,:,i), kfromto, nrots, ldim, smpd, synth_img)
            call coeff_ptcls(i)%copy(avg_img)
            call coeff_ptcls(i)%add(synth_img)
        end do
        write(logfhandle,'(A,I8,A,I8,A,I8,A,I8,A,I8,A,I8)') 'Graph coefficient denoising: n=', nptcls, &
            ' modes=', mode_max, ' rank=', rank_keep, ' pftsz=', pftsz, ' kmin=', kfromto(1), ' kmax=', kfromto(2)
        call flush(logfhandle)
        if( synth_img%exists() ) call synth_img%kill
        call avg_img%kill
        deallocate(pfts, pfts_hat, pft_tmp, mode_coeff, mode_hat)
    end subroutine so2_coeffproj_denoise

    subroutine se2_sync_coeffproj_denoise(params, imgs, avg, graph, coeff_ptcls, rank_keep_override)
        type(parameters),         intent(in)  :: params
        class(image),             intent(in)  :: imgs(:)
        real,                     intent(in)  :: avg(:)
        type(diffmap_graph), target, intent(in) :: graph
        type(image), allocatable, intent(out) :: coeff_ptcls(:)
        integer, optional,        intent(in)  :: rank_keep_override
        type(image), allocatable :: imgs_sync(:), den_sync(:)
        type(image) :: avg_img, tmp_img
        real, allocatable :: phi(:), sx(:), sy(:), rmat_rot(:,:,:)
        integer :: nptcls, ldim(3), i
        real :: smpd, angle_deg
        nptcls = size(imgs)
        if( nptcls < 2 ) return
        if( graph%n /= nptcls .or. .not. graph%has_theta() .or. .not. graph%has_shift() ) return
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        if( size(avg) /= product(ldim) ) return
        call synchronize_graph_potential(graph, 'theta', phi)
        call synchronize_graph_potential(graph, 'shift_x', sx)
        call synchronize_graph_potential(graph, 'shift_y', sy)
        allocate(rmat_rot(ldim(1),ldim(2),1), source=0.)
        imgs_sync = copy_imgarr(imgs)
        call avg_img%new(ldim, smpd)
        call tmp_img%new(ldim, smpd)
        call avg_img%unserialize(avg)
        do i = 1,nptcls
            call tmp_img%copy(imgs(i))
            call tmp_img%subtr(avg_img)
            call tmp_img%fft()
            call tmp_img%shift2Dserial([-sx(i), -sy(i)])
            call tmp_img%ifft()
            angle_deg = rad2deg(phi(i))
            call tmp_img%rtsq_serial(angle_deg, 0., 0., rmat_rot)
            call tmp_img%set_rmat(rmat_rot, .false.)
            call imgs_sync(i)%copy(avg_img)
            call imgs_sync(i)%add(tmp_img)
        end do
        call so2_coeffproj_denoise(params, imgs_sync, avg, graph, den_sync, force_zero_theta=.true., &
                                   rank_keep_override=rank_keep_override)
        if( .not. allocated(den_sync) )then
            call dealloc_imgarr(imgs_sync)
            call avg_img%kill
            call tmp_img%kill
            deallocate(phi, sx, sy, rmat_rot)
            return
        endif
        coeff_ptcls = copy_imgarr(imgs)
        do i = 1,nptcls
            call tmp_img%copy(den_sync(i))
            call tmp_img%subtr(avg_img)
            angle_deg = -rad2deg(phi(i))
            call tmp_img%rtsq_serial(angle_deg, 0., 0., rmat_rot)
            call tmp_img%set_rmat(rmat_rot, .false.)
            call tmp_img%fft()
            call tmp_img%shift2Dserial([sx(i), sy(i)])
            call tmp_img%ifft()
            call coeff_ptcls(i)%copy(avg_img)
            call coeff_ptcls(i)%add(tmp_img)
        end do
        write(logfhandle,'(A,I8,A,I8,A,F8.3,A,F8.3)') 'SE2 graph coefficient denoising: n=', nptcls, &
            ' directed_edges=', graph%nnz, ' phi_span_deg=', rad2deg(maxval(phi) - minval(phi)), &
            ' shift_span_px=', max(maxval(sx)-minval(sx), maxval(sy)-minval(sy))
        call flush(logfhandle)
        call dealloc_imgarr(imgs_sync)
        call dealloc_imgarr(den_sync)
        call avg_img%kill
        call tmp_img%kill
        deallocate(phi, sx, sy, rmat_rot)
    end subroutine se2_sync_coeffproj_denoise

    subroutine synchronize_graph_potential(graph, field, potential)
        type(diffmap_graph), intent(in) :: graph
        character(len=*),    intent(in) :: field
        real, allocatable,   intent(out) :: potential(:)
        logical, allocatable :: known(:)
        integer :: seed, i, j, p, nknown
        allocate(potential(graph%n), source=0.)
        allocate(known(graph%n), source=.false.)
        known = .false.
        seed = 1
        do while( seed <= graph%n )
            if( .not. known(seed) )then
                known(seed) = .true.
                do
                    nknown = count(known)
                    do i = 1,graph%n
                        if( .not. known(i) ) cycle
                        do p = graph%rowptr(i), graph%rowptr(i+1) - 1
                            j = graph%colind(p)
                            if( known(j) ) cycle
                            potential(j) = potential(i) + edge_delta(graph, p, field)
                            known(j) = .true.
                        end do
                    end do
                    if( count(known) == nknown ) exit
                end do
            endif
            seed = seed + 1
        end do
        deallocate(known)
    end subroutine synchronize_graph_potential

    real function edge_delta(graph, p, field) result(delta)
        type(diffmap_graph), intent(in) :: graph
        integer,             intent(in) :: p
        character(len=*),    intent(in) :: field
        select case(trim(field))
            case('theta')
                delta = graph%theta(p)
            case('shift_x')
                delta = graph%shift_x(p)
            case('shift_y')
                delta = graph%shift_y(p)
            case DEFAULT
                delta = 0.
        end select
    end function edge_delta

    subroutine prepare_steerable_pft_params(params, ldim, smpd, params_pft)
        type(parameters), intent(in)  :: params
        integer,          intent(in)  :: ldim(3)
        real,             intent(in)  :: smpd
        type(parameters), intent(out) :: params_pft
        params_pft             = params
        params_pft%objfun      = 'cc'
        params_pft%cc_objfun   = OBJFUN_CC
        params_pft%ctf         = 'no'
        params_pft%ldim        = ldim
        params_pft%box         = ldim(1)
        params_pft%box_crop    = ldim(1)
        params_pft%smpd        = smpd
        params_pft%smpd_crop   = smpd
        if( params_pft%pftsz < 1 ) params_pft%pftsz = magic_pftsz((real(ldim(1)) - COSMSKHALFWIDTH) / 2.)
    end subroutine prepare_steerable_pft_params

    subroutine project_pft_angular_mode(graph, mode, mode_coeff, mode_hat, rank_keep, zero_theta)
        type(diffmap_graph), target, intent(in) :: graph
        integer,     intent(in)  :: mode, rank_keep
        complex(sp), intent(in)  :: mode_coeff(:,:)
        complex(sp), intent(out) :: mode_hat(:,:)
        logical,     intent(in)  :: zero_theta
        type(graph_mode_context) :: ctx
        real, allocatable :: evals(:), evecs(:,:), y(:,:), yhat(:,:), coeff_r(:)
        complex(sp), allocatable :: coeff_c(:)
        integer :: n, n2, nk, nkeep, a, i, info, max_basis
        n  = size(mode_coeff,1)
        nk = size(mode_coeff,2)
        mode_hat = cmplx(0.,0.,kind=sp)
        if( n < 1 .or. nk < 1 ) return
        if( mode == 0 .or. zero_theta )then
            nkeep = min(max(1, rank_keep), n)
            allocate(evals(nkeep), evecs(n,nkeep), coeff_c(nk))
            max_basis = min(n, max(160, 8 * nkeep + 80))
            call sparse_eigh(graph_matvec, graph, n, nkeep, evals, evecs, tol=1.e-5, max_basis=max_basis, info=info)
            do a = 1,nkeep
                coeff_c = cmplx(0.,0.,kind=sp)
                do i = 1,n
                    coeff_c = coeff_c + real(evecs(i,a), kind=sp) * mode_coeff(i,:)
                end do
                do i = 1,n
                    mode_hat(i,:) = mode_hat(i,:) + real(evecs(i,a), kind=sp) * coeff_c
                end do
            end do
            deallocate(evals, evecs, coeff_c)
        else
            n2 = 2 * n
            nkeep = min(max(2, 2 * rank_keep), n2)
            ctx%graph => graph
            ctx%mode  = mode
            allocate(evals(nkeep), evecs(n2,nkeep), y(n2,nk), yhat(n2,nk), coeff_r(nk))
            max_basis = min(n2, max(160, 8 * nkeep + 80))
            call sparse_eigh(so2_graph_matvec, ctx, n2, nkeep, evals, evecs, tol=1.e-5, max_basis=max_basis, info=info)
            do i = 1,n
                y(i,:)   = real(mode_coeff(i,:))
                y(n+i,:) = aimag(mode_coeff(i,:))
            end do
            yhat = 0.
            do a = 1,nkeep
                coeff_r = 0.
                do i = 1,n2
                    coeff_r = coeff_r + evecs(i,a) * y(i,:)
                end do
                do i = 1,n2
                    yhat(i,:) = yhat(i,:) + evecs(i,a) * coeff_r
                end do
            end do
            do i = 1,n
                mode_hat(i,:) = cmplx(yhat(i,:), yhat(n+i,:), kind=sp)
            end do
            deallocate(evals, evecs, y, yhat, coeff_r)
        endif
    end subroutine project_pft_angular_mode

    subroutine so2_graph_matvec(ctx_any, x, y)
        class(*), intent(in)  :: ctx_any
        real,     intent(in)  :: x(:)
        real,     intent(out) :: y(:)
        type(diffmap_graph), pointer :: graph
        real :: w, phase, c, s
        integer :: n, i, j, p, mode
        select type(ctx => ctx_any)
        type is (graph_mode_context)
            graph => ctx%graph
            mode = ctx%mode
        class default
            THROW_HARD('invalid SO2 graph matvec context')
        end select
        n = graph%n
        if( size(x) /= 2*n .or. size(y) /= 2*n ) THROW_HARD('SO2 graph matvec shape mismatch')
        y = 0.
        do i = 1,n
            do p = graph%rowptr(i), graph%rowptr(i+1) - 1
                j = graph%colind(p)
                if( allocated(graph%wnorm) )then
                    w = graph%wnorm(p)
                else
                    w = graph%w(p)
                endif
                phase = real(mode) * graph%theta(p)
                c = cos(phase)
                s = sin(phase)
                y(i)   = y(i)   + w * ( c * x(j) - s * x(n+j) )
                y(n+i) = y(n+i) + w * ( s * x(j) + c * x(n+j) )
            end do
        end do
    end subroutine so2_graph_matvec

    subroutine extract_pft_angular_mode(pfts, kfromto, nrots, mode, mode_coeff)
        complex(sp), intent(in)  :: pfts(:,:,:)
        integer,     intent(in)  :: kfromto(2), nrots, mode
        complex(sp), intent(out) :: mode_coeff(:,:)
        complex(sp) :: acc, phase, val
        real :: theta, phase_arg, parity
        integer :: nptcls, pftsz, nk, i, irot, k
        nptcls = size(pfts,3)
        pftsz  = size(pfts,1)
        nk     = kfromto(2) - kfromto(1) + 1
        mode_coeff = cmplx(0.,0.,kind=sp)
        parity = merge(1., -1., mod(abs(mode), 2) == 0)
        do i = 1,nptcls
            do k = 1,nk
                acc = cmplx(0.,0.,kind=sp)
                do irot = 1,pftsz
                    theta     = twopi * real(irot - 1) / real(nrots)
                    phase_arg = -real(mode) * theta
                    phase     = cmplx(cos(phase_arg), sin(phase_arg), kind=sp)
                    val       = pfts(irot,k,i)
                    acc       = acc + val * phase + conjg(val) * phase * real(parity, kind=sp)
                end do
                mode_coeff(i,k) = acc / real(nrots, kind=sp)
            end do
        end do
    end subroutine extract_pft_angular_mode

    subroutine accumulate_pft_angular_mode(pfts_hat, kfromto, nrots, mode, mode_hat)
        complex(sp), intent(inout) :: pfts_hat(:,:,:)
        integer,     intent(in)    :: kfromto(2), nrots, mode
        complex(sp), intent(in)    :: mode_hat(:,:)
        complex(sp) :: phase
        real :: theta, phase_arg
        integer :: nptcls, pftsz, nk, i, irot, k
        nptcls = size(pfts_hat,3)
        pftsz  = size(pfts_hat,1)
        nk     = kfromto(2) - kfromto(1) + 1
        do i = 1,nptcls
            do k = 1,nk
                do irot = 1,pftsz
                    theta     = twopi * real(irot - 1) / real(nrots)
                    phase_arg = real(mode) * theta
                    phase     = cmplx(cos(phase_arg), sin(phase_arg), kind=sp)
                    pfts_hat(irot,k,i) = pfts_hat(irot,k,i) + mode_hat(i,k) * phase
                end do
            end do
        end do
    end subroutine accumulate_pft_angular_mode

    subroutine synthesize_pft_residual(pft_half, kfromto, nrots, ldim, smpd, img_out)
        complex(sp), intent(in)    :: pft_half(:,:)
        integer,     intent(in)    :: kfromto(2), nrots, ldim(3)
        real,        intent(in)    :: smpd
        type(image), intent(inout) :: img_out
        complex, allocatable :: cmat(:,:,:)
        real,    allocatable :: rho(:,:,:)
        real :: theta, x, y
        integer :: ashape(3), lims(3,2), pftsz, nk, irot, k, krad, h, l
        if( img_out%exists() ) call img_out%kill
        call img_out%new(ldim, smpd)
        call img_out%zero_and_flag_ft
        ashape = img_out%get_array_shape()
        allocate(rho(ashape(1),ashape(2),ashape(3)), source=0.)
        lims  = img_out%loop_lims(3)
        pftsz = size(pft_half,1)
        nk    = kfromto(2) - kfromto(1) + 1
        do irot = 1,pftsz
            theta = twopi * real(irot - 1) / real(nrots)
            do k = 1,nk
                krad = kfromto(1) + k - 1
                x =  sin(theta) * real(krad)
                y = -cos(theta) * real(krad)
                call scatter_fourier_sample(img_out, rho, lims, x, y, pft_half(irot,k))
            end do
        end do
        cmat = img_out%get_cmat()
        do l = 1,ashape(3)
            do k = 1,ashape(2)
                do h = 1,ashape(1)
                    if( rho(h,k,l) > 1.e-6 )then
                        cmat(h,k,l) = cmat(h,k,l) / rho(h,k,l)
                    else
                        cmat(h,k,l) = cmplx(0.,0.)
                    endif
                end do
            end do
        end do
        call img_out%set_cmat(cmat)
        call img_out%ifft()
        deallocate(cmat, rho)
    end subroutine synthesize_pft_residual

    subroutine scatter_fourier_sample(img, rho, lims, x, y, comp)
        type(image), intent(inout) :: img
        real,        intent(inout) :: rho(:,:,:)
        integer,     intent(in)    :: lims(3,2)
        real,        intent(in)    :: x, y
        complex(sp), intent(in)    :: comp
        complex :: comp_here
        integer :: h0, h1, k0, k1, h, k, ih, ik, phys(3)
        real :: dh, dk, wh(2), wk(2), w
        if( .not. ieee_is_finite(real(comp)) .or. .not. ieee_is_finite(aimag(comp)) ) return
        h0 = floor(x)
        k0 = floor(y)
        h1 = h0 + 1
        k1 = k0 + 1
        dh = x - real(h0)
        dk = y - real(k0)
        wh = [1. - dh, dh]
        wk = [1. - dk, dk]
        do ih = 1,2
            h = merge(h0, h1, ih == 1)
            if( h < lims(1,1) .or. h > lims(1,2) ) cycle
            do ik = 1,2
                k = merge(k0, k1, ik == 1)
                if( k < lims(2,1) .or. k > lims(2,2) ) cycle
                w = wh(ih) * wk(ik)
                if( w <= 1.e-6 ) cycle
                comp_here = cmplx(real(comp) * w, aimag(comp) * w)
                call img%add([h,k,0], comp_here, phys_out=phys)
                rho(phys(1),phys(2),phys(3)) = rho(phys(1),phys(2),phys(3)) + w
            end do
        end do
    end subroutine scatter_fourier_sample

end module simple_diff_map_denoise
