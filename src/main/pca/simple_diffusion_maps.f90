!@descr: diffusion maps embedding for nonlinear class splitting
module simple_diffusion_maps
use simple_core_module_api
!$ use omp_lib
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_image,                  only: image
use simple_imgarr_utils,           only: copy_imgarr, dealloc_imgarr
use simple_linalg,                 only: eigh, matinv
use simple_parameters,             only: parameters
use simple_polarft_calc,           only: polarft_calc
use simple_srch_sort_loc,          only: hpsort
use simple_stat,                   only: median, pearsn
implicit none
#include "simple_local_flags.inc"

private
public :: diffusion_map_embedder
public :: steerable_diffusion_map_embedder
public :: diffusion_map_generate_subavg_vectors
public :: steerable_transport_denoise
public :: steerable_coeffproj_denoise

type :: diffusion_map_embedder
    integer :: ndiff = 4  !< 0 requests automatic dimensionality selection from the diffusion spectrum
    integer :: k_nn  = 10
    logical :: keep_preimage_model = .false.
    real    :: decoder_ridge       = 1.e-3
    logical :: fit_joint_decoder   = .false.
    integer :: joint_decoder_iters = 80
    real    :: joint_decoder_weight = 1.
    real    :: joint_decoder_tol   = 1.e-4
    real, allocatable :: train_coords(:,:), preimage_targets(:,:), decoder_mat(:,:), joint_decoder_mat(:,:)
contains
    procedure :: set_params
    procedure :: set_preimage_params
    procedure :: embed
    procedure :: decode
    procedure :: joint_decode
    procedure :: trim_preimage_model
    procedure :: kill_preimage_model
end type diffusion_map_embedder

type :: steerable_diffusion_map_embedder
    integer :: ndiff  = 4
    integer :: k_nn   = 10
    integer :: nmodes = 4
contains
    procedure :: set_params => steerable_set_params
    procedure :: embed      => steerable_embed
end type steerable_diffusion_map_embedder

contains

    subroutine set_params(self, ndiff, k_nn)
        class(diffusion_map_embedder), intent(inout) :: self
        integer,                       intent(in)    :: ndiff, k_nn
        self%ndiff = max(0, ndiff)
        self%k_nn  = max(2, k_nn)
    end subroutine set_params

    subroutine set_preimage_params(self, keep_model, decoder_ridge, joint_decoder, &
                                   joint_decoder_iters, joint_decoder_weight, joint_decoder_tol)
        class(diffusion_map_embedder), intent(inout) :: self
        logical, optional,             intent(in)    :: keep_model
        logical, optional,             intent(in)    :: joint_decoder
        real,    optional,             intent(in)    :: decoder_ridge
        real,    optional,             intent(in)    :: joint_decoder_weight, joint_decoder_tol
        integer, optional,             intent(in)    :: joint_decoder_iters
        if( present(keep_model)     ) self%keep_preimage_model = keep_model
        if( present(decoder_ridge)  ) self%decoder_ridge       = max(decoder_ridge, 0.)
        if( present(joint_decoder)  ) self%fit_joint_decoder   = joint_decoder
        if( present(joint_decoder_iters) ) self%joint_decoder_iters = max(0, joint_decoder_iters)
        if( present(joint_decoder_weight) ) self%joint_decoder_weight = max(joint_decoder_weight, 0.)
        if( present(joint_decoder_tol) ) self%joint_decoder_tol = max(joint_decoder_tol, 1.e-8)
    end subroutine set_preimage_params

    subroutine embed(self, pcavecs, coords, eigvals, preimage_targets)
        class(diffusion_map_embedder), intent(inout) :: self
        real,                          intent(in)    :: pcavecs(:,:)
        real, allocatable,             intent(out)   :: coords(:,:)
        real, allocatable, optional,   intent(out)   :: eigvals(:)
        real, optional,                 intent(in)    :: preimage_targets(:,:)
        real, allocatable :: d2(:,:), aff(:,:), evals(:), evecs(:,:), kth_d2(:), deg(:), norm_aff(:,:)
        real, allocatable :: diff_evals(:)
        real :: eps, scale
        integer :: nptcls, npix, i, j, k, ndiff_used, ndiff_scan, k_used, nev
        integer(int64) :: t0, t1, trate, t_total0
        npix   = size(pcavecs, 1)
        nptcls = size(pcavecs, 2)
        if( self%keep_preimage_model ) call self%kill_preimage_model()
        if( nptcls < 3 )then
            allocate(coords(1, nptcls), source=0.)
            return
        endif
        if( self%ndiff <= 0 )then
            ndiff_scan = min(16, max(1, nptcls - 2))
        else
            ndiff_scan = min(max(1, self%ndiff), max(1, nptcls - 2))
        endif
        ndiff_used = ndiff_scan
        k_used     = min(max(2, self%k_nn), nptcls - 1)
        call system_clock(t_total0, trate)
        write(logfhandle,'(A,I8,A,I8,A,I8,A,I8)') 'Diffusion maps start: N=', nptcls, ' D=', npix, ' ndiff=', ndiff_scan, ' k_nn=', k_used
        call flush(logfhandle)
        allocate(d2(nptcls,nptcls), aff(nptcls,nptcls), kth_d2(nptcls), deg(nptcls), source=0.)
        call system_clock(t0)
        !$omp parallel default(shared) private(i,j,k,scale)
        !$omp do schedule(dynamic)
        do i = 1, nptcls - 1
            do j = i + 1, nptcls
                scale = 0.
                do k = 1, npix
                    scale = scale + (pcavecs(k,i) - pcavecs(k,j))**2
                end do
                d2(i,j) = scale
                d2(j,i) = d2(i,j)
            end do
        end do
        !$omp end do
        !$omp do schedule(static)
        do i = 1, nptcls
            call kth_distance_for_point(d2(i,:), i, k_used, kth_d2(i))
        end do
        !$omp end do
        !$omp end parallel
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A)') 'Diffusion maps pairwise distances: ', real(t1-t0)/real(trate), ' s'
        call flush(logfhandle)
        call system_clock(t0)
        eps = median(kth_d2)
        if( eps < DTINY ) eps = max(sum(kth_d2) / real(max(size(kth_d2),1)), 1.e-6)
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A,ES10.3)') 'Diffusion maps local scale selection: ', real(t1-t0)/real(trate), ' s; eps=', eps
        call flush(logfhandle)
        call system_clock(t0)
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do i = 1, nptcls
            do j = 1, nptcls
                if( i == j ) cycle
                if( d2(i,j) <= kth_d2(i) .or. d2(i,j) <= kth_d2(j) )then
                    aff(i,j) = exp(-d2(i,j) / eps)
                endif
            end do
        end do
        !$omp end parallel do
        aff = 0.5 * (aff + transpose(aff))
        deg = sum(aff, dim=2)
        do i = 1, nptcls
            if( deg(i) > DTINY ) cycle
            j = nearest_neighbor_index(d2(i,:), i)
            aff(i,j) = 1.
            aff(j,i) = 1.
        end do
        deg = sum(aff, dim=2)
        allocate(norm_aff(nptcls,nptcls), source=0.)
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do i = 1, nptcls
            do j = 1, nptcls
                if( aff(i,j) <= DTINY ) cycle
                norm_aff(i,j) = aff(i,j) / sqrt(max(deg(i), DTINY) * max(deg(j), DTINY))
            end do
        end do
        !$omp end parallel do
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A)') 'Diffusion maps graph/normalization: ', real(t1-t0)/real(trate), ' s'
        call flush(logfhandle)
        nev = ndiff_scan + 1
        allocate(evals(nptcls), evecs(nptcls,nptcls))
        call system_clock(t0)
        call eigh(nptcls, norm_aff, nev, evals, evecs)
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A,I8)') 'Diffusion maps eigensolve: ', real(t1-t0)/real(trate), ' s; nev=', nev
        call flush(logfhandle)
        call system_clock(t0)
        allocate(diff_evals(ndiff_scan), source=0.)
        do k = 1, ndiff_scan
            diff_evals(k) = evals(nev - k)
        end do
        if( self%ndiff <= 0 ) ndiff_used = auto_ndiff_from_eigengap(diff_evals)
        allocate(coords(ndiff_used, nptcls), source=0.)
        if( present(eigvals) ) allocate(eigvals(ndiff_scan), source=diff_evals)
        !$omp parallel do default(shared) private(i,k,j) schedule(static)
        do k = 1, ndiff_used
            j = nev - k
            do i = 1, nptcls
                coords(k,i) = evals(j) * evecs(i,j)
            end do
        end do
        !$omp end parallel do
        if( self%keep_preimage_model )then
            call normalize_coords(coords)
            if( present(preimage_targets) )then
                call store_preimage_model(self, preimage_targets, coords)
            else
                call store_preimage_model(self, pcavecs, coords)
            endif
        else
            call normalize_coords(coords)
        endif
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A,I8)') 'Diffusion maps embedding build: ', real(t1-t0)/real(trate), ' s; dims=', ndiff_used
        call flush(logfhandle)
        write(logfhandle,'(A,F8.3,A,I8,A,ES10.3,A,ES10.3)') 'Diffusion maps total: ', real(t1-t_total0)/real(trate), &
            ' s; dims=', ndiff_used, ' lambda1=', diff_evals(1), ' lambda_last=', diff_evals(ndiff_scan)
        call flush(logfhandle)
        deallocate(d2, aff, kth_d2, deg, norm_aff, evals, evecs, diff_evals)
    end subroutine embed

    subroutine store_preimage_model(self, targets, coords)
        class(diffusion_map_embedder), intent(inout) :: self
        real,                          intent(in)    :: targets(:,:), coords(:,:)
        if( size(targets,2) /= size(coords,2) ) THROW_HARD('preimage target particle-count mismatch; diffusion map embed')
        allocate(self%train_coords(size(coords,1),size(coords,2)), source=coords)
        allocate(self%preimage_targets(size(targets,1),size(targets,2)), source=targets)
        call fit_latent_decoder(self, coords, targets)
        if( self%fit_joint_decoder ) call fit_joint_latent_decoder(self, coords, targets)
    end subroutine store_preimage_model

    subroutine trim_preimage_model(self, nkeep)
        class(diffusion_map_embedder), intent(inout) :: self
        integer,                       intent(in)    :: nkeep
        real, allocatable :: coords_tmp(:,:)
        integer :: n
        if( .not. allocated(self%train_coords) ) return
        if( .not. allocated(self%preimage_targets) ) return
        n = min(max(1, nkeep), size(self%train_coords,1))
        if( n == size(self%train_coords,1) ) return
        allocate(coords_tmp(n,size(self%train_coords,2)), source=self%train_coords(1:n,:))
        call move_alloc(coords_tmp, self%train_coords)
        if( allocated(self%decoder_mat) ) deallocate(self%decoder_mat)
        if( allocated(self%joint_decoder_mat) ) deallocate(self%joint_decoder_mat)
        call fit_latent_decoder(self, self%train_coords, self%preimage_targets)
        if( self%fit_joint_decoder ) call fit_joint_latent_decoder(self, self%train_coords, self%preimage_targets)
    end subroutine trim_preimage_model

    subroutine fit_latent_decoder(self, coords, targets)
        class(diffusion_map_embedder), intent(inout) :: self
        real,                          intent(in)    :: coords(:,:), targets(:,:)
        real, allocatable :: design(:,:), gram(:,:), invgram(:,:), rhs(:,:)
        integer :: ndim, nptcls, nfeat, i, err
        ndim   = size(coords,1)
        nptcls = size(coords,2)
        nfeat  = ndim + 1
        allocate(design(nfeat,nptcls), gram(nfeat,nfeat), invgram(nfeat,nfeat), source=0.)
        design(1,:) = 1.
        design(2:nfeat,:) = coords
        gram = matmul(design, transpose(design))
        gram(1,1) = gram(1,1) + self%decoder_ridge * 1.e-3
        do i = 2, nfeat
            gram(i,i) = gram(i,i) + self%decoder_ridge
        end do
        call matinv(gram, invgram, nfeat, err)
        if( err == -1 )then
            invgram = 0.
            do i = 1, nfeat
                invgram(i,i) = 1. / max(gram(i,i), 1.e-6)
            end do
        endif
        rhs = matmul(targets, transpose(design))
        self%decoder_mat = matmul(rhs, invgram)
        deallocate(design, gram, invgram, rhs)
    end subroutine fit_latent_decoder

    subroutine fit_image_encoder(self, targets, coords, encoder_mat)
        class(diffusion_map_embedder), intent(in)  :: self
        real,                          intent(in)  :: targets(:,:), coords(:,:)
        real, allocatable,             intent(out) :: encoder_mat(:,:)
        real, allocatable :: loadings(:,:), gram(:,:), invgram(:,:)
        integer :: npix, ndim, i, err
        if( .not. allocated(self%decoder_mat) ) THROW_HARD('diffusion decoder not fitted; fit_image_encoder')
        npix   = size(targets,1)
        ndim   = size(coords,1)
        if( size(self%decoder_mat,1) /= npix .or. size(self%decoder_mat,2) /= ndim + 1 )then
            THROW_HARD('diffusion decoder/encoder shape mismatch; fit_image_encoder')
        endif
        allocate(loadings(npix,ndim), source=self%decoder_mat(:,2:ndim+1))
        allocate(gram(ndim,ndim), invgram(ndim,ndim), source=0.)
        gram = matmul(transpose(loadings), loadings)
        do i = 1, ndim
            gram(i,i) = gram(i,i) + self%decoder_ridge
        end do
        call matinv(gram, invgram, ndim, err)
        if( err == -1 )then
            invgram = 0.
            do i = 1, ndim
                invgram(i,i) = 1. / max(gram(i,i), 1.e-6)
            end do
        endif
        allocate(encoder_mat(ndim,npix+1), source=0.)
        encoder_mat(:,2:npix+1) = matmul(invgram, transpose(loadings))
        encoder_mat(:,1) = -matmul(encoder_mat(:,2:npix+1), self%decoder_mat(:,1))
        deallocate(loadings, gram, invgram)
    end subroutine fit_image_encoder

    subroutine fit_joint_latent_decoder(self, coords, targets)
        class(diffusion_map_embedder), intent(inout) :: self
        real,                          intent(in)    :: coords(:,:), targets(:,:)
        real, allocatable :: design(:,:), encoder_mat(:,:), xhat(:,:), zhat(:,:), img_resid(:,:), lat_resid(:,:)
        real, allocatable :: grad_x(:,:), grad(:,:)
        real :: loss, prev_loss, step, grad_norm, weight
        integer :: ndim, nptcls, nfeat, iter
        if( .not. allocated(self%decoder_mat) ) return
        ndim   = size(coords,1)
        nptcls = size(coords,2)
        nfeat  = ndim + 1
        weight = self%joint_decoder_weight
        call fit_image_encoder(self, targets, coords, encoder_mat)
        allocate(design(nfeat,nptcls), source=0.)
        design(1,:)       = 1.
        design(2:nfeat,:) = coords
        if( allocated(self%joint_decoder_mat) ) deallocate(self%joint_decoder_mat)
        allocate(self%joint_decoder_mat(size(self%decoder_mat,1),size(self%decoder_mat,2)), source=self%decoder_mat)
        allocate(xhat(size(targets,1),nptcls), zhat(ndim,nptcls), img_resid(size(targets,1),nptcls), &
                 lat_resid(ndim,nptcls), grad_x(size(targets,1),nptcls), grad(size(self%joint_decoder_mat,1),nfeat), source=0.)
        prev_loss = huge(prev_loss)
        do iter = 1, self%joint_decoder_iters
            xhat = matmul(self%joint_decoder_mat, design)
            zhat = spread(encoder_mat(:,1), 2, nptcls) + matmul(encoder_mat(:,2:size(encoder_mat,2)), xhat)
            if( .not. all(ieee_is_finite(xhat)) .or. .not. all(ieee_is_finite(zhat)) ) exit
            img_resid = xhat - targets
            lat_resid = zhat - coords
            loss = 0.5 * sum(img_resid**2) / real(max(1, nptcls)) + &
                   0.5 * weight * sum(lat_resid**2) / real(max(1, nptcls))
            if( .not. ieee_is_finite(loss) ) exit
            if( iter > 1 .and. abs(prev_loss - loss) <= self%joint_decoder_tol * max(1., abs(prev_loss)) ) exit
            prev_loss = loss
            grad_x = img_resid + weight * matmul(transpose(encoder_mat(:,2:size(encoder_mat,2))), lat_resid)
            grad = matmul(grad_x, transpose(design)) / real(max(1, nptcls))
            grad(:,1) = grad(:,1) + self%decoder_ridge * 1.e-3 * self%joint_decoder_mat(:,1)
            grad(:,2:nfeat) = grad(:,2:nfeat) + self%decoder_ridge * self%joint_decoder_mat(:,2:nfeat)
            grad_norm = sqrt(sum(grad**2) / real(max(1, size(grad))))
            if( .not. ieee_is_finite(grad_norm) .or. grad_norm <= self%joint_decoder_tol ) exit
            step = 0.005 / sqrt(real(iter))
            self%joint_decoder_mat = self%joint_decoder_mat - step * grad
        end do
        deallocate(design, encoder_mat, xhat, zhat, img_resid, lat_resid, grad_x, grad)
    end subroutine fit_joint_latent_decoder

    subroutine decode(self, z, x)
        class(diffusion_map_embedder), intent(in)  :: self
        real,                          intent(in)  :: z(:)
        real, allocatable,             intent(out) :: x(:)
        integer :: ndim
        if( .not. allocated(self%decoder_mat) ) THROW_HARD('diffusion preimage decoder is not fitted')
        ndim = size(self%decoder_mat,2) - 1
        if( size(z) /= ndim ) THROW_HARD('latent coordinate size mismatch; diffusion decode')
        allocate(x(size(self%decoder_mat,1)), source=self%decoder_mat(:,1))
        x = x + matmul(self%decoder_mat(:,2:ndim+1), z)
    end subroutine decode

    subroutine joint_decode(self, z, x)
        class(diffusion_map_embedder), intent(in)  :: self
        real,                          intent(in)  :: z(:)
        real, allocatable,             intent(out) :: x(:)
        integer :: ndim
        if( .not. allocated(self%joint_decoder_mat) )then
            call self%decode(z, x)
            return
        endif
        ndim = size(self%joint_decoder_mat,2) - 1
        if( size(z) /= ndim ) THROW_HARD('latent coordinate size mismatch; diffusion joint_decode')
        allocate(x(size(self%joint_decoder_mat,1)), source=self%joint_decoder_mat(:,1))
        x = x + matmul(self%joint_decoder_mat(:,2:ndim+1), z)
    end subroutine joint_decode

    subroutine kill_preimage_model(self)
        class(diffusion_map_embedder), intent(inout) :: self
        if( allocated(self%train_coords)  ) deallocate(self%train_coords)
        if( allocated(self%preimage_targets) ) deallocate(self%preimage_targets)
        if( allocated(self%decoder_mat)   ) deallocate(self%decoder_mat)
        if( allocated(self%joint_decoder_mat) ) deallocate(self%joint_decoder_mat)
    end subroutine kill_preimage_model

    integer function auto_ndiff_from_eigengap(eigvals) result(ndiff)
        real, intent(in) :: eigvals(:)
        real :: best_gap, gap
        integer :: i, maxcand
        maxcand = size(eigvals)
        if( maxcand <= 1 )then
            ndiff = 1
            return
        endif
        best_gap = -huge(best_gap)
        ndiff = min(4, maxcand)
        do i = 1, maxcand - 1
            gap = eigvals(i) - eigvals(i+1)
            if( gap > best_gap )then
                best_gap = gap
                ndiff = i
            endif
        end do
        ! Avoid collapsing the embedding to one coordinate unless the spectrum is tiny.
        ndiff = max(min(2, maxcand), min(maxcand, ndiff))
    end function auto_ndiff_from_eigengap

    subroutine normalize_coords(coords, coord_mu, coord_sigma)
        real, intent(inout) :: coords(:,:)
        real, allocatable, optional, intent(out) :: coord_mu(:), coord_sigma(:)
        real :: mu, sigma
        integer :: i
        if( present(coord_mu) )then
            if( allocated(coord_mu) ) deallocate(coord_mu)
            allocate(coord_mu(size(coords,1)), source=0.)
        endif
        if( present(coord_sigma) )then
            if( allocated(coord_sigma) ) deallocate(coord_sigma)
            allocate(coord_sigma(size(coords,1)), source=1.)
        endif
        do i = 1, size(coords,1)
            mu = sum(coords(i,:)) / real(size(coords,2))
            sigma = sqrt(sum((coords(i,:) - mu)**2) / real(max(size(coords,2)-1, 1)))
            if( sigma < 1.e-6 ) sigma = 1.0
            if( present(coord_mu)    ) coord_mu(i)    = mu
            if( present(coord_sigma) ) coord_sigma(i) = sigma
            coords(i,:) = (coords(i,:) - mu) / sigma
        end do
    end subroutine normalize_coords

    integer function nearest_neighbor_index(row, self_idx) result(ind)
        real,    intent(in) :: row(:)
        integer, intent(in) :: self_idx
        real :: best
        integer :: i
        ind  = 1
        best = huge(best)
        do i = 1, size(row)
            if( i == self_idx ) cycle
            if( row(i) < best )then
                best = row(i)
                ind = i
            endif
        end do
    end function nearest_neighbor_index

    subroutine kth_distance_for_point(d2row, self_idx, k_used, kth)
        real,    intent(in)  :: d2row(:)
        integer, intent(in)  :: self_idx, k_used
        real,    intent(out) :: kth
        real :: work(max(1,size(d2row)-1))
        integer :: j, k
        k = 0
        do j = 1, size(d2row)
            if( j == self_idx ) cycle
            k = k + 1
            work(k) = d2row(j)
        end do
        call hpsort(work(1:k))
        kth = work(min(max(1, k_used), k))
    end subroutine kth_distance_for_point

    subroutine kth_distance_with_self(d2row, k_used, kth)
        real,    intent(in)  :: d2row(:)
        integer, intent(in)  :: k_used
        real,    intent(out) :: kth
        real :: work(size(d2row))
        work = d2row
        call hpsort(work)
        kth = work(min(max(1, k_used), size(work)))
    end subroutine kth_distance_with_self

    subroutine diffusion_map_generate_subavg_vectors(params, context_id, coords, pcavecs, avg, labels, nlabels, &
                                                     subavg_vecs, l_generated)
        type(parameters),  intent(in)  :: params
        integer,           intent(in)  :: context_id, nlabels
        real,              intent(in)  :: coords(:,:), pcavecs(:,:), avg(:)
        integer,           intent(in)  :: labels(:)
        real, allocatable, intent(out) :: subavg_vecs(:,:)
        logical,           intent(out) :: l_generated
        real, allocatable :: d2(:,:), kth_d2(:)
        integer, allocatable :: pops(:)
        real :: tau, w, wsum, scale
        integer :: nptcls, ndim, npix, k_used, i, j, lab, nedges, nfinite
        l_generated = .false.
        nptcls = size(coords, 2)
        ndim   = size(coords, 1)
        npix   = size(avg)
        if( nptcls < 2 .or. ndim < 1 .or. nlabels < 1 ) return
        if( size(pcavecs, 1) /= npix .or. size(pcavecs, 2) /= nptcls .or. size(labels) /= nptcls )then
            write(logfhandle,'(A,I8)') 'Diffusion map vector preimage skipped: shape mismatch context=', context_id
            call flush(logfhandle)
            return
        endif
        if( .not. all(ieee_is_finite(coords)) )then
            write(logfhandle,'(A,I8)') 'Diffusion map vector preimage skipped: nonfinite coords context=', context_id
            call flush(logfhandle)
            return
        endif
        k_used = min(max(2, params%k_nn), nptcls)
        allocate(d2(nptcls,nptcls), kth_d2(nptcls), source=0.)
        allocate(pops(nlabels), source=0)
        allocate(subavg_vecs(npix,nlabels), source=0.)
        do i = 1, nptcls - 1
            do j = i + 1, nptcls
                d2(i,j) = sum((coords(1:ndim,i) - coords(1:ndim,j))**2)
                d2(j,i) = d2(i,j)
            end do
        end do
        do i = 1, nptcls
            call kth_distance_with_self(d2(i,:), k_used, kth_d2(i))
        end do
        tau = median(kth_d2)
        if( .not. ieee_is_finite(tau) .or. tau < DTINY )then
            nfinite = count(ieee_is_finite(kth_d2))
            if( nfinite > 0 ) tau = max(sum(kth_d2, mask=ieee_is_finite(kth_d2)) / real(nfinite), 1.e-6)
        endif
        if( .not. ieee_is_finite(tau) .or. tau < DTINY )then
            write(logfhandle,'(A,I8)') 'Diffusion map vector preimage skipped: invalid bandwidth context=', context_id
            call flush(logfhandle)
            deallocate(d2, kth_d2, pops, subavg_vecs)
            return
        endif
        nedges = 0
        do i = 1, nptcls
            lab = labels(i)
            if( lab < 1 .or. lab > nlabels ) cycle
            wsum = 0.
            do j = 1, nptcls
                if( d2(i,j) > kth_d2(i) ) cycle
                w = exp(-d2(i,j) / tau)
                if( w <= DTINY .or. .not. ieee_is_finite(w) ) cycle
                wsum = wsum + w
            end do
            pops(lab) = pops(lab) + 1
            if( wsum <= DTINY )then
                subavg_vecs(:,lab) = subavg_vecs(:,lab) + pcavecs(:,i)
            else
                do j = 1, nptcls
                    if( d2(i,j) > kth_d2(i) ) cycle
                    w = exp(-d2(i,j) / tau)
                    if( w <= DTINY .or. .not. ieee_is_finite(w) ) cycle
                    scale = w / wsum
                    subavg_vecs(:,lab) = subavg_vecs(:,lab) + scale * pcavecs(:,j)
                    nedges = nedges + 1
                end do
            endif
        end do
        do lab = 1, nlabels
            if( pops(lab) > 0 )then
                subavg_vecs(:,lab) = avg + subavg_vecs(:,lab) / real(pops(lab))
            else
                subavg_vecs(:,lab) = avg
            endif
        end do
        l_generated = .true.
        write(logfhandle,'(A,I8,A,I8,A,I8,A,ES10.3)') 'Diffusion map vector preimage: context=', context_id, &
            ' nlabels=', nlabels, ' directed_edges=', nedges, ' tau=', tau
        call flush(logfhandle)
        deallocate(d2, kth_d2, pops)
    end subroutine diffusion_map_generate_subavg_vectors

    subroutine steerable_set_params(self, ndiff, k_nn, nmodes)
        class(steerable_diffusion_map_embedder), intent(inout) :: self
        integer,                                 intent(in)    :: ndiff, k_nn, nmodes
        self%ndiff  = max(0, ndiff)
        self%k_nn   = max(2, k_nn)
        self%nmodes = max(0, nmodes)
    end subroutine steerable_set_params

    subroutine steerable_embed(self, params, imgs, coords, eigvals, graph_aff, graph_theta)
        class(steerable_diffusion_map_embedder), intent(inout) :: self
        type(parameters),                        intent(in)    :: params
        class(image),                            intent(in)    :: imgs(:)
        real, allocatable,                       intent(out)   :: coords(:,:)
        real, allocatable, optional,             intent(out)   :: eigvals(:)
        real, allocatable, optional,             intent(out)   :: graph_aff(:,:), graph_theta(:,:)
        type(image), allocatable :: imgs_work(:)
        real, allocatable        :: d2(:,:), aff(:,:), kth_d2(:), deg(:), norm_aff(:,:), theta(:,:)
        real, allocatable        :: inpl_corrs(:,:), cand_coords(:,:), cand_eigvals(:), out_eigvals(:)
        type(parameters)    :: params_pft
        type(polarft_calc)  :: pftc
        real                :: corr, eps, smpd
        integer             :: nptcls, ndiff_scan, k_used, ldim(3), kfromto(2), pdim_srch(3), nrots, nthr_use
        integer             :: i, j, loc, ncand, ithr
        integer(int64)      :: t0, t1, trate, t_total0
        nptcls = size(imgs)
        if( nptcls < 3 )then
            allocate(coords(1, nptcls), source=0.)
            if( present(eigvals) ) allocate(eigvals(1), source=0.)
            if( present(graph_aff)   ) allocate(graph_aff(nptcls,nptcls), source=0.)
            if( present(graph_theta) ) allocate(graph_theta(nptcls,nptcls), source=0.)
            return
        endif
        if( self%ndiff <= 0 )then
            ndiff_scan = min(16, max(1, nptcls - 2))
        else
            ndiff_scan = min(max(1, self%ndiff), max(1, nptcls - 2))
        endif
        k_used = min(max(2, self%k_nn), nptcls - 1)
        ldim   = imgs(1)%get_ldim()
        smpd   = imgs(1)%get_smpd()
        call prepare_steerable_pft_params(params, ldim, smpd, params_pft)
        nthr_use = max(1, params_pft%nthr)
        kfromto(1) = max(2, calc_fourier_index(params_pft%hp, params_pft%box, params_pft%smpd))
        kfromto(2) =        calc_fourier_index(params_pft%lp, params_pft%box, params_pft%smpd)
        if( kfromto(2) < kfromto(1) )then
            write(logfhandle,'(A,I8,A,I8)') 'Steerable diffusion maps invalid Fourier limits: kmin=', &
                kfromto(1), ' kmax=', kfromto(2)
            call flush(logfhandle)
            allocate(coords(1, nptcls), source=0.)
            if( present(eigvals) ) allocate(eigvals(1), source=0.)
            if( present(graph_aff)   ) allocate(graph_aff(nptcls,nptcls), source=0.)
            if( present(graph_theta) ) allocate(graph_theta(nptcls,nptcls), source=0.)
            return
        endif
        call system_clock(t_total0, trate)
        write(logfhandle,'(A,I8,A,I8,A,I8,A,I8,A,I4,A,I8)') 'Steerable diffusion maps start: N=', nptcls, &
            ' box=', ldim(1), ' ndiff=', ndiff_scan, ' k_nn=', k_used, ' nmodes=', self%nmodes, ' pftsz=', params_pft%pftsz
        call flush(logfhandle)
        allocate(d2(nptcls,nptcls), aff(nptcls,nptcls), kth_d2(nptcls), deg(nptcls), theta(nptcls,nptcls), source=0.)
        call system_clock(t0)
        imgs_work = copy_imgarr(imgs)
        call pftc%new(params_pft, nptcls, [1,nptcls], kfromto)
        pdim_srch = pftc%get_pdim_srch()
        call imgs_work(1)%memoize4polarize(pdim_srch)
        do i = 1, nptcls
            call imgs_work(i)%fft()
            call pftc%polarize_ref_pft( imgs_work(i), i, iseven=.true., pdim=pdim_srch, oversamp=.false.)
            call pftc%polarize_ptcl_pft(imgs_work(i), i,               pdim=pdim_srch, oversamp=.false.)
            call imgs_work(i)%ifft()
        end do
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        nrots = pftc%get_nrots()
        allocate(inpl_corrs(nrots,nthr_use))
        ithr = 1
        !$omp parallel private(i,j,loc,corr,ithr) default(shared) num_threads(nthr_use) proc_bind(close)
        !$omp do schedule(dynamic)
        do i = 1, nptcls - 1
            do j = i + 1, nptcls
                call pftc%gen_objfun_vals(j, i, [0.,0.], inpl_corrs(:,ithr))
                loc  = maxloc(inpl_corrs(:,ithr), dim=1)
                corr = inpl_corrs(loc,ithr)
                if( .not. ieee_is_finite(corr) ) corr = 0.
                corr = min(1., max(-1., corr))
                d2(i,j) = max(0., 1. - corr)
                d2(j,i) = d2(i,j)
                theta(i,j) = deg2rad(pftc%get_rot(loc))
                theta(j,i) = -theta(i,j)
            end do
        end do
        !$omp end do
        !$omp end parallel
        deallocate(inpl_corrs)
        call pftc%kill
        call dealloc_imgarr(imgs_work)
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A)') 'Steerable diffusion maps PFT pairwise rotations: ', real(t1-t0)/real(trate), ' s'
        call flush(logfhandle)
        call system_clock(t0)
        !$omp parallel do default(shared) private(i) schedule(static)
        do i = 1, nptcls
            call kth_distance_for_point(d2(i,:), i, k_used, kth_d2(i))
        end do
        !$omp end parallel do
        eps = median(kth_d2)
        if( eps < DTINY ) eps = max(sum(kth_d2) / real(max(size(kth_d2),1)), 1.e-6)
        do i = 1, nptcls
            do j = 1, nptcls
                if( i == j ) cycle
                if( d2(i,j) <= kth_d2(i) .or. d2(i,j) <= kth_d2(j) ) aff(i,j) = exp(-d2(i,j) / eps)
            end do
        end do
        aff = 0.5 * (aff + transpose(aff))
        deg = sum(aff, dim=2)
        do i = 1, nptcls
            if( deg(i) > DTINY ) cycle
            j = nearest_neighbor_index(d2(i,:), i)
            aff(i,j) = 1.
            aff(j,i) = 1.
        end do
        deg = sum(aff, dim=2)
        allocate(norm_aff(nptcls,nptcls), source=0.)
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do i = 1, nptcls
            do j = 1, nptcls
                if( aff(i,j) <= DTINY ) cycle
                norm_aff(i,j) = aff(i,j) / sqrt(max(deg(i), DTINY) * max(deg(j), DTINY))
            end do
        end do
        !$omp end parallel do
        if( present(graph_aff)   ) allocate(graph_aff(nptcls,nptcls), source=aff)
        if( present(graph_theta) ) allocate(graph_theta(nptcls,nptcls), source=theta)
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A,ES10.3)') 'Steerable diffusion maps graph/normalization: ', &
            real(t1-t0)/real(trate), ' s; eps=', eps
        call flush(logfhandle)
        call system_clock(t0)
        allocate(cand_coords(nptcls, max(1, ndiff_scan * (self%nmodes + 1))), &
                 cand_eigvals(max(1, ndiff_scan * (self%nmodes + 1))), source=0.)
        call collect_steerable_candidates(norm_aff, theta, self%nmodes, ndiff_scan, cand_coords, cand_eigvals, ncand)
        call select_output_candidates(cand_coords, cand_eigvals, ncand, ndiff_scan, coords, out_eigvals)
        call normalize_coords(coords)
        if( present(eigvals) ) allocate(eigvals(size(out_eigvals)), source=out_eigvals)
        call system_clock(t1)
        write(logfhandle,'(A,F8.3,A,I8,A,I8,A,ES10.3)') 'Steerable diffusion maps eigensolve/embedding: ', &
            real(t1-t0)/real(trate), ' s; candidates=', ncand, ' dims=', size(coords,1), ' lambda1=', out_eigvals(1)
        write(logfhandle,'(A,F8.3,A,I8)') 'Steerable diffusion maps total: ', &
            real(t1-t_total0)/real(trate), ' s; dims=', size(coords,1)
        call flush(logfhandle)
        deallocate(d2, aff, kth_d2, deg, norm_aff, theta, cand_coords, cand_eigvals, out_eigvals)
    end subroutine steerable_embed

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

    subroutine collect_steerable_candidates(norm_aff, theta, nmodes, ndiff_scan, cand_coords, cand_eigvals, ncand)
        real,    intent(in)    :: norm_aff(:,:), theta(:,:)
        integer, intent(in)    :: nmodes, ndiff_scan
        real,    intent(inout) :: cand_coords(:,:), cand_eigvals(:)
        integer, intent(out)   :: ncand
        real, allocatable :: mat(:,:), evals(:), evecs(:,:), feat(:)
        integer :: n, n2, nev, mode, idx, k
        n = size(norm_aff,1)
        n2 = 2 * n
        ncand = 0
        allocate(feat(n))
        nev = min(ndiff_scan + 1, n)
        allocate(mat(n,n), evals(n), evecs(n,n))
        mat = norm_aff
        call eigh(n, mat, nev, evals, evecs)
        do k = 1, min(ndiff_scan, nev - 1)
            idx = nev - k
            feat = evals(idx) * evecs(:,idx)
            call add_steerable_candidate(feat, evals(idx), cand_coords, cand_eigvals, ncand)
        end do
        deallocate(mat, evals, evecs)
        do mode = 1, nmodes
            if( ncand >= size(cand_eigvals) ) exit
            nev = min(2 * ndiff_scan + 2, n2)
            allocate(mat(n2,n2), evals(n2), evecs(n2,n2))
            call build_connection_block(norm_aff, theta, mode, mat)
            call eigh(n2, mat, nev, evals, evecs)
            do idx = nev, 1, -1
                feat = evals(idx) * sqrt(evecs(1:n,idx)**2 + evecs(n+1:n2,idx)**2)
                call add_steerable_candidate(feat, evals(idx), cand_coords, cand_eigvals, ncand)
                if( ncand >= size(cand_eigvals) ) exit
            end do
            deallocate(mat, evals, evecs)
        end do
        deallocate(feat)
    end subroutine collect_steerable_candidates

    subroutine build_connection_block(norm_aff, theta, mode, conn_block)
        real,    intent(in)  :: norm_aff(:,:), theta(:,:)
        integer, intent(in)  :: mode
        real,    intent(out) :: conn_block(:,:)
        real :: a, b, phase
        integer :: n, i, j
        n = size(norm_aff,1)
        conn_block = 0.
        do i = 1, n
            do j = 1, n
                if( norm_aff(i,j) <= DTINY ) cycle
                phase = real(mode) * theta(i,j)
                a = norm_aff(i,j) * cos(phase)
                b = norm_aff(i,j) * sin(phase)
                conn_block(i,   j)   =  a
                conn_block(i,   n+j) = -b
                conn_block(n+i, j)   =  b
                conn_block(n+i, n+j) =  a
            end do
        end do
    end subroutine build_connection_block

    subroutine add_steerable_candidate(feat, eigval, cand_coords, cand_eigvals, ncand)
        real,    intent(in)    :: feat(:), eigval
        real,    intent(inout) :: cand_coords(:,:), cand_eigvals(:)
        integer, intent(inout) :: ncand
        real :: mu, sigma, corr
        integer :: i
        if( ncand >= size(cand_eigvals) ) return
        mu = sum(feat) / real(size(feat))
        sigma = sqrt(sum((feat - mu)**2) / real(max(size(feat)-1, 1)))
        if( sigma < 1.e-6 ) return
        do i = 1, ncand
            corr = pearsn(feat, cand_coords(:,i))
            if( abs(corr) > 0.999 ) return
        end do
        ncand = ncand + 1
        cand_coords(:,ncand) = feat
        cand_eigvals(ncand)  = eigval
    end subroutine add_steerable_candidate

    subroutine select_output_candidates(cand_coords, cand_eigvals, ncand, ndiff_scan, coords, eigvals)
        real,                 intent(in)  :: cand_coords(:,:), cand_eigvals(:)
        integer,              intent(in)  :: ncand, ndiff_scan
        real, allocatable,    intent(out) :: coords(:,:), eigvals(:)
        logical, allocatable :: used(:)
        real :: best
        integer :: n, nkeep, k, i, best_i
        n = size(cand_coords,1)
        if( ncand < 1 )then
            allocate(coords(1,n), source=0.)
            allocate(eigvals(1), source=0.)
            return
        endif
        nkeep = min(ndiff_scan, ncand)
        allocate(coords(nkeep,n), source=0.)
        allocate(eigvals(nkeep), source=0.)
        allocate(used(ncand), source=.false.)
        do k = 1, nkeep
            best   = -huge(best)
            best_i = 0
            do i = 1, ncand
                if( used(i) ) cycle
                if( cand_eigvals(i) > best )then
                    best   = cand_eigvals(i)
                    best_i = i
                endif
            end do
            if( best_i == 0 ) exit
            used(best_i) = .true.
            coords(k,:) = cand_coords(:,best_i)
            eigvals(k)  = cand_eigvals(best_i)
        end do
        deallocate(used)
    end subroutine select_output_candidates

    subroutine steerable_transport_denoise(params, imgs, avg, graph_aff, graph_theta, den_ptcls)
        type(parameters),         intent(in)  :: params
        class(image),             intent(in)  :: imgs(:)
        real,                     intent(in)  :: avg(:), graph_aff(:,:), graph_theta(:,:)
        type(image), allocatable, intent(out) :: den_ptcls(:)
        type(image) :: avg_img, acc_img, tmp_img
        real, allocatable :: rmat_rot(:,:,:)
        real :: w, wsum, mean_wsum, angle_deg, denom
        integer :: nptcls, i, j, ldim(3), nedges
        nptcls = size(imgs)
        if( nptcls < 2 ) return
        if( size(graph_aff,1) /= nptcls .or. size(graph_aff,2) /= nptcls .or. &
            size(graph_theta,1) /= nptcls .or. size(graph_theta,2) /= nptcls )then
            write(logfhandle,'(A)') 'Steerable transport denoising skipped: reason=graph shape mismatch'
            call flush(logfhandle)
            return
        endif
        nedges = count(graph_aff > DTINY)
        if( nedges < 1 )then
            write(logfhandle,'(A)') 'Steerable transport denoising skipped: reason=empty steerable graph'
            call flush(logfhandle)
            return
        endif
        ldim = imgs(1)%get_ldim()
        if( size(avg) /= product(ldim) )then
            write(logfhandle,'(A,I8,A,I8)') 'Steerable transport denoising skipped: avg_size=', size(avg), &
                ' expected=', product(ldim)
            call flush(logfhandle)
            return
        endif
        den_ptcls = copy_imgarr(imgs)
        call avg_img%new(ldim, imgs(1)%get_smpd())
        call acc_img%new(ldim, imgs(1)%get_smpd())
        call tmp_img%new(ldim, imgs(1)%get_smpd())
        call avg_img%unserialize(avg)
        allocate(rmat_rot(ldim(1),ldim(2),1), source=0.)
        mean_wsum = 0.
        do i = 1, nptcls
            call acc_img%copy(imgs(i))
            call acc_img%subtr(avg_img)
            wsum = 1.
            do j = 1, nptcls
                if( i == j ) cycle
                w = graph_aff(i,j)
                if( w <= DTINY ) cycle
                if( .not. ieee_is_finite(w) ) cycle
                call tmp_img%copy(imgs(j))
                call tmp_img%subtr(avg_img)
                ! theta(j,i) transports neighbor j into target i's in-plane gauge.
                angle_deg = rad2deg(graph_theta(j,i))
                call tmp_img%rtsq_serial(angle_deg, 0., 0., rmat_rot)
                call tmp_img%set_rmat(rmat_rot, .false.)
                call acc_img%add(tmp_img, w)
                wsum = wsum + w
            end do
            denom = max(wsum, real(DTINY))
            call acc_img%div(denom)
            call den_ptcls(i)%copy(avg_img)
            call den_ptcls(i)%add(acc_img)
            mean_wsum = mean_wsum + wsum
        end do
        mean_wsum = mean_wsum / real(nptcls)

        write(logfhandle,'(A,I8,A,I8,A,F8.3,A,I8)') &
            'Steerable transport denoising: n=', nptcls, ' directed_edges=', nedges, &
            ' mean_wsum=', mean_wsum, ' k_nn=', params%k_nn
        call flush(logfhandle)

        deallocate(rmat_rot)
        call avg_img%kill
        call acc_img%kill
        call tmp_img%kill
    end subroutine steerable_transport_denoise

    subroutine steerable_coeffproj_denoise(params, imgs, avg, graph_aff, graph_theta, coeff_ptcls)
        type(parameters),         intent(in)  :: params
        class(image),             intent(in)  :: imgs(:)
        real,                     intent(in)  :: avg(:), graph_aff(:,:), graph_theta(:,:)
        type(image), allocatable, intent(out) :: coeff_ptcls(:)
        type(parameters) :: params_pft
        type(polarft_calc) :: pftc
        type(image), allocatable :: imgs_work(:)
        type(image) :: avg_img, synth_img
        complex(sp), allocatable :: pfts(:,:,:), pfts_hat(:,:,:), pft_tmp(:,:)
        complex(sp), allocatable :: mode_coeff(:,:), mode_hat(:,:)
        real, allocatable :: deg(:), norm_aff(:,:)
        integer :: nptcls, ldim(3), kfromto(2), pdim_srch(3), pftsz, nrots, nk
        integer :: i, j, k, mode, mode_max, rank_keep, nedges
        real    :: smpd
        nptcls = size(imgs)
        if( nptcls < 2 ) return
        if( size(graph_aff,1) /= nptcls .or. size(graph_aff,2) /= nptcls .or. &
            size(graph_theta,1) /= nptcls .or. size(graph_theta,2) /= nptcls )then
            write(logfhandle,'(A)') 'Steerable coefficient denoising skipped: reason=graph shape mismatch'
            call flush(logfhandle)
            return
        endif
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        if( size(avg) /= product(ldim) )then
            write(logfhandle,'(A,I8,A,I8)') 'Steerable coefficient denoising skipped: avg_size=', size(avg), &
                ' expected=', product(ldim)
            call flush(logfhandle)
            return
        endif
        allocate(deg(nptcls), norm_aff(nptcls,nptcls), source=0.)
        deg = sum(graph_aff, dim=2)
        do i = 1, nptcls
            do j = 1, nptcls
                if( graph_aff(i,j) <= DTINY ) cycle
                if( .not. ieee_is_finite(graph_aff(i,j)) ) cycle
                norm_aff(i,j) = graph_aff(i,j) / sqrt(max(deg(i), DTINY) * max(deg(j), DTINY))
            end do
        end do
        nedges = count(norm_aff > DTINY)
        if( nedges < 1 )then
            write(logfhandle,'(A)') 'Steerable coefficient denoising skipped: reason=empty steerable graph'
            call flush(logfhandle)
            deallocate(deg, norm_aff)
            return
        endif
        call prepare_steerable_pft_params(params, ldim, smpd, params_pft)
        kfromto(1) = max(2, calc_fourier_index(params_pft%hp, params_pft%box, params_pft%smpd))
        kfromto(2) =        calc_fourier_index(params_pft%lp, params_pft%box, params_pft%smpd)
        if( kfromto(2) < kfromto(1) )then
            write(logfhandle,'(A,I8,A,I8)') 'Steerable coefficient denoising skipped: invalid Fourier limits kmin=', &
                kfromto(1), ' kmax=', kfromto(2)
            call flush(logfhandle)
            deallocate(deg, norm_aff)
            return
        endif
        call avg_img%new(ldim, smpd)
        call avg_img%unserialize(avg)
        call pftc%new(params_pft, 1, [1,nptcls], kfromto)
        pdim_srch = pftc%get_pdim_srch()
        pftsz     = pftc%get_pftsz()
        nrots     = pftc%get_nrots()
        nk        = kfromto(2) - kfromto(1) + 1
        mode_max  = min(max(0, params%steerable_nmodes), max(0, pftsz - 1))
        rank_keep = min(max(2, params%k_nn), nptcls)
        allocate(pfts(pftsz,nk,nptcls), pfts_hat(pftsz,nk,nptcls))
        pfts     = cmplx(0.,0.,kind=sp)
        pfts_hat = cmplx(0.,0.,kind=sp)
        allocate(pft_tmp(pftsz,kfromto(1):kfromto(2)))
        imgs_work = copy_imgarr(imgs)
        call imgs_work(1)%memoize4polarize(pdim_srch)
        do i = 1, nptcls
            call imgs_work(i)%subtr(avg_img)
            call imgs_work(i)%fft()
            call pftc%polarize_ptcl_pft(imgs_work(i), i, pdim_srch, oversamp=.false.)
            call pftc%get_ptcl_pft(i, pft_tmp)
            do k = 1, nk
                pfts(:,k,i) = pft_tmp(:,kfromto(1) + k - 1)
            end do
        end do
        call pftc%kill
        call dealloc_imgarr(imgs_work)
        allocate(mode_coeff(nptcls,nk), mode_hat(nptcls,nk))
        pfts_hat = cmplx(0.,0.,kind=sp)
        do mode = 0, mode_max
            call extract_pft_angular_mode(pfts, kfromto, nrots, mode, mode_coeff)
            call project_pft_angular_mode(norm_aff, graph_theta, mode, mode_coeff, mode_hat, rank_keep)
            call accumulate_pft_angular_mode(pfts_hat, kfromto, nrots, mode, mode_hat)
            if( mode > 0 )then
                call extract_pft_angular_mode(pfts, kfromto, nrots, -mode, mode_coeff)
                call project_pft_angular_mode(norm_aff, graph_theta, -mode, mode_coeff, mode_hat, rank_keep)
                call accumulate_pft_angular_mode(pfts_hat, kfromto, nrots, -mode, mode_hat)
            endif
        end do
        coeff_ptcls = copy_imgarr(imgs)
        do i = 1, nptcls
            call synthesize_pft_residual(pfts_hat(:,:,i), kfromto, nrots, ldim, smpd, synth_img)
            call coeff_ptcls(i)%copy(avg_img)
            call coeff_ptcls(i)%add(synth_img)
        end do
        write(logfhandle,'(A,I8,A,I8,A,I8,A,I8,A,I8,A,I8)') &
            'Steerable coefficient denoising: n=', nptcls, ' modes=', mode_max, &
            ' rank=', rank_keep, ' pftsz=', pftsz, ' kmin=', kfromto(1), ' kmax=', kfromto(2)
        call flush(logfhandle)
        if( synth_img%exists() ) call synth_img%kill
        call avg_img%kill
        deallocate(deg, norm_aff, pfts, pfts_hat, pft_tmp, mode_coeff, mode_hat)
    end subroutine steerable_coeffproj_denoise

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
        do i = 1, nptcls
            do k = 1, nk
                acc = cmplx(0.,0.,kind=sp)
                do irot = 1, pftsz
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
        do i = 1, nptcls
            do k = 1, nk
                do irot = 1, pftsz
                    theta     = twopi * real(irot - 1) / real(nrots)
                    phase_arg = real(mode) * theta
                    phase     = cmplx(cos(phase_arg), sin(phase_arg), kind=sp)
                    pfts_hat(irot,k,i) = pfts_hat(irot,k,i) + mode_hat(i,k) * phase
                end do
            end do
        end do
    end subroutine accumulate_pft_angular_mode

    subroutine project_pft_angular_mode(norm_aff, theta, mode, mode_coeff, mode_hat, rank_keep)
        real,        intent(in)  :: norm_aff(:,:), theta(:,:)
        integer,     intent(in)  :: mode, rank_keep
        complex(sp), intent(in)  :: mode_coeff(:,:)
        complex(sp), intent(out) :: mode_hat(:,:)
        real, allocatable :: mat(:,:), evals(:), evecs(:,:), y(:,:), yhat(:,:), coeff_r(:)
        complex(sp), allocatable :: coeff_c(:)
        integer :: n, n2, nk, nkeep, i, a
        n  = size(mode_coeff,1)
        nk = size(mode_coeff,2)
        mode_hat = cmplx(0.,0.,kind=sp)
        if( n < 1 .or. nk < 1 ) return
        if( mode == 0 )then
            nkeep = min(max(1, rank_keep), n)
            allocate(mat(n,n), evals(n), evecs(n,n), coeff_c(nk))
            mat = norm_aff
            call eigh(n, mat, nkeep, evals, evecs)
            do a = 1, nkeep
                coeff_c = cmplx(0.,0.,kind=sp)
                do i = 1, n
                    coeff_c = coeff_c + real(evecs(i,a), kind=sp) * mode_coeff(i,:)
                end do
                do i = 1, n
                    mode_hat(i,:) = mode_hat(i,:) + real(evecs(i,a), kind=sp) * coeff_c
                end do
            end do
            deallocate(mat, evals, evecs, coeff_c)
        else
            n2    = 2 * n
            nkeep = min(max(2, 2 * rank_keep), n2)
            allocate(mat(n2,n2), evals(n2), evecs(n2,n2), y(n2,nk), yhat(n2,nk), coeff_r(nk))
            call build_connection_block(norm_aff, theta, mode, mat)
            call eigh(n2, mat, nkeep, evals, evecs)
            do i = 1, n
                y(i,:)   = real(mode_coeff(i,:))
                y(n+i,:) = aimag(mode_coeff(i,:))
            end do
            yhat = 0.
            do a = 1, nkeep
                coeff_r = 0.
                do i = 1, n2
                    coeff_r = coeff_r + evecs(i,a) * y(i,:)
                end do
                do i = 1, n2
                    yhat(i,:) = yhat(i,:) + evecs(i,a) * coeff_r
                end do
            end do
            do i = 1, n
                mode_hat(i,:) = cmplx(yhat(i,:), yhat(n+i,:), kind=sp)
            end do
            deallocate(mat, evals, evecs, y, yhat, coeff_r)
        endif
    end subroutine project_pft_angular_mode

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
        do irot = 1, pftsz
            theta = twopi * real(irot - 1) / real(nrots)
            do k = 1, nk
                krad = kfromto(1) + k - 1
                x =  sin(theta) * real(krad)
                y = -cos(theta) * real(krad)
                call scatter_fourier_sample(img_out, rho, lims, x, y, pft_half(irot,k))
            end do
        end do
        cmat = img_out%get_cmat()
        do l = 1, ashape(3)
            do k = 1, ashape(2)
                do h = 1, ashape(1)
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
        real    :: dh, dk, wh(2), wk(2), w
        if( .not. ieee_is_finite(real(comp)) .or. .not. ieee_is_finite(aimag(comp)) ) return
        h0 = floor(x)
        k0 = floor(y)
        h1 = h0 + 1
        k1 = k0 + 1
        dh = x - real(h0)
        dk = y - real(k0)
        wh = [1. - dh, dh]
        wk = [1. - dk, dk]
        do ih = 1, 2
            h = merge(h0, h1, ih == 1)
            if( h < lims(1,1) .or. h > lims(1,2) ) cycle
            do ik = 1, 2
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

end module simple_diffusion_maps