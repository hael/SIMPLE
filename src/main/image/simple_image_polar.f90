!@descr: polar 2D Fourier transform generation by convolution interpolation (gridding)
submodule (simple_image) simple_image_polar
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine memoize4polarize( self, pdim, instrfun_img )
        class(image), intent(in) :: self
        integer,      intent(in) :: pdim(3)
        class(image), optional, intent(inout) :: instrfun_img
        type(kbinterpol)  :: kbwin
        real, allocatable :: w(:,:)
        real :: xpd, ypd, dang, ang
        integer :: win(2,2), lims(2,3), i, k, l, m, pf, iwinsz, wdim
        pf = STRIDE_GRID_PAD_FAC
        if( allocated(mem_polweights_mat) ) deallocate(mem_polweights_mat)
        if( allocated(mem_polcyc1_mat)    ) deallocate(mem_polcyc1_mat)
        if( allocated(mem_polcyc2_mat)    ) deallocate(mem_polcyc2_mat)
        mem_poldim = pdim
        ! IMPORTANT: lims from PADDED grid because self is padded
        lims = transpose(self%loop_lims(3))
        ! KB kernel for padded-grid interpolation
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        wdim   = kbwin%get_wdim()
        iwinsz = ceiling(kbwin%get_winsz() - 0.5)
        mem_polwdim = wdim
        mem_polwlen = wdim*wdim
        allocate( mem_polcyc1_mat(  1:wdim,        1:mem_poldim(1), mem_poldim(2):mem_poldim(3)), &
                mem_polcyc2_mat(    1:wdim,        1:mem_poldim(1), mem_poldim(2):mem_poldim(3)), &
                mem_polweights_mat( 1:mem_polwlen, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)), &
                w(1:wdim,1:wdim) )
        dang = twopi / real(2 * mem_poldim(1))
        !$omp parallel do collapse(2) schedule(static) private(i,k,ang,xpd,ypd,win,w,l,m) default(shared) proc_bind(close)
        do i = 1, mem_poldim(1)
            do k = mem_poldim(2), mem_poldim(3)
                ang = real(i-1) * dang
                ! Sample location in PADDED logical units
                xpd =  sin(ang) * real(pf*k)
                ypd = -cos(ang) * real(pf*k)
                call sqwin_2d(xpd, ypd, kbwin%get_winsz(), win)
                ! Store cyclic neighbor indices (PADDED logical)
                do l = 1, wdim
                    mem_polcyc1_mat(l,i,k) = cyci_1d(lims(:,1), win(1,1) + l - 1)
                    mem_polcyc2_mat(l,i,k) = cyci_1d(lims(:,2), win(2,1) + l - 1)
                end do
                ! Separable weights (PADDED geometry)
                w = 1.0
                do l = 1, wdim
                    w(l,:) = w(l,:) * kbwin%apod(real(win(1,1)+l-1) - xpd)
                    w(:,l) = w(:,l) * kbwin%apod(real(win(2,1)+l-1) - ypd)
                end do
                mem_polweights_mat(:,i,k) = reshape(w, (/mem_polwlen/))
                mem_polweights_mat(:,i,k) = mem_polweights_mat(:,i,k) / sum(w)
            end do
        end do
        !$omp end parallel do
        deallocate(w)
    end subroutine memoize4polarize

    ! keep serial
    module subroutine polarize( self, pft, mask )
        class(image),      intent(in)    :: self     !< image instance to polarize
        complex,           intent(inout) :: pft(mem_poldim(1),mem_poldim(2):mem_poldim(3)) !< polarft to be filled
        logical, optional, intent(in)    :: mask(:)  !< interpolation mask, all .false. set to CMPLX_ZERO
        complex(kind=c_float_complex) :: acc, fcomp
        logical :: h_negative
        integer :: i, k, l, m, ind, h_val, k_val, phys1, phys2, ithr
        ! interpolate
        !$OMP SIMD COLLAPSE(2) PRIVATE(i,k,acc,ind,m,l,h_val,k_val,phys1,phys2,h_negative,fcomp)
        do k=mem_poldim(2),mem_poldim(3)
            do i=1,mem_poldim(1)
                acc = CMPLX_ZERO
                ind = 0
                do m = 1, mem_polwdim
                    k_val = mem_polcyc2_mat(m,i,k)
                    do l = 1, mem_polwdim
                        ind         = ind + 1
                        ! Get h and k values
                        h_val       = mem_polcyc1_mat(l,i,k)
                        ! Branch-free indexing computation
                        h_negative  = (h_val < 0)
                        phys1       = merge(-h_val, h_val, h_negative) + 1
                        phys2       = merge(-k_val, k_val, h_negative) + 1 + merge(self%ldim(2), 0, merge(-k_val, k_val, h_negative) < 0)
                        ! Fetch complex value
                        fcomp       = merge(conjg(self%cmat(phys1,phys2,1)), self%cmat(phys1,phys2,1), h_negative)
                        ! accumulate dot product
                        acc         = acc + mem_polweights_mat(ind,i,k) * fcomp
                    end do
                end do
                pft(i,k) = acc
            end do
        end do
        !$OMP END SIMD
        if( present(mask) )then
            ! band masking
            !$OMP SIMD
            do k=mem_poldim(2),mem_poldim(3)
                if( .not.mask(k) ) pft(:,k) = CMPLX_ZERO
            enddo
            !$OMP END SIMD
        endif
    end subroutine polarize

    module subroutine polarize_oversamp( self, pft, mask )
        class(image),      intent(in)    :: self
        complex,           intent(inout) :: pft(mem_poldim(1),mem_poldim(2):mem_poldim(3))
        logical, optional, intent(in)    :: mask(:)
        complex(kind=c_float_complex) :: acc, fcomp
        real    :: padding_factor_scaling
        logical :: h_negative
        integer :: padded_box
        integer :: i, k, l, m, ind
        integer :: h_val_pd, k_val_pd
        integer :: h_abs, k_eff
        integer :: phys1p, phys2p
        padded_box = self%ldim(1)
        padding_factor_scaling = real(STRIDE_GRID_PAD_FAC**2)
        !$OMP SIMD COLLAPSE(2) PRIVATE(i,k,acc,ind,m,l,h_val_pd,k_val_pd,h_abs,k_eff,phys1p,phys2p,h_negative,fcomp)
        do k = mem_poldim(2), mem_poldim(3)
            do i = 1, mem_poldim(1)
                acc = CMPLX_ZERO
                ind = 0
                do m = 1, mem_polwdim
                    k_val_pd = mem_polcyc2_mat(m,i,k)
                    do l = 1, mem_polwdim
                        ind = ind + 1
                        h_val_pd = mem_polcyc1_mat(l,i,k)
                        h_negative = (h_val_pd < 0)
                        h_abs = abs(h_val_pd)
                        k_eff = merge(-k_val_pd, k_val_pd, h_negative)
                        phys1p = h_abs + 1
                        phys2p = k_eff + 1
                        if (k_eff < 0) phys2p = phys2p + padded_box
                        fcomp = merge(conjg(self%cmat(phys1p,phys2p,1)), self%cmat(phys1p,phys2p,1), h_negative)
                        acc = acc + mem_polweights_mat(ind,i,k) * fcomp
                    end do
                end do
                pft(i,k) = acc * padding_factor_scaling
            end do
        end do
        !$OMP END SIMD
        if( present(mask) )then
            !$OMP SIMD
            do k = mem_poldim(2), mem_poldim(3)
                if( .not.mask(k) ) pft(:,k) = CMPLX_ZERO
            end do
            !$OMP END SIMD
        endif
    end subroutine polarize_oversamp

end submodule simple_image_polar
