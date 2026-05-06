!@descr: polar 2D Fourier transform generation by convolution interpolation (gridding)
submodule (simple_image) simple_image_polar
implicit none
#include "simple_local_flags.inc"

integer, parameter :: POLMEM_NORMAL = 1, POLMEM_OVERSAMP = 2

contains

    !> \brief  initialises the image polarizer
    module subroutine memoize4polarize( self, pdim )
        class(image), intent(in) :: self    !< instance
        integer,      intent(in) :: pdim(3) !< pftsz,kfrom,kto
        type(kbinterpol)  :: kbwin          !< KB kernel object
        real, allocatable :: wx(:), wy(:), sinang(:), negcosang(:)
        real              :: x, y, dang
        integer           :: win_lo(2), lims(2,3), i, k, l, m, ind, wdim, iwinsz
        if( polar_memo_valid(self, pdim, POLMEM_NORMAL) ) return
        call clear_polar_memo
        mem_poldim   = pdim
        mem_polldim  = self%ldim
        mem_polmode  = POLMEM_NORMAL
        lims         = transpose(self%loop_lims(3)) ! fortran layered memory
        kbwin        = kbinterpol(KBWINSZ, 1.0)     ! no oversampling, since self is not padded
        wdim         = kbwin%get_wdim()
        iwinsz       = ceiling(kbwin%get_winsz() - 0.5)
        mem_polwdim  = wdim
        mem_polwlen  = wdim**2
        allocate( mem_polph_mat(      1:mem_polwlen, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)),&
                  &mem_polk_mat(      1:mem_polwlen, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)),&
                  &mem_polconjg_mat(  1:mem_polwlen, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)),&
                  &mem_polweights_mat(1:mem_polwlen, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)),&
                  &sinang(1:mem_poldim(1)), negcosang(1:mem_poldim(1)))
        dang = twopi / real(2 * mem_poldim(1))
        do i = 1, mem_poldim(1)
            sinang(i)    =  sin(real(i-1) * dang)
            negcosang(i) = -cos(real(i-1) * dang)
        end do
        !$omp parallel default(shared) private(i,k,x,y,win_lo,l,m,ind,wx,wy) proc_bind(close)
        allocate(wx(1:wdim), wy(1:wdim))
        !$omp do collapse(2) schedule(static)
        do i = 1, mem_poldim(1)
            do k = mem_poldim(2), mem_poldim(3)
                x = sinang(i)    * real(k)
                y = negcosang(i) * real(k)
                win_lo = nint([x,y]) - iwinsz
                do l = 1, wdim
                    wx(l) = kbwin%apod(real(win_lo(1)+l-1) - x)
                    wy(l) = kbwin%apod(real(win_lo(2)+l-1) - y)
                end do
                wx = wx * (1.0 / sum(wx))
                wy = wy * (1.0 / sum(wy))
                ind = 0
                do m = 1, wdim
                    do l = 1, wdim
                        ind = ind + 1
                        call memoize_polar_point(ind, i, k, wx(l) * wy(m),&
                            &cyci_1d(lims(:,1), win_lo(1)+l-1),&
                            &cyci_1d(lims(:,2), win_lo(2)+m-1), self%ldim(2))
                    end do
                end do
            enddo
        enddo
        !$omp end do
        deallocate(wx, wy)
        !$omp end parallel
        deallocate(sinang, negcosang)
    end subroutine memoize4polarize

    module subroutine memoize4polarize_oversamp( self, pdim )
        class(image), intent(in) :: self
        integer,      intent(in) :: pdim(3)
        type(kbinterpol)  :: kbwin
        real, allocatable :: wx(:), wy(:), sinang(:), negcosang(:)
        real    :: xpd, ypd, dang
        integer :: win_lo(2), lims(2,3), i, k, l, m, ind, pf, wdim, iwinsz
        if( polar_memo_valid(self, pdim, POLMEM_OVERSAMP) ) return
        call clear_polar_memo
        mem_poldim  = pdim
        mem_polldim = self%ldim
        mem_polmode = POLMEM_OVERSAMP
        pf          = OSMPL_PAD_FAC
        ! IMPORTANT: lims from PADDED grid because self is padded
        lims = transpose(self%loop_lims(3))
        ! KB kernel for padded-grid interpolation
        kbwin       = kbinterpol(KBWINSZ, KBALPHA)
        wdim        = kbwin%get_wdim()
        iwinsz      = ceiling(kbwin%get_winsz() - 0.5)
        mem_polwdim = wdim
        mem_polwlen = wdim*wdim
        allocate( mem_polph_mat(      1:mem_polwlen, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)),&
                  &mem_polk_mat(     1:mem_polwlen, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)),&
                  &mem_polconjg_mat( 1:mem_polwlen, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)),&
                  &mem_polweights_mat(1:mem_polwlen, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)),&
                  &sinang(1:mem_poldim(1)), negcosang(1:mem_poldim(1)))
        dang = twopi / real(2 * mem_poldim(1))
        do i = 1, mem_poldim(1)
            sinang(i)    =  sin(real(i-1) * dang)
            negcosang(i) = -cos(real(i-1) * dang)
        end do
        !$omp parallel default(shared) private(i,k,xpd,ypd,win_lo,l,m,ind,wx,wy) proc_bind(close)
        allocate(wx(1:wdim), wy(1:wdim))
        !$omp do collapse(2) schedule(static)
        do i = 1, mem_poldim(1)
            do k = mem_poldim(2), mem_poldim(3)
                xpd = sinang(i)    * real(pf*k)
                ypd = negcosang(i) * real(pf*k)
                win_lo = nint([xpd,ypd]) - iwinsz
                do l = 1, wdim
                    wx(l) = kbwin%apod_fast(real(win_lo(1)+l-1) - xpd)
                    wy(l) = kbwin%apod_fast(real(win_lo(2)+l-1) - ypd)
                end do
                wx = wx * (1.0 / sum(wx))
                wy = wy * (1.0 / sum(wy))
                ind = 0
                do m = 1, wdim
                    do l = 1, wdim
                        ind = ind + 1
                        call memoize_polar_point(ind, i, k, wx(l) * wy(m),&
                            &cyci_1d(lims(:,1), win_lo(1)+l-1),&
                            &cyci_1d(lims(:,2), win_lo(2)+m-1), self%ldim(1))
                    end do
                end do
            end do
        end do
        !$omp end do
        deallocate(wx, wy)
        !$omp end parallel
        deallocate(sinang, negcosang)
    end subroutine memoize4polarize_oversamp

    ! keep serial
    module subroutine polarize( self, pft )
        class(image),      intent(in)    :: self     !< image instance to polarize
        complex,           intent(inout) :: pft(mem_poldim(1),mem_poldim(2):mem_poldim(3)) !< polarft to be filled
        complex(kind=c_float_complex) :: acc, fcomp
        integer :: i, k, ind
        ! interpolate
        !$OMP SIMD COLLAPSE(2) PRIVATE(i,k,acc,ind,fcomp)
        do k = mem_poldim(2), mem_poldim(3)
            do i = 1, mem_poldim(1)
                acc = CMPLX_ZERO
                do ind = 1, mem_polwlen
                    fcomp = self%cmat(mem_polph_mat(ind,i,k), mem_polk_mat(ind,i,k), 1)
                    if( mem_polconjg_mat(ind,i,k) ) fcomp = conjg(fcomp)
                    acc = acc + mem_polweights_mat(ind,i,k) * fcomp
                end do
                pft(i,k) = acc
            end do
        end do
        !$OMP END SIMD
    end subroutine polarize

    module subroutine polarize_oversamp( self, pft )
        class(image),      intent(in)    :: self
        complex,           intent(inout) :: pft(mem_poldim(1),mem_poldim(2):mem_poldim(3))
        complex(kind=c_float_complex) :: acc, fcomp
        real    :: padding_factor_scaling
        integer :: i, k, ind
        padding_factor_scaling = real(OSMPL_PAD_FAC**2)
        !$OMP SIMD COLLAPSE(2) PRIVATE(i,k,acc,ind,fcomp)
        do k = mem_poldim(2), mem_poldim(3)
            do i = 1, mem_poldim(1)
                acc = CMPLX_ZERO
                do ind = 1, mem_polwlen
                    fcomp = self%cmat(mem_polph_mat(ind,i,k), mem_polk_mat(ind,i,k), 1)
                    if( mem_polconjg_mat(ind,i,k) ) fcomp = conjg(fcomp)
                    acc = acc + mem_polweights_mat(ind,i,k) * fcomp
                end do
                pft(i,k) = acc * padding_factor_scaling
            end do
        end do
        !$OMP END SIMD
    end subroutine polarize_oversamp

    pure logical function polar_memo_valid( self, pdim, mode )
        class(image), intent(in) :: self
        integer,      intent(in) :: pdim(3), mode
        polar_memo_valid = allocated(mem_polweights_mat) .and. allocated(mem_polph_mat) .and.&
            &allocated(mem_polk_mat) .and. allocated(mem_polconjg_mat)
        if( .not. polar_memo_valid ) return
        polar_memo_valid = all(mem_poldim == pdim) .and. all(mem_polldim == self%ldim) .and. mem_polmode == mode
    end function polar_memo_valid

    subroutine clear_polar_memo
        if( allocated(mem_polweights_mat) ) deallocate(mem_polweights_mat)
        if( allocated(mem_polph_mat)      ) deallocate(mem_polph_mat)
        if( allocated(mem_polk_mat)       ) deallocate(mem_polk_mat)
        if( allocated(mem_polconjg_mat)   ) deallocate(mem_polconjg_mat)
        mem_poldim  = 0
        mem_polldim = 0
        mem_polmode = 0
        mem_polwdim = 0
        mem_polwlen = 0
    end subroutine clear_polar_memo

    subroutine memoize_polar_point( ind, i, k, weight, h_val, k_val, kdim )
        integer, intent(in) :: ind, i, k, h_val, k_val, kdim
        real,    intent(in) :: weight
        integer :: k_eff
        mem_polweights_mat(ind,i,k) = weight
        mem_polconjg_mat(ind,i,k)   = h_val < 0
        mem_polph_mat(ind,i,k)      = abs(h_val) + 1
        k_eff = merge(-k_val, k_val, mem_polconjg_mat(ind,i,k))
        mem_polk_mat(ind,i,k) = k_eff + 1
        if( k_eff < 0 ) mem_polk_mat(ind,i,k) = mem_polk_mat(ind,i,k) + kdim
    end subroutine memoize_polar_point

end submodule simple_image_polar
