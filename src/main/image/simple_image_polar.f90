!@descr: polar 2D Fourier transform generation by convolution interpolation (gridding)
submodule (simple_image) simple_image_polar
implicit none
#include "simple_local_flags.inc"

contains

    !> \brief  initialises the image polarizer
    module subroutine memoize4polarize( self, pdim, instrfun_img )
        use simple_gridding, only: gen_instrfun_img
        class(image),           intent(in)    :: self         !< instance
        integer,                intent(in)    :: pdim(3)      !< pftsz,kfrom,kto
        class(image), optional, intent(inout) :: instrfun_img !< instrument function
        type(kbinterpol)  :: kbwin                            !< KB kernel  object
        real, allocatable :: w(:,:)
        real              :: x, y, d1, d2, dang, ang
        integer           :: win(2,2), lims(2,3), i, k, l, cnt, f1, f2
        if( allocated(mem_polweights_mat) ) deallocate(mem_polweights_mat)
        if( allocated(mem_polcyc1_mat)    ) deallocate(mem_polcyc1_mat)
        if( allocated(mem_polcyc2_mat)    ) deallocate(mem_polcyc2_mat)
        mem_poldim   = pdim
        lims         = transpose(self%loop_lims(3)) ! fortran layered memory
        kbwin        = kbinterpol(KBWINSZ, 1.0)     ! no oversampling
        mem_polwdim  = kbwin%get_wdim()
        mem_polwlen  = mem_polwdim**2
        allocate( mem_polcyc1_mat(    1:mem_polwdim, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)),&
                  &mem_polcyc2_mat(   1:mem_polwdim, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)),&
                  &mem_polweights_mat(1:mem_polwlen, 1:mem_poldim(1), mem_poldim(2):mem_poldim(3)),&
                  &w(1:mem_polwdim,1:mem_polwdim))
        ! cartesian to polar with Kaiser-Bessel
        dang = twopi/real(2 * mem_poldim(1))
        !$omp parallel do collapse(2) schedule(static) private(i,k,ang,x,y,cnt,l,w,win) default(shared) proc_bind(close)
        do i=1,mem_poldim(1)
            do k=mem_poldim(2),mem_poldim(3)
                ang = real(i-1) * dang
                ! polar coordinates
                x =  sin(ang)*real(k) ! x-coordinate
                y = -cos(ang)*real(k) ! y-coordinate
                call sqwin_2d(x, y, kbwin%get_winsz(), win)
                w   = 1.
                cnt = 0
                do l=1,mem_polwdim
                    cnt = cnt + 1
                    ! interpolation weights
                    w(l,:) = w(l,:) * kbwin%apod(real(win(1,1)+l-1)-x)
                    w(:,l) = w(:,l) * kbwin%apod(real(win(2,1)+l-1)-y)
                    ! cyclic addresses
                    mem_polcyc1_mat(cnt, i, k) = cyci_1d(lims(:,1), win(1,1)+l-1)
                    mem_polcyc2_mat(cnt, i, k) = cyci_1d(lims(:,2), win(2,1)+l-1)
                end do
                mem_polweights_mat(:,i,k) = reshape(w,(/mem_polwlen/))
                mem_polweights_mat(:,i,k) = mem_polweights_mat(:,i,k) / sum(w)
            enddo
        enddo
        !$omp end parallel do
        ! memoize instrument function
        if( present(instrfun_img) )then
            call instrfun_img%new(self%ldim, self%smpd, self%wthreads)
            call gen_instrfun_img(instrfun_img, kbwin=kbwin)
        endif
        deallocate(w)
    end subroutine memoize4polarize

    ! keep serial
    module subroutine polarize( self, pft, mask )
        class(image),      intent(in)    :: self     !< image instance to polarize
        complex,           intent(inout) :: pft(mem_poldim(1),mem_poldim(2):mem_poldim(3)) !< polarft to be filled
        logical, optional, intent(in)    :: mask(:)  !< interpolation mask, all .false. set to CMPLX_ZERO
        complex(kind=c_float_complex) :: acc, fcomp
        logical :: h_negative
        integer :: i, k, l, m, ind, h_val, k_val, phys1, phys2
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

end submodule simple_image_polar
