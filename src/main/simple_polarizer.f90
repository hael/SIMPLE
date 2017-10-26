! polar 2D Fourier transform generation by convolution interpolation (gridding)

module simple_polarizer
!$ use omp_lib
!$ use omp_lib_kinds
#include "simple_lib.f08"
    
use simple_kbinterpol, only: kbinterpol
use simple_image,      only: image
use simple_params,     only: params
implicit none

public :: polarizer
private

type, extends(image) :: polarizer
    private
    type(kbinterpol)      :: kbwin                 !< window function object
    real,    allocatable  :: polweights_mat(:,:,:) !< polar weights matrix for the image to polar transformer
    integer, allocatable  :: polcyc1_mat(:,:,:)    !< image cyclic adresses for the image to polar transformer
    integer, allocatable  :: polcyc2_mat(:,:,:)    !< image cyclic adresses for the image to polar transformer
    integer               :: wdim = 0              !< dimension of K-B window
    integer               :: pdim(3) = 0           ! Polar-FT matrix dimensions
    complex, allocatable  :: pft(:,:), comps(:,:)
  contains
    procedure :: init_polarizer
    procedure :: polarize
    procedure :: kill_polarizer
end type polarizer

contains

    ! IMAGE TO POLAR FT TRANSFORMER

    !> \brief  initialises the image polarizer
    subroutine init_polarizer( self, pftcc, alpha )
        use simple_math,             only: sqwin_2d, cyci_1d
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(polarizer),        intent(inout) :: self   !< projector instance
        class(polarft_corrcalc), intent(inout) :: pftcc  !< polarft_corrcalc object to be filled
        real,                    intent(in)    :: alpha  !< oversampling factor
        real, allocatable :: w(:,:)
        real              :: loc(2)
        integer           :: win(2,2), lims(3,2)
        integer           :: i, k, l, wlen, cnt
        if( .not. pftcc%exists() ) stop 'polarft_corrcalc object needs to be created; init_imgpolarizer; simple_projector'
        call self%kill_polarizer
        self%kbwin = kbinterpol(KBWINSZ, alpha)
        self%wdim  = self%kbwin%get_wdim()
        wlen       = self%wdim**2
        self%pdim  = pftcc%get_pdim()
        lims       = self%loop_lims(3)
        allocate( self%polcyc1_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wdim),&
                  &self%polcyc2_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wdim),&
                  &self%polweights_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:wlen),&
                  &w(1:self%wdim,1:self%wdim), self%comps(1:self%wdim,1:self%wdim),&
                  &self%pft(self%pdim(1),self%pdim(2):self%pdim(3)), stat=alloc_stat)
        allocchk('in simple_projector :: init_imgpolarizer')
        !$omp parallel do collapse(2) schedule(static) default(shared)&
        !$omp private(i,k,l,w,loc,cnt,win) proc_bind(close)
        do i=1,self%pdim(1)
            do k=self%pdim(2),self%pdim(3)
                ! polar coordinates
                loc = pftcc%get_coord(i,k)
                call sqwin_2d(loc(1), loc(2), KBWINSZ, win)
                w   = 1.
                cnt = 0
                do l=1,self%wdim
                    cnt = cnt + 1
                    ! interpolation weights
                    w(l,:) = w(l,:) * self%kbwin%apod( real(win(1,1)+l-1)-loc(1) )
                    w(:,l) = w(:,l) * self%kbwin%apod( real(win(2,1)+l-1)-loc(2) )
                    ! cyclic addresses
                    self%polcyc1_mat(i, k, cnt) = cyci_1d(lims(1,:), win(1,1)+l-1)
                    self%polcyc2_mat(i, k, cnt) = cyci_1d(lims(2,:), win(2,1)+l-1)
                end do
                self%polweights_mat(i,k,:) = reshape(w,(/wlen/))
            enddo
        enddo
        !$omp end parallel do
        deallocate(w)
    end subroutine init_polarizer

    !> \brief  creates the polar Fourier transform
    subroutine polarize( self, pftcc, img_ind, isptcl, iseven )
        use simple_math, only: sqwin_2d, cyci_1d
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(polarizer),        intent(inout) :: self    !< projector instance
        class(polarft_corrcalc), intent(inout) :: pftcc   !< polarft_corrcalc object to be filled
        integer,                 intent(in)    :: img_ind !< image index
        logical,                 intent(in)    :: isptcl  !< is ptcl (or reference)
        logical,                 intent(in)    :: iseven  !< is even (or odd)
        integer :: logi(3), phys(3), i, k, l, m, vecdim, addr_l
        vecdim   = self%wdim**2
        !!$omp parallel do collapse(2) schedule(static) default(shared)&
        !!$omp private(i,k,l,m,logi,phys,comps,addr_l) proc_bind(close)
        do i=1,self%pdim(1)
            do k=self%pdim(2),self%pdim(3)
                do l=1,self%wdim
                    addr_l = self%polcyc1_mat(i,k,l)
                    do m=1,self%wdim
                        logi = [addr_l,self%polcyc2_mat(i,k,m),0]
                        phys = self%comp_addr_phys(logi)
                        self%comps(l,m) = self%get_fcomp(logi,phys)
                    enddo
                enddo
                self%pft(i,k) = dot_product(self%polweights_mat(i,k,:), reshape(self%comps,(/vecdim/)))
            end do
        end do
        !!$omp end parallel do
        if( isptcl )then
            call pftcc%set_ptcl_pft(img_ind, self%pft)
        else
            call pftcc%set_ref_pft(img_ind, self%pft, iseven)
        endif
    end subroutine polarize

    ! DESTRUCTOR

    !>  \brief  is a destructor of impolarizer
    subroutine kill_polarizer( self )
        class(polarizer), intent(inout) :: self !< projector instance
        if( allocated(self%polweights_mat) ) deallocate(self%polweights_mat)
        if( allocated(self%polcyc1_mat)    ) deallocate(self%polcyc1_mat)
        if( allocated(self%polcyc2_mat)    ) deallocate(self%polcyc2_mat)
        if( allocated(self%pft)            ) deallocate(self%pft)
        if( allocated(self%comps)          ) deallocate(self%comps)
    end subroutine kill_polarizer

end module simple_polarizer
