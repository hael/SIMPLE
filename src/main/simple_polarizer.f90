! polar 2D Fourier transform generation by convolution interpolation (gridding)
module simple_polarizer
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: polarizer
private
#include "simple_local_flags.inc"

complex, parameter :: CMPLX_ZERO = cmplx(0.,0.)

type, extends(image) :: polarizer
    private
    type(image)           :: instrfun_img          !< weights to divide with in real-space
    complex, allocatable  :: pft(:,:)              !< Polar-FT matrix
    complex, allocatable  :: comps(:,:)            !< pre-allocated for performance kernel Fourier components
    real,    allocatable  :: polweights_mat(:,:,:) !< polar weights matrix for the image to polar transformer
    real,    allocatable  :: gridweights(:,:)      !< weights to divide with in real-space
    integer, allocatable  :: polcyc1_mat(:,:,:)    !< image cyclic adresses for the image to polar transformer
    integer, allocatable  :: polcyc2_mat(:,:,:)    !< image cyclic adresses for the image to polar transformer
    integer               :: wdim = 0              !< dimension of K-B window
    integer               :: wlen = 0              !< dimension squared of K-B window
    integer               :: pdim(3) = 0           !< Polar-FT matrix dimensions
  contains
    procedure :: init_polarizer
    procedure :: copy_polarizer
    procedure :: div_by_instrfun
    procedure :: polarize
    procedure :: kill_polarizer
end type polarizer

contains

    ! IMAGE TO POLAR FT TRANSFORMER

    !> \brief  initialises the image polarizer
    subroutine init_polarizer( self, pftcc, alpha )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_kbinterpol,       only: kbinterpol
        use simple_gridding,         only: gen_instrfun_img
        class(polarizer),        intent(inout) :: self   !< projector instance
        class(polarft_corrcalc), intent(inout) :: pftcc  !< polarft_corrcalc object to be filled
        real,                    intent(in)    :: alpha  !< oversampling factor
        type(kbinterpol)  :: kbwin                 !< window function object
        real, allocatable :: w(:,:)
        real              :: loc(2)
        integer           :: win(2,2), lims(3,2), i, k, l, cnt
        if( .not. pftcc%exists() ) THROW_HARD('polarft_corrcalc object needs to be created; init_polarizer')
        call self%kill_polarizer
        self%pdim  = pftcc%get_pdim()
        kbwin      = kbinterpol(KBWINSZ, alpha)
        self%wdim  = kbwin%get_wdim()
        self%wlen  = self%wdim**2
        lims       = self%loop_lims(3)
        allocate( self%polcyc1_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wdim),&
                  &self%polcyc2_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wdim),&
                  &self%polweights_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wlen),&
                  &w(1:self%wdim,1:self%wdim), self%comps(1:self%wdim,1:self%wdim),&
                  &self%pft(self%pdim(1),self%pdim(2):self%pdim(3)), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('in simple_projector :: init_imgpolarizer',alloc_stat)
        ! instrument function
        call self%instrfun_img%new(self%get_ldim(), self%get_smpd())
        call gen_instrfun_img(self%instrfun_img, kbwin)
        ! cartesian to polar
        !$omp parallel do collapse(2) schedule(static) private(i,k,l,w,loc,cnt,win) default(shared) proc_bind(close)
        do i=1,self%pdim(1)
            do k=self%pdim(2),self%pdim(3)
                ! polar coordinates
                loc = pftcc%get_coord(i,k)
                call sqwin_2d(loc(1), loc(2), kbwin%get_winsz(), win)
                w   = 1.
                cnt = 0
                do l=1,self%wdim
                    cnt = cnt + 1
                    ! interpolation weights
                    w(l,:) = w(l,:) * kbwin%apod( real(win(1,1)+l-1)-loc(1) )
                    w(:,l) = w(:,l) * kbwin%apod( real(win(2,1)+l-1)-loc(2) )
                    ! cyclic addresses
                    self%polcyc1_mat(i, k, cnt) = cyci_1d(lims(1,:), win(1,1)+l-1)
                    self%polcyc2_mat(i, k, cnt) = cyci_1d(lims(2,:), win(2,1)+l-1)
                end do
                self%polweights_mat(i,k,:) = reshape(w,(/self%wlen/))
            enddo
        enddo
        !$omp end parallel do
        deallocate(w)
    end subroutine init_polarizer

    subroutine copy_polarizer(self, self_in)
        class(polarizer), intent(inout) :: self    !< projector instance
        class(polarizer), intent(inout) :: self_in
        call self%kill_polarizer
        self%pdim = self_in%pdim
        self%wdim = self_in%wdim
        self%wlen = self_in%wlen
        allocate( self%polcyc1_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wdim),    source=self_in%polcyc1_mat )
        allocate( self%polcyc2_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wdim),    source=self_in%polcyc2_mat )
        allocate( self%polweights_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wlen), source=self_in%polweights_mat )
        allocate( self%comps(1:self%wdim,1:self%wdim),                                         source=CMPLX_ZERO )
        allocate( self%pft(self%pdim(1),self%pdim(2):self%pdim(3)),                            source=CMPLX_ZERO )
        call self%instrfun_img%copy(self_in%instrfun_img)
    end subroutine copy_polarizer

    !> \brief  divide by gridding weights in real-space prior to FFT &
    !>         change to polar coordinate system, keep serial
    subroutine div_by_instrfun( self )
        class(polarizer), intent(inout) :: self
        if( self%is_ft() ) THROW_HARD('Polarizer must be in real-space in div_by_gridweights')
        if( .not.(self.eqdims.self%instrfun_img) ) THROW_HARD('Incompatible dimensions in div_by_gridweights')
        call self%div(self%instrfun_img)
    end subroutine div_by_instrfun

    !> \brief  creates the polar Fourier transform
    !!         KEEP THIS ROUTINE SERIAL
    subroutine polarize( self, pftcc, img_ind, isptcl, iseven, mask )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(polarizer),        intent(inout) :: self    !< projector instance
        class(polarft_corrcalc), intent(inout) :: pftcc   !< polarft_corrcalc object to be filled
        integer,                 intent(in)    :: img_ind !< image index
        logical,                 intent(in)    :: isptcl  !< is ptcl (or reference)
        logical,                 intent(in)    :: iseven  !< is even (or odd)
        logical,                 intent(in)    :: mask(:) !< interpolation mask, all .false. set to CMPLX_ZERO
        integer :: logi(3), phys(3), i, k, l, m, addr_l
        do i=1,self%pdim(1)
            do k=self%pdim(2),self%pdim(3)
                if( mask(k) )then
                    do l=1,self%wdim
                        addr_l = self%polcyc1_mat(i,k,l)
                        do m=1,self%wdim
                            logi = [addr_l,self%polcyc2_mat(i,k,m),0]
                            phys = self%comp_addr_phys(logi)
                            self%comps(l,m) = self%get_fcomp(logi,phys)
                        enddo
                    enddo
                    self%pft(i,k) = dot_product(self%polweights_mat(i,k,:), reshape(self%comps,(/self%wlen/)))
                else
                    self%pft(i,k) = CMPLX_ZERO
                endif
            end do
        end do
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
        call self%instrfun_img%kill
    end subroutine kill_polarizer

end module simple_polarizer
