! polar 2D Fourier transform generation by convolution interpolation (gridding)
module simple_polarizer
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image
use simple_parameters, only: params_glob
implicit none

public :: polarizer
private
#include "simple_local_flags.inc"

complex, parameter :: CMPLX_ZERO = cmplx(0.,0.)

type, extends(image) :: polarizer
    private
    type(image)           :: instrfun_img          !< weights to divide with in real-space
    real,    allocatable  :: polweights_mat(:,:,:) !< polar weights matrix for the image to polar transformer
    real,    allocatable  :: gridweights(:,:)      !< weights to divide with in real-space
    integer, allocatable  :: polcyc1_mat(:,:,:)    !< image cyclic adresses for the image to polar transformer
    integer, allocatable  :: polcyc2_mat(:,:,:)    !< image cyclic adresses for the image to polar transformer
    integer               :: wdim = 0              !< dimension of K-B window
    integer               :: wlen = 0              !< dimension squared of K-B window
    integer               :: pdim(3) = 0           !< Polar-FT matrix dimensions
    logical               :: polarizer_exists = .false.
  contains
    procedure          :: init_polarizer
    procedure, private :: copy_polarizer ! is private as it should not be required anymore
    procedure          :: div_by_instrfun
    procedure, private :: polarize_1, polarize_2
    generic            :: polarize => polarize_1, polarize_2
    procedure          :: polarizer_initialized
    procedure          :: kill_polarizer
end type polarizer

contains

    ! IMAGE TO POLAR FT TRANSFORMER

    !> \brief  initialises the image polarizer
    subroutine init_polarizer( self, pftcc, alpha )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_gridding,         only: gen_instrfun_img
        class(polarizer),        intent(inout) :: self   !< projector instance
        class(polarft_corrcalc), intent(inout) :: pftcc  !< polarft_corrcalc object to be filled
        real,                    intent(in)    :: alpha  !< oversampling factor
        type(kbinterpol)  :: kbwin                 !< window function object
        real, allocatable :: w(:,:)
        real              :: loc(2), d1, d2
        integer           :: win(2,2), lims(2,3), i, k, l, cnt, f1, f2
        logical           :: normalize_weights
        if( .not. pftcc%exists() ) THROW_HARD('polarft_corrcalc object needs to be created; init_polarizer')
        call self%kill_polarizer
        self%pdim = pftcc%get_pdim()
        lims      = transpose(self%loop_lims(3)) ! fortran layered memory
        select case(trim(params_glob%interpfun))
        case('kb')
            kbwin      = kbinterpol(KBWINSZ, alpha)
            self%wdim  = kbwin%get_wdim()
        case('linear')
            self%wdim  = 2
        case DEFAULT
            THROW_HARD('Unsupported interpolation function: '//trim(params_glob%interpfun))
        end select
        self%wlen = self%wdim**2
        allocate( self%polcyc1_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wdim),&
                  &self%polcyc2_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wdim),&
                  &self%polweights_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wlen),&
                  &w(1:self%wdim,1:self%wdim))
        ! instrument function
        normalize_weights = trim(params_glob%interpfun) == 'kb'
        if( params_glob%gridding.eq.'yes' )then
            normalize_weights = .false.
            call self%instrfun_img%new(self%get_ldim(), self%get_smpd())
            call gen_instrfun_img(self%instrfun_img, params_glob%interpfun, kbwin=kbwin)
        endif
        ! cartesian to polar
        if( trim(params_glob%interpfun) == 'kb' )then
            ! Kaiser-Bessel
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
                        self%polcyc1_mat(i, k, cnt) = cyci_1d(lims(:,1), win(1,1)+l-1)
                        self%polcyc2_mat(i, k, cnt) = cyci_1d(lims(:,2), win(2,1)+l-1)
                    end do
                    self%polweights_mat(i,k,:) = reshape(w,(/self%wlen/))
                    if( normalize_weights ) self%polweights_mat(i,k,:) = self%polweights_mat(i,k,:) / sum(w)
                enddo
            enddo
            !$omp end parallel do
        else
            ! Bi-linear
            !$omp parallel do collapse(2) schedule(static) private(i,k,w,loc,f1,f2,d1,d2) default(shared) proc_bind(close)
            do i=1,self%pdim(1)
                do k=self%pdim(2),self%pdim(3)
                    loc = pftcc%get_coord(i,k)
                    f1  = int(floor(loc(1)))
                    d1  = loc(1) - f1
                    f2  = int(floor(loc(2)))
                    d2  = loc(2) - f2
                    self%polcyc1_mat(i, k,  1)  = cyci_1d(lims(:,1), f1)
                    self%polcyc1_mat(i, k,  2)  = cyci_1d(lims(:,1), f1+1)
                    self%polcyc2_mat(i, k,  1)  = cyci_1d(lims(:,2), f2)
                    self%polcyc2_mat(i, k,  2)  = cyci_1d(lims(:,2), f2+1)
                    w      = 1.
                    w(1,:) = w(1,:) * (1.0-d1)
                    w(:,1) = w(:,1) * (1.0-d2)
                    w(2,:) = w(2,:) * d1
                    w(:,2) = w(:,2) * d2
                    self%polweights_mat(i,k,:) = reshape(w,(/self%wlen/))
                enddo
            enddo
            !$omp end parallel do
        endif
        deallocate(w)
        self%polarizer_exists = .true.
    end subroutine init_polarizer

    subroutine copy_polarizer(self, self_in)
        class(polarizer), intent(inout) :: self    !< projector instance
        class(polarizer), intent(in)    :: self_in
        if( .not.self_in%polarizer_exists ) THROW_HARD('Polarizer instance has not bee initialized!')
        call self%kill_polarizer
        self%pdim = self_in%pdim
        self%wdim = self_in%wdim
        self%wlen = self_in%wlen
        allocate( self%polcyc1_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wdim),    source=self_in%polcyc1_mat )
        allocate( self%polcyc2_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wdim),    source=self_in%polcyc2_mat )
        allocate( self%polweights_mat(1:self%pdim(1), self%pdim(2):self%pdim(3), 1:self%wlen), source=self_in%polweights_mat )
        if( self_in%instrfun_img%exists() ) call self%instrfun_img%copy(self_in%instrfun_img)
    end subroutine copy_polarizer

    !> \brief  divide by gridding weights in real-space prior to FFT &
    !>         change to polar coordinate system, keep serial
    subroutine div_by_instrfun( self, img )
        class(polarizer), intent(inout) :: self
        class(image),     intent(inout) :: img
        if( img%is_ft() ) THROW_HARD('Image must be in real-space in div_by_instrfun')
        if( .not.(self.eqdims.self%instrfun_img) ) THROW_HARD('Incompatible dimensions in div_by_instrfun 1')
        if( .not.(self.eqdims.img) ) THROW_HARD('Incompatible dimensions in div_by_instrfun 2')
        call img%div(self%instrfun_img)
    end subroutine div_by_instrfun

    !> \brief  creates the polar Fourier transform from self
    !!         KEEP THIS ROUTINE SERIAL
    subroutine polarize_1( self, pftcc, img_ind, isptcl, iseven, mask )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(polarizer),        intent(in)    :: self    !< projector instance
        class(polarft_corrcalc), intent(inout) :: pftcc   !< polarft_corrcalc object to be filled
        integer,                 intent(in)    :: img_ind !< image index
        logical,                 intent(in)    :: isptcl  !< is ptcl (or reference)
        logical,                 intent(in)    :: iseven  !< is even (or odd)
        logical, optional,       intent(in)    :: mask(:) !< interpolation mask, all .false. set to CMPLX_ZERO
        call self%polarize_2(pftcc, self, img_ind, isptcl, iseven, mask )
    end subroutine polarize_1

    !> \brief  creates the polar Fourier transform from a given image
    !!         KEEP THIS ROUTINE SERIAL
    subroutine polarize_2( self, pftcc, img, img_ind, isptcl, iseven, mask )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(polarizer),        intent(in)    :: self    !< projector instance
        class(polarft_corrcalc), intent(inout) :: pftcc   !< polarft_corrcalc object to be filled
        class(image),            intent(in)    :: img
        integer,                 intent(in)    :: img_ind !< image index
        logical,                 intent(in)    :: isptcl  !< is ptcl (or reference)
        logical,                 intent(in)    :: iseven  !< is even (or odd)
        logical, optional,       intent(in)    :: mask(:) !< interpolation mask, all .false. set to CMPLX_ZERO
        complex, pointer :: pft(:,:)
        complex :: fcomps(1:self%wdim,1:self%wdim)
        integer :: i, k, l, m, addr_l
        ! get temporary pft matrix
        call pftcc%get_work_pft_ptr(pft)
        ! interpolate
        do i=1,self%pdim(1)
            do k=self%pdim(2),self%pdim(3)
                do l=1,self%wdim
                    addr_l = self%polcyc1_mat(i,k,l)
                    do m=1,self%wdim
                        fcomps(l,m) = img%get_fcomp2D(addr_l,self%polcyc2_mat(i,k,m))
                    enddo
                enddo
                pft(i,k) = dot_product(self%polweights_mat(i,k,:), reshape(fcomps,(/self%wlen/)))
            end do
        end do
        if( present(mask) )then
            ! band masking
            do k=self%pdim(2),self%pdim(3)
                if( .not.mask(k) ) pft(:,k) = CMPLX_ZERO
            enddo
        endif
        if( isptcl )then
            call pftcc%set_ptcl_pft(img_ind, pft)
        else
            call pftcc%set_ref_pft(img_ind, pft, iseven)
        endif
        nullify(pft)
    end subroutine polarize_2

    logical function polarizer_initialized( self )
        class(polarizer), intent(in) :: self
        polarizer_initialized = self%polarizer_exists
    end function polarizer_initialized

    ! DESTRUCTOR

    !>  \brief  is a destructor of impolarizer
    subroutine kill_polarizer( self )
        class(polarizer), intent(inout) :: self !< projector instance
        if( allocated(self%polweights_mat) ) deallocate(self%polweights_mat)
        if( allocated(self%polcyc1_mat)    ) deallocate(self%polcyc1_mat)
        if( allocated(self%polcyc2_mat)    ) deallocate(self%polcyc2_mat)
        call self%instrfun_img%kill
        self%polarizer_exists = .false.
    end subroutine kill_polarizer

end module simple_polarizer
