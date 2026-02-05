!@descr: re-arranged Fourier transformed image for fast 3D gridding
module simple_fplane
use simple_core_module_api
use simple_memoize_ft_maps
use simple_image,         only: image
use simple_parameters,    only: params_glob
use simple_euclid_sigma2, only: eucl_sigma2_glob
use simple_ctf,           only: ctf
implicit none

public :: fplane
private
#include "simple_local_flags.inc"

type :: fplane
    complex, allocatable :: cmplx_plane(:,:)      !< On output image pre-multiplied by CTF
    real,    allocatable :: ctfsq_plane(:,:)      !< On output CTF normalization
    integer              :: frlims_crop(3,2) = 0  !< Redundant Fourier cropped limits
    real                 :: shconst(3)       = 0. !< memoized constants for origin shifting
    integer              :: nyq_crop         = 0  !< cropped Nyqvist Fourier index
    logical              :: exists           = .false.
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! GETTERS
    procedure :: does_exist
    ! SETTERS
    procedure :: gen_planes
    ! MODIFIERS
    procedure :: zero
    ! DESTRUCTOR
    procedure :: kill
end type fplane

contains

    ! CONSTRUCTORS

    subroutine new( self, img )
        use simple_ftiter, only: ftiter
        class(fplane), intent(inout) :: self
        class(image),  intent(in)    :: img
        type(ftiter) :: fiterator
        call self%kill
        ! shift is with respect to the original image dimensions
        self%shconst = img%get_shconst()
        ! cropped Fourier limits & dimensions
        call fiterator%new([params_glob%box_crop, params_glob%box_crop, 1], params_glob%smpd_crop)
        self%frlims_crop = fiterator%loop_lims(3)
        self%nyq_crop    = fiterator%get_lfny(1)
        allocate(self%cmplx_plane(self%frlims_crop(1,1):self%frlims_crop(1,2),self%frlims_crop(2,1):self%frlims_crop(2,2)),&
        &self%ctfsq_plane(self%frlims_crop(1,1):self%frlims_crop(1,2),self%frlims_crop(2,1):self%frlims_crop(2,2)))
        call self%zero
        self%exists = .true.
    end subroutine new

    logical pure function does_exist( self )
        class(fplane), intent(in) :: self
        does_exist = self%exists
    end function does_exist

    ! example usage of upsample_sigma2
    ! real, allocatable :: sigma2_noise(:)
    ! integer           :: kfromto(2), nyq_croppd
    ! nyq_croppd = padded_img%get_lfny(1)
    ! allocate(sigma2_noise(0:nyq_croppd), source=0.0)
    ! kfromto = eucl_sigma2_glob%get_kfromto()
    ! call upsample_sigma2(kfromto(1), self%nyq_crop,&
    !   &eucl_sigma2_glob%sigma2_noise(kfromto(1):self%nyq_crop, iptcl), nyq_croppd, sigma2_noise)

    !> Produces shifted, CTF multiplied fourier & CTF-squared planes
    subroutine gen_planes( self, img, ctfvars, shift, iptcl )
        class(fplane),                  intent(inout) :: self
        class(image),                   intent(in)    :: img
        class(ctfparams),               intent(in)    :: ctfvars
        real,                           intent(in)    :: shift(2)
        integer,                        intent(in)    :: iptcl
        type(ctf)                :: tfun
        real,    allocatable     :: sigma2_noise(:) !< Noise power spectrum for ML regularization
        complex(c_float_complex) :: c, w1, w2, ph0, ph_h, ph_k
        real(dp) :: pshift(2)
        real     :: tval, tvalsq, add_phshift
        integer  :: physh, physk, sigma2_kfromto(2), h,k,shell
        logical  :: l_ctf, l_flip
        ! CTF
        l_ctf = ctfvars%ctfflag /= CTFFLAG_NO
        if( l_ctf )then
            l_flip = ctfvars%ctfflag == CTFFLAG_FLIP
            tfun   = ctf(ctfvars%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
            call tfun%init(ctfvars%dfx, ctfvars%dfy, ctfvars%angast)
            add_phshift = merge(ctfvars%phshift, 0.0, ctfvars%l_phaseplate)
        endif
        if( params_glob%l_ml_reg )then
            ! the object will be used to prep image for reconstruction
            ! otherwise the following allocations are not done to reduce memory usage
            allocate(sigma2_noise(1:self%nyq_crop), source=0.0)
            sigma2_kfromto(1) = lbound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_kfromto(2) = ubound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_noise(sigma2_kfromto(1):self%nyq_crop) = eucl_sigma2_glob%sigma2_noise(sigma2_kfromto(1):self%nyq_crop,iptcl)
        end if
        ! Shift precomputation
        pshift = real(-shift * self%shconst(1:2),dp)
        ! phase increments exp(i*sh1), exp(i*sh2)
        w1 = cmplx( real(cos(pshift(1)), c_float), real(sin(pshift(1)), c_float), kind=c_float_complex )
        w2 = cmplx( real(cos(pshift(2)), c_float), real(sin(pshift(2)), c_float), kind=c_float_complex )
        ! starting phase for h: exp(i*hmin*sh1) & k: exp(i*kmin*sh2)
        ph0  = cmplx( real(cos(real(self%frlims_crop(1,1),dp)*pshift(1)), c_float), &
                    real(sin(real(self%frlims_crop(1,1),dp)*pshift(1)), c_float), kind=c_float_complex )
        ph_k = cmplx( real(cos(real(self%frlims_crop(2,1),dp)*pshift(2)), c_float), &
                    real(sin(real(self%frlims_crop(2,1),dp)*pshift(2)), c_float), kind=c_float_complex )
        ! prep slice for k in [N/2;0]
        do k = self%frlims_crop(2,1),0
            ph_h = ph0
            do h = self%frlims_crop(1,1),self%frlims_crop(1,2)
                shell = nint(sqrt(real(h*h + k*k)))
                if( shell > self%nyq_crop )then
                    c      = cmplx(0.,0.)
                    tvalsq = 0.
                else
                    ! Retrieve component & shift
                    physh = ft_map_phys_addrh(h,k)
                    physk = ft_map_phys_addrk(h,k)
                    c    = merge(conjg(img%get_cmat_at(physh,physk,1)), img%get_cmat_at(physh,physk,1), h<0) * (ph_k * ph_h)
                    ! CTF
                    if( l_ctf )then
                        tval   = tfun%eval(ft_map_spafreqsq(h,k), ft_map_astigang(h,k), add_phshift)
                        tvalsq = tval * tval
                        c      = merge(abs(tval)*c, tval*c, l_flip)
                    else
                        tvalsq = 1.0
                    endif
                    ! sigma2 weighing
                    if( params_glob%l_ml_reg ) then
                        if(shell >= sigma2_kfromto(1))then
                            c      = c      / sigma2_noise(shell)
                            tvalsq = tvalsq / sigma2_noise(shell)
                        else
                            c      = c      / sigma2_noise(sigma2_kfromto(1))
                            tvalsq = tvalsq / sigma2_noise(sigma2_kfromto(1))
                        endif
                    endif
                endif
                self%cmplx_plane(h,k) = c
                self%ctfsq_plane(h,k) = tvalsq
                ph_h = ph_h * w1    ! phase shift along h
            enddo
            ph_k = ph_k * w2        ! phase shift along k
        enddo
        ! prep slice for k in [1;N/2-1] with Friedel symmetry
        do k = 1,self%frlims_crop(2,2)
            do h = self%frlims_crop(1,1),self%frlims_crop(1,2)
                shell = nint(sqrt(real(h*h + k*k)))
                if( shell > self%nyq_crop )then
                    self%cmplx_plane(h,k) = cmplx(0.,0.)
                    self%ctfsq_plane(h,k) = 0.
                else
                    self%cmplx_plane(h,k) = conjg(self%cmplx_plane(-h,-k))
                    self%ctfsq_plane(h,k) = self%ctfsq_plane(-h,-k)
                endif
            enddo
        enddo
        if( allocated(sigma2_noise) ) deallocate(sigma2_noise)
    end subroutine gen_planes

    elemental subroutine zero( self )
        class(fplane), intent(inout) :: self
        self%cmplx_plane  = cmplx(0.,0.)
        self%ctfsq_plane  = 0.
    end subroutine zero

    !>  \brief  is a destructor
    elemental subroutine kill( self )
        class(fplane), intent(inout) :: self !< this instance
        if( allocated(self%cmplx_plane)  ) deallocate(self%cmplx_plane)
        if( allocated(self%ctfsq_plane)  ) deallocate(self%ctfsq_plane)
        self%exists   = .false.
    end subroutine kill

end module simple_fplane
