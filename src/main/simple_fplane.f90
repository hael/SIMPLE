! 3D reconstruction from projections using convolution interpolation (gridding)
module simple_fplane
!$ use omp_lib
include 'simple_lib.f08'
use simple_image,         only: image
use simple_parameters,    only: params_glob
use simple_euclid_sigma2, only: eucl_sigma2_glob
use simple_ctf,           only: ctf
implicit none

public :: fplane
private
#include "simple_local_flags.inc"

type :: fplane
    private
    complex, allocatable, public :: cmplx_plane(:,:)             !< On output image pre-multiplied by CTF
    real,    allocatable, public :: ctfsq_plane(:,:)             !< On output CTF normalization
    real,    allocatable         :: ctf_ang(:,:)                 !< CTF effective defocus
    integer,              public :: ldim(3)       = 0            !< dimensions of original image
    integer,              public :: frlims_crop(3,2) = 0         !< Redundant Fourier cropped limits
    integer,              public :: ldim_crop(3)  = 0            !< dimensions of cropped image
    real,                 public :: shconst(3)    = 0.           !< memoized constants for origin shifting
    integer,              public :: nyq_crop      = 0            !< cropped Nyqvist Fourier index
    logical                      :: exists        = .false.      !< Volta phaseplate images or not
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! GETTERS
    procedure :: does_exist
    ! SETTERS
    procedure :: gen_planes
    ! DESTRUCTOR
    procedure :: kill
end type fplane

contains

    ! CONSTRUCTORS

    subroutine new( self, img )
        use simple_ftiter, only: ftiter
        class(fplane),     intent(inout) :: self
        class(image),      intent(in)    :: img
        type(ftiter) :: fiterator
        integer      :: h, k
        call self%kill
        ! Original image dimension
        self%ldim    = img%get_ldim()
        ! shift is with respect to the original image dimensions
        self%shconst = img%get_shconst()
        ! cropped Fourier limits & dimensions
        self%ldim_crop    = [params_glob%box_crop, params_glob%box_crop, 1]
        call fiterator%new(self%ldim_crop, params_glob%smpd_crop)
        self%frlims_crop  = fiterator%loop_lims(3)
        self%nyq_crop     = fiterator%get_lfny(1)
        ! allocations
        allocate(self%cmplx_plane(self%frlims_crop(1,1):self%frlims_crop(1,2),self%frlims_crop(2,1):self%frlims_crop(2,2)),&
                &self%ctfsq_plane(self%frlims_crop(1,1):self%frlims_crop(1,2),self%frlims_crop(2,1):self%frlims_crop(2,2)),&
                self%ctf_ang(self%frlims_crop(1,1):self%frlims_crop(1,2), self%frlims_crop(2,1):self%frlims_crop(2,2)))
        self%cmplx_plane = cmplx(0.,0.)
        self%ctfsq_plane = 0.
        ! CTF pre-calculations
        do k=self%frlims_crop(2,1),self%frlims_crop(2,2)
            do h=self%frlims_crop(1,1),self%frlims_crop(1,2)
                self%ctf_ang(h,k) = atan2(real(k), real(h))
            enddo
        enddo
        self%exists = .true.
    end subroutine new

    logical pure function does_exist( self )
        class(fplane), intent(in) :: self
        does_exist = self%exists
    end function does_exist

    !> Produces CTF multiplied fourier & CTF-squared planes
    subroutine gen_planes( self, img, ctfvars, iptcl)
        class(fplane),    intent(inout) :: self
        class(image),     intent(in)    :: img
        class(ctfparams), intent(in)    :: ctfvars
        integer, optional,intent(in)    :: iptcl
        type(ctf) :: tfun
        complex   :: c
        real      :: invldim(2),inv(2),tval,tvalsq,sqSpatFreq,add_phshift
        integer   :: sigma2_kfromto(2), h,k,sh
        logical   :: use_sigmas
        if( ctfvars%ctfflag /= CTFFLAG_NO )then
            tfun = ctf(ctfvars%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
            call tfun%init(ctfvars%dfx, ctfvars%dfy, ctfvars%angast)
            invldim = 1./real(self%ldim(1:2))
            add_phshift = 0.0
            if( ctfvars%l_phaseplate ) add_phshift = ctfvars%phshift
        endif
        use_sigmas = present(iptcl) .and. params_glob%l_ml_reg
        if( use_sigmas )then
            sigma2_kfromto(1) = lbound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_kfromto(2) = ubound(eucl_sigma2_glob%sigma2_noise,1)
        end if
        !$omp parallel do collapse(2) default(shared) schedule(static) proc_bind(close)&
        !$omp private(h,k,sh,c,tval,tvalsq,inv,sqSpatFreq)
        do h = self%frlims_crop(1,1),self%frlims_crop(1,2)
            do k = self%frlims_crop(2,1),self%frlims_crop(2,2)
                sh = nint(sqrt(real(h*h + k*k)))
                if( sh > self%nyq_crop )then
                    c = cmplx(0.,0.)
                    tvalsq = 0.
                else
                    ! CTF
                    if( ctfvars%ctfflag /= CTFFLAG_NO )then
                        inv        = real([h,k]) * invldim
                        sqSpatFreq = dot_product(inv,inv)
                        tval       = tfun%eval(sqSpatFreq, self%ctf_ang(h,k), add_phshift, .not.params_glob%l_wiener_part)
                        tvalsq     = tval * tval
                    else
                        tval = 1.0
                        tvalsq = tval
                    endif
                    if( ctfvars%ctfflag == CTFFLAG_FLIP ) tval = abs(tval)
                    ! CTF pre-multiplied Fourier component
                    c = tval * img%get_fcomp2D(h,k)
                    ! sigma2 weighing
                    if( use_sigmas) then
                        if(sh >= sigma2_kfromto(1))then
                            c      = c      / eucl_sigma2_glob%sigma2_noise(sh,iptcl)
                            tvalsq = tvalsq / eucl_sigma2_glob%sigma2_noise(sh,iptcl)
                        else
                            c      = c      / eucl_sigma2_glob%sigma2_noise(sigma2_kfromto(1),iptcl)
                            tvalsq = tvalsq / eucl_sigma2_glob%sigma2_noise(sigma2_kfromto(1),iptcl)
                        endif
                    endif
                endif
                self%cmplx_plane(h,k) = c
                self%ctfsq_plane(h,k) = tvalsq
            enddo
        enddo
        !$omp end parallel do
    end subroutine gen_planes

    !>  \brief  is a destructor
    subroutine kill( self )
        class(fplane), intent(inout) :: self !< this instance
        if(allocated(self%cmplx_plane)    ) deallocate(self%cmplx_plane)
        if(allocated(self%ctfsq_plane)    ) deallocate(self%ctfsq_plane)
        if(allocated(self%ctf_ang)        ) deallocate(self%ctf_ang)
        self%exists = .false.
    end subroutine kill

end module simple_fplane
