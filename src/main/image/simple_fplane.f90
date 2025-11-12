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
    type(kbinterpol),     public :: kbwin                        !< window function object
    complex, allocatable, public :: cmplx_plane(:,:)             !< On output image pre-multiplied by CTF
    real,    allocatable, public :: ctfsq_plane(:,:)             !< On output CTF normalization
    real,    allocatable         :: ctf_ang(:,:)                 !< CTF effective defocus
    real(dp),allocatable         :: hcos(:), hsin(:)             !< For fast shift
    real,    allocatable         :: sigma2_noise(:)              !< Noise power
    integer,              public :: ldim(3)       = 0            !< dimensions of original image
    integer,              public :: frlims_crop(3,2)   = 0       !< Redundant Fourier cropped limits
    integer,              public :: frlims_croppd(3,2) = 0       !< Redundant Fourier cropped padded limits
    integer,              public :: ldim_crop(3)  = 0            !< dimensions of cropped image
    real,                 public :: shconst(3)    = 0.           !< memoized constants for origin shifting
    integer,              public :: nyq_crop      = 0            !< cropped Nyqvist Fourier index
    integer,              public :: nyq_croppd    = 0            !< cropped Nyqvist Fourier index
    real                         :: winsz         = RECWINSZ     !< window half-width
    real                         :: alpha         = KBALPHA      !< oversampling ratio
    integer                      :: wdim          = 0            !< dim of interpolation matrix
    logical                      :: genplane      = .true.
    logical,              public :: padded        = .false.      !< Whether the resulting planes are padded
    logical                      :: exists        = .false.      !< Volta phaseplate images or not
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! GETTERS
    procedure :: does_exist
    procedure :: convert2img
    ! SETTERS
    procedure :: gen_planes
    procedure :: gen_planes_pad
    ! MODIFIERS
    procedure :: zero
    procedure :: neg
    ! DESTRUCTOR
    procedure :: kill
end type fplane

contains

    ! CONSTRUCTORS

    subroutine new( self, img, pad, genplane )
        use simple_ftiter, only: ftiter
        class(fplane),     intent(inout) :: self
        class(image),      intent(in)    :: img
        logical, optional, intent(in)    :: pad, genplane
        type(ftiter) :: fiterator
        integer      :: h, k
        call self%kill
        self%padded = .false.
        if( present(pad) ) self%padded = pad
        self%genplane = .true.
        if( present(genplane) ) self%genplane = genplane
        ! Original image dimension
        self%ldim    = img%get_ldim()
        ! shift is with respect to the original image dimensions
        self%shconst = img%get_shconst()
        ! cropped Fourier limits & dimensions
        self%ldim_crop    = [params_glob%box_crop, params_glob%box_crop, 1]
        call fiterator%new(self%ldim_crop, params_glob%smpd_crop)
        self%frlims_crop = fiterator%loop_lims(3)
        self%nyq_crop    = fiterator%get_lfny(1)
        if( self%padded )then
            ! the resulting Fourier slices will be padded
            call fiterator%new([params_glob%box_croppd, params_glob%box_croppd, 1], params_glob%smpd_crop)
            self%frlims_croppd = fiterator%loop_lims(3)
            self%nyq_croppd    = fiterator%get_lfny(1)
            ! intepolation window
            self%winsz = params_glob%winsz
            self%alpha = params_glob%alpha
            self%kbwin = kbinterpol(self%winsz,self%alpha)
            self%wdim  = self%kbwin%get_wdim()
            ! allocations
            allocate(self%cmplx_plane(self%frlims_croppd(1,1):self%frlims_croppd(1,2),self%frlims_croppd(2,1):self%frlims_croppd(2,2)),&
            &self%ctfsq_plane(self%frlims_croppd(1,1):self%frlims_croppd(1,2),self%frlims_croppd(2,1):self%frlims_croppd(2,2)))
        else
            self%frlims_croppd = 0
            self%nyq_croppd    = 0
            allocate(self%cmplx_plane(self%frlims_crop(1,1):self%frlims_crop(1,2),self%frlims_crop(2,1):self%frlims_crop(2,2)),&
            &self%ctfsq_plane(self%frlims_crop(1,1):self%frlims_crop(1,2),self%frlims_crop(2,1):self%frlims_crop(2,2)))
        endif
        call self%zero
        if( self%genplane )then
            ! the object will be used to prep image for reconstruction
            ! otherwise the following allocations are not done to reduce memory usage
            allocate(self%ctf_ang(self%frlims_crop(1,1):self%frlims_crop(1,2), self%frlims_crop(2,1):self%frlims_crop(2,2)),&
            &self%hcos(self%frlims_crop(1,1):self%frlims_crop(1,2)), self%hsin(self%frlims_crop(1,1):self%frlims_crop(1,2)),&
            &self%sigma2_noise(1:self%nyq_crop))
            ! CTF pre-calculations
            do k=self%frlims_crop(2,1),self%frlims_crop(2,2)
                do h=self%frlims_crop(1,1),self%frlims_crop(1,2)
                    self%ctf_ang(h,k) = atan2(real(k), real(h))
                enddo
            enddo
            self%sigma2_noise = 0.
        endif
        self%exists = .true.
    end subroutine new

    logical pure function does_exist( self )
        class(fplane), intent(in) :: self
        does_exist = self%exists
    end function does_exist

    !> Produces shifted, CTF multiplied fourier & CTF-squared planes
    subroutine gen_planes( self, img, ctfvars, shift, iptcl )
        class(fplane),                  intent(inout) :: self
        class(image),                   intent(in)    :: img
        class(ctfparams),               intent(in)    :: ctfvars
        real,                           intent(in)    :: shift(2)
        integer,                        intent(in)    :: iptcl
        type(ctf) :: tfun
        complex   :: c, s
        real(dp)  :: pshift(2), arg, kcos,ksin
        real      :: invldim(2),inv(2),tval,tvalsq,sqSpatFreq,add_phshift
        integer   :: sigma2_kfromto(2), h,k,sh
        if( ctfvars%ctfflag /= CTFFLAG_NO )then
            tfun = ctf(ctfvars%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
            call tfun%init(ctfvars%dfx, ctfvars%dfy, ctfvars%angast)
            invldim = 1./real(self%ldim(1:2))
            add_phshift = 0.0
            if( ctfvars%l_phaseplate ) add_phshift = ctfvars%phshift
        endif
        pshift = real(-shift * self%shconst(1:2),dp)
        do h = self%frlims_crop(1,1),self%frlims_crop(1,2)
            arg = real(h,dp)*pshift(1)
            self%hcos(h) = dcos(arg)
            self%hsin(h) = dsin(arg)
        enddo
        if( params_glob%l_ml_reg )then
            sigma2_kfromto(1) = lbound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_kfromto(2) = ubound(eucl_sigma2_glob%sigma2_noise,1)
            self%sigma2_noise(sigma2_kfromto(1):self%nyq_crop) = eucl_sigma2_glob%sigma2_noise(sigma2_kfromto(1):self%nyq_crop,iptcl)
        end if
        ! prep slice for k in [N/2;0]
        do k = self%frlims_crop(2,1),0
            arg  = real(k,dp)*pshift(2)
            kcos = dcos(arg)
            ksin = dsin(arg)
            do h = self%frlims_crop(1,1),self%frlims_crop(1,2)
                sh = nint(sqrt(real(h*h + k*k)))
                if( sh > self%nyq_crop )then
                    c      = cmplx(0.,0.)
                    tvalsq = 0.
                else
                    ! Shift
                    s = cmplx(kcos*self%hcos(h)-ksin*self%hsin(h), kcos*self%hsin(h)+ksin*self%hcos(h),sp)
                    c = img%get_fcomp2D(h,k) * s
                    ! CTF
                    if( ctfvars%ctfflag /= CTFFLAG_NO )then
                        inv        = real([h,k]) * invldim
                        sqSpatFreq = dot_product(inv,inv)
                        tval       = tfun%eval(sqSpatFreq, self%ctf_ang(h,k), add_phshift, .not.params_glob%l_wiener_part)
                        tvalsq     = tval * tval
                        if( ctfvars%ctfflag == CTFFLAG_FLIP ) tval = abs(tval)
                        c          = tval * c
                    else
                        tvalsq = 1.0
                    endif
                    ! sigma2 weighing
                    if( params_glob%l_ml_reg ) then
                        if(sh >= sigma2_kfromto(1))then
                            c      = c      / self%sigma2_noise(sh)
                            tvalsq = tvalsq / self%sigma2_noise(sh)
                        else
                            c      = c      / self%sigma2_noise(sigma2_kfromto(1))
                            tvalsq = tvalsq / self%sigma2_noise(sigma2_kfromto(1))
                        endif
                    endif
                endif
                self%cmplx_plane(h,k) = c
                self%ctfsq_plane(h,k) = tvalsq
            enddo
        enddo
        ! prep slice for k in [1;N/2-1] with Friedel symmetry
        do k = 1,self%frlims_crop(2,2)
            do h = self%frlims_crop(1,1),self%frlims_crop(1,2)
                sh = nint(sqrt(real(h*h + k*k)))
                if( sh > self%nyq_crop )then
                    self%cmplx_plane(h,k) = cmplx(0.,0.)
                    self%ctfsq_plane(h,k) = 0.
                else
                    self%cmplx_plane(h,k) = conjg(self%cmplx_plane(-h,-k))
                    self%ctfsq_plane(h,k) = self%ctfsq_plane(-h,-k)
                endif
            enddo
        enddo
    end subroutine gen_planes

    !> Produces padded shifted, rotated, CTF multiplied fourier & CTF-squared planes
    subroutine gen_planes_pad( self, img, ctfvars, shift, e3, iptcl, linear )
        class(fplane),                  intent(inout) :: self
        class(image),                   intent(inout) :: img
        class(ctfparams),               intent(in)    :: ctfvars
        real,                           intent(in)    :: shift(2), e3
        integer,                        intent(in)    :: iptcl
        logical,                        intent(in)    :: linear
        type(ctf) :: tfun
        complex   :: c, s
        real(dp)  :: pshift(2), arg, kcos,ksin
        real      :: kbw(self%wdim,self%wdim), w(2,2),rmat(2,2),loc(2),inv(2),d(2)
        real      :: tval,tvalsq,sqSpatFreq,add_phshift, dh,dk
        integer   :: win(2,2),sigma2_kfromto(2), i,j,h,k,hh,kk,sh, iwinsz, fh,fk, hlim
        if( .not.self%padded ) THROW_HARD('gen_planes_pad only for use with padding!')
        ! CTF
        if( ctfvars%ctfflag /= CTFFLAG_NO )then
            tfun = ctf(ctfvars%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
            call tfun%init(ctfvars%dfx, ctfvars%dfy, ctfvars%angast)
            add_phshift = 0.0
            if( ctfvars%l_phaseplate ) add_phshift = ctfvars%phshift
        endif
        ! Shift precomputation
        pshift = real(-shift * self%shconst(1:2),dp)
        do h = self%frlims_crop(1,1),self%frlims_crop(1,2)
            arg = real(h,dp)*pshift(1)
            self%hcos(h) = dcos(arg)
            self%hsin(h) = dsin(arg)
        enddo
        ! rotation & scale
        call rotmat2D(e3, rmat)
        rmat = self%alpha * rmat
        ! sigma2
        if( params_glob%l_ml_reg )then
            sigma2_kfromto(1) = lbound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_kfromto(2) = ubound(eucl_sigma2_glob%sigma2_noise,1)
            self%sigma2_noise(sigma2_kfromto(1):self%nyq_crop) = eucl_sigma2_glob%sigma2_noise(sigma2_kfromto(1):self%nyq_crop,iptcl)
        end if
        ! Interpolation
        self%cmplx_plane = cmplx(0.,0.)
        self%ctfsq_plane = 0.
        ! Cycling within non-redundant components:
        ! k in [-N/2,-1], h in [-N/2;N/2]
        ! k = 0, h in [-N/2;-1]
        ! h=0, k=0 is treated after the loop
        if( linear )then
            ! Bilinear interpolation
            do k = self%frlims_crop(2,1),0
                arg  = real(k,dp)*pshift(2)
                kcos = dcos(arg)
                ksin = dsin(arg)
                if( k == 0 )then
                    hlim = -1
                else
                    hlim = self%frlims_crop(1,2)
                endif
                do h = self%frlims_crop(1,1),hlim
                    sh = nint(sqrt(real(h*h + k*k)))
                    if( sh > self%nyq_crop ) cycle
                    ! Shift
                    s = cmplx(kcos*self%hcos(h)-ksin*self%hsin(h), kcos*self%hsin(h)+ksin*self%hcos(h),sp)
                    c = img%get_fcomp2D(h,k) * s
                    ! CTF
                    if( ctfvars%ctfflag /= CTFFLAG_NO )then
                        inv        = real([h,k]) / real(self%ldim(1:2))
                        sqSpatFreq = dot_product(inv,inv)
                        tval       = tfun%eval(sqSpatFreq, self%ctf_ang(h,k), add_phshift, .not.params_glob%l_wiener_part)
                        tvalsq     = tval * tval
                        if( ctfvars%ctfflag == CTFFLAG_FLIP ) tval = abs(tval)
                        c          = tval * c
                    else
                        tvalsq = 1.0
                    endif
                    ! sigma2 weighing
                    if( params_glob%l_ml_reg ) then
                        if(sh >= sigma2_kfromto(1))then
                            c      = c      / self%sigma2_noise(sh)
                            tvalsq = tvalsq / self%sigma2_noise(sh)
                        else
                            c      = c      / self%sigma2_noise(sigma2_kfromto(1))
                            tvalsq = tvalsq / self%sigma2_noise(sigma2_kfromto(1))
                        endif
                    endif
                    ! rotation
                    loc = matmul(real([h,k]), rmat)
                    ! bottom left corner
                    fh  = floor(loc(1))
                    fk  = floor(loc(2))
                    ! kernel
                    dh = loc(1) - real(fh)
                    dk = loc(2) - real(fk)
                    w(1,:) = 1.-dh
                    w(2,:) = dh
                    w(:,1) = w(:,1)*(1.-dk)
                    w(:,2) = w(:,2)*dk
                    ! interpolation
                    i = 0
                    do hh = fh,fh+1
                        i = i+1
                        if( abs(hh) > self%nyq_croppd )cycle
                        j = 0
                        do kk = fk,fk+1
                            j = j+1
                            if( (kk >= self%frlims_croppd(2,1)) .and. (kk <= self%frlims_croppd(2,2)) )then
                                self%cmplx_plane(hh,kk) = self%cmplx_plane(hh,kk) + w(i,j) * c
                                self%ctfsq_plane(hh,kk) = self%ctfsq_plane(hh,kk) + w(i,j) * tvalsq
                            endif
                            ! Friedel symmetric
                            ! hh is in [-N/2;N/2] so no need to check for bounds of -hh again
                            if( (-kk >= self%frlims_croppd(2,1)) .and. (-kk <= self%frlims_croppd(2,2)) )then
                                self%cmplx_plane(-hh,-kk) = self%cmplx_plane(-hh,-kk) + w(i,j) * conjg(c)
                                self%ctfsq_plane(-hh,-kk) = self%ctfsq_plane(-hh,-kk) + w(i,j) * tvalsq
                            endif
                        enddo
                    enddo
                enddo
            enddo
        else
            ! KB interpolation
            iwinsz = ceiling(self%winsz - 0.5)
            do k = self%frlims_crop(2,1),0
                arg  = real(k,dp)*pshift(2)
                kcos = dcos(arg)
                ksin = dsin(arg)
                if( k == 0 )then
                    hlim = -1 ! h=0, k=0 is treated after the loop
                else
                    hlim = self%frlims_crop(1,2)
                endif
                do h = self%frlims_crop(1,1),hlim
                    sh = nint(sqrt(real(h*h + k*k)))
                    if( sh > self%nyq_crop ) cycle
                    ! Shift
                    s = cmplx(kcos*self%hcos(h)-ksin*self%hsin(h), kcos*self%hsin(h)+ksin*self%hcos(h),sp)
                    c = img%get_fcomp2D(h,k) * s
                    ! CTF
                    if( ctfvars%ctfflag /= CTFFLAG_NO )then
                        inv        = real([h,k]) / real(self%ldim(1:2))
                        sqSpatFreq = dot_product(inv,inv)
                        tval       = tfun%eval(sqSpatFreq, self%ctf_ang(h,k), add_phshift, .not.params_glob%l_wiener_part)
                        tvalsq     = tval * tval
                        if( ctfvars%ctfflag == CTFFLAG_FLIP ) tval = abs(tval)
                        c          = tval * c
                    else
                        tvalsq = 1.0
                    endif
                    ! sigma2 weighing
                    if( params_glob%l_ml_reg ) then
                        if(sh >= sigma2_kfromto(1))then
                            c      = c      / self%sigma2_noise(sh)
                            tvalsq = tvalsq / self%sigma2_noise(sh)
                        else
                            c      = c      / self%sigma2_noise(sigma2_kfromto(1))
                            tvalsq = tvalsq / self%sigma2_noise(sigma2_kfromto(1))
                        endif
                    endif
                    ! rotation
                    loc = matmul(real([h,k]), rmat)
                    ! window
                    win(1,:) = nint(loc)
                    win(2,:) = win(1,:) + iwinsz
                    win(1,:) = win(1,:) - iwinsz
                    ! kernel
                    kbw = 1.
                    do i = 1,self%wdim
                        d = real(win(1,:) + i - 1) - loc
                        kbw(i,:) = kbw(i,:) * self%kbwin%apod(d(1))
                        kbw(:,i) = kbw(:,i) * self%kbwin%apod(d(2))
                    enddo
                    kbw = kbw / sum(kbw)
                    ! interpolation
                    i = 0
                    do hh = win(1,1),win(2,1)
                        i = i+1
                        if( abs(hh) > self%nyq_croppd )cycle
                        j = 0
                        do kk = win(1,2),win(2,2)
                            j = j+1
                            if( (kk >= self%frlims_croppd(2,1)) .and. (kk <= self%frlims_croppd(2,2)) )then
                                self%cmplx_plane(hh,kk) = self%cmplx_plane(hh,kk) + kbw(i,j) * c
                                self%ctfsq_plane(hh,kk) = self%ctfsq_plane(hh,kk) + kbw(i,j) * tvalsq
                            endif
                            ! Friedel symmetric
                            ! hh is in [-N/2;N/2] so no need to check for bounds of -hh
                            if( (-kk >= self%frlims_croppd(2,1)) .and. (-kk <= self%frlims_croppd(2,2)) )then
                                self%cmplx_plane(-hh,-kk) = self%cmplx_plane(-hh,-kk) + kbw(i,j) * conjg(c)
                                self%ctfsq_plane(-hh,-kk) = self%ctfsq_plane(-hh,-kk) + kbw(i,j) * tvalsq
                            endif
                        enddo
                    enddo
                enddo
            enddo
        endif
        ! Interpolation free DC
        c = img%get_fcomp2D(0,0)
        if( ctfvars%ctfflag /= CTFFLAG_NO )then
            tval       = tfun%eval(0., 0., add_phshift, .not.params_glob%l_wiener_part)
            tvalsq     = tval * tval
            if( ctfvars%ctfflag == CTFFLAG_FLIP ) tval = abs(tval)
            c          = tval * c
        else
            tvalsq = 1.0
        endif
        if( params_glob%l_ml_reg ) then
            c      = c      / self%sigma2_noise(sigma2_kfromto(1))
            tvalsq = tvalsq / self%sigma2_noise(sigma2_kfromto(1))
        endif
        self%cmplx_plane(0,0) = c
        self%ctfsq_plane(0,0) = tvalsq
    end subroutine gen_planes_pad

    subroutine convert2img( self, fcimg, ctfsqimg )
        class(fplane), intent(in)    :: self
        class(image),  intent(inout) :: fcimg, ctfsqimg
        complex, pointer :: pcmat(:,:,:)
        integer :: c
        if( self%padded )then
            call fcimg%new([params_glob%box_croppd,params_glob%box_croppd,1], params_glob%smpd_crop)
            call ctfsqimg%new([params_glob%box_croppd,params_glob%box_croppd,1], params_glob%smpd_crop)
            c = params_glob%box_croppd/2+1
            call fcimg%set_ft(.true.)
            call fcimg%get_cmat_ptr(pcmat)
            pcmat(1:self%frlims_croppd(1,2),c+self%frlims_croppd(2,1):c+self%frlims_croppd(2,2),1) =&
                &self%cmplx_plane(0:self%frlims_croppd(1,2)-1,:)
            call fcimg%shift_phorig
            call ctfsqimg%set_ft(.true.)
            call ctfsqimg%get_cmat_ptr(pcmat)
            pcmat(1:self%frlims_croppd(1,2),c+self%frlims_croppd(2,1):c+self%frlims_croppd(2,2),1) =&
                &cmplx(self%ctfsq_plane(0:self%frlims_croppd(1,2)-1,:),0.)
        else
            call fcimg%new(self%ldim, params_glob%smpd)
            call ctfsqimg%new(self%ldim, params_glob%smpd)
            c = self%ldim(1)/2+1
            call fcimg%set_ft(.true.)
            call fcimg%get_cmat_ptr(pcmat)
            pcmat(1:self%frlims_crop(1,2),c+self%frlims_crop(2,1):c+self%frlims_crop(2,2),1) = self%cmplx_plane(0:self%frlims_crop(1,2)-1,:)
            call fcimg%shift_phorig
            call ctfsqimg%set_ft(.true.)
            call ctfsqimg%get_cmat_ptr(pcmat)
            pcmat(1:self%frlims_crop(1,2),c+self%frlims_crop(2,1):c+self%frlims_crop(2,2),1) = cmplx(self%ctfsq_plane(0:self%frlims_crop(1,2)-1,:),0.)
        endif
        call ctfsqimg%shift_phorig
        nullify(pcmat)
    end subroutine convert2img

    elemental subroutine zero( self )
        class(fplane), intent(inout) :: self
        self%cmplx_plane  = cmplx(0.,0.)
        self%ctfsq_plane  = 0.
    end subroutine zero

    subroutine neg( self )
        class(fplane),    intent(inout) :: self
        self%cmplx_plane = -self%cmplx_plane
        self%ctfsq_plane = -self%ctfsq_plane
    end subroutine neg

    !>  \brief  is a destructor
    elemental subroutine kill( self )
        class(fplane), intent(inout) :: self !< this instance
        if( allocated(self%cmplx_plane)  ) deallocate(self%cmplx_plane)
        if( allocated(self%ctfsq_plane)  ) deallocate(self%ctfsq_plane)
        if( allocated(self%ctf_ang)      ) deallocate(self%ctf_ang)
        if( allocated(self%hcos)         ) deallocate(self%hcos)
        if( allocated(self%hsin)         ) deallocate(self%hsin)
        if( allocated(self%sigma2_noise) ) deallocate(self%sigma2_noise)
        self%padded   = .false.
        self%genplane = .true.
        self%exists   = .false.
    end subroutine kill

end module simple_fplane
