!@descr: one thread's multislice image of a particle from its rotated cube's spectrum, its linearisation and adjoint
! (the sweep and its checkpoints: doc/policies/3D/wave_model_policy.md, section 3)
module simple_wave_multislice
use, intrinsic :: iso_c_binding, only: c_associated, c_f_pointer, c_float_complex, c_int, c_ptr, c_size_t, c_null_ptr
use simple_defs,         only: PI
use simple_type_defs,    only: ctfparams
use simple_error,        only: simple_exception
use simple_fftw3,        only: fftwf_alloc_complex, fftwf_free, fftwf_plan_many_dft, fftwf_execute_dft, &
    &fftwf_destroy_plan, FFTW_FORWARD, FFTW_BACKWARD, FFTW_ESTIMATE
use simple_wave_field,   only: wave_field, wave_freq
use simple_wave_imaging, only: wave_imaging, WAVE_MSLICE
implicit none
private
#include "simple_local_flags.inc"

public :: wave_multislice

!> The multislice image of one particle from its rotated cube's native spectrum C(h,k,l): one transform along l
!! gives the sections, transmitted and propagated on the fine grid, then wave_imaging's lens and detector; the adjoint
!! sweeps back through segments of ceil(sqrt(nat)) sections from checkpoints (wave_model_policy.md, section 3)
type :: wave_multislice
    private
    integer :: nat = 0, R = 0, rl = 0, n = 0, osf = 0, kappa = 0, kseg = 0, nseg = 0
    complex :: c = cmplx(0.,0.)
    type(wave_imaging)   :: wimg                   !< the lens and detector of the exit wave
    type(wave_field)     :: fpsi                   !< the wave
    type(wave_field)     :: fw                     !< work field (transmissions, tangents, cotangents)
    complex, allocatable :: prop(:,:)              !< one-section propagator on the fine grid, FFT order
    complex, allocatable :: pmid(:,:)              !< exit to midplane
    complex, allocatable :: ck(:,:,:)              !< the wave entering each segment
    complex, allocatable :: segp(:,:,:)            !< a segment's waves after each transmission
    complex, allocatable :: segt(:,:,:)            !< a segment's transmissions
    complex, allocatable :: dpsi(:,:), tt(:,:), sj(:,:)
    integer, allocatable :: order(:)               !< the sections in the beam's order
    ! the sections' spectra along l and z: zs(1:nat,-R:R,-R:R), transformed along the first dimension
    type(c_ptr) :: zs_ptr = c_null_ptr, zd_ptr = c_null_ptr, plan_bwd = c_null_ptr, plan_fwd = c_null_ptr
    complex(c_float_complex), pointer, contiguous :: zs(:,:,:) => null()   !< of the point forward stored
    complex(c_float_complex), pointer, contiguous :: zd(:,:,:) => null()   !< of a direction or of the cotangents
    logical :: l_optics = .false., l_point = .false.
  contains
    procedure :: new
    procedure :: set_optics
    procedure :: forward
    procedure :: linearize
    procedure :: adjoint
    procedure :: get_exit_wave
    procedure :: get_n
    procedure :: kill
    procedure, private :: cube_to_sections
    procedure, private :: sections_to_cube_t
    procedure, private :: section_field
end type wave_multislice

integer, parameter :: ZS_POINT = 1, ZS_DIR = 2   !< which spectra: the point's (zs) or the direction's (zd)

contains

    !> nat the native box (even), osf the oversampling of the fine grid
    subroutine new( self, nat, osf )
        class(wave_multislice), intent(inout) :: self
        integer,                intent(in)    :: nat, osf
        complex(c_float_complex), pointer :: zin(:), zout(:)
        integer(c_size_t) :: nz
        integer(c_int)    :: dims(1), emb(1)
        call self%kill
        if( nat < 4 .or. mod(nat,2) /= 0 ) THROW_HARD('the native box must be even; wave_multislice%new')
        if( osf < 1 ) THROW_HARD('the oversampling must be at least 1; wave_multislice%new')
        self%nat  = nat
        self%R    = nat / 2
        self%rl   = self%R - 1
        self%osf  = osf
        self%n    = osf * nat
        self%kseg = ceiling(sqrt(real(nat)))
        self%nseg = (nat + self%kseg - 1) / self%kseg
        call self%wimg%new(nat, osf, WAVE_MSLICE)
        call self%fpsi%new(self%n)
        call self%fw%new(self%n)
        allocate(self%prop(self%n,self%n), self%pmid(self%n,self%n), self%dpsi(self%n,self%n), self%tt(self%n,self%n))
        allocate(self%ck(self%n,self%n,self%nseg), self%segp(self%n,self%n,self%kseg), self%segt(self%n,self%n,self%kseg))
        allocate(self%sj(-self%R:self%R,-self%R:self%R), self%order(nat))
        ! the sections' spectra, transformed along l (contiguous) for every (h,k): FFTW-aligned, one plan per sign
        nz = int(nat, c_size_t) * int(2*self%R+1, c_size_t)**2
        self%zs_ptr = fftwf_alloc_complex(nz)
        self%zd_ptr = fftwf_alloc_complex(nz)
        if( .not. c_associated(self%zs_ptr) .or. .not. c_associated(self%zd_ptr) ) &
            &THROW_HARD('FFTW allocation failed; wave_multislice%new')
        call c_f_pointer(self%zs_ptr, self%zs, [nat, 2*self%R+1, 2*self%R+1])
        call c_f_pointer(self%zd_ptr, self%zd, [nat, 2*self%R+1, 2*self%R+1])
        self%zs(1:,-self%R:,-self%R:) => self%zs
        self%zd(1:,-self%R:,-self%R:) => self%zd
        self%zs = cmplx(0.,0.)
        self%zd = cmplx(0.,0.)
        dims = int(nat, c_int)
        emb  = int(nat, c_int)
        ! in place: two pointers to the one buffer, as wave_field's cmat and cout
        call c_f_pointer(self%zs_ptr, zin,  [nz])
        call c_f_pointer(self%zs_ptr, zout, [nz])
        !$omp critical
        self%plan_bwd = fftwf_plan_many_dft(1_c_int, dims, int((2*self%R+1)**2, c_int), zin, emb, 1_c_int, &
            &int(nat, c_int), zout, emb, 1_c_int, int(nat, c_int), FFTW_BACKWARD, FFTW_ESTIMATE)
        self%plan_fwd = fftwf_plan_many_dft(1_c_int, dims, int((2*self%R+1)**2, c_int), zin, emb, 1_c_int, &
            &int(nat, c_int), zout, emb, 1_c_int, int(nat, c_int), FFTW_FORWARD, FFTW_ESTIMATE)
        !$omp end critical
        if( .not. c_associated(self%plan_bwd) .or. .not. c_associated(self%plan_fwd) ) &
            &THROW_HARD('FFTW plan creation failed; wave_multislice%new')
    end subroutine new

    !> The particle's lens, scattering coefficient and propagation: ctfparms at the native sampling, bfac the
    !! B-factor (A^2) of the lens envelope, kappa -1 or +1 the side the beam enters (ewald=neg, pos)
    subroutine set_optics( self, ctfparms, bfac, kappa )
        class(wave_multislice), intent(inout) :: self
        type(ctfparams),        intent(in)    :: ctfparms
        real,                   intent(in)    :: bfac
        integer,                intent(in)    :: kappa
        real    :: wl, fsq, ph, dmid
        integer :: i, j, hf, kf, js
        if( self%n == 0 ) THROW_HARD('new has not been called; wave_multislice%set_optics')
        if( abs(kappa) /= 1 ) THROW_HARD('multislice propagates through the sections: ewald must be neg or pos')
        call self%wimg%set_optics(ctfparms, bfac, 0.)
        self%c     = self%wimg%get_c()
        self%kappa = kappa
        wl         = self%wimg%get_wl()
        do js = 1, self%nat
            if( kappa == -1 )then
                self%order(js) = js
            else
                self%order(js) = self%nat - js + 1
            endif
        end do
        ! the midplane, section nat/2+1, from the exit: the last section the beam crossed
        dmid = -real(self%R - 1)
        if( kappa == 1 ) dmid = -real(self%R)
        do j = 1, self%n
            kf = wave_freq(j, self%n)
            do i = 1, self%n
                hf  = wave_freq(i, self%n)
                ! the fine grid's frequencies in native cycles per pixel, the section one native pixel thick
                fsq = real(hf*hf + kf*kf) / real(self%nat)**2
                ph  = -PI * wl * fsq
                self%prop(i,j) = cmplx(cos(ph), sin(ph))
                ph  = -PI * wl * dmid * fsq
                self%pmid(i,j) = cmplx(cos(ph), sin(ph))
            end do
        end do
        self%l_optics = .true.
        self%l_point  = .false.
    end subroutine set_optics

    !> ihat, the native DFT of the relative intensity, of the cube spectrum cube; stores the point
    subroutine forward( self, cube, ihat )
        class(wave_multislice), intent(inout) :: self
        complex,                intent(in)    :: cube(-self%R:self%R,-self%R:self%R,-self%rl:self%rl)
        complex,                intent(out)   :: ihat(-self%R:self%R,-self%R:self%R)
        complex, pointer :: psi(:,:), w(:,:)
        integer :: js
        if( .not. self%l_optics ) THROW_HARD('set_optics has not been called; wave_multislice%forward')
        psi => self%fpsi%ptr()
        w   => self%fw%ptr()
        call self%cube_to_sections(cube, ZS_POINT)
        psi = cmplx(1.,0.)
        do js = 1, self%nat
            if( mod(js-1, self%kseg) == 0 ) self%ck(:,:,(js-1)/self%kseg+1) = psi
            call self%section_field(ZS_POINT, self%order(js))
            psi = psi * exp(self%c * w)
            if( js < self%nat ) call propagate(self%fpsi, self%prop, .false.)
        end do
        call propagate(self%fpsi, self%pmid, .false.)
        call self%wimg%forward_wave(psi, ihat)
        self%l_point = .true.
    end subroutine forward

    !> d ihat = J dcube at the point forward stored (the tangent propagated beside the wave)
    subroutine linearize( self, dcube, dihat )
        class(wave_multislice), intent(inout) :: self
        complex,                intent(in)    :: dcube(-self%R:self%R,-self%R:self%R,-self%rl:self%rl)
        complex,                intent(out)   :: dihat(-self%R:self%R,-self%R:self%R)
        complex, pointer :: psi(:,:), w(:,:)
        integer :: js
        if( .not. self%l_point ) THROW_HARD('forward has not been called; wave_multislice%linearize')
        psi => self%fpsi%ptr()
        w   => self%fw%ptr()
        call self%cube_to_sections(dcube, ZS_DIR)
        psi       = cmplx(1.,0.)
        self%dpsi = cmplx(0.,0.)
        do js = 1, self%nat
            call self%section_field(ZS_POINT, self%order(js))
            self%tt = exp(self%c * w)
            call self%section_field(ZS_DIR, self%order(js))
            ! d(psi t) = t dpsi + psi t c du
            self%dpsi = self%tt * (self%dpsi + psi * self%c * w)
            psi       = psi * self%tt
            if( js < self%nat )then
                call propagate(self%fpsi, self%prop, .false.)
                w = self%dpsi
                call propagate(self%fw, self%prop, .false.)
                self%dpsi = w
            endif
        end do
        w = self%dpsi
        call propagate(self%fw, self%pmid, .false.)
        self%dpsi = w
        call self%wimg%linearize_wave(self%dpsi, dihat)
    end subroutine linearize

    !> dcube_adj = J^T q at the point forward stored, under Re<a,b>
    subroutine adjoint( self, q, dcube_adj )
        class(wave_multislice), intent(inout) :: self
        complex,                intent(in)    :: q(-self%R:self%R,-self%R:self%R)
        complex,                intent(out)   :: dcube_adj(-self%R:self%R,-self%R:self%R,-self%rl:self%rl)
        complex, pointer :: psi(:,:), w(:,:)
        integer :: sg, ja, jb, js, ii
        if( .not. self%l_point ) THROW_HARD('forward has not been called; wave_multislice%adjoint')
        psi => self%fpsi%ptr()
        w   => self%fw%ptr()
        ! the cotangent of the exit wave, carried back to the last section (the transpose of the midplane step)
        call self%wimg%adjoint_wave(q, self%dpsi)
        w = self%dpsi
        call propagate(self%fw, self%pmid, .true.)
        self%dpsi = w
        self%zd   = cmplx(0.,0.)
        do sg = self%nseg, 1, -1
            ja = (sg - 1) * self%kseg + 1
            jb = min(sg * self%kseg, self%nat)
            ! the segment's waves and transmissions again, from its checkpoint
            psi = self%ck(:,:,sg)
            do js = ja, jb
                ii = js - ja + 1
                call self%section_field(ZS_POINT, self%order(js))
                self%segt(:,:,ii) = exp(self%c * w)
                psi = psi * self%segt(:,:,ii)
                self%segp(:,:,ii) = psi
                if( js < jb ) call propagate(self%fpsi, self%prop, .false.)
            end do
            ! back through the segment: dpsi holds the cotangent of the wave after section js's propagation
            do js = jb, ja, -1
                ii = js - ja + 1
                if( js < self%nat )then
                    w = self%dpsi
                    call propagate(self%fw, self%prop, .true.)
                    self%dpsi = w
                endif
                ! psi+ = psi- t, t = exp(c u): d/du = conj(c) conj(psi+) lambda, lambda- = conj(t) lambda
                w = conjg(self%c) * conjg(self%segp(:,:,ii)) * self%dpsi
                self%dpsi = conjg(self%segt(:,:,ii)) * self%dpsi
                ! (bwd o embed)^T = embed_t o bwd^H = embed_t o fwd / n^2
                call self%fw%fwd
                call self%fw%embed_t(self%sj, self%R)
                self%zd(self%order(js),:,:) = self%sj / real(self%n)**2
            end do
        end do
        call self%sections_to_cube_t(dcube_adj)
    end subroutine adjoint

    !> z(j,h,k) = s_j(h,k) = (1/nat) sum_l cube(h,k,l) exp(+2 pi i l z_j / nat), z_j = j - 1 - nat/2, into zs
    !! (ZS_POINT) or zd (ZS_DIR): the l axis stored at index mod(l,nat)+1 with the factor (-1)^l of the shift by
    !! nat/2, transformed backward
    subroutine cube_to_sections( self, cube, which )
        class(wave_multislice), intent(inout) :: self
        complex,                intent(in)    :: cube(-self%R:self%R,-self%R:self%R,-self%rl:self%rl)
        integer,                intent(in)    :: which
        complex(c_float_complex), pointer, contiguous :: z(:,:,:)
        integer :: h, k, l
        if( which == ZS_POINT )then
            z => self%zs
        else
            z => self%zd
        endif
        z = cmplx(0.,0.)
        do k = -self%R, self%R
            do h = -self%R, self%R
                do l = -self%rl, self%rl
                    if( h*h + k*k + l*l > self%R**2 ) cycle
                    z(modulo(l, self%nat)+1,h,k) = real(1 - 2*modulo(l,2)) * cube(h,k,l) / real(self%nat)
                end do
            end do
        end do
        call fftwf_execute_dft(self%plan_bwd, z, z)
    end subroutine cube_to_sections

    !> the transpose of cube_to_sections under Re<a,b>, from the section cotangents in zd: forward along z, the
    !! factor (-1)^l / nat, the sphere
    subroutine sections_to_cube_t( self, cube )
        class(wave_multislice), intent(inout) :: self
        complex,                intent(out)   :: cube(-self%R:self%R,-self%R:self%R,-self%rl:self%rl)
        integer :: h, k, l
        call fftwf_execute_dft(self%plan_fwd, self%zd, self%zd)
        cube = cmplx(0.,0.)
        do l = -self%rl, self%rl
            do k = -self%R, self%R
                do h = -self%R, self%R
                    if( h*h + k*k + l*l > self%R**2 ) cycle
                    cube(h,k,l) = real(1 - 2*modulo(l,2)) * self%zd(modulo(l, self%nat)+1,h,k) / real(self%nat)
                end do
            end do
        end do
    end subroutine sections_to_cube_t

    !> fw receives u_j = bwd(embed(s_j)) of section j of zs (ZS_POINT) or zd (ZS_DIR)
    subroutine section_field( self, which, j )
        class(wave_multislice), intent(inout) :: self
        integer,                intent(in)    :: which, j
        if( which == ZS_POINT )then
            self%sj = self%zs(j,:,:)
        else
            self%sj = self%zd(j,:,:)
        endif
        call self%fw%embed(self%sj, self%R)
        call self%fw%bwd
    end subroutine section_field

    !> field <- bwd(p fwd(field)), or with conj(p) (the transpose: bwd and fwd are mutual transposes up to n^2)
    subroutine propagate( field, p, l_conj )
        class(wave_field), intent(inout) :: field
        complex,           intent(in)    :: p(:,:)
        logical,           intent(in)    :: l_conj
        complex, pointer :: f(:,:)
        f => field%ptr()
        call field%fwd
        if( l_conj )then
            f = f * conjg(p)
        else
            f = f * p
        endif
        call field%bwd
    end subroutine propagate

    !> the exit wave at the midplane of the point forward stored (fine grid, origin at the first point)
    subroutine get_exit_wave( self, psi_out )
        class(wave_multislice), intent(in)  :: self
        complex,                intent(out) :: psi_out(self%n,self%n)
        if( .not. self%l_point ) THROW_HARD('forward has not been called; wave_multislice%get_exit_wave')
        call self%wimg%get_exit_wave(psi_out)
    end subroutine get_exit_wave

    pure integer function get_n( self )
        class(wave_multislice), intent(in) :: self
        get_n = self%n
    end function get_n

    !> idempotent
    subroutine kill( self )
        class(wave_multislice), intent(inout) :: self
        if( c_associated(self%plan_bwd) )then
            !$omp critical
            call fftwf_destroy_plan(self%plan_bwd)
            !$omp end critical
        endif
        if( c_associated(self%plan_fwd) )then
            !$omp critical
            call fftwf_destroy_plan(self%plan_fwd)
            !$omp end critical
        endif
        if( c_associated(self%zs_ptr) ) call fftwf_free(self%zs_ptr)
        if( c_associated(self%zd_ptr) ) call fftwf_free(self%zd_ptr)
        nullify(self%zs, self%zd)
        self%zs_ptr   = c_null_ptr
        self%zd_ptr   = c_null_ptr
        self%plan_bwd = c_null_ptr
        self%plan_fwd = c_null_ptr
        call self%wimg%kill
        call self%fpsi%kill
        call self%fw%kill
        if( allocated(self%prop)  ) deallocate(self%prop)
        if( allocated(self%pmid)  ) deallocate(self%pmid)
        if( allocated(self%ck)    ) deallocate(self%ck)
        if( allocated(self%segp)  ) deallocate(self%segp)
        if( allocated(self%segt)  ) deallocate(self%segt)
        if( allocated(self%dpsi)  ) deallocate(self%dpsi)
        if( allocated(self%tt)    ) deallocate(self%tt)
        if( allocated(self%sj)    ) deallocate(self%sj)
        if( allocated(self%order) ) deallocate(self%order)
        self%nat = 0; self%R = 0; self%rl = 0; self%n = 0; self%osf = 0; self%kappa = 0
        self%kseg = 0; self%nseg = 0
        self%c = cmplx(0.,0.)
        self%l_optics = .false.
        self%l_point  = .false.
    end subroutine kill

end module simple_wave_multislice
