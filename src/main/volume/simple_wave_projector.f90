!@descr: the wave operator of a particle set (projection, first Born, Rytov, multislice) on the PCG lattice
! (doc/implementation_notes/planned/ghostbuster_operators.md, sections 3.1 to 3.5, 3.7 and 3.8)
module simple_wave_projector
use simple_defs,              only: dp, PI
use simple_type_defs,         only: ctfparams, ctfvars, CTFFLAG_FLIP
use simple_error,             only: simple_exception
use simple_ori,               only: ori
use simple_oris,              only: oris
use simple_sym,               only: sym
use simple_image,             only: image
use simple_ctf,               only: ctf
use simple_cartesian_fourier, only: center_embed_real3d
use simple_gridding,          only: kb_stencil_centered_crop_inv_envelope_1d
use simple_pcg_lattice,       only: pcg_lattice
use simple_wave_imaging,      only: wave_imaging, WAVE_PROJECTION, WAVE_FIRST_BORN, WAVE_RYTOV, WAVE_MSLICE
use simple_wave_multislice,   only: wave_multislice
!$ use omp_lib,               only: omp_get_max_threads, omp_get_thread_num, omp_in_parallel
implicit none
private
#include "simple_local_flags.inc"

public :: wave_projector

!> The image of a potential V (volts, native box) predicted by a single-scatter wave model, for every
!! particle of a selection, its derivative and transpose:
!!
!!   v -> v / env (KB deapodisation) -> pad, transform (SIMPLE's fft, origin at the box centre, 1/boxpd^3)
!!   S_i(h,k) = (1/nsym) sum_g gather(padf [h, k, kz] R_i S_g),  kz = curv_i (h^2+k^2), curv_i = ewald_curv
!!   s_i = nat^3 S_i exp(+2 pi i (h sx + k sy) / nat)         (the shift on the model side, as PCG's)
!!   p_i = -g ihat_i(s_i) / nat^2                              (wave_imaging; SIMPLE's particle units)
!!
!! g the gain sqrt(dose a^2) that maps relative intensity onto contrast-inverted, noise-normalised
!! particles (section 3.4); projection samples the central section, the curved models the paraboloid.
!! The misfit of a particle is 1/2 sum_(h,k) w_i(h,k) |y_i - p_i|^2, w_i = w_i inv_sig2(shell), zero at
!! the origin. The volume-to-plane map and its transpose (scatter at R_i S_g, Hermitian fold, inverse
!! transform, 1/boxpd^3, deapodisation) are exact transposes, so the gradient is the true gradient in
!! volts. Each thread runs the 2D chain of one particle with its own wave_imaging; the gradient planes of a
!! batch are scattered afterwards in particle order with the lattice's coloured scatter, and losses are
!! summed per particle in double precision in particle order, so results do not depend on the threads.
!! Multislice (section 3.7) gathers the particle's rotated cube instead of a section,
!!   C_i(h,k,l) = nat^3 (1/nsym) sum_g gather(padf [h, k, l] R_i S_g) exp(+2 pi i (h sx + k sy) / nat),
!! and runs wave_multislice on it; the transpose scatters the cube. Its gradient passes go in chunks of one
!! particle per thread, each chunk's cubes scattered in particle order before the next.
!! The operator runs its own parallel regions and is not called from inside another one.
type, extends(pcg_lattice) :: wave_projector
    private
    integer :: model = 0
    integer :: kappa = 0          !< -1, 0, +1 for ewald=neg, no, pos (projection: 0)
    integer :: osf   = 2          !< oversampling of the image formation
    real    :: gain  = 1.0        !< g = sqrt(dose a^2)
    real    :: bfac  = 0.0        !< B-factor of the lens envelope (A^2)
    integer :: fft_nthreads = 1
    type(image)          :: wimg                  !< persistent padded work image (keeps its FFTW plans)
    complex, allocatable :: cmat_x(:,:,:)         !< packed transform of the current volume
    complex, allocatable :: cmat_d(:,:,:)         !< packed transform of the search direction
    real,    allocatable :: invenv(:,:,:)         !< KB deapodisation of the native box
    integer              :: nsym = 1
    real,    allocatable :: symmats(:,:,:)        !< (3,3,nsym), symmats(:,:,1) = I
    ! ---- particles ----
    integer                      :: nptcls = 0
    real,            allocatable :: rotmats(:,:,:)   !< (3,3,nptcls)
    real,            allocatable :: shifts(:,:)      !< (2,nptcls), native pixels
    real,            allocatable :: curvs(:)         !< (nptcls) paraboloid coefficient
    type(ctfparams), allocatable :: ctfparms(:)      !< (nptcls), at the native sampling
    real,            allocatable :: inv_sig2(:)      !< (0:R) 1/sigma2 of the selection's half (1 without)
    real,            allocatable :: wts(:)           !< (nptcls) particle weights
    integer,         allocatable :: wshell(:,:)      !< (h,k) resolution shell, for the loss weights
    type(wave_imaging), allocatable :: imgs(:)       !< one per thread (the single-scatter models)
    type(wave_multislice), allocatable :: ms(:)      !< one per thread (multislice)
    ! ---- gradient accumulator ----
    complex, allocatable :: gacc(:,:,:)
    logical :: l_vol = .false., l_dir = .false., l_grad = .false., exists = .false.
  contains
    procedure :: new
    procedure :: set_sym
    procedure :: prep_particles
    procedure :: set_volume
    procedure :: set_direction
    procedure :: begin_grad
    procedure :: get_grad
    procedure :: predict
    procedure :: loss_grad
    procedure :: adjoint_planes
    procedure :: curvature
    procedure :: linearize_planes
    procedure :: line_losses
    procedure :: gain_terms
    procedure :: get_exit_wave
    procedure :: get_nptcls
    procedure :: get_model
    procedure, private :: transform
    procedure, private :: spectrum
    procedure, private :: spectrum_cube
    procedure, private :: prep_optics
    procedure, private :: ensure_imgs
    procedure, private :: pass
    procedure :: kill
end type wave_projector

integer, parameter :: PASS_LOSS_GRAD = 1, PASS_ADJOINT = 2, PASS_CURVATURE = 3, PASS_PREDICT = 4, PASS_LINEAR = 5
integer, parameter :: PASS_GAIN = 6

contains

    !> box, smpd: the native box and sampling; model one of WAVE_*; kappa -1, 0, +1 (ewald neg, no, pos;
    !! forced to 0 for projection); gain g = sqrt(dose a^2); osf the oversampling (2 unless measuring);
    !! bfac the lens B-factor (A^2, 0 by default)
    subroutine new( self, box, smpd, model, kappa, gain, osf, bfac, fft_nthreads, winsz )
        class(wave_projector), intent(inout) :: self
        integer,               intent(in)    :: box, model, kappa
        real,                  intent(in)    :: smpd, gain
        integer, optional,     intent(in)    :: osf, fft_nthreads
        real,    optional,     intent(in)    :: bfac, winsz
        real, allocatable :: inv1d(:)
        integer :: i, j, k, h
        call self%kill
        select case(model)
            case(WAVE_PROJECTION, WAVE_FIRST_BORN, WAVE_RYTOV, WAVE_MSLICE)
            case DEFAULT
                THROW_HARD('unknown wave model; wave_projector%new')
        end select
        if( abs(kappa) > 1 ) THROW_HARD('kappa must be -1, 0 or 1; wave_projector%new')
        if( model == WAVE_MSLICE .and. kappa == 0 ) THROW_HARD('multislice propagates through the box: kappa -1 or 1')
        if( gain <= 0. ) THROW_HARD('the gain must be positive; wave_projector%new')
        if( box < 2 .or. mod(box,2) /= 0 ) THROW_HARD('the box must be even; wave_projector%new')
        self%model = model
        self%kappa = kappa
        if( model == WAVE_PROJECTION ) self%kappa = 0
        self%gain  = gain
        self%osf   = 2
        if( present(osf) ) self%osf = osf
        self%bfac  = 0.
        if( present(bfac) ) self%bfac = bfac
        self%fft_nthreads = 1
        if( present(fft_nthreads) ) self%fft_nthreads = max(1, fft_nthreads)
        call self%new_lattice(box, smpd, self%fft_nthreads, winsz)
        call self%wimg%new([self%boxpd,self%boxpd,self%boxpd], smpd, &
            &wthreads=self%fft_nthreads > 1, fft_nthreads=self%fft_nthreads)
        ! the KB deapodisation of the native box, as the PCG backend's (separable, centred crop)
        call kb_stencil_centered_crop_inv_envelope_1d(self%kbwin, self%boxpd, self%box, inv1d)
        allocate(self%invenv(box,box,box))
        do k = 1, box
            do j = 1, box
                do i = 1, box
                    self%invenv(i,j,k) = inv1d(i) * inv1d(j) * inv1d(k)
                end do
            end do
        end do
        ! the shell of each plane point, for the loss weights
        allocate(self%wshell(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                self%wshell(h,k) = min(nint(sqrt(real(h*h + k*k))), self%Rnat)
            end do
        end do
        self%nsym = 1
        allocate(self%symmats(3,3,1), source=0.)
        self%symmats(1,1,1) = 1.; self%symmats(2,2,1) = 1.; self%symmats(3,3,1) = 1.
        allocate(self%gacc(self%lims3(1,1):self%lims3(1,2), self%lims3(2,1):self%lims3(2,2), &
            &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        self%exists = .true.
    end subroutine new

    !> point-group replication: the model sees the symmetrised volume, the mean of the sections at R_i S_g
    subroutine set_sym( self, pgrpsyms )
        class(wave_projector), intent(inout) :: self
        class(sym),            intent(in)    :: pgrpsyms
        integer :: g
        if( allocated(self%symmats) ) deallocate(self%symmats)
        self%nsym = pgrpsyms%get_nsym()
        allocate(self%symmats(3,3,self%nsym), source=0.)
        do g = 1, self%nsym
            call pgrpsyms%get_sym_rmat(g, self%symmats(:,:,g))
        end do
    end subroutine set_sym

    !> the particles' orientations, shifts (native pixels) and CTF (at the native sampling) from os, as
    !! rec_obs_oris gives them; inv_sig2(0:) the selection's noise spectrum at the native shells, weights
    !! one per particle
    subroutine prep_particles( self, os, inv_sig2, weights )
        class(wave_projector), intent(inout) :: self
        class(oris),           intent(inout) :: os
        real, optional,        intent(in)    :: inv_sig2(0:)
        real, optional,        intent(in)    :: weights(:)
        type(ori)     :: o
        type(ctf)     :: tfun
        type(ctfvars) :: vars
        integer :: i, rtop
        if( .not. self%exists ) THROW_HARD('new has not been called; wave_projector%prep_particles')
        self%nptcls = os%get_noris()
        if( allocated(self%rotmats)  ) deallocate(self%rotmats)
        if( allocated(self%shifts)   ) deallocate(self%shifts)
        if( allocated(self%curvs)    ) deallocate(self%curvs)
        if( allocated(self%ctfparms) ) deallocate(self%ctfparms)
        if( allocated(self%inv_sig2) ) deallocate(self%inv_sig2)
        if( allocated(self%wts)      ) deallocate(self%wts)
        allocate(self%rotmats(3,3,self%nptcls), self%shifts(2,self%nptcls), self%curvs(self%nptcls), &
            &self%ctfparms(self%nptcls))
        allocate(self%inv_sig2(0:self%Rnat), source=1.)
        allocate(self%wts(self%nptcls), source=1.)
        call o%new(.true.)
        do i = 1, self%nptcls
            call os%get_ori(i, o)
            self%rotmats(:,:,i) = o%get_mat()
            self%shifts(:,i)    = o%get_2Dshift()
            self%ctfparms(i)    = o%get_ctfvars()
            if( abs(self%ctfparms(i)%smpd - self%smpd) > 1.e-4 * self%smpd )then
                THROW_HARD('the particles are not at the operator sampling; wave_projector%prep_particles')
            endif
            ! what the image formation refuses, refused before its parallel regions
            if( self%ctfparms(i)%ctfflag == CTFFLAG_FLIP ) THROW_HARD('phase-flipped particles have no wave model')
            if( self%ctfparms(i)%fraca < 0. .or. self%ctfparms(i)%fraca >= 1. ) THROW_HARD('fraca must be in [0,1)')
            if( self%kappa == 0 .or. self%model == WAVE_MSLICE )then
                ! no section: multislice gathers the cube
                self%curvs(i) = 0.
            else
                tfun = ctf(self%ctfparms(i)%smpd, self%ctfparms(i)%kv, self%ctfparms(i)%cs, self%ctfparms(i)%fraca)
                vars = tfun%get_ctfvars(0.)
                self%curvs(i) = self%ewald_curv(self%kappa, vars%wl)
            endif
        end do
        call o%kill
        if( present(inv_sig2) )then
            rtop = min(self%Rnat, ubound(inv_sig2,1))
            self%inv_sig2(0:rtop) = inv_sig2(0:rtop)
        endif
        if( present(weights) )then
            if( size(weights) /= self%nptcls ) THROW_HARD('one weight per particle; wave_projector%prep_particles')
            self%wts = weights
        endif
    end subroutine prep_particles

    !> the current point: the potential in volts on the native box
    subroutine set_volume( self, v )
        class(wave_projector), intent(inout) :: self
        real,                  intent(in)    :: v(self%box,self%box,self%box)
        if( .not. self%exists ) THROW_HARD('new has not been called; wave_projector%set_volume')
        call self%transform(v, self%cmat_x)
        self%l_vol = .true.
    end subroutine set_volume

    !> the search direction d (volts) for curvature and line_losses
    subroutine set_direction( self, d )
        class(wave_projector), intent(inout) :: self
        real,                  intent(in)    :: d(self%box,self%box,self%box)
        if( .not. self%exists ) THROW_HARD('new has not been called; wave_projector%set_direction')
        call self%transform(d, self%cmat_d)
        self%l_dir = .true.
    end subroutine set_direction

    subroutine transform( self, v, cmat )
        class(wave_projector), intent(inout) :: self
        real,                  intent(in)    :: v(self%box,self%box,self%box)
        complex, allocatable,  intent(inout) :: cmat(:,:,:)
        call self%wimg%set_rmat(center_embed_real3d(v * self%invenv, self%boxpd), .false.)
        call self%wimg%fft()
        cmat = self%wimg%get_cmat()
    end subroutine transform

    !> zeroes the gradient accumulator
    subroutine begin_grad( self )
        class(wave_projector), intent(inout) :: self
        self%gacc   = cmplx(0.,0.)
        self%l_grad = .true.
    end subroutine begin_grad

    !> the accumulated gradient in volts: the transpose of the volume-to-plane map applied to the
    !! scattered planes (Hermitian fold, inverse transform, 1/boxpd^3, deapodisation)
    subroutine get_grad( self, g )
        class(wave_projector), intent(inout) :: self
        real,                  intent(out)   :: g(self%box,self%box,self%box)
        if( .not. self%l_grad ) THROW_HARD('begin_grad has not been called; wave_projector%get_grad')
        g = self%pcg_lattice%fold_lattice_accum(self%gacc, self%wimg, .true.)
        g = g * self%invenv / real(self%boxpd)**3
    end subroutine get_grad

    !> s_i of the current volume (l_dir = .false.) or of the direction (.true.), on the native disk
    subroutine spectrum( self, i, l_dir, s )
        class(wave_projector), intent(in)    :: self
        integer,               intent(in)    :: i
        logical,               intent(in)    :: l_dir
        complex,               intent(out)   :: s(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2))
        complex, allocatable :: plane(:,:)
        real    :: rot(3,3), arg, scale
        integer :: g, h, k
        allocate(plane(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        s = cmplx(0.,0.)
        do g = 1, self%nsym
            rot = matmul(self%rotmats(:,:,i), self%symmats(:,:,g))
            if( l_dir )then
                call self%gather_plane(self%cmat_d, rot, self%curvs(i), plane)
            else
                call self%gather_plane(self%cmat_x, rot, self%curvs(i), plane)
            endif
            s = s + plane
        end do
        scale = real(self%box)**3 / real(self%nsym)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp ) cycle
                arg    = 2.0*PI * (real(h)*self%shifts(1,i) + real(k)*self%shifts(2,i)) / real(self%box)
                s(h,k) = scale * s(h,k) * cmplx(cos(arg), sin(arg))
            end do
        end do
        deallocate(plane)
    end subroutine spectrum

    !> C_i of the current volume (l_dir = .false.) or of the direction (.true.): the rotated cube's spectrum
    subroutine spectrum_cube( self, i, l_dir, cube )
        class(wave_projector), intent(in)  :: self
        integer,               intent(in)  :: i
        logical,               intent(in)  :: l_dir
        complex,               intent(out) :: cube(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2),&
                                                  &-(self%Rnat-1):self%Rnat-1)
        complex, allocatable :: part(:,:,:)
        complex :: ph
        real    :: rot(3,3), arg, scale
        integer :: g, h, k
        ! the first rotation straight into the cube, a second buffer only for a point group
        rot = matmul(self%rotmats(:,:,i), self%symmats(:,:,1))
        if( l_dir )then
            call self%gather_cube(self%cmat_d, rot, cube)
        else
            call self%gather_cube(self%cmat_x, rot, cube)
        endif
        if( self%nsym > 1 ) allocate(part, mold=cube)
        do g = 2, self%nsym
            rot = matmul(self%rotmats(:,:,i), self%symmats(:,:,g))
            if( l_dir )then
                call self%gather_cube(self%cmat_d, rot, part)
            else
                call self%gather_cube(self%cmat_x, rot, part)
            endif
            cube = cube + part
        end do
        scale = real(self%box)**3 / real(self%nsym)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp ) cycle
                arg = 2.0*PI * (real(h)*self%shifts(1,i) + real(k)*self%shifts(2,i)) / real(self%box)
                ph  = scale * cmplx(cos(arg), sin(arg))
                cube(h,k,:) = ph * cube(h,k,:)
            end do
        end do
        if( allocated(part) ) deallocate(part)
    end subroutine spectrum_cube

    subroutine prep_optics( self, ithr, i )
        class(wave_projector), intent(inout) :: self
        integer,               intent(in)    :: ithr, i
        real :: depth
        ! Rytov on the paraboloid exponentiates at the exit face, half the box downstream (the specimen
        ! fills the box); without curvature the model has no depth
        if( self%model == WAVE_MSLICE )then
            call self%ms(ithr)%set_optics(self%ctfparms(i), self%bfac, self%kappa)
            return
        endif
        depth = 0.
        if( self%model == WAVE_RYTOV .and. self%kappa /= 0 ) depth = real(self%box) / 2.
        call self%imgs(ithr)%set_optics(self%ctfparms(i), self%bfac, depth)
    end subroutine prep_optics

    !> one wave_imaging per thread, sized at the point of use; every threaded entry passes here first, so the
    !! check that no caller's parallel region encloses the operator's own lives here
    subroutine ensure_imgs( self )
        class(wave_projector), intent(inout) :: self
        integer :: nthr, t
        !$ if( omp_in_parallel() ) THROW_HARD('wave_projector runs its own parallel regions, not inside another')
        nthr = 1
        !$ nthr = omp_get_max_threads()
        if( self%model == WAVE_MSLICE )then
            if( allocated(self%ms) )then
                if( size(self%ms) >= nthr ) return
                do t = 1, size(self%ms)
                    call self%ms(t)%kill
                end do
                deallocate(self%ms)
            endif
            allocate(self%ms(nthr))
            do t = 1, nthr
                call self%ms(t)%new(self%box, self%osf)
            end do
            return
        endif
        if( allocated(self%imgs) )then
            if( size(self%imgs) >= nthr ) return
            do t = 1, size(self%imgs)
                call self%imgs(t)%kill
            end do
            deallocate(self%imgs)
        endif
        allocate(self%imgs(nthr))
        do t = 1, nthr
            call self%imgs(t)%new(self%box, self%osf, self%model)
        end do
    end subroutine ensure_imgs

    !> p_i of particles ifrom..ifrom+nb-1 of the selection at the current volume
    subroutine predict( self, nb, ifrom, planes )
        class(wave_projector), intent(inout) :: self
        integer,               intent(in)    :: nb, ifrom
        complex,               intent(out)   :: planes(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2),nb)
        real(dp) :: dummy(nb)
        call self%pass(PASS_PREDICT, nb, ifrom, dummy, pout=planes)
    end subroutine predict

    !> the misfit of the batch at the current volume, its gradient planes scattered into the accumulator
    real(dp) function loss_grad( self, y_batch, nb, ifrom ) result( loss )
        class(wave_projector), intent(inout) :: self
        integer,               intent(in)    :: nb, ifrom
        complex,               intent(in)    :: y_batch(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2),nb)
        real(dp) :: losses(nb)
        integer  :: ib
        call self%pass(PASS_LOSS_GRAD, nb, ifrom, losses, yin=y_batch)
        loss = 0.0_dp
        do ib = 1, nb
            loss = loss + losses(ib)
        end do
    end function loss_grad

    !> J_i^T q_i for particles ifrom..ifrom+nb-1 at the current volume (q in the units of p, real inner
    !! product), scattered into the accumulator
    subroutine adjoint_planes( self, q_batch, nb, ifrom )
        class(wave_projector), intent(inout) :: self
        integer,               intent(in)    :: nb, ifrom
        complex,               intent(in)    :: q_batch(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2),nb)
        real(dp) :: dummy(nb)
        call self%pass(PASS_ADJOINT, nb, ifrom, dummy, yin=q_batch)
    end subroutine adjoint_planes

    !> the Gauss-Newton curvature of the batch's misfit along the direction: sum_i sum w_i |J_i d|^2
    real(dp) function curvature( self, nb, ifrom ) result( curv )
        class(wave_projector), intent(inout) :: self
        integer,               intent(in)    :: nb, ifrom
        real(dp) :: curvs(nb)
        integer  :: ib
        call self%pass(PASS_CURVATURE, nb, ifrom, curvs)
        curv = 0.0_dp
        do ib = 1, nb
            curv = curv + curvs(ib)
        end do
    end function curvature

    !> J_i d of particles ifrom..ifrom+nb-1: the change of p_i along the direction, at the current volume
    subroutine linearize_planes( self, nb, ifrom, planes )
        class(wave_projector), intent(inout) :: self
        integer,               intent(in)    :: nb, ifrom
        complex,               intent(out)   :: planes(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2),nb)
        real(dp) :: dummy(nb)
        call self%pass(PASS_LINEAR, nb, ifrom, dummy, pout=planes)
    end subroutine linearize_planes

    !> the batch's misfit at x + t(j) d for every step t(j): the sections are linear in the volume, so
    !! s(x + t d) = s(x) + t s(d) without a new transform
    subroutine line_losses( self, y_batch, nb, ifrom, t, losses )
        class(wave_projector), intent(inout) :: self
        integer,               intent(in)    :: nb, ifrom
        complex,               intent(in)    :: y_batch(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2),nb)
        real(dp),              intent(in)    :: t(:)
        real(dp),              intent(out)   :: losses(:)
        complex, allocatable :: sx(:,:), sd(:,:), st(:,:), ihat(:,:), cx(:,:,:), cd(:,:,:), ct(:,:,:)
        real(dp), allocatable :: lossmat(:,:)
        integer :: ib, i, ithr, j
        if( .not. (self%l_vol .and. self%l_dir) ) THROW_HARD('set_volume and set_direction first; line_losses')
        if( size(losses) /= size(t) ) THROW_HARD('one loss per step; line_losses')
        if( nb < 1 ) return
        if( ifrom < 1 .or. ifrom + nb - 1 > self%nptcls ) THROW_HARD('batch outside the particle selection; line_losses')
        call self%ensure_imgs
        allocate(lossmat(size(t),nb))
        !$omp parallel do default(shared) private(ib,i,ithr,j,sx,sd,st,ihat,cx,cd,ct) schedule(static) proc_bind(close)
        do ib = 1, nb
            i    = ifrom + ib - 1
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            allocate(ihat(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
            call self%prep_optics(ithr, i)
            if( self%model == WAVE_MSLICE )then
                allocate(cx(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2), -(self%Rnat-1):self%Rnat-1))
                allocate(cd, ct, mold=cx)
                call self%spectrum_cube(i, .false., cx)
                call self%spectrum_cube(i, .true.,  cd)
                do j = 1, size(t)
                    ct = cx + real(t(j)) * cd
                    call self%ms(ithr)%forward(ct, ihat)
                    lossmat(j,ib) = misfit(self, i, y_batch(:,:,ib), ihat)
                end do
                deallocate(cx, cd, ct)
            else
                allocate(sx, sd, st, mold=ihat)
                call self%spectrum(i, .false., sx)
                call self%spectrum(i, .true.,  sd)
                do j = 1, size(t)
                    st = sx + real(t(j)) * sd
                    call self%imgs(ithr)%forward(st, ihat)
                    lossmat(j,ib) = misfit(self, i, y_batch(:,:,ib), ihat)
                end do
                deallocate(sx, sd, st)
            endif
            deallocate(ihat)
        end do
        !$omp end parallel do
        do j = 1, size(t)
            losses(j) = 0.0_dp
            do ib = 1, nb
                losses(j) = losses(j) + lossmat(j,ib)
            end do
        end do
    end subroutine line_losses

    !> The terms of the one-scalar fit of the data's nonlinear part (section 3.4 of the plan): with p_i the
    !! prediction at the current volume and l_i = J_i(0) v its linear part (the weak-phase image of the same
    !! volume), q_i = p_i - l_i; num(ib) = sum w Re(conj(y_i - l_i) q_i) and den(ib) = sum w |q_i|^2. The data
    !! fix g v through the linear part, so y - l = (g / g_data) q, and the gain the data prefer is
    !! g sum(den) / sum(num)
    subroutine gain_terms( self, y_batch, nb, ifrom, num, den )
        class(wave_projector), intent(inout) :: self
        integer,               intent(in)    :: nb, ifrom
        complex,               intent(in)    :: y_batch(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2),nb)
        real(dp),              intent(out)   :: num(nb), den(nb)
        call self%pass(PASS_GAIN, nb, ifrom, num, yin=y_batch, vals2=den)
    end subroutine gain_terms

    !> the fine-grid exit wave psi of particle i at the current volume (real space, origin at its first point)
    subroutine get_exit_wave( self, i, psi )
        class(wave_projector), intent(inout) :: self
        integer,               intent(in)    :: i
        complex,               intent(out)   :: psi(:,:)
        complex, allocatable :: s(:,:), ihat(:,:), cube(:,:,:)
        if( .not. self%l_vol ) THROW_HARD('set_volume has not been called; get_exit_wave')
        call self%ensure_imgs
        if( i < 1 .or. i > self%nptcls ) THROW_HARD('particle outside the selection; get_exit_wave')
        if( size(psi,1) /= self%osf*self%box .or. size(psi,2) /= self%osf*self%box )then
            THROW_HARD('exit wave array size; get_exit_wave')
        endif
        allocate(s(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        allocate(ihat, mold=s)
        call self%prep_optics(1, i)
        if( self%model == WAVE_MSLICE )then
            allocate(cube(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2), -(self%Rnat-1):self%Rnat-1))
            call self%spectrum_cube(i, .false., cube)
            call self%ms(1)%forward(cube, ihat)
            call self%ms(1)%get_exit_wave(psi)
        else
            call self%spectrum(i, .false., s)
            call self%imgs(1)%forward(s, ihat)
            call self%imgs(1)%get_exit_wave(psi)
        endif
    end subroutine get_exit_wave

    !> one pass over a batch: every particle's chain on its thread, then (for the gradient passes) the
    !! scatter of its plane or cube in particle order. yin holds y (loss_grad) or q (adjoint), pout receives p
    !! (predict) or J d (linearize_planes). The single-scatter models run the batch as one chunk; multislice in
    !! chunks of one particle per thread, which bounds the gradient cubes held for the scatter
    subroutine pass( self, mode, nb, ifrom, vals, yin, pout, vals2 )
        class(wave_projector), intent(inout) :: self
        integer,               intent(in)    :: mode, nb, ifrom
        real(dp),              intent(out)   :: vals(nb)
        real(dp), optional,    intent(out)   :: vals2(nb)
        complex, optional,     intent(in)    :: yin(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2),nb)
        complex, optional,     intent(out)   :: pout(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2),nb)
        complex, allocatable :: s(:,:), ihat(:,:), q(:,:), sadj(:,:), ds(:,:), dihat(:,:)
        complex, allocatable :: cx(:,:,:), cd(:,:,:)
        complex, allocatable :: gplanes(:,:,:), gcubes(:,:,:,:)
        real    :: rots(3,3,self%nsym), arg, pscale, scale
        integer :: ib, ib0, ib1, i, ithr, h, k, g, nchunk, rl
        logical :: l_ms, l_grad
        if( .not. self%l_vol ) THROW_HARD('set_volume has not been called; wave_projector pass')
        if( nb < 1 ) return
        if( ifrom < 1 .or. ifrom + nb - 1 > self%nptcls ) THROW_HARD('batch outside the particle selection')
        if( (mode == PASS_CURVATURE .or. mode == PASS_LINEAR) .and. .not. self%l_dir )then
            THROW_HARD('set_direction has not been called; curvature, linearize_planes')
        endif
        if( (mode == PASS_LOSS_GRAD .or. mode == PASS_ADJOINT) .and. .not. self%l_grad )then
            THROW_HARD('begin_grad has not been called; wave_projector pass')
        endif
        if( mode == PASS_GAIN .and. .not. (present(yin) .and. present(vals2)) )then
            THROW_HARD('gain_terms passes its planes and both sums; wave_projector pass')
        endif
        if( (mode == PASS_LOSS_GRAD .or. mode == PASS_ADJOINT) .and. .not. present(yin) )then
            THROW_HARD('loss_grad and adjoint_planes pass their planes; wave_projector pass')
        endif
        if( (mode == PASS_PREDICT .or. mode == PASS_LINEAR) .and. .not. present(pout) )then
            THROW_HARD('predict and linearize_planes receive their planes; wave_projector pass')
        endif
        call self%ensure_imgs
        l_ms   = self%model == WAVE_MSLICE
        l_grad = mode == PASS_LOSS_GRAD .or. mode == PASS_ADJOINT
        rl     = self%Rnat - 1
        pscale = -self%gain / real(self%box)**2
        scale  = real(self%box)**3 / real(self%nsym)
        nchunk = nb
        if( l_ms )then
            nchunk = 1
            !$ nchunk = omp_get_max_threads()
            nchunk = min(nb, nchunk)
        endif
        if( l_grad )then
            if( l_ms )then
                allocate(gcubes(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2), -rl:rl, nchunk))
            else
                allocate(gplanes(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2), nb))
            endif
        endif
        vals = 0.0_dp
        if( present(vals2) ) vals2 = 0.0_dp
        do ib0 = 1, nb, nchunk
            ib1 = min(nb, ib0 + nchunk - 1)
            !$omp parallel do default(shared) private(ib,i,ithr,h,k,arg,s,ihat,q,sadj,ds,dihat,cx,cd) &
            !$omp schedule(static) proc_bind(close)
            do ib = ib0, ib1
                i    = ifrom + ib - 1
                ithr = 1
                !$ ithr = omp_get_thread_num() + 1
                allocate(s(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
                allocate(ihat, q, sadj, ds, dihat, mold=s)
                call self%prep_optics(ithr, i)
                if( l_ms )then
                    allocate(cx(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2), -rl:rl))
                    ! the direction's cube only where a pass uses one
                    if( mode == PASS_LINEAR .or. mode == PASS_CURVATURE .or. mode == PASS_GAIN ) allocate(cd, mold=cx)
                    call self%spectrum_cube(i, .false., cx)
                    call self%ms(ithr)%forward(cx, ihat)
                else
                    call self%spectrum(i, .false., s)
                    call self%imgs(ithr)%forward(s, ihat)
                endif
                select case(mode)
                    case(PASS_PREDICT)
                        pout(:,:,ib) = pscale * ihat
                    case(PASS_LOSS_GRAD, PASS_ADJOINT)
                        if( mode == PASS_LOSS_GRAD )then
                            vals(ib) = misfit(self, i, yin(:,:,ib), ihat)
                            ! dL/dp = -w (y - p), in the units of p
                            q = cmplx(0.,0.)
                            do k = self%lims2(2,1), self%lims2(2,2)
                                do h = self%lims2(1,1), self%lims2(1,2)
                                    if( h*h + k*k > self%sqlp ) cycle
                                    q(h,k) = -weight(self, i, h, k) * (yin(h,k,ib) - pscale * ihat(h,k))
                                end do
                            end do
                        else
                            q = yin(:,:,ib)
                        endif
                        ! p = pscale ihat: J_p^T q = J_ihat^T (pscale q)
                        q = pscale * q
                        ! the transpose of the spectrum's nat^3/nsym exp(+i phi): conjugate phase, the same scale
                        if( l_ms )then
                            ! straight into the chunk's gradient cube, the phase and scale applied in place
                            call self%ms(ithr)%adjoint(q, gcubes(:,:,:,ib-ib0+1))
                            do k = self%lims2(2,1), self%lims2(2,2)
                                do h = self%lims2(1,1), self%lims2(1,2)
                                    if( h*h + k*k > self%sqlp )then
                                        gcubes(h,k,:,ib-ib0+1) = cmplx(0.,0.)
                                        cycle
                                    endif
                                    arg = 2.0*PI * (real(h)*self%shifts(1,i) + real(k)*self%shifts(2,i)) / real(self%box)
                                    gcubes(h,k,:,ib-ib0+1) = (scale * cmplx(cos(arg), -sin(arg))) * gcubes(h,k,:,ib-ib0+1)
                                end do
                            end do
                        else
                            call self%imgs(ithr)%adjoint(q, sadj)
                            gplanes(:,:,ib) = cmplx(0.,0.)
                            do k = self%lims2(2,1), self%lims2(2,2)
                                do h = self%lims2(1,1), self%lims2(1,2)
                                    if( h*h + k*k > self%sqlp ) cycle
                                    arg = 2.0*PI * (real(h)*self%shifts(1,i) + real(k)*self%shifts(2,i)) / real(self%box)
                                    gplanes(h,k,ib) = scale * sadj(h,k) * cmplx(cos(arg), -sin(arg))
                                end do
                            end do
                        endif
                    case(PASS_LINEAR, PASS_CURVATURE)
                        if( l_ms )then
                            call self%spectrum_cube(i, .true., cd)
                            call self%ms(ithr)%linearize(cd, dihat)
                        else
                            call self%spectrum(i, .true., ds)
                            call self%imgs(ithr)%linearize(ds, dihat)
                        endif
                        if( mode == PASS_LINEAR )then
                            pout(:,:,ib) = pscale * dihat
                        else
                            vals(ib) = 0.0_dp
                            do k = self%lims2(2,1), self%lims2(2,2)
                                do h = self%lims2(1,1), self%lims2(1,2)
                                    if( h*h + k*k > self%sqlp ) cycle
                                    vals(ib) = vals(ib) + real(weight(self, i, h, k), dp) * &
                                        &real(pscale,dp)**2 * real(abs(dihat(h,k)),dp)**2
                                end do
                            end do
                        endif
                    case(PASS_GAIN)
                        ! the linear part J(0) s: the chain linearised at zero
                        if( l_ms )then
                            cd = cmplx(0.,0.)
                            call self%ms(ithr)%forward(cd, sadj)
                            call self%ms(ithr)%linearize(cx, dihat)
                        else
                            ds = cmplx(0.,0.)
                            call self%imgs(ithr)%forward(ds, sadj)
                            call self%imgs(ithr)%linearize(s, dihat)
                        endif
                        do k = self%lims2(2,1), self%lims2(2,2)
                            do h = self%lims2(1,1), self%lims2(1,2)
                                if( h*h + k*k > self%sqlp ) cycle
                                q(h,k)    = pscale * (ihat(h,k) - dihat(h,k))
                                vals(ib)  = vals(ib) + real(weight(self, i, h, k), dp) * &
                                    &real(real(conjg(yin(h,k,ib) - pscale * dihat(h,k)) * q(h,k)), dp)
                                vals2(ib) = vals2(ib) + real(weight(self, i, h, k), dp) * real(abs(q(h,k)),dp)**2
                            end do
                        end do
                end select
                deallocate(s, ihat, q, sadj, ds, dihat)
                if( allocated(cx) ) deallocate(cx)
                if( allocated(cd) ) deallocate(cd)
            end do
            !$omp end parallel do
            if( l_grad )then
                ! the transposes of the gathers, in particle order: G^T at every R_i S_g
                do ib = ib0, ib1
                    i = ifrom + ib - 1
                    do g = 1, self%nsym
                        rots(:,:,g) = matmul(self%rotmats(:,:,i), self%symmats(:,:,g))
                    end do
                    if( l_ms )then
                        call self%scatter_cube(gcubes(:,:,:,ib-ib0+1), rots, self%gacc)
                    else
                        call self%scatter_plane(gplanes(:,:,ib), rots, self%curvs(i), self%gacc)
                    endif
                end do
            endif
        end do
        if( allocated(gplanes) ) deallocate(gplanes)
        if( allocated(gcubes)  ) deallocate(gcubes)
    end subroutine pass

    !> the loss weight of plane point (h,k) of particle i: inv_sig2(shell) w_i, zero at the origin
    pure real function weight( self, i, h, k )
        class(wave_projector), intent(in) :: self
        integer,               intent(in) :: i, h, k
        if( h == 0 .and. k == 0 )then
            weight = 0.
        else
            weight = self%inv_sig2(self%wshell(h,k)) * self%wts(i)
        endif
    end function weight

    !> 1/2 sum w |y - p|^2 in double precision, p = -g ihat / nat^2
    real(dp) function misfit( self, i, y, ihat )
        class(wave_projector), intent(in) :: self
        integer,               intent(in) :: i
        complex,               intent(in) :: y(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2))
        complex,               intent(in) :: ihat(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2))
        real    :: pscale
        integer :: h, k
        pscale = -self%gain / real(self%box)**2
        misfit = 0.0_dp
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp ) cycle
                misfit = misfit + real(weight(self, i, h, k), dp) * real(abs(y(h,k) - pscale * ihat(h,k)),dp)**2
            end do
        end do
        misfit = 0.5_dp * misfit
    end function misfit

    pure integer function get_nptcls( self )
        class(wave_projector), intent(in) :: self
        get_nptcls = self%nptcls
    end function get_nptcls

    pure integer function get_model( self )
        class(wave_projector), intent(in) :: self
        get_model = self%model
    end function get_model

    !> idempotent
    subroutine kill( self )
        class(wave_projector), intent(inout) :: self
        integer :: t
        if( allocated(self%imgs) )then
            do t = 1, size(self%imgs)
                call self%imgs(t)%kill
            end do
            deallocate(self%imgs)
        endif
        if( allocated(self%ms) )then
            do t = 1, size(self%ms)
                call self%ms(t)%kill
            end do
            deallocate(self%ms)
        endif
        if( self%exists ) call self%wimg%kill
        if( allocated(self%cmat_x)   ) deallocate(self%cmat_x)
        if( allocated(self%cmat_d)   ) deallocate(self%cmat_d)
        if( allocated(self%invenv)   ) deallocate(self%invenv)
        if( allocated(self%symmats)  ) deallocate(self%symmats)
        if( allocated(self%rotmats)  ) deallocate(self%rotmats)
        if( allocated(self%shifts)   ) deallocate(self%shifts)
        if( allocated(self%curvs)    ) deallocate(self%curvs)
        if( allocated(self%ctfparms) ) deallocate(self%ctfparms)
        if( allocated(self%inv_sig2) ) deallocate(self%inv_sig2)
        if( allocated(self%wts)      ) deallocate(self%wts)
        if( allocated(self%wshell)   ) deallocate(self%wshell)
        if( allocated(self%gacc)     ) deallocate(self%gacc)
        call self%kill_lattice
        self%nptcls = 0; self%nsym = 1; self%model = 0; self%kappa = 0
        self%osf = 2; self%gain = 1.0; self%bfac = 0.0; self%fft_nthreads = 1
        self%l_vol = .false.; self%l_dir = .false.; self%l_grad = .false.; self%exists = .false.
    end subroutine kill

end module simple_wave_projector
