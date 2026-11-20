!@descr: continuous view refinement on polar sections: B-spline sections cut from cartft_calc's references, their view derivatives, the per-thread section cache and the five-variable cost of a drawn pose
! The polar side of polar_continuous_view.md. A section is cut on polarft_calc's polar grid from cartft_calc's B-spline
! references (gather_ref) at the frame F = T(w) R_ref, T(w) = exp([wx, wy, 0]x) a tilt about the drawn reference's in-plane
! axes; the pose is Rz(psi) T(w) R_ref. The cost's variables are the joint shift search's (shift x, shift y, rotation
! index) and the two tilts (rad), the tilt columns the exact derivatives of the section in that chart. polarft_calc
! correlates the supplied section; cartft_calc and polarft_calc do not import each other.
module simple_polarft_cont_calc
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api, only: dp, sp, simple_exception
use simple_polarft_calc,    only: polarft_calc
use simple_cartft_calc,     only: cartft_calc
!$ use omp_lib,             only: omp_get_thread_num, omp_get_max_threads
implicit none
private
public :: polarft_cont_calc
public :: run_all_polarft_cont_calc_tests

#include "simple_local_flags.inc"

!> one thread's drawn reference and the section it last cut: reused while the view variables stay unchanged
type :: cont_sect_cache
    complex(sp), allocatable :: sect(:,:)        !< (pftsz, kfrom:kto) the section at the cached view
    complex(sp), allocatable :: dsect(:,:,:)     !< (pftsz, kfrom:kto, 2) its derivatives in the two tilts
    real(dp)                 :: rref(3,3) = 0.d0 !< the drawn reference's rotation
    real(dp)                 :: omega(2)  = 0.d0 !< the tilts of the cached section
    integer                  :: state     = 0
    logical                  :: iseven    = .true.
    logical                  :: l_ref     = .false. !< a reference is set
    logical                  :: l_sect    = .false. !< sect and dsect hold the section at omega
    integer(int64)           :: nsections = 0_int64 !< sections cut since new
    integer(int64)           :: nsolves   = 0_int64 !< drawn references set since new (one per solve)
    integer(int64)           :: nevals    = 0_int64 !< cost evaluations since new
end type cont_sect_cache

type :: polarft_cont_calc
    private
    class(polarft_calc), pointer :: pftc => null()
    class(cartft_calc),  pointer :: cftc => null()
    integer                      :: pftsz      = 0
    integer                      :: kfromto(2) = 0
    real(sp),        allocatable :: q(:,:,:)     !< (2, pftsz, kfrom:kto) polar sample positions in the section frame
    type(cont_sect_cache), allocatable :: cache(:) !< one per thread of the team at new
    logical                      :: exists = .false.
contains
    procedure          :: new
    procedure          :: kill
    procedure          :: set_view_ref
    procedure, private :: section
    procedure          :: eval
    procedure          :: sigma_contrib
    procedure          :: get_counts
    procedure, nopass  :: view_frame
    procedure, private :: thread_cache
    procedure, private :: ensure_section
    procedure, private :: ensure_ref
end type polarft_cont_calc

interface
    !> white-box tests (submodule simple_polarft_cont_calc_tester, built only with BUILD_TESTS=ON)
    module subroutine run_all_polarft_cont_calc_tests()
    end subroutine run_all_polarft_cont_calc_tests
end interface

contains

    !> On the polar grid of pftc (its half ring and searched shells) and the references of cftc (made on pftc's box);
    !! one cache per thread of the team in use. Both must outlive this object, renewed with pftc; kills a previous one.
    subroutine new( self, pftc, cftc )
        class(polarft_cont_calc),    intent(inout) :: self
        class(polarft_calc), target, intent(inout) :: pftc
        class(cartft_calc),  target, intent(inout) :: cftc
        real(sp), allocatable :: tmp(:,:)
        integer :: nthr, ithr, k
        call self%kill
        self%pftc    => pftc
        self%cftc    => cftc
        self%pftsz   =  pftc%get_pftsz()
        self%kfromto =  pftc%get_kfromto()
        if( self%kfromto(2) > cftc%get_box()/2 ) THROW_HARD('polar shells beyond the Cartesian references')
        allocate(tmp(self%kfromto(1):self%kfromto(2),self%pftsz), self%q(2,self%pftsz,self%kfromto(1):self%kfromto(2)))
        call pftc%get_polar_coords(1, tmp)
        do k = self%kfromto(1), self%kfromto(2)
            self%q(1,:,k) = tmp(k,:)
        enddo
        call pftc%get_polar_coords(2, tmp)
        do k = self%kfromto(1), self%kfromto(2)
            self%q(2,:,k) = tmp(k,:)
        enddo
        deallocate(tmp)
        nthr = 1
        !$ nthr = omp_get_max_threads()
        allocate(self%cache(nthr))
        do ithr = 1, nthr
            allocate(self%cache(ithr)%sect(self%pftsz,self%kfromto(1):self%kfromto(2)), &
                &self%cache(ithr)%dsect(self%pftsz,self%kfromto(1):self%kfromto(2),2))
        enddo
        self%exists = .true.
    end subroutine new

    subroutine kill( self )
        class(polarft_cont_calc), intent(inout) :: self
        if( allocated(self%cache) ) deallocate(self%cache)
        if( allocated(self%q)     ) deallocate(self%q)
        self%pftc    => null()
        self%cftc    => null()
        self%pftsz   = 0
        self%kfromto = 0
        self%exists  = .false.
    end subroutine kill

    !> The drawn reference of the calling thread's solve (state, half, rotation); empties the thread's section cache
    subroutine set_view_ref( self, state, iseven, rref )
        class(polarft_cont_calc), intent(inout) :: self
        integer,                  intent(in)    :: state
        logical,                  intent(in)    :: iseven
        real(dp),                 intent(in)    :: rref(3,3)
        integer :: ithr
        ithr = self%thread_cache()
        self%cache(ithr)%state  = state
        self%cache(ithr)%iseven = iseven
        self%cache(ithr)%rref   = rref
        self%cache(ithr)%l_ref  = .true.
        self%cache(ithr)%l_sect = .false.
        self%cache(ithr)%nsolves = self%cache(ithr)%nsolves + 1_int64
    end subroutine set_view_ref

    !> The section and its normal derivative dV/dn = grad V . n, n = F^T z, at the frame F (the tester's point query,
    !! T1 and T2), in the native scale of a polar reference section and per native lattice unit
    subroutine section( self, state, iseven, frame, sect, dnv )
        class(polarft_cont_calc), intent(in)  :: self
        integer,                  intent(in)  :: state
        logical,                  intent(in)  :: iseven
        real(dp),                 intent(in)  :: frame(3,3)
        complex(sp),              intent(out) :: sect(self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(sp),              intent(out) :: dnv(self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(sp), allocatable :: dsect(:,:,:)
        real(dp) :: v(3,2)
        allocate(dsect(self%pftsz,self%kfromto(1):self%kfromto(2),2))
        v = 0.d0
        call cut(self%cftc, self%q, state, iseven, frame, v, sect, dsect, dnv)
    end subroutine section

    !> The cost of the calling thread's drawn reference at x = (shift x, shift y, rotation index, tilt x, tilt y): the raw
    !! Euclidean loss f, its gradient and, when curv is present, its Gauss-Newton matrix in the same order. A section is
    !! cut only when the tilts differ from the cached ones.
    subroutine eval( self, iptcl, x, f, grad, curv )
        class(polarft_cont_calc), intent(inout) :: self
        integer,                  intent(in)    :: iptcl
        real(dp),                 intent(in)    :: x(5)
        real(dp),                 intent(out)   :: f, grad(5)
        real(dp), optional,       intent(out)   :: curv(5,5)
        integer :: ithr
        ithr = self%thread_cache()
        call self%ensure_section(ithr, iptcl, x(4:5))
        self%cache(ithr)%nevals = self%cache(ithr)%nevals + 1_int64
        call self%pftc%gen_sect_euclid_grad_at_angle(iptcl, self%cache(ithr)%sect, self%cache(ithr)%dsect, x(1:2), x(3), &
            &f, grad, curv)
    end subroutine eval

    !> The per-shell residual of the reference (state, half, rotation rref) at x (the sigma2 contribution of a view
    !! solve's point): the section at the tilts x(4:5), shifted by x(1:2), at the rotation x(3) rounded to its grid
    !! index; the calling thread's cache takes the reference when another solve ran on the thread since
    subroutine sigma_contrib( self, iptcl, state, iseven, rref, x, contrib )
        class(polarft_cont_calc), intent(inout) :: self
        integer,                  intent(in)    :: iptcl, state
        logical,                  intent(in)    :: iseven
        real(dp),                 intent(in)    :: rref(3,3), x(5)
        real(sp),                 intent(out)   :: contrib(self%kfromto(1):self%kfromto(2))
        integer :: ithr, irot, nrots
        ithr  = self%thread_cache()
        call self%ensure_ref(ithr, state, iseven, rref)
        call self%ensure_section(ithr, iptcl, x(4:5))
        nrots = self%pftc%get_nrots()
        irot  = modulo(nint(x(3)) - 1, nrots) + 1
        call self%pftc%gen_sigma_contrib_sect(iptcl, self%cache(ithr)%sect, real(x(1:2),sp), irot, contrib)
    end subroutine sigma_contrib

    !> solves (drawn references set), sections cut and cost evaluations since new, over every thread
    subroutine get_counts( self, nsolves, nsections, nevals )
        class(polarft_cont_calc), intent(in)  :: self
        integer(int64),           intent(out) :: nsolves, nsections, nevals
        nsolves = 0_int64; nsections = 0_int64; nevals = 0_int64
        if( .not. allocated(self%cache) ) return
        nsolves   = sum(self%cache(:)%nsolves)
        nsections = sum(self%cache(:)%nsections)
        nevals    = sum(self%cache(:)%nevals)
    end subroutine get_counts

    !> the frame F = T(omega) R_ref at which the section of tilts omega is cut
    function view_frame( rref, omega ) result( frame )
        real(dp), intent(in) :: rref(3,3), omega(2)
        real(dp) :: frame(3,3), v(3,2)
        call view_frame_jacobian(rref, omega, frame, v)
    end function view_frame

    ! PRIVATE

    !> the calling thread's section at the tilts omega of its drawn reference, cut unless cached
    subroutine ensure_section( self, ithr, iptcl, omega )
        class(polarft_cont_calc), intent(inout) :: self
        integer,                  intent(in)    :: ithr, iptcl
        real(dp),                 intent(in)    :: omega(2)
        real(dp) :: frame(3,3), v(3,2)
        associate( c => self%cache(ithr) )
            if( .not. c%l_ref ) THROW_HARD('polarft_cont_calc evaluated without a drawn reference')
            if( self%pftc%get_pftsz() /= self%pftsz .or. any(self%pftc%get_kfromto() /= self%kfromto) )then
                THROW_HARD('polar grid changed after polarft_cont_calc%new')
            endif
            if( self%pftc%ptcl_iseven(iptcl) .neqv. c%iseven ) THROW_HARD('drawn reference and particle of different halves')
            if( .not. c%l_sect .or. any(omega /= c%omega) )then
                call view_frame_jacobian(c%rref, omega, frame, v)
                call cut(self%cftc, self%q, c%state, c%iseven, frame, v, c%sect, c%dsect)
                c%omega     = omega
                c%l_sect    = .true.
                c%nsections = c%nsections + 1_int64
            endif
        end associate
    end subroutine ensure_section

    !> the calling thread's cache on the reference (state, half, rotation) unless it holds it already; no solve counted
    subroutine ensure_ref( self, ithr, state, iseven, rref )
        class(polarft_cont_calc), intent(inout) :: self
        integer,                  intent(in)    :: ithr, state
        logical,                  intent(in)    :: iseven
        real(dp),                 intent(in)    :: rref(3,3)
        associate( c => self%cache(ithr) )
            if( c%l_ref .and. c%state == state .and. (c%iseven .eqv. iseven) )then
                if( all(c%rref == rref) ) return
            endif
            c%state  = state
            c%iseven = iseven
            c%rref   = rref
            c%l_ref  = .true.
            c%l_sect = .false.
        end associate
    end subroutine ensure_ref

    integer function thread_cache( self ) result( ithr )
        class(polarft_cont_calc), intent(in) :: self
        if( .not. self%exists ) THROW_HARD('polarft_cont_calc used before new')
        ithr = 1
        !$ ithr = omp_get_thread_num() + 1
        if( ithr > size(self%cache) ) THROW_HARD('polarft_cont_calc thread beyond the team at new')
    end function thread_cache

    !> The section at the frame F and its derivatives along the rotation directions v(:,a) (left increments F ->
    !! exp([v_a]x) F, in the section frame): with q the polar sample in the section frame and g = F grad V its gradient
    !! in that frame, dS = -g . (v_a x q); dnv (optional) is g(3)
    subroutine cut( cftc, q, state, iseven, frame, v, sect, dsect, dnv )
        class(cartft_calc),    intent(in)  :: cftc
        real(sp),              intent(in)  :: q(:,:,:)
        integer,               intent(in)  :: state
        logical,               intent(in)  :: iseven
        real(dp),              intent(in)  :: frame(3,3), v(3,2)
        complex(sp),           intent(out) :: sect(:,:), dsect(:,:,:)
        complex(sp), optional, intent(out) :: dnv(:,:)
        complex(sp) :: val, grad(3), gs(3)
        real(sp)    :: f(3,3), vs(3,2), loc(3), q1, q2
        integer     :: k, irot, a
        f  = real(frame, sp)
        vs = real(v, sp)
        ! q, sect, dsect and dnv share the shape (:, pftsz, nk) of the searched shells, from 1
        do k = 1, size(q,3)
            do irot = 1, size(q,2)
                q1  = q(1,irot,k)
                q2  = q(2,irot,k)
                loc = q1*f(1,:) + q2*f(2,:)
                call cftc%gather_ref(state, iseven, loc, val, grad)
                gs  = f(:,1)*grad(1) + f(:,2)*grad(2) + f(:,3)*grad(3)
                sect(irot,k) = val
                do a = 1, 2
                    dsect(irot,k,a) = -(vs(3,a)*(q1*gs(2) - q2*gs(1)) + (vs(1,a)*q2 - vs(2,a)*q1)*gs(3))
                enddo
                if( present(dnv) ) dnv(irot,k) = gs(3)
            enddo
        enddo
    end subroutine cut

    !> F = exp([w]x) R_ref for w = (omega, 0) and the left-Jacobian columns v(:,a) = J_l(w) e_a: a change d of omega
    !! moves the frame by exp([J_l d]x) to first order, J_l = I + (1 - cos t)/t^2 [w]x + (t - sin t)/t^3 [w]x^2
    pure subroutine view_frame_jacobian( rref, omega, frame, v )
        real(dp), intent(in)  :: rref(3,3), omega(2)
        real(dp), intent(out) :: frame(3,3), v(3,2)
        real(dp) :: w(3), wx(3,3), wx2(3,3), t, t2, s1, c1, a1, e(3,3), jl(3,3)
        integer  :: i
        w  = [omega(1), omega(2), 0.d0]
        wx = reshape([0.d0, w(3), -w(2), -w(3), 0.d0, w(1), w(2), -w(1), 0.d0], [3,3])
        wx2 = matmul(wx, wx)
        t2 = sum(w*w)
        t  = sqrt(t2)
        if( t < 1.d-4 )then
            s1 = 1.d0 - t2/6.d0              ! sin t/t
            c1 = 0.5d0 - t2/24.d0            ! (1 - cos t)/t^2
            a1 = 1.d0/6.d0 - t2/120.d0       ! (t - sin t)/t^3
        else
            s1 = sin(t)/t
            c1 = (1.d0 - cos(t))/t2
            a1 = (t - sin(t))/(t2*t)
        endif
        e  = s1*wx + c1*wx2
        jl = c1*wx + a1*wx2
        do i = 1, 3
            e(i,i)  = e(i,i)  + 1.d0
            jl(i,i) = jl(i,i) + 1.d0
        enddo
        frame = matmul(e, rref)
        v     = jl(:,1:2)
    end subroutine view_frame_jacobian

end module simple_polarft_cont_calc
