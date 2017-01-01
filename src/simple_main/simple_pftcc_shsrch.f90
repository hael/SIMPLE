module simple_pftcc_shsrch
use simple_opt_factory,      only: opt_factory
use simple_opt_spec,         only: opt_spec
use simple_optimizer,        only: optimizer
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_defs
implicit none

public :: pftcc_shsrch_init, pftcc_shsrch_set_indices, pftcc_shsrch_minimize, pftcc_shsrch_get_nevals,&
pftcc_shsrch_cost
private

type(opt_factory)                :: ofac                 !< optimizer factory
type(opt_spec)                   :: ospec                !< optimizer specification object
class(optimizer), pointer        :: nlopt      =>null()  !< pointer to nonlinear optimizer
class(polarft_corrcalc), pointer :: pftcc_ptr  =>null()  !< pointer to pftcc object
integer                          :: reference  = 0       !< reference pft
integer                          :: particle   = 0       !< particle pft
integer                          :: rot        = 1       !< in-plane rotation
integer                          :: ldim(3)    = [0,0,0] !< logical dimension of Cartesian image
real                             :: rotmat(2,2)= 0.      !< rotation matrix for checking limits
real                             :: maxshift   = 0.      !< maximal shift
logical                          :: shbarr     = .true.  !< shift barrier constraint or not
integer, parameter               :: NRESTARTS  =  5      !< simplex restarts (randomized bounds)


contains

    subroutine pftcc_shsrch_init( pftcc, lims, shbarrier )
        class(polarft_corrcalc), intent(in), target :: pftcc
        real,             intent(in)                :: lims(2,2)
        character(len=*), intent(in), optional      :: shbarrier
        integer :: nnrestarts
        ! flag the barrier constraint
        shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) shbarr = .false.
        endif
        ! make optimizer spec
        call ospec%specify('simplex', 2, ftol=1e-4, gtol=1e-4, limits=lims, nrestarts=NRESTARTS)
        ! set optimizer cost function
        call ospec%set_costfun(pftcc_shsrch_cost)
        ! generate optimizer object with the factory
        call ofac%new(ospec, nlopt)
        ! set pointer to corrcalc object
        pftcc_ptr => pftcc
        ! get logical dimension
        ldim = pftcc_ptr%get_ldim()
        ! set maxshift
        maxshift = real(maxval(ldim))/2.
        ! rotmat init
        call init_rotmat
    end subroutine pftcc_shsrch_init
    
    subroutine pftcc_shsrch_set_indices( ref, ptcl, r )
        integer, intent(in) :: ref, ptcl, r
        reference = ref 
        particle  = ptcl
        rot       = r
    end subroutine pftcc_shsrch_set_indices

    function pftcc_shsrch_cost( vec, D ) result( cost )
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        real    :: vec_here(2)    ! current set of values
        real    :: rotvec_here(2) ! current set of values rotated to frame of reference
        real    :: cost
        vec_here = vec
        if( abs(vec(1)) < 1e-6 ) vec_here(1) = 0.
        if( abs(vec(2)) < 1e-6 ) vec_here(2) = 0.
        rotvec_here = matmul(vec_here,rotmat)
        if( shbarr )then
            if( rotvec_here(1) < ospec%limits(1,1) .or. rotvec_here(1) > ospec%limits(1,2) )then
                cost = 1.
                return
            else if( rotvec_here(2) < ospec%limits(2,1) .or. rotvec_here(2) > ospec%limits(2,2) )then
                cost = 1.
                return
            endif
        endif
        cost = -pftcc_ptr%corr(reference, particle, rot, vec_here)
    end function pftcc_shsrch_cost
    
    function pftcc_shsrch_get_nevals() result( nevals )
        integer :: nevals
        nevals = ospec%nevals
    end function pftcc_shsrch_get_nevals
    
    function pftcc_shsrch_minimize( ) result( cxy )
        use simple_math, only: rotmat2d
        real :: cxy(3), cost, cost_init
        ! set rotmat for boudary checking and final rotation
        rotmat = rotmat2d( pftcc_ptr%get_rot(rot) )
        ! minimisation
        ospec%x = 0.
        ospec%nevals = 0
        cost_init = pftcc_shsrch_cost( ospec%x, ospec%ndim )
        call nlopt%minimize(ospec, cost)
        if( cost < cost_init )then
            cxy(1)  = -cost ! correlation
            ! rotate the shift vector to the frame of reference
            cxy(2:) = ospec%x ! shift
            cxy(2:) = matmul(cxy(2:),rotmat)
            if( any(cxy(2:) > maxshift) .or. any(cxy(2:) < -maxshift) )then
                cxy(1)  = -1.
                cxy(2:) = 0.
            endif
        else
             cxy(1)  = -cost_init ! correlation
             cxy(2:) = 0.
        endif
        call init_rotmat ! clean exit
    end function pftcc_shsrch_minimize

    subroutine init_rotmat()
        ! identity matrix
        rotmat      = 0.
        rotmat(1,1) = 1.
        rotmat(2,2) = 1.
    end subroutine init_rotmat
    
end module simple_pftcc_shsrch
