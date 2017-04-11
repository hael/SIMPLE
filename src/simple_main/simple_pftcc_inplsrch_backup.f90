module simple_pftcc_inplsrch
use simple_opt_factory,      only: opt_factory
use simple_opt_spec,         only: opt_spec
use simple_optimizer,        only: optimizer
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_defs
implicit none

public :: pftcc_inplsrch_init, pftcc_inplsrch_set_indices, pftcc_inplsrch_minimize, pftcc_inplsrch_get_nevals
private

type(opt_factory)                :: ofac                   !< optimizer factory
type(opt_spec)                   :: ospec                  !< optimizer specification object
class(optimizer), pointer        :: nlopt       => null()  !< pointer to nonlinear optimizer
class(polarft_corrcalc), pointer :: pftcc_ptr   => null()  !< pointer to pftcc object
integer                          :: reference   =  0       !< reference pft
integer                          :: particle    =  0       !< particle pft
integer                          :: ldim(3)     =  [0,0,0] !< logical dimension of Cartesian image
real                             :: rotmat(2,2) =  0.      !< rotation matrix for checking limits
real                             :: maxshift    =  0.      !< maximal shift
logical                          :: shbarr      =  .true.  !< shift barrier constraint or not
integer, parameter               :: NRESTARTS_DEFAULT =  5 !< simplex restarts (randomized bounds)

contains

    subroutine pftcc_inplsrch_init( pftcc, shlims, shbarrier, nrestarts )
        class(polarft_corrcalc), target, intent(in) :: pftcc
        real,                            intent(in) :: shlims(2,2)
        character(len=*), optional,      intent(in) :: shbarrier
        integer,          optional,      intent(in) :: nrestarts
        real :: lims(3,2)
        integer :: nnrestarts
        ! flag the barrier constraint
        shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) shbarr = .false.
        endif
        nnrestarts = NRESTARTS_DEFAULT
        if( present(nrestarts) ) nnrestarts = nrestarts 
        ! make optimizer spec
        lims(1,1)  = 0.
        lims(1,2)  = 360.
        lims(2:,:) = shlims
        call ospec%specify('simplex', 3, ftol=1e-4, gtol=1e-4, limits=lims, nrestarts=nnrestarts)
        ! set optimizer cost function
        call ospec%set_costfun(pftcc_inplsrch_cost)
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
    end subroutine pftcc_inplsrch_init
    
    subroutine pftcc_inplsrch_set_indices( ref, ptcl )
        integer, intent(in) :: ref, ptcl
        reference = ref 
        particle  = ptcl
    end subroutine pftcc_inplsrch_set_indices

    function pftcc_inplsrch_cost( vec, D ) result( cost )
        use simple_math, only: rotmat2d, enforce_cyclic_limit
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        real    :: vec_here(3)    ! current set of values
        real    :: rotvec_here(2) ! current set of shift values rotated to frame of reference
        real    :: cost
        integer :: rot
        vec_here = vec
        ! check so that the in-plane rotation is within the limit
        call enforce_cyclic_limit(vec_here(1), 360.)
        ! zero small shifts
        if( abs(vec(2)) < 1e-6 ) vec_here(2) = 0.
        if( abs(vec(3)) < 1e-6 ) vec_here(3) = 0.
        ! set rotmat for boudary checking and final rotation
        rotmat      = rotmat2d( vec_here(1) )
        rotvec_here = matmul(vec_here(2:3),rotmat)
        ! check boundary 
        if( shbarr )then
            if( rotvec_here(1) < ospec%limits(2,1) .or. rotvec_here(1) > ospec%limits(2,2) )then
                cost = 1.
                return
            else if( rotvec_here(2) < ospec%limits(3,1) .or. rotvec_here(2) > ospec%limits(3,2) )then
                cost = 1.
                return
            endif
        endif
        rot  =  pftcc_ptr%get_roind(vec_here(1))
        cost = -pftcc_ptr%corr(reference, particle, rot, vec_here(2:3))
    end function pftcc_inplsrch_cost
    
    function pftcc_inplsrch_get_nevals() result( nevals )
        integer :: nevals
        nevals = ospec%nevals
    end function pftcc_inplsrch_get_nevals
    
    function pftcc_inplsrch_minimize( irot, shvec, rxy ) result( crxy )
        use simple_math, only: rotmat2d, enforce_cyclic_limit
        use simple_rnd,  only: ran3
        integer, intent(in), optional :: irot
        real,    intent(in), optional :: shvec(2)
        real,    intent(in), optional :: rxy(3)
        logical :: irot_here, shvec_here, rxy_here
        real    :: crxy(4)
        irot_here  = .false.
        shvec_here = .false.
        rxy_here   = .false.
        if( present(irot)  ) irot_here  = .true.
        if( present(shvec) ) shvec_here = .true.
        if( present(rxy)   ) rxy_here   = .true.
        if( irot_here  .and. rxy_here ) stop 'irot & rxy cannot be simultaneously present; pftcc_inplsrch_minimize'
        if( shvec_here .and. rxy_here ) stop 'shvec_here & rxy_here cannot be simultaneously present; pftcc_inplsrch_minimize'
        ! initialisation
        ospec%x      = 0.
        ospec%x(1)   = ran3()*360.
        ospec%nevals = 0
        if( rxy_here )then
            ospec%x(1) = rxy(1)
            ospec%x(2) = rxy(2)
            ospec%x(3) = rxy(3)
        endif
        if( irot_here ) ospec%x(1) = pftcc_ptr%get_rot(irot)
        if( shvec_here )then
            ospec%x(1) = shvec(1)
            ospec%x(2) = shvec(2)
        endif
        ! minimisation
        call nlopt%minimize(ospec, crxy(1))
        crxy(1)  = -crxy(1)   ! correlation
        crxy(2)  = ospec%x(1) ! in-plane rotation
        ! check so that the in-plane rotation is within the limit
        call enforce_cyclic_limit(crxy(2), 360.)
        ! set rotmat for boudary checking and final rotation
        rotmat   = rotmat2d( crxy(2) )
        ! rotate the shift vector to the frame of reference
        crxy(3:) = ospec%x(2:) ! shift
        crxy(3:) = matmul(crxy(3:), rotmat)
        if( any(crxy(3:) > maxshift) .or. any(crxy(3:) < -maxshift) )then
            crxy(1)  = -1.
            crxy(3:) = 0.
        endif
        call init_rotmat ! clean exit
    end function pftcc_inplsrch_minimize

    subroutine init_rotmat
        ! identity matrix
        rotmat      = 0.
        rotmat(1,1) = 1.
        rotmat(2,2) = 1.
    end subroutine init_rotmat
    
end module simple_pftcc_inplsrch
