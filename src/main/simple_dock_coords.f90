! Cartesian coordinates-coordiantes alignment based on rmsd
module simple_dock_coords
include 'simple_lib.f08'
use simple_atoms,       only: atoms
use simple_opt_spec,    only: opt_spec
use simple_opt_simplex, only: opt_simplex
use simple_ori
implicit none

public :: dock_coords_init, dock_coords_minimize
private

logical, parameter :: INI_W_DISCR_SRCH = .false.

type(opt_spec)        :: ospec                             !< optimizer specification object
type(opt_simplex)     :: nlopt                             !< optimizer object
integer               :: nrestarts = 3                     !< simplex restarts (randomized bounds)
real                  :: lims(6,2)                         !< variable bounds
type(atoms)           :: aref                              !< reference atoms
type(atoms)           :: atarg                             !< target atoms (subjected to shift and rotation)
real, allocatable     :: ref_coords(:,:), targ_coords(:,:) !< set of coordinates of the referece/target
integer               :: nref, ntarg                       !< number of points in the reference/target

contains

    subroutine dock_coords_init( a_ref, a_targ, trs_in, nrestarts_in )
        class(atoms),      intent(inout) :: a_ref, a_targ
        real,              intent(in)    :: trs_in ! for shift search
        integer, optional, intent(in)    :: nrestarts_in
        integer :: i
        real    :: r_com(3), t_com(3)
        nrestarts = 3
        if( present(nrestarts_in) ) nrestarts = nrestarts_in
        call aref%copy(a_ref)
        call atarg%copy(a_targ)
        nref  = aref%get_n()
        ntarg = atarg%get_n()
        ! calc center of mass of the ref and translate both the ref and target by the opposite of that.
        r_com = aref%get_geom_center()
        ! t_com = atarg%get_geom_center()
        call aref%translate(-r_com)
        call atarg%translate(-r_com)
        ! same center of mass
        ! call atarg%translate()
        ! fill in the coordinates
        allocate(ref_coords(3,nref), targ_coords(3,ntarg), source = 0.)
        do i = 1, nref
            ref_coords(:,i) = aref%get_coord(i)
        enddo
        do i = 1, ntarg
            targ_coords(:,i) = atarg%get_coord(i)
        enddo

        lims = 0.
        lims(1,2)   = 359.99
        lims(2,2)   = 180.
        lims(3,2)   = 359.99      ! lims for eul angle
        lims(4:6,1) =  - trs_in   ! lims for translation
        lims(4:6,2) =    trs_in
        print *, 'lims: '
        call vis_mat(lims)

        call ospec%specify('de', 6, ftol=1e-4,&
            &gtol=1e-4, limits=lims, nrestarts=nrestarts, maxits=200) ! 6 is thedimensionality
        call ospec%set_costfun(dock_coords_srch_costfun)
        call nlopt%new(ospec)
    end subroutine dock_coords_init

    ! The cost is atoms rmsd when rotating second set of points with angles angles
    ! and translating it with the vector vec
    function dock_coords_srch_costfun( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)    !angles, shift
        real                    :: cost, rmat(3,3)
        rmat = euler2m(vec(1:3))
        cost = rmsd() ! rmat > rotation, vec(2:6) > shift
        print *, 'Rotation: ', vec(1:3), 'Translation: ', vec(4:6), 'Cost: ', cost
    contains
        !     ! TO COMPLETE
        function rmsd() result(cost)
            real                    :: cost
            type(atoms) :: arott
            integer :: nmin, nmax, location(1)
            integer :: i, cnt1
            real,    allocatable :: dist(:), dist_sq(:)
            logical, allocatable :: mask(:)
            call arott%copy(atarg)
            call arott%translate(-vec(4:6)) ! - for convention
            call arott%rotate(rmat)
            ! set coords targ_coords with the translated ones
            do i = 1, ntarg
                targ_coords(:,i) = arott%get_coord(i)
            enddo
            call arott%kill
            ! Rotate and translate a_targ in a copy
            if(nref <= ntarg) then
                nmin = nref
                nmax = ntarg
                allocate(dist(nmax), dist_sq(nmax), source = 0.)
                allocate(mask(nmin), source = .true.)
                cnt1 = 0
                do i = 1, nmax !compare based on target             ! TO FIX FOR WHEN THEY HAVE DIFFERENT NB OF ATOMS
                    if(cnt1+1 <= nmin) then ! just nmin couples, starting from 0
                        dist(i) = pixels_dist(targ_coords(:,i),ref_coords(:,:),'min',mask,location, keep_zero=.true.)
                        cnt1 = cnt1 + 1 ! number of couples
                        mask(location(1)) = .false. ! not to consider the same atom more than once
                        dist_sq(i) = dist(i)**2 !formula wants them square, could improve performance here
                    endif
                enddo
            else
                nmin = ntarg
                nmax = nref
                allocate(dist(nmax), dist_sq(nmax), source = 0.)
                allocate(mask(nmin), source = .true.)
                cnt1 = 0
                do i = 1, nmax !compare based on reference
                    if(cnt1+1 <= nmin) then ! just nmin couples, starting from 0
                        dist(i) = pixels_dist(ref_coords(:,i),targ_coords(:,:),'min',mask,location, keep_zero = .true.)
                        cnt1 = cnt1 + 1
                        mask(location(1)) = .false. ! not to consider the same atom more than once
                        dist_sq(i) = dist(i)**2     ! formula wants them square
                    endif
                enddo
            endif
            !RMSD
            cost = sqrt(sum(dist_sq)/real(nmin))
        end function rmsd
    end function dock_coords_srch_costfun

    function dock_coords_minimize( ) result( rot_trans )
        real :: cost_init, cost, rot_trans(7), cost_best
        real :: shvec(3), shvec_best(3),ang(3), ang_best(3), xsh, ysh, zsh, xang, yang, zang, vec(6)
        class(*), pointer :: fun_self => null()
        shvec_best = 0.
        ang_best   = 0.
        ! if( INI_W_DISCR_SRCH )then
        !     ! discrete search to start-off with (half-pixel resolution for shift, 1 degree for angle)
        !     cost_best = 100. ! to be sure to go into the loop
        !     xang = lims(1,1)
        !     do while(xang <= lims(1,2))
        !         yang = lims(2,1)
        !         do while(yang <= lims(2,2))
        !             zang = lims(3,1)
        !             do while(zang <= lims(3,2))
        !                 xsh = lims(4,1)
        !                 do while( xsh <= lims(4,2) )
        !                     ysh = lims(5,1)
        !                     do while( ysh <= lims(5,2) )
        !                         zsh = lims(6,1)
        !                         do while( zsh <= lims(6,2) )
        !                             shvec = [xsh,ysh,zsh]
        !                             ang   = [xang,yang,zang]
        !                             vec(1:3) = ang
        !                             vec(4:6) = shvec
        !                             cost  = dock_coords_srch_costfun(fun_self, vec, 6)
        !                             if( cost < cost_best )then
        !                                 shvec_best = shvec
        !                                 ang_best   = ang
        !                                 cost_best  = cost
        !                             endif
        !                             zsh = zsh + 1.
        !                         end do
        !                         ysh = ysh + 1.
        !                     end do
        !                     xsh = xsh + 1.
        !                 end do
        !                 zang = zang + 2. ! stepsz two degrees
        !             enddo
        !             yang = yang + 2.
        !         enddo
        !         xang = xang + 2.
        !     enddo
        ! endif
        ! refinement with simplex
        ospec%x(1:3)  = ang_best
        ospec%x(4:6)  = shvec_best ! assumed that vol is shifted to previous centre
        cost_init = dock_coords_srch_costfun(fun_self, ospec%x, 6)
        print *, 'cost_init', cost_init
        call nlopt%minimize(ospec, fun_self, cost)
        print *, 'cost, cost_init: ', cost, cost_init
        if( cost < cost_init )then
            rot_trans(1)   = cost    ! rmsd
            rot_trans(2:7) = ospec%x ! angles, shift
        else
            rot_trans(1)   = cost_init
            rot_trans(2:7) =  0.
        endif
    end function dock_coords_minimize

end module simple_dock_coords
