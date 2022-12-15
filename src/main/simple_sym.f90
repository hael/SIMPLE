! defines protein point-group symmetries
module simple_sym
! include 'simple_lib.f08'
use simple_ori
use simple_oris
use simple_math
use simple_rnd
implicit none

public :: sym, sym_tester
private
#include "simple_local_flags.inc"

type sym
    private
    type(oris)                    :: e_sym                    !< symmetry eulers
    character(len=3), allocatable :: subgrps(:)               !< subgroups
    real                          :: eullims(3,2) = 0.        !< euler angles limits (asymetric unit, including mirrors)
    real                          :: eullims_nomirr(3,2) = 0. !< euler angles limits (asymetric unit, excluding mirrors)
    integer                       :: n=1                      !< # symmetry ops
    integer                       :: ncsym=1                  !< num of order C sym
    integer                       :: ndsym=1                  !< num of order D sym
    integer                       :: t_or_o=0                 !< tetahedral or octahedral symmetry
    character(len=3)              :: pgrp='c1'                !< point-group symmetry
    logical                       :: c_or_d=.false.           !< c- or d-symmetry
  contains
    ! constructor
    procedure          :: new
    ! getters
    procedure          :: get_eullims
    procedure          :: get_nsym
    procedure          :: get_pgrp
    procedure          :: get_symori
    procedure          :: get_nsubgrp
    procedure          :: get_subgrp
    procedure          :: has_subgrp
    procedure          :: get_subgrp_descr
    procedure          :: get_all_subgrps_descr
    ! setters
    procedure, private :: set_subgrps
    ! modifiers
    procedure, private :: apply_1, apply_2
    generic            :: apply => apply_1, apply_2
    procedure, private :: apply2all
    procedure          :: apply_sym_with_shift
    procedure          :: find_closest_proj
    procedure          :: nearest_proj_neighbors
    procedure          :: rot_to_asym
    procedure          :: rotall_to_asym
    procedure          :: symrandomize
    procedure          :: build_refspiral
    procedure, private :: rnd_euler_1
    procedure, private :: rnd_euler_2
    procedure, private :: rnd_euler_3
    generic            :: rnd_euler => rnd_euler_1, rnd_euler_2, rnd_euler_3
    ! calculators
    procedure, private :: sym_dists_1, sym_dists_2
    generic            :: sym_dists => sym_dists_1, sym_dists_2
    procedure          :: nearest_sym_neighbors
    ! utils
    procedure, private :: build_eullims
    procedure, private :: make_c_and_d
    procedure, private :: make_t
    procedure, private :: make_o
    procedure, private :: make_i_spider
    procedure, private :: make_i_relion
    procedure, private :: within_asymunit
    procedure, private :: within_platonic_asymunit
    ! destructor
    procedure          :: kill
end type sym

interface sym
    module procedure constructor
end interface sym

integer, parameter          :: ntet   = 12 ! # tetahedral symmetry operations
integer, parameter          :: noct   = 24 ! # octahedral symmetry operations
integer, parameter          :: nico   = 60 ! # icosahedral symmetry operations
double precision, parameter :: delta2 = 180.d0
double precision, parameter :: delta3 = 120.d0
double precision, parameter :: delta5 = 72.d0
double precision, parameter :: deltan = 36.d0
double precision, parameter :: gamma  = 54.735610317245346d0
real,             parameter :: tet_theta_lim = 70.52877936550931
real,             parameter :: oct_theta_lim = real(gamma)
real,             parameter :: ico_theta_lim = 37.37736814064969

contains

    !>  \brief  is a constructor
    function constructor( pgrp, icorelion ) result( self )
        character(len=*), intent(in)  :: pgrp        !< sym group string
        logical, optional, intent(in) :: icorelion
        type(sym)                     :: self
        call self%new(pgrp, icorelion)
    end function constructor

    !>  \brief  is a constructor
    subroutine new( self, pgrp, icorelion )
        class(sym),        intent(inout) :: self
        character(len=*),  intent(in)    :: pgrp         !< sym group string
        logical, optional, intent(in)    :: icorelion
        call self%kill
        self%c_or_d = .false.
        self%n      = 1
        self%ncsym  = 1
        self%ndsym  = 1
        self%pgrp   = pgrp
        self%t_or_o = 0
        if(pgrp(1:1).eq.'c' .or. pgrp(1:1).eq.'C')then
            if( self%pgrp(1:1).eq.'C' ) self%pgrp(1:1) = 'c'
            self%c_or_d = .true.
            read(pgrp(2:),'(I2)') self%ncsym
            self%n = self%ncsym
            call self%e_sym%new(self%n, is_ptcl=.false.)
            call self%make_c_and_d
        else if(pgrp(1:1).eq.'d' .or. pgrp(1:1).eq.'D')then
            if( self%pgrp(1:1).eq.'D' ) self%pgrp(1:1) = 'd'
            self%c_or_d = .true.
            self%ndsym  = 2
            read(pgrp(2:),'(I2)') self%ncsym
            self%n = self%ncsym*self%ndsym
            call self%e_sym%new(self%n, is_ptcl=.false.)
            call self%make_c_and_d
        else if( pgrp(1:1).eq.'t' .or. pgrp(1:1).eq.'T' )then
            if( self%pgrp(1:1).eq.'T' ) self%pgrp(1:1) = 't'
            self%t_or_o = 1
            self%n      = ntet
            call self%e_sym%new(self%n, is_ptcl=.false.)
            call self%make_t
        else if( pgrp(1:1).eq.'o' .or. pgrp(1:1).eq.'O' )then
            if( self%pgrp(1:1).eq.'O' ) self%pgrp(1:1) = 'o'
            self%t_or_o = 3
            self%n      = noct
            call self%e_sym%new(self%n, is_ptcl=.false.)
            call self%make_o
        else if( pgrp(1:1).eq.'i' .or. pgrp(1:1).eq.'I' )then
            if( self%pgrp(1:1).eq.'I' ) self%pgrp(1:1) = 'i'
            self%n = nico
            call self%e_sym%new(self%n, is_ptcl=.false.)
            if( present(icorelion) )then
                if( icorelion )then
                    call self%make_i_relion
                else
                    call self%make_i_spider
                endif
            else
                call self%make_i_spider
            endif
        else
            write(logfhandle,'(a)') 'symmetry not supported; new; simple_sym', pgrp
            stop
        endif
        if(trim(self%pgrp).ne.'d1')call self%set_subgrps
        call self%build_eullims
    end subroutine new

    !>  \brief  builds the search range for the point-group
    subroutine build_eullims( self )
        class(sym), intent(inout) :: self
        self%eullims(:,1) = 0.
        self%eullims(:,2) = 360.
        self%eullims(2,2) = 180.
        self%eullims_nomirr = self%eullims
        if( self%pgrp(1:1).eq.'c' )then
            if( self%n > 1 )then
                self%eullims(1,2)  = 360./real(self%ncsym)
            endif
            self%eullims_nomirr      = self%eullims
            self%eullims_nomirr(2,2) = 90.
        else if( self%pgrp(1:1).eq.'d' .and. self%ncsym > 1 )then
            self%eullims(1,2)      = 360./real(self%ncsym)
            self%eullims(2,2)      = 90.
            self%eullims_nomirr      = self%eullims
            if(is_even(self%n/2))then
                self%eullims_nomirr(1,2) = 360./(2.*real(self%ncsym))
            else
                self%eullims_nomirr(1,1) = 360./(4.*real(self%ncsym))
                self%eullims_nomirr(1,2) = 360.*3./(4.*real(self%ncsym))
            endif
        else if( self%pgrp(1:1).eq.'t' )then
            self%eullims(1,2)   = 120.
            self%eullims(2,2)   = tet_theta_lim
            self%eullims_nomirr = self%eullims
        else if( self%pgrp(1:1).eq.'o' )then
            self%eullims(1,2)        = 90.
            self%eullims(2,2)        = oct_theta_lim
            self%eullims_nomirr      = self%eullims
            self%eullims_nomirr(1,2) = 45.
        else if( self%pgrp(1:1).eq.'i' )then
            self%eullims(1,2)        = real(delta5)
            self%eullims(2,2)        = ico_theta_lim
            self%eullims_nomirr      = self%eullims
            self%eullims_nomirr(1,2) = 36.
        endif
    end subroutine build_eullims

    !>  \brief  returns the search range for the point-group
    function get_eullims( self ) result( eullims )
        class(sym), intent(inout) :: self !< this instance
        real                      :: eullims(3,2) !< 3-axis search range
        eullims = self%eullims
    end function get_eullims

    !>  \brief  is a getter
    pure function get_nsym( self ) result( n )
        class(sym), intent(in) :: self !< this instance
        integer :: n
        n = self%n
    end function get_nsym

    !>  \brief  is a getter
    pure function get_pgrp( self ) result( pgrp_str )
        class(sym), intent(in) :: self !< this instance
        character(len=3) :: pgrp_str  !< sym group string
        pgrp_str = self%pgrp
    end function get_pgrp

    !>  \brief  is a getter
    function get_nsubgrp( self )result( n )
        class(sym) :: self
        integer :: n      !< size of sym groups
        n = size(self%subgrps)
    end function get_nsubgrp

    !>  \brief  is a getter
    function get_subgrp( self, i )result( symobj )
        class(sym),intent(in) :: self
        type(sym) :: symobj
        integer   :: i, n
        n = size(self%subgrps)
        if( (i>n).or.(i<1) )then
            write(logfhandle,*)'Index out of bonds on simple_sym; get_subgroup'
            stop
        endif
        symobj = sym(self%subgrps(i))
    end function get_subgrp

    function has_subgrp( self, subgrp ) result( has )
        class(sym),       intent(in) :: self
        character(len=*), intent(in) :: subgrp
        integer :: i, n
        logical :: has
        has = .false.
        n   = size(self%subgrps)
        if( n < 1 ) return
        do i=1,n
            if( trim(self%subgrps(i)) .eq. trim(subgrp) )then
                has = .true.
                return
            endif
        end do
    end function has_subgrp

    !>  \brief  is a getter
    function get_subgrp_descr( self, i )result( str_descr )
        class(sym),intent(in) :: self
        character(len=3)      :: str_descr
        integer   :: i, n
        n = size(self%subgrps)
        if( (i>n).or.(i<1) )then
            write(logfhandle,*)'Index out of bonds on simple_sym; get_subgrp_descr'
            stop
        endif
        str_descr = self%subgrps(i)
    end function get_subgrp_descr

    function get_all_subgrps_descr( self )result( str_descr )
        class(sym),intent(in) :: self
        character(len=:), allocatable :: str_descr, str_tmp
        integer :: i, n
        n = size(self%subgrps)
        allocate(str_descr, source=trim(self%subgrps(1)))
        if( n == 1 ) return
        do i=2,n
            str_tmp   = str_descr
            str_descr = str_tmp//' '//trim(self%subgrps(i))
            deallocate(str_tmp)
        end do
    end function get_all_subgrps_descr

    !>  \brief  is a symmetry adaptor
    subroutine apply_1( self, e_in, symop, e_sym )
        class(sym), intent(in)    :: self
        class(ori), intent(in)    :: e_in  !< orientation object
        integer,    intent(in)    :: symop !< index of sym operation
        type(ori),  intent(inout) :: e_sym
        type(ori)                 :: e_symop, e_tmp
        e_sym = e_in ! transfer of parameters
        call self%e_sym%get_ori(symop, e_symop)
        call e_symop%compose(e_in, e_tmp)
        call e_sym%set_euler(e_tmp%get_euler())
        call e_symop%kill
        call e_tmp%kill
    end subroutine apply_1

    !>  \brief  is a symmetry adaptor, operates on euler triplet
    subroutine apply_2( self, euls_in, symop, euls_sym )
        class(sym), intent(in)    :: self
        real,       intent(in)    :: euls_in(3)
        integer,    intent(in)    :: symop !< index of sym operation
        real,       intent(inout) :: euls_sym(3)
        real :: euls_symop(3)
        euls_symop = self%e_sym%get_euler(symop)
        call euler_compose(euls_symop, euls_in, euls_sym)
    end subroutine apply_2

    !>  \brief  is a symmetry adaptor
    subroutine apply2all( self, e_in, symop )
        class(sym), intent(inout) :: self
        class(oris), intent(inout):: e_in  !< orientations object
        integer,    intent(in)    :: symop !< index of sym operation
        type(ori)                 :: o, o2
        integer :: i
        do i=1,e_in%get_noris()
            call e_in%get_ori(i, o)
            call self%apply(o, symop, o2)
            call e_in%set_ori(i, o2)
        enddo
        call o%kill
        call o2%kill
    end subroutine apply2all

    ! !>  \brief  rotates any orientation to the asymmetric unit
    subroutine rot_to_asym( self, osym )
        class(sym), intent(inout) :: self
        class(ori), intent(inout) :: osym !< orientation
        type(ori) :: oasym
        integer   :: nsym
        if( self%within_asymunit(osym, incl_mirr=.true.) )then
            ! already in asymetric unit
        else
            do nsym=2,self%n ! nsym=1 is the identity operator
                call self%apply(osym, nsym, oasym)
                if( self%within_asymunit(oasym, incl_mirr=.true.) ) exit
            enddo
            osym = oasym
        endif
        call oasym%kill
    end subroutine rot_to_asym

    !>  \brief  rotates orientations to the asymmetric unit
    subroutine rotall_to_asym( self, osyms )
        class(sym),  intent(inout) :: self
        class(oris), intent(inout) :: osyms     !< orientations
        type(ori) :: o
        integer   :: i
        do i = 1, osyms%get_noris()
            call osyms%get_ori(i, o)
            call self%rot_to_asym(o)
            call osyms%set_ori(i, o)
        enddo
        call o%kill
    end subroutine rotall_to_asym

    !>  \brief  determines euler distance and corresponding symmetrized
    !>          orientation between two orientations
    subroutine sym_dists_1( self, oref, oasym, osym, euldist, inplrotdist )
        class(sym), intent(inout) :: self
        class(ori), intent(in)    :: oref        !< is the untouched reference
        class(ori), intent(in)    :: oasym       !< is the orientation determined within assymetric unit
        class(ori), intent(inout) :: osym        !< is the orientatiom that minimises the distance to oref
        real,       intent(out)   :: euldist     !< Euler distance
        real,       intent(out)   :: inplrotdist !< in-plane rotational distance
        type(ori) :: o
        real      :: dist
        integer   :: isym
        euldist     = oasym.euldist.oref
        inplrotdist = oasym.inplrotdist.oref
        osym        = oasym
        if( self%n == 1 )then
            ! C1: nothing to do
        else
            do isym=2,self%n
                call self%apply(oasym, isym, o)
                dist = o.euldist.oref
                if( dist < euldist )then
                    euldist = dist
                    osym    = o
                endif
            enddo
            inplrotdist = osym.inplrotdist.oref
        endif
        euldist     = rad2deg( euldist )
        inplrotdist = rad2deg( inplrotdist )
        call o%kill
    end subroutine sym_dists_1

        !>  \brief  determines euler distance and corresponding symmetrized
    !>          orientation between two orientations
    subroutine sym_dists_2( self, euls_ref, euls_asym, euls_sym, euldist, inplrotdist )
        class(sym), intent(in)    :: self
        real,       intent(in)    :: euls_ref(3), euls_asym(3)
        real,       intent(inout) :: euls_sym(3)
        real,       intent(out)   :: euldist     !< Euler distance
        real,       intent(out)   :: inplrotdist !< in-plane rotational distance
        real    :: euls(3), dist
        integer :: isym
        euldist     = euler_dist(euls_asym, euls_ref)
        inplrotdist = euler_inplrotdist(euls_asym, euls_ref)
        euls_sym    = euls_asym
        if( self%n == 1 )then
            ! C1: nothing to do
        else
            do isym=2,self%n
                call self%apply(euls_asym, isym, euls)
                dist = euler_dist(euls, euls_ref)
                if( dist < euldist )then
                    euldist = dist
                    euls_sym  = euls
                endif
            enddo
            inplrotdist = euler_inplrotdist(euls_sym, euls_ref)
        endif
        euldist     = rad2deg( euldist )
        inplrotdist = rad2deg( inplrotdist )
    end subroutine sym_dists_2

    !>  \brief  is a getter
    subroutine get_symori( self, symop, e_sym )
        class(sym), intent(inout) :: self
        integer,    intent(in)    :: symop !< symmetry orientation
        type(ori),  intent(inout) :: e_sym
        call self%e_sym%get_ori(symop, e_sym)
    end subroutine get_symori

    !>  \brief  whether or not an orientation falls within the asymetric unit excluding mirror
    logical function within_asymunit( self, o, incl_mirr )
        class(sym),        intent(inout) :: self
        class(ori),        intent(in)    :: o !< symmetry orientation
        logical, optional, intent(in)    :: incl_mirr
        real    :: e1, e2
        logical :: incl_mirr_here
        within_asymunit = .false.
        incl_mirr_here  = .true.
        if(present(incl_mirr))incl_mirr_here = incl_mirr
        e1 = o%e1get()
        if( incl_mirr_here )then
            if(e1 <  self%eullims(1,1) )return
            if(e1 >= self%eullims(1,2) )return
            e2 = o%e2get()
            if(e2 <  self%eullims(2,1) )return
            if(e2 >= self%eullims(2,2) )return
        else
            if(e1 <  self%eullims_nomirr(1,1) )return
            if(e1 >= self%eullims_nomirr(1,2) )return
            e2 = o%e2get()
            if(e2 <  self%eullims_nomirr(2,1) )return
            if(e2 >= self%eullims_nomirr(2,2) )return
        endif
        within_asymunit = self%within_platonic_asymunit(e1,e2,incl_mirr_here)
    end function within_asymunit

    !>  Modified from SPARX
    logical function within_platonic_asymunit(self, e1, e2, incl_mirr)
        class(sym), intent(inout) :: self
        real,       intent(in)    :: e1,e2
        logical,    intent(in)    :: incl_mirr
        real :: e1_here
        within_platonic_asymunit = .true.
        if( self%c_or_d )return
        e1_here = min(e1, self%eullims(1,2)-e1)
        select case(self%n)
            case(24)    ! o
                within_platonic_asymunit = e2 < e2_bound(self%eullims(1,2),45.0,oct_theta_lim)
            case(60)    ! i
                within_platonic_asymunit = e2 < e2_bound(self%eullims(1,2),31.717474411461005,ico_theta_lim)
            case(12)    ! t
                if(e2 > e2_bound(self%eullims(1,2),oct_theta_lim,tet_theta_lim))then
                    within_platonic_asymunit = .false.
                else
                    if(.not.incl_mirr)then
                        within_platonic_asymunit = e2 < e2_upper_bound(self%eullims(1,2),real(gamma),tet_theta_lim)
                    endif
                endif
        end select

        contains

            real function e2_bound( ang1, ang2, ang3 )
                real, intent(in) :: ang1, ang2, ang3
                e2_bound = sin(deg2rad(ang1/2.-e1_here)) / tan(deg2rad(ang2))
                e2_bound = e2_bound + sin(deg2rad(e1_here)) / tan(deg2rad(ang3))
                e2_bound = e2_bound / sin(deg2rad(ang1/2.))
                e2_bound = rad2deg(atan(1./e2_bound))
            end function e2_bound

            real function e2_upper_bound( ang1, ang2, ang3 )
                real, intent(in) :: ang1, ang2, ang3
                e2_upper_bound = sin(deg2rad(ang1/2.-e1_here)) / tan(deg2rad(ang2))
                e2_upper_bound = e2_upper_bound + sin(deg2rad(e1_here)) / tan(deg2rad(ang3/2.))
                e2_upper_bound = e2_upper_bound / sin(deg2rad(ang1/2.))
                e2_upper_bound = rad2deg(atan(1./e2_upper_bound))
            end function e2_upper_bound

    end function within_platonic_asymunit

    !>  apply symmetry orientations with shift
    subroutine apply_sym_with_shift( self, os, symaxis_ori, shvec, state )
        class(sym),        intent(inout) :: self
        class(oris),       intent(inout) :: os          !< output orientations
        class(ori),        intent(in)    :: symaxis_ori !< input orientations
        real,              intent(in)    :: shvec(3)    !< shift vector
        integer, optional, intent(in)    :: state       !< current state
        type(ori) :: o
        integer   :: i, s
        if( present(state) )then
            do i=1,os%get_noris()
                s = nint(os%get(i, 'state'))
                if(s .ne. state) cycle
                call os%map3dshift22d(i, shvec)
                ! transposed rotation to get the correct sign on rotation
                ! the old fetching versus inserting issue
                call os%rot_transp(i,symaxis_ori)
                call os%get_ori(i, o)
                call self%rot_to_asym(o)
                call os%set_ori(i, o)
            end do
        else
            call os%map3dshift22d( shvec )
            ! transposed rotation to get the correct sign on rotation
            ! the old fetching versus inserting issue
            call os%rot_transp(symaxis_ori)
            call self%rotall_to_asym(os)
        endif
        call o%kill
    end subroutine apply_sym_with_shift

    !>  \brief  to find the closest matching projection direction
    !! KEEP THIS ROUTINE SERIAL
    function find_closest_proj( self, os_asym_unit, o_in ) result( closest )
        class(sym),  intent(inout) :: self
        class(oris), intent(in)    :: os_asym_unit
        class(ori),  intent(in)    :: o_in
        real      :: dists(os_asym_unit%get_noris())
        integer   :: closest, i, isym
        type(ori) :: oasym, osym
        if( trim(self%pgrp).eq.'c1' )then
            closest = os_asym_unit%find_closest_proj(o_in)
        else
            dists = pi
            do i=1,os_asym_unit%get_noris()
                call os_asym_unit%get_ori(i, oasym)
                do isym = 1, self%n
                    if(isym == 1)then
                        call osym%copy(oasym)
                    else
                        call self%apply(oasym, isym, osym)
                    endif
                    dists(i) = min(dists(i), osym.euldist.o_in)
                end do
            end do
            closest = minloc(dists, dim=1)
        endif
        call oasym%kill
        call osym%kill
    end function find_closest_proj

    !>  \brief  is for retrieving nearest neighbors in symmetric cases
    !! the policy here is based solely on angular distance and initialization of lnns is
    !! deferred to the calling unit, so that we can add additional neighborhoods on top of
    !! of each other to create more complex search spaces
    subroutine nearest_proj_neighbors( self, os_asym_unit, o, euldist_thres, lnns )
        class(sym), intent(inout) :: self
        type(oris), intent(inout) :: os_asym_unit !< sampled orientations from assymetric unit, eg from spiral with symmetry
        class(ori), intent(in)    :: o
        real,       intent(in)    :: euldist_thres ! in degrees
        logical,    intent(inout) :: lnns(os_asym_unit%get_noris())
        real      :: dists(os_asym_unit%get_noris()), euldist_thres_rad
        integer   :: i, isym
        type(ori) :: oasym, osym
        if( trim(self%pgrp).eq.'c1' )then
            call os_asym_unit%nearest_proj_neighbors(o, euldist_thres, lnns)
        else
            euldist_thres_rad = deg2rad(euldist_thres)
            dists = pi
            do i=1,os_asym_unit%get_noris()
                call os_asym_unit%get_ori(i, oasym)
                do isym = 1, self%n
                    if(isym == 1)then
                        call osym%copy(oasym)
                    else
                        call self%apply(oasym, isym, osym)
                    endif
                    dists(i) = min(dists(i), osym.euldist.o)
                end do
            end do
            where(dists <= euldist_thres_rad) lnns = .true.
        endif
        call oasym%kill
        call osym%kill
    end subroutine nearest_proj_neighbors

    !>  \brief  is for retrieving nearest symmetry neighbours in an assymetric set of projection directions
    subroutine nearest_sym_neighbors( self, asym_os, nnmat )
        class(sym),           intent(inout) :: self
        type(oris),           intent(inout) :: asym_os    !< C1 projection directions
        integer, allocatable, intent(inout) :: nnmat(:,:) !< nearest-neighbour matrix
        real,    allocatable :: dists(:)
        type(ori)  :: oasym, osym, oj
        integer    :: i, j, n_os, isym, loc(1), nsym
        if( trim(self%pgrp).eq.'c1' )then
            call asym_os%nearest_proj_neighbors(1, nnmat)
        else
            n_os = asym_os%get_noris()
            nsym = self%n
            allocate( nnmat(n_os,nsym), dists(n_os) )
            do i = 1, n_os
                call asym_os%get_ori(i, oasym)
                do isym = 1, nsym
                    if( isym == 1 )then
                        osym = oasym
                    else
                        call self%apply(oasym, isym, osym)
                    endif
                    do j = 1, n_os
                        call asym_os%get_ori(j, oj)
                        dists(j) = osym.euldist.oj
                    enddo
                    loc = minloc(dists)
                    nnmat(i,isym) = loc(1)
                enddo
            enddo
            deallocate(dists)
        endif
        call oasym%kill
        call osym%kill
        call oj%kill
    end subroutine nearest_sym_neighbors

    !>  \brief  is for randomizing over the symmetry operations of the point-group
    subroutine symrandomize( self, os_asym_unit )
        class(sym), intent(inout) :: self
        type(oris), intent(inout) :: os_asym_unit !< sampled orientations from assymetric unit, eg from spiral with symmetry
        integer   :: n_os, nsym, symop, i
        type(ori) :: oasym, osym
        if( trim(self%pgrp).eq.'c1' )then
            THROW_HARD('nothing to do, pgrp=c1; symrandomize')
        endif
        n_os = os_asym_unit%get_noris()
        nsym = self%n
        do i = 1, n_os
            call os_asym_unit%get_ori(i, oasym)
            symop = irnd_uni(nsym)
            if( symop > 1 )then
                call self%apply(oasym, symop, osym)
                call os_asym_unit%set_ori(i,osym)
            endif
        end do
        call oasym%kill
        call osym%kill
    end subroutine symrandomize

    !>  \brief  is for building a spiral INCLUDING mirror projections
    subroutine build_refspiral( self, os )
        class(sym), intent(inout) :: self
        type(oris), intent(inout) :: os
        type(oris) :: tmp, os_nomirr
        integer    :: cnt, i, n, nprojs, lim, nos, nos_nomirr
        logical, allocatable :: avail(:)
        nos = os%get_noris()
        if(is_odd(nos))THROW_HARD('odd number of projections directions not supported; build_refspiral')
        nos_nomirr = nos/2
        call os%new(nos, is_ptcl=.false.)
        call os_nomirr%new(nos_nomirr, is_ptcl=.false.)
        n = self%n * nos
        call gen_c1
        nprojs = count(avail)
        if( nprojs < nos_nomirr )then
            ! under sampling
            n = n + nint(real(nos_nomirr)/2.)
            call gen_c1
            nprojs = count(avail)
        endif
        if( nprojs > nos_nomirr )then
            ! over sampling
            lim = max(2, floor(real(nprojs)/real(nprojs-nos_nomirr+1)) )
            cnt  = 0
            do i=1,n
                if(.not.avail(i))cycle
                cnt  = cnt + 1
                if(cnt == lim) then
                    avail(i) = .false.
                    cnt      = 0
                    nprojs   = nprojs-1
                    if( nprojs == nos_nomirr )exit
                endif
            enddo
        endif
        ! copy un-mirrored asu
        cnt = 0
        do i=1,n
            if( .not.avail(i) )cycle
            cnt = cnt + 1
            if(cnt > nos_nomirr)exit
            call os_nomirr%set_euler(cnt, tmp%get_euler(i))
        enddo
        ! stash un-mirrored
        do i=1,nos_nomirr
            call os%set_euler(i, os_nomirr%get_euler(i))
        enddo
        ! mirror
        call os_nomirr%mirror2d
        ! now apply symmetry to mirrored projection directions
        ! such that they fall within the mirrored asymmetric unit
        if( self%n > 1 )then
            call self%rotall_to_asym(os_nomirr)
        endif
        do i=1,nos_nomirr
            call os_nomirr%e3set(i,0.)
        enddo
        ! append
        cnt = 0
        do i=nos_nomirr+1,nos
            cnt = cnt + 1
            call os%set_euler(i, os_nomirr%get_euler(cnt))
        enddo
        ! cleanup
        call os_nomirr%kill
        call tmp%kill
        deallocate(avail)

        contains

            subroutine gen_c1
                integer   :: j, north_pole_ind
                type(ori) :: o, o2, o_tmp, north_pole
                real      :: min_dist, dist
                if( allocated(avail) )deallocate(avail)
                allocate(avail(n), source=.false.)
                call tmp%new(n, is_ptcl=.false.)
                call tmp%spiral
                call north_pole%new(is_ptcl=.false.)
                call north_pole%set_euler([0.,0.,0.])
                north_pole_ind = 0
                min_dist = PI
                do j=1,n
                    call tmp%get_ori(j, o)
                    dist = north_pole.euldist.o
                    if(dist < 0.001 )then !in radians
                        ! the north pole shall always be present
                        avail(j) = .true.
                        north_pole_ind = j
                    else
                        call tmp%get_ori(j, o2)
                        avail(j) = self%within_asymunit(o2, incl_mirr=.false.)
                        if(avail(j)) min_dist = min(dist,min_dist)
                    endif
                end do
                select case(self%pgrp(1:1))
                case('d','i','o')
                    ! this is to jitter the north pole such that
                    ! its mirror symmetric is not the north pole
                    if(north_pole_ind>0)then
                        min_dist = min(0.5,rad2deg(min_dist)/5.)
                        call o%set_euler([min_dist,min_dist,0.])
                        call o%compose(north_pole, o_tmp)
                        o_tmp = o
                        call tmp%set_euler(north_pole_ind,o%get_euler())
                    endif
                case DEFAULT
                    ! all is good
                end select
                call o%kill
                call o2%kill
                call o_tmp%kill
                call north_pole%kill
            end subroutine gen_c1

    end subroutine build_refspiral

    subroutine rnd_euler_1( self, osym )
        class(sym), intent(inout) :: self
        class(ori), intent(inout) :: osym
        call osym%rnd_euler(self%eullims)
    end subroutine rnd_euler_1

    subroutine rnd_euler_2( self, o_prev, athres, osym )
        class(sym),     intent(inout) :: self
        class(ori),     intent(in)    :: o_prev
        real,           intent(in)    :: athres
        class(ori),     intent(inout) :: osym
        call osym%rnd_euler(o_prev, athres, self%eullims)
    end subroutine rnd_euler_2

    subroutine rnd_euler_3( self, o_prev, athres_proj, athres_inpl, osym )
        class(sym),     intent(inout) :: self
        class(ori),     intent(in)    :: o_prev
        real,           intent(in)    :: athres_proj, athres_inpl
        class(ori),     intent(inout) :: osym
        call osym%rnd_euler(o_prev, athres_proj, athres_inpl, self%eullims)
    end subroutine rnd_euler_3

    !>  \brief  SPIDER code for making c and d symmetries
    subroutine make_c_and_d( self )
        class(sym), intent(inout) :: self
        double precision :: delta, degree
        double precision, dimension(3,3) :: a,b,g
        integer   :: i,j,cnt
        delta = 360.d0/dble(self%ncsym)
        cnt = 0
        do i=0,1
           degree = dble(i)*delta2
           a = matcreate(0, degree)
           do j=0,self%ncsym-1
              cnt = cnt+1
              degree = dble(j)*delta
              b = matcreate(1, degree)
              g = matmul(a,b) ! this is the c-symmetry
              call self%e_sym%set_euler(cnt,real(matextract(g)))
           end do
           if(self%pgrp(1:1).ne.'d' .and. self%pgrp(1:1).ne.'D') return
        end do
    end subroutine make_c_and_d

    !>  \brief  hardcoded euler angles taken from SPARX
    subroutine make_o( self )
        class(sym), intent(inout) :: self
        integer   :: i,j,cnt
        double precision :: psi, phi
        cnt = 0
        do i = 1,4
            phi = dble(i-1) * 90.d0
            cnt = cnt + 1
            call self%e_sym%set_euler(cnt,real([0.d0, 0.d0, phi]))
            do j = 1,4
                psi = dble(j-1) * 90.d0
                cnt = cnt + 1
                call self%e_sym%set_euler(cnt,real([psi, 90.d0, phi]))
            end do
            cnt = cnt + 1
            call self%e_sym%set_euler(cnt,real([0.d0, 180.d0, phi]))
        enddo
    end subroutine make_o

    !>  \brief  SPIDER code for making tetahedral symmetry
    !>          tetrahedral, with 3axis align w/z axis, point on +ve x axis
    subroutine make_t( self )
        class(sym), intent(inout) :: self
        double precision :: tester,dt,degree,psi,theta,phi
        integer :: i,j,cnt
        ! value from (90 -dihedral angle) + 90 =? 109.47
        tester = 1.d0/3.d0
        ! degree = 180.0 - (acos(tester) * (180.0/pi))
        dt = max(-1.0d0,min(1.0d0,tester))
        degree = 180.d0-(dacos(dt)*(180.d0/dpi))
        cnt = 0
        do i=0,2
            cnt = cnt+1
            psi   = 0.d0
            theta = 0.d0
            phi   = dble(i)*delta3
            call self%e_sym%set_euler(cnt,real([psi,theta,phi]))
            psi = phi
            do j=0,2
                cnt   = cnt+1
                phi   = 60.d0+dble(j)*delta3
                theta = degree
                call self%e_sym%set_euler(cnt,real([psi,theta,phi]))
            end do
        end do
    end subroutine make_t

    !>  \brief  SPIDER code for making icosahedral symmetry
    subroutine make_i_spider( self )
        class(sym), intent(inout) :: self
        double precision :: psi, theta, phi
        integer   :: i, j, cnt
        cnt = 0
        do i=0,1
            do j=0,4
                cnt   = cnt + 1
                psi   = 0.d0
                theta = dble(i) * delta2
                phi   = dble(j) * delta5
                call self%e_sym%set_euler(cnt, real([psi,theta,phi]))
            end do
        end do
        theta = 63.43494882292201d0
        do i=0,4
            do j=0,4
                cnt = cnt + 1
                psi = dble(i) * delta5
                phi = dble(j) * delta5 + deltan
                call self%e_sym%set_euler(cnt, real([psi,theta,phi]))
            end do
        end do
        theta = 116.56505117707799d0
        do i=0,4
            do j=0,4
                cnt = cnt + 1
                psi = dble(i) * delta5 + deltan
                phi = dble(j) * delta5
                call self%e_sym%set_euler(cnt, real([psi,theta,phi]))
            end do
        end do
    end subroutine make_i_spider

    subroutine make_i_relion( self )
        class(sym), intent(inout) :: self
        real :: rmat(3,3)
        type(ori) :: o
        ! Symmetry operation: 1
        rmat(1,1) = 1
        rmat(1,2) = 0
        rmat(1,3) = 0
        rmat(2,1) = 0
        rmat(2,2) = 1
        rmat(2,3) = 0
        rmat(3,1) = 0
        rmat(3,2) = 0
        rmat(3,3) = 1
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(1, o)

        ! Symmetry operation: 2
        rmat(1,1) = -1
        rmat(1,2) = 0
        rmat(1,3) = 0
        rmat(2,1) = 0
        rmat(2,2) = -1
        rmat(2,3) = 0
        rmat(3,1) = 0
        rmat(3,2) = 0
        rmat(3,3) = 1
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(2, o)

        ! Symmetry operation: 3
        rmat(1,1) = 0.5
        rmat(1,2) = 0.80902
        rmat(1,3) = 0.30902
        rmat(2,1) = -0.80902
        rmat(2,2) = 0.30902
        rmat(2,3) = 0.5
        rmat(3,1) = 0.30902
        rmat(3,2) = -0.5
        rmat(3,3) = 0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(3, o)

        ! Symmetry operation: 4
        rmat(1,1) = -0.30902
        rmat(1,2) = 0.5
        rmat(1,3) = 0.80902
        rmat(2,1) = -0.5
        rmat(2,2) = -0.80902
        rmat(2,3) = 0.30902
        rmat(3,1) = 0.80902
        rmat(3,2) = -0.30902
        rmat(3,3) = 0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(4, o)

        ! Symmetry operation: 5
        rmat(1,1) = -0.30902
        rmat(1,2) = -0.5
        rmat(1,3) = 0.80902
        rmat(2,1) = 0.5
        rmat(2,2) = -0.80902
        rmat(2,3) = -0.30902
        rmat(3,1) = 0.80902
        rmat(3,2) = 0.30902
        rmat(3,3) = 0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(5, o)

        ! Symmetry operation: 6
        rmat(1,1) = 0.5
        rmat(1,2) = -0.80902
        rmat(1,3) = 0.30902
        rmat(2,1) = 0.80902
        rmat(2,2) = 0.30902
        rmat(2,3) = -0.5
        rmat(3,1) = 0.30902
        rmat(3,2) = 0.5
        rmat(3,3) = 0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(6, o)

        ! Symmetry operation: 7
        rmat(1,1) = -0.5
        rmat(1,2) = 0.80902
        rmat(1,3) = -0.30902
        rmat(2,1) = -0.80902
        rmat(2,2) = -0.30902
        rmat(2,3) = 0.5
        rmat(3,1) = 0.30902
        rmat(3,2) = 0.5
        rmat(3,3) = 0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(7, o)

        ! Symmetry operation: 8
        rmat(1,1) = -0.5
        rmat(1,2) = -0.80902
        rmat(1,3) = 0.30902
        rmat(2,1) = 0.80902
        rmat(2,2) = -0.30902
        rmat(2,3) = 0.5
        rmat(3,1) = -0.30902
        rmat(3,2) = 0.5
        rmat(3,3) = 0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(8, o)

        ! Symmetry operation: 9
        rmat(1,1) = -0.5
        rmat(1,2) = -0.80902
        rmat(1,3) = -0.30902
        rmat(2,1) = 0.80902
        rmat(2,2) = -0.30902
        rmat(2,3) = -0.5
        rmat(3,1) = 0.30902
        rmat(3,2) = -0.5
        rmat(3,3) = 0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(9, o)

        ! Symmetry operation: 10
        rmat(1,1) = 0.30902
        rmat(1,2) = -0.5
        rmat(1,3) = -0.80902
        rmat(2,1) = 0.5
        rmat(2,2) = 0.80902
        rmat(2,3) = -0.30902
        rmat(3,1) = 0.80902
        rmat(3,2) = -0.30902
        rmat(3,3) = 0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(10, o)

        ! Symmetry operation: 11
        rmat(1,1) = 0.30902
        rmat(1,2) = -0.5
        rmat(1,3) = 0.80902
        rmat(2,1) = 0.5
        rmat(2,2) = 0.80902
        rmat(2,3) = 0.30902
        rmat(3,1) = -0.80902
        rmat(3,2) = 0.30902
        rmat(3,3) = 0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(11, o)

        ! Symmetry operation: 12
        rmat(1,1) = 0.30902
        rmat(1,2) = 0.5
        rmat(1,3) = -0.80902
        rmat(2,1) = -0.5
        rmat(2,2) = 0.80902
        rmat(2,3) = 0.30902
        rmat(3,1) = 0.80902
        rmat(3,2) = 0.30902
        rmat(3,3) = 0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(12, o)

        ! Symmetry operation: 13
        rmat(1,1) = 0.30902
        rmat(1,2) = 0.5
        rmat(1,3) = 0.80902
        rmat(2,1) = -0.5
        rmat(2,2) = 0.80902
        rmat(2,3) = -0.30902
        rmat(3,1) = -0.80902
        rmat(3,2) = -0.30902
        rmat(3,3) = 0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(13, o)

        ! Symmetry operation: 14
        rmat(1,1) = -0.5
        rmat(1,2) = 0.80902
        rmat(1,3) = 0.30902
        rmat(2,1) = -0.80902
        rmat(2,2) = -0.30902
        rmat(2,3) = -0.5
        rmat(3,1) = -0.30902
        rmat(3,2) = -0.5
        rmat(3,3) = 0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(14, o)

        ! Symmetry operation: 15
        rmat(1,1) = -0.80902
        rmat(1,2) = 0.30902
        rmat(1,3) = 0.5
        rmat(2,1) = 0.30902
        rmat(2,2) = -0.5
        rmat(2,3) = 0.80902
        rmat(3,1) = 0.5
        rmat(3,2) = 0.80902
        rmat(3,3) = 0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(15, o)

        ! Symmetry operation: 16
        rmat(1,1) = 0
        rmat(1,2) = 0
        rmat(1,3) = 1
        rmat(2,1) = 1
        rmat(2,2) = 0
        rmat(2,3) = 0
        rmat(3,1) = 0
        rmat(3,2) = 1
        rmat(3,3) = 0
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(16, o)

        ! Symmetry operation: 17
        rmat(1,1) = 0.80902
        rmat(1,2) = 0.30902
        rmat(1,3) = 0.5
        rmat(2,1) = 0.30902
        rmat(2,2) = 0.5
        rmat(2,3) = -0.80902
        rmat(3,1) = -0.5
        rmat(3,2) = 0.80902
        rmat(3,3) = 0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(17, o)

        ! Symmetry operation: 18
        rmat(1,1) = 0.5
        rmat(1,2) = 0.80902
        rmat(1,3) = -0.30902
        rmat(2,1) = -0.80902
        rmat(2,2) = 0.30902
        rmat(2,3) = -0.5
        rmat(3,1) = -0.30902
        rmat(3,2) = 0.5
        rmat(3,3) = 0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(18, o)

        ! Symmetry operation: 19
        rmat(1,1) = 0.5
        rmat(1,2) = -0.80902
        rmat(1,3) = -0.30902
        rmat(2,1) = 0.80902
        rmat(2,2) = 0.30902
        rmat(2,3) = 0.5
        rmat(3,1) = -0.30902
        rmat(3,2) = -0.5
        rmat(3,3) = 0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(19, o)

        ! Symmetry operation: 20
        rmat(1,1) = 0
        rmat(1,2) = 1
        rmat(1,3) = 0
        rmat(2,1) = 0
        rmat(2,2) = 0
        rmat(2,3) = 1
        rmat(3,1) = 1
        rmat(3,2) = 0
        rmat(3,3) = 0
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(20, o)

        ! Symmetry operation: 21
        rmat(1,1) = 0.80902
        rmat(1,2) = 0.30902
        rmat(1,3) = -0.5
        rmat(2,1) = 0.30902
        rmat(2,2) = 0.5
        rmat(2,3) = 0.80902
        rmat(3,1) = 0.5
        rmat(3,2) = -0.80902
        rmat(3,3) = 0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(21, o)

        ! Symmetry operation: 22
        rmat(1,1) = 0.80902
        rmat(1,2) = -0.30902
        rmat(1,3) = 0.5
        rmat(2,1) = -0.30902
        rmat(2,2) = 0.5
        rmat(2,3) = 0.80902
        rmat(3,1) = -0.5
        rmat(3,2) = -0.80902
        rmat(3,3) = 0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(22, o)

        ! Symmetry operation: 23
        rmat(1,1) = 0
        rmat(1,2) = 0
        rmat(1,3) = 1
        rmat(2,1) = -1
        rmat(2,2) = 0
        rmat(2,3) = 0
        rmat(3,1) = 0
        rmat(3,2) = -1
        rmat(3,3) = 0
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(23, o)

        ! Symmetry operation: 24
        rmat(1,1) = -0.80902
        rmat(1,2) = -0.30902
        rmat(1,3) = 0.5
        rmat(2,1) = -0.30902
        rmat(2,2) = -0.5
        rmat(2,3) = -0.80902
        rmat(3,1) = 0.5
        rmat(3,2) = -0.80902
        rmat(3,3) = 0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(24, o)

        ! Symmetry operation: 25
        rmat(1,1) = -0.30902
        rmat(1,2) = 0.5
        rmat(1,3) = -0.80902
        rmat(2,1) = -0.5
        rmat(2,2) = -0.80902
        rmat(2,3) = -0.30902
        rmat(3,1) = -0.80902
        rmat(3,2) = 0.30902
        rmat(3,3) = 0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(25, o)

        ! Symmetry operation: 26
        rmat(1,1) = 0.80902
        rmat(1,2) = -0.30902
        rmat(1,3) = -0.5
        rmat(2,1) = -0.30902
        rmat(2,2) = 0.5
        rmat(2,3) = -0.80902
        rmat(3,1) = 0.5
        rmat(3,2) = 0.80902
        rmat(3,3) = 0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(26, o)

        ! Symmetry operation: 27
        rmat(1,1) = 0.80902
        rmat(1,2) = 0.30902
        rmat(1,3) = 0.5
        rmat(2,1) = -0.30902
        rmat(2,2) = -0.5
        rmat(2,3) = 0.80902
        rmat(3,1) = 0.5
        rmat(3,2) = -0.80902
        rmat(3,3) = -0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(27, o)

        ! Symmetry operation: 28
        rmat(1,1) = 0.30902
        rmat(1,2) = -0.5
        rmat(1,3) = 0.80902
        rmat(2,1) = -0.5
        rmat(2,2) = -0.80902
        rmat(2,3) = -0.30902
        rmat(3,1) = 0.80902
        rmat(3,2) = -0.30902
        rmat(3,3) = -0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(28, o)

        ! Symmetry operation: 29
        rmat(1,1) = 0
        rmat(1,2) = -1
        rmat(1,3) = 0
        rmat(2,1) = 0
        rmat(2,2) = 0
        rmat(2,3) = -1
        rmat(3,1) = 1
        rmat(3,2) = 0
        rmat(3,3) = 0
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(29, o)

        ! Symmetry operation: 30
        rmat(1,1) = -0.80902
        rmat(1,2) = -0.30902
        rmat(1,3) = -0.5
        rmat(2,1) = 0.30902
        rmat(2,2) = 0.5
        rmat(2,3) = -0.80902
        rmat(3,1) = 0.5
        rmat(3,2) = -0.80902
        rmat(3,3) = -0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(30, o)

        ! Symmetry operation: 31
        rmat(1,1) = -0.80902
        rmat(1,2) = 0.30902
        rmat(1,3) = -0.5
        rmat(2,1) = 0.30902
        rmat(2,2) = -0.5
        rmat(2,3) = -0.80902
        rmat(3,1) = -0.5
        rmat(3,2) = -0.80902
        rmat(3,3) = 0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(31, o)

        ! Symmetry operation: 32
        rmat(1,1) = -0.30902
        rmat(1,2) = -0.5
        rmat(1,3) = -0.80902
        rmat(2,1) = 0.5
        rmat(2,2) = -0.80902
        rmat(2,3) = 0.30902
        rmat(3,1) = -0.80902
        rmat(3,2) = -0.30902
        rmat(3,3) = 0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(32, o)

        ! Symmetry operation: 33
        rmat(1,1) = 0
        rmat(1,2) = 0
        rmat(1,3) = -1
        rmat(2,1) = -1
        rmat(2,2) = 0
        rmat(2,3) = 0
        rmat(3,1) = 0
        rmat(3,2) = 1
        rmat(3,3) = 0
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(33, o)

        ! Symmetry operation: 34
        rmat(1,1) = -0.80902
        rmat(1,2) = -0.30902
        rmat(1,3) = -0.5
        rmat(2,1) = -0.30902
        rmat(2,2) = -0.5
        rmat(2,3) = 0.80902
        rmat(3,1) = -0.5
        rmat(3,2) = 0.80902
        rmat(3,3) = 0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(34, o)

        ! Symmetry operation: 35
        rmat(1,1) = -0.80902
        rmat(1,2) = -0.30902
        rmat(1,3) = 0.5
        rmat(2,1) = 0.30902
        rmat(2,2) = 0.5
        rmat(2,3) = 0.80902
        rmat(3,1) = -0.5
        rmat(3,2) = 0.80902
        rmat(3,3) = -0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(35, o)

        ! Symmetry operation: 36
        rmat(1,1) = 0.30902
        rmat(1,2) = 0.5
        rmat(1,3) = 0.80902
        rmat(2,1) = 0.5
        rmat(2,2) = -0.80902
        rmat(2,3) = 0.30902
        rmat(3,1) = 0.80902
        rmat(3,2) = 0.30902
        rmat(3,3) = -0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(36, o)

        ! Symmetry operation: 37
        rmat(1,1) = 0.80902
        rmat(1,2) = -0.30902
        rmat(1,3) = 0.5
        rmat(2,1) = 0.30902
        rmat(2,2) = -0.5
        rmat(2,3) = -0.80902
        rmat(3,1) = 0.5
        rmat(3,2) = 0.80902
        rmat(3,3) = -0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(37, o)

        ! Symmetry operation: 38
        rmat(1,1) = -0.30902
        rmat(1,2) = -0.5
        rmat(1,3) = -0.80902
        rmat(2,1) = -0.5
        rmat(2,2) = 0.80902
        rmat(2,3) = -0.30902
        rmat(3,1) = 0.80902
        rmat(3,2) = 0.30902
        rmat(3,3) = -0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(38, o)

        ! Symmetry operation: 39
        rmat(1,1) = -0.80902
        rmat(1,2) = 0.30902
        rmat(1,3) = -0.5
        rmat(2,1) = -0.30902
        rmat(2,2) = 0.5
        rmat(2,3) = 0.80902
        rmat(3,1) = 0.5
        rmat(3,2) = 0.80902
        rmat(3,3) = -0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(39, o)

        ! Symmetry operation: 40
        rmat(1,1) = -0.30902
        rmat(1,2) = 0.5
        rmat(1,3) = -0.80902
        rmat(2,1) = 0.5
        rmat(2,2) = 0.80902
        rmat(2,3) = 0.30902
        rmat(3,1) = 0.80902
        rmat(3,2) = -0.30902
        rmat(3,3) = -0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(40, o)

        ! Symmetry operation: 41
        rmat(1,1) = 0
        rmat(1,2) = 0
        rmat(1,3) = -1
        rmat(2,1) = 1
        rmat(2,2) = 0
        rmat(2,3) = 0
        rmat(3,1) = 0
        rmat(3,2) = -1
        rmat(3,3) = 0
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(41, o)

        ! Symmetry operation: 42
        rmat(1,1) = 0
        rmat(1,2) = -1
        rmat(1,3) = 0
        rmat(2,1) = 0
        rmat(2,2) = 0
        rmat(2,3) = 1
        rmat(3,1) = -1
        rmat(3,2) = 0
        rmat(3,3) = 0
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(42, o)

        ! Symmetry operation: 43
        rmat(1,1) = -0.30902
        rmat(1,2) = -0.5
        rmat(1,3) = 0.80902
        rmat(2,1) = -0.5
        rmat(2,2) = 0.80902
        rmat(2,3) = 0.30902
        rmat(3,1) = -0.80902
        rmat(3,2) = -0.30902
        rmat(3,3) = -0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(43, o)

        ! Symmetry operation: 44
        rmat(1,1) = -0.80902
        rmat(1,2) = 0.30902
        rmat(1,3) = 0.5
        rmat(2,1) = -0.30902
        rmat(2,2) = 0.5
        rmat(2,3) = -0.80902
        rmat(3,1) = -0.5
        rmat(3,2) = -0.80902
        rmat(3,3) = -0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(44, o)

        ! Symmetry operation: 45
        rmat(1,1) = -0.30902
        rmat(1,2) = 0.5
        rmat(1,3) = 0.80902
        rmat(2,1) = 0.5
        rmat(2,2) = 0.80902
        rmat(2,3) = -0.30902
        rmat(3,1) = -0.80902
        rmat(3,2) = 0.30902
        rmat(3,3) = -0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(45, o)

        ! Symmetry operation: 46
        rmat(1,1) = 0
        rmat(1,2) = 1
        rmat(1,3) = 0
        rmat(2,1) = 0
        rmat(2,2) = 0
        rmat(2,3) = -1
        rmat(3,1) = -1
        rmat(3,2) = 0
        rmat(3,3) = 0
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(46, o)

        ! Symmetry operation: 47
        rmat(1,1) = -0.5
        rmat(1,2) = 0.80902
        rmat(1,3) = 0.30902
        rmat(2,1) = 0.80902
        rmat(2,2) = 0.30902
        rmat(2,3) = 0.5
        rmat(3,1) = 0.30902
        rmat(3,2) = 0.5
        rmat(3,3) = -0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(47, o)

        ! Symmetry operation: 48
        rmat(1,1) = 0.5
        rmat(1,2) = 0.80902
        rmat(1,3) = 0.30902
        rmat(2,1) = 0.80902
        rmat(2,2) = -0.30902
        rmat(2,3) = -0.5
        rmat(3,1) = -0.30902
        rmat(3,2) = 0.5
        rmat(3,3) = -0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(48, o)

        ! Symmetry operation: 49
        rmat(1,1) = 0.80902
        rmat(1,2) = 0.30902
        rmat(1,3) = -0.5
        rmat(2,1) = -0.30902
        rmat(2,2) = -0.5
        rmat(2,3) = -0.80902
        rmat(3,1) = -0.5
        rmat(3,2) = 0.80902
        rmat(3,3) = -0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(49, o)

        ! Symmetry operation: 50
        rmat(1,1) = 0.5
        rmat(1,2) = 0.80902
        rmat(1,3) = -0.30902
        rmat(2,1) = 0.80902
        rmat(2,2) = -0.30902
        rmat(2,3) = 0.5
        rmat(3,1) = 0.30902
        rmat(3,2) = -0.5
        rmat(3,3) = -0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(50, o)

        ! Symmetry operation: 51
        rmat(1,1) = 0.80902
        rmat(1,2) = -0.30902
        rmat(1,3) = -0.5
        rmat(2,1) = 0.30902
        rmat(2,2) = -0.5
        rmat(2,3) = 0.80902
        rmat(3,1) = -0.5
        rmat(3,2) = -0.80902
        rmat(3,3) = -0.30902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(51, o)

        ! Symmetry operation: 52
        rmat(1,1) = 0.5
        rmat(1,2) = -0.80902
        rmat(1,3) = 0.30902
        rmat(2,1) = -0.80902
        rmat(2,2) = -0.30902
        rmat(2,3) = 0.5
        rmat(3,1) = -0.30902
        rmat(3,2) = -0.5
        rmat(3,3) = -0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(52, o)

        ! Symmetry operation: 53
        rmat(1,1) = -0.5
        rmat(1,2) = -0.80902
        rmat(1,3) = 0.30902
        rmat(2,1) = -0.80902
        rmat(2,2) = 0.30902
        rmat(2,3) = -0.5
        rmat(3,1) = 0.30902
        rmat(3,2) = -0.5
        rmat(3,3) = -0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(53, o)

        ! Symmetry operation: 54
        rmat(1,1) = 1
        rmat(1,2) = 0
        rmat(1,3) = 0
        rmat(2,1) = 0
        rmat(2,2) = -1
        rmat(2,3) = 0
        rmat(3,1) = 0
        rmat(3,2) = 0
        rmat(3,3) = -1
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(54, o)

        ! Symmetry operation: 55
        rmat(1,1) = 0.5
        rmat(1,2) = -0.80902
        rmat(1,3) = -0.30902
        rmat(2,1) = -0.80902
        rmat(2,2) = -0.30902
        rmat(2,3) = -0.5
        rmat(3,1) = 0.30902
        rmat(3,2) = 0.5
        rmat(3,3) = -0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(55, o)

        ! Symmetry operation: 56
        rmat(1,1) = 0.30902
        rmat(1,2) = -0.5
        rmat(1,3) = -0.80902
        rmat(2,1) = -0.5
        rmat(2,2) = -0.80902
        rmat(2,3) = 0.30902
        rmat(3,1) = -0.80902
        rmat(3,2) = 0.30902
        rmat(3,3) = -0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(56, o)

        ! Symmetry operation: 57
        rmat(1,1) = -0.5
        rmat(1,2) = -0.80902
        rmat(1,3) = -0.30902
        rmat(2,1) = -0.80902
        rmat(2,2) = 0.30902
        rmat(2,3) = 0.5
        rmat(3,1) = -0.30902
        rmat(3,2) = 0.5
        rmat(3,3) = -0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(57, o)

        ! Symmetry operation: 58
        rmat(1,1) = -0.5
        rmat(1,2) = 0.80902
        rmat(1,3) = -0.30902
        rmat(2,1) = 0.80902
        rmat(2,2) = 0.30902
        rmat(2,3) = -0.5
        rmat(3,1) = -0.30902
        rmat(3,2) = -0.5
        rmat(3,3) = -0.80902
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(58, o)

        ! Symmetry operation: 59
        rmat(1,1) = 0.30902
        rmat(1,2) = 0.5
        rmat(1,3) = -0.80902
        rmat(2,1) = 0.5
        rmat(2,2) = -0.80902
        rmat(2,3) = -0.30902
        rmat(3,1) = -0.80902
        rmat(3,2) = -0.30902
        rmat(3,3) = -0.5
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(59, o)

        ! Symmetry operation: 60
        rmat(1,1) = -1
        rmat(1,2) = 0
        rmat(1,3) = 0
        rmat(2,1) = 0
        rmat(2,2) = 1
        rmat(2,3) = 0
        rmat(3,1) = 0
        rmat(3,2) = 0
        rmat(3,3) = -1
        call o%ori_from_rotmat(rmat, is_ptcl=.false.)
        call self%e_sym%set_ori(60, o)
    end subroutine make_i_relion

    !>  \brief Sets the array of subgroups (character identifier) including itself
    subroutine set_subgrps( self )
        class(sym), intent(inout) :: self
        real                      :: symorders(self%n)
        integer                   :: i, cnt, inds(self%n)
        character(len=1)          :: pgrp
        character(len=3)          :: pgrps(self%n)
        pgrp = self%pgrp(1:1)
        cnt  = 0
        if( pgrp.eq.'c' .or. pgrp.eq.'d' )then
            symorders = 0
            inds = (/(i,i=1,self%n)/)
            if( pgrp.eq.'c' )then
                if( is_even(self%n) )then
                    call getevensym('c', self%n)
                else
                    cnt            = cnt+1
                    pgrps(cnt)     = self%pgrp
                    symorders(cnt) = real(self%n)
                endif
                if(self%n > 1)then
                    cnt            = cnt+1
                    pgrps(cnt)     = 'c1'
                    symorders(cnt) = 1.
                endif
            else if( pgrp.eq.'d' )then
                if( is_even(self%n/2) )then
                    call getevensym('c', self%n/2)
                    call getevensym('d', self%n/2)
                else
                    cnt            = cnt + 1
                    pgrps(cnt)     = self%pgrp
                    symorders(cnt) = real(self%n)
                    cnt            = cnt + 1
                    pgrps(cnt)     = fmtsymstr('c', self%n/2)
                    symorders(cnt) = real(self%n/2)
                    cnt            = cnt + 1
                    pgrps(cnt)     = 'c2'
                    symorders(cnt) = 2.
                endif
                if(self%n > 1)then
                    cnt            = cnt + 1
                    pgrps(cnt)     = 'c1'
                    symorders(cnt) = 1.
                endif
            endif
            call hpsort(symorders, inds)
            allocate( self%subgrps(cnt) )
            do i=1,cnt
                self%subgrps(i) = pgrps(inds(self%n-i+1))
            enddo
        else
            if( pgrp.eq.'t' )then
                cnt = 5
                allocate(self%subgrps(cnt))
                self%subgrps(5) = 'c1'
                self%subgrps(4) = 'c2'
                self%subgrps(3) = 'c3'
                self%subgrps(2) = 'd2'
                self%subgrps(1) = 't'
            else if( pgrp.eq.'o' )then
                cnt = 9
                allocate(self%subgrps(cnt))
                self%subgrps(9) = 'c1'
                self%subgrps(8) = 'c2'
                self%subgrps(7) = 'c3'
                self%subgrps(6) = 'c4'
                self%subgrps(5) = 'd2'
                self%subgrps(4) = 'd3'
                self%subgrps(3) = 'd4'
                self%subgrps(2) = 't'
                self%subgrps(1) = 'o'
            else if( pgrp.eq.'i' )then
                cnt = 9
                allocate(self%subgrps(cnt))
                self%subgrps(9) = 'c1'
                self%subgrps(8) = 'c2'
                self%subgrps(7) = 'c3'
                self%subgrps(6) = 'c5'
                self%subgrps(5) = 'd2'
                self%subgrps(4) = 'd3'
                self%subgrps(3) = 'd5'
                self%subgrps(2) = 't'
                self%subgrps(1) = 'i'
            endif
        endif

        contains

            !> get even symmetry
            subroutine getevensym( cstr, o )
                integer,intent(in)           :: o     !< order
                character(len=1), intent(in) :: cstr  !< sym type
                integer          :: i
                do i=2,o
                    if( (mod(o,i) .eq. 0) )then
                        cnt = cnt + 1
                        pgrps(cnt) = fmtsymstr(cstr, i)
                        if(cstr .eq. 'c')then
                            symorders(cnt) = real(i)
                        else
                            symorders(cnt) = real(2*i)
                        endif
                    endif
                enddo
            end subroutine getevensym

            !> format symmetry string, symmetry type, one char options  c, d, t or i
            function fmtsymstr( symtype, iord ) result( ostr )
                character(len=1), intent(in)  :: symtype
                integer, intent(in)           :: iord    !< order
                character(len=2)              :: ord     !< temp
                character(len=3)              :: ostr    !< formatted output
                write(ord,'(I2)') iord
                write(ostr,'(A1,A2)') symtype, adjustl(ord)
            end function fmtsymstr

    end subroutine set_subgrps

    !>  \brief  is a destructor
    subroutine kill( self )
        class(sym) :: self
        if( allocated(self%subgrps) )deallocate( self%subgrps )
        call self%e_sym%kill
    end subroutine kill

    ! PRIVATE STUFF

    !>  \brief  from SPIDER, creates a rotation matrix around either x or z assumes
    !!          rotation ccw looking towards 0 from +axis accepts 2 arguments,
    !!          0=x or 1=z (first reg.)and rot angle in deg.(second reg.)
    function matcreate( inxorz, indegr ) result( newmat )
        integer, intent(in)          :: inxorz  !< input axis
        double precision, intent(in) :: indegr  !< input rotation angle
        double precision             :: newmat(3,3) !< output rot. matrix
        double precision             :: inrad
        ! create blank rotation matrix
        newmat = 0.d0
        ! change input of degrees to radians for fortran
        inrad = indegr*(dpi/180.d0)
        ! for x rot matrix, place 1 and add cos&sin values
        if( inxorz .eq. 0 )then
            newmat(1,1) = 1.d0
            newmat(2,2) = cos(inrad)
            newmat(3,3) = cos(inrad)
            newmat(3,2) = sin(inrad)
            newmat(2,3) = -sin(inrad)
        ! for z rot matrix, do as above, but for z
        elseif( inxorz .eq. 1 )then
            newmat(3,3) = 1.d0
            newmat(1,1) = cos(inrad)
            newmat(2,2) = cos(inrad)
            newmat(2,1) = sin(inrad)
            newmat(1,2) = -sin(inrad)
        else
            THROW_HARD('Unsupported matrix spec; matcreate')
        endif
    end function matcreate

    !> \brief  from SPIDER, used to calculate the angles SPIDER expects from rot. matrix.
    !!  rotmat assumes sin(theta) is positive (0-180 deg),
    !! \return euls(3) is returned in the SPIDER convention psi, theta, phi
    function matextract( rotmat ) result( euls )
        double precision, intent(inout) :: rotmat(3,3) !<  rotmat rotational matrix
        double precision :: euls(3),radtha,sintha
        double precision :: radone,radtwo,dt
        ! calculate euls(2) from lower/right corner
        ! radtha = acos(rotmat(3,3))
        dt      = max(-1.0d0,min(1.0d0,rotmat(3,3)))
        radtha  = dacos(dt)
        euls(2) = radtha*(180.d0/dpi)
        sintha  = sin(radtha)
        ! close enough test, set corner to -1
        if(abs(1.d0-(rotmat(3,3)/(-1.d0))).lt.(1.e-6))then
            rotmat(3,3) = -1.d0
        endif
        ! close enough test, set corner to 1
        if (abs(1.d0-(rotmat(3,3)/(1.d0))).lt.(1.e-6))then
            rotmat(3,3) = 1.d0
        endif
        ! special case of euls(2) rotation/ y rotaion = 180 or 0
        ! if we do this, need only one z rotation
        if( is_equal(rotmat(3,3),1.d0) )then
            euls(1) = 0.d0
            ! find euls(3), if sin=-, switch sign all in radians
            dt     = max(-1.0d0, min(1.0d0,rotmat(1,1)))
            radone = dacos(dt)
            radtwo = rotmat(2,1)
            if(radtwo < 0.d0)then
                euls(3) = 2.d0*dpi-radone
                euls(3) = euls(3)*(180.d0/dpi)
            else
                euls(3) = radone*(180.d0/dpi)
            endif
        else
            ! normal case of all three rotations
            ! find euls(1), if sin(euls(1)) =- then switch around
            dt     = -rotmat(3,1)/sintha
            dt     = max(-1.0d0,min(1.0d0,dt))
            radone = dacos(dt)
            radtwo = (rotmat(3,2)/sintha)
            if(radtwo < 0.d0)then
                euls(1) = 2.d0*dpi-radone
                euls(1) = euls(1)*(180.d0/dpi)
            else
                euls(1) = radone*(180.d0/dpi)
            endif
            ! find euls(3), similar to before
            dt     = rotmat(1,3)/sintha
            dt     = max(-1.0d0,min(1.0d0,dt))
            radone = dacos(dt)
            radtwo = rotmat(2,3)/sintha
            if(radtwo < 0.d0)then
                euls(3) = 2.d0*dpi-radone
                euls(3) = euls(3)*(180.d0/dpi)
            else
                euls(3) = radone*(180.d0/dpi)
            endif
        endif
        ! catch to change 360 euls(3) to 0
        if(abs(1.d0-(euls(3)/(360.d0))).lt.(1.e-2))then
            euls(3) = 0.d0
        endif
        ! catch to change really small euls(1) to 0, for oct.
        if(abs(1.d0-((euls(1)+1.d0)/1.d0)).lt.(1.e-4))then
            euls(1) = 0.d0
        endif
        ! catch to change really small euls(3) to 0, for oct.
        if(abs(1.d0-((euls(3)+1.)/1.d0)).lt.(1.e-4))then
            euls(3) = 0.d0
        endif
    end function matextract

    subroutine sym_tester(pgrp)
        character(len=*), intent(in) :: pgrp
        type(sym)  :: se, tmp
        type(ori)  :: o, oj, north_pole
        type(oris) :: os
        integer    :: i,j,isym,n
        logical    :: found
        write(logfhandle,'(A,A)')'>>>'
        write(logfhandle,'(A,A)')'>>> POINT GROUP: ', trim(pgrp)
        call se%new(pgrp)
        write(logfhandle,'(A)')'>>> SYMMETRY SUB-GROUPS'
        do isym=1, se%get_nsubgrp()
            tmp = se%get_subgrp(isym)
            write(logfhandle,'(I3,1X,A3)') isym, tmp%get_pgrp()
        enddo
        write(logfhandle,'(A,2F8.3)')'>>> ANGULAR RANGE PHI  :', se%eullims_nomirr(1,:)
        write(logfhandle,'(A,2F8.3)')'>>> ANGULAR RANGE THETA:', se%eullims_nomirr(2,:)
        write(logfhandle,'(A,2F8.3)')'>>> ANGULAR RANGE INCL MIRROR PHI  :', se%eullims(1,:)
        write(logfhandle,'(A,2F8.3)')'>>> ANGULAR RANGE INCL MIRROR THETA:', se%eullims(2,:)
        write(logfhandle,'(A)')'>>> SPIRAL'
        call os%new(1000, is_ptcl=.false.)
        call se%build_refspiral(os)
        call os%write(pgrp//'.txt')
        call os%write2bild(pgrp//'.bild')
        ! redundancy
        n = 0
        call oj%new(is_ptcl=.false.)
        do i=1,os%get_noris()-1
            call os%get_ori(i, o)
            do j=i+1,os%get_noris()
                call os%get_ori(j, oj)
                if( rad2deg(o.euldist.oj) < 0.001 )then
                    n=n+1
                    write(logfhandle,*)i,j,o%get_euler(),oj%get_euler()
                endif
            enddo
        enddo
        if(n==0)then
            write(logfhandle,'(A)')'>>> SPIRAL REDUNDANCY PASSED'
        else
            write(logfhandle,'(A)')'>>> SPIRAL REDUNDANCY FAILED'
            write(logfhandle,*)n
        endif
        ! north pole
        call north_pole%new(is_ptcl=.false.)
        call north_pole%set_euler([0.,0.,0.])
        found = .false.
        do i=1,os%get_noris()
            call os%get_ori(i, o)
            if( (o.euldist.north_pole) < 0.01 )then
                found = .true.
                exit
            endif
        enddo
        if(found)then
            write(logfhandle,'(A)')'>>> NORTH POLE PASSED'
        else
            write(logfhandle,'(A)')'>>> NORTH POLE FAILED'
        endif
        write(logfhandle,'(A)')'>>> SYMMETRY OPERATORS EULER ANGLES:'
        do isym=1, se%get_nsym()
            call se%e_sym%print_(isym)
        enddo
        ! cleanup
        call se%kill
        call os%kill
        call tmp%kill
        call o%kill
        call oj%kill
    end subroutine sym_tester

end module simple_sym
