! defines protein point-group symmetries

module simple_sym
#include "simple_lib.f08"
use simple_oris,   only: oris
use simple_ori,   only: ori
implicit none

public :: sym
private

type sym
    private
    type(oris)                    :: e_sym                 !< symmetry eulers
    character(len=3), allocatable :: subgrps(:)            !< subgroups
    real                          :: eullims(3,2)= 0.      !< euler angles limits (asymetric unit)
    integer                       :: n=1                   !< nr of symmetry ops
    integer                       :: ncsym=1               !< num of order C sym
    integer                       :: ndsym=1               !< num of order D sym
    integer                       :: t_or_o=0              !< tetahedral or octahedral symmetry
    character(len=3)              :: pgrp='c1'             !< point-group symmetry
    logical                       :: c_or_d=.false.        !< c- or d-symmetry
  contains
    procedure          :: new
    procedure          :: srchrange
    procedure          :: which
    procedure          :: get_nsym
    procedure          :: get_pgrp
    procedure          :: apply
    procedure          :: apply2all
    procedure          :: rot_to_asym
    procedure          :: rotall_to_asym
    procedure          :: sym_dists
    procedure          :: get_symori
    procedure          :: get_nsubgrp
    procedure          :: get_subgrp
    procedure          :: get_all_subgrps
    procedure          :: within_asymunit
    procedure          :: apply_sym_with_shift
    procedure          :: nearest_proj_neighbors
    procedure          :: nearest_sym_neighbors
    procedure          :: symrandomize
    procedure, private :: build_srchrange
    procedure, private :: make_c_and_d
    procedure, private :: make_t
    procedure, private :: make_o
    procedure, private :: make_i
    procedure, private :: set_subgrps
    procedure, private :: get_all_cd_subgrps
    procedure :: kill
end type sym

interface sym
    module procedure constructor
end interface sym

integer, parameter          :: ntet=12 ! number of tetahedral symmetry operations
integer, parameter          :: noct=24 ! number of octahedral symmetry operations
integer, parameter          :: nico=60 ! number of icosahedral symmetry operations
double precision, parameter :: delta2 = 180.d0
double precision, parameter :: delta3 = 120.d0
double precision, parameter :: delta5 = 72.d0
double precision, parameter :: alpha = 58.282524d0
double precision, parameter :: beta  = 20.905157d0
double precision, parameter :: gamma = 54.735611d0

contains

    !>  \brief  is a constructor
    function constructor( pgrp) result( self )
        character(len=*), intent(in) :: pgrp !< sym group string
        type(sym)                    :: self
        call self%new(pgrp)
    end function constructor

    !>  \brief  is a constructor
    subroutine new( self, pgrp )
        class(sym), intent(inout)    :: self !< this instance
        character(len=*), intent(in) :: pgrp !< sym group string
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
            call self%e_sym%new(self%n)
            call self%make_c_and_d
        else if(pgrp(1:1).eq.'d' .or. pgrp(1:1).eq.'D')then
            if( self%pgrp(1:1).eq.'D' ) self%pgrp(1:1) = 'd'
            self%c_or_d = .true.
            self%ndsym = 2
            read(pgrp(2:),'(I2)') self%ncsym
            self%n = self%ncsym*self%ndsym
            call self%e_sym%new(self%n)
            call self%make_c_and_d
        else if( pgrp(1:1).eq.'t' .or. pgrp(1:1).eq.'T' )then
            if( self%pgrp(1:1).eq.'T' ) self%pgrp(1:1) = 't'
            self%t_or_o = 1
            self%n = ntet
            call self%e_sym%new(self%n)
            call self%make_t
        else if( pgrp(1:1).eq.'o' .or. pgrp(1:1).eq.'O' )then
            if( self%pgrp(1:1).eq.'O' ) self%pgrp(1:1) = 'o'
            self%t_or_o = 3
            self%n = noct
            call self%e_sym%new(self%n)
            call self%make_o
        else if( pgrp(1:1).eq.'i' .or. pgrp(1:1).eq.'I' )then
            if( self%pgrp(1:1).eq.'I' ) self%pgrp(1:1) = 'i'
            self%n = nico
            call self%e_sym%new(self%n)
            call self%make_i
        else
            write(*,'(a)') 'symmetry not supported; new; simple_sym', pgrp
            stop
        endif
        call self%e_sym%swape1e3
        call self%set_subgrps
        self%eullims = self%build_srchrange()
    end subroutine new

    !>  \brief  builds the search range for the point-group
    function build_srchrange( self ) result( eullims )
        class(sym), intent(inout) :: self !< this instance
        real                      :: eullims(3,2) !< 3-axis search range
        eullims(:,1) = 0.
        eullims(:,2) = 360.
        eullims(2,2) = 180.
        if( self%pgrp(1:1).eq.'c' .and. self%ncsym > 1 )then
            eullims(1,2) = 360./real(self%ncsym)
        else if( self%pgrp(1:1).eq.'d' .and. self%ncsym > 1 )then
            eullims(1,2) = 360./real(self%ncsym)
            eullims(2,2) = 90.
        else if( self%pgrp(1:1).eq.'t' )then
            eullims(1,2) = 180.
            eullims(2,2) = 54.7
        else if( self%pgrp(1:1).eq.'o' )then
            eullims(1,2) = 90.
            eullims(2,2) = 54.7
        else if( self%pgrp(1:1).eq.'i' )then
            eullims(1,2) = 180.
            eullims(2,2) = 31.7
        endif
    end function build_srchrange

    !>  \brief  returns the search range for the point-group
    function srchrange( self ) result( eullims )
        class(sym), intent(inout) :: self !< this instance
        real                      :: eullims(3,2) !< 3-axis search range
        eullims = self%eullims
    end function srchrange

    !>  \brief  to check which point-group symmetry
    pure function which( self ) result( pgrp )
        class(sym), intent(in) :: self !< this instance
        character(len=3) :: pgrp        !< sym group string
        pgrp = self%pgrp
    end function which

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
            write(*,*)'Index out of bonds on simple_sym; get_subgroup'
            stop
        endif
        symobj = sym(self%subgrps(i))
    end function get_subgrp

    !>  get_all_cd_subgrps
    subroutine get_all_cd_subgrps( self, subgrps)
        use simple_math,                   only: is_even
        class(sym), intent(inout)             :: self
        integer                               :: i, cnt
        character(len=1)                      :: pgrp
        character(len=3),allocatable          :: pgrps(:), subgrps(:)
        allocate( pgrps(self%n), stat=alloc_stat )
        allocchk( 'get_all_cd_subgrps; simple_sym; 1')
        pgrp = self%pgrp(1:1)
        cnt  = 0
        if( pgrp=='c' )then
            do i=2,self%n
                cnt = cnt+1
                pgrps(cnt) = fmtsymstr(pgrp, i)
            enddo
        else if( pgrp=='d' )then
            do i=2,self%n/2
                cnt = cnt+1
                pgrps(cnt) = fmtsymstr('c', i)
                cnt = cnt+1
                pgrps(cnt) = fmtsymstr(pgrp, i)
            enddo
        else
            write(*,*)'Unsupported point-group; simple_sym; get_all_cd_subgrps'
            stop
        endif
        allocate( subgrps(cnt), stat=alloc_stat )
        allocchk( 'get_all_cd_subgrps; simple_sym; 2')
        do i=1,cnt
            subgrps(i) = pgrps(i)
        enddo
        deallocate(pgrps)

        contains

            function fmtsymstr( symtype, iord )result( ostr )
                integer, intent(in)           :: iord     !< order index
                character(len=1), intent(in)  :: symtype  !< symmetry
                character(len=2)              :: ord
                character(len=3)              :: ostr     !< formatted output string
                write(ord,'(I2)') iord
                write(ostr,'(A1,A2)') symtype, adjustl(ord)
            end function fmtsymstr

    end subroutine get_all_cd_subgrps

    !>  \brief  Returns array of all symmetry subgroups in c and d
    function get_all_subgrps( self )result( subgrps )
        class(sym), intent(inout)             :: self
        character(len=3),allocatable          :: subgrps(:)
        character(len=1)                      :: pgrp
        pgrp = self%pgrp(1:1)
        if( (pgrp=='c').or.(pgrp=='d') )then
            call self%get_all_cd_subgrps( subgrps )
        else if( pgrp=='t' )then
            allocate( subgrps(1), stat=alloc_stat )
            allocchk( 'get_all_subgrps; simple_sym; 1')
            subgrps(1)  = 't'
        else if( pgrp=='o' )then
            allocate( subgrps(2), stat=alloc_stat )
            allocchk( 'get_all_subgrps; simple_sym; 2')
            subgrps(1)  = 't'
            subgrps(2)  = 'o'
        else if( pgrp=='i' )then
            allocate( subgrps(3), stat=alloc_stat )
            allocchk( 'get_all_subgrps; simple_sym; 3')
            subgrps(1)  = 't'
            subgrps(2)  = 'o'
            subgrps(3)  = 'i'
        endif
    end function get_all_subgrps

    !>  \brief  is a symmetry adaptor
    function apply( self, e_in, symop ) result( e_sym )
        class(sym), intent(inout) :: self
        class(ori), intent(in)    :: e_in  !< orientation object
        integer,    intent(in)    :: symop !< index of sym operation
        type(ori)                 :: e_sym, e_symop, e_tmp
        e_sym   = e_in ! transfer of parameters
        e_symop = self%e_sym%get_ori(symop)
        e_tmp   = e_symop.compose.e_in
        call e_sym%set_euler(e_tmp%get_euler())
    end function apply

    !>  \brief  rotates any orientation to the asymmetric unit
    subroutine rot_to_asym( self, osym )
        class(sym), intent(inout) :: self
        class(ori), intent(inout) :: osym !< orientation
        type(ori) :: oasym
        integer   :: nsym
        if( self%within_asymunit(osym) )then
            ! already in asymetric unit
        else
            do nsym=2,self%n     ! nsym=1 is the identity operator
                oasym = self%apply(osym, nsym)
                if( self%within_asymunit(oasym) )exit
            enddo
            osym = oasym
        endif
    end subroutine rot_to_asym

    !>  \brief  rotates orientations to the asymmetric unit
    subroutine rotall_to_asym( self, osyms )
        class(sym),  intent(inout) :: self     
        class(oris), intent(inout) :: osyms     !< orientations
        type(ori) :: o
        integer   :: i
        do i = 1, osyms%get_noris()
            o = osyms%get_ori(i)
            call self%rot_to_asym(o)
            call osyms%set_ori(i, o)
        enddo
    end subroutine rotall_to_asym

    !>  \brief  determines euler distance and corresponding symmetrized
    !>          orientation between two orientations
    subroutine sym_dists( self, oref, oasym, euldist, inplrotdist )
        use simple_math, only: rad2deg
        class(sym), intent(inout) :: self
        class(ori), intent(in)    :: oref   !< is the untouched reference
        class(ori), intent(in)    :: oasym  !< is outputted as symmetrized
        real,       intent(out)   :: euldist !< Euler distance
        real,       intent(out)   :: inplrotdist !< in-plane rotational distance
        type(ori) :: o, osym
        real      :: dist
        integer   :: isym
        euldist     = oasym.euldist.oref
        inplrotdist = oasym.inplrotdist.oref
        if( self%n == 1 )then
            ! C1: nothing to do
        else
            osym = oasym
            do isym = 2, self%n
                o    = self%apply(oasym, isym)
                dist = o.euldist.oref
                if( dist < euldist )then
                    euldist = dist
                    osym = o
                endif
            enddo
            inplrotdist = osym.inplrotdist.oref
        endif
        euldist     = rad2deg( euldist )
        inplrotdist = rad2deg( inplrotdist )
    end subroutine sym_dists

    !>  \brief  is a getter 
    function get_symori( self, symop ) result( e_sym )
        class(sym), intent(inout) :: self
        integer, intent(in)       :: symop !< symmetry orientation
        type(ori) :: e_sym
        e_sym = self%e_sym%get_ori(symop)
    end function get_symori

    !>  \brief  is a symmetry adaptor
    subroutine apply2all( self, e_in )
        class(sym),  intent(inout) :: self
        class(oris), intent(inout) :: e_in !< symmetry orientations
        type(ori)                  :: orientation
        integer                    :: j, cnt
        cnt = 0
        do j=1,e_in%get_noris()
            cnt = cnt+1
            orientation = e_in%get_ori(j)
            call e_in%set_ori(j, self%apply(orientation, cnt))
            if( cnt == self%n ) cnt = 0
        end do
    end subroutine apply2all

    !>  \brief  whether or not an orientation falls within the asymetric unit
    function within_asymunit( self, o )result( is_within )
        class(sym), intent(inout) :: self
        class(ori), intent(in)    :: o !< symmetry orientation
        logical :: is_within
        real    :: euls(3)
        euls = o%get_euler()
        is_within = .false.
        if( euls(1)<self%eullims(1,1) )return
        if( euls(1)>=self%eullims(1,2) )return
        if( euls(2)<self%eullims(2,1) )return
        if( euls(2)>=self%eullims(2,2) )return
        if( euls(3)<self%eullims(3,1) )return
        if( euls(3)>=self%eullims(3,2) )return
        is_within = .true.
    end function within_asymunit
    
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
                call os%rot(i,symaxis_ori)
                o = os%get_ori(i)
                call self%rot_to_asym(o)
                call os%set_ori(i, o)
            end do
        else
            call os%map3dshift22d( shvec )
            call os%rot(symaxis_ori)
            call self%rotall_to_asym(os)
        endif
    end subroutine apply_sym_with_shift

    !>  \brief  is for retrieving nearest neighbors in symmetric cases 
    subroutine nearest_proj_neighbors( self, os_asym_unit, k, nnmat )
        use simple_math, only: hpsort
        class(sym),           intent(inout) :: self
        type(oris),           intent(inout) :: os_asym_unit !< sampled orientations from assymetric unit, eg from spiral with symmetry
        integer,              intent(inout) :: k
        integer, allocatable, intent(inout) :: nnmat(:,:) !< nearest-neighbour matrix
        real,    allocatable :: dists(:)
        integer, allocatable :: inds(:)
        type(ori)  :: oasym, osym, oj
        integer    :: i, j, n_os, isym
        if( trim(self%pgrp).eq.'c1' )then
            call os_asym_unit%nearest_proj_neighbors(k, nnmat)
        else
            n_os = os_asym_unit%get_noris()
            allocate( nnmat(n_os,k), dists(n_os), inds(n_os), stat=alloc_stat )
            allocchk("In: nearest_proj_neighbors; simple_sym")
            do i = 1, n_os
                dists = pi
                oasym = os_asym_unit%get_ori(i)
                do isym = 1, self%n
                    if(isym == 1)then
                        osym = oasym
                    else
                        osym = self%apply(oasym, isym)
                    endif
                    do j = 1, n_os
                        if( j == i )then
                            dists(j) = 0.
                        else
                            oj = os_asym_unit%get_ori(j)
                            dists(j) = min(dists(j), osym.euldist.oj)
                        endif
                    enddo
                enddo
                inds = (/(j,j=1,n_os)/)
                call hpsort(n_os, dists, inds)
                nnmat(i,:) = inds(:k)
            enddo
            deallocate(inds, dists)
        endif
    end subroutine nearest_proj_neighbors

    !>  \brief  is for retrieving nearest symmetry neighbours in an assymetric set of projection directions
    subroutine nearest_sym_neighbors( self, asym_os, nnmat )
        use simple_math, only: hpsort
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
            allocate( nnmat(n_os,nsym), dists(n_os), stat=alloc_stat )
            allocchk("In: nearest_proj_neighbors; simple_sym")
            do i = 1, n_os
                oasym = asym_os%get_ori(i)
                do isym = 1, nsym
                    if( isym == 1 )then
                        osym = oasym
                    else
                        osym = self%apply(oasym, isym)
                    endif
                    do j = 1, n_os
                        oj = asym_os%get_ori(j)
                        dists(j) = osym.euldist.oj
                    enddo
                    loc = minloc(dists)
                    nnmat(i,isym) = loc(1)
                enddo
            enddo
            deallocate(dists)
        endif
    end subroutine nearest_sym_neighbors

    !>  \brief  is for randomizing over the symmetry operations of the point-group
    subroutine symrandomize( self, os_asym_unit )
        class(sym), intent(inout) :: self
        type(oris), intent(inout) :: os_asym_unit !< sampled orientations from assymetric unit, eg from spiral with symmetry
        integer   :: n_os, nsym, symop, i
        type(ori) :: oasym, osym
        if( trim(self%pgrp).eq.'c1' )then
            stop 'nothing to do, pgrp=c1; simple_sym :: symrandomize'
        endif
        n_os = os_asym_unit%get_noris()
        nsym = self%n
        do i = 1, n_os
            oasym = os_asym_unit%get_ori(i)
            symop = irnd_uni(nsym)
            if( symop > 1 )then
                osym = self%apply(oasym, symop)
                call os_asym_unit%set_ori(i,osym)
            endif
        end do
    end subroutine symrandomize

    !>  \brief  SPIDER code for making c and d symmetries
    subroutine make_c_and_d( self )
        class(sym), intent(inout) :: self
        double precision :: delta, degree
        double precision, dimension(3,3) :: a,b,g
        integer   :: i,j,cnt
        delta = 360.d0/dble(self%ncsym)
        cnt = 0
        do i=0,1
           degree = i*delta2
           a = matcreate(0, degree)
           do j=0,self%ncsym-1
              cnt = cnt+1
                  degree = j*delta
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
            phi = (i-1) * 90.d0
            cnt = cnt + 1
            call self%e_sym%set_euler(cnt,real([0.d0, 0.d0, phi]))
            do j = 1,4
                psi = (j-1) * 90.d0
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
            phi   = i*delta3
            call self%e_sym%set_euler(cnt,real([psi,theta,phi]))
            psi = phi
            do j=0,2
                cnt = cnt+1
                phi = 60.d0+j*delta3
                theta = degree
                call self%e_sym%set_euler(cnt,real([psi,theta,phi]))
            end do
        end do
    end subroutine make_t

    !>  \brief  SPIDER code for making icosahedral symmetry
    subroutine make_i( self )
        class(sym), intent(inout) :: self
        double precision :: deltan, psi, theta, phi
        integer   :: i, j, cnt
        deltan = 36.0
        cnt = 0
        do i=0,1
            do j=0,4
                cnt = cnt+1
                psi = 0.d0
                theta = i*delta2
                phi = j*delta5
                call self%e_sym%set_euler(cnt, real([psi,theta,phi]))
            end do
        end do
        theta = 63.4349d0
        do i=0,4
            do j=0,4
                cnt = cnt+1
                psi = i*delta5
                phi = j*delta5+deltan
                call self%e_sym%set_euler(cnt, real([psi,theta,phi]))
            end do
        end do
        theta = 116.5651
        do i=0,4
            do j=0,4
                cnt = cnt+1
                psi = i*delta5+deltan
                phi = j*delta5
                call self%e_sym%set_euler(cnt, real([psi,theta,phi]))
            end do
        end do
    end subroutine make_i

    !>  \brief Sets the array of subgroups (character identifier) including itself
    subroutine set_subgrps( self )
        use simple_math, only: is_even
        class(sym), intent(inout)     :: self
        integer                       :: i, cnt
        character(len=1)              :: pgrp
        character(len=3), allocatable :: pgrps(:)
        allocate( pgrps(self%n), stat=alloc_stat )
        pgrp = self%pgrp(1:1)
        cnt  = 0
        if( pgrp.eq.'c' )then
            if( is_even(self%n) )then
                call getevensym('c', self%n)
            else
                cnt        = cnt+1
                pgrps(cnt) = self%pgrp
            endif
        else if( pgrp.eq.'d' )then
            if( is_even(self%n/2) )then
                call getevensym('c', self%n/2)
                call getevensym('d', self%n/2)
            else
                cnt        = cnt+1
                pgrps(cnt) = self%pgrp
                cnt        = cnt+1
                pgrps(cnt) = fmtsymstr('c', self%n/2)
                cnt        = cnt+1
                pgrps(cnt) = 'c2'
            endif
        else if( pgrp.eq.'t' )then
            cnt      = 4
            pgrps(1) = 'c2'
            pgrps(2) = 'c3'
            pgrps(3) = 'd2'
            pgrps(4) = 't'
        else if( pgrp.eq.'o' )then
            cnt      = 8
            pgrps(1) = 'c2'
            pgrps(2) = 'c3'
            pgrps(3) = 'c4'
            pgrps(4) = 'd2'
            pgrps(5) = 'd3'
            pgrps(6) = 'd4'
            pgrps(7) = 't'
            pgrps(8) = 'o'
        else if( pgrp.eq.'i' )then
            cnt      = 8
            pgrps(1) = 'c2'
            pgrps(2) = 'c3'
            pgrps(3) = 'c5'
            pgrps(4) = 'd2'
            pgrps(5) = 'd3'
            pgrps(6) = 'd5'
            pgrps(7) = 't'
            pgrps(8) = 'i'
        endif
        if( allocated(self%subgrps) )deallocate( self%subgrps )
        allocate( self%subgrps(cnt), stat=alloc_stat )
        do i=1,cnt
            self%subgrps(i) = pgrps(i)
        enddo
        deallocate(pgrps)

        contains
            !> get even symmetry
            subroutine getevensym( cstr, o )
                integer,intent(in)           :: o     !< order
                character(len=1), intent(in) :: cstr  !< sym type 
                integer          :: i
                do i=2,o
                    if( (mod( o, i).eq.0) )then
                        cnt = cnt+1
                        pgrps(cnt) = fmtsymstr(cstr, i)
                    endif
                enddo
            end subroutine getevensym
            !> format symmetry string
            !! \param symtype   symmetry type, one char options  c, d, t or i
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
            stop 'Unsupported matrix spec; matcreate; simple_ori'
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
        if(abs(1.-(rotmat(3,3)/(-1.))).lt.(1.e-6))then
            rotmat(3,3) = -1.d0
        endif
        ! close enough test, set corner to 1
        if (abs(1.-(rotmat(3,3)/(1.))).lt.(1.e-6))then
            rotmat(3,3) = 1.d0
        endif
        ! special case of euls(2) rotation/ y rotaion = 180 or 0
        ! if we do this, need only one z rotation
        if((dabs(rotmat(3,3))-1.d0)<DTINY)then
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
        if(abs(1.-(euls(3)/(360.))).lt.(1.e-2))then
            euls(3) = 0.d0
        endif
        ! catch to change really small euls(1) to 0, for oct.
        if(abs(1.-((euls(1)+1.)/1.)).lt.(1.e-4))then
            euls(1) = 0.d0
        endif
        ! catch to change really small euls(3) to 0, for oct.
        if(abs(1.-((euls(3)+1.)/1.)).lt.(1.e-4))then
            euls(3) = 0.d0
        endif
    end function matextract

end module simple_sym
