! an orientation
module simple_ori
!include 'simple_lib.f08'
use simple_defs
use simple_error, only: simple_stop
use simple_hash,  only: hash
use simple_chash, only: chash
use simple_math,  only: myacos, deg2rad, rad2deg, make_transfmat, rotmat2d
use simple_rnd,   only: ran3
use simple_stat,  only: pearsn
implicit none

public :: ori, test_ori
private

!>  orientation parameter stuct and operations
type :: ori
    private
    type(hash)  :: htab              !< hash table for the parameters
    type(chash) :: chtab             !< hash table for the filenames etc.
    logical     :: existence=.false. !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new_ori
    procedure          :: new => new_ori
    procedure          :: ori_from_rotmat
    ! SETTERS
    procedure          :: reject
    procedure          :: assign_ori
    generic            :: assignment(=) => assign_ori
    procedure          :: copy_ori
    procedure          :: copy => copy_ori
    procedure          :: delete_entry
    procedure          :: set_euler
    procedure          :: e1set
    procedure          :: e2set
    procedure          :: e3set
    procedure          :: swape1e3
    procedure          :: set_shift
    procedure, private :: set_1
    procedure, private :: set_2
    generic            :: set => set_1, set_2
    procedure          :: rnd_euler_1
    procedure          :: rnd_euler_2
    generic            :: rnd_euler => rnd_euler_1, rnd_euler_2
    procedure          :: rnd_ori
    procedure          :: rnd_inpl
    procedure          :: rnd_shift
    procedure          :: str2ori
    ! GETTERS
    procedure          :: exists
    procedure          :: get_euler
    procedure          :: e1get
    procedure          :: e2get
    procedure          :: e3get
    procedure          :: get_normal
    procedure          :: get_mat
    procedure          :: get
    procedure          :: get_static
    procedure, private :: getter_1
    procedure, private :: getter_2
    generic            :: getter => getter_1, getter_2
    procedure          :: get_2Dshift
    procedure          :: get_3Dshift
    procedure          :: get_state
    procedure          :: hash_size
    procedure          :: isthere
    procedure          :: ischar
    procedure          :: isstatezero
    procedure          :: ori2str
    procedure          :: ori2strlen_trim
    procedure          :: ori2chash
    ! PRINTING & I/O
    procedure          :: print_mat
    procedure          :: print_ori
    procedure          :: write
    procedure          :: read
    ! CALCULATORS
    procedure          :: round_shifts
    procedure, private :: shift
    procedure, private :: compeuler
    procedure          :: compose3d2d
    procedure          :: map3dshift22d
    procedure          :: mirror3d
    procedure          :: mirror2d
    procedure, private :: geodesic_dist
    procedure, private :: geodesic_dist_scaled
    generic            :: operator(.geod.)   => geodesic_dist
    generic            :: operator(.geodsc.) => geodesic_dist_scaled
    procedure, private :: euldist
    procedure, private :: inpldist
    procedure, private :: inplrotdist
    generic            :: operator(.compose.)     => compeuler
    generic            :: operator(.euldist.)     => euldist
    generic            :: operator(.inpldist.)    => inpldist
    generic            :: operator(.inplrotdist.) => inplrotdist
    ! DESTRUCTORS
    procedure          :: kill_chash
    procedure          :: kill
end type ori

interface ori
    module procedure constructor
end interface

real, parameter :: zvec(3) = [0.,0.,1.]

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    function constructor( ) result( self )
        type(ori) :: self
        call self%new_ori
    end function constructor

    !>  \brief  is a constructor
    subroutine new_ori( self )
        class(ori), intent(inout) :: self
        call self%kill
        self%htab  = hash()
        self%chtab = chash()
        self%existence = .true.
    end subroutine new_ori

    !>  \brief  is a parameterized constructor
    subroutine ori_from_rotmat( self, rotmat )
        class(ori), intent(inout) :: self
        real, intent(in)          :: rotmat(3,3) !< rotation matrix
        real :: euls(3)
        call self%new
        euls = m2euler(rotmat)
        call self%set_euler(euls)
    end subroutine ori_from_rotmat

    ! SETTERS

    !>  \brief  sets parameters for particle rejection
    subroutine reject( self )
        class(ori), intent(inout) :: self
        call self%set_euler([0., 0., 0.])
        call self%set_shift([0., 0.])
        call self%htab%set('state', 0.)
        if( self%isthere('corr') )     call self%htab%set('corr',     -1.)
        if( self%isthere('specscore') )call self%htab%set('specscore', 0.)
        if( self%isthere('eo') )       call self%htab%set('eo', -1.)
    end subroutine reject

    !>  \brief  is a polymorphic assigner
    subroutine assign_ori( self_out, self_in )
        type(ori), intent(in)     :: self_in
        class(ori), intent(inout) :: self_out
        call self_out%copy(self_in)
    end subroutine assign_ori

    !>  \brief  is a polymorphic copier
    subroutine copy_ori( self_out, self_in )
        class(ori), intent(in)    :: self_in
        class(ori), intent(inout) :: self_out
        self_out%htab   = self_in%htab
        self_out%chtab  = self_in%chtab
    end subroutine copy_ori

    !>  \brief  is a setter
    subroutine delete_entry( self, key )
        class(ori),       intent(inout) :: self
        character(len=*), intent(in)    :: key
        if( self%htab%isthere(key) )then
            call self%htab%delete(key)
        else
            call self%chtab%delete(key)
        endif
    end subroutine delete_entry

    !>  \brief  is a setter
    subroutine set_euler( self, euls )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: euls(3) !< Euler angle
        real :: euls_here(3), rmat(3,3)
        rmat = euler2m(euls)
        euls_here = m2euler(rmat)
        call self%htab%set('e1',euls_here(1))
        call self%htab%set('e2',euls_here(2))
        call self%htab%set('e3',euls_here(3))
    end subroutine set_euler

    !>  \brief  is a setter
    subroutine e1set( self, e1 )
        class(ori), intent(inout) :: self
        real, intent(in)          :: e1     !< Euler angle
        call self%htab%set('e1', e1)
    end subroutine e1set

    !>  \brief  is a setter
    subroutine e2set( self, e2 )
        class(ori), intent(inout) :: self
        real, intent(in)          :: e2      !< Euler angle
        call self%htab%set('e2', e2)
    end subroutine e2set

    !>  \brief  is a setter
    subroutine e3set( self, e3 )
        class(ori), intent(inout) :: self
        real, intent(in)          :: e3      !< Euler angle
        call self%htab%set('e3', e3)
    end subroutine e3set

    !>  \brief  is a setter
    subroutine swape1e3( self )
        class(ori), intent(inout) :: self
        real :: e1, e3
        e1 = self%htab%get('e1')
        e3 = self%htab%get('e3')
        call self%htab%set('e1', e3)
        call self%htab%set('e3', e1)
    end subroutine swape1e3

    !>  \brief  is a setter
    subroutine set_shift( self, shvec )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: shvec(2) !< shift vector
        call self%htab%set( 'x', shvec(1) )
        call self%htab%set( 'y', shvec(2) )
    end subroutine set_shift

    !>  \brief  is a setter
    subroutine set_1( self, key, val )
        class(ori),       intent(inout) :: self
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
        select case(key)
            case('e1')
                call self%e1set(val)
            case('e2')
                call self%e2set(val)
            case('e3')
                call self%e3set(val)
            case DEFAULT
                call self%htab%set(key, val)
        end select
    end subroutine set_1

    !>  \brief  is a setter
    subroutine set_2( self, key, val )
        class(ori),       intent(inout) :: self
        character(len=*), intent(in)    :: key, val
        call self%chtab%set(key, val)
    end subroutine set_2

    !>  \brief  for generating a random Euler angle
    subroutine rnd_euler_1( self, eullims )
        class(ori),     intent(inout) :: self
        real, optional, intent(inout) :: eullims(3,2)!< Euler angles
        logical :: found
        real    :: euls(3), rmat(3,3)
        if( present(eullims) )then
            found = .false.
            do while( .not.found )
                call rnd_romat(rmat)
                euls = m2euler(rmat)
                if( euls(1) >= eullims(1,2) )cycle
                if( euls(2) >= eullims(2,2) )cycle
                if( euls(3) >= eullims(3,2) )cycle
                if( euls(1) < eullims(1,1) )cycle
                if( euls(2) < eullims(2,1) )cycle
                if( euls(3) < eullims(3,1) )cycle
                found = .true.
            end do
        else
            call rnd_romat(rmat)
            euls = m2euler(rmat)
        endif
        call self%set_euler(euls)
    end subroutine rnd_euler_1

    !>  \brief  for generating a random Euler angle neighbour to o_prev
    subroutine rnd_euler_2( self, o_prev, athres, proj )
        class(ori),        intent(inout) :: self   !< instance
        class(ori),        intent(in)    :: o_prev !< template ori
        real,              intent(in)    :: athres !< angle threshold in degrees
        logical, optional, intent(in)    :: proj   !< projection only thres or not
        real    :: athres_rad, dist
        logical :: pproj
        pproj = .false.
        if( present(proj) ) pproj = proj
        athres_rad = deg2rad(athres)
        dist = 2.*athres_rad
        do while( dist > athres_rad )
            call self%rnd_euler_1
            if( pproj )then
                dist = self.euldist.o_prev
            else
                dist = self.geodsc.o_prev
            endif
        end do
    end subroutine rnd_euler_2

    !>  \brief  for generating random ori
    subroutine rnd_ori( self, trs, eullims )
        class(ori),     intent(inout) :: self
        real, optional, intent(in)    :: trs         !< threshold
        real, optional, intent(inout) :: eullims(3,2)!< Euler angles
        real :: x, y, z
        if( present(trs) )then
            if( abs(trs) < 1e-3 )then
                x = 0.
                y = 0.
                z = 0.
            else
                x = ran3()*2*trs-trs
                y = ran3()*2*trs-trs
                z = ran3()*2*trs-trs
            endif
            call self%set('x', x)
            call self%set('y', y)
            call self%set('z', z)
        endif
        if( present( eullims ) )then
            call self%rnd_euler( eullims )
        else
            call self%rnd_euler
        endif
    end subroutine rnd_ori

    !>  \brief  for generating random in-plane parameters
    subroutine rnd_inpl( self, trs )
        class(ori), intent(inout)  :: self
        real, intent(in), optional :: trs         !< threshold
        if( present(trs) )call self%rnd_shift(trs)
        call self%e3set(ran3()*359.99)
    end subroutine rnd_inpl

    !>  \brief  for generating random in-plane parameters
    subroutine rnd_shift( self, trs )
        class(ori), intent(inout)  :: self
        real,       intent(in)     :: trs         !< threshold
        real :: x, y
        if( abs(trs) < 1e-3 )then
            x = 0.
            y = 0.
        else
            x = ran3()*2*trs-trs
            y = ran3()*2*trs-trs
        endif
        call self%set('x', x)
        call self%set('y', y)
    end subroutine rnd_shift

    !>  \brief  reads all orientation info (by line) into the hash-tables
    subroutine str2ori( self, line )
        use simple_sauron, only: sauron_line_parser
        class(ori),       intent(inout) :: self
        character(len=*), intent(inout) :: line
        logical :: isthere(3)
        call self%new
        call sauron_line_parser(line, self%htab, self%chtab)
        isthere(1) = self%htab%isthere('e1')
        isthere(2) = self%htab%isthere('e2')
        isthere(3) = self%htab%isthere('e3')
        if( any(isthere) )&
          &call self%set_euler([self%htab%get('e1'),self%htab%get('e2'),self%htab%get('e3')])
    end subroutine str2ori

    ! GETTERS

     !>  \brief  is a getter
    pure function exists( self ) result( t )
        class(ori), intent(in) :: self
        logical :: t
        t = self%existence
    end function exists

    !>  \brief  is a getter
    pure function get_euler( self ) result( euls )
        class(ori), intent(in) :: self
        real :: euls(3)
        euls(1) = self%htab%get('e1')
        euls(2) = self%htab%get('e2')
        euls(3) = self%htab%get('e3')
    end function get_euler

    !>  \brief  is a getter
    pure function e1get( self ) result( e1 )
        class(ori), intent(in) :: self
        real :: e1
        e1 = self%htab%get('e1')
    end function e1get

    !>  \brief  is a getter
    pure function e2get( self ) result( e2 )
        class(ori), intent(in) :: self
        real :: e2
        e2 = self%htab%get('e2')
    end function e2get

    !>  \brief  is a getter
    pure function e3get( self ) result( e3 )
        class(ori), intent(in) :: self
        real :: e3
        e3 = self%htab%get('e3')
    end function e3get

    !>  \brief  is a getter
    pure function get_normal( self ) result( normal )
        class(ori), intent(in) :: self
        real :: normal(3), rmat(3,3), euls(3)
        euls = self%get_euler()
        rmat   = euler2m(euls)
        normal = matmul(zvec, rmat)
    end function get_normal

    !>  \brief  is a getter
    pure function get_mat( self ) result( mat )
        class(ori), intent(in) :: self
        real, dimension(3,3) :: mat
        real :: euls(3)
        euls = self%get_euler()
        mat  = euler2m(euls)
    end function get_mat

    !>  \brief  is a getter
    function get( self, key ) result( val )
        class(ori), intent(inout)    :: self
        character(len=*), intent(in) :: key
        real :: val
        val = self%htab%get(key)
    end function get

    !>  \brief  is a getter with fixed length return string
    function get_static( self, key )result( val )
        class(ori),           intent(inout) :: self
        character(len=*),     intent(in)    :: key
        character(len=STDLEN)               :: val
        val = trim(self%chtab%get_static(key))
    end function get_static

    !>  \brief  is a getter
    subroutine getter_1( self, key, val )
        class(ori),                    intent(inout) :: self
        character(len=*),              intent(in)    :: key
        character(len=:), allocatable, intent(inout) :: val
        if( allocated(val) ) deallocate(val)
        val = self%chtab%get(key)
    end subroutine getter_1

    !>  \brief  is a getter
    subroutine getter_2( self, key, val )
        class(ori),       intent(inout) :: self
        character(len=*), intent(in)    :: key
        real,             intent(inout) :: val
        val = self%htab%get(key)
    end subroutine getter_2

    !>  \brief  is a getter
    function get_2Dshift( self ) result( vec )
        class(ori), intent(inout) :: self
        real :: vec(2)
        vec(1) = self%htab%get('x')
        vec(2) = self%htab%get('y')
    end function get_2Dshift

    !>  \brief  is a getter
    function get_3Dshift( self ) result( vec )
        class(ori), intent(inout) :: self
        real :: vec(3)
        vec(1) = self%htab%get('x')
        vec(2) = self%htab%get('y')
        vec(2) = self%htab%get('z')
    end function get_3Dshift

    !>  \brief  is a getter
    integer function get_state( self )
        class(ori), intent(inout) :: self
        get_state = nint(self%htab%get('state'))
    end function get_state

    !>  \brief  returns size of hash
    function hash_size( self ) result( sz )
        class(ori), intent(in) :: self
        integer :: sz
        sz = self%htab%size_of()
    end function hash_size

    !>  \brief  check for presence of key in the ori hash
    function isthere( self, key ) result( found )
        class(ori),        intent(inout) :: self
        character(len=*),  intent(in)    :: key
        logical :: hash_found, chash_found, found
        hash_found  = self%htab%isthere(key)
        chash_found = self%chtab%isthere(key)
        found = .false.
        if( hash_found .or. chash_found ) found = .true.
    end function isthere

    !>  \brief  test for character key
    function ischar( self, key ) result( is )
        class(ori),        intent(inout) :: self
        character(len=*),  intent(in)    :: key
        logical :: is
        is = self%chtab%isthere(key)
    end function ischar

    !>  \brief  whether state is zero
    logical function isstatezero( self )
        class(ori),       intent(inout) :: self
        isstatezero = (self%get_state() == 0)
    end function isstatezero

    !>  \brief  joins the hashes into a string that represent the ori
    function ori2str( self ) result( str )
        class(ori), intent(inout) :: self
        character(len=:), allocatable :: str, str_chtab, str_htab
        integer :: sz_chash, sz_hash
        sz_chash = self%chtab%size_of()
        sz_hash  = self%htab%size_of()
        if( sz_chash > 0 ) str_chtab = self%chtab%chash2str()
        if( sz_hash  > 0 ) str_htab  = self%htab%hash2str()
        if( sz_chash > 0 .and. sz_hash > 0 )then
            allocate( str, source=str_chtab//' '//str_htab)
        else if( sz_hash > 0 )then
            allocate( str, source=str_htab)
        else if( sz_chash > 0 )then
            allocate( str, source=str_chtab)
        endif
    end function ori2str

    function ori2strlen_trim( self ) result( len )
        class(ori), intent(inout) :: self
        character(len=:), allocatable :: str
        integer :: len
        str = self%ori2str()
        len = len_trim(str)
    end function ori2strlen_trim

    function ori2chash( self ) result( ch )
        use simple_chash, only: chash
        class(ori), intent(in) :: self
        type(chash) :: ch
        ch = self%chtab
        call self%htab%push2chash(ch, as_ints=.true.)
    end function ori2chash

    !<  \brief  to print the rotation matrix
    subroutine print_mat( self )
        class(ori), intent(inout) :: self
        real :: euls(3), rmat(3,3)
        euls = self%get_euler()
        rmat = euler2m(euls)
        write(*,*) rmat(1,1), rmat(1,2), rmat(1,3), &
            &      rmat(2,1), rmat(2,2), rmat(2,3), &
            &      rmat(3,1), rmat(3,2), rmat(3,3)
    end subroutine print_mat

    !>  \brief prints ori data
    subroutine print_ori( self )
        class(ori), intent(inout) :: self
        call self%htab%print()
        call self%chtab%print_key_val_pairs
    end subroutine print_ori

    !>  \brief  writes orientation info
    subroutine write( self, fhandle )
        class(ori), intent(inout) :: self
        integer,    intent(in)    :: fhandle
        character(len=:), allocatable :: str
        str = self%ori2str()
        if( allocated(str) ) write(fhandle,'(a)') str
    end subroutine write

    !>  \brief  reads all orientation info (by line) into the hash-tables
    subroutine read( self, fhandle )
        use simple_sauron, only: sauron_line_parser
        class(ori), intent(inout) :: self
        integer,    intent(in)    :: fhandle
        character(len=2048) :: line
        logical :: isthere(3)
        read(fhandle,fmt='(A)') line
        call sauron_line_parser( line, self%htab, self%chtab )
        isthere(1) = self%htab%isthere('e1')
        isthere(2) = self%htab%isthere('e2')
        isthere(3) = self%htab%isthere('e3')
        if( any(isthere) )&
          &call self%set_euler([self%htab%get('e1'),self%htab%get('e2'),self%htab%get('e3')])
    end subroutine read

    ! CALCULATORS

    !>  \brief  rounds the origin shifts
    subroutine round_shifts( self )
        class(ori), intent(inout) :: self
        call self%set('x', real(nint(self%get('x'))))
        call self%set('y', real(nint(self%get('y'))))
    end subroutine round_shifts

    !>  \brief  corrects the Euler angle bounds
    subroutine shift( self )
        class(ori), intent(inout) :: self
        logical :: doshift
        real    :: euls(3), rmat(3,3)
        doshift = .false.
        euls = self%get_euler()
        if( euls(1) < 0. .or. euls(1) > 360. ) doshift = .true.
        if( euls(2) < 0. .or. euls(2) > 180. ) doshift = .true.
        if( euls(3) < 0. .or. euls(3) > 360. ) doshift = .true.
        if( doshift )then
            rmat = euler2m(euls)
            euls = m2euler(rmat)
            call self%set_euler(euls)
        endif
    end subroutine shift

    !>  \brief  is for composing Euler angles
    !! multiplication of two rotation matrices commute (n>2 not)
    !! the composed rotation matrix is constructed
    !! \param self1,self2 ori class type rotational matrices
    function compeuler( self1, self2 ) result( self_out )
        class(ori), intent(in) :: self1, self2
        type(ori) :: self_out
        real      :: euls(3), euls1(3), euls2(3), rmat(3,3), rmat1(3,3), rmat2(3,3)
        call self_out%new_ori
        euls1 = self1%get_euler()
        euls2 = self2%get_euler()
        rmat1 = euler2m(euls1)
        rmat2 = euler2m(euls2)
        rmat = matmul(rmat2,rmat1)  ! multiplication of two rotation matrices commute (n>2 not)
        ! the composed rotation matrix is constructed
        ! convert to euler
        euls = m2euler(rmat)
        call self_out%set_euler(euls)
    end function compeuler

    !>  combines 3d and 2d oris and flags for ptcl mirroring
    !! \param ori3d,ori2d,o_out ori class type rotational matrices
    subroutine compose3d2d( ori3d, ori2d, o_out )
        class(ori), intent(inout) :: ori3d, ori2d, o_out
        real                      :: ori3dx,ori3dy,x,y,e3,euls(3)
        real                      :: R3d(3,3), R2d(3,3), R(3,3)
        logical                   :: mirr
        ! inpl from 3D
        ori3dx = ori3d%get('x')
        ori3dy = ori3d%get('y')
        euls   = ori3d%get_euler()
        mirr   = .false.
        if( nint(ori3d%get('mirr')) == 1 )then
            mirr    = .true.
            euls(3) = -euls(3)
            ori3dy  = -ori3dy
        endif
        R3d = make_transfmat(euls(3), ori3dx, ori3dy)
        x  = ori2d%get('x')
        y  = ori2d%get('y')
        e3 = ori2d%e3get()
        ! inpls composition and determination
        R2d = make_transfmat(e3, x, y)
        R   = matmul(R2d,R3d)
        call transfmat2inpls(R, e3, x, y)
        ! mirror for southern hemisphere
        if( mirr )then
            call o_out%set('mirr', 1.)
            e3 = -e3
            y  = -y
        endif
        call o_out%set('x',x)
        call o_out%set('y',y)
        euls(3) = e3
        call o_out%set_euler(euls)

        contains

            !>  extracts in-plane parameters from transformation matrix
            subroutine transfmat2inpls( R, psi, tx, ty )
                real,intent(inout) :: psi,tx,ty
                real,intent(in)    :: R(3,3)
                psi = rad2deg( myacos( R(1,1) ))
                if( R(1,2)<0. )psi=360.-psi
                tx  = R(1,3)
                ty  = R(2,3)
            end subroutine transfmat2inpls

    end subroutine compose3d2d

    subroutine map3dshift22d( self, sh3d )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: sh3d(3) !< 3D shift
        real :: old_x,old_y,x,y,cx,cy
        real :: u(3),v(3),shift(3),euls(3)
        real :: phi,cosphi,sinphi,rmat(3,3)
        euls   = self%get_euler()
        phi    = deg2rad(euls(3))
        cosphi = cos( phi )
        sinphi = sin( phi )
        old_x  = self%get('x')
        old_y  = self%get('y')
        ! 3D Shift rotated with respect to projdir
        rmat  = euler2m(euls)
        shift = matmul(rmat,sh3d)
        ! Projection onto xy plane
        shift = shift - dot_product(shift,zvec) * zvec
        ! e3-composed xy plane unit vectors
        u = [  cosphi,sinphi,0. ]
        v = [ -sinphi,cosphi,0. ]
        ! composed shift clockwise components
        cx =  old_x * cosphi + old_y * sinphi + dot_product(u,shift)
        cy = -old_x * sinphi + old_y * cosphi + dot_product(v,shift)
        ! Identification
        x = cx * cosphi - cy * sinphi
        y = cx * sinphi + cy * cosphi
        ! the end
        call self%set('x', x)
        call self%set('y', y)
    end subroutine map3dshift22d

    !>  \brief  generates the opposite hand of an Euler angle
    !!          so that a set of Euler angles transformed by
    !!          this operation changes the handedness of the volume
    subroutine mirror3d( self )
        class(ori), intent(inout) :: self
        real :: mirr_mat(3,3), rmat(3,3), euls(3)
        mirr_mat = 0.
        mirr_mat(1,1) = 1.
        mirr_mat(2,2) = 1.
        mirr_mat(3,3) = -1.
        euls = self%get_euler()
        rmat = euler2m(euls)
        rmat = matmul(mirr_mat,matmul(rmat,mirr_mat))
        ! the mirrored rotation matrix is constructed
        ! convert to euler
        euls = m2euler(rmat)
        call self%set_euler(euls)
    end subroutine mirror3d

    !>  \brief  generates the mirror of the projection
    !!          can be equivalently accomplished by mirror('y')
    !!          the image after projection
    subroutine mirror2d( self )
        class(ori), intent(inout) :: self
        type(ori) :: old
        real      :: euls(3), rmat(3,3), euls_mirr(3)
        old = self
        euls = self%get_euler()
        euls_mirr(1) = euls(1)
        euls_mirr(2) = 180. + euls(2)
        euls_mirr(3) = 180. - euls(3)
        rmat = euler2m(euls_mirr)
        ! the mirrored rotation matrix is constructed
        ! convert to euler
        euls = m2euler(rmat)
        call self%set_euler(euls)
    end subroutine mirror2d

    ! GEODESIC DISTANCE METRIC

    ! Here we seek to define a bi-invariant metric that respect the topology on SO(3), i.e.
    ! we will define ta metric that attemts to find the rotation required to bring
    ! R1 in register with R2. Our goal is to find R such that R1=matmul(R,R2), thus
    ! R=matmul(R1,transpose(R2)).
    !
    ! Formally:
    !
    ! dist: SO(3)xSO(3)->R+
    !
    ! This metric is sensitive to in-plane rotation

    !>  \brief  this metric is measuring the frobenius deviation from the identity matrix .in.[0,2*sqrt(2)]
    !!          Larochelle, P.M., Murray, A.P., Angeles, J., A distance metric for finite
    !!          sets of rigid-body displacement in the polar decomposition. ASME J. Mech. Des.
    !!          129, 883--886 (2007)
    !! \param self1,self2 ori class type rotational matrices
    pure real function geodesic_dist( self1, self2 )
        class(ori), intent(in) :: self1, self2
        real :: Imat(3,3), sumsq, diffmat(3,3)
        real :: euls1(3), euls2(3), rmat1(3,3), rmat2(3,3)
        Imat      = 0.
        Imat(1,1) = 1.
        Imat(2,2) = 1.
        Imat(3,3) = 1.
        euls1 = self1%get_euler()
        euls2 = self2%get_euler()
        rmat1 = euler2m(euls1)
        rmat2 = euler2m(euls2)
        diffmat = Imat - matmul(rmat1,transpose(rmat2))
        sumsq   = sum(diffmat*diffmat)
        if( sumsq > 0.0001 )then
            geodesic_dist = sqrt(sumsq)
        else
            geodesic_dist = 0.
        endif
    end function geodesic_dist

    !>  \brief  this is the same metric as above, scaled to [0,pi]
    !! \param self1,self2 ori class type rotational matrices
    pure real function geodesic_dist_scaled( self1, self2 )
        class(ori), intent(in) :: self1, self2
        real, parameter :: old_max = 2.*sqrt(2.)
        geodesic_dist_scaled = self1%geodesic_dist(self2)*(pi/old_max)
    end function geodesic_dist_scaled

    ! CLASSIC DISTANCE METRICS

    !>  \brief  calculates the distance (in radians) btw two Euler angles
    !! \param self1,self2 ori class type rotational matrices
    pure function euldist( self1, self2 ) result( dist )
        class(ori), intent(in) :: self1, self2
        real :: dist, normal1(3), normal2(3)
        normal1 = self1%get_normal()
        normal2 = self2%get_normal()
        dist = myacos(dot_product(normal1,normal2))
    end function euldist

    !>  \brief  calculates the distance (in radians) btw all df:s
    !! \param self1,self2 ori class type rotaional matrices
    pure function inpldist( self1, self2 ) result( dist )
        class(ori), intent(in) :: self1, self2
        real :: dist
        dist = ((self1.inplrotdist.self2)+(self1.euldist.self2))/2.
    end function inpldist

    !>  \brief  calculates the distance (in radians) btw the in-plane rotations
    !! \param self1,self2 ori class type rotaional matrices
    pure function inplrotdist( self1, self2 ) result( dist )
        class(ori), intent(in) :: self1, self2
        real :: dist, mat(2,2), u(2), x1(2), x2(2)
        u(1) = 0.
        u(2) = 1.
        mat  = rotmat2d(self1%e3get())
        x1   = matmul(u,mat)
        mat  = rotmat2d(self2%e3get())
        x2   = matmul(u,mat)
        dist = myacos(dot_product(x1,x2))
    end function inplrotdist

    ! DESTRUCTORS

    !>  \brief  is a destructor
    subroutine kill_chash( self )
        class(ori), intent(inout) :: self
        call self%chtab%kill
    end subroutine kill_chash

    !>  \brief  is a destructor
    subroutine kill( self )
        class(ori), intent(inout) :: self
        if( self%existence )then
            call self%htab%kill
            call self%chtab%kill
            self%existence = .false.
        endif
    end subroutine kill

    ! PRIVATE STUFF

    !>  \brief  makes a rotation matrix from a Spider format Euler triplet
    pure function euler2m( euls ) result( r )
        real, intent(in)     :: euls(3)
        real, dimension(3,3) :: r1, r2, r3, r, tmp
        r1 = rotmat(euls(1),3) ! rotation around z
        r2 = rotmat(euls(2),2) ! tilt
        r3 = rotmat(euls(3),3) ! rotation around z
        ! order of multiplication is r3r2r1
        tmp = matmul(r3,r2)
        r   = matmul(tmp,r1)
    end function euler2m

    !> \brief  calculates the derivatives of the rotation matrices w.r.t. one Euler angle
    subroutine euler2dm( euls, drmat )
        real, intent(in)     :: euls(3)
        real, intent(out)    :: drmat(3,3,3)
        real, dimension(3,3) :: r1, r2, r3, dr1, dr2, dr3
        r1  =  rotmat(euls(1),3) ! rotation around z
        r2  =  rotmat(euls(2),2) ! tilt
        r3  =  rotmat(euls(3),3) ! rotation around z
        dr1 = drotmat(euls(1),3) ! derivative of r1 w.r.t. e1
        dr2 = drotmat(euls(2),2) ! derivative of r2 w.r.t. e2
        dr3 = drotmat(euls(3),3) ! derivative of r3 w.r.t. e3
        drmat(:,:,1) = matmul(matmul(r3,r2),dr1)
        drmat(:,:,2) = matmul(matmul(r3,dr2),r1)
        drmat(:,:,3) = matmul(matmul(dr3,r2),r1)
    end subroutine euler2dm

    !>  \brief  returns a Spider format Euler triplet, given a rotation matrix
    pure function m2euler( mat ) result ( euls )
        real, intent(in), dimension(3,3) :: mat(3,3)
        real, dimension(3,3):: tmp, anticomp1z, anticomp2y, get_euler3eul
        real :: imagevec(3), absxy, mceul1deg, mceul2deg, phitmp, euls(3)
        ! get_euler image of (0,0,1)
        imagevec = matmul( zvec, mat )
        ! extract eul1 from imagevec:
        absxy = sqrt(imagevec(1)**2+imagevec(2)**2)
        if(absxy < 0.0000001) then
            ! normal parallel to z, phi undefined
            euls(1) = 0.
        else
            phitmp = myacos(imagevec(1)/absxy)
            if (imagevec(2) >= 0.) then
                euls(1) = phitmp
            else
                euls(1) = -phitmp
            endif
        endif
        ! extract the theta
        euls(2) = myacos(imagevec(3))
        ! set_eulerup the rotation matrices to retrieve the resulting psi angle:
        mceul1deg  = rad2deg( -euls(1) )
        mceul2deg  = rad2deg( -euls(2) )
        anticomp1z = rotmat( mceul1deg,3 )
        anticomp2y = rotmat( mceul2deg,2 )
        tmp = matmul( anticomp1z,anticomp2y  )
        get_euler3eul = matmul( mat,tmp )
        ! get_euler3eul is a rotation around z, 1st element is cosangle
        if( get_euler3eul( 1,2 ) > 0. ) then
            euls(3) = myacos(get_euler3eul(1,1))
        else
            euls(3) = -myacos(get_euler3eul(1,1))
        endif
        euls(3)=rad2deg( euls(3) )
        euls(2)=rad2deg( euls(2) )
        euls(1)=rad2deg( euls(1) )
        do while( euls(3) < 0. )
            euls(3)=euls(3)+360.
        end do
        do while( euls(1) < 0. )
            euls(1)=euls(1)+360.
        end do
    end function m2euler

    !>  \brief  returns the rotation matrix for _ang_ degrees of rotation
    !! around x,y or z for _choice_ = _1_,_2_ or _3_
    pure function rotmat( ang, choice ) result( r )
        real, intent(in)           :: ang
        integer, intent(in)        :: choice
        real :: r(3,3)
        real :: ang_in_rad
        ang_in_rad = ang*pi/180.
        if ( choice == 1 ) then
            r( 1,1 ) = 1.
            r( 1,2 ) = 0.
            r( 1,3 ) = 0.
            r( 2,1 ) = 0.
            r( 2,2 ) = cos( ang_in_rad )
            r( 2,3 ) =-sin( ang_in_rad )
            r( 3,1 ) = 0.
            r( 3,2 ) = sin( ang_in_rad )
            r( 3,3 ) = cos( ang_in_rad )
        elseif ( choice == 2 ) then
            r( 1,1 ) = cos( ang_in_rad )
            r( 1,2 ) = 0.
            r( 1,3 ) = -sin( ang_in_rad )
            r( 2,1 ) = 0.
            r( 2,2 ) = 1.
            r( 2,3 ) = 0.
            r( 3,1 ) = sin( ang_in_rad )
            r( 3,2 ) = 0.
            r( 3,3 ) = cos( ang_in_rad )
        elseif ( choice == 3 ) then
            r( 1,1 ) = cos( ang_in_rad )
            r( 1,2 ) = sin( ang_in_rad )
            r( 1,3 ) = 0.
            r( 2,1 ) = -sin( ang_in_rad )
            r( 2,2 ) = cos( ang_in_rad )
            r( 2,3 ) = 0.
            r( 3,1 ) = 0.
            r( 3,2 ) = 0.
            r( 3,3 ) = 1.
            ! beware of the signs:z-rot is really negative
        endif
    end function rotmat

    !>  \brief  returns the derivative of the rotation matrix for _ang_ degrees of rotation
    !! around x,y or z for _choice_ = _1_,_2_ or _3_
    pure function drotmat( ang, choice ) result( r )
        real, intent(in)           :: ang
        integer, intent(in)        :: choice
        real :: r(3,3)
        real :: ang_in_rad
        ang_in_rad = ang*pi/180.
        if ( choice == 1 ) then
            r( 1,1 ) = 0.
            r( 1,2 ) = 0.
            r( 1,3 ) = 0.
            r( 2,1 ) = 0.
            r( 2,2 ) = -sin( ang_in_rad )
            r( 2,3 ) = -cos( ang_in_rad )
            r( 3,1 ) = 0.
            r( 3,2 ) =  cos( ang_in_rad )
            r( 3,3 ) = -sin( ang_in_rad )
        elseif ( choice == 2 ) then
            r( 1,1 ) = -sin( ang_in_rad )
            r( 1,2 ) = 0.
            r( 1,3 ) = -cos( ang_in_rad )
            r( 2,1 ) = 0.
            r( 2,2 ) = 0.
            r( 2,3 ) = 0.
            r( 3,1 ) =  cos( ang_in_rad )
            r( 3,2 ) = 0.
            r( 3,3 ) = -sin( ang_in_rad )
        elseif ( choice == 3 ) then
            r( 1,1 ) = -sin( ang_in_rad )
            r( 1,2 ) =  cos( ang_in_rad )
            r( 1,3 ) = 0.
            r( 2,1 ) = -cos( ang_in_rad )
            r( 2,2 ) = -sin( ang_in_rad )
            r( 2,3 ) = 0.
            r( 3,1 ) = 0.
            r( 3,2 ) = 0.
            r( 3,3 ) = 0.
            ! beware of the signs:z-rot is really negative
        endif
    end function drotmat

    !>  \brief  for generating a random rotation matrix
    subroutine rnd_romat( rmat )
        ! Fast Random Rotation Matrices, Arvo, Graphics Gems III, 1992
        real, intent(out) :: rmat(3,3)
        real :: theta, phi, z, vx, vy, vz
        real :: r, st, ct, sx, sy
        ! init
        theta = ran3()*TWOPI
        phi   = ran3()*TWOPI
        z     = ran3()*2.
        ! V
        r  = sqrt( z )
        vx = r*sin( phi )
        vy = r*cos( phi )
        vz = sqrt( 2.-z )
        ! S=Vt*R; sz=vz
        st = sin( theta )
        ct = cos( theta )
        sx = vx*ct - vy*st
        sy = vx*st + vy*ct
        ! M
        rmat(1,1) = vx * sx - ct
        rmat(1,2) = vx * sy - st
        rmat(1,3) = vx * vz
        rmat(2,1) = vy * sx + st
        rmat(2,2) = vy * sy - ct
        rmat(2,3) = vy * vz
        rmat(3,1) = vz * sx
        rmat(3,2) = vz * sy
        rmat(3,3) = 1. - z
    end subroutine rnd_romat

    ! UNIT TEST

    !>  \brief  is the unit test for the ori class
    subroutine test_ori
    ! this test only tests the Euler part of ori,
    ! the rest is tested in the oris class
        type(ori) :: e1, e2, e3
        real :: euls(3), normal(3), mat(3,3), normal2(3), dist
        logical :: passed
        call e1%new_ori
        call e2%new_ori
        call e3%new_ori
        write(*,'(a)') '**info(simple_ori_unit_test: testing all functionality'
        passed = .false.
        call e1%set_euler([1.,2.,3.])
        euls = e1%get_euler()
        if( abs(euls(1)-1.+euls(2)-2.+euls(3)-3.) < 0.001 )then
            passed = .true.
        endif
        if( .not. passed ) call simple_stop('Euler assignment/getters corrupt!')
        passed = .false.
        call e1%e1set(99.)
        call e1%e2set(98.)
        call e1%e3set(97.)
        euls = e1%get_euler()
        if( abs(euls(1)-99.+euls(2)-98.+euls(3)-97.) < 0.001 ) passed = .true.
        if( .not. passed ) call simple_stop('Euler e-setters corrupt!')
        passed = .false.
        call e1%set_euler([0.,0.,0.])
        normal = e1%get_normal()
        if( abs(normal(1))+abs(normal(2))+abs(normal(3)-1.)<3.*TINY ) passed = .true.
        if( .not. passed ) call simple_stop('Euler normal derivation corrupt!')
        passed = .false.
        mat = e1%get_mat()
        if( abs(mat(1,1)-1.)+abs(mat(2,2)-1.)+abs(mat(3,3)-1.)<3.*TINY )then
            mat(1,1) = 0.
            mat(2,2) = 0.
            mat(3,3) = 0.
            if( abs(sum(mat)) < TINY ) passed = .true.
        endif
        if( .not. passed ) call simple_stop('Euler rotation matrix derivation corrupt!')
        passed = .false.
        call e2%set_euler([20.,20.,20.])
        e3 = e2
        normal = e2%get_normal()
        normal2 = e3%get_normal()
        if( pearsn(normal,normal2) > 0.99 ) passed = .true.
        passed = .false.
        dist = rad2deg(e1.euldist.e2) ! IFORT CANNOT DEAL WITH THE OPERATORS HERE
        if( dist < 20.00001 .and. dist > 19.99999 ) passed = .true.
        passed = .false.
        e3 = e1%compeuler(e2)         ! IFORT CANNOT DEAL WITH THE OPERATORS HERE
        euls = e3%get_euler()
        if( euls(1) < 20.0001 .and. euls(1) > 19.9999 .and.&
            euls(2) < 20.0001 .and. euls(2) > 19.9999 .and.&
            euls(3) < 20.0001 .and. euls(3) > 19.9999 ) passed = .true.
        if( .not. passed ) call simple_stop('Euler composer corrupt!')
        write(*,'(a)') 'SIMPLE_ORI_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ori

end module simple_ori
