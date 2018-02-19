! an orientation
module simple_ori
#include "simple_lib.f08"
use simple_hash,   only: hash
use simple_chash,  only: chash
implicit none

public :: ori, test_ori, test_ori_dists
private

!>  orientation parameter stuct and operations
type :: ori
    private
    real        :: euls(3)=0.        !< Euler angle
    real        :: normal(3)=0.      !< Fourier plane normal
    real        :: rmat(3,3)=0.      !< rotation matrix
    type(hash)  :: htab              !< hash table for the parameters
    type(chash) :: chtab             !< hash table for the filenames etc.
    logical     :: existence=.false. !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new_ori
    procedure          :: new => new_ori
    procedure          :: new_ori_clean
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
    procedure, private :: rnd_romat
    procedure          :: rnd_euler_1
    procedure          :: rnd_euler_2
    generic            :: rnd_euler => rnd_euler_1, rnd_euler_2
    procedure          :: rnd_ori
    procedure          :: rnd_inpl
    procedure          :: rnd_shift
    procedure          :: revshsgn
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
    procedure, private :: getter_1
    procedure, private :: getter_2
    generic            :: getter => getter_1, getter_2
    procedure          :: get_2Dshift
    procedure          :: get_3Dshift
    procedure          :: get_state
    procedure          :: hash_size
    procedure          :: hash_keys
    procedure          :: hash_vals
    procedure          :: chash_size
    procedure          :: chash_nmax
    procedure          :: isthere
    procedure          :: isstatezero
    procedure          :: isevenodd
    procedure          :: iseven
    procedure          :: isodd
    procedure          :: key_is_real
    procedure          :: ori2str
    procedure          :: ori2strlen_trim
    ! PRINTING & I/O
    procedure          :: print_mat
    procedure          :: print_ori
    procedure          :: write
    procedure          :: read
    ! CALCULATORS
    procedure          :: round_shifts
    procedure          :: mul_shifts
    procedure, private :: shift
    procedure, private :: compeuler
    procedure          :: compose3d2d
    procedure          :: map3dshift22d
    procedure          :: mirror3d
    procedure          :: mirror2d
    procedure, private :: geodesic_dist
    procedure, private :: geodesic_dist_scaled
    procedure          :: oripair_diverse
    procedure          :: oripair_diverse_projdir
    procedure          :: ori_generator
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

real,       parameter :: zvec(3)  = [0.,0.,1.]
integer,    parameter :: NNAMES   = 10
class(ori), pointer   :: class_self1=>null(), class_self2=>null(), class_self3=>null()
real                  :: angthres = 0.

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
        call self%set_euler([0.,0.,0.])
        call self%htab%set('e1',0.)
        call self%htab%set('e2',0.)
        call self%htab%set('e3',0.)
        call self%htab%set('x',0.)
        call self%htab%set('y',0.)
        call self%htab%set('dist',180.)
        call self%htab%set('state',1.)
        call self%htab%set('state_balance',1.)
        call self%htab%set('frac',0.)
        call self%htab%set('eo',-1.) ! -1. is default (low-pass set); 0. for even; 1. for odd
        self%chtab = chash(NNAMES)
        self%existence = .true.
    end subroutine new_ori

    !>  \brief  is a constructor
    subroutine new_ori_clean( self )
        class(ori), intent(inout) :: self
        self%chtab = chash(NNAMES)
        self%existence = .true.
    end subroutine new_ori_clean

    !>  \brief  is a parameterized constructor
    subroutine ori_from_rotmat( self, rotmat )
        class(ori), intent(inout) :: self
        real, intent(in)          :: rotmat(3,3) !< rotation matrix
        call self%new_ori_clean
        self%rmat = rotmat
        self%euls = m2euler(self%rmat)
        call self%htab%set('e1',self%euls(1))
        call self%htab%set('e2',self%euls(2))
        call self%htab%set('e3',self%euls(3))
        self%normal = matmul(zvec, self%rmat)
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
        self_out%euls   = self_in%euls
        self_out%normal = self_in%normal
        self_out%rmat   = self_in%rmat
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
    subroutine set_euler( self, euls  )
        class(ori), intent(inout) :: self
        real, intent(in)          :: euls(3)!< Euler angle
        self%rmat   = euler2m(euls(1),euls(2),euls(3))
        self%euls   = m2euler(self%rmat)
        call self%htab%set('e1',self%euls(1))
        call self%htab%set('e2',self%euls(2))
        call self%htab%set('e3',self%euls(3))
        self%normal = matmul(zvec, self%rmat)
    end subroutine set_euler

    !>  \brief  is a setter
    subroutine e1set( self, e1 )
        class(ori), intent(inout) :: self
        real, intent(in)          :: e1     !< Euler angle
        call self%set_euler([e1,self%euls(2),self%euls(3)])
    end subroutine e1set

    !>  \brief  is a setter
    subroutine e2set( self, e2 )
        class(ori), intent(inout) :: self
        real, intent(in)          :: e2      !< Euler angle
        call self%set_euler([self%euls(1),e2,self%euls(3)])
    end subroutine e2set

    !>  \brief  is a setter
    subroutine e3set( self, e3 )
        class(ori), intent(inout) :: self
        real, intent(in)          :: e3      !< Euler angle
        call self%set_euler([self%euls(1),self%euls(2),e3])
    end subroutine e3set

    !>  \brief  is a setter
    subroutine swape1e3( self )
        class(ori), intent(inout) :: self
        real :: e
        e = self%euls(1)
        self%euls(1) = self%euls(3)
        self%euls(3) = e
        call self%set_euler(self%euls)
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

    !>  \brief  for generating a random rotation matrix
    subroutine rnd_romat( self )
        ! Fast Random Rotation Matrices, Arvo, Graphics Gems III, 1992
        use simple_rnd, only: ran3
        class(ori), intent(inout) :: self
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
        self%rmat(1,1) = vx * sx - ct
        self%rmat(1,2) = vx * sy - st
        self%rmat(1,3) = vx * vz
        self%rmat(2,1) = vy * sx + st
        self%rmat(2,2) = vy * sy - ct
        self%rmat(2,3) = vy * vz
        self%rmat(3,1) = vz * sx
        self%rmat(3,2) = vz * sy
        self%rmat(3,3) = 1. - z
    end subroutine rnd_romat

    !>  \brief  for generating a random Euler angle
    subroutine rnd_euler_1( self, eullims )
        class(ori),     intent(inout) :: self
        real, optional, intent(inout) :: eullims(3,2)!< Euler angles
        logical :: found
        if( present(eullims) )then
            found = .false.
            do while( .not.found )
                call self%rnd_romat
                self%euls = m2euler(self%rmat)
                if( self%euls(1) >= eullims(1,2) )cycle
                if( self%euls(2) >= eullims(2,2) )cycle
                if( self%euls(3) >= eullims(3,2) )cycle
                if( self%euls(1) < eullims(1,1) )cycle
                if( self%euls(2) < eullims(2,1) )cycle
                if( self%euls(3) < eullims(3,1) )cycle
                found = .true.
            end do
        else
            call self%rnd_romat
            self%euls = m2euler(self%rmat)
        endif
        call self%set('e1',self%euls(1))
        call self%set('e2',self%euls(2))
        call self%set('e3',self%euls(3))
        self%normal = matmul(zvec, self%rmat)
    end subroutine rnd_euler_1

    !>  \brief  for generating a random Euler angle neighbour to o_prev
    subroutine rnd_euler_2( self, o_prev, athres, proj )
        use simple_math, only: deg2rad
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
        use simple_rnd, only: ran3
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
        use simple_rnd, only: ran3
        class(ori), intent(inout)  :: self
        real, intent(in), optional :: trs         !< threshold
        if( present(trs) )call self%rnd_shift(trs)
        call self%e3set(ran3()*359.99)
    end subroutine rnd_inpl

    !>  \brief  for generating random in-plane parameters
    subroutine rnd_shift( self, trs )
        use simple_rnd, only: ran3
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

    !>  \brief  for reversing the shift signs (to fit convention)
    subroutine revshsgn( self )
        class(ori), intent(inout) :: self
        real                      :: x, y
        x = self%get('x')
        y = self%get('y')
        call self%set('x', -x)
        call self%set('y', -y)
    end subroutine revshsgn

    !>  \brief  reads all orientation info (by line) into the hash-tables
    subroutine str2ori( self, line )
        use simple_sauron, only: sauron_line_parser
        class(ori),       intent(inout) :: self
        character(len=*), intent(inout) :: line
        logical :: isthere(3)
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
        euls = self%euls
    end function get_euler

    !>  \brief  is a getter
    pure function e1get( self ) result( e1 )
        class(ori), intent(in) :: self
        real :: e1
        e1 = self%euls(1)
    end function e1get

    !>  \brief  is a getter
    pure function e2get( self ) result( e2 )
        class(ori), intent(in) :: self
        real :: e2
        e2 = self%euls(2)
    end function e2get

    !>  \brief  is a getter
    pure function e3get( self ) result( e3 )
        class(ori), intent(in) :: self
        real :: e3
        e3 = self%euls(3)
    end function e3get

    !>  \brief  is a getter
    pure function get_normal( self ) result( normal )
        class(ori), intent(in) :: self
        real :: normal(3)
        normal = self%normal
    end function get_normal

    !>  \brief  is a getter
    pure function get_mat( self ) result( mat )
        class(ori), intent(in) :: self
        real, dimension(3,3) :: mat
        mat = self%rmat
    end function get_mat

    !>  \brief  is a getter
    function get( self, key ) result( val )
        class(ori), intent(inout)    :: self
        character(len=*), intent(in) :: key
        real :: val
        val = self%htab%get(key)
    end function get

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
        sz = self%htab%size_of_hash()
    end function hash_size

    !>  \brief  returns the keys of the hash
    function hash_keys( self ) result( keys )
        class(ori), intent(inout) :: self
        character(len=32), allocatable :: keys(:)
        keys = self%htab%get_keys()
    end function hash_keys

    !>  \brief  returns the keys of the hash
    function hash_vals( self ) result( vals )
        class(ori), intent(inout) :: self
        real(kind=4), allocatable :: vals(:)
        vals = self%htab%get_vals()
    end function hash_vals

    !>  \brief  returns size of chash
    function chash_size( self ) result( sz )
        class(ori), intent(in) :: self
        integer :: sz
        sz = self%chtab%size_of_chash()
    end function chash_size

    !>  \brief  returns size of chash
    function chash_nmax( self ) result( nmax )
        class(ori), intent(in) :: self
        integer :: nmax
        nmax = self%chtab%get_nmax()
    end function chash_nmax

    !>  \brief  check for presence of key in the ori hash
    function isthere( self, key ) result( found )
        class(ori),       intent(inout) :: self
        character(len=*), intent(in)    :: key
        logical :: hash_found, chash_found, found
        hash_found  = self%htab%isthere(key)
        chash_found = self%chtab%isthere(key)
        found = .false.
        if( hash_found .or. chash_found ) found = .true.
    end function isthere

    !>  \brief  whether state is zero
    logical function isstatezero( self )
        class(ori),       intent(inout) :: self
        isstatezero = (self%get_state() == 0)
    end function isstatezero

    !>  \brief  whether orientation has been atributed an even/odd partition
    logical function isevenodd( self )
        class(ori),       intent(inout) :: self
        if( self%isthere('eo') )then
            isevenodd = self%htab%get('eo') > -.5
        else
            isevenodd = .false.
        endif
    end function isevenodd

    !>  \brief  whether orientation is part of the even partition
    logical function iseven( self )
        class(ori), intent(inout) :: self
        real :: val
        val = self%htab%get('eo')
        iseven = (val > -0.5) .and. (val < 0.5)
    end function iseven

    !>  \brief  whether orientation is part of the odd partition
    logical function isodd( self )
        class(ori),       intent(inout) :: self
        isodd = self%htab%get('eo') > 0.5
    end function isodd

    !>  \brief  is for checking whether key maps to a real value or not
    function key_is_real( self, key ) result( is )
        class(ori),       intent(inout) :: self
        character(len=*), intent(in)    :: key
        logical :: hash_found, chash_found, is
        hash_found  = self%htab%isthere(key)
        chash_found = self%chtab%isthere(key)
        is = .false.
        if( hash_found )then
            is = .true.
            if( chash_found ) stop 'ERROR, ambigous keys; simple_ori :: key_is_real'
        endif
    end function key_is_real

    !>  \brief  joins the hashes into a string that represent the ori
    function ori2str( self ) result( str )
        class(ori), intent(inout) :: self
        character(len=:), allocatable :: str, str_chtab, str_htab
        integer :: sz_chash, sz_hash
        sz_chash = self%chtab%size_of_chash()
        sz_hash  = self%htab%size_of_hash()
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

    !<  \brief  to print the rotation matrix
    subroutine print_mat( self )
        class(ori), intent(inout) :: self
        write(*,*) self%rmat(1,1), self%rmat(1,2), self%rmat(1,3), &
            &      self%rmat(2,1), self%rmat(2,2), self%rmat(2,3), &
            &      self%rmat(3,1), self%rmat(3,2), self%rmat(3,3)
    end subroutine print_mat

    !>  \brief prints oris data based on the existence of intent(in) character
    !! strings whose variable names are equivalent to the variables contained in
    !! the align derived type. The input character variables are printed first,
    !! preceeding printing of their associated value (key-value style for simple
    !! parsing of data)
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

    !>  \brief  rounds the origin shifts
    subroutine mul_shifts( self, mul )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: mul
        call self%set('x', mul * self%get('x'))
        call self%set('y', mul * self%get('y'))
    end subroutine mul_shifts

    !>  \brief  corrects the Euler angle bounds
    subroutine shift( self )
        class(ori), intent(inout) :: self
        logical :: doshift
        doshift = .false.
        if( self%euls(1) < 0. .or. self%euls(1) > 360. ) doshift = .true.
        if( self%euls(2) < 0. .or. self%euls(2) > 180. ) doshift = .true.
        if( self%euls(3) < 0. .or. self%euls(3) > 360. ) doshift = .true.
        if( doshift ) self%euls = m2euler(self%rmat)
    end subroutine shift

    !>  \brief  is for composing Euler angles
    !! multiplication of two rotation matrices commute (n>2 not)
    !! the composed rotation matrix is constructed
    !! \param self1,self2 ori class type rotational matrices
    function compeuler( self1, self2 ) result( self_out )
        class(ori), intent(in) :: self1, self2
        type(ori) :: self_out
        call self_out%new_ori
        self_out%rmat = matmul(self2%rmat,self1%rmat)  ! multiplication of two rotation matrices commute (n>2 not)
        ! the composed rotation matrix is constructed
        ! convert to euler
        self_out%euls = m2euler(self_out%rmat)
        call self_out%set('e1',self_out%euls(1))
        call self_out%set('e2',self_out%euls(2))
        call self_out%set('e3',self_out%euls(3))
        self_out%normal = matmul(zvec, self_out%rmat)
    end function compeuler

    !>  combines 3d and 2d oris and flags for ptcl mirroring
    !! \param ori3d,ori2d,o_out ori class type rotational matrices
    subroutine compose3d2d( ori3d, ori2d, o_out )
        use simple_math, only: make_transfmat, transfmat2inpls
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
    end subroutine compose3d2d

    subroutine map3dshift22d( self, sh3d )
        use simple_math, only: deg2rad
        class(ori), intent(inout) :: self
        real, intent(in)          :: sh3d(3) !< 3D shift
        real                      :: old_x,old_y,x,y,cx,cy
        real                      :: u(3),v(3),shift(3)
        real                      :: phi,cosphi,sinphi
        phi    = deg2rad( self%e3get() )
        cosphi = cos( phi )
        sinphi = sin( phi )
        old_x  = self%get('x')
        old_y  = self%get('y')
        ! 3D Shift rotated with respect to projdir
        shift = matmul( self%rmat,sh3d )
        ! Projection onto xy plane
        shift = shift-dot_product( shift,zvec )*zvec
        ! e3-composed xy plane unit vectors
        u = [  cosphi,sinphi,0. ]
        v = [ -sinphi,cosphi,0. ]
        ! composed shift clockwise components
        cx =  old_x*cosphi+old_y*sinphi+dot_product( u,shift )
        cy = -old_x*sinphi+old_y*cosphi+dot_product( v,shift )
        ! Identification
        x = cx*cosphi-cy*sinphi
        y = cx*sinphi+cy*cosphi
        ! the end
        call self%set( 'x',x )
        call self%set( 'y',y )
    end subroutine map3dshift22d

    !>  \brief  generates the opposite hand of an Euler angle
    !!          so that a set of Euler angles transformed by
    !!          this operation changes the handedness of the volume
    subroutine mirror3d( self )
        class(ori), intent(inout) :: self
        real :: mirr_mat(3,3)
        mirr_mat = 0.
        mirr_mat(1,1) = 1.
        mirr_mat(2,2) = 1.
        mirr_mat(3,3) = -1.
        self%rmat = matmul(mirr_mat,matmul(self%rmat,mirr_mat))
        ! the mirrored rotation matrix is constructed
        ! convert to euler
        self%euls = m2euler(self%rmat)
        call self%set('e1',self%euls(1))
        call self%set('e2',self%euls(2))
        call self%set('e3',self%euls(3))
        self%normal = matmul(zvec, self%rmat)
    end subroutine mirror3d

    !>  \brief  generates the mirror of the projection
    !!          can be equivalently accomplished by mirror('y')
    !!          the image after projection
    subroutine mirror2d( self )
        use simple_math, only: rad2deg
        class(ori), intent(inout) :: self
        type(ori) :: old
        old = self
        self%rmat = euler2m(self%euls(1), 180.+self%euls(2), 180.-self%euls(3))
        ! the mirrored rotation matrix is constructed
        ! convert to euler
        self%euls = m2euler(self%rmat)
        call self%set('e1',self%euls(1))
        call self%set('e2',self%euls(2))
        call self%set('e3',self%euls(3))
        self%normal = matmul(zvec, self%rmat)
    end subroutine mirror2d

    !! \param self1,self2 ori class type rotational matrices
    subroutine oripair_diverse( self1, self2 )
        use simple_opt_simplex, only: opt_simplex
        use simple_opt_spec,    only: opt_spec
        class(ori), target, intent(inout) :: self1, self2
        real, parameter                   :: TOL=1e-6
        type(opt_spec)                    :: ospec
        type(opt_simplex)                 :: opt
        real                              :: dist, lims(3,2)
        class(*), pointer                 :: fun_self
        if( .not. self1%exists() ) call self1%new_ori
        if( .not. self2%exists() ) call self2%new_ori
        call self1%rnd_euler ! fixed
        call self2%rnd_euler ! changing
        ! associate class pointers
        class_self1 => self1
        class_self2 => self2
        lims(1,1) = 0.
        lims(1,2) = 359.99
        lims(2,1) = 0.
        lims(2,2) = 180.
        lims(3,1) = 0.
        lims(3,2) = 359.99
        call ospec%specify('simplex', 3, ftol=TOL, limits=lims)
        ospec%x(1) = self2%e1get()
        ospec%x(2) = self2%e2get()
        ospec%x(3) = self2%e3get()
        call ospec%set_costfun(costfun_diverse_1)
        call opt%new(ospec)
        call opt%minimize(ospec, fun_self, dist)
        call self2%set_euler(ospec%x)
        call ospec%kill
        call opt%kill
    end subroutine oripair_diverse

    subroutine oripair_diverse_projdir( self1, self2, rangthres )
        use simple_opt_simplex, only: opt_simplex
        use simple_opt_spec,    only: opt_spec
        class(ori), target, intent(inout) :: self1, self2
        real,               intent(in)    :: rangthres
        real, parameter                   :: TOL=1e-6
        type(opt_spec)                    :: ospec
        type(opt_simplex)                 :: opt
        real                              :: e3,dist,lims(2,2)
        class(*), pointer                 :: fun_self => null()
        if( .not. self1%exists() )stop 'need a fixed ori! simple_oris::oripair_diverse_projdir'
        if( .not. self2%exists() )call self2%new_ori
        angthres = rangthres
        e3 = self1%e3get()
        call self2%rnd_euler_2( self1, angthres, .true. )
        ! associate class pointers
        class_self1 => self1
        class_self2 => self2
        lims(1,1) = 0.
        lims(1,2) = 359.99
        lims(2,1) = 0.
        lims(2,2) = 180.
        call ospec%specify('simplex', 2, ftol=TOL, limits=lims)
        ospec%x(1) = self2%e1get()
        ospec%x(2) = self2%e2get()
        call ospec%set_costfun(costfun_diverse_projdir)
        call opt%new(ospec)
        call opt%minimize(ospec, fun_self, dist)
        call self2%set_euler([ospec%x(1),ospec%x(2),e3])
        call ospec%kill
        call opt%kill
        angthres = 0.
    end subroutine oripair_diverse_projdir

    !>  \brief creates a spatial median or maximally diverse rotation matrix
    !!  if mode='median' this function creates a spatial median rotation matrix
    !!          that is "in between" the two inputted rotation matrices in a geodesic
    !!          distance sense and if mode='diverse' it creates the rotation matrix that
    !!          is maximally diverse with respect to the inputted two
    !! \param mode 'median' or 'diverse'
    !! \param self1,self2 ori class type rotaional matrices
    function ori_generator( self1, self2, mode ) result( oout )
        use simple_opt_simplex, only: opt_simplex
        use simple_opt_spec,    only: opt_spec
        class(ori), target, intent(in) :: self1, self2
        character(len=*),   intent(in) :: mode
        real,      parameter           :: TOL=1e-6
        type(ori)                      :: oout
        type(ori), target              :: otst
        type(opt_spec)                 :: ospec
        type(opt_simplex)              :: opt
        real                           :: dist, lims(3,2)
        class(*), pointer              :: fun_self => null()
        call oout%new_ori
        call otst%new_ori
        call otst%rnd_euler
        ! set class pointers
        class_self1 => self1
        class_self2 => self2
        class_self3 => otst
        ! init
        lims(1,1) = 0.
        lims(1,2) = 359.99
        lims(2,1) = 0.
        lims(2,2) = 180.
        lims(3,1) = 0.
        lims(3,2) = 359.99
        call ospec%specify('simplex', 3, ftol=TOL, limits=lims)
        ospec%x(1) = otst%e1get()
        ospec%x(2) = otst%e2get()
        ospec%x(3) = otst%e3get()
        call opt%new(ospec)
        select case(mode)
            case('median')
                call ospec%set_costfun(costfun_median)
            case('diverse')
                call ospec%set_costfun(costfun_diverse_2)
            case DEFAULT
                stop 'unsupported mode inputted; simple_ori :: ori_generator'
        end select
        call opt%minimize(ospec, fun_self, dist)
        call oout%set_euler(ospec%x)
        call ospec%kill
        call opt%kill
    end function ori_generator

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
    !!          129, 883–886 (2007)
    !! \param self1,self2 ori class type rotational matrices
    pure real function geodesic_dist( self1, self2 )
        class(ori), intent(in) :: self1, self2
        real :: Imat(3,3), sumsq, diffmat(3,3)
        Imat      = 0.
        Imat(1,1) = 1.
        Imat(2,2) = 1.
        Imat(3,3) = 1.
        diffmat = Imat-matmul(self1%rmat,transpose(self2%rmat))
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
        use simple_math, only: myacos
        class(ori), intent(in) :: self1, self2
        real :: dist
        dist = myacos(dot_product(self1%normal,self2%normal))
    end function euldist

    !>  \brief  calculates the distance (in radians) btw all df:s
    !! \param self1,self2 ori class type rotaional matrices
    pure function inpldist( self1, self2 ) result( dist )
        use simple_math, only: myacos, rotmat2d
        class(ori), intent(in) :: self1, self2
        real :: dist
        dist = ((self1.inplrotdist.self2)+(self1.euldist.self2))/2.
    end function inpldist

    !>  \brief  calculates the distance (in radians) btw the in-plane rotations
    !! \param self1,self2 ori class type rotaional matrices
    pure function inplrotdist( self1, self2 ) result( dist )
        use simple_math, only: myacos, rotmat2d
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
            call self%chtab%kill
            self%existence = .false.
        endif
    end subroutine kill

    ! PRIVATE STUFF

    !>  \brief  makes a rotation matrix from a Spider format Euler triplet
    !! \param e1,e2,e3 Euler triplet
    pure function euler2m( e1, e2, e3 ) result( r )
        real, intent(in)     :: e1, e2, e3
        real, dimension(3,3) :: r1, r2, r3, r, tmp
        r1 = rotmat(e1,3) ! rotation around z
        r2 = rotmat(e2,2) ! tilt
        r3 = rotmat(e3,3) ! rotation around z
        ! order of multiplication is r3r2r1
        tmp = matmul(r3,r2)
        r = matmul(tmp,r1)
    end function euler2m

    !>  \brief  calculates the normal of a Fourier plane
    pure function euler2vec( phi, theta, psi ) result( normal )
        real, intent(in) :: phi, theta, psi
        real :: normal(3), mat(3,3)
        mat = euler2m( phi, theta, psi )
        normal = matmul(zvec,mat)
    end function euler2vec

    !>  \brief  returns a Spider format Euler triplet, given a rotation matrix
    pure function m2euler( mat ) result ( euls )
        use simple_math, only: rad2deg, myacos
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

    function costfun_diverse_1( fun_self, vec, D ) result( dist )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: dist
        call class_self2%set_euler(vec)
        ! negation because we try to maximize the distance
        dist = -class_self1%geodesic_dist(class_self2)
    end function costfun_diverse_1

    function costfun_diverse_2( fun_self, vec, D ) result( dist )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: d1, d2, dist
        call class_self3%set_euler(vec)
        d1   = class_self3%geodesic_dist(class_self1)
        d2   = class_self3%geodesic_dist(class_self2)
        dist = -d1-d2
    end function costfun_diverse_2

    !>  \brief  uses an angular threshold
    function costfun_diverse_projdir( fun_self, vec, D ) result( dist )
        use simple_math, only: rad2deg
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: dist
        call class_self2%e1set(vec(1))
        call class_self2%e2set(vec(2))
        dist = class_self1.euldist.class_self2
        if(rad2deg(dist) > angthres) dist=0.
        ! negation because we try to maximize the distance
        dist = -dist
    end function costfun_diverse_projdir

    function costfun_median( fun_self, vec, D ) result( dist )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: d1, d2, dist
        call class_self3%set_euler(vec)
        d1   = class_self3%geodesic_dist(class_self1)
        d2   = class_self3%geodesic_dist(class_self2)
        dist = abs(d1-d2)
    end function costfun_median

    ! UNIT TESTS

    subroutine test_ori_dists
        use simple_math, only: rad2deg
        type(ori) :: o1, o2, omed
        real    :: frob_lim, dist, adist
        integer :: itst
        integer, parameter :: NTESTS=5000
        frob_lim = 2.*sqrt(2.)
        call o1%new_ori
        call o2%new_ori
        do itst=1,NTESTS
            call o1%rnd_euler
            call o2%rnd_euler
            dist = o1%geodesic_dist(o2)
            if( dist > frob_lim )then
                print *, 'ERROR dist > lim, dist/lim: ', dist, frob_lim
            endif
        end do
        adist = 0.
        do itst=1,NTESTS
            call o1%rnd_euler
            call o2%rnd_euler
            omed  = o1%ori_generator(o2, 'median')
            dist  = abs(omed%geodesic_dist(o1)-omed%geodesic_dist(o2))
            adist = adist+dist
        end do
        print *, 'average distance for spatial medians: ', adist/real(NTESTS)
        adist = 0.
        do itst=1,NTESTS
            call o1%rnd_euler
            call o2%rnd_euler
            omed  = o1%ori_generator(o2, 'diverse')
            dist  = (omed%geodesic_dist(o1)+omed%geodesic_dist(o2))/2.
            adist = adist+dist
        end do
        adist = adist/real(NTESTS)
        print *, 'average distance for spatial extremes: ', adist, ' maxdist: ',&
        frob_lim, ' frac_of_max: ', 100*(adist/frob_lim), ' %'
        write(*,'(a)') 'SIMPLE_ORI_DISTS_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ori_dists

    !>  \brief  is the unit test for the ori class
    subroutine test_ori
    ! this test only tests the Euler part of ori,
    ! the rest is tested in the oris class
        use simple_stat, only: pearsn
        use simple_math, only: rad2deg
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
        if( .not. passed ) stop 'Euler assignment/getters corrupt!'
        passed = .false.
        call e1%e1set(99.)
        call e1%e2set(98.)
        call e1%e3set(97.)
        euls = e1%get_euler()
        if( abs(euls(1)-99.+euls(2)-98.+euls(3)-97.) < 0.001 ) passed = .true.
        if( .not. passed ) stop 'Euler e-setters corrupt!'
        passed = .false.
        call e1%set_euler([0.,0.,0.])
        normal = e1%get_normal()
        if( abs(normal(1))+abs(normal(2))+abs(normal(3)-1.)<3.*TINY ) passed = .true.
        if( .not. passed ) stop 'Euler normal derivation corrupt!'
        passed = .false.
        mat = e1%get_mat()
        if( abs(mat(1,1)-1.)+abs(mat(2,2)-1.)+abs(mat(3,3)-1.)<3.*TINY )then
            mat(1,1) = 0.
            mat(2,2) = 0.
            mat(3,3) = 0.
            if( abs(sum(mat)) < TINY ) passed = .true.
        endif
        if( .not. passed ) stop 'Euler rotation matrix derivation corrupt!'
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
        if( .not. passed ) stop 'Euler composer corrupt!'
        write(*,'(a)') 'SIMPLE_ORI_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ori

end module simple_ori
