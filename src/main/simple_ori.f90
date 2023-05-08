! an orientation
module simple_ori
use simple_defs
use simple_defs_ori
use simple_hash
use simple_chash
use simple_rnd
use simple_is_check_assert
use simple_strings
use simple_math
use simple_stat
use simple_nrtxtfile
implicit none

public :: ori, test_ori
public :: euler2m, m2euler, euler_dist, euler_inplrotdist, euler_compose, euler_mirror
public :: geodesic_frobdev
private
#include "simple_local_flags.inc"

type :: ori
    private
    real        :: pparms(N_PTCL_ORIPARAMS) = 0. !< hardcoded per-particle parameters (see simple_defs_ori)
    type(hash)  :: htab                          !< hash table for dynamic parameters
    type(chash) :: chtab                         !< hash table for filenames etc.
    logical     :: is_ptcl   = .false.           !< to indicate whether the info is of particle kind
    logical     :: existence = .false.           !< to indicate existence
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
    procedure          :: append_ori
    procedure          :: delete_entry
    procedure          :: delete_2Dclustering
    procedure          :: delete_3Dalignment
    procedure          :: transfer_2Dparams
    procedure          :: transfer_3Dparams
    procedure          :: set_euler
    procedure          :: e1set
    procedure          :: e2set
    procedure          :: e3set
    procedure          :: swape1e3
    procedure          :: set_shift
    procedure          :: set_dfx, set_dfy
    procedure, private :: set_1
    procedure, private :: set_2
    generic            :: set => set_1, set_2
    procedure, private :: rnd_euler_1
    procedure, private :: rnd_euler_2
    procedure, private :: rnd_euler_3
    procedure, private :: rnd_euler_4
    generic            :: rnd_euler => rnd_euler_1, rnd_euler_2, rnd_euler_3, rnd_euler_4
    procedure          :: rnd_ori
    procedure          :: rnd_inpl
    procedure          :: rnd_shift
    procedure          :: str2ori
    procedure          :: set_boxfile
    ! GETTERS
    procedure          :: exists
    procedure          :: is_particle
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
    procedure          :: set_state
    procedure          :: get_class
    procedure          :: get_dfx, get_dfy
    procedure          :: isthere
    procedure          :: ischar
    procedure          :: isstatezero
    procedure          :: has_been_searched
    procedure          :: pparms2str
    procedure          :: pparms_strlen
    procedure          :: ori2str
    procedure          :: ori2prec
    procedure          :: prec2ori
    procedure          :: ori_strlen_trim
    procedure          :: ori2chash
    procedure          :: chash2ori
    procedure          :: get_ctfvars
    procedure          :: get_axis_angle
    procedure          :: set_ctfvars
    procedure          :: get_keys
    ! PRINTING & I/O
    procedure          :: print_mat
    procedure          :: print_ori
    procedure          :: write
    procedure          :: write2bild
    procedure          :: read
    ! CALCULATORS
    procedure          :: round_shifts
    procedure, private :: compeuler
    procedure          :: compose3d2d
    procedure          :: map3dshift22d
    procedure          :: mirror3d
    procedure          :: mirror2d
    procedure          :: transp
    procedure          :: geodesic_dist_frobdev
     generic           :: operator(.geod.)        => geodesic_dist_frobdev
    procedure          :: geodesic_dist_trace
    procedure, private :: euldist
    procedure, private :: inplrotdist
    procedure          :: compose                 => compeuler
    generic            :: operator(.euldist.)     => euldist
    generic            :: operator(.inplrotdist.) => inplrotdist
    ! DESTRUCTORS
    procedure          :: reset_pparms
    procedure          :: kill_hash
    procedure          :: kill_chash
    procedure          :: kill
end type ori

interface ori
    module procedure constructor
end interface

real, parameter :: zvec(3) = [0.,0.,1.]

contains

    ! CONSTRUCTORS

    function constructor( is_ptcl ) result( self )
        logical, intent(in) :: is_ptcl
        type(ori) :: self
        call self%new_ori(is_ptcl)
    end function constructor

    subroutine new_ori( self, is_ptcl )
        class(ori), intent(inout) :: self
        logical,    intent(in)    :: is_ptcl
        call self%kill
        self%is_ptcl = is_ptcl
        if( self%is_ptcl )then
            self%htab  = hash(0)
            self%chtab = chash(0)
        else
            self%htab  = hash()
            self%chtab = chash()
        endif
        self%existence = .true.
    end subroutine new_ori

    !>  \brief  is a parameterized constructor
    subroutine ori_from_rotmat( self, rotmat, is_ptcl )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: rotmat(3,3) !< rotation matrix
        logical,    intent(in)    :: is_ptcl
        real :: euls(3)
        call self%new(is_ptcl)
        euls = m2euler(rotmat)
        call self%set_euler(euls)
    end subroutine ori_from_rotmat

    ! SETTERS

    !>  \brief  sets parameters for particle rejection
    subroutine reject( self )
        class(ori), intent(inout) :: self
        call self%set_euler([0., 0., 0.])
        call self%set_shift([0., 0.])
        if( self%is_ptcl )then
            self%pparms(I_STATE)     =  0.
            self%pparms(I_CORR)      = -1.
            self%pparms(I_SPECSCORE) =  0.
            self%pparms(I_EO)        = -1.
        else
            call self%htab%set('state', 0.)
            if( self%isthere('corr') )     call self%htab%set('corr',     -1.)
            if( self%isthere('specscore') )call self%htab%set('specscore', 0.)
            if( self%isthere('eo') )       call self%htab%set('eo', -1.)
        endif
    end subroutine reject

    subroutine assign_ori( self_out, self_in )
        type(ori), intent(in)     :: self_in
        class(ori), intent(inout) :: self_out
        call self_out%copy(self_in)
    end subroutine assign_ori

    subroutine copy_ori( self_out, self_in )
        class(ori), intent(in)     :: self_in
        class(ori), intent(inout)  :: self_out
        character(len=KEYLEN)      :: key
        type(str4arr), allocatable :: keys(:)
        integer :: sz, i, ind
        self_out%pparms = self_in%pparms
        self_out%htab   = self_in%htab
        self_out%chtab  = self_in%chtab
        if( self_out%is_ptcl.neqv.self_in%is_ptcl )then
            if( self_in%is_ptcl )then
                do i=1,N_PTCL_ORIPARAMS
                    if( oriparam_isthere(i, self_in%pparms(i)) )then
                        key = get_oriparam_flag(i)
                        call self_out%set(trim(key), self_in%pparms(i))
                    endif
                end do
            else ! self_out is ptcl, self_in is not
                sz   = self_in%htab%size_of()
                keys = self_in%htab%get_keys()
                do i=1,sz
                    ind = get_oriparam_ind(trim(keys(i)%str))
                    if( ind /= 0 )then
                        self_out%pparms(ind) = self_in%get(trim(keys(i)%str))
                    endif
                end do
            endif
        endif
    end subroutine copy_ori

    subroutine append_ori( self_out, self_in )
        class(ori), intent(in)     :: self_in
        class(ori), intent(inout)  :: self_out
        type(str4arr), allocatable :: keys(:)
        integer :: sz, i
        if(.not. self_out%is_ptcl .and. .not. self_in%is_ptcl) then
            sz   = self_in%htab%size_of()
            keys = self_in%htab%get_keys()
            do i=1,sz
                call self_out%set(trim(keys(i)%str), self_in%get(trim(keys(i)%str)))
            end do
            sz   = self_in%chtab%size_of()
            do i=1,sz
                call self_out%set(trim(self_in%chtab%get_key(i)), trim(self_in%chtab%get(i)))
            end do
        end if
    end subroutine append_ori

    subroutine delete_entry( self, key )
        class(ori),       intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer :: ind
        if( self%htab%isthere(key) )then
            call self%htab%delete(key)
        else
            call self%chtab%delete(key)
        endif
        ind = get_oriparam_ind(key)
        if( ind /= 0 ) self%pparms(ind) = 0. ! default value on init
    end subroutine delete_entry

    subroutine delete_2Dclustering( self, keepshifts, keepcls )
        class(ori),        intent(inout) :: self
        logical, optional, intent(in)    :: keepshifts, keepcls
        logical :: kkeepshifts, kkeepcls
        kkeepshifts = .true.
        if( present(keepshifts) ) kkeepshifts = keepshifts
        kkeepcls = .false.
        if( present(keepcls) ) kkeepcls = keepcls
        if( self%is_ptcl )then
            if( .not. kkeepcls ) self%pparms(I_CLASS) = 0.
            self%pparms(I_E3)    = 0.
            if( .not. kkeepshifts )then
                self%pparms(I_X) = 0.
                self%pparms(I_Y) = 0.
            endif
            self%pparms(I_CORR) = 0.
            self%pparms(I_FRAC) = 0.
        else
            if( .not. kkeepcls ) call self%htab%delete('class')
            call self%htab%delete('e3')
            if( .not. kkeepshifts )then
                call self%htab%delete('x')
                call self%htab%delete('y')
            endif
            call self%htab%delete('corr')
            call self%htab%delete('frac')
        endif
    end subroutine delete_2Dclustering

    subroutine delete_3Dalignment( self, keepshifts )
        class(ori),        intent(inout) :: self
        logical, optional, intent(in)    :: keepshifts
        logical :: kkeepshifts
        kkeepshifts = .false.
        if( present(keepshifts) ) kkeepshifts = keepshifts
        if( self%is_ptcl )then
            self%pparms(I_PROJ) = 0.
            self%pparms(I_E1)   = 0.
            self%pparms(I_E2)   = 0.
            self%pparms(I_E3)   = 0.
            if( .not. kkeepshifts )then
                self%pparms(I_X) = 0.
                self%pparms(I_Y) = 0.
            endif
            self%pparms(I_CORR) =  0.
            self%pparms(I_FRAC) =  0.
        else
            call self%htab%delete('proj')
            call self%htab%delete('e1')
            call self%htab%delete('e2')
            call self%htab%delete('e3')
            if( .not. kkeepshifts )then
                call self%htab%delete('x')
                call self%htab%delete('y')
            endif
            call self%htab%delete('corr')
            call self%htab%delete('frac')
        endif
    end subroutine delete_3Dalignment

    subroutine transfer_2Dparams( self_out, self_in )
        class(ori), intent(inout) :: self_out
        class(ori), intent(in)    :: self_in
        if( self_out%is_ptcl.neqv.self_in%is_ptcl ) THROW_HARD('non-conforming types (is_ptcl)')
        if( self_out%is_ptcl )then
            self_out%pparms(I_CLASS)     = self_in%pparms(I_CLASS)
            self_out%pparms(I_CORR)      = self_in%pparms(I_CORR)
            self_out%pparms(I_FRAC)      = self_in%pparms(I_FRAC)
            self_out%pparms(I_SPECSCORE) = self_in%pparms(I_SPECSCORE)
            self_out%pparms(I_UPDATECNT) = self_in%pparms(I_UPDATECNT)
            self_out%pparms(I_W)         = self_in%pparms(I_W)
            self_out%pparms(I_EO)        = self_in%pparms(I_EO)
        else
            call self_out%htab%set('class',    self_in%htab%get('class'))
            call self_out%htab%set('corr',     self_in%htab%get('corr'))
            call self_out%htab%set('frac',     self_in%htab%get('frac'))
            call self_out%htab%set('specscore',self_in%htab%get('specscore'))
            call self_out%htab%set('updatecnt',self_in%htab%get('updatecnt'))
            call self_out%htab%set('w',        self_in%htab%get('w'))
            call self_out%htab%set('eo',       self_in%htab%get('eo'))
        endif
        call self_out%set_euler(self_in%get_euler())
        call self_out%set_shift(self_in%get_2Dshift())
    end subroutine transfer_2Dparams

    subroutine transfer_3Dparams( self_out, self_in )
        class(ori), intent(inout) :: self_out
        class(ori), intent(in)    :: self_in
        if( self_out%is_ptcl.neqv.self_in%is_ptcl ) THROW_HARD('non-conforming types (is_ptcl)')
        if( self_out%is_ptcl )then
            self_out%pparms(I_PROJ)      = self_in%pparms(I_PROJ)
            self_out%pparms(I_CORR)      = self_in%pparms(I_CORR)
            self_out%pparms(I_FRAC)      = self_in%pparms(I_FRAC)
            self_out%pparms(I_SPECSCORE) = self_in%pparms(I_SPECSCORE)
            self_out%pparms(I_UPDATECNT) = self_in%pparms(I_UPDATECNT)
            self_out%pparms(I_W)         = self_in%pparms(I_W)
            self_out%pparms(I_EO)        = self_in%pparms(I_EO)
        else
            call self_out%htab%set('proj',     self_in%htab%get('proj'))
            call self_out%htab%set('corr',     self_in%htab%get('corr'))
            call self_out%htab%set('frac',     self_in%htab%get('frac'))
            call self_out%htab%set('specscore',self_in%htab%get('specscore'))
            call self_out%htab%set('updatecnt',self_in%htab%get('updatecnt'))
            call self_out%htab%set('w',        self_in%htab%get('w'))
            call self_out%htab%set('eo',       self_in%htab%get('eo'))
        endif
        call self_out%set_euler(self_in%get_euler())
        call self_out%set_shift(self_in%get_2Dshift())
    end subroutine transfer_3Dparams

    subroutine set_euler( self, euls )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: euls(3)
        real :: euls_here(3), rmat(3,3)
        rmat      = euler2m(euls)
        euls_here = m2euler(rmat)
        if( self%is_ptcl )then
            self%pparms(I_E1:I_E3) = euls_here
        else
            call self%htab%set('e1',euls_here(1))
            call self%htab%set('e2',euls_here(2))
            call self%htab%set('e3',euls_here(3))
        endif
    end subroutine set_euler

    subroutine e1set( self, e1 )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: e1
        if( self%is_ptcl )then
            self%pparms(I_E1) = e1
        else
            call self%htab%set('e1', e1)
        endif
    end subroutine e1set

    subroutine e2set( self, e2 )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: e2
        if( self%is_ptcl )then
            self%pparms(I_E2) = e2
        else
            call self%htab%set('e2', e2)
        endif
    end subroutine e2set

    subroutine e3set( self, e3 )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: e3
        if( self%is_ptcl )then
            self%pparms(I_E3) = e3
        else
            call self%htab%set('e3', e3)
        endif
    end subroutine e3set

    subroutine swape1e3( self )
        class(ori), intent(inout) :: self
        real :: e1, e3
        if( self%is_ptcl )then
            e1 = self%pparms(I_E1)
            e3 = self%pparms(I_E3)
            self%pparms(I_E1) = e3
            self%pparms(I_E3) = e1
        else
            e1 = self%htab%get('e1')
            e3 = self%htab%get('e3')
            call self%htab%set('e1', e3)
            call self%htab%set('e3', e1)
        endif
    end subroutine swape1e3

    subroutine set_shift( self, shvec )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: shvec(2)
        if( self%is_ptcl )then
            self%pparms(I_X) = shvec(1)
            self%pparms(I_Y) = shvec(2)
        else
            call self%htab%set('x', shvec(1))
            call self%htab%set('y', shvec(2))
        endif
    end subroutine set_shift

    subroutine set_dfx( self, dfx)
        class(ori), intent(inout) :: self
        real,       intent(in)    :: dfx
        if( self%is_ptcl )then
            self%pparms(I_DFX) = dfx
        else
            call self%htab%set('dfx', dfx)
        endif
    end subroutine set_dfx

    subroutine set_dfy( self, dfy)
        class(ori), intent(inout) :: self
        real,       intent(in)    :: dfy
        if( self%is_ptcl )then
            self%pparms(I_DFY) = dfy
        else
            call self%htab%set('dfy', dfy)
        endif
    end subroutine set_dfy

    subroutine set_1( self, key, val )
        class(ori),       intent(inout) :: self
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
        integer :: ind
        if( self%is_ptcl )then
            ind = get_oriparam_ind(key)
            if( ind == 0 )then
                call self%htab%set(key, val)
            else
                self%pparms(ind) = val
            endif
        else
            call self%htab%set(key, val)
        endif
    end subroutine set_1

    subroutine set_2( self, key, val )
        class(ori),       intent(inout) :: self
        character(len=*), intent(in)    :: key, val
        call self%chtab%set(key, val)
    end subroutine set_2

    !>  \brief  for generating a random Euler angle
    subroutine rnd_euler_1( self )
        class(ori), intent(inout) :: self
        real    :: euls(3), rmat(3,3)
        call rnd_romat(rmat)
        euls = m2euler(rmat)
        call self%set_euler(euls)
    end subroutine rnd_euler_1

    !>  \brief  for generating a random Euler angle
    subroutine rnd_euler_2( self, eullims )
        class(ori), intent(inout) :: self
        real,       intent(inout) :: eullims(3,2) !< Euler angle limits
        logical :: found
        real    :: euls(3), rmat(3,3)
        found = .false.
        do while( .not.found )
            call rnd_romat(rmat)
            euls = m2euler(rmat)
            if( euls(1) >= eullims(1,2) ) cycle
            if( euls(2) >= eullims(2,2) ) cycle
            if( euls(3) >= eullims(3,2) ) cycle
            if( euls(1) <  eullims(1,1) ) cycle
            if( euls(2) <  eullims(2,1) ) cycle
            if( euls(3) <  eullims(3,1) ) cycle
            found = .true.
        end do
        call self%set_euler(euls)
    end subroutine rnd_euler_2

    !>  \brief  for generating a random Euler angle neighbour to o_prev
    subroutine rnd_euler_3( self, o_prev, athres, eullims )
        class(ori), intent(inout) :: self         !< instance
        class(ori), intent(in)    :: o_prev       !< template ori
        real,       intent(in)    :: athres       !< Euler angle threshold in degrees
        real,       intent(inout) :: eullims(3,2) !< Euler angle limits
        real :: euls(3)
        euls(1) = o_prev%e1get()
        euls(2) = o_prev%e2get()
        euls(3) = o_prev%e3get()
        euls(1) = euls(1) + 2. * (ran3() - 0.5) * athres
        euls(2) = euls(2) + 2. * (ran3() - 0.5) * athres
        euls(3) = euls(3) + 2. * (ran3() - 0.5) * athres
        euls(1) = max(eullims(1,1), euls(1))
        euls(2) = max(eullims(2,1), euls(2))
        euls(3) = max(eullims(3,1), euls(3))
        euls(1) = min(eullims(1,2), euls(1))
        euls(2) = min(eullims(2,2), euls(2))
        euls(3) = min(eullims(3,2), euls(3))
        call self%set_euler(euls)
    end subroutine rnd_euler_3

    subroutine rnd_euler_4( self, o_prev, athres_proj, athres_inpl, eullims )
        class(ori), intent(inout) :: self                     !< instance
        class(ori), intent(in)    :: o_prev                   !< template ori
        real,       intent(in)    :: athres_proj, athres_inpl !< Euler angle thresholds in degrees
        real,       intent(inout) :: eullims(3,2)             !< Euler angle limits
        real :: euls(3)
        euls(1) = o_prev%e1get()
        euls(2) = o_prev%e2get()
        euls(3) = o_prev%e3get()
        euls(1) = euls(1) + 2. * (ran3() - 0.5) * athres_proj
        euls(2) = euls(2) + 2. * (ran3() - 0.5) * athres_proj
        euls(3) = euls(3) + 2. * (ran3() - 0.5) * athres_inpl
        euls(1) = max(eullims(1,1), euls(1))
        euls(2) = max(eullims(2,1), euls(2))
        euls(3) = max(eullims(3,1), euls(3))
        euls(1) = min(eullims(1,2), euls(1))
        euls(2) = min(eullims(2,2), euls(2))
        euls(3) = min(eullims(3,2), euls(3))
        call self%set_euler(euls)
    end subroutine rnd_euler_4

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
        if( present(eullims) )then
            call self%rnd_euler( eullims )
        else
            call self%rnd_euler
        endif
    end subroutine rnd_ori

    subroutine rnd_inpl( self, trs )
        class(ori), intent(inout)  :: self
        real, intent(in), optional :: trs !< threshold
        if( present(trs) ) call self%rnd_shift(trs)
        call self%e3set(ran3()*359.99)
    end subroutine rnd_inpl

    subroutine rnd_shift( self, trs )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: trs !< threshold
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

    subroutine str2ori( self, line, is_ptcl )
        use simple_sauron, only: sauron_ori_parser
        class(ori),       intent(inout) :: self
        character(len=*), intent(inout) :: line
        logical,          intent(in)    :: is_ptcl
        call self%new(is_ptcl)
        if( is_ptcl )then
            call sauron_ori_parser(line, self%pparms, self%htab, self%chtab)
        else
            call sauron_ori_parser(line, self%htab, self%chtab)
        endif
    end subroutine str2ori

    subroutine set_boxfile( self, boxfname, nptcls )
        class(ori),        intent(inout) :: self
        character(len=*),  intent(in)    :: boxfname
        integer, optional, intent(in)    :: nptcls
        type(nrtxtfile) :: boxfile
        integer         :: nptcls_here
        if( present(nptcls) )then
            nptcls_here = nptcls
            if( nptcls_here == 0 )then
                call self%set('nptcls',  0.)
                return
            endif
        else
            call boxfile%new(boxfname, 1)
            nptcls_here = boxfile%get_ndatalines()
            call boxfile%kill
        endif
        call self%set('boxfile', trim(boxfname))
        call self%set('nptcls',  real(nptcls_here))
    end subroutine set_boxfile

    ! GETTERS

    pure function exists( self ) result( t )
        class(ori), intent(in) :: self
        logical :: t
        t = self%existence
    end function exists

    pure function is_particle( self ) result( t )
        class(ori), intent(in) :: self
        logical :: t
        t = self%is_ptcl
    end function is_particle

    pure function get_euler( self ) result( euls )
        class(ori), intent(in) :: self
        real :: euls(3)
        if( self%is_ptcl )then
            euls(1) = self%pparms(I_E1)
            euls(2) = self%pparms(I_E2)
            euls(3) = self%pparms(I_E3)
        else
            euls(1) = self%htab%get('e1')
            euls(2) = self%htab%get('e2')
            euls(3) = self%htab%get('e3')
        endif
    end function get_euler

    pure function e1get( self ) result( e1 )
        class(ori), intent(in) :: self
        real :: e1
        if( self%is_ptcl )then
            e1 = self%pparms(I_E1)
        else
            e1 = self%htab%get('e1')
        endif
    end function e1get

    pure function e2get( self ) result( e2 )
        class(ori), intent(in) :: self
        real :: e2
        if( self%is_ptcl )then
            e2 = self%pparms(I_E2)
        else
            e2 = self%htab%get('e2')
        endif
    end function e2get

    pure function e3get( self ) result( e3 )
        class(ori), intent(in) :: self
        real :: e3
        if( self%is_ptcl )then
            e3 = self%pparms(I_E3)
        else
            e3 = self%htab%get('e3')
        endif
    end function e3get

    pure function get_normal( self ) result( normal )
        class(ori), intent(in) :: self
        real :: normal(3), rmat(3,3), euls(3)
        euls   = self%get_euler()
        rmat   = euler2m(euls)
        normal = matmul(zvec, rmat)
    end function get_normal

    pure function get_mat( self ) result( mat )
        class(ori), intent(in) :: self
        real, dimension(3,3) :: mat
        real :: euls(3)
        euls = self%get_euler()
        mat  = euler2m(euls)
    end function get_mat

    pure function get( self, key ) result( val )
        class(ori),       intent(in) :: self
        character(len=*), intent(in) :: key
        integer :: ind
        real    :: val
        if( self%is_ptcl )then
            ind = get_oriparam_ind(key)
            if( ind == 0 )then
                val = self%htab%get(key)
            else
                val = self%pparms(ind)
            endif
        else
            val = self%htab%get(key)
        endif
    end function get

    !>  \brief  is a getter with fixed length return string
    function get_static( self, key )result( val )
        class(ori),       intent(in) :: self
        character(len=*), intent(in) :: key
        character(len=STDLEN)        :: val
        val = trim(self%chtab%get_static(key))
    end function get_static

    subroutine getter_1( self, key, val )
        class(ori),                    intent(inout) :: self
        character(len=*),              intent(in)    :: key
        character(len=:), allocatable, intent(inout) :: val
        if( allocated(val) ) deallocate(val)
        val = self%chtab%get(key)
    end subroutine getter_1

    subroutine getter_2( self, key, val )
        class(ori),       intent(in)    :: self
        character(len=*), intent(in)    :: key
        real,             intent(inout) :: val
        integer :: ind
        if( self%is_ptcl )then
            ind = get_oriparam_ind(key)
            if( ind == 0 )then
                val = self%htab%get(key)
            else
                val = self%pparms(ind)
            endif
        else
            val = self%htab%get(key)
        endif
    end subroutine getter_2

    pure function get_2Dshift( self ) result( vec )
        class(ori), intent(in) :: self
        real :: vec(2)
        if( self%is_ptcl )then
            vec(1) = self%pparms(I_X)
            vec(2) = self%pparms(I_Y)
        else
            vec(1) = self%htab%get('x')
            vec(2) = self%htab%get('y')
        endif
    end function get_2Dshift

    pure function get_3Dshift( self ) result( vec )
        class(ori), intent(in) :: self
        real :: vec(3)
        if( self%is_ptcl )then
            vec(1) = self%pparms(I_X)
            vec(2) = self%pparms(I_Y)
            vec(3) = self%htab%get('z')
        else
            vec(1) = self%htab%get('x')
            vec(2) = self%htab%get('y')
            vec(3) = self%htab%get('z')
        endif
    end function get_3Dshift

    pure integer function get_state( self )
        class(ori), intent(in) :: self
        if( self%is_ptcl )then
            get_state = nint(self%pparms(I_STATE))
        else
            get_state = nint(self%htab%get('state'))
        endif
    end function get_state

    subroutine set_state( self, state )
        class(ori), intent(inout) :: self
        integer,    intent(in)    :: state
        if( self%is_ptcl )then
            self%pparms(I_STATE) = real(state)
        else
            call self%set('state', real(state))
        endif
    end subroutine set_state

    pure integer function get_class( self )
        class(ori), intent(in) :: self
        if( self%is_ptcl )then
            get_class = nint(self%pparms(I_CLASS))
        else
            get_class = nint(self%htab%get('class'))
        endif
    end function get_class

    pure real function get_dfx( self )
        class(ori), intent(in) :: self
        if( self%is_ptcl )then
            get_dfx = self%pparms(I_DFX)
        else
            get_dfx = self%htab%get('dfx')
        endif
    end function get_dfx

    pure real function get_dfy( self )
        class(ori), intent(in) :: self
        if( self%is_ptcl )then
            get_dfy = self%pparms(I_DFY)
        else
            get_dfy = self%htab%get('dfy')
        endif
    end function get_dfy

    pure function isthere( self, key ) result( found )
        class(ori),        intent(in) :: self
        character(len=*),  intent(in) :: key
        logical :: hash_found, chash_found, found
        integer :: ind
        found = .false.
        if( self%is_ptcl )then
            ind = get_oriparam_ind(key)
            if( ind /= 0 )then
                if( oriparam_isthere(ind, self%pparms(ind)) )then
                    found = .true.
                    return
                endif
            endif
            hash_found  = self%htab%isthere(key)
            chash_found = self%chtab%isthere(key)
            if( hash_found .or. chash_found ) found = .true.
        else
            hash_found  = self%htab%isthere(key)
            chash_found = self%chtab%isthere(key)
            if( hash_found .or. chash_found ) found = .true.
        endif
    end function isthere

    logical function ischar( self, key )
        class(ori),       intent(in) :: self
        character(len=*), intent(in) :: key
        ischar = self%chtab%isthere(key)
    end function ischar

    logical function isstatezero( self )
        class(ori), intent(in) :: self
        isstatezero = (self%get_state() == 0)
    end function isstatezero

    !>  \brief  check wether the orientation has any typical search parameter
    logical function has_been_searched( self )
        class(ori), intent(in) :: self
        has_been_searched = .true.
        if (.not. is_zero(self%e1get())      ) return
        if (.not. is_zero(self%e2get())      ) return
        if (.not. is_zero(self%e3get())      ) return
        if (.not. is_zero(self%get('corr'))  ) return
        ! removed shifts from here as we may want to keep shifts from ini3D
        ! but restart the refinement (not initialize with the class oris)
        ! if(.not. is_zero(self%get('x'))     )return
        ! if(.not. is_zero(self%get('y'))     )return
        has_been_searched = .false.
    end function has_been_searched

    pure function pparms2str( self ) result( str )
        class(ori), intent(in)     :: self
        character(len=LONGSTRLEN)  :: str
        character(len=XLONGSTRLEN) :: tmpstr
        integer :: i, cnt
        str = repeat(' ',LONGSTRLEN)
        write(tmpstr,*)(trim(get_oriparam_flag(i)),'=',self%pparms(i),'/', i=1,N_PTCL_ORIPARAMS)
        cnt = 0
        do i=1,len_trim(tmpstr)
            if( tmpstr(i:i) == ' ' ) cycle
            cnt = cnt + 1
            str(cnt:cnt) = tmpstr(i:i)
        enddo
        do i=4,cnt
            if( str(i:i) == '/' ) str(i:i) = ' '
        enddo
        str = trim(str)
    end function pparms2str

    pure integer function pparms_strlen( self )
        class(ori), intent(in)     :: self
        character(len=XLONGSTRLEN) :: tmpstr
        integer :: i, n
        pparms_strlen = 0
        n = N_PTCL_ORIPARAMS
        write(tmpstr,*)(trim(get_oriparam_flag(i)), self%pparms(i), i=1,n)
        do i=1,len_trim(tmpstr)
            if( tmpstr(i:i) == ' ' ) cycle
            pparms_strlen = pparms_strlen + 1
        enddo
        pparms_strlen = pparms_strlen + n     ! for '=' separator
        pparms_strlen = pparms_strlen + n - 1 ! for ' ' separator
    end function pparms_strlen

    pure function ori2str( self ) result( str )
        class(ori), intent(in)        :: self
        character(len=:), allocatable :: str
        character(len=LONGSTRLEN)     :: str_pparms, str_htab
        character(len=XLONGSTRLEN)    :: str_chtab
        integer :: sz_chash, sz_hash
        sz_chash = self%chtab%size_of()
        sz_hash  = self%htab%size_of()
        if( sz_chash > 0 ) str_chtab = self%chtab%chash2str()
        if( sz_hash  > 0 ) str_htab  = self%htab%hash2str()
        if( self%is_ptcl )then
            str_pparms = self%pparms2str()
            if( sz_chash > 0 .and. sz_hash > 0 )then
                allocate( str, source=trim(str_chtab)//' '//trim(str_pparms)//' '//trim(str_htab))
            else if( sz_hash > 0 )then
                allocate( str, source=trim(str_pparms)//' '//trim(str_htab))
            else if( sz_chash > 0 )then
                allocate( str, source=trim(str_chtab)//' '//trim(str_pparms))
            else
                allocate( str, source=trim(str_pparms))
            endif
        else
            if( sz_chash > 0 .and. sz_hash > 0 )then
                allocate( str, source=trim(str_chtab)//' '//trim(str_htab))
            else if( sz_hash > 0 )then
                allocate( str, source=trim(str_htab))
            else if( sz_chash > 0 )then
                allocate( str, source=trim(str_chtab))
            endif
        endif
    end function ori2str

    pure subroutine ori2prec( self, prec )
        class(ori), intent(in)    :: self
        real,       intent(inout) :: prec(N_PTCL_ORIPARAMS)
        prec = self%pparms
    end subroutine ori2prec

    pure subroutine prec2ori( self, prec )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: prec(N_PTCL_ORIPARAMS)
        self%pparms = prec
    end subroutine prec2ori

    pure integer function ori_strlen_trim( self )
        class(ori), intent(in) :: self
        integer :: chashlen, pprmslen, hashlen
        chashlen = self%chtab%chash_strlen()
        hashlen  = self%htab%hash_strlen()
        ori_strlen_trim = chashlen + hashlen
        if( chashlen > 0 .and. hashlen > 0 )then
            ori_strlen_trim = ori_strlen_trim + 1            ! + 1 for ' ' separator
        endif
        if( self%is_ptcl )then
            pprmslen = self%pparms_strlen()
            ori_strlen_trim = ori_strlen_trim + pprmslen + 1 ! + 1 for ' ' separator
        endif
    end function ori_strlen_trim

    ! used by qsys_env, pparms omitted deliberatly
    function ori2chash( self ) result( ch )
        class(ori), intent(in) :: self
        type(chash) :: ch
        type(hash) :: h
        integer :: i
        ch = self%chtab
        h  = self%htab
        if( h%size_of() >= 1 )then
            do i=1,h%size_of()
                call ch%push(h%get_str(i), int2str(int(h%get_value_at(i))))
            end do
        endif
    end function ori2chash

    subroutine chash2ori( self, ch )
        class(ori),   intent(inout)   :: self
        class(chash), intent(in)      :: ch
        character(len=:), allocatable :: line
        line = ch%chash2str()
        call self%str2ori(line, is_ptcl=.false.)
    end subroutine chash2ori

    !>  \brief  converts euler angles into BILD Chimera readable format
    subroutine write2bild( self, fnr )
        class(ori), intent(inout) :: self
        integer,    intent(in)    :: fnr
        real :: xyz_start(3), xyz_end(3), radius
        xyz_start = self%get_normal()
        xyz_end   = 1.1 * xyz_start
        radius    = 0.05
        write(fnr,'(A,3F7.3,F7.3,F7.3,F7.3,F6.3)')'.cylinder ', xyz_start, xyz_end, radius
    end subroutine write2bild

    function get_ctfvars( self ) result( ctfvars )
        class(ori), intent(in) :: self
        type(ctfparams)        :: ctfvars
        character(len=STDLEN)  :: ctfstr, phplate
        if( self%isthere('ctf') )then
            ctfstr = self%get_static('ctf')
            ctfvars%ctfflag = CTFFLAG_YES
            select case( trim(ctfstr) )
                case('no')
                    ctfvars%ctfflag = CTFFLAG_NO
                case('yes')
                    ctfvars%ctfflag = CTFFLAG_YES
                case('flip')
                    ctfvars%ctfflag = CTFFLAG_FLIP
            end select
        endif
        ctfvars%l_phaseplate = .false.
        if( self%isthere('phaseplate') )then
            phplate = self%get_static('phaseplate')
            ctfvars%l_phaseplate = phplate .eq. 'yes'
        endif
        ctfvars%smpd    = self%get('smpd')
        ctfvars%cs      = self%get('cs')
        ctfvars%kv      = self%get('kv')
        ctfvars%fraca   = self%get('fraca')
        ctfvars%dfx     = self%get_dfx()
        ctfvars%dfy     = self%get_dfy()
        ctfvars%angast  = self%get('angast')
        ctfvars%phshift = self%get('phshift')
    end function get_ctfvars

    subroutine get_axis_angle( self, vec, angle )
        class(ori), intent(in)  :: self
        real,       intent(out) :: vec(3), angle
        real :: c1, c2, c3, s1, s2, s3, denom, euls(3), c1c2, s1s2
        euls   = self%get_euler()
        c1     = cos(euls(2)/2.)
        c2     = cos(euls(3)/2.)
        c3     = cos(euls(1)/2.)
        s1     = sin(euls(2)/2.)
        s2     = sin(euls(3)/2.)
        s3     = sin(euls(1)/2.)
        c1c2   = c1*c2
        s1s2   = s1*s2
        vec(1) =  c1c2*s3 +  s1s2*c3
        vec(2) = s1*c2*c3 + c1*s2*s3
        vec(3) = c1*s2*c3 - s1*c2*s3
        angle  = 2 * acos(c1c2*c3 - s1s2*s3)
        denom  = norm2(vec)
        if( denom < TINY )then
            vec = [1., 0., 0.]
        else
            vec = vec/denom
        endif
    end subroutine get_axis_angle

    subroutine set_ctfvars( self, ctfvars )
        class(ori),       intent(inout) :: self
        class(ctfparams), intent(in)    :: ctfvars
        call self%set('smpd', ctfvars%smpd)
        select case( ctfvars%ctfflag )
            case(CTFFLAG_NO)
                call self%set('ctf', 'no')
            case(CTFFLAG_YES)
                call self%set('ctf', 'yes')
            case(CTFFLAG_FLIP)
                call self%set('ctf', 'flip')
        end select
        call self%set('cs',    ctfvars%cs)
        call self%set('kv',    ctfvars%kv)
        call self%set('fraca', ctfvars%fraca)
        if( ctfvars%l_phaseplate )then
            call self%set('phaseplate', 'yes')
            call self%set('phshift', ctfvars%phshift)
        else
            call self%set('phaseplate', 'no')
        endif
        call self%set_dfx(ctfvars%dfx)
        call self%set_dfy(ctfvars%dfy)
        call self%set('angast', ctfvars%angast)
    end subroutine set_ctfvars

    function get_keys( self ) result (keys)
        class(ori),                 intent(in)  :: self
        type(str4arr),              allocatable :: hkeys(:)
        character(len=XLONGSTRLEN), allocatable :: keys(:)
        integer :: sz_chash, sz_hash, i
        sz_chash = self%chtab%size_of()
        sz_hash  = self%htab%size_of()
        allocate(keys(sz_chash + sz_hash))
        hkeys = self%htab%get_keys()
        do i = 1,sz_chash
            keys(i) = trim(adjustl(self%chtab%get_key(i)))
        end do
        do i=1, sz_hash
            keys(i + sz_chash) = trim(adjustl(hkeys(i)%str))
        end do
        if(allocated(hkeys)) deallocate(hkeys)
    end function get_keys

    !<  \brief  to print the rotation matrix
    subroutine print_mat( self )
        class(ori), intent(inout) :: self
        real :: euls(3), rmat(3,3)
        euls = self%get_euler()
        rmat = euler2m(euls)
        write(logfhandle,*) rmat(1,1), rmat(1,2), rmat(1,3), &
                     &      rmat(2,1), rmat(2,2), rmat(2,3), &
                     &      rmat(3,1), rmat(3,2), rmat(3,3)
    end subroutine print_mat

    subroutine print_ori( self )
        class(ori), intent(in) :: self
        character(len=KEYLEN) :: flag
        integer :: i
        if( self%is_ptcl )then
            do i=1,N_PTCL_ORIPARAMS-1,1
                flag = get_oriparam_flag(i)
                write(logfhandle,"(1X,A,A)", advance="no") trim(flag), '='
                write(logfhandle,"(A)", advance="no") trim(real2str(self%pparms(i)))
            end do
            flag = get_oriparam_flag(N_PTCL_ORIPARAMS)
            write(logfhandle,"(1X,A,A)", advance="no") trim(flag), '='
            write(logfhandle,"(A)") trim(real2str(self%pparms(N_PTCL_ORIPARAMS)))
        endif
        call self%htab%print()
        call self%chtab%print_key_val_pairs(logfhandle)
    end subroutine print_ori

    subroutine write( self, fhandle )
        class(ori), intent(inout) :: self
        integer,    intent(in)    :: fhandle
        character(len=:), allocatable :: str
        str = self%ori2str()
        if( allocated(str) ) write(fhandle,'(a)') str
    end subroutine write

    subroutine read( self, fhandle )
        class(ori), intent(inout) :: self
        integer,    intent(in)    :: fhandle
        character(len=2048) :: line
        logical :: is_ptcl
        read(fhandle,fmt='(A)') line
        is_ptcl = self%is_ptcl
        call self%str2ori(line, is_ptcl)
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
    subroutine compeuler( self1, self2, self_out )
        class(ori), intent(in)    :: self1, self2
        type(ori),  intent(inout) :: self_out
        real :: euls(3), euls1(3), euls2(3), rmat(3,3), rmat1(3,3), rmat2(3,3)
        call self_out%new_ori(self1%is_ptcl)
        euls1 = self1%get_euler()
        euls2 = self2%get_euler()
        rmat1 = euler2m(euls1)
        rmat2 = euler2m(euls2)
        rmat = matmul(rmat2,rmat1)  ! multiplication of two rotation matrices commute (n>2 not)
        ! the composed rotation matrix is constructed
        ! convert to euler
        euls = m2euler(rmat)
        call self_out%set_euler(euls)
    end subroutine compeuler

    !>  combines 3d and 2d oris and flags for ptcl mirroring
    !! \param ori3d,ori2d,o_out ori class type rotational matrices
    subroutine compose3d2d( ori3d, ori2d, o_out )
        class(ori), intent(inout) :: ori3d, ori2d, o_out
        real    :: ori3dx, ori3dy, x, y, e3, euls(3)
        real    :: R3d(3,3), R2d(3,3), R(3,3)
        logical :: mirr
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
        x   = ori2d%get('x')
        y   = ori2d%get('y')
        e3  = ori2d%e3get()
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
        call o_out%set('x', x)
        call o_out%set('y', y)
        euls(3) = e3
        call o_out%set_euler(euls)

        contains

            !>  extracts in-plane parameters from transformation matrix
            subroutine transfmat2inpls( R, psi, tx, ty )
                real,intent(inout) :: psi, tx, ty
                real,intent(in)    :: R(3,3)
                psi = rad2deg( myacos( R(1,1) ))
                if( R(1,2) < 0. ) psi = 360. - psi
                tx = R(1,3)
                ty = R(2,3)
            end subroutine transfmat2inpls

    end subroutine compose3d2d

    subroutine map3dshift22d( self, sh3d )
        class(ori), intent(inout) :: self
        real,       intent(in)    :: sh3d(3) !< 3D shift
        real :: old_x, old_y, x, y, cx, cy
        real :: u(3), v(3), shift(3), euls(3)
        real :: phi, cosphi, sinphi, rmat(3,3)
        euls   = self%get_euler()
        phi    = deg2rad(euls(3))
        cosphi = cos( phi )
        sinphi = sin( phi )
        old_x  = self%get('x')
        old_y  = self%get('y')
        ! 3D Shift rotated with respect to projdir
        rmat  = euler2m(euls)
        shift = matmul(rmat, sh3d)
        ! Projection onto xy plane
        shift = shift - dot_product(shift, zvec) * zvec
        ! e3-composed xy plane unit vectors
        u = [  cosphi,sinphi,0. ]
        v = [ -sinphi,cosphi,0. ]
        ! composed shift clockwise components
        cx =  old_x * cosphi + old_y * sinphi + dot_product(u, shift)
        cy = -old_x * sinphi + old_y * cosphi + dot_product(v, shift)
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
        mirr_mat      =  0.
        mirr_mat(1,1) =  1.
        mirr_mat(2,2) =  1.
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
        real :: euls(3), euls_mirr(3)
        euls = self%get_euler()
        call euler_mirror(euls, euls_mirr)
        call self%set_euler(euls_mirr)
    end subroutine mirror2d

    subroutine transp( self )
        class(ori), intent(inout) :: self
        real :: euls(3), rmat(3,3)
        euls = self%get_euler()
        rmat = euler2m(euls)
        rmat = transpose(rmat)
        euls = m2euler(rmat)
        call self%set_euler(euls)
    end subroutine transp

    ! GEODESIC DISTANCE METRIC

    !! \param self1,self2 ori class type rotational matrices
    pure real function geodesic_dist_frobdev( self1, self2 )
        class(ori), intent(in) :: self1, self2
        geodesic_dist_frobdev = geodesic_frobdev(self1%get_euler(), self2%get_euler())
    end function geodesic_dist_frobdev

    ! computing the angle of the difference rotation, see http://www.boris-belousov.net/2016/12/01/quat-dist/
    pure function geodesic_dist_trace( self1, self2 ) result(angle)
        class(ori), intent(in) :: self1, self2
        real :: angle, mat_diff(3, 3), arg_tr, mat1(3,3), mat2(3,3), euls1(3), euls2(3)
        euls1 = self1%get_euler()
        euls2 = self2%get_euler()
        mat1  = euler2m(euls1)
        mat2  = euler2m(euls2)
        mat_diff = matmul(mat1, transpose(mat2))
        arg_tr   = (trace(mat_diff) - 1.)/2.
        if( arg_tr < -1 .or. arg_tr > 1)then
            angle = huge(arg_tr)
        else
            angle = acos(arg_tr)
        endif
    end function geodesic_dist_trace

    ! CLASSIC DISTANCE METRICS

    !>  \brief  calculates the distance (in radians) btw two Euler angles
    !! \param self1,self2 ori class type rotational matrices
    pure function euldist( self1, self2 ) result( dist )
        class(ori), intent(in) :: self1, self2
        real :: dist, normal1(3), normal2(3)
        normal1 = self1%get_normal()
        normal2 = self2%get_normal()
        dist = myacos(dot_product(normal1, normal2))
    end function euldist

    !>  \brief  calculates the distance (in radians) btw the in-plane rotations
    !! \param self1,self2 ori class type rotaional matrices
    pure function inplrotdist( self1, self2 ) result( dist )
        class(ori), intent(in) :: self1, self2
        real :: dist, mat(2,2), u(2), x1(2), x2(2)
        u(1) = 0.
        u(2) = 1.
        call rotmat2d(self1%e3get(), mat)
        x1   = matmul(u,mat)
        call rotmat2d(self2%e3get(), mat)
        x2   = matmul(u,mat)
        dist = myacos(dot_product(x1,x2))
    end function inplrotdist

    ! DESTRUCTORS

    subroutine reset_pparms( self )
        class(ori), intent(inout) :: self
        self%pparms = 0.
    end subroutine reset_pparms

    subroutine kill_hash( self )
        class(ori), intent(inout) :: self
        call self%htab%kill
    end subroutine kill_hash

    subroutine kill_chash( self )
        class(ori), intent(inout) :: self
        call self%chtab%kill
    end subroutine kill_chash

    subroutine kill( self )
        class(ori), intent(inout) :: self
        if( self%existence )then
            self%pparms = 0.
            call self%htab%kill
            call self%chtab%kill
            self%is_ptcl   = .false.
            self%existence = .false.
        endif
    end subroutine kill

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
        absxy = sqrt(imagevec(1)**2.+imagevec(2)**2.)
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

    !>  in-plane parameters to 3x3 transformation matrix
    function make_transfmat( psi, tx, ty )result( R )
        real,intent(in) :: psi,tx,ty
        real            :: R(3,3),radpsi,cospsi,sinpsi
        radpsi = deg2rad( psi )
        cospsi = cos( radpsi )
        sinpsi = sin( radpsi )
        R      = 0.
        R(1,1) = cospsi
        R(2,2) = cospsi
        R(1,2) = sinpsi
        R(2,1) = -sinpsi
        R(1,3) = tx
        R(2,3) = ty
        R(3,3) = 1.
    end function make_transfmat

    !>  \brief  returns the rotation matrix for _ang_ degrees of rotation
    !! around x,y or z for _choice_ = _1_,_2_ or _3_
    pure function rotmat( ang, choice ) result( r )
        real, intent(in)           :: ang
        integer, intent(in)        :: choice
        real :: r(3,3)
        real :: ang_in_rad, cosang, sinang
        ang_in_rad = deg2rad(ang)
        cosang = cos( ang_in_rad )
        sinang = sin( ang_in_rad )
        if ( choice == 1 ) then
            r( 1,1 ) = 1.
            r( 1,2 ) = 0.
            r( 1,3 ) = 0.
            r( 2,1 ) = 0.
            r( 2,2 ) = cosang
            r( 2,3 ) =-sinang
            r( 3,1 ) = 0.
            r( 3,2 ) = sinang
            r( 3,3 ) = cosang
        elseif ( choice == 2 ) then
            r( 1,1 ) = cosang
            r( 1,2 ) = 0.
            r( 1,3 ) = -sinang
            r( 2,1 ) = 0.
            r( 2,2 ) = 1.
            r( 2,3 ) = 0.
            r( 3,1 ) = sinang
            r( 3,2 ) = 0.
            r( 3,3 ) = cosang
        elseif ( choice == 3 ) then
            r( 1,1 ) = cosang
            r( 1,2 ) = sinang
            r( 1,3 ) = 0.
            r( 2,1 ) = -sinang
            r( 2,2 ) = cosang
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
        real :: ang_in_rad, cosang, sinang
        ang_in_rad = deg2rad(ang)
        cosang = cos( ang_in_rad )
        sinang = sin( ang_in_rad )
        if ( choice == 1 ) then
            r( 1,1 ) = 0.
            r( 1,2 ) = 0.
            r( 1,3 ) = 0.
            r( 2,1 ) = 0.
            r( 2,2 ) = -sinang
            r( 2,3 ) = -cosang
            r( 3,1 ) = 0.
            r( 3,2 ) =  cosang
            r( 3,3 ) = -sinang
        elseif ( choice == 2 ) then
            r( 1,1 ) = -sinang
            r( 1,2 ) = 0.
            r( 1,3 ) = -cosang
            r( 2,1 ) = 0.
            r( 2,2 ) = 0.
            r( 2,3 ) = 0.
            r( 3,1 ) =  cosang
            r( 3,2 ) = 0.
            r( 3,3 ) = -sinang
        elseif ( choice == 3 ) then
            r( 1,1 ) = -sinang
            r( 1,2 ) =  cosang
            r( 1,3 ) = 0.
            r( 2,1 ) = -cosang
            r( 2,2 ) = -sinang
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

    !>  \brief  calculates the distance (in radians) btw two Euler triplets
    pure real function euler_dist( euls1, euls2 )
        real, intent(in) :: euls1(3), euls2(3)
        real             :: normal1(3), normal2(3)
        normal1 = euler_normal(euls1)
        normal2 = euler_normal(euls2)
        euler_dist = myacos(dot_product(normal1, normal2))
    end function euler_dist

    !>  \brief  calculates the normal vector of a Euler triplet
    pure function euler_normal( euls ) result( normal )
        real, intent(in) :: euls(3)
        real             :: normal(3), rmat(3,3)
        rmat   = euler2m(euls)
        normal = matmul(zvec, rmat)
    end function euler_normal

    !>  \brief  calculates the in-plane distance of a euler triplet (so just psi) in radians
    pure real function euler_inplrotdist( euls1, euls2 )
        real, intent(in) :: euls1(3), euls2(3)
        real, parameter  :: u(2) = [0.0, 1.0]
        real             :: mat(2,2), x1(2), x2(2)
        call rotmat2d(euls1(3), mat)
        x1   = matmul(u,mat)
        call rotmat2d(euls2(3), mat)
        x2   = matmul(u,mat)
        euler_inplrotdist = myacos(dot_product(x1,x2))
    end function euler_inplrotdist

    !>  \brief  is for composing Euler triplets
    pure subroutine euler_compose( euls1, euls2, euls_out )
        real, intent(in)    :: euls1(3), euls2(3)
        real, intent(inout) :: euls_out(3)
        real                :: rmat(3,3), rmat1(3,3), rmat2(3,3)
        rmat1 = euler2m(euls1)
        rmat2 = euler2m(euls2)
        ! the composed rotation matrix is constructed
        rmat = matmul(rmat2,rmat1)  ! multiplication of two rotation matrices commute (n>2 not)
        ! convert to euler
        euls_out = m2euler(rmat)
    end subroutine euler_compose

    !>  \brief  Generates mirror euler angles
    pure subroutine euler_mirror( euls, euls_out)
        real, intent(in)    :: euls(3)
        real, intent(inout) :: euls_out(3)
        real :: rmat(3,3)
        euls_out(1) = euls(1)
        euls_out(2) = 180. + euls(2)
        euls_out(3) = 180. - euls(3)
        ! the mirrored rotation matrix is constructed
        rmat = euler2m(euls_out)
        ! convert to euler
        euls_out = m2euler(rmat)
    end subroutine euler_mirror

    !>  \brief  this metric is measuring the frobenius deviation from the identity matrix .in.[0,2*sqrt(2)]
    !!          Larochelle, P.M., Murray, A.P., Angeles, J., A distance metric for finite
    !!          sets of rigid-body displacement in the polar decomposition. ASME J. Mech. Des.
    !!          129, 883--886 (2007)
    pure real function geodesic_frobdev( euls1, euls2 )
        real, intent(in) :: euls1(3), euls2(3)
        real :: Imat(3,3), sumsq, diffmat(3,3)
        real :: rmat1(3,3), rmat2(3,3)
        Imat      = 0.
        Imat(1,1) = 1.
        Imat(2,2) = 1.
        Imat(3,3) = 1.
        rmat1 = euler2m(euls1)
        rmat2 = euler2m(euls2)
        diffmat = Imat - matmul(rmat1, transpose(rmat2))
        sumsq   = sum(diffmat * diffmat)
        if( sumsq > 1e-6 )then
            geodesic_frobdev = sqrt(sumsq)
        else
            geodesic_frobdev = 0.
        endif
    end function geodesic_frobdev

    ! UNIT TEST

    !>  \brief  is the unit test for the ori class
    subroutine test_ori
    ! this test only tests the Euler part of ori,
    ! the rest is tested in the oris class
        type(ori) :: e1, e2, e3
        real :: euls(3), normal(3), mat(3,3), normal2(3), dist
        logical :: passed
        call e1%new_ori(.false.)
        call e2%new_ori(.false.)
        call e3%new_ori(.false.)
        write(logfhandle,'(a)') '**info(simple_ori_unit_test: testing all functionality'
        passed = .false.
        call e1%set_euler([1.,2.,3.])
        euls = e1%get_euler()
        if( abs(euls(1)-1.+euls(2)-2.+euls(3)-3.) < 0.001 )then
            passed = .true.
        endif
        if( .not. passed ) THROW_HARD('Euler assignment/getters corrupt!')
        passed = .false.
        call e1%e1set(99.)
        call e1%e2set(98.)
        call e1%e3set(97.)
        euls = e1%get_euler()
        if( abs(euls(1)-99.+euls(2)-98.+euls(3)-97.) < 0.001 ) passed = .true.
        if( .not. passed ) THROW_HARD('Euler e-setters corrupt!')
        passed = .false.
        call e1%set_euler([0.,0.,0.])
        normal = e1%get_normal()
        if( abs(normal(1))+abs(normal(2))+abs(normal(3)-1.)<3.*TINY ) passed = .true.
        if( .not. passed ) THROW_HARD('Euler normal derivation corrupt!')
        passed = .false.
        mat = e1%get_mat()
        if( abs(mat(1,1)-1.)+abs(mat(2,2)-1.)+abs(mat(3,3)-1.)<3.*TINY )then
            mat(1,1) = 0.
            mat(2,2) = 0.
            mat(3,3) = 0.
            if( abs(sum(mat)) < TINY ) passed = .true.
        endif
        if( .not. passed ) THROW_HARD('Euler rotation matrix derivation corrupt!')
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
        call e1%compeuler(e2, e3)         ! IFORT CANNOT DEAL WITH THE OPERATORS HERE
        euls = e3%get_euler()
        if( euls(1) < 20.0001 .and. euls(1) > 19.9999 .and.&
            euls(2) < 20.0001 .and. euls(2) > 19.9999 .and.&
            euls(3) < 20.0001 .and. euls(3) > 19.9999 ) passed = .true.
        if( .not. passed ) THROW_HARD('Euler composer corrupt!')
        write(logfhandle,'(a)') 'SIMPLE_ORI_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
        call e1%kill
        call e2%kill
        call e3%kill
    end subroutine test_ori

end module simple_ori
