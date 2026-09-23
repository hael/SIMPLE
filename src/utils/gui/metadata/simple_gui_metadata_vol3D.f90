!@descr: GUI metadata type for a single 3D volume entry (product paths + stats).
!==============================================================================
! MODULE: simple_gui_metadata_vol3D
!
! PURPOSE:
!   Extends gui_metadata_base with fields that describe one reconstructed
!   3D volume for GUI display:
!     reprojpath    — absolute path to the reprojections image
!     volpath       — absolute path to the raw reconstructed volume
!     pprocpath     — absolute path to the postprocessed volume (_pproc)
!     lppath        — absolute path to the low-pass filtered volume (_lp)
!     pprocmirrpath — absolute path to the mirrored postprocessed volume (_pproc_mirr)
!     state         — state index (multistate reconstructions)
!     box           — volume box size in pixels
!     smpd          — pixel size (Angstroms)
!     res0143       — FSC=0.143 resolution estimate (Angstroms, optional)
!     res05         — FSC=0.5 resolution estimate (Angstroms, optional)
!     pop           — particle population count contributing to this volume (optional)
!     fsc_invres    — FSC curve x-axis, 1/resolution (up to 1000 points, optional)
!     fsc_corr      — FSC curve correlation values, same length as fsc_invres (optional)
!     oridist       — orientation distribution 2D histogram, azimuth (-180..180) x
!                     elevation (-90..90), 5 degree bins (72 x 36 bins, optional)
!     oridistpath   — absolute path to the orientation-distribution histogram image (optional)
!   Provides set/get for all fields and a jsonise override that emits all
!   mandatory fields plus the optional res0143, res05, pop, fsc curve,
!   oridist histogram and oridistpath when set.
!
! DEPENDENCIES:
!   json_module, simple_defs, simple_string, simple_error,
!   simple_gui_metadata_base, simple_gui_metadata_types
!==============================================================================
module simple_gui_metadata_vol3D
use json_module,               only: json_core, json_value
use simple_defs,               only: LONGSTRLEN
use simple_error,              only: simple_exception
use simple_string,             only: string
use simple_gui_metadata_base,  only: gui_metadata_base
use simple_gui_metadata_types, only: GUI_METADATA_VOL3D_TYPE

implicit none

public :: gui_metadata_vol3D
private
#include "simple_local_flags.inc"

integer, parameter :: MAX_FSC_VOL3D    = 1000
integer, parameter :: ORIDIST_BINWIDTH = 5
integer, parameter :: ORIDIST_NBINS_X  = 360 / ORIDIST_BINWIDTH  ! azimuth,   -180..180
integer, parameter :: ORIDIST_NBINS_Y  = 180 / ORIDIST_BINWIDTH  ! elevation, -90..90

type, extends( gui_metadata_base ) :: gui_metadata_vol3D
  private
  character(len=LONGSTRLEN) :: reprojpath    = ''    ! absolute path to the reprojections image
  character(len=LONGSTRLEN) :: volpath       = ''    ! absolute path to the raw reconstructed volume
  character(len=LONGSTRLEN) :: pprocpath     = ''    ! absolute path to the postprocessed volume (_pproc)
  character(len=LONGSTRLEN) :: lppath        = ''    ! absolute path to the low-pass filtered volume (_lp)
  character(len=LONGSTRLEN) :: pprocmirrpath = ''    ! absolute path to the mirrored postprocessed volume (_pproc_mirr)
  character(len=LONGSTRLEN) :: oridistpath   = ''    ! absolute path to the orientation-distribution histogram image; empty means unset
  integer                   :: state      = 1     ! state index (multistate reconstructions)
  integer                   :: box        = 0     ! volume box size in pixels
  real                      :: smpd       = 0.0   ! pixel size (Angstroms)
  real                      :: res0143    = 0.0   ! FSC=0.143 resolution estimate (Angstroms); valid only when l_res0143
  real                      :: res05      = 0.0   ! FSC=0.5 resolution estimate (Angstroms); valid only when l_res05
  integer                   :: pop        = 0     ! particle population count; valid only when l_pop
  logical                   :: l_res0143  = .false.
  logical                   :: l_res05    = .false.
  logical                   :: l_pop      = .false.
  real                      :: fsc_invres(MAX_FSC_VOL3D) = 0.0  ! FSC curve x-axis, 1/resolution
  real                      :: fsc_corr(MAX_FSC_VOL3D)   = 0.0  ! FSC curve correlation values
  integer                   :: n_fsc      = 0     ! number of populated FSC points
  integer                   :: oridist(ORIDIST_NBINS_X, ORIDIST_NBINS_Y) = 0  ! orientation distribution histogram
  logical                   :: l_oridist  = .false.
  integer                   :: i          = 1     ! index of this entry within the current batch (IPC routing)
  integer                   :: i_max      = 1     ! total entries in the current batch (IPC routing)
contains
  procedure :: set
  procedure :: set_fsc
  procedure :: set_oridist
  procedure :: get
  procedure :: get_fsc
  procedure :: get_oridist
  procedure :: get_state
  procedure :: get_i
  procedure :: get_i_max
  procedure :: jsonise => jsonise_override
end type gui_metadata_vol3D

contains

  !---------------- setters ----------------

  ! Set all volume fields and mark the object as assigned.
  ! res0143, res05, pop and oridistpath are optional; omitting them clears the
  ! corresponding value/flag.
  ! i and i_max are required IPC routing fields (batch index / batch size).
  subroutine set( self, reprojpath, volpath, pprocpath, lppath, pprocmirrpath, state, box, smpd, i, i_max, res0143, res05, pop, oridistpath )
    class(gui_metadata_vol3D), intent(inout) :: self
    type(string),               intent(in)    :: reprojpath, volpath, pprocpath, lppath, pprocmirrpath
    integer,                    intent(in)    :: state, box
    real,                       intent(in)    :: smpd
    integer,                    intent(in)    :: i, i_max
    real,    optional,          intent(in)    :: res0143, res05
    integer, optional,          intent(in)    :: pop
    type(string), optional,     intent(in)    :: oridistpath
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned    = .true.
    self%reprojpath    = reprojpath%to_char()
    self%volpath       = volpath%to_char()
    self%pprocpath     = pprocpath%to_char()
    self%lppath        = lppath%to_char()
    self%pprocmirrpath = pprocmirrpath%to_char()
    if( present(oridistpath) ) self%oridistpath = oridistpath%to_char()
    self%state      = state
    self%box        = box
    self%smpd       = smpd
    self%i          = i
    self%i_max      = i_max
    self%l_res0143  = present(res0143)
    if( self%l_res0143 ) self%res0143 = res0143
    self%l_res05    = present(res05)
    if( self%l_res05 ) self%res05 = res05
    self%l_pop      = present(pop)
    if( self%l_pop ) self%pop = pop
  end subroutine set

  ! Store the FSC curve. invres and corr must be the same length and not
  ! exceed MAX_FSC_VOL3D points.
  subroutine set_fsc( self, invres, corr )
    class(gui_metadata_vol3D), intent(inout) :: self
    real,                       intent(in)    :: invres(:), corr(:)
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( size(invres) /= size(corr) ) THROW_HARD('invres and corr differ in size')
    if( size(invres) > MAX_FSC_VOL3D ) THROW_HARD('fsc curve exceeds maximum size')
    self%l_assigned = .true.
    self%n_fsc      = size(invres)
    self%fsc_invres(1:self%n_fsc) = invres
    self%fsc_corr(1:self%n_fsc)   = corr
  end subroutine set_fsc

  ! Store the orientation distribution histogram. hist must be shaped
  ! (ORIDIST_NBINS_X, ORIDIST_NBINS_Y) (azimuth x elevation, 5 degree bins).
  subroutine set_oridist( self, hist )
    class(gui_metadata_vol3D), intent(inout) :: self
    integer,                    intent(in)    :: hist(ORIDIST_NBINS_X, ORIDIST_NBINS_Y)
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned = .true.
    self%l_oridist  = .true.
    self%oridist    = hist
  end subroutine set_oridist

  !---------------- getters ----------------

  ! Return all fields; result is .true. if the object has been assigned.
  ! res0143, res05, pop and oridistpath are set only when the corresponding
  ! optional was supplied to set() (oridistpath returns '' otherwise).
  function get( self, reprojpath, volpath, pprocpath, lppath, pprocmirrpath, state, box, smpd, res0143, res05, pop, oridistpath ) result( l_assigned )
    class(gui_metadata_vol3D), intent(in)  :: self
    type(string),               intent(out) :: reprojpath, volpath, pprocpath, lppath, pprocmirrpath
    integer,                    intent(out) :: state, box
    real,                       intent(out) :: smpd
    real,    optional,          intent(out) :: res0143, res05
    integer, optional,          intent(out) :: pop
    type(string), optional,     intent(out) :: oridistpath
    logical                                 :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned    = self%l_assigned
    reprojpath    = trim(self%reprojpath)
    volpath       = trim(self%volpath)
    pprocpath     = trim(self%pprocpath)
    lppath        = trim(self%lppath)
    pprocmirrpath = trim(self%pprocmirrpath)
    if( present(oridistpath) ) oridistpath = trim(self%oridistpath)
    state      = self%state
    box        = self%box
    smpd       = self%smpd
    if( present(res0143) .and. self%l_res0143 ) res0143 = self%res0143
    if( present(res05)   .and. self%l_res05   ) res05   = self%res05
    if( present(pop) .and. self%l_pop ) pop = self%pop
  end function get

  ! Return the FSC curve; result is .true. if a curve has been set.
  function get_fsc( self, invres, corr ) result( l_set )
    class(gui_metadata_vol3D), intent(in)  :: self
    real,          allocatable, intent(out) :: invres(:), corr(:)
    logical                                 :: l_set
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_set = self%n_fsc > 0
    allocate(invres(self%n_fsc), corr(self%n_fsc))
    invres = self%fsc_invres(1:self%n_fsc)
    corr   = self%fsc_corr(1:self%n_fsc)
  end function get_fsc

  ! Return the orientation distribution histogram; result is .true. if set.
  function get_oridist( self, hist ) result( l_set )
    class(gui_metadata_vol3D), intent(in)  :: self
    integer,                    intent(out) :: hist(ORIDIST_NBINS_X, ORIDIST_NBINS_Y)
    logical                                 :: l_set
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_set = self%l_oridist
    hist  = self%oridist
  end function get_oridist

  ! Return the state index.
  function get_state( self ) result( state )
    class(gui_metadata_vol3D), intent(in) :: self
    integer                                :: state
    state = self%state
  end function get_state

  ! Return the index of this entry within the current IPC batch.
  function get_i( self ) result( i )
    class(gui_metadata_vol3D), intent(in) :: self
    integer                                :: i
    i = self%i
  end function get_i

  ! Return the total number of entries in the current IPC batch.
  function get_i_max( self ) result( i_max )
    class(gui_metadata_vol3D), intent(in) :: self
    integer                                :: i_max
    i_max = self%i_max
  end function get_i_max

  !---------------- serialisation ----------------

  ! Emit all mandatory fields plus optional res0143/res05/pop/fsc curve/oridist as a JSON object.
  ! Returns a null pointer when the object has not been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_vol3D), intent(inout) :: self
    type(json_core)                           :: json
    type(json_value),          pointer       :: json_ptr, json_invres_ptr, json_corr_ptr
    type(json_value),          pointer       :: json_oridist_ptr, json_row_ptr
    integer                                   :: i_fsc, ix, iy
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( .not.self%l_assigned ) then
      nullify(json_ptr)
      return
    endif
    call json%create_object(json_ptr, '')
    call json%add(json_ptr, "reprojpath",    trim(self%reprojpath))
    call json%add(json_ptr, "volpath",       trim(self%volpath))
    call json%add(json_ptr, "pprocpath",     trim(self%pprocpath))
    call json%add(json_ptr, "lppath",        trim(self%lppath))
    call json%add(json_ptr, "pprocmirrpath", trim(self%pprocmirrpath))
    if( len_trim(self%oridistpath) > 0 ) call json%add(json_ptr, "oridistpath", trim(self%oridistpath))
    call json%add(json_ptr, "state",   self%state)
    call json%add(json_ptr, "box",     self%box)
    call json%add(json_ptr, "smpd",    dble(self%smpd))
    if( self%l_res0143 ) call json%add(json_ptr, "res0143", dble(self%res0143))
    if( self%l_res05   ) call json%add(json_ptr, "res05",   dble(self%res05))
    if( self%l_pop ) call json%add(json_ptr, "pop", self%pop)
    if( self%n_fsc > 0 ) then
      call json%create_array(json_invres_ptr, 'fsc_invres')
      call json%create_array(json_corr_ptr,   'fsc_corr')
      do i_fsc = 1, self%n_fsc
        call json%add(json_invres_ptr, '', dble(self%fsc_invres(i_fsc)))
        call json%add(json_corr_ptr,   '', dble(self%fsc_corr(i_fsc)))
      end do
      call json%add(json_ptr, json_invres_ptr)
      call json%add(json_ptr, json_corr_ptr)
    endif
    if( self%l_oridist ) then
      call json%create_array(json_oridist_ptr, 'oridist')
      do ix = 1, ORIDIST_NBINS_X
        call json%create_array(json_row_ptr, '')
        do iy = 1, ORIDIST_NBINS_Y
          call json%add(json_row_ptr, '', self%oridist(ix, iy))
        end do
        call json%add(json_oridist_ptr, json_row_ptr)
      end do
      call json%add(json_ptr, json_oridist_ptr)
    endif
  end function jsonise_override

end module simple_gui_metadata_vol3D
