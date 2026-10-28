!@descr: the frame-quality feature bank (names, families, directions, descriptions), its constants, and the per-trajectory table
module single_frame_quality_types
use simple_error, only: simple_exception
implicit none
private
#include "simple_local_flags.inc"

public :: frame_quality_table, fq_feature_name, fq_feature_family, fq_feature_direction, fq_feature_description
public :: FQ_NFEATS, FQ_NLOC, FQ_SPOT_HALFWIN, FQ_CLIP, FQ_WINSZ_DEFAULT, FQ_FLOOR_K, FQ_GUARD_WEIGHT, FQ_GUARD_SEP
public :: FQ_DISAGREE, FQ_MSKDIAM_FAC, FQ_BAND_LOW_RES
public :: I_LATTICE_PEAK, I_BAND_MAX, I_PEAK_MEDIAN, I_SPOT_COUNT, I_SPOT_ARC, I_SPOT_DIAM, I_SPOT_SNR
public :: I_BAND_PRESENCE, I_NEFF

integer, parameter :: FQ_NFEATS = 9
integer, parameter :: I_LATTICE_PEAK = 1, I_BAND_MAX = 2, I_PEAK_MEDIAN = 3, I_SPOT_COUNT = 4, I_SPOT_ARC = 5
integer, parameter :: I_SPOT_DIAM = 6, I_SPOT_SNR = 7, I_BAND_PRESENCE = 8, I_NEFF = 9
integer, parameter :: FQ_NLOC          = 64     !< local maxima kept per frame, strongest first (spot_count)
integer, parameter :: FQ_SPOT_HALFWIN  = 7      !< the 15 x 15 window about the strongest spot
integer, parameter :: FQ_WINSZ_DEFAULT = 10     !< peak_median half-window, frames
real,    parameter :: FQ_CLIP          = 8.     !< clip of the normalized features (model input)
real,    parameter :: FQ_FLOOR_K       = 4.     !< nano_floor threshold and the spot_count threshold, floor sd
real,    parameter :: FQ_GUARD_WEIGHT  = 0.05   !< the floor must own at least this fraction of the frames
real,    parameter :: FQ_GUARD_SEP     = 4.     !< and lie at least this many floor sd below the upper component
real,    parameter :: FQ_DISAGREE      = 8.     !< label-disagreement margin on peak_median, floor sd
real,    parameter :: FQ_MSKDIAM_FAC   = 1.5    !< mask diameter = this times the particle diameter
real,    parameter :: FQ_BAND_LOW_RES  = 3.     !< the lattice band runs from this resolution (A) to Nyquist

type :: fq_feature_def
    character(len=16)  :: name        = ''
    character(len=16)  :: direction   = ''
    character(len=8)   :: family      = ''
    character(len=128) :: description = ''
end type fq_feature_def

type(fq_feature_def), parameter :: FEATURE_DEFS(FQ_NFEATS) = [ &
    fq_feature_def('lattice_peak',  'higher_is_better', 'peak',    &
        'log10 of the largest smoothed band power over its ring mean'), &
    fq_feature_def('band_max',      'higher_is_better', 'peak',    &
        'log10 of the largest smoothed band power'), &
    fq_feature_def('peak_median',   'higher_is_better', 'peak',    &
        'running median of lattice_peak over 2 winsz + 1 frames'), &
    fq_feature_def('spot_count',    'higher_is_better', 'lattice', &
        'local maxima of the smoothed ratio above the floor threshold'), &
    fq_feature_def('spot_arc',      'lower_is_better',  'lattice', &
        'log10 of the strongest spot extent across over along its radius'), &
    fq_feature_def('spot_diam',     'none',             'shape',   &
        'projected-sphere diameter (A) fitted to the strongest spot'), &
    fq_feature_def('spot_snr',      'higher_is_better', 'shape',   &
        'log10 of the fitted spot amplitude over its background'), &
    fq_feature_def('band_presence', 'higher_is_better', 'energy',  &
        'log10 of band-passed variance in the central disc over the annulus'), &
    fq_feature_def('neff',          'none',             'denoise', &
        'ln of the denoiser Kish population of the frame') ]

!> one trajectory: raw features, the floor, the normalized features, the model's verdict; one row per frame
type :: frame_quality_table
    integer :: n = 0
    integer :: box = 0
    real    :: smpd = 0.
    integer :: winsz = FQ_WINSZ_DEFAULT
    integer :: sh_lo = 0, sh_hi = 0
    logical :: has_floor = .false.
    real    :: floor_mean = 0., floor_sd = 0., floor_weight = 0., upper_mean = 0., upper_sd = 0.
    real    :: sf_floor_mean = 0., sf_floor_sd = 0. !< the floor of the single-frame lattice_peak (wider than the 5-median's)
    real    :: bm_floor_mean = 0., bm_floor_sd = 0.
    real    :: spot_threshold = 0.          !< log10 ratio: floor mean + FQ_FLOOR_K floor sd
    real    :: white_noise = 0.             !< lattice_peak of white noise for this box, band and mask
    real    :: diam = 0., mskdiam = 0.      !< particle diameter (A) from the accepted frames, mask diameter
    logical, allocatable :: invalid(:)      !< invalid pixels: hard rejected
    logical, allocatable :: floor_frames(:) !< owned by the floor component
    logical, allocatable :: accept(:)
    integer, allocatable :: spot_shell(:)   !< shell of the strongest spot (0: none)
    real,    allocatable :: raw(:,:)        !< (n, FQ_NFEATS) the feature values
    real,    allocatable :: z(:,:)          !< (n, FQ_NFEATS) normalized, unclipped
    real,    allocatable :: prob(:)
    real,    allocatable :: locmax(:,:)     !< (n, FQ_NLOC) smoothed ratios of the strongest local maxima
contains
    procedure :: new
    procedure :: zclip
    procedure :: kill
end type frame_quality_table

contains

    subroutine new( self, n, box, smpd )
        class(frame_quality_table), intent(inout) :: self
        integer,                    intent(in)    :: n, box
        real,                       intent(in)    :: smpd
        call self%kill
        if( n < 1 ) THROW_HARD('no frames; frame_quality_table%new')
        self%n    = n
        self%box  = box
        self%smpd = smpd
        allocate(self%invalid(n), self%floor_frames(n), self%accept(n), source=.false.)
        allocate(self%raw(n,FQ_NFEATS), self%z(n,FQ_NFEATS), self%locmax(n,FQ_NLOC), source=0.)
        allocate(self%prob(n), source=0.)
        allocate(self%spot_shell(n), source=0)
    end subroutine new

    !> the model input of frame i: normalized features clipped to FQ_CLIP
    function zclip( self, i ) result( zc )
        class(frame_quality_table), intent(in) :: self
        integer,                    intent(in) :: i
        real :: zc(FQ_NFEATS)
        zc = max(-FQ_CLIP, min(FQ_CLIP, self%z(i,:)))
    end function zclip

    subroutine kill( self )
        class(frame_quality_table), intent(inout) :: self
        if( allocated(self%invalid)      ) deallocate(self%invalid)
        if( allocated(self%floor_frames) ) deallocate(self%floor_frames)
        if( allocated(self%accept)       ) deallocate(self%accept)
        if( allocated(self%raw)          ) deallocate(self%raw)
        if( allocated(self%z)            ) deallocate(self%z)
        if( allocated(self%prob)         ) deallocate(self%prob)
        if( allocated(self%locmax)       ) deallocate(self%locmax)
        if( allocated(self%spot_shell)   ) deallocate(self%spot_shell)
        self%n = 0
        self%has_floor = .false.
    end subroutine kill

    function fq_feature_name( i ) result( name )
        integer, intent(in) :: i
        character(len=16) :: name
        if( i < 1 .or. i > FQ_NFEATS ) THROW_HARD('invalid feature index; fq_feature_name')
        name = FEATURE_DEFS(i)%name
    end function fq_feature_name

    function fq_feature_family( i ) result( family )
        integer, intent(in) :: i
        character(len=8) :: family
        if( i < 1 .or. i > FQ_NFEATS ) THROW_HARD('invalid feature index; fq_feature_family')
        family = FEATURE_DEFS(i)%family
    end function fq_feature_family

    function fq_feature_direction( i ) result( direction )
        integer, intent(in) :: i
        character(len=16) :: direction
        if( i < 1 .or. i > FQ_NFEATS ) THROW_HARD('invalid feature index; fq_feature_direction')
        direction = FEATURE_DEFS(i)%direction
    end function fq_feature_direction

    function fq_feature_description( i ) result( description )
        integer, intent(in) :: i
        character(len=128) :: description
        if( i < 1 .or. i > FQ_NFEATS ) THROW_HARD('invalid feature index; fq_feature_description')
        description = FEATURE_DEFS(i)%description
    end function fq_feature_description

end module single_frame_quality_types
