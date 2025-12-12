module simple_type_defs

! CTF flag type
enum, bind(c)
    enumerator :: ENUM_CTFFLAG = -1
    enumerator :: CTFFLAG_NO = 0, CTFFLAG_YES = 1, CTFFLAG_FLIP = 2
end enum

! Objective function type
enum, bind(c)
    enumerator :: ENUM_OBJFUN= -1
    enumerator :: OBJFUN_CC = 0, OBJFUN_EUCLID = 1
end enum

! oritype enumeration
enum, bind(c)
    enumerator :: ENUM_ORISEG  = 0
    enumerator :: MIC_SEG      = 1,  STK_SEG     = 2,  PTCL2D_SEG  = 3, CLS2D_SEG  = 4
    enumerator :: CLS3D_SEG    = 5,  PTCL3D_SEG  = 6,  OUT_SEG     = 7, OPTICS_SEG = 8
    enumerator :: PROJINFO_SEG = 11, JOBPROC_SEG = 12, COMPENV_SEG = 13
end enum
integer(kind=kind(ENUM_ORISEG)), parameter :: GENERIC_SEG = PTCL3D_SEG

! weighting methods enumeration
enum, bind(c)
    enumerator :: ENUM_WCRIT        = 0
    enumerator :: RANK_SUM_CRIT     = 1
    enumerator :: RANK_CEN_CRIT     = 2
    enumerator :: RANK_EXP_CRIT     = 3
    enumerator :: RANK_INV_CRIT     = 4
    enumerator :: CORRW_CRIT        = 5
    enumerator :: CORRW_ZSCORE_CRIT = 6
end enum

! export (to STAR) type enumeration
enum, bind(c)
    enumerator :: ENUM_STARTYPE = 0
    enumerator :: MOV_STAR = 1, MIC_STAR = 2, STK_STAR = 3, PTCL_STAR = 4, CLSAVG_STAR = 5, OUT_STAR = 6
    enumerator :: PROJINFO_STAR = 7
end enum
integer(kind=kind(ENUM_STARTYPE)), parameter :: GENERIC_STAR = PTCL_STAR

! type for CTF parameters
type ctfparams
    integer(kind(ENUM_CTFFLAG)) :: ctfflag = CTFFLAG_YES
    real    :: smpd    = 0.
    real    :: kv      = 0.
    real    :: cs      = 0.
    real    :: fraca   = 0.
    real    :: dfx     = 0.
    real    :: dfy     = 0.
    real    :: angast  = 0.
    real    :: phshift = 0.
    logical :: l_phaseplate = .false. !< image obtained with Volta phaseplate
end type ctfparams

type stats_struct
    real :: avg  = 0.
    real :: med  = 0.
    real :: sdev = 0.
    real :: maxv = 0.
    real :: minv = 0.
end type stats_struct

type inpl_struct
    real    :: e3 = 0., x = 0., y = 0., corr = 0.
    logical :: l_mirr = .false.
    integer :: find_fsc05 = 0
end type inpl_struct

type clust_inpl
    type(inpl_struct), allocatable :: params(:)
end type clust_inpl

 type clust_info
    integer          :: nptcls      = 0
    integer          :: pop         = 0
    integer          :: good_bad    = 0
    integer          :: scoreclust  = 0
    type(clust_inpl) :: algninfo
    real             :: res         = 0.
    real             :: euclid      = 0.
    real             :: homogeneity = 0.
    real             :: resscore    = 0.
    real             :: clustscore  = 0.
    real             :: jointscore  = 0.
    real             :: icescore    = 0.
end type clust_info

! type for particle reference relation in eul_prob_tab
type ptcl_ref
    integer :: pind = 0, iproj = 0, inpl = 0, istate=0
    real    :: dist = 0., x = 0., y = 0.
    logical :: has_sh = .false.
end type ptcl_ref

type lp_crop_inf
    real    :: lp=0., smpd_crop=0., scale=1., trslim=0., frc_crit=0.
    integer :: box_crop=0
    logical :: l_autoscale=.false., l_lpset=.false.
end type lp_crop_inf

type class_sample
    integer :: clsind = 0, pop = 0, nsample = 0
    integer, allocatable :: pinds(:)
    real,    allocatable :: ccs(:)
end type class_sample

end module simple_type_defs
