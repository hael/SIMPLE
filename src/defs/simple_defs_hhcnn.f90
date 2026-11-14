!@descr: constants of the polar students on the HHCNN (the pose and class students): their grids and candidate caps, the teacher's rows and the training designs, the held-out calibration, the shift registration, and the learn, predict and validation programs
! The network's own constants (CNN_*) and the trainer's step settings stay private to simple_hhcnn and
! simple_polar_student. doc/algorithms/pose_student.md and doc/algorithms/class_student.md say what each controls.
module simple_defs_hhcnn
implicit none
private

public :: NSPACE_STUDENT, KMAX_STUDENT, KMAX_CLASS
public :: STUDENT_KMAX_ROW, STUDENT_WIDTHS, STUDENT_WIDEST_MIN_ROWS, STUDENT_MAX_EPOCHS, STUDENT_PATIENCE_EPOCHS
public :: CLASS_STUDENT_NDESIGNS
public :: STUDENT_SHIFT_AUG_NSIG, STUDENT_SHIFT_SIGMA, STUDENT_TRS, CLASS_STUDENT_TEACHER_ITER
public :: STUDENT_HOLDOUT_EVERY, STUDENT_CONF_FRACTION, STUDENT_CONF_NONE, STUDENT_NBOOT, STUDENT_DESIGN_QUANTILE
public :: STUDENT_SMPD_TOL, STUDENT_MSKDIAM_TOL
public :: STUDENT_SHIFT_REG_HALFWIDTH, STUDENT_SHIFT_REG_MAXITS
public :: POSE_SHIFT_CHECK_NPTCLS, POSE_SHIFT_CHECK_MAGS, POSE_SHIFT_CHECK_NDIR, POSE_SHIFT_CHECK_NGROUPS
public :: STUDENT_LP_VALIDATE, POSE_AGREE_DEG, POSE_TWO_STEPS_DEG
public :: CLASS_DENSE_SHIFT_SIGMA, CLASS_DENSE_PROB_BETA, CLASS_MINPOP_FRC
public :: STUDENT_TRAIN_PROJ, STUDENT_ALL_PROJ, STUDENT_UNSEEN_PROJ
public :: POSE_STUDENT_MODEL_FNAME, POSE_STUDENT_PRED_TABLE, VALIDATE_POSE_STUDENT_REPORT
public :: CLASS_STUDENT_MODEL_FNAME, CLASS_STUDENT_PRED_TABLE, CLASS_STUDENT_CAVGS, VALIDATE_CLASS_STUDENT_REPORT
public :: VALIDATE_CLASS_STUDENT_DETAILS

! ---- the students' classes and candidates
integer, parameter :: NSPACE_STUDENT = 1000 !< the pose student's direction grid
integer, parameter :: KMAX_STUDENT   = 64   !< the pose student's cone cap: a particle whose 99% needs more directions takes the dense search
integer, parameter :: KMAX_CLASS     = 16   !< the class student's candidate cap: a particle whose 99% needs more classes is searched densely

! ---- the teacher's rows and the training
integer, parameter :: STUDENT_KMAX_ROW        = 256   !< the row record's kept entries, the tail's mass beside them
integer, parameter :: STUDENT_WIDTHS(3)       = [16, 32, 64] !< the designs, cheapest first: the first block's channels
integer, parameter :: STUDENT_WIDEST_MIN_ROWS = 10000 !< the widest design is tried from this many training rows: below, it rarely beats the next and costs two and a half times it
integer, parameter :: CLASS_STUDENT_NDESIGNS  = 1     !< the class student's designs tried, cheapest first: in 2D a wider one has not paid for its time (doc/algorithms/class_student.md); raise it to restore the choice
integer, parameter :: STUDENT_MAX_EPOCHS      = 30    !< the training budget in epochs (the best held-out epoch comes at 3 to 7 on the fixed sets)
integer, parameter :: STUDENT_PATIENCE_EPOCHS = 3     !< epochs without a held-out improvement before the training stops (improvements came at most 2 epochs apart)
real,    parameter :: STUDENT_SHIFT_AUG_NSIG  = 3.    !< the augmentation's disk: this many prior widths (the dynamic centring's disk)
real,    parameter :: STUDENT_SHIFT_SIGMA     = 2.    !< the training's shift-prior width (A): the teacher's particles are centred (10 A taught the 2D student nothing)
real,    parameter :: STUDENT_TRS             = 5.    !< the teacher's shift refinement box (pixels)
integer, parameter :: CLASS_STUDENT_TEACHER_ITER = 1  !< the class averages prepared as a first iteration's: every class polarized, none re-centred

! ---- the held-out calibration and the design choice
integer, parameter :: STUDENT_HOLDOUT_EVERY = 10      !< every tenth of the selected particles' micrographs (stacks) is held out, every tenth selected particle when they lie on fewer stacks
real,    parameter :: STUDENT_CONF_FRACTION = 0.95    !< the confidence threshold: above it the teacher's answer lies among the candidates for this fraction of the held-out particles
real,    parameter :: STUDENT_CONF_NONE     = 1.00001 !< the threshold when no prefix qualifies: above every confidence, so that the student seeds nobody
integer, parameter :: STUDENT_NBOOT         = 1000    !< resamples by micrograph of the paired difference between two designs
real,    parameter :: STUDENT_DESIGN_QUANTILE = 0.05    !< a wider design is kept when this quantile of its resampled recovery gain is above zero

! ---- a model's geometry against a project's
real, parameter :: STUDENT_SMPD_TOL    = 1.e-4 !< a model's sampling (A/pixel) equals the project's within the header's five decimals
real, parameter :: STUDENT_MSKDIAM_TOL = 1.e-2 !< a model's mask diameter (A) equals the project's within the header's three decimals

! ---- the shift registration of a prediction
real,    parameter :: STUDENT_SHIFT_REG_HALFWIDTH = 25. !< its box: +/- this (A) about the stored shift
integer, parameter :: STUDENT_SHIFT_REG_MAXITS    = 50  !< its L-BFGS-B evaluations after the coarse grid

! ---- predict_poses' shift check (ruling R2 of the pose student)
integer, parameter :: POSE_SHIFT_CHECK_NPTCLS  = 256           !< the first selected particles predicted again
real,    parameter :: POSE_SHIFT_CHECK_MAGS(3) = [2., 5., 10.] !< shifted through the calculator's phase ramps by these pixels
integer, parameter :: POSE_SHIFT_CHECK_NDIR    = 4             !< in this many directions (0, 90, 180, 270 degrees)
integer, parameter :: POSE_SHIFT_CHECK_NGROUPS = 1             !< the readout map as one feature group (ruling R3's readout replaced the pooled features)

! ---- the validation programs
real,    parameter :: STUDENT_LP_VALIDATE     = 15.  !< the students' band in validation (A), the pose student's size-curve study's
real,    parameter :: POSE_AGREE_DEG          = 6.   !< one step of the pose student's grid
real,    parameter :: POSE_TWO_STEPS_DEG      = 12.  !< two steps of the pose student's grid
real,    parameter :: CLASS_DENSE_SHIFT_SIGMA = 10.  !< the dense reference pass's prior (A): center_ptcls's, whose 3-sigma disk covers the picking error
real,    parameter :: CLASS_DENSE_PROB_BETA   = 100. !< the dense pass draws at this factor: each row's mode, the teacher's class, as the held-out rows are judged
integer, parameter :: CLASS_MINPOP_FRC        = 10   !< a class's average is compared when both routes hold this many particles

! ---- file names
character(len=*), parameter :: STUDENT_TRAIN_PROJ            = 'training_set.simple'
character(len=*), parameter :: STUDENT_ALL_PROJ              = 'all_particles.simple'
character(len=*), parameter :: STUDENT_UNSEEN_PROJ           = 'unseen_particles.simple'
character(len=*), parameter :: POSE_STUDENT_MODEL_FNAME      = 'pose_student.txt'
character(len=*), parameter :: POSE_STUDENT_PRED_TABLE       = 'pose_student_pred.txt'
character(len=*), parameter :: VALIDATE_POSE_STUDENT_REPORT  = 'validate_pose_student_report.txt'
character(len=*), parameter :: CLASS_STUDENT_MODEL_FNAME     = 'class_student.txt'
character(len=*), parameter :: CLASS_STUDENT_PRED_TABLE      = 'class_student_pred.txt'
character(len=*), parameter :: CLASS_STUDENT_CAVGS           = 'cavgs_student.mrc'
character(len=*), parameter :: VALIDATE_CLASS_STUDENT_REPORT = 'validate_class_student_report.txt'
character(len=*), parameter :: VALIDATE_CLASS_STUDENT_DETAILS = 'validate_class_student_details.txt'

end module simple_defs_hhcnn
