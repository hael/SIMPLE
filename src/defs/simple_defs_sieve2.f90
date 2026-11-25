!@descr: constants of the program particle_sieving_2: its seed, batches, greedy rounds, junk threshold and student sample, and its folder and file names
! The plan's decisions (doc/implementation_notes/planned/particle_sieving_2.md) say what each controls.
module simple_defs_sieve2
use simple_defs, only: dp, LPLIM3ITERBOUND
implicit none
private

public :: SIEVE2_NMICS, SIEVE2_SEED_NPTCLS, SIEVE2_NCLS, SIEVE2_BATCH_NPTCLS, SIEVE2_ROUND_NBATCHES, SIEVE2_ROUND_OLD_FRAC
public :: SIEVE2_TRAIN_NPTCLS, SIEVE2_JUNK_THRES, SIEVE2_LPSTOP, SIEVE2_ROUND_ITER0, SIEVE2_FREEZE_CORR, SIEVE2_FREEZE_NROUNDS
public :: SIEVE2_SEED_DIR, SIEVE2_BATCH_DIR, SIEVE2_ROUND_DIR, SIEVE2_LEARN_DIR, SIEVE2_LABELS_FBODY
public :: SIEVE2_REPORT, SIEVE2_STUDENT_TABLE

! ---- the workflow
integer,  parameter :: SIEVE2_NMICS          = 100    !< micrographs of the seed solve2D at most (nmics)
integer,  parameter :: SIEVE2_SEED_NPTCLS    = 10000  !< the seed takes no further micrograph once it holds as many
integer,  parameter :: SIEVE2_NCLS           = 100    !< classes (ncls)
integer,  parameter :: SIEVE2_BATCH_NPTCLS   = 10000  !< selected particles a batch holds at least, the tail excepted
integer,  parameter :: SIEVE2_ROUND_NBATCHES = 5      !< batches between greedy rounds
real,     parameter :: SIEVE2_ROUND_OLD_FRAC = 0.25   !< share of the earlier particles a round searches, at least the new count
integer,  parameter :: SIEVE2_TRAIN_NPTCLS   = 10000  !< the class student's training sample
real(dp), parameter :: SIEVE2_JUNK_THRES     = 0.9d0  !< a populated class is junk when the CNN's p_junk exceeds it
real,     parameter :: SIEVE2_LPSTOP         = 10.    !< the alignment's finest low-pass limit (A, lpstop)
!> the references freeze once SIEVE2_FREEZE_NROUNDS rounds running moved the class averages less than this (weighted
!! correlation with the previous ones); later rounds search the new particles only and keep the references
real,     parameter :: SIEVE2_FREEZE_CORR    = 0.99
integer,  parameter :: SIEVE2_FREEZE_NROUNDS = 2
!> the refine2D iteration of the first greedy round: past the in-plane shift search's onset (3) and the band from the
!! class FRCs (LPLIM1ITERBOUND), as solve2D's terminal pass is
integer,  parameter :: SIEVE2_ROUND_ITER0    = LPLIM3ITERBOUND + 1

! ---- the folders of the steps and their files, in the execution directory; a step's copy of the work project keeps
!      its file name, so that side files tied to the particle layout (the sigma2 store) stay valid
character(len=*), parameter :: SIEVE2_SEED_DIR     = 'seed'
character(len=*), parameter :: SIEVE2_BATCH_DIR    = 'batch_'     !< + the batch, three digits
character(len=*), parameter :: SIEVE2_ROUND_DIR    = 'round_'     !< + the round, two digits
character(len=*), parameter :: SIEVE2_LEARN_DIR    = 'learn_'     !< + the rounds done before it, two digits
character(len=*), parameter :: SIEVE2_LABELS_FBODY = 'labels_'    !< + seed or round_<nn>, .txt: the CNN's labels
character(len=*), parameter :: SIEVE2_REPORT       = 'particle_sieving_2_report.txt'
character(len=*), parameter :: SIEVE2_STUDENT_TABLE = 'student_by_class.txt' !< the student's agreement by class, ranked

end module simple_defs_sieve2
