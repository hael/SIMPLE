!@descr: flex_pca stage protocol: the stage and fit identifiers a round carries, the typed stage request, the mod-4 half rule
!!
!! What a distributable phase asks of the rounds object is a `flex_stage_request`: which stage
!! body every worker runs, over which fit, at which global iteration and budget, and how many
!! fits the part carries. The identifiers travel to the workers through job_descr under
!! registered keys (stage, pcafit, which_iter, maxits, nfits); the request is the master-side
!! value they are taken from. The mod-4 half rule lives here because the paired engine and the
!! two-job half harness must partition a selection identically by construction.
module simple_flex_pca_stages
implicit none
private

public :: flex_stage_request
public :: PCA_STAGE_EMBED, PCA_STAGE_STATES, PCA_STAGE_PROBE, PCA_STAGE_POLISH
public :: FLEX_FIT_ALL, FLEX_FIT_A, FLEX_FIT_B
public :: flex_pca_half_of

! Stage selector carried to the worker in params%stage. 1..3 were the moment estimator's SNR /
! column / reduced-solve rounds; that path is gone and the numbering is kept so an old part file
! cannot be mistaken for a current one.
integer, parameter :: PCA_STAGE_EMBED  = 4
integer, parameter :: PCA_STAGE_STATES = 5
! One qsys round per probe EM iteration: the basis changes every iteration, so workers are
! re-launched against the master's refreshed flex_pca_pc*.mrc rather than looping locally.
integer, parameter :: PCA_STAGE_PROBE  = 6
! The joint (polish) fit after the paired merge: the same E-step pass over ALL particles against
! the merged basis under the polished namespace. A stage, not a string stamp: the stage says what
! to compute and which basis namespace to load.
integer, parameter :: PCA_STAGE_POLISH = 7

! Which fit a round serves, carried to the worker in params%pcafit beside the stage: the stage
! says what to compute, this says over which particles.
integer, parameter :: FLEX_FIT_ALL = 0
integer, parameter :: FLEX_FIT_A   = 1
integer, parameter :: FLEX_FIT_B   = 2

!> One round's control state. Constructed with keywords at the call site; every field has the
!! value a round without that piece of state always carried.
type :: flex_stage_request
    integer :: stage      = 0            !< PCA_STAGE_*
    integer :: which_iter = 0            !< the global EM iteration (0 when the stage has none)
    integer :: maxits     = 0            !< the iteration budget (0 when the stage has none)
    integer :: nfits      = 1            !< 1 single fit, 2 paired halves in one part
    character(len=64) :: label = ''      !< what the log calls the round
end type flex_stage_request

contains

    !> The ONE mod-4 split rule, shared by the two-job pcafit harness (validate_covariance_inputs)
    !! and the paired engine's driver -- so the two instruments partition the selection identically
    !! by construction. Pairing 1 (default) puts row residues {0,1} in half A; pairing 3 puts
    !! {0,3}. Both pair one even-row residue with one odd-row residue, so each half keeps both
    !! internal e/o classes under the row-alternating project eo split. Pairing 2 ({0,2}|{1,3})
    !! is eo-degenerate by construction and is REFUSED where the pairing is validated --
    !! never silently mapped here.
    pure integer function flex_pca_half_of( pind, vpair ) result( ifit )
        integer, intent(in) :: pind, vpair
        integer :: r4, a2
        a2 = 1
        if( vpair == 3 ) a2 = 3
        r4 = mod(pind, 4)
        if( r4 == 0 .or. r4 == a2 )then
            ifit = FLEX_FIT_A
        else
            ifit = FLEX_FIT_B
        endif
    end function flex_pca_half_of

end module simple_flex_pca_stages
