!@descr: the pose student: the direction grid's mirror pairs as the classes of the polar student's core, the folded cone, the read-out pose, the point group, and the model file's grid keys
! The core (simple_polar_student) holds the network, the features, the band, the training and the model file's common
! part; this class adds what is 3D: the grid, its mirror pairs and their frames, the point group, the cone on the folded
! sphere, the full pose a readout gives, the route's bound and the map the rows came from.
module simple_pose_student
use, intrinsic :: iso_fortran_env, only: int64
use simple_defs,          only: sp, dp, logfhandle
use simple_error,         only: simple_exception
use simple_string,        only: string
use simple_oris,          only: oris
use simple_sym,           only: sym
use simple_polar_feats,   only: polar_feat
use simple_polar_student, only: polar_student, polar_student_header, student_train_profile, polar_student_targets
use simple_math,          only: hpsort
use simple_defs_hhcnn,    only: KMAX_STUDENT
implicit none

public :: pose_student, student_pred
public :: POSE_STUDENT_VERSION, run_all_pose_student_tests
private
#include "simple_local_flags.inc"

integer,  parameter :: POSE_STUDENT_VERSION = 4     !< 4: the weights a file name beside the header, the threshold written to round-trip
real(sp), parameter :: CONE_MASS            = 0.99  !< the cone holds this share of the posterior
character(len=*), parameter :: MODEL_TITLE  = 'pose student model'

!> the record of one prediction on the folded sphere (ruling R1: a direction and its mirror partner, the antipode,
!! are one class to the Friedel-symmetric student, so their masses are summed before the cone is read): the mode,
!! the cone (the smallest set of folded directions holding CONE_MASS, at most KMAX_STUDENT), its mass and its
!! angular radius about the mode (the spread on the folded sphere). mode and cone(:) name the partner of larger mass
type :: student_pred
    integer  :: mode   = 0     !< the direction of largest folded mass (its heavier partner)
    real(sp) :: pmax   = 0.    !< the folded mass of the mode
    integer  :: k99    = 0     !< folded directions holding CONE_MASS
    integer  :: ncone  = 0     !< the cone's size (k99 capped)
    real(sp) :: conf   = 0.    !< the folded mass inside the cone
    real(sp) :: radius = 0.    !< the cone's angular radius about the mode on the folded sphere (degrees)
    integer  :: cone(KMAX_STUDENT) = 0
end type student_pred

type :: pose_student
    private
    type(polar_student), allocatable :: core   !< the network, the features, the band, the training, the common model file
    integer :: nspace = 0
    character(len=:), allocatable :: pgrp
    type(oris), allocatable :: grid
    integer,    allocatable :: mirr(:)           !< the grid's mirror partner per direction (the antipode's class)
    integer,    allocatable :: rep(:)            !< (nspace/2) the network's classes: the unmirrored direction of each mirror pair
    integer,    allocatable :: pair_of(:)        !< (nspace) the class of each direction, the pair it belongs to
    real(sp),   allocatable :: frame(:)          !< (nspace/2) each pair's mirror frame angle (the partner's psi, degrees)
    integer                 :: nclass = 0
    type(sym),  allocatable :: symop
    real(sp)       :: cone_bound = 0.   !< the route's bound on the folded cone's radius (degrees; ruling R1: the dense search's cell), recorded with the threshold it was calibrated with
    integer(int64) :: vol1_hash  = 0_int64
    logical :: exists = .false.
  contains
    procedure :: new
    procedure :: init_readout
    procedure :: train
    procedure :: prepare_predict
    procedure :: predict_feat
    procedure :: predict_tensor
    procedure :: get_features
    procedure :: get_nfeat
    procedure :: cone_of
    procedure :: pose_of
    procedure :: fold_index
    procedure :: mirror_of
    procedure, private :: heavier
    procedure :: write
    procedure :: read
    procedure :: get_pftsz
    procedure :: get_kfromto
    procedure :: get_nspace
    procedure :: get_lp
    procedure :: get_pgrp
    procedure :: get_box
    procedure :: get_box_crop
    procedure :: get_smpd
    procedure :: get_mskdiam
    procedure :: get_with_ctf
    procedure :: get_weights
    procedure :: get_grid
    procedure :: grid_euler
    procedure :: grid_closest
    procedure :: grid_spacing
    procedure :: angle_between
    procedure :: get_record
    procedure :: get_train_profile
    procedure :: set_provenance
    procedure :: set_threshold
    procedure :: get_threshold
    procedure :: set_cone_bound
    procedure :: get_cone_bound
    procedure :: kill
end type pose_student

interface
    module subroutine run_all_pose_student_tests()
    end subroutine run_all_pose_student_tests
end interface

contains

    !> A student for training: the calculator's geometry (pftsz over the half ring, the band kfromto, the box the
    !! shifts are measured in and the cropped box), the band's low-pass lp, the sampling, the mask, the point group,
    !! the grid's size, the shell weights wk the rows were whitened with (frozen), whether the CTF is applied, the
    !! network's width, and the seed its weights start from
    subroutine new( self, pftsz, kfromto, box, box_crop, smpd, mskdiam, pgrp, nspace, lp, wk, with_ctf, width, seed )
        class(pose_student), intent(inout) :: self
        integer,             intent(in)    :: pftsz, kfromto(2), box, box_crop, nspace, width
        real(sp),            intent(in)    :: smpd, mskdiam, lp
        character(len=*),    intent(in)    :: pgrp
        real(sp),            intent(in)    :: wk(kfromto(1):kfromto(2),2)
        logical,             intent(in)    :: with_ctf
        integer(int64),      intent(in)    :: seed
        integer :: k
        integer, allocatable :: half_partners(:)
        real,    allocatable :: half_frames(:)
        call self%kill
        self%pgrp   = trim(pgrp)
        self%nspace = nspace
        allocate(self%grid, self%symop)
        call self%symop%new(trim(pgrp))
        call self%grid%new(nspace, is_ptcl=.false.)
        call self%symop%build_refspiral(self%grid)
        ! the spiral's mirror pairs (the second half mirrors the first): the folded sphere's classes
        allocate(self%mirr(nspace), source=0)
        do k = 1,nspace
            self%mirr(k) = nint(self%grid%get(k, 'mirr'))
            if( self%mirr(k) < 1 .or. self%mirr(k) > nspace .or. self%mirr(k) == k ) THROW_HARD('the grid carries no mirror partner; pose_student%new')
        enddo
        ! the network classifies the mirror pairs (the spiral's unmirrored half), invariant to the mirror; the search
        ! tells a direction from its mirror (prob_neigh's cone about both)
        call self%grid%get_mirror_pairs(half_partners, half_frames)
        self%nclass = size(half_partners)
        allocate(self%rep(self%nclass), self%pair_of(nspace), self%frame(self%nclass))
        do k = 1,self%nclass
            if( self%mirr(k) /= half_partners(k) ) THROW_HARD('the grid''s mirror pairs disagree; pose_student%new')
            self%rep(k)                  = k
            self%pair_of(k)              = k
            self%pair_of(half_partners(k)) = k
            self%frame(k)                = half_frames(k)
        enddo
        allocate(self%core)
        call self%core%new(pftsz, kfromto, box, box_crop, smpd, mskdiam, lp, wk, with_ctf, self%nclass, .true., width, seed)
        self%exists = .true.
    end subroutine new

    !> the readout's templates from the classes' reprojections (ruling R3): xrefs(pftsz_net, nk, 2, nspace) the tensors
    !! of the grid's references through the features (no CTF: the ring-averaged CTF squared one); done on every new
    !! network before its training
    subroutine init_readout( self, xrefs )
        class(pose_student), intent(inout) :: self
        real(sp),            intent(in)    :: xrefs(:,:,:,:)
        if( .not. self%exists ) THROW_HARD('student not constructed; pose_student%init_readout')
        call self%core%init_templates(xrefs(:,:,:,self%rep))
    end subroutine init_readout

    !> The training on the teacher's rows of n particles (the core's training, simple_polar_student): x(pftsz_net, nk,
    !! 2, n) the packed tensors, ctf2(nk, n), side(n), dists(nspace, n) the rows' losses and temps(n) their
    !! temperatures, whose posterior over the grid folded onto the mirror pairs is each row's soft target
    subroutine train( self, n, x, ctf2, side, dists, temps, held, prob_beta, shift_aug_px, max_epochs, patience, &
        &nthr, log_unit )
        class(pose_student), intent(inout) :: self
        integer,             intent(in)    :: n, max_epochs, patience, nthr, log_unit
        real(sp),            intent(in)    :: x(:,:,:,:), ctf2(:,:), side(n)
        logical,             intent(in)    :: held(n)
        real(sp),            intent(in)    :: dists(self%nspace, n), temps(n), prob_beta, shift_aug_px
        real(sp), allocatable :: targets(:,:)
        real(sp) :: tdir(self%nspace)
        integer  :: isample, idir
        if( .not. self%exists ) THROW_HARD('student not constructed; pose_student%train')
        allocate(targets(self%nclass, n))
        !$omp parallel do default(shared) private(isample,idir,tdir) schedule(static) proc_bind(close)
        do isample = 1,n
            call polar_student_targets(self%nspace, dists(:,isample), temps(isample), tdir)
            targets(:,isample) = 0.
            do idir = 1,self%nspace
                targets(self%pair_of(idir),isample) = targets(self%pair_of(idir),isample) + tdir(idir)
            enddo
        enddo
        !$omp end parallel do
        call self%core%train(n, x, ctf2, side, targets, held, prob_beta, shift_aug_px, max_epochs, patience, nthr, log_unit)
        deallocate(targets)
    end subroutine train

    !> the workspaces of a team of nthr threads for predict_feat
    subroutine prepare_predict( self, nthr )
        class(pose_student), intent(inout) :: self
        integer,             intent(in)    :: nthr
        if( .not. self%exists ) THROW_HARD('student not constructed; pose_student%prepare_predict')
        call self%core%prepare_predict(nthr)
    end subroutine prepare_predict

    !> the posterior over the grid of one feature (built with this student's geometry and weights) on thread ithr
    subroutine predict_feat( self, ithr, feat, prob )
        class(pose_student), intent(inout) :: self
        integer,             intent(in)    :: ithr
        type(polar_feat),    intent(in)    :: feat
        real(sp),            intent(out)   :: prob(self%nspace)
        call self%predict_tensor(ithr, feat%x, feat%ctf2, feat%logpow, prob)
    end subroutine predict_feat

    !> the posterior of a packed tensor with its ring-averaged CTF squared and raw side scalar (the training rows'
    !! standardisation applied by the core)
    subroutine predict_tensor( self, ithr, x, ctf2, side, prob )
        class(pose_student), intent(inout) :: self
        integer,             intent(in)    :: ithr
        real(sp),            intent(in)    :: x(:,:,:), ctf2(:), side
        real(sp),            intent(out)   :: prob(self%nspace)
        real(sp) :: pc(self%nclass)
        call self%core%predict_tensor(ithr, x, ctf2, side, pc)
        ! a pair's mass on its unmirrored direction; the search tests that direction and its mirror
        prob = 0.
        prob(self%rep) = pc
    end subroutine predict_tensor

    !> the readout map of thread ithr's last prediction, flattened (ruling R2's shift check reads it)
    subroutine get_features( self, ithr, f )
        class(pose_student), intent(in)  :: self
        integer,             intent(in)  :: ithr
        real(sp),            intent(out) :: f(:)
        call self%core%get_features(ithr, f)
    end subroutine get_features

    !> the readout map's size
    integer function get_nfeat( self )
        class(pose_student), intent(in) :: self
        get_nfeat = self%core%get_nfeat()
    end function get_nfeat

    subroutine cone_of( self, prob, pred )
        class(pose_student), intent(in)  :: self
        real(sp),            intent(in)  :: prob(self%nspace)
        type(student_pred),  intent(out) :: pred
        real(sp) :: wf(self%nspace), euls_mode(3), d
        integer  :: inds(self%nspace), i, j, k, p, npair, imode
        real(dp) :: csum
        wf = 0.
        do i = 1,self%nspace
            p     = self%fold_index(i)
            wf(p) = wf(p) + prob(i)
        enddo
        npair = 0
        do i = 1,self%nspace
            if( self%fold_index(i) /= i ) cycle
            npair       = npair + 1
            inds(npair) = i
        enddo
        call sort_pairs
        imode     = inds(npair)
        pred%mode = self%heavier(imode, prob)
        pred%pmax = wf(imode)
        csum      = 0.d0
        pred%k99  = 0
        do i = npair,1,-1
            csum     = csum + real(wf(inds(i)),dp)
            pred%k99 = pred%k99 + 1
            if( csum >= real(CONE_MASS,dp) ) exit
        enddo
        pred%ncone = min(pred%k99, KMAX_STUDENT)
        pred%cone  = 0
        pred%conf  = 0.
        do j = 1,pred%ncone
            k            = inds(npair - j + 1)
            pred%cone(j) = self%heavier(k, prob)
            pred%conf    = pred%conf + wf(k)
        enddo
        euls_mode   = self%grid%get_euler(pred%mode)
        pred%radius = 0.
        do j = 2,pred%ncone
            d = min(self%angle_between(pred%cone(j), euls_mode), self%angle_between(self%mirr(pred%cone(j)), euls_mode))
            pred%radius = max(pred%radius, d)
        enddo

    contains

        !> the folded masses of the pairs' representatives sorted ascending with their indices
        subroutine sort_pairs
            real(sp) :: wsort(npair)
            integer  :: ii
            do ii = 1,npair
                wsort(ii) = wf(inds(ii))
            enddo
            call hpsort(wsort, inds(1:npair))
        end subroutine sort_pairs

    end subroutine cone_of

    !> the full pose the readout of thread ithr's last prediction gives a direction's pair, without any search: the
    !! readout turns the sample's map by ang from the template (the unmirrored direction's reprojection at e3 = 0), so
    !! the unmirrored direction with e3 = 360 - ang, as the search's rotation index; scored by the reflected map, the
    !! mirror partner with e3 = ang + the pair's frame angle (the reflection that mirror_ref_bank applies)
    subroutine pose_of( self, ithr, iproj, euls, iproj_pose )
        class(pose_student), intent(in)  :: self
        integer,             intent(in)  :: ithr, iproj
        real(sp),            intent(out) :: euls(3)
        integer,             intent(out) :: iproj_pose   !< the grid direction of the pose (the pair's member)
        real(sp) :: ang
        integer  :: ic
        logical  :: refl
        ic = self%pair_of(iproj)
        call self%core%class_rotation(ithr, ic, ang, refl)
        if( refl )then
            iproj_pose = self%mirr(self%rep(ic))
            euls       = self%grid%get_euler(iproj_pose)
            euls(3)    = modulo(ang + self%frame(ic), 360.)
        else
            iproj_pose = self%rep(ic)
            euls       = self%grid%get_euler(iproj_pose)
            euls(3)    = modulo(360. - ang, 360.)
        endif
    end subroutine pose_of

    !> the folded sphere's class of a direction: the smaller index of the direction and its mirror partner
    pure integer function fold_index( self, iproj )
        class(pose_student), intent(in) :: self
        integer,             intent(in) :: iproj
        fold_index = min(iproj, self%mirr(iproj))
    end function fold_index

    !> the mirror partner of a direction on the grid (the antipode's class under the point group)
    pure integer function mirror_of( self, iproj )
        class(pose_student), intent(in) :: self
        integer,             intent(in) :: iproj
        mirror_of = self%mirr(iproj)
    end function mirror_of

    !> of a direction and its mirror partner, the one of larger mass
    pure integer function heavier( self, iproj, prob )
        class(pose_student), intent(in) :: self
        integer,             intent(in) :: iproj
        real(sp),            intent(in) :: prob(self%nspace)
        heavier = iproj
        if( prob(self%mirr(iproj)) > prob(iproj) ) heavier = self%mirr(iproj)
    end function heavier

    ! ===== THE MODEL FILE =====

    !> the header fname (key=value) and the weights beside it (the core's model file), with the grid's keys and the
    !! caller's after the core's
    subroutine write( self, fname, extra_keys, extra_vals )
        class(pose_student),        intent(in) :: self
        class(string),              intent(in) :: fname
        character(len=*), optional, intent(in) :: extra_keys(:), extra_vals(:)
        integer, parameter :: NOWN = 6
        character(len=64)   :: keys(NOWN)
        character(len=1024) :: vals(NOWN)
        character(len=64),   allocatable :: allkeys(:)
        character(len=1024), allocatable :: allvals(:)
        integer :: nextra
        if( .not. self%exists ) THROW_HARD('student not constructed; pose_student%write')
        keys(1) = 'nspace';         write(vals(1),'(I0)')    self%nspace
        keys(2) = 'pgrp';           vals(2) = self%pgrp
        keys(3) = 'kmax_student';   write(vals(3),'(I0)')    KMAX_STUDENT
        keys(4) = 'cone_mass';      write(vals(4),'(F6.3)')  CONE_MASS
        keys(5) = 'cone_bound_deg'; write(vals(5),'(ES17.9)') self%cone_bound
        keys(6) = 'vol1_hash';      write(vals(6),'(I0)')    self%vol1_hash
        nextra = 0
        if( present(extra_keys) ) nextra = size(extra_keys)
        allocate(allkeys(NOWN + nextra), allvals(NOWN + nextra))
        allkeys(1:NOWN) = keys
        allvals(1:NOWN) = vals
        if( nextra > 0 )then
            allkeys(NOWN+1:) = extra_keys
            allvals(NOWN+1:) = extra_vals
        endif
        call self%core%write(fname, MODEL_TITLE, POSE_STUDENT_VERSION, allkeys, allvals)
        deallocate(allkeys, allvals)
    end subroutine write

    !> the student from its header and the weights the header names beside it; the weights must be the network the
    !! header describes
    subroutine read( self, fname )
        class(pose_student), intent(inout) :: self
        class(string),       intent(in)    :: fname
        type(polar_student_header) :: hdr
        character(len=:), allocatable :: pgrp, key, val
        real(sp)       :: cone_bound
        integer(int64) :: vol1_hash
        integer        :: nspace, i
        call hdr%read(fname)
        if( hdr%version /= POSE_STUDENT_VERSION ) THROW_HARD('model_version mismatch; pose_student%read')
        nspace     = 0
        pgrp       = 'c1'
        cone_bound = 0.
        vol1_hash  = 0_int64
        do i = 1,hdr%nextra
            key = hdr%keys(i)%to_char()
            val = hdr%vals(i)%to_char()
            select case(key)
                case('nspace');         read(val,*) nspace
                case('pgrp');           pgrp = val
                case('cone_bound_deg'); read(val,*) cone_bound
                case('vol1_hash');      read(val,*) vol1_hash
                case default
            end select
        enddo
        if( nspace < 2 ) THROW_HARD('incomplete geometry in the model header; pose_student%read')
        call self%new(hdr%pftsz, hdr%kfromto, hdr%box, hdr%box_crop, hdr%smpd, hdr%mskdiam, pgrp, nspace, hdr%lp, &
            &hdr%wk, hdr%with_ctf, hdr%width, hdr%seed)
        call self%core%load(hdr, fname)
        self%cone_bound = cone_bound
        self%vol1_hash  = vol1_hash
        call hdr%kill
    end subroutine read

    ! ===== GETTERS =====

    pure integer function get_pftsz( self )
        class(pose_student), intent(in) :: self
        get_pftsz = self%core%get_pftsz()
    end function get_pftsz

    pure function get_kfromto( self ) result( kfromto )
        class(pose_student), intent(in) :: self
        integer :: kfromto(2)
        kfromto = self%core%get_kfromto()
    end function get_kfromto

    pure integer function get_nspace( self )
        class(pose_student), intent(in) :: self
        get_nspace = self%nspace
    end function get_nspace

    pure real(sp) function get_lp( self )
        class(pose_student), intent(in) :: self
        get_lp = self%core%get_lp()
    end function get_lp

    function get_pgrp( self ) result( pgrp )
        class(pose_student), intent(in) :: self
        character(len=:), allocatable :: pgrp
        pgrp = self%pgrp
    end function get_pgrp

    pure integer function get_box( self )
        class(pose_student), intent(in) :: self
        get_box = self%core%get_box()
    end function get_box

    pure integer function get_box_crop( self )
        class(pose_student), intent(in) :: self
        get_box_crop = self%core%get_box_crop()
    end function get_box_crop

    pure real(sp) function get_smpd( self )
        class(pose_student), intent(in) :: self
        get_smpd = self%core%get_smpd()
    end function get_smpd

    pure real(sp) function get_mskdiam( self )
        class(pose_student), intent(in) :: self
        get_mskdiam = self%core%get_mskdiam()
    end function get_mskdiam

    pure logical function get_with_ctf( self )
        class(pose_student), intent(in) :: self
        get_with_ctf = self%core%get_with_ctf()
    end function get_with_ctf

    subroutine get_weights( self, wk )
        class(pose_student),   intent(in)  :: self
        real(sp), allocatable, intent(out) :: wk(:,:)
        call self%core%get_weights(wk)
    end subroutine get_weights

    !> the direction grid
    subroutine get_grid( self, grid )
        class(pose_student), intent(in)    :: self
        type(oris),          intent(inout) :: grid
        grid = self%grid
    end subroutine get_grid

    function grid_euler( self, iproj ) result( euls )
        class(pose_student), intent(in) :: self
        integer,             intent(in) :: iproj
        real(sp) :: euls(3)
        euls = self%grid%get_euler(iproj)
    end function grid_euler

    !> the grid direction nearest an orientation
    integer function grid_closest( self, o )
        use simple_ori, only: ori
        class(pose_student), intent(in) :: self
        class(ori),          intent(in) :: o
        grid_closest = self%grid%find_closest_proj(o)
    end function grid_closest

    !> the grid's spacing: the angle (degrees) from the first direction to its nearest neighbour across the point group
    real(sp) function grid_spacing( self )
        class(pose_student), intent(in) :: self
        real(sp) :: euls1(3), d
        integer  :: j
        euls1 = self%grid%get_euler(1)
        grid_spacing = 180.
        do j = 2,self%nspace
            d = self%angle_between(j, euls1)
            if( d > 1.e-3 ) grid_spacing = min(grid_spacing, d)
        enddo
    end function grid_spacing

    !> the angle (degrees) between a grid direction and an orientation across the point group
    real(sp) function angle_between( self, iproj, euls_ref )
        class(pose_student), intent(in) :: self
        integer,             intent(in) :: iproj
        real(sp),            intent(in) :: euls_ref(3)
        real(sp) :: euls(3), euls_sym(3), inpldist
        euls = self%grid%get_euler(iproj)
        call self%symop%sym_dists(euls, euls_ref, euls_sym, angle_between, inpldist)
    end function angle_between

    !> the confidence threshold set on the held-out rows: above it the teacher's direction lies inside the cone in
    !! conf_fraction of the held-out particles; heldout_seeded the held-out fraction above it with a cone within the
    !! cap, heldout_recovered the in-cone fraction among those
    subroutine set_threshold( self, conf_threshold, conf_fraction, heldout_seeded, heldout_recovered )
        class(pose_student), intent(inout) :: self
        real(sp),            intent(in)    :: conf_threshold, conf_fraction, heldout_seeded, heldout_recovered
        call self%core%set_threshold(conf_threshold, conf_fraction, heldout_seeded, heldout_recovered)
    end subroutine set_threshold

    pure subroutine get_threshold( self, conf_threshold, conf_fraction, heldout_seeded, heldout_recovered )
        class(pose_student), intent(in)  :: self
        real(sp),            intent(out) :: conf_threshold, conf_fraction, heldout_seeded, heldout_recovered
        call self%core%get_threshold(conf_threshold, conf_fraction, heldout_seeded, heldout_recovered)
    end subroutine get_threshold

    !> the route's bound on the folded cone's radius (degrees), set with the threshold it is calibrated with
    subroutine set_cone_bound( self, cone_bound )
        class(pose_student), intent(inout) :: self
        real(sp),            intent(in)    :: cone_bound
        self%cone_bound = cone_bound
    end subroutine set_cone_bound

    pure real(sp) function get_cone_bound( self )
        class(pose_student), intent(in) :: self
        get_cone_bound = self%cone_bound
    end function get_cone_bound

    !> the teacher's provenance, recorded in the model file: the core's and the map the rows came from
    subroutine set_provenance( self, teacher_project, selection_hash, vol1_hash, holdout_rule )
        class(pose_student), intent(inout) :: self
        character(len=*),    intent(in)    :: teacher_project, holdout_rule
        integer(int64),      intent(in)    :: selection_hash, vol1_hash
        call self%core%set_provenance(teacher_project, selection_hash, holdout_rule)
        self%vol1_hash = vol1_hash
    end subroutine set_provenance

    !> the training record for the caller's reports
    subroutine get_record( self, nptcls_train, nheldout, epochs_run, best_epoch, loss_train, loss_held, mode_held, &
        &time_train_s, nweights )
        class(pose_student), intent(in)  :: self
        integer,             intent(out) :: nptcls_train, nheldout, epochs_run, best_epoch, nweights
        real(sp),            intent(out) :: loss_train, loss_held, mode_held, time_train_s
        call self%core%get_record(nptcls_train, nheldout, epochs_run, best_epoch, loss_train, loss_held, mode_held, &
            &time_train_s, nweights)
    end subroutine get_record

    !> the cost record of the last training
    function get_train_profile( self ) result( prof )
        class(pose_student), intent(in) :: self
        type(student_train_profile) :: prof
        prof = self%core%get_train_profile()
    end function get_train_profile

    subroutine kill( self )
        class(pose_student), intent(inout) :: self
        if( allocated(self%core) )then
            call self%core%kill
            deallocate(self%core)
        endif
        if( allocated(self%grid) )then
            call self%grid%kill
            deallocate(self%grid)
        endif
        if( allocated(self%mirr)    ) deallocate(self%mirr)
        if( allocated(self%rep)     ) deallocate(self%rep, self%pair_of, self%frame)
        if( allocated(self%symop) )then
            call self%symop%kill
            deallocate(self%symop)
        endif
        if( allocated(self%pgrp) ) deallocate(self%pgrp)
        self%nclass     = 0
        self%nspace     = 0
        self%cone_bound = 0.
        self%vol1_hash  = 0_int64
        self%exists     = .false.
    end subroutine kill

end module simple_pose_student
