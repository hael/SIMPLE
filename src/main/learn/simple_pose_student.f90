!@descr: the pose student: the Hermitian network with the direction grid, the band and the frozen shell weights it was trained with, its training on the teacher's rows, its prediction (the mode, the cone, the confidence) and its model file
! The model file is a key=value text header beside the weights' binary file (doc/algorithms/pose_student.md); the
! header names the weights by their file name alone, resolved against the header's directory, and a reader refuses
! weights of another network and another grid, band or geometry.
module simple_pose_student
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64
use simple_defs,         only: sp, dp, logfhandle
use simple_error,        only: simple_exception
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_fileio,       only: fopen, fclose, fileiochk, basename, get_fpath
use simple_syslib,       only: file_exists
use simple_timer,        only: tic, toc, timer_int_kind
use simple_oris,         only: oris
use simple_sym,          only: sym
use simple_hhcnn,        only: hhcnn, hhcnn_profile
use simple_polar_feats,  only: polar_feats, polar_feat
use simple_math,         only: hpsort
!$ use omp_lib,          only: omp_get_max_threads, omp_get_thread_num
implicit none

public :: pose_student, student_pred, pose_student_targets, student_train_profile
public :: NSPACE_STUDENT, KMAX_STUDENT, POSE_STUDENT_VERSION, run_all_pose_student_tests
private
#include "simple_local_flags.inc"

integer,  parameter :: NSPACE_STUDENT       = 1000  !< the student's direction grid
integer,  parameter :: KMAX_STUDENT         = 64    !< the cone's cap: a particle whose 99% needs more directions takes the dense search
integer,  parameter :: POSE_STUDENT_VERSION = 4     !< 4: the weights a file name beside the header, the threshold written to round-trip
integer,  parameter :: RECORD_CHUNK         = 1024  !< the header reader's unit: a record of any length is read in pieces of this many characters
integer,  parameter :: BATCH_MAX            = 256   !< the batch: a sixteenth of the training rows within these bounds
integer,  parameter :: BATCH_MIN            = 16
integer,  parameter :: BATCH_FRAC           = 16
real(sp), parameter :: CONE_MASS            = 0.99  !< the cone holds this share of the posterior
real(sp), parameter :: LEARNING_RATE        = 1.e-3
real(sp), parameter :: WEIGHT_DECAY         = 1.e-5

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

!> the cost record of the last training: wall seconds of the batch assembly (the shift augmentation, the folded
!! targets), the network's own record (the gradients, the optimiser steps, the held-out evaluations) and the epochs
type :: student_train_profile
    real(dp)            :: rt_assemble = 0.d0, rt_total = 0.d0
    integer             :: epochs = 0, batch_size = 0, nbatches_epoch = 0
    type(hhcnn_profile) :: net
end type student_train_profile

type :: pose_student
    private
    ! the geometry, the band, the grid
    integer  :: pftsz = 0, nrots = 0, nk = 0, kfromto(2) = 0, nspace = 0, box = 0, box_crop = 0, width = 0
    real(sp) :: lp = 0., smpd = 0., mskdiam = 0.
    logical  :: with_ctf = .true.
    character(len=:), allocatable :: pgrp
    real(sp), allocatable :: wk(:,:)             !< (kfromto(1):kfromto(2), 2) the frozen shell weights per half
    type(oris), allocatable :: grid
    integer,    allocatable :: mirr(:)           !< the grid's mirror partner per direction (the antipode's class)
    integer,    allocatable :: rep(:)            !< (nspace/2) the network's classes: the unmirrored direction of each mirror pair
    integer,    allocatable :: pair_of(:)        !< (nspace) the class of each direction, the pair it belongs to
    real(sp),   allocatable :: frame(:)          !< (nspace/2) each pair's mirror frame angle (the partner's psi, degrees)
    integer                 :: nclass = 0
    type(sym),  allocatable :: symop
    type(polar_feats)       :: feats
    type(hhcnn), allocatable :: net
    ! the teacher and the training record
    real(sp)       :: prob_beta     = 0.
    real(sp)       :: shift_aug_px  = 0.
    real(sp)       :: side_mean     = 0., side_sd = 1.  !< the side scalar's standardisation over the training rows
    integer(int64) :: seed          = 0_int64
    integer        :: nptcls_train  = 0, nheldout = 0, batch_size = 0, epochs_run = 0, best_epoch = 0, nweights = 0
    real(sp)       :: loss_train    = 0., loss_held = 0., mode_held = 0., time_train_s = 0.
    real(dp)       :: rt_assemble   = 0.d0   !< the last training's batch assembly (wall seconds)
    real(sp)       :: conf_threshold = 0., conf_fraction = 0., heldout_seeded = 0., heldout_recovered = 0. !< the confidence threshold (0: none), the declared in-cone fraction above it, the held-out fraction seeded and, of the seeded, recovered
    real(sp)       :: cone_bound    = 0.   !< the route's bound on the folded cone's radius (degrees; ruling R1: the dense search's cell), recorded with the threshold it was calibrated with
    integer(int64) :: vol1_hash     = 0_int64, selection_hash = 0_int64
    character(len=:), allocatable :: teacher_project, holdout_rule
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
    ! the model file's checks the tester drives on their own (bound so that the submodule links them)
    procedure, private, nopass :: header_weights_path
    procedure, private, nopass :: weights_header_mismatch
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
        self%pftsz    = pftsz
        self%nrots    = 2 * pftsz
        self%kfromto  = kfromto
        self%nk       = kfromto(2) - kfromto(1) + 1
        self%box      = box
        self%box_crop = box_crop
        self%smpd     = smpd
        self%mskdiam  = mskdiam
        self%pgrp     = trim(pgrp)
        self%nspace   = nspace
        self%lp       = lp
        self%with_ctf = with_ctf
        self%width    = width
        self%seed     = seed
        allocate(self%wk(kfromto(1):kfromto(2),2), source=wk)
        call self%feats%new(pftsz, kfromto, box, wk)
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
        allocate(self%net)
        call self%net%new(pftsz, self%nk, self%nclass, width, seed)
        self%nweights = self%net%get_nweights()
        self%exists   = .true.
    end subroutine new

    !> the readout's templates from the classes' reprojections (ruling R3): xrefs(pftsz, nk, 2, nspace) the tensors
    !! of the grid's references through the features (no CTF: the ring-averaged CTF squared one); done on every new
    !! network before its training
    subroutine init_readout( self, xrefs )
        class(pose_student), intent(inout) :: self
        real(sp),            intent(in)    :: xrefs(self%pftsz, self%nk, 2, self%nspace)
        real(sp) :: g(self%nk)
        if( .not. self%exists ) THROW_HARD('student not constructed; pose_student%init_readout')
        g = 1.
        call self%net%prepare_work(1)
        call self%net%init_templates(xrefs(:,:,:,self%rep), g)
    end subroutine init_readout

    !> the soft targets of a row: the posterior exp(-temp (L - L_min)) over the valid losses, normalised
    subroutine pose_student_targets( n, dists, temp, t )
        integer,  intent(in)  :: n
        real(sp), intent(in)  :: dists(n), temp
        real(sp), intent(out) :: t(n)
        real(dp) :: z, arg
        real(sp) :: dmin, invalid
        integer  :: i, imin
        invalid = 0.1 * huge(invalid)
        dmin    = huge(dmin)
        imin    = 0
        t       = 0.
        do i = 1,n
            if( ieee_is_finite(dists(i)) .and. dists(i) < invalid )then
                if( dists(i) < dmin )then
                    dmin = dists(i)
                    imin = i
                endif
            endif
        enddo
        if( imin == 0 ) return
        if( .not. (temp > 0.) )then
            t(imin) = 1.
            return
        endif
        z = 0.d0
        do i = 1,n
            if( ieee_is_finite(dists(i)) .and. dists(i) < invalid )then
                arg = real(temp,dp) * real(dists(i) - dmin,dp)
                if( arg < 80.d0 )then
                    t(i) = real(exp(-arg))
                    z    = z + real(t(i),dp)
                endif
            endif
        enddo
        if( z > 0.d0 ) t = real(real(t,dp) / z)
    end subroutine pose_student_targets

    !> two uniforms in (0,1) keyed on (seed, epoch, sample) through a Lehmer generator warmed over the key
    subroutine keyed_uniforms( seed, epoch, isample, u1, u2 )
        integer(int64), intent(in)  :: seed
        integer,        intent(in)  :: epoch, isample
        real(dp),       intent(out) :: u1, u2
        integer(int64), parameter :: M = 2147483647_int64, A = 48271_int64
        integer(int64) :: s
        integer :: k
        s = modulo(seed + 7919_int64 * int(epoch,int64) + 104729_int64 * int(isample,int64), M - 1_int64) + 1_int64
        do k = 1,4
            s = modulo(s * A, M)
            if( s == 0_int64 ) s = 1_int64
        enddo
        u1 = real(s,dp) / real(M,dp)
        s  = modulo(s * A, M)
        if( s == 0_int64 ) s = 1_int64
        u2 = real(s,dp) / real(M,dp)
    end subroutine keyed_uniforms

    !> a shift drawn uniformly in the disk of radius r (pixels)
    subroutine disk_shift( seed, epoch, isample, r, shift )
        integer(int64), intent(in)  :: seed
        integer,        intent(in)  :: epoch, isample
        real(sp),       intent(in)  :: r
        real(sp),       intent(out) :: shift(2)
        real(dp) :: u1, u2, rho, phi
        call keyed_uniforms(seed, epoch, isample, u1, u2)
        rho   = real(r,dp) * sqrt(u1)
        phi   = 6.283185307179586d0 * u2
        shift = [real(rho * cos(phi)), real(rho * sin(phi))]
    end subroutine disk_shift

    !> A permutation of 1..n keyed on (seed, epoch): a Fisher-Yates shuffle on the keyed generator
    subroutine keyed_permutation( seed, epoch, n, perm )
        integer(int64), intent(in)  :: seed
        integer,        intent(in)  :: epoch, n
        integer,        intent(out) :: perm(n)
        integer(int64), parameter :: M = 2147483647_int64, A = 48271_int64
        integer(int64) :: s
        integer :: i, j, tmp
        s = modulo(seed + 15485863_int64 * int(epoch,int64), M - 1_int64) + 1_int64
        do i = 1,n
            perm(i) = i
        enddo
        do i = n,2,-1
            s = modulo(s * A, M)
            if( s == 0_int64 ) s = 1_int64
            j = 1 + int(real(s,dp) / real(M,dp) * real(i,dp))
            if( j > i ) j = i
            tmp     = perm(i)
            perm(i) = perm(j)
            perm(j) = tmp
        enddo
    end subroutine keyed_permutation

    !> The training on the teacher's rows of n particles: x(pftsz, nk, 2, n) the packed tensors, ctf2(nk, n) the
    !! ring-averaged CTF squared, side(n), dists(nspace, n) the rows' losses and temps(n)
    !! their temperatures (the tempering factor prob_beta times the particle's own: the targets), held(n) the
    !! held-out rows (early stopping). Every epoch sees every training row shifted by a keyed draw in the disk
    !! of shift_aug_px; the held-out rows take one keyed shift each. Adam at LEARNING_RATE with WEIGHT_DECAY on batches
    !! of a sixteenth of the training rows within BATCH_MIN and BATCH_MAX, at most max_epochs epochs, stopped when
    !! the held-out loss has not improved for patience epochs; the best epoch's weights are kept
    subroutine train( self, n, x, ctf2, side, dists, temps, held, prob_beta, shift_aug_px, max_epochs, patience, &
        &nthr, log_unit )
        class(pose_student), intent(inout) :: self
        integer,             intent(in)    :: n, max_epochs, patience, nthr, log_unit
        real(sp),            intent(in)    :: x(self%pftsz, self%nk, 2, n), ctf2(self%nk, n), side(n)
        logical,             intent(in)    :: held(n)
        real(sp),            intent(in)    :: dists(self%nspace, n), temps(n), prob_beta, shift_aug_px
        type(hhcnn), allocatable :: best
        real(sp), allocatable :: xb(:,:,:,:), gb(:,:), sb(:), tb(:,:)
        integer,  allocatable :: train_inds(:), held_inds(:), perm(:)
        real(sp) :: loss, fm, loss_held, fm_held, best_loss
        real(dp) :: loss_sum, fm_sum
        integer  :: ntr, nhe, epoch, ib, nb, i0, i1, nbatch, i, j, since_best, nthr_here, bsz
        integer(timer_int_kind) :: t0, ta
        if( .not. self%exists ) THROW_HARD('student not constructed; pose_student%train')
        self%rt_assemble = 0.d0
        call self%net%reset_profile
        if( n < 2 ) THROW_HARD('too few rows to train on; pose_student%train')
        if( max_epochs < 1 .or. patience < 1 ) THROW_HARD('the training budget is at least one epoch; pose_student%train')
        self%prob_beta    = prob_beta
        self%shift_aug_px = shift_aug_px
        ntr = count(.not. held)
        nhe = count(held)
        if( ntr < 1 ) THROW_HARD('no training row; pose_student%train')
        if( nhe < 1 ) THROW_HARD('no held-out row: the early stopping needs them; pose_student%train')
        allocate(train_inds(ntr), held_inds(nhe), perm(ntr))
        j = 0
        do i = 1,n
            if( held(i) ) cycle
            j = j + 1
            train_inds(j) = i
        enddo
        j = 0
        do i = 1,n
            if( .not. held(i) ) cycle
            j = j + 1
            held_inds(j) = i
        enddo
        self%nptcls_train = ntr
        self%nheldout     = nhe
        ! the side scalar standardised by the training rows' mean and spread
        self%side_mean = real(sum(real(side(train_inds),dp)) / real(ntr,dp))
        self%side_sd   = real(sqrt(sum((real(side(train_inds),dp) - real(self%side_mean,dp))**2) / real(ntr,dp)))
        if( .not. (self%side_sd > 0.) ) self%side_sd = 1.
        nthr_here = nthr
        !$ nthr_here = omp_get_max_threads()
        call self%net%prepare_work(nthr_here)
        allocate(xb(self%pftsz, self%nk, 2, BATCH_MAX), gb(self%nk, BATCH_MAX), sb(BATCH_MAX), tb(self%nclass, BATCH_MAX))
        allocate(best)
        call best%new(self%pftsz, self%nk, self%nclass, self%width, self%seed)
        best_loss  = huge(best_loss)
        since_best = 0
        bsz        = min(BATCH_MAX, max(BATCH_MIN, ntr / BATCH_FRAC))
        nb         = (ntr + bsz - 1) / bsz
        self%batch_size = bsz
        if( log_unit > 0 ) write(log_unit,'(A,I0,A,I0,A,I0,A,I0,A)') 'training budget: ', max_epochs, ' epochs in batches of ', &
            &bsz, ' (', nb, ' per epoch), stopped after ', patience, ' epochs without a held-out improvement'
        t0 = tic()
        do epoch = 1,max_epochs
            call keyed_permutation(self%seed, epoch, ntr, perm)
            loss_sum = 0.d0
            fm_sum   = 0.d0
            do ib = 1,nb
                i0     = (ib-1) * bsz + 1
                i1     = min(ntr, ib * bsz)
                nbatch = i1 - i0 + 1
                ta = tic()
                call assemble_batch(nbatch, train_inds(perm(i0:i1)), epoch)
                self%rt_assemble = self%rt_assemble + real(toc(ta),dp)
                call self%net%train_batch(nbatch, xb, gb, sb, tb, LEARNING_RATE, WEIGHT_DECAY, loss, fm)
                loss_sum = loss_sum + real(loss,dp) * real(nbatch,dp)
                fm_sum   = fm_sum   + real(fm,dp)   * real(nbatch,dp)
            enddo
            self%loss_train = real(loss_sum / real(ntr,dp))
            self%epochs_run = epoch
            call evaluate_held(loss_held, fm_held)
            if( log_unit > 0 )then
                write(log_unit,'(A,I3,A,F9.5,A,F7.4,A,F9.5,A,F7.4,A,F9.1,A)') 'epoch ', epoch, ': train loss ', self%loss_train, &
                    &' mode ', real(fm_sum / real(ntr,dp)), ' | held-out loss ', loss_held, ' mode ', fm_held, ' | ', &
                    &real(toc(t0)), ' s'
                call flush(log_unit)
            endif
            if( loss_held < best_loss )then
                best_loss       = loss_held
                since_best      = 0
                self%best_epoch = epoch
                self%loss_held  = loss_held
                self%mode_held  = fm_held
                call best%copy_weights_from(self%net)
            else
                since_best = since_best + 1
                if( since_best >= patience ) exit
            endif
        enddo
        call self%net%copy_weights_from(best)
        self%time_train_s = real(toc(t0))
        call best%kill
        deallocate(best, xb, gb, sb, tb, train_inds, held_inds, perm)

    contains

        !> the batch's tensors: each row shifted by its keyed draw, its weights times CTF squared, its targets
        subroutine assemble_batch( nb_here, inds, epoch_here )
            integer, intent(in) :: nb_here, inds(nb_here), epoch_here
            integer  :: ii, isample, idir
            real(sp) :: sh(2), tdir(self%nspace)
            !$omp parallel do default(shared) private(ii,isample,sh,idir,tdir) schedule(static) proc_bind(close)
            do ii = 1,nb_here
                isample = inds(ii)
                if( self%shift_aug_px > 0. )then
                    call disk_shift(self%seed, epoch_here, isample, self%shift_aug_px, sh)
                    call self%feats%shift_feat(x(:,:,:,isample), sh, xb(:,:,:,ii))
                else
                    xb(:,:,:,ii) = x(:,:,:,isample)
                endif
                gb(:,ii) = ctf2(:,isample)
                sb(ii)   = (side(isample) - self%side_mean) / self%side_sd
                call pose_student_targets(self%nspace, dists(:,isample), temps(isample), tdir)
                tb(1:self%nclass,ii) = 0.
                do idir = 1,self%nspace
                    tb(self%pair_of(idir),ii) = tb(self%pair_of(idir),ii) + tdir(idir)
                enddo
            enddo
            !$omp end parallel do
        end subroutine assemble_batch

        !> the held-out rows' mean loss and mode agreement, each row at its own keyed shift (epoch 0)
        subroutine evaluate_held( loss_mean, frac_mode )
            real(sp), intent(out) :: loss_mean, frac_mode
            real(dp) :: lsum, fsum
            integer  :: jb, j0, j1, nbh, nbb
            lsum = 0.d0
            fsum = 0.d0
            nbh  = (nhe + BATCH_MAX - 1) / BATCH_MAX
            do jb = 1,nbh
                j0  = (jb-1) * BATCH_MAX + 1
                j1  = min(nhe, jb * BATCH_MAX)
                nbb = j1 - j0 + 1
                ta = tic()
                call assemble_batch(nbb, held_inds(j0:j1), 0)
                self%rt_assemble = self%rt_assemble + real(toc(ta),dp)
                call self%net%eval_batch(nbb, xb, gb, sb, tb, loss, fm)
                lsum = lsum + real(loss,dp) * real(nbb,dp)
                fsum = fsum + real(fm,dp)   * real(nbb,dp)
            enddo
            loss_mean = real(lsum / real(nhe,dp))
            frac_mode = real(fsum / real(nhe,dp))
        end subroutine evaluate_held

    end subroutine train

    !> the workspaces of a team of nthr threads for predict_feat
    subroutine prepare_predict( self, nthr )
        class(pose_student), intent(inout) :: self
        integer,             intent(in)    :: nthr
        if( .not. self%exists ) THROW_HARD('student not constructed; pose_student%prepare_predict')
        call self%net%prepare_work(nthr)
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
    !! standardisation applied here)
    subroutine predict_tensor( self, ithr, x, ctf2, side, prob )
        class(pose_student), intent(inout) :: self
        integer,             intent(in)    :: ithr
        real(sp),            intent(in)    :: x(self%pftsz, self%nk, 2), ctf2(self%nk), side
        real(sp),            intent(out)   :: prob(self%nspace)
        real(sp) :: pc(self%nclass)
        call self%net%predict(ithr, x, ctf2, (side - self%side_mean) / self%side_sd, pc)
        ! a pair's mass on its unmirrored direction; the search tests that direction and its mirror
        prob = 0.
        prob(self%rep) = pc
    end subroutine predict_tensor

    !> the readout map of thread ithr's last prediction, flattened (ruling R2's shift check reads it)
    subroutine get_features( self, ithr, f )
        class(pose_student), intent(in)  :: self
        integer,             intent(in)  :: ithr
        real(sp),            intent(out) :: f(:)
        call self%net%get_features(ithr, f)
    end subroutine get_features

    !> the readout map's size
    integer function get_nfeat( self )
        class(pose_student), intent(in) :: self
        integer :: pftsz, nk, nclasses, width
        call self%net%get_shape(pftsz, nk, nclasses, width, get_nfeat)
    end function get_nfeat

    !> the mode, the cone and the confidence on the folded sphere (ruling R1): the mass of a direction and of its
    !! mirror partner summed, the mode and the smallest set holding CONE_MASS taken over the folded directions (capped
    !! at KMAX_STUDENT), the cone's radius the set's spread about the mode where a member's distance is the smaller of
    !! its own and its partner's, across the point group; mode and cone(:) name each pair's heavier partner
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
        call self%net%get_class_rotation(ithr, ic, ang, refl)
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

    !> the header fname (key=value) and the weights beside it (<base>_weights.bin), named in the header by their file
    !! name alone
    subroutine write( self, fname, extra_keys, extra_vals )
        class(pose_student),        intent(in) :: self
        class(string),              intent(in) :: fname
        character(len=*), optional, intent(in) :: extra_keys(:), extra_vals(:)
        type(string) :: wfile, wname
        integer :: funit, io_stat, k, i
        if( .not. self%exists ) THROW_HARD('student not constructed; pose_student%write')
        wfile = weights_fname(fname)
        wname = basename(wfile)
        call self%net%write_weights(wfile)
        call fopen(funit, fname, status='replace', action='write', iostat=io_stat)
        call fileiochk('pose_student%write; open '//fname%to_char(), io_stat)
        write(funit,'(A)')        '# pose student model'
        write(funit,'(A,I0)')     'model_version=', POSE_STUDENT_VERSION
        write(funit,'(A,A)')      'weights_file=', wname%to_char()
        write(funit,'(A,I0)')     'nspace=', self%nspace
        write(funit,'(A,A)')      'pgrp=', self%pgrp
        write(funit,'(A,I0)')     'pftsz=', self%pftsz
        write(funit,'(A,I0)')     'nrots=', self%nrots
        write(funit,'(A,I0,A,I0)') 'kfromto=', self%kfromto(1), ',', self%kfromto(2)
        write(funit,'(A,F10.4)')  'lp=', self%lp
        write(funit,'(A,I0)')     'box=', self%box
        write(funit,'(A,I0)')     'box_crop=', self%box_crop
        write(funit,'(A,F10.5)')  'smpd=', self%smpd
        write(funit,'(A,F10.3)')  'mskdiam=', self%mskdiam
        write(funit,'(A,L1)')     'with_ctf=', self%with_ctf
        write(funit,'(A,I0)')     'width=', self%width
        write(funit,'(A,I0)')     'nweights=', self%nweights
        write(funit,'(A,I0)')     'seed=', self%seed
        write(funit,'(A,F8.4)')   'prob_beta=', self%prob_beta
        write(funit,'(A,F8.3)')   'shift_aug_px=', self%shift_aug_px
        write(funit,'(A,ES17.9)') 'side_mean=', self%side_mean
        write(funit,'(A,ES17.9)') 'side_sd=', self%side_sd
        write(funit,'(A,I0)')     'kmax_student=', KMAX_STUDENT
        write(funit,'(A,F6.3)')   'cone_mass=', CONE_MASS
        write(funit,'(A,I0)')     'nptcls_train=', self%nptcls_train
        write(funit,'(A,I0)')     'nheldout=', self%nheldout
        if( allocated(self%holdout_rule) ) write(funit,'(A,A)') 'holdout_rule=', self%holdout_rule
        write(funit,'(A,I0)')     'batch_size=', self%batch_size
        write(funit,'(A,I0)')     'epochs_run=', self%epochs_run
        write(funit,'(A,I0)')     'best_epoch=', self%best_epoch
        write(funit,'(A,ES14.6)') 'loss_train=', self%loss_train
        write(funit,'(A,ES14.6)') 'loss_heldout=', self%loss_held
        write(funit,'(A,F8.5)')   'mode_heldout=', self%mode_held
        write(funit,'(A,F10.2)')  'time_train_s=', self%time_train_s
        write(funit,'(A,ES17.9)') 'conf_threshold=', self%conf_threshold
        write(funit,'(A,F8.5)')   'conf_fraction=', self%conf_fraction
        write(funit,'(A,F8.5)')   'heldout_seeded=', self%heldout_seeded
        write(funit,'(A,F8.5)')   'heldout_recovered=', self%heldout_recovered
        write(funit,'(A,ES17.9)') 'cone_bound_deg=', self%cone_bound
        if( allocated(self%teacher_project) ) write(funit,'(A,A)') 'teacher_project=', self%teacher_project
        write(funit,'(A,I0)')     'teacher_selection_hash=', self%selection_hash
        write(funit,'(A,I0)')     'vol1_hash=', self%vol1_hash
        write(funit,'(A)', advance='no') 'wk_even='
        do k = self%kfromto(1),self%kfromto(2)
            if( k > self%kfromto(1) ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES17.9)', advance='no') self%wk(k,1)
        enddo
        write(funit,*)
        write(funit,'(A)', advance='no') 'wk_odd='
        do k = self%kfromto(1),self%kfromto(2)
            if( k > self%kfromto(1) ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES17.9)', advance='no') self%wk(k,2)
        enddo
        write(funit,*)
        if( present(extra_keys) )then
            do i = 1,size(extra_keys)
                write(funit,'(A,A,A)') trim(extra_keys(i)), '=', trim(extra_vals(i))
            enddo
        endif
        call fclose(funit)
        call wfile%kill
        call wname%kill
    end subroutine write

    !> the weights file's name beside the header: <base>_weights.bin
    function weights_fname( fname ) result( wfile )
        class(string), intent(in) :: fname
        type(string) :: wfile
        character(len=:), allocatable :: s
        integer :: n
        s = fname%to_char()
        n = len_trim(s)
        if( n > 4 )then
            if( s(n-3:n) == '.txt' ) s = s(1:n-4)
        endif
        wfile = trim(s)//'_weights.bin'
    end function weights_fname

    !> the weights file a header fname names (wname, a file name without a directory), in the header's directory;
    !! msg says why the name is refused, empty otherwise
    subroutine header_weights_path( fname, wname, wfile, msg )
        class(string),                 intent(in)  :: fname
        character(len=*),              intent(in)  :: wname
        type(string),                  intent(out) :: wfile
        character(len=:), allocatable, intent(out) :: msg
        type(string) :: dir
        msg = ''
        if( len_trim(wname) == 0 )then
            msg = 'the model header names no weights file; pose_student%read'
        else if( scan(trim(wname), '/\') > 0 )then
            msg = 'the model header must name its weights by a file name beside it; pose_student%read'
        else
            dir   = get_fpath(fname)
            wfile = dir%to_char()//trim(adjustl(wname))
            call dir%kill
        endif
    end subroutine header_weights_path

    !> empty when the network read from a weights file is the one a header describes (its half ring, rings, classes,
    !! width, seed and weight count), else what differs
    function weights_header_mismatch( net, pftsz, nk, nclass, width, seed, nweights ) result( msg )
        type(hhcnn),    intent(in) :: net
        integer,        intent(in) :: pftsz, nk, nclass, width, nweights
        integer(int64), intent(in) :: seed
        character(len=:), allocatable :: msg
        integer :: pftsz_w, nk_w, nclass_w, width_w, nfeat_w
        call net%get_shape(pftsz_w, nk_w, nclass_w, width_w, nfeat_w)
        msg = ''
        if( pftsz_w /= pftsz .or. nk_w /= nk .or. nclass_w /= nclass .or. width_w /= width )then
            msg = 'the weights file holds a network of another shape than its header''s; pose_student%read'
        else if( net%get_seed() /= seed )then
            msg = 'the weights file holds a network of another seed than its header''s; pose_student%read'
        else if( net%get_nweights() /= nweights )then
            msg = 'the weights file holds another weight count than its header''s; pose_student%read'
        endif
    end function weights_header_mismatch

    !> one record of a text file, of any length (the header's shell weights grow with the band)
    subroutine read_record( funit, line, io_stat )
        integer,                       intent(in)  :: funit
        character(len=:), allocatable, intent(out) :: line
        integer,                       intent(out) :: io_stat
        character(len=RECORD_CHUNK) :: chunk
        integer :: got
        line = ''
        do
            read(funit, '(A)', advance='no', size=got, iostat=io_stat) chunk
            line = line//chunk(1:got)
            if( io_stat /= 0 ) exit
        enddo
        if( is_iostat_eor(io_stat) ) io_stat = 0
        if( is_iostat_end(io_stat) .and. len(line) > 0 ) io_stat = 0
    end subroutine read_record

    !> the student from its header and the weights the header names beside it; the weights must be the network the
    !! header describes
    subroutine read( self, fname )
        class(pose_student), intent(inout) :: self
        class(string),       intent(in)    :: fname
        character(len=:), allocatable :: line, key, val, msg
        type(string) :: wfile
        integer  :: funit, io_stat, ieq, version, pftsz, nrots, kfromto(2), nspace, box, box_crop, width, k
        integer  :: nptcls_train, nheldout, batch_size, epochs_run, best_epoch, pass, nweights
        real(sp) :: lp, smpd, mskdiam, prob_beta, shift_aug_px, loss_train, loss_held, mode_held, time_train_s
        real(sp) :: conf_threshold, conf_fraction, heldout_seeded, heldout_recovered, cone_bound
        real(sp) :: side_mean, side_sd
        real(sp), allocatable :: wk(:,:)
        integer(int64) :: seed, vol1_hash, selection_hash
        logical  :: with_ctf, l_wk_even, l_wk_odd
        character(len=:), allocatable :: pgrp, holdout_rule, teacher_project, wfile_name
        if( .not. file_exists(fname) ) THROW_HARD('no model file '//fname%to_char()//'; pose_student%read')
        version = -1
        pftsz = 0; nrots = 0; kfromto = 0; nspace = 0; box = 0; box_crop = 0; width = 0
        lp = 0.; smpd = 0.; mskdiam = 0.; prob_beta = 0.; shift_aug_px = 0.; with_ctf = .true.
        side_mean = 0.; side_sd = 1.
        seed = 0_int64; vol1_hash = 0_int64; selection_hash = 0_int64
        nptcls_train = 0; nheldout = 0; batch_size = 0; epochs_run = 0; best_epoch = 0; nweights = 0
        loss_train = 0.; loss_held = 0.; mode_held = 0.; time_train_s = 0.
        conf_threshold = 0.; conf_fraction = 0.; heldout_seeded = 0.; heldout_recovered = 0.; cone_bound = 0.
        pgrp = 'c1'; holdout_rule = ''; teacher_project = ''; wfile_name = ''
        l_wk_even = .false.
        l_wk_odd  = .false.
        ! two passes: the geometry first (the weights' array needs the band), then the spectra
        do pass = 1,2
            call fopen(funit, fname, status='old', action='read', iostat=io_stat)
            call fileiochk('pose_student%read; open '//fname%to_char(), io_stat)
            do
                call read_record(funit, line, io_stat)
                if( io_stat /= 0 ) exit
                if( len_trim(line) == 0 ) cycle
                if( line(1:1) == '#' ) cycle
                ieq = index(line, '=')
                if( ieq < 2 ) cycle
                key = trim(adjustl(line(1:ieq-1)))
                val = trim(adjustl(line(ieq+1:)))
                if( pass == 1 )then
                    select case(key)
                        case('model_version');  read(val,*) version
                        case('weights_file');   wfile_name = val
                        case('nspace');         read(val,*) nspace
                        case('pgrp');           pgrp = val
                        case('pftsz');          read(val,*) pftsz
                        case('nrots');          read(val,*) nrots
                        case('kfromto');        read(val,*) kfromto(1), kfromto(2)
                        case('lp');             read(val,*) lp
                        case('box');            read(val,*) box
                        case('box_crop');       read(val,*) box_crop
                        case('smpd');           read(val,*) smpd
                        case('mskdiam');        read(val,*) mskdiam
                        case('with_ctf');       read(val,*) with_ctf
                        case('width');          read(val,*) width
                        case('nweights');       read(val,*) nweights
                        case('seed');           read(val,*) seed
                        case('prob_beta');      read(val,*) prob_beta
                        case('shift_aug_px');   read(val,*) shift_aug_px
                        case('side_mean');      read(val,*) side_mean
                        case('side_sd');        read(val,*) side_sd
                        case('nptcls_train');   read(val,*) nptcls_train
                        case('nheldout');       read(val,*) nheldout
                        case('holdout_rule');   holdout_rule = val
                        case('batch_size');     read(val,*) batch_size
                        case('epochs_run');     read(val,*) epochs_run
                        case('best_epoch');     read(val,*) best_epoch
                        case('loss_train');     read(val,*) loss_train
                        case('loss_heldout');   read(val,*) loss_held
                        case('mode_heldout');   read(val,*) mode_held
                        case('time_train_s');   read(val,*) time_train_s
                        case('conf_threshold'); read(val,*) conf_threshold
                        case('conf_fraction');  read(val,*) conf_fraction
                        case('heldout_seeded'); read(val,*) heldout_seeded
                        case('heldout_recovered'); read(val,*) heldout_recovered
                        case('cone_bound_deg'); read(val,*) cone_bound
                        case('teacher_project'); teacher_project = val
                        case('teacher_selection_hash'); read(val,*) selection_hash
                        case('vol1_hash');      read(val,*) vol1_hash
                        case default
                    end select
                else
                    select case(key)
                        case('wk_even')
                            read(val,*) (wk(k,1), k = kfromto(1),kfromto(2))
                            l_wk_even = .true.
                        case('wk_odd')
                            read(val,*) (wk(k,2), k = kfromto(1),kfromto(2))
                            l_wk_odd = .true.
                        case default
                    end select
                endif
            enddo
            call fclose(funit)
            if( pass == 1 )then
                if( version /= POSE_STUDENT_VERSION ) THROW_HARD('model_version mismatch; pose_student%read')
                if( pftsz < 2 .or. nspace < 2 .or. kfromto(1) < 1 .or. kfromto(2) < kfromto(1) .or. box < 2 ) &
                    &THROW_HARD('incomplete geometry in the model header; pose_student%read')
                allocate(wk(kfromto(1):kfromto(2),2), source=0.)
            endif
        enddo
        if( .not. (l_wk_even .and. l_wk_odd) ) THROW_HARD('the model header lacks its shell weights; pose_student%read')
        call header_weights_path(fname, wfile_name, wfile, msg)
        if( len(msg) > 0 ) THROW_HARD(msg)
        call self%new(pftsz, kfromto, box, box_crop, smpd, mskdiam, pgrp, nspace, lp, wk, with_ctf, width, seed)
        call self%net%read_weights(wfile)
        msg = weights_header_mismatch(self%net, pftsz, self%nk, self%nclass, width, seed, nweights)
        if( len(msg) > 0 ) THROW_HARD(msg)
        self%prob_beta      = prob_beta
        self%shift_aug_px   = shift_aug_px
        self%side_mean      = side_mean
        self%side_sd        = side_sd
        self%nptcls_train   = nptcls_train
        self%nheldout       = nheldout
        self%holdout_rule   = holdout_rule
        self%batch_size     = batch_size
        self%epochs_run     = epochs_run
        self%best_epoch     = best_epoch
        self%loss_train     = loss_train
        self%loss_held      = loss_held
        self%mode_held      = mode_held
        self%time_train_s   = time_train_s
        self%conf_threshold = conf_threshold
        self%conf_fraction  = conf_fraction
        self%heldout_seeded = heldout_seeded
        self%heldout_recovered = heldout_recovered
        self%cone_bound     = cone_bound
        self%teacher_project = teacher_project
        self%selection_hash = selection_hash
        self%vol1_hash      = vol1_hash
        call wfile%kill
        deallocate(wk)
    end subroutine read

    ! ===== GETTERS =====

    pure integer function get_pftsz( self )
        class(pose_student), intent(in) :: self
        get_pftsz = self%pftsz
    end function get_pftsz

    pure function get_kfromto( self ) result( kfromto )
        class(pose_student), intent(in) :: self
        integer :: kfromto(2)
        kfromto = self%kfromto
    end function get_kfromto

    pure integer function get_nspace( self )
        class(pose_student), intent(in) :: self
        get_nspace = self%nspace
    end function get_nspace

    pure real(sp) function get_lp( self )
        class(pose_student), intent(in) :: self
        get_lp = self%lp
    end function get_lp

    function get_pgrp( self ) result( pgrp )
        class(pose_student), intent(in) :: self
        character(len=:), allocatable :: pgrp
        pgrp = self%pgrp
    end function get_pgrp

    pure integer function get_box( self )
        class(pose_student), intent(in) :: self
        get_box = self%box
    end function get_box

    pure integer function get_box_crop( self )
        class(pose_student), intent(in) :: self
        get_box_crop = self%box_crop
    end function get_box_crop

    pure real(sp) function get_smpd( self )
        class(pose_student), intent(in) :: self
        get_smpd = self%smpd
    end function get_smpd

    pure real(sp) function get_mskdiam( self )
        class(pose_student), intent(in) :: self
        get_mskdiam = self%mskdiam
    end function get_mskdiam

    pure logical function get_with_ctf( self )
        class(pose_student), intent(in) :: self
        get_with_ctf = self%with_ctf
    end function get_with_ctf

    subroutine get_weights( self, wk )
        class(pose_student),   intent(in)  :: self
        real(sp), allocatable, intent(out) :: wk(:,:)
        allocate(wk(self%kfromto(1):self%kfromto(2),2), source=self%wk)
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

    !> the teacher's provenance, recorded in the model file
    !> the confidence threshold set on the held-out rows: above it the teacher's direction lies inside the cone in
    !! conf_fraction of the held-out particles; heldout_seeded the held-out fraction above it with a cone within the
    !! cap, heldout_recovered the in-cone fraction among those
    subroutine set_threshold( self, conf_threshold, conf_fraction, heldout_seeded, heldout_recovered )
        class(pose_student), intent(inout) :: self
        real(sp),            intent(in)    :: conf_threshold, conf_fraction, heldout_seeded, heldout_recovered
        self%conf_threshold    = conf_threshold
        self%conf_fraction     = conf_fraction
        self%heldout_seeded    = heldout_seeded
        self%heldout_recovered = heldout_recovered
    end subroutine set_threshold

    pure subroutine get_threshold( self, conf_threshold, conf_fraction, heldout_seeded, heldout_recovered )
        class(pose_student), intent(in)  :: self
        real(sp),            intent(out) :: conf_threshold, conf_fraction, heldout_seeded, heldout_recovered
        conf_threshold    = self%conf_threshold
        conf_fraction     = self%conf_fraction
        heldout_seeded    = self%heldout_seeded
        heldout_recovered = self%heldout_recovered
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

    subroutine set_provenance( self, teacher_project, selection_hash, vol1_hash, holdout_rule )
        class(pose_student), intent(inout) :: self
        character(len=*),    intent(in)    :: teacher_project, holdout_rule
        integer(int64),      intent(in)    :: selection_hash, vol1_hash
        self%teacher_project = trim(teacher_project)
        self%selection_hash  = selection_hash
        self%vol1_hash       = vol1_hash
        self%holdout_rule    = trim(holdout_rule)
    end subroutine set_provenance

    !> the training record for the caller's reports
    subroutine get_record( self, nptcls_train, nheldout, epochs_run, best_epoch, loss_train, loss_held, mode_held, &
        &time_train_s, nweights )
        class(pose_student), intent(in)  :: self
        integer,             intent(out) :: nptcls_train, nheldout, epochs_run, best_epoch, nweights
        real(sp),            intent(out) :: loss_train, loss_held, mode_held, time_train_s
        nptcls_train = self%nptcls_train
        nheldout     = self%nheldout
        epochs_run   = self%epochs_run
        best_epoch   = self%best_epoch
        loss_train   = self%loss_train
        loss_held    = self%loss_held
        mode_held    = self%mode_held
        time_train_s = self%time_train_s
        nweights     = self%nweights
    end subroutine get_record

    !> the cost record of the last training (the network's since its construction or the training's start)
    function get_train_profile( self ) result( prof )
        class(pose_student), intent(in) :: self
        type(student_train_profile) :: prof
        prof%rt_assemble    = self%rt_assemble
        prof%rt_total       = real(self%time_train_s,dp)
        prof%epochs         = self%epochs_run
        prof%batch_size     = self%batch_size
        prof%nbatches_epoch = 0
        if( self%batch_size > 0 ) prof%nbatches_epoch = (self%nptcls_train + self%batch_size - 1) / self%batch_size
        prof%net            = self%net%get_profile()
    end function get_train_profile

    subroutine kill( self )
        class(pose_student), intent(inout) :: self
        if( allocated(self%net) )then
            call self%net%kill
            deallocate(self%net)
        endif
        if( allocated(self%grid) )then
            call self%grid%kill
            deallocate(self%grid)
        endif
        if( allocated(self%mirr)    ) deallocate(self%mirr)
        if( allocated(self%rep)     ) deallocate(self%rep, self%pair_of, self%frame)
        self%nclass = 0
        self%cone_bound = 0.
        if( allocated(self%symop) )then
            call self%symop%kill
            deallocate(self%symop)
        endif
        call self%feats%kill
        if( allocated(self%wk)   ) deallocate(self%wk)
        if( allocated(self%pgrp) ) deallocate(self%pgrp)
        if( allocated(self%holdout_rule)    ) deallocate(self%holdout_rule)
        if( allocated(self%teacher_project) ) deallocate(self%teacher_project)
        self%pftsz    = 0
        self%nrots    = 0
        self%nk       = 0
        self%kfromto  = 0
        self%nspace   = 0
        self%box      = 0
        self%box_crop = 0
        self%width    = 0
        self%exists   = .false.
    end subroutine kill

end module simple_pose_student
