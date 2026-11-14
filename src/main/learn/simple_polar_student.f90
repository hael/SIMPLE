!@descr: the polar student's core shared by the pose and the class students: the Hermitian network over the features of the calculator's half ring, the band, the frozen shell weights, the training on soft class targets, the prediction, the threshold record, the provenance and the model file's common header and weights
! A student's classes are its owner's (the pose student's grid's mirror pairs, the class student's class averages); the
! core sees nclass classes and their templates. The model file is a key=value text header beside the weights' binary
! file (doc/algorithms/pose_student.md): the core writes and reads its own keys, the owner adds and parses its own.
module simple_polar_student
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64
use simple_defs,         only: sp, dp
use simple_error,        only: simple_exception
use simple_string,       only: string
use simple_fileio,       only: fopen, fclose, fileiochk, basename, get_fpath
use simple_syslib,       only: file_exists
use simple_timer,        only: tic, toc, timer_int_kind
use simple_hhcnn,        only: hhcnn, hhcnn_profile
use simple_polar_feats,  only: polar_feats, polar_feat
!$ use omp_lib,          only: omp_get_max_threads
implicit none

public :: polar_student, polar_student_header, student_train_profile, polar_student_targets
public :: run_all_polar_student_tests
private
#include "simple_local_flags.inc"

integer,  parameter :: RECORD_CHUNK  = 1024  !< the header reader's unit: a record of any length is read in pieces of this many characters
integer,  parameter :: BATCH_MAX     = 256   !< the batch: a sixteenth of the training rows within these bounds
integer,  parameter :: BATCH_MIN     = 16
integer,  parameter :: BATCH_FRAC    = 16
integer,  parameter :: MAX_EXTRA     = 64    !< the header keys the core leaves to its owner
real(sp), parameter :: LEARNING_RATE = 1.e-3
real(sp), parameter :: WEIGHT_DECAY  = 1.e-5

!> the cost record of the last training: wall seconds of the batch assembly (the shift augmentation, the targets),
!! the network's own record (the gradients, the optimiser steps, the held-out evaluations) and the epochs
type :: student_train_profile
    real(dp)            :: rt_assemble = 0.d0, rt_total = 0.d0
    integer             :: epochs = 0, batch_size = 0, nbatches_epoch = 0
    type(hhcnn_profile) :: net
end type student_train_profile

!> a model header's values: the core's keys, the shell weights, and the keys it does not know, for its owner
type :: polar_student_header
    integer  :: version = -1, pftsz = 0, nrots = 0, kfromto(2) = 0, box = 0, box_crop = 0, width = 0, nweights = 0
    integer  :: nptcls_train = 0, nheldout = 0, batch_size = 0, epochs_run = 0, best_epoch = 0, nextra = 0
    real(sp) :: lp = 0., smpd = 0., mskdiam = 0., prob_beta = 0., shift_aug_px = 0., side_mean = 0., side_sd = 1.
    real(sp) :: loss_train = 0., loss_held = 0., mode_held = 0., time_train_s = 0.
    real(sp) :: conf_threshold = 0., conf_fraction = 0., heldout_seeded = 0., heldout_recovered = 0.
    integer(int64) :: seed = 0_int64, selection_hash = 0_int64
    logical  :: with_ctf = .true.
    character(len=:), allocatable :: weights_file, holdout_rule, teacher_project
    real(sp),     allocatable :: wk(:,:)
    type(string), allocatable :: keys(:), vals(:)
  contains
    procedure :: read => header_read
    procedure :: kill => header_kill
end type polar_student_header

type :: polar_student
    private
    ! the geometry, the band, the classes
    integer  :: pftsz = 0, pftsz_net = 0, nrots = 0, nk = 0, kfromto(2) = 0, box = 0, box_crop = 0, width = 0
    integer  :: nclass = 0, nweights = 0
    real(sp) :: lp = 0., smpd = 0., mskdiam = 0.
    logical  :: with_ctf = .true.
    logical  :: l_refl   = .true.               !< the classes score the reflected map too (the pose student's mirror pairs)
    real(sp), allocatable :: wk(:,:)             !< (kfromto(1):kfromto(2), 2) the frozen shell weights per half
    type(polar_feats)        :: feats
    type(hhcnn), allocatable :: net
    ! the teacher and the training record
    real(sp)       :: prob_beta     = 0.
    real(sp)       :: shift_aug_px  = 0.
    real(sp)       :: side_mean     = 0., side_sd = 1.  !< the side scalar's standardisation over the training rows
    integer(int64) :: seed          = 0_int64
    integer        :: nptcls_train  = 0, nheldout = 0, batch_size = 0, epochs_run = 0, best_epoch = 0
    real(sp)       :: loss_train    = 0., loss_held = 0., mode_held = 0., time_train_s = 0.
    real(dp)       :: rt_assemble   = 0.d0   !< the last training's batch assembly (wall seconds)
    real(sp)       :: conf_threshold = 0., conf_fraction = 0., heldout_seeded = 0., heldout_recovered = 0. !< the confidence threshold (0: none), the declared fraction above it, the held-out fraction seeded and, of the seeded, recovered
    integer(int64) :: selection_hash = 0_int64
    character(len=:), allocatable :: teacher_project, holdout_rule
    logical :: exists = .false.
  contains
    procedure :: new
    procedure :: init_templates
    procedure :: train
    procedure :: prepare_predict
    procedure :: predict_tensor
    procedure :: get_features
    procedure :: get_nfeat
    procedure :: class_rotation
    procedure :: write
    procedure :: load
    procedure :: get_pftsz
    procedure :: get_pftsz_net
    procedure :: get_kfromto
    procedure :: get_nk
    procedure :: get_nclass
    procedure :: get_lp
    procedure :: get_box
    procedure :: get_box_crop
    procedure :: get_smpd
    procedure :: get_mskdiam
    procedure :: get_with_ctf
    procedure :: get_weights
    procedure :: get_seed
    procedure :: get_record
    procedure :: get_train_profile
    procedure :: shift_feat
    procedure :: set_provenance
    procedure :: set_threshold
    procedure :: get_threshold
    procedure :: kill
    ! the model file's checks the tester drives on their own (bound so that the submodule links them)
    procedure, private, nopass :: header_weights_path
    procedure, private, nopass :: weights_header_mismatch
end type polar_student

interface
    module subroutine run_all_polar_student_tests()
    end subroutine run_all_polar_student_tests
end interface

contains

    !> A core for training: the calculator's geometry (half ring, band, box and crop), sampling, mask, low-pass, the
    !! frozen shell weights wk, the CTF treatment, the classes, whether they score the reflected map, the network's
    !! width and seed
    subroutine new( self, pftsz, kfromto, box, box_crop, smpd, mskdiam, lp, wk, with_ctf, nclass, reflect, width, seed )
        class(polar_student), intent(inout) :: self
        integer,              intent(in)    :: pftsz, kfromto(2), box, box_crop, nclass, width
        real(sp),             intent(in)    :: smpd, mskdiam, lp
        real(sp),             intent(in)    :: wk(kfromto(1):kfromto(2),2)
        logical,              intent(in)    :: with_ctf, reflect
        integer(int64),       intent(in)    :: seed
        call self%kill
        self%pftsz    = pftsz
        self%kfromto  = kfromto
        self%nk       = kfromto(2) - kfromto(1) + 1
        self%box      = box
        self%box_crop = box_crop
        self%smpd     = smpd
        self%mskdiam  = mskdiam
        self%lp       = lp
        self%with_ctf = with_ctf
        self%nclass   = nclass
        self%l_refl   = reflect
        self%width    = width
        self%seed     = seed
        allocate(self%wk(kfromto(1):kfromto(2),2), source=wk)
        call self%feats%new(pftsz, kfromto, box, wk)
        ! the network works on the features' half ring: the calculator's capped at the band's sampling limit
        self%pftsz_net = self%feats%get_pftsz()
        self%nrots     = 2 * self%pftsz_net
        allocate(self%net)
        call self%net%new(self%pftsz_net, self%nk, nclass, width, seed, reflect)
        self%nweights = self%net%get_nweights()
        self%exists   = .true.
    end subroutine new

    !> the readout's templates from the classes' references: xtmpl(pftsz_net, nk, 2, nclass) their tensors through the
    !! features (no CTF: the ring-averaged CTF squared one); done on every new network before its training
    subroutine init_templates( self, xtmpl )
        class(polar_student), intent(inout) :: self
        real(sp),             intent(in)    :: xtmpl(self%pftsz_net, self%nk, 2, self%nclass)
        real(sp) :: g(self%nk)
        if( .not. self%exists ) THROW_HARD('student not constructed; polar_student%init_templates')
        g = 1.
        call self%net%prepare_work(1)
        call self%net%init_templates(xtmpl, g)
    end subroutine init_templates

    !> the soft targets of a row: the posterior exp(-temp (L - L_min)) over the valid losses, normalised
    subroutine polar_student_targets( n, dists, temp, t )
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
    end subroutine polar_student_targets

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

    !> The training on n rows' tensors and soft class targets, held(n) the early stopping's rows: every epoch shifts each
    !! training row by a keyed draw in the disk of shift_aug_px; Adam on batches until the held-out loss has not
    !! improved for patience epochs (at most max_epochs), the best epoch's weights kept
    subroutine train( self, n, x, ctf2, side, targets, held, prob_beta, shift_aug_px, max_epochs, patience, nthr, &
        &log_unit )
        class(polar_student), intent(inout) :: self
        integer,              intent(in)    :: n, max_epochs, patience, nthr, log_unit
        real(sp),             intent(in)    :: x(self%pftsz_net, self%nk, 2, n), ctf2(self%nk, n), side(n)
        real(sp),             intent(in)    :: targets(self%nclass, n), prob_beta, shift_aug_px
        logical,              intent(in)    :: held(n)
        type(hhcnn), allocatable :: best
        real(sp), allocatable :: xb(:,:,:,:), gb(:,:), sb(:), tb(:,:)
        integer,  allocatable :: train_inds(:), held_inds(:), perm(:)
        real(sp) :: loss, fm, loss_held, fm_held, best_loss
        real(dp) :: loss_sum, fm_sum
        integer  :: ntr, nhe, epoch, ib, nb, i0, i1, nbatch, i, j, since_best, nthr_here, bsz
        integer(timer_int_kind) :: t0, ta
        if( .not. self%exists ) THROW_HARD('student not constructed; polar_student%train')
        self%rt_assemble = 0.d0
        call self%net%reset_profile
        if( n < 2 ) THROW_HARD('too few rows to train on; polar_student%train')
        if( max_epochs < 1 .or. patience < 1 ) THROW_HARD('the training budget is at least one epoch; polar_student%train')
        self%prob_beta    = prob_beta
        self%shift_aug_px = shift_aug_px
        ntr = count(.not. held)
        nhe = count(held)
        if( ntr < 1 ) THROW_HARD('no training row; polar_student%train')
        if( nhe < 1 ) THROW_HARD('no held-out row: the early stopping needs them; polar_student%train')
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
        allocate(xb(self%pftsz_net, self%nk, 2, BATCH_MAX), gb(self%nk, BATCH_MAX), sb(BATCH_MAX), tb(self%nclass, BATCH_MAX))
        allocate(best)
        call best%new(self%pftsz_net, self%nk, self%nclass, self%width, self%seed, self%l_refl)
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
            integer  :: ii, isample
            real(sp) :: sh(2)
            !$omp parallel do default(shared) private(ii,isample,sh) schedule(static) proc_bind(close)
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
                tb(1:self%nclass,ii) = targets(:,isample)
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

    !> the workspaces of a team of nthr threads for predict_tensor
    subroutine prepare_predict( self, nthr )
        class(polar_student), intent(inout) :: self
        integer,              intent(in)    :: nthr
        if( .not. self%exists ) THROW_HARD('student not constructed; polar_student%prepare_predict')
        call self%net%prepare_work(nthr)
    end subroutine prepare_predict

    !> the class posterior of a packed tensor with its ring-averaged CTF squared and raw side scalar (the training
    !! rows' standardisation applied here), on thread ithr
    subroutine predict_tensor( self, ithr, x, ctf2, side, pc )
        class(polar_student), intent(inout) :: self
        integer,              intent(in)    :: ithr
        real(sp),             intent(in)    :: x(self%pftsz_net, self%nk, 2), ctf2(self%nk), side
        real(sp),             intent(out)   :: pc(self%nclass)
        call self%net%predict(ithr, x, ctf2, (side - self%side_mean) / self%side_sd, pc)
    end subroutine predict_tensor

    !> the readout map of thread ithr's last prediction, flattened
    subroutine get_features( self, ithr, f )
        class(polar_student), intent(in)  :: self
        integer,              intent(in)  :: ithr
        real(sp),             intent(out) :: f(:)
        call self%net%get_features(ithr, f)
    end subroutine get_features

    !> the readout map's size
    integer function get_nfeat( self )
        class(polar_student), intent(in) :: self
        integer :: pftsz, nk, nclasses, width
        call self%net%get_shape(pftsz, nk, nclasses, width, get_nfeat)
    end function get_nfeat

    !> the readout's turn (degrees) of class ic's template in thread ithr's last prediction, and whether the reflected
    !! map scored it
    subroutine class_rotation( self, ithr, ic, ang, refl )
        class(polar_student), intent(in)  :: self
        integer,              intent(in)  :: ithr, ic
        real(sp),             intent(out) :: ang
        logical,              intent(out) :: refl
        call self%net%get_class_rotation(ithr, ic, ang, refl)
    end subroutine class_rotation

    !> the packed tensor shifted by shift (pixels of the box) through the features' phase ramp
    subroutine shift_feat( self, x, shift, x_sh )
        class(polar_student), intent(in)  :: self
        real(sp),             intent(in)  :: x(self%pftsz_net, self%nk, 2), shift(2)
        real(sp),             intent(out) :: x_sh(self%pftsz_net, self%nk, 2)
        call self%feats%shift_feat(x, shift, x_sh)
    end subroutine shift_feat

    ! ===== THE MODEL FILE =====

    !> the header fname (key=value) and the weights beside it (<base>_weights.bin), named in the header by their file
    !! name alone: the title line, the version, the core's keys and shell weights, then the owner's keys
    subroutine write( self, fname, title, version, extra_keys, extra_vals )
        class(polar_student),       intent(in) :: self
        class(string),              intent(in) :: fname
        character(len=*),           intent(in) :: title
        integer,                    intent(in) :: version
        character(len=*), optional, intent(in) :: extra_keys(:), extra_vals(:)
        type(string) :: wfile, wname
        integer :: funit, io_stat, k, i
        if( .not. self%exists ) THROW_HARD('student not constructed; polar_student%write')
        wfile = weights_fname(fname)
        wname = basename(wfile)
        call self%net%write_weights(wfile)
        call fopen(funit, fname, status='replace', action='write', iostat=io_stat)
        call fileiochk('polar_student%write; open '//fname%to_char(), io_stat)
        write(funit,'(A)')        '# '//trim(title)
        write(funit,'(A,I0)')     'model_version=', version
        write(funit,'(A,A)')      'weights_file=', wname%to_char()
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
        if( allocated(self%teacher_project) ) write(funit,'(A,A)') 'teacher_project=', self%teacher_project
        write(funit,'(A,I0)')     'teacher_selection_hash=', self%selection_hash
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
            msg = 'the model header names no weights file; polar_student%load'
        else if( scan(trim(wname), '/\') > 0 )then
            msg = 'the model header must name its weights by a file name beside it; polar_student%load'
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
            msg = 'the weights file holds a network of another shape than its header''s; polar_student%load'
        else if( net%get_seed() /= seed )then
            msg = 'the weights file holds a network of another seed than its header''s; polar_student%load'
        else if( net%get_nweights() /= nweights )then
            msg = 'the weights file holds another weight count than its header''s; polar_student%load'
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

    !> A model header's values: the core's keys in a first pass (the shell weights' array needs the band), the shell
    !! weights in a second; the keys the core does not know are kept, in order, for the header's owner
    subroutine header_read( self, fname )
        class(polar_student_header), intent(inout) :: self
        class(string),               intent(in)    :: fname
        character(len=:), allocatable :: line, key, val
        integer :: funit, io_stat, ieq, k, pass
        logical :: l_wk_even, l_wk_odd
        call self%kill
        if( .not. file_exists(fname) ) THROW_HARD('no model file '//fname%to_char()//'; polar_student_header%read')
        self%weights_file = ''
        self%holdout_rule = ''
        self%teacher_project = ''
        allocate(self%keys(MAX_EXTRA), self%vals(MAX_EXTRA))
        l_wk_even = .false.
        l_wk_odd  = .false.
        do pass = 1,2
            call fopen(funit, fname, status='old', action='read', iostat=io_stat)
            call fileiochk('polar_student_header%read; open '//fname%to_char(), io_stat)
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
                        case('model_version');  read(val,*) self%version
                        case('weights_file');   self%weights_file = val
                        case('pftsz');          read(val,*) self%pftsz
                        case('nrots');          read(val,*) self%nrots
                        case('kfromto');        read(val,*) self%kfromto(1), self%kfromto(2)
                        case('lp');             read(val,*) self%lp
                        case('box');            read(val,*) self%box
                        case('box_crop');       read(val,*) self%box_crop
                        case('smpd');           read(val,*) self%smpd
                        case('mskdiam');        read(val,*) self%mskdiam
                        case('with_ctf');       read(val,*) self%with_ctf
                        case('width');          read(val,*) self%width
                        case('nweights');       read(val,*) self%nweights
                        case('seed');           read(val,*) self%seed
                        case('prob_beta');      read(val,*) self%prob_beta
                        case('shift_aug_px');   read(val,*) self%shift_aug_px
                        case('side_mean');      read(val,*) self%side_mean
                        case('side_sd');        read(val,*) self%side_sd
                        case('nptcls_train');   read(val,*) self%nptcls_train
                        case('nheldout');       read(val,*) self%nheldout
                        case('holdout_rule');   self%holdout_rule = val
                        case('batch_size');     read(val,*) self%batch_size
                        case('epochs_run');     read(val,*) self%epochs_run
                        case('best_epoch');     read(val,*) self%best_epoch
                        case('loss_train');     read(val,*) self%loss_train
                        case('loss_heldout');   read(val,*) self%loss_held
                        case('mode_heldout');   read(val,*) self%mode_held
                        case('time_train_s');   read(val,*) self%time_train_s
                        case('conf_threshold'); read(val,*) self%conf_threshold
                        case('conf_fraction');  read(val,*) self%conf_fraction
                        case('heldout_seeded'); read(val,*) self%heldout_seeded
                        case('heldout_recovered'); read(val,*) self%heldout_recovered
                        case('teacher_project'); self%teacher_project = val
                        case('teacher_selection_hash'); read(val,*) self%selection_hash
                        case('wk_even', 'wk_odd')
                        case default
                            if( self%nextra >= MAX_EXTRA ) THROW_HARD('too many keys in the model header; polar_student_header%read')
                            self%nextra = self%nextra + 1
                            self%keys(self%nextra) = key
                            self%vals(self%nextra) = val
                    end select
                else
                    select case(key)
                        case('wk_even')
                            read(val,*) (self%wk(k,1), k = self%kfromto(1),self%kfromto(2))
                            l_wk_even = .true.
                        case('wk_odd')
                            read(val,*) (self%wk(k,2), k = self%kfromto(1),self%kfromto(2))
                            l_wk_odd = .true.
                        case default
                    end select
                endif
            enddo
            call fclose(funit)
            if( pass == 1 )then
                if( self%pftsz < 2 .or. self%kfromto(1) < 1 .or. self%kfromto(2) < self%kfromto(1) .or. self%box < 2 ) &
                    &THROW_HARD('incomplete geometry in the model header; polar_student_header%read')
                allocate(self%wk(self%kfromto(1):self%kfromto(2),2), source=0.)
            endif
        enddo
        if( .not. (l_wk_even .and. l_wk_odd) ) THROW_HARD('the model header lacks its shell weights; polar_student_header%read')
    end subroutine header_read

    subroutine header_kill( self )
        class(polar_student_header), intent(inout) :: self
        if( allocated(self%wk)   ) deallocate(self%wk)
        if( allocated(self%keys) ) deallocate(self%keys, self%vals)
        if( allocated(self%weights_file)    ) deallocate(self%weights_file)
        if( allocated(self%holdout_rule)    ) deallocate(self%holdout_rule)
        if( allocated(self%teacher_project) ) deallocate(self%teacher_project)
        self%nextra  = 0
        self%version = -1
    end subroutine header_kill

    !> after new with a header's geometry: the weights the header fname names beside it, which must be the network the
    !! header describes, and the header's records
    subroutine load( self, hdr, fname )
        class(polar_student),       intent(inout) :: self
        type(polar_student_header), intent(in)    :: hdr
        class(string),              intent(in)    :: fname
        character(len=:), allocatable :: msg
        type(string) :: wfile
        if( .not. self%exists ) THROW_HARD('student not constructed; polar_student%load')
        call header_weights_path(fname, hdr%weights_file, wfile, msg)
        if( len(msg) > 0 ) THROW_HARD(msg)
        call self%net%read_weights(wfile)
        msg = weights_header_mismatch(self%net, self%pftsz_net, self%nk, self%nclass, hdr%width, hdr%seed, hdr%nweights)
        if( len(msg) > 0 ) THROW_HARD(msg)
        self%prob_beta         = hdr%prob_beta
        self%shift_aug_px      = hdr%shift_aug_px
        self%side_mean         = hdr%side_mean
        self%side_sd           = hdr%side_sd
        self%nptcls_train      = hdr%nptcls_train
        self%nheldout          = hdr%nheldout
        self%holdout_rule      = hdr%holdout_rule
        self%batch_size        = hdr%batch_size
        self%epochs_run        = hdr%epochs_run
        self%best_epoch        = hdr%best_epoch
        self%loss_train        = hdr%loss_train
        self%loss_held         = hdr%loss_held
        self%mode_held         = hdr%mode_held
        self%time_train_s      = hdr%time_train_s
        self%conf_threshold    = hdr%conf_threshold
        self%conf_fraction     = hdr%conf_fraction
        self%heldout_seeded    = hdr%heldout_seeded
        self%heldout_recovered = hdr%heldout_recovered
        self%teacher_project   = hdr%teacher_project
        self%selection_hash    = hdr%selection_hash
        call wfile%kill
    end subroutine load

    ! ===== GETTERS AND SETTERS =====

    pure integer function get_pftsz( self )
        class(polar_student), intent(in) :: self
        get_pftsz = self%pftsz
    end function get_pftsz

    pure integer function get_pftsz_net( self )
        class(polar_student), intent(in) :: self
        get_pftsz_net = self%pftsz_net
    end function get_pftsz_net

    pure function get_kfromto( self ) result( kfromto )
        class(polar_student), intent(in) :: self
        integer :: kfromto(2)
        kfromto = self%kfromto
    end function get_kfromto

    pure integer function get_nk( self )
        class(polar_student), intent(in) :: self
        get_nk = self%nk
    end function get_nk

    pure integer function get_nclass( self )
        class(polar_student), intent(in) :: self
        get_nclass = self%nclass
    end function get_nclass

    pure real(sp) function get_lp( self )
        class(polar_student), intent(in) :: self
        get_lp = self%lp
    end function get_lp

    pure integer function get_box( self )
        class(polar_student), intent(in) :: self
        get_box = self%box
    end function get_box

    pure integer function get_box_crop( self )
        class(polar_student), intent(in) :: self
        get_box_crop = self%box_crop
    end function get_box_crop

    pure real(sp) function get_smpd( self )
        class(polar_student), intent(in) :: self
        get_smpd = self%smpd
    end function get_smpd

    pure real(sp) function get_mskdiam( self )
        class(polar_student), intent(in) :: self
        get_mskdiam = self%mskdiam
    end function get_mskdiam

    pure logical function get_with_ctf( self )
        class(polar_student), intent(in) :: self
        get_with_ctf = self%with_ctf
    end function get_with_ctf

    subroutine get_weights( self, wk )
        class(polar_student),  intent(in)  :: self
        real(sp), allocatable, intent(out) :: wk(:,:)
        allocate(wk(self%kfromto(1):self%kfromto(2),2), source=self%wk)
    end subroutine get_weights

    pure integer(int64) function get_seed( self )
        class(polar_student), intent(in) :: self
        get_seed = self%seed
    end function get_seed

    !> the training record for the caller's reports
    subroutine get_record( self, nptcls_train, nheldout, epochs_run, best_epoch, loss_train, loss_held, mode_held, &
        &time_train_s, nweights )
        class(polar_student), intent(in)  :: self
        integer,              intent(out) :: nptcls_train, nheldout, epochs_run, best_epoch, nweights
        real(sp),             intent(out) :: loss_train, loss_held, mode_held, time_train_s
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
        class(polar_student), intent(in) :: self
        type(student_train_profile) :: prof
        prof%rt_assemble    = self%rt_assemble
        prof%rt_total       = real(self%time_train_s,dp)
        prof%epochs         = self%epochs_run
        prof%batch_size     = self%batch_size
        prof%nbatches_epoch = 0
        if( self%batch_size > 0 ) prof%nbatches_epoch = (self%nptcls_train + self%batch_size - 1) / self%batch_size
        prof%net            = self%net%get_profile()
    end function get_train_profile

    !> the teacher's provenance, recorded in the model file
    subroutine set_provenance( self, teacher_project, selection_hash, holdout_rule )
        class(polar_student), intent(inout) :: self
        character(len=*),     intent(in)    :: teacher_project, holdout_rule
        integer(int64),       intent(in)    :: selection_hash
        self%teacher_project = trim(teacher_project)
        self%selection_hash  = selection_hash
        self%holdout_rule    = trim(holdout_rule)
    end subroutine set_provenance

    !> the confidence threshold set on the held-out rows: above it the teacher's class lies among the candidates in
    !! conf_fraction of the held-out particles; heldout_seeded the held-out fraction above it within the owner's
    !! bounds, heldout_recovered the recovered fraction among those
    subroutine set_threshold( self, conf_threshold, conf_fraction, heldout_seeded, heldout_recovered )
        class(polar_student), intent(inout) :: self
        real(sp),             intent(in)    :: conf_threshold, conf_fraction, heldout_seeded, heldout_recovered
        self%conf_threshold    = conf_threshold
        self%conf_fraction     = conf_fraction
        self%heldout_seeded    = heldout_seeded
        self%heldout_recovered = heldout_recovered
    end subroutine set_threshold

    pure subroutine get_threshold( self, conf_threshold, conf_fraction, heldout_seeded, heldout_recovered )
        class(polar_student), intent(in)  :: self
        real(sp),             intent(out) :: conf_threshold, conf_fraction, heldout_seeded, heldout_recovered
        conf_threshold    = self%conf_threshold
        conf_fraction     = self%conf_fraction
        heldout_seeded    = self%heldout_seeded
        heldout_recovered = self%heldout_recovered
    end subroutine get_threshold

    subroutine kill( self )
        class(polar_student), intent(inout) :: self
        if( allocated(self%net) )then
            call self%net%kill
            deallocate(self%net)
        endif
        call self%feats%kill
        if( allocated(self%wk) ) deallocate(self%wk)
        if( allocated(self%holdout_rule)    ) deallocate(self%holdout_rule)
        if( allocated(self%teacher_project) ) deallocate(self%teacher_project)
        self%pftsz     = 0
        self%pftsz_net = 0
        self%nrots     = 0
        self%nk        = 0
        self%kfromto   = 0
        self%box       = 0
        self%box_crop  = 0
        self%width     = 0
        self%nclass    = 0
        self%nweights  = 0
        self%exists    = .false.
    end subroutine kill

end module simple_polar_student
