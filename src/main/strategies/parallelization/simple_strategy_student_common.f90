!@descr: the students' shared strategy helpers: the held-out cohort, the resampled design comparison, the confidence threshold, a model's geometry against a project's, a part's range of the selected particles, the prediction's distributed master and its parts' sigma2, a file's hash
! Stateless routines the learn and predict programs of the pose and the class students share (doc/algorithms/pose_student.md).
module simple_strategy_student_common
use, intrinsic :: iso_fortran_env, only: int32, int64
use simple_defs,        only: sp, dp, logfhandle
use simple_timer,       only: tic, toc, timer_int_kind
use simple_error,       only: simple_exception
use simple_string,      only: string
use simple_fileio,      only: fopen, fclose, fileiochk
use simple_syslib,      only: file_exists
use simple_math,        only: hpsort
use simple_polar_feats, only: polar_feats_hash_ints
use simple_defs_hhcnn,  only: STUDENT_HOLDOUT_EVERY, STUDENT_CONF_NONE, STUDENT_SMPD_TOL, STUDENT_MSKDIAM_TOL, &
    &STUDENT_NBOOT, STUDENT_DESIGN_QUANTILE
implicit none

public :: student_holdout_split, student_boot_quantile, student_conf_threshold, student_geometry_mismatch
public :: student_part_pinds, student_predict_parts, student_sigma2_for_pass, concat_part_tables
public :: file_hash, run_all_student_strategy_tests
public :: student_resampled_gain, choose_student_design, student_design_trial, student_design_keep, fname_base
private
#include "simple_local_flags.inc"


interface
    module subroutine run_all_student_strategy_tests()
    end subroutine run_all_student_strategy_tests
end interface

abstract interface
    !> the owner's trial of design iw: built, trained, and judged on the held-out rows, recov the rows it recovers
    subroutine student_design_trial( iw, recov )
        integer, intent(in)    :: iw
        logical, intent(inout) :: recov(:)
    end subroutine student_design_trial
    !> the owner's record of design iw as the kept model
    subroutine student_design_keep( iw )
        integer, intent(in) :: iw
    end subroutine student_design_keep
end interface

contains

    !> The held-out cohort of n selected particles on stacks stkinds: with STUDENT_HOLDOUT_EVERY or more distinct stacks every
    !! STUDENT_HOLDOUT_EVERY-th of them in order of stkind (group_of the stack), otherwise every STUDENT_HOLDOUT_EVERY-th particle
    !! (group_of the particle); msg says why a cohort would be empty (too few particles), empty otherwise
    subroutine student_holdout_split( n, stkinds, held, group_of, l_by_stack, msg )
        integer,                       intent(in)  :: n, stkinds(n)
        logical,                       intent(out) :: held(n)
        integer,                       intent(out) :: group_of(n)
        logical,                       intent(out) :: l_by_stack
        character(len=:), allocatable, intent(out) :: msg
        integer, allocatable :: rank(:)
        integer :: i, ng, smax
        held       = .false.
        group_of   = 0
        l_by_stack = .false.
        msg        = ''
        if( n < 1 )then
            msg = 'no selected particle to split; student_holdout_split'
            return
        endif
        if( any(stkinds < 1) )then
            msg = 'a selected particle has no stack; student_holdout_split'
            return
        endif
        ! the selected stacks ranked in order of stkind
        smax = maxval(stkinds)
        allocate(rank(smax), source=0)
        do i = 1,n
            rank(stkinds(i)) = 1
        enddo
        ng = 0
        do i = 1,smax
            if( rank(i) > 0 )then
                ng      = ng + 1
                rank(i) = ng
            endif
        enddo
        l_by_stack = ng >= STUDENT_HOLDOUT_EVERY
        do i = 1,n
            if( l_by_stack )then
                group_of(i) = stkinds(i)
                held(i)     = mod(rank(stkinds(i)), STUDENT_HOLDOUT_EVERY) == 0
            else
                group_of(i) = i
                held(i)     = mod(i, STUDENT_HOLDOUT_EVERY) == 0
            endif
        enddo
        deallocate(rank)
        if( count(held) < 1 .or. count(.not. held) < 1 ) msg = 'too few selected particles to hold out every tenth '//&
            &'(at least ten are needed); student_holdout_split'
    end subroutine student_holdout_split

    !> the STUDENT_DESIGN_QUANTILE of STUDENT_NBOOT means of the ngroups gains resampled with replacement by a Lehmer generator seeded
    !! from seed; zero groups give zero, no evidence of a gain
    function student_boot_quantile( gain, ngroups, seed ) result( q )
        real,           intent(in) :: gain(:)
        integer,        intent(in) :: ngroups
        integer(int64), intent(in) :: seed
        real     :: q, boot_means(STUDENT_NBOOT)
        real(dp) :: gsum
        integer(int64) :: lehmer
        integer  :: ib, gg, iq
        q = 0.
        if( ngroups < 1 ) return
        lehmer = modulo(seed, 2147483646_int64) + 1_int64
        do ib = 1,STUDENT_NBOOT
            gsum = 0.d0
            do gg = 1,ngroups
                lehmer = modulo(48271_int64 * lehmer, 2147483647_int64)
                gsum   = gsum + real(gain(1 + int(modulo(lehmer, int(ngroups,int64)))),dp)
            enddo
            boot_means(ib) = real(gsum / real(ngroups,dp))
        enddo
        call hpsort(boot_means)
        iq = max(1, min(STUDENT_NBOOT, nint(STUDENT_DESIGN_QUANTILE * real(STUDENT_NBOOT))))
        q  = boot_means(iq)
    end function student_boot_quantile

    !> The paired difference of held-out recovery between a candidate design and the kept one: per micrograph
    !! (group_held, numbered here in order of appearance) the mean difference of the rows' recovery, resampled by
    !! micrograph from seed; the STUDENT_DESIGN_QUANTILE of the resampled means (student_boot_quantile)
    real function student_resampled_gain( nheld, group_held, recov_new, recov_old, seed ) result( q )
        integer,        intent(in) :: nheld, group_held(nheld)
        logical,        intent(in) :: recov_new(nheld), recov_old(nheld)
        integer(int64), intent(in) :: seed
        integer, allocatable :: groups(:), gfirst(:), gsize(:)
        real,    allocatable :: gain(:)
        integer :: jj, gg, ig, ngroups
        allocate(groups(max(1,nheld)), gfirst(max(1,nheld)), source=0)
        ngroups = 0
        do jj = 1,nheld
            gg = 0
            do ig = 1,ngroups
                if( group_held(gfirst(ig)) == group_held(jj) )then
                    gg = ig
                    exit
                endif
            enddo
            if( gg == 0 )then
                ngroups = ngroups + 1
                groups(jj) = ngroups
                gfirst(ngroups) = jj
            else
                groups(jj) = gg
            endif
        enddo
        allocate(gain(max(1,ngroups)), source=0.)
        allocate(gsize(max(1,ngroups)), source=0)
        do jj = 1,nheld
            gg = groups(jj)
            gsize(gg) = gsize(gg) + 1
            gain(gg)  = gain(gg) + merge(1., 0., recov_new(jj)) - merge(1., 0., recov_old(jj))
        enddo
        do gg = 1,ngroups
            gain(gg) = gain(gg) / real(max(1,gsize(gg)))
        enddo
        q = student_boot_quantile(gain, ngroups, seed)
        deallocate(groups, gfirst, gsize, gain)
    end function student_resampled_gain

    !> The design choice, cheapest first: design 1 is kept; each wider one replaces it when the STUDENT_DESIGN_QUANTILE of its
    !! resampled recovery gain is above zero, and the first that does not ends the trials (at most ntry); recov the rows
    !! each design recovers, iw_best the kept one, gain_q and rt_boot each comparison's quantile and cost
    subroutine choose_student_design( label, widths, ntry, nheld, group_held, seed, try_design, keep_design, recov, &
        &iw_best, ntried, gain_q, rt_boot )
        character(len=*),               intent(in)    :: label
        integer,                        intent(in)    :: widths(:), ntry, nheld, group_held(nheld)
        integer(int64),                 intent(in)    :: seed
        procedure(student_design_trial)               :: try_design
        procedure(student_design_keep)                :: keep_design
        logical,                        intent(inout) :: recov(:,:)
        integer,                        intent(out)   :: iw_best, ntried
        real,                           intent(inout) :: gain_q(:)
        real(timer_int_kind),           intent(inout) :: rt_boot(:)
        integer(timer_int_kind) :: tb
        real    :: recov_new, recov_best
        integer :: iw
        iw_best = 0
        ntried  = 0
        do iw = 1,ntry
            ntried = iw
            call try_design(iw, recov(:,iw))
            if( iw == 1 )then
                iw_best = 1
                call keep_design(1)
                cycle
            endif
            tb = tic()
            gain_q(iw)  = student_resampled_gain(nheld, group_held, recov(1:nheld,iw), recov(1:nheld,iw_best), seed)
            rt_boot(iw) = toc(tb)
            recov_new   = real(count(recov(1:nheld,iw)))      / real(nheld)
            recov_best  = real(count(recov(1:nheld,iw_best))) / real(nheld)
            write(logfhandle,'(A,I0,A,I0,A,F8.5,A,F5.2,A,F8.5)') '>>> '//trim(label)//' CHOICE: width ', widths(iw), &
                &' against ', widths(iw_best), ': recovery gain ', recov_new - recov_best, &
                &', its resampled quantile ', STUDENT_DESIGN_QUANTILE, ' = ', gain_q(iw)
            if( gain_q(iw) > 0. )then
                iw_best = iw
                call keep_design(iw)
            else
                exit
            endif
        enddo
    end subroutine choose_student_design

    !> The confidence threshold of n held-out rows: over the eligible rows in decreasing confidence, the longest prefix
    !! that ends with a group of equal confidences and whose recovered share is at least frac; threshold its last
    !! confidence (STUDENT_CONF_NONE when none qualifies), nsel its rows, nrec the recovered among them
    subroutine student_conf_threshold( n, conf, eligible, recovered, frac, threshold, nsel, nrec )
        integer, intent(in)  :: n
        real,    intent(in)  :: conf(n), frac
        logical, intent(in)  :: eligible(n), recovered(n)
        real,    intent(out) :: threshold
        integer, intent(out) :: nsel, nrec
        real,    allocatable :: negconf(:)
        integer, allocatable :: order(:)
        integer :: ne, i, k, nq, nr
        threshold = STUDENT_CONF_NONE
        nsel      = 0
        nrec      = 0
        ne        = count(eligible)
        if( ne < 1 ) return
        allocate(negconf(ne), order(ne))
        k = 0
        do i = 1,n
            if( .not. eligible(i) ) cycle
            k = k + 1
            order(k)   = i
            negconf(k) = -conf(i)
        enddo
        call hpsort(negconf, order)
        nq = 0
        nr = 0
        do k = 1,ne
            nq = nq + 1
            if( recovered(order(k)) ) nr = nr + 1
            ! a prefix that ends inside a group of equal confidences is no threshold: conf >= it takes the whole group
            if( k < ne )then
                if( conf(order(k+1)) >= conf(order(k)) ) cycle
            endif
            if( real(nr) >= frac * real(nq) )then
                threshold = conf(order(k))
                nsel      = nq
                nrec      = nr
            endif
        enddo
        deallocate(negconf, order)
    end subroutine student_conf_threshold

    !> empty when a model's geometry is the project's: the half ring, the box, the crop and the CTF treatment exactly,
    !! the sampling within STUDENT_SMPD_TOL and the mask within STUDENT_MSKDIAM_TOL (the header's decimals); else what differs
    function student_geometry_mismatch( pftsz_m, box_m, box_crop_m, smpd_m, mskdiam_m, ctf_m, pftsz, box, box_crop, &
        &smpd, mskdiam, ctf ) result( msg )
        integer, intent(in) :: pftsz_m, box_m, box_crop_m, pftsz, box, box_crop
        real,    intent(in) :: smpd_m, mskdiam_m, smpd, mskdiam
        logical, intent(in) :: ctf_m, ctf
        character(len=:), allocatable :: msg
        msg = ''
        if( pftsz_m /= pftsz )then
            msg = 'the model''s half ring differs from this project''s (mskdiam, box); student_geometry_mismatch'
        else if( box_m /= box )then
            msg = 'the model''s box differs from this project''s; student_geometry_mismatch'
        else if( box_crop_m /= box_crop )then
            msg = 'the model''s crop (box_crop) differs from this project''s; student_geometry_mismatch'
        else if( abs(smpd_m - smpd) > STUDENT_SMPD_TOL )then
            msg = 'the model''s sampling (smpd) differs from this project''s; student_geometry_mismatch'
        else if( abs(mskdiam_m - mskdiam) > STUDENT_MSKDIAM_TOL )then
            msg = 'the model''s mask (mskdiam) differs from this project''s; student_geometry_mismatch'
        else if( ctf_m .neqv. ctf )then
            msg = 'the model''s CTF treatment differs from this project''s; student_geometry_mismatch'
        endif
    end function student_geometry_mismatch

    !> the selected particles pinds_sel of part part of nparts within its range fromp..top of a field of nptcls
    !! particles (possibly none); msg says why the part or its range is refused, empty otherwise
    subroutine student_part_pinds( pinds_sel, part, nparts, fromp, top, nptcls, pinds, msg )
        integer,                       intent(in)  :: pinds_sel(:), part, nparts, fromp, top, nptcls
        integer, allocatable,          intent(out) :: pinds(:)
        character(len=:), allocatable, intent(out) :: msg
        msg = ''
        if( nparts < 1 .or. part < 1 .or. part > nparts )then
            msg = 'a part is one of nparts: 1 <= part <= nparts; student_part_pinds'
            allocate(pinds(0))
        else if( fromp < 1 .or. top < fromp .or. top > nptcls )then
            msg = 'a part''s range is 1 <= fromp <= top <= the particles; student_part_pinds'
            allocate(pinds(0))
        else
            pinds = pack(pinds_sel, pinds_sel >= fromp .and. pinds_sel <= top)
        endif
    end subroutine student_part_pinds

    !> The distributed master of a student's prediction: the parts scheduled over the project's rows (balanced over
    !! the selected ones) with the command line as their job, their records merged into the loaded field and their
    !! tables concatenated in part order into <base>_pred.txt. The caller writes the field; part_bases are the parts'
    !! file bases (<base>_part<i>), for the caller's own part files
    subroutine student_predict_parts( params, build, cline, base, part_bases )
        use simple_parameters, only: parameters
        use simple_builder,    only: builder
        use simple_cmdline,    only: cmdline
        use simple_chash,      only: chash
        use simple_qsys_env,   only: qsys_env
        use simple_qsys_funs,  only: qsys_cleanup
        use simple_defs,       only: L_USE_SLURM_ARR
        use simple_defs_fname, only: ALGN_FBODY
        use simple_syslib,     only: del_file
        use simple_string_utils, only: int2str
        class(parameters),         intent(inout) :: params
        class(builder),            intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
        class(string),             intent(in)    :: base
        type(string), allocatable, intent(out)   :: part_bases(:)
        type(qsys_env) :: qenv
        type(chash)    :: job_descr
        type(string), allocatable :: tables(:)
        integer :: ipart
        write(logfhandle,'(A,I0,A)') '>>> STUDENT PREDICTION: ', params%nparts, ' parts'
        call qenv%new(params, params%nparts, l_active=build%spproj_field%included())
        call cline%gen_job_descr(job_descr)
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR, &
            &extra_params=params)
        call build%spproj%merge_algndocs(params%nptcls, params%nparts, params%oritype, ALGN_FBODY, write=.false.)
        allocate(part_bases(params%nparts), tables(params%nparts))
        do ipart = 1,params%nparts
            part_bases(ipart) = base%to_char()//'_part'//int2str(ipart)
            tables(ipart)     = part_bases(ipart)%to_char()//'_pred.txt'
        enddo
        call concat_part_tables(string(base%to_char()//'_pred.txt'), tables)
        do ipart = 1,params%nparts
            call del_file(tables(ipart))
        enddo
        call qsys_cleanup(params)
        call qenv%kill
        call job_descr%kill
    end subroutine student_predict_parts

    !> The canonical sigma2 of a pass that registers shifts: made when missing (ensure_sigma2_for_iteration) by a
    !! shared-memory run or a master; a distributed part only checks it, since its master made it
    subroutine student_sigma2_for_pass( cline, projfile, oritype, iter, box, smpd, label, l_rebuilt )
        use simple_cmdline,          only: cmdline
        use simple_defs,             only: l_distr_worker_glob
        use simple_sigma2_bootstrap, only: ensure_sigma2_for_iteration, sigma2_estimate_available
        class(cmdline),   intent(in)  :: cline
        class(string),    intent(in)  :: projfile
        character(len=*), intent(in)  :: oritype, label
        integer,          intent(in)  :: iter, box
        real,             intent(in)  :: smpd
        logical,          intent(out) :: l_rebuilt
        l_rebuilt = .false.
        if( l_distr_worker_glob )then
            if( .not. sigma2_estimate_available(projfile, oritype, box, smpd) ) &
                &THROW_HARD('no consumable sigma2 for a part: its master makes it; student_sigma2_for_pass')
        else
            call ensure_sigma2_for_iteration(cline, projfile, oritype, iter, box, smpd, label, l_rebuilt)
        endif
    end subroutine student_sigma2_for_pass

    !> the text tables part_fnames joined into fname_out byte for byte: the first whole, the others without their
    !! header line
    subroutine concat_part_tables( fname_out, part_fnames )
        class(string), intent(in) :: fname_out, part_fnames(:)
        character(len=:), allocatable :: buf
        integer(int64) :: nbytes
        integer :: ipart, uin, uout, ios, i0
        open(newunit=uout, file=fname_out%to_char(), access='stream', form='unformatted', status='replace', &
            &action='write', iostat=ios)
        if( ios /= 0 ) THROW_HARD('could not open '//fname_out%to_char()//'; concat_part_tables')
        do ipart = 1,size(part_fnames)
            open(newunit=uin, file=part_fnames(ipart)%to_char(), access='stream', form='unformatted', status='old', &
                &action='read', iostat=ios)
            if( ios /= 0 ) THROW_HARD('a part''s table is missing: '//part_fnames(ipart)%to_char())
            inquire(unit=uin, size=nbytes)
            if( allocated(buf) ) deallocate(buf)
            allocate(character(len=nbytes) :: buf)
            if( nbytes > 0 ) read(uin) buf
            close(uin)
            i0 = 1
            if( ipart > 1 ) i0 = index(buf, new_line('a')) + 1
            if( i0 > 1 .or. ipart == 1 )then
                if( i0 <= len(buf) ) write(uout) buf(i0:)
            endif
        enddo
        close(uout)
    end subroutine concat_part_tables

    ! ===== HELPERS =====

    !> the hash of a file's bytes, read as 32-bit words (the references' provenance)
    function file_hash( fname ) result( h )
        class(string), intent(in) :: fname
        integer(int32), allocatable :: words(:)
        integer(int64) :: h, nbytes
        integer        :: funit, io_stat, nwords
        if( .not. file_exists(fname) ) THROW_HARD('no such file: '//fname%to_char()//'; file_hash')
        call fopen(funit, fname, 'old', 'read', iostat=io_stat, access='stream', form='unformatted')
        call fileiochk('file_hash; open '//fname%to_char(), io_stat)
        inquire(unit=funit, size=nbytes)
        nwords = int(nbytes / 4_int64)
        allocate(words(max(1,nwords)), source=0_int32)
        if( nwords > 0 ) read(funit) words(1:nwords)
        call fclose(funit)
        h = polar_feats_hash_ints(words)
        deallocate(words)
    end function file_hash

    !> fname without a trailing .txt
    function fname_base( fname ) result( base )
        class(string), intent(in) :: fname
        type(string) :: base
        character(len=:), allocatable :: s
        integer :: n
        s = fname%to_char()
        n = len_trim(s)
        if( n > 4 )then
            if( s(n-3:n) == '.txt' ) s = s(1:n-4)
        endif
        base = trim(s)
    end function fname_base

end module simple_strategy_student_common
