!@descr: unit tests for the carried 2D class sums and the worker contributions (simple_cavg_sums)
! Covers the file contract (round trip, refusals, publication through a temporary name), the
! owner-side reduction and blend (independent of the number of parts, centering shift applied
! once) and the crop-box upsample of one carried set, on small deterministic arrays.
module simple_cavg_sums_tester
use simple_core_module_api
use simple_test_utils
use simple_cavg_sums, only: cavg_sums, cavg_contrib_fname, CAVG_SUMS_STATE, CAVG_SUMS_CONTRIB, &
                           &CAVG_SUMS_OK, CAVG_SUMS_MISSING, CAVG_SUMS_CORRUPT, CAVG_SUMS_MISMATCH
implicit none
private
public :: run_all_cavg_sums_tests

integer, parameter :: NCLS = 3
integer, parameter :: BOX  = 16
real,    parameter :: SMPD = 2.5
real,    parameter :: TOL  = 1.e-5

contains

    subroutine run_all_cavg_sums_tests()
        write(*,'(A)') '**** running all class-average carry-over tests ****'
        call test_cavg_sums_roundtrip()
        call test_cavg_sums_refusals()
        call test_owner_blend_split()
        call test_cavg_sums_pad()
    end subroutine run_all_cavg_sums_tests

    ! deterministic pseudo-random sums (a linear congruential sequence, no global RNG state)
    subroutine fill( sums, seed )
        type(cavg_sums), intent(inout) :: sums
        integer,         intent(in)    :: seed
        complex(sp) :: ce(fdim(BOX),BOX,NCLS), co(fdim(BOX),BOX,NCLS)
        real        :: te(fdim(BOX),BOX,NCLS), to(fdim(BOX),BOX,NCLS)
        integer(kind=8) :: state
        integer :: i, j, k
        state = int(seed, 8)
        do k = 1, NCLS
            do j = 1, BOX
                do i = 1, fdim(BOX)
                    ce(i,j,k) = cmplx(next(), next(), kind=sp)
                    co(i,j,k) = cmplx(next(), next(), kind=sp)
                    te(i,j,k) = 1. + next()
                    to(i,j,k) = 1. + next()
                enddo
            enddo
        enddo
        call sums%set_sums(ce, co, te, to)

    contains

        real function next()
            state = modulo(state * 1103515245_8 + 12345_8, 2147483648_8)
            next  = real(state) / 2147483648.
        end function next

    end subroutine fill

    ! largest absolute difference over the four accumulators
    real function maxdiff( a, b )
        type(cavg_sums), intent(in) :: a, b
        complex(sp) :: ae(fdim(BOX),BOX,NCLS), ao(fdim(BOX),BOX,NCLS), be(fdim(BOX),BOX,NCLS), bo(fdim(BOX),BOX,NCLS)
        real        :: ate(fdim(BOX),BOX,NCLS), ato(fdim(BOX),BOX,NCLS), bte(fdim(BOX),BOX,NCLS), bto(fdim(BOX),BOX,NCLS)
        call a%get_sums(ae, ao, ate, ato)
        call b%get_sums(be, bo, bte, bto)
        maxdiff = max(maxval(abs(ae - be)), maxval(abs(ao - bo)), maxval(abs(ate - bte)), maxval(abs(ato - bto)))
    end function maxdiff

    subroutine test_cavg_sums_roundtrip()
        type(cavg_sums)      :: a, b
        type(string)         :: fname
        real,    allocatable :: mrep(:), offsets(:,:)
        integer, allocatable :: pops(:,:)
        integer :: status
        write(*,'(A)') 'test_cavg_sums_roundtrip'
        fname = 'cavg_sums_test_state.bin'
        call a%new(CAVG_SUMS_STATE, NCLS, BOX, SMPD)
        call fill(a, 7)
        call a%set_mrep([10., 0., 3.5])
        call a%write(fname)
        call assert_false(file_exists(fname//'.tmp'), 'state round trip: no temporary file left')
        call b%read(fname, CAVG_SUMS_STATE, status)
        call assert_int(CAVG_SUMS_OK, status, 'state round trip: read status')
        call assert_true(b%matches(NCLS, BOX, SMPD), 'state round trip: geometry')
        call assert_real(0., maxdiff(a, b), 0., 'state round trip: sums bit identical')
        call b%get_mrep(mrep)
        call assert_true(all(mrep == [10., 0., 3.5]), 'state round trip: M(c)')
        call del_file(fname)
        ! contribution: offsets, populations and the carry-over flag
        fname = cavg_contrib_fname(3)
        call assert_true(fname%to_char() == 'cavg_contrib_part3.bin', 'contribution name has no numlen padding')
        call a%new(CAVG_SUMS_CONTRIB, NCLS, BOX, SMPD, part=3)
        call fill(a, 11)
        call a%set_contrib_meta(reshape([1.,-2., 0.,0., 0.5,0.25], [2,NCLS]), reshape([4,5, 0,0, 1,2], [2,NCLS]), .true.)
        call a%write(fname)
        call b%read(fname, CAVG_SUMS_CONTRIB, status)
        call assert_int(CAVG_SUMS_OK, status, 'contribution round trip: read status')
        call assert_int(3, b%get_part(), 'contribution round trip: part')
        call assert_true(b%get_l_frac(), 'contribution round trip: carry-over flag')
        call b%get_offsets(offsets)
        call b%get_eo_pops(pops)
        call assert_true(all(offsets == reshape([1.,-2., 0.,0., 0.5,0.25], [2,NCLS])), 'contribution round trip: offsets')
        call assert_true(all(pops == reshape([4,5, 0,0, 1,2], [2,NCLS])), 'contribution round trip: populations')
        call assert_real(0., maxdiff(a, b), 0., 'contribution round trip: sums bit identical')
        call del_file(fname)
        call a%kill
        call b%kill
    end subroutine test_cavg_sums_roundtrip

    ! Every unusable set is refused with a status that selects a full update
    subroutine test_cavg_sums_refusals()
        type(cavg_sums)  :: a, b
        type(string)     :: fname
        character(len=1), allocatable :: bytes(:)
        integer(kind=8)  :: fsize
        integer :: status, funit, ios
        write(*,'(A)') 'test_cavg_sums_refusals'
        fname = 'cavg_sums_test_refuse.bin'
        call del_file(fname)
        call b%read(fname, CAVG_SUMS_STATE, status)
        call assert_int(CAVG_SUMS_MISSING, status, 'refusal: missing set')
        ! a leftover temporary file is never read
        call a%new(CAVG_SUMS_STATE, NCLS, BOX, SMPD)
        call fill(a, 3)
        call a%write(fname)
        call simple_rename(fname, fname//'.tmp')
        call b%read(fname, CAVG_SUMS_STATE, status)
        call assert_int(CAVG_SUMS_MISSING, status, 'refusal: only a temporary file (write never completed)')
        ! a later publication replaces the leftover
        call a%write(fname)
        call assert_false(file_exists(fname//'.tmp'), 'refusal: publication removes a leftover temporary file')
        ! another kind
        call b%read(fname, CAVG_SUMS_CONTRIB, status)
        call assert_int(CAVG_SUMS_MISMATCH, status, 'refusal: a state is not a contribution')
        ! geometry disagreeing with the run
        call b%read(fname, CAVG_SUMS_STATE, status)
        call assert_false(b%matches(NCLS + 1, BOX, SMPD), 'refusal: class count')
        call assert_false(b%matches(NCLS, BOX + 2, SMPD), 'refusal: box_crop')
        call assert_false(b%matches(NCLS, BOX, 1.1 * SMPD), 'refusal: smpd_crop')
        ! truncated (a half-written set published by another tool)
        inquire(file=fname%to_char(), size=fsize)
        allocate(bytes(fsize - 100))
        open(newunit=funit, file=fname%to_char(), access='stream', status='old', action='read', iostat=ios)
        read(funit) bytes
        close(funit)
        open(newunit=funit, file=fname%to_char(), access='stream', status='replace', action='write', iostat=ios)
        write(funit) bytes
        close(funit)
        call b%read(fname, CAVG_SUMS_STATE, status)
        call assert_int(CAVG_SUMS_CORRUPT, status, 'refusal: truncated set')
        ! foreign content with the right length
        bytes = 'x'
        open(newunit=funit, file=fname%to_char(), access='stream', status='replace', action='write', iostat=ios)
        write(funit) bytes
        close(funit)
        call b%read(fname, CAVG_SUMS_STATE, status)
        call assert_int(CAVG_SUMS_CORRUPT, status, 'refusal: foreign file')
        call del_file(fname)
        call a%kill
        call b%kill
    end subroutine test_cavg_sums_refusals

    ! The owner's result depends on the current sums only through their total: splitting them
    ! into 1, 2 or 4 worker contributions, written, read back and reduced in ascending part order,
    ! gives the single-process result up to summation order. The centering shift is applied once
    ! to the previous set: a result with the shift applied twice differs.
    subroutine test_owner_blend_split()
        real, parameter :: S(NCLS) = [1., 1., 2.]
        real, parameter :: W(NCLS) = [0.9, 1.25, 0.]
        real, parameter :: OFFSETS(2,NCLS) = reshape([1.5,-0.5, 0.,0., -2.,3.], [2,NCLS])
        integer, parameter :: NSPLIT(3) = [1, 2, 4]
        type(cavg_sums) :: cur, prev, ref, total, part_sums, twice
        complex(sp) :: ce(fdim(BOX),BOX,NCLS), co(fdim(BOX),BOX,NCLS)
        real        :: te(fdim(BOX),BOX,NCLS), to(fdim(BOX),BOX,NCLS)
        real        :: frac
        integer     :: isplit, ipart, icls, status, nparts
        write(*,'(A)') 'test_owner_blend_split'
        call cur%new(CAVG_SUMS_STATE, NCLS, BOX, SMPD)
        call fill(cur, 21)
        call prev%new(CAVG_SUMS_STATE, NCLS, BOX, SMPD)
        call fill(prev, 42)
        call prev%write(string('cavg_sums_test_prev.bin'))
        ! single-process reference
        call ref%new(CAVG_SUMS_STATE, NCLS, BOX, SMPD)
        call ref%set_sums(zero_c(), zero_c(), zero_r(), zero_r())
        call ref%accumulate(cur)
        call shifted_prev(1, prev)
        call ref%blend(prev, S, W)
        call cur%get_sums(ce, co, te, to)
        do isplit = 1, size(NSPLIT)
            nparts = NSPLIT(isplit)
            do ipart = 1, nparts
                ! unequal shares that sum to one
                frac = real(ipart) / real(nparts * (nparts + 1) / 2)
                call part_sums%new(CAVG_SUMS_CONTRIB, NCLS, BOX, SMPD, part=ipart)
                call part_sums%set_sums(frac * ce, frac * co, frac * te, frac * to)
                call part_sums%write(cavg_contrib_fname(ipart))
            enddo
            call total%new(CAVG_SUMS_STATE, NCLS, BOX, SMPD)
            do ipart = 1, nparts
                call part_sums%read(cavg_contrib_fname(ipart), CAVG_SUMS_CONTRIB, status)
                call assert_int(CAVG_SUMS_OK, status, 'owner blend: contribution read')
                call total%accumulate(part_sums)
                call del_file(cavg_contrib_fname(ipart))
            enddo
            call shifted_prev(1, prev)
            call total%blend(prev, S, W)
            call assert_real(0., maxdiff(total, ref), TOL, 'owner blend: independent of the number of parts')
        enddo
        ! the shift is applied once
        call twice%new(CAVG_SUMS_STATE, NCLS, BOX, SMPD)
        call twice%accumulate(cur)
        call shifted_prev(2, prev)
        call twice%blend(prev, S, W)
        call assert_true(maxdiff(twice, ref) > 100. * TOL, 'owner blend: centering shift applied once')
        ! shifting there and back restores the previous set
        call shifted_prev(0, prev)
        do icls = 1, NCLS
            call prev%shift_class(icls,  OFFSETS(:,icls))
            call prev%shift_class(icls, -OFFSETS(:,icls))
        enddo
        call shifted_prev(0, twice)
        call assert_real(0., maxdiff(prev, twice), TOL, 'shift and inverse shift restore the sums')
        call del_file(string('cavg_sums_test_prev.bin'))
        call cur%kill; call prev%kill; call ref%kill; call total%kill; call part_sums%kill; call twice%kill

    contains

        ! the previous set as read by the owner, shifted ntimes by the centering offsets
        subroutine shifted_prev( ntimes, sums )
            integer,         intent(in)    :: ntimes
            type(cavg_sums), intent(inout) :: sums
            integer :: i, c, st
            call sums%read(string('cavg_sums_test_prev.bin'), CAVG_SUMS_STATE, st)
            do i = 1, ntimes
                do c = 1, NCLS
                    call sums%shift_class(c, OFFSETS(:,c))
                enddo
            enddo
        end subroutine shifted_prev

        function zero_c() result( z )
            complex(sp) :: z(fdim(BOX),BOX,NCLS)
            z = cmplx(0.,0.,kind=sp)
        end function zero_c

        function zero_r() result( z )
            real :: z(fdim(BOX),BOX,NCLS)
            z = 0.
        end function zero_r

    end subroutine test_owner_blend_split

    ! The pool's crop-box upsample pads the one carried set once; padding is linear, so it
    ! agrees with the legacy per-part padding summed over parts. Padding keeps every logical
    ! Fourier component, zero-fills the new ones, and keeps M(c).
    subroutine test_cavg_sums_pad()
        integer, parameter :: BOX_PD = 24
        type(cavg_sums) :: a, b, total, padsum
        complex(sp) :: ce(fdim(BOX),BOX,NCLS), co(fdim(BOX),BOX,NCLS)
        real        :: te(fdim(BOX),BOX,NCLS), to(fdim(BOX),BOX,NCLS)
        complex(sp) :: pe(fdim(BOX_PD),BOX_PD,NCLS), po(fdim(BOX_PD),BOX_PD,NCLS)
        real        :: pte(fdim(BOX_PD),BOX_PD,NCLS), pto(fdim(BOX_PD),BOX_PD,NCLS)
        complex(sp) :: qe(fdim(BOX_PD),BOX_PD,NCLS), qo(fdim(BOX_PD),BOX_PD,NCLS)
        real        :: qte(fdim(BOX_PD),BOX_PD,NCLS), qto(fdim(BOX_PD),BOX_PD,NCLS)
        real, allocatable :: mrep(:)
        real    :: smpd_pd
        integer :: nonzero
        write(*,'(A)') 'test_cavg_sums_pad'
        smpd_pd = SMPD * real(BOX) / real(BOX_PD)
        call a%new(CAVG_SUMS_STATE, NCLS, BOX, SMPD)
        call fill(a, 5)
        call b%new(CAVG_SUMS_STATE, NCLS, BOX, SMPD)
        call fill(b, 6)
        call total%new(CAVG_SUMS_STATE, NCLS, BOX, SMPD)
        call total%accumulate(a)
        call total%accumulate(b)
        call total%set_mrep([1., 2., 3.])
        ! one set padded once
        call total%pad_to(BOX_PD, smpd_pd)
        call assert_true(total%matches(NCLS, BOX_PD, smpd_pd), 'pad: new geometry')
        call total%get_mrep(mrep)
        call assert_true(all(mrep == [1., 2., 3.]), 'pad: M(c) unchanged')
        ! per-part padding summed over parts
        call a%pad_to(BOX_PD, smpd_pd)
        call b%pad_to(BOX_PD, smpd_pd)
        call padsum%new(CAVG_SUMS_STATE, NCLS, BOX_PD, smpd_pd)
        call padsum%accumulate(a)
        call padsum%accumulate(b)
        call total%get_sums(pe, po, pte, pto)
        call padsum%get_sums(qe, qo, qte, qto)
        call assert_real(0., max(maxval(abs(pe - qe)), maxval(abs(pte - qte)), maxval(abs(po - qo))), TOL, &
            &'pad: one padded set equals the per-part padding summed')
        ! components are moved, not changed: the padded CTF^2 total equals the original total
        call a%new(CAVG_SUMS_STATE, NCLS, BOX, SMPD)
        call fill(a, 5)
        call a%get_sums(ce, co, te, to)
        call a%pad_to(BOX_PD, smpd_pd)
        call a%get_sums(pe, po, pte, pto)
        call assert_real(sum(te), sum(pte), 1.e-3, 'pad: CTF^2 mass preserved')
        nonzero = count(pte /= 0.)
        call assert_int(count(te /= 0.), nonzero, 'pad: new components are zero')
        call a%kill; call b%kill; call total%kill; call padsum%kill
    end subroutine test_cavg_sums_pad

end module simple_cavg_sums_tester
