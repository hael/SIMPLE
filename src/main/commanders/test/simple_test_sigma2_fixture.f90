!@descr: sigma2 store fixtures shared by testers: a committed unit store and a header-word patch (simple_sigma2_state_file)
module simple_test_sigma2_fixture
use, intrinsic :: iso_fortran_env, only: int64, real32
use simple_syslib,            only: del_file
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_init_header, sigma2_state_create_candidate,&
                                   &sigma2_state_write_particles, sigma2_state_write_global, sigma2_state_publish,&
                                   &SIGMA2_PROV_PSPEC
use simple_test_utils,        only: assert_true, assert_int
implicit none
private
public :: write_test_sigma2_store, patch_sigma2_header_word

contains

    !> overwrite header word iword of the store fname (the header carries no checksum)
    subroutine patch_sigma2_header_word( fname, iword, value )
        character(len=*), intent(in) :: fname
        integer,        intent(in) :: iword
        integer(int64), intent(in) :: value
        integer(int64) :: words(32)
        integer :: funit, io_stat
        open(newunit=funit, file=fname, access='stream', form='unformatted', status='old', &
            &action='readwrite', iostat=io_stat)
        call assert_true(io_stat == 0, 'store can be opened for patching')
        read(funit, pos=17, iostat=io_stat) words
        words(iword) = value
        write(funit, pos=17, iostat=io_stat) words
        close(funit)
        call assert_true(io_stat == 0, 'store header can be patched')
    end subroutine patch_sigma2_header_word

    !> A committed unit store of nptcls records over the shells 1..box/2 (fdim(box)-1): the fixture
    !! testers read their residual records from (read_part), as a matcher reads them
    subroutine write_test_sigma2_store( fname, nptcls, box )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: nptcls, box
        type(sigma2_state_header) :: header
        real(real32), allocatable :: records(:,:)
        real(real32) :: globals(box/2,2)
        character(len=128) :: candidate, message
        integer :: status
        candidate = trim(fname)//'.next'
        call del_file(fname)
        call del_file(candidate)
        call sigma2_state_init_header(header, 1, box/2, nptcls, box, 1.0, 1_int64, 4242_int64, SIGMA2_PROV_PSPEC)
        call sigma2_state_create_candidate(candidate, header, status, message)
        call assert_int(0, status, 'sigma2 state operation succeeds: '//trim(message))
        allocate(records(box/2,nptcls), source=1.0_real32)
        call sigma2_state_write_particles(candidate, 1, records, status, message)
        call assert_int(0, status, 'sigma2 state operation succeeds: '//trim(message))
        globals = 1.0_real32
        call sigma2_state_write_global(candidate, globals, status, message)
        call assert_int(0, status, 'sigma2 state operation succeeds: '//trim(message))
        call sigma2_state_publish(candidate, fname, status, message)
        call assert_int(0, status, 'sigma2 state operation succeeds: '//trim(message))
    end subroutine write_test_sigma2_store

end module simple_test_sigma2_fixture
