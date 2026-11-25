!@descr: phase-randomized FSC correction, input isolation and fixed-seed thread repeatability tests
module simple_fsc_tester
use iso_fortran_env,      only: int64
use ieee_arithmetic,      only: ieee_is_finite
use simple_string,        only: string
use simple_image,         only: image
use simple_fsc,           only: phase_rand_fsc
use simple_defs,          only: dp, PI
use simple_estimate_ssnr, only: get_resolution
use simple_test_utils,    only: assert_true, assert_real, assert_int, enter_fixture, leave_fixture, &
    &tests_failed, set_fixed_seed, offset_seed, fix_test_seed
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads, omp_get_num_threads, &
!$     omp_get_thread_num, omp_get_dynamic, omp_set_dynamic
implicit none
private
public :: run_all_fsc_tests
integer, parameter :: BOX=48, NSH=BOX/2, FIXED_SEED=20261124, TEST_THREADS=3
real,    parameter :: TOL=2.e-6

type :: fsc_curves
    real, allocatable :: corrected(:), masked(:), randomized(:), unmasked(:)
end type fsc_curves

contains

    subroutine run_all_fsc_tests()
        class(image), allocatable :: even, odd, mask
        type(fsc_curves)          :: serial_reference
        type(string)              :: cwd_saved, fixture_root
        integer                   :: previous_threads, nfail_before
        logical                   :: previous_dynamic
        real, allocatable         :: before_even(:,:,:), before_odd(:,:,:), before_mask(:,:,:)
        real, pointer             :: pixels(:,:,:)
        write(*,'(A)') '**** running all FSC tests ****'
        write(*,'(A,1X,I0)') 'FSC_TEST_SEED',offset_seed(FIXED_SEED)
        previous_threads=1
        previous_dynamic=.false.
        !$ previous_threads=omp_get_max_threads()
        !$ previous_dynamic=omp_get_dynamic()
        !$ call omp_set_dynamic(.false.)
        nfail_before=tests_failed
        call enter_fixture('tmp_fsc_phase',cwd_saved,fixture_root)
        allocate(even,odd,mask)
        call make_fixture(even,odd,mask)
        call even%get_rmat_ptr(pixels)
        before_even=pixels
        call odd%get_rmat_ptr(pixels)
        before_odd=pixels
        call mask%get_rmat_ptr(pixels)
        before_mask=pixels

        call test_serial_fsc(even,odd,mask,serial_reference)
        call test_seed_sensitivity(even,odd,mask,serial_reference)
        call test_threaded_fsc(even,odd,mask,serial_reference)
        call test_fresh_realizations(even,odd,mask)
        call test_input_isolation(even,odd,mask,before_even,before_odd,before_mask)
        call test_no_crossing_fallback(even,mask)
        call even%kill
        call odd%kill
        call mask%kill
        deallocate(even,odd,mask)
        !$ call omp_set_num_threads(previous_threads)
        !$ call omp_set_dynamic(previous_dynamic)
        call fix_test_seed()
        call leave_fixture(cwd_saved,fixture_root,nfail_before)
    end subroutine run_all_fsc_tests

    subroutine test_serial_fsc(even,odd,mask,a)
        class(image), intent(in)      :: even,odd
        class(image), intent(inout)   :: mask
        type(fsc_curves), intent(out) :: a
        type(fsc_curves)              :: b
        integer                       :: k,onset
        real                          :: expected,denominator
        write(*,'(A)') 'test_serial_fsc'
        !$ call omp_set_num_threads(1)
        call seeded_curves(even,odd,mask,FIXED_SEED,.false.,a)
        call seeded_curves(even,odd,mask,FIXED_SEED,.false.,b)
        call report_repeat('serial',1,a,b)
        call assert_repeat(a,b,'serial fixed seed')
        call assert_reference_curves(even,odd,mask,a)
        call assert_true(all(ieee_is_finite(a%corrected)), 'corrected FSC is finite')
        onset=0
        do k=2,NSH
            if(a%unmasked(k)<0.8)then
                onset=k
                exit
            endif
        enddo
        call assert_true(onset>=2.and.onset<=NSH-2,'fixture exercises the phase-randomization branch')
        if(onset>=2.and.onset<=NSH-2)then
            do k=1,NSH
                expected=a%masked(k)
                if(k>=onset+2)then
                    denominator=1.-a%randomized(k)
                    expected=0.
                    if(a%masked(k)>=a%randomized(k).and.denominator>tiny(1.)) &
                        expected=max(0.,min(1.,(a%masked(k)-a%randomized(k))/denominator))
                endif
                call assert_real(expected,a%corrected(k),TOL,'FSC follows solvent-correction equation')
            enddo
        endif
    end subroutine test_serial_fsc

    subroutine test_seed_sensitivity(even,odd,mask,a)
        class(image), intent(in)     :: even,odd
        class(image), intent(inout)  :: mask
        type(fsc_curves), intent(in) :: a
        type(fsc_curves)             :: b
        write(*,'(A)') 'test_seed_sensitivity'
        call seeded_curves(even,odd,mask,FIXED_SEED+1,.false.,b)
        call report_repeat('changed_seed',1,a,b)
        call assert_real(0.,maxval(abs(a%unmasked-b%unmasked)),TOL,'seed leaves unmasked FSC unchanged')
        call assert_real(0.,maxval(abs(a%masked-b%masked)),TOL,'seed leaves masked FSC unchanged')
        call assert_true(maxval(abs(a%randomized-b%randomized))>1.e-5,'changed seed changes random-phase estimate')
    end subroutine test_seed_sensitivity

    subroutine test_threaded_fsc(even,odd,mask,serial_reference)
        class(image), intent(in)     :: even,odd
        class(image), intent(inout)  :: mask
        type(fsc_curves), intent(in) :: serial_reference
        type(fsc_curves)             :: a,b
        integer                      :: actual
        write(*,'(A)') 'test_threaded_fsc'
        !$ call omp_set_num_threads(TEST_THREADS)
        actual=1
        !$omp parallel default(shared) num_threads(TEST_THREADS)
        !$omp single
        !$ actual=omp_get_num_threads()
        !$omp end single
        !$omp end parallel
        !$ call assert_int(TEST_THREADS,actual,'FSC repeatability exercises the requested OpenMP team')
        call seeded_curves(even,odd,mask,FIXED_SEED,.false.,a)
        call seeded_curves(even,odd,mask,FIXED_SEED,.false.,b)
        call report_repeat('ordinary_reseed',actual,a,b)
        call assert_repeat(a,b,'parallel fixed seed')
        call assert_repeat(serial_reference,a,'keyed phase thread-count invariance')
        call seeded_curves(even,odd,mask,FIXED_SEED,.true.,a)
        call seeded_curves(even,odd,mask,FIXED_SEED,.true.,b)
        call report_repeat('all_streams_reset',actual,a,b)
        call assert_repeat(a,b,'explicit worker-stream control')
    end subroutine test_threaded_fsc

    subroutine test_fresh_realizations(even,odd,mask)
        class(image), intent(in)    :: even,odd
        class(image), intent(inout) :: mask
        type(fsc_curves)            :: a,b
        write(*,'(A)') 'test_fresh_realizations'
        call set_fixed_seed(FIXED_SEED,propagate=.true.)
        call phase_rand_fsc(even,odd,mask,1,NSH,a%corrected,a%masked,a%randomized,a%unmasked)
        call phase_rand_fsc(even,odd,mask,1,NSH,b%corrected,b%masked,b%randomized,b%unmasked)
        call assert_true(maxval(abs(a%randomized-b%randomized))>1.e-5,'unreset consecutive calls keep fresh realizations')
    end subroutine test_fresh_realizations

    subroutine test_input_isolation(even,odd,mask,before_even,before_odd,before_mask)
        class(image), intent(in) :: even,odd,mask
        real, intent(in)         :: before_even(:,:,:),before_odd(:,:,:),before_mask(:,:,:)
        real, pointer            :: pixels(:,:,:)
        write(*,'(A)') 'test_input_isolation'
        call even%get_rmat_ptr(pixels)
        call assert_real(0.,maxval(abs(pixels-before_even)),0.,'even input remains unchanged')
        call odd%get_rmat_ptr(pixels)
        call assert_real(0.,maxval(abs(pixels-before_odd)),0.,'odd input remains unchanged')
        call mask%get_rmat_ptr(pixels)
        call assert_real(0.,maxval(abs(pixels-before_mask)),0.,'mask pixels remain unchanged')
    end subroutine test_input_isolation

    subroutine test_no_crossing_fallback(even,mask)
        class(image), intent(in)    :: even
        class(image), intent(inout) :: mask
        type(fsc_curves)            :: a
        write(*,'(A)') 'test_no_crossing_fallback'
        !$ call omp_set_num_threads(1)
        call seeded_curves(even,even,mask,FIXED_SEED,.false.,a)
        call assert_true(all(a%unmasked>0.99),'identical half maps never cross the onset threshold')
        call assert_real(0.,maxval(abs(a%corrected-a%unmasked)),TOL,'no-crossing fallback uses unmasked FSC')
        call assert_real(0.,maxval(abs(a%randomized)),0.,'no-crossing fallback never randomizes')
    end subroutine test_no_crossing_fallback

    subroutine make_fixture(even,odd,mask)
        class(image), intent(inout) :: even,odd,mask
        real, pointer               :: ep(:,:,:),op(:,:,:),mp(:,:,:)
        real                        :: x,y,z,r,signal
        integer                     :: i,j,k
        call even%new([BOX,BOX,BOX],1.,wthreads=.false.)
        call odd%new([BOX,BOX,BOX],1.,wthreads=.false.)
        call mask%new([BOX,BOX,BOX],1.,wthreads=.false.)
        call even%get_rmat_ptr(ep)
        call odd%get_rmat_ptr(op)
        call mask%get_rmat_ptr(mp)
        do k=1,BOX
            do j=1,BOX
                do i=1,BOX
                    x=real(i)-0.5*real(BOX+1)
                    y=real(j)-0.5*real(BOX+1)
                    z=real(k)-0.5*real(BOX+1)
                    r=sqrt(x*x+y*y+z*z)
                    signal=2.*exp(-(x*x+y*y+z*z)/36.)+exp(-((x-6.)**2+(y+3.)**2+z*z)/12.) &
                        +2.*exp(-((x+5.)**2+(y-4.)**2+(z+2.)**2)/4.)
                    ep(i,j,k)=signal+0.35*fixture_noise(i,j,k,17)
                    op(i,j,k)=signal+0.35*fixture_noise(i,j,k,7919)
                    mp(i,j,k)=0.
                    if(r<13.)then
                        mp(i,j,k)=1.
                    else if(r<17.)then
                        mp(i,j,k)=0.5*(1.+cos(PI*(r-13.)/4.))
                    endif
                enddo
            enddo
        enddo
    end subroutine make_fixture

    pure real function fixture_noise(i,j,k,offset) result(value)
        integer, intent(in) :: i,j,k,offset
        integer(int64)      :: n
        n=modulo(104729_int64*i+13007_int64*j+8191_int64*k+16381_int64*offset,2147483647_int64)
        n=modulo(n*n+104729_int64,2147483647_int64)
        n=modulo(n*48271_int64,2147483647_int64)
        value=real(real(n,dp)/2147483647._dp)-0.5
    end function fixture_noise

    subroutine seeded_curves(even,odd,mask,seed,reset_workers,curves)
        class(image), intent(in)      :: even,odd
        class(image), intent(inout)   :: mask
        integer, intent(in)           :: seed
        logical, intent(in)           :: reset_workers
        type(fsc_curves), intent(out) :: curves
        call set_fixed_seed(seed,propagate=.true.)
        if(reset_workers)call reset_worker_streams()
        call phase_rand_fsc(even,odd,mask,1,NSH,curves%corrected,curves%masked,curves%randomized,curves%unmasked)
    end subroutine seeded_curves

    !> A diagnostic control only: assign each intrinsic worker stream a known state.
    subroutine reset_worker_streams()
        integer, allocatable :: seed(:)
        integer              :: n,i,tid
        integer(int64)       :: worker_base
        worker_base=offset_seed(12345)
        !$omp parallel default(shared) private(seed,n,i,tid)
        tid=0
        !$ tid=omp_get_thread_num()
        call random_seed(size=n)
        allocate(seed(n))
        do i=1,n
            seed(i)=int(modulo(worker_base+104729_int64*tid+37_int64*i,int(huge(0),int64)))
        enddo
        call random_seed(put=seed)
        deallocate(seed)
        !$omp end parallel
    end subroutine reset_worker_streams

    subroutine assert_reference_curves(even,odd,mask,curves)
        class(image), intent(in)     :: even,odd
        class(image), intent(inout)  :: mask
        type(fsc_curves), intent(in) :: curves
        class(image), allocatable    :: e,o
        real                         :: expected(NSH)
        allocate(e,o)
        call e%copy(even)
        call o%copy(odd)
        call e%fft
        call o%fft
        call e%fsc(o,expected)
        call assert_real(0.,maxval(abs(expected-curves%unmasked)),TOL,'unmasked FSC matches direct half-map comparison')
        call e%ifft
        call o%ifft
        call e%zero_env_background(mask)
        call o%zero_env_background(mask)
        call e%mul(mask)
        call o%mul(mask)
        call e%fft
        call o%fft
        call e%fsc(o,expected)
        call assert_real(0.,maxval(abs(expected-curves%masked)),TOL,'masked FSC matches independent direct comparison')
        call e%kill
        call o%kill
        deallocate(e,o)
    end subroutine assert_reference_curves

    subroutine assert_repeat(a,b,label)
        type(fsc_curves), intent(in) :: a,b
        character(len=*), intent(in) :: label
        call assert_real(0.,maxval(abs(a%unmasked-b%unmasked)),TOL,label//': unmasked')
        call assert_real(0.,maxval(abs(a%masked-b%masked)),TOL,label//': masked')
        call assert_real(0.,maxval(abs(a%randomized-b%randomized)),TOL,label//': randomized')
        call assert_real(0.,maxval(abs(a%corrected-b%corrected)),TOL,label//': corrected')
    end subroutine assert_repeat

    subroutine report_repeat(label,threads,a,b)
        character(len=*), intent(in) :: label
        integer, intent(in)          :: threads
        type(fsc_curves), intent(in) :: a,b
        real                         :: res(NSH),r05,r143,s05,s143
        integer                      :: k
        do k=1,NSH
            res(k)=real(BOX)/real(k)
        enddo
        call get_resolution(a%corrected,res,r05,r143)
        call get_resolution(b%corrected,res,s05,s143)
        write(*,'(A,1X,A,1X,I0,4(1X,ES12.4))')'FSC_REPRO',label,threads, &
            maxval(abs(a%unmasked-b%unmasked)),maxval(abs(a%masked-b%masked)), &
            maxval(abs(a%randomized-b%randomized)),maxval(abs(a%corrected-b%corrected))
        write(*,'(A,1X,A,4(1X,F10.5))')'FSC_RESOLUTION',label,r05,s05,r143,s143
        do k=1,NSH
            write(*,'(A,1X,A,1X,I0,1X,I0,8(1X,ES16.8))')'FSC_SHELL',label,threads,k, &
                a%unmasked(k),b%unmasked(k),a%masked(k),b%masked(k), &
                a%randomized(k),b%randomized(k),a%corrected(k),b%corrected(k)
        enddo
    end subroutine report_repeat
end module simple_fsc_tester
