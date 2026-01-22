program simple_test_msk_routines
    use simple_core_module_api
    use simple_image
    use simple_imgarr_utils, only: write_imgarr
    implicit none
    type(image) :: img
    type(image), allocatable :: stk(:)
    integer :: i, nimgs
    integer :: ldim2(3), ldim3(3)
    real    :: smpd
    smpd  = 3.0
    nimgs = 8
    ! -----------------------------
    ! POSITIVE TESTS (2D serial)
    ! -----------------------------
    ldim2 = [256,256,1]
    call img%new(ldim2, smpd)
    call img%ran()
    call img%memoize_mask_coords
    call img%mask2D_soft(ldim2(1)/3.0)
    call img%mask2D_softavg(ldim2(1)/3.0)
    call img%mask2D_hard(ldim2(1)/3.0)
    write(*,*) "OK: 2D serial"

    ! -----------------------------
    ! POSITIVE TESTS (2D parallel)
    ! -----------------------------
    call unmemoize_mask_coords()
    allocate(stk(nimgs))
    do i = 1, nimgs
        call stk(i)%new(ldim2, smpd)
        call stk(i)%ran()
    end do

    ! Ensure memoization is established outside parallel
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask2D_soft(ldim2(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 2D parallel soft"
    call write_imgarr(stk, string('mask2D_soft_openmp.mrcs'))

    do i = 1, nimgs
        call stk(i)%new(ldim2, smpd)
        call stk(i)%ran()
    end do
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask2D_softavg(ldim2(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 2D parallel softavg"
    call write_imgarr(stk, string('mask2D_softavg_openmp.mrc'))

    do i = 1, nimgs
        call stk(i)%new(ldim2, smpd)
        call stk(i)%ran()
    end do
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask2D_hard(ldim2(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 2D parallel hard"
    call write_imgarr(stk, string('mask2D_hard_openmp.mrc'))

    ! -----------------------------
    ! POSITIVE TESTS (3D serial)
    ! -----------------------------
    ldim3 = [64,64,64]
    call img%new(ldim3, smpd)
    call img%ran()
    call img%memoize_mask_coords
    call img%mask3D_soft(ldim3(1)/3.0)
    call img%mask3D_softavg(ldim3(1)/3.0)
    call img%mask3D_hard(ldim3(1)/3.0)
    write(*,*) "OK: 3D serial"

    ! -----------------------------
    ! POSITIVE TESTS (3D parallel)
    ! -----------------------------
    call unmemoize_mask_coords()
    do i = 1, nimgs
        call stk(i)%new(ldim3, smpd)
        call stk(i)%ran()
    end do

    ! Pre-memoize outside parallel so 3D calls don't attempt memoize inside parallel
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask3D_soft(ldim3(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 3D parallel soft"

    do i = 1, nimgs
        call stk(i)%new(ldim3, smpd)
        call stk(i)%ran()
    end do
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask3D_softavg(ldim3(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 3D parallel softavg"

    do i = 1, nimgs
        call stk(i)%new(ldim3, smpd)
        call stk(i)%ran()
    end do
    call stk(1)%memoize_mask_coords
    !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
    do i = 1, nimgs
        call stk(i)%mask3D_hard(ldim3(1)/3.0)
    end do
    !$omp end parallel do
    write(*,*) "OK: 3D parallel hard"

    write(*,*) "ALL MASK TESTS PASSED"

end program simple_test_msk_routines





