!@descr: unit tests of the pcg_lattice plane gather and scatter: adjoint, Ewald curvature, thread identity
! box 24 (padded 48). The scatter is the transpose of the gather to 1e-5 relative, flat and curved, with
! point-group replication. The curved gather sits at the paraboloid's location padf [h, k, kz] R to 1e-5
! (against the KB interpolation at the location the formula gives, identity and random rotation). A
! point displaced by z0 = 5 along the beam:
! relative to the flat plane it carries the phase exp(s 2 pi i kz z0 / box) to 0.05 rad, s the sign an
! in-plane displacement shows, at kz up to 2.1 voxels. Gather and scatter are bitwise identical at
! 1, 11, 16 and 24 threads, each run twice. The cube of multislice: its scatter is the transpose of its gather
! to 1e-5 with point-group replication, its plane l = 0 is the flat plane gather bitwise, and both are bitwise
! identical at 1, 11, 16 and 24 threads, twice each.
module simple_pcg_lattice_tester
use simple_defs,        only: dp, PI
use simple_image,       only: image
use simple_pcg_lattice, only: pcg_lattice
use simple_cartesian_fourier, only: gather_packed_window
use simple_test_utils,  only: assert_true, fix_test_seed
!$ use omp_lib,          only: omp_get_max_threads, omp_set_num_threads
implicit none
private
public :: run_all_pcg_lattice_tests

integer, parameter :: BOX = 24
real,    parameter :: SMPD = 1.0
real,    parameter :: WL_TEST = 0.7     !< a wavelength in pixels large enough to curve the paraboloid by 2.1 voxels

contains

    subroutine run_all_pcg_lattice_tests()
        write(*,'(A)') '**** running all PCG lattice gather/scatter tests ****'
        call fix_test_seed()
        call test_adjoint()
        call test_curvature_phase()
        call test_thread_identity()
        call test_cube()
    end subroutine run_all_pcg_lattice_tests

    !> <G x, q> = <x, G^T q> for flat and curved planes, two replicated rotations
    subroutine test_adjoint()
        class(pcg_lattice), allocatable :: lat
        complex, allocatable :: cmat(:,:,:), accum(:,:,:), plane(:,:), q(:,:)
        real     :: rots(3,3,2), curv
        real(dp) :: lhs, rhs
        integer  :: ic, g
        allocate(lat)
        call lat%new_lattice(BOX, SMPD)
        cmat = random_lattice(lat)
        rots(:,:,1) = random_rotation()
        rots(:,:,2) = random_rotation()
        allocate(plane(lat%lims2(1,1):lat%lims2(1,2), lat%lims2(2,1):lat%lims2(2,2)))
        q = random_plane(lat)
        do ic = 1, 3
            select case(ic)
                case(1); curv = 0.
                case(2); curv = lat%ewald_curv(-1, WL_TEST)
                case(3); curv = lat%ewald_curv( 1, WL_TEST)
            end select
            lhs = 0.0_dp
            do g = 1, 2
                call lat%gather_plane(cmat, rots(:,:,g), curv, plane)
                lhs = lhs + plane_dot(lat, plane, q)
            end do
            allocate(accum(lat%lims3(1,1):lat%lims3(1,2), lat%lims3(2,1):lat%lims3(2,2), &
                &lat%lims3(3,1):lat%lims3(3,2)), source=cmplx(0.,0.))
            call lat%scatter_plane(q, rots, curv, accum)
            rhs = lattice_dot(lat, cmat, accum)
            deallocate(accum)
            call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), &
                &'the plane scatter is the transpose of the plane gather, flat and on the paraboloid')
        end do
        call lat%kill_lattice
        deallocate(lat)
    end subroutine test_adjoint

    !> a point at (0,0,z0) seen on the paraboloid. Location: the gather equals the KB interpolation of the
    !! lattice at padf*[h, k, kappa wl (h^2+k^2) / (2 box)], computed here from the formula, to 1e-5
    !! (rounding only). Phase: relative to the flat plane it is exp(s 2 pi i kz z0 / box), s the sign of an
    !! in-plane displacement, to 0.05 rad; the KB interpolation of this off-centre point errs by up to 1%
    !! between lattice points, while a wrong sign errs by twice the phase (up to 5.5 rad here) and a 2%
    !! error in the size by 0.055 rad at the rim
    subroutine test_curvature_phase()
        integer, parameter :: Z0 = 5, X0 = 5
        class(pcg_lattice), allocatable :: lat
        complex, allocatable :: cz(:,:,:), cx(:,:,:), c0(:,:,:), pflat(:,:), pcurv(:,:), px(:,:), p0(:,:)
        real    :: ident(3,3), rot(3,3), curv, kz, s, loc(3), w(15,15,15), dphase, err_loc, err_phase
        complex :: expected
        integer :: h, k, kappa, i0(3)
        allocate(lat)
        call lat%new_lattice(BOX, SMPD)
        if( lat%wdim > size(w,1) ) error stop 'test_curvature_phase: KB window larger than the work stencil'
        ident = 0.; ident(1,1) = 1.; ident(2,2) = 1.; ident(3,3) = 1.
        c0 = point_lattice(lat, [0, 0, 0])
        cx = point_lattice(lat, [X0, 0, 0])
        cz = point_lattice(lat, [0, 0, Z0])
        allocate(pflat(lat%lims2(1,1):lat%lims2(1,2), lat%lims2(2,1):lat%lims2(2,2)))
        allocate(pcurv, px, p0, mold=pflat)
        ! the sign of the transform: the x-displaced point in the flat plane, at (1,0)
        call lat%gather_plane(c0, ident, 0., p0)
        call lat%gather_plane(cx, ident, 0., px)
        s = atan2(aimag(px(1,0)/p0(1,0)), real(px(1,0)/p0(1,0))) / (2.0*PI*real(X0)/real(BOX))
        call assert_true(abs(abs(s) - 1.0) < 1.0e-3, 'an in-plane displacement is a unit phase ramp')
        s = sign(1.0, s)
        call lat%gather_plane(cz, ident, 0., pflat)
        err_loc   = 0.
        err_phase = 0.
        do kappa = -1, 1, 2
            curv = lat%ewald_curv(kappa, WL_TEST)
            call lat%gather_plane(cz, ident, curv, pcurv)
            do k = lat%lims2(2,1), lat%lims2(2,2)
                do h = lat%lims2(1,1), lat%lims2(1,2)
                    if( h*h + k*k > lat%sqlp ) cycle
                    ! the location from the formula, independently of the lattice's own evaluation
                    kz  = real(kappa) * WL_TEST * real(h*h + k*k) / (2.0 * real(BOX))
                    loc = real(lat%padf) * [real(h), real(k), kz]
                    i0  = nint(loc) - lat%iwinsz
                    call lat%kbwin%apod_mat_3d_fast(loc, lat%iwinsz, lat%wdim, w(:lat%wdim,:lat%wdim,:lat%wdim))
                    expected = lat%padsc * gather_packed_window(cz, lbound(lat%wrap,1), lat%wrap, i0, &
                        &w(:lat%wdim,:lat%wdim,:lat%wdim))
                    err_loc = max(err_loc, abs(pcurv(h,k) - expected) / abs(pflat(h,k)))
                    ! the phase the paraboloid adds to the flat plane
                    dphase = atan2(aimag(pcurv(h,k)/pflat(h,k)), real(pcurv(h,k)/pflat(h,k))) &
                        &- s * 2.0 * PI * kz * real(Z0) / real(BOX)
                    dphase = modulo(dphase + PI, 2.0*PI) - PI
                    err_phase = max(err_phase, abs(dphase))
                end do
            end do
        end do
        ! the location under a random rotation: padf [h, k, kz] R, rows times the rotation as for the flat plane
        rot  = random_rotation()
        curv = lat%ewald_curv(-1, WL_TEST)
        call lat%gather_plane(cx, rot, curv, pcurv)
        do k = lat%lims2(2,1), lat%lims2(2,2)
            do h = lat%lims2(1,1), lat%lims2(1,2)
                if( h*h + k*k > lat%sqlp ) cycle
                kz  = -WL_TEST * real(h*h + k*k) / (2.0 * real(BOX))
                loc = real(lat%padf) * matmul([real(h), real(k), kz], rot)
                i0  = nint(loc) - lat%iwinsz
                call lat%kbwin%apod_mat_3d_fast(loc, lat%iwinsz, lat%wdim, w(:lat%wdim,:lat%wdim,:lat%wdim))
                expected = lat%padsc * gather_packed_window(cx, lbound(lat%wrap,1), lat%wrap, i0, &
                    &w(:lat%wdim,:lat%wdim,:lat%wdim))
                err_loc = max(err_loc, abs(pcurv(h,k) - expected) / maxval(abs(pcurv)))
            end do
        end do
        write(*,'(A,ES10.3,A,ES10.3,A)') '    paraboloid: location error ', err_loc, ' (relative), phase error ', &
            &err_phase, ' rad'
        call assert_true(err_loc <= 1.0e-5, &
            &'the paraboloid points sit at padf [h, k, kappa wl (h^2+k^2) / (2 box)] R, flat and rotated')
        call assert_true(err_phase <= 0.05, 'the paraboloid adds the phase of its height, with its sign')
        call assert_true(abs(lat%ewald_curv(0, WL_TEST)) == 0., 'ewald=no is the central section')
        call lat%kill_lattice
        deallocate(lat)
    end subroutine test_curvature_phase

    !> gather and scatter on the paraboloid at 1, 11, 16 and 24 threads, twice each: bitwise equal
    subroutine test_thread_identity()
        integer, parameter :: NTHR(4) = [1, 11, 16, 24]
        class(pcg_lattice), allocatable :: lat
        complex, allocatable :: cmat(:,:,:), q(:,:), plane(:,:), plane_ref(:,:), accum(:,:,:), accum_ref(:,:,:)
        real    :: rots(3,3,3), curv
        integer :: it, irep, nthr0
        logical :: l_same
        allocate(lat)
        call lat%new_lattice(BOX, SMPD)
        cmat = random_lattice(lat)
        q    = random_plane(lat)
        rots(:,:,1) = random_rotation()
        rots(:,:,2) = random_rotation()
        rots(:,:,3) = random_rotation()
        curv = lat%ewald_curv(-1, WL_TEST)
        allocate(plane(lat%lims2(1,1):lat%lims2(1,2), lat%lims2(2,1):lat%lims2(2,2)))
        allocate(plane_ref, mold=plane)
        allocate(accum(lat%lims3(1,1):lat%lims3(1,2), lat%lims3(2,1):lat%lims3(2,2), &
            &lat%lims3(3,1):lat%lims3(3,2)))
        allocate(accum_ref, mold=accum)
        nthr0 = 1
        !$ nthr0 = omp_get_max_threads()
        l_same = .true.
        do it = 1, size(NTHR)
            do irep = 1, 2
                !$ call omp_set_num_threads(NTHR(it))
                call lat%gather_plane(cmat, rots(:,:,1), curv, plane)
                accum = cmplx(0.,0.)
                call lat%scatter_plane(q, rots, curv, accum)
                if( it == 1 .and. irep == 1 )then
                    plane_ref = plane
                    accum_ref = accum
                else
                    l_same = l_same .and. all(plane == plane_ref) .and. all(accum == accum_ref)
                endif
            end do
        end do
        !$ call omp_set_num_threads(nthr0)
        call assert_true(l_same, 'gather and scatter on the paraboloid are identical at 1, 11, 16 and 24 threads')
        call lat%kill_lattice
        deallocate(lat)
    end subroutine test_thread_identity

    !> the cube of multislice: adjoint, central plane, thread identity
    subroutine test_cube()
        integer, parameter :: NTHR(4) = [1, 11, 16, 24]
        class(pcg_lattice), allocatable :: lat
        complex, allocatable :: cmat(:,:,:), cube(:,:,:), q(:,:,:), cube_ref(:,:,:), plane(:,:), accum(:,:,:), &
            &accum_ref(:,:,:)
        real     :: rots(3,3,2), re
        real(dp) :: lhs, rhs
        integer  :: h, k, l, g, rl, it, irep, nthr0
        logical  :: l_same
        allocate(lat)
        call lat%new_lattice(BOX, SMPD)
        cmat = random_lattice(lat)
        rots(:,:,1) = random_rotation()
        rots(:,:,2) = random_rotation()
        rl = lat%lims2(1,2) - 1
        allocate(cube(lat%lims2(1,1):lat%lims2(1,2), lat%lims2(2,1):lat%lims2(2,2), -rl:rl))
        allocate(q, cube_ref, mold=cube)
        do l = -rl, rl
            do k = lat%lims2(2,1), lat%lims2(2,2)
                do h = lat%lims2(1,1), lat%lims2(1,2)
                    call random_number(re)
                    q(h,k,l) = cmplx(re - 0.5, 0.)
                    call random_number(re)
                    q(h,k,l) = q(h,k,l) + cmplx(0., re - 0.5)
                    if( h*h + k*k + l*l > lat%sqlp ) q(h,k,l) = cmplx(0.,0.)
                end do
            end do
        end do
        allocate(accum(lat%lims3(1,1):lat%lims3(1,2), lat%lims3(2,1):lat%lims3(2,2), &
            &lat%lims3(3,1):lat%lims3(3,2)), source=cmplx(0.,0.))
        allocate(accum_ref, mold=accum)
        lhs = 0.0_dp
        do g = 1, 2
            call lat%gather_cube(cmat, rots(:,:,g), cube)
            lhs = lhs + sum(real(cube,dp)*real(q,dp) + real(aimag(cube),dp)*real(aimag(q),dp))
        end do
        call lat%scatter_cube(q, rots, accum)
        rhs = lattice_dot(lat, cmat, accum)
        call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), &
            &'the cube scatter is the transpose of the cube gather')
        ! the central plane of the cube is the flat plane
        allocate(plane(lat%lims2(1,1):lat%lims2(1,2), lat%lims2(2,1):lat%lims2(2,2)))
        call lat%gather_cube(cmat, rots(:,:,1), cube)
        call lat%gather_plane(cmat, rots(:,:,1), 0., plane)
        call assert_true(all(cube(:,:,0) == plane), 'the plane l = 0 of the cube is the flat plane gather, bitwise')
        ! threads
        nthr0 = 1
        !$ nthr0 = omp_get_max_threads()
        l_same = .true.
        do it = 1, size(NTHR)
            do irep = 1, 2
                !$ call omp_set_num_threads(NTHR(it))
                call lat%gather_cube(cmat, rots(:,:,2), cube)
                accum = cmplx(0.,0.)
                call lat%scatter_cube(q, rots, accum)
                if( it == 1 .and. irep == 1 )then
                    cube_ref  = cube
                    accum_ref = accum
                else
                    l_same = l_same .and. all(cube == cube_ref) .and. all(accum == accum_ref)
                endif
            end do
        end do
        !$ call omp_set_num_threads(nthr0)
        call assert_true(l_same, 'the cube gather and scatter are identical at 1, 11, 16 and 24 threads')
        call lat%kill_lattice
        deallocate(lat)
    end subroutine test_cube

    ! ---- helpers ----

    !> the packed transform of a random real volume on the padded lattice
    function random_lattice( lat ) result( cmat )
        class(pcg_lattice), intent(in) :: lat
        complex, allocatable :: cmat(:,:,:)
        type(image) :: vol
        real, allocatable :: r(:,:,:)
        allocate(r(lat%boxpd,lat%boxpd,lat%boxpd))
        call random_number(r)
        call vol%new([lat%boxpd,lat%boxpd,lat%boxpd], SMPD)
        call vol%set_rmat(r - 0.5, .false.)
        call vol%fft
        cmat = vol%get_cmat()
        call vol%kill
    end function random_lattice

    !> the packed transform of a unit point at voxel offset off from the padded lattice's centre
    function point_lattice( lat, off ) result( cmat )
        class(pcg_lattice), intent(in) :: lat
        integer,            intent(in) :: off(3)
        complex, allocatable :: cmat(:,:,:)
        type(image) :: vol
        real, allocatable :: r(:,:,:)
        integer :: c
        allocate(r(lat%boxpd,lat%boxpd,lat%boxpd), source=0.)
        c = lat%boxpd/2 + 1
        r(c+off(1), c+off(2), c+off(3)) = 1.
        call vol%new([lat%boxpd,lat%boxpd,lat%boxpd], SMPD)
        call vol%set_rmat(r, .false.)
        call vol%fft
        cmat = vol%get_cmat()
        call vol%kill
    end function point_lattice

    function random_plane( lat ) result( q )
        class(pcg_lattice), intent(in) :: lat
        complex, allocatable :: q(:,:)
        real, allocatable :: re(:,:), im(:,:)
        integer :: h, k
        allocate(re(lat%lims2(1,1):lat%lims2(1,2), lat%lims2(2,1):lat%lims2(2,2)))
        allocate(im, mold=re)
        call random_number(re)
        call random_number(im)
        ! allocated with the disk's bounds: an expression would give q lower bounds of 1
        allocate(q(lat%lims2(1,1):lat%lims2(1,2), lat%lims2(2,1):lat%lims2(2,2)))
        q(:,:) = cmplx(re - 0.5, im - 0.5)
        do k = lat%lims2(2,1), lat%lims2(2,2)
            do h = lat%lims2(1,1), lat%lims2(1,2)
                if( h*h + k*k > lat%sqlp ) q(h,k) = cmplx(0.,0.)
            end do
        end do
    end function random_plane

    function random_rotation() result( rot )
        real :: rot(3,3), a(3), ca(3), sa(3), rz1(3,3), ry(3,3), rz2(3,3)
        call random_number(a)
        a  = a * [2.0*PI, PI, 2.0*PI]
        ca = cos(a); sa = sin(a)
        rz1 = reshape([ca(1), -sa(1), 0., sa(1), ca(1), 0., 0., 0., 1.], [3,3])
        ry  = reshape([ca(2), 0., sa(2), 0., 1., 0., -sa(2), 0., ca(2)], [3,3])
        rz2 = reshape([ca(3), -sa(3), 0., sa(3), ca(3), 0., 0., 0., 1.], [3,3])
        rot = matmul(rz2, matmul(ry, rz1))
    end function random_rotation

    real(dp) function plane_dot( lat, a, b )
        class(pcg_lattice), intent(in) :: lat
        complex,            intent(in) :: a(lat%lims2(1,1):,lat%lims2(2,1):), b(lat%lims2(1,1):,lat%lims2(2,1):)
        plane_dot = sum(real(conjg(a)*b, dp))
    end function plane_dot

    !> Re <x_hat, accum> over the full-range lattice, x_hat unpacked from the half-complex cmat by
    !! Friedel symmetry (an independent reading of the packed layout)
    real(dp) function lattice_dot( lat, cmat, accum )
        class(pcg_lattice), intent(in) :: lat
        complex,            intent(in) :: cmat(:,:,:)
        complex,            intent(in) :: accum(lat%lims3(1,1):,lat%lims3(2,1):,lat%lims3(3,1):)
        complex :: xv
        integer :: h, k, m, n
        n = lat%boxpd
        lattice_dot = 0.0_dp
        do m = lat%lims3(3,1), lat%lims3(3,2)
            do k = lat%lims3(2,1), lat%lims3(2,2)
                do h = lat%lims3(1,1), lat%lims3(1,2)
                    if( h >= 0 )then
                        xv = cmat(h+1, modulo(k,n)+1, modulo(m,n)+1)
                    else
                        xv = conjg(cmat(-h+1, modulo(-k,n)+1, modulo(-m,n)+1))
                    endif
                    lattice_dot = lattice_dot + real(conjg(xv)*accum(h,k,m), dp)
                end do
            end do
        end do
    end function lattice_dot

end module simple_pcg_lattice_tester
