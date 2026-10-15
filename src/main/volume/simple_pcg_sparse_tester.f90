!@descr: library tests of PCG reconstruction from sparse and uneven view sets against gridding (simple_reconstructor_pcg)
! An atomic phantom with strong Nyquist content is reconstructed from a single-axis tilt series, a preferred
! orientation and a dense control, through the production PCG path (raw accumulation, kernel operator, sphere
! support, two fixed iterations) and through production gridding. Unmeasured Fourier voxels must not be
! amplified: PCG must track gridding and keep the power beyond Nyquist small.
module simple_pcg_sparse_tester
use simple_defs,              only: dp, PI, OSMPL_PAD_FAC
use simple_image,             only: image
use simple_fplane_grid,       only: fplane_grid
use simple_ori,               only: ori
use simple_oris,              only: oris
use simple_reconstructor_pcg, only: reconstructor_pcg, PCG_OP_KERNEL, PCG_LAMBDA
use simple_reconstructor,     only: reconstructor
use simple_parameters,        only: parameters
use simple_sp_project,        only: sp_project
use simple_sym,               only: sym
use simple_memoize_ft_maps,   only: forget_ft_maps, memoize_ft_maps
use simple_type_defs,         only: CTFFLAG_NO, ctfparams, fplane_type
use simple_test_utils,        only: assert_true
use ieee_arithmetic,          only: ieee_is_finite
implicit none
private
public :: run_all_pcg_sparse_tests

integer,  parameter :: BOX               = 32
real,     parameter :: SMPD              = 1.0
real,     parameter :: ATOM_SIGMA        = 0.55   !< px: the atoms carry strong Nyquist content, as AET data do
real,     parameter :: ATOM_RAD          = 9.6    !< px: atoms lie within this radius
real,     parameter :: ATOM_DMIN         = 2.2    !< px: minimum atom spacing
integer,  parameter :: NATOMS            = 75
real,     parameter :: MSKRAD            = 13.
real,     parameter :: CORR_RAD          = 11.    !< px: correlation with the truth inside this sphere
integer,  parameter :: MAXITS            = 2      !< the production budget
real(dp), parameter :: SPARSE_CORR_SLACK = 0.05_dp
real(dp), parameter :: DENSE_CORR_SLACK  = 0.02_dp
real(dp), parameter :: BEYOND_NYQ_MAX    = 0.03_dp

contains

    subroutine run_all_pcg_sparse_tests()
        write(*,'(A)') '**** running all pcg_sparse tests ****'
        call test_atom_fixture()
        call test_view_set('tilt',  'single-axis tilt, 34 views from -65 to +64 degrees', SPARSE_CORR_SLACK)
        call test_view_set('cone',  'preferred orientation, 200 views within 20 degrees',  SPARSE_CORR_SLACK)
        call test_view_set('dense', 'dense control, 400 spiral views',                     DENSE_CORR_SLACK)
    end subroutine run_all_pcg_sparse_tests

    ! ---- fixture ---------------------------------------------------------------

    !> deterministic atom positions (3-D additive recurrence, no RNG) inside ATOM_RAD, ATOM_DMIN apart
    subroutine build_atoms( pos )
        real, allocatable, intent(out) :: pos(:,:)
        real(dp), parameter :: G = 1.2207440846057595_dp   ! the root of x**4 = x + 1
        real(dp) :: a(3), u(3)
        real     :: p(3), acc(3,NATOMS)
        integer  :: n, na
        a  = [1._dp/G, 1._dp/G**2, 1._dp/G**3]
        na = 0
        n  = 0
        do while( na < NATOMS .and. n < 100000 )
            n = n + 1
            u = modulo(0.5_dp + real(n,dp)*a, 1._dp)
            p = real((2._dp*u - 1._dp) * real(ATOM_RAD,dp))
            if( norm2(p) > ATOM_RAD ) cycle
            if( na > 0 )then
                if( any(norm2(acc(:,1:na) - spread(p,2,na), dim=1) < ATOM_DMIN) ) cycle
            endif
            na = na + 1
            acc(:,na) = p
        end do
        allocate(pos(3,na), source=acc(:,1:na))
    end subroutine build_atoms

    !> Gaussian atoms sampled at the voxel centres, origin at voxel BOX/2+1
    subroutine build_truth( pos, vol )
        real,              intent(in)  :: pos(:,:)
        real, allocatable, intent(out) :: vol(:,:,:)
        real    :: x(3), c
        integer :: i, j, k, ia
        c = real(BOX/2 + 1)
        allocate(vol(BOX,BOX,BOX), source=0.)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    x = [real(i), real(j), real(k)] - c
                    do ia = 1, size(pos,2)
                        vol(i,j,k) = vol(i,j,k) + exp(-sum((x - pos(:,ia))**2) / (2. * ATOM_SIGMA**2))
                    end do
                end do
            end do
        end do
    end subroutine build_truth

    !> exact line integrals of the Gaussian atoms: rows 1 and 2 of the rotation are the image x and y axes
    !! (loc = [h,k,0] R), so an atom at p projects to (p.R1, p.R2); no inverse crime with the operator
    subroutine project_atoms( pos, rot, img )
        real, intent(in)  :: pos(:,:), rot(3,3)
        real, intent(out) :: img(BOX,BOX)
        real, allocatable :: xy(:,:)
        real    :: c, amp
        integer :: i, j, ia
        c   = real(BOX/2 + 1)
        amp = sqrt(2. * PI) * ATOM_SIGMA
        allocate(xy(2,size(pos,2)))
        xy(1,:) = matmul(rot(1,:), pos)
        xy(2,:) = matmul(rot(2,:), pos)
        img = 0.
        do j = 1, BOX
            do i = 1, BOX
                do ia = 1, size(pos,2)
                    img(i,j) = img(i,j) + amp * exp(-((real(i) - c - xy(1,ia))**2 + (real(j) - c - xy(2,ia))**2) &
                        &/ (2. * ATOM_SIGMA**2))
                end do
            end do
        end do
    end subroutine project_atoms

    subroutine build_view_set( kind, os )
        character(len=*), intent(in)    :: kind
        type(oris),       intent(inout) :: os
        integer :: i, n
        select case( kind )
            case( 'tilt' )
                n = 34
                call os%new(n, .false.)
                do i = 1, n
                    call os%set_euler(i, [0., -65. + real(i-1) * 129. / real(n-1), 0.])
                end do
            case( 'cone' )
                n = 200
                call os%new(n, .false.)
                do i = 1, n
                    call os%set_euler(i, [modulo(real(i-1) * 137.508, 360.), 20. * sqrt((real(i) - 0.5) / real(n)), &
                        &modulo(real(i-1) * 29., 360.)])
                end do
            case( 'dense' )
                call os%new(400, .false.)
                call os%spiral()
            case default
                error stop 'build_view_set: unknown view set'
        end select
    end subroutine build_view_set

    ! ---- reconstructions -------------------------------------------------------

    !> the production PCG path: raw accumulation, sphere support, kernel operator, MAXITS fixed iterations from zero
    subroutine reconstruct_pcg( os, images, vol, niters )
        type(oris),        intent(inout) :: os
        real,              intent(in)    :: images(:,:,:)
        real, allocatable, intent(out)   :: vol(:,:,:)
        integer,           intent(out)   :: niters
        class(reconstructor_pcg), allocatable :: pcgop
        complex, allocatable :: planes(:,:,:)
        type(image) :: proj
        integer     :: i, n, lims2(2,2)
        n = os%get_noris()
        allocate(pcgop)
        call pcgop%new(BOX, SMPD, PCG_LAMBDA)
        lims2 = pcgop%get_lims2()
        allocate(planes(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),n))
        call proj%new([BOX,BOX,1], SMPD, wthreads=.false.)
        do i = 1, n
            call proj%set_rmat(images(:,:,i:i), .false.)
            call proj%fft()
            planes(:,:,i) = pcgop%extract_native_plane(proj)
        end do
        call proj%kill()
        call pcgop%set_window_sphere(MSKRAD)
        call pcgop%prep_particles(os, use_ctf=.false.)
        call pcgop%begin_accum
        call pcgop%accumulate_batch(planes, n, 1)
        call pcgop%end_accum(.true.)
        call pcgop%set_op_mode(PCG_OP_KERNEL)
        allocate(vol(BOX,BOX,BOX), source=0.)
        call pcgop%solve_accum(vol, maxits=MAXITS, rtol=0., niters=niters)
        call pcgop%kill()
        deallocate(pcgop)
    end subroutine reconstruct_pcg

    !> production gridding (KB insertion, no CTF), shipped under the soft window the PCG map carries
    subroutine reconstruct_gridding( os, images, vol )
        type(oris),        intent(inout) :: os
        real,              intent(in)    :: images(:,:,:)
        real, allocatable, intent(out)   :: vol(:,:,:)
        class(parameters),    allocatable, target :: params
        class(reconstructor), allocatable :: gridder
        class(sp_project),    allocatable :: project
        type(ctfparams)   :: ctfparms
        type(fplane_type) :: fplane
        type(fplane_grid) :: grid
        type(image)       :: observation, observation_padded, restored
        type(ori)         :: o
        type(sym)         :: c1sym
        integer :: i
        allocate(params, gridder, project)
        params%box        = BOX
        params%box_crop   = BOX
        params%box_croppd = OSMPL_PAD_FAC * BOX
        params%smpd_crop  = SMPD
        params%nstates    = 1
        params%numlen     = 1
        params%oritype    = 'cls3D'   ! class-volume field: no project dependency
        call gridder%new_accumulator(params, project, expand=.true., wthreads=.false.)
        call gridder%set_ft(.true.)
        call restored%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call observation%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call observation_padded%new([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD, wthreads=.false.)
        call grid%new(OSMPL_PAD_FAC*BOX)
        call o%new(.false.)
        call c1sym%new('c1')
        call memoize_ft_maps([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD)
        ctfparms%smpd    = SMPD
        ctfparms%ctfflag = CTFFLAG_NO
        do i = 1, os%get_noris()
            call os%get_ori(i, o)
            call observation%set_rmat(images(:,:,i:i), .false.)
            call observation%pad(observation_padded, backgr=0., antialiasing=.false.)
            call observation_padded%fft()
            call observation_padded%gen_fplane4rec(SMPD, ctfparms, [0.,0.], fplane, grid)
            call gridder%insert_plane_oversamp(c1sym, o, fplane)
        end do
        call gridder%compress_exp()
        call gridder%restore_final(restored)
        call restored%mask3D_soft(MSKRAD, backgr=0.)
        vol = restored%get_rmat()
        if( allocated(fplane%cmplx_plane)    ) deallocate(fplane%cmplx_plane)
        if( allocated(fplane%ctfsq_plane)    ) deallocate(fplane%ctfsq_plane)
        if( allocated(fplane%transfer_plane) ) deallocate(fplane%transfer_plane)
        call forget_ft_maps()
        call c1sym%kill()
        call o%kill()
        call restored%kill()
        call observation_padded%kill()
        call observation%kill()
        call gridder%kill()
        call project%kill()
        deallocate(gridder, project, params)
    end subroutine reconstruct_gridding

    ! ---- measures --------------------------------------------------------------

    !> correlation inside the CORR_RAD sphere about voxel BOX/2+1
    real(dp) function corr_sphere( a, b )
        real, intent(in) :: a(:,:,:), b(:,:,:)
        logical, allocatable :: inside(:,:,:)
        real(dp) :: am, bm, n, sab, saa, sbb
        real     :: c
        integer  :: i, j, k
        c = real(BOX/2 + 1)
        allocate(inside(BOX,BOX,BOX))
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    inside(i,j,k) = (real(i)-c)**2 + (real(j)-c)**2 + (real(k)-c)**2 <= CORR_RAD**2
                end do
            end do
        end do
        n   = real(count(inside),dp)
        am  = sum(real(a,dp), mask=inside) / n
        bm  = sum(real(b,dp), mask=inside) / n
        sab = sum((real(a,dp)-am) * (real(b,dp)-bm), mask=inside)
        saa = sum((real(a,dp)-am)**2, mask=inside)
        sbb = sum((real(b,dp)-bm)**2, mask=inside)
        corr_sphere = sab / max(sqrt(saa*sbb), tiny(1._dp))
    end function corr_sphere

    !> fraction of the Fourier power beyond the Nyquist sphere, where no view has data
    real(dp) function beyond_nyquist_fraction( vol )
        real, intent(in) :: vol(:,:,:)
        type(image) :: img
        complex, allocatable :: cmat(:,:,:)
        real(dp) :: w, p, total, beyond
        integer  :: cdim(3), logi(3), i, j, k, nyq
        nyq = BOX / 2
        call img%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call img%set_rmat(vol, .false.)
        call img%fft()
        cmat = img%get_cmat()
        cdim = img%get_array_shape()
        total  = 0._dp
        beyond = 0._dp
        do k = 1, cdim(3)
            do j = 1, cdim(2)
                do i = 1, cdim(1)
                    logi = img%comp_addr_logi(i,j,k)
                    ! the packed half holds h >= 0; h = 0 and h = Nyquist have no stored mate
                    w = 2._dp
                    if( logi(1) == 0 .or. logi(1) == nyq ) w = 1._dp
                    p = w * real(abs(cmat(i,j,k)),dp)**2
                    total = total + p
                    if( sum(logi**2) > nyq*nyq ) beyond = beyond + p
                end do
            end do
        end do
        call img%kill()
        beyond_nyquist_fraction = beyond / max(total, tiny(1._dp))
    end function beyond_nyquist_fraction

    !> |centre voxel| over max |map| inside the correlation sphere: the solver-made spike reads ~1
    real function centre_spike( vol )
        real, intent(in) :: vol(:,:,:)
        integer :: c
        c = BOX/2 + 1
        centre_spike = abs(vol(c,c,c)) / max(maxval(abs(vol(c-nint(CORR_RAD):c+nint(CORR_RAD), &
            &c-nint(CORR_RAD):c+nint(CORR_RAD),c-nint(CORR_RAD):c+nint(CORR_RAD)))), tiny(1.))
    end function centre_spike

    ! ---- tests -----------------------------------------------------------------

    !> the phantom is reproducible, holds the requested atoms and has Nyquist content to lose
    subroutine test_atom_fixture()
        real, allocatable :: pos(:,:), again(:,:), truth(:,:,:)
        write(*,'(A)') 'test_atom_fixture'
        call build_atoms(pos)
        call build_atoms(again)
        call assert_true(size(pos,2) == NATOMS, 'the recurrence places all atoms')
        call assert_true(all(pos == again), 'the atom positions are bit-identical when regenerated')
        call build_truth(pos, truth)
        call assert_true(all(ieee_is_finite(truth)) .and. minval(truth) >= 0., 'the atomic phantom is finite and non-negative')
        ! the phantom itself is sampled on the lattice, so its own power beyond Nyquist is aliasing
        write(*,'(A,F8.5)') '  phantom power beyond Nyquist: ', beyond_nyquist_fraction(truth)
    end subroutine test_atom_fixture

    !> one view set through both backends; PCG must track gridding and keep the Nyquist sphere clean
    subroutine test_view_set( kind, label, corr_slack )
        character(len=*), intent(in) :: kind, label
        real(dp),         intent(in) :: corr_slack
        type(oris) :: os
        real, allocatable :: pos(:,:), truth(:,:,:), images(:,:,:), vpcg(:,:,:), vgrid(:,:,:)
        real(dp) :: corr_pcg, corr_grid, beyond_pcg, beyond_grid
        integer  :: i, n, niters
        write(*,'(A)') 'test_view_set: '//trim(label)
        call build_atoms(pos)
        call build_truth(pos, truth)
        call build_view_set(kind, os)
        n = os%get_noris()
        allocate(images(BOX,BOX,n))
        do i = 1, n
            call project_atoms(pos, os%get_mat(i), images(:,:,i))
        end do
        call reconstruct_pcg(os, images, vpcg, niters)
        call reconstruct_gridding(os, images, vgrid)
        call assert_true(niters == MAXITS, 'the PCG solve runs the production budget of two iterations')
        call assert_true(all(ieee_is_finite(vpcg)) .and. all(ieee_is_finite(vgrid)), 'both maps are finite')
        corr_pcg    = corr_sphere(vpcg,  truth)
        corr_grid   = corr_sphere(vgrid, truth)
        beyond_pcg  = beyond_nyquist_fraction(vpcg)
        beyond_grid = beyond_nyquist_fraction(vgrid)
        write(*,'(A,2F8.4)') '  correlation with the truth, PCG / gridding: ', corr_pcg, corr_grid
        write(*,'(A,2F8.5)') '  power beyond Nyquist, PCG / gridding:       ', beyond_pcg, beyond_grid
        write(*,'(A,2F8.3)') '  centre / max inside the sphere, PCG / grid: ', centre_spike(vpcg), centre_spike(vgrid)
        call assert_true(corr_pcg >= corr_grid - corr_slack, trim(kind)//': PCG correlates with the truth as well as gridding')
        call assert_true(beyond_pcg <= BEYOND_NYQ_MAX, trim(kind)//': PCG keeps the power beyond Nyquist under 3%')
        call os%kill()
    end subroutine test_view_set

end module simple_pcg_sparse_tester
