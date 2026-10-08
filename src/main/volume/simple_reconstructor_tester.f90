!@descr: unit tests of the weighted gridding insertion (simple_reconstructor)
! Linearity: a plane inserted with weight w equals the plane inserted with weight one and the
! accumulators scaled by w, to round-off; weight one is the unweighted insertion bit for bit.
module simple_reconstructor_tester
use simple_defs,            only: OSMPL_PAD_FAC
use simple_image,           only: image
use simple_ori,             only: ori
use simple_oris,            only: oris
use simple_parameters,      only: parameters
use simple_sp_project,      only: sp_project
use simple_sym,             only: sym
use simple_memoize_ft_maps, only: forget_ft_maps, memoize_ft_maps
use simple_type_defs,       only: CTFFLAG_NO, ctfparams, fplane_type
use simple_reconstructor,   only: reconstructor
use simple_test_utils
implicit none
private
public :: run_all_reconstructor_tests

integer, parameter :: BOX     = 16
integer, parameter :: NPLANES = 12
integer, parameter :: SEED    = 20261007
real,    parameter :: SMPD    = 2.0
real,    parameter :: RELTOL  = 1.e-6

contains

    subroutine run_all_reconstructor_tests()
        write(*,'(A)') '**** running all reconstructor tests ****'
        call test_weighted_insertion_linearity()
    end subroutine run_all_reconstructor_tests

    subroutine test_weighted_insertion_linearity()
        real, parameter :: WEIGHTS(3) = [0.37, 0.05, 0.81]
        type(parameters), target :: params
        type(sp_project)         :: project
        type(reconstructor)      :: rw, rs, rone, rplain
        type(oris)               :: poses
        type(ori)                :: o
        type(sym)                :: c1sym
        type(ctfparams)          :: ctfparms
        type(fplane_type)        :: fplane
        type(image)              :: obs, obs_pad
        complex, allocatable     :: cw(:,:,:), cs(:,:,:), cone(:,:,:), cplain(:,:,:)
        real,    allocatable     :: img(:,:,:), dw(:,:,:), ds(:,:,:), done(:,:,:), dplain(:,:,:)
        real    :: r, cerr, derr
        integer :: i, j, k, iw
        write(*,'(A)') 'test_weighted_insertion_linearity'
        params%box        = BOX
        params%box_crop   = BOX
        params%box_croppd = OSMPL_PAD_FAC * BOX
        params%smpd_crop  = SMPD
        params%nstates    = 1
        params%numlen     = 1
        params%oritype    = 'cls3D'
        call poses%new(NPLANES, .false.)
        call poses%spiral()
        call obs%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call obs_pad%new([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD, wthreads=.false.)
        call o%new(.false.)
        call c1sym%new('c1')
        call memoize_ft_maps([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD)
        ctfparms%smpd    = SMPD
        ctfparms%ctfflag = CTFFLAG_NO
        call set_fixed_seed(SEED)
        allocate(img(BOX,BOX,1))
        do iw = 1, size(WEIGHTS)
            call rw%new_accumulator(params, project, expand=.true., wthreads=.false.)
            call rs%new_accumulator(params, project, expand=.true., wthreads=.false.)
            call rone%new_accumulator(params, project, expand=.true., wthreads=.false.)
            call rplain%new_accumulator(params, project, expand=.true., wthreads=.false.)
            do i = 1, NPLANES
                do k = 1, BOX
                    do j = 1, BOX
                        call random_number(r)
                        img(j,k,1) = r - 0.5
                    enddo
                enddo
                call poses%get_ori(i, o)
                call obs%set_rmat(img, .false.)
                call obs%pad(obs_pad, backgr=0., antialiasing=.false.)
                call obs_pad%fft()
                call obs_pad%gen_fplane4rec([0,BOX/2], SMPD, ctfparms, [0.,0.], fplane)
                call rw%insert_plane_oversamp(c1sym, o, fplane, w=WEIGHTS(iw))
                call rs%insert_plane_oversamp(c1sym, o, fplane)
                call rone%insert_plane_oversamp(c1sym, o, fplane, w=1.0)
                call rplain%insert_plane_oversamp(c1sym, o, fplane)
            enddo
            ! one plane set per weight: scaling the unweighted sums by w must equal weighting each plane
            call rs%apply_weight(WEIGHTS(iw))
            call rw%compress_exp()
            call rs%compress_exp()
            call rone%compress_exp()
            call rplain%compress_exp()
            cw = rw%get_cmat()
            cs = rs%get_cmat()
            call rw%get_rho_copy(dw)
            call rs%get_rho_copy(ds)
            cerr = maxval(abs(cw - cs)) / max(maxval(abs(cs)), tiny(1.))
            derr = maxval(abs(dw - ds)) / max(maxval(abs(ds)), tiny(1.))
            call assert_true(cerr < RELTOL .and. derr < RELTOL, &
                &'insertion with weight w equals weight one scaled by w (data and density)')
            cone   = rone%get_cmat()
            cplain = rplain%get_cmat()
            call rone%get_rho_copy(done)
            call rplain%get_rho_copy(dplain)
            call assert_true(all(cone == cplain) .and. all(done == dplain), &
                &'insertion with weight one equals the unweighted insertion bit for bit')
            call rw%kill
            call rs%kill
            call rone%kill
            call rplain%kill
        enddo
        call forget_ft_maps()
        call c1sym%kill()
        call o%kill()
        call poses%kill()
        call obs%kill()
        call obs_pad%kill()
    end subroutine test_weighted_insertion_linearity

end module simple_reconstructor_tester
