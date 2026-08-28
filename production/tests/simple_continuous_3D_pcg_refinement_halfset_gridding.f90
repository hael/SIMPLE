module continuous_3D_pcg_refinement_halfset_gridding
use continuous_3D_pcg_refinement_halfset_support, only: HALFSET_SMPD
use continuous_3D_pcg_refinement_test_helpers, only: BOX => TRUTH_VOLUME_BOX
use simple_defs, only: OSMPL_PAD_FAC
use simple_image, only: image
use simple_memoize_ft_maps, only: forget_ft_maps, memoize_ft_maps
use simple_ori, only: ori
use simple_oris, only: oris
use simple_parameters, only: parameters
use simple_reconstructor, only: reconstructor
use simple_sp_project, only: sp_project
use simple_sym, only: sym
use simple_type_defs, only: CTFFLAG_NO, ctfparams, fplane_type
implicit none
private

public :: reconstruct_half_conventionally

contains

!> Reconstruct one half with SIMPLE's conventional gridding path for an A/B baseline.
subroutine reconstruct_half_conventionally(orientations, observations, reconstruction)
    type(oris),        intent(inout) :: orientations
    real,              intent(in)    :: observations(:,:,:)
    real, allocatable, intent(out)   :: reconstruction(:,:,:)
    type(ctfparams) :: ctfparms
    type(fplane_type) :: fplane
    type(image) :: observation, observation_padded, restored
    type(ori) :: orientation
    type(parameters), target :: params
    type(reconstructor) :: gridder
    type(sp_project) :: project
    type(sym) :: c1sym
    integer :: i

    params%box        = BOX
    params%box_crop   = BOX
    params%box_croppd = OSMPL_PAD_FAC * BOX
    params%smpd_crop  = HALFSET_SMPD
    params%nstates    = 1
    params%numlen     = 1
    ! A class-volume field gives this CTF-free in-memory control no project dependency.
    params%oritype    = 'cls3D'

    call gridder%new_accumulator(params, project, expand=.true., wthreads=.false.)
    call gridder%set_ft(.true.)
    call restored%new([BOX,BOX,BOX], HALFSET_SMPD, wthreads=.false.)
    call observation%new([BOX,BOX,1], HALFSET_SMPD, wthreads=.false.)
    call observation_padded%new([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], &
        &HALFSET_SMPD, wthreads=.false.)
    call orientation%new(.false.)
    call c1sym%new('c1')
    call memoize_ft_maps([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], HALFSET_SMPD)
    ctfparms%smpd = HALFSET_SMPD
    ctfparms%ctfflag = CTFFLAG_NO

    do i = 1, orientations%get_noris()
        call orientations%get_ori(i, orientation)
        call observation%set_rmat(observations(:,:,i:i), .false.)
        call observation%pad(observation_padded, backgr=0., antialiasing=.false.)
        call observation_padded%fft()
        call observation_padded%gen_fplane4rec([0,BOX/2], HALFSET_SMPD, ctfparms, &
            &[0.,0.], fplane)
        ! Production gridding: padded Fourier plane -> KB insertion and sampling density.
        call gridder%insert_plane_oversamp(c1sym, orientation, fplane)
    enddo

    call gridder%compress_exp()
    ! Production convention: no box division, native-lattice KB deapodization
    ! (matches one-half gridding finalization after drop_legacy_box_division)
    call gridder%restore_final(restored)
    reconstruction = restored%get_rmat()

    if( allocated(fplane%cmplx_plane) ) deallocate(fplane%cmplx_plane)
    if( allocated(fplane%ctfsq_plane) ) deallocate(fplane%ctfsq_plane)
    if( allocated(fplane%transfer_plane) ) deallocate(fplane%transfer_plane)
    call forget_ft_maps()
    call c1sym%kill()
    call orientation%kill()
    call restored%kill()
    call observation_padded%kill()
    call observation%kill()
    call gridder%kill()
    call project%kill()
end subroutine reconstruct_half_conventionally

end module continuous_3D_pcg_refinement_halfset_gridding
