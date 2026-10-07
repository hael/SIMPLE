!@descr: common 3-D delivery of reconstructed flex_pca state maps.
!! Owns eo-FSC, filtering/masking, naming, publication and the project update.
module simple_flex_pca_state_delivery
use simple_core_module_api, only: cosmskhalfwidth, fdim, get_resolution, logfhandle, simple_exception, string
use simple_defs_flex,                only: FLEX_FSC_SIGNAL_THRESHOLD
use simple_builder,                  only: builder
use simple_parameters,               only: parameters
use simple_image,                    only: image
use simple_estimate_ssnr,            only: get_resolution
use simple_flex_pca_util,            only: flex_pca_write_state
use simple_flex_pca_project_gateway, only: prepare_project_fsc_lowpass_filters, publish_state_volume, write_out_segment
use simple_flex_pca_states_backend,  only: flex_state_maps, flex_state_delivery_policy
implicit none

public :: flex_state_delivery
private
#include "simple_local_flags.inc"

!> One delivery per state round. Per-state eo-FSC delivery: the even and odd maps of one state
!! carry identical kernel weights, so their FSC measures that state's own information content
!! (kernel regression borrows strength across the dataset, so a small state can resolve beyond its
!! bin count). The project FSC is only a fallback, where the backend's policy allows it.
type :: flex_state_delivery
    type(flex_state_delivery_policy) :: policy
    logical :: l_fuse = .false.
    integer :: box_rec = 0, filtsz = 0
    real    :: smpd_rec = 0., mskrad = 0.
    type(string) :: outvol_bak, outvol_even, outvol_odd, state_vol_fname
    real,    allocatable :: lowpass_filters(:,:), fsc_eo(:), filt_merged(:)
    logical, allocatable :: has_lowpass_filter(:)
    integer, allocatable :: lowpass_source_state(:)
  contains
    procedure :: new     => delivery_new
    procedure :: deliver => delivery_deliver
    procedure :: finish  => delivery_finish
    procedure, private :: write_view
    procedure :: kill    => delivery_kill
end type flex_state_delivery

contains

    subroutine delivery_new( self, params, policy, nstates, l_fuse, box_rec, smpd_rec, outvol_even, outvol_odd )
        class(flex_state_delivery),       intent(inout) :: self
        class(parameters),                intent(in)    :: params
        type(flex_state_delivery_policy), intent(in)    :: policy
        integer,                          intent(in)    :: nstates, box_rec
        logical,                          intent(in)    :: l_fuse
        real,                             intent(in)    :: smpd_rec
        type(string), optional,           intent(in)    :: outvol_even, outvol_odd
        call self%kill
        self%policy   = policy
        self%l_fuse   = l_fuse
        self%box_rec  = box_rec
        self%smpd_rec = smpd_rec
        self%filtsz   = max(1, fdim(box_rec) - 1)
        ! delivery mask at mskdiam, capped at the box edge; a box/2 mask lets solvent noise
        ! dominate the state FSC
        self%mskrad     = min(real(box_rec/2) - COSMSKHALFWIDTH - 1., 0.5*params%mskdiam/smpd_rec)
        self%outvol_bak = params%outvol
        if( present(outvol_even) ) self%outvol_even = outvol_even
        if( present(outvol_odd) )  self%outvol_odd  = outvol_odd
        if( l_fuse )then
            if( .not. (present(outvol_even) .and. present(outvol_odd)) ) &
                &THROW_HARD('flex halfset state delivery needs outvol_even/outvol_odd')
        endif
        allocate(self%fsc_eo(self%filtsz), self%filt_merged(self%filtsz), source=0.)
        if( self%policy%l_project_fsc_fallback )then
            call prepare_project_fsc_lowpass_filters(params, self%filtsz, nstates, self%lowpass_filters, &
                &self%has_lowpass_filter, self%lowpass_source_state)
        else
            allocate(self%lowpass_filters(self%filtsz,nstates), source=0.)
            allocate(self%has_lowpass_filter(nstates), source=.false.)
            allocate(self%lowpass_source_state(nstates), source=0)
        endif
    end subroutine delivery_new

    !> Deliver one state: with halves, the eo-FSC sets the 8th-order Butterworth low-pass at
    !! eo-FSC(0.143) for the combined map and both halves. A single-set state takes the
    !! project-FSC fallback where the policy allows it.
    subroutine delivery_deliver( self, params, build, state, maps )
        class(flex_state_delivery), intent(inout) :: self
        class(parameters),          intent(inout) :: params
        class(builder),             intent(inout) :: build
        integer,                    intent(in)    :: state
        type(flex_state_maps),      intent(inout) :: maps
        type(image) :: state_img, msk_e, msk_o
        real, allocatable :: res_arr(:)
        real    :: kc_lp, fsc05, fsc0143
        integer :: k_lp, iv
        logical :: l_eo_fsc
        if( maps%l_halves )then
            ! per-state eo FSC (on masked copies, background zeroed, where the policy masks)
            call msk_e%copy(maps%even)
            call msk_o%copy(maps%odd)
            if( self%policy%l_mask )then
                call msk_e%zero_background
                call msk_e%mask3D_soft(self%mskrad, backgr=0.)
                call msk_o%zero_background
                call msk_o%mask3D_soft(self%mskrad, backgr=0.)
            endif
            call msk_e%fft
            call msk_o%fft
            call msk_e%fsc(msk_o, self%fsc_eo)
            call msk_e%kill
            call msk_o%kill
            l_eo_fsc = any(self%fsc_eo > real(FLEX_FSC_SIGNAL_THRESHOLD))
            if( l_eo_fsc )then
                res_arr = maps%even%get_res()
                call get_resolution(self%fsc_eo, res_arr, fsc05, fsc0143)
                kc_lp = real(self%box_rec) * self%smpd_rec / fsc0143
                do k_lp = 1, self%filtsz
                    self%filt_merged(k_lp) = 1.0 / (1.0 + (real(k_lp)/max(kc_lp,1.0))**8)
                end do
                write(logfhandle,'(A,I3,A,F7.2,A,F7.2,A)') '>>> FLEX STATE'//trim(self%policy%tag)//' eo-FSC state=', &
                    &state,'  res(0.143)=',fsc0143,' A  res(0.5)=',fsc05,' A -- low-pass at the state eo-FSC(0.143) applied'
                deallocate(res_arr)
            else
                if( self%policy%l_project_fsc_fallback )then
                    write(logfhandle,'(A,I3,A)') '>>> FLEX STATE'//trim(self%policy%tag)//' eo-FSC state=',state, &
                        &'  unmeasurable (no shell above 0.143); project-FSC low-pass if the project has one, else unfiltered'
                else
                    write(logfhandle,'(A,I3,A)') '>>> FLEX STATE'//trim(self%policy%tag)//' eo-FSC state=',state, &
                        &'  unmeasurable (no shell above 0.143); delivered unfiltered'
                endif
            endif
            call flush(logfhandle)
            do iv = 1, 3
                select case(iv)
                    case(1); call state_img%copy(maps%combined); params%outvol = self%outvol_bak
                    case(2); call state_img%copy(maps%even);     params%outvol = self%outvol_even
                    case(3); call state_img%copy(maps%odd);      params%outvol = self%outvol_odd
                end select
                if( l_eo_fsc )then
                    call state_img%apply_filter(self%filt_merged)
                else if( self%policy%l_project_fsc_fallback .and. self%has_lowpass_filter(state) )then
                    call state_img%apply_filter(self%lowpass_filters(:,state))
                endif
                call self%write_view(params, build, state_img, state, iv == 1)
            end do
        else
            params%outvol = self%outvol_bak
            call state_img%copy(maps%combined)
            if( self%policy%l_project_fsc_fallback .and. self%has_lowpass_filter(state) )then
                call state_img%apply_filter(self%lowpass_filters(:,state))
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE applied project-FSC low-pass filter to state=',state, &
                    &' using_source_state=',self%lowpass_source_state(state)
            else
                write(logfhandle,'(A,I3,A)') '>>> FLEX STATE'//trim(self%policy%tag)//' state=', state, &
                    &'  single-set map delivered unfiltered'
            endif
            call self%write_view(params, build, state_img, state, .true.)
        endif
    end subroutine delivery_deliver

    !> Mask (per the policy), write under params%outvol's naming, and publish the combined view.
    !! Background removal + soft spherical mask: each state carries a different total kernel
    !! weight, so without both the states differ by a baseline offset and solvent noise rather
    !! than by conformation. Radius is capped at the broadest soft-maskable sphere in the box.
    subroutine write_view( self, params, build, state_img, state, l_publish )
        class(flex_state_delivery), intent(inout) :: self
        class(parameters),          intent(in)    :: params
        class(builder),             intent(inout) :: build
        type(image),                intent(inout) :: state_img
        integer,                    intent(in)    :: state
        logical,                    intent(in)    :: l_publish
        if( self%policy%l_mask )then
            call state_img%zero_background
            call state_img%mask3D_soft(self%mskrad, backgr=0.)
        endif
        call flex_pca_write_state(params%outvol, state_img, state, self%state_vol_fname)
        if( l_publish ) call publish_state_volume(build, self%state_vol_fname, state_img%get_smpd(), state, state_img%get_box())
        call state_img%kill
    end subroutine write_view

    !> After every state: the project's out segment, and params%outvol restored.
    subroutine delivery_finish( self, params, build )
        class(flex_state_delivery), intent(inout) :: self
        class(parameters),          intent(inout) :: params
        class(builder),             intent(inout) :: build
        call write_out_segment(build, params%projfile)
        params%outvol = self%outvol_bak
    end subroutine delivery_finish

    subroutine delivery_kill( self )
        class(flex_state_delivery), intent(inout) :: self
        call self%outvol_bak%kill
        call self%outvol_even%kill
        call self%outvol_odd%kill
        call self%state_vol_fname%kill
        if( allocated(self%lowpass_filters) )      deallocate(self%lowpass_filters)
        if( allocated(self%has_lowpass_filter) )   deallocate(self%has_lowpass_filter)
        if( allocated(self%lowpass_source_state) ) deallocate(self%lowpass_source_state)
        if( allocated(self%fsc_eo) )      deallocate(self%fsc_eo)
        if( allocated(self%filt_merged) ) deallocate(self%filt_merged)
        self%l_fuse = .false.; self%box_rec = 0; self%filtsz = 0; self%smpd_rec = 0.; self%mskrad = 0.
    end subroutine delivery_kill

end module simple_flex_pca_state_delivery
