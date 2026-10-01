!@descr: flex_pca: the gridding state-reconstruction backend -- kernel-weighted backprojection of every state in one pass
module simple_flex_pca_states_gridding
use simple_core_module_api
use simple_builder,         only: builder
use simple_parameters,      only: parameters
use simple_image,           only: image
use simple_reconstructor,   only: reconstructor
use simple_gridding,        only: prep3D_inv_kbenvelope4mul
use simple_matcher_3Drec,   only: init_rec, prep_imgs4rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: insert_planes_oversamp_multi_scaled_batch
use simple_flex_pca_rounds,         only: flex_pca_rounds
use simple_flex_pca_run_types,      only: flex_run_settings
use simple_flex_pca_state_parts,    only: flex_state_part_fbody, flex_pca_rho_part_name
use simple_flex_pca_states_backend, only: flex_states_backend, flex_state_maps, flex_state_delivery_policy, &
    &flex_rec_box, flex_rec_smpd
implicit none

public :: flex_states_gridding
private
#include "simple_local_flags.inc"

!> Each state is a weighted backprojection of all particles. With halves the combined, even and
!! odd maps come from one pass: insertion is linear in the weights and every nonlinear
!! finalisation runs after compress_exp, so combined = even + odd. Each particle is inserted
!! once, into its own halfset.
type, extends(flex_states_backend) :: flex_states_gridding
    type(reconstructor), allocatable :: state_recs(:)  !< even (or single-set) accumulators
    type(reconstructor), allocatable :: recs_o(:)      !< odd accumulators (l_fuse)
    type(image) :: gridcorr_img
    logical     :: l_reduced = .false.   !< accumulators hold compressed global sums (the parts' reduce)
  contains
    procedure :: begin                          => gridding_begin
    procedure :: accumulate_local_or_write_part => gridding_accumulate
    procedure :: fold_parts                     => gridding_fold_parts
    procedure :: finalize_maps                  => gridding_finalize_maps
    procedure :: delivery_policy                => gridding_delivery_policy
    procedure :: kill                           => gridding_kill
end type flex_states_gridding

contains

    subroutine gridding_begin( self, params, build, rounds, pinds, state_weights, nstates, l_fuse, l_floor_rho, box_rec, smpd_rec )
        class(flex_states_gridding), intent(inout) :: self
        class(parameters),           intent(inout) :: params
        class(builder),              intent(inout) :: build
        class(flex_pca_rounds),      intent(inout) :: rounds
        integer,                     intent(in)    :: pinds(:), nstates, box_rec
        real,                        intent(in)    :: state_weights(:,:), smpd_rec
        logical,                     intent(in)    :: l_fuse, l_floor_rho
        integer :: state
        call self%kill
        call self%set_selection(pinds, state_weights, nstates, l_fuse, l_floor_rho, box_rec, smpd_rec)
        ! the delivered state maps are always low-passed at their own eo-FSC(0.143) resolution:
        ! a poorly determined state must look poorly determined
        write(logfhandle,'(A)') '>>> FLEX_PCA state maps delivered under a per-state low-pass at each state''s own &
            &eo-FSC(0.143) resolution'
        allocate(self%state_recs(nstates))
        do state=1,nstates
            call init_state_reconstructor(params,build,self%state_recs(state))
        end do
        if( l_fuse )then
            allocate(self%recs_o(nstates))
            do state=1,nstates
                call init_state_reconstructor(params,build,self%recs_o(state))
            end do
        endif
        self%gridcorr_img = prep3D_inv_kbenvelope4mul([box_rec,box_rec,box_rec], smpd_rec)
        self%l_reduced    = .false.
    end subroutine gridding_begin

    !> One pass over the selection into the expanded accumulators; a worker then compresses and
    !! writes its per-(state, half) parts.
    subroutine gridding_accumulate( self, params, build, rounds )
        class(flex_states_gridding), intent(inout) :: self
        class(parameters),           intent(inout) :: params
        class(builder),              intent(inout) :: build
        class(flex_pca_rounds),      intent(inout) :: rounds
        type(fplane_type), allocatable :: fpls(:)
        type(ori) :: orientation
        ! batch buffers for insert_planes_oversamp_multi_scaled_batch; bvalid_c/bvalid_o are disjoint
        type(ori),            allocatable :: borientations(:)
        real(dp),             allocatable :: bscales(:,:)
        logical,              allocatable :: bvalid_c(:), bvalid_o(:)
        type(string) :: pf
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, state, eo_i
        real    :: smpd_crop_bak
        integer(timer_int_kind) :: t_ins, t_sec
        real(timer_int_kind) :: sec_read, sec_prep, sec_ins
        call init_rec(params,build,MAXIMGBATCHSZ,fpls)
        call prepimgbatch(params,build,MAXIMGBATCHSZ)
        ! prep_imgs4rec builds planes at params%smpd_crop, which fixes the CTF frequencies; the maps
        ! live on the box_rec/smpd_rec lattice, so smpd_crop is pointed at smpd_rec for the batch loop
        ! and restored afterwards (prep_imgs4rec is the only consumer in between; no-op unless decoupled).
        smpd_crop_bak    = params%smpd_crop
        params%smpd_crop = self%smpd_rec
        allocate(borientations(MAXIMGBATCHSZ), bscales(self%nstates,MAXIMGBATCHSZ))
        allocate(bvalid_c(MAXIMGBATCHSZ), bvalid_o(MAXIMGBATCHSZ))
        t_ins = tic()
        sec_read = 0.; sec_prep = 0.; sec_ins = 0.
        do ibatch=1,size(self%pinds),MAXIMGBATCHSZ
            batchlims=[ibatch,min(size(self%pinds),ibatch+MAXIMGBATCHSZ-1)]
            batchsz=batchlims(2)-batchlims(1)+1
            t_sec = tic()
            call discrete_read_imgbatch(params,build,size(self%pinds),self%pinds,batchlims)
            sec_read = sec_read + toc(t_sec)
            t_sec = tic()
            call prep_imgs4rec(params,build,batchsz,build%imgbatch(:batchsz), &
                &self%pinds(batchlims(1):batchlims(2)),fpls(:batchsz))
            sec_prep = sec_prep + toc(t_sec)
            ! gather serially (get_ori/get_eo are not guaranteed thread-safe), then one parallel
            ! region per target; under l_fuse each halfset gets its own mask and call
            do i=1,batchsz
                iptcl=self%pinds(batchlims(1)+i-1)
                call build%spproj_field%get_ori(iptcl,orientation)
                bscales(:,i) = real(self%state_weights(batchlims(1)+i-1,:),dp)
                bvalid_c(i)  = .not. orientation%isstatezero()
                bvalid_o(i)  = .false.
                if( bvalid_c(i) .and. self%l_fuse )then
                    ! one insertion per particle, into its own halfset
                    eo_i = build%spproj_field%get_eo(iptcl)
                    if( eo_i == 1 )then
                        bvalid_o(i) = .true.
                        bvalid_c(i) = .false.
                    endif
                endif
                call borientations(i)%copy(orientation)
            end do
            t_sec = tic()
            if( self%l_fuse .and. any(bvalid_o(:batchsz)) )then
                call insert_planes_oversamp_multi_scaled_batch(self%recs_o, build%pgrpsyms, &
                    &borientations(:batchsz), fpls(:batchsz), bscales(:,:batchsz), &
                    &bscales(:,:batchsz), bvalid_o(:batchsz), batchsz)
            endif
            if( any(bvalid_c(:batchsz)) )then
                call insert_planes_oversamp_multi_scaled_batch(self%state_recs, build%pgrpsyms, &
                    &borientations(:batchsz), fpls(:batchsz), bscales(:,:batchsz), &
                    &bscales(:,:batchsz), bvalid_c(:batchsz), batchsz)
            endif
            sec_ins = sec_ins + toc(t_sec)
        end do
        write(logfhandle,'(A,F8.1)') '>>> FLEX_PCA STATEREC read+prep+insert seconds=', toc(t_ins)
        write(logfhandle,'(A,F7.1,A,F7.1,A,F7.1)') '>>> FLEX_PCA STATEREC SPLIT (seconds): read=', &
            &sec_read,'  prep=',sec_prep,'  insert=',sec_ins
        call flush(logfhandle)
        do i = 1, MAXIMGBATCHSZ
            call borientations(i)%kill
        end do
        deallocate(borientations, bscales, bvalid_c, bvalid_o)
        call orientation%kill
        params%smpd_crop = smpd_crop_bak
        call cleanup_rec_buffers(build,fpls)
        if( rounds%is_worker() )then
            do state=1,self%nstates
                call self%state_recs(state)%compress_exp
                pf = flex_state_part_fbody(params, params%part, state, 0)
                call self%state_recs(state)%write(pf//MRC_EXT, del_if_exists=.true.)
                call self%state_recs(state)%write_rho(flex_pca_rho_part_name(pf))
                call pf%kill
                if( self%l_fuse )then
                    call self%recs_o(state)%compress_exp
                    pf = flex_state_part_fbody(params, params%part, state, 1)
                    call self%recs_o(state)%write(pf//MRC_EXT, del_if_exists=.true.)
                    call self%recs_o(state)%write_rho(flex_pca_rho_part_name(pf))
                    call pf%kill
                endif
            end do
        endif
    end subroutine gridding_accumulate

    !> Distributed master: sum the compressed partials of every part into the accumulators; every
    !! nonlinear finalisation then runs once on the global sums.
    subroutine gridding_fold_parts( self, params, build, rounds )
        class(flex_states_gridding), intent(inout) :: self
        class(parameters),           intent(inout) :: params
        class(builder),              intent(inout) :: build
        class(flex_pca_rounds),      intent(inout) :: rounds
        type(reconstructor) :: rec_read
        type(string) :: pf
        integer :: ipart, state, eo_i
        integer(timer_int_kind) :: t_red
        t_red = tic()
        call init_state_reconstructor(params,build,rec_read)
        do ipart = 1, rounds%nparts()
            do state = 1, self%nstates
                ! on a split round each part carries both halfsets; reduce each into its own
                ! accumulator so combined = even + odd below sums two populated halves
                do eo_i = 0, merge(1, 0, self%l_fuse)
                    pf = flex_state_part_fbody(params, ipart, state, eo_i)
                    if( .not. file_exists(pf//MRC_EXT) ) THROW_HARD('missing states part: '//pf%to_char())
                    call rec_read%read(pf//MRC_EXT)
                    call rec_read%read_rho(flex_pca_rho_part_name(pf))
                    if( eo_i == 1 )then
                        call self%recs_o(state)%sum_reduce(rec_read)
                    else
                        call self%state_recs(state)%sum_reduce(rec_read)
                    endif
                    call del_file(pf//MRC_EXT)
                    call del_file(flex_pca_rho_part_name(pf))
                    call pf%kill
                end do
            end do
        end do
        call rec_read%dealloc_rho; call rec_read%kill
        write(logfhandle,'(A,I0,A,F8.1)') '>>> FLEX_PCA reduced states parts=', &
            &rounds%nparts(),' seconds=',toc(t_red)
        call flush(logfhandle)
        self%l_reduced = .true.
    end subroutine gridding_fold_parts

    !> The views of one state (destructive on its accumulators): with halves the combined map is
    !! the finalisation of the summed even+odd accumulators.
    subroutine gridding_finalize_maps( self, params, build, rounds, state, maps )
        class(flex_states_gridding), intent(inout) :: self
        class(parameters),           intent(inout) :: params
        class(builder),              intent(inout) :: build
        class(flex_pca_rounds),      intent(inout) :: rounds
        integer,                     intent(in)    :: state
        type(flex_state_maps),       intent(inout) :: maps
        type(reconstructor) :: rec_c
        call maps%kill
        if( .not. self%l_reduced )then
            call self%state_recs(state)%compress_exp
            if( self%l_fuse ) call self%recs_o(state)%compress_exp
        endif
        if( self%l_fuse )then
            call init_state_reconstructor(params,build,rec_c)
            call rec_c%sum_reduce(self%state_recs(state))
            call rec_c%sum_reduce(self%recs_o(state))
            call finalize_state_rec(self%state_recs(state), self%gridcorr_img, self%l_floor_rho, maps%even)
            call finalize_state_rec(self%recs_o(state),     self%gridcorr_img, self%l_floor_rho, maps%odd)
            call finalize_state_rec(rec_c,                  self%gridcorr_img, self%l_floor_rho, maps%combined)
            call rec_c%dealloc_rho; call rec_c%kill
            call self%recs_o(state)%dealloc_rho; call self%recs_o(state)%kill
            maps%l_halves = .true.
        else
            ! kernel weights are mostly near zero, so rho is small and an unfloored divide
            ! amplifies noise wherever occupancy is low (the opt-in floor inside finalize_state_rec)
            call finalize_state_rec(self%state_recs(state), self%gridcorr_img, self%l_floor_rho, maps%combined)
            maps%l_halves = .false.
        endif
        call self%state_recs(state)%dealloc_rho; call self%state_recs(state)%kill
    end subroutine gridding_finalize_maps

    !> The gridding delivery: low-pass at the state's own eo-FSC(0.143), background removal + soft
    !! mask on every delivered view, the project FSC low-pass as the fallback. The run settings'
    !! eofilt/filt switches are the PCG backend's; here the values are fixed (user decision 2026-09-16).
    function gridding_delivery_policy( self, cfg ) result( policy )
        class(flex_states_gridding), intent(in) :: self
        type(flex_run_settings),     intent(in) :: cfg
        type(flex_state_delivery_policy) :: policy
        policy%l_state_eofilt = .false.
        policy%l_state_filt   = .true.
        policy%l_mask         = .true.
        policy%l_project_fsc_fallback = .true.
        policy%tag = ''
    end function gridding_delivery_policy

    subroutine gridding_kill( self )
        class(flex_states_gridding), intent(inout) :: self
        integer :: state
        if( allocated(self%state_recs) )then
            do state = 1, size(self%state_recs)
                call self%state_recs(state)%dealloc_rho; call self%state_recs(state)%kill
            end do
            deallocate(self%state_recs)
        endif
        if( allocated(self%recs_o) )then
            do state = 1, size(self%recs_o)
                call self%recs_o(state)%dealloc_rho; call self%recs_o(state)%kill
            end do
            deallocate(self%recs_o)
        endif
        call self%gridcorr_img%kill
        self%l_reduced = .false.
        call self%kill_selection
    end subroutine gridding_kill

    !> One state view: rho floor (opt-in), density correction, ifft, gridding correction -> image.
    !! Destructive on the reconstructor's Fourier state.
    subroutine finalize_state_rec( rec, gridcorr_img, l_floor_rho, img )
        type(reconstructor), intent(inout) :: rec
        type(image),         intent(in)    :: gridcorr_img
        logical,             intent(in)    :: l_floor_rho
        type(image),         intent(inout) :: img
        if( l_floor_rho ) call rec%floor_rho_shellwise
        call rec%sampl_dens_correct
        call rec%ifft
        call rec%mul(gridcorr_img)
        call img%copy(rec)
    end subroutine finalize_state_rec

    subroutine init_state_reconstructor( params, build, state_rec )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        type(reconstructor), intent(inout) :: state_rec
        integer :: box_rec
        box_rec = flex_rec_box(params)
        call state_rec%new([box_rec,box_rec,box_rec],flex_rec_smpd(params))
        call state_rec%alloc_rho(params,build%spproj,expand=.true.)
        call state_rec%reset
        call state_rec%reset_exp
    end subroutine init_state_reconstructor

end module simple_flex_pca_states_gridding
