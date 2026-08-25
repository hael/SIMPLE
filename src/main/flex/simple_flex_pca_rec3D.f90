!@descr: kernel-weighted state reconstruction for flex_pca
module simple_flex_pca_rec3D
use simple_core_module_api
use simple_builder,                only: builder
use simple_sp_project,             only: sp_project
use simple_gridding,               only: prep3D_inv_kbenvelope4mul
use simple_image,                  only: image
use simple_matcher_3Drec,          only: init_rec, prep_imgs4rec, cleanup_rec_buffers
use simple_matcher_ptcl_io,        only: discrete_read_imgbatch, discrete_read_imgbatch_source, prepimgbatch
use simple_parameters,             only: parameters
use simple_flex_pca_distr, only: flex_pca_is_master, flex_pca_is_worker, flex_pca_nparts, &
    &flex_pca_run_stage, PCA_STAGE_STATES
use simple_flex_pca_parts, only: write_state_weights_round
use simple_flex_reconstructor_latent_ops, only: insert_planes_oversamp_multi_scaled_batch
use simple_reconstructor,          only: reconstructor
implicit none
private
#include "simple_local_flags.inc"

public :: reconstruct_flex_weighted_states
public :: flex_rec_box, flex_rec_smpd

contains

    !> Direct manifold pre-image reconstruction with kernel weights.
    !! Medoid-selected manifold descriptors are converted to soft particle
    !! weights upstream; this routine reconstructs each state from weighted
    !! contributions of all particles rather than from a global residual basis.
    !> With outvol_even/outvol_odd present this performs the COMBINED, EVEN and ODD reconstructions in
    !! ONE pass instead of three: plane insertion is linear in the weights and every nonlinear
    !! finalisation runs after compress_exp, so the combined accumulator is exactly even+odd. Each
    !! particle is inserted once, into its own halfset. Distributed, it collapses three qsys rounds into one.
    subroutine reconstruct_flex_weighted_states( params, build, pinds, state_weights, nstates, fsc_projfile, &
        &floor_rho, outvol_even, outvol_odd, split_eo )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nstates
        real,              intent(in)    :: state_weights(:,:)
        type(string), optional, intent(in) :: fsc_projfile
        ! RELION-style shellwise rho floor before the divide. OFF by default so the diffusion-map
        ! callers keep their previous behaviour exactly; flex_pca opts in. See the note at the
        ! floor call below for why the kernel-weighted path needs it.
        logical,      optional, intent(in) :: floor_rho
        type(string), optional, intent(in) :: outvol_even, outvol_odd
        !! Worker-side entry point for the halfset split: the master decides by supplying the two output
        !! names, but a worker has none (it writes part files), so it is told through the round-weights
        !! table. Both routes must set l_fuse identically or master and workers disagree on the halves.
        logical,      optional, intent(in) :: split_eo
        logical :: l_floor_rho, l_reduced, l_fuse
        type(reconstructor), allocatable :: recs_o(:), recs_c(:), cur(:)
        type(string) :: outvol_bak, state_vol_fname
        integer :: eo_i, iview, nview
        type(reconstructor), allocatable :: state_recs(:)
        type(fplane_type), allocatable :: fpls(:)
        type(image) :: gridcorr_img, state_img
        type(ori) :: orientation
        ! batch buffers for insert_planes_oversamp_multi_scaled_batch; bvalid_c/bvalid_o are disjoint
        type(ori),            allocatable :: borientations(:)
        real(dp),             allocatable :: bscales(:,:)
        logical,              allocatable :: bvalid_c(:), bvalid_o(:)
        real(dp), allocatable :: scales(:)
        real, allocatable :: lowpass_filters(:,:)
        logical, allocatable :: has_lowpass_filter(:)
        integer, allocatable :: lowpass_source_state(:)
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, state, box_rec
        real    :: smpd_rec, smpd_crop_bak
        l_floor_rho = .false.
        if( present(floor_rho) ) l_floor_rho = floor_rho
        l_fuse = present(outvol_even) .and. present(outvol_odd)
        if( present(split_eo) ) l_fuse = l_fuse .or. split_eo
        if( size(pinds)<1 .or. nstates<1 ) THROW_HARD('invalid flex weighted state reconstruction dimensions')
        if( any(shape(state_weights)/=[size(pinds),nstates]) ) THROW_HARD('flex weighted state table mismatch')
        box_rec  = flex_rec_box(params)
        smpd_rec = flex_rec_smpd(params)
        if( box_rec /= params%box_crop )then
            write(logfhandle,'(A,I0,A,F6.3,A,I0,A,F6.3,A)') '>>> FLEX STATE RECONSTRUCTION decoupled box: rec box=',box_rec, &
                &' smpd=',smpd_rec,' A (covariance box=',params%box_crop,' smpd=',params%smpd_crop,' A)'
        endif
        allocate(state_recs(nstates),scales(nstates))
        call prepare_project_fsc_lowpass_filters(params,build,nstates,lowpass_filters,has_lowpass_filter,lowpass_source_state, &
            &fsc_projfile)
        do state=1,nstates
            call init_state_reconstructor(params,build,state_recs(state))
        end do
        if( l_fuse )then
            allocate(recs_o(nstates))
            do state=1,nstates
                call init_state_reconstructor(params,build,recs_o(state))
            end do
        endif
        call init_rec(params,build,MAXIMGBATCHSZ,fpls)
        call prepimgbatch(params,build,MAXIMGBATCHSZ)
        ! prep_imgs4rec builds the Fourier planes at params%smpd_crop, which fixes the physical
        ! frequency each plane sample carries and hence the CTF evaluated there. These state maps
        ! are reconstructed on the box_rec/smpd_rec lattice, so the planes must be built at
        ! smpd_rec or the CTF lands on the wrong frequencies. Point smpd_crop at the reconstruction
        ! sampling for the batch loop and restore it afterwards. Safe because prep_imgs4rec is the
        ! only consumer of smpd_crop in between: init_state_reconstructor sizes from
        ! flex_rec_box/flex_rec_smpd and init_rec uses boxpd/smpd. A no-op by default, since
        ! box_rec defaults to box_crop and therefore smpd_rec == smpd_crop.
        ! DISTRIBUTED: the master ships this round's weight table, fans the particle range out and sums
        ! the compressed partial reconstructions. Every nonlinear finalisation below -- the density
        ! floor, sampl_dens_correct, zero_background, mask3D_soft -- then runs ONCE on the global sums,
        ! which is what makes the reduction equivalent to the in-process accumulation.
        l_reduced = .false.
        if( flex_pca_is_master() )then
            call write_state_weights_round(pinds, state_weights, size(pinds), nstates, l_fuse)
            call flex_pca_run_stage(PCA_STAGE_STATES, 'state reconstruction')
            block
                type(reconstructor) :: rec_read
                type(string) :: pf
                integer :: ipart
                integer(timer_int_kind) :: t_red
                t_red = tic()
                call init_state_reconstructor(params,build,rec_read)
                do ipart = 1, flex_pca_nparts()
                    do state = 1, nstates
                        ! on a split round each part carries BOTH halfsets; reduce each into its own
                        ! accumulator so combined = even + odd below is a sum of two populated halves
                        do eo_i = 0, merge(1, 0, l_fuse)
                            pf = flex_state_part_fbody(params, ipart, state, eo_i)
                            if( .not. file_exists(pf//MRC_EXT) ) THROW_HARD('missing states part: '//pf%to_char())
                            call rec_read%read(pf//MRC_EXT)
                            call rec_read%read_rho(string('rho_')//pf//MRC_EXT)
                            if( eo_i == 1 )then
                                call recs_o(state)%sum_reduce(rec_read)
                            else
                                call state_recs(state)%sum_reduce(rec_read)
                            endif
                            call del_file(pf//MRC_EXT)
                            call del_file(string('rho_')//pf//MRC_EXT)
                            call pf%kill
                        end do
                    end do
                end do
                call rec_read%dealloc_rho; call rec_read%kill
                write(logfhandle,'(A,I0,A,F8.1)') '>>> FLEX_PCA reduced states parts=', &
                    &flex_pca_nparts(),' seconds=',toc(t_red)
                call flush(logfhandle)
            end block
            l_reduced = .true.
            goto 300
        endif
        smpd_crop_bak    = params%smpd_crop
        params%smpd_crop = smpd_rec
        allocate(borientations(MAXIMGBATCHSZ), bscales(nstates,MAXIMGBATCHSZ))
        allocate(bvalid_c(MAXIMGBATCHSZ), bvalid_o(MAXIMGBATCHSZ))
        do ibatch=1,size(pinds),MAXIMGBATCHSZ
            batchlims=[ibatch,min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
            batchsz=batchlims(2)-batchlims(1)+1
            if( params%l_ptcl_src_den )then
                call discrete_read_imgbatch_source(params,build,'den',batchsz, &
                    &pinds(batchlims(1):batchlims(2)),[1,batchsz],build%imgbatch(:batchsz))
            else
                call discrete_read_imgbatch(params,build,size(pinds),pinds,batchlims)
            endif
            call prep_imgs4rec(params,build,batchsz,build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)),fpls(:batchsz))
            ! gather serially (get_ori/get_eo are not guaranteed thread-safe), then one parallel
            ! region per target; under l_fuse each halfset gets its own mask and call
            do i=1,batchsz
                iptcl=pinds(batchlims(1)+i-1)
                call build%spproj_field%get_ori(iptcl,orientation)
                bscales(:,i) = real(state_weights(batchlims(1)+i-1,:),dp)
                bvalid_c(i)  = .not. orientation%isstatezero()
                bvalid_o(i)  = .false.
                if( bvalid_c(i) .and. l_fuse )then
                    ! one insertion per particle, into the halfset it belongs to; the union is the
                    ! combined map and each half is its own deliverable
                    eo_i = build%spproj_field%get_eo(iptcl)
                    if( eo_i == 1 )then
                        bvalid_o(i) = .true.
                        bvalid_c(i) = .false.
                    endif
                endif
                call borientations(i)%copy(orientation)
            end do
            if( l_fuse .and. any(bvalid_o(:batchsz)) )then
                call insert_planes_oversamp_multi_scaled_batch(recs_o, build%pgrpsyms, &
                    &borientations(:batchsz), fpls(:batchsz), bscales(:,:batchsz), &
                    &bscales(:,:batchsz), bvalid_o(:batchsz), batchsz)
            endif
            if( any(bvalid_c(:batchsz)) )then
                call insert_planes_oversamp_multi_scaled_batch(state_recs, build%pgrpsyms, &
                    &borientations(:batchsz), fpls(:batchsz), bscales(:,:batchsz), &
                    &bscales(:,:batchsz), bvalid_c(:batchsz), batchsz)
            endif
        end do
        do i = 1, MAXIMGBATCHSZ
            call borientations(i)%kill
        end do
        deallocate(borientations, bscales, bvalid_c, bvalid_o)
        call orientation%kill
        params%smpd_crop = smpd_crop_bak
        call cleanup_rec_buffers(build,fpls)
        if( flex_pca_is_worker() )then
            block
                type(string) :: pf
                do state=1,nstates
                    call state_recs(state)%compress_exp
                    pf = flex_state_part_fbody(params, params%part, state, 0)
                    call state_recs(state)%write(pf//MRC_EXT, del_if_exists=.true.)
                    call state_recs(state)%write_rho(string('rho_')//pf//MRC_EXT)
                    call pf%kill
                    if( l_fuse )then
                        call recs_o(state)%compress_exp
                        pf = flex_state_part_fbody(params, params%part, state, 1)
                        call recs_o(state)%write(pf//MRC_EXT, del_if_exists=.true.)
                        call recs_o(state)%write_rho(string('rho_')//pf//MRC_EXT)
                        call pf%kill
                    endif
                end do
            end block
            do state=1,nstates
                call state_recs(state)%dealloc_rho; call state_recs(state)%kill
                if( l_fuse )then
                    call recs_o(state)%dealloc_rho; call recs_o(state)%kill
                endif
            end do
            deallocate(state_recs,scales,lowpass_filters,has_lowpass_filter,lowpass_source_state)
            if( l_fuse ) deallocate(recs_o)
            return
        endif
300     continue
        gridcorr_img=prep3D_inv_kbenvelope4mul([box_rec,box_rec,box_rec], smpd_rec)
        outvol_bak = params%outvol
        nview = 1
        if( l_fuse )then
            ! Build the combined accumulator BEFORE any finalisation, because finalisation is
            ! destructive (compress_exp, ifft, masking). combined = even + odd exactly.
            do state=1,nstates
                if( .not. l_reduced )then
                    call state_recs(state)%compress_exp
                    call recs_o(state)%compress_exp
                endif
            end do
            allocate(recs_c(nstates))
            do state=1,nstates
                call init_state_reconstructor(params,build,recs_c(state))
                call recs_c(state)%sum_reduce(state_recs(state))
                call recs_c(state)%sum_reduce(recs_o(state))
            end do
            l_reduced = .true.
            nview     = 3
        endif
        do iview = 1, nview
            if( l_fuse )then
                select case(iview)
                    case(1); call move_alloc(recs_c,     cur); params%outvol = outvol_bak
                    case(2); call move_alloc(state_recs, cur); params%outvol = outvol_even
                    case(3); call move_alloc(recs_o,     cur); params%outvol = outvol_odd
                end select
            else
                call move_alloc(state_recs, cur)
            endif
            do state=1,nstates
                if( .not. l_reduced ) call cur(state)%compress_exp
                ! Kernel weights live in [0,1] with most near zero, so rho is small and highly
                ! variable here and an unfloored divide amplifies noise wherever occupancy is low.
                ! Opt-in: the platform reconstructor is unchanged, and so is the diffusion-map path.
                if( l_floor_rho ) call cur(state)%floor_rho_shellwise
                call cur(state)%sampl_dens_correct
                call cur(state)%ifft
                call cur(state)%mul(gridcorr_img)
                call state_img%copy(cur(state))
                if( has_lowpass_filter(state) )then
                    call state_img%apply_filter(lowpass_filters(:,state))
                    write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE applied project-FSC low-pass filter to state=',state, &
                        &' using_source_state=',lowpass_source_state(state)
                endif
                ! Background removal + soft spherical mask. Without both the states look like one smeared map:
                ! each carries a different total kernel weight, so the difference between two states is
                ! dominated by a constant baseline offset rather than by the conformational change
                ! (zero_background reads the level off the box faces, shape-preserving), and solvent noise
                ! dominates any unmasked comparison. Radius convention is the platform's (see reconstructor_eo).
                call state_img%zero_background
                call state_img%mask3D_soft(real(box_rec/2) - COSMSKHALFWIDTH - 1., backgr=0.)
                call write_state(params, state_img, state, state_vol_fname)
                ! Add merged volume only to project
                if( iview == 1 ) then
                    call build%spproj%add_vol2os_out(state_vol_fname, state_img%get_smpd(), state, 'vol_flex',&
                        &box=state_img%get_box())
                endif
                ! clean up
                call state_img%kill
                call cur(state)%dealloc_rho
                call cur(state)%kill
            end do
        end do
        call build%spproj%write_segment_inside('out', params%projfile)
        params%outvol = outvol_bak
        ! clean up
        call gridcorr_img%kill
            call state_vol_fname%kill
        deallocate(scales,lowpass_filters,has_lowpass_filter,lowpass_source_state)
    end subroutine reconstruct_flex_weighted_states

    !> Part-file body for one worker's partial reconstruction of one state. eo selects the halfset
    !! accumulator on a split round: 0 is the even/single accumulator, 1 the odd one. A split round
    !! writes both per state, so the two must not collide on disk.
    function flex_state_part_fbody( params, part, state, eo ) result( fbody )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: part, state, eo
        type(string) :: fbody
        fbody = string('flex_pca_statepart')//int2str_pad(part,max(1,params%numlen))// &
            &'_'//int2str_pad(state,2)
        if( eo == 1 ) fbody = fbody//'_o'
    end function flex_state_part_fbody

    subroutine prepare_project_fsc_lowpass_filters( params, build, nstates, lowpass_filters, has_filter, source_state, &
        &fsc_projfile )
        class(parameters), intent(in) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in) :: nstates
        real, allocatable, intent(out) :: lowpass_filters(:,:)
        logical, allocatable, intent(out) :: has_filter(:)
        integer, allocatable, intent(out) :: source_state(:)
        type(string), optional, intent(in) :: fsc_projfile
        type(sp_project) :: spproj
        type(string) :: fsc_fname, imgkind_here, proj_for_fsc
        real, allocatable :: fsc(:)
        integer :: filtsz, state, fsc_box, i, state1_fsc_count
        logical :: out_loaded
        ! sized to the DELIVERED map, which is box_rec (== box_crop unless decoupled)
        filtsz=fdim(flex_rec_box(params))-1
        allocate(lowpass_filters(filtsz,nstates),has_filter(nstates),source_state(nstates))
        lowpass_filters=0.
        has_filter=.false.
        source_state=0
        if( filtsz<1 ) return
        proj_for_fsc=params%projfile
        if( present(fsc_projfile) )then
            if( len_trim(fsc_projfile%to_char())>0 ) proj_for_fsc=fsc_projfile
        endif
        if( .not.file_exists(proj_for_fsc) )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass filtering skipped: projfile not found'
            call proj_for_fsc%kill
            return
        endif
        call spproj%read_segment('out',proj_for_fsc)
        out_loaded=spproj%os_out%get_noris()>0
        if( .not.out_loaded )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass filtering skipped: empty out segment'
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        state1_fsc_count=0
        do i=1,spproj%os_out%get_noris()
            if( .not.spproj%os_out%isthere(i,'imgkind') ) cycle
            call spproj%os_out%getter(i,'imgkind',imgkind_here)
            if( imgkind_here%to_char()/='fsc' ) cycle
            if( spproj%os_out%get_state(i)==1 ) state1_fsc_count=state1_fsc_count+1
        end do
        if( state1_fsc_count/=1 )then
            write(logfhandle,'(A,I0)') '>>> FLEX PRE-IMAGE project-FSC low-pass filtering skipped: expected exactly one state=1 FSC in out segment, found=', &
                &state1_fsc_count
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        call spproj%get_fsc(1,fsc_fname,fsc_box)
        if( .not.file_exists(fsc_fname) )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass filtering skipped: state=1 FSC file missing'
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        fsc=file2rarr(fsc_fname)
        if( size(fsc)/=filtsz )then
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE project-FSC low-pass filtering skipped: state=1 FSC size mismatch; fsc_nyq=', &
                &size(fsc),' model_nyq=',filtsz
            deallocate(fsc)
            call fsc_fname%kill
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        do state=1,nstates
            call fsc2optlp_sub(filtsz,fsc,lowpass_filters(:,state),merged=.false.)
            has_filter(state)=any(lowpass_filters(:,state)>0.)
            source_state(state)=1
        end do
        deallocate(fsc)
        call fsc_fname%kill
        call imgkind_here%kill
        call proj_for_fsc%kill
        call spproj%kill
    end subroutine prepare_project_fsc_lowpass_filters

    !> Box/sampling of the DELIVERED state maps.  Decoupled from box_crop because the covariance and
    !! the embedding are low-frequency objects whose column FSC dies well short of Nyquist, whereas
    !! the state maps are ordinary backprojections of the same particles and carry signal beyond it.
    !! With box_rec==box_crop this is a no-op.  The embedding never sees the extra band, so the added
    !! shells cannot be selected on themselves and this cannot manufacture structure from noise.
    pure integer function flex_rec_box( params ) result( box_rec )
        class(parameters), intent(in) :: params
        box_rec = params%box_crop
        if( params%box_rec >= 1 ) box_rec = params%box_rec
    end function flex_rec_box

    pure real function flex_rec_smpd( params ) result( smpd_rec )
        class(parameters), intent(in) :: params
        smpd_rec = params%smpd_crop
        if( params%box_rec >= 1 .and. params%smpd_rec > 0. ) smpd_rec = params%smpd_rec
    end function flex_rec_smpd

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

    subroutine write_state( params, img, state, vol_fname )
        class(parameters), intent(in)    :: params
        class(image),      intent(inout) :: img
        integer,           intent(in)    :: state
        class(string),     intent(inout) :: vol_fname
        type(string) :: prefix, ext
        character(len=:), allocatable :: stem
        character(len=3) :: tag
        if( state==1 )then
            vol_fname = params%outvol
        else
            ext=fname2ext(params%outvol)
            prefix=get_fbody(params%outvol,ext)
            stem=prefix%to_char()
            if( len_trim(stem)>4 )then
                if( stem(len_trim(stem)-3:len_trim(stem))=='_001' ) stem=stem(:len_trim(stem)-4)
            endif
            prefix=string(stem)
            write(tag,'(I3.3)') state
            vol_fname = prefix//'_'//tag//MRC_EXT
        endif
        call img%write(vol_fname,del_if_exists=.true.)
        write(logfhandle,'(A,I0,A,A)') '>>> FLEX DIFFMAP NYSTROM PRE-IMAGE ',state,': ',vol_fname%to_char()
        call prefix%kill
        call ext%kill
    end subroutine write_state

end module simple_flex_pca_rec3D
