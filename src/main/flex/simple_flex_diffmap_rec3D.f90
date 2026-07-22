!@descr: projection-aware residual Nyström pre-images for flex_analysis
module simple_flex_diffmap_rec3D
use simple_core_module_api
use simple_builder,                only: builder
use simple_gridding,               only: prep3D_inv_instrfun4mul
use simple_image,                  only: image
use simple_matcher_3Drec,          only: init_rec, prep_imgs4rec, cleanup_rec_buffers
use simple_matcher_ptcl_io,        only: discrete_read_imgbatch, discrete_read_imgbatch_source, prepimgbatch
use simple_parameters,             only: parameters
use simple_flex_projected_latent_model, only: update_basis_from_latents, write_mstep_stats_part_file, &
    &update_basis_from_mstep_stats_part_files, cleanup_planes, test_projected_model_plane_preparation
use simple_reconstructor,          only: reconstructor
implicit none
private
#include "simple_local_flags.inc"

public :: reconstruct_flex_diffmap_states, write_flex_diffmap_rec_parts, reduce_flex_diffmap_rec_parts
public :: cleanup_flex_diffmap_rec_parts
public :: test_fake_preimage_against_reconstruct3D

character(len=*), parameter :: TEST_FAKE_VOL = 'flex_fake_preimage_test.mrc'
character(len=*), parameter :: TEST_REF_VOL  = 'flex_reconstruct3D_reference_test.mrc'
character(len=*), parameter :: TEST_DIFF_VOL = 'flex_fake_preimage_difference_test.mrc'
character(len=*), parameter :: TEST_FSC_FILE = 'flex_fake_preimage_fsc_test.txt'
character(len=*), parameter :: TEST_METRICS_FILE = 'flex_fake_preimage_metrics_test.txt'
real(dp), parameter :: TEST_MIN_CC = 0.995d0
real(dp), parameter :: TEST_MAX_RAW_REL_L2 = 0.10d0
real(dp), parameter :: TEST_MAX_SCALED_REL_L2 = 0.05d0

contains

    !> Data-driven reconstruction identity test.  With one constant latent
    !! coordinate z_i=1 and target coefficient 1, the projection-aware model
    !! must satisfy mean + residual == the ordinary reconstruct3D solution.
    !! Output-side masking and filtering are disabled so this is an algebraic
    !! reconstruction test rather than a policy/postprocessing comparison.
    subroutine test_fake_preimage_against_reconstruct3D( params, build, pinds )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:)
        type(parameters) :: test_params
        type(reconstructor) :: reference_rec
        type(image) :: fake_img, diff_img, reference_ft, fake_ft
        real, allocatable :: z(:,:), target_coeffs(:,:), ref_data(:,:,:), fake_data(:,:,:)
        real, allocatable :: fsc(:), resolutions(:)
        real(dp) :: dot_rf, norm_ref2, norm_fake2, scale_fake, raw_rel_l2, scaled_rel_l2, cc
        integer :: funit, k
        if( size(pinds)<3 ) THROW_HARD('fake pre-image reconstruction test requires at least three particles')
        test_params = params
        test_params%lp       = 0.
        test_params%msk_crop = 0.
        test_params%ml_reg   = 'no'
        test_params%l_ml_reg = .false.
        test_params%outvol   = TEST_FAKE_VOL
        allocate(z(size(pinds),1),source=1.0)
        allocate(target_coeffs(1,1),source=1.0)
        write(logfhandle,'(A,I0)') '>>> FLEX PRE-IMAGE IDENTITY TEST PARTICLES: ',size(pinds)
        call test_projected_model_plane_preparation(test_params,build,pinds)
        write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE IDENTITY TEST: ORDINARY RECONSTRUCT3D REFERENCE'
        call reconstruct3D_reference(test_params,build,pinds,reference_rec)
        call reference_rec%write(string(TEST_REF_VOL),del_if_exists=.true.)
        write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE IDENTITY TEST: CONSTANT z=1 PRE-IMAGE'
        call reconstruct_flex_diffmap_states(test_params,build,pinds,z,target_coeffs,1)
        call fake_img%read(string(TEST_FAKE_VOL))
        call diff_img%copy(fake_img)
        call diff_img%subtr(reference_rec)
        call diff_img%write(string(TEST_DIFF_VOL),del_if_exists=.true.)
        ref_data  = reference_rec%get_rmat()
        fake_data = fake_img%get_rmat()
        dot_rf     = sum(real(ref_data,dp)*real(fake_data,dp))
        norm_ref2  = sum(real(ref_data,dp)*real(ref_data,dp))
        norm_fake2 = sum(real(fake_data,dp)*real(fake_data,dp))
        if( norm_ref2<=DTINY .or. norm_fake2<=DTINY ) &
            &THROW_HARD('fake pre-image reconstruction test produced an empty volume')
        scale_fake = dot_rf/norm_fake2
        raw_rel_l2 = sqrt(sum((real(fake_data,dp)-real(ref_data,dp))**2)/norm_ref2)
        scaled_rel_l2 = sqrt(sum((scale_fake*real(fake_data,dp)-real(ref_data,dp))**2)/norm_ref2)
        cc = reference_rec%corr(fake_img)
        call reference_ft%copy(reference_rec)
        call fake_ft%copy(fake_img)
        call reference_ft%fft
        call fake_ft%fft
        allocate(fsc(fdim(params%box_crop)-1),source=0.)
        call reference_ft%fsc(fake_ft,fsc)
        resolutions=get_resarr(params%box_crop,params%smpd_crop)
        open(newunit=funit,file=TEST_FSC_FILE,status='replace',action='write')
        write(funit,'(A)') '# shell resolution_A correlation'
        do k=1,min(size(fsc),size(resolutions))
            write(funit,'(I8,1X,F12.5,1X,F12.7)') k,resolutions(k),fsc(k)
        end do
        close(funit)
        open(newunit=funit,file=TEST_METRICS_FILE,status='replace',action='write')
        write(funit,'(A,I0)') 'particles=',size(pinds)
        write(funit,'(A,ES16.8)') 'fourier_cc=',cc
        write(funit,'(A,ES16.8)') 'optimal_fake_scale=',scale_fake
        write(funit,'(A,ES16.8)') 'relative_l2_raw=',raw_rel_l2
        write(funit,'(A,ES16.8)') 'relative_l2_after_scale=',scaled_rel_l2
        write(funit,'(A,ES16.8)') 'minimum_required_cc=',TEST_MIN_CC
        write(funit,'(A,ES16.8)') 'maximum_allowed_raw_relative_l2=',TEST_MAX_RAW_REL_L2
        write(funit,'(A,ES16.8)') 'maximum_allowed_scaled_relative_l2=',TEST_MAX_SCALED_REL_L2
        close(funit)
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
            &'>>> FLEX PRE-IMAGE IDENTITY METRICS cc=',cc,' optimal_scale=',scale_fake, &
            &' raw_relative_l2=',raw_rel_l2,' scaled_relative_l2=',scaled_rel_l2
        write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE IDENTITY OUTPUTS: '//TEST_REF_VOL//', '//TEST_FAKE_VOL//', '//TEST_DIFF_VOL
        call reference_rec%dealloc_rho
        call reference_rec%kill
        call fake_img%kill
        call diff_img%kill
        call reference_ft%kill
        call fake_ft%kill
        deallocate(z,target_coeffs,ref_data,fake_data,fsc,resolutions)
        if( cc<TEST_MIN_CC .or. raw_rel_l2>TEST_MAX_RAW_REL_L2 .or. &
            &scaled_rel_l2>TEST_MAX_SCALED_REL_L2 .or. scale_fake<=0.d0 )then
            THROW_HARD('constant flex pre-image does not reproduce reconstruct3D; inspect flex_fake_preimage_*_test outputs')
        endif
        write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE IDENTITY TEST PASSED'
    end subroutine test_fake_preimage_against_reconstruct3D

    subroutine reconstruct3D_reference( params, build, pinds, recvol )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:)
        type(reconstructor), intent(inout) :: recvol
        type(fplane_type), allocatable :: fpls(:)
        type(image) :: gridcorr_img
        type(ori) :: orientation
        integer :: batchlims(2), batchsz, ibatch, i, iptcl
        call init_basis_reconstructor(params,build,recvol)
        call init_rec(params,build,MAXIMGBATCHSZ,fpls)
        call prepimgbatch(params,build,MAXIMGBATCHSZ)
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
            do i=1,batchsz
                iptcl=pinds(batchlims(1)+i-1)
                call build%spproj_field%get_ori(iptcl,orientation)
                if( orientation%isstatezero() ) cycle
                call recvol%insert_plane_oversamp(build%pgrpsyms,orientation,fpls(i))
            end do
        end do
        call orientation%kill
        call cleanup_rec_buffers(build,fpls)
        call recvol%compress_exp
        call recvol%sampl_dens_correct
        call recvol%ifft
        call recvol%div(real(params%box))
        gridcorr_img=prep3D_inv_instrfun4mul([params%box_crop,params%box_crop,params%box_crop], &
            &OSMPL_PAD_FAC*[params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
        call recvol%mul(gridcorr_img)
        call gridcorr_img%kill
    end subroutine reconstruct3D_reference

    subroutine reconstruct_flex_diffmap_states( params, build, pinds, z, target_coeffs, nstates )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer, intent(in) :: pinds(:), nstates
        real, intent(in) :: z(:,:), target_coeffs(:,:)
        type(parameters) :: model_params
        type(reconstructor) :: mean_rec
        type(reconstructor), allocatable :: basis_recs(:)
        type(fplane_type), allocatable :: fpls(:)
        real(dp), allocatable :: z_dp(:,:)
        integer :: ncomp, q
        call validate_model_tables(pinds,z,target_coeffs,nstates,ncomp)
        call prepare_unfiltered_model_params(params,model_params)
        allocate(z_dp(size(z,1),ncomp),source=real(z,dp))
        allocate(basis_recs(ncomp))
        call init_mean_reconstructor(model_params,build,mean_rec)
        do q=1,ncomp
            call init_basis_reconstructor(model_params,build,basis_recs(q))
        end do
        call update_basis_from_latents(model_params,build,mean_rec,basis_recs,z_dp,pinds,size(pinds),ncomp,fpls, &
            &log_label='FLEX_DIFFMAP')
        call write_preimage_states(params,basis_recs,target_coeffs,nstates,ncomp)
        call cleanup_planes(fpls)
        call cleanup_reconstructors(mean_rec,basis_recs)
        deallocate(z_dp)
    end subroutine reconstruct_flex_diffmap_states

    subroutine write_flex_diffmap_rec_parts( params, build, pinds, z, ncomp, part )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        integer, intent(in) :: pinds(:), ncomp, part
        real, intent(in) :: z(:,:)
        type(parameters) :: model_params
        type(reconstructor) :: mean_rec
        type(fplane_type), allocatable :: fpls(:)
        real(dp), allocatable :: z_dp(:,:)
        type(string) :: fname
        if( size(pinds)<1 .or. part<1 .or. ncomp<1 ) THROW_HARD('invalid flex residual-model worker assignment')
        if( size(z,1)/=size(pinds) .or. size(z,2)/=ncomp ) THROW_HARD('flex residual-model worker coordinate mismatch')
        call prepare_unfiltered_model_params(params,model_params)
        allocate(z_dp(size(z,1),ncomp),source=real(z,dp))
        call init_mean_reconstructor(model_params,build,mean_rec)
        fname=flex_part_fname(part,params%numlen)
        call write_mstep_stats_part_file(model_params,build,mean_rec,z_dp,pinds,size(pinds),ncomp,[1,size(pinds)], &
            &fname,fpls,log_label='FLEX_DIFFMAP_WORKER')
        call cleanup_planes(fpls)
        call mean_rec%dealloc_rho
        call mean_rec%kill
        call fname%kill
        deallocate(z_dp)
    end subroutine write_flex_diffmap_rec_parts

    subroutine reduce_flex_diffmap_rec_parts( params, build, nparts, ncomp, target_coeffs, nstates )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        integer, intent(in) :: nparts, ncomp, nstates
        real, intent(in) :: target_coeffs(:,:)
        type(parameters) :: model_params
        type(reconstructor), allocatable :: basis_recs(:)
        type(string), allocatable :: part_fnames(:)
        integer :: ipart, q
        if( nparts<1 .or. ncomp<1 .or. nstates<1 ) THROW_HARD('invalid distributed flex residual-model reduction')
        if( any(shape(target_coeffs)/=[nstates,ncomp]) ) THROW_HARD('flex target coefficient shape mismatch')
        call prepare_unfiltered_model_params(params,model_params)
        allocate(basis_recs(ncomp),part_fnames(nparts))
        do q=1,ncomp
            call init_basis_reconstructor(model_params,build,basis_recs(q))
        end do
        do ipart=1,nparts
            part_fnames(ipart)=flex_part_fname(ipart,params%numlen)
        end do
        call update_basis_from_mstep_stats_part_files(model_params,basis_recs,ncomp,part_fnames,nparts, &
            &log_label='FLEX_DIFFMAP')
        call write_preimage_states(params,basis_recs,target_coeffs,nstates,ncomp)
        do q=1,ncomp
            call basis_recs(q)%dealloc_rho
            call basis_recs(q)%kill
        end do
        do ipart=1,nparts
            call part_fnames(ipart)%kill
        end do
        deallocate(basis_recs,part_fnames)
    end subroutine reduce_flex_diffmap_rec_parts

    subroutine cleanup_flex_diffmap_rec_parts( nparts, numlen )
        integer, intent(in) :: nparts, numlen
        type(string) :: fname
        integer :: ipart
        do ipart=1,nparts
            fname=flex_part_fname(ipart,numlen)
            call del_file(fname)
            call fname%kill
        end do
    end subroutine cleanup_flex_diffmap_rec_parts

    subroutine validate_model_tables( pinds, z, target_coeffs, nstates, ncomp )
        integer, intent(in) :: pinds(:), nstates
        real, intent(in) :: z(:,:), target_coeffs(:,:)
        integer, intent(out) :: ncomp
        ncomp=size(z,2)
        if( size(pinds)<1 .or. ncomp<1 .or. nstates<1 ) THROW_HARD('invalid flex residual-model dimensions')
        if( size(z,1)/=size(pinds) ) THROW_HARD('flex residual-model particle coordinate mismatch')
        if( any(shape(target_coeffs)/=[nstates,ncomp]) ) THROW_HARD('flex residual-model target mismatch')
    end subroutine validate_model_tables

    subroutine prepare_unfiltered_model_params( params, model_params )
        class(parameters), intent(in) :: params
        type(parameters), intent(out) :: model_params
        model_params=params
        model_params%lp=0.
    end subroutine prepare_unfiltered_model_params

    subroutine init_mean_reconstructor( params, build, mean_rec )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        call mean_rec%read_and_crop(params%vols(1),params%smpd,params%box_crop,params%smpd_crop)
        call mean_rec%fft
        call mean_rec%alloc_rho(params,build%spproj,expand=.true.)
        call mean_rec%expand_exp
    end subroutine init_mean_reconstructor

    subroutine init_basis_reconstructor( params, build, basis_rec )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        type(reconstructor), intent(inout) :: basis_rec
        call basis_rec%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
        call basis_rec%alloc_rho(params,build%spproj,expand=.true.)
        call basis_rec%reset
        call basis_rec%reset_exp
    end subroutine init_basis_reconstructor

    subroutine write_preimage_states( params, basis_recs, target_coeffs, nstates, ncomp )
        class(parameters), intent(in) :: params
        integer, intent(in) :: nstates, ncomp
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        real, intent(in) :: target_coeffs(nstates,ncomp)
        type(image) :: mean_img, state_img
        integer :: state, q
        call mean_img%read_and_crop(params%vols(1),params%smpd,params%box_crop,params%smpd_crop)
        call mean_img%fft
        do state=1,nstates
            call state_img%copy(mean_img)
            do q=1,ncomp
                if( abs(target_coeffs(state,q))>real(DTINY) ) &
                    &call state_img%add(basis_recs(q),target_coeffs(state,q))
            end do
            call state_img%ifft
            call write_state(params,state_img,state)
            call state_img%kill
        end do
        call mean_img%kill
    end subroutine write_preimage_states

    subroutine write_state( params, img, state )
        class(parameters), intent(in) :: params
        type(image), intent(inout) :: img
        integer, intent(in) :: state
        type(string) :: fname, prefix, ext
        character(len=:), allocatable :: stem
        character(len=3) :: tag
        if( state==1 )then
            fname=params%outvol
        else
            ext=fname2ext(params%outvol)
            prefix=get_fbody(params%outvol,ext)
            stem=prefix%to_char()
            if( len_trim(stem)>4 )then
                if( stem(len_trim(stem)-3:len_trim(stem))=='_001' ) stem=stem(:len_trim(stem)-4)
            endif
            prefix=string(stem)
            write(tag,'(I3.3)') state
            fname=prefix//'_'//tag//MRC_EXT
        endif
        call img%write(fname,del_if_exists=.true.)
        write(logfhandle,'(A,I0,A,A)') '>>> FLEX DIFFMAP NYSTROM PRE-IMAGE ',state,': ',fname%to_char()
        call fname%kill
        call prefix%kill
        call ext%kill
    end subroutine write_state

    subroutine cleanup_reconstructors( mean_rec, basis_recs )
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), allocatable, intent(inout) :: basis_recs(:)
        integer :: q
        call mean_rec%dealloc_rho
        call mean_rec%kill
        if( allocated(basis_recs) )then
            do q=1,size(basis_recs)
                call basis_recs(q)%dealloc_rho
                call basis_recs(q)%kill
            end do
            deallocate(basis_recs)
        endif
    end subroutine cleanup_reconstructors

    function flex_part_fname( part, numlen ) result(fname)
        integer, intent(in) :: part, numlen
        type(string) :: fname
        fname=string('flex_diffmap_mstep_part')//int2str_pad(part,numlen)//'.bin'
    end function flex_part_fname

end module simple_flex_diffmap_rec3D
