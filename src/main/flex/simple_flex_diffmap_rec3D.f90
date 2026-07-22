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
    &update_basis_from_mstep_stats_part_files, cleanup_planes, test_projected_model_plane_preparation, &
    &prep_imgs4projected_model, solve_coupled_basis_exp
use simple_flex_reconstructor_latent_ops, only: project_fplane_mean
use simple_flex_reconstructor_latent_ops, only: insert_planes_oversamp_coupled_batch_scaled
use simple_reconstructor,          only: reconstructor
implicit none
private
#include "simple_local_flags.inc"

public :: reconstruct_flex_diffmap_states, write_flex_diffmap_rec_parts, reduce_flex_diffmap_rec_parts
public :: cleanup_flex_diffmap_rec_parts
public :: test_fake_preimage_against_reconstruct3D
public :: canonicalize_flex_preimage_coordinates

character(len=*), parameter :: TEST_FAKE_VOL = 'flex_fake_preimage_test.mrc'
character(len=*), parameter :: TEST_REF_VOL  = 'flex_reconstruct3D_reference_test.mrc'
character(len=*), parameter :: TEST_DIFF_VOL = 'flex_fake_preimage_difference_test.mrc'
character(len=*), parameter :: TEST_STANDARD_RESIDUAL_VOL = 'flex_standard_residual_identity_test.mrc'
character(len=*), parameter :: TEST_STANDARD_RESIDUAL_DIFF_VOL = 'flex_standard_residual_difference_test.mrc'
character(len=*), parameter :: TEST_FSC_FILE = 'flex_fake_preimage_fsc_test.txt'
character(len=*), parameter :: TEST_METRICS_FILE = 'flex_fake_preimage_metrics_test.txt'
real(dp), parameter :: TEST_MIN_CC = 0.995d0
real(dp), parameter :: TEST_MAX_RAW_REL_L2 = 0.10d0
real(dp), parameter :: TEST_MAX_SCALED_REL_L2 = 0.05d0
! The coupled and ordinary paths accumulate the same single-precision
! Fourier samples in different valid orders.  Their solved grids agree at
! approximately 2.3e-3 on the production identity data; retain margin for
! this arithmetic difference while keeping the map-level check much stricter.
real(dp), parameter :: TEST_MAX_COUPLED_SOLVE_REL_ERR = 5.0d-3

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
        type(image) :: fake_img, diff_img, reference_ft, fake_ft, standard_residual_img, standard_residual_diff_img
        real, allocatable :: z(:,:), target_coeffs(:,:), ref_data(:,:,:), fake_data(:,:,:), standard_residual_data(:,:,:)
        real, allocatable :: fsc(:), resolutions(:)
        real(dp) :: dot_rf, norm_ref2, norm_fake2, scale_fake, raw_rel_l2, scaled_rel_l2, cc
        real(dp) :: dot_rs, norm_standard_residual2, scale_standard_residual, standard_residual_raw_rel_l2
        real(dp) :: standard_residual_scaled_rel_l2, standard_residual_cc
        real(dp) :: dot_fc, scale_fake_to_standard_residual, fake_to_standard_residual_raw_rel_l2
        real(dp) :: fake_to_standard_residual_scaled_rel_l2, fake_to_standard_residual_cc
        real(dp) :: coupled_rhs_relative_error, coupled_rho_relative_error
        real(dp) :: coupled_solution_relative_error
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
        ref_data = reference_rec%get_rmat()
        norm_ref2 = sum(real(ref_data,dp)*real(ref_data,dp))
        if( norm_ref2<=DTINY ) THROW_HARD('ordinary reconstruct3D reference produced an empty volume')
        write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE IDENTITY TEST: STANDARD RESIDUAL CONTROL'
        call reconstruct3D_standard_residual_control(test_params,build,pinds,standard_residual_img, &
            &coupled_rhs_relative_error,coupled_rho_relative_error,coupled_solution_relative_error)
        call standard_residual_img%write(string(TEST_STANDARD_RESIDUAL_VOL),del_if_exists=.true.)
        call standard_residual_diff_img%copy(standard_residual_img)
        call standard_residual_diff_img%subtr(reference_rec)
        call standard_residual_diff_img%write(string(TEST_STANDARD_RESIDUAL_DIFF_VOL),del_if_exists=.true.)
        standard_residual_data = standard_residual_img%get_rmat()
        dot_rs = sum(real(ref_data,dp)*real(standard_residual_data,dp))
        norm_standard_residual2 = sum(real(standard_residual_data,dp)*real(standard_residual_data,dp))
        if( norm_standard_residual2<=DTINY ) THROW_HARD('standard flex residual control produced an empty volume')
        scale_standard_residual = dot_rs/norm_standard_residual2
        standard_residual_raw_rel_l2 = sqrt(sum((real(standard_residual_data,dp)-real(ref_data,dp))**2)/norm_ref2)
        standard_residual_scaled_rel_l2 = sqrt(sum((scale_standard_residual*real(standard_residual_data,dp)- &
            &real(ref_data,dp))**2)/norm_ref2)
        standard_residual_cc = reference_rec%corr(standard_residual_img)
        write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE IDENTITY TEST: CONSTANT z=1 PRE-IMAGE'
        call reconstruct_flex_diffmap_states(test_params,build,pinds,z,target_coeffs,1)
        call fake_img%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
        call fake_img%read(string(TEST_FAKE_VOL))
        call diff_img%copy(fake_img)
        call diff_img%subtr(reference_rec)
        call diff_img%write(string(TEST_DIFF_VOL),del_if_exists=.true.)
        fake_data = fake_img%get_rmat()
        dot_rf     = sum(real(ref_data,dp)*real(fake_data,dp))
        norm_fake2 = sum(real(fake_data,dp)*real(fake_data,dp))
        if( norm_ref2<=DTINY .or. norm_fake2<=DTINY ) &
            &THROW_HARD('fake pre-image reconstruction test produced an empty volume')
        scale_fake = dot_rf/norm_fake2
        raw_rel_l2 = sqrt(sum((real(fake_data,dp)-real(ref_data,dp))**2)/norm_ref2)
        scaled_rel_l2 = sqrt(sum((scale_fake*real(fake_data,dp)-real(ref_data,dp))**2)/norm_ref2)
        cc = reference_rec%corr(fake_img)
        dot_fc = sum(real(standard_residual_data,dp)*real(fake_data,dp))
        scale_fake_to_standard_residual = dot_fc/norm_fake2
        fake_to_standard_residual_raw_rel_l2 = sqrt(sum((real(fake_data,dp)- &
            &real(standard_residual_data,dp))**2)/norm_standard_residual2)
        fake_to_standard_residual_scaled_rel_l2 = sqrt(sum((scale_fake_to_standard_residual*real(fake_data,dp)- &
            &real(standard_residual_data,dp))**2)/norm_standard_residual2)
        fake_to_standard_residual_cc = standard_residual_img%corr(fake_img)
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
        write(funit,'(A,ES16.8)') 'standard_residual_cc=',standard_residual_cc
        write(funit,'(A,ES16.8)') 'standard_residual_optimal_scale=',scale_standard_residual
        write(funit,'(A,ES16.8)') 'standard_residual_relative_l2_raw=',standard_residual_raw_rel_l2
        write(funit,'(A,ES16.8)') 'standard_residual_relative_l2_after_scale=',standard_residual_scaled_rel_l2
        write(funit,'(A,ES16.8)') 'coupled_raw_rhs_relative_error=',coupled_rhs_relative_error
        write(funit,'(A,ES16.8)') 'coupled_raw_rho_relative_error=',coupled_rho_relative_error
        write(funit,'(A,ES16.8)') 'coupled_solution_relative_error=',coupled_solution_relative_error
        write(funit,'(A,ES16.8)') 'coupled_to_standard_residual_cc=',fake_to_standard_residual_cc
        write(funit,'(A,ES16.8)') 'coupled_to_standard_residual_optimal_scale=',scale_fake_to_standard_residual
        write(funit,'(A,ES16.8)') 'coupled_to_standard_residual_relative_l2_raw=',fake_to_standard_residual_raw_rel_l2
        write(funit,'(A,ES16.8)') 'coupled_to_standard_residual_relative_l2_after_scale=', &
            &fake_to_standard_residual_scaled_rel_l2
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
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
            &'>>> FLEX STANDARD RESIDUAL CONTROL cc=',standard_residual_cc, &
            &' optimal_scale=',scale_standard_residual,' raw_relative_l2=',standard_residual_raw_rel_l2, &
            &' scaled_relative_l2=',standard_residual_scaled_rel_l2
        write(logfhandle,'(A,ES12.4,A,ES12.4)') '>>> FLEX RAW COUPLED INSERTION rhs_relative_error=', &
            &coupled_rhs_relative_error,' rho_relative_error=',coupled_rho_relative_error
        write(logfhandle,'(A,ES12.4)') '>>> FLEX COUPLED SOLVE relative_error=',coupled_solution_relative_error
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
            &'>>> FLEX COUPLED-TO-STANDARD-RESIDUAL cc=',fake_to_standard_residual_cc, &
            &' optimal_scale=',scale_fake_to_standard_residual,' raw_relative_l2=', &
            &fake_to_standard_residual_raw_rel_l2,' scaled_relative_l2=',fake_to_standard_residual_scaled_rel_l2
        write(logfhandle,'(A,5(1X,A))') '>>> FLEX PRE-IMAGE IDENTITY OUTPUTS:', TEST_REF_VOL, &
            &TEST_STANDARD_RESIDUAL_VOL, TEST_STANDARD_RESIDUAL_DIFF_VOL, TEST_FAKE_VOL, TEST_DIFF_VOL
        call reference_rec%dealloc_rho
        call reference_rec%kill
        call fake_img%kill
        call diff_img%kill
        call standard_residual_img%kill
        call standard_residual_diff_img%kill
        call reference_ft%kill
        call fake_ft%kill
        deallocate(z,target_coeffs,ref_data,fake_data,standard_residual_data,fsc,resolutions)
        if( coupled_rhs_relative_error>1.d-3 .or. coupled_rho_relative_error>1.d-3 )then
            THROW_HARD('coupled batch insertion differs from ordinary reconstruction insertion')
        endif
        if( coupled_solution_relative_error>TEST_MAX_COUPLED_SOLVE_REL_ERR )then
            THROW_HARD('coupled solve differs from ordinary sampling-density correction')
        endif
        if( fake_to_standard_residual_cc<TEST_MIN_CC .or. &
            &fake_to_standard_residual_raw_rel_l2>TEST_MAX_RAW_REL_L2 .or. &
            &fake_to_standard_residual_scaled_rel_l2>TEST_MAX_SCALED_REL_L2 .or. &
            &scale_fake_to_standard_residual<=0.d0 )then
            THROW_HARD('coupled flex pre-image does not reproduce the standard residual control; inspect flex_fake_preimage_*_test outputs')
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

    !> Reconstruct the residual against the supplied mean through the ordinary
    !! refine3D plane/insertion path.  For a constant latent z=1, adding this
    !! residual back to the mean must reproduce reconstruct3D.  This control
    !! deliberately does not use the coupled latent accumulator or solver.
    subroutine reconstruct3D_standard_residual_control( params, build, pinds, outvol, coupled_rhs_relative_error, &
        &coupled_rho_relative_error, coupled_solution_relative_error )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:)
        type(image),       intent(inout) :: outvol
        real(dp),          intent(out)   :: coupled_rhs_relative_error, coupled_rho_relative_error
        real(dp),          intent(out)   :: coupled_solution_relative_error
        type(reconstructor) :: residual_rec, mean_rec, coupled_rec(1)
        type(fplane_type), allocatable :: standard_fpls(:), model_fpls(:), mean_fpls(:)
        type(image) :: gridcorr_img, mean_img
        type(ori), allocatable :: orientations(:)
        real(dp), allocatable :: unit_z(:,:), unit_second(:,:,:)
        real, allocatable :: coupled_rho_cross(:,:,:,:)
        complex, allocatable :: standard_cmat(:,:,:), coupled_cmat(:,:,:)
        logical, allocatable :: valid(:)
        integer :: batchlims(2), batchsz, ibatch, i, iptcl
        real(dp) :: rhs_ref2, rhs_err2, rho_ref2, rho_err2, solution_ref2, solution_err2
        call init_basis_reconstructor(params,build,residual_rec)
        call init_basis_reconstructor(params,build,coupled_rec(1))
        call init_mean_reconstructor(params,build,mean_rec)
        call init_rec(params,build,MAXIMGBATCHSZ,standard_fpls)
        allocate(model_fpls(MAXIMGBATCHSZ),mean_fpls(1),orientations(MAXIMGBATCHSZ),valid(MAXIMGBATCHSZ))
        allocate(unit_z(1,MAXIMGBATCHSZ),unit_second(1,1,MAXIMGBATCHSZ),source=1.d0)
        allocate(coupled_rho_cross(1,size(coupled_rec(1)%cmat_exp,1),size(coupled_rec(1)%cmat_exp,2), &
            &size(coupled_rec(1)%cmat_exp,3)),source=0.)
        coupled_rhs_relative_error=0.d0
        coupled_rho_relative_error=0.d0
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
                &pinds(batchlims(1):batchlims(2)),standard_fpls(:batchsz))
            ! The preparation routines modify their images.  Re-read before
            ! generating the model/operator plane used for the subtraction.
            if( params%l_ptcl_src_den )then
                call discrete_read_imgbatch_source(params,build,'den',batchsz, &
                    &pinds(batchlims(1):batchlims(2)),[1,batchsz],build%imgbatch(:batchsz))
            else
                call discrete_read_imgbatch(params,build,size(pinds),pinds,batchlims)
            endif
            call prep_imgs4projected_model(params,build,batchsz,build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)),model_fpls(:batchsz))
            valid(:batchsz)=.false.
            do i=1,batchsz
                iptcl=pinds(batchlims(1)+i-1)
                call orientations(i)%kill
                call build%spproj_field%get_ori(iptcl,orientations(i))
                if( orientations(i)%isstatezero() ) cycle
                if( .not.allocated(model_fpls(i)%transfer_plane) ) &
                    &THROW_HARD('standard residual control lacks a forward transfer plane')
                call project_fplane_mean(mean_rec,orientations(i),model_fpls(i),mean_fpls(1),apply_ctf_amp=.true.)
                standard_fpls(i)%cmplx_plane = standard_fpls(i)%cmplx_plane - &
                    &conjg(model_fpls(i)%transfer_plane)*mean_fpls(1)%cmplx_plane
                model_fpls(i)%cmplx_plane = model_fpls(i)%cmplx_plane - mean_fpls(1)%cmplx_plane
                call residual_rec%insert_plane_oversamp(build%pgrpsyms,orientations(i),standard_fpls(i))
                valid(i)=.true.
            end do
            call insert_planes_oversamp_coupled_batch_scaled(coupled_rec,coupled_rho_cross,build%pgrpsyms, &
                &orientations,model_fpls,unit_z,unit_second,valid,batchsz)
        end do
        rhs_ref2=sum(real(residual_rec%cmat_exp*conjg(residual_rec%cmat_exp),dp))
        rhs_err2=sum(real((coupled_rec(1)%cmat_exp-residual_rec%cmat_exp)* &
            &conjg(coupled_rec(1)%cmat_exp-residual_rec%cmat_exp),dp))
        rho_ref2=sum(real(residual_rec%rho_exp,dp)*real(residual_rec%rho_exp,dp))
        rho_err2=sum((real(coupled_rho_cross(1,:,:,:),dp)-real(residual_rec%rho_exp,dp))**2)
        coupled_rhs_relative_error=sqrt(rhs_err2/max(rhs_ref2,DTINY))
        coupled_rho_relative_error=sqrt(rho_err2/max(rho_ref2,DTINY))
        call solve_coupled_basis_exp(coupled_rec,coupled_rho_cross,1)
        call residual_rec%compress_exp
        call residual_rec%sampl_dens_correct
        call coupled_rec(1)%compress_exp
        standard_cmat=residual_rec%get_cmat()
        coupled_cmat=coupled_rec(1)%get_cmat()
        solution_ref2=sum(real(standard_cmat*conjg(standard_cmat),dp))
        solution_err2=sum(real((coupled_cmat-standard_cmat)*conjg(coupled_cmat-standard_cmat),dp))
        coupled_solution_relative_error=sqrt(solution_err2/max(solution_ref2,DTINY))
        deallocate(standard_cmat,coupled_cmat)
        do i=1,size(orientations)
            call orientations(i)%kill
        end do
        call cleanup_planes(model_fpls)
        call cleanup_planes(mean_fpls)
        call cleanup_rec_buffers(build,standard_fpls)
        call residual_rec%ifft
        call residual_rec%div(real(params%box))
        gridcorr_img=prep3D_inv_instrfun4mul([params%box_crop,params%box_crop,params%box_crop], &
            &OSMPL_PAD_FAC*[params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
        call residual_rec%mul(gridcorr_img)
        call gridcorr_img%kill
        call mean_img%read_and_crop(params%vols(1),params%smpd,params%box_crop,params%smpd_crop)
        call residual_rec%add(mean_img)
        call outvol%copy(residual_rec)
        call mean_img%kill
        call residual_rec%dealloc_rho
        call residual_rec%kill
        call coupled_rec(1)%dealloc_rho
        call coupled_rec(1)%kill
        call mean_rec%dealloc_rho
        call mean_rec%kill
        deallocate(orientations,unit_z,unit_second,coupled_rho_cross,valid)
    end subroutine reconstruct3D_standard_residual_control

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

    !> Build a full-rank, zero-mean-preserving-free coordinate whitening for
    !! the pre-image A/B diagnostic.  This is deliberately *not* the generic
    !! projected-model canonicalization: the diffusion residual model has no
    !! intercept mode, so subtracting a column mean would alter its span.
    !!
    !! The returned coordinates obey z_canonical = z_raw * transform.  A
    !! state coefficient must therefore obey transform * phi_canonical =
    !! phi_raw.  Applying the inverse transform to targets is essential: using
    !! the same transform on both tables would compare different states.
    subroutine canonicalize_flex_preimage_coordinates( z_raw, target_raw, z_canonical, target_canonical, &
        &transform, transform_condition, gram_error, target_error )
        real, intent(in)  :: z_raw(:,:), target_raw(:,:)
        real, intent(out) :: z_canonical(size(z_raw,1),size(z_raw,2))
        real, intent(out) :: target_canonical(size(target_raw,1),size(target_raw,2))
        real(dp), intent(out) :: transform(size(z_raw,2),size(z_raw,2))
        real(dp), intent(out) :: transform_condition, gram_error, target_error
        real(dp), allocatable :: zwork(:,:)
        real(dp) :: coeff, normq, target_work(size(z_raw,2)), target_sol(size(z_raw,2))
        integer :: nptcls, ncomp, q, r, state
        if( size(z_raw,1)<1 .or. size(z_raw,2)<1 ) THROW_HARD('empty flex pre-image coordinate table')
        if( size(target_raw,2)/=size(z_raw,2) ) THROW_HARD('flex pre-image target/component mismatch')
        nptcls=size(z_raw,1)
        ncomp=size(z_raw,2)
        allocate(zwork(nptcls,ncomp),source=real(z_raw,dp))
        transform=0.d0
        do q=1,ncomp
            transform(q,q)=1.d0
            do r=1,q-1
                coeff=dot_product(zwork(:,r),zwork(:,q))/max(dot_product(zwork(:,r),zwork(:,r)),DTINY)
                zwork(:,q)=zwork(:,q)-coeff*zwork(:,r)
                transform(:,q)=transform(:,q)-coeff*transform(:,r)
            end do
            normq=sqrt(dot_product(zwork(:,q),zwork(:,q)))
            if( normq<=sqrt(DTINY)*sqrt(real(nptcls,dp)) )then
                THROW_HARD('flex pre-image canonicalization is rank deficient; reduce neigs')
            endif
            zwork(:,q)=sqrt(real(nptcls,dp))*zwork(:,q)/normq
            transform(:,q)=sqrt(real(nptcls,dp))*transform(:,q)/normq
        end do
        transform_condition=maxval(abs(transform))/max(minval(abs([(transform(q,q),q=1,ncomp)])),DTINY)
        target_error=0.d0
        do state=1,size(target_raw,1)
            target_work=real(target_raw(state,:),dp)
            call solve_upper_triangular(transform,target_work,target_sol,ncomp)
            target_canonical(state,:)=real(target_sol)
            target_error=max(target_error,maxval(abs(matmul(transform,real(target_canonical(state,:),dp))-target_work)))
        end do
        z_canonical=real(zwork)
        gram_error=maxval(abs(matmul(transpose(real(z_canonical,dp)),real(z_canonical,dp))/real(nptcls,dp)- &
            &identity_matrix(ncomp)))
        deallocate(zwork)
    contains
        function identity_matrix(n) result(identity)
            integer, intent(in) :: n
            real(dp) :: identity(n,n)
            integer :: iq
            identity=0.d0
            do iq=1,n
                identity(iq,iq)=1.d0
            end do
        end function identity_matrix

        subroutine solve_upper_triangular( upper, rhs, sol, n )
            integer, intent(in) :: n
            real(dp), intent(in) :: upper(n,n), rhs(n)
            real(dp), intent(out) :: sol(n)
            real(dp) :: sumv
            integer :: iq, ir
            sol=0.d0
            do iq=n,1,-1
                sumv=rhs(iq)
                do ir=iq+1,n
                    sumv=sumv-upper(iq,ir)*sol(ir)
                end do
                if( abs(upper(iq,iq))<=DTINY ) THROW_HARD('singular flex pre-image canonical transform')
                sol(iq)=sumv/upper(iq,iq)
            end do
        end subroutine solve_upper_triangular
    end subroutine canonicalize_flex_preimage_coordinates

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
