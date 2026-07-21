!@descr: projection-aware residual Nyström pre-images for flex_eigenvol
module simple_flex_diffmap_rec3D
use simple_core_module_api
use simple_builder,                only: builder
use simple_image,                  only: image
use simple_parameters,             only: parameters
use simple_projected_latent_model, only: update_basis_from_latents, write_mstep_stats_part_file, &
    &update_basis_from_mstep_stats_part_files, cleanup_planes
use simple_reconstructor,          only: reconstructor
implicit none
private
#include "simple_local_flags.inc"

public :: reconstruct_flex_diffmap_states, write_flex_diffmap_rec_parts, reduce_flex_diffmap_rec_parts
public :: cleanup_flex_diffmap_rec_parts

contains

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
