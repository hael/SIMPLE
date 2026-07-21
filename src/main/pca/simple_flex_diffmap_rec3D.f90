!@descr: Hermitian 3D reconstruction of diffusion-map coordinates
module simple_flex_diffmap_rec3D
use simple_core_module_api
use simple_builder,          only: builder
use simple_image,            only: image
use simple_matcher_3Drec,    only: init_rec
use simple_matcher_ptcl_io,  only: discrete_read_imgbatch, discrete_read_imgbatch_source, prepimgbatch, killimgbatch
use simple_memoize_ft_maps,  only: memoize_ft_maps, forget_ft_maps
use simple_ori,              only: ori
use simple_parameters,       only: parameters
use simple_reconstructor,    only: reconstructor
use simple_reconstructor_latent_ops, only: insert_plane_oversamp_multi_scaled, project_fplane_mean
implicit none
private
#include "simple_local_flags.inc"

public :: reconstruct_flex_diffmap_modes, write_flex_diffmap_rec_parts, reduce_flex_diffmap_rec_parts
public :: cleanup_flex_diffmap_rec_parts

integer, parameter :: FLEX_REC_MODE_BLOCK_MAX=4

contains

    subroutine reconstruct_flex_diffmap_modes( params, build, pinds, coords, nmodes )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nmodes
        real,              intent(in)    :: coords(:,:)
        type(reconstructor), allocatable :: mode_recs(:)
        integer :: q,qfirst,qlast,nblock,blocksz
        blocksz=reconstruction_block_size(params,build,nmodes)
        do qfirst=1,nmodes,blocksz
            qlast=min(nmodes,qfirst+blocksz-1)
            nblock=qlast-qfirst+1
            call accumulate_flex_modes(params,build,pinds,coords(:,qfirst:qlast),nblock,1,size(pinds),mode_recs)
            do q=1,nblock
                call finalize_mode(params, mode_recs(q))
                call write_mode(params, mode_recs(q), qfirst+q-1)
                call mode_recs(q)%dealloc_rho
                call mode_recs(q)%kill
            end do
            deallocate(mode_recs)
        end do
    end subroutine reconstruct_flex_diffmap_modes

    subroutine write_flex_diffmap_rec_parts( params, build, pinds, coords, nmodes, part )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        integer, intent(in) :: pinds(:),nmodes,part
        real, intent(in) :: coords(:,:)
        type(reconstructor), allocatable :: recs(:)
        type(string) :: vol_fname, rho_fname
        integer :: q,qfirst,qlast,nblock,blocksz
        if( size(pinds)<1 .or. part<1 ) THROW_HARD('invalid flex reconstruction worker assignment')
        blocksz=reconstruction_block_size(params,build,nmodes)
        do qfirst=1,nmodes,blocksz
            qlast=min(nmodes,qfirst+blocksz-1)
            nblock=qlast-qfirst+1
            call accumulate_flex_modes(params,build,pinds,coords(:,qfirst:qlast),nblock,1,size(pinds),recs)
            do q=1,nblock
                call recs(q)%compress_exp
                call flex_part_fnames(qfirst+q-1,part,params%numlen,vol_fname,rho_fname)
                call recs(q)%write(vol_fname,del_if_exists=.true.)
                call recs(q)%write_rho(rho_fname)
                call recs(q)%dealloc_rho; call recs(q)%kill
                call vol_fname%kill; call rho_fname%kill
            end do
            deallocate(recs)
        end do
    end subroutine write_flex_diffmap_rec_parts

    subroutine reduce_flex_diffmap_rec_parts( params, build, nparts, nmodes )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        integer, intent(in) :: nparts, nmodes
        type(reconstructor), allocatable :: recs(:)
        type(reconstructor) :: rec_read
        type(string) :: vol_fname, rho_fname
        integer :: q,iq,qfirst,qlast,nblock,blocksz,ipart
        if( nparts<1 .or. nmodes<1 ) THROW_HARD('invalid distributed flex reconstruction reduction')
        blocksz=reconstruction_block_size(params,build,nmodes)
        do qfirst=1,nmodes,blocksz
            qlast=min(nmodes,qfirst+blocksz-1)
            nblock=qlast-qfirst+1
            allocate(recs(nblock))
            do q=1,nblock; call init_reduced_mode(params,build,recs(q)); end do
            call init_reduced_mode(params,build,rec_read)
            do ipart=1,nparts
                do iq=qfirst,qlast
                    call flex_part_fnames(iq,ipart,params%numlen,vol_fname,rho_fname)
                    if( .not.file_exists(vol_fname) .or. .not.file_exists(rho_fname) )then
                        THROW_HARD('missing flex reconstruction volume/rho partial')
                    endif
                    call rec_read%read(vol_fname)
                    call rec_read%read_rho(rho_fname)
                    q=iq-qfirst+1
                    call recs(q)%sum_reduce(rec_read)
                    call vol_fname%kill; call rho_fname%kill
                end do
            end do
            call rec_read%dealloc_rho; call rec_read%kill
            do q=1,nblock
                call finalize_mode(params,recs(q)); call write_mode(params,recs(q),qfirst+q-1)
                call recs(q)%dealloc_rho; call recs(q)%kill
            end do
            deallocate(recs)
        end do
    end subroutine reduce_flex_diffmap_rec_parts

    subroutine cleanup_flex_diffmap_rec_parts( nparts, nmodes, numlen )
        integer, intent(in) :: nparts,nmodes,numlen
        type(string) :: vol_fname,rho_fname
        integer :: ipart,q
        do ipart=1,nparts
            do q=1,nmodes
                call flex_part_fnames(q,ipart,numlen,vol_fname,rho_fname)
                call del_file(vol_fname); call del_file(rho_fname)
                call vol_fname%kill; call rho_fname%kill
            end do
        end do
    end subroutine cleanup_flex_diffmap_rec_parts

    integer function reconstruction_block_size( params, build, nmodes, lb, ub ) result(blocksz)
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nmodes
        integer, optional, intent(out)   :: lb(3),ub(3)
        type(reconstructor) :: probe
        integer(kind=8) :: bytes_per_mode
        call init_mode(params,build,probe)
        if( present(lb) ) lb=lbound(probe%cmat_exp)
        if( present(ub) ) ub=ubound(probe%cmat_exp)
        bytes_per_mode=int(size(probe%cmat_exp),8)*int(storage_size(probe%cmat_exp)/8,8) &
            &+int(size(probe%rho_exp),8)*int(storage_size(probe%rho_exp)/8,8)
        call probe%dealloc_rho
        call probe%kill
        blocksz=min(nmodes,FLEX_REC_MODE_BLOCK_MAX)
        write(logfhandle,'(A,I0,A,F10.1,A,F10.1)') '>>> FLEX DIFFMAP reconstruction_mode_block=',blocksz, &
            &' per_mode_MiB=',real(bytes_per_mode,dp)/(1024._dp**2), &
            &' modes_plus_mean_MiB=',real(blocksz+1,dp)*real(bytes_per_mode,dp)/(1024._dp**2)
        call flush(logfhandle)
    end function reconstruction_block_size

    subroutine accumulate_flex_modes( params, build, pinds, coords, nmodes, first, last, mode_recs )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        integer, intent(in) :: pinds(:),nmodes,first,last
        real, intent(in) :: coords(:,:)
        type(reconstructor), allocatable, intent(out) :: mode_recs(:)
        type(reconstructor) :: mean_rec
        type(fplane_type), allocatable :: fpls(:)
        type(fplane_type) :: mean_fpl
        type(ori) :: orientation
        real(dp), allocatable :: data_scales(:),density_scales(:)
        integer :: ibatch,batch_end,batchsz,i,row,q
        if( nmodes<1 .or. first<1 .or. last>size(pinds) .or. first>last ) THROW_HARD('invalid flex reconstruction range')
        if( size(coords,1)/=size(pinds) .or. size(coords,2)<nmodes ) THROW_HARD('flex reconstruction coordinate mismatch')
        call init_mean(params,build,mean_rec)
        allocate(mode_recs(nmodes),data_scales(nmodes),density_scales(nmodes))
        do q=1,nmodes; call init_mode(params,build,mode_recs(q)); end do
        call init_rec(params,build,MAXIMGBATCHSZ,fpls); call prepimgbatch(params,build,MAXIMGBATCHSZ)
        do ibatch=first,last,MAXIMGBATCHSZ
            batch_end=min(last,ibatch+MAXIMGBATCHSZ-1); batchsz=batch_end-ibatch+1
            call read_particles(params,build,size(pinds),pinds,[ibatch,batch_end],batchsz)
            call prepare_fplanes(params,build,batchsz,pinds(ibatch:batch_end),fpls(:batchsz))
            do i=1,batchsz
                row=ibatch+i-1; call build%spproj_field%get_ori(pinds(row),orientation)
                if( orientation%isstatezero() ) cycle
                call project_fplane_mean(mean_rec,orientation,fpls(i),mean_fpl,apply_ctf_amp=.true.)
                fpls(i)%cmplx_plane=fpls(i)%cmplx_plane-mean_fpl%cmplx_plane
                data_scales=real(coords(row,1:nmodes),dp); density_scales=data_scales*data_scales
                call insert_plane_oversamp_multi_scaled(mode_recs,build%pgrpsyms,orientation,fpls(i),data_scales,density_scales)
            end do
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX DIFFMAP RECONSTRUCTION: ',batch_end,' / ',last
        end do
        call orientation%kill; call cleanup_planes(fpls); call killimgbatch(build); call forget_ft_maps
        call cleanup_plane(mean_fpl); call mean_rec%dealloc_rho; call mean_rec%kill
        deallocate(data_scales,density_scales)
    end subroutine accumulate_flex_modes

    subroutine init_mean( params, build, rec )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        type(reconstructor), intent(inout) :: rec
        call rec%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
        call rec%fft
        call rec%alloc_rho(params, build%spproj, expand=.true.)
        call rec%expand_exp
    end subroutine init_mean

    subroutine init_mode( params, build, rec )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        type(reconstructor), intent(inout) :: rec
        call rec%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        call rec%alloc_rho(params, build%spproj, expand=.true.)
        call rec%reset
        call rec%reset_exp
    end subroutine init_mode

    subroutine init_reduced_mode( params, build, rec )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        type(reconstructor), intent(inout) :: rec
        call rec%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        call rec%alloc_rho(params, build%spproj, expand=.false.)
        call rec%reset
    end subroutine init_reduced_mode

    subroutine read_particles( params, build, nptcls, pinds, lims, batchsz )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer, intent(in) :: nptcls, pinds(nptcls), lims(2), batchsz
        if( trim(params%ptcl_src) == 'raw' )then
            call discrete_read_imgbatch(params, build, nptcls, pinds, lims)
        else
            call discrete_read_imgbatch_source(params, build, trim(params%ptcl_src), nptcls, pinds, lims, &
                &build%imgbatch(:batchsz))
        endif
    end subroutine read_particles

    subroutine prepare_fplanes( params, build, nptcls, pinds, fplanes )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls, pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        type(ctfparams) :: ctfparms(nthr_glob)
        real :: shift(2)
        integer :: i, iptcl, ithr, kfromto(2)
        call memoize_ft_maps([params%boxpd,params%boxpd,1],params%smpd)
        kfromto = flex_kfromto(params)
        !$omp parallel do default(shared) private(i,iptcl,ithr,shift) schedule(static) proc_bind(close)
        do i = 1,nptcls
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds(i)
            call build%imgbatch(i)%norm_noise_taper_edge_pad_fft(build%lmsk, build%img_pad_heap(ithr))
            ctfparms(ithr) = build%spproj%get_ctfparams('ptcl3D',iptcl)
            shift = build%spproj_field%get_2Dshift(iptcl)
            call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr), shift, &
                &fplanes(i), store_transfer=.true., observation_model=.true.)
            if( fplanes(i)%nyq > 0 ) fplanes(i)%nyq = min(fplanes(i)%nyq, OSMPL_PAD_FAC*kfromto(2))
        end do
        !$omp end parallel do
    end subroutine prepare_fplanes

    function flex_kfromto( params ) result(kfromto)
        class(parameters), intent(in) :: params
        integer :: kfromto(2), kto_full
        real :: dstep
        kto_full = max(1,fdim(params%box_crop)-1)
        kfromto = [1,kto_full]
        if( params%lp > 2.*params%smpd_crop + TINY )then
            dstep = real(max(1,params%box_crop-1))*params%smpd_crop
            kfromto(2) = max(1,min(kto_full,int(dstep/params%lp)))
        endif
    end function flex_kfromto

    subroutine finalize_mode( params, rec )
        class(parameters), intent(in)    :: params
        type(reconstructor), intent(inout) :: rec
        if( allocated(rec%cmat_exp) ) call rec%compress_exp
        call rec%sampl_dens_correct
        call rec%ifft
        call rec%div(real(params%box))
        if( params%msk_crop > TINY ) call rec%mask3D_soft(params%msk_crop,backgr=0.)
        if( params%lp > 2.*params%smpd_crop + TINY ) call rec%bp(0.,params%lp)
        if( params%msk_crop > TINY ) call rec%mask3D_soft(params%msk_crop,backgr=0.)
    end subroutine finalize_mode

    subroutine write_mode( params, rec, q )
        class(parameters), intent(in)    :: params
        type(reconstructor), intent(inout) :: rec
        integer, intent(in) :: q
        type(string) :: fname, prefix, ext
        character(len=:), allocatable :: stem
        character(len=3) :: tag
        if( q == 1 )then
            fname = params%outvol
        else
            ext = fname2ext(params%outvol)
            prefix = get_fbody(params%outvol,ext)
            stem = prefix%to_char()
            if( len_trim(stem)>4 )then
                if( stem(len_trim(stem)-3:len_trim(stem))=='_001' ) stem=stem(:len_trim(stem)-4)
            endif
            prefix=string(trim(stem))
            write(tag,'(I3.3)') q
            fname = prefix//'_'//tag//MRC_EXT
        endif
        call rec%write(fname,del_if_exists=.true.)
        write(logfhandle,'(A,I0,A,A)') '>>> FLEX DIFFMAP MODE ',q,': ',fname%to_char()
        call fname%kill
        call prefix%kill
        call ext%kill
    end subroutine write_mode

    subroutine flex_part_fnames( mode, part, numlen, vol_fname, rho_fname )
        integer, intent(in) :: mode,part,numlen
        type(string), intent(out) :: vol_fname,rho_fname
        vol_fname = string('flex_diffmap_mode_')//int2str_pad(mode,3)//'_part'//int2str_pad(part,numlen)//MRC_EXT
        rho_fname = string('rho_')//vol_fname
    end subroutine flex_part_fnames

    subroutine cleanup_plane( fpl )
        type(fplane_type), intent(inout) :: fpl
        if( allocated(fpl%cmplx_plane) ) deallocate(fpl%cmplx_plane)
        if( allocated(fpl%ctfsq_plane) ) deallocate(fpl%ctfsq_plane)
        if( allocated(fpl%transfer_plane) ) deallocate(fpl%transfer_plane)
    end subroutine cleanup_plane

    subroutine cleanup_planes( fpls )
        type(fplane_type), allocatable, intent(inout) :: fpls(:)
        integer :: i
        if( .not. allocated(fpls) ) return
        do i = 1,size(fpls)
            call cleanup_plane(fpls(i))
        end do
        deallocate(fpls)
    end subroutine cleanup_planes

end module simple_flex_diffmap_rec3D
