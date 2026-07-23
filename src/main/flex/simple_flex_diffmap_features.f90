!@descr: registered residual feature preparation for flex_analysis diffusion maps
module simple_flex_diffmap_features
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,          only: builder
use simple_classaverager,    only: transform_ptcls
use simple_ctf,              only: ctf
use simple_image,            only: image
use simple_imgarr_utils,     only: dealloc_imgarr
use simple_math_ft,          only: upsample_sigma2
use simple_memoize_ft_maps,  only: memoize_ft_maps, forget_ft_maps
use simple_ori,              only: ori
use simple_parameters,       only: parameters
use simple_simple_volinterp, only: reproject
use simple_sp_project,       only: sp_project
use simple_stack_io,         only: stack_io
implicit none
private
#include "simple_local_flags.inc"

public :: prepare_flex_diffmap_features, prepare_flex_diffmap_feature_part
public :: assemble_flex_diffmap_feature_parts, read_flex_diffmap_feature_parts
public :: flex_projection_directions, write_flex_mean_projection_stack
public :: map_flex_registered_to_native_project

contains

    !> Map the registered-frame flex model back onto native particle images.
    !! Registered images are created by applying the inverse in-plane
    !! ptcl3D transform to native images.  The registered project therefore
    !! holds the canonical view (zero in-plane transform), while this output
    !! restores the native ptcl3D in-plane transform onto that view.  ptcl2D
    !! is deliberately neither read nor modified here.
    subroutine map_flex_registered_to_native_project( native_spproj, registered_project, registered_pinds, &
        &native_project, native_pinds, manifest_fname, native_project_name )
        type(sp_project), intent(inout) :: native_spproj
        type(string), intent(in) :: registered_project
        integer, intent(in) :: registered_pinds(:)
        type(string), intent(out) :: native_project
        integer, allocatable, intent(out) :: native_pinds(:)
        character(len=*), optional, intent(in) :: manifest_fname
        character(len=*), optional, intent(in) :: native_project_name
        type(sp_project) :: registered_spproj,outproj
        type(ori) :: native_3d,registered_3d,mapped_3d
        logical, allocatable :: mapped_native(:)
        integer :: nraw,nregistered
        integer :: i,ireg,inative
        integer :: native_stkind,native_indstk,registered_stkind,registered_indstk
        integer :: u
        logical :: write_manifest
        real, parameter :: FRAME_TOL=1.e-4

        if( size(registered_pinds)<1 ) THROW_HARD('cannot map an empty flex particle selection to native images')
        call registered_spproj%read(registered_project)
        nraw=native_spproj%os_ptcl3D%get_noris()
        nregistered=registered_spproj%os_ptcl3D%get_noris()
        if( nraw<1 .or. nregistered/=nraw ) &
            &THROW_HARD('flex registered and native projects require identical ptcl3D row counts')
        allocate(native_pinds(size(registered_pinds)),source=0)
        allocate(mapped_native(nraw),source=.false.)
        do i=1,size(registered_pinds)
            ireg=registered_pinds(i)
            if( ireg<1 .or. ireg>nregistered ) THROW_HARD('registered flex particle row outside project')
            native_pinds(i)=ireg
        end do

        if( present(native_project_name) )then
            native_project=native_project_name
        else
            native_project='flex_native_model.simple'
        endif
        call del_file(native_project)
        call outproj%copy(native_spproj)
        call outproj%update_projinfo(native_project)
        do inative=1,nraw
            call outproj%os_ptcl3D%set(inative,'state',0)
        end do
        write_manifest=present(manifest_fname)
        if( write_manifest )then
            open(newunit=u,file=manifest_fname,status='replace',action='write')
            write(u,'(A)') '# mapping=ptcl3D_row'
            write(u,'(A)') '# model_row registered_row native_row native_stkind native_indstk '// &
                &'registered_stkind registered_indstk canonical_e1 canonical_e2 native_e3 native_x native_y '// &
                &'mapped_e1 mapped_e2 mapped_e3 mapped_x mapped_y'
        endif
        do i=1,size(registered_pinds)
            ireg=registered_pinds(i)
            inative=native_pinds(i)
            if( mapped_native(inative) ) THROW_HARD('duplicate registered-to-native flex particle mapping')
            if( registered_spproj%os_ptcl3D%get_state(ireg)<=0 ) &
                &THROW_HARD('inactive registered particle selected for flex native mapping')
            if( native_spproj%os_ptcl3D%get_state(inative)<=0 ) &
                &THROW_HARD('registered flex particle maps to an inactive native particle')
            call native_spproj%os_ptcl3D%get_ori(inative,native_3d)
            call registered_spproj%os_ptcl3D%get_ori(ireg,registered_3d)
            if( abs(registered_3d%e3get())>FRAME_TOL .or. &
                &any(abs(registered_3d%get_2Dshift())>FRAME_TOL) ) &
                &THROW_HARD('registered flex project does not contain canonical ptcl3D in-plane geometry')
            if( abs(registered_3d%e1get()-native_3d%e1get())>FRAME_TOL .or. &
                &abs(registered_3d%e2get()-native_3d%e2get())>FRAME_TOL ) &
                &THROW_HARD('registered and native flex projects disagree on the ptcl3D viewing direction')
            ! The second argument supplies only e3/x/y.  Its e1/e2 values are
            ! irrelevant to compose3d2d, so this restores precisely the
            ! ptcl3D transform removed by transform_ptcls.
            call registered_3d%compose3d2d(native_3d,mapped_3d)
            call outproj%os_ptcl3D%e1set(inative,mapped_3d%e1get())
            call outproj%os_ptcl3D%e2set(inative,mapped_3d%e2get())
            call outproj%os_ptcl3D%e3set(inative,mapped_3d%e3get())
            call outproj%os_ptcl3D%set_shift(inative,mapped_3d%get_2Dshift())
            call outproj%os_ptcl3D%set(inative,'state',registered_spproj%os_ptcl3D%get_state(ireg))
            if( registered_spproj%os_ptcl3D%isthere(ireg,'proj') ) &
                &call outproj%os_ptcl3D%set(inative,'proj',registered_spproj%os_ptcl3D%get_int(ireg,'proj'))
            call native_spproj%map_ptcl_ind2stk_ind('ptcl3D',inative,native_stkind,native_indstk)
            call registered_spproj%map_ptcl_ind2stk_ind('ptcl3D',ireg,registered_stkind,registered_indstk)
            if( native_stkind<1 .or. native_indstk<1 .or. registered_stkind<1 .or. registered_indstk<1 ) &
                &THROW_HARD('invalid stack address in registered-to-native flex particle mapping')
            if( write_manifest ) write(u,*) i,ireg,inative,native_stkind,native_indstk,registered_stkind, &
                &registered_indstk,registered_3d%e1get(),registered_3d%e2get(),native_3d%e3get(), &
                &native_3d%get_2Dshift(),mapped_3d%e1get(),mapped_3d%e2get(),mapped_3d%e3get(), &
                &mapped_3d%get_2Dshift()
            mapped_native(inative)=.true.
        end do
        if( write_manifest ) close(u)
        call outproj%write(native_project)
        write(logfhandle,'(A,I0,A,A)') '>>> FLEX REGISTERED-TO-NATIVE MAP VERIFIED: ',size(native_pinds), &
            &' mapping=ptcl3D_row project=',native_project%to_char()
        call native_3d%kill; call registered_3d%kill; call mapped_3d%kill
        call registered_spproj%kill; call outproj%kill
        deallocate(mapped_native)

    end subroutine map_flex_registered_to_native_project

    subroutine prepare_flex_diffmap_features( params, build, pinds, features, proj_ids, proj_dirs, &
        &registered_stack, registered_project, part, retain_features )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:)
        real, allocatable, intent(out)   :: features(:,:), proj_dirs(:,:)
        integer, allocatable, intent(out):: proj_ids(:)
        type(string), intent(out)        :: registered_stack, registered_project
        integer, optional, intent(in)    :: part
        logical, optional, intent(in)    :: retain_features
        type(sp_project) :: registered_spproj
        type(stack_io)   :: registered_io, residual_io
        type(image)      :: mean_vol
        type(image), allocatable :: mean_reprojs(:), registered(:), registered_crops(:), mean_ctfs(:)
        type(ctfparams)  :: ctfparms
        type(ctf)        :: tfun
        type(ori)        :: ptcl_ori
        real, allocatable :: vec(:)
        integer, allocatable :: batch_pinds(:), p2row(:), batches(:,:)
        logical, allocatable :: processed(:), used_proj(:)
        real :: e3
        integer :: nptcls, nall, nproj, npix, batchsz_max, nbatches
        type(string) :: residual_stack,map_fname
        integer :: i, j, iptcl, iproj, row, ithr, map_unit, ctf_mode, ibatch, batch_start, batch_end,part_here
        logical :: do_phaseflip,l_part,l_retain
        nptcls = size(pinds)
        if( nptcls < 1 ) THROW_HARD('flex_analysis feature preparation received no particles')
        l_part=present(part)
        l_retain=.true.
        if( present(retain_features) ) l_retain=retain_features
        part_here=0
        if( l_part ) part_here=part
        nall = build%spproj%os_ptcl3D%get_noris()
        if( nall < 1 ) THROW_HARD('flex_analysis project has no ptcl3D records')
        allocate(p2row(nall), source=0)
        allocate(processed(nall), source=.false.)
        nproj = build%eulspace%get_noris()
        if( nproj < 1 ) THROW_HARD('flex_analysis requires the builder projection-direction grid')
        allocate(proj_ids(nptcls),source=0)
        do i = 1,nptcls
            if( pinds(i) < 1 .or. pinds(i) > nall ) THROW_HARD('invalid active particle index in flex feature preparation')
            if( p2row(pinds(i)) /= 0 ) THROW_HARD('duplicate active particle index in flex feature preparation')
            p2row(pinds(i)) = i
            ! Stored proj values may have been assigned on a different
            ! nspace during refinement.  The graph must instead use bins on
            ! this run's explicit eulspace grid.
            call build%spproj_field%get_ori(pinds(i),ptcl_ori)
            proj_ids(i)=build%eulspace%find_closest_proj(ptcl_ori)
            if( proj_ids(i)<1 .or. proj_ids(i)>nproj ) &
                &THROW_HARD('active particle projection index outside explicit nspace grid')
        end do
        npix = params%box_crop * params%box_crop
        if( l_retain ) allocate(features(npix,nptcls),source=0.)
        if( l_part )then
            registered_stack = string('flex_registered_particles_part')//int2str_pad(part_here,params%numlen)//'.mrcs'
            residual_stack   = string('flex_residual_features_part')//int2str_pad(part_here,params%numlen)//'.mrcs'
            map_fname        = string('flex_registered_particle_map_part')//int2str_pad(part_here,params%numlen)//TXT_EXT
            registered_project = ''
        else
            registered_stack   = 'flex_registered_particles.mrcs'
            map_fname          = 'flex_registered_particle_map.txt'
            registered_project = 'flex_registered_particles.simple'
        endif
        call del_file(registered_stack)
        if( l_part ) call del_file(residual_stack)
        if( .not.l_part ) call del_file(registered_project)
        call del_file(map_fname)
        open(newunit=map_unit, file=map_fname%to_char(), status='replace', action='write')
        write(map_unit,'(A)') '# feature_row raw_particle_index projection_index registered_stack_index'

        if( l_part )then
            call read_mean_projection_stack(params,nproj,mean_reprojs)
        else
            call mean_vol%read_and_crop(params%vols(1),params%smpd,params%box_crop,params%smpd_crop)
            if( params%msk_crop>TINY ) call mean_vol%mask3D_soft(params%msk_crop,backgr=0.)
            mean_reprojs=reproject(mean_vol,build%eulspace,top=nproj)
        endif
        if( l_retain ) call flex_projection_directions(build,proj_dirs)
        allocate(used_proj(nproj),source=.false.)
        do i=1,nptcls
            used_proj(proj_ids(i))=.true.
        end do
        !$omp parallel do default(shared) schedule(static) proc_bind(close)
        do iproj=1,nproj
            if( used_proj(iproj) ) call mean_reprojs(iproj)%fft
        end do
        !$omp end parallel do
        allocate(mean_ctfs(max(1,params%nthr)),registered_crops(min(nptcls,MAXIMGBATCHSZ)))
        !$omp parallel do default(shared) schedule(static) proc_bind(close)
        do i=1,size(mean_ctfs)
            call mean_ctfs(i)%new([params%box_crop,params%box_crop,1],params%smpd_crop,wthreads=.false.)
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) schedule(static) proc_bind(close)
        do i=1,size(registered_crops)
            call registered_crops(i)%new([params%box_crop,params%box_crop,1],params%smpd_crop,wthreads=.false.)
        end do
        !$omp end parallel do
        call mean_ctfs(1)%memoize_mask_coords
        ctf_mode = build%spproj%get_ctfflag_type('ptcl3D',pinds(1))
        select case(ctf_mode)
            case(CTFFLAG_YES)
                do_phaseflip = .true.
            case(CTFFLAG_FLIP)
                do_phaseflip = .false.
            case DEFAULT
                THROW_HARD('flex_analysis requires ctf=yes or ctf=flip particle images')
        end select
        do i=2,nptcls
            if( build%spproj%get_ctfflag_type('ptcl3D',pinds(i)) /= ctf_mode )then
                THROW_HARD('flex_analysis does not support mixed ctf=yes and ctf=flip particle stacks')
            endif
        end do

        batchsz_max = min(nptcls, min(MAXIMGBATCHSZ,max(1,params%nthr) * BATCHTHRSZ))
        nbatches    = ceiling(real(nptcls) / real(batchsz_max))
        batches     = split_nobjs_even(nptcls, nbatches)
        call registered_io%open(registered_stack, params%smpd_crop, 'write', box=params%box_crop)
        if( l_part ) call residual_io%open(residual_stack, params%smpd_crop, 'write', box=params%box_crop)
        do ibatch = 1,nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            call transform_ptcls(params, build, build%spproj, 'ptcl3D', 0, registered, batch_pinds, &
                &phflip=do_phaseflip, pinds_in=pinds(batch_start:batch_end))
            if( .not. allocated(batch_pinds) .or. .not. allocated(registered) )then
                THROW_HARD('batch registration returned no particles in flex feature preparation')
            endif
            if( size(batch_pinds) /= batch_end - batch_start + 1 .or. size(registered) /= size(batch_pinds) )then
                THROW_HARD('batch registration size mismatch in flex feature preparation')
            endif
            do j = 1,size(batch_pinds)
                row   = batch_start + j - 1
                iptcl = batch_pinds(j)
                if( iptcl /= pinds(row) .or. p2row(iptcl) /= row )then
                    THROW_HARD('batch registration changed active particle order')
                endif
                if( processed(iptcl) ) THROW_HARD('active flex particle was registered twice')
            end do
            ! transform_ptcls uses and clears the global Fourier maps for its
            ! native and padded boxes, so restore the crop-box maps here.
            call memoize_ft_maps([params%box_crop,params%box_crop,1],params%smpd_crop)
            !$omp parallel do default(shared) schedule(static) proc_bind(close)
            do j = 1,size(batch_pinds)
                call registered(j)%fft
                call registered(j)%clip(registered_crops(j))
                call registered_crops(j)%ifft
            end do
            !$omp end parallel do
            call dealloc_imgarr(registered)
            do j = 1,size(batch_pinds)
                row   = batch_start + j - 1
                iptcl = batch_pinds(j)
                iproj = proj_ids(row)
                call registered_io%write(row,registered_crops(j))
                write(map_unit,'(I10,1X,I10,1X,I8,1X,I10)') row, iptcl, iproj, row
            end do
            !$omp parallel do default(shared) private(ithr,row,iptcl,iproj,ctfparms,e3,tfun,vec) &
            !$omp& schedule(static) proc_bind(close)
            do j = 1,size(batch_pinds)
                ithr  = omp_get_thread_num()+1
                row   = batch_start+j-1
                iptcl = batch_pinds(j)
                iproj = proj_ids(row)
                call mean_ctfs(ithr)%copy_fast(mean_reprojs(iproj))
                ctfparms = build%spproj%get_ctfparams('ptcl3D',iptcl)
                ! transform_ptcls rotates the image by -e3.  The astigmatism
                ! axes therefore rotate into the registered frame by +e3.
                e3 = build%spproj%os_ptcl3D%e3get(iptcl)
                ctfparms%angast = registered_angast(ctfparms%angast,e3)
                tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call mean_ctfs(ithr)%apply_ctf(tfun,'abs',ctfparms)
                call registered_crops(j)%fft
                call registered_crops(j)%subtr(mean_ctfs(ithr))
                if( params%lp > 2. * params%smpd_crop + TINY ) call registered_crops(j)%bp(0.,params%lp)
                call apply_sigma_shell_whitening(registered_crops(j),params,build,iptcl)
                call registered_crops(j)%ifft
                call registered_crops(j)%mask2D_softavg(params%msk_crop)
                if( l_retain )then
                    vec=registered_crops(j)%serialize()
                    if( size(vec)/=npix ) THROW_HARD('registered residual feature size mismatch')
                    features(:,row)=vec
                    deallocate(vec)
                endif
                processed(iptcl) = .true.
            end do
            !$omp end parallel do
            if( l_part )then
                do j=1,size(batch_pinds)
                    row=batch_start+j-1
                    call residual_io%write(row,registered_crops(j))
                end do
            endif
            call forget_ft_maps
            if( allocated(batch_pinds) ) deallocate(batch_pinds)
        end do
        call registered_io%close
        if( l_part ) call residual_io%close
        close(map_unit)
        if( any(proj_ids == 0) .or. any(.not. processed(pinds)) )then
            THROW_HARD('not every active particle received a registered residual feature')
        endif
        if( l_retain )then
            if( .not.all(ieee_is_finite(features)) ) THROW_HARD('nonfinite registered residual feature')
        endif
        if( .not.l_part )then
            call write_registered_project(params,build%spproj,pinds,proj_ids,registered_stack,registered_project,registered_spproj)
            call registered_spproj%kill
        endif
        call mean_vol%kill
        if( allocated(mean_reprojs) ) call dealloc_imgarr(mean_reprojs)
        call dealloc_imgarr(mean_ctfs)
        call dealloc_imgarr(registered_crops)
        call ptcl_ori%kill
        deallocate(p2row, processed, used_proj, batches)
        call residual_stack%kill; call map_fname%kill
    end subroutine prepare_flex_diffmap_features

    subroutine prepare_flex_diffmap_feature_part( params, build, pinds, part )
        class(parameters), intent(in) :: params
        class(builder), intent(inout) :: build
        integer, intent(in) :: pinds(:),part
        real, allocatable :: features(:,:),proj_dirs(:,:)
        integer, allocatable :: proj_ids(:)
        type(string) :: registered_stack,registered_project
        call prepare_flex_diffmap_features(params,build,pinds,features,proj_ids,proj_dirs,registered_stack, &
            &registered_project,part=part,retain_features=.false.)
        if( allocated(features) ) deallocate(features)
        deallocate(proj_ids)
        if( allocated(proj_dirs) ) deallocate(proj_dirs)
        call registered_stack%kill; call registered_project%kill
    end subroutine prepare_flex_diffmap_feature_part

    subroutine write_flex_mean_projection_stack( params, build )
        class(parameters), intent(in) :: params
        class(builder), intent(inout) :: build
        type(image) :: mean_vol
        type(image), allocatable :: mean_reprojs(:)
        type(stack_io) :: proj_io
        type(string) :: fname
        integer :: iproj,nproj,nproj_written,ldim(3)
        nproj=build%eulspace%get_noris()
        if( nproj<1 ) THROW_HARD('cannot prepare flex mean projections without eulspace')
        fname='flex_mean_projections.mrcs'; call del_file(fname)
        call mean_vol%read_and_crop(params%vols(1),params%smpd,params%box_crop,params%smpd_crop)
        if( params%msk_crop>TINY ) call mean_vol%mask3D_soft(params%msk_crop,backgr=0.)
        mean_reprojs=reproject(mean_vol,build%eulspace,top=nproj)
        call proj_io%open(fname,params%smpd_crop,'write',box=params%box_crop)
        do iproj=1,nproj
            call proj_io%write(iproj,mean_reprojs(iproj))
        end do
        call proj_io%close; call mean_vol%kill
        if( .not.file_exists(fname) ) THROW_HARD('flex mean projection stack was not written')
        call find_ldim_nptcls(fname,ldim,nproj_written)
        if( nproj_written/=nproj ) THROW_HARD('flex mean projection stack count mismatch after write')
        write(logfhandle,'(A,I0)') '>>> FLEX mean projection stack ready: projections=',nproj
        call flush(logfhandle)
        call dealloc_imgarr(mean_reprojs); call fname%kill
    end subroutine write_flex_mean_projection_stack

    subroutine read_mean_projection_stack( params, nproj, mean_reprojs )
        class(parameters), intent(in) :: params
        integer, intent(in) :: nproj
        type(image), allocatable, intent(out) :: mean_reprojs(:)
        type(stack_io) :: proj_io
        type(string) :: fname
        integer :: iproj
        fname='flex_mean_projections.mrcs'
        allocate(mean_reprojs(nproj))
        call proj_io%open(fname,params%smpd_crop,'read',bufsz=min(nproj,BUFSZ_DEFAULT))
        if( proj_io%get_nptcls()/=nproj ) THROW_HARD('flex mean projection stack count mismatch')
        do iproj=1,nproj
            call mean_reprojs(iproj)%new([params%box_crop,params%box_crop,1],params%smpd_crop,wthreads=.false.)
            call proj_io%read(iproj,mean_reprojs(iproj))
        end do
        call proj_io%close; call fname%kill
    end subroutine read_mean_projection_stack

    subroutine flex_projection_directions( build, proj_dirs )
        class(builder), intent(inout) :: build
        real, allocatable, intent(out) :: proj_dirs(:,:)
        type(ori) :: o
        real :: rotmat(3,3)
        integer :: iproj,nproj
        nproj=build%eulspace%get_noris()
        if( nproj<1 ) THROW_HARD('flex_analysis requires the builder projection-direction grid')
        allocate(proj_dirs(3,nproj))
        do iproj=1,nproj
            call build%eulspace%get_ori(iproj,o)
            call o%e3set(0.)
            rotmat=o%get_mat()
            proj_dirs(:,iproj)=rotmat(3,:)
        end do
        call o%kill
    end subroutine flex_projection_directions

    subroutine read_flex_diffmap_feature_parts( params, nparts, features, proj_ids, build, pinds )
        class(parameters), intent(in) :: params
        integer, intent(in) :: nparts
        real, allocatable, intent(out) :: features(:,:)
        integer, allocatable, intent(out) :: proj_ids(:)
        class(builder), optional, intent(inout) :: build
        integer, optional, intent(in) :: pinds(:)
        type(stack_io) :: residual_io
        type(image) :: residual
        type(ori) :: ptcl_ori
        type(string) :: stack_fname,map_fname
        real, allocatable :: rmat(:,:,:)
        integer, allocatable :: counts(:)
        integer :: ipart,nlocal,total,row,local_row,u,ios,map_local,pind,iproj,stack_ind,ldim(3)
        integer :: graph_step,graph_side,ix,iy,ixhi,iyhi,ifeat
        character(len=XLONGSTRLEN) :: line
        integer, parameter :: GRAPH_FEATURE_MAX_BOX=96
        integer, parameter :: GRAPH_FEATURE_BUFSZ=8
        allocate(counts(nparts),source=0)
        total=0
        do ipart=1,nparts
            stack_fname=string('flex_residual_features_part')//int2str_pad(ipart,params%numlen)//'.mrcs'
            call find_ldim_nptcls(stack_fname,ldim,nlocal)
            if( any(ldim(1:2)/=params%box_crop) ) THROW_HARD('flex residual feature part dimensions mismatch')
            counts(ipart)=nlocal; total=total+nlocal
            call stack_fname%kill
        end do
        if( total<1 ) THROW_HARD('distributed flex feature preparation produced no residuals')
        if( present(build) .neqv. present(pinds) ) THROW_HARD('flex graph particle metadata must be supplied together')
        if( present(pinds) ) then
            if( size(pinds)/=total ) THROW_HARD('flex graph particle metadata size mismatch')
        endif
        ! A graph worker needs every particle as a candidate, but retaining
        ! every crop-box pixel makes the distributed feature table several
        ! gigabytes.  Average-pool the already low-pass residuals onto a
        ! bounded grid before assembling them.  This is an anti-aliased
        ! compact graph descriptor; particle rows and angular gating remain
        ! unchanged.
        graph_step=max(1,(params%box_crop+GRAPH_FEATURE_MAX_BOX-1)/GRAPH_FEATURE_MAX_BOX)
        graph_side=(params%box_crop+graph_step-1)/graph_step
        allocate(features(graph_side*graph_side,total),source=0.)
        allocate(proj_ids(total),source=0)
        call residual%new([params%box_crop,params%box_crop,1],params%smpd_crop,wthreads=.false.)
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX GRAPH compact features raw_box=',params%box_crop, &
            &' compact_box=',graph_side,' values_per_particle=',graph_side*graph_side
        call flush(logfhandle)
        row=0
        do ipart=1,nparts
            stack_fname=string('flex_residual_features_part')//int2str_pad(ipart,params%numlen)//'.mrcs'
            map_fname=string('flex_registered_particle_map_part')//int2str_pad(ipart,params%numlen)//TXT_EXT
            ! The graph descriptor is assembled one image at a time; a large
            ! stack buffer only duplicates hundreds of megabytes alongside
            ! the compact global feature table.
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX GRAPH loading feature part=',ipart,' rows=',counts(ipart)
            call flush(logfhandle)
            call residual_io%open(stack_fname,params%smpd_crop,'read',bufsz=min(counts(ipart),GRAPH_FEATURE_BUFSZ))
            write(logfhandle,'(A,I0)') '>>> FLEX GRAPH feature stack open part=',ipart
            call flush(logfhandle)
            local_row=0
            if( present(build) ) then
                do while( local_row<counts(ipart) )
                    local_row=local_row+1
                    row=row+1
                    if( local_row==1 ) then
                        write(logfhandle,'(A)') '>>> FLEX GRAPH feature first-row orientation lookup'
                        call flush(logfhandle)
                    endif
                    call build%spproj_field%get_ori(pinds(row),ptcl_ori)
                    if( local_row==1 ) then
                        write(logfhandle,'(A)') '>>> FLEX GRAPH feature first-row Euler-bin lookup'
                        call flush(logfhandle)
                    endif
                    iproj=build%eulspace%find_closest_proj(ptcl_ori)
                    if( iproj<1 .or. iproj>build%eulspace%get_noris() ) &
                        &THROW_HARD('flex graph projection index outside explicit eulspace grid')
                    if( local_row==1 ) then
                        write(logfhandle,'(A)') '>>> FLEX GRAPH feature first-row residual read'
                        call flush(logfhandle)
                    endif
                    call residual_io%read(local_row,residual)
                    if( local_row==1 ) then
                        write(logfhandle,'(A)') '>>> FLEX GRAPH feature first-row residual copied'
                        call flush(logfhandle)
                    endif
                    rmat=residual%get_rmat()
                    if( any(shape(rmat)/=[params%box_crop,params%box_crop,1]) ) &
                        &THROW_HARD('flex residual feature image dimensions mismatch')
                    ifeat=0
                    do iy=1,params%box_crop,graph_step
                        iyhi=min(params%box_crop,iy+graph_step-1)
                        do ix=1,params%box_crop,graph_step
                            ixhi=min(params%box_crop,ix+graph_step-1)
                            ifeat=ifeat+1
                            features(ifeat,row)=sum(rmat(ix:ixhi,iy:iyhi,1))/real((ixhi-ix+1)*(iyhi-iy+1))
                        end do
                    end do
                    if( ifeat/=size(features,1) ) THROW_HARD('invalid compact flex feature size')
                    proj_ids(row)=iproj
                    deallocate(rmat)
                    if( mod(local_row,1024)==0 ) then
                        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX GRAPH feature progress part=',ipart, &
                            &' local_row=',local_row,' global_row=',row
                        call flush(logfhandle)
                    endif
                end do
            else
                open(newunit=u,file=map_fname%to_char(),status='old',action='read',iostat=ios)
                call fileiochk('opening flex feature part map '//map_fname%to_char(),ios)
                do
                    read(u,'(A)',iostat=ios) line
                    if( ios/=0 ) exit
                    if( len_trim(line)==0 .or. line(1:1)=='#' ) cycle
                    read(line,*) map_local,pind,iproj,stack_ind
                local_row=local_row+1; row=row+1
                if( map_local/=local_row .or. stack_ind/=local_row ) THROW_HARD('invalid flex feature part map ordering')
                call residual_io%read(local_row,residual)
                rmat=residual%get_rmat()
                if( any(shape(rmat)/=[params%box_crop,params%box_crop,1]) ) &
                    &THROW_HARD('flex residual feature image dimensions mismatch')
                ifeat=0
                do iy=1,params%box_crop,graph_step
                    iyhi=min(params%box_crop,iy+graph_step-1)
                    do ix=1,params%box_crop,graph_step
                        ixhi=min(params%box_crop,ix+graph_step-1)
                        ifeat=ifeat+1
                        features(ifeat,row)=sum(rmat(ix:ixhi,iy:iyhi,1))/real((ixhi-ix+1)*(iyhi-iy+1))
                    end do
                end do
                if( ifeat/=size(features,1) ) THROW_HARD('invalid compact flex feature size')
                proj_ids(row)=iproj
                deallocate(rmat)
                if( mod(local_row,1024)==0 ) then
                    write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX GRAPH feature progress part=',ipart, &
                        &' local_row=',local_row,' global_row=',row
                    call flush(logfhandle)
                endif
                end do
                close(u)
            endif
            call residual_io%close
            if( local_row/=counts(ipart) ) THROW_HARD('flex feature part map/stack count mismatch')
            write(logfhandle,'(A,I0)') '>>> FLEX GRAPH feature part complete=',ipart
            call flush(logfhandle)
            call stack_fname%kill; call map_fname%kill
        end do
        call residual%kill
        if( row/=total .or. any(proj_ids<1) ) THROW_HARD('incomplete distributed flex feature assembly')
        if( .not.all(ieee_is_finite(features)) ) THROW_HARD('nonfinite distributed flex feature')
        deallocate(counts)
    end subroutine read_flex_diffmap_feature_parts

    subroutine assemble_flex_diffmap_feature_parts( params, spproj, pinds, nparts, parts, project_fname )
        class(parameters), intent(in) :: params
        type(sp_project), intent(inout) :: spproj
        integer, intent(in) :: pinds(:),nparts,parts(nparts,2)
        type(string), intent(out) :: project_fname
        type(sp_project) :: outproj
        type(string) :: stack_fname,frcs_fname,map_fname
        logical, allocatable :: active(:)
        integer, allocatable :: part_proj_ids(:)
        character(len=XLONGSTRLEN) :: line
        real :: e3,angast
        integer :: ipart,i,iptcl,nptcls,nlocal,ldim(3),stkind,indstk,map_unit
        integer :: map_read_unit,ios,map_local,map_pind,map_proj,map_stack,nread
        project_fname='flex_registered_particles.simple'
        call del_file(project_fname)
        nptcls=spproj%os_ptcl3D%get_noris()
        allocate(active(nptcls),source=.false.); active(pinds)=.true.
        allocate(part_proj_ids(size(pinds)),source=0)
        do ipart=1,nparts
            map_fname=string('flex_registered_particle_map_part')//int2str_pad(ipart,params%numlen)//TXT_EXT
            open(newunit=map_read_unit,file=map_fname%to_char(),status='old',action='read',iostat=ios)
            call fileiochk('opening registered flex feature map '//map_fname%to_char(),ios)
            nread=0
            do
                read(map_read_unit,'(A)',iostat=ios) line
                if( ios/=0 ) exit
                if( len_trim(line)==0 .or. line(1:1)=='#' ) cycle
                read(line,*) map_local,map_pind,map_proj,map_stack
                nread=nread+1
                i=parts(ipart,1)+nread-1
                if( i>parts(ipart,2) .or. map_local/=nread .or. map_stack/=nread ) &
                    &THROW_HARD('invalid registered flex feature-map ordering')
                if( map_pind/=pinds(i) .or. map_proj<1 .or. map_proj>params%nspace ) &
                    &THROW_HARD('registered flex feature-map content mismatch')
                part_proj_ids(i)=map_proj
            end do
            close(map_read_unit); call map_fname%kill
            if( nread/=parts(ipart,2)-parts(ipart,1)+1 ) THROW_HARD('incomplete registered flex feature map')
        end do
        if( any(part_proj_ids<1) ) THROW_HARD('distributed registered flex project has invalid projections')
        call outproj%copy(spproj); call outproj%update_projinfo(project_fname)
        if( outproj%os_stk%get_noris()>0 ) call outproj%os_stk%kill()
        call outproj%os_stk%new(nparts,is_ptcl=.false.)
        do ipart=1,nparts
            stack_fname=string('flex_registered_particles_part')//int2str_pad(ipart,params%numlen)//'.mrcs'
            call find_ldim_nptcls(stack_fname,ldim,nlocal)
            if( nlocal/=parts(ipart,2)-parts(ipart,1)+1 .or. any(ldim(1:2)/=params%box_crop) ) &
                &THROW_HARD('registered flex feature part stack mismatch')
            call outproj%os_stk%set(ipart,'stk',simple_abspath(stack_fname))
            call outproj%os_stk%set(ipart,'box',params%box_crop)
            call outproj%os_stk%set(ipart,'nptcls',nlocal)
            call outproj%os_stk%set(ipart,'nptcls_stk',nlocal)
            call outproj%os_stk%set(ipart,'fromp',parts(ipart,1))
            call outproj%os_stk%set(ipart,'top',parts(ipart,2))
            call outproj%os_stk%set(ipart,'stkkind','single')
            call outproj%os_stk%set(ipart,'imgkind','ptcl')
            call outproj%os_stk%set(ipart,'smpd',params%smpd_crop)
            call outproj%os_stk%set(ipart,'kv',params%kv); call outproj%os_stk%set(ipart,'cs',params%cs)
            call outproj%os_stk%set(ipart,'fraca',params%fraca); call outproj%os_stk%set(ipart,'phshift',0.)
            call outproj%os_stk%set(ipart,'ctf','flip'); call outproj%os_stk%set(ipart,'state',1)
            call stack_fname%kill
        end do
        do iptcl=1,nptcls
            call outproj%os_ptcl3D%set(iptcl,'stkind',1); call outproj%os_ptcl3D%set(iptcl,'indstk',1)
            if( .not.active(iptcl) ) call outproj%os_ptcl3D%set(iptcl,'state',0)
        end do
        do ipart=1,nparts
            do i=parts(ipart,1),parts(ipart,2)
                iptcl=pinds(i)
                call outproj%os_ptcl3D%set(iptcl,'stkind',ipart)
                call outproj%os_ptcl3D%set(iptcl,'indstk',i-parts(ipart,1)+1)
                call outproj%os_ptcl3D%set(iptcl,'proj',part_proj_ids(i))
                if( outproj%os_ptcl3D%isthere(iptcl,'angast') )then
                    e3=spproj%os_ptcl3D%e3get(iptcl); angast=spproj%os_ptcl3D%get(iptcl,'angast')
                    call outproj%os_ptcl3D%set(iptcl,'angast',registered_angast(angast,e3))
                endif
            end do
        end do
        call outproj%os_ptcl3D%zero_inpl()
        do i=1,size(pinds)
            call outproj%map_ptcl_ind2stk_ind('ptcl3D',pinds(i),stkind,indstk)
            if( stkind<1 .or. indstk<1 ) THROW_HARD('distributed registered project mapping mismatch')
        end do
        if( outproj%os_out%get_noris()>0 ) call outproj%os_out%kill()
        call spproj%get_frcs(frcs_fname,'frc3D',fail=.false.)
        if( file_exists(frcs_fname) ) call outproj%add_frcs2os_out(frcs_fname,'frc3D')
        if( frcs_fname%is_allocated() ) call frcs_fname%kill
        call outproj%write(project_fname); call outproj%kill
        call del_file('flex_registered_particle_map.txt')
        open(newunit=map_unit,file='flex_registered_particle_map.txt',status='replace',action='write')
        write(map_unit,'(A)') '# feature_row raw_particle_index projection_index registered_stack_index stack_part'
        do ipart=1,nparts
            do i=parts(ipart,1),parts(ipart,2)
                iptcl=pinds(i)
                write(map_unit,'(I10,1X,I10,1X,I8,1X,I10,1X,I8)') i,iptcl, &
                    &part_proj_ids(i),i-parts(ipart,1)+1,ipart
            end do
        end do
        close(map_unit)
        deallocate(active,part_proj_ids)
    end subroutine assemble_flex_diffmap_feature_parts

    subroutine write_registered_project( params, spproj, pinds, proj_ids, stack_fname, project_fname, outproj )
        class(parameters), intent(in)    :: params
        type(sp_project),  intent(inout) :: spproj
        integer,           intent(in)    :: pinds(:),proj_ids(:)
        type(string),      intent(in)    :: stack_fname, project_fname
        type(sp_project),  intent(inout) :: outproj
        type(string) :: frcs_fname
        logical, allocatable :: active(:)
        real :: e3, angast
        integer :: iptcl, i, nptcls, nactive, nstack, ldim_stack(3), stkind, indstk
        nptcls = spproj%os_ptcl3D%get_noris()
        nactive = size(pinds)
        if( size(proj_ids)/=nactive ) THROW_HARD('registered flex project projection-map size mismatch')
        allocate(active(nptcls), source=.false.)
        active(pinds) = .true.
        call find_ldim_nptcls(stack_fname, ldim_stack, nstack)
        if( nstack /= nactive ) THROW_HARD('registered flex stack count does not match active particle count')
        if( any(ldim_stack(1:2) /= params%box_crop) ) THROW_HARD('registered flex stack has unexpected dimensions')
        call outproj%copy(spproj)
        call outproj%update_projinfo(project_fname)
        if( outproj%os_stk%get_noris() > 0 ) call outproj%os_stk%kill()
        call outproj%os_stk%new(1, is_ptcl=.false.)
        call outproj%os_stk%set(1, 'stk',        simple_abspath(stack_fname))
        call outproj%os_stk%set(1, 'box',        params%box_crop)
        call outproj%os_stk%set(1, 'nptcls',     nactive)
        call outproj%os_stk%set(1, 'nptcls_stk', nactive)
        call outproj%os_stk%set(1, 'fromp',      1)
        call outproj%os_stk%set(1, 'top',        nactive)
        call outproj%os_stk%set(1, 'stkkind',    'single')
        call outproj%os_stk%set(1, 'imgkind',    'ptcl')
        call outproj%os_stk%set(1, 'smpd',       params%smpd_crop)
        call outproj%os_stk%set(1, 'kv',         params%kv)
        call outproj%os_stk%set(1, 'cs',         params%cs)
        call outproj%os_stk%set(1, 'fraca',      params%fraca)
        call outproj%os_stk%set(1, 'phshift',    0.)
        call outproj%os_stk%set(1, 'ctf',        'flip')
        call outproj%os_stk%set(1, 'state',      1)
        do iptcl = 1,nptcls
            call outproj%os_ptcl3D%set(iptcl,'stkind',1)
            call outproj%os_ptcl3D%set(iptcl,'indstk',1)
            if( .not. active(iptcl) ) call outproj%os_ptcl3D%set(iptcl,'state',0)
        end do
        do i = 1,nactive
            call outproj%os_ptcl3D%set(pinds(i),'indstk',i)
            call outproj%os_ptcl3D%set(pinds(i),'proj',proj_ids(i))
            if( outproj%os_ptcl3D%isthere(pinds(i),'angast') )then
                e3     = spproj%os_ptcl3D%e3get(pinds(i))
                angast = spproj%os_ptcl3D%get(pinds(i),'angast')
                call outproj%os_ptcl3D%set(pinds(i),'angast',registered_angast(angast,e3))
            endif
        end do
        call outproj%os_ptcl3D%zero_inpl()
        do i=1,nactive
            call outproj%map_ptcl_ind2stk_ind('ptcl3D',pinds(i),stkind,indstk)
            if( stkind/=1 .or. indstk/=i ) THROW_HARD('registered ptcl3D stack mapping mismatch')
        end do
        if( outproj%os_out%get_noris() > 0 ) call outproj%os_out%kill()
        call spproj%get_frcs(frcs_fname, 'frc3D', fail=.false.)
        if( file_exists(frcs_fname%to_char()) ) call outproj%add_frcs2os_out(frcs_fname, 'frc3D')
        if( frcs_fname%is_allocated() ) call frcs_fname%kill
        call outproj%write(project_fname)
        deallocate(active)
    end subroutine write_registered_project

    subroutine apply_sigma_shell_whitening( residual_img, params, build, iptcl )
        class(image),      intent(inout) :: residual_img
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: iptcl
        real, allocatable :: sigma_shells(:), sigma_filter(:)
        integer :: kfromto(2), nyq, sigma_nyq, sh
        real :: sigma_floor
        if( .not. allocated(build%esig%sigma2_noise) ) return
        if( iptcl < lbound(build%esig%sigma2_noise,2) .or. iptcl > ubound(build%esig%sigma2_noise,2) ) return
        kfromto = build%esig%get_kfromto()
        nyq = residual_img%get_nyq()
        if( nyq < 1 ) return
        sigma_nyq = fdim(params%box_crop) - 1
        allocate(sigma_shells(0:nyq), sigma_filter(nyq))
        call upsample_sigma2(kfromto(1), sigma_nyq, build%esig%sigma2_noise(kfromto(1):kfromto(2),iptcl), nyq, sigma_shells)
        sigma_floor = max(TINY, minval(sigma_shells(1:nyq), mask=sigma_shells(1:nyq) > TINY))
        do sh=1,nyq
            sigma_filter(sh) = 1. / sqrt(max(sigma_shells(sh), sigma_floor))
        end do
        call residual_img%apply_filter(sigma_filter)
        deallocate(sigma_shells, sigma_filter)
    end subroutine apply_sigma_shell_whitening

    pure elemental real function registered_angast( angast, e3 ) result(canonical_angast)
        real, intent(in) :: angast,e3
        canonical_angast=modulo(angast+e3,180.)
    end function registered_angast

end module simple_flex_diffmap_features
