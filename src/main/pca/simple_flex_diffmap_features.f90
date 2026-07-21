!@descr: registered residual feature preparation for flex_eigenvol diffusion maps
module simple_flex_diffmap_features
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,          only: builder
use simple_classaverager,    only: transform_ptcls
use simple_ctf,              only: ctf
use simple_image,            only: image
use simple_imgarr_utils,     only: dealloc_imgarr
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

contains

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
        type(image)      :: mean_vol, mean_ctf, reg_crop, residual
        type(image), allocatable :: mean_reprojs(:), registered(:)
        type(ctfparams)  :: ctfparms
        type(ctf)        :: tfun
        real, allocatable :: vec(:)
        integer, allocatable :: batch_pinds(:), p2row(:), batches(:,:)
        logical, allocatable :: processed(:)
        real :: e3
        integer :: nptcls, nall, nproj, npix, batchsz_max, nbatches
        type(string) :: residual_stack,map_fname
        integer :: i, j, iptcl, iproj, row, map_unit, ctf_mode, ibatch, batch_start, batch_end,part_here
        logical :: do_phaseflip,l_part,l_retain
        nptcls = size(pinds)
        if( nptcls < 1 ) THROW_HARD('flex_eigenvol feature preparation received no particles')
        l_part=present(part)
        l_retain=.true.
        if( present(retain_features) ) l_retain=retain_features
        part_here=0
        if( l_part ) part_here=part
        nall = build%spproj%os_ptcl3D%get_noris()
        if( nall < 1 ) THROW_HARD('flex_eigenvol project has no ptcl3D records')
        allocate(p2row(nall), source=0)
        allocate(processed(nall), source=.false.)
        nproj = build%eulspace%get_noris()
        if( nproj < 1 ) THROW_HARD('flex_eigenvol requires the builder projection-direction grid')
        allocate(proj_ids(nptcls),source=0)
        do i = 1,nptcls
            if( pinds(i) < 1 .or. pinds(i) > nall ) THROW_HARD('invalid active particle index in flex feature preparation')
            if( p2row(pinds(i)) /= 0 ) THROW_HARD('duplicate active particle index in flex feature preparation')
            p2row(pinds(i)) = i
            proj_ids(i)=build%spproj%os_ptcl3D%get_int(pinds(i),'proj')
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
        call mean_ctf%new([params%box_crop,params%box_crop,1], params%smpd_crop, wthreads=.false.)
        call reg_crop%new([params%box_crop,params%box_crop,1], params%smpd_crop, wthreads=.false.)
        call residual%new([params%box_crop,params%box_crop,1], params%smpd_crop, wthreads=.false.)
        ctf_mode = build%spproj%get_ctfflag_type('ptcl3D',pinds(1))
        select case(ctf_mode)
            case(CTFFLAG_YES)
                do_phaseflip = .true.
            case(CTFFLAG_FLIP)
                do_phaseflip = .false.
            case DEFAULT
                THROW_HARD('flex_eigenvol requires ctf=yes or ctf=flip particle images')
        end select
        do i=2,nptcls
            if( build%spproj%get_ctfflag_type('ptcl3D',pinds(i)) /= ctf_mode )then
                THROW_HARD('flex_eigenvol does not support mixed ctf=yes and ctf=flip particle stacks')
            endif
        end do

        batchsz_max = min(nptcls, max(1,params%nthr) * BATCHTHRSZ)
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
            call memoize_ft_maps([params%box_crop,params%box_crop,1],params%smpd_crop)
            do j = 1,size(batch_pinds)
                row   = batch_start + j - 1
                iptcl = batch_pinds(j)
                if( iptcl /= pinds(row) .or. p2row(iptcl) /= row )then
                    THROW_HARD('batch registration changed active particle order')
                endif
                if( processed(iptcl) ) THROW_HARD('active flex particle was registered twice')
                iproj = proj_ids(row)
                call registered(j)%fft
                call registered(j)%clip(reg_crop)
                call reg_crop%ifft
                call registered_io%write(row, reg_crop)
                write(map_unit,'(I10,1X,I10,1X,I8,1X,I10)') row, iptcl, iproj, row
                call mean_ctf%copy_fast(mean_reprojs(iproj))
                call mean_ctf%fft
                ctfparms = build%spproj%get_ctfparams('ptcl3D',iptcl)
                ! transform_ptcls rotates the image by -e3.  The astigmatism
                ! axes therefore rotate into the registered frame by +e3.
                e3 = build%spproj%os_ptcl3D%e3get(iptcl)
                ctfparms%angast = registered_angast(ctfparms%angast,e3)
                tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call mean_ctf%apply_ctf(tfun, 'abs', ctfparms)
                call reg_crop%fft
                residual = reg_crop - mean_ctf
                call residual%ifft
                if( params%lp > 2. * params%smpd_crop + TINY ) call residual%bp(0., params%lp)
                call residual%mask2D_softavg(params%msk_crop)
                if( l_part ) call residual_io%write(row,residual)
                if( l_retain )then
                    vec=residual%serialize()
                    if( size(vec)/=npix ) THROW_HARD('registered residual feature size mismatch')
                    features(:,row)=vec
                    deallocate(vec)
                endif
                processed(iptcl) = .true.
            end do
            call forget_ft_maps
            if( allocated(registered) ) call dealloc_imgarr(registered)
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
        call mean_ctf%kill
        call reg_crop%kill
        call residual%kill
        deallocate(p2row, processed, batches)
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
        integer :: iproj,nproj
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
        if( nproj<1 ) THROW_HARD('flex_eigenvol requires the builder projection-direction grid')
        allocate(proj_dirs(3,nproj))
        do iproj=1,nproj
            call build%eulspace%get_ori(iproj,o)
            call o%e3set(0.)
            rotmat=o%get_mat()
            proj_dirs(:,iproj)=rotmat(3,:)
        end do
        call o%kill
    end subroutine flex_projection_directions

    subroutine read_flex_diffmap_feature_parts( params, nparts, features, proj_ids )
        class(parameters), intent(in) :: params
        integer, intent(in) :: nparts
        real, allocatable, intent(out) :: features(:,:)
        integer, allocatable, intent(out) :: proj_ids(:)
        type(stack_io) :: residual_io
        type(image) :: residual
        type(string) :: stack_fname,map_fname
        real, allocatable :: vec(:)
        integer, allocatable :: counts(:)
        integer :: ipart,nlocal,total,row,local_row,u,ios,map_local,pind,iproj,stack_ind,ldim(3)
        character(len=XLONGSTRLEN) :: line
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
        allocate(features(params%box_crop*params%box_crop,total),source=0.)
        allocate(proj_ids(total),source=0)
        row=0
        do ipart=1,nparts
            stack_fname=string('flex_residual_features_part')//int2str_pad(ipart,params%numlen)//'.mrcs'
            map_fname=string('flex_registered_particle_map_part')//int2str_pad(ipart,params%numlen)//TXT_EXT
            call residual_io%open(stack_fname,params%smpd_crop,'read',bufsz=min(counts(ipart),BUFSZ_DEFAULT))
            open(newunit=u,file=map_fname%to_char(),status='old',action='read',iostat=ios)
            call fileiochk('opening flex feature part map '//map_fname%to_char(),ios)
            local_row=0
            do
                read(u,'(A)',iostat=ios) line
                if( ios/=0 ) exit
                if( len_trim(line)==0 .or. line(1:1)=='#' ) cycle
                read(line,*) map_local,pind,iproj,stack_ind
                local_row=local_row+1; row=row+1
                if( map_local/=local_row .or. stack_ind/=local_row ) THROW_HARD('invalid flex feature part map ordering')
                call residual_io%read(local_row,residual)
                vec=residual%serialize()
                if( size(vec)/=size(features,1) ) THROW_HARD('flex residual feature size mismatch')
                features(:,row)=vec; proj_ids(row)=iproj
                deallocate(vec)
            end do
            close(u); call residual_io%close
            if( local_row/=counts(ipart) ) THROW_HARD('flex feature part map/stack count mismatch')
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
        if( outproj%os_ptcl2D%get_noris()==nptcls )then
            do iptcl=1,nptcls
                call outproj%os_ptcl2D%set(iptcl,'stkind',1); call outproj%os_ptcl2D%set(iptcl,'indstk',1)
                if( .not.active(iptcl) ) call outproj%os_ptcl2D%set(iptcl,'state',0)
            end do
            do ipart=1,nparts
                do i=parts(ipart,1),parts(ipart,2)
                    iptcl=pinds(i)
                    call outproj%os_ptcl2D%set(iptcl,'stkind',ipart)
                    call outproj%os_ptcl2D%set(iptcl,'indstk',i-parts(ipart,1)+1)
                    if( outproj%os_ptcl2D%isthere(iptcl,'angast') )then
                        e3=spproj%os_ptcl3D%e3get(iptcl); angast=spproj%os_ptcl2D%get(iptcl,'angast')
                        call outproj%os_ptcl2D%set(iptcl,'angast',registered_angast(angast,e3))
                    endif
                end do
            end do
            call outproj%os_ptcl2D%zero_inpl()
        endif
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
        if( outproj%os_ptcl2D%get_noris() == nptcls )then
            do iptcl = 1,nptcls
                call outproj%os_ptcl2D%set(iptcl,'stkind',1)
                call outproj%os_ptcl2D%set(iptcl,'indstk',1)
                if( .not. active(iptcl) ) call outproj%os_ptcl2D%set(iptcl,'state',0)
            end do
            do i = 1,nactive
                call outproj%os_ptcl2D%set(pinds(i),'indstk',i)
                if( outproj%os_ptcl2D%isthere(pinds(i),'angast') )then
                    e3     = spproj%os_ptcl3D%e3get(pinds(i))
                    angast = spproj%os_ptcl2D%get(pinds(i),'angast')
                    call outproj%os_ptcl2D%set(pinds(i),'angast',registered_angast(angast,e3))
                endif
            end do
            call outproj%os_ptcl2D%zero_inpl()
        endif
        do i=1,nactive
            call outproj%map_ptcl_ind2stk_ind('ptcl3D',pinds(i),stkind,indstk)
            if( stkind/=1 .or. indstk/=i ) THROW_HARD('registered ptcl3D stack mapping mismatch')
            if( outproj%os_ptcl2D%get_noris()==nptcls )then
                call outproj%map_ptcl_ind2stk_ind('ptcl2D',pinds(i),stkind,indstk)
                if( stkind/=1 .or. indstk/=i ) THROW_HARD('registered ptcl2D stack mapping mismatch')
            endif
        end do
        if( outproj%os_out%get_noris() > 0 ) call outproj%os_out%kill()
        call spproj%get_frcs(frcs_fname, 'frc3D', fail=.false.)
        if( file_exists(frcs_fname%to_char()) ) call outproj%add_frcs2os_out(frcs_fname, 'frc3D')
        if( frcs_fname%is_allocated() ) call frcs_fname%kill
        call outproj%write(project_fname)
        deallocate(active)
    end subroutine write_registered_project

    pure elemental real function registered_angast( angast, e3 ) result(canonical_angast)
        real, intent(in) :: angast,e3
        canonical_angast=modulo(angast+e3,180.)
    end function registered_angast

end module simple_flex_diffmap_features
