!@descr: registered residual feature preparation for flex_eigenvol diffusion maps
module simple_flex_diffmap_features
use, intrinsic :: iso_c_binding, only: c_float
use simple_core_module_api
use simple_builder,          only: builder
use simple_classaverager,    only: transform_ptcls
use simple_ctf,              only: ctf
use simple_image,            only: image
use simple_imgfile,          only: imgfile
use simple_imgarr_utils,     only: dealloc_imgarr
use simple_memoize_ft_maps,  only: memoize_ft_maps, forget_ft_maps
use simple_ori,              only: ori
use simple_parameters,       only: parameters
use simple_projector,        only: projector
use simple_sp_project,       only: sp_project
implicit none
private
#include "simple_local_flags.inc"

public :: prepare_flex_diffmap_features

contains

    subroutine prepare_flex_diffmap_features( params, build, pinds, features, proj_ids, proj_dirs, &
        &registered_stack, registered_project )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:)
        real, allocatable, intent(out)   :: features(:,:), proj_dirs(:,:)
        integer, allocatable, intent(out):: proj_ids(:)
        type(string), intent(out)        :: registered_stack, registered_project
        type(sp_project) :: registered_spproj
        type(projector)  :: mean_vol_pad
        type(image)      :: mean_vol, proj_pad, mean_img, mean_ctf, reg_crop, residual
        type(image), allocatable :: registered(:)
        type(ori)        :: o
        type(ctfparams)  :: ctfparms
        type(ctf)        :: tfun
        real, allocatable :: vec(:)
        integer, allocatable :: group_pinds(:), p2row(:)
        real :: rotmat(3,3)
        integer :: nptcls, nall, nproj, npix, boxpd, ldim_pd(3)
        integer :: i, j, iptcl, iproj, row, map_unit, ctf_mode
        logical :: do_phaseflip
        nptcls = size(pinds)
        if( nptcls < 3 ) THROW_HARD('flex_eigenvol requires at least three active particles')
        nall = build%spproj%os_ptcl3D%get_noris()
        if( nall < 1 ) THROW_HARD('flex_eigenvol project has no ptcl3D records')
        allocate(p2row(nall), source=0)
        nproj = 0
        do i = 1,nptcls
            if( pinds(i) < 1 .or. pinds(i) > nall ) THROW_HARD('invalid active particle index in flex feature preparation')
            p2row(pinds(i)) = i
            nproj = max(nproj, build%spproj%os_ptcl3D%get_int(pinds(i),'proj'))
        end do
        if( nproj < 1 ) THROW_HARD('active flex particles require positive proj indices')
        npix = params%box_crop * params%box_crop
        allocate(features(npix,nptcls), source=0.)
        allocate(proj_ids(nptcls), source=0)
        allocate(proj_dirs(3,nproj), source=0.)
        registered_stack   = string('flex_registered_particles.mrcs')
        registered_project = string('flex_registered_particles.simple')
        call del_file(registered_stack)
        call del_file(registered_project)
        call del_file('flex_registered_particle_map.txt')
        open(newunit=map_unit, file='flex_registered_particle_map.txt', status='replace', action='write')
        write(map_unit,'(A)') '# feature_row raw_particle_index projection_index'

        call mean_vol%read_and_crop(params%vols(1), params%smpd, params%box_crop, params%smpd_crop)
        boxpd   = 2 * round2even(KBALPHA * real(params%box_crop / 2))
        ldim_pd = [boxpd,boxpd,boxpd]
        call mean_vol_pad%new(ldim_pd, params%smpd_crop)
        call mean_vol%pad(mean_vol_pad)
        call mean_vol_pad%fft
        call mean_vol_pad%expand_cmat(params%box_crop)
        call proj_pad%new([boxpd,boxpd,1], params%smpd_crop, wthreads=.false.)
        call mean_img%new([params%box_crop,params%box_crop,1], params%smpd_crop, wthreads=.false.)
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

        do iproj = 1,nproj
            call transform_ptcls(params, build, build%spproj, 'ptcl3D', iproj, registered, group_pinds, &
                &phflip=do_phaseflip)
            if( .not. allocated(group_pinds) ) cycle
            if( size(group_pinds) == 0 ) cycle
            iptcl = group_pinds(1)
            call build%spproj%os_ptcl3D%get_ori(iptcl, o)
            call o%e3set(0.)
            rotmat = o%get_mat()
            proj_dirs(:,iproj) = rotmat(3,:)
            call mean_vol_pad%fproject_serial(o, proj_pad)
            call proj_pad%ifft
            call proj_pad%clip(mean_img)
            call mean_img%fft
            call memoize_ft_maps([params%box_crop,params%box_crop,1],params%smpd_crop)
            do j = 1,size(group_pinds)
                iptcl = group_pinds(j)
                if( iptcl < 1 .or. iptcl > nall ) cycle
                row = p2row(iptcl)
                if( row == 0 ) cycle
                call registered(j)%fft
                call registered(j)%clip(reg_crop)
                call reg_crop%ifft
                call write_registered_stack_image(registered_stack, iptcl, reg_crop)
                write(map_unit,'(I10,1X,I10,1X,I8)') row, iptcl, iproj
                call mean_ctf%copy_fast(mean_img)
                ctfparms = build%spproj%get_ctfparams('ptcl3D',iptcl)
                tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call mean_ctf%apply_ctf(tfun, 'abs', ctfparms)
                call reg_crop%fft
                residual = reg_crop - mean_ctf
                call residual%ifft
                call residual%mask2D_softavg(params%msk_crop)
                vec = residual%serialize()
                if( size(vec) /= npix ) THROW_HARD('registered residual feature size mismatch')
                features(:,row) = vec
                proj_ids(row) = iproj
                deallocate(vec)
            end do
            call forget_ft_maps
            call o%kill
            if( allocated(registered) ) call dealloc_imgarr(registered)
            if( allocated(group_pinds) ) deallocate(group_pinds)
        end do
        close(map_unit)
        if( any(proj_ids == 0) ) THROW_HARD('not every active particle received a registered residual feature')
        call reg_crop%zero_and_unflag_ft
        do iptcl=1,nall
            if( p2row(iptcl) == 0 ) call write_registered_stack_image(registered_stack,iptcl,reg_crop)
        end do
        call write_registered_project(params, build%spproj, registered_stack, registered_project, registered_spproj)
        call registered_spproj%kill
        call mean_vol_pad%kill_expanded
        call mean_vol_pad%kill
        call mean_vol%kill
        call proj_pad%kill
        call mean_img%kill
        call mean_ctf%kill
        call reg_crop%kill
        call residual%kill
        deallocate(p2row)
    end subroutine prepare_flex_diffmap_features

    subroutine write_registered_stack_image( fname, indstk, img )
        type(string), intent(in)    :: fname
        integer,      intent(in)    :: indstk
        type(image),  intent(inout) :: img
        type(imgfile) :: ioimg
        real(kind=c_float), pointer :: rmat_ptr(:,:,:) => null()
        integer :: ldim(3)
        ldim = img%get_ldim()
        call ioimg%open(fname, ldim, img%get_smpd(), formatchar='M', readhead=file_exists(fname))
        call img%get_rmat_ptr(rmat_ptr)
        call ioimg%wmrcSlices(indstk, indstk, rmat_ptr, ldim, img%is_ft())
        call ioimg%close
        rmat_ptr => null()
    end subroutine write_registered_stack_image

    subroutine write_registered_project( params, spproj, stack_fname, project_fname, outproj )
        class(parameters), intent(in)    :: params
        type(sp_project),  intent(inout) :: spproj
        type(string),      intent(in)    :: stack_fname, project_fname
        type(sp_project),  intent(inout) :: outproj
        integer :: iptcl, nptcls
        nptcls = spproj%os_ptcl3D%get_noris()
        call outproj%copy(spproj)
        call outproj%update_projinfo(project_fname)
        if( outproj%os_stk%get_noris() > 0 ) call outproj%os_stk%kill()
        call outproj%os_stk%new(1, is_ptcl=.false.)
        call outproj%os_stk%set(1, 'stk',        simple_abspath(stack_fname))
        call outproj%os_stk%set(1, 'box',        params%box_crop)
        call outproj%os_stk%set(1, 'nptcls',     nptcls)
        call outproj%os_stk%set(1, 'nptcls_stk', nptcls)
        call outproj%os_stk%set(1, 'fromp',      1)
        call outproj%os_stk%set(1, 'top',        nptcls)
        call outproj%os_stk%set(1, 'stkkind',    'single')
        call outproj%os_stk%set(1, 'imgkind',    'ptcl')
        call outproj%os_stk%set(1, 'smpd',       params%smpd_crop)
        call outproj%os_stk%set(1, 'ctf',        'flip')
        call outproj%os_stk%set(1, 'state',      1)
        do iptcl = 1,nptcls
            call outproj%os_ptcl3D%set(iptcl,'stkind',1)
            call outproj%os_ptcl3D%set(iptcl,'indstk',iptcl)
        end do
        call outproj%os_ptcl3D%zero_inpl()
        if( outproj%os_ptcl2D%get_noris() == nptcls )then
            do iptcl = 1,nptcls
                call outproj%os_ptcl2D%set(iptcl,'stkind',1)
                call outproj%os_ptcl2D%set(iptcl,'indstk',iptcl)
            end do
            call outproj%os_ptcl2D%zero_inpl()
        endif
        if( outproj%os_out%get_noris() > 0 ) call outproj%os_out%kill()
        call outproj%write(project_fname)
    end subroutine write_registered_project

end module simple_flex_diffmap_features
