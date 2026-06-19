module simple_diffmap_denoise_project_strategy
use simple_core_module_api
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_cmdline,            only: cmdline
use simple_sp_project,         only: sp_project
use simple_image,              only: image
use simple_image_msk,          only: automask2D
use simple_classaverager,      only: transform_ptcls
use simple_imgfile,            only: imgfile
implicit none

public :: diffmap_denoise_project_strategy
public :: diffmap_denoise_project_shmem_strategy
public :: create_diffmap_denoise_project_strategy

private
#include "simple_local_flags.inc"

type, abstract :: diffmap_denoise_project_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
end type diffmap_denoise_project_strategy

type, extends(diffmap_denoise_project_strategy) :: diffmap_denoise_project_shmem_strategy
contains
    procedure :: initialize   => shmem_initialize
    procedure :: execute      => shmem_execute
    procedure :: finalize_run => shmem_finalize_run
    procedure :: cleanup      => shmem_cleanup
end type diffmap_denoise_project_shmem_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: diffmap_denoise_project_strategy, parameters, builder, cmdline
        class(diffmap_denoise_project_strategy), intent(inout) :: self
        type(parameters),                        intent(inout) :: params
        type(builder),                           intent(inout) :: build
        class(cmdline),                          intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, build, cline)
        import :: diffmap_denoise_project_strategy, parameters, builder, cmdline
        class(diffmap_denoise_project_strategy), intent(inout) :: self
        type(parameters),                        intent(inout) :: params
        type(builder),                           intent(inout) :: build
        class(cmdline),                          intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, build, cline)
        import :: diffmap_denoise_project_strategy, parameters, builder, cmdline
        class(diffmap_denoise_project_strategy), intent(inout) :: self
        type(parameters),                        intent(in)    :: params
        type(builder),                           intent(inout) :: build
        class(cmdline),                          intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params)
        import :: diffmap_denoise_project_strategy, parameters
        class(diffmap_denoise_project_strategy), intent(inout) :: self
        type(parameters),                        intent(in)    :: params
    end subroutine cleanup_interface
end interface

contains

    function create_diffmap_denoise_project_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(diffmap_denoise_project_strategy), allocatable :: strategy
        integer :: nparts
        nparts = 1
        if( cline%defined('nparts') ) nparts = max(1, cline%get_iarg('nparts'))
        if( nparts > 1 .or. cline%defined('part') )then
            THROW_HARD('diffmap_denoise_project currently supports shared-memory execution only')
        endif
        allocate(diffmap_denoise_project_shmem_strategy :: strategy)
        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DIFFMAP_DENOISE_PROJECT (SHARED-MEMORY)'
    end function create_diffmap_denoise_project_strategy

    subroutine apply_defaults(cline)
        use simple_default_clines, only: set_automask2D_defaults
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',    'yes')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype',  'ptcl2D')
        if( .not. cline%defined('pca_mode') ) call cline%set('pca_mode', 'diffusion_maps')
        if( .not. cline%defined('graph')    ) call cline%set('graph',    'euc')
        if( .not. cline%defined('steering') ) call cline%set('steering', 'none')
        call set_automask2D_defaults(cline)
    end subroutine apply_defaults

    subroutine shmem_initialize(self, params, build, cline)
        class(diffmap_denoise_project_shmem_strategy), intent(inout) :: self
        type(parameters),                              intent(inout) :: params
        type(builder),                                 intent(inout) :: build
        class(cmdline),                                intent(inout) :: cline
        call apply_defaults(cline)
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
    end subroutine shmem_initialize

    subroutine shmem_execute(self, params, build, cline)
        class(diffmap_denoise_project_shmem_strategy), intent(inout) :: self
        type(parameters),                              intent(inout) :: params
        type(builder),                                 intent(inout) :: build
        class(cmdline),                                intent(inout) :: cline
        type(sp_project) :: spproj
        call spproj%read(params%projfile)
        call run_diffmap_denoise_project(params, build, spproj)
        call spproj%kill
    end subroutine shmem_execute

    subroutine shmem_finalize_run(self, params, build, cline)
        class(diffmap_denoise_project_shmem_strategy), intent(inout) :: self
        type(parameters),                              intent(in)    :: params
        type(builder),                                 intent(inout) :: build
        class(cmdline),                                intent(inout) :: cline
    end subroutine shmem_finalize_run

    subroutine shmem_cleanup(self, params)
        class(diffmap_denoise_project_shmem_strategy), intent(inout) :: self
        type(parameters),                              intent(in)    :: params
    end subroutine shmem_cleanup

    subroutine run_diffmap_denoise_project(params, build, spproj)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(sp_project), intent(inout) :: spproj
        type(sp_project) :: outproj
        type(string), allocatable :: raw_stks(:), den_stks(:)
        integer, allocatable :: cls_inds(:), cls_pops(:)
        logical, allocatable :: processed(:)
        logical :: l_phflip
        integer :: icls, nptcls, nprocessed
        call validate_diffmap_denoise_project(params, spproj, cls_inds, cls_pops, l_phflip)
        nptcls = spproj%os_ptcl2D%get_noris()
        allocate(processed(nptcls), source=.false.)
        call prepare_diffmap_particle_stacks(spproj, raw_stks, den_stks)
        nprocessed = 0
        do icls = 1, size(cls_inds)
            write(logfhandle,'(A,I8,A,I8,A,I8)') 'Diffmap denoise project: class=', cls_inds(icls), &
                ' class_index=', icls, ' nptcls=', cls_pops(icls)
            call flush(logfhandle)
            call write_diffmap_denoised_class(params, build, spproj, cls_inds(icls), l_phflip, &
                raw_stks, den_stks, processed, nprocessed)
        end do
        if( any(.not. processed) )then
            write(logfhandle,'(A,I10)') 'unprocessed particles: ', count(.not. processed)
            THROW_HARD('diffmap_denoise_project requires every ptcl2D particle to have a valid class assignment')
        endif
        call write_diffmap_project(params, spproj, raw_stks, den_stks, outproj)
        write(logfhandle,'(A,A)') 'Diffmap denoise project written: ', params%projfile%to_char()
        write(logfhandle,'(A,I10)') 'Diffmap denoise project particles: ', nprocessed
        call flush(logfhandle)
        call outproj%kill
        if( allocated(raw_stks) ) deallocate(raw_stks)
        if( allocated(den_stks) ) deallocate(den_stks)
        if( allocated(cls_inds) ) deallocate(cls_inds)
        if( allocated(cls_pops) ) deallocate(cls_pops)
        if( allocated(processed) ) deallocate(processed)
    end subroutine run_diffmap_denoise_project

    subroutine validate_diffmap_denoise_project(params, spproj, cls_inds, cls_pops, l_phflip)
        type(parameters), intent(in)    :: params
        type(sp_project), intent(inout) :: spproj
        integer, allocatable, intent(out) :: cls_inds(:), cls_pops(:)
        logical, intent(out) :: l_phflip
        type(string) :: ctfstr
        integer, allocatable :: pinds(:)
        logical, allocatable :: seen_slots(:)
        integer :: iptcl, istk, icls, nptcls, nstks, cls, stkind, indstk, slot, ctf_mode
        if( trim(params%oritype) /= 'ptcl2D' )then
            THROW_HARD('diffmap_denoise_project supports oritype=ptcl2D input only')
        endif
        if( trim(params%pca_mode) /= 'diffusion_maps' )then
            THROW_HARD('diffmap_denoise_project supports pca_mode=diffusion_maps only')
        endif
        select case(trim(params%graph))
            case('euc','euclidean')
            case DEFAULT
                THROW_HARD('diffmap_denoise_project supports non-steerable Euclidean diffusion maps only; use graph=euc')
        end select
        if( trim(params%steering) /= 'none' )then
            THROW_HARD('diffmap_denoise_project supports non-steerable diffusion maps only; use steering=none')
        endif
        nptcls = spproj%os_ptcl2D%get_noris()
        if( nptcls < 1 ) THROW_HARD('empty ptcl2D field; diffmap_denoise_project')
        if( spproj%os_ptcl3D%get_noris() /= nptcls )then
            THROW_HARD('ptcl2D/ptcl3D particle counts differ; diffmap_denoise_project')
        endif
        nstks = spproj%os_stk%get_noris()
        if( nstks < 1 ) THROW_HARD('empty stack field; diffmap_denoise_project')
        ctf_mode = 0
        do istk = 1, nstks
            if( .not. spproj%os_stk%isthere(istk, 'ctf') ) THROW_HARD('stack ctf key missing; diffmap_denoise_project')
            call spproj%os_stk%getter(istk, 'ctf', ctfstr)
            select case(trim(ctfstr%to_char()))
                case('yes')
                    if( ctf_mode == 0 ) ctf_mode = CTFFLAG_YES
                    if( ctf_mode /= CTFFLAG_YES ) THROW_HARD('mixed ctf=yes/ctf=flip stacks are not supported; diffmap_denoise_project')
                case('flip')
                    if( ctf_mode == 0 ) ctf_mode = CTFFLAG_FLIP
                    if( ctf_mode /= CTFFLAG_FLIP ) THROW_HARD('mixed ctf=yes/ctf=flip stacks are not supported; diffmap_denoise_project')
                case DEFAULT
                    THROW_HARD('diffmap_denoise_project requires ctf=yes or ctf=flip input stacks')
            end select
        enddo
        l_phflip = ctf_mode == CTFFLAG_YES
        allocate(seen_slots(nptcls), source=.false.)
        do iptcl = 1, nptcls
            if( spproj%os_ptcl2D%get_state(iptcl) <= 0 )then
                THROW_HARD('diffmap_denoise_project requires all ptcl2D particles to be active')
            endif
            if( .not. spproj%os_ptcl2D%isthere(iptcl, 'class') )then
                THROW_HARD('ptcl2D class assignment missing; diffmap_denoise_project')
            endif
            cls = spproj%os_ptcl2D%get_int(iptcl, 'class')
            if( cls <= 0 ) THROW_HARD('ptcl2D class assignment must be positive; diffmap_denoise_project')
            call spproj%map_ptcl_ind2stk_ind('ptcl2D', iptcl, stkind, indstk)
            if( stkind < 1 .or. stkind > nstks ) THROW_HARD('invalid particle stkind; diffmap_denoise_project')
            if( indstk < 1 .or. indstk > spproj%os_stk%get_int(stkind, 'nptcls_stk') )then
                THROW_HARD('invalid particle indstk; diffmap_denoise_project')
            endif
            slot = spproj%os_stk%get_fromp(stkind) + indstk - 1
            if( slot < 1 .or. slot > nptcls ) THROW_HARD('invalid particle stack slot; diffmap_denoise_project')
            if( seen_slots(slot) ) THROW_HARD('duplicate particle stack slot; diffmap_denoise_project')
            seen_slots(slot) = .true.
        enddo
        if( any(.not. seen_slots) ) THROW_HARD('particle stack slots do not cover all particles; diffmap_denoise_project')
        cls_inds = spproj%os_ptcl2D%get_label_inds('class')
        cls_inds = pack(cls_inds, mask=cls_inds > 0)
        if( size(cls_inds) < 1 ) THROW_HARD('No positive ptcl2D classes found; diffmap_denoise_project')
        allocate(cls_pops(size(cls_inds)), source=0)
        do icls = 1, size(cls_inds)
            call spproj%os_ptcl2D%get_pinds(cls_inds(icls), 'class', pinds)
            if( allocated(pinds) )then
                cls_pops(icls) = size(pinds)
                deallocate(pinds)
            endif
            if( cls_pops(icls) < 3 )then
                write(logfhandle,'(A,I8,A,I8)') 'class=', cls_inds(icls), ' pop=', cls_pops(icls)
                THROW_HARD('diffmap_denoise_project requires at least 3 particles per class')
            endif
        end do
        if( sum(cls_pops) /= nptcls ) THROW_HARD('class populations do not cover all particles; diffmap_denoise_project')
        if( allocated(seen_slots) ) deallocate(seen_slots)
        call ctfstr%kill
    end subroutine validate_diffmap_denoise_project

    subroutine prepare_diffmap_particle_stacks(spproj, raw_stks, den_stks)
        type(sp_project), intent(inout) :: spproj
        type(string), allocatable, intent(out) :: raw_stks(:), den_stks(:)
        type(string) :: raw_dir, den_dir, raw_tab, den_tab
        integer :: istk, nstks, numlen
        nstks = spproj%os_stk%get_noris()
        allocate(raw_stks(nstks), den_stks(nstks))
        raw_dir = 'diffmap_raw_stacks'
        den_dir = 'diffmap_denoised_stacks'
        call simple_mkdir(raw_dir)
        call simple_mkdir(den_dir)
        numlen = max(1, len_trim(int2str(nstks)))
        do istk = 1, nstks
            raw_stks(istk) = filepath(raw_dir, string('diffmap_raw_stk')//int2str_pad(istk, numlen)//STK_EXT)
            den_stks(istk) = filepath(den_dir, string('diffmap_den_stk')//int2str_pad(istk, numlen)//STK_EXT)
            call del_file(raw_stks(istk))
            call del_file(den_stks(istk))
        end do
        raw_tab = 'diffmap_raw_stktab.txt'
        den_tab = 'diffmap_denoised_stktab.txt'
        call write_filetable(raw_tab, raw_stks)
        call write_filetable(den_tab, den_stks)
        call raw_dir%kill
        call den_dir%kill
        call raw_tab%kill
        call den_tab%kill
    end subroutine prepare_diffmap_particle_stacks

    subroutine write_diffmap_denoised_class(params, build, spproj, cls_id, l_phflip, raw_stks, den_stks, &
                                            processed, nprocessed)
        use simple_diff_map_graphs,  only: diffmap_graph, build_cls_split_graph
        use simple_diff_map_denoise, only: graph_coeffproj_denoise
        use simple_imgarr_utils,     only: dealloc_imgarr, copy_imgarr
        use simple_imgproc,          only: make_pcavecs
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: cls_id
        logical,          intent(in)    :: l_phflip
        type(string),     intent(in)    :: raw_stks(:), den_stks(:)
        logical,          intent(inout) :: processed(:)
        integer,          intent(inout) :: nprocessed
        type(parameters) :: params_mask
        type(image), allocatable :: imgs(:), imgs_ppca(:), class_mask(:), den_ptcls(:)
        type(image) :: cavg_raw
        type(diffmap_graph) :: graph
        real, allocatable :: avg(:), avg_ptcls(:), pcavecs(:,:)
        real, allocatable :: class_diams(:), class_shifts(:,:)
        integer, allocatable :: pinds(:)
        integer :: nptcls, npix, j, class_ldim(3), stkind, indstk
        real :: class_moldiam, class_mskdiam, class_mskrad, sdev_noise
        call transform_ptcls(params, build, spproj, 'ptcl2D', cls_id, imgs, pinds, phflip=l_phflip, cavg=cavg_raw)
        if( .not. allocated(imgs) ) THROW_HARD('transform_ptcls returned no images; diffmap_denoise_project')
        nptcls = size(imgs)
        if( nptcls < 3 ) THROW_HARD('at least 3 particles per class required; diffmap_denoise_project')
        imgs_ppca = copy_imgarr(imgs)
        call imgs_ppca(1)%memoize_mask_coords()
        allocate(class_mask(1))
        call class_mask(1)%copy(cavg_raw)
        params_mask%ngrow   = params%ngrow
        params_mask%winsz   = params%winsz
        params_mask%edge    = params%edge
        params_mask%amsklp  = params%amsklp
        params_mask%automsk = params%automsk
        params_mask%part    = params%part
        class_ldim          = cavg_raw%get_ldim()
        params_mask%box     = class_ldim(1)
        params_mask%smpd    = cavg_raw%get_smpd()
        params_mask%msk     = real(class_ldim(1) / 2) - COSMSKHALFWIDTH
        call automask2D(params_mask, class_mask, params_mask%ngrow, nint(params_mask%winsz), params_mask%edge, &
                        class_diams, class_shifts)
        class_moldiam = params_mask%smpd * real(min(round2even(class_diams(1) / params_mask%smpd + &
            2. * COSMSKHALFWIDTH), class_ldim(1)))
        class_mskdiam = class_moldiam * MSK_EXP_FAC
        class_mskrad  = min(real(class_ldim(1) / 2) - COSMSKHALFWIDTH - 1., 0.5 * class_mskdiam / params_mask%smpd)
        do j = 1, nptcls
            call imgs_ppca(j)%norm_noise(build%lmsk, sdev_noise)
            call imgs_ppca(j)%mask2D_softavg(class_mskrad)
        end do
        if( trim(params%pre_norm) == 'yes' )then
            do j = 1, nptcls
                call imgs_ppca(j)%norm
            end do
        endif
        call make_pcavecs(imgs_ppca, npix, avg, pcavecs, transp=.false.)
        call build_cls_split_graph(params, spproj, pinds, pcavecs, imgs_ppca, graph)
        if( graph%n /= nptcls ) THROW_HARD('diffusion-map graph size mismatch; diffmap_denoise_project')
        if( trim(graph%metric) /= 'euclidean' .or. trim(graph%steering) /= 'none' )then
            THROW_HARD('diffmap_denoise_project expected a non-steerable Euclidean graph')
        endif
        avg_ptcls = cavg_raw%serialize()
        call graph_coeffproj_denoise(params, imgs, avg_ptcls, graph, den_ptcls)
        if( .not. allocated(den_ptcls) ) THROW_HARD('diffusion-map generative denoising failed; diffmap_denoise_project')
        do j = 1, nptcls
            if( processed(pinds(j)) ) THROW_HARD('particle was processed twice; diffmap_denoise_project')
            call spproj%map_ptcl_ind2stk_ind('ptcl2D', pinds(j), stkind, indstk)
            if( stkind < 1 .or. stkind > size(raw_stks) ) THROW_HARD('invalid output stack index; diffmap_denoise_project')
            call write_diffmap_stack_image(raw_stks(stkind), indstk, imgs(j))
            call write_diffmap_stack_image(den_stks(stkind), indstk, den_ptcls(j))
            processed(pinds(j)) = .true.
            nprocessed = nprocessed + 1
        end do
        call cavg_raw%kill
        call graph%kill()
        if( allocated(imgs) ) call dealloc_imgarr(imgs)
        if( allocated(imgs_ppca) ) call dealloc_imgarr(imgs_ppca)
        if( allocated(class_mask) ) call dealloc_imgarr(class_mask)
        if( allocated(den_ptcls) ) call dealloc_imgarr(den_ptcls)
        if( allocated(pinds) ) deallocate(pinds)
        if( allocated(avg) ) deallocate(avg)
        if( allocated(avg_ptcls) ) deallocate(avg_ptcls)
        if( allocated(pcavecs) ) deallocate(pcavecs)
        if( allocated(class_diams) ) deallocate(class_diams)
        if( allocated(class_shifts) ) deallocate(class_shifts)
    end subroutine write_diffmap_denoised_class

    subroutine write_diffmap_stack_image(fname, indstk, img)
        type(string), intent(in)    :: fname
        integer,      intent(in)    :: indstk
        type(image),  intent(inout) :: img
        type(imgfile) :: ioimg
        real(kind=c_float), pointer :: rmat_ptr(:,:,:) => null()
        integer :: ldim(3)
        if( indstk < 1 ) THROW_HARD('invalid stack-local output index; diffmap_denoise_project')
        if( .not. img%is_2d() ) THROW_HARD('2D particle image expected; diffmap_denoise_project')
        ldim = img%get_ldim()
        ldim(3) = 1
        call ioimg%open(fname, ldim, img%get_smpd(), formatchar='M', readhead=file_exists(fname))
        call img%get_rmat_ptr(rmat_ptr)
        call ioimg%wmrcSlices(indstk, indstk, rmat_ptr, ldim, img%is_ft())
        call ioimg%close
        rmat_ptr => null()
    end subroutine write_diffmap_stack_image

    subroutine write_diffmap_project(params, spproj, raw_stks, den_stks, outproj)
        type(parameters), intent(in)    :: params
        type(sp_project), intent(inout) :: spproj
        type(string),     intent(in)    :: raw_stks(:), den_stks(:)
        type(sp_project), intent(inout) :: outproj
        integer :: istk, nstks, nraw, nden, ldim_raw(3), ldim_den(3), expected_nptcls
        nstks = spproj%os_stk%get_noris()
        if( size(raw_stks) /= nstks .or. size(den_stks) /= nstks )then
            THROW_HARD('stack-table size mismatch; diffmap_denoise_project')
        endif
        call outproj%copy(spproj)
        call outproj%update_projinfo(params%projfile)
        do istk = 1, nstks
            call find_ldim_nptcls(raw_stks(istk), ldim_raw, nraw)
            call find_ldim_nptcls(den_stks(istk), ldim_den, nden)
            expected_nptcls = spproj%os_stk%get_int(istk, 'nptcls_stk')
            if( nraw /= expected_nptcls .or. nden /= expected_nptcls )then
                THROW_HARD('output stack count does not match input stack row; diffmap_denoise_project')
            endif
            if( .not. all(ldim_raw(1:2) == ldim_den(1:2)) )then
                THROW_HARD('raw/denoised output stack dimensions differ; diffmap_denoise_project')
            endif
            call outproj%os_stk%set(istk, 'stk',        simple_abspath(raw_stks(istk)))
            call outproj%os_stk%set(istk, 'stk_den',    simple_abspath(den_stks(istk)))
            call outproj%os_stk%set(istk, 'box',        ldim_raw(1))
            call outproj%os_stk%set(istk, 'nptcls',     nraw)
            call outproj%os_stk%set(istk, 'nptcls_stk', nraw)
            call outproj%os_stk%set(istk, 'ctf',        'flip')
            call outproj%os_stk%set(istk, 'imgkind',    'ptcl')
        end do
        call outproj%os_ptcl2D%zero_inpl()
        call outproj%os_ptcl3D%zero_inpl()
        if( outproj%os_out%get_noris() > 0 ) call outproj%os_out%kill()
        call outproj%write(params%projfile)
    end subroutine write_diffmap_project

end module simple_diffmap_denoise_project_strategy
