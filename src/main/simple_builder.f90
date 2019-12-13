! centralised builder (the main object constructor in SIMPLE)
module simple_builder
include 'simple_lib.f08'
use simple_image,            only: image
use simple_binimage,         only: binimage
use simple_sp_project,       only: sp_project
use simple_oris,             only: oris
use simple_reconstructor,    only: reconstructor
use simple_reconstructor_eo, only: reconstructor_eo
use simple_sym,              only: sym
use simple_projector,        only: projector
use simple_polarizer,        only: polarizer
use simple_masker,           only: masker
use simple_projection_frcs,  only: projection_frcs
use simple_parameters,       only: parameters
use simple_cmdline,          only: cmdline
implicit none

public :: builder, build_glob
private
#include "simple_local_flags.inc"

type :: builder
    ! GENERAL TOOLBOX
    type(sp_project)                    :: spproj                 !< centralised single-particle project meta-data handler
    class(oris), pointer                :: spproj_field => null() !< pointer to field in spproj
    type(oris)                          :: eulspace, eulspace_red !< discrete spaces(red for reduced)
    type(sym)                           :: pgrpsyms               !< symmetry elements object
    type(image)                         :: img                    !< individual image/projector objects
    type(polarizer)                     :: img_match              !< -"-
    type(image)                         :: img_pad                !< -"-
    type(image)                         :: img_tmp                !< -"-
    type(image)                         :: img_msk                !< -"-
    type(image)                         :: img_copy               !< -"-
    type(projector)                     :: vol, vol_odd
    type(image)                         :: vol2                   !< -"-
    type(masker)                        :: mskimg                 !< mask image
    type(projection_frcs)               :: projfrcs               !< projection FRC's used in the anisotropic Wiener filter
    type(projection_frcs)               :: projpssnrs             !<
    type(image),            allocatable :: imgbatch(:)            !< batch of images
    ! RECONSTRUCTION TOOLBOX
    type(reconstructor_eo)              :: eorecvol               !< object for eo reconstruction
    ! STRATEGY3D TOOLBOX
    type(reconstructor_eo), allocatable :: eorecvols(:)           !< array of volumes for eo-reconstruction
    real,                   allocatable :: fsc(:,:)               !< Fourier Shell Correlation
    integer,                allocatable :: nnmat(:,:)             !< matrix with nearest neighbor indices
    logical,                allocatable :: lmsk(:,:,:)            !< logical circular 2D mask
    logical,                allocatable :: l_resmsk(:)            !< logical resolution mask
    ! PRIVATE EXISTENCE VARIABLES
    logical, private                    :: general_tbox_exists    = .false.
    logical, private                    :: cluster_tbox_exists    = .false.
    logical, private                    :: rec_tbox_exists        = .false.
    logical, private                    :: eo_rec_tbox_exists     = .false.
    logical, private                    :: strategy3D_tbox_exists = .false.
    logical, private                    :: strategy2D_tbox_exists = .false.
    logical, private                    :: extremal3D_tbox_exists = .false.
  contains
    ! HIGH-LEVEL BUILDERS
    procedure                           :: init_params_and_build_spproj
    procedure                           :: init_params_and_build_general_tbox
    procedure                           :: init_params_and_build_strategy2D_tbox
    procedure                           :: init_params_and_build_strategy3D_tbox
    ! LOW-LEVEL BUILDERS
    procedure, private                  :: build_spproj
    procedure                           :: build_general_tbox
    procedure                           :: kill_general_tbox
    procedure                           :: build_rec_eo_tbox
    procedure, private                  :: kill_rec_eo_tbox
    procedure                           :: build_strategy3D_tbox
    procedure, private                  :: kill_strategy3D_tbox
    procedure, private                  :: build_strategy2D_tbox
    procedure                           :: kill_strategy2D_tbox
end type builder

class(builder), pointer :: build_glob  => null()

contains

    ! HIGH-LEVEL BUILDERS

    subroutine init_params_and_build_spproj( self, cline, params )
        class(builder),    target, intent(inout)  :: self
        class(cmdline),    intent(inout)          :: cline
        class(parameters), intent(inout)          :: params
        call params%new(cline)
        call self%build_spproj(params, cline)
        build_glob => self
    end subroutine init_params_and_build_spproj

    subroutine init_params_and_build_general_tbox( self, cline, params, do3d, boxmatch_off )
        class(builder),    target, intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        class(parameters),         intent(inout) :: params
        logical,         optional, intent(in)    :: do3d, boxmatch_off
        logical :: bboxmatch_off
        bboxmatch_off = .false.
        if( present(boxmatch_off) ) bboxmatch_off = boxmatch_off
        call params%new(cline)
        if( bboxmatch_off ) params%boxmatch = params%box
        call self%build_general_tbox(params, cline, do3d=do3d)
        build_glob => self
    end subroutine init_params_and_build_general_tbox

    subroutine init_params_and_build_strategy2D_tbox( self, cline, params )
        class(builder),    target, intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        class(parameters),         intent(inout) :: params
        call params%new(cline)
        call self%build_general_tbox(params, cline, do3d=.false.)
        call self%build_strategy2D_tbox(params)
        build_glob => self
    end subroutine init_params_and_build_strategy2D_tbox

    subroutine init_params_and_build_strategy3D_tbox( self, cline, params )
        class(builder),    target,  intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(parameters),          intent(inout) :: params
        call params%new(cline)
        call self%build_general_tbox(params, cline, do3d=.true.)
        call self%build_strategy3D_tbox(params)
        build_glob => self
    end subroutine init_params_and_build_strategy3D_tbox

    ! LOW-LEVEL BUILDERS

    subroutine build_spproj( self, params, cline )
        use simple_binoris_io,     only: binread_ctfparams_state_eo, binread_oritab
        use simple_user_interface, only: get_prg_ptr
        class(builder), target, intent(inout) :: self
        class(parameters),      intent(inout) :: params
        class(cmdline),         intent(inout) :: cline
        logical :: read_spproj
        ! create object for orientations
        ! b%a is now a pointer to a field in b%spproj
        select case(params%spproj_iseg)
            case(MIC_SEG)
                call self%spproj%os_mic%new(params%nptcls)
                self%spproj_field => self%spproj%os_mic
            case(STK_SEG)
                call self%spproj%os_stk%new(params%nptcls)
                self%spproj_field => self%spproj%os_stk
            case(PTCL2D_SEG)
                call self%spproj%os_ptcl2D%new(params%nptcls)
                self%spproj_field => self%spproj%os_ptcl2D
            case(CLS2D_SEG)
                call self%spproj%os_cls2D%new(params%nptcls)
                self%spproj_field => self%spproj%os_cls2D
            case(CLS3D_SEG)
                call self%spproj%os_cls3D%new(params%nptcls)
                self%spproj_field => self%spproj%os_cls3D
            case(PTCL3D_SEG)
                call self%spproj%os_ptcl3D%new(params%nptcls)
                self%spproj_field => self%spproj%os_ptcl3D
            case DEFAULT
                ! using ptcl3D as the generic segment
                call self%spproj%os_ptcl3D%new(params%nptcls)
                self%spproj_field => self%spproj%os_ptcl3D
        end select
        ! read project file
        read_spproj = .false.
        if( params%sp_required )then
            read_spproj = .true.
        else
            if( cline%defined('projfile') )read_spproj = .true.
        endif
        if( read_spproj )then
            call self%spproj%read(params%projfile)
            ! update cwd of project (in case the params class changed exec dir)
            call self%spproj%projinfo%set(1, 'cwd', trim(params%cwd))
        else
            ! we need the oritab to override the deftab in order not to loose parameters
            if( params%deftab /= '' ) call binread_ctfparams_state_eo(params%deftab,  self%spproj, self%spproj_field, [1,params%nptcls])
            if( params%oritab /= '' ) call binread_oritab(params%oritab,              self%spproj, self%spproj_field, [1,params%nptcls])
        endif
        if( .not. associated(build_glob) ) build_glob => self
        write(logfhandle,'(A)') '>>> DONE BUILDING SP PROJECT'
    end subroutine build_spproj

    subroutine build_general_tbox( self, params, cline, do3d )
        class(builder), target, intent(inout) :: self
        class(parameters),      intent(inout) :: params
        class(cmdline),         intent(inout) :: cline
        logical, optional,      intent(in)    :: do3d
        type(image) :: mskimg
        integer     :: lfny, cyc_lims(3,2)
        logical     :: ddo3d
        call self%kill_general_tbox
        ddo3d = .true.
        if( present(do3d) ) ddo3d = do3d
        ! set up symmetry functionality
        call self%pgrpsyms%new(trim(params%pgrp))
        params%nsym    = self%pgrpsyms%get_nsym()
        params%eullims = self%pgrpsyms%get_eullims()
        ! build spproj
        call self%build_spproj(params, cline)
        ! states exception
        if( self%spproj_field%get_n('state') > 1 )then
            if( .not. cline%defined('nstates') )then
                THROW_WARN('your input doc has multiple states but NSTATES is not given')
            endif
        endif
        ! generate discrete projection direction spaces
        if( ddo3d )then
            call self%eulspace%new(params%nspace)
            call self%pgrpsyms%build_refspiral(self%eulspace)
            call self%eulspace_red%new(NSPACE_REDUCED)
            call self%pgrpsyms%build_refspiral(self%eulspace_red)
        endif
        if( params%box > 0 )then
            ! build image objects
            ! box-sized ones
            call self%img%new([params%box,params%box,1],params%smpd,                 wthreads=.false.)
            call self%img_match%new([params%boxmatch,params%boxmatch,1],params%smpd, wthreads=.false.)
            call self%img_copy%new([params%box,params%box,1],params%smpd,  wthreads=.false.)
            ! for thread safety in the image class
            call self%img%construct_thread_safe_tmp_imgs(params%nthr)
            ! boxmatch-sized ones
            call self%img_tmp%new([params%boxmatch,params%boxmatch,1],params%smpd,   wthreads=.false.)
            call self%img_msk%new([params%boxmatch,params%boxmatch,1],params%smpd,   wthreads=.false.)
            call self%mskimg%new ([params%boxmatch,params%boxmatch,1],params%smpd,   wthreads=.false.)
            ! boxpd-sized ones
            call self%img_pad%new([params%boxpd,params%boxpd,1],params%smpd)
            if( ddo3d )then
                call self%vol%new([params%box,params%box,params%box], params%smpd)
                call self%vol2%new([params%box,params%box,params%box], params%smpd)
            endif
            ! build arrays
            lfny       = self%img%get_lfny(1)
            cyc_lims   = self%img_pad%loop_lims(3)
            allocate( self%fsc(params%nstates,lfny), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk("In: build_general_tbox; simple_builder, 1", alloc_stat)
            self%fsc  = 0.
            ! set default amsklp
            if( .not. cline%defined('amsklp') .and. cline%defined('lp') )then
                params%amsklp = self%img%get_lp(self%img%get_find(params%lp)-2)
            endif
            ! generate logical circular 2D mask
            call mskimg%disc([params%box,params%box,1], params%smpd, params%msk, self%lmsk)
            call mskimg%kill
        endif
        if( params%projstats .eq. 'yes' )then
            if( .not. self%spproj_field%isthere('proj') ) call self%spproj_field%set_projs(self%eulspace)
        endif
        ! resolution mask for correlation calculation (omitting shells corresponding to the graphene signal if params%l_graphene = .true.)
        if( params%boxmatch > 0 ) then
            self%l_resmsk = calc_graphene_mask(params%boxmatch, params%smpd)
            if( .not. params%l_graphene ) self%l_resmsk = .true.
        end if
        ! associate global build pointer
        if( .not. associated(build_glob) ) build_glob => self
        self%general_tbox_exists = .true.
        write(logfhandle,'(A)') '>>> DONE BUILDING GENERAL TOOLBOX'
    end subroutine build_general_tbox

    subroutine kill_general_tbox( self )
        class(builder), intent(inout)  :: self
        if( self%general_tbox_exists )then
            call self%pgrpsyms%kill
            call self%spproj_field%kill
            self%spproj_field => null()
            call self%spproj%kill
            call self%eulspace%kill
            call self%eulspace_red%kill
            call self%img%kill
            call self%img_match%kill_polarizer
            call self%img_match%kill
            call self%img_copy%kill
            call self%img_tmp%kill
            call self%img_msk%kill
            call self%img_pad%kill
            call self%vol%kill_expanded
            call self%vol%kill
            call self%vol2%kill
            call self%mskimg%kill
            if( allocated(self%fsc)  )     deallocate(self%fsc)
            if( allocated(self%lmsk) )     deallocate(self%lmsk)
            if( allocated(self%l_resmsk) ) deallocate(self%l_resmsk)
            self%general_tbox_exists = .false.
        endif
    end subroutine kill_general_tbox

    subroutine build_rec_eo_tbox( self, params )
        class(builder), target, intent(inout) :: self
        class(parameters),      intent(inout) :: params
        call self%kill_rec_eo_tbox
        call self%eorecvol%new(self%spproj)
        if( .not. self%spproj_field%isthere('proj') ) call self%spproj_field%set_projs(self%eulspace)
        call self%projfrcs%new(NSPACE_REDUCED, params%box, params%smpd, params%nstates)
        if( .not. associated(build_glob) ) build_glob => self
        self%eo_rec_tbox_exists = .true.
        write(logfhandle,'(A)') '>>> DONE BUILDING EO RECONSTRUCTION TOOLBOX'
    end subroutine build_rec_eo_tbox

    subroutine kill_rec_eo_tbox( self )
        class(builder), intent(inout) :: self
        if( self%eo_rec_tbox_exists )then
            call self%eorecvol%kill
            call self%projfrcs%kill
            self%eo_rec_tbox_exists = .false.
        endif
    end subroutine kill_rec_eo_tbox

    subroutine build_strategy2D_tbox( self, params )
        class(builder), target, intent(inout) :: self
        class(parameters),      intent(inout) :: params
        call self%kill_strategy2D_tbox
        if( params%neigh.eq.'yes' )then
            if( self%spproj%os_cls3D%get_noris() == params%ncls )then
                allocate(self%nnmat(params%ncls,params%nnn))
                call self%spproj%os_cls3D%nearest_proj_neighbors(params%nnn, self%nnmat)
            else
               THROW_HARD('size of os_cls3D segment of spproj does not conform with # clusters (ncls)')
            endif
        endif
        call self%projfrcs%new(params%ncls, params%box, params%smpd, params%nstates)
        call self%projpssnrs%new(params%ncls, params%box, params%smpd, params%nstates)
        if( .not. associated(build_glob) ) build_glob => self
        self%strategy2D_tbox_exists = .true.
        write(logfhandle,'(A)') '>>> DONE BUILDING STRATEGY2D TOOLBOX'
    end subroutine build_strategy2D_tbox

    subroutine kill_strategy2D_tbox( self )
        class(builder), intent(inout) :: self
        if( self%strategy2D_tbox_exists )then
            if( allocated(self%nnmat) ) deallocate(self%nnmat)
            call self%projfrcs%kill
            call self%projpssnrs%kill
            self%strategy2D_tbox_exists = .false.
        endif
    end subroutine kill_strategy2D_tbox

    subroutine build_strategy3D_tbox( self, params )
        class(builder), target, intent(inout) :: self
        class(parameters),      intent(inout) :: params
        call self%kill_strategy3D_tbox
        allocate( self%eorecvols(params%nstates), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('build_strategy3D_tbox; simple_builder, 1', alloc_stat)
        if( str_has_substr(params%refine, 'neigh') )then
            allocate(self%nnmat(params%npeaks,params%nnn))
        else if( params%neigh .eq. 'yes' )then
            allocate(self%nnmat(params%nspace,params%nnn))
            call self%pgrpsyms%nearest_proj_neighbors(self%eulspace, params%nnn, self%nnmat)
        endif
        if( .not. self%spproj_field%isthere('proj') ) call self%spproj_field%set_projs(self%eulspace)
        if( params%clsfrcs.eq.'yes' )then
            ! done later
        else
            call self%projfrcs%new(NSPACE_REDUCED, params%box, params%smpd, params%nstates)
        endif
        if( .not. associated(build_glob) ) build_glob => self
        self%strategy3D_tbox_exists = .true.
        write(logfhandle,'(A)') '>>> DONE BUILDING STRATEGY3D TOOLBOX'
    end subroutine build_strategy3D_tbox

    subroutine kill_strategy3D_tbox( self )
        class(builder), intent(inout) :: self
        integer :: i
        if( self%strategy3D_tbox_exists )then
            if( allocated(self%eorecvols) )then
                do i=1,size(self%eorecvols)
                    call self%eorecvols(i)%kill
                end do
                deallocate(self%eorecvols)
            endif
            if( allocated(self%nnmat) ) deallocate(self%nnmat)
            call self%projfrcs%kill
            self%strategy3D_tbox_exists = .false.
        endif
    end subroutine kill_strategy3D_tbox

end module simple_builder
