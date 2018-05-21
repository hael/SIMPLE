! centralised builder (the main object constructor in SIMPLE)
module simple_builder
include 'simple_lib.f08'
use simple_comlin,           only: comlin
use simple_image,            only: image
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
    type(image),            allocatable :: imgbatch(:)            !< batch of images
    ! COMMON LINES TOOLBOX
    type(image),            allocatable :: imgs(:)                !< images (all should be read in)
    type(image),            allocatable :: imgs_sym(:)            !< images (all should be read in)
    type(comlin)                        :: clins                  !< common lines data structure
    type(image),            allocatable :: ref_imgs(:,:)          !< array of reference images
    ! RECONSTRUCTION TOOLBOX
    type(reconstructor_eo)              :: eorecvol               !< object for eo reconstruction
    type(reconstructor)                 :: recvol                 !< object for reconstruction
    ! STRATEGY3D TOOLBOX
    type(reconstructor),    allocatable :: recvols(:)             !< array of volumes for reconstruction
    type(reconstructor_eo), allocatable :: eorecvols(:)           !< array of volumes for eo-reconstruction
    real,                   allocatable :: fsc(:,:)               !< Fourier Shell Correlation
    integer,                allocatable :: nnmat(:,:)             !< matrix with nearest neighbor indices
    integer,                allocatable :: grid_projs(:)          !< projection directions for coarse grid search
    ! PRIVATE EXISTENCE VARIABLES
    logical, private                    :: general_tbox_exists    = .false.
    logical, private                    :: cluster_tbox_exists    = .false.
    logical, private                    :: comlin_tbox_exists     = .false.
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
    procedure                           :: build_comlin_tbox
    procedure, private                  :: kill_comlin_tbox
    procedure                           :: build_rec_tbox
    procedure, private                  :: kill_rec_tbox
    procedure                           :: build_rec_eo_tbox
    procedure, private                  :: kill_rec_eo_tbox
    procedure, private                  :: build_strategy3D_tbox
    procedure, private                  :: kill_strategy3D_tbox
    procedure, private                  :: build_strategy2D_tbox
    procedure, private                  :: kill_strategy2D_tbox
    procedure                           :: build_extremal3D_tbox
    procedure, private                  :: kill_extremal3D_tbox
end type builder

class(builder), pointer :: build_glob  => null()

contains

    ! HIGH-LEVEL BUILDERS

    subroutine init_params_and_build_spproj( self, cline, params )
        class(builder),    target, intent(inout)  :: self
        class(cmdline),    intent(inout)          :: cline
        class(parameters), target,  intent(inout) :: params
        call params%new(cline)
        call self%build_spproj(params, cline)
        build_glob => self
    end subroutine init_params_and_build_spproj

    subroutine init_params_and_build_general_tbox( self, cline, params, do3d, boxmatch_off )
        class(builder),    target, intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        class(parameters), target, intent(inout) :: params
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
        class(parameters), target, intent(inout) :: params
        call params%new(cline)
        call self%build_general_tbox(params, cline, do3d=.false.)
        call self%build_strategy2D_tbox(params)
        build_glob => self
    end subroutine init_params_and_build_strategy2D_tbox

    subroutine init_params_and_build_strategy3D_tbox( self, cline, params )
        class(builder),    target,  intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(parameters), target,  intent(inout) :: params
        call params%new(cline)
        call self%build_general_tbox(params, cline, do3d=.true.)
        call self%build_strategy2D_tbox(params)
        build_glob => self
    end subroutine init_params_and_build_strategy3D_tbox

    ! LOW-LEVEL BUILDERS

    subroutine build_spproj( self, params, cline )
        use simple_binoris_io,     only: binread_ctfparams_state_eo, binread_oritab
        use simple_user_interface, only: get_prg_ptr
        class(builder), target, intent(inout) :: self
        class(parameters),      intent(inout) :: params
        class(cmdline),         intent(inout) :: cline
        logical ::  metadata_read
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
        ! read from project file
        metadata_read = .false.
        if( cline%defined('projfile') )then
            if( .not. file_exists(params%projfile) )then
                write(*,*) 'project file not in exec dir: ', trim(params%projfile)
                stop 'ERROR! build :: build_spproj'
            endif
            metadata_read = .true.
        else if( associated(params%ptr2prg) )then
            if( params%ptr2prg%requires_sp_project() ) metadata_read = .true.
        endif
        if( metadata_read )then
            call self%spproj%read(params%projfile)
            ! update cwd of project (in case the params class changed exec dir)
            call self%spproj%projinfo%set(1, 'cwd', trim(params%cwd))
        else
            ! we need the oritab to override the deftab in order not to loose parameters
            if( params%deftab /= '' ) call binread_ctfparams_state_eo(params%deftab,  self%spproj, self%spproj_field, [1,params%nptcls])
            if( params%oritab /= '' ) call binread_oritab(params%oritab,              self%spproj, self%spproj_field, [1,params%nptcls])
            DebugPrint 'read deftab'
        endif
        write(*,'(A)') '>>> DONE BUILDING SP PROJECT'
    end subroutine build_spproj

    subroutine build_general_tbox( self, params, cline, do3d )
        class(builder),    intent(inout) :: self
        class(parameters), intent(inout) :: params
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: do3d
        integer :: lfny,  lfny_match, cyc_lims(3,2)
        logical :: ddo3d, fforce_ctf
        call self%kill_general_tbox
        ddo3d = .true.
        if( present(do3d) ) ddo3d = do3d
        ! set up symmetry functionality
        call self%pgrpsyms%new(trim(params%pgrp))
        params%nsym    = self%pgrpsyms%get_nsym()
        params%eullims = self%pgrpsyms%srchrange()
        DebugPrint   'did setup symmetry functionality'
        ! build spproj
        call self%build_spproj(params, cline)
        ! states exception
        if( self%spproj_field%get_n('state') > 1 )then
            if( .not. cline%defined('nstates') )then
                write(*,'(a)') 'WARNING, your input doc has multiple states but NSTATES is not given'
            endif
        endif
        DebugPrint 'created & filled object for orientations'
        ! generate discrete projection direction spaces
        call self%eulspace%new(params%nspace)
        call self%eulspace%spiral(params%nsym, params%eullims)
        call self%eulspace_red%new(NSPACE_REDUCED)
        call self%eulspace_red%spiral(params%nsym, params%eullims)
        ! create angular subspace
        self%grid_projs = self%eulspace%create_proj_subspace(NPDIRS_SUBSPACE, params%nsym, params%eullims)
        DebugPrint 'generated discrete projection direction space'
        if( params%box > 0 )then
            ! build image objects
            ! box-sized ones
            call self%img%new([params%box,params%box,1],params%smpd,                 wthreads=.false.)
            call self%img_match%new([params%boxmatch,params%boxmatch,1],params%smpd, wthreads=.false.)
            call self%img_copy%new([params%box,params%box,1],params%smpd,  wthreads=.false.)
            DebugPrint   'did build box-sized image objects'
            ! for thread safety in the image class
            call self%img%construct_thread_safe_tmp_imgs(params%nthr)
            ! boxmatch-sized ones
            call self%img_tmp%new([params%boxmatch,params%boxmatch,1],params%smpd,   wthreads=.false.)
            call self%img_msk%new([params%boxmatch,params%boxmatch,1],params%smpd,   wthreads=.false.)
            call self%mskimg%new([params%boxmatch, params%boxmatch, 1],params%smpd,  wthreads=.false.)
            DebugPrint  'did build boxmatch-sized image objects'
            ! boxpd-sized ones
            call self%img_pad%new([params%boxpd,params%boxpd,1],params%smpd)
            if( ddo3d )then
                call self%vol%new([params%box,params%box,params%box], params%smpd)
                call self%vol2%new([params%box,params%box,params%box], params%smpd)
            endif
            DebugPrint  'did build boxpd-sized image objects'
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
            DebugPrint   'did set default values'
        endif
        if( params%projstats .eq. 'yes' )then
            if( .not. self%spproj_field%isthere('proj') ) call self%spproj_field%set_projs(self%eulspace)
        endif
        write(*,'(A)') '>>> DONE BUILDING GENERAL TOOLBOX'
        self%general_tbox_exists = .true.
    end subroutine build_general_tbox

    subroutine kill_general_tbox( self )
        class(builder), intent(inout)  :: self
        if( self%general_tbox_exists )then
            call self%pgrpsyms%kill
            call self%spproj_field%kill
            self%spproj_field => null()
            call self%spproj%kill
            call self%eulspace%kill
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
            if( allocated(self%fsc)  ) deallocate(self%fsc)
            self%general_tbox_exists = .false.
        endif
    end subroutine kill_general_tbox

    subroutine build_comlin_tbox( self, params )
        class(builder),    intent(inout) :: self
        class(parameters), intent(inout) :: params
        integer :: i
        call self%kill_comlin_tbox
        if( params%pgrp /= 'c1' )then ! set up symmetry functionality
            DebugPrint 'build_comlin_tbox:  set up symmetry functionality'
            ! make object for symmetrized orientations
            call self%spproj_field%symmetrize(params%nsym)
            DebugPrint 'build_comlin_tbox:  symmetrized orientations obj made'
            allocate( self%imgs_sym(1:params%nsym*params%nptcls), self%ref_imgs(params%nstates,params%nspace), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk( 'build_comlin_tbox; simple_builder, 1', alloc_stat)
            DebugPrint 'build_comlin_tbox: allocated imgs '
            do i=1,params%nptcls*params%nsym
               call self%imgs_sym(i)%new([params%box,params%box,1],params%smpd)
            end do
            DebugPrint 'build_comlin_tbox: sym ptcls created '
            self%clins = comlin(self%spproj_field, self%imgs_sym, params%lp)
            DebugPrint 'build_comlin_tbox: comlin called '
        else ! set up assymetrical common lines-based alignment functionality
            DebugPrint 'build_comlin_tbox:  set up assymetrical common lines mode'
            allocate( self%imgs(1:params%nptcls), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk( 'build_comlin_tbox; simple_builder, 2', alloc_stat)
            DebugPrint 'build_comlin_tbox: allocated imgs '
            do i=1,params%nptcls
                call self%imgs(i)%new([params%box,params%box,1],params%smpd)
            end do
            DebugPrint 'build_comlin_tbox: sym ptcls created '
            self%clins = comlin( self%spproj_field, self%imgs, params%lp )
            DebugPrint 'build_comlin_tbox: comlin called '
        endif
        write(*,'(A)') '>>> DONE BUILDING COMLIN TOOLBOX'
        self%comlin_tbox_exists = .true.
    end subroutine build_comlin_tbox

    subroutine kill_comlin_tbox( self )
        class(builder), intent(inout) :: self
        integer :: i,j
        if( self%comlin_tbox_exists )then
            call self%spproj_field%kill
            if( allocated(self%imgs_sym) )then
                do i=1,size(self%imgs_sym)
                    call self%imgs_sym(i)%kill
                end do
                deallocate(self%imgs_sym)
            endif
            if( allocated(self%ref_imgs) )then
                do i=1,size(self%ref_imgs,1)
                    do j=1,size(self%ref_imgs,2)
                        call self%ref_imgs(i,j)%kill
                    end do
                end do
                deallocate(self%ref_imgs)
            endif
            if( allocated(self%imgs) )then
                do i=1,size(self%imgs)
                    call self%imgs(i)%kill
                end do
                deallocate(self%imgs)
            endif
            call self%clins%kill
            self%comlin_tbox_exists = .false.
        endif
    end subroutine kill_comlin_tbox

    subroutine build_rec_tbox( self, params )
        class(builder),    intent(inout) :: self
        class(parameters), intent(inout) :: params
        call self%kill_rec_tbox
        call self%recvol%new([params%boxpd,params%boxpd,params%boxpd],params%smpd)
        call self%recvol%alloc_rho( self%spproj)
        if( .not. self%spproj_field%isthere('proj') ) call self%spproj_field%set_projs(self%eulspace)
        write(*,'(A)') '>>> DONE BUILDING RECONSTRUCTION TOOLBOX'
        self%rec_tbox_exists = .true.
    end subroutine build_rec_tbox

    subroutine kill_rec_tbox( self )
        class(builder), intent(inout) :: self
        if( self%rec_tbox_exists )then
            call self%recvol%dealloc_rho
            call self%recvol%kill
            self%rec_tbox_exists = .false.
        endif
    end subroutine kill_rec_tbox

    subroutine build_rec_eo_tbox( self, params )
        class(builder),    intent(inout) :: self
        class(parameters), intent(inout) :: params
        call self%kill_rec_eo_tbox
        call self%eorecvol%new(self%spproj)
        if( .not. self%spproj_field%isthere('proj') ) call self%spproj_field%set_projs(self%eulspace)
        call self%projfrcs%new(NSPACE_REDUCED, params%box, params%smpd, params%nstates)
        write(*,'(A)') '>>> DONE BUILDING EO RECONSTRUCTION TOOLBOX'
        self%eo_rec_tbox_exists = .true.
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
        class(builder),    intent(inout) :: self
        class(parameters), intent(inout) :: params
        call self%kill_strategy2D_tbox
        if( params%neigh.eq.'yes' )then
            if( self%spproj%os_cls3D%get_noris() == params%ncls )then
                call self%spproj%os_cls3D%nearest_proj_neighbors(params%nnn, self%nnmat)
            else
               call simple_stop('simple_builder::build_strategy2D_tbox size of os_cls3D segment of spproj does not conform with # clusters (ncls)')
            endif
        endif
        call self%projfrcs%new(params%ncls, params%box, params%smpd, params%nstates)
        write(*,'(A)') '>>> DONE BUILDING STRATEGY2D TOOLBOX'
        self%strategy2D_tbox_exists = .true.
    end subroutine build_strategy2D_tbox

    subroutine kill_strategy2D_tbox( self )
        class(builder), intent(inout) :: self
        if( self%strategy2D_tbox_exists )then
            call self%projfrcs%kill
            self%strategy2D_tbox_exists = .false.
        endif
    end subroutine kill_strategy2D_tbox

    subroutine build_strategy3D_tbox( self, params )
        class(builder),    intent(inout) :: self
        class(parameters), intent(inout) :: params
        integer :: nnn
        call self%kill_strategy3D_tbox
        if( params%eo .ne. 'no' )then
            allocate( self%eorecvols(params%nstates), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('build_strategy3D_tbox; simple_builder, 1', alloc_stat)
        else
            allocate( self%recvols(params%nstates), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('build_strategy3D_tbox; simple_builder, 2', alloc_stat)
        endif
        if( params%neigh.eq.'yes' ) then
            nnn = params%nnn
            call self%pgrpsyms%nearest_proj_neighbors(self%eulspace, nnn, self%nnmat)
        endif
        if( .not. self%spproj_field%isthere('proj') ) call self%spproj_field%set_projs(self%eulspace)
        call self%projfrcs%new(NSPACE_REDUCED, params%box, params%smpd, params%nstates)
        write(*,'(A)') '>>> DONE BUILDING STRATEGY3D TOOLBOX'
        self%strategy3D_tbox_exists = .true.
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
            if( allocated(self%recvols) )then
                do i=1,size(self%recvols)
                    call self%recvols(i)%dealloc_rho
                    call self%recvols(i)%kill
                end do
                deallocate(self%recvols)
            endif
            if( allocated(self%nnmat) ) deallocate(self%nnmat)
            call self%projfrcs%kill
            self%strategy3D_tbox_exists = .false.
        endif
    end subroutine kill_strategy3D_tbox

    subroutine build_extremal3D_tbox( self, params )
        class(builder),    intent(inout) :: self
        class(parameters), intent(inout) :: params
        call self%kill_extremal3D_tbox
        allocate( self%recvols(1), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('build_strategy3D_tbox; simple_builder, 2', alloc_stat)
        call self%recvols(1)%new([params%boxpd,params%boxpd,params%boxpd],params%smpd)
        call self%recvols(1)%alloc_rho(self%spproj)
        if( .not. self%spproj_field%isthere('proj') ) call self%spproj_field%set_projs(self%eulspace)
        write(*,'(A)') '>>> DONE BUILDING EXTREMAL3D TOOLBOX'
        self%extremal3D_tbox_exists = .true.
    end subroutine build_extremal3D_tbox

    subroutine kill_extremal3D_tbox( self )
        class(builder), intent(inout) :: self
        if( self%extremal3D_tbox_exists )then
            call self%recvols(1)%dealloc_rho
            call self%recvols(1)%kill
            deallocate(self%recvols)
            self%extremal3D_tbox_exists = .false.
        endif
    end subroutine kill_extremal3D_tbox

end module simple_builder
