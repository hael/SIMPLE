! centralised builder (the main object constructor in SIMPLE)
module simple_build
include 'simple_lib.f08'
use simple_cmdline,          only: cmdline
use simple_params,           only: params
use simple_comlin,           only: comlin
use simple_image,            only: image
use simple_sp_project,       only: sp_project
use simple_oris,             only: oris
use simple_reconstructor,    only: reconstructor
use simple_reconstructor_eo, only: reconstructor_eo
use simple_sym,              only: sym
use simple_convergence,      only: convergence
use simple_projector,        only: projector
use simple_polarizer,        only: polarizer
use simple_masker,           only: masker
use simple_projection_frcs,  only: projection_frcs
implicit none

public :: build
private
#include "simple_local_flags.inc"

type :: build
    ! GENERAL TOOLBOX
    type(sp_project)                    :: spproj             !< centralised single-particle project meta-data handler
    class(oris), pointer                :: a => null()        !< pointer to field in spproj
    type(oris)                          :: e, e_bal           !< discrete spaces
    type(sym)                           :: se                 !< symmetry elements object
    type(convergence)                   :: conv               !< object for convergence checking of the PRIME2D/3D approaches
    type(image)                         :: img                !< individual image/projector objects
    type(polarizer)                     :: img_match          !< -"-
    type(image)                         :: img_pad            !< -"-
    type(image)                         :: img_tmp            !< -"-
    type(image)                         :: img_msk            !< -"-
    type(image)                         :: img_copy           !< -"-
    type(projector)                     :: vol, vol_odd
    type(image)                         :: vol2               !< -"-
    type(masker)                        :: mskimg             !< mask image
    type(projection_frcs)               :: projfrcs           !< projection FRC's used in the anisotropic Wiener filter
    type(image),            allocatable :: imgbatch(:)        !< batch of images
    ! COMMON LINES TOOLBOX
    type(image),            allocatable :: imgs(:)            !< images (all should be read in)
    type(image),            allocatable :: imgs_sym(:)        !< images (all should be read in)
    type(comlin)                        :: clins              !< common lines data structure
    type(image),            allocatable :: ref_imgs(:,:)      !< array of reference images
    ! RECONSTRUCTION TOOLBOX
    type(reconstructor_eo)              :: eorecvol           !< object for eo reconstruction
    type(reconstructor)                 :: recvol             !< object for reconstruction
    ! PRIME TOOLBOX
    type(reconstructor),    allocatable :: recvols(:)         !< array of volumes for reconstruction
    type(reconstructor_eo), allocatable :: eorecvols(:)       !< array of volumes for eo-reconstruction
    real,                   allocatable :: fsc(:,:)           !< Fourier Shell Correlation
    integer,                allocatable :: nnmat(:,:)         !< matrix with nearest neighbor indices
    integer,                allocatable :: grid_projs(:)      !< projection directions for coarse grid search
    ! PRIVATE EXISTENCE VARIABLES
    logical, private                    :: general_tbox_exists          = .false.
    logical, private                    :: cluster_tbox_exists          = .false.
    logical, private                    :: comlin_tbox_exists           = .false.
    logical, private                    :: rec_tbox_exists              = .false.
    logical, private                    :: eo_rec_tbox_exists           = .false.
    logical, private                    :: hadamard_prime3D_tbox_exists = .false.
    logical, private                    :: hadamard_prime2D_tbox_exists = .false.
    logical, private                    :: extremal3D_tbox_exists       = .false.
    logical, private                    :: read_features_exists         = .false.
  contains
    procedure                           :: build_spproj
    procedure                           :: build_general_tbox
    procedure                           :: kill_general_tbox
    procedure                           :: build_comlin_tbox
    procedure                           :: kill_comlin_tbox
    procedure                           :: build_rec_tbox
    procedure                           :: kill_rec_tbox
    procedure                           :: build_rec_eo_tbox
    procedure                           :: kill_rec_eo_tbox
    procedure                           :: build_strategy3D_tbox
    procedure                           :: kill_strategy3D_tbox
    procedure                           :: build_strategy2D_tbox
    procedure                           :: kill_strategy2D_tbox
    procedure                           :: build_extremal3D_tbox
    procedure                           :: kill_extremal3D_tbox
    procedure                           :: raise_hard_ctf_exception
end type build

contains

    !> \brief  constructs the sp project (part of general builder)
    subroutine build_spproj( self, p, cline, nooritab )
        use simple_binoris_io, only: binread_ctfparams_state_eo, binread_oritab
        use simple_user_interface,   only: get_prg_ptr
        class(build), target, intent(inout) :: self
        class(params),        intent(inout) :: p
        class(cmdline),       intent(inout) :: cline
        logical, optional,    intent(in)    :: nooritab
        logical ::  metadata_read
        ! create object for orientations
        ! b%a is now a pointer to a field in b%spproj
        select case(p%spproj_a_seg)
            case(MIC_SEG)
                call self%spproj%os_mic%new_clean(p%nptcls)
                self%a => self%spproj%os_mic
            case(STK_SEG)
                call self%spproj%os_stk%new_clean(p%nptcls)
                self%a => self%spproj%os_stk
            case(PTCL2D_SEG)
                call self%spproj%os_ptcl2D%new_clean(p%nptcls)
                self%a => self%spproj%os_ptcl2D
            case(CLS2D_SEG)
                call self%spproj%os_cls2D%new_clean(p%nptcls)
                self%a => self%spproj%os_cls2D
            case(CLS3D_SEG)
                call self%spproj%os_cls3D%new_clean(p%nptcls)
                self%a => self%spproj%os_cls3D
            case(PTCL3D_SEG)
                call self%spproj%os_ptcl3D%new_clean(p%nptcls)
                self%a => self%spproj%os_ptcl3D
            case DEFAULT
                ! using ptcl3D as the generic segment
                call self%spproj%os_ptcl3D%new_clean(p%nptcls)
                self%a => self%spproj%os_ptcl3D
        end select
        ! read from project file
        metadata_read = .false.
        if( cline%defined('projfile') )then
            if( .not. file_exists(p%projfile) )then
                write(*,*) 'project file not in exec dir: ', trim(p%projfile)
                stop 'ERROR! build :: build_spproj'
            endif
            metadata_read = .true.
        else if( associated(p%ptr2prg) )then
            if( p%ptr2prg%requires_sp_project() ) metadata_read = .true.
        endif
        if( metadata_read )then
            call self%spproj%read(p%projfile)
            ! ctf planning
            select case(trim(p%oritype))
                case('mic')
                    ! nothing to do
                case DEFAULT
                    p%tfplan%mode = self%spproj%get_ctfmode(p%oritype)
            end select
            ! update cwd of project (in case the params class changed exec dir)
            call self%spproj%projinfo%set(1, 'cwd', trim(p%cwd))
        else
            ! revert to oldschool logic
            if( present(nooritab) )then
                call self%a%spiral(p%nsym, p%eullims)
            else
                ! we need the oritab to override the deftab in order not to loose parameters
                if( p%deftab /= '' ) call binread_ctfparams_state_eo(p%deftab,  self%spproj, self%a, [1,p%nptcls])
                if( p%oritab /= '' ) call binread_oritab(p%oritab,              self%spproj, self%a, [1,p%nptcls])
                DebugPrint 'read deftab'
            endif
            ! ctf planning
            if( self%a%isthere('dfx') .and. self%a%isthere('dfy'))then
                p%tfplan%mode = 'astig'
            else if( self%a%isthere('dfx') )then
                p%tfplan%mode = 'noastig'
            else
                p%tfplan%mode = 'no'
            endif
            if( p%tfplan%flag .ne. 'no' .and. p%tfplan%mode .eq. 'no' )then
                write(*,'(a)') 'WARNING! It looks like you want to do Wiener restoration (p%ctf .ne. no)'
                write(*,'(a)') 'but your input orientation table lacks defocus values'
            endif
        endif
        write(*,'(A)') '>>> DONE BUILDING SP PROJECT'
    end subroutine build_spproj

    !> \brief  constructs the general toolbox
    subroutine build_general_tbox( self, p, cline, do3d, nooritab, force_ctf )
        class(build), target, intent(inout) :: self
        class(params),        intent(inout) :: p
        class(cmdline),       intent(inout) :: cline
        logical, optional,    intent(in)    :: do3d, nooritab, force_ctf
        integer :: lfny,  lfny_match, cyc_lims(3,2)
        logical :: ddo3d, fforce_ctf
        call self%kill_general_tbox
        ddo3d = .true.
        if( present(do3d) ) ddo3d = do3d
        fforce_ctf = .false.
        if( present(force_ctf) ) fforce_ctf = force_ctf
        ! seed the random number generator
        call seed_rnd
        DebugPrint   'seeded random number generator'
        ! set up symmetry functionality
        call self%se%new(trim(p%pgrp))
        p%nsym    = self%se%get_nsym()
        p%eullims = self%se%srchrange()
        DebugPrint   'did setup symmetry functionality'
        ! build spproj
        call self%build_spproj(p, cline, nooritab)
        ! states exception
        if( self%a%get_n('state') > 1 )then
            if( .not. cline%defined('nstates') )then
                write(*,'(a)') 'WARNING, your input doc has multiple states but NSTATES is not given'
            endif
        endif
        DebugPrint 'created & filled object for orientations'
        if( fforce_ctf ) call self%raise_hard_ctf_exception(p)
        ! generate discrete projection direction spaces
        call self%e%new_clean( p%nspace )
        call self%e%spiral( p%nsym, p%eullims )
        call self%e_bal%new_clean(NSPACE_BALANCE)
        call self%e_bal%spiral( p%nsym, p%eullims )
        ! create angular subspace
        self%grid_projs = self%e%create_proj_subspace(NPDIRS_SUBSPACE, p%nsym, p%eullims)
        ! store angular resolution of search space (ares)
        p%ares = self%e%find_angres()
        DebugPrint 'generated discrete projection direction space'
        if( p%box > 0 )then
            ! build image objects
            ! box-sized ones
            call self%img%new([p%box,p%box,1],p%smpd,                 wthreads=.false.)
            call self%img_match%new([p%boxmatch,p%boxmatch,1],p%smpd, wthreads=.false.)
            call self%img_copy%new([p%box,p%box,1],p%smpd,  wthreads=.false.)
            DebugPrint   'did build box-sized image objects'
            ! for thread safety in the image class
            call self%img%construct_thread_safe_tmp_imgs(p%nthr)
            ! boxmatch-sized ones
            call self%img_tmp%new([p%boxmatch,p%boxmatch,1],p%smpd,   wthreads=.false.)
            call self%img_msk%new([p%boxmatch,p%boxmatch,1],p%smpd,   wthreads=.false.)
            call self%mskimg%new([p%boxmatch, p%boxmatch, 1],p%smpd,  wthreads=.false.)
            DebugPrint  'did build boxmatch-sized image objects'
            ! boxpd-sized ones
            call self%img_pad%new([p%boxpd,p%boxpd,1],p%smpd)
            if( ddo3d )then
                call self%vol%new([p%box,p%box,p%box], p%smpd)
                call self%vol2%new([p%box,p%box,p%box], p%smpd)
            endif
            DebugPrint  'did build boxpd-sized image objects'
            ! build arrays
            lfny       = self%img%get_lfny(1)
            cyc_lims   = self%img_pad%loop_lims(3)
            allocate( self%fsc(p%nstates,lfny), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk("In: build_general_tbox; simple_build, 1", alloc_stat)
            self%fsc  = 0.
            ! set default amsklp
            if( .not. cline%defined('amsklp') .and. cline%defined('lp') )then
                p%amsklp = self%img%get_lp(self%img%get_find(p%lp)-2)
            endif
            DebugPrint   'did set default values'
        endif
        self%conv = convergence(self%a, p, cline)
        if( p%projstats .eq. 'yes' )then
            if( .not. self%a%isthere('proj') ) call self%a%set_projs(self%e)
        endif
        write(*,'(A)') '>>> DONE BUILDING GENERAL TOOLBOX'
        self%general_tbox_exists = .true.
    end subroutine build_general_tbox

    !> \brief  destructs the general toolbox
    subroutine kill_general_tbox( self )
        class(build), intent(inout)  :: self
        if( self%general_tbox_exists )then
            call self%conv%kill
            call self%se%kill
            call self%a%kill
            call self%e%kill
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

    !> \brief  constructs the common lines toolbox
    subroutine build_comlin_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        integer :: i
        call self%kill_comlin_tbox
        if( p%pgrp /= 'c1' )then ! set up symmetry functionality
            DebugPrint 'build_comlin_tbox:  set up symmetry functionality'
            ! make object for symmetrized orientations
            call self%a%symmetrize(p%nsym)
            DebugPrint 'build_comlin_tbox:  symmetrized orientations obj made'
            allocate( self%imgs_sym(1:p%nsym*p%nptcls), self%ref_imgs(p%nstates,p%nspace), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk( 'build_comlin_tbox; simple_build, 1', alloc_stat)
            DebugPrint 'build_comlin_tbox: allocated imgs '
            do i=1,p%nptcls*p%nsym
               call self%imgs_sym(i)%new([p%box,p%box,1],p%smpd)
            end do
            DebugPrint 'build_comlin_tbox: sym ptcls created '
            self%clins = comlin(self%a, self%imgs_sym, p%lp)
            DebugPrint 'build_comlin_tbox: comlin called '
        else ! set up assymetrical common lines-based alignment functionality
            DebugPrint 'build_comlin_tbox:  set up assymetrical common lines mode'
            allocate( self%imgs(1:p%nptcls), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk( 'build_comlin_tbox; simple_build, 2', alloc_stat)
            DebugPrint 'build_comlin_tbox: allocated imgs '
            do i=1,p%nptcls
                call self%imgs(i)%new([p%box,p%box,1],p%smpd)
            end do
            DebugPrint 'build_comlin_tbox: sym ptcls created '
            self%clins = comlin( self%a, self%imgs, p%lp )
            DebugPrint 'build_comlin_tbox: comlin called '
        endif
        write(*,'(A)') '>>> DONE BUILDING COMLIN TOOLBOX'
        self%comlin_tbox_exists = .true.
    end subroutine build_comlin_tbox

    !> \brief  destructs the common lines toolbox
    subroutine kill_comlin_tbox( self )
        class(build), intent(inout) :: self
        integer :: i,j
        if( self%comlin_tbox_exists )then
            call self%a%kill
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

    !> \brief  constructs the reconstruction toolbox
    subroutine build_rec_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        call self%kill_rec_tbox
        call self%raise_hard_ctf_exception(p)
        call self%recvol%new([p%boxpd,p%boxpd,p%boxpd],p%smpd)
        call self%recvol%alloc_rho(p, self%spproj)
        if( .not. self%a%isthere('proj') ) call self%a%set_projs(self%e)
        write(*,'(A)') '>>> DONE BUILDING RECONSTRUCTION TOOLBOX'
        self%rec_tbox_exists = .true.
    end subroutine build_rec_tbox

    !> \brief  destructs the reconstruction toolbox
    subroutine kill_rec_tbox( self )
        class(build), intent(inout) :: self
        if( self%rec_tbox_exists )then
            call self%recvol%dealloc_rho
            call self%recvol%kill
            self%rec_tbox_exists = .false.
        endif
    end subroutine kill_rec_tbox

    !> \brief  constructs the eo reconstruction toolbox
    subroutine build_rec_eo_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        call self%kill_rec_eo_tbox
        call self%raise_hard_ctf_exception(p)
        call self%eorecvol%new(p, self%spproj)
        if( .not. self%a%isthere('proj') ) call self%a%set_projs(self%e)
        call self%projfrcs%new(NSPACE_BALANCE, p%box, p%smpd, p%nstates)
        write(*,'(A)') '>>> DONE BUILDING EO RECONSTRUCTION TOOLBOX'
        self%eo_rec_tbox_exists = .true.
    end subroutine build_rec_eo_tbox

    !> \brief  destructs the eo reconstruction toolbox
    subroutine kill_rec_eo_tbox( self )
        class(build), intent(inout) :: self
        if( self%eo_rec_tbox_exists )then
            call self%eorecvol%kill
            call self%projfrcs%kill
            self%eo_rec_tbox_exists = .false.
        endif
    end subroutine kill_rec_eo_tbox

    !> \brief  constructs the prime2D toolbox
    subroutine build_strategy2D_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(inout) :: p
        call self%kill_strategy2D_tbox
        call self%raise_hard_ctf_exception(p)
        if( p%neigh.eq.'yes' )then
            if( self%spproj%os_cls3D%get_noris() == p%ncls )then
                call self%spproj%os_cls3D%nearest_proj_neighbors(p%nnn, self%nnmat)
            else
               call simple_stop('simple_build::build_hadamard_prime2D_tbox size of os_cls3D segment of spproj does not conform with # clusters (ncls)')
            endif
        endif
        call self%projfrcs%new(p%ncls, p%box, p%smpd, p%nstates)
        write(*,'(A)') '>>> DONE BUILDING HADAMARD PRIME2D TOOLBOX'
        self%hadamard_prime2D_tbox_exists = .true.
    end subroutine build_strategy2D_tbox

    !> \brief  destructs the prime2D toolbox
    subroutine kill_strategy2D_tbox( self )
        class(build), intent(inout) :: self
        if( self%hadamard_prime2D_tbox_exists )then
            call self%projfrcs%kill
            self%hadamard_prime2D_tbox_exists = .false.
        endif
    end subroutine kill_strategy2D_tbox

    !> \brief  constructs the prime3D toolbox
    subroutine build_strategy3D_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        integer :: nnn
        call self%kill_strategy3D_tbox
        call self%raise_hard_ctf_exception(p)
        if( p%eo .ne. 'no' )then
            allocate( self%eorecvols(p%nstates), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('build_strategy3D_tbox; simple_build, 1', alloc_stat)
        else
            allocate( self%recvols(p%nstates), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('build_strategy3D_tbox; simple_build, 2', alloc_stat)
        endif
        if( p%neigh.eq.'yes' ) then
            nnn = p%nnn
            call self%se%nearest_proj_neighbors(self%e, nnn, self%nnmat)
        endif
        if( .not. self%a%isthere('proj') ) call self%a%set_projs(self%e)
        call self%projfrcs%new(NSPACE_BALANCE, p%box, p%smpd, p%nstates)
        write(*,'(A)') '>>> DONE BUILDING HADAMARD PRIME3D TOOLBOX'
        self%hadamard_prime3D_tbox_exists = .true.
    end subroutine build_strategy3D_tbox

    !> \brief  destructs the prime3D toolbox
    subroutine kill_strategy3D_tbox( self )
        class(build), intent(inout) :: self
        integer :: i
        if( self%hadamard_prime3D_tbox_exists )then
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
            self%hadamard_prime3D_tbox_exists = .false.
        endif
    end subroutine kill_strategy3D_tbox

    !> \brief  constructs the extremal3D toolbox
    subroutine build_extremal3D_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        call self%kill_extremal3D_tbox
        call self%raise_hard_ctf_exception(p)
        allocate( self%recvols(1), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('build_strategy3D_tbox; simple_build, 2', alloc_stat)
        call self%recvols(1)%new([p%boxpd,p%boxpd,p%boxpd],p%smpd)
        call self%recvols(1)%alloc_rho(p, self%spproj)
        if( .not. self%a%isthere('proj') ) call self%a%set_projs(self%e)
        write(*,'(A)') '>>> DONE BUILDING EXTREMAL3D TOOLBOX'
        self%extremal3D_tbox_exists = .true.
    end subroutine build_extremal3D_tbox

    !> \brief  destructs the toolbox for continuous refinement
    subroutine kill_extremal3D_tbox( self )
        class(build), intent(inout) :: self
        if( self%extremal3D_tbox_exists )then
            call self%recvols(1)%dealloc_rho
            call self%recvols(1)%kill
            deallocate(self%recvols)
            self%extremal3D_tbox_exists = .false.
        endif
    end subroutine kill_extremal3D_tbox

    !> \brief  fall-over if CTF params are missing
    subroutine raise_hard_ctf_exception( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        if( p%tfplan%flag.ne.'no' )then
            select case(p%spproj_a_seg)
                case(MIC_SEG)
                    !
                case(STK_SEG)
                    if( .not. self%a%isthere('cs') ) write(*,*) 'ERROR! ctf .ne. no and input doc lacks cs'
                    if( .not. self%a%isthere('kv') ) write(*,*) 'ERROR! ctf .ne. no and input doc lacks kv'
                    if( .not. self%a%isthere('fraca') ) write(*,*) 'ERROR! ctf .ne. no and input doc lacks fraca'
                case(PTCL2D_SEG,PTCL3D_SEG)
                    if( .not. self%a%isthere('dfx') ) write(*,*) 'ERROR! ctf .ne. no and input doc lacks defocus values'
                case(CLS2D_SEG)
                    !
                case(CLS3D_SEG)
                    !
                case DEFAULT
                    !
            end select
        endif
        if( p%tfplan%l_phaseplate )then
            if( .not. self%a%isthere('phshift') )then
                write(*,*) 'ERROR! l_phaseplate = .true. and input doc lacks radian phshift'
                stop
            endif
        endif
    end subroutine raise_hard_ctf_exception

end module simple_build
