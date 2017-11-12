! centralised builder (the main object constructor in SIMPLE)

module simple_build
#include "simple_lib.f08"
!!import classes
use simple_cmdline,          only: cmdline
use simple_comlin,           only: comlin
use simple_image,            only: image
use simple_oris,             only: oris
use simple_reconstructor,    only: reconstructor
use simple_eo_reconstructor, only: eo_reconstructor
use simple_params,           only: params
use simple_sym,              only: sym
use simple_opt_spec,         only: opt_spec
use simple_convergence,      only: convergence
use simple_projector,        only: projector
use simple_polarizer,        only: polarizer
use simple_masker,           only: masker
use simple_projection_frcs,  only: projection_frcs
use simple_ran_tabu,         only: ran_tabu
!! import functions
use simple_timer,            only: tic, toc, timer_int_kind
use simple_binoris_io,       only: binread_ctfparams_state_eo, binread_oritab
implicit none

public :: build, test_build
private
#include "simple_local_flags.inc"
integer(timer_int_kind) :: tbuild

type :: build
    ! GENERAL TOOLBOX
    type(oris)                          :: a, e, e_bal        !< aligndata, discrete spaces
    type(sym)                           :: se                 !< symmetry elements object
    type(convergence)                   :: conv               !< object for convergence checking of the PRIME2D/3D approaches
    type(image)                         :: img                !< individual image/projector objects 
    type(polarizer)                     :: img_match          !< -"-
    type(image)                         :: img_pad            !< -"-
    type(image)                         :: img_tmp            !< -"-
    type(image)                         :: img_msk            !< -"-
    type(image)                         :: img_copy           !< -"-
    type(projector)                     :: vol                !< -"-
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
    type(eo_reconstructor)              :: eorecvol           !< object for eo reconstruction
    type(reconstructor)                 :: recvol             !< object for reconstruction
    ! PRIME TOOLBOX
    type(reconstructor),    allocatable :: recvols(:)         !< array of volumes for reconstruction
    type(eo_reconstructor), allocatable :: eorecvols(:)       !< array of volumes for eo-reconstruction
    complex(sp),            allocatable :: cmat(:,:)          !< 2D matrix for reading FTs prepped 4 cgrid
    real,                   allocatable :: fsc(:,:)           !< Fourier Shell Correlation
    integer,                allocatable :: nnmat(:,:)         !< matrix with nearest neighbor indices
    integer,                allocatable :: pbatch(:)          !< particle index batch
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
    procedure                           :: build_general_tbox
    procedure                           :: kill_general_tbox
    procedure                           :: build_comlin_tbox
    procedure                           :: kill_comlin_tbox
    procedure                           :: build_rec_tbox
    procedure                           :: kill_rec_tbox
    procedure                           :: build_eo_rec_tbox
    procedure                           :: kill_eo_rec_tbox
    procedure                           :: build_hadamard_prime3D_tbox
    procedure                           :: kill_hadamard_prime3D_tbox
    procedure                           :: build_hadamard_prime2D_tbox
    procedure                           :: kill_hadamard_prime2D_tbox
    procedure                           :: build_extremal3D_tbox
    procedure                           :: kill_extremal3D_tbox
    procedure                           :: raise_hard_ctf_exception
end type build

contains

    !> \brief  constructs the general toolbox
    subroutine build_general_tbox( self, p, cline, do3d, nooritab, force_ctf )
        class(build),      intent(inout) :: self
        class(params),     intent(inout) :: p
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: do3d, nooritab, force_ctf
        type(ran_tabu) :: rt
        integer        :: lfny, partsz, lfny_match, cyc_lims(3,2)
        logical        :: ddo3d, fforce_ctf
        verbose=.false.
        if(verbose.or.global_verbose) tbuild=tic()
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
        ! create object for orientations
        call self%a%new(p%nptcls)
        if( present(nooritab) )then
            call self%a%spiral(p%nsym, p%eullims)
        else
            ! we need the oritab to override the deftab in order not to loose parameters
            if( p%deftab /= '' ) call binread_ctfparams_state_eo(p%deftab, self%a, [1,p%nptcls])
            if( p%oritab /= '' )then
                if( .not. cline%defined('nstates') )then
                    if( p%nstates > 1 )then
                        print *,'Multiple states detected, please input the NSTATES key'
                        stop
                    endif
                    call binread_oritab(p%oritab, self%a, [1,p%nptcls], p%nstates)
                else
                    call binread_oritab(p%oritab, self%a, [1,p%nptcls])
                endif
            endif
        endif
        if( self%a%get_n('state') > 1 )then
            if( .not. cline%defined('nstates') )then
                write(*,'(a)') 'WARNING, your input doc has multiple states but NSTATES is not given'
            endif
        endif
        DebugPrint 'created & filled object for orientations'
        DebugPrint 'read deftab'
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
        DebugPrint 'did set number of dimensions and ctfmode'
        if( fforce_ctf ) call self%raise_hard_ctf_exception(p)
        ! generate discrete projection direction spaces
        call self%e%new( p%nspace )
        call self%e%spiral( p%nsym, p%eullims )
        call self%e_bal%new(NSPACE_BALANCE)
        call self%e_bal%spiral( p%nsym, p%eullims )
        self%grid_projs = self%e%create_proj_subspace(p%nsub, p%nsym, p%eullims )
        DebugPrint 'generated discrete projection direction space'
        if( p%box > 0 )then
            ! build image objects
            ! box-sized ones
            call self%img%new([p%box,p%box,1],p%smpd,                 wthreads=.false.)
            call self%img_match%new([p%boxmatch,p%boxmatch,1],p%smpd, wthreads=.false.)
            call self%img_copy%new([p%boxmatch,p%boxmatch,1],p%smpd,  wthreads=.false.)
            DebugPrint   'did build box-sized image objects'
            ! for thread safety in the image class
            call self%img%construct_thread_safe_tmp_imgs(p%nthr)
            ! boxmatch-sized ones
            call self%img_tmp%new([p%boxmatch,p%boxmatch,1],p%smpd,   wthreads=.false.)
            call self%img_msk%new([p%boxmatch,p%boxmatch,1],p%smpd,   wthreads=.false.)
            call self%mskimg%new([p%boxmatch, p%boxmatch, 1],p%smpd,  wthreads=.false.)
            DebugPrint  'did build boxmatch-sized image objects'
            ! boxpd-sized ones
            call self%img_pad%new([p%boxpd,p%boxpd,1],p%smpd,         wthreads=.false.)
            if( ddo3d )then
                call self%vol%new([p%box,p%box,p%box], p%smpd)
                call self%vol2%new([p%box,p%box,p%box], p%smpd)
            endif
            DebugPrint  'did build boxpd-sized image objects'
            ! build arrays
            lfny       = self%img%get_lfny(1)
            lfny_match = self%img_match%get_lfny(1)
            cyc_lims   = self%img_pad%loop_lims(3)
            allocate( self%fsc(p%nstates,lfny),&
                &self%cmat(cyc_lims(1,1):cyc_lims(1,2),cyc_lims(2,1):cyc_lims(2,2)), stat=alloc_stat )
            allocchk("In: build_general_tbox; simple_build, 1")
            self%fsc  = 0.
            self%cmat = cmplx(0.,0.)
            ! set record-lenght for direct access I/O
            inquire(iolength=p%recl_cgrid) self%cmat
            ! set default amsklp
            if( .not. cline%defined('amsklp') .and. cline%defined('lp') )then
                p%amsklp = self%img%get_lp(self%img%get_find(p%lp)-2)
            endif
            DebugPrint   'did set default values'
        endif
        ! build convergence checker
        self%conv = convergence(self%a, p, cline)
        ! generate random particle batch
        if( cline%defined('batchfrac') )then
            ! allocate index array
            partsz    = p%top - p%fromp + 1
            p%batchsz = nint(real(partsz) * p%batchfrac)
            allocate(self%pbatch(p%batchsz), stat=alloc_stat)
            allocchk("In: build_general_tbox; simple_build, 2")
            ! select random indices
            rt = ran_tabu(partsz)
            call rt%ne_ran_iarr(self%pbatch)
            ! shift the indicies
            self%pbatch = self%pbatch + p%fromp - 1
            call rt%kill
        endif
        if( p%projstats .eq. 'yes' )then
            if( .not. self%a%isthere('proj') ) call self%a%set_projs(self%e)
        endif
        if (verbose.or.global_verbose)then
            write(*,'(A,1x,1ES20.5)') '>>> DONE BUILDING GENERAL TOOLBOX             time (s)', toc(tbuild)
        else
            write(*,'(A)') '>>> DONE BUILDING GENERAL TOOLBOX'
        endif
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
            if( allocated(self%cmat) ) deallocate(self%cmat)
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
            ! make object for symmetrized orientations
            call self%a%symmetrize(p%nsym)
            allocate( self%imgs_sym(1:p%nsym*p%nptcls), self%ref_imgs(p%nstates,p%nspace), stat=alloc_stat )
            allocchk( 'build_comlin_tbox; simple_build, 1')
            do i=1,p%nptcls*p%nsym
                call self%imgs_sym(i)%new([p%box,p%box,1],p%smpd)
            end do
            self%clins = comlin(self%a, self%imgs_sym, p%lp)
        else ! set up assymetrical common lines-based alignment functionality
            allocate( self%imgs(1:p%nptcls), stat=alloc_stat )
            allocchk( 'build_comlin_tbox; simple_build, 2')
            do i=1,p%nptcls
                call self%imgs(i)%new([p%box,p%box,1],p%smpd)
            end do
            self%clins = comlin( self%a, self%imgs, p%lp )
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
        verbose=.false.
        if(verbose.or.global_verbose) tbuild=tic()
        call self%kill_rec_tbox
        call self%raise_hard_ctf_exception(p)
        call self%recvol%new([p%boxpd,p%boxpd,p%boxpd],p%smpd)
        call self%recvol%alloc_rho(p)
        if( .not. self%a%isthere('proj') ) call self%a%set_projs(self%e)
        if (verbose.or.global_verbose)then
            write(*,'(A,1x,1ES20.5)') '>>> DONE BUILDING RECONSTRUCTION TOOLBOX      time (s)', toc(tbuild)
        else
            write(*,'(A)') '>>> DONE BUILDING RECONSTRUCTION TOOLBOX'
        endif
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
    subroutine build_eo_rec_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        if(verbose.or.global_verbose) tbuild=tic()
        call self%kill_eo_rec_tbox
        call self%raise_hard_ctf_exception(p)
        call self%eorecvol%new(p)
        if( .not. self%a%isthere('proj') ) call self%a%set_projs(self%e)
        call self%projfrcs%new(NSPACE_BALANCE, p%box, p%smpd, p%nstates)
        if (verbose.or.global_verbose)then
            write(*,'(A,1x,1ES20.5)') '>>> DONE BUILDING EO RECONSTRUCTION TOOLBOX   time (s)', toc(tbuild)
        else
            write(*,'(A)') '>>> DONE BUILDING EO RECONSTRUCTION TOOLBOX'
        endif
        self%eo_rec_tbox_exists = .true.
    end subroutine build_eo_rec_tbox

    !> \brief  destructs the eo reconstruction toolbox
    subroutine kill_eo_rec_tbox( self )
        class(build), intent(inout) :: self
        if( self%eo_rec_tbox_exists )then
            call self%eorecvol%kill
            call self%projfrcs%kill
            self%eo_rec_tbox_exists = .false.
        endif
    end subroutine kill_eo_rec_tbox

    !> \brief  constructs the prime2D toolbox
    subroutine build_hadamard_prime2D_tbox( self, p )
        use simple_binoris_io, only: binread_oritab
        class(build),  intent(inout) :: self
        class(params), intent(inout) :: p
        type(oris) :: os
        if(verbose.or.global_verbose) tbuild=tic()
        call self%kill_hadamard_prime2D_tbox
        call self%raise_hard_ctf_exception(p)
        if( str_has_substr(p%refine,'neigh') )then
            if( file_exists(p%oritab3D) )then
                call os%new(p%ncls)
                call binread_oritab(p%oritab3D, os, [1,p%ncls])
                call os%nearest_proj_neighbors(p%nnn, self%nnmat)
                call os%kill
            else
                stop 'need oritab3D input for prime2D refine=neigh mode; simple_build :: build_hadamard_prime2D_tbox'
            endif
        endif
        ! build projection frcs
        call self%projfrcs%new(p%ncls, p%box, p%smpd, p%nstates)
        if (verbose.or.global_verbose)then
            write(*,'(A,1x,1ES20.5)') '>>> DONE BUILDING HADAMARD PRIME2D TOOLBOX    time (s)', toc(tbuild)
        else
            write(*,'(A)') '>>> DONE BUILDING HADAMARD PRIME2D TOOLBOX'
        endif
        self%hadamard_prime2D_tbox_exists = .true.
    end subroutine build_hadamard_prime2D_tbox

    !> \brief  destructs the prime2D toolbox
    subroutine kill_hadamard_prime2D_tbox( self )
        class(build), intent(inout) :: self
        if( self%hadamard_prime2D_tbox_exists )then
            call self%projfrcs%kill
            self%hadamard_prime2D_tbox_exists = .false.
        endif
    end subroutine kill_hadamard_prime2D_tbox

    !> \brief  constructs the prime3D toolbox
    subroutine build_hadamard_prime3D_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        integer :: nnn
        if(verbose.or.global_verbose) tbuild=tic()
        call self%kill_hadamard_prime3D_tbox
        call self%raise_hard_ctf_exception(p)
        ! reconstruction objects
        if( p%eo .ne. 'no' )then
            allocate( self%eorecvols(p%nstates), stat=alloc_stat )
            allocchk('build_hadamard_prime3D_tbox; simple_build, 1')
        else
            allocate( self%recvols(p%nstates), stat=alloc_stat )
            allocchk('build_hadamard_prime3D_tbox; simple_build, 2')
        endif
        if( str_has_substr(p%refine,'neigh') .or. trim(p%refine).eq.'states' )then
            nnn = p%nnn
            call self%se%nearest_proj_neighbors(self%e, nnn, self%nnmat)
        endif
        if( .not. self%a%isthere('proj') ) call self%a%set_projs(self%e)
        call self%projfrcs%new(NSPACE_BALANCE, p%box, p%smpd, p%nstates)
        if (verbose.or.global_verbose)then
            write(*,'(A,1x,1ES20.5)') '>>> DONE BUILDING HADAMARD PRIME3D TOOLBOX    time (s)', toc(tbuild)
        else
            write(*,'(A)') '>>> DONE BUILDING HADAMARD PRIME3D TOOLBOX'
        endif
        self%hadamard_prime3D_tbox_exists = .true.
    end subroutine build_hadamard_prime3D_tbox

    !> \brief  destructs the prime3D toolbox
    subroutine kill_hadamard_prime3D_tbox( self )
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
    end subroutine kill_hadamard_prime3D_tbox

    !> \brief  constructs the extremal3D toolbox
    subroutine build_extremal3D_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        if(verbose.or.global_verbose) tbuild=tic()
        call self%kill_extremal3D_tbox
        call self%raise_hard_ctf_exception(p)
        allocate( self%recvols(1), stat=alloc_stat )
        allocchk('build_hadamard_prime3D_tbox; simple_build, 2')
        call self%recvols(1)%new([p%boxpd,p%boxpd,p%boxpd],p%smpd)
        call self%recvols(1)%alloc_rho(p)
        if( .not. self%a%isthere('proj') ) call self%a%set_projs(self%e)
        if (verbose.or.global_verbose)then
            write(*,'(A,1x,1ES20.5)') '>>> DONE BUILDING EXTREMAL3D TOOLBOX          time (s)', toc(tbuild)
        else
            write(*,'(A)') '>>> DONE BUILDING EXTREMAL3D TOOLBOX'
        end if
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
        logical :: params_present(4)
        if( p%tfplan%flag.ne.'no' )then
            params_present(1) = self%a%isthere('kv')
            params_present(2) = self%a%isthere('cs')
            params_present(3) = self%a%isthere('fraca')
            params_present(4) = self%a%isthere('dfx')
            if( all(params_present) )then
                ! alles ok
            else
                if( .not. params_present(1) ) write(*,*) 'ERROR! ctf .ne. no and input doc lacks kv'
                if( .not. params_present(2) ) write(*,*) 'ERROR! ctf .ne. no and input doc lacks cs'
                if( .not. params_present(3) ) write(*,*) 'ERROR! ctf .ne. no and input doc lacks fraca'
                if( .not. params_present(4) ) write(*,*) 'ERROR! ctf .ne. no and input doc lacks defocus'
                stop
            endif
        endif
        if( p%tfplan%l_phaseplate )then
            if( .not. self%a%isthere('phshift') )then
                write(*,*) 'ERROR! l_phaseplate = .true. and input doc lacks radian phshift'
                stop
            endif
        endif
    end subroutine raise_hard_ctf_exception

    ! UNIT TEST

    !> \brief  build unit test
    subroutine test_build
        type(build)   :: myb
        type(cmdline) :: mycline_static, mycline_varying
        type(params)  :: myp
        write(*,'(a)') '**info(simple_build_unit_test): testing the different building options'
        ! setup command line
        call mycline_static%set('box',      100.)
        call mycline_static%set('msk',       40.)
        call mycline_static%set('smpd',       2.)
        call mycline_static%set('nvars',     40.)
        call mycline_static%set('nptcls', 10000.)
        call mycline_static%set('ncls',     100.)
        call mycline_static%set('nstates',    2.)
        write(*,'(a)') '**info(simple_build_unit_test): generated command line'
        ! 12 cases to test
        ! case 1:  refine=no, pgrp=c1, eo=yes
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=no, pgrp=c1, eo=yes'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 1 passed'
        ! case 2:  refine=no, pgrp=c1, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=no, pgrp=c1, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 2 passed'
        ! case 3:  refine=no, pgrp=c2, eo=yes
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=no, pgrp=c2, eo=yes'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 3 passed'
        ! case 4:  refine=no, pgrp=c2, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=no, pgrp=c2, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 4 passed'
        ! case 5:  refine=neigh, pgrp=c1, eo=yes
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c1, eo=yes'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 5 passed'
        ! case 6:  refine=neigh, pgrp=c1, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c1, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 6 passed'
        ! case 7:  refine=neigh, pgrp=c2, eo=yes
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c2, eo=yes'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 7 passed'
        ! case 8:  refine=neigh, pgrp=c2, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c2, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 8 passed'
        ! case 9:  refine=neigh, pgrp=c1, eo=yes
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c1, eo=yes'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 9 passed'
        ! case 10: refine=neigh, pgrp=c1, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c1, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        write(*,'(a)') '**info(simple_build_unit_test): case 10 passed'
        ! case 11: refine=neigh, pgrp=c2, eo=yes
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c2, eo=yes'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 11 passed'
        ! case 12: refine=neigh, pgrp=c2, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c2, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 12 passed'
        write(*,'(a)') 'SIMPLE_BUILD_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'

      contains

          subroutine tester
              call myb%build_general_tbox(myp, mycline_varying, do3d=.true., nooritab=.true.)
              call myb%build_general_tbox(myp, mycline_varying, do3d=.true., nooritab=.true.)
              call myb%kill_general_tbox
              write(*,'(a)') 'build_general_tbox passed'
              call myb%build_comlin_tbox(myp)
              call myb%build_comlin_tbox(myp)
              call myb%kill_comlin_tbox
              write(*,'(a)') 'build_comlin_tbox passed'
              call myb%build_rec_tbox(myp)
              call myb%build_rec_tbox(myp)
              call myb%kill_rec_tbox
              write(*,'(a)') 'build_rec_tbox passed'
              call myb%build_eo_rec_tbox(myp)
              call myb%build_eo_rec_tbox(myp)
              call myb%kill_eo_rec_tbox
              write(*,'(a)') 'build_eo_rec_tbox passed'
              call myb%build_hadamard_prime3D_tbox(myp)
              call myb%build_hadamard_prime3D_tbox(myp)
              call myb%kill_hadamard_prime3D_tbox
              write(*,'(a)') 'build_hadamard_prime3D_tbox passed'
              call myb%build_hadamard_prime2D_tbox(myp)
              call myb%build_hadamard_prime2D_tbox(myp)
              call myb%kill_hadamard_prime2D_tbox
              write(*,'(a)') 'build_hadamard_prime2D_tbox passed'
              call myb%build_extremal3D_tbox(myp)
              call myb%build_extremal3D_tbox(myp)
              call myb%kill_extremal3D_tbox
              write(*,'(a)') 'build_extremal3D_tbox passed'
          end subroutine tester

    end subroutine test_build

end module simple_build
