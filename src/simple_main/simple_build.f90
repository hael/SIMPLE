!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> simple_build is the builder class for the methods in _SIMPLE_. Access is global in the using unit
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution or modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund, 2009-06-11.
!
!==Changes are documented below
!
!* deugged and incorporated in the _SIMPLE_ library, HE 2009-06-25
!* reshaped according to the new simple_params class that will deal with all global parameters, HE 2011-08-18
!* rewritten with the new language features, HE 2012-06-18
!
module simple_build
use simple_defs
use simple_cmdline,             only: cmdline
use simple_comlin,              only: comlin
use simple_image,               only: image
use simple_centre_clust,        only: centre_clust
use simple_oris,                only: oris
use simple_pair_dtab,           only: pair_dtab
use simple_ppca,                only: ppca
use simple_reconstructor,       only: reconstructor
use simple_eo_reconstructor,    only: eo_reconstructor
use simple_params,              only: params
use simple_sym,                 only: sym
use simple_opt_spec,            only: opt_spec
use simple_convergence,         only: convergence
use simple_convergence_perptcl, only: convergence_perptcl
use simple_jiffys,              only: alloc_err
use simple_projector,           only: projector
use simple_polarizer,           only: polarizer
use simple_masker,              only: masker
use simple_filehandling         ! use all in there
implicit none

public :: build, test_build
private
#include "simple_local_flags.inc"

type build
    ! GENERAL TOOLBOX
    type(oris)                          :: a, e               !< aligndata, discrete space
    type(sym)                           :: se                 !< symmetry elements object
    type(convergence)                   :: conv               !< object for convergence checking of the PRIME2D/3D approaches
    type(convergence_perptcl)           :: ppconv             !< per-particle convergence checking object
    type(image)                         :: img                !< individual image objects
    type(polarizer)                     :: img_match          !< -"-
    type(image)                         :: img_pad            !< -"-
    type(image)                         :: img_tmp            !< -"-
    type(image)                         :: img_msk            !< -"-
    type(image)                         :: img_copy           !< -"-
    type(projector)                     :: vol                !< -"-
    type(projector)                     :: vol_pad            !< -"-
    type(masker)                        :: mskimg             !< mask image
    type(masker)                        :: mskvol             !< mask volume
    ! CLUSTER TOOLBOX
    type(ppca)                          :: pca                !< 4 probabilistic pca
    type(centre_clust)                  :: cenclust           !< centre-based clustering object
    real, allocatable                   :: features(:,:)      !< features for clustering
    ! COMMON LINES TOOLBOX
    type(image), allocatable            :: imgs(:)            !< images (all should be read in)
    type(image), allocatable            :: imgs_sym(:)        !< images (all should be read in)
    type(comlin)                        :: clins              !< common lines data structure
    type(image), allocatable            :: ref_imgs(:,:)      !< array of reference images
    ! RECONSTRUCTION TOOLBOX
    type(eo_reconstructor)              :: eorecvol           !< object for eo reconstruction
    type(reconstructor)                 :: recvol             !< object for reconstruction
    ! PRIME TOOLBOX
    type(image),            allocatable :: cavgs(:)           !< class averages (Wiener normalised references)
    type(image),            allocatable :: ctfsqsums(:)       !< CTF**2 sums for Wiener normalisation
    type(projector),        allocatable :: refvols(:)         !< reference volumes for quasi-continuous search
    type(reconstructor),    allocatable :: recvols(:)         !< array of volumes for reconstruction
    type(eo_reconstructor), allocatable :: eorecvols(:)       !< array of volumes for eo-reconstruction
    real,    allocatable                :: fsc(:,:)           !< Fourier shell correlation
    integer, allocatable                :: nnmat(:,:)         !< matrix with nearest neighbor indices
    integer, allocatable                :: pbatch(:)          !< particle index batch
    integer, allocatable                :: grid_projs(:)      !< projection directions for coarse grid search
    ! PRIVATE EXISTENCE VARIABLES
    logical, private                    :: general_tbox_exists          = .false.
    logical, private                    :: cluster_tbox_exists          = .false.
    logical, private                    :: comlin_tbox_exists           = .false.
    logical, private                    :: rec_tbox_exists              = .false.
    logical, private                    :: eo_rec_tbox_exists           = .false.
    logical, private                    :: hadamard_prime3D_tbox_exists = .false.
    logical, private                    :: hadamard_prime2D_tbox_exists = .false.
    logical, private                    :: cont3D_tbox_exists           = .false.
    logical, private                    :: extremal3D_tbox_exists       = .false.
    logical, private                    :: read_features_exists         = .false.
  contains
    procedure                           :: build_general_tbox
    procedure                           :: kill_general_tbox
    procedure                           :: build_cluster_tbox
    procedure                           :: kill_cluster_tbox
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
    procedure                           :: build_cont3D_tbox
    procedure                           :: kill_cont3D_tbox
    procedure                           :: build_extremal3D_tbox
    procedure                           :: kill_extremal3D_tbox
    procedure                           :: read_features
    procedure                           :: raise_hard_ctf_exception
end type build

contains

    !> \brief  constructs the general toolbox
    subroutine build_general_tbox( self, p, cline, do3d, nooritab, force_ctf )
        use simple_ran_tabu, only: ran_tabu
        use simple_math,     only: nvoxfind, rad2deg
        use simple_rnd,      only: seed_rnd
        use simple_strings,  only: str_has_substr
        class(build),      intent(inout) :: self
        class(params),     intent(inout) :: p
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: do3d, nooritab, force_ctf
        type(ran_tabu) :: rt
        integer        :: alloc_stat, lfny, partsz, lfny_match
        real           :: slask(3)
        logical        :: err, ddo3d, fforce_ctf
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
            if( p%deftab /= '' ) call self%a%read(p%deftab)
            if( p%oritab /= '' )then
                if( .not. cline%defined('nstates') )then
                    call self%a%read(p%oritab, p%nstates)
                    if( p%nstates>1 .and. str_has_substr(p%refine,'shc') )then
                        print *,'Multiple states detected, please input the NSTATES key'
                        stop
                    endif
                else
                    call self%a%read(p%oritab)
                endif
                if( self%a%get_noris() > 1 )then
                    call self%a%stats('corr', slask(1), slask(2), slask(3), err)
                    if( err )then
                    else
                        if( p%frac < 0.99 ) call self%a%calc_hard_ptcl_weights(p%var, bystate=.true.)
                    endif
                endif
            endif
        endif
        DebugPrint   'created & filled object for orientations'
        DebugPrint   'read deftab'
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
        DebugPrint   'did set number of dimensions and ctfmode'
        if( fforce_ctf ) call self%raise_hard_ctf_exception(p)
        ! generate discrete projection direction spaces
        call self%e%new( p%nspace )
        call self%e%spiral( p%nsym, p%eullims )
        self%grid_projs = self%e%create_proj_subspace(p%nsub, p%nsym, p%eullims )
        DebugPrint 'generated discrete projection direction space'
        if( p%box > 0 )then
            ! build image objects
            ! box-sized ones
            call self%img%new([p%box,p%box,1],p%smpd,p%imgkind)
            call self%img_match%new([p%boxmatch,p%boxmatch,1],p%smpd,p%imgkind)
            call self%img_copy%new([p%boxmatch,p%boxmatch,1],p%smpd,p%imgkind)
            DebugPrint   'did build box-sized image objects'
            ! boxmatch-sized ones
            call self%img_tmp%new([p%boxmatch,p%boxmatch,1],p%smpd,p%imgkind)
            call self%img_msk%new([p%boxmatch,p%boxmatch,1],p%smpd,p%imgkind)
            call self%mskimg%new([p%boxmatch, p%boxmatch, 1],p%smpd,p%imgkind)
            DebugPrint  'did build boxmatch-sized image objects'
            ! boxpd-sized ones
            call self%img_pad%new([p%boxpd,p%boxpd,1],p%smpd,p%imgkind)
            if( ddo3d )then
                call self%vol%new([p%box,p%box,p%box], p%smpd, p%imgkind)
                call self%vol_pad%new([p%boxpd,p%boxpd,p%boxpd],p%smpd,p%imgkind)
            endif
            DebugPrint  'did build boxpd-sized image objects'
            ! build arrays
            lfny = self%img%get_lfny(1)
            lfny_match = self%img_match%get_lfny(1)         
            allocate( self%fsc(p%nstates,lfny), source=0., stat=alloc_stat )
            call alloc_err("In: build_general_tbox; simple_build, 1", alloc_stat)
            ! set default amsklp
            if( .not. cline%defined('amsklp') .and. cline%defined('lp') )then
                p%amsklp = self%img%get_lp(self%img%get_find(p%lp)-2)
            endif
            DebugPrint   'did set default values'
        endif
        ! build convergence checkers
        self%conv   = convergence(self%a, p, cline)
        self%ppconv = convergence_perptcl(self%a, p, cline)
        ! generate random particle batch
        if( cline%defined('batchfrac') )then
            ! allocate index array
            partsz    = p%top - p%fromp + 1
            p%batchsz = nint(real(partsz) * p%batchfrac)
            allocate(self%pbatch(p%batchsz), stat=alloc_stat)
            call alloc_err("In: build_general_tbox; simple_build, 2", alloc_stat)
            ! select random indices
            rt = ran_tabu(partsz)
            call rt%ne_ran_iarr(self%pbatch)
            ! shift the indicies
            self%pbatch = self%pbatch + p%fromp - 1
            call rt%kill
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
            call self%mskvol%kill_masker
            call self%mskvol%kill
            call self%mskimg%kill_masker
            call self%mskimg%kill
            call self%vol_pad%kill_expanded
            call self%vol_pad%kill
            if( allocated(self%fsc) )deallocate(self%fsc)
            self%general_tbox_exists = .false.
        endif
    end subroutine kill_general_tbox

    !> \brief  constructs the cluster toolbox
    subroutine build_cluster_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(inout) :: p
        type(image)                  :: img
        integer                      :: alloc_stat
        call self%kill_cluster_tbox
        call img%new([p%box,p%box,1], p%smpd, p%imgkind)
        p%ncomps = img%get_npix(p%msk)
        call img%kill
        DebugPrint   'ncomps (npixels): ', p%ncomps
        DebugPrint   'nvars (nfeatures): ', p%nvars
        call self%pca%new(p%nptcls, p%ncomps, p%nvars) !< Dmat not pre-allocated (possible bug... /ME)
        call self%cenclust%new(self%features, self%a, p%nptcls, p%nvars, p%ncls)
        allocate( self%features(p%nptcls,p%nvars), stat=alloc_stat )
        call alloc_err('build_cluster_toolbox', alloc_stat)
        write(*,'(A)') '>>> DONE BUILDING CLUSTER TOOLBOX'
        self%cluster_tbox_exists = .true.
    end subroutine build_cluster_tbox

    !> \brief  destructs the cluster toolbox
    subroutine kill_cluster_tbox( self )
        class(build), intent(inout) :: self
        if( self%cluster_tbox_exists )then
            call self%pca%kill
            call self%cenclust%kill
            deallocate( self%features )
            self%cluster_tbox_exists = .false.
        endif
    end subroutine kill_cluster_tbox

    !> \brief  constructs the common lines toolbox
    subroutine build_comlin_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        integer :: alloc_stat, i
        call self%kill_comlin_tbox
        if( p%pgrp /= 'c1' )then ! set up symmetry functionality
            ! make object for symmetrized orientations
            call self%a%symmetrize(p%nsym)
            allocate( self%imgs_sym(1:p%nsym*p%nptcls), self%ref_imgs(p%nstates,p%nspace), stat=alloc_stat )
            call alloc_err( 'build_comlin_tbox; simple_build, 1', alloc_stat )
            do i=1,p%nptcls*p%nsym
                call self%imgs_sym(i)%new([p%box,p%box,1],p%smpd,p%imgkind)
            end do
            self%clins = comlin(self%a, self%imgs_sym, p%lp)
        else ! set up assymetrical common lines-based alignment functionality
            allocate( self%imgs(1:p%nptcls), stat=alloc_stat )
            call alloc_err( 'build_comlin_tbox; simple_build, 2', alloc_stat )
            do i=1,p%nptcls
                call self%imgs(i)%new([p%box,p%box,1],p%smpd,p%imgkind)
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
        call self%kill_rec_tbox
        call self%raise_hard_ctf_exception(p)
        call self%recvol%new([p%boxpd,p%boxpd,p%boxpd],p%smpd,p%imgkind)
        call self%recvol%alloc_rho(p)
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
    subroutine build_eo_rec_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        call self%kill_eo_rec_tbox
        call self%raise_hard_ctf_exception(p)
        call self%eorecvol%new(p)
        write(*,'(A)') '>>> DONE BUILDING EO RECONSTRUCTION TOOLBOX'
        self%eo_rec_tbox_exists = .true.
    end subroutine build_eo_rec_tbox

    !> \brief  destructs the eo reconstruction toolbox
    subroutine kill_eo_rec_tbox( self )
        class(build), intent(inout) :: self
        if( self%eo_rec_tbox_exists )then
            call self%eorecvol%kill
            self%eo_rec_tbox_exists = .false.
        endif
    end subroutine kill_eo_rec_tbox

    !> \brief  constructs the prime2D toolbox
    subroutine build_hadamard_prime2D_tbox( self, p )
        use simple_strings, only: str_has_substr
        class(build),  intent(inout) :: self
        class(params), intent(inout) :: p
        type(oris) :: os
        integer    :: icls, alloc_stat
        call self%kill_hadamard_prime2D_tbox
        call self%raise_hard_ctf_exception(p)
        allocate( self%cavgs(p%ncls), self%ctfsqsums(p%ncls), stat=alloc_stat )
        call alloc_err('build_hadamard_prime2D_tbox; simple_build, 1', alloc_stat)
        do icls=1,p%ncls
            call self%cavgs(icls)%new([p%box,p%box,1],p%smpd,p%imgkind)
            call self%ctfsqsums(icls)%new([p%box,p%box,1],p%smpd,p%imgkind)
        end do
        if( str_has_substr(p%refine,'neigh') )then
            if( file_exists(p%oritab3D) )then
                call os%new(p%ncls)
                call os%read(p%oritab3D)
                call os%nearest_neighbors(p%nnn, self%nnmat)
                call os%kill
            else
                stop 'need oritab3D input for prime2D refine=neigh mode; simple_build :: build_hadamard_prime2D_tbox'
            endif
        endif
        write(*,'(A)') '>>> DONE BUILDING HADAMARD PRIME2D TOOLBOX'
        self%hadamard_prime2D_tbox_exists = .true.
    end subroutine build_hadamard_prime2D_tbox

    !> \brief  destructs the prime2D toolbox
    subroutine kill_hadamard_prime2D_tbox( self )
        class(build), intent(inout) :: self
        integer :: i
        if( self%hadamard_prime2D_tbox_exists )then
            do i=1,size(self%cavgs)
                call self%cavgs(i)%kill
                call self%ctfsqsums(i)%kill
            end do
            deallocate(self%cavgs, self%ctfsqsums)
            self%hadamard_prime2D_tbox_exists = .false.
        endif
    end subroutine kill_hadamard_prime2D_tbox

    !> \brief  constructs the prime3D toolbox
    subroutine build_hadamard_prime3D_tbox( self, p )
        use simple_strings, only: str_has_substr
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        integer :: s, alloc_stat, nnn
        call self%kill_hadamard_prime3D_tbox
        call self%raise_hard_ctf_exception(p)
        ! reconstruction objects
        if( p%eo .eq. 'yes' )then
            allocate( self%eorecvols(p%nstates), stat=alloc_stat )
            call alloc_err('build_hadamard_prime3D_tbox; simple_build, 1', alloc_stat)
        else
            allocate( self%recvols(p%nstates), stat=alloc_stat )
            call alloc_err('build_hadamard_prime3D_tbox; simple_build, 2', alloc_stat)
        endif
        if( str_has_substr(p%refine,'neigh') )then
            nnn = p%nnn
            call self%se%nearest_neighbors(self%e, nnn, self%nnmat)
        endif
        write(*,'(A)') '>>> DONE BUILDING HADAMARD PRIME3D TOOLBOX'
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
            self%hadamard_prime3D_tbox_exists = .false.
        endif
    end subroutine kill_hadamard_prime3D_tbox

    !> \brief  constructs the toolbox for continuous refinement
    subroutine build_cont3D_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        integer :: s, alloc_stat
        call self%kill_cont3D_tbox
        call self%raise_hard_ctf_exception(p)
        if( p%norec .eq. 'yes' )then
            ! no reconstruction objects needed
        else
            if( p%eo .eq. 'yes' )then
                allocate( self%eorecvols(p%nstates), stat=alloc_stat )
                call alloc_err('build_hadamard_prime3D_tbox; simple_build, 1', alloc_stat)
            else
                allocate( self%recvols(p%nstates), stat=alloc_stat )
                call alloc_err('build_hadamard_prime3D_tbox; simple_build, 2', alloc_stat)
            endif
        endif
        allocate( self%refvols(p%nstates), stat=alloc_stat)
        call alloc_err('build_cont3D_tbox; simple_build, 3', alloc_stat)
        do s=1,p%nstates
            call self%refvols(s)%new([p%boxmatch,p%boxmatch,p%boxmatch],p%smpd,p%imgkind)
        end do
        write(*,'(A)') '>>> DONE BUILDING CONT3D TOOLBOX'
        self%cont3D_tbox_exists = .true.
    end subroutine build_cont3D_tbox

    !> \brief  destructs the toolbox for continuous refinement
    subroutine kill_cont3D_tbox( self )
        class(build), intent(inout) :: self
        integer :: i
        if( self%cont3D_tbox_exists )then
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
            if( allocated(self%refvols) )then
                do i=1,size(self%refvols)
                    call self%refvols(i)%kill_expanded
                    call self%refvols(i)%kill
                end do
                deallocate(self%refvols)
            endif
            self%cont3D_tbox_exists = .false.
        endif
    end subroutine kill_cont3D_tbox

    !> \brief  constructs the extremal3D toolbox
    subroutine build_extremal3D_tbox( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        integer :: alloc_stat
        call self%kill_extremal3D_tbox
        call self%raise_hard_ctf_exception(p)
        allocate( self%recvols(1), stat=alloc_stat )
        call alloc_err('build_hadamard_prime3D_tbox; simple_build, 2', alloc_stat)
        call self%recvols(1)%new([p%boxpd,p%boxpd,p%boxpd],p%smpd,p%imgkind)
        call self%recvols(1)%alloc_rho(p)
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

    !>  \brief  for reading feature vectors from disk
    subroutine read_features( self, p )
        class(build),  intent(inout) :: self
        class(params), intent(in)    :: p
        integer :: k, file_stat, funit, recsz
        inquire( iolength=recsz ) self%features(1,:)
        funit = get_fileunit()
        open(unit=funit, status='old', action='read', file=p%featstk,&
        access='direct', form='unformatted', recl=recsz, iostat=file_stat)
        call fopen_err('build_read_features', file_stat)
        ! read features from disk
        do k=1,p%nptcls
            read(funit, rec=k) self%features(k,:)
        end do
        close(unit=funit)
    end subroutine read_features

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
        myp = params(mycline_varying, checkdistr=.false.)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 1 passed'
        ! case 2:  refine=no, pgrp=c1, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=no, pgrp=c1, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying, checkdistr=.false.)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 2 passed'
        ! case 3:  refine=no, pgrp=c2, eo=yes
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=no, pgrp=c2, eo=yes'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying, checkdistr=.false.)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 3 passed'
        ! case 4:  refine=no, pgrp=c2, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=no, pgrp=c2, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying, checkdistr=.false.)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 4 passed'
        ! case 5:  refine=neigh, pgrp=c1, eo=yes
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c1, eo=yes'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying, checkdistr=.false.)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 5 passed'
        ! case 6:  refine=neigh, pgrp=c1, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c1, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying, checkdistr=.false.)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 6 passed'
        ! case 7:  refine=neigh, pgrp=c2, eo=yes
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c2, eo=yes'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying, checkdistr=.false.)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 7 passed'
        ! case 8:  refine=neigh, pgrp=c2, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c2, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying, checkdistr=.false.)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 8 passed'
        ! case 9:  refine=neigh, pgrp=c1, eo=yes
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c1, eo=yes'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying, checkdistr=.false.)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 9 passed'
        ! case 10: refine=neigh, pgrp=c1, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c1, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying, checkdistr=.false.)
        write(*,'(a)') '**info(simple_build_unit_test): case 10 passed'
        ! case 11: refine=neigh, pgrp=c2, eo=yes
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c2, eo=yes'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying, checkdistr=.false.)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 11 passed'
        ! case 12: refine=neigh, pgrp=c2, eo=no
        write(*,'(a)') '**info(simple_build_unit_test): testing case: refine=neigh, pgrp=c2, eo=no'
        mycline_varying = mycline_static
        call mycline_varying%set('refine', 'no')
        call mycline_varying%set('pgrp',   'c1')
        call mycline_varying%set('eo',    'yes')
        myp = params(mycline_varying, checkdistr=.false.)
        call tester
        write(*,'(a)') '**info(simple_build_unit_test): case 12 passed'
        write(*,'(a)') 'SIMPLE_BUILD_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'

      contains

          subroutine tester
              call myb%build_general_tbox(myp, mycline_varying, do3d=.true., nooritab=.true.)
              call myb%build_general_tbox(myp, mycline_varying, do3d=.true., nooritab=.true.)
              call myb%kill_general_tbox
              write(*,'(a)') 'build_general_tbox passed'
              call myb%build_cluster_tbox(myp)
              call myb%build_cluster_tbox(myp)
              call myb%kill_cluster_tbox
              write(*,'(a)') 'build_cluster_tbox passed'
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
              call myb%build_cont3D_tbox(myp)
              call myb%build_cont3D_tbox(myp)
              call myb%kill_cont3D_tbox
              write(*,'(a)') 'build_cont3D_tbox passed'
              call myb%build_extremal3D_tbox(myp)
              call myb%build_extremal3D_tbox(myp)
              call myb%kill_extremal3D_tbox
              write(*,'(a)') 'build_extremal3D_tbox passed'
          end subroutine tester

    end subroutine test_build

end module simple_build
