module simple_nanoparticle
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,                  only: image
use simple_binimage,               only: binimage
use simple_atoms,                  only: atoms, get_Z_and_radius_from_name
use simple_lattice_fitting,        only: fit_lattice, find_radius_for_coord_number
use simple_np_coordination_number, only: run_coord_number_analysis
implicit none

public :: nanoparticle
private
#include "simple_local_flags.inc"

! module global constants
integer, parameter :: N_THRESH      = 20      ! number of thresholds for binarization
logical, parameter :: DEBUG         = .false. ! for debugging purposes
logical, parameter :: GENERATE_FIGS = .false. ! for figures generation
integer, parameter :: SOFT_EDGE     = 6
integer, parameter :: N_DISCRET     = 1000

! container for per-atom statistics
type :: atom_stats
    integer :: cc_ind          = 0  ! index of the connected component        **filled-in**
    integer :: size            = 0  ! number of voxels in connected component **filled-in**
    integer :: cn_std          = 0  ! standard coordination number            **filled-in**
    integer :: loc_ldist(3)         ! for indentific of the vxl that determins the longest dim of the atom
    real    :: bondl           = 0. ! (dists)
    real    :: cn_gen          = 0. ! generalized coordination number         **filled-in**
    real    :: center(3)            ! atom center (centers)                   **filled-in**
    real    :: aspect_ratio    = 0. ! (ratios)                                **filled-in**
    real    :: polar_angle     = 0. ! polarization angle (ang_var)
    real    :: diam            = 0. ! atom diameter
    real    :: avg_int         = 0. ! average grey level intensity across the connected component **filled-in**
    real    :: max_int         = 0. ! maximum            -"-                                      **filled-in**
    real    :: sdev_int        = 0. ! standard deviation -"-                                      **filled-in**
    real    :: dist_from_NPcen = 0. ! distance from the centre of mass of the nanoparticle        **filled-in**
    real    :: strain          = 0. ! tensile strain in %
end type atom_stats

type :: nanoparticle
    private
    type(atoms)    :: centers_pdb
    type(image)    :: img,img_raw
    type(binimage) :: img_bin, img_cc
    integer        :: ldim(3)            = 0  ! logical dimension of image             **filled-in**
    integer        :: n_cc               = 0  ! number of atoms (connected components) **filled-in**
    real           :: smpd               = 0. ! sampling distance                      **filled-in**
    real           :: NPcen(3)           = 0. ! coordinates of the center of mass of the nanoparticle **filled-in**
    real           :: NPdiam             = 0. ! diameter of the nanoparticle           **filled-in**
    real           :: avg_bondl          = 0. ! average bond length in A
    real           :: max_bondl          = 0. ! maximum bond length in A
    real           :: min_bondl          = 0. ! minimum bond length in A
    real           :: sdev_bondl         = 0. ! standard deviation of  bond length in A
    real           :: med_bondl          = 0. ! median bond length in A
    real           :: theoretical_radius = 0. ! theoretical atom radius               **filled-in**
    real           :: net_dipole(3)      = 0. ! sum of all the directions of polarization
    type(atom_stats), allocatable :: atominfo(:)
    character(len=2)      :: element     = ' '
    character(len=4)      :: atom_name   = '    '
    character(len=STDLEN) :: partname    = '' ! fname
    character(len=STDLEN) :: fbody       = '' ! fbody
  contains
    ! constructor
    procedure          :: new => new_nanoparticle
    ! utils
    procedure, private :: atominfo2centers
    procedure, private :: atominfo2centers_A
    ! getters/setters
    procedure          :: get_ldim
    procedure          :: get_natoms
    procedure          :: set_img
    procedure          :: set_partname
    procedure          :: set_atomic_coords
    ! segmentation and statistics
    procedure, private :: binarize
    procedure, private :: find_centers
    procedure, private :: masscen
    procedure, private :: est_aspect_ratios
    procedure, private :: calc_aspect_ratio
    procedure, private :: discard_outliers
    procedure, private :: discard_outliers_cn
    procedure, private :: validate_atomic_positions
    procedure, private :: distances_distribution
    ! phase correlation
    procedure, private :: phasecorr_one_atom
    ! clustering
    procedure, private :: cluster_atom_intensity
    procedure          :: cluster_atom_maxint
    procedure          :: cluster_atom_intint
    procedure, private :: search_polarization
    procedure          :: cluster_ang
    procedure          :: cluster_ar
    procedure          :: cluster_bondl
    procedure          :: geometry_analysis
    ! execution
    procedure          :: identify_atomic_pos
    procedure          :: identify_atomic_pos_thresh
    procedure          :: fillin_atominfo
    procedure          :: radial_dependent_stats
    procedure          :: atoms_stats
    procedure          :: make_soft_mask
    ! visualization and output
    procedure          :: write_centers
    procedure          :: print_atoms
    ! others
    procedure          :: update_self_ncc
    procedure, private :: keep_atomic_pos_at_radius
    procedure, private :: center_on_atom
    procedure          :: mask
    ! kill
    procedure          :: kill => kill_nanoparticle
end type nanoparticle

contains

    subroutine new_nanoparticle(self, fname, cline_smpd, element)
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        real,                intent(in)    :: cline_smpd
        character(len=2),    intent(inout) :: element
        integer :: nptcls
        integer :: Z ! atomic number
        real    :: smpd
        call self%kill
        call self%set_partname(fname)
        self%fbody     = get_fbody(trim(basename(fname)), trim(fname2ext(fname)))
        self%smpd      = cline_smpd
        self%atom_name = ' '//element
        self%element   = element
        call get_Z_and_radius_from_name(self%atom_name, Z, self%theoretical_radius)
        call find_ldim_nptcls(self%partname, self%ldim, nptcls, smpd)
        call self%img%new(self%ldim, self%smpd)
        call self%img_raw%new(self%ldim, self%smpd)
        call self%img_bin%new_bimg(self%ldim, self%smpd)
        call self%img_bin%new(self%ldim, self%smpd)
        call self%img%read(fname)
        call self%img_raw%read(fname)
    end subroutine new_nanoparticle

    function atominfo2centers( self ) result( centers )
        class(nanoparticle), intent(in) :: self
        real, allocatable :: centers(:,:)
        integer :: sz, i
        sz = size(self%atominfo)
        allocate(centers(3,sz), source=0.)
        do i = 1, sz
            centers(:,i) = self%atominfo(i)%center(:)
        end do
    end function atominfo2centers

    function atominfo2centers_A( self ) result( centers_A )
        class(nanoparticle), intent(in) :: self
        real, allocatable :: centers_A(:,:)
        integer :: sz, i
        sz = size(self%atominfo)
        allocate(centers_A(3,sz), source=0.)
        do i = 1, sz
            centers_A(:,i) = (self%atominfo(i)%center(:) - 1.) * self%smpd
        end do
    end function atominfo2centers_A

    subroutine get_ldim(self,ldim)
        class(nanoparticle), intent(in)  :: self
        integer,             intent(out) :: ldim(3)
        ldim = self%img%get_ldim()
    end subroutine get_ldim

    function get_natoms(self) result(n)
       class(nanoparticle), intent(inout)  :: self
       integer :: n
       call self%img_cc%get_nccs(n)
    end function get_natoms

    ! set one of the images of the nanoparticle type
    subroutine set_img( self, imgfile, which )
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: imgfile
        character(len=*),    intent(in)    :: which
        select case(which)
            case('img')
                call self%img%new(self%ldim, self%smpd)
                call self%img%read(imgfile)
            case('img_bin')
                call self%img_bin%new_bimg(self%ldim, self%smpd)
                call self%img_bin%read_bimg(imgfile)
            case('img_cc')
                call self%img_cc%new_bimg(self%ldim, self%smpd)
                call self%img_cc%read_bimg(imgfile)
            case DEFAULT
                THROW_HARD('Wrong input parameter img type; set_img')
        end select
    end subroutine set_img

    ! set particle name (name of the input volume)
    subroutine set_partname( self, name )
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: name
        self%partname = name
    end subroutine set_partname

    subroutine update_self_ncc(self, img_cc)
        class(nanoparticle),      intent(inout) :: self
        type(binimage), optional, intent(inout) :: img_cc
        if( present(img_cc) )then
            call img_cc%get_nccs(self%n_cc)
        else
            call self%img_cc%get_nccs(self%n_cc)
        endif
    end subroutine update_self_ncc

    subroutine print_atoms( self )
        class(nanoparticle), intent(inout) :: self
        type(binimage)       :: img_atom
        integer, allocatable :: imat(:,:,:), imat_atom(:,:,:)
        integer :: i
        call img_atom%copy_bimg(self%img_cc)
        allocate(imat_atom(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
        call img_atom%get_imat(imat)
        do i = 1, maxval(imat)
            where(imat == i)
                imat_atom = 1
            elsewhere
                imat_atom = 0
            endwhere
            call img_atom%set_imat(imat_atom)
            call img_atom%write_bimg('Atom'//trim(int2str(i))//'.mrc')
        enddo
        deallocate(imat,imat_atom)
        call img_atom%kill_bimg
    end subroutine print_atoms

    ! Find the centers coordinates of the atoms in the particle
    ! and save it in the global variable centers.
    ! If coords is present, it saves it also in coords.
    subroutine find_centers(self, img_bin, img_cc, coords)
       class(nanoparticle),         intent(inout) :: self
       type(binimage), optional,    intent(inout) :: img_bin, img_cc
       real, optional, allocatable, intent(out)   :: coords(:,:)
       real,        pointer :: rmat_raw(:,:,:)
       integer, allocatable :: imat_cc_in(:,:,:)
       logical, allocatable :: mask(:,:,:)
       integer :: i, ii, jj, kk
       real    :: m(3), sum_mass
       ! sanity check
       if(present(img_bin) .and. .not. present(img_cc)) THROW_HARD('img_bin and img_cc have to be both present in input')
       ! global variables allocation
       if( allocated(self%atominfo) ) deallocate(self%atominfo)
       allocate(self%atominfo(self%n_cc))
       if(present(img_cc)) then
           call img_cc%get_imat(imat_cc_in)
       else
           call self%img_cc%get_imat(imat_cc_in)
       endif
       call self%img_raw%get_rmat_ptr(rmat_raw)
       allocate(mask(self%ldim(1),self%ldim(2),self%ldim(3)), source = .true.)
       !$omp parallel do default(shared) private(i,ii,jj,kk,mask,m,sum_mass) schedule(static) proc_bind(close)
       do i=1,self%n_cc
           mask = .true.
           where( imat_cc_in /= i ) mask = .false.
           m        = 0.
           sum_mass = 0.
           do ii = 1, self%ldim(1)
               do jj = 1, self%ldim(2)
                   do kk = 1, self%ldim(3)
                       if(mask(ii,jj,kk)) then
                           m = m + real([ii,jj,kk]) * rmat_raw(ii,jj,kk)
                           sum_mass = sum_mass + rmat_raw(ii,jj,kk)
                       endif
                   enddo
               enddo
           enddo
           self%atominfo(i)%center(:) = m/sum_mass
       enddo
       !$omp end parallel do
       ! saving centers coordinates, optional
       if( present(coords) )then
           allocate(coords(3,self%n_cc))
           do i=1,self%n_cc
               coords(:,i) = self%atominfo(i)%center(:)
           end do
       endif
    end subroutine find_centers

    subroutine write_centers(self, fname, coords)
        class(nanoparticle),        intent(inout) :: self
        character(len=*), optional, intent(in)    :: fname
        real,             optional, intent(in)    :: coords(:,:)
        integer :: cc
        if(present(coords)) then
            call self%centers_pdb%new(size(coords, dim = 2), dummy=.true.)
            do cc=1,size(coords, dim = 2)
                call self%centers_pdb%set_name(cc,self%atom_name)
                call self%centers_pdb%set_element(cc,self%element)
                call self%centers_pdb%set_coord(cc,(coords(:,cc)-1.)*self%smpd)
            enddo
        else
            call self%centers_pdb%new(self%n_cc, dummy=.true.)
            do cc=1,self%n_cc
                call self%centers_pdb%set_name(cc,self%atom_name)
                call self%centers_pdb%set_element(cc,self%element)
                call self%centers_pdb%set_coord(cc,(self%atominfo(cc)%center(:)-1.)*self%smpd)
                call self%centers_pdb%set_beta(cc,self%atominfo(cc)%cn_gen) ! use generalised coordination number
                call self%centers_pdb%set_resnum(cc,cc)
            enddo
        endif
        if(present(fname)) then
            call self%centers_pdb%writepdb(fname)
        else
            call self%centers_pdb%writepdb(trim(self%fbody)//'_atom_centers')
        endif
    end subroutine write_centers

    ! This subroutine sets the atom positions to be the ones in the inputted PDB file.
    subroutine set_atomic_coords(self, pdb_file)
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: pdb_file
        type(atoms) :: a
        integer     :: i, N
        if(fname2ext(pdb_file) .ne. 'pdb') THROW_HARD('Inputted filename has to be PDB file!; set_atomic_coords')
        if( allocated(self%atominfo) ) deallocate(self%atominfo)
        call a%new(pdb_file)
        N = a%get_n() ! number of atoms
        allocate(self%atominfo(N))
        do i = 1, N
            self%atominfo(i)%center(:) = a%get_coord(i)/self%smpd + 1.
        enddo
        self%n_cc = N
        call a%kill
    end subroutine set_atomic_coords

    ! calc the avg of the centers coords
    function masscen( self ) result( m )
        class(nanoparticle), intent(inout) :: self
        real    :: m(3) ! mass center coords
        integer :: i
        m = [0.,0.,0.]
        do i = 1, self%n_cc
            m = m + self%atominfo(i)%center(:)
        end do
        m = m / real(self%n_cc)
    end function masscen

    ! FORMULA: phasecorr = ifft2(fft2(field).*conj(fft2(reference)));
    subroutine phasecorr_one_atom(self, out_img)
        use simple_segmentation
        class(nanoparticle), intent(inout) :: self
        type(image),         intent(inout) :: out_img  ! where to save the correlation values
        type(image) :: one_atom, img_copy, phasecorr
        type(atoms) :: atom
        real :: cutoff
        call img_copy%new(self%ldim, self%smpd)
        call img_copy%copy(self%img)
        call phasecorr%new(self%ldim, self%smpd)
        call phasecorr%set_ft(.true.)
        call one_atom%new(self%ldim,self%smpd)
        cutoff = 8.*self%smpd
        call atom%new(1)
        call atom%set_element(1,self%element)
        call atom%set_coord(1,self%smpd*(real(self%ldim)/2.)) !DO NOT NEED THE +1
        call atom%convolve(one_atom, cutoff)
        call one_atom%fft()
        call img_copy%fft()
        call img_copy%phase_corr(one_atom,phasecorr,1.)
        if( GENERATE_FIGS ) call phasecorr%write(trim(self%fbody)//'CorrFiltered.mrc')
        ! Save Output
        call out_img%kill() ! if it exists already
        call out_img%new(self%ldim, self%smpd)
        call out_img%fft()
        call out_img%copy(phasecorr)
        call one_atom%kill()
        call phasecorr%kill()
    end subroutine phasecorr_one_atom

    ! This subrotuine takes in input a nanoparticle and
    ! binarizes it by thresholding. The gray level histogram is split
    ! in 20 parts, which corrispond to 20 possible threshold.
    ! Among those threshold, the selected one is the for which
    ! tha correlation between the raw map and a simulated distribution
    ! obtained with that threshold reaches the maximum value.
    subroutine binarize( self )
        class(nanoparticle), intent(inout) :: self
        type(binimage)       :: img_bin_thresh(N_THRESH/2-1)
        type(binimage)       :: img_ccs_thresh(N_THRESH/2-1)
        type(image)          :: pc
        type(atoms)          :: atom
        type(image)          :: simulated_distrib
        integer, allocatable :: imat_t(:,:,:)
        real,    allocatable :: x_mat(:)  !vectorization of the volume
        real,    allocatable :: coords(:,:)
        real,    allocatable :: rmat(:,:,:)
        integer :: i, t
        real    :: otsu_thresh
        real    :: x_thresh(N_THRESH/2-1)
        real    :: step, maximum, mm(2)
        call otsu_nano(self%img,otsu_thresh) ! find initial threshold
        write(logfhandle,*) '****binarization, init'
        rmat = self%img%get_rmat()
        x_mat = pack(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) >= 0.)
        allocate(imat_t(self%ldim(1), self%ldim(2), self%ldim(3)), source = 0)
        step = (maxval(x_mat)-otsu_thresh )/real(N_THRESH)
        deallocate(x_mat)
        call simulated_distrib%new(self%ldim,self%smpd)
        call pc%new(self%ldim,self%smpd)
        call self%img%fft ! for pc calculation
        t = 1
        maximum = 0.
        do i = 1, N_THRESH/2-1
            call progress(i, N_THRESH/2-1)
            call img_bin_thresh(i)%new_bimg(self%ldim, self%smpd)
            if(i == 1) then
                x_thresh(i) = otsu_thresh
            else
                x_thresh(i) = x_thresh(i-1) + step
            endif
            where(rmat > x_thresh(i))
                imat_t = 1
            elsewhere
                imat_t = 0
            endwhere
            ! Generate binary image and cc image
            call img_bin_thresh(i)%set_imat(imat_t)
            call img_bin_thresh(i)%find_ccs(img_ccs_thresh(i))
            ! Find atom centers in the generated distributions
            call self%update_self_ncc(img_ccs_thresh(i)) ! self%n_cc is needed in find_centers
            call self%find_centers(img_bin_thresh(i), img_ccs_thresh(i), coords)
            ! Generate a simulated distribution based on those center
            call self%write_centers('centers_'//trim(int2str(i))//'_iteration', coords)
            call atom%new          ('centers_'//trim(int2str(i))//'_iteration.pdb')
            call atom%convolve(simulated_distrib, cutoff = 8.*self%smpd)
            call del_file('centers_'//trim(int2str(i))//'_iteration.pdb')
            if(DEBUG) call simulated_distrib%write('simulated_'//trim(int2str(i))//'_iteration.mrc')
            call atom%kill
            ! Take care of Fourier status, for phase_corr calculation
            call simulated_distrib%fft
            call pc%zero_and_flag_ft
            ! Correlation volume generation
            call self%img%phase_corr(simulated_distrib, pc, lp=3.)
            ! Calculation and update of the maxval the correlation reaches
            mm = pc%minmax()
            if(mm(2) > maximum) then
                maximum = mm(2)
                t       = i
            endif
            call simulated_distrib%set_ft(.false.)
        enddo
        call self%img%ifft ! To remove
        write(logfhandle,*) 'Selected threshold: ', x_thresh(t)
        ! Update img_bin and img_cc
        call self%img_bin%copy_bimg(img_bin_thresh(t))
        if(GENERATE_FIGS) call self%img_bin%write_bimg(trim(self%fbody)//'SelectedThreshold.mrc')
        call self%img_cc%copy_bimg(img_ccs_thresh(t))
        call self%update_self_ncc()
        call self%find_centers()
        do i = 1,  N_THRESH/2-1
            call img_bin_thresh(i)%kill_bimg
            call img_ccs_thresh(i)%kill_bimg
        enddo
        ! deallocate and kill
        if(allocated(rmat))   deallocate(rmat)
        if(allocated(imat_t)) deallocate(imat_t)
        if(allocated(coords)) deallocate(coords)
        call simulated_distrib%kill
        write(logfhandle,*) '****binarization, completed'
    contains
        !Otsu binarization for nanoparticle maps
        !It considers the gray level value just in the positive range.
        !It doesn't threshold the map. It just returns the ideal threshold.
        ! This is based on the implementation of 1D otsu
         subroutine otsu_nano(img, scaled_thresh)
             use simple_math, only : otsu
             type(image),    intent(inout) :: img
             real,           intent(out)   :: scaled_thresh !returns the threshold in the correct range
             real, pointer     :: rmat(:,:,:)
             real, allocatable :: x(:)
             call img%get_rmat_ptr(rmat)
             x = pack(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) > 0.)
             call otsu(x, scaled_thresh)
         end subroutine otsu_nano
    end subroutine binarize

    ! This subroutine discards outliers that resisted binarization.
    ! If generalised is true
    ! It calculates the generalised coordination number (cn_gen) of each atom and discards
    ! the atoms with cn_gen < cn_thresh
    ! If generalised is false
    ! It calculates the standard coordination number (cn) of each atom and discards
    ! the atoms with cn_std < cn_thresh
    ! It modifies the img_bin and img_cc instances deleting the
    ! identified outliers.
    subroutine discard_outliers_cn(self, cn_thresh, generalised)
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: cn_thresh   ! threshold for discard outliers based on coordination number
        logical,             intent(in)    :: generalised ! use cn_gen or cn standard?
        integer, allocatable :: imat_bin(:,:,:), imat_cc(:,:,:)
        logical, allocatable :: mask(:)
        real, allocatable    :: centers_A(:,:) ! coordinates of the atoms in ANGSTROMS
        real    :: radius  ! radius of the sphere to consider for cn calculation
        real    :: a(3)    ! lattice parameter
        real    :: cn_gen(self%n_cc)
        integer :: cn(self%n_cc), cc, n_discard, filnum, filnum1, filnum2, i
        write(logfhandle, *) '****outliers discarding cn, init'
        centers_A = self%atominfo2centers_A()
        ! In order to find radius for cn calculation, FCC Lattice has to be fit.
        if(DEBUG) then
            call fopen(filnum, file='CentersAngCN.txt')
            do i = 1, size(centers_A,2)
                write(filnum,*) centers_A(:,i)
            enddo
            call fclose(filnum)
            call fopen(filnum1, file='CentersPxlsCN.txt')
            do i = 1, size(self%atominfo)
                write(filnum1,*) self%atominfo(i)%center(:)
            enddo
            call fclose(filnum1)
        endif
        call fit_lattice(centers_A,a)
        call find_radius_for_coord_number(a,radius)
        if( DEBUG ) write(logfhandle,*) 'Radius for coord nr calculation ', radius
        call run_coord_number_analysis(centers_A,radius,cn,cn_gen)
        if( DEBUG )then
            call fopen(filnum2, file='CoordNumberGenBeforeOutliers.txt')
            write(filnum2,*) cn_gen
            call fclose(filnum2)
        endif
        allocate(mask(self%n_cc), source = .true.)
        if( generalised .eqv. .true. )then
            if( DEBUG )then
                write(logfhandle, *) 'Before outliers discarding cn_gen is'
                write(logfhandle, *)  cn_gen
            endif
            write(logfhandle, *) 'Using GENERALISED cn'
            where(cn_gen < cn_thresh) mask = .false. ! false where atom has to be discarded
            n_discard = count(cn_gen < cn_thresh)
        else
            if(DEBUG) then
                write(logfhandle, *) 'Before outliers discarding cn standard is'
                write(logfhandle, *)  cn
            endif
            write(logfhandle, *) 'Using STANDARD cn'
            where(cn < cn_thresh) mask = .false. ! false where atom has to be discarded
            n_discard = count(cn < cn_thresh)
        endif
        write(logfhandle, *) 'Numbers of atoms discarded because of low cn ', n_discard
        call self%img_cc%get_imat(imat_cc)
        call self%img_bin%get_imat(imat_bin)
        ! Removing outliers from the binary image and the connected components image
        if( generalised .eqv. .true. )then
            do cc = 1, self%n_cc
                if(cn_gen(cc)<cn_thresh) then
                    where(imat_cc == cc) imat_bin = 0
                endif
            enddo
        else
            do cc = 1, self%n_cc
                if(cn(cc)<cn_thresh) then
                    where(imat_cc == cc) imat_bin = 0
                endif
            enddo
        endif
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%find_ccs(self%img_cc)
        ! update number of connected components
        call self%img_cc%get_nccs(self%n_cc)
        call self%find_centers()
        deallocate(centers_A)
        centers_A = self%atominfo2centers_A()
        call run_coord_number_analysis(centers_A,radius,self%atominfo(:)%cn_std,self%atominfo(:)%cn_gen)
        ! ATTENTION: you will see low coord numbers because they are UPDATED, after elimination
        ! of the atoms with low cn. It is like this in order to be consistent with the figure.
        if(DEBUG) then
           write(logfhandle, *) 'After outliers discarding cn is'
           write(logfhandle, *)  self%atominfo(:)%cn_std
           write(logfhandle, *) 'And generalised cn is'
           write(logfhandle, *)  self%atominfo(:)%cn_gen
        endif
        write(logfhandle, *) '****outliers discarding cn, completed'
    end subroutine discard_outliers_cn

    ! This subroutine discards outliers that resisted binarization.
    ! It calculates the generalised coordination number of each atom and discards the bottom
    ! PERCENT_DISCARD% of the atoms according to the contact score.
    ! It modifies the img_bin and img_cc instances deleting the
    ! identified outliers.
    subroutine discard_outliers( self )
        class(nanoparticle), intent(inout) :: self
        integer, allocatable :: imat_bin(:,:,:), imat_cc(:,:,:)
        real,    parameter   :: PERCENT_DISCARD = 5./100.
        logical, allocatable :: mask(:)
        real,    allocatable :: centers_A(:,:) ! coordinates of the atoms in ANGSTROMS
        real    :: cn_gen(self%n_cc), a(3), dist
        integer :: cn(self%n_cc)
        real    :: radius ! radius of the sphere to consider
        integer :: cc, loc(1), n_discard, filnum1, filnum, i, label(1)
        write(logfhandle, *) '****outliers discarding default, init'
        centers_A = self%atominfo2centers_A() ! atom centers in A
        ! In order to find radius for cn calculation, FCC Lattice has to be fit.
        if( DEBUG )then
            call fopen(filnum, file='CentersAngCN.txt')
            do i = 1, size(centers_A,2)
                write(filnum,*) centers_A(:,i)
            enddo
            call fclose(filnum)
            call fopen(filnum1, file='CentersPxlsCN.txt')
            do i = 1, size(self%atominfo)
                write(filnum1,*) self%atominfo(i)%center(:)
            enddo
            call fclose(filnum1)
        endif
        call fit_lattice(centers_A,a)
        call find_radius_for_coord_number(a,radius)
        if( DEBUG ) write(logfhandle,*) 'Radius for coord nr calculation ', radius
        call run_coord_number_analysis(centers_A, radius, cn, cn_gen)
        allocate(mask(self%n_cc), source = .true.)
        n_discard = nint(self%n_cc*PERCENT_DISCARD)
        write(logfhandle, *) 'Numbers of atoms discarded because outliers ', n_discard
        call self%img_cc%get_imat(imat_cc)
        call self%img_bin%get_imat(imat_bin)
        ! Removing outliers from the binary image and the connected components image
        mask(1:self%n_cc) = .true. ! reset
        do cc = 1, n_discard
            label = minloc(cn_gen, mask)
            where(imat_cc == label(1))
               imat_bin = 0
            endwhere
            mask(label(1)) = .false.
        enddo
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%find_ccs(self%img_cc)
        ! update number of connected components
        call self%img_cc%get_nccs(self%n_cc)
        call self%find_centers()
        deallocate(centers_A)
        centers_A = self%atominfo2centers_A()
        call run_coord_number_analysis(centers_A, radius, self%atominfo(:)%cn_std, self%atominfo(:)%cn_gen)
        ! ATTENTION: you will see low coord numbers because they are UPDATED, after elimination
        ! of the atoms with low cn. It is like this in order to be consistent with the figure.
        if( DEBUG )then
           write(logfhandle, *) 'After outliers discarding cn is'
           write(logfhandle, *)  self%atominfo(:)%cn_std
           write(logfhandle, *) 'And generalised cn is'
           write(logfhandle, *)  self%atominfo(:)%cn_gen
        endif
        write(logfhandle, *) '****outliers discarding default, completed'
    end subroutine discard_outliers

    ! This subrouitne validates the identified atomic positions
    subroutine validate_atomic_positions( self )
        class(nanoparticle), intent(inout) :: self
        real,    allocatable :: x(:)
        real,    pointer     :: rmat_pc(:,:,:)
        integer, allocatable :: imat(:,:,:), imat_cc(:,:,:), imat_bin(:,:,:)
        integer, parameter   :: RANK_THRESH = 4
        integer :: icc, cnt
        integer :: rank, m(1)
        real    :: new_centers(3,2*self%n_cc)               ! will pack it afterwards if it has too many elements
        integer :: new_coordination_number(2*self%n_cc)     ! will pack it afterwards if it has too many elements
        real    :: new_coordination_number_gen(2*self%n_cc) ! will pack it afterwards if it has too many elements
        real    :: pc
        call self%img%get_rmat_ptr(rmat_pc) ! now img contains the phase correlation
        call self%img_cc%get_imat(imat_cc)  ! to pass to the subroutine split_atoms
        allocate(imat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), source = imat_cc)
        cnt = 0
        ! Remember to update the centers
        do icc =1, self%n_cc ! for each cc check if the center corresponds with the local max of the phase corr
            pc = rmat_pc(nint(self%atominfo(icc)%center(1)),nint(self%atominfo(icc)%center(2)),nint(self%atominfo(icc)%center(3)))
            ! calculate the rank
            x = pack(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), imat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) == icc)
            call hpsort(x)
            m(:) = minloc(abs(x - pc))
            rank = size(x) - m(1)
            deallocate(x)
            if( rank > RANK_THRESH )then
                call split_atom(new_centers,cnt)
            else
                cnt = cnt + 1 ! new number of centers derived from splitting
                new_centers(:,cnt)               = self%atominfo(icc)%center(:)
                new_coordination_number(cnt)     = self%atominfo(icc)%cn_std
                new_coordination_number_gen(cnt) = self%atominfo(icc)%cn_gen
            endif
        enddo
        deallocate(self%atominfo)
        self%n_cc = cnt ! update
        allocate(self%atominfo(cnt))
        ! update centers
        do icc = 1, cnt
            self%atominfo(icc)%center(:) = new_centers(:,icc)
        enddo
        ! update contact scores
        do icc = 1, cnt
            self%atominfo(icc)%cn_std = new_coordination_number(icc)
            self%atominfo(icc)%cn_gen = new_coordination_number_gen(icc)
        enddo
        call self%img_bin%get_imat(imat_bin)
        if( DEBUG ) call self%img_bin%write_bimg(trim(self%fbody)//'BINbeforeValidation.mrc')
        ! update binary image
        where( imat_cc > 0 )
            imat_bin = 1
        elsewhere
            imat_bin = 0
        endwhere
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%update_img_rmat()
        call self%img_bin%write_bimg(trim(self%fbody)//'BIN.mrc')
        call self%img_bin%find_ccs(self%img_cc)
        call self%img_cc%write_bimg(trim(self%fbody)//'CC.mrc')
        ! update number of ccs
        call self%update_self_ncc(self%img_cc)
        ! update and write centers
        call self%find_centers()
        call self%write_centers()

    contains

        subroutine split_atom(new_centers,cnt)
            real,    intent(inout) :: new_centers(:,:)  ! updated coordinates of the centers
            integer, intent(inout) :: cnt               ! atom counter, to update the center coords
            integer :: new_center1(3), new_center2(3), new_center3(3)
            integer :: i, j, k
            logical :: mask(self%ldim(1),self%ldim(2),self%ldim(3)) ! false in the layer of connection of the atom to be split
            mask = .false.  ! initialization
            ! Identify first new center
            new_center1 = maxloc(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), imat == icc)
            cnt = cnt + 1
            new_centers(:,cnt)               = real(new_center1)
            new_coordination_number(cnt)     = self%atominfo(icc)%cn_std
            new_coordination_number_gen(cnt) = self%atominfo(icc)%cn_gen
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    do k = 1, self%ldim(3)
                        if(((real(i-new_center1(1)))**2 + (real(j-new_center1(2)))**2 + &
                        &   (real(k-new_center1(3)))**2)*self%smpd  <=  (0.9*self%theoretical_radius)**2) then
                            if(imat(i,j,k) == icc) mask(i,j,k) = .true.
                        endif
                    enddo
                enddo
            enddo
            ! Second likely center.
            new_center2 = maxloc(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), (imat == icc) .and. .not. mask)
            if( any(new_center2 > 0) )then ! if anything was found
                ! Validate second center (check if it's 2 merged atoms, or one pointy one)
                if(sum(real(new_center2-new_center1)**2.)*self%smpd <= 2.*self%theoretical_radius) then
                    ! Set the merged cc back to 0
                    where(imat_cc == icc .and. (.not.mask) ) imat_cc = 0
                    return
                else
                    cnt                                = cnt + 1
                    new_centers(:,cnt)                 = real(new_center2)
                    new_coordination_number(cnt)       = self%atominfo(icc)%cn_std + 1
                    new_coordination_number(cnt-1)     = new_coordination_number(cnt)
                    new_coordination_number_gen(cnt)   = self%atominfo(icc)%cn_gen + 1.
                    new_coordination_number_gen(cnt-1) = new_coordination_number_gen(cnt)
                    ! In the case two merged atoms, build the second atom
                    do i = 1, self%ldim(1)
                        do j = 1, self%ldim(2)
                            do k = 1, self%ldim(3)
                                if(((real(i-new_center2(1)))**2 + (real(j-new_center2(2)))**2 + (real(k-new_center2(3)))**2)*self%smpd < &
                                & (0.9*self%theoretical_radius)**2) then
                                    if(imat(i,j,k) == icc)   mask(i,j,k) = .true.
                                endif
                            enddo
                        enddo
                    enddo
                endif
            endif
            ! Third likely center.
            new_center3 = maxloc(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), &
            & (imat == icc) .and. .not. mask)
            if(any(new_center3 > 0)) then ! if anything was found
                ! Validate third center
                if(sum(real(new_center3-new_center1)**2.)*self%smpd <= 2.*self%theoretical_radius .or. &
                &  sum(real(new_center3-new_center2)**2.)*self%smpd <= 2.*self%theoretical_radius ) then
                    ! Set the merged cc back to 0
                    where( imat_cc == icc .and. (.not.mask) ) imat_cc = 0
                    return
                else
                    cnt = cnt + 1
                    new_centers(:,cnt) = real(new_center3)
                    new_coordination_number(cnt)   = self%atominfo(icc)%cn_std + 2
                    new_coordination_number(cnt-1) = new_coordination_number(cnt)
                    new_coordination_number(cnt-2) = new_coordination_number(cnt)
                    new_coordination_number_gen(cnt)   = self%atominfo(icc)%cn_gen + 2.
                    new_coordination_number_gen(cnt-1) = new_coordination_number_gen(cnt)
                    new_coordination_number_gen(cnt-2) = new_coordination_number_gen(cnt)
                    ! In the case two merged atoms, build the second atom
                    do i = 1, self%ldim(1)
                        do j = 1, self%ldim(2)
                            do k = 1, self%ldim(3)
                                if(((real(i-new_center3(1)))**2 + (real(j-new_center3(2)))**2 + &
                                &(real(k-new_center3(3)))**2)*self%smpd < (0.9*self%theoretical_radius)**2) then ! a little smaller to be sure
                                    if(imat(i,j,k) == icc)   mask(i,j,k) = .true.
                                endif
                            enddo
                        enddo
                    enddo
                endif
            endif
            ! Set the merged cc back to 0
            where(imat_cc == icc .and. (.not.mask) ) imat_cc = 0
            call self%img_cc%set_imat(imat_cc)
        end subroutine split_atom

    end subroutine validate_atomic_positions

    subroutine fillin_atominfo( self )
        class(nanoparticle), intent(inout) :: self
        logical, allocatable :: mask(:,:,:)
        real,    allocatable :: centers_A(:,:), tmpcens(:,:)
        real,    pointer     :: rmat(:,:,:)
        integer, allocatable :: imat_cc(:,:,:)
        real    :: tmp_diam, a(3), radius_cn
        logical :: cc_mask(self%n_cc)
        integer :: i, j, k, cc
        ! calc cn and cn_gen
        centers_A = self%atominfo2centers_A()
        call fit_lattice(centers_A,a)
        call find_radius_for_coord_number(a,radius_cn)
        call run_coord_number_analysis(centers_A,radius_cn,self%atominfo(:)%cn_std,self%atominfo(:)%cn_gen)
        deallocate(centers_A)
        ! calc NPdiam & NPcen
        tmpcens     = self%atominfo2centers()
        cc_mask     = .true.
        self%NPdiam = 0.
        do i = 1, self%n_cc
            tmp_diam = pixels_dist(self%atominfo(i)%center(:), tmpcens, 'max', cc_mask)
            if( tmp_diam > self%NPdiam ) self%NPdiam = tmp_diam
        enddo
        deallocate(tmpcens)
        self%NPdiam = self%NPdiam * self%smpd ! in A
        self%NPcen  = self%masscen()
        ! extract atominfo
        allocate(mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), source = .false.)
        call self%img%get_rmat_ptr(rmat)
        call self%img_cc%get_imat(imat_cc)
        ! type :: atom_stats
        !     integer :: cc_ind          = 0  ! index of the connected component        **filled-in**
        !     integer :: size            = 0  ! number of voxels in connected component **filled-in**
        !     integer :: cn_std          = 0  ! standard coordination number            **filled-in**
        !     integer :: loc_ldist(3)         ! for indentific of the vxl that determins the longest dim of the atom
        !     real    :: bondl           = 0. ! (dists)
        !     real    :: cn_gen          = 0. ! generalized coordination number         **filled-in**
        !     real    :: center(3)            ! atom center (centers)                   **filled-in**
        !     real    :: aspect_ratio    = 0. ! (ratios)                                **filled-in**
        !     real    :: polar_angle     = 0. ! polarization angle (ang_var)
        !     real    :: diam            = 0. ! atom diameter
        !     real    :: avg_int         = 0. ! average grey level intensity across the connected component **filled-in**
        !     real    :: max_int         = 0. ! maximum            -"-                                      **filled-in**
        !     real    :: sdev_int        = 0. ! standard deviation -"-                                      **filled-in**
        !     real    :: dist_from_NPcen = 0. ! distance from the centre of mass of the nanoparticle        **filled-in**
        !     real    :: strain          = 0. ! tensile strain in %
        ! end type atom_stats
        do cc = 1, self%n_cc
            ! index of the connected component
            self%atominfo(cc)%cc_ind = cc
            ! number of voxels in connected component
            where( imat_cc == cc ) mask = .true.
            self%atominfo(cc)%size = count(mask)
            ! distance from the centre of mass of the nanoparticle
            self%atominfo(cc)%dist_from_NPcen = euclid(self%atominfo(cc)%center(:), self%NPcen) * self%smpd
            ! atom aspect ratio
            call self%calc_aspect_ratio(cc, self%atominfo(cc)%aspect_ratio, lld=.false., print_ar=.false.)
            ! maximum grey level intensity across the connected component
            self%atominfo(cc)%max_int = maxval(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), mask)
            ! average grey level intensity across the connected component
            self%atominfo(cc)%avg_int = sum(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), mask)
            self%atominfo(cc)%avg_int = self%atominfo(cc)%avg_int / real(count(mask))
            ! standard deviation of the grey level intensity across the connected component
            self%atominfo(cc)%sdev_int = 0.
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    do k = 1, self%ldim(3)
                        if( mask(i,j,k) ) self%atominfo(cc)%sdev_int = &
                        &self%atominfo(cc)%sdev_int + (rmat(i,j,k) - self%atominfo(cc)%avg_int)**2.
                    enddo
                enddo
            enddo
            if( count(mask) > 1 )then
                self%atominfo(cc)%sdev_int = sqrt(self%atominfo(cc)%sdev_int / real(count(mask)-1))
            else ! atom composed by one voxel
                self%atominfo(cc)%sdev_int = 0.
            endif
            mask = .false. ! reset mask
        end do
        deallocate(mask, imat_cc)
        call self%distances_distribution(file=.false.) ! across the whole nano
   end subroutine fillin_atominfo

    ! This subroutine calculates some basic stats in the nanoparticle.
    ! In particular it reports the diameter of the nanocrystal
    !(the longest distance between any two atoms) and the overall
    ! atomic density (nr of atoms per volume unit A^3) and the
    ! atomic density with radial dependency
    ! Moreover, this subroutine calculates the histogram of the within-atoms
    ! distances distribution within the nanoparticle nano.
    ! To each atom the distance assigned is the min distance
    ! to the other atoms. There is a threshold (3.*self%theoretical_radius) for
    ! outliers discarding.
    ! If coords in input, then it considers just the atom-to-atom
    ! distance between the atoms with centers in coords.
    ! At the same time, it calculated statistics in a radial-dependent
    ! way, and save the result on a file.
    subroutine radial_dependent_stats(self,min_rad, max_rad, step, cn_min, cn_max, cn)
        class(nanoparticle), intent(inout) :: self
        real,                intent(in)    :: min_rad, max_rad, step
        integer, optional,   intent(in)    :: cn_min, cn_max
        integer, optional,   intent(in)    :: cn
        type(atoms)          :: radial_atoms_just,  radial_atoms_all
        type(image)          :: simulated_density
        logical, allocatable :: mask(:,:,:)
        real,    allocatable :: coords(:,:) !coordinates of the centers of the atoms according to radial distances
        real,    allocatable :: max_intensity(:), avg_intensity(:), stdev_intensity(:), ratios(:), centers_A(:,:), tmpcens(:,:)
        real,    pointer     :: rmat(:,:,:)
        integer, allocatable :: imat_cc(:,:,:)
        real    :: nano_diameter, temp_diameter !for the diameter of the nanoparticle
        real    :: m(3)    !mass center of c3 map
        real    :: d       !distance atoms from the center
        real    :: radius  !radius of the sphere to consider
        real    :: avg_int, max_int, stdev_int, cutoff, a(3), radius_cn
        real    :: volume ! to calculate the atomic density
        integer :: cnt, cnt_just, cnt_all, filnum, io_stat, filnum2
        integer :: nsteps, i, j, k, l, cc
        logical :: cn_mask(self%n_cc)
        ! Calculate cn and cn_gen
        centers_A = self%atominfo2centers_A()
        call fit_lattice(centers_A,a)
        call find_radius_for_coord_number(a,radius_cn)
        call run_coord_number_analysis(centers_A,radius_cn,self%atominfo(:)%cn_std,self%atominfo(:)%cn_gen)
        deallocate(centers_A)
        ! Mask for elimination of atoms with std cn outside selected range
        if(present(cn_min)) then
            do cc = 1, self%n_cc
                if( self%atominfo(cc)%cn_std > cn_max .or. self%atominfo(cc)%cn_std < cn_min ) cn_mask(cc) = .false.
            enddo
        else
            cn_mask(:) = .true.
        endif
        nano_diameter = 0.
        tmpcens = self%atominfo2centers()
        do i = 1, self%n_cc
            temp_diameter = pixels_dist(self%atominfo(i)%center(:), tmpcens, 'max', cn_mask)
            if(temp_diameter > nano_diameter) nano_diameter = temp_diameter
        enddo
        deallocate(tmpcens)
        nano_diameter = nano_diameter*self%smpd ! in A
        volume = 4./3.*pi*(nano_diameter/2.*self%smpd)**3
        cnt_all  = 0
        cutoff = 8.*self%smpd
        ! min_step and max_step is in A
        self%n_cc = size(self%atominfo)
        if(present(cn_min)) then
          write(logfhandle, *) '****radial atom-to-atom distances estimation with cn std range, init'
        else
          write(logfhandle, *) '****radial atom-to-atom distances estimation, init'
        endif
        if(abs(max_rad-min_rad) > TINY) then
          nsteps = floor((max_rad-min_rad)/step)+1
        else
          nsteps = 1 ! minrad and maxrad coincide
        endif
        m = self%masscen()
        ! Report statistics in dedicated directory
        call fopen(filnum, file='RadialStat.txt', iostat=io_stat)
        call fopen(filnum2, file='SimpleStat.txt', iostat=io_stat)
        write(unit = filnum2, fmt = '(a,f6.3)')  'Nanoparticle Diameter: ', nano_diameter
        write(unit = filnum2, fmt = '(a,f6.3)')  'Overall density:       ', real(self%n_cc)/volume
        allocate(mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), source = .false.)
        call simulated_density%new(self%ldim,self%smpd)
        if(nsteps ==1) then ! one fixed radius
          write(logfhandle, *) 'Fixed radius of ', min_rad
          radius = min_rad
          cnt_all  = 0
          do cc = 1, self%n_cc
              d = euclid(self%atominfo(cc)%center(:), m)*self%smpd
              ! Count nb of atoms in the sphere of radius radius
             if(d<=radius .and. cn_mask(cc)) then
                 cnt_all = cnt_all+1
             endif
          enddo
          call radial_atoms_all%new (cnt_all,  dummy=.true.)
          if(radius < nano_diameter/2.)  then
              volume = 4./3.*pi*(radius*self%smpd)**3
          else
              volume = 4./3.*pi*(nano_diameter/2.*self%smpd)**3
          endif
          write(unit = filnum2, fmt = '(a,f4.1)')  'Radius:   ', radius
          write(unit = filnum2, fmt = '(a,f6.3)')  'Density: ', real(cnt_all)/volume
          cnt_all  = 0
          ! Save coords
          do cc = 1, self%n_cc
              d = euclid(self%atominfo(cc)%center(:), m)*self%smpd
             if(d<=radius .and. cn_mask(cc)) then
                 cnt_all = cnt_all+1
                 call radial_atoms_all%set_name(cnt_all,self%atom_name)
                 call radial_atoms_all%set_element(cnt_all,self%element)
                 call radial_atoms_all%set_coord(cnt_all,(self%atominfo(cc)%center(:)-1.)*self%smpd)
             endif
          enddo
          call radial_atoms_all%writepdb (trim(int2str(nint(radius)))//'radial_atoms_all')
          ! Generate distribution based on atomic position
          call radial_atoms_all%convolve(simulated_density, cutoff)
          call simulated_density%write(trim(int2str(nint(radius))//'density_atoms_all.mrc'))
          ! Estimation of avg distance and stdev among atoms in radial dependent shells
          write(unit = filnum, fmt = '(a)')   '*********************************************************'
          write(unit = filnum, fmt = '(a,a)') 'Estimation of atom-to-atom statistics in shell of radius ', trim(int2str(nint(radius)))
          allocate(coords(3,cnt_all), source = 0.)
          allocate(max_intensity(cnt_all),avg_intensity(cnt_all),stdev_intensity(cnt_all), source = 0. )
          allocate(ratios(cnt_all), source = 0. )
          cnt = 0
          call self%img%get_rmat_ptr(rmat)
          call self%img_cc%get_imat(imat_cc)
          do cc = 1, size(self%atominfo)
              d = euclid(self%atominfo(cc)%center(:), m)*self%smpd
              if(d<=radius .and. cn_mask(cc)) then
                cnt = cnt + 1
                coords(:3,cnt) = self%atominfo(cc)%center(:)
                call self%calc_aspect_ratio(cc, ratios(cnt),lld=.false., print_ar=.false.)
                where( imat_cc == cc ) mask = .true.
                max_intensity(cnt) = maxval(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
                avg_intensity(cnt) = sum(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
                avg_intensity(cnt) = avg_intensity(cnt)/real(count(mask))
                do i = 1, self%ldim(1)
                    do j = 1, self%ldim(2)
                        do k = 1, self%ldim(3)
                            if(mask(i,j,k)) stdev_intensity(cnt) = stdev_intensity(cnt) + (rmat(i,j,k)-avg_intensity(cnt))**2
                        enddo
                    enddo
                enddo
                if(count(mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))) > 1) then
                    stdev_intensity(cnt) = sqrt(stdev_intensity(cnt)/real(count(mask)-1))
                else ! atom composed by one voxel
                    stdev_intensity(cnt) = 0.
                endif
                write(unit = filnum, fmt = '(a,i3,a,f9.5,a,f9.5,a,f9.5,a,f9.5)')  'ATOM ', cc,' maxval ', max_intensity(cnt), '   avg ', avg_intensity(cnt), '   stdev ', stdev_intensity(cnt),' aspect ratio ', ratios(cnt)
                mask = .false. !Reset
              endif
          enddo
          avg_int   = sum(avg_intensity)/real(cnt)
          max_int   = maxval(max_intensity)
          stdev_int = 0.
          do i = 1, cnt
              stdev_int = stdev_int + (avg_intensity(i)-avg_int)**2
          enddo
          stdev_int = sqrt(stdev_int/real(cnt-1))
          write(unit = filnum, fmt = '(a,f9.5)') 'Maxval  int shell :', max_int
          write(unit = filnum, fmt = '(a,f9.5)') 'Average int shell :', avg_int
          write(unit = filnum, fmt = '(a,f9.5)') 'Stdev   int shell :', stdev_int
          call self%distances_distribution(file=.false., coords=coords)
          deallocate(coords)
          deallocate(avg_intensity)
          deallocate(max_intensity)
          deallocate(stdev_intensity)
          deallocate(ratios)
          call radial_atoms_all%kill
        else
          do l = 1, nsteps
              radius = min_rad+real(l-1)*step
              cnt_all  = 0
              cnt_just = 0
              do cc = 1, self%n_cc
                  d = euclid(self%atominfo(cc)%center(:), m)*self%smpd
                  ! Count nb of atoms in the sphere of radius radius
                 if(d<=radius .and. cn_mask(cc)) then
                     cnt_all = cnt_all+1
                     if(l == 1) then
                         cnt_just = cnt_all
                     else
                         if(d>radius-step .and. d<=radius) cnt_just = cnt_just+1
                     endif
                 endif
              enddo
              if(cnt_all == 0) then ! first selected radius, no atoms belonging to it
                write(logfhandle, *) 'WARNING! No atoms in radius ', nint(radius)
                cycle
              endif
              call radial_atoms_just%new(cnt_just, dummy=.true.)
              call radial_atoms_all%new (cnt_all,  dummy=.true.)
              if(radius < nano_diameter/2.)  then
                  volume = 4./3.*pi*(radius*self%smpd)**3
              else
                  volume = 4./3.*pi*(nano_diameter/2.*self%smpd)**3
              endif
              write(unit = filnum2, fmt = '(a,f4.1,a,f6.3)')  'Radius: ', radius, ' Density: ', real(cnt_all)/volume
              cnt_all  = 0
              cnt_just = 0
              ! Save coords
              do cc = 1, self%n_cc
                  d = euclid(self%atominfo(cc)%center(:), m)*self%smpd
                 if(d<=radius .and. cn_mask(cc)) then
                     cnt_all = cnt_all+1
                     call radial_atoms_all%set_name(cnt_all,self%atom_name)
                     call radial_atoms_all%set_element(cnt_all,self%element)
                     call radial_atoms_all%set_coord(cnt_all,(self%atominfo(cc)%center(:)-1.)*self%smpd)
                     if(l == 1) then
                         cnt_just= cnt_just+1
                         call radial_atoms_just%set_name(cnt_just,self%atom_name)
                         call radial_atoms_just%set_element(cnt_just,self%element)
                         call radial_atoms_just%set_coord(cnt_just,(self%atominfo(cc)%center(:)-1.)*self%smpd)
                     endif
                     if(d>(radius-step) .and. l > 1) then
                         cnt_just = cnt_just+1
                         call radial_atoms_just%set_name(cnt_just,self%atom_name)
                         call radial_atoms_just%set_element(cnt_just,self%element)
                         call radial_atoms_just%set_coord(cnt_just,(self%atominfo(cc)%center(:)-1.)*self%smpd)
                     endif
                 endif
              enddo
              call radial_atoms_just%writepdb(trim(int2str(nint(radius)))//'radial_atoms_just')
              call radial_atoms_all%writepdb (trim(int2str(nint(radius)))//'radial_atoms_all')
              ! Generate distribution based on atomic position
              call radial_atoms_just%convolve(simulated_density, cutoff)
              call simulated_density%write(trim(int2str(nint(radius))//'density_atoms_just.mrc'))
              call radial_atoms_all%convolve(simulated_density, cutoff)
              call simulated_density%write(trim(int2str(nint(radius))//'density_atoms_all.mrc'))
              ! Estimation of avg distance and stdev among atoms in radial dependent shells
              write(unit = filnum, fmt = '(a)')   '*********************************************************'
              write(unit = filnum, fmt = '(a,a)') 'Estimation of atom-to-atom statistics in shell of radius ', trim(int2str(nint(radius)))
              allocate(coords(3,cnt_all), source = 0.)
              allocate(max_intensity(cnt_all),avg_intensity(cnt_all),stdev_intensity(cnt_all), source = 0. )
              allocate(ratios(cnt_all), source = 0. )
              cnt = 0
              call self%img%get_rmat_ptr(rmat)
              call self%img_cc%get_imat(imat_cc)
              do cc = 1, size(self%atominfo)
                  d = euclid(self%atominfo(cc)%center(:), m)*self%smpd
                  if(d<=radius .and. cn_mask(cc)) then
                      if(l == 1) then
                          cnt = cnt + 1
                          coords(:3,cnt) = self%atominfo(cc)%center(:)
                          call self%calc_aspect_ratio(cc, ratios(cnt),lld=.false., print_ar=.false.)
                          where( imat_cc == cc ) mask = .true.
                          max_intensity(cnt) = maxval(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
                          avg_intensity(cnt) = sum(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
                          avg_intensity(cnt) = avg_intensity(cnt)/real(count(mask))
                          do i = 1, self%ldim(1)
                              do j = 1, self%ldim(2)
                                  do k = 1, self%ldim(3)
                                      if(mask(i,j,k)) stdev_intensity(cnt) = stdev_intensity(cnt) + (rmat(i,j,k)-avg_intensity(cnt))**2
                                  enddo
                              enddo
                          enddo
                          if(count(mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))) > 1) then
                              stdev_intensity(cnt) = sqrt(stdev_intensity(cnt)/real(count(mask)-1))
                          else ! atom composed by one voxel
                              stdev_intensity(cnt) = 0.
                          endif
                          write(unit = filnum, fmt = '(a,i3,a,f9.5,a,f9.5,a,f9.5,a,f9.5)')  'ATOM ', cc,' maxval ', max_intensity(cnt), '   avg ', avg_intensity(cnt), '   stdev ', stdev_intensity(cnt),' aspect ratio ', ratios(cnt)
                          mask = .false. !Reset
                      elseif(d>(radius-step)) then
                          cnt = cnt + 1
                          coords(:3,cnt) = self%atominfo(cc)%center(:)
                          call self%calc_aspect_ratio(cc, ratios(cnt),lld=.false., print_ar=.false.)
                          where( imat_cc == cc ) mask = .true.
                          max_intensity(cnt) = maxval(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
                          avg_intensity(cnt) = sum(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
                          avg_intensity(cnt) = avg_intensity(cnt)/real(count(mask))
                          do i = 1, self%ldim(1)
                              do j = 1, self%ldim(2)
                                  do k = 1, self%ldim(3)
                                      if(mask(i,j,k)) stdev_intensity(cnt) = stdev_intensity(cnt) + (rmat(i,j,k)-avg_intensity(cnt))**2
                                  enddo
                              enddo
                          enddo
                          if(count(mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))) > 1) then
                              stdev_intensity(cnt) = sqrt(stdev_intensity(cnt)/real(count(mask)-1))
                          else ! atom composed by one voxel
                              stdev_intensity(cnt) = 0.
                          endif
                          write(unit = filnum, fmt = '(a,i3,a,f9.5,a,f9.5,a,f9.5,a,f9.5)')  'ATOM ', cc,' maxval ', max_intensity(cnt), '   avg ', avg_intensity(cnt), '   stdev ', stdev_intensity(cnt), ' aspect ratio ', ratios(cnt)
                          mask = .false. !Reset
                      endif
                  endif
              enddo
              avg_int   = sum(avg_intensity)/real(cnt)
              max_int   = maxval(max_intensity)
              stdev_int = 0.
              do i = 1, cnt
                  stdev_int = stdev_int + (avg_intensity(i)-avg_int)**2
              enddo
              stdev_int = sqrt(stdev_int/real(cnt-1))
              write(unit = filnum, fmt = '(a,f9.5)') 'Maxval  int shell :', max_int
              write(unit = filnum, fmt = '(a,f9.5)') 'Average int shell :', avg_int
              write(unit = filnum, fmt = '(a,f9.5)') 'Stdev   int shell :', stdev_int
              call self%distances_distribution(file=.false., coords=coords)
              deallocate(coords)
              deallocate(avg_intensity)
              deallocate(max_intensity)
              deallocate(stdev_intensity)
              deallocate(ratios)
              call radial_atoms_all%kill
              call radial_atoms_just%kill
          enddo
        endif
        call fclose(filnum)
        call fclose(filnum2)
        deallocate(mask, imat_cc)
        call simulated_density%kill
        call self%distances_distribution(file=.false.) ! across the whole nano
        if(present(cn_min)) then
          write(logfhandle,*) 'Calculating stats for atoms with std cn in range ', cn_min, '-', cn_max
          if(present(cn)) then
            call self%atoms_stats(.true.,cn_min=cn_min,cn_max=cn_max,cn=cn)
          else
            call self%atoms_stats(.true.,cn_min=cn_min,cn_max=cn_max)
          endif
        else
          write(logfhandle,*) 'Calculating stats for all the atoms '
          if(present(cn)) then
            write(logfhandle,*) 'Fixing cn to ', cn, ' for dipole calculation'
            call self%atoms_stats(.true.,cn=cn)
          else
            call self%atoms_stats(.true.)
          endif
        endif
        ! Kill
        call self%kill
        if(present(cn_min)) then
          write(logfhandle, *) '****radial atom-to-atom distances estimation with cn std range, completed'
        else
          write(logfhandle, *) '****radial atom-to-atom distances estimation, completed'
        endif
   end subroutine radial_dependent_stats

   ! This subroutine calculates some statistics (min,max,avg,stdev)
   ! in the intensity gray level value of the nanoparticle
   ! map in each atom. It is likely that these statistics
   ! are going to be able to distinguish between the different
   ! atom compositions in heterogeneous nanoparticles.
   subroutine atoms_stats(self,print_file,cn_min,cn_max,cn)
       class(nanoparticle), intent(inout) :: self
       logical,             intent(in)    :: print_file    ! whether to print on a file all the calc stats
       integer, optional,   intent(in)    :: cn_min,cn_max ! min/max cn standard for calculation of stats in range
       integer, optional,   intent(in)    :: cn            ! given fixed std cn for dipole calculation
       real,    pointer     :: rmat(:,:,:), rmat_corr(:,:,:)
       logical, allocatable :: mask(:,:,:), cn_mask(:)
       integer, allocatable :: imat_cc(:,:,:), sz(:)
       real,    allocatable :: ar(:),diameter(:),avg_gr(:),stdev_gr(:)
       type(image)          :: phasecorr
       integer :: i, j, k,n_atom, filnum, io_stat, min_cn, max_cn
       real    :: max_intensity(self%n_cc), avg_intensity(self%n_cc), max_corr(self%n_cc), int_corr(self%n_cc)
       real    :: stdev_intensity(self%n_cc), radii(self%n_cc), int_intensity(self%n_cc) ! integrated intensity
       real    :: avg_int, stdev_int, max_int, radius
       real    :: a(3),m(3) ! center of mass of the nanoparticle
       write(logfhandle,*)'**atoms intensity statistics calculations init'
       call self%phasecorr_one_atom(phasecorr)
       call phasecorr%get_rmat_ptr(rmat_corr)
       call self%img%get_rmat_ptr(rmat)
       call self%img_cc%get_imat(imat_cc)
       allocate(mask(self%ldim(1),self%ldim(2),self%ldim(3)), source = .false.)
       self%n_cc = maxval(imat_cc) ! update number of connected components
       m = self%masscen()
       ! Generate mask to eliminate atoms outside the selected cn range for stats calculation
       if(present(cn_min)) then
           allocate(cn_mask(self%n_cc), source = .true.)
           do n_atom = 1, self%n_cc
               if( self%atominfo(n_atom)%cn_std > cn_max .or. self%atominfo(n_atom)%cn_std < cn_min ) cn_mask(n_atom) = .false.
           enddo
       else
           allocate(cn_mask(self%n_cc), source = .true.) ! consider all the atoms
       endif
       ! initialise
       max_intensity(:)   = 0.
       avg_intensity(:)   = 0.
       int_intensity(:)   = 0.
       stdev_intensity(:) = 0.
       ! Size of the atoms
       sz = self%img_cc%size_ccs()
       ! Aspect Ratios and Diameters
       call self%est_aspect_ratios(print_ar=.false.,ar=ar,diameter=diameter)
       ! Polarization angles
       if(present(cn)) then
         call self%search_polarization(output_files=.true.,cn=cn)
       else
         call self%search_polarization(output_files=.true.)
       endif
       ! Write on a file
       if(print_file) then
          call fopen(filnum, file ='AllAtomsStats.txt', iostat=io_stat)
          write(unit = filnum, fmt = '(a)') '______________________  LEGEND  ______________________'
          write(unit = filnum, fmt = '(a)') 'diameter: diameter of the atom'
          write(unit = filnum, fmt = '(a)') 'sz:       number of vxls that compose the atom'
          write(unit = filnum, fmt = '(a)') 'radius:   radius corresponding to the location of the atom'
          write(unit = filnum, fmt = '(a)') 'maxval:   maximum intensity value across the atom'
          write(unit = filnum, fmt = '(a)') 'maxcorr:  maximum correlation value across the atom'
          write(unit = filnum, fmt = '(a)') 'intval:   integrated intensity across the atom'
          write(unit = filnum, fmt = '(a)') 'intcorr:  integrated correlation across the atom'
          write(unit = filnum, fmt = '(a)') 'avgval:   average intensity across the atom'
          write(unit = filnum, fmt = '(a)') 'stdevval: standard deviation of the intensities across the atom'
          write(unit = filnum, fmt = '(a)') 'cn_std:   standard coordination number'
          write(unit = filnum, fmt = '(a)') 'cn_gen:   generalised coordination number'
          write(unit = filnum, fmt = '(a)') 'ar:       aspect ratio'
          write(unit = filnum, fmt = '(a)') 'polang:   polarization angle'
          write(unit = filnum, fmt = '(a)') '______________________END LEGEND______________________'
          write(unit = filnum, fmt = '(a)') '  '
       endif
       do n_atom = 1, self%n_cc
           where(imat_cc == n_atom) mask = .true.
           max_intensity(n_atom) = maxval(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
           int_intensity(n_atom) = sum(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
           max_corr(n_atom)      = maxval(rmat_corr(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
           int_corr(n_atom)      = sum(rmat_corr(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
           avg_intensity(n_atom) = int_intensity(n_atom)/real(count(mask))
           radii(n_atom)         = euclid(self%atominfo(n_atom)%center(:), m)*self%smpd
           do i = 1, self%ldim(1)
               do j = 1, self%ldim(2)
                   do k = 1, self%ldim(3)
                       if(mask(i,j,k)) stdev_intensity(n_atom) = stdev_intensity(n_atom) + (rmat(i,j,k)-avg_intensity(n_atom))**2
                   enddo
               enddo
           enddo
           if(count(mask) > 1) then
               stdev_intensity(n_atom) = sqrt(stdev_intensity(n_atom)/real(count(mask)-1))
           else ! atom composed by one voxel
               stdev_intensity(n_atom) = 0.
           endif
            if(print_file .and. cn_mask(n_atom)) then
              write(unit = filnum, fmt = '(a,i3,a)') '______________________ATOM # ', n_atom, ' ____________________________'
              write(unit = filnum, fmt = '(a,f9.5)') ' diameter ', diameter(n_atom)
              write(unit = filnum, fmt = '(a,i3)')   ' sz       ', sz(n_atom)
              write(unit = filnum, fmt = '(a,f9.5)') ' radius   ', radii(n_atom)
              write(unit = filnum, fmt = '(a,f9.5)') ' maxval   ', max_intensity(n_atom)
              write(unit = filnum, fmt = '(a,f9.5)') ' maxcorr  ', max_corr(n_atom)
              write(unit = filnum, fmt = '(a,f9.5)') ' intint   ', int_intensity(n_atom)
              write(unit = filnum, fmt = '(a,f9.5)') ' intcorr  ', int_corr(n_atom)
              write(unit = filnum, fmt = '(a,f9.5)') ' avgval   ', avg_intensity(n_atom)
              write(unit = filnum, fmt = '(a,f9.5)') ' stdevval ', stdev_intensity(n_atom)
              write(unit = filnum, fmt = '(a,i3)')   ' cn_std   ', self%atominfo(n_atom)%cn_std
              write(unit = filnum, fmt = '(a,f9.5)') ' cn_gen   ', self%atominfo(n_atom)%cn_gen
              write(unit = filnum, fmt = '(a,f9.5)') ' ar       ', ar(n_atom)
              write(unit = filnum, fmt = '(a,f9.5)') ' polang   ', self%atominfo(n_atom)%polar_angle
            endif
           mask = .false. !Reset
       enddo
       avg_int   = sum(avg_intensity, cn_mask)/real(count(cn_mask))
       max_int   = maxval(max_intensity, cn_mask)
       stdev_int = 0.
       do n_atom= 1, self%n_cc
           if(cn_mask(n_atom))stdev_int = stdev_int + (avg_intensity(n_atom)-avg_int)**2
       enddo
       stdev_int = sqrt(stdev_int/real(count(cn_mask)-1))
       if(print_file) write(unit = filnum, fmt = '(a,f9.5,a,f9.5,a,f9.5)') 'maxval_general ', max_int, '   avg_general ', avg_int, '   stdev_general ', stdev_int
       if(print_file) call fclose(filnum)
       ! print one by one the stats anyway
       ! call fopen(filnum, file='../'//'RadialPos.csv', iostat=io_stat)
       ! Avg intensity value and std deviation within each cn group.
       min_cn = minval(self%atominfo(:)%cn_std, cn_mask)
       max_cn = maxval(self%atominfo(:)%cn_std, cn_mask)
       allocate(avg_gr(max_cn-min_cn+1),stdev_gr(max_cn-min_cn+1),source=0.)
       call cn_analysis(avg_gr, stdev_gr)
       ! Report on a file
       call fopen(filnum, file='RadialPos.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'radius'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(real2str(radii(n_atom)))
       end do
       call fclose(filnum)
       call fopen(filnum, file='MaxIntensity.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'maxint'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(real2str(max_intensity(n_atom)))
       end do
       call fclose(filnum)
       call fopen(filnum, file='MaxCorr.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'maxcorr'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(real2str(max_corr(n_atom)))
       end do
       call fclose(filnum)
       call fopen(filnum, file='AvgIntensity.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'avgint'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(real2str(avg_intensity(n_atom)))
       end do
       call fclose(filnum)
       call fopen(filnum, file='IntIntensity.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'intint'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(real2str(int_intensity(n_atom)))
       end do
       call fclose(filnum)
       call fopen(filnum, file='IntCorr.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'intcorr'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(real2str(int_corr(n_atom)))
       end do
       call fclose(filnum)
       call fopen(filnum, file='Size.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'sz'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(int2str(sz(n_atom)))
       end do
       call fclose(filnum)
       call fopen(filnum, file='AspectRatio.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'ar'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(real2str(ar(n_atom)))
       end do
       call fclose(filnum)
       call fopen(filnum, file='Diameter.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'diam'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(real2str(diameter(n_atom)))
       end do
       call fclose(filnum)
       call fopen(filnum, file='PolarizationAngle.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'polar_ang'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(real2str(self%atominfo(n_atom)%polar_angle))
       end do
       call fclose(filnum)
       call fopen(filnum, file='CnStd.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'cn_std'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(int2str(self%atominfo(n_atom)%cn_std))
       end do
       call fclose(filnum)
       call fopen(filnum, file='CnGen.csv', action='readwrite', iostat=io_stat)
       write (filnum,*) 'cn_gen'
       do n_atom = 1, self%n_cc
           if(cn_mask(n_atom)) write (filnum,'(A)', advance='yes') trim(real2str(self%atominfo(n_atom)%cn_gen))
       end do
       call fclose(filnum)
       deallocate(sz,ar,diameter)
       call phasecorr%kill()
       write(logfhandle,*)'**atoms intensity statistics calculations completed'
     contains

       subroutine cn_analysis(avg_gr,stdev_gr)
         real, intent(inout) :: avg_gr(:)   ! avg intensity value per cn group
         real, intent(inout) :: stdev_gr(:) ! standard deviation per cn group
         integer :: i,cn,cnt,cnt_gr,filnum1,filnum2,filnum3
         cnt_gr = 0
         do cn = min_cn,max_cn
           cnt = 0  ! number of atoms with coord number = cn
           cnt_gr = cnt_gr + 1 ! number of groups of coordination number
           do i = 1, self%n_cc
             if( cn_mask(i) .and. self%atominfo(i)%cn_std == cn) then
               cnt    = cnt + 1
               avg_gr(cnt_gr) = avg_gr(cnt_gr) + max_intensity(i)
             endif
           enddo
           avg_gr(cnt_gr) = avg_gr(cnt_gr)/real(cnt)
         enddo
         ! standard deviation of the max intensity values in the coordination number group
         call fopen(filnum1, file='AvgMaxIntCnGroup.csv', iostat=io_stat)
         call fopen(filnum2, file='StdevMaxIntCnGroup.csv', iostat=io_stat)
         call fopen(filnum3, file='CnGroup.csv', iostat=io_stat)
         write (filnum1,*) 'avg_max_int'
         write (filnum2,*) 'stdev_max_int'
         write (filnum3,*) 'cn'
         cnt_gr = 0
         do cn = min_cn,max_cn
           cnt = 0  ! number of atoms with coord number = cn
           cnt_gr = cnt_gr + 1 ! number of groups of coordination number
           do i = 1, self%n_cc
             if(cn_mask(i) .and. self%atominfo(i)%cn_std == cn) then
               cnt    = cnt + 1
               stdev_gr(cnt_gr) = stdev_gr(cnt_gr) + (max_intensity(i)-avg_gr(cnt_gr))**2
             endif
           enddo
           stdev_gr(cnt_gr) = sqrt(stdev_gr(cnt_gr)/real(cnt))
           print *, 'cn: ', cn, 'nb of elements in the group ', cnt, 'avg_gr ', avg_gr(cnt_gr), 'cnt_gr ', cnt_gr, 'stdev_gr ', stdev_gr(cnt_gr)
           write (filnum1,'(A)', advance='yes') trim(real2str(avg_gr(cnt_gr)))
           write (filnum2,'(A)', advance='yes') trim(real2str(stdev_gr(cnt_gr)))
           write (filnum3,'(A)', advance='yes') trim(int2str(cn))
         enddo
         call fclose(filnum1)
         call fclose(filnum2)
         call fclose(filnum3)
       end subroutine cn_analysis
   end subroutine atoms_stats

   subroutine distances_distribution( self, file, coords, volume )
       class(nanoparticle), intent(inout) :: self
       logical,             intent(in)    :: file ! output on a file
       real,    optional,   intent(in)    :: coords(:,:)
       integer, optional,   intent(in)    :: volume
       real, allocatable :: dist(:), tmpcens(:,:)
       real    :: stdev, med
       integer :: i, n_discard, filnum, io_stat
       logical :: mask(self%n_cc)
       ! Initialisations
       mask       = .true.
       stdev      = 0.
       n_discard  = 0
       if( present(coords) )then
           allocate(dist(size(coords,dim=2)), source = 0.)
           tmpcens = self%atominfo2centers()
           do i = 1, size(coords,dim=2)
               dist(i) =  pixels_dist(coords(:,i), tmpcens, 'min', mask=mask)      ! Use all the atoms
               mask(:) = .true. ! restore
               ! Discard outliers
               if( dist(i) * self%smpd > 3. * self%theoretical_radius ) then       ! maximum interatomic distance
                   dist(i)   = 0.
                   n_discard = n_discard + 1
               else if( dist(i) * self%smpd < 1.5 * self%theoretical_radius ) then ! minimum interatomic distance
                   dist(i)   = 0.
                   n_discard = n_discard + 1
               endif
           enddo
           self%avg_bondl = (sum(dist) / real(size(coords,dim=2) - n_discard)) * self%smpd
           self%max_bondl = maxval(dist) * self%smpd
           self%min_bondl = minval(dist) * self%smpd
           do i = 1, size(coords,dim=2)
               if( dist(i) * self%smpd <= 3. * self%theoretical_radius) stdev = stdev + (dist(i) - self%avg_bondl)**2.
           enddo
           if( size(coords,dim=2) - 1 - n_discard > 0 )then
               stdev = sqrt(stdev / real(size(coords,dim=2) - 1 - n_discard))
           else
               stdev = 0.
           endif
           med             = median(dist)
           self%sdev_bondl = stdev * self%smpd
           self%med_bondl  = med   * self%smpd
           if( file )then
               call fopen(filnum, file='DistancesDistr.txt', iostat=io_stat)
               if( present(volume) )then
                   write(unit = filnum, fmt = '(a,a,a,f6.3,a)') 'Average dist atoms vol ', trim(int2str(volume)),':', self%avg_bondl,  ' A'
                   write(unit = filnum, fmt = '(a,a,a,f6.3,a)') 'StDev   dist atoms vol ', trim(int2str(volume)),':', self%sdev_bondl, ' A'
                   write(unit = filnum, fmt = '(a,a,a,f6.3,a)') 'Median  dist atoms vol ', trim(int2str(volume)),':', self%med_bondl,  ' A'
               else
                   write(unit = filnum, fmt = '(a,f6.3,a)') 'Average dist atoms: ', self%avg_bondl,  ' A'
                   write(unit = filnum, fmt = '(a,f6.3,a)') 'StDev   dist atoms: ', self%sdev_bondl, ' A'
                   write(unit = filnum, fmt = '(a,f6.3,a)') 'Median  dist atoms: ', self%med_bondl,  ' A'
               endif
               call fclose(filnum)
           endif
           deallocate(dist, tmpcens)
       else
           tmpcens = self%atominfo2centers()
           do i = 1, size(self%atominfo)
               self%atominfo(i)%bondl = pixels_dist(self%atominfo(i)%center(:), tmpcens, 'min', mask=mask) ! Use all the atoms
               mask(:) = .true. ! restore
               ! Discard outliers
               if( self%atominfo(i)%bondl * self%smpd > 3. * self%theoretical_radius )then                 ! maximum interatomic distance
                   self%atominfo(i)%bondl = 0.
                   n_discard              = n_discard + 1
               else if( self%atominfo(i)%bondl * self%smpd < 1.5 * self%theoretical_radius ) then          ! minimum interatomic distance
                   self%atominfo(i)%bondl = 0.
                   n_discard              = n_discard + 1
               endif
           enddo
           self%avg_bondl = (sum(self%atominfo(:)%bondl)/real(size(self%atominfo)-n_discard)) * self%smpd
           self%max_bondl = maxval(self%atominfo(:)%bondl) * self%smpd
           self%min_bondl = minval(self%atominfo(:)%bondl) * self%smpd
           do i = 1, size(self%atominfo)
               if(self%atominfo(i)%bondl * self%smpd <= 3. * self%theoretical_radius)&
               &stdev = stdev + (self%atominfo(i)%bondl-self%avg_bondl)**2.
           enddo
           stdev = sqrt(stdev/real(size(self%atominfo)-1-n_discard))
           med   = median(self%atominfo(:)%bondl)
           self%sdev_bondl = stdev * self%smpd
           self%med_bondl  = med   * self%smpd
           if(file) then
               call fopen(filnum, file='DistancesDistr.txt', iostat=io_stat)
               write(unit = filnum, fmt = '(a,f6.3,a)') 'Average dist atoms: ', self%avg_bondl, ' A'
               write(unit = filnum, fmt = '(a,f6.3,a)') 'StDev   dist atoms: ', self%sdev_bondl, ' A'
               write(unit = filnum, fmt = '(a,f6.3,a)') 'Median  dist atoms: ', self%med_bondl, ' A'
               call fclose(filnum)
           endif
           deallocate(tmpcens)
       endif
   end subroutine distances_distribution

   subroutine est_aspect_ratios(self, print_ar, ar, diameter)
       class(nanoparticle), intent(inout) :: self
       logical, optional,   intent(in)    :: print_ar !print longest/shortest dim and ratio
       real, allocatable, optional, intent(inout) :: ar(:), diameter(:)
       integer, allocatable :: imat_cc(:,:,:)
       integer              :: label, filnum, io_stat
       real                 :: avg_diameter, median_diameter, min_diameter, max_diameter, stdev_diameter
       call self%img_cc%get_imat(imat_cc)
       self%n_cc = maxval(imat_cc) ! update number of connected components
       if( .not. allocated(self%atominfo) )then
           allocate(self%atominfo(self%n_cc))
       else
           if( size(self%atominfo) /= self%n_cc )then
               deallocate(self%atominfo)
               allocate(self%atominfo(self%n_cc))
           endif
       endif
       call self%find_centers()
       do label = 1, self%n_cc
           call self%calc_aspect_ratio(label, self%atominfo(label)%aspect_ratio, lld=.true., ld=self%atominfo(label)%diam, print_ar=print_ar)
       enddo
       self%atominfo(:)%diam = 2.*self%atominfo(:)%diam ! radius --> diameter
       min_diameter = minval(self%atominfo(:)%diam)
       max_diameter = maxval(self%atominfo(:)%diam)
       median_diameter = median(self%atominfo(:)%diam)
       avg_diameter    = sum(self%atominfo(:)%diam)/real(self%n_cc)
       stdev_diameter = 0.
       do label = 1, self%n_cc
           stdev_diameter = stdev_diameter + (avg_diameter-self%atominfo(label)%diam)**2
       enddo
       stdev_diameter = sqrt(stdev_diameter/real(self%n_cc-1))
       if(present(ar) .and. present(diameter))then
         allocate(ar(self%n_cc),       source = self%atominfo(:)%aspect_ratio)
         allocate(diameter(self%n_cc), source = self%atominfo(:)%diam)
       endif
       if(DEBUG) then
           write(logfhandle,*) 'minimum  value diameter ', min_diameter, 'A'
           write(logfhandle,*) 'maximum  value diameter ', max_diameter, 'A'
           write(logfhandle,*) 'median   value diameter ', median_diameter, 'A'
           write(logfhandle,*) 'average  value diameter ', avg_diameter, 'A'
           write(logfhandle,*) 'stdev    value diameter ', stdev_diameter, 'A'
       endif
       if(print_ar) then
         call fopen(filnum, file='AspectRatio.csv', iostat=io_stat)
         write (filnum,*) 'ar'
         do label = 1, size(self%atominfo)
           write (filnum,'(A)', advance='yes') trim(real2str(self%atominfo(label)%aspect_ratio))
         end do
         call fclose(filnum)
       endif
       write(logfhandle,*)'**aspect ratio calculations completed'
   end subroutine est_aspect_ratios

   ! This subroutine takes in input a connected component (cc) image
   ! and the label of one of its ccs and calculates the aspect ratio of the cc,
   ! defined as the ratio of the width and the height.
   ! The idea behind this is that the center of the cc is calculated,
   ! then everything is deleted except the borders of the cc. Finally,
   ! in order to calculate the width and the height, the min/max
   ! distances between the center and the borders are calculated. The
   ! aspect ratio is the ratio of those 2 distances.
   subroutine calc_aspect_ratio( self, label, ratio, lld, ld, print_ar )
       class(nanoparticle), intent(inout) :: self
       integer,             intent(in)    :: label
       real,                intent(out)   :: ratio
       logical,             intent(in)    :: lld ! fill in self%loc_longest_dim
       real   , optional,   intent(out)   :: ld  ! longest dist
       logical, optional,   intent(in)    :: print_ar
       integer, allocatable :: pos(:,:)
       integer, allocatable :: imat_cc(:,:,:)
       logical, allocatable :: border(:,:,:)
       logical, allocatable :: mask_dist(:) ! for min and max dist calculation
       integer :: location(1)               ! location of vxls of the atom farthest from its center
       real    :: shortest_dist, longest_dist
       call self%img_cc%get_imat(imat_cc)
       call self%img_cc%border_mask(border, label, .true.) ! use 4neigh instead of 8neigh
       where(border .eqv. .true.)
           imat_cc = 1
       elsewhere
           imat_cc = 0
       endwhere
       call get_pixel_pos( imat_cc, pos ) ! pxls positions of the shell
       allocate(mask_dist(size(pos, dim=2)), source = .true.)
       shortest_dist = pixels_dist(self%atominfo(label)%center(:), real(pos), 'min', mask_dist, location)
       if(size(pos,2) == 1) then ! if the connected component has size 1 (just 1 vxl)
           shortest_dist = 0.
           longest_dist  = shortest_dist
           ratio         = 1.
           if( lld ) self%atominfo(label)%loc_ldist(:) = nint(self%atominfo(label)%center(:))
           if( present(print_ar) )then
               if( print_ar .eqv. .true. ) then
                   write(logfhandle,*) 'ATOM #          ', label
                   write(logfhandle,*) 'shortest dist = ', shortest_dist
                   write(logfhandle,*) 'longest  dist = ', longest_dist
                   write(logfhandle,*) 'RATIO         = ', ratio
               endif
           endif
           return
       else
           longest_dist  = pixels_dist(self%atominfo(label)%center(:), real(pos),'max', mask_dist, location)
           if( lld ) self%atominfo(label)%loc_ldist(:) =  pos(:,location(1))
       endif
       if( abs(longest_dist) > TINY .and. size(pos,2) > 1 )then
           ratio = shortest_dist/longest_dist
       else
            ratio = 0.
            if( DEBUG ) write(logfhandle,*) 'cc ', label, 'LONGEST DIST = 0'
       endif
       longest_dist  = longest_dist  * self%smpd  ! in A
       shortest_dist = shortest_dist * self%smpd
       if( present(print_ar) )then
           if( print_ar .eqv. .true. )then
               write(logfhandle,*) 'ATOM #          ', label
               write(logfhandle,*) 'shortest dist = ', shortest_dist
               write(logfhandle,*) 'longest  dist = ', longest_dist
               write(logfhandle,*) 'RATIO         = ', ratio
           endif
       endif
       if( present(ld) ) ld = longest_dist
       deallocate(imat_cc, border, pos, mask_dist)
   end subroutine calc_aspect_ratio

    ! This subroutine clusters the atoms with respect to the maximum intensity
    ! or the integrated density (according to the values contained in feature)
    ! using kmeans algorithm with 2 classes. The initial guess fo the centers
    ! is intentionally biased. It supposes there are two distinguished classes
    ! with different avgs (proved with simulated data).
    subroutine cluster_atom_intensity( self, feature )
        use gnufor2
        use simple_nanoML, only : nanoML
        class(nanoparticle), intent(inout) :: self
        real,                intent(inout) :: feature(:)
        integer, parameter   :: MAX_IT = 50 ! maximum number of iterations for
        integer, allocatable :: imat_cc(:,:,:)
        real, pointer        :: rmat1(:,:,:), rmat2(:,:,:)
        type(image)          :: class1, class2
        type(nanoML)         :: emfit
        real    :: centers_kmeans(2) ! output of k-means
        real    :: avgs(2),  vars(2), gammas(self%n_cc,2) ! output of ML
        integer :: i, cnt1, cnt2, filnum, io_stat
        feature = feature/maxval(feature)*10. ! normalise with maxval*10
        call hist(feature, 20)
        write(logfhandle,*) '****clustering wrt maximum intensity, init'
        ! Report clusters on images in dedicated directory
        call class1%new(self%ldim, self%smpd)
        call class2%new(self%ldim, self%smpd)
        call class1%get_rmat_ptr(rmat1)
        call class2%get_rmat_ptr(rmat2)
        call self%img_cc%get_imat(imat_cc)
        call fopen(filnum, file='ClusterIntensities.txt', iostat=io_stat)
        ! kmeans
        call emfit%kmeans_biased2classes(feature, centers_kmeans)
        write(filnum,*) 'centers_kmeans', centers_kmeans
        ! ML, fit
        call emfit%new(self%n_cc,2)
        call emfit%set_data(feature)
        call emfit%fit(MAX_IT,centers_kmeans)
        avgs   = emfit%get_avgs()
        vars   = emfit%get_vars()
        gammas = emfit%get_gammas()
        write(filnum,*) 'AVG/VAR 1:', avgs(1), vars(1)
        write(filnum,*) 'AVG/VAR 2:', avgs(2), vars(2)
        cnt2 = 0
        do i = 1, self%n_cc
            if( (avgs(1) - feature(i))**2. < (avgs(2) - feature(i))**2. ) then
                cnt2 = cnt2 + 1
                write(filnum,*) 'connected component #', i, 'belongs to class 1 with probability', max(gammas(i,1),gammas(i,2))
              else
                write(filnum,*) 'connected component #', i, 'belongs to class 2 with probability', max(gammas(i,1),gammas(i,2))
            endif
        enddo
        cnt1 = count((avgs(1) - feature)**2. <  (avgs(2) - feature)**2. )
        cnt2 = count((avgs(1) - feature)**2. >= (avgs(2) - feature)**2. )
        do i = 1, self%n_cc
            if( (avgs(1) - feature(i))**2. < (avgs(2) - feature(i))**2. ) then
                where( imat_cc == i ) rmat1(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 1.
            else
                where( imat_cc == i ) rmat2(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 1.
            endif
        enddo
        call class1%write('Class1.mrc')
        call class2%write('Class2.mrc')
        write(filnum,*)  'Class1 contains ', cnt1 ,' atoms'
        write(filnum,*)  'Class2 contains ', cnt2 ,' atoms'
        write(filnum,*)  'Total ', cnt1+cnt2, ' atoms'
        call fclose(filnum)
        call class1%kill
        call class2%kill
        deallocate(imat_cc)
        ! come back to root directory
        write(logfhandle,*) '****clustering wrt maximum intensity, completed'

    contains

        subroutine initialise_centers(data,cen1,cen2)
            real, intent(inout) :: cen1,cen2
            real, intent(inout) :: data(:)
            !>   rheapsort from numerical recepies (largest last)
            call hpsort(data)
            cen1 = sum(data(1:self%n_cc/2))/real(self%n_cc/2)
            cen2 = sum(data(self%n_cc/2+1:size(data)))/real(self%n_cc/2)
        end subroutine initialise_centers

        subroutine update_centers(cen1,cen2,converged,val1,val2)
            real,    intent(inout) :: cen1,cen2
            logical, intent(inout) :: converged
            integer, intent(inout) :: val1, val2
            integer :: i
            integer :: cnt1, cnt2
            real    :: sum1, sum2
            real :: cen1_new, cen2_new
            sum1 = 0.
            cnt1 = 0
            do i=1,self%n_cc
                if( (cen1-feature(i))**2. < (cen2-feature(i))**2. )then
                    cnt1 = cnt1 + 1 ! number of elements in cluster 1
                    sum1 = sum1 + feature(i)
                endif
            end do
            cnt2 = self%n_cc - cnt1       ! number of elements in cluster 2
            sum2 = sum(feature)- sum1
            cen1_new = sum1 / real(cnt1)
            cen2_new = sum2 / real(cnt2)
            if(abs(cen1_new - cen1) < TINY .and. abs(cen2_new - cen2) < TINY) then
                converged = .true.
            else
                converged = .false.
            endif
            ! update
            cen1 = cen1_new
            cen2 = cen2_new
            ! assign values to the centers
            if( cen1 > cen2 )then
                val1           = 1
                val2           = 0
            else
                val1           = 0
                val2           = 1
            endif
        end subroutine update_centers

    end subroutine cluster_atom_intensity

    ! Cluster the atoms wrt the maximum intensity
    ! k-means 2 classes + ML
    subroutine cluster_atom_maxint( self )
        class(nanoparticle), intent(inout) :: self
        real    :: max_intensity(self%n_cc)
        integer :: i, io_stat, filnum
        call fopen(filnum, file='MaxIntensity.csv', action='readwrite', iostat=io_stat)
        if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file MaxIntensity.csv Did you run atoms_stats?; cluster_atom_maxint')
        endif
        read(filnum,*) ! first line is variable name
        do i = 1, self%n_cc
            read(filnum,*) max_intensity(i)
        enddo
        call fclose(filnum)
        call self%cluster_atom_intensity(max_intensity)
    end subroutine cluster_atom_maxint

    ! Cluster the atoms wrt the integrated density
    ! k-means 2 classes + ML
    subroutine cluster_atom_intint( self )
        class(nanoparticle), intent(inout) :: self
        real    :: int_intensity(self%n_cc)
        integer :: i, io_stat, filnum
        call fopen(filnum, file='IntIntensity.csv', action='readwrite', iostat=io_stat)
        if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file IntIntensity.csv Did you run atoms_stats?; cluster_atom_intint')
        endif
        read(filnum,*) ! first line is variable name
        do i = 1, self%n_cc
            read(filnum,*) int_intensity(i)
        enddo
        call fclose(filnum)
        call self%cluster_atom_intensity(int_intensity)
    end subroutine cluster_atom_intint

    ! Polarization search. The considered angle is
    ! the angle between the vector [0,0,1] and the direction of the
    ! longest dim of the atom.
    ! This is done for each of the atoms and then there is a search
    ! of the similarity between the calculated angles.
    subroutine search_polarization(self, output_files, cn)
        class(nanoparticle), intent(inout) :: self
        logical,             intent(in)    :: output_files
        integer, optional,   intent(in)    :: cn ! calculate net polarization for given std cn
        real    :: loc_ld_real(3,self%n_cc)
        real    :: vec_fixed(3)                  ! vector indentifying z direction
        real    :: m_adjusted(3), radius, a(3)
        integer :: k, i, filnum, io_stat
        integer, allocatable :: sz(:)
        real,    allocatable :: centers_A(:,:), tmpcens(:,:)
        logical :: cn_mask(self%n_cc)
        ! Calculate cn and cn_gen
        centers_A = self%atominfo2centers_A()
        call fit_lattice(centers_A,a)
        call find_radius_for_coord_number(a,radius)
        call run_coord_number_analysis(centers_A,radius,self%atominfo(:)%cn_std,self%atominfo(:)%cn_gen)
        deallocate(centers_A)
        ! Generate mask for cn
        if( present(cn) )then
            where(self%atominfo(:)%cn_std .ne. cn)
                cn_mask = .false.
            elsewhere
                cn_mask = .true.
            endwhere
            write(logfhandle,*) count(cn_mask), ' atoms out of ', self%n_cc, ' have cn ', cn
        else
            cn_mask(:) = .true.
        endif
        self%atominfo(:)%polar_angle = 0.
        sz = self%img_cc%size_ccs()
        ! bring vector back to center
        do i = 1, self%n_cc
            loc_ld_real(:3,i) = real(self%atominfo(i)%loc_ldist(:))- self%atominfo(i)%center(:)
        enddo
        ! consider fixed vector [0,0,1] (z direction)
        vec_fixed(1) = 0.
        vec_fixed(2) = 0.
        vec_fixed(3) = 1.
        write(logfhandle,*)'>>>>>>>>>>>>>>>> calculating angles wrt the vector [0,0,1]'
        self%net_dipole(1:3) = 0. ! inizialization
        do k = 1, self%n_cc
            if(sz(k) > 2) then
                self%atominfo(k)%polar_angle = ang3D_vecs(vec_fixed(:),loc_ld_real(:,k))
                if( cn_mask(k) ) self%net_dipole(:) = self%net_dipole(:) + loc_ld_real(:,k)
            else
                self%atominfo(k)%polar_angle = 0. ! If the cc is too small it doesn't make sense
            endif
            if( DEBUG ) write(logfhandle,*) 'ATOM ', k, 'angle between direction longest dim and vec [0,0,1] ', self%atominfo(k)%polar_angle
        enddo
        if(count(cn_mask) == 0) then
            write(logfhandle,*) 'No atoms with cn ', cn
            m_adjusted = 0.
            self%net_dipole = 0.
        else
            tmpcens = self%atominfo2centers()
            m_adjusted = sum(tmpcens(:,:)-1.,dim=2)*self%smpd/real(self%n_cc)
            self%net_dipole = (self%net_dipole-1.)*self%smpd + m_adjusted
            deallocate(tmpcens)
        endif
        if( output_files )then
            if( present(cn) )then
                write(logfhandle, *) 'Generating output dipole file for fixed cn ', cn
                call fopen(filnum, file='NetDipole'//int2str(cn)//'.bild', iostat=io_stat)
            else
                call fopen(filnum, file='NetDipole.bild', iostat=io_stat)
            endif
            write (filnum,'(a6,6i4)') trim('.arrow '), nint(m_adjusted), nint(self%net_dipole)
            call fclose(filnum)
        endif

        contains

            ! This subroutine takes in input 2 3D vectors, centered in the origin
            ! and it gives as an output the angle between them, IN DEGREES.
            function ang3D_vecs( vec1, vec2 ) result( ang )
                real, intent(inout) :: vec1(3), vec2(3)
                real :: ang        ! output angle
                real :: ang_rad    ! angle in radians
                real :: mod1, mod2
                real :: dot_prod
                mod1 = sqrt(vec1(1)**2+vec1(2)**2+vec1(3)**2)
                mod2 = sqrt(vec2(1)**2+vec2(2)**2+vec2(3)**2)
                ! normalise
                vec1 = vec1/mod1
                vec2 = vec2/mod2
                ! dot product
                dot_prod = (vec1(1)*vec2(1))+(vec1(2)*vec2(2))+(vec1(3)*vec2(3))
                ! sanity check
                if(dot_prod > 1. .or. dot_prod < -1.) then
                    THROW_WARN('Out of the domain of definition of arccos; search_polarization :: ang3D_vecs')
                    ang_rad = 0.
                else
                    ang_rad = acos(dot_prod)
                endif
                if(DEBUG) then
                    write(logfhandle,*)'>>>>>>>>>>>>>>>>>>>>>>'
                    write(logfhandle,*)'mod_1     = ', mod1
                    write(logfhandle,*)'mod_2     = ', mod2
                    write(logfhandle,*)'dot_prod  = ', dot_prod
                    write(logfhandle,*)'ang in radians', acos(dot_prod)
                endif
                ang = rad2deg(ang_rad)
            end function ang3D_vecs

    end subroutine search_polarization

    ! Cluster wrt to the angle between the direction of the longest dimension and a fixed vector
    subroutine cluster_ang(self, thresh)
        class(nanoparticle), intent(inout) :: self
        real,    intent(in)  :: thresh ! threshold for class definition, user inputted
        real,    allocatable :: stdev_within(:), centroids(:)
        integer, allocatable :: labels(:), populations(:)
        integer              :: i, j, ncls, dim, filnum, io_stat
        real                 :: avg, stdev
        type(binimage)       :: img_1clss
        integer, allocatable :: imat_cc(:,:,:), imat_1clss(:,:,:)
        character(len=4)     :: str_thres
        ! Preparing for clustering
        ! Aspect ratios and loc_longest dist estimation
        call self%est_aspect_ratios(print_ar=.false.) ! it's needed to identify the dir of longest dim
        call self%search_polarization(output_files =.false.)
        dim = size(self%atominfo)
        allocate(labels(dim), source = 0)
        call hac_1d(self%atominfo(:)%polar_angle,thresh,labels,centroids,populations)
        ! Stats calculations
        ncls = maxval(labels)
        allocate(stdev_within(ncls), source = 0.)
        !stdev within the same class
        do i = 1, ncls
            do j = 1, dim
                if(labels(j) == i) stdev_within(i) = stdev_within(i) + (self%atominfo(j)%polar_angle-centroids(i))**2
            enddo
        enddo
        where (populations>1)
            stdev_within = sqrt(stdev_within/(real(populations)-1.))
        elsewhere
            stdev_within = 0.
        endwhere
        ! avg and stdev among different classes
        avg = 0.
        do i = 1, ncls
            avg = avg + centroids(i)
        enddo
        avg = avg/real(ncls)
        stdev = 0.
        if(ncls>1) then
          do i = 1, ncls
              stdev = stdev + (centroids(i) - avg)**2
          enddo
          stdev = sqrt(stdev/(real(ncls)-1.))
        endif
        ! Output on a file
        str_thres = trim(real2str(thresh))
        call fopen(filnum, file='ClusterAngThresh'//str_thres//'.txt', iostat=io_stat)
        write(unit = filnum,fmt ='(a,i2,a,f6.2)') 'NR OF IDENTIFIED CLUSTERS:', ncls, ' SELECTED THRESHOLD: ',  thresh
        write(unit = filnum,fmt ='(a)') 'CLASSIFICATION '
        do i = 1, dim
          write(unit = filnum,fmt ='(a,i3,a,f6.2,a,i3)') 'Atom #: ', i, '; data (deg): ', self%atominfo(i)%polar_angle, '; class: ', labels(i)
        enddo
        write(unit = filnum,fmt ='(a)') 'CLASS STATISTICS '
        do i = 1, ncls
            write(unit = filnum,fmt ='(a,i3,a,i3,a,f6.2,a,f6.2,a)') 'class: ', i, '; cardinality: ', populations(i), '; centroid: ', centroids(i), ' degrees; stdev within the class: ', stdev_within(i), ' degrees'
        enddo
        write(unit = filnum,fmt ='(a,f6.2,a,f6.2,a)') 'AVG among the classes: ', avg, ' degrees; STDEV among the classes: ', stdev, ' degrees'
        call fclose(filnum)
        call fopen(filnum, file='Ang.csv', iostat=io_stat)
        write(filnum,*) 'ang'
        do i  = 1, self%n_cc
          write(filnum,*) self%atominfo(i)%polar_angle
        enddo
        call fclose(filnum)
        !Generate one figure for each class
        if(GENERATE_FIGS) then
          call img_1clss%new_bimg(self%ldim, self%smpd)
          call self%img_cc%get_imat(imat_cc)
          allocate(imat_1clss(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
          do i = 1, ncls
              imat_1clss = 0
              do j = 1, dim
                  if(labels(j) == i) then
                      where(imat_cc == j) imat_1clss = i
                  endif
              enddo
              call img_1clss%set_imat(imat_1clss)
              call img_1clss%write_bimg('AngClass'//int2str(i)//'.mrc')
          enddo
          call img_1clss%kill_bimg
        endif
        if(allocated(stdev_within)) deallocate(stdev_within)
        deallocate(centroids, labels, populations)
      end subroutine cluster_ang

  ! Cluster the atoms wrt to the aspect ratio
  ! Additionally, calculate the median std cn
  ! (coordination number) in each class
  subroutine cluster_ar(self, thresh)
      class(nanoparticle), intent(inout) :: self
      real,    intent(in)  :: thresh ! threshold for class definition, user inputted
      real,    allocatable :: centroids(:)
      integer, allocatable :: labels(:), populations(:)
      real,    allocatable :: stdev_within(:), centers_A(:,:), cn_cls(:), median_cn(:)
      integer              :: i, j, ncls, dim, filnum, io_stat, cnt
      real                 :: avg, stdev, radius_cn, a(3)
      type(binimage)       :: img_1clss
      integer, allocatable :: imat_cc(:,:,:), imat_1clss(:,:,:)
      character(len=4)     :: str_thres
      ! Calculate cn and cn_gen
      centers_A = self%atominfo2centers_A()
      call fit_lattice(centers_A,a)
      call find_radius_for_coord_number(a,radius_cn)
      call run_coord_number_analysis(centers_A,radius_cn,self%atominfo(:)%cn_std,self%atominfo(:)%cn_gen)
      deallocate(centers_A)
      if(thresh > 1. .or. thresh < 0.) THROW_HARD('Invalid input threshold! AR is in [0,1]; cluster_ar')
      ! Preparing for clustering
      ! Aspect ratios calculations
      call self%est_aspect_ratios(print_ar=.false.) ! it's needed to identify the dir of longest dim
      dim = size(self%atominfo)
      allocate(labels(dim), source = 0)
      ! classify
      call hac_1d(self%atominfo(:)%aspect_ratio,thresh,labels,centroids,populations)
      ! Stats calculations
      ncls = maxval(labels)
      allocate(stdev_within(ncls), source = 0.)
      ! stdev within the same class
      do i = 1, ncls
          do j = 1, dim
              if(labels(j) == i) stdev_within(i) = stdev_within(i) + (self%atominfo(j)%aspect_ratio-centroids(i))**2
          enddo
      enddo
      where (populations>1)
          stdev_within = sqrt(stdev_within/(real(populations)-1.))
      elsewhere
          stdev_within = 0.
      endwhere
      ! avg and stdev among different classes
      avg = 0.
      do i = 1, ncls
          avg = avg + centroids(i)
      enddo
      avg = avg/real(ncls)
      stdev = 0.
      if(ncls>1) then
        do i = 1, ncls
            stdev = stdev + (centroids(i) - avg)**2
        enddo
        stdev = sqrt(stdev/(real(ncls)-1.))
      endif
      ! median std cn on each class
      allocate(median_cn(ncls), source = 0.)
      do i = 1, ncls
        if(allocated(cn_cls)) deallocate(cn_cls)
        allocate(cn_cls(populations(i)), source = 0.)
        cnt = 0
        do j = 1, dim
            if(labels(j) == i) then
              cnt = cnt + 1
              cn_cls(cnt) = self%atominfo(j)%cn_std
            endif
        enddo
        median_cn(i) = median(cn_cls)
      enddo
      ! Output on a file
      str_thres = trim(real2str(thresh))
      call fopen(filnum, file='ClusterARThresh'//str_thres//'.txt', iostat=io_stat)
      write(unit = filnum,fmt ='(a,i2,a,f6.2)') 'NR OF IDENTIFIED CLUSTERS:', ncls, ' SELECTED THRESHOLD: ',  thresh
      write(unit = filnum,fmt ='(a)') 'CLASSIFICATION '
      do i = 1, dim
        write(unit = filnum,fmt ='(a,i3,a,f6.2,a,i3,a,f6.2)') 'Atom #: ', i, '; data (adimensional): ', self%atominfo(i)%aspect_ratio, '; class: ', labels(i)
      enddo
      write(unit = filnum,fmt ='(a)') 'CLASS STATISTICS '
      do i = 1, ncls
        write(unit = filnum,fmt ='(a,i3,a,i3,a,f6.2,a,f6.2,a,f6.2)') 'class: ', i, '; cardinality: ', populations(i), '; centroid: ', centroids(i), '; stdev within the class: ', stdev_within(i), '; median std cn: ', median_cn(i)
      enddo
      write(unit = filnum,fmt ='(a,f6.2,a,f6.2)') 'AVG among the classes: ', avg, '; STDEV among the classes: ', stdev
      call fclose(filnum)
      ! Generate file that contains the calculater AR, that can be read afterwards
      call fopen(filnum, file='AspectRatio.csv', iostat=io_stat)
      write(filnum,*) 'ar'
      do i  = 1, self%n_cc
        write(filnum,*) self%atominfo(i)%aspect_ratio
      enddo
      call fclose(filnum)
      ! Generate one figure for each class
      if(GENERATE_FIGS) then
        call img_1clss%new_bimg(self%ldim, self%smpd)
        call self%img_cc%get_imat(imat_cc)
        allocate(imat_1clss(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
        do i = 1, ncls
            imat_1clss = 0
            do j = 1, dim
                if(labels(j) == i) then
                    where(imat_cc == j) imat_1clss = i
                endif
            enddo
            call img_1clss%set_imat(imat_1clss)
            call img_1clss%write_bimg('ArClass'//int2str(i)//'.mrc')
        enddo
        call img_1clss%kill_bimg
      endif
      if(allocated(stdev_within)) deallocate(stdev_within)
      deallocate(centroids, labels, populations)
    end subroutine cluster_ar

    ! Cluster the atoms wrt to the interatomic distances
    subroutine cluster_bondl(self, thresh)
        class(nanoparticle), intent(inout) :: self
        real,    intent(in)  :: thresh ! threshold for class definition, user inputted
        real,    allocatable :: centroids(:)
        integer, allocatable :: labels(:), populations(:)
        real,    allocatable :: stdev_within(:), avg_dist_cog(:)
        integer              :: i, j, ncls, dim, filnum, io_stat
        real                 :: avg, stdev, cog(3)
        type(binimage)       :: img_1clss
        integer, allocatable :: imat_cc(:,:,:), imat_1clss(:,:,:)
        character(len=4)     :: str_thres
        ! Preparing for clustering
        call self%distances_distribution(file=.true.)
        dim = size(self%atominfo)
        allocate(labels(dim), source = 0)
        self%atominfo(:)%bondl = self%atominfo(:)%bondl * self%smpd   ! pass from pixels to A
        ! classify
        call hac_1d(self%atominfo(:)%bondl, thresh, labels, centroids, populations)
        ! Stats calculations
        ncls = maxval(labels)
        allocate(stdev_within(ncls), source = 0.)
        allocate(avg_dist_cog(ncls), source = 0.)
        cog = self%masscen()
        ! stdev within the same class
        ! avg dist to the center of gravity of each class
        do i = 1, ncls
            do j = 1, dim
                if(labels(j) == i) then
                   stdev_within(i) = stdev_within(i) + (self%atominfo(j)%bondl - centroids(i))**2.
                   avg_dist_cog(i) = avg_dist_cog(i) + euclid(cog,self%atominfo(j)%center(:))
                endif
            enddo
        enddo
        avg_dist_cog = (avg_dist_cog*self%smpd)/real(populations) ! in A
        where (populations>1)
            stdev_within = sqrt(stdev_within/(real(populations)-1.))
        elsewhere
            stdev_within = 0.
        endwhere
        ! avg and stdev among different classes
        avg = 0.
        do i = 1, ncls
            avg = avg + centroids(i)
        enddo
        avg = avg/real(ncls)
        stdev = 0.
        if(ncls>1) then
          do i = 1, ncls
              stdev = stdev + (centroids(i) - avg)**2
          enddo
          stdev = sqrt(stdev/(real(ncls)-1.))
        endif
        ! Output on a file
        str_thres = trim(real2str(thresh))
        call fopen(filnum, file='ClusterInterDistThresh'//str_thres//'.txt', iostat=io_stat)
        write(unit = filnum,fmt ='(a,i2,a,f6.2)') 'NR OF IDENTIFIED CLUSTERS:', ncls, ' SELECTED THRESHOLD: ',  thresh
        write(unit = filnum,fmt ='(a)') 'CLASSIFICATION '
        do i = 1, dim
          write(unit = filnum,fmt ='(a,i3,a,f6.2,a,i3)') 'Atom #: ', i, '; data (A): ', self%atominfo(i)%bondl, '; class: ', labels(i)
        enddo
        write(unit = filnum,fmt ='(a)') 'CLASS STATISTICS '
        do i = 1, ncls
          write(unit = filnum,fmt ='(a,i3,a,i3,a,f6.2,a,f6.2,a)') 'class: ', i, '; cardinality: ', populations(i), '; centroid: ', centroids(i), ' A; stdev within the class: ', stdev_within(i), ' A'
        enddo
        do i = 1, ncls
          write(unit = filnum,fmt ='(a,i3,a,f6.2,a)') 'class: ', i, '; average distance to the center of gravity: ', avg_dist_cog(i), ' A'
        enddo
        write(unit = filnum,fmt ='(a,f6.2,a,f6.2,a)') 'AVG among the classes: ', avg, ' A; STDEV among the classes: ', stdev, ' A'
        call fclose(filnum)
        call fopen(filnum, file='Dist.csv', iostat=io_stat)
        write(filnum,*) 'dist'
        do i  = 1, self%n_cc
          write(filnum,*) self%atominfo(i)%bondl
        enddo
        call fclose(filnum)
        ! Generate one figure for each class
        if(GENERATE_FIGS) then
          call img_1clss%new_bimg(self%ldim, self%smpd)
          call self%img_cc%get_imat(imat_cc)
          allocate(imat_1clss(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
          do i = 1, ncls
              imat_1clss = 0
              do j = 1, dim
                  if(labels(j) == i) then
                      where(imat_cc == j) imat_1clss = i
                  endif
              enddo
              call img_1clss%set_imat(imat_1clss)
              call img_1clss%write_bimg('DistClass'//int2str(i)//'.mrc')
          enddo
          call img_1clss%kill_bimg
        endif
        if(allocated(stdev_within)) deallocate(stdev_within)
        deallocate(centroids, labels, populations)
      end subroutine cluster_bondl

    subroutine make_soft_mask(self) !change the name
        class(nanoparticle), intent(inout) :: self
        type(binimage) :: simulated_density
        type(image)    :: img_cos
        type(atoms)    :: atomic_pos
        call simulated_density%new_bimg(self%ldim, self%smpd)
        call atomic_pos%new(trim(self%fbody)//'_atom_centers.pdb')
        call atomic_pos%convolve(simulated_density, cutoff=8.*self%smpd)
        call simulated_density%grow_bins(nint(0.5*self%theoretical_radius/self%smpd)+1)
        call simulated_density%cos_edge(SOFT_EDGE, img_cos)
        call img_cos%write(trim(self%fbody)//'SoftMask.mrc')
        !kill
        call img_cos%kill
        call simulated_density%kill_bimg
        call atomic_pos%kill
    end subroutine make_soft_mask

    ! Detect atoms. User does NOT input threshold for binarization..
    ! User might have inputted threshold for outliers
    ! removal based on contact score.
    subroutine identify_atomic_pos(self, cn_thresh, cn_type )
      class(nanoparticle), intent(inout) :: self
      integer,             intent(in)    :: cn_thresh
      character(len=6),    intent(in)    :: cn_type ! use generalised cn or standard?
      ! Phase correlations approach
      call self%phasecorr_one_atom(self%img)
      ! Nanoparticle binarization
      call self%binarize()
      ! Outliers discarding
      if(cn_thresh > 0) then
        if(cn_type .eq. 'cn_gen') then
          call self%discard_outliers_cn(cn_thresh, .true.)
        elseif(cn_type .eq. 'cn_std') then
          call self%discard_outliers_cn(cn_thresh, .false.)
        else
          write(logfhandle,*) 'Invalid cn_type, proceeding with cn_gen'
          call self%discard_outliers_cn(cn_thresh, .true.)
        endif
      else
        call self%discard_outliers()
      endif
      if(GENERATE_FIGS) call self%img_bin%write(trim(self%fbody)//'AfterOutliersRemoval.mrc')
      ! Validation of the selected atomic positions
       call self%validate_atomic_positions()
       if(GENERATE_FIGS) call self%img_bin%write(trim(self%fbody)//'AfterAPValidation.mrc')
    end subroutine identify_atomic_pos

    ! Detect atoms. User DOES input threshold for binarization.
    ! No phasecorrelation filter is applied.
    ! User might have inputted threshold for outliers
    ! removal based on contact score.
    subroutine identify_atomic_pos_thresh(self, thresh, cn_thresh, cn_type)
      class(nanoparticle), intent(inout) :: self
      real,                intent(in)    :: thresh
      integer,             intent(in)    :: cn_thresh
      character(len=2),    intent(in)    :: cn_type ! use generalised or standard coord number?
      ! Nanoparticle binarization
      call self%img%binarize(thres=thresh,self_out=self%img_bin)
      call self%img_bin%set_imat()
      ! Find connected components
      call self%img_bin%find_ccs(self%img_cc)
      call self%update_self_ncc
      ! Find atom centers
      call self%find_centers(self%img_bin, self%img_cc)
      ! Outliers discarding
      if(cn_thresh>0) then
        if(cn_type .eq. 'cn_gen') then
          call self%discard_outliers_cn(cn_thresh, .true.)
        elseif(cn_type .eq. 'cn_std') then
          call self%discard_outliers_cn(cn_thresh, .false.)
        else
          write(logfhandle,*) 'Invalid cn_type, proceeding with cn_gen'
          call self%discard_outliers_cn(cn_thresh, .true.)
        endif
      else
        call self%discard_outliers()
      endif
      ! Validation of the selected atomic positions
       call self%validate_atomic_positions()
    end subroutine identify_atomic_pos_thresh

    subroutine keep_atomic_pos_at_radius(self, radius, element, fname)
      class(nanoparticle), intent(inout) :: self
      real,                intent(in)    :: radius
      character(len=2),    intent(inout) :: element
      character(len=*),    intent(in)    :: fname !name of the pdb file where to write the coords
      character(len=4) :: atom_name
      integer     :: n_cc,cnt
      real        :: m(3) !coords of the center of mass of the nano
      real        :: d    !distance of the atom from the center of mass
      type(atoms) :: radial_atom
      element = upperCase(element)
      select case(element)
      case('PT')
          atom_name     = ' PT '
          ! thoretical radius is already set
      case('PD')
          atom_name     = ' PD '
      case('FE')
          atom_name     = 'FE'
      case('AU')
          self%atom_name     = ' AU '
      case default
          THROW_HARD('Unknown atom element; keep_atomic_pos_at_radius')
     end select
      cnt = 0
      m = self%masscen()
      ! Count nb of atoms in the selected radius
      do n_cc = 1, self%n_cc
            d = euclid(self%atominfo(n_cc)%center(:), m)*self%smpd
           if(d<=radius) then
               cnt = cnt+1
           endif
      enddo
      call radial_atom%new (cnt, dummy=.true.)
      cnt = 0
      ! Fill in radial_atom with the atomic positions
      do n_cc = 1, self%n_cc
          d = euclid(self%atominfo(n_cc)%center(:), m)*self%smpd
          if(d<=radius) then
             cnt = cnt+1
             call radial_atom%set_element(cnt,element)
             call radial_atom%set_name(cnt,atom_name)
             call radial_atom%set_coord(cnt,(self%atominfo(n_cc)%center(:)-1.)*self%smpd)
           endif
      enddo
      call radial_atom%writepdb(fname)
    end subroutine keep_atomic_pos_at_radius

    ! Translate the identified atomic positions so that the center of mass
    ! of the nanoparticle coincides with its closest atom
    subroutine center_on_atom(self, pdbfile_in, pdbfile_out)
      class(nanoparticle), intent(inout) :: self
      character(len=*),    intent(in)    :: pdbfile_in
      character(len=*),    intent(inout) :: pdbfile_out
      type(atoms) :: atom_centers
      real        :: m(3), vec(3), d, d_before
      integer     :: i
      call atom_centers%new(pdbfile_in)
      m(:)  = self%masscen()
      d_before = huge(d_before)
      do i = 1, self%n_cc
        d = euclid(m,self%atominfo(i)%center(:))
        if(d < d_before) then
          vec(:)   = m(:) - self%atominfo(i)%center(:)
          d_before = d
        endif
      enddo
      do i = 1, self%n_cc
        self%atominfo(i)%center(:) = self%atominfo(i)%center(:) + vec
        call atom_centers%set_coord(i,(self%atominfo(i)%center(:)-1.)*self%smpd)
      enddo
      call atom_centers%writePDB(pdbfile_out)
      call atom_centers%kill
    end subroutine center_on_atom

    subroutine mask(self, msk_rad)
      class(nanoparticle), intent(inout) :: self
      real,                intent(in)    :: msk_rad
      call self%img%mask(msk_rad, 'soft')
      if(DEBUG) call self%img%write(trim(self%fbody)//'_masked.mrc')
    end subroutine mask

    subroutine geometry_analysis(self, pdbfile2, thresh)
      use simple_math, only : plane_from_points
      class(nanoparticle),        intent(inout) :: self
      character(len=*),           intent(in)    :: pdbfile2  ! atomic pos of the 2 (or 3) selected atoms
      real,             optional, intent(in)    :: thresh    ! for belonging/not belonging to the plane/column
      type(atoms)           :: init_atoms, final_atoms
      type(binimage)        :: img_out
      character(len=3)      :: aux_var
      integer, allocatable  :: imat_cc(:,:,:), imat(:,:,:)
      real,    allocatable  :: line(:,:), plane(:,:,:),points(:,:),distances_totheplane(:),pointsTrans(:,:)
      real,    allocatable  :: radii(:),max_intensity(:),w(:),v(:,:),d(:),distances_totheline(:)
      integer :: i, n, t, s, filnum, io_stat, cnt_intersect, cnt
      logical :: flag(self%n_cc)
      real    :: atom1(3), atom2(3), atom3(3), dir_1(3), dir_2(3), vec(3), m(3), dist_plane, dist_line
      real    :: t_vec(N_DISCRET), s_vec(N_DISCRET), denominator, prod(3), centroid(3), tthresh
      if(present(thresh)) then
          tthresh = thresh
      else
          tthresh = 0.9*self%theoretical_radius
      endif
      call init_atoms%new(pdbfile2)
      n = init_atoms%get_n()
      if(n < 2 .or. n > 3 ) THROW_HARD('Inputted pdb file contains the wrong number of atoms!; geometry_analysis')
      do i = 1, N_DISCRET/2
          t_vec(i) = -real(i)/10.
      enddo
      t_vec(N_DISCRET/2+1:N_DISCRET) = -t_vec(1:N_DISCRET/2)
      s_vec(:) = t_vec(:)
      call self%img_cc%get_imat(imat_cc)
      call img_out%new(self%ldim, self%smpd)
      allocate(imat(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
      flag(:) = .false.  ! initialization
      if(n == 2) then
          write(logfhandle,*)'COLUMN IDENTIFICATION, INITIATION'
          allocate(line(3, N_DISCRET), source = 0.)
          atom1(:) = init_atoms%get_coord(1)/self%smpd + 1.
          atom2(:) = init_atoms%get_coord(2)/self%smpd + 1.
          dir_1 = atom1-atom2
          do t = 1, N_DISCRET
            line(1,t) = atom1(1) + t_vec(t)* dir_1(1)
            line(2,t) = atom1(2) + t_vec(t)* dir_1(2)
            line(3,t) = atom1(3) + t_vec(t)* dir_1(3)
          enddo
          ! calculate how many atoms does the line intersect and flag them
          do i = 1, self%n_cc
            do t = 1, N_DISCRET
              dist_line = euclid(self%atominfo(i)%center(:),line(:3,t))
              if(dist_line*self%smpd <= tthresh) then ! it intersects atoms i
                flag(i) = .true. !flags also itself
              endif
            enddo
          enddo
          ! generate volume for visualisation
          imat          = 0
          cnt_intersect = 0
          call img_out%new_bimg(self%ldim, self%smpd)
          call final_atoms%new(count(flag), dummy=.true.)
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt_intersect = cnt_intersect + 1
              call final_atoms%set_name(cnt_intersect,self%atom_name)
              call final_atoms%set_element(cnt_intersect,self%element)
              call final_atoms%set_coord(cnt_intersect,(self%atominfo(i)%center(:)-1.)*self%smpd)
              where(imat_cc == i) imat = 1
            endif
          enddo
          call img_out%set_imat(imat)
          call img_out%write_bimg('ImageColumn.mrc')
          call final_atoms%writePDB('AtomColumn')
          call final_atoms%kill
          ! Find the line that best fits the atoms
          allocate(points(3,count(flag)), source = 0.)
          cnt = 0
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt = cnt + 1
              points(:3,cnt) = self%atominfo(i)%center(:)
            endif
          enddo
          ! svd fit (https://au.mathworks.com/matlabcentral/answers/424591-3d-best-fit-line)
          allocate(pointsTrans(count(flag),3), source = 0.) ! because svdcmp modifies its input
          ! translate
          centroid = sum(points, dim = 2)/real(count(flag)) ! it belongs to the line
          do i = 1, count(flag)
            pointsTrans(i,:3) = points(:3,i) - centroid(:3)
          enddo
          allocate(w(3), v(3,3), source = 0.)
          allocate(d(3), source = 0.)
          call svdcmp(pointsTrans,w,v)
          d = v(:,1)
          write(logfhandle,*) 'Directional vector of the line', d
          ! line
          ! line(1,t) = centroid(1) + t_vec(t)* d(1)
          ! line(2,t) = centroid(2) + t_vec(t)* d(2)
          ! line(3,t) = centroid(3) + t_vec(t)* d(3)
          ! calculate the distance to the points from the identified line
          allocate(distances_totheline(cnt), source = 0.)
          allocate(radii(cnt), source = 0.) ! which radius is the atom center belonging to
          denominator = sqrt(d(1)**2+d(2)**2+d(3)**2)
          m = self%masscen()
          cnt = count(flag)
          do i = 1, cnt
            vec  = centroid(:3)-points(:3,i)
            prod = cross(vec,d)
            distances_totheline(i) = sqrt(prod(1)**2+prod(2)**2+prod(3)**2)/denominator
            radii(i) = euclid(points(:,i), m)*self%smpd
          enddo
          distances_totheline = distances_totheline*self%smpd ! go to A
          call fopen(filnum, file='Radii.csv', iostat=io_stat)
          write (filnum,*) 'r'
          do i = 1, cnt
            write (filnum,'(A)', advance='yes') trim(real2str(radii(i)))
          end do
          call fclose(filnum)
          call fopen(filnum, file='DistancesToTheLine.csv',iostat=io_stat)
          write (filnum,*) 'd'
          do i = 1, cnt
            write (filnum,'(A)', advance='yes') trim(real2str(distances_totheline(i)))
          end do
          call fclose(filnum)
          ! Read ratios, ang_var, dists, max_intensity
          if(.not. allocated(self%atominfo)) allocate(self%atominfo(self%n_cc))
          call fopen(filnum, file='../'//'../'//'AspectRatio.csv', action='readwrite',iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file AspectRatio.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*) ! first line is variable name
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%aspect_ratio
          enddo
          call fclose(filnum)
          call fopen(filnum, file='../'//'../'//'Ang.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file Ang.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%polar_angle
          enddo
          call fclose(filnum)
          call fopen(filnum, file='../'//'../'//'Dist.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file Dist.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%bondl
          enddo
          call fclose(filnum)
          allocate(max_intensity(self%n_cc))
          call fopen(filnum, file='../'//'../'//'MaxIntensity.csv', action='readwrite', iostat=io_stat)
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) max_intensity(i)
          enddo
          call fclose(filnum)
          ! Output on Excel file all the stats on the atoms belonging to the plane
          cnt = 0
          call fopen(filnum, file='AtomColumnInfo.txt',iostat=io_stat)
          write (filnum,*) '        Atom #    ','   Ar        ','       Dist     ','         Ang     ','         MaxInt     '
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt = cnt + 1
              write (filnum,*) i, '   ', self%atominfo(i)%aspect_ratio, self%atominfo(i)%bondl, self%atominfo(i)%polar_angle, max_intensity(i)
            endif
          end do
          call fclose(filnum)
      elseif(n == 3) then
          write(logfhandle,*)'PLANE IDENTIFICATION, INITIATION'
          atom1(:) = init_atoms%get_coord(1)/self%smpd + 1.
          atom2(:) = init_atoms%get_coord(2)/self%smpd + 1.
          atom3(:) = init_atoms%get_coord(3)/self%smpd + 1.
          dir_1 = atom1-atom2
          dir_2 = atom1-atom3
          allocate(plane(3, N_DISCRET, N_DISCRET), source = 0.)
          do t = 1, N_DISCRET
            do s = 1, N_DISCRET
              plane(1,t,s) = atom1(1) + t_vec(t)* dir_1(1) + s_vec(s)* dir_2(1)
              plane(2,t,s) = atom1(2) + t_vec(t)* dir_1(2) + s_vec(s)* dir_2(2)
              plane(3,t,s) = atom1(3) + t_vec(t)* dir_1(3) + s_vec(s)* dir_2(3)
            enddo
          enddo
          ! calculate how many atoms does the plane intersect and flag them
          do i = 1, self%n_cc
            do t = 1, N_DISCRET
              do s = 1, N_DISCRET
                dist_plane = euclid(self%atominfo(i)%center(:),plane(:3,t,s))
                if(dist_plane*self%smpd <= tthresh) then ! it intersects atoms i
                  flag(i) = .true. !flags also itself
                endif
              enddo
            enddo
          enddo
          ! generate volume for visualisation
          ! reset
          imat    = 0
          cnt_intersect = 0
          call img_out%new_bimg(self%ldim, self%smpd)
          call final_atoms%new(count(flag), dummy=.true.)
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt_intersect = cnt_intersect + 1
              call final_atoms%set_name(cnt_intersect,self%atom_name)
              call final_atoms%set_element(cnt_intersect,self%element)
              call final_atoms%set_coord(cnt_intersect,(self%atominfo(i)%center(:)-1.)*self%smpd)
              where(imat_cc == i) imat = 1
            endif
          enddo
          call img_out%set_imat(imat)
          call img_out%write_bimg('ImagePlane.mrc')
          call final_atoms%writePDB('AtomPlane')
          call final_atoms%kill
          allocate(points(3, count(flag)), source = 0.)
          m = self%masscen()
          cnt = 0
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt = cnt + 1
              points(:3,cnt) = self%atominfo(i)%center(:)-m(:)
            endif
          enddo
          vec = plane_from_points(points)
          allocate(distances_totheplane(cnt), source = 0.)
          allocate(radii(cnt), source = 0.) ! which radius is the atom center belonging to
          cnt = 0
          denominator = sqrt(vec(1)**2+vec(2)**2+1.)
          write(logfhandle,*) 'Normal vector: [', vec(1), ',', vec(2), ',', -1., ']'
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt = cnt + 1
              ! formula for distance of a point to a plane
              distances_totheplane(cnt) = abs(vec(1)*points(1,cnt)+vec(2)*points(2,cnt)-points(3,cnt)+vec(3))/denominator
              radii(cnt) = euclid(self%atominfo(i)%center(:), m)*self%smpd
            endif
          enddo
          distances_totheplane = (distances_totheplane)*self%smpd
          call fopen(filnum, file='Radii.csv', action='readwrite', iostat=io_stat)
          write (filnum,*) 'r'
          do i = 1, cnt
            write (filnum,'(A)', advance='yes') trim(real2str(radii(i)))
          end do
          call fclose(filnum)
          call fopen(filnum, file='DistancesToThePlane.csv', action='readwrite',iostat=io_stat)
          write (filnum,*) 'd'
          do i = 1, cnt
            write (filnum,'(A)', advance='yes') trim(real2str(distances_totheplane(i)))
          end do
          call fclose(filnum)
          ! Read ratios, ang_var, dists, max_intensity
          if(.not. allocated(self%atominfo)) allocate(self%atominfo(self%n_cc))
          call fopen(filnum, file='../'//'../'//'AspectRatio.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file AspectRatio.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*)  ! first line is variable name
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%aspect_ratio
          enddo
          call fclose(filnum)
          call fopen(filnum, file='../'//'../'//'Ang.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file Ang.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%polar_angle
          enddo
          call fclose(filnum)
          call fopen(filnum, file='../'//'../'//'Dist.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file Dist.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%bondl
          enddo
          call fclose(filnum)
          allocate(max_intensity(self%n_cc))
          call fopen(filnum, file='../'//'../'//'MaxIntensity.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file MaxIntensity.csv Did you run atoms_stats?; geometry_analysis')
          endif
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) max_intensity(i)
          enddo
          call fclose(filnum)
          ! Output on txt file all the stats on the atoms belonging to the plane
          cnt = 0
          call fopen(filnum, file='AtomPlaneInfo.txt',action='write',iostat=io_stat)
          write (filnum,*) '        Atom #    ','   Ar        ','       Dist     ','         Ang     ','      MaxInt     '
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt = cnt + 1
              write (filnum,*) i, '   ', self%atominfo(i)%aspect_ratio, self%atominfo(i)%bondl, self%atominfo(i)%polar_angle, max_intensity(i)
            endif
          end do
          call fclose(filnum)
      endif
      call img_out%kill_bimg
      call init_atoms%kill
      if(allocated(line))  deallocate(line)
      if(allocated(plane)) deallocate(plane)
      deallocate(imat_cc, imat)

    contains

      ! Compute the cross product of 2 3D real vectors
      function cross(a, b) result(c)
          real, intent(in) :: a(3),b(3)
          real :: c(3)
          c(1) = a(2) * b(3) - a(3) * b(2)
          c(2) = a(3) * b(1) - a(1) * b(3)
          c(3) = a(1) * b(2) - a(2) * b(1)
      end function cross

    end subroutine geometry_analysis

    subroutine kill_nanoparticle(self)
        class(nanoparticle), intent(inout) :: self
        self%ldim(3)         = 0
        self%smpd            = 0.
        self%NPcen(3)        = 0.
        self%avg_bondl = 0.
        self%n_cc            = 0
        call self%img%kill()
        call self%img_raw%kill
        call self%img_bin%kill_bimg()
        call self%img_cc%kill_bimg()
        call self%centers_pdb%kill
        if( allocated(self%atominfo) ) deallocate(self%atominfo)
    end subroutine kill_nanoparticle

end module simple_nanoparticle
