module simple_nanoparticle
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,     only : image
use simple_binimage,  only : binimage
use simple_atoms,     only : atoms
implicit none

public :: nanoparticle
private
#include "simple_local_flags.inc"

! module global constants
integer, parameter :: N_THRESH   = 20      ! number of thresholds for binarization
logical, parameter :: DEBUG      = .false. ! for debugging purposes
integer, parameter :: SOFT_EDGE  = 6
integer, parameter :: N_DISCRET  = 1000

type :: nanoparticle
    private
    type(atoms)    :: centers_pdb
    type(image)    :: img,img_raw
    type(binimage) :: img_bin, img_cc
    integer        :: ldim(3)            = 0
    real           :: smpd               = 0.
    real           :: nanop_mass_cen(3)  = 0. !coordinates of the center of mass of the nanoparticle
    real           :: avg_dist_atoms     = 0.
    real           :: theoretical_radius = 0. !theoretical atom radius
    integer        :: n_cc               = 0  !number of atoms (connected components)
    real,    allocatable  :: centers(:,:)
    real,    allocatable  :: ratios(:)
    real,    allocatable  :: ang_var(:)
    real,    allocatable  :: dists(:)
    integer, allocatable  :: loc_longest_dist(:,:)   ! for indentific of the vxl that determins the longest dim of the atom
    integer, allocatable  :: contact_scores(:)
    character(len=2)      :: element     = ' '
    character(len=4)      :: atom_name   = '    '
    character(len=STDLEN) :: partname    = ''     ! fname
    character(len=STDLEN) :: fbody       = ''     ! fbody
  contains
    ! constructor
    procedure          :: new => new_nanoparticle
    ! getters/setters
    procedure          :: get_ldim
    procedure          :: get_natoms
    procedure          :: set_img
    procedure          :: set_partname
    procedure          :: set_atomic_coords
    ! segmentation and statistics
    procedure          :: binarize
    procedure          :: find_centers
    procedure          :: nanopart_masscen
    procedure          :: aspect_ratios_estimation
    procedure          :: calc_aspect_ratio
    procedure          :: discard_outliers
    procedure          :: atom_intensity_stats
    procedure          :: validate_atomic_positions
    procedure          :: distances_distribution
    ! phase correlation
    procedure          :: phasecorrelation_nano_gaussian
    ! clustering
    procedure, private :: hierarc_clust
    procedure, private :: cluster_atom_intensity
    procedure          :: cluster_ang
    procedure          :: cluster_ar
    procedure          :: cluster_interdist
    procedure          :: cluster_atom_maxint
    procedure          :: cluster_atom_intint
    procedure          :: search_polarization
    procedure          :: geometry_analysis
    ! execution
    procedure          :: identify_atomic_pos
    procedure          :: radial_dependent_stats
    procedure          :: atoms_rmsd
    procedure          :: make_soft_mask
    ! visualization and output
    procedure          :: write_centers
    procedure          :: print_asym_unit
    procedure          :: print_atoms
    ! others
    procedure          :: update_self_ncc
    procedure          :: keep_atomic_pos_at_radius
    procedure          :: center_on_atom
    ! kill
    procedure          :: kill => kill_nanoparticle
end type nanoparticle

contains

    subroutine new_nanoparticle(self, fname, cline_smpd, element)
        class(nanoparticle),        intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        real,                       intent(in)    :: cline_smpd
        character(len=2), optional, intent(inout) :: element
        integer :: nptcls
        real    :: smpd
        call self%kill
        call self%set_partname(fname)
        self%fbody = get_fbody(trim(basename(fname)), trim(fname2ext(fname)))
        self%smpd  = cline_smpd
        ! CHIARA, the below should be in a defs module atom_constants or similar, discuss with C# please
        if(.not. present(element)) then
            self%element   = 'PT'  !default is pt
            self%atom_name = ' PT '
            self%theoretical_radius = 1.1
        else
            element = upperCase(element)
            select case(element)
                case('PT')
                    self%element       = 'PT'
                    self%atom_name     = ' PT '
                    self%theoretical_radius = 1.1
                    ! thoretical radius is already set
                case('PD')
                    self%element       = 'PD'
                    self%atom_name     = ' PD '
                    self%theoretical_radius = 1.12
                case('FE')
                    self%element       = 'FE'
                    self%atom_name     = ' FE '
                    self%theoretical_radius = 1.02
                case('AU')
                    self%element       = 'AU'
                    self%atom_name     = ' AU '
                    self%theoretical_radius = 1.23
                case default
                    THROW_HARD('Unknown atom element; new_nanoparticle')
           end select
        endif
        !self%sigma = 0.8*theoretical_radius/(2.*sqrt(2.*log(2.))*self%smpd) !0.8 not to have it too big (avoid connecting atoms)
        call find_ldim_nptcls(self%partname,  self%ldim, nptcls, smpd)
        call self%img%new(self%ldim, self%smpd)
        call self%img_raw%new(self%ldim, self%smpd)
        call self%img_bin%new(self%ldim, self%smpd)
        call self%img%read(fname)
        call self%img_raw%read(fname)
    end subroutine new_nanoparticle

     subroutine get_ldim(self,ldim) !TO REMOVE??
       class(nanoparticle), intent(in)  :: self
       integer,             intent(out) :: ldim(3)
       ldim = self%img%get_ldim()
     end subroutine get_ldim

     function get_natoms(self) result(n)
       class(nanoparticle), intent(inout)  :: self
       integer :: n
       call self%img_cc%get_nccs(n)
     end function get_natoms

     !set one of the images of the nanoparticle type
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
        integer, allocatable :: imat_cc(:,:,:)
        if(present(img_cc)) then
          call img_cc%get_nccs(self%n_cc)
        else
          call self%img_cc%get_nccs(self%n_cc)
        endif
    end subroutine update_self_ncc

    ! This subroutine has visualization purpose only.
    ! It prints out the 3 asymmetric units in a nanoparticle
    ! with c3 symmetry.
    subroutine print_asym_unit(self,nb_asym_units)
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: nb_asym_units
        real, allocatable :: rmat1(:,:,:), rmat2(:,:,:), rmat3(:,:,:)
        type(image)   :: img_asym
        integer       :: i, j
        integer       :: x, y
        real :: vec(2), vec1(2)
        real :: ang1, ang2
        rmat1 = self%img%get_rmat()
        rmat2 = self%img%get_rmat()
        rmat3 = self%img%get_rmat()
        x = self%ldim(1)/2 !x-coords of the center
        y = self%ldim(2)/2 !y-coords of the center
        vec1(:)= [0.,1.]   !fixed 2D vector
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                vec(:) = real([i-x,j-y]) !translate the vector in the origin
                ang1 = ang2D_vecs(vec1,vec)
                if(ang1 > 360./real(nb_asym_units) .or. i <= x) then
                    rmat1(i,j,:) = 0.    !set to 0 the all line
                endif
                if(ang1 > 360./real(nb_asym_units) .or. i > x) then
                    rmat2(i,j,:) = 0.    !set to 0 the all line
                endif
                if(ang1 < 360./real(nb_asym_units)) rmat3(i,j,:) = 0.
             enddo
        enddo
        call img_asym%new(self%ldim, self%smpd)
        call img_asym%set_rmat(rmat1)
        call img_asym%write(trim(self%fbody)//'FirstAsymUnit.mrc')
        call img_asym%set_rmat(rmat2)
        call img_asym%write(trim(self%fbody)//'SecondAsymUnit.mrc')
        call img_asym%set_rmat(rmat3)
        call img_asym%write(trim(self%fbody)//'ThirdAsymUnit.mrc')
        call img_asym%kill
        deallocate(rmat1,rmat2,rmat3)
    end subroutine print_asym_unit

    subroutine print_atoms(self)
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
       class(nanoparticle),      intent(inout) :: self
       type(binimage), optional,     intent(inout) :: img_bin, img_cc
       real, optional, allocatable,  intent(out)   :: coords(:,:)
       real,        pointer :: rmat_raw(:,:,:)
       integer, allocatable :: imat_cc_in(:,:,:)
       logical, allocatable :: mask(:,:,:)
       integer :: i, ii, jj, kk
       real    :: m(3), sum_mass
       ! sanity check
       if(present(img_bin) .and. .not. present(img_cc)) THROW_HARD('img_bin and img_cc have to be both present in input')
       ! global variables allocation
       if(allocated(self%centers)) deallocate(self%centers)
       allocate(self%centers(3,self%n_cc),source = 0.)
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
           self%centers(:,i) = m/sum_mass
       enddo
       !$omp end parallel do
       ! saving centers coordinates, optional
       if(present(coords)) allocate(coords(3,self%n_cc), source = self%centers)
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
                call self%centers_pdb%set_coord(cc,(self%centers(:,cc)-1.)*self%smpd)
                call self%centers_pdb%set_beta(cc,real(self%contact_scores(cc)))
                call self%centers_pdb%set_resnum(cc,cc)
            enddo
        endif
        if(present(fname)) then
            call self%centers_pdb%writepdb(fname)
        else
            call self%centers_pdb%writepdb(trim(self%fbody)//'_atom_centers')
        endif
    end subroutine write_centers

    ! This subroutine sets the atom positions to be the ones
    ! indicated in the inputted PDB file.
    subroutine set_atomic_coords(self, pdb_file)
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: pdb_file
        type(atoms) :: a
        integer     :: n_atoms, N
        if(fname2ext(pdb_file) .ne. 'pdb') THROW_HARD('Inputted filename has to be PDB file!; set_atomic_coords')
        if(allocated(self%centers)) deallocate(self%centers)
        call a%new(pdb_file)
        N = a%get_n() ! number of atoms
        allocate(self%centers(3,N), source = 0.)
        do n_atoms = 1, N
            self%centers(:,n_atoms) = a%get_coord(n_atoms)/self%smpd + 1.
        enddo
        self%n_cc = N
        call a%kill
    end subroutine set_atomic_coords

    ! calc the avg of the centers coords
     function nanopart_masscen(self) result(m)
         class(nanoparticle), intent(inout) :: self
         real    :: m(3)  !mass center coords
         integer :: i, j, k
          m = sum(self%centers(:,:), dim=2) /real(self%n_cc)
     end function nanopart_masscen

    ! This subroutine takes in input 2 2D vectors, centered in the origin
    ! and it gives as an output the angle between them, IN DEGREES.
    function ang2D_vecs(vec1, vec2) result(ang)
        real, intent(inout) :: vec1(2), vec2(2)
        real :: ang        !output angle
        real :: ang_rad    !angle in radians
        real :: mod1, mod2, dot_prod
        mod1 = sqrt(vec1(1)**2+vec1(2)**2)
        mod2 = sqrt(vec2(1)**2+vec2(2)**2)
        ! normalise
        vec1 = vec1/mod1
        vec2 = vec2/mod2
        ! dot product
        dot_prod = vec1(1)*vec2(1)+vec1(2)*vec2(2)
        ! sanity check
        if(dot_prod > 1. .or. dot_prod< -1.) then
            THROW_WARN('Out of the domain of definition of arccos; ang2D_vecs')
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
        !output angle
        ang = rad2deg(ang_rad)
    end function ang2D_vecs

    ! This subroutine takes in input 2 3D vectors, centered in the origin
    ! and it gives as an output the angle between them, IN DEGREES.
    function ang3D_vecs(vec1, vec2) result(ang)
        real, intent(inout) :: vec1(3), vec2(3)
        real :: ang        !output angle
        real :: ang_rad    !angle in radians
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
            THROW_WARN('Out of the domain of definition of arccos; ang3D_vecs')
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

    ! FORMULA: phasecorr = ifft2(fft2(field).*conj(fft2(reference)));
    subroutine phasecorrelation_nano_gaussian(self)
        use simple_segmentation
        class(nanoparticle), intent(inout) :: self
        type(image) :: one_atom, phasecorr
        type(atoms) :: atom
        real :: cutoff, o_t
        call phasecorr%new(self%ldim, self%smpd)
        call phasecorr%set_ft(.true.)
        call one_atom%new(self%ldim,self%smpd)
        cutoff = 8.*self%smpd
        call atom%new(1)
        call atom%set_element(1,self%element)
        call atom%set_coord(1,self%smpd*(real(self%ldim)/2.)) !DO NOT NEED THE +1
        call atom%convolve(one_atom, cutoff)
        call one_atom%fft()
        call self%img%fft()
        call self%img%phase_corr(one_atom,phasecorr,1.)
        call phasecorr%write(trim(self%fbody)//'CorrFiltered.mrc')
        call self%img%copy(phasecorr)
        call one_atom%kill()
        call phasecorr%kill()
    end subroutine phasecorrelation_nano_gaussian

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
        integer :: i, cc, t
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

    ! This subroutine discard outliers that resisted binarization.
    ! It calculates the contact score of each atom and discards the bottom
    ! PERCENT_DISCARD% of the atoms according to the contact score.
    ! It modifies the img_bin and img_cc instances deleting the
    ! identified outliers.
    ! Contact score: fix a radius and an atom A. Count the number N of atoms
    ! in the sphere centered in A of that radius. Define contact score of A
    ! equal to N. Do it for each atom.
    subroutine discard_outliers(self, cs_thresh)
        class(nanoparticle), intent(inout) :: self
        integer, optional,       intent(in)    :: cs_thresh ! threshold for discard outliers based on contact score
        integer, allocatable :: imat_bin(:,:,:), imat_cc(:,:,:)
        integer, parameter   :: PERCENT_DISCARD = 5/100
        integer, allocatable :: contact_scores(:)
        logical, allocatable :: mask(:)
        real    :: dist
        real    :: radius  !radius of the sphere to consider
        integer :: cc
        integer :: cnt     !contact_score
        integer :: loc(1)
        integer :: n_discard
        integer :: label(1)
        write(logfhandle, *) '****outliers discarding, init'
        ! Update number of connected components?
        ! Outliers removal using contact score
        radius = 2.*(2.*self%theoretical_radius)/self%smpd ! In pixels
        allocate(mask(self%n_cc), source = .true.)
        allocate(contact_scores(self%n_cc), source = 0)
        do cc = 1, self%n_cc !fix the atom
            dist = 0.
            mask = .true.
            cnt = 0
            mask(cc)  = .false.
            do while(dist < radius)
                dist = pixels_dist(self%centers(:,cc), self%centers,'min', mask, loc)
                mask(loc) = .false.
                cnt = cnt + 1
            enddo
            contact_scores(cc) = cnt - 1 !-1 because while loop counts one extra before exiting
        enddo
        if(present(cs_thresh)) then
          n_discard = count(contact_scores<cs_thresh) !unnecessary
        else
          n_discard = self%n_cc*PERCENT_DISCARD
        endif
        allocate(self%contact_scores(self%n_cc-n_discard), source = 0)
        call self%img_cc%get_imat(imat_cc)
        call self%img_bin%get_imat(imat_bin)
        cnt = 0
        ! Removing outliers from the binary image and the connected components image
        if(present(cs_thresh)) then
          do cc = 1, self%n_cc
            if(contact_scores(cc)<cs_thresh) then
              where(imat_cc == cc) imat_bin = 0
            else
              cnt = cnt + 1
              self%contact_scores(cnt) = contact_scores(cc)
            endif
          enddo
        else
          mask(1:self%n_cc) = .true. ! reset
          do cc = 1, n_discard
              label = minloc(contact_scores, mask)
              where(imat_cc == label(1))
                 imat_bin = 0
              endwhere
              mask(label(1)) = .false.
          enddo
          do cc = 1, self%n_cc ! fill in self%contact_scores
            if(mask(cc)) then
              cnt = cnt + 1
              self%contact_scores(cnt) = contact_scores(cc)
            endif
          enddo
        endif
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%find_ccs(self%img_cc)
        ! update number of connected components
        call self%img_cc%get_nccs(self%n_cc)
        call self%find_centers()
        write(logfhandle, *) '****outliers discarding, completed'
    end subroutine discard_outliers

    ! This subrouitne validates the identified atomic positions
    subroutine validate_atomic_positions(self)
        class(nanoparticle), intent(inout) :: self
        real,    allocatable :: x(:)
        real,    pointer     :: rmat_pc(:,:,:)
        integer, allocatable :: imat(:,:,:), imat_cc(:,:,:), imat_bin(:,:,:)
        integer, parameter   :: RANK_THRESH = 4
        integer :: n_cc, cnt
        integer :: rank, m(1)
        real    :: new_centers(3,2*self%n_cc)   !will pack it afterwards if it has too many elements
        real    :: new_contact_scores(2*self%n_cc)   !will pack it afterwards if it has too many elements
        real    :: pc
        call self%img%get_rmat_ptr(rmat_pc)     !now img contains the phase correlation
        call self%img_cc%get_imat(imat_cc)        !to pass to the subroutine split_atoms
        allocate(imat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), source = imat_cc)
        cnt = 0
        ! Remember to update the centers
        do n_cc =1, self%n_cc !for each cc check if the center corresponds with the local max of the phase corr
            pc = rmat_pc(nint(self%centers(1,n_cc)),nint(self%centers(2,n_cc)),nint(self%centers(3,n_cc)))
            !calculate the rank
            x = pack(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), imat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) == n_cc)
            call hpsort(x)
            m(:) = minloc(abs(x-pc))
            rank = size(x)-m(1)
            deallocate(x)
            if(rank > RANK_THRESH) then
                call split_atom(new_centers,cnt)
            else
                cnt = cnt + 1 !new number of centers deriving from splitting
                new_centers(:,cnt)      = self%centers(:,n_cc)
                new_contact_scores(cnt) = self%contact_scores(n_cc)
            endif
        enddo
        deallocate(self%centers)
        self%n_cc = cnt !update
        allocate(self%centers(3,cnt), source = 0.)
        ! update centers
        do n_cc =1, cnt
            self%centers(:,n_cc) = new_centers(:,n_cc)
        enddo
        deallocate(self%contact_scores)
        allocate(self%contact_scores(cnt), source = 0)
        ! update contact scores
        do n_cc =1, cnt
            self%contact_scores(n_cc) = new_contact_scores(n_cc)
        enddo
        call self%img_bin%get_imat(imat_bin)
        if(DEBUG) call self%img_bin%write_bimg(trim(self%fbody)//'BINbeforeValidation.mrc')
        ! update binary image
        where(imat_cc > 0)
            imat_bin = 1
        elsewhere
            imat_bin = 0
        endwhere
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%write_bimg(trim(self%fbody)//'BIN.mrc')
        call self%img_bin%find_ccs(self%img_cc)
        call self%img_cc%write_bimg(trim(self%fbody)//'CC.mrc')
        ! update number of ccs
        call self%update_self_ncc()
        ! update and write centers
        call self%find_centers()
        call self%write_centers()
    contains
        subroutine split_atom(new_centers,cnt)
            real,    intent(inout) :: new_centers(:,:)  !updated coordinates of the centers
            integer, intent(inout) :: cnt               !atom counter, to u pdate the center coords
            integer :: new_center1(3),new_center2(3),new_center3(3)
            integer :: i, j, k
            logical :: mask(self%ldim(1),self%ldim(2),self%ldim(3)) !false in the layer of connection of the atom to be split
            mask = .false.  !initialization
            ! Identify first new center
            new_center1 = maxloc(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), imat == n_cc)
            cnt = cnt + 1
            new_centers(:,cnt) = real(new_center1)
            new_contact_scores(cnt) = self%contact_scores(n_cc)
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    do k = 1, self%ldim(3)
                        if(((real(i-new_center1(1)))**2 + (real(j-new_center1(2)))**2 + &
                        &   (real(k-new_center1(3)))**2)*self%smpd  <=  (0.9*self%theoretical_radius)**2) then
                            if(imat(i,j,k) == n_cc) mask(i,j,k) = .true.
                        endif
                    enddo
                enddo
            enddo
            ! Second likely center.
                new_center2 = maxloc(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), &
                              & (imat == n_cc) .and. .not. mask)
                if(any(new_center2 > 0)) then !if anything was found
                    !Validate second center (check if it's 2 merged atoms, or one pointy one)
                    if(sum(real(new_center2-new_center1)**2.)*self%smpd <= 2.*self%theoretical_radius) then
                        ! Set the merged cc back to 0
                        where(imat_cc == n_cc .and. (.not.mask) ) imat_cc = 0
                        return
                    else
                      cnt = cnt + 1
                      new_centers(:,cnt) = real(new_center2)
                      new_contact_scores(cnt)   = self%contact_scores(n_cc) + 1
                      new_contact_scores(cnt-1) = new_contact_scores(cnt)
                      !In the case two merged atoms, build the second atom
                      do i = 1, self%ldim(1)
                          do j = 1, self%ldim(2)
                              do k = 1, self%ldim(3)
                                  if(((real(i-new_center2(1)))**2 + (real(j-new_center2(2)))**2 + (real(k-new_center2(3)))**2)*self%smpd < &
                                    & (0.9*self%theoretical_radius)**2) then
                                      if(imat(i,j,k) == n_cc)   mask(i,j,k) = .true.
                                  endif
                              enddo
                          enddo
                      enddo
                    endif
                  endif
                  ! Third likely center.
                      new_center3 = maxloc(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), &
                      & (imat == n_cc) .and. .not. mask)
                      if(any(new_center3 > 0)) then !if anything was found
                          !Validate third center
                          if(sum(real(new_center3-new_center1)**2.)*self%smpd <= 2.*self%theoretical_radius .or. &
                          &  sum(real(new_center3-new_center2)**2.)*self%smpd <= 2.*self%theoretical_radius ) then
                              ! Set the merged cc back to 0
                              where( imat_cc == n_cc .and. (.not.mask) ) imat_cc = 0
                              return
                          else
                            cnt = cnt + 1
                            new_centers(:,cnt) = real(new_center3)
                            new_contact_scores(cnt)   = self%contact_scores(n_cc) + 2
                            new_contact_scores(cnt-1) = new_contact_scores(cnt)
                            new_contact_scores(cnt-2) = new_contact_scores(cnt)
                            !In the case two merged atoms, build the second atom
                            do i = 1, self%ldim(1)
                                do j = 1, self%ldim(2)
                                    do k = 1, self%ldim(3)
                                        if(((real(i-new_center3(1)))**2 + (real(j-new_center3(2)))**2 + (real(k-new_center3(3)))**2)*self%smpd < (0.9*self%theoretical_radius)**2) then ! a little smaller to be sure
                                            if(imat(i,j,k) == n_cc)   mask(i,j,k) = .true.
                                        endif
                                    enddo
                                enddo
                            enddo
                          endif
                        endif
            ! Set the merged cc back to 0
            where(imat_cc == n_cc .and. (.not.mask) ) imat_cc = 0
            call self%img_cc%set_imat(imat_cc)
        end subroutine split_atom
    end subroutine validate_atomic_positions

    ! This subroutine calculates the histogram of the within-atoms
    ! distances distribution within the nanoparticle nano.
    ! To each atom the distance assigned is the min distance
    ! to the other atoms. There is a threshold (3.*self%theoretical_radius) for
    ! outliers discarding.
    ! If coords in input, then it considers just the atom-to-atom
    ! distance between the atoms with centers in coords.
    ! At the same time, it calculated statistics in a radial-dependent
    ! way, and save the result on a file.
    subroutine radial_dependent_stats(self,min_rad, max_rad, step)
        class(nanoparticle), intent(inout) :: self
        real,                intent(in)    :: min_rad, max_rad, step
        type(atoms)          :: radial_atoms_just,  radial_atoms_all
        type(image)          :: simulated_density
        logical, allocatable :: mask(:,:,:)
        real,    allocatable :: coords(:,:) !coordinates of the centers of the atoms according to radial distances
        real,    allocatable :: max_intensity(:), avg_intensity(:), stdev_intensity(:), ratios(:)
        real,    pointer     :: rmat(:,:,:)
        integer, allocatable :: imat_cc(:,:,:)
        real    :: m(3)    !mass center of c3 map
        real    :: d       !distance atoms from the center
        real    :: radius  !radius of the sphere to consider
        real    :: avg_int, max_int, stdev_int, cutoff
        integer :: cnt, cnt_just, cnt_all, filnum, io_stat
        integer :: nsteps, i, j, k, l, cc
        cutoff = 8.*self%smpd
        ! min_step and max_step is in A
        self%n_cc = size(self%centers, dim = 2)
        write(logfhandle, *) '****radial atom-to-atom distances estimation, init'
        nsteps = floor((max_rad-min_rad)/step)+1
        m        = self%nanopart_masscen()
        ! Report statistics in dedicated directory
        call simple_mkdir('./'//'/RadialDependentStat',errmsg="simple_nanoparticles :: radial_dependent_stats, simple_mkdir; ")
        call simple_chdir('./'//'/RadialDependentStat',errmsg="simple_nanoparticles :: radial_dependent_stats, simple_chdir; ")
        call fopen(filnum, file='RadialDependentStat.txt', iostat=io_stat)
        allocate(mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), source = .false.)
        call simulated_density%new(self%ldim,self%smpd)
        do l = 1, nsteps
            radius = min_rad+real(l-1)*step
            cnt_all  = 0
            cnt_just = 0
            do cc = 1, self%n_cc
                d = euclid(self%centers(:,cc), m)*self%smpd
                ! Count nb of atoms in the sphere of radius radius
               if(d<=radius) then
                   cnt_all = cnt_all+1
                   if(l == 1) then
                       cnt_just = cnt_all
                   else
                       if(d>radius-step .and. d<=radius) cnt_just = cnt_just+1
                   endif
               endif
            enddo
            call radial_atoms_just%new(cnt_just, dummy=.true.)
            call radial_atoms_all%new (cnt_all,  dummy=.true.)
            cnt_all  = 0
            cnt_just = 0
            ! Save coords
            do cc = 1, self%n_cc
                d = euclid(self%centers(:,cc), m)*self%smpd
               if(d<=radius) then
                   cnt_all = cnt_all+1
                   call radial_atoms_all%set_name(cnt_all,self%atom_name)
                   call radial_atoms_all%set_element(cnt_all,self%element)
                   call radial_atoms_all%set_coord(cnt_all,(self%centers(:,cc)-1.)*self%smpd)
                   if(l == 1) then
                       cnt_just= cnt_just+1
                       call radial_atoms_just%set_name(cnt_just,self%atom_name)
                       call radial_atoms_just%set_element(cnt_just,self%element)
                       call radial_atoms_just%set_coord(cnt_just,(self%centers(:,cc)-1.)*self%smpd)
                   endif
                   if(d>(radius-step) .and. l > 1) then
                       cnt_just = cnt_just+1
                       call radial_atoms_just%set_name(cnt_just,self%atom_name)
                       call radial_atoms_just%set_element(cnt_just,self%element)
                       call radial_atoms_just%set_coord(cnt_just,(self%centers(:,cc)-1.)*self%smpd)
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
            do cc = 1, size(self%centers, dim = 2)
                d = euclid(self%centers(:,cc), m)*self%smpd
                if(d<=radius) then
                    if(l == 1) then
                        cnt = cnt + 1
                        coords(:3,cnt) = self%centers(:,cc)
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
                        coords(:3,cnt) = self%centers(:,cc)
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
            call self%distances_distribution(coords)
            deallocate(coords)
            deallocate(avg_intensity)
            deallocate(max_intensity)
            deallocate(stdev_intensity)
            deallocate(ratios)
            call radial_atoms_all%kill
            call radial_atoms_just%kill
        enddo
        call fclose(filnum)
        deallocate(mask, imat_cc)
        call simulated_density%kill
        call self%distances_distribution() ! across the whole nano
        call self%atom_intensity_stats(.false.)
        call self%kill
        write(logfhandle, *) '****radial atom-to-atom distances estimation, completed'
   end subroutine radial_dependent_stats

   subroutine distances_distribution(self,coords,volume)
       class(nanoparticle), intent(inout) :: self
       real,    optional,   intent(in)    :: coords(:,:)
       integer, optional,   intent(in)    :: volume
       real, allocatable :: dist(:)
       real    :: stdev, med
       integer :: i, j, n_discard, filnum, io_stat
       logical :: mask(self%n_cc)
       ! Initialisations
       mask       = .true.
       stdev      = 0.
       n_discard  = 0
       if(present(coords)) then
           call fopen(filnum, file='DistancesDistr.txt', iostat=io_stat)
           allocate(dist(size(coords,dim=2)), source = 0.)
           do i = 1, size(coords,dim=2)
               dist(i) =  pixels_dist(coords(:,i), self%centers(:,:), 'min', mask=mask) !I have to use all the atoms when
               mask(:) = .true. ! restore
               !Discard outliers
               if(dist(i)*self%smpd > 3.*self%theoretical_radius ) then      !maximum interatomic distance
                   dist(i) = 0.
                   n_discard = n_discard + 1
               else if(dist(i)*self%smpd < 1.5*self%theoretical_radius ) then ! minimum interatomic distance
                   dist(i) = 0.
                   n_discard = n_discard + 1
               endif
           enddo
           self%avg_dist_atoms = sum(dist)/real(size(coords,dim=2)-n_discard)
           do i = 1, size(coords,dim=2)
               if(dist(i)*self%smpd <=3.*self%theoretical_radius) stdev = stdev + (dist(i)-self%avg_dist_atoms)**2
           enddo
           stdev = sqrt(stdev/real(size(coords,dim=2)-1-n_discard))
           med = median(dist)
           if(present(volume)) then
               write(unit = filnum, fmt = '(a,a,a,f6.3,a)') 'Average dist atoms vol ', trim(int2str(volume)),':', self%avg_dist_atoms*self%smpd, ' A'
               write(unit = filnum, fmt = '(a,a,a,f6.3,a)') 'StDev   dist atoms vol ', trim(int2str(volume)),':', stdev*self%smpd, ' A'
               write(unit = filnum, fmt = '(a,a,a,f6.3,a)') 'Median  dist atoms vol ', trim(int2str(volume)),':', med*self%smpd, ' A'
           else
               write(unit = filnum, fmt = '(a,f6.3,a)') 'Average dist atoms: ', self%avg_dist_atoms*self%smpd, ' A'
               write(unit = filnum, fmt = '(a,f6.3,a)') 'StDev   dist atoms: ', stdev*self%smpd, ' A'
               write(unit = filnum, fmt = '(a,f6.3,a)') 'Median  dist atoms: ', med*self%smpd, ' A'
           endif
           deallocate(dist)
       else
           call fopen(filnum, file='DistancesDistr.txt', iostat=io_stat)
           allocate(self%dists(size(self%centers, dim = 2)), source = 0.)
           do i = 1, size(self%centers, dim = 2)
               self%dists(i) =  pixels_dist(self%centers(:,i), self%centers(:,:), 'min', mask=mask) !Use all the atoms
               mask(:) = .true. ! restore
               !Discard outliers
               if(self%dists(i)*self%smpd > 3.*self%theoretical_radius ) then
                   self%dists(i) = 0.
                   n_discard = n_discard + 1
               else if(self%dists(i)*self%smpd < 1.5*self%theoretical_radius ) then
                   self%dists(i) = 0.
                   n_discard = n_discard + 1
               endif
           enddo
           self%avg_dist_atoms = sum(self%dists)/real(size(self%centers, dim = 2)-n_discard)
           do i = 1, size(self%centers, dim = 2)
               if(self%dists(i)*self%smpd <=3.*self%theoretical_radius) stdev = stdev + (self%dists(i)-self%avg_dist_atoms)**2
           enddo
           stdev = sqrt(stdev/real(size(self%centers, dim = 2)-1-n_discard))
           med = median(self%dists)
           write(unit = filnum, fmt = '(a,f6.3,a)') 'Average dist atoms: ', self%avg_dist_atoms*self%smpd, ' A'
           write(unit = filnum, fmt = '(a,f6.3,a)') 'StDev   dist atoms: ', stdev*self%smpd, ' A'
           write(unit = filnum, fmt = '(a,f6.3,a)') 'Median  dist atoms: ', med*self%smpd, ' A'
       endif
       call fclose(filnum)
   end subroutine distances_distribution

   subroutine aspect_ratios_estimation(self, print_ar)
       class(nanoparticle), intent(inout) :: self
       logical, optional,   intent(in)    :: print_ar !print longest/shortest dim and ratio
       real,    allocatable :: longest_dist(:)
       integer, allocatable :: imat_cc(:,:,:)
       integer              :: label, filnum, io_stat
       real                 :: avg_diameter, median_diameter, min_diameter, max_diameter, stdev_diameter
       call self%img_cc%get_imat(imat_cc)
       self%n_cc = maxval(imat_cc) ! update number of connected components
       allocate(self%ratios (self%n_cc),             source = 0.)
       allocate(longest_dist(self%n_cc),             source = 0.)
       allocate(self%loc_longest_dist(3,self%n_cc),  source = 0 )
       call self%find_centers() !TO KEEP
       do label = 1, self%n_cc
           call self%calc_aspect_ratio(label, self%ratios(label),lld=.true., ld=longest_dist(label), print_ar=print_ar)
       enddo
       longest_dist = 2.*longest_dist ! radius --> diameter
       min_diameter = minval(longest_dist(1:self%n_cc))
       max_diameter = maxval(longest_dist(1:self%n_cc))
       median_diameter = median(longest_dist(1:self%n_cc))
       avg_diameter    = sum(longest_dist(1:self%n_cc))/real(self%n_cc)
       stdev_diameter = 0.
       do label = 1, self%n_cc
           stdev_diameter = stdev_diameter + (avg_diameter-longest_dist(label))**2
       enddo
       stdev_diameter = sqrt(stdev_diameter/real(self%n_cc-1))
       if(DEBUG) then
           write(logfhandle,*) 'minimum  value diameter ', min_diameter, 'A'
           write(logfhandle,*) 'maximum  value diameter ', max_diameter, 'A'
           write(logfhandle,*) 'median   value diameter ', median_diameter, 'A'
           write(logfhandle,*) 'average  value diameter ', avg_diameter, 'A'
           write(logfhandle,*) 'stdev    value diameter ', stdev_diameter, 'A'
       endif
       call fopen(filnum, file='AspectRatio.csv', iostat=io_stat)
       write (filnum,*) 'AR'
       do label = 1, size(self%ratios)
         write (filnum,'(A)', advance='yes') trim(real2str(self%ratios(label)))
       end do
       call fclose(filnum)
       deallocate(longest_dist)
       write(logfhandle,*)'**aspect ratio calculations completed'
   end subroutine aspect_ratios_estimation

   ! This subroutine takes in input a connected component (cc) image
   ! and the label of one of its ccs and calculates its aspect ratio, which
   ! is defined as the ratio of the width and the height.
   ! The idea behind this is that the center of the cc is calculated,
   ! than everything is deleted except the borders of the cc. Finally,
   ! in order to calculate the width and the height, the min/max
   ! distances between the center and the borders are calculated. The
   ! aspect ratio is the ratio of those 2 distances.
   subroutine calc_aspect_ratio(self,label,ratio,lld, ld,print_ar)
       class(nanoparticle), intent(inout) :: self
       integer,             intent(in)    :: label
       real,                intent(out)   :: ratio
       logical,             intent(in)    :: lld ! fill in self%loc_longest_dim
       real   , optional,   intent(out)   :: ld  ! longest dist
       logical, optional,   intent(in)    :: print_ar
       integer, allocatable :: pos(:,:)
       integer, allocatable :: imat_cc(:,:,:)
       logical, allocatable :: border(:,:,:)
       logical, allocatable :: mask_dist(:) !for min and max dist calculation
       integer :: location(1) !location of the farest vxls of the atom from its center
       integer :: i
       real    :: shortest_dist, longest_dist
       call self%img_cc%get_imat(imat_cc)
       call self%img_cc%border_mask(border, label, .true.) !use 4neigh instead of 8neigh
       where(border .eqv. .true.)
           imat_cc = 1
       elsewhere
           imat_cc = 0
       endwhere
       call get_pixel_pos(imat_cc,pos)   !pxls positions of the shell
       if(allocated(mask_dist)) deallocate(mask_dist)
       allocate(mask_dist(size(pos, dim = 2)), source = .true. )
       shortest_dist = pixels_dist(self%centers(:,label), real(pos),'min', mask_dist, location)
       if(size(pos,2) == 1) then !if the connected component has size 1 (just 1 vxl)
           shortest_dist = 0.
           longest_dist  = shortest_dist
           ratio = 1.
           if(lld) self%loc_longest_dist(:3, label) = self%centers(:,label)
           if(present(print_ar) .and. (print_ar .eqv. .true.)) then
                write(logfhandle,*) 'ATOM #          ', label
                write(logfhandle,*) 'shortest dist = ', shortest_dist
                write(logfhandle,*) 'longest  dist = ', longest_dist
                write(logfhandle,*) 'RATIO         = ', ratio
           endif
           return
       else
           longest_dist  = pixels_dist(self%centers(:,label), real(pos),'max', mask_dist, location)
           if(lld) self%loc_longest_dist(:3, label) =  pos(:3,location(1))
       endif
       if(abs(longest_dist) > TINY .and. size(pos,2) > 1) then
           ratio = shortest_dist/longest_dist
       else
            ratio = 0.
            if(DEBUG) write(logfhandle,*) 'cc ', label, 'LONGEST DIST = 0'
       endif
       longest_dist  = longest_dist*self%smpd  !in A
       shortest_dist = shortest_dist*self%smpd
      if(present(print_ar) .and. (print_ar .eqv. .true.)) then
           write(logfhandle,*) 'ATOM #          ', label
           write(logfhandle,*) 'shortest dist = ', shortest_dist
           write(logfhandle,*) 'longest  dist = ', longest_dist
           write(logfhandle,*) 'RATIO         = ', ratio
      endif
      if(present(ld)) ld=longest_dist
      deallocate(imat_cc, border, pos, mask_dist)
   end subroutine calc_aspect_ratio

    ! This subroutine calculates some statistics (min,max,avg,stdev)
    ! in the intensity gray level value of the nanoparticle
    ! map in each atom. It is likely that these statistics
    ! are going to be able to distinguish between the different
    ! atom compositions in heterogeneous nanoparticles.
    subroutine atom_intensity_stats(self, print_file)
        class(nanoparticle), intent(inout) :: self
        logical,             intent(in)    :: print_file ! whether to print on a file all the calc stats
        real,    pointer     :: rmat(:,:,:)
        logical, allocatable :: mask(:,:,:)
        integer, allocatable :: imat_cc(:,:,:), sz(:)
        integer :: i, j, k,n_atom, filnum, io_stat
        real    :: max_intensity(self%n_cc), avg_intensity(self%n_cc)
        real    :: stdev_intensity(self%n_cc), radii(self%n_cc), int_intensity(self%n_cc) ! integrated intensity
        real    :: avg_int, stdev_int, max_int
        real    :: m(3) ! center of mass of the nanoparticle
        write(logfhandle,*)'**atoms intensity statistics calculations init'
        call self%img%get_rmat_ptr(rmat)
        call self%img_cc%get_imat(imat_cc)
        allocate(mask(self%ldim(1),self%ldim(2),self%ldim(3)), source = .false.)
        self%n_cc = maxval(imat_cc) ! update number of connected components
        m = self%nanopart_masscen()
        ! initialise
        max_intensity(:)   = 0.
        avg_intensity(:)   = 0.
        int_intensity(:)   = 0.
        stdev_intensity(:) = 0.
        ! Write on a file
        if(print_file) call fopen(filnum, file ='IntensityStats.txt', iostat=io_stat)
        do n_atom = 1, self%n_cc
            where(imat_cc == n_atom) mask = .true.
            max_intensity(n_atom) = maxval(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
            int_intensity(n_atom) = sum(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
            avg_intensity(n_atom) = int_intensity(n_atom)/real(count(mask))
            radii(n_atom)         = euclid(self%centers(:,n_atom), m)*self%smpd
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
             if(print_file) write(unit = filnum, fmt = '(a,i3,a,f9.5,a,f9.5,a,f9.5)') &
             & 'ATOM # ', n_atom, 'maxval ', max_intensity(n_atom),'   avg ', avg_intensity(n_atom), '   stdev ', stdev_intensity(n_atom), ' integrated ', int_intensity(n_atom)
            mask = .false. !Reset
        enddo
        avg_int   = sum(avg_intensity)/real(self%n_cc)
        max_int   = maxval(max_intensity)
        stdev_int = 0.
        do n_atom= 1, self%n_cc
            stdev_int = stdev_int + (avg_intensity(n_atom)-avg_int)**2
        enddo
        stdev_int = sqrt(stdev_int/real(self%n_cc-1))
        if(print_file) write(unit = filnum, fmt = '(a,f9.5,a,f9.5,a,f9.5)') 'maxval_general ', max_int, '   avg_general ', avg_int, '   stdev_general ', stdev_int
        if(print_file) call fclose(filnum)
        ! calculate sz of connected components, just for output on a file
        sz = self%img_cc%size_ccs()
        ! print one by one the stats anyway
        call fopen(filnum, file='../'//'RadialPos.csv', iostat=io_stat)
        write (filnum,*) 'radius'
        do n_atom = 1, self%n_cc
            write (filnum,'(A)', advance='yes') trim(real2str(radii(n_atom)))
        end do
        call fclose(filnum)
        call fopen(filnum, file='../'//'MaxIntensity.csv', iostat=io_stat)
        write (filnum,*) 'maxint'
        do n_atom = 1, self%n_cc
            write (filnum,'(A)', advance='yes') trim(real2str(max_intensity(n_atom)))
        end do
        call fclose(filnum)
        call fopen(filnum, file='../'//'AvgIntensity.csv', iostat=io_stat)
        write (filnum,*) 'avgint'
        do n_atom = 1, self%n_cc
            write (filnum,'(A)', advance='yes') trim(real2str(avg_intensity(n_atom)))
        end do
        call fclose(filnum)
        call fopen(filnum, file='../'//'IntIntensity.csv', iostat=io_stat)
        write (filnum,*) 'intint'
        do n_atom = 1, self%n_cc
            write (filnum,'(A)', advance='yes') trim(real2str(int_intensity(n_atom)))
        end do
        call fclose(filnum)
        call fopen(filnum, file='../'//'Size.csv', iostat=io_stat)
        write (filnum,*) 'sz'
        do n_atom = 1, self%n_cc
            write (filnum,'(A)', advance='yes') trim(int2str(sz(n_atom)))
        end do
        call fclose(filnum)
        deallocate(sz)
        write(logfhandle,*)'**atoms intensity statistics calculations completed'
    end subroutine atom_intensity_stats

    ! This subroutine clusters the atoms with respect to the maximum intensity
    ! or the integrated density (according to the values contained in feature)
    ! using kmean algorithm for 2 classes. The initial guess fo the centers
    ! is intentionally biased. It supposes there are two distinguished classes
    ! with different avgs (proved with simulated data).
    subroutine cluster_atom_intensity(self, feature)
        use gnufor2
        use simple_nanoML, only : nanoML
        class(nanoparticle), intent(inout) :: self
        real,                intent(inout) :: feature(:)
        integer, parameter   :: MAX_IT = 50 ! maximum number of iterations for
        integer, allocatable :: imat_cc(:,:,:)
        real, pointer        :: rmat1(:,:,:), rmat2(:,:,:)
        type(image)          :: class1, class2
        type(nanoML)         :: emfit
        real    :: centers_kmeans(2) !output of k-means
        real    :: avgs(2),  vars(2), gammas(self%n_cc,2) !output of ML
        integer :: i, cnt1, cnt2, filnum, io_stat
        feature = feature/maxval(feature)*10. !normalise with maxval 10
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
        avgs = emfit%get_avgs()
        vars = emfit%get_vars()
        gammas  = emfit%get_gammas()
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
    subroutine cluster_atom_maxint(self)
      class(nanoparticle), intent(inout) :: self
      real    :: max_intensity(self%n_cc)
      integer :: i, io_stat, filnum
      call fopen(filnum, file='MaxIntensity.csv', iostat=io_stat)
      if( io_stat .ne. 0 ) then
        THROW_HARD('Unable to read file MaxIntensity.csv Did you run radial_dependent_stats?; cluster_atom_maxint')
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
    subroutine cluster_atom_intint(self)
      class(nanoparticle), intent(inout) :: self
      real    :: int_intensity(self%n_cc)
      integer :: i, io_stat, filnum
      call fopen(filnum, file='IntIntensity.csv', iostat=io_stat)
      if( io_stat .ne. 0 ) then
        THROW_HARD('Unable to read file IntIntensity.csv Did you run radial_dependent_stats?; cluster_atom_intint')
      endif
      read(filnum,*) ! first line is variable name
      do i = 1, self%n_cc
        read(filnum,*) int_intensity(i)
      enddo
      call fclose(filnum)
      call self%cluster_atom_intensity(int_intensity)
    end subroutine cluster_atom_intint

    ! Polarization search via angle variance. The considered angle is
    ! the angle between the vector [0,0,1] and the direction of the
    ! longest dim of the atom.
    ! This is done for each of the atoms and then there is a search
    ! of the similarity between the calculated angles.
    subroutine search_polarization(self)
        class(nanoparticle), intent(inout) :: self
        real    :: loc_ld_real(3,self%n_cc)
        real    :: vec_fixed(3)     !vector indentifying z direction
        real    :: theta, mod_1, mod_2, dot_prod
        integer :: k, i, filnum, io_stat
        integer, allocatable :: sz(:)
        if(allocated(self%ang_var)) deallocate(self%ang_var)
           allocate (self%ang_var(self%n_cc), source = 0.)
        sz = self%img_cc%size_ccs()
        !bring vector back to center
        do i = 1, self%n_cc
            loc_ld_real(:3,i) = real(self%loc_longest_dist(:3,i))- self%centers(:3,i)
        enddo
        !consider fixed vector [0,0,1] (z direction)
        vec_fixed(1) = 0.
        vec_fixed(2) = 0.
        vec_fixed(3) = 1.
        write(logfhandle,*)'>>>>>>>>>>>>>>>> calculating angles wrt the vector [0,0,1]'
        do k = 1, self%n_cc
            if(sz(k) > 2) then
                self%ang_var(k) = ang3D_vecs(vec_fixed(:),loc_ld_real(:,k))
            else
                self%ang_var(k) = 0. ! If the cc is too small it doesn't make sense
            endif
            if(DEBUG) write(logfhandle,*) 'ATOM ', k, 'angle between direction longest dim and vec [0,0,1] ', self%ang_var(k)
        enddo
        call fopen(filnum, file='AnglesLongestDims.csv', iostat=io_stat)
        write (filnum,*) 'ang'
        do k = 1, self%n_cc
            write (filnum,'(A)', advance='yes') trim(real2str(self%ang_var(k)))
        end do
        call fclose(filnum)
    end subroutine search_polarization

    ! Source https://www.mathworks.com/help/stats/hierarchical-clustering.html#bq_679x-10
    subroutine hierarc_clust(self,vec,thresh,labels,centroids, populations)
      class(nanoparticle), intent(inout) :: self
      real,    intent(in)  :: vec(:)       ! input vector
      real,    intent(in)  :: thresh       ! threshold for class merging
      integer, intent(out) :: labels(:)    ! labels of the elements in vec
      real,    allocatable, intent(out) :: centroids(:) ! centroids of the classes
      integer, allocatable, intent(out) :: populations(:) ! number of elements belonging to each class
      real,    allocatable :: mat(:,:)        ! matrix that contains the distances
      logical, allocatable :: mask(:), outliers(:)
      real    :: d
      integer :: N, i, j, index(1), loc1(1), loc2(1), cnt, ncls, filnum, io_stat
      N = size(vec) ! number of data points
      ! 1) calc all the couples of distances, using euclid dist
      allocate(mat(N,N), source = 0.)
      do i = 1, N-1
        do j = i+1, N
          mat(i,j) = abs(vec(i)-vec(j))
          mat(j,i) = mat(i,j)
        enddo
      enddo
      ! 2) Generate binary clusters
      allocate(mask(N),source = .true. )
      allocate(outliers(N), source = .false.)
      ncls = 0
      do i = 1, N ! fin the class for vec(i)
        if(mask(i)) then ! if it's not already been clustered
          mask(i) = .false.
          ! find the index of the couple
          d = minval(mat(i,:), mask)
          index(:) = minloc(mat(i,:), mask)
          ncls = ncls + 1
          ! assign lavels
          labels(i) = ncls
          if(d <= thresh) then ! if it's not an outlier (it has a couple)
            labels(index(1)) = labels(i)
            mask(index(1)) = .false. ! index(1) has already been clustered
          else
            outliers(i) = .true.
          endif
        endif
      enddo
      ! 3) Calculate centroids
      allocate(centroids(ncls), source = 0.)
      mask = .true. ! reset
      do i = 1, ncls
        ! find the first member of the class
        loc1(:) = minloc(abs(labels-i))
        if(.not. outliers(loc1(1))) then
          mask(loc1(1)) = .false.
          ! find the second member of the class
          loc2(:) = minloc(abs(labels-i), mask)
          mask(loc2(1)) = .false.
          centroids(i) = (vec(loc1(1)) + vec(loc2(1)))/2.
        else ! the class has just one element
          loc1(:) = minloc(abs(labels-i))
          centroids(i) = vec(loc1(1))
          mask(loc1(1)) = .false.
        endif
      enddo
      mask  = .true. ! reset
      ! 4) merge classes
      do i = 1, ncls-1
        do j = i+1, ncls
          if(abs(centroids(i)-centroids(j)) <= thresh) then ! merge classes
            ! change label to class j
            where(labels == j) labels = i
          endif
        enddo
      enddo
      ! 5) Reoder labels
      cnt = 0
      do i = 1, ncls
         if(any(labels== i)) then !there is a class labelled i
            cnt = cnt + 1                  !increasin order cc
            where(labels == i) labels = cnt
          endif
      enddo
      ! 6) recalculate centroids
      deallocate(centroids)
      ncls = maxval(labels) ! the nb of classes is maxval(labels)
      allocate(centroids(ncls),   source = 0.)
      allocate(populations(ncls), source = 0 )
      mask = .true. ! reset
      do i = 1, ncls
        populations(i) = count(labels == i) ! counts the nb of elements in the class
        ! find the all the cnt member of the class and update the centroids
        do j = 1, populations(i)
          loc1(:) = minloc(abs(labels-i), mask)
          mask(loc1(1)) = .false. ! do not consider this member of the class anymore
          centroids(i) = centroids(i)+ vec(loc1(1))
        enddo
        centroids(i) = centroids(i)/real(populations(i))
      enddo
    end subroutine hierarc_clust

    ! Cluster wrt to the angle between the direction of the longest dimension and a fixed vector
    subroutine cluster_ang(self, thresh)
        class(nanoparticle), intent(inout) :: self
        real,    intent(in)  :: thresh ! threshold for class definition, user inputted
        real,    allocatable :: stdev_within(:), centroids(:)
        integer, allocatable :: labels(:), populations(:)
        integer              :: i, j, ncls, dim, filnum, io_stat, cnt
        real                 :: avg, stdev
        ! Preparing for clustering
        ! Aspect ratios and loc_longest dist estimation
        call self%aspect_ratios_estimation(print_ar=.false.) ! it's needed to identify the dir of longest dim
        call self%search_polarization()
        dim = size(self%ang_var)
        allocate(labels(dim), source = 0)
        call self%hierarc_clust(self%ang_var,thresh,labels,centroids,populations)
        ! Stats calculations
        ncls = maxval(labels)
        allocate(stdev_within(ncls), source = 0.)
        !stdev within the same class
        do i = 1, ncls
            do j = 1, dim
                if(labels(j) == i) stdev_within(i) = stdev_within(i) + (self%ang_var(j)-centroids(i))**2
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
        call fopen(filnum, file='ClusterAngThresh'//trim(real2str(thresh))//'.txt', iostat=io_stat)
        write(unit = filnum,fmt ='(a,i2,a,f6.2)') 'NR OF IDENTIFIED CLUSTERS:', ncls, ' SELECTED THRESHOLD: ',  thresh
        write(unit = filnum,fmt ='(a)') 'CLASSIFICATION '
        do i = 1, dim
          write(unit = filnum,fmt ='(a,i3,a,f6.2,a,i3)') 'Atom #: ', i, '; data (deg): ', self%ang_var(i), '; class: ', labels(i)
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
          write(filnum,*) self%ang_var(i)
        enddo
        call fclose(filnum)
        if(allocated(stdev_within)) deallocate(stdev_within)
        deallocate(centroids, labels, populations)
      end subroutine cluster_ang

  ! Cluster the atoms wrt to the aspect ratio
  subroutine cluster_ar(self, thresh)
      class(nanoparticle), intent(inout) :: self
      real,    intent(in)  :: thresh ! threshold for class definition, user inputted
      real,    allocatable :: centroids(:)
      integer, allocatable :: labels(:), populations(:)
      real,    allocatable :: stdev_within(:)
      integer              :: i, j, ncls, dim, filnum, io_stat, cnt
      real                 :: avg, stdev
      if(thresh > 1. .or. thresh < 0.) THROW_HARD('Invalid input threshold! AR is in [0,1]; cluster_ar')
      ! Preparing for clustering
      ! Aspect ratios calculations
      call self%aspect_ratios_estimation(print_ar=.false.) ! it's needed to identify the dir of longest dim
      dim = size(self%ratios)
      allocate(labels(dim), source = 0)
      ! classify
      call self%hierarc_clust(self%ratios,thresh,labels,centroids,populations)
      ! Stats calculations
      ncls = maxval(labels)
      allocate(stdev_within(ncls), source = 0.)
      ! stdev within the same class
      do i = 1, ncls
          do j = 1, dim
              if(labels(j) == i) stdev_within(i) = stdev_within(i) + (self%ratios(j)-centroids(i))**2
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
      call fopen(filnum, file='ClusterARThresh'//trim(real2str(thresh))//'.txt', iostat=io_stat)
      write(unit = filnum,fmt ='(a,i2,a,f6.2)') 'NR OF IDENTIFIED CLUSTERS:', ncls, ' SELECTED THRESHOLD: ',  thresh
      write(unit = filnum,fmt ='(a)') 'CLASSIFICATION '
      do i = 1, dim
        write(unit = filnum,fmt ='(a,i3,a,f6.2,a,i3,a,f6.2)') 'Atom #: ', i, '; data (adimensional): ', self%ratios(i), '; class: ', labels(i)
      enddo
      write(unit = filnum,fmt ='(a)') 'CLASS STATISTICS '
      do i = 1, ncls
        write(unit = filnum,fmt ='(a,i3,a,i3,a,f6.2,a,f6.2)') 'class: ', i, '; cardinality: ', populations(i), '; centroid: ', centroids(i), '; stdev within the class: ', stdev_within(i)
      enddo
      write(unit = filnum,fmt ='(a,f6.2,a,f6.2)') 'AVG among the classes: ', avg, '; STDEV among the classes: ', stdev
      call fclose(filnum)
      ! Generate file that contains the calculater AR, that can be read afterwards
      call fopen(filnum, file='Ar.csv', iostat=io_stat)
      write(filnum,*) 'ar'
      do i  = 1, self%n_cc
        write(filnum,*) self%ratios(i)
      enddo
      call fclose(filnum)
      if(allocated(stdev_within)) deallocate(stdev_within)
      deallocate(centroids, labels, populations)
    end subroutine cluster_ar

    ! Cluster the atoms wrt to the interatomic distances
    subroutine cluster_interdist(self, thresh)
        class(nanoparticle), intent(inout) :: self
        real,    intent(in)  :: thresh ! threshold for class definition, user inputted
        real,    allocatable :: centroids(:)
        integer, allocatable :: labels(:), populations(:)
        real,    allocatable :: stdev_within(:)
        integer              :: i, j, ncls, dim, filnum, io_stat, cnt
        real                 :: avg, stdev
        ! Preparing for clustering
        ! need to recalculate self%dists
        call self%distances_distribution()
        ! Aspect ratios calculations
        dim = size(self%dists)
        ! pass from pixels to A
        self%dists = self%dists*self%smpd
        allocate(labels(dim), source = 0)
        ! classify
        call self%hierarc_clust(self%dists,thresh,labels,centroids,populations)
        ! Stats calculations
        ncls = maxval(labels)
        allocate(stdev_within(ncls), source = 0.)
        ! stdev within the same class
        do i = 1, ncls
            do j = 1, dim
                if(labels(j) == i) stdev_within(i) = stdev_within(i) + (self%dists(j)-centroids(i))**2
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
        call fopen(filnum, file='ClusterInterDistThresh'//trim(real2str(thresh))//'.txt', iostat=io_stat)
        write(unit = filnum,fmt ='(a,i2,a,f6.2)') 'NR OF IDENTIFIED CLUSTERS:', ncls, ' SELECTED THRESHOLD: ',  thresh
        write(unit = filnum,fmt ='(a)') 'CLASSIFICATION '
        do i = 1, dim
          write(unit = filnum,fmt ='(a,i3,a,f6.2,a,i3)') 'Atom #: ', i, '; data (A): ', self%dists(i), '; class: ', labels(i)
        enddo
        write(unit = filnum,fmt ='(a)') 'CLASS STATISTICS '
        do i = 1, ncls
          write(unit = filnum,fmt ='(a,i3,a,i3,a,f6.2,a,f6.2,a)') 'class: ', i, '; cardinality: ', populations(i), '; centroid: ', centroids(i), ' A; stdev within the class: ', stdev_within(i), ' A'
        enddo
        write(unit = filnum,fmt ='(a,f6.2,a,f6.2,a)') 'AVG among the classes: ', avg, ' A; STDEV among the classes: ', stdev, ' A'
        call fclose(filnum)
        call fopen(filnum, file='Dist.csv', iostat=io_stat)
        write(filnum,*) 'dist'
        do i  = 1, self%n_cc
          write(filnum,*) self%dists(i)
        enddo
        call fclose(filnum)
        if(allocated(stdev_within)) deallocate(stdev_within)
        deallocate(centroids, labels, populations)
      end subroutine cluster_interdist



    subroutine make_soft_mask(self) !change the name
        class(nanoparticle), intent(inout) :: self
        type(binimage) :: simulated_density
        type(image)   :: img_cos
        type(atoms)   :: atomic_pos
        call simulated_density%new_bimg(self%ldim, self%smpd)
        call atomic_pos%new(trim(self%fbody)//'_atom_centers.pdb')
        call atomic_pos%convolve(simulated_density, cutoff=8.*self%smpd)
        call simulated_density%grow_bins(nint(0.5*self%theoretical_radius/self%smpd)+1)
        call simulated_density%cos_edge(img_cos,SOFT_EDGE)
        call img_cos%write(trim(self%fbody)//'SoftMask.mrc')
        !kill
        call img_cos%kill
        call simulated_density%kill_bimg
        call atomic_pos%kill
    end subroutine make_soft_mask

    subroutine atoms_rmsd(nano1,nano2,r)
        use simple_atoms, only : atoms
        use gnufor2
        class(nanoparticle), intent(inout) :: nano1, nano2 !nanoparticles to compare
        real, optional,      intent(out)   :: r   !rmsd calculated
        logical, allocatable :: mask(:)
        real,    allocatable :: dist(:), dist_sq(:), dist_no_zero(:), dist_close(:)
        real,    parameter   :: CLOSE_THRESH = 1. ! 1 Amstrong
        integer :: location(1), i, j
        integer :: N_min !min{nb atoms in nano1, nb atoms in nano2}
        integer :: N_max !max{nb atoms in nano1, nb atoms in nano2}
        integer :: cnt,cnt2,cnt3,filnum,io_stat
        real    :: sum, rmsd, coord(3)
        real    :: avg, stdev, m(3), tmp_max, d !for statistics calculation
        type(atoms) :: centers_coupled1, centers_coupled2 !visualization purposes
        type(atoms) :: centers_close1, centers_close2
        type(atoms) :: couples1
        call fopen(filnum, file='CompareAtomicModels.txt', iostat=io_stat)
        write(unit = filnum, fmt = '(a)') '>>>>>>>>>   COMPARE NANO   >>>>>>>>>'
        write(unit = filnum, fmt = '(a)') ''
        write(unit = filnum, fmt = '(a)') 'Comparing atomic models of particles'
        write(unit = filnum, fmt = '(a,a)') trim(nano1%fbody), ' ---> vol1'
        write(unit = filnum, fmt = '(a)') 'and'
        write(unit = filnum, fmt = '(a,a)') trim(nano2%fbody), ' ---> vol2'
        write(unit = filnum, fmt = '(a)')  '>>>>>>>>>VOLUME COMPARISION>>>>>>>>'
        ! If they don't have the same nb of atoms
        if(size(nano1%centers, dim = 2) <= size(nano2%centers, dim = 2)) then
            N_min = size(nano1%centers, dim = 2)
            N_max = size(nano2%centers, dim = 2)
            write(unit = filnum, fmt = '(a,i3)') 'NB atoms in vol1   ', N_min
            write(unit = filnum, fmt = '(a,i3)') 'NB atoms in vol2   ', N_max
            write(unit = filnum, fmt = '(a,i3,a)') '--->', abs(N_max-N_min), ' atoms do NOT correspond'
            call centers_coupled1%new(N_max, dummy=.true.)
            call centers_coupled2%new(N_max, dummy=.true.)
            call centers_close1%new  (N_max, dummy=.true.)
            call centers_close2%new  (N_max, dummy=.true.)
            call couples1%new        (N_max, dummy=.true.)
            allocate(dist(N_max), dist_sq(N_max), source = 0.)
            allocate(dist_close(N_max), source = 0.) ! there are going to be unused entry of the vector
            allocate(mask(N_min), source = .true.)
            cnt  = 0
            cnt2 = 0
            cnt3 = 0
            do i = 1, N_max !compare based on centers2
                dist(i) = pixels_dist(nano2%centers(:,i),nano1%centers(:,:),'min',mask,location, keep_zero=.true.)
                if(dist(i)*nano2%smpd > 2.*nano2%theoretical_radius) then
                    dist(i) = 0. !it means there is no correspondent atom in the other nano
                    cnt = cnt + 1  !to discard them in the rmsd calculation
                    call couples1%set_coord(i,(nano2%centers(:,location(1))-1.)*nano2%smpd)
                    ! remove the atoms from the pdb file
                    call centers_coupled1%set_occupancy(i,0.)
                    call centers_coupled2%set_occupancy(i,0.)
                    call centers_close1%set_occupancy(i,0.)
                    call centers_close2%set_occupancy(i,0.)
                    call couples1%set_occupancy(i,0.)
                elseif(dist(i)*nano2%smpd < CLOSE_THRESH) then
                    cnt3 = cnt3 + 1
                    dist_close(i) = dist(i)**2
                    call centers_close2%set_coord(i,(nano2%centers(:,i)-1.)*nano2%smpd)
                    call centers_close1%set_coord(i,(nano1%centers(:,location(1))-1.)*nano1%smpd)
                    ! remove the atoms from the pdb file
                    call centers_coupled1%set_occupancy(i,0.)
                    call centers_coupled2%set_occupancy(i,0.)
                    call couples1%set_occupancy(i,0.)
                elseif(dist(i)*nano2%smpd > CLOSE_THRESH .and. dist(i)*nano2%smpd<=2.*nano2%theoretical_radius ) then  !to save the atoms which correspond with a precision in the range [0,220] pm
                    cnt2 = cnt2 + 1
                    call centers_coupled2%set_coord(i,(nano2%centers(:,i)-1.)*nano2%smpd)
                    call centers_coupled1%set_coord(i,(nano1%centers(:,location(1))-1.)*nano1%smpd)
                    call couples1%set_coord(i,(nano2%centers(:,location(1))-1.)*nano2%smpd)
                    mask(location(1)) = .false. ! not to consider the same atom more than once
                    ! remove the atoms from the pdb file
                    call centers_close1%set_occupancy(i,0.)
                    call centers_close2%set_occupancy(i,0.)
                endif
                dist_sq(i) = dist(i)**2 !formula wants them square, could improve performance here
                if(DEBUG) then
                     write(logfhandle,*) 'ATOM', i,'coords: ', nano2%centers(:,i), 'coupled with '
                     write(logfhandle,*) '    ',location, 'coordinates: ', nano1%centers(:,location(1)), 'DIST^2= ', dist_sq(i), 'DIST = ', dist(i)
                endif
            enddo
        else
            N_min = size(nano2%centers, dim = 2)
            N_max = size(nano1%centers, dim = 2)
            write(unit = filnum, fmt = '(a,i3)') 'NB atoms in vol1   ', N_max
            write(unit = filnum, fmt = '(a,i3)') 'NB atoms in vol2   ', N_min
            write(unit = filnum, fmt = '(a,i3,a)') '--->', abs(N_max-N_min), ' atoms do NOT correspond'
            allocate(dist(N_max), dist_sq(N_max), source = 0.)
            allocate(dist_close(N_max), source = 0.) ! there are going to be unused entry of the vector
            call centers_coupled1%new(N_max, dummy=.true.)
            call centers_coupled2%new(N_max, dummy=.true.)
            call centers_close1%new  (N_max, dummy=.true.)
            call centers_close2%new  (N_max, dummy=.true.)
            call couples1%new        (N_max, dummy=.true.)
            allocate(mask(N_min), source = .true.)
            cnt  = 0
            cnt2 = 0
            cnt3 = 0
            do i = 1, N_max !compare based on centers1
                dist(i) = pixels_dist(nano1%centers(:,i),nano2%centers(:,:),'min',mask,location, keep_zero = .true.)
                if(dist(i)*nano2%smpd > 2.*nano2%theoretical_radius) then
                    dist(i) = 0.
                    cnt = cnt + 1
                    call couples1%set_coord(i,(nano1%centers(:,location(1))-1.)*nano1%smpd)
                    ! remove the atoms from the pdb file
                    call centers_coupled1%set_occupancy(i,0.)
                    call centers_coupled2%set_occupancy(i,0.)
                    call centers_close1%set_occupancy(i,0.)
                    call centers_close2%set_occupancy(i,0.)
                    call couples1%set_occupancy(i,0.)
                elseif(dist(i)*nano2%smpd <= CLOSE_THRESH) then
                    cnt3 = cnt3 + 1
                    dist_close(i) = dist(i)**2
                    call centers_close1%set_coord(i,(nano1%centers(:,i)-1.)*nano1%smpd)
                    call centers_close2%set_coord(i,(nano2%centers(:,location(1))-1.)*nano2%smpd)
                    ! remove the atoms from the pdb file
                    call centers_coupled1%set_occupancy(i,0.)
                    call centers_coupled2%set_occupancy(i,0.)
                    call couples1%set_occupancy(i,0.)
                elseif(dist(i)*nano2%smpd > CLOSE_THRESH .and. dist(i)*nano2%smpd<=2.*nano2%theoretical_radius ) then  !to save the atoms which correspond with a precision in the range [0,2*theoretical_radius] pm
                    cnt2 = cnt2 + 1
                    call centers_coupled1%set_coord(i,(nano1%centers(:,i)-1.)*nano1%smpd)
                    call centers_coupled2%set_coord(i,(nano2%centers(:,location(1))-1.)*nano2%smpd)
                    call couples1%set_coord(i,(nano1%centers(:,location(1))-1.)*nano1%smpd)
                    mask(location(1)) = .false. ! not to consider the same atom more than once
                    ! remove the atoms from the pdb file
                    call centers_close1%set_occupancy(i,0.)
                    call centers_close2%set_occupancy(i,0.)
                endif
                dist_sq(i) = dist(i)**2 !formula wants them square
                if(DEBUG) then
                    write(logfhandle,*) 'ATOM', i,'coordinates: ', nano1%centers(:,i), 'coupled with '
                    write(logfhandle,*) '    ',location, 'coordinates: ', nano2%centers(:,location(1)), 'DIST^2= ', dist_sq(i), 'DIST = ', dist(i)
                endif
            enddo
        endif
        write(unit = filnum, fmt = '(i3,a,i2,a)')      cnt3,' atoms correspond within        1 A. (', cnt3*100/N_min, '% of the atoms )'
        write(unit = filnum, fmt = '(i3,a,f3.1,a,i2,a)')  cnt2,' atoms correspond within  1 - ',2.*nano2%theoretical_radius,' A. (', cnt2*100/N_min, '% of the atoms )'
        write(unit = filnum, fmt = '(i3,a,f3.1,a,i2,a)')  cnt-(N_max-N_min),' atoms have error bigger than ',2.*nano2%theoretical_radius,' A. (',(cnt-N_max+N_min)*100/N_min, '% of the atoms )' !remove the extra atoms
        ! remove unused atoms from the pdb file
        do i = 1, N_max
            coord(:) = centers_close1%get_coord(i)
            if(coord(1)<TINY .and. coord(2)<TINY .and. coord(3)<TINY) call centers_close1%set_occupancy(i,0.)
            coord(:) = centers_close2%get_coord(i)
            if(coord(1)<TINY .and. coord(2)<TINY .and. coord(3)<TINY) call centers_close2%set_occupancy(i,0.)
            coord(:) = centers_coupled1%get_coord(i)
            if(coord(1)<TINY .and. coord(2)<TINY .and. coord(3)<TINY) call centers_coupled1%set_occupancy(i,0.)
            coord(:) = centers_coupled2%get_coord(i)
            if(coord(1)<TINY .and. coord(2)<TINY .and. coord(3)<TINY) call centers_coupled2%set_occupancy(i,0.)
        enddo
        call centers_close1%writepdb  (trim(nano1%fbody)//'_atom_close_couples')
        call centers_close2%writepdb  (trim(nano2%fbody)//'_atom_close_couples')
        call centers_coupled1%writepdb(trim(nano1%fbody)//'_atom_couples')
        call centers_coupled2%writepdb(trim(nano2%fbody)//'_atom_couples')
        call couples1%writepdb('extra_atoms')
        !Avg dist and stdev symmetry breaking atoms from the center
        !Max dist atoms from the center
        avg     = 0.
        cnt3    = 0
        m(:)    = nano1%nanopart_masscen()
        m(:)    = (m(:)-1.)*nano1%smpd
        tmp_max = 0.
        do i = 1, N_max
            coord(:) = centers_coupled1%get_coord(i)
            if(coord(1)>TINY .and. coord(2)>TINY .and. coord(3)>TINY) then
                cnt3 = cnt3 + 1
                avg = avg + euclid(coord,m)
                if(euclid(coord,m) > tmp_max) tmp_max = euclid(coord,m)
            endif
        enddo
        avg   = avg/real(cnt3)
        cnt3  = 0
        stdev = 0.
        do i = 1, N_max
            coord(:) = centers_coupled1%get_coord(i)
            if(coord(1)>TINY .and. coord(2)>TINY .and. coord(3)>TINY) then
                cnt3 = cnt3 + 1
                stdev = stdev + (euclid(coord,m)-avg)**2
            endif
        enddo
        stdev = sqrt(stdev/(real(cnt3)-1.))
        write(unit = filnum, fmt = '(a)')       '--->    IN VOL1    <---'
        write(unit = filnum, fmt = '(a,f6.3,a)')'AVG     DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', avg, ' A'
        write(unit = filnum, fmt = '(a,f6.3,a)')'STDEV   DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', stdev, ' A'
        write(unit = filnum, fmt = '(a,f6.3,a)')'MAX     DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', tmp_max, ' A'
        tmp_max = 0. ! reset
        do i = 1, size(nano1%centers, dim = 2)
            coord(:) = (nano1%centers(:,i)-1.)*nano1%smpd
                d =  euclid(coord,m)
                if(d > tmp_max) tmp_max = d
        enddo
        write(unit = filnum, fmt = '(a,f6.3,a)')'MAX     DIST ATOMS IN VOL1 3D RECONSTRUCT. TO THE CENTER: ', tmp_max, ' A'
        avg     = 0.
        cnt3    = 0
        m(:)    = nano2%nanopart_masscen()
        m(:)    = (m(:)-1.)*nano2%smpd
        tmp_max = 0.
        do i = 1, N_max
            coord(:) = centers_coupled2%get_coord(i)
            if(coord(1)>TINY .and. coord(2)>TINY .and. coord(3)>TINY) then
                cnt3 = cnt3 + 1
                avg = avg + euclid(coord,m)
                if(euclid(coord,m) > tmp_max) tmp_max = euclid(coord,m)
            endif
        enddo
        avg   = avg/real(cnt3)
        cnt3  = 0
        stdev = 0.
        do i = 1, N_max
            coord(:) = centers_coupled2%get_coord(i)
            if(coord(1)>TINY .and. coord(2)>TINY .and. coord(3)>TINY) then
                cnt3 = cnt3 + 1
                stdev = stdev + (euclid(coord,m)-avg)**2
            endif
        enddo
        stdev = sqrt(stdev/(real(cnt3)-1.))
        write(unit = filnum, fmt = '(a)')       '--->    IN VOL2    <---'
        write(unit = filnum, fmt = '(a,f6.3,a)')'AVG     DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', avg, ' A'
        write(unit = filnum, fmt = '(a,f6.3,a)')'STDEV   DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', stdev, ' A'
        write(unit = filnum, fmt = '(a,f6.3,a)')'MAX     DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', tmp_max, ' A'
        tmp_max = 0. ! reset
        do i = 1, size(nano2%centers, dim = 2)
            coord(:) = (nano2%centers(:,i)-1.)*nano2%smpd
                d =  euclid(coord,m)
                if(d > tmp_max) tmp_max = d
        enddo
        write(unit = filnum, fmt = '(a,f6.3,a)')'MAX     DIST ATOMS IN VOL2 3D RECONSTRUCT. TO THE CENTER: ', tmp_max, ' A'
        ! kill atoms instances
        call centers_close1%kill
        call centers_close2%kill
        call centers_coupled1%kill
        call centers_coupled2%kill
        call couples1%kill
        !RMSD
        rmsd = sqrt(sum(dist_sq)/real(count(dist_sq > TINY)))
        write(unit = filnum, fmt = '(a,f6.3,a)') 'RMSD CALCULATED CONSIDERING ALL ATOMS = ', rmsd*nano1%smpd, ' A'
        write(unit = filnum, fmt = '(a,f6.3,a)') 'RMSD ATOMS THAT CORRESPOND WITHIN 1 A = ', (sqrt(sum(dist_close)/real(count(dist_close > TINY))))*nano1%smpd, ' A'
        call fclose(filnum)
        dist_no_zero = pack(dist, dist>TINY)
        dist_no_zero = dist_no_zero*nano1%smpd ! report distances in Amstrongs
        call hist(dist_no_zero, 50)
        !For CSV files
        call fopen(filnum, file='RMSDhist.csv')
        write (filnum,*) 'r'
        do i = 1, size(dist_no_zero)
            write (filnum,'(A)', advance='yes') trim(real2str(dist_no_zero(i)))
        end do
        call fclose(filnum)
        if(present(r)) r=rmsd
        deallocate(dist, dist_sq, dist_no_zero, mask)
    end subroutine atoms_rmsd

    ! Check for symmetry at different radii.
    ! radii: 5, 7, 9, 12 A
    ! The idea is to identify the atomic positions
    ! contained in a certain radius, generate a distribution
    ! based on that and check for symmetry.
    subroutine identify_atomic_pos(self, cs_thresh)
      class(nanoparticle), intent(inout) :: self
      integer, optional,       intent(in)    :: cs_thresh
      ! Phase correlations approach
      call self%phasecorrelation_nano_gaussian()
      ! Nanoparticle binarization
      call self%binarize()
      ! Outliers discarding
      if(present(cs_thresh)) then
        call self%discard_outliers(cs_thresh)
      else
        call self%discard_outliers()
      endif
      ! Validation of the selected atomic positions
       call self%validate_atomic_positions()
    end subroutine identify_atomic_pos

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
      m = self%nanopart_masscen()
      ! Count nb of atoms in the selected radius
      do n_cc = 1, self%n_cc
            d = euclid(self%centers(:,n_cc), m)*self%smpd
           if(d<=radius) then
               cnt = cnt+1
           endif
      enddo
      call radial_atom%new (cnt, dummy=.true.)
      cnt = 0
      ! Fill in radial_atom with the atomic positions
      do n_cc = 1, self%n_cc
          d = euclid(self%centers(:,n_cc), m)*self%smpd
          if(d<=radius) then
             cnt = cnt+1
             call radial_atom%set_element(cnt,element)
             call radial_atom%set_name(cnt,atom_name)
             call radial_atom%set_coord(cnt,(self%centers(:,n_cc)-1.)*self%smpd)
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
      real        :: m(3), mm(3), vec(3), d, d_before
      integer     :: i
      !!!!!!!!!!! DEBUG !!!!!!!!!!
      type(atoms) :: center_of_mass
      integer :: index
      real :: coords(3)
      call atom_centers%new(pdbfile_in)
      m(:)  = self%nanopart_masscen()
      d_before = huge(d_before)
      do i = 1, self%n_cc
        d = euclid(m,self%centers(:,i))
        if(d < d_before) then
          vec(:)   = m(:) - self%centers(:,i)
          d_before = d
          index = i
        endif
      enddo
      print *, 'self%centers(:, index) + vec: ', self%centers(:, index) + vec, 'm(:): ', m(:)
      print *, 'self%smpd:', self%smpd
      vec = (vec-1.)*self%smpd
      m   = (m-1.)*self%smpd
      coords = atom_centers%get_coord(index)
      print *, 'coords + vec: ', coords + vec, '(m(:)-1.)*self%smpd: ', m
      call atom_centers%translate(vec)
      call atom_centers%writePDB(pdbfile_out)
      !!!!!!!!!
      call center_of_mass%new(1)
      call center_of_mass%set_coord(1,m)
      call center_of_mass%writePDB('CenterOfMass')
      call center_of_mass%kill
      !!!!!!!!!
      call atom_centers%kill
    end subroutine center_on_atom

    subroutine geometry_analysis(self, pdbfile2, outfile)
      use simple_math, only : plane_from_points
      class(nanoparticle),        intent(inout) :: self
      character(len=*),           intent(in)    :: pdbfile2  ! atomic pos of the 2 (or 3) selected atoms
      character(len=*), optional, intent(in)    :: outfile
      type(atoms)           :: init_atoms, final_atoms
      type(binimage)        :: img_out
      integer, allocatable  :: imat_cc(:,:,:), imat(:,:,:)
      real,    allocatable  :: line(:,:), plane(:,:,:),points(:,:), pointsCopy(:,:),distances_totheplane(:)
      real,    allocatable  :: radii(:), max_intensity(:), int_intensity(:), w(:), v(:,:), d(:)
      integer :: i, j, n, t, s, filnum, io_stat, cnt_intersect, cnt, n_cc
      logical :: flag(self%n_cc)
      real    :: atom1(3), atom2(3), atom3(3), dir_1(3), dir_2(3), vec(3), m(3), dist_plane, dist_line
      real    :: t_vec(N_DISCRET), s_vec(N_DISCRET), denominator, t1, t2
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
              dist_line = euclid(self%centers(:3,i),line(:3,t))
              if(dist_line*self%smpd <= 0.9*self%theoretical_radius) then ! it intersects atoms i
                   flag(i) = .true. !flags also itself
               endif
            enddo
        enddo
        write(logfhandle, *) 'generating volumes for visualization'
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
               call final_atoms%set_coord(cnt_intersect,(self%centers(:,i)-1.)*self%smpd)
               where(imat_cc == i) imat = 1
           endif
       enddo
       call img_out%set_imat(imat)
       if(present(outfile)) then
           call img_out%write_bimg(outfile)
       else
         call img_out%write_bimg('ImageColumn.mrc')
       endif
       call final_atoms%writePDB('AtomColumn')
       call final_atoms%kill
       ! Find the plane that best fits the atoms belonging to the line
       allocate(points(3, count(flag)), source = 0.)
       m = self%nanopart_masscen()
       cnt = 0
       do i = 1, self%n_cc
           if(flag(i)) then
               cnt = cnt + 1
               points(:3,cnt) = self%centers(:3,i)-m(:)
           endif
       enddo
       !!!!!! SVDFIT!!!!!!! TO COMPLETE
       ! allocate(pointsCopy(3,count(flag)), source = points) ! because svdcmp modifies its input
       ! allocate(w(3), v(3,3), source = 0.)
       ! allocate(d(count(flag)), source = 0.)
       ! call svdcmp(points,w,v)
       ! print *, 'u'
       ! ! call vis_mat(points)
       ! print *, 'w', w
       ! print *, 'v'
       ! call vis_mat(v)
       ! d = A(:,2)
       ! print *, 'directional vector of the line', d


       ! line
       ! line(1,t) = m(1) + t_vec(t)* d(1)
       ! line(2,t) = m(2) + t_vec(t)* d(2)
       ! line(3,t) = m(3) + t_vec(t)* d(3)

       ! t   = matmul(d,copyA)
       ! t1  = minval(t)
       ! t2  = maxval(t)
       ! print *, 't', t
       ! print *, 't1', t1
       ! print *, 't2', t2

       stop
       !!!!!!!!!!!!!!!!!!!!!!
       ! TO FIX
       vec = plane_from_points(points)
       allocate(distances_totheplane(cnt), source = 0.)
       allocate(radii(cnt), source = 0.) ! which radius is the atom center belonging to
       cnt = 0
       denominator = sqrt(vec(1)*2+vec(2)*2+1.)
       write(logfhandle,*) 'Directional vector: [', -1., ',', -1., ',', -vec(1)-vec(2), ']'
       ! Read ratios, ang_var, dists, max_intensity
       if(.not. allocated(self%ratios)) allocate(self%ratios(self%n_cc))
       call fopen(filnum, file='../'//'Ar.csv', iostat=io_stat)
       if( io_stat .ne. 0 ) then
         THROW_HARD('Unable to read file Ar.csv Did you run cluster_analysis?; geometry_analysis')
       endif
       read(filnum,*) ! first line is variable name
       do i = 1, self%n_cc
         read(filnum,*) self%ratios(i)
       enddo
       call fclose(filnum)
       if(.not. allocated(self%ang_var)) allocate(self%ang_var(self%n_cc))
       call fopen(filnum, file='../'//'Ang.csv', iostat=io_stat)
       if( io_stat .ne. 0 ) then
         THROW_HARD('Unable to read file Ang.csv Did you run cluster_analysis?; geometry_analysis')
       endif
       read(filnum,*)
       do i = 1, self%n_cc
         read(filnum,*) self%ang_var(i)
       enddo
       call fclose(filnum)
       if(.not. allocated(self%dists)) allocate(self%dists(self%n_cc))
       call fopen(filnum, file='../'//'Dist.csv', iostat=io_stat)
       if( io_stat .ne. 0 ) then
         THROW_HARD('Unable to read file Dist.csv Did you run cluster_analysis?; geometry_analysis')
       endif
       read(filnum,*)
       do i = 1, self%n_cc
         read(filnum,*) self%dists(i)
       enddo
       call fclose(filnum)
       allocate(max_intensity(self%n_cc))
       call fopen(filnum, file='../'//'MaxIntensity.csv', iostat=io_stat)
       if( io_stat .ne. 0 ) then
         THROW_HARD('Unable to read file MaxIntensity.csv Did you run radial_dependent_stats?; geometry_analysis')
       endif
       read(filnum,*)
       do i = 1, self%n_cc
         read(filnum,*) max_intensity(i)
       enddo
       call fclose(filnum)
       allocate(int_intensity(self%n_cc))
       call fopen(filnum, file='../'//'IntIntensity.csv', iostat=io_stat)
       if( io_stat .ne. 0 ) then
         THROW_HARD('Unable to read file IntIntensity.csv Did you run radial_dependent_stats?; geometry_analysis')
       endif
       read(filnum,*)
       do i = 1, self%n_cc
         read(filnum,*) int_intensity(i)
       enddo
       call fclose(filnum)
       ! Output on Excel file all the stats on the atoms belonging to the plane
       cnt = 0
       call fopen(filnum, file='AtomColumnInfo.xls',iostat=io_stat)
       write (filnum,*) '        Atom #    ','   Ar        ','       Dist     ','         Ang     ','         MaxInt     ', '          IntInt     '
       do i = 1, self%n_cc
         if(flag(i)) then
           cnt = cnt + 1
             write (filnum,*) i, '   ', self%ratios(i), self%dists(i), self%ang_var(i), max_intensity(i), int_intensity(i)
         endif
       end do
       call fclose(filnum)
       ! do i = 1, self%n_cc
       !     if(flag(i)) then
       !         cnt = cnt + 1
       !         ! formula for distance of a point to a plane
       !         distances_totheplane(cnt) = abs(vec(1)*points(1,cnt)+vec(2)*points(2,cnt)-points(3,cnt)+vec(3))/denominator
       !         radii(cnt) = euclid(self%centers(:,i), m)*self%smpd
       !     endif
       ! enddo
       ! distances_totheplane = (distances_totheplane)*self%smpd
       ! call fopen(filnum, file='Radii.csv', iostat=io_stat)
       ! write (filnum,*) 'r'
       ! do i = 1, cnt
       !     write (filnum,'(A)', advance='yes') trim(real2str(radii(i)))
       ! end do
       ! call fclose(filnum)
       ! call fopen(filnum, file='DistancesToThePlane.csv',iostat=io_stat)
       ! write (filnum,*) 'd'
       ! do i = 1, cnt
       !     write (filnum,'(A)', advance='yes') trim(real2str(distances_totheplane(i)))
       ! end do
       ! call fclose(filnum)
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
                dist_plane = euclid(self%centers(:3,i),plane(:3,t,s))
                if(dist_plane*self%smpd <= 0.9*self%theoretical_radius) then ! it intersects atoms i
                     flag(i) = .true. !flags also itself
                 endif
              enddo
            enddo
        enddo
        write(logfhandle, *) 'generating volumes for visualization'
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
                call final_atoms%set_coord(cnt_intersect,(self%centers(:,i)-1.)*self%smpd)
                where(imat_cc == i) imat = 1
            endif
        enddo
        call img_out%set_imat(imat)
        if(present(outfile)) then
            call img_out%write_bimg(outfile)
        else
          call img_out%write_bimg('ImagePlane.mrc')
        endif
        call final_atoms%writePDB('AtomPlane')
        call final_atoms%kill
        allocate(points(3, count(flag)), source = 0.)
        m = self%nanopart_masscen()
        cnt = 0
        do i = 1, self%n_cc
            if(flag(i)) then
                cnt = cnt + 1
                points(:3,cnt) = self%centers(:3,i)-m(:)
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
                radii(cnt) = euclid(self%centers(:,i), m)*self%smpd
            endif
        enddo
        distances_totheplane = (distances_totheplane)*self%smpd
        call fopen(filnum, file='Radii.csv', iostat=io_stat)
        write (filnum,*) 'r'
        do i = 1, cnt
            write (filnum,'(A)', advance='yes') trim(real2str(radii(i)))
        end do
        call fclose(filnum)
        call fopen(filnum, file='DistancesToThePlane.csv',iostat=io_stat)
        write (filnum,*) 'd'
        do i = 1, cnt
            write (filnum,'(A)', advance='yes') trim(real2str(distances_totheplane(i)))
        end do
        call fclose(filnum)
        ! Read ratios, ang_var, dists, max_intensity
        if(.not. allocated(self%ratios)) allocate(self%ratios(self%n_cc))
        call fopen(filnum, file='../'//'Ar.csv', iostat=io_stat)
        if( io_stat .ne. 0 ) then
          THROW_HARD('Unable to read file Ar.csv Did you run cluster_analysis?; geometry_analysis')
        endif
        read(filnum,*) ! first line is variable name
        do i = 1, self%n_cc
          read(filnum,*) self%ratios(i)
        enddo
        call fclose(filnum)
        if(.not. allocated(self%ang_var)) allocate(self%ang_var(self%n_cc))
        call fopen(filnum, file='../'//'Ang.csv', iostat=io_stat)
        if( io_stat .ne. 0 ) then
          THROW_HARD('Unable to read file Ang.csv Did you run cluster_analysis?; geometry_analysis')
        endif
        read(filnum,*)
        do i = 1, self%n_cc
          read(filnum,*) self%ang_var(i)
        enddo
        call fclose(filnum)
        if(.not. allocated(self%dists)) allocate(self%dists(self%n_cc))
        call fopen(filnum, file='../'//'Dist.csv', iostat=io_stat)
        if( io_stat .ne. 0 ) then
          THROW_HARD('Unable to read file Dist.csv Did you run cluster_analysis?; geometry_analysis')
        endif
        read(filnum,*)
        do i = 1, self%n_cc
          read(filnum,*) self%dists(i)
        enddo
        call fclose(filnum)
        allocate(max_intensity(self%n_cc))
        call fopen(filnum, file='../'//'MaxIntensity.csv', iostat=io_stat)
        if( io_stat .ne. 0 ) then
          THROW_HARD('Unable to read file MaxIntensity.csv Did you run radial_dependent_stats?; geometry_analysis')
        endif
        read(filnum,*)
        do i = 1, self%n_cc
          read(filnum,*) max_intensity(i)
        enddo
        call fclose(filnum)
        allocate(int_intensity(self%n_cc))
        call fopen(filnum, file='../'//'IntIntensity.csv', iostat=io_stat)
        if( io_stat .ne. 0 ) then
          THROW_HARD('Unable to read file IntIntensity.csv Did you run radial_dependent_stats?; geometry_analysis')
        endif
        read(filnum,*)
        do i = 1, self%n_cc
          read(filnum,*) int_intensity(i)
        enddo
        call fclose(filnum)
        ! Output on Excel file all the stats on the atoms belonging to the plane
        cnt = 0
        call fopen(filnum, file='AtomPlaneInfo.xls',iostat=io_stat)
        write (filnum,*) '        Atom #    ','   Ar        ','       Dist     ','         Ang     ','      MaxInt     ','       IntInt     '
        do i = 1, self%n_cc
          if(flag(i)) then
            cnt = cnt + 1
              write (filnum,*) i, '   ', self%ratios(i), self%dists(i), self%ang_var(i), max_intensity(i), int_intensity(i)
          endif
        end do
        call fclose(filnum)
      endif
      call img_out%kill_bimg
      call init_atoms%kill
      if(allocated(line))  deallocate(line)
      if(allocated(plane)) deallocate(plane)
      deallocate(imat_cc, imat)
    end subroutine geometry_analysis

    subroutine kill_nanoparticle(self)
        class(nanoparticle), intent(inout) :: self
        self%ldim(3)           = 0
        self%smpd              = 0.
        self%nanop_mass_cen(3) = 0.
        self%avg_dist_atoms    = 0.
        self%n_cc              = 0
        call self%img%kill()
        call self%img_raw%kill
        call self%img_bin%kill_bimg()
        call self%img_cc%kill_bimg()
        call self%centers_pdb%kill
        if(allocated(self%centers))          deallocate(self%centers)
        if(allocated(self%ratios))           deallocate(self%ratios)
        if(allocated(self%dists))            deallocate(self%dists)
        if(allocated(self%loc_longest_dist)) deallocate(self%loc_longest_dist)
        if(allocated(self%contact_scores))   deallocate(self%contact_scores)
        if(allocated(self%ang_var))          deallocate(self%ang_var)
    end subroutine kill_nanoparticle

end module simple_nanoparticle
