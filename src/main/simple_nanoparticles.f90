!USAGE 1: simple_exec prg=detect_atoms smpd=0.358 vol1='particle1.mrc'
!USAGE 2: simple_exec prg=compare_nano vol1=vol_c1_aligned2_c3axis.mrc vol2=particle1_c3.mrc
module simple_nanoparticles_mod
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,         only : image
use simple_atoms,         only : atoms
implicit none

private
public :: nanoparticle

#include "simple_local_flags.inc"

! module global constants
integer, parameter :: N_THRESH    = 20    !number of thresholds for binarization
logical, parameter :: DEBUG_HERE  = .false.!for debugging purposes
integer, parameter :: SOFT_EDGE   = 6

type :: nanoparticle
    private
    type(atoms) :: centers_pdb
    type(image) :: img, img_bin, img_cc, img_raw
    integer     :: ldim(3)            = 0
    real        :: smpd               = 0.
    real        :: nanop_mass_cen(3)  = 0.!coordinates of the center of mass of the nanoparticle
    real        :: avg_dist_atoms     = 0.
    real        :: theoretical_radius = 0.!theoretical atom radius
    integer     :: n_cc               = 0 !number of atoms (connected components
    real,    allocatable  :: centers(:,:)
    real,    allocatable  :: ratios(:)
    real,    allocatable  :: ang_var(:)
    real,    allocatable  :: dists(:)
    integer, allocatable  :: loc_longest_dist(:,:)   !for indentific of the vxl that determins the longest dim of the atom
    character(len=2)      :: element     = ' '
    character(len=4)      :: atom_name   = '    '
    character(len=STDLEN) :: partname    = ''     !fname
    character(len=STDLEN) :: fbody       = ''     !fbody
    character(len=STDLEN) :: output_dir  = ''

  contains
    ! constructor
    procedure          :: new => new_nanoparticle
    ! getters/setters
    procedure          :: get_img
    procedure          :: get_ldim
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
    procedure, private :: affprop_cluster_ar
    procedure, private :: affprop_cluster_ang
    procedure, private :: affprop_cluster_dist_distr
    procedure          :: cluster_atom_intensity
    procedure          :: search_polarization
    procedure          :: identify_atom_columns
    procedure          :: identify_atom_planes
    ! execution
    procedure          :: identify_atomic_pos
    procedure          :: radial_dependent_stats
    procedure          :: cluster
    procedure          :: atoms_rmsd
    procedure          :: make_soft_mask
    ! visualization and output
    procedure          :: write_centers
    procedure          :: print_asym_unit
    ! others
    procedure, private :: update_self_ncc
    procedure          :: keep_atomic_pos_at_radius
    ! kill
    procedure          :: kill => kill_nanoparticle
end type nanoparticle

contains

    !constructor
    subroutine new_nanoparticle(self, fname, cline_smpd, element)
        use simple_syslib
        class(nanoparticle),      intent(inout) :: self
        character(len=*),         intent(in)    :: fname
        real,                     intent(in)    :: cline_smpd
        character(len=2),optional,intent(in)    :: element
        integer :: nptcls
        real    :: smpd
        call self%kill
        call simple_getcwd(self%output_dir)
        call self%set_partname(fname)
        self%fbody = get_fbody(trim(basename(fname)), trim(fname2ext(fname)))
        self%smpd  = cline_smpd
        if(.not. present(element)) then
            self%element   = 'pt'  !default is pt
            self%atom_name = ' pt '
            self%theoretical_radius = 1.1
        else
            select case(element)
                case('pt')
                    self%element       = 'pt'
                    self%atom_name     = ' pt '
                    self%theoretical_radius = 1.1
                    ! thoretical radius is already set
                case('pd')
                    self%element       = 'pd'
                    self%atom_name     = ' pd '
                    self%theoretical_radius = 1.12
                case('fe')
                    self%element       = 'fe'
                    self%atom_name     = ' fe '
                    self%theoretical_radius = 1.02
                case('au')
                    self%element       = 'au'
                    self%atom_name     = ' au '
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

     !get one of the images of the nanoparticle type
     subroutine get_img( self, img, which )
         class(nanoparticle), intent(in)    :: self
         type(image),         intent(inout) :: img
         character(len=*),    intent(in)    :: which
         integer :: ldim(3)
         ldim = img%get_ldim()
         if(self%ldim(1) .ne. ldim(1) .or. self%ldim(2) .ne. ldim(2) .or. self%ldim(3) .ne. ldim(3)) THROW_HARD('Wrong dimension in input img; get_img')
         select case(which)
             case('img')
                 img = self%img
             case('img_bin')
                 img = self%img_bin
             case('img_cc')
                 img = self%img_cc
             case DEFAULT
                 THROW_HARD('Wrong input parameter img type; get_img')
         end select
     end subroutine get_img

     subroutine get_ldim(self,ldim)
       class(nanoparticle), intent(in)  :: self
       integer,             intent(out) :: ldim(3)
       ldim = self%img%get_ldim()
     end subroutine get_ldim

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
                 call self%img_bin%new(self%ldim, self%smpd)
                 call self%img_bin%read(imgfile)
             case('img_cc')
                 call self%img_cc%new(self%ldim, self%smpd)
                 call self%img_cc%read(imgfile)
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
        class(nanoparticle),   intent(inout) :: self
        type(image), optional, intent(in)    :: img_cc
        real, pointer :: rmat_cc(:,:,:)
        if(present(img_cc)) then
            call img_cc%get_rmat_ptr(rmat_cc)
        else
            call self%img_cc%get_rmat_ptr(rmat_cc)
        endif
        self%n_cc = nint(maxval(rmat_cc))
        rmat_cc => null()
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

    ! Find the centers coordinates of the atoms in the particle
    ! and save it in the global variable centers.
    ! If coords is present, it saves it also in coords.
    subroutine find_centers(self, img_bin, img_cc, coords)
       class(nanoparticle),          intent(inout) :: self
       type(image), optional,        intent(inout) :: img_bin, img_cc
       real, optional, allocatable,  intent(out)   :: coords(:,:)
       real,        pointer :: rmat_cc_in(:,:,:), rmat_raw(:,:,:)
       logical, allocatable :: mask(:,:,:)
       integer :: i, ii, jj, kk
       real    :: m(3), sum_mass
       ! sanity check
       if(present(img_bin) .and. .not. present(img_cc)) THROW_HARD('img_bin and img_cc have to be both present in input')
       ! global variables allocation
       if(allocated(self%centers)) deallocate(self%centers)
       allocate(self%centers(3,self%n_cc),source = 0.)
       if(present(img_cc)) then
           call img_cc%get_rmat_ptr(rmat_cc_in)
       else
           call self%img_cc%get_rmat_ptr(rmat_cc_in)
       endif
       call self%img_raw%get_rmat_ptr(rmat_raw)
       allocate(mask(self%ldim(1),self%ldim(2),self%ldim(3)), source = .true.)
       !$omp parallel do default(shared) private(i,ii,jj,kk,mask,m,sum_mass) schedule(static) proc_bind(close)
       do i=1,self%n_cc
           mask = .true.
           where(abs(rmat_cc_in(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))-real(i)) > TINY) mask = .false.
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
         ! H&C
         ! m = 0.
         ! !$omp parallel do schedule(static) reduction(+:m) private(i)
         ! do i = 1, self%n_cc
         !      m = m + 1.*self%centers(:,i)
         ! enddo
         ! !$omp end parallel do
         !
         ! m = m/real(self%n_cc)
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
        if(DEBUG_HERE) then
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
        if(DEBUG_HERE) then
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
        type(image)       :: img_bin_thresh(N_THRESH/2-1)
        type(image)       :: img_ccs_thresh(N_THRESH/2-1)
        type(image)       :: pc
        type(atoms)       :: atom
        type(image)       :: simulated_distrib
        integer, allocatable :: imat_t(:,:,:)
        real,    allocatable :: x_mat(:)  !vectorization of the volume
        real,    allocatable :: coords(:,:)
        real,    allocatable :: rmat(:,:,:)
        integer :: i, cc, t
        real    :: otsu_thresh
        real    :: x_thresh(N_THRESH/2-1)
        real    :: step, maximum, mm(2)
        call otsu_nano(self%img,otsu_thresh)
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
            call img_bin_thresh(i)%new(self%ldim, self%smpd)
            if(i == 1) then
                x_thresh(i) = otsu_thresh
            else
                x_thresh(i) = x_thresh(i-1) + step
            endif
            where(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) > x_thresh(i))
                imat_t = 1
            elsewhere
                imat_t = 0
            endwhere
            ! Generate binary image and cc image
            call img_bin_thresh(i)%set_rmat(real(imat_t))
            call img_bin_thresh(i)%find_connected_comps(img_ccs_thresh(i))
            ! Find atom centers in the generated distributions
            call self%update_self_ncc(img_ccs_thresh(i)) ! self%n_cc is needed in find_centers
            call self%find_centers(img_bin_thresh(i), img_ccs_thresh(i), coords)
            ! Generate a simulated distribution based on those center
            call self%write_centers('centers_'//trim(int2str(i))//'_iteration', coords)
            call atom%new          ('centers_'//trim(int2str(i))//'_iteration.pdb')
            call atom%convolve(simulated_distrib, cutoff = 8.*self%smpd)
            call del_file('centers_'//trim(int2str(i))//'_iteration.pdb')
            if(DEBUG_HERE) call simulated_distrib%write('simulated_'//trim(int2str(i))//'_iteration.mrc')
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
        call self%img_bin%copy(img_bin_thresh(t))
        call self%img_cc%copy(img_ccs_thresh(t))
        do i = 1,  N_THRESH/2-1
            call img_bin_thresh(i)%kill
            call img_ccs_thresh(i)%kill
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
    subroutine discard_outliers(self)
        class(nanoparticle), intent(inout) :: self
        real, pointer      :: rmat(:,:,:), rmat_cc(:,:,:)
        integer, parameter :: PERCENT_DISCARD = 5/100
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
        n_discard = self%n_cc*PERCENT_DISCARD
        mask(1:self%n_cc) = .true.
        call self%img_cc%get_rmat_ptr(rmat_cc)
        call self%img_bin%get_rmat_ptr(rmat)
        ! Removing outliers from the binary image and the connected components image
        !omp parallel do schedule(static) private(cc,label,rmat,mask) !not sure about mask
        do cc = 1, n_discard
            label = minloc(contact_scores, mask)
            where(abs(rmat_cc-real(label(1))) < TINY)
               rmat = 0.   ! discard atom corresponding to label_discard(i)
            endwhere
            mask(label(1)) = .false.
            self%centers(:,label) = -1. !identify to discard them
        enddo
        !omp end parallel do
        call self%img_bin%find_connected_comps(self%img_cc)
        call self%img_cc%get_rmat_ptr(rmat_cc)
        self%n_cc = nint(maxval(rmat_cc)) !update
        call self%find_centers()
        !root folder
        call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: discard_outliers, simple_chdir; ")
        write(logfhandle, *) '****outliers discarding, completed'
    end subroutine discard_outliers

    ! This subrouitne validates the identified atomic positions
    subroutine validate_atomic_positions(self)
        class(nanoparticle), intent(inout) :: self
        integer, allocatable :: imat(:,:,:)
        real,    allocatable :: x(:)
        real,    pointer     :: rmat_pc(:,:,:), rmat_cc(:,:,:), rmat_bin(:,:,:)
        integer, parameter   :: RANK_THRESH = 4
        integer :: n_cc, cnt
        integer :: rank, m(1)
        real    :: new_centers(3,2*self%n_cc)   !will pack it afterwards if it has too many elements
        real    :: pc
        call self%img%get_rmat_ptr(rmat_pc)    !now img contains the phase correlation
        imat = nint(self%img_cc%get_rmat())
        call self%img_cc%get_rmat_ptr(rmat_cc) !to pass to the subroutine split_atoms
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
                call split_atom(rmat_pc,rmat_cc,imat,n_cc,new_centers,cnt)
            else
                cnt = cnt + 1 !new number of centers deriving from splitting
                new_centers(:,cnt) = self%centers(:,n_cc)
            endif
        enddo
        deallocate(self%centers)
        self%n_cc = cnt !update
        allocate(self%centers(3,cnt), source = 0.)
        !update centers
        do n_cc =1, cnt
            self%centers(:,n_cc) = new_centers(:,n_cc)
        enddo
        call self%img_bin%get_rmat_ptr(rmat_bin)
        call self%img_bin%write(trim(self%fbody)//'BINbeforeValidation.mrc') !if(DEBUG_HERE)
        ! Update binary image
        where(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) > 0.)
            rmat_bin(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 1.
        elsewhere
            rmat_bin(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 0.
        endwhere
        call self%img_bin%write(trim(self%fbody)//'BIN.mrc')
        call self%img_bin%find_connected_comps(self%img_cc)
        call self%img_cc%write(trim(self%fbody)//'CC.mrc')
        ! Update number of ccs
        call self%update_self_ncc()
        ! Update and write centers
        call self%find_centers()
        call self%write_centers()
    contains
        subroutine split_atom(rmat_pc, rmat_cc, imat,n_cc,new_centers,cnt)
            real,    intent(in)    :: rmat_pc(:,:,:)    !correlation matrix
            real,    intent(inout) :: rmat_cc(:,:,:)    !ccs matrix, variable
            integer, intent(in)    :: imat(:,:,:)       !ccs matrix, fixed
            integer, intent(in)    :: n_cc              !label of the cc to split
            real,    intent(inout) :: new_centers(:,:)  !updated coordinates of the centers
            integer, intent(inout) :: cnt               !atom counter, to u pdate the center coords
            integer :: new_center1(3),new_center2(3),new_center3(3)
            integer :: i, j, k
            logical :: mask(self%ldim(1),self%ldim(2),self%ldim(3)) !false in the layer of connection of the atom to be split
            mask = .false.  !initialization
            ! Identify first new center
            new_center1 = maxloc(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), imat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) == n_cc)
            cnt = cnt + 1
            new_centers(:,cnt) = real(new_center1)
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
                & (imat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) == n_cc) .and. .not. mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)))
                if(any(new_center2 > 0)) then !if anything was found
                    !Validate second center (check if it's 2 merged atoms, or one pointy one)

                    ! READABILITY
                    ! sum(real(new_center2-new_center1)**2.)

                    if(sqrt((real(new_center2(1)-new_center1(1))**2+real(new_center2(2)-new_center1(2))**2+real(new_center2(3)-new_center1(3))**2)*self%smpd) <= 2.*self%theoretical_radius) then
                        ! Set the merged cc back to 0
                        where((abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))-real(n_cc))<TINY) .and. (.not.mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)))) rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 0.
                        return
                    else
                      cnt = cnt + 1
                      new_centers(:,cnt) = real(new_center2)
                      !In the case two merged atoms, build the second atom
                      do i = 1, self%ldim(1)
                          do j = 1, self%ldim(2)
                              do k = 1, self%ldim(3)
                                  if(((real(i-new_center2(1)))**2 + (real(j-new_center2(2)))**2 + (real(k-new_center2(3)))**2)*self%smpd < (0.9*self%theoretical_radius)**2) then
                                      if(imat(i,j,k) == n_cc)   mask(i,j,k) = .true.
                                  endif
                              enddo
                          enddo
                      enddo
                    endif
                  endif
                  ! Third likely center.
                      new_center3 = maxloc(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), &
                      & (imat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) == n_cc) .and. .not. mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)))
                      if(any(new_center3 > 0)) then !if anything was found
                          !Validate third center
                          if(sqrt((real(new_center3(1)-new_center1(1))**2+real(new_center3(2)-new_center1(2))**2+real(new_center3(3)-new_center1(3))**2)*self%smpd) <= 2.*self%theoretical_radius .or. &
                          &  sqrt((real(new_center3(1)-new_center2(1))**2+real(new_center3(2)-new_center2(2))**2+real(new_center3(3)-new_center2(3))**2)*self%smpd) <= 2.*self%theoretical_radius ) then
                              ! Set the merged cc back to 0
                              where((abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))-real(n_cc))<TINY) .and. (.not.mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)))) rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 0.
                              return
                          else
                            cnt = cnt + 1
                            new_centers(:,cnt) = real(new_center3)
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
            where((abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))-real(n_cc))<TINY) .and. (.not.mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)))) rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 0.
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
        type(atoms)        :: radial_atoms_just,  radial_atoms_all
        logical, allocatable :: mask(:,:,:)
        real,    allocatable :: coords(:,:) !coordinates of the centers of the atoms according to radial distances
        real,    allocatable :: max_intensity(:), avg_intensity(:), stdev_intensity(:), ratios(:)
        real,    pointer     :: rmat(:,:,:), rmat_cc(:,:,:)
        real    :: m(3)    !mass center of c3 map
        real    :: d       !distance atoms from the center
        real    :: radius  !radius of the sphere to consider
        integer :: cnt, cnt_just, cnt_all
        integer :: nsteps, i, j, k, l, cc
        real    :: avg_int, max_int, stdev_int
        ! min_step and max_step is in A
        self%n_cc = size(self%centers, dim = 2)
        write(logfhandle, *) '****radial atom-to-atom distances estimation, init'
        nsteps = floor((max_rad-min_rad)/step)+1
        m        = self%nanopart_masscen()
        ! Report statistics in dedicated directory
        call simple_mkdir(trim(self%output_dir)//'/DistDistribANDRadialDependentStat',errmsg="simple_nanoparticles :: radial_dependent_stats, simple_mkdir; ")
        call simple_chdir(trim(self%output_dir)//'/DistDistribANDRadialDependentStat',errmsg="simple_nanoparticles :: radial_dependent_stats, simple_chdir; ")
        open(11, file='RadialDependentStat')
        allocate(mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), source = .false.)
        ! allocate(max_intensity(self%n_cc),avg_intensity(self%n_cc),stdev_intensity(self%n_cc), source= 0.)
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
            ! Estimation of avg distance and stdev among atoms in radial dependent shells
            write(unit = 11, fmt = '(a)')   '*********************************************************'
            write(unit = 11, fmt = '(a,a)') 'Estimation of atom-to-atom statistics in shell of radius ', trim(int2str(nint(radius)))
            allocate(coords(3,cnt_all), source = 0.)
            allocate(max_intensity(cnt_all),avg_intensity(cnt_all),stdev_intensity(cnt_all), source = 0. )
            allocate(ratios(cnt_all), source = 0. )
            cnt = 0
            call self%img%get_rmat_ptr(rmat)
            call self%img_cc%get_rmat_ptr(rmat_cc)
            do cc = 1, size(self%centers, dim = 2)
                d = euclid(self%centers(:,cc), m)*self%smpd
                if(d<=radius) then
                    if(l == 1) then
                        cnt = cnt + 1
                        coords(:3,cnt) = self%centers(:,cc)
                        call self%calc_aspect_ratio(cc, ratios(cnt),lld=.false., print_ar=.false.)
                        where(abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) - real(cc)) < TINY) mask = .true.
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
                        write(unit = 11, fmt = '(a,i3,a,f9.5,a,f9.5,a,f9.5,a,f9.5)')  'ATOM ', cc,' maxval ', max_intensity(cnt), '   avg ', avg_intensity(cnt), '   stdev ', stdev_intensity(cnt),' aspect ratio ', ratios(cnt)
                        mask = .false. !Reset
                    elseif(d>(radius-step)) then
                        cnt = cnt + 1
                        coords(:3,cnt) = self%centers(:,cc)
                        call self%calc_aspect_ratio(cc, ratios(cnt),lld=.false., print_ar=.false.)
                        where(abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) - real(cc)) < TINY) mask = .true.
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
                        write(unit = 11, fmt = '(a,i3,a,f9.5,a,f9.5,a,f9.5,a,f9.5)')  'ATOM ', cc,' maxval ', max_intensity(cnt), '   avg ', avg_intensity(cnt), '   stdev ', stdev_intensity(cnt), ' aspect ratio ', ratios(cnt)
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
            write(unit = 11, fmt = '(a,f9.5)') 'Maxval  int shell :', max_int
            write(unit = 11, fmt = '(a,f9.5)') 'Average int shell :', avg_int
            write(unit = 11, fmt = '(a,f9.5)') 'Stdev   int shell :', stdev_int
            call self%distances_distribution(coords)
            deallocate(coords)
            deallocate(avg_intensity)
            deallocate(max_intensity)
            deallocate(stdev_intensity)
            deallocate(ratios)
            call radial_atoms_all%kill
            call radial_atoms_just%kill
        enddo
        close(11)
        deallocate(mask)
        call self%distances_distribution() ! across the whole nano
        ! Come back to root directory
        ! call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: radial_dependent_stats, simple_chdir; ")
        call self%kill
        write(logfhandle, *) '****radial atom-to-atom distances estimation, completed'
   end subroutine radial_dependent_stats

   subroutine distances_distribution(self,coords,volume)
       class(nanoparticle), intent(inout) :: self
       real,    optional,   intent(in)    :: coords(:,:)
       integer, optional,   intent(in)    :: volume
       real, allocatable :: dist(:)
       real    :: stdev, med
       integer :: i, j, n_discard
       logical :: mask(self%n_cc)
       ! Initialisations
       mask       = .true.
       stdev      = 0.
       n_discard  = 0
       if(present(coords)) then
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
                   print *, 'atom ',i, 'dist(i)', dist(i), 'INSERT THRESHOLD FOR MIN DIST BETWEEN ATOMS'
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
               write(unit = 11, fmt = '(a,a,a,f6.3,a)') 'Average dist atoms vol ', trim(int2str(volume)),':', self%avg_dist_atoms*self%smpd, ' A'
               write(unit = 11, fmt = '(a,a,a,f6.3,a)') 'StDev   dist atoms vol ', trim(int2str(volume)),':', stdev*self%smpd, ' A'
               write(unit = 11, fmt = '(a,a,a,f6.3,a)') 'Median  dist atoms vol ', trim(int2str(volume)),':', med*self%smpd, ' A'
           else
               write(unit = 11, fmt = '(a,f6.3,a)') 'Average dist atoms: ', self%avg_dist_atoms*self%smpd, ' A'
               write(unit = 11, fmt = '(a,f6.3,a)') 'StDev   dist atoms: ', stdev*self%smpd, ' A'
               write(unit = 11, fmt = '(a,f6.3,a)') 'Median  dist atoms: ', med*self%smpd, ' A'
           endif
           deallocate(dist)
       else
           open(15, file='DistancesDistr')
           allocate(self%dists(size(self%centers, dim = 2)), source = 0.)
           do i = 1, size(self%centers, dim = 2)
               self%dists(i) =  pixels_dist(self%centers(:,i), self%centers(:,:), 'min', mask=mask) !Use all the atoms
               mask(:) = .true. ! restore
               !Discard outliers
               if(self%dists(i)*self%smpd > 3.*self%theoretical_radius ) then
                   self%dists(i) = 0.
                   n_discard = n_discard + 1
               else if(self%dists(i)*self%smpd < 1.5*self%theoretical_radius ) then
                   print *, 'atom ',i, 'dist(i)', self%dists(i), 'INSERT THRESHOLD FOR MIN DIST BETWEEN ATOMS'
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
           write(unit = 15, fmt = '(a,f6.3,a)') 'Average dist atoms: ', self%avg_dist_atoms*self%smpd, ' A'
           write(unit = 15, fmt = '(a,f6.3,a)') 'StDev   dist atoms: ', stdev*self%smpd, ' A'
           write(unit = 15, fmt = '(a,f6.3,a)') 'Median  dist atoms: ', med*self%smpd, ' A'
           !write on pdb file
           call self%centers_pdb%new(self%n_cc, dummy=.true.)
           do i = 1,self%n_cc
               call self%centers_pdb%set_element(i,self%element)
               call self%centers_pdb%set_name(i,self%atom_name)
               call self%centers_pdb%set_beta(i,self%dists(i)*self%smpd)
           enddo
           call self%centers_pdb%writepdb('DistancesDistr')
       endif
       close(15)
   end subroutine distances_distribution

   subroutine aspect_ratios_estimation(self, print_ar)
       use gnufor2
       class(nanoparticle), intent(inout) :: self
       logical, optional,   intent(in)    :: print_ar !print longest/shortest dim and ratio
       real, pointer     :: rmat_cc(:,:,:)
       real, allocatable :: longest_dist(:)
       integer       :: label
       real :: avg_diameter, median_diameter, min_diameter, max_diameter, stdev_diameter
       call self%img_cc%get_rmat_ptr(rmat_cc)
       self%n_cc = int(maxval(rmat_cc))
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
       if(DEBUG_HERE) then
           write(logfhandle,*) 'minimum  value diameter ', min_diameter, 'A'
           write(logfhandle,*) 'maximum  value diameter ', max_diameter, 'A'
           write(logfhandle,*) 'median   value diameter ', median_diameter, 'A'
           write(logfhandle,*) 'average  value diameter ', avg_diameter, 'A'
           write(logfhandle,*) 'stdev    value diameter ', stdev_diameter, 'A'
       endif
       call hist(self%ratios, 20)
       ! To dump some of the analysis on aspect ratios on file compatible
       ! with Matlab.
        ! open(123, file='AspectRatios')
        ! write (123,*) 'AR=[...'
        ! do label = 1, size(self%ratios)
        !     write (123,'(A)', advance='no') trim(real2str(self%ratios(label)))
        !     if(label < size(self%ratios)) write (123,'(A)', advance='no') ', '
        ! end do
        ! write (123,*) '];'
        ! close(123)
        deallocate(longest_dist)
        write(logfhandle,*)'**aspect ratios calculations completed'
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
       imat_cc = int(self%img_cc%get_rmat())
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
            if(DEBUG_HERE) write(logfhandle,*) 'cc ', label, 'LONGEST DIST = 0'
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
    subroutine atom_intensity_stats(self,mmax_intensity)
        class(nanoparticle), intent(inout) :: self
        real, optional,      intent(inout) :: mmax_intensity(self%n_cc)
        logical, allocatable :: mask(:,:,:)
        real,    pointer     :: rmat(:,:,:)
        real,    pointer     :: rmat_cc(:,:,:)
        integer :: n_atom
        integer :: i, j, k
        real    :: max_intensity(self%n_cc), avg_intensity(self%n_cc), stdev_intensity(self%n_cc), radii(self%n_cc)
        real    :: avg_int, stdev_int, max_int
        real    :: m(3) ! center of mass of the nanoparticle
        write(logfhandle,*)'**atoms intensity statistics calculations init'
        call self%img%get_rmat_ptr(rmat)
        call self%img_cc%get_rmat_ptr(rmat_cc)
        allocate(mask(self%ldim(1),self%ldim(2),self%ldim(3)), source = .false.)
        self%n_cc = nint(maxval(rmat_cc))
        m = self%nanopart_masscen()
        ! initialise
        max_intensity(:)   = 0.
        avg_intensity(:)   = 0.
        stdev_intensity(:) = 0.
        ! Write on a file
        open(13, file = trim(self%output_dir)//'/IntensityStats.txt')
        do n_atom = 1, self%n_cc
            write(unit = 13, fmt = '(a,i3)') 'ATOM ', n_atom
            where(abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) - real(n_atom)) < TINY) mask = .true.
            max_intensity(n_atom) = maxval(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
            avg_intensity(n_atom) = sum(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
            avg_intensity(n_atom) = avg_intensity(n_atom)/real(count(mask))
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
            write(unit = 13, fmt = '(a,f9.5,a,f9.5,a,f9.5)') 'maxval ', max_intensity(n_atom), '   avg ', avg_intensity(n_atom), '   stdev ', stdev_intensity(n_atom)
            mask = .false. !Reset
        enddo
        if(present(mmax_intensity)) mmax_intensity = max_intensity
        avg_int   = sum(avg_intensity)/real(self%n_cc)
        max_int   = maxval(max_intensity)
        stdev_int = 0.
        do n_atom= 1, self%n_cc
            stdev_int = stdev_int + (avg_intensity(n_atom)-avg_int)**2
        enddo
        stdev_int = sqrt(stdev_int/real(self%n_cc-1))
        write(unit = 13, fmt = '(a,f9.5,a,f9.5,a,f9.5)') 'maxval_general ', max_int, '   avg_general ', avg_int, '   stdev_general ', stdev_int
        close(13)
        ! Write on a Matlab compatible file the output
        ! to enable scatterplots
        open(119, file='RadialPosCCs')
        write (119,*) 'rccs=[...'
        do n_atom = 1, self%n_cc
            write (119,'(A)', advance='no') trim(real2str(radii(n_atom)))
            if(n_atom < self%n_cc) write (119,'(A)', advance='no') ', '
        end do
        write (119,*) '];'
        close(119)
        open(119, file='MaxIntensityCCs')
        write (119,*) 'mccs=[...'
        do n_atom = 1, self%n_cc
            write (119,'(A)', advance='no') trim(real2str(max_intensity(n_atom)))
            if(n_atom < self%n_cc) write (119,'(A)', advance='no') ', '
        end do
        write (119,*) '];'
        open(119, file='AvgIntensityCCs')
        write (119,*) 'accs=[...'
        do n_atom = 1, self%n_cc
            write (119,'(A)', advance='no') trim(real2str(avg_intensity(n_atom)))
            if(n_atom < self%n_cc) write (119,'(A)', advance='no') ', '
        end do
        write (119,*) '];'
        close(119)
        write(logfhandle,*)'**atoms intensity statistics calculations completed'
    end subroutine atom_intensity_stats

    ! This subroutine clusters the atoms with respect to the maximum intensity
    ! using kmean algorithm for 2 classes. The initial guess fo the centers
    ! is intentionally biased. It supposes there are two distinguished classes
    ! with different avgs (proved with simulated data).
    subroutine cluster_atom_intensity(self, max_intensity)
        use gnufor2
        use simple_nanoML, only : nanoML
        use simple_math  ! for norm calculation
        class(nanoparticle), intent(inout) :: self
        real,                intent(inout) :: max_intensity(:)
        integer, parameter :: MAX_IT = 50 ! maximum number of iterations for
        type(image)        :: class1, class2
        real, pointer      :: rmat1(:,:,:), rmat2(:,:,:), rmat_cc(:,:,:)
        type(nanoML)       :: emfit
        real    :: centers_kmeans(2) !output of k-means
        real    :: avgs(2),  vars(2), gammas(self%n_cc,2) !output of ML
        integer :: i, cnt1, cnt2
        max_intensity = max_intensity/maxval(max_intensity)*10. !normalise with maxval 10
        call hist(max_intensity, 20)
        write(logfhandle,*) '****clustering wrt maximum intensity, init'
        ! Report clusters on images in dedicated directory
        call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: cluster_atom_intensity, simple_chdir1; ")
        call simple_mkdir(trim(self%output_dir)//'/ClusterAtomsIntensities',errmsg="simple_nanoparticles :: cluster_atom_intensity, simple_mkdir; ")
        call simple_chdir(trim(self%output_dir)//'/ClusterAtomsIntensities',errmsg="simple_nanoparticles :: cluster_atom_intensity, simple_chdir; ")
        call class1%new(self%ldim, self%smpd)
        call class2%new(self%ldim, self%smpd)
        call class1%get_rmat_ptr(rmat1)
        call class2%get_rmat_ptr(rmat2)
        call self%img_cc%get_rmat_ptr(rmat_cc)
        ! kmeans
        call emfit%kmeans_biased2classes(max_intensity, centers_kmeans)
        print *, 'centers_kmeans', centers_kmeans
        ! ML, fit
        call emfit%new(self%n_cc,2)
        call emfit%set_data(max_intensity)
        call emfit%fit(MAX_IT,centers_kmeans)
        avgs = emfit%get_avgs()
        vars = emfit%get_vars()
        gammas  = emfit%get_gammas()
        write(*,*) 'AVG/VAR 1:', avgs(1), vars(1)
        write(*,*) 'AVG/VAR 2:', avgs(2), vars(2)
        cnt2 = 0
        do i = 1, self%n_cc
            if( (avgs(1) - max_intensity(i))**2. < (avgs(2) - max_intensity(i))**2. ) then
                cnt2 = cnt2 + 1
                print *, 'connected component #', i, 'belongs to class 1 with probability', max(gammas(i,1),gammas(i,2))
            endif
        enddo
        ! assign
        ! avg shoudl be the output of ML
        ! where( (avg(1) - max_intensity)**2. < (avg(2) - max_intensity)**2. )
        !     classes = val1
        ! elsewhere
        !     classes = val2
        ! end where
        cnt1 = count((avgs(1) - max_intensity)**2. <  (avgs(2) - max_intensity)**2. )
        cnt2 = count((avgs(1) - max_intensity)**2. >= (avgs(2) - max_intensity)**2. )
        do i = 1, self%n_cc
            if( (avgs(1) - max_intensity(i))**2. < (avgs(2) - max_intensity(i))**2. ) then
                where(abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))-real(i)) < TINY) rmat1(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 1.
            else
                where(abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))-real(i)) < TINY) rmat2(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 1.
            endif
        enddo
        call class1%write('Class1.mrc')
        call class2%write('Class2.mrc')
        print *, 'Class1 contains ', cnt1 ,' atoms'
        print *, 'Class2 contains ', cnt2 ,' atoms'
        print *, 'Total ', cnt1+cnt2, ' atoms'
        call class1%kill
        call class2%kill
        ! come back to root directory
        call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: cluster_atom_intensity, simple_chdir; ")
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
                if( (cen1-max_intensity(i))**2. < (cen2-max_intensity(i))**2. )then
                    cnt1 = cnt1 + 1 ! number of elements in cluster 1
                    sum1 = sum1 + max_intensity(i)
                endif
            end do
            cnt2 = self%n_cc - cnt1       ! number of elements in cluster 2
            sum2 = sum(max_intensity)- sum1
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
        integer :: k, i
        integer, allocatable :: sz(:)
        if(allocated(self%ang_var)) deallocate(self%ang_var)
           allocate (self%ang_var(self%n_cc), source = 0.)
        sz = self%img_cc%size_connected_comps()
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
            if(DEBUG_HERE) write(logfhandle,*) 'ATOM ', k, 'angle between direction longest dim and vec [0,0,1] ', self%ang_var(k)
        enddo
        ! Matlab compatible File
        ! open(129, file='AnglesLongestDims')
        ! write (129,*) 'ang=[...'
        ! do k = 1, self%n_cc
        !     write (129,'(A)', advance='no') trim(int2str(k))
        !     write (129,'(A)', advance='no') ', '
        !     write (129,'(A)', advance='no') trim(real2str(self%ang_var(k)))
        !     write (129,'(A)')'; ...'
        ! end do
        ! write (129,*) '];'
        ! close(129)
    end subroutine search_polarization

    ! Affinity propagation clustering based on agles wrt vec [0,0,1].
    subroutine affprop_cluster_ang(self)
        use simple_aff_prop
        class(nanoparticle), intent(inout) :: self
        type(aff_prop)           :: ap_nano
        type(image), allocatable :: img_clusters(:)
        real, pointer            :: rmat_cc(:,:,:)
        real,    allocatable     :: simmat(:,:)
        real,    allocatable     :: avg_within(:), stdev_within(:)
        integer, allocatable     :: centers_ap(:), labels_ap(:)
        integer, allocatable     :: imat_onecls(:,:,:,:)
        integer                  :: i, j, ncls, nerr
        integer                  :: dim
        real                     :: simsum
        real                     :: avg, stdev
        integer, allocatable     :: cnt(:)
        dim = size(self%ang_var)!self%n_cc
        call self%centers_pdb%new(dim, dummy = .true. )
        allocate(simmat(dim, dim), source = 0.)
        write(logfhandle,*) '****clustering wrt direction longest dim, init'
        do i=1,dim-1
            do j=i+1,dim
                simmat(i,j) = -sqrt((self%ang_var(i)-self%ang_var(j))**2)
                simmat(j,i) = simmat(i,j)
            end do
        end do
        call ap_nano%new(dim, simmat)
        call ap_nano%propagate(centers_ap, labels_ap, simsum)
        ncls = size(centers_ap)
        ! Report clusters on images in dedicated directory
        call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: affprop_cluster_ang, simple_chdir1; ")
        call simple_mkdir(trim(self%output_dir)//'/ClusterAngLongestDim',errmsg="simple_nanoparticles :: affprop_cluster_ang, simple_mkdir; ")
        call simple_chdir(trim(self%output_dir)//'/ClusterAngLongestDim',errmsg="simple_nanoparticles :: affprop_cluster_ang, simple_chdir; ")
        call self%img_cc%get_rmat_ptr(rmat_cc)
        allocate(imat_onecls(self%ldim(1),self%ldim(2),self%ldim(3), ncls), source = 0)
        allocate(img_clusters(ncls))
        do i = 1, dim
            do j = 1, ncls
                call img_clusters(j)%new(self%ldim,self%smpd)
                if(labels_ap(i) == j) then
                    where(abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))-real(i))<TINY) imat_onecls(:,:,:,j) = 1
                endif
            enddo
        enddo
        do j = 1, ncls
            call img_clusters(j)%set_rmat(real(imat_onecls(:,:,:,j)))
            call img_clusters(j)%write(int2str(j)//'ClusterAng.mrc')
            call img_clusters(j)%kill
        enddo
        if(allocated(img_clusters)) deallocate(img_clusters)
        if(allocated(imat_onecls))  deallocate(imat_onecls)
        open(111, file='ClusterAngLongestDim')
        write(unit = 111,fmt ='(a,i2)') 'NR OF CLUSTERS FOUND ANG:', ncls
        do i = 1, dim
            call self%centers_pdb%set_coord(i,(self%centers(:,i)-1.)*self%smpd)
            call self%centers_pdb%set_element(i,self%element)
            call self%centers_pdb%set_beta(i, real(labels_ap(i)))
        enddo
        write(unit = 111,fmt ='(a)') 'CENTERS ANG'
        do i=1,ncls
            write(unit = 111,fmt ='(f7.3)') self%ang_var(centers_ap(i))
        end do
        !standard deviation within each center
        allocate(  avg_within(ncls), source = 0.)
        allocate(stdev_within(ncls), source = 0.)
        allocate(cnt(ncls), source = 0)
        do i = 1, ncls
            do j = 1, dim
                if(labels_ap(j) == i) then
                    cnt(i) = cnt(i) + 1 !cnt is how many atoms I have in each class
                    avg_within(i) = avg_within(i) + self%ang_var(j)
            endif
            enddo
        enddo
        avg_within = avg_within/cnt
        write(unit = 111,fmt ='(a,i3)')     'number of atoms in each class = ', cnt(1)
        do i =2, ncls
            write(unit = 111,fmt ='(a,i3)') '                                ', cnt(i)
        enddo
        write(unit = 111,fmt ='(a,i3)')     'total atoms =                   ', sum(cnt)
        do i = 1, ncls
            do j = 1, dim
                if(labels_ap(j) == i) stdev_within(i) = stdev_within(i) + (self%ang_var(j)-avg_within(i))**2
            enddo
        enddo
        where (cnt>1)
            stdev_within = stdev_within/(real(cnt)-1.)
        elsewhere
            stdev_within = 0.
        endwhere
        stdev_within = sqrt(stdev_within)
        write(unit = 111,fmt ='(a,f5.3)')  'stdev_within = ', stdev_within(1)
        do i = 2, ncls
            write(unit = 111,fmt ='(a,f5.3)')   '               ', stdev_within(i)
        enddo
        !standard deviation among different centers
        avg = 0.
        do i = 1, ncls
            avg = avg + self%ang_var(centers_ap(i))
        enddo
        avg = avg/real(ncls)
        stdev = 0.
        do i = 1, ncls
            stdev = stdev + (self%ang_var(centers_ap(i)) - avg)**2
        enddo
        stdev = stdev/(real(ncls)-1.)
        stdev = sqrt(stdev)
        write(unit = 111,fmt ='(a,f7.3)') 'standard deviation among centers: ', stdev
        call self%centers_pdb%writepdb('ClusterAngLongestDim')
        close(111)
        write(logfhandle,*) '****clustering wrt direction longest dim, completed'
        deallocate(simmat, labels_ap, centers_ap)
        if(allocated(cnt))          deallocate(cnt)
        if(allocated(avg_within))   deallocate(avg_within)
        if(allocated(stdev_within)) deallocate(stdev_within)
    end subroutine affprop_cluster_ang

    !Affinity propagation clustering based on aspect ratios
    subroutine affprop_cluster_ar(self)
        use simple_aff_prop
        class(nanoparticle), intent(inout) :: self
        type(aff_prop)           :: ap_nano
        type(image), allocatable :: img_clusters(:)
        real, pointer            :: rmat_cc(:,:,:)
        real,    allocatable     :: simmat(:,:)
        integer, allocatable     :: imat_onecls(:,:,:,:)
        integer, allocatable     :: centers_ap(:), labels_ap(:)
        real                     :: simsum
        integer                  :: i, j, ncls, nerr
        integer                  :: dim
        dim = size(self%ratios) !self%n_cc
        call self%centers_pdb%new(dim, dummy = .true. )
        allocate(simmat(dim, dim), source = 0.)
        write(logfhandle,*) '****clustering wrt aspect ratios, init'
        do i=1,dim-1
            do j=i+1,dim
                simmat(i,j) = -sqrt((self%ratios(i)-self%ratios(j))**2)
                simmat(j,i) = simmat(i,j)
            end do
        end do
        call ap_nano%new(dim, simmat)
        call ap_nano%propagate(centers_ap, labels_ap, simsum)
        ncls = size(centers_ap)
        ! Report clusters on images in dedicated directory
        call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: affprop_cluster_ar, simple_chdir1; ")
        call simple_mkdir(trim(self%output_dir)//'/ClusterAspectRatio',errmsg="simple_nanoparticles :: affprop_cluster_ar, simple_mkdir; ")
        call simple_chdir(trim(self%output_dir)//'/ClusterAspectRatio',errmsg="simple_nanoparticles :: affprop_cluster_ar, simple_chdir; ")
        call self%img_cc%get_rmat_ptr(rmat_cc)
        allocate(imat_onecls(self%ldim(1),self%ldim(2),self%ldim(3), ncls), source = 0)
        allocate(img_clusters(ncls))
        do i = 1, dim
            do j = 1, ncls
                call img_clusters(j)%new(self%ldim,self%smpd)
                if(labels_ap(i) == j) then
                        where(abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))-real(i))<TINY) imat_onecls(:,:,:,j) = 1
                endif
            enddo
        enddo
        do j = 1, ncls
            call img_clusters(j)%set_rmat(real(imat_onecls(:,:,:,j)))
            call img_clusters(j)%write(int2str(j)//'ClusterAr.mrc')
            call img_clusters(j)%kill
        enddo
        if(allocated(img_clusters)) deallocate(img_clusters)
        if(allocated(imat_onecls))  deallocate(imat_onecls)
        open(119, file='ClusterAspectRatio')
        write(unit = 119,fmt ='(a,i2)') 'NR OF CLUSTERS FOUND AR:', ncls
        do i = 1, dim
            call self%centers_pdb%set_coord(i,(self%centers(:,i)-1.)*self%smpd)
            call self%centers_pdb%set_element(i,self%element)
            call self%centers_pdb%set_beta(i, real(labels_ap(i)))
        enddo
        write(unit = 119, fmt = '(a)') 'CENTERS AR'
        do i=1,size(centers_ap)
            write(unit = 119, fmt = '(f6.3)') self%ratios(centers_ap(i))
        end do
        call self%centers_pdb%writepdb('ClustersAspectRatio')
        close(119)
        write(logfhandle,*) '****clustering wrt aspect ratios, completed'
        deallocate(simmat, labels_ap, centers_ap)
    end subroutine affprop_cluster_ar

    ! Affinity propagation clustering based on teh distribuition of
    ! the atom distances within the nanoparticle.
    ! I want to study the dustribuition of the atom distances within the nanoparticle.
    ! For example, I would expect to have atoms close to eachother in the center
    ! and far apart in the surface.
    subroutine affprop_cluster_dist_distr(self)
        use simple_aff_prop
        class(nanoparticle), intent(inout) :: self
        type(aff_prop)       :: ap_nano
        type(image), allocatable :: img_clusters(:)
        real, pointer            :: rmat_cc(:,:,:)
        real, allocatable        :: simmat(:,:)
        real                     :: simsum
        integer, allocatable     :: centers_ap(:), labels_ap(:)
        integer, allocatable     :: centers_merged(:), labels_merged(:)
        integer                  :: i, j, ncls, nerr
        integer                  :: dim
        integer, allocatable     :: imat_onecls(:,:,:,:)
        dim = size(self%dists)
        call self%centers_pdb%new(dim, dummy = .true. )
        allocate(simmat(dim, dim), source = 0.)
        write(logfhandle,*) '****clustering wrt distances distribution, init'
        do i=1,dim-1
            do j=i+1,dim
                simmat(i,j) = -sqrt((self%dists(i)-self%dists(j))**2)
                simmat(j,i) = simmat(i,j)
            end do
        end do
        call ap_nano%new(dim, simmat)
        call ap_nano%propagate(centers_ap, labels_ap, simsum)
        ncls = size(centers_ap)
        call self%img_cc%get_rmat_ptr(rmat_cc)
        allocate(imat_onecls(self%ldim(1),self%ldim(2),self%ldim(3), ncls), source = 0)
        allocate(img_clusters(ncls))
        do i = 1, dim
            do j = 1, ncls
                call img_clusters(j)%new(self%ldim,self%smpd)
                if(labels_ap(i) == j) then
                    where(abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))-real(i))<TINY) imat_onecls(:,:,:,j) = 1
                endif
            enddo
        enddo
        do j = 1, ncls
            call img_clusters(j)%set_rmat(real(imat_onecls(:,:,:,j)))
            call img_clusters(j)%write(int2str(j)//'ClusterDistDistr.mrc')
            call img_clusters(j)%kill
        enddo
        if(allocated(img_clusters)) deallocate(img_clusters)
        if(allocated(imat_onecls))  deallocate(imat_onecls)
        open(125, file = 'ClusterDistDistr')
        write(unit = 125, fmt = '(a,i2)') 'NR OF CLUSTERS FOUND DISTS DISTR:', ncls
        do i = 1, dim
            call self%centers_pdb%set_coord(i,(self%centers(:,i)-1.)*self%smpd)
            call self%centers_pdb%set_element(i,self%element)
            call self%centers_pdb%set_beta(i, real(labels_ap(i)))
        enddo
        write(unit = 125, fmt = '(a)') 'CENTERS DISTS DISTR'
        do i=1,ncls
            write(unit = 125, fmt = '(f6.3)') self%dists(centers_ap(i))
        end do
        call self%centers_pdb%writepdb('ClusterDistDistr')
        close(125)
        write(logfhandle,*) '****clustering wrt distances distribution, completed'
        deallocate(simmat, labels_ap, centers_ap)
    end subroutine affprop_cluster_dist_distr

    subroutine cluster(self, biatomic)
        class(nanoparticle), intent(inout) :: self
        character(len=3),       intent(in) :: biatomic
        real,    allocatable :: max_intensity(:)
        integer, allocatable :: classes(:)
        real, pointer :: rmat_cc(:,:,:)
        call self%img_cc%get_rmat_ptr(rmat_cc)
        self%n_cc = nint(maxval(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))))
        allocate(max_intensity(self%n_cc), classes(self%n_cc))
        rmat_cc => null()
        ! PREPARING FOR CLUSTERING
        ! Aspect ratios calculations
        call self%aspect_ratios_estimation(print_ar=.false.)
        ! Polarization search
        call self%search_polarization()
        ! CLUSTERING
        ! clustering wrt to angle longest_dim-vector z=[0,0,1]
        call self%affprop_cluster_ang()
        ! clustering wrt aspect ratio
        call self%affprop_cluster_ar()
        ! clustering wrt interatomic distances distribution
        call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: cluster, simple_chdir1; ")
        call simple_mkdir(trim(self%output_dir)//'/ClusterDistDistr',errmsg="simple_nanoparticles :: cluster, simple_mkdir; ")
        call simple_chdir(trim(self%output_dir)//'/ClusterDistDistr',errmsg="simple_nanoparticles :: cluster, simple_chdir; ")
        ! need to recalculate self%dists
        call self%distances_distribution() ! calc distances distrib across the whole nano
        call self%affprop_cluster_dist_distr()
        ! come back to root directory
        call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: cluster, simple_chdir1; ")
        ! performs clustering with respect to max intensity in the
        ! gray values of each atoms, with 2 fixed classes by using
        ! a Maximum Likelihood algorithm initialised with kmeans,
        ! only of the nanoparticle map is heterogeneous.
        if(biatomic .eq. 'yes') then
        ! Kmeans clustering (2 classes) wrt atom intensities
            call self%atom_intensity_stats(max_intensity)
            call self%cluster_atom_intensity(max_intensity)
        endif
        call self%kill
    end subroutine cluster

    ! This subrotuine indentifies the 'rows' of atoms in the z direction.
    ! The idea is to look at the projection of the nanoparticle on the
    ! xy plane and identify all the atoms whose center has (almost, see
    ! MAX_DIST_CENTERS) the same x and y coords.
    ! The inputs are:
    ! -) centers, coordinates of the centers of mass of the atoms;
    ! For how it is built, self%centers(:3,i) contains the coords of the
    ! center of mass of the i-th cc.
    subroutine make_soft_mask(self) !change the name
        class(nanoparticle), intent(inout) :: self
        call self%img_bin%grow_bins(nint(0.5*self%avg_dist_atoms)+1)
        call self%img_bin%cos_edge(SOFT_EDGE)
        call self%img_bin%write(trim(self%fbody)//'SoftMask.mrc')
    end subroutine make_soft_mask

    subroutine atoms_rmsd(nano1,nano2,r)
        use simple_atoms, only : atoms
        use gnufor2
        class(nanoparticle), intent(inout) :: nano1, nano2 !nanoparticles to compare
        real, optional,      intent(out)   :: r   !rmsd calculated
        logical, allocatable :: mask(:)
        real,    allocatable :: dist(:), dist_sq(:), dist_no_zero(:), dist_close(:)
        real,    parameter   :: CLOSE_THRESH = 1. ! 1 Amstrong
        integer :: location(1)
        integer :: i, j
        integer :: N_min !min{nb atoms in nano1, nb atoms in nano2}
        integer :: N_max !max{nb atoms in nano1, nb atoms in nano2}
        integer :: cnt,cnt2,cnt3
        real    :: sum, rmsd
        real    :: coord(3)
        real    :: avg, stdev, m(3), tmp_max, d !for statistics calculation
        type(atoms) :: centers_coupled1, centers_coupled2 !visualization purposes
        type(atoms) :: centers_close1, centers_close2
        type(atoms) :: couples1
        open(121, file='CompareAtomicModels')
        write(unit = 121, fmt = '(a)') '>>>>>>>>>   COMPARE NANO   >>>>>>>>>'
        write(unit = 121, fmt = '(a)') ''
        write(unit = 121, fmt = '(a)') 'Comparing atomic models of particles'
        write(unit = 121, fmt = '(a,a)') trim(nano1%fbody), ' ---> vol1'
        write(unit = 121, fmt = '(a)') 'and'
        write(unit = 121, fmt = '(a,a)') trim(nano2%fbody), ' ---> vol2'
        write(unit = 121, fmt = '(a)')  '>>>>>>>>>VOLUME COMPARISION>>>>>>>>'
        ! If they don't have the same nb of atoms
        if(size(nano1%centers, dim = 2) <= size(nano2%centers, dim = 2)) then
            N_min = size(nano1%centers, dim = 2)
            N_max = size(nano2%centers, dim = 2)
            write(unit = 121, fmt = '(a,i3)') 'NB atoms in vol1   ', N_min
            write(unit = 121, fmt = '(a,i3)') 'NB atoms in vol2   ', N_max
            write(unit = 121, fmt = '(a,i3,a)') '--->', abs(N_max-N_min), ' atoms do NOT correspond'
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
                if(DEBUG_HERE) then
                     write(logfhandle,*) 'ATOM', i,'coords: ', nano2%centers(:,i), 'coupled with '
                     write(logfhandle,*) '    ',location, 'coordinates: ', nano1%centers(:,location(1)), 'DIST^2= ', dist_sq(i), 'DIST = ', dist(i)
                endif
            enddo
        else
            N_min = size(nano2%centers, dim = 2)
            N_max = size(nano1%centers, dim = 2)
            write(unit = 121, fmt = '(a,i3)') 'NB atoms in vol1   ', N_max
            write(unit = 121, fmt = '(a,i3)') 'NB atoms in vol2   ', N_min
            write(unit = 121, fmt = '(a,i3,a)') '--->', abs(N_max-N_min), ' atoms do NOT correspond'
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
                if(DEBUG_HERE) then
                    write(logfhandle,*) 'ATOM', i,'coordinates: ', nano1%centers(:,i), 'coupled with '
                    write(logfhandle,*) '    ',location, 'coordinates: ', nano2%centers(:,location(1)), 'DIST^2= ', dist_sq(i), 'DIST = ', dist(i)
                endif
            enddo
        endif
        write(unit = 121, fmt = '(i3,a,i2,a)')      cnt3,' atoms correspond within        1 A. (', cnt3*100/N_min, '% of the atoms )'
        write(unit = 121, fmt = '(i3,a,f3.1,a,i2,a)')  cnt2,' atoms correspond within  1 - ',2.*nano2%theoretical_radius,' A. (', cnt2*100/N_min, '% of the atoms )'
        write(unit = 121, fmt = '(i3,a,f3.1,a,i2,a)')  cnt-(N_max-N_min),' atoms have error bigger than ',2.*nano2%theoretical_radius,' A. (',(cnt-N_max+N_min)*100/N_min, '% of the atoms )' !remove the extra atoms
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
        write(unit = 121, fmt = '(a)')       '--->    IN VOL1    <---'
        write(unit = 121, fmt = '(a,f6.3,a)')'AVG     DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', avg, ' A'
        write(unit = 121, fmt = '(a,f6.3,a)')'STDEV   DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', stdev, ' A'
        write(unit = 121, fmt = '(a,f6.3,a)')'MAX     DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', tmp_max, ' A'
        tmp_max = 0. ! reset
        do i = 1, size(nano1%centers, dim = 2)
            coord(:) = (nano1%centers(:,i)-1.)*nano1%smpd
                d =  euclid(coord,m)
                if(d > tmp_max) tmp_max = d
        enddo
        write(unit = 121, fmt = '(a,f6.3,a)')'MAX     DIST ATOMS IN VOL1 3D RECONSTRUCT. TO THE CENTER: ', tmp_max, ' A'
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
        write(unit = 121, fmt = '(a)')       '--->    IN VOL2    <---'
        write(unit = 121, fmt = '(a,f6.3,a)')'AVG     DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', avg, ' A'
        write(unit = 121, fmt = '(a,f6.3,a)')'STDEV   DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', stdev, ' A'
        write(unit = 121, fmt = '(a,f6.3,a)')'MAX     DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', tmp_max, ' A'
        tmp_max = 0. ! reset
        do i = 1, size(nano2%centers, dim = 2)
            coord(:) = (nano2%centers(:,i)-1.)*nano2%smpd
                d =  euclid(coord,m)
                if(d > tmp_max) tmp_max = d
        enddo
        write(unit = 121, fmt = '(a,f6.3,a)')'MAX     DIST ATOMS IN VOL2 3D RECONSTRUCT. TO THE CENTER: ', tmp_max, ' A'
        ! kill atoms instances
        call centers_close1%kill
        call centers_close2%kill
        call centers_coupled1%kill
        call centers_coupled2%kill
        call couples1%kill
        !RMSD
        rmsd = sqrt(sum(dist_sq)/real(count(dist_sq > TINY)))
        write(unit = 121, fmt = '(a,f6.3,a)') 'RMSD CALCULATED CONSIDERING ALL ATOMS = ', rmsd*nano1%smpd, ' A'
        write(unit = 121, fmt = '(a,f6.3,a)') 'RMSD ATOMS THAT CORRESPOND WITHIN 1 A = ', (sqrt(sum(dist_close)/real(count(dist_close > TINY))))*nano1%smpd, ' A'
        dist_no_zero = pack(dist, dist>TINY)
        dist_no_zero = dist_no_zero*nano1%smpd ! report distances in Amstrongs
        call hist(dist_no_zero, 50)
        !For SCV files
        open(117, file='RMSDhist')
        write (117,*) 'r'
        do i = 1, size(dist_no_zero)
            write (117,'(A)', advance='yes') trim(real2str(dist_no_zero(i)))
        end do
        close(117)
        if(present(r)) r=rmsd
        close(121)
        deallocate(dist, dist_sq, dist_no_zero, mask)
    end subroutine atoms_rmsd

    ! Check for symmetry at different radii.
    ! radii: 5, 7, 9, 12 A
    ! The idea is to identify the atomic positions
    ! contained in a certain radius, generate a distribution
    ! based on that and check for symmetry.
    subroutine identify_atomic_pos(self)
        class(nanoparticle), intent(inout) :: self
        ! Phase correlations approach
        call self%phasecorrelation_nano_gaussian()
        ! Nanoparticle binarization
        call self%binarize()
        ! Outliers discarding
        call self%discard_outliers()
        ! Validation of the selected atomic positions
         call self%validate_atomic_positions()
    end subroutine identify_atomic_pos

    subroutine keep_atomic_pos_at_radius(self, radius, element, fname)
      class(nanoparticle), intent(inout) :: self
      real,                intent(in)    :: radius
      character(len=2),    intent(in)    :: element
      character(len=*),    intent(in)    :: fname !name of the pdb file where to write the coords
      character(len=4) :: atom_name
      integer     :: n_cc,cnt
      real        :: m(3) !coords of the center of mass of the nano
      real        :: d    !distance of the atom from the center of mass
      type(atoms) :: radial_atom
      ! Make subroutine for set_atom_name()
      select case(element)
      case('pt')
          atom_name     = ' pt '
          ! thoretical radius is already set
      case('pd')
          atom_name     = ' pd '
      case('fe')
          atom_name     = 'fe'
      case('au')
          self%atom_name     = ' au '
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

    subroutine identify_atom_columns(self)
        class(nanoparticle), intent(inout) :: self
        real, pointer      :: rmat_cc(:,:,:), rmat_col(:,:,:)
        type(image)        :: img_col
        integer, parameter :: N_DISCRET = 1000
        integer, parameter :: N_LINES   = 6    ! number of lines to consider for each atom (exhaustive search) (12???)
        integer :: n_cc, loc(1), i, j
        integer :: x, y, z, t
        integer :: cnt_intersect
        logical :: mask(self%n_cc),mask_line(N_DISCRET)
        logical :: flag_matrix(N_LINES,self%n_cc,self%n_cc)
        real    :: dir(3), dist
        real    :: line(3, N_DISCRET), t_vec(N_DISCRET)
        do i = 1, N_DISCRET/2
            t_vec(i) = -real(i)/10.
        enddo
        t_vec(N_DISCRET/2+1:N_DISCRET) = -t_vec(1:N_DISCRET/2)
        x = 1
        y = 2
        z = 3 ! just for notation
        call self%img_cc%get_rmat_ptr(rmat_cc)
        call img_col%new(self%ldim, self%smpd)
        call img_col%get_rmat_ptr(rmat_col)
        flag_matrix = .false.  ! initialization
        do n_cc = 1, self%n_cc ! for each atom
            mask   = .true.
            do j = 1, N_LINES        ! multiple lines for each atom
                cnt_intersect = 0  ! how many atoms does the line intersect
                ! find nearest neighbour
                dist   = pixels_dist(self%centers(:,n_cc), self%centers,'min', mask, loc)
                mask(loc(1)) = .false.                                 ! next time pick other neigh
                if(dist*self%smpd >  4.*self%theoretical_radius) cycle ! disregard if too far
                ! direction of the line
                dir    = self%centers(:,n_cc) - self%centers(:,loc(1))
                ! line
                do t = 1, N_DISCRET
                    line(x,t) = self%centers(1,n_cc) + t_vec(t)* dir(1)
                    line(y,t) = self%centers(2,n_cc) + t_vec(t)* dir(2)
                    line(z,t) = self%centers(3,n_cc) + t_vec(t)* dir(3)
                enddo
                ! calculate how many atoms does the line intersect and flag them
                do i = 1, self%n_cc
                    mask_line = .true.
                    dist   = pixels_dist(self%centers(:3,i), line(:3,:),'min', mask_line)
                    if(dist*self%smpd <= 0.9*self%theoretical_radius) then ! it intersects atoms i
                         cnt_intersect = cnt_intersect + 1
                         flag_matrix(j,n_cc,i) = .true. !flags also itself
                     endif
                enddo
            enddo
        enddo
        ! remove redundancies
        do j = 1, N_LINES
            do t = 1, N_LINES
                if(j /= t) then
                    do n_cc = 1, self%n_cc     ! fix the first row
                        do i = 1, self%n_cc   ! fix the second row
                            if(all(flag_matrix(j,n_cc,:) .eqv. flag_matrix(t,i,:))) then
                                flag_matrix(t,i,:)= .false. ! unflag them
                            endif
                        enddo
                    enddo
                else
                    do n_cc = 1, self%n_cc-1     ! fix the first row
                        do i = n_cc+1, self%n_cc   ! fix the second row
                            if(all(flag_matrix(j,n_cc,:) .eqv. flag_matrix(t,i,:))) then
                                flag_matrix(t,i,:)= .false. ! unflag them
                            endif
                        enddo
                    enddo
                endif
            enddo
        enddo
        ! generate images for visualisation
        do j = 1, N_LINES
            do n_cc = 1, self%n_cc
                ! reset
                rmat_col(:,:,:) = 0.
                cnt_intersect   = 0
                do i = 1, self%n_cc
                    if(flag_matrix(j,n_cc,i)) then
                        cnt_intersect = cnt_intersect + 1
                        where(abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) - real(i)) < TINY) rmat_col(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 1.
                    endif
                enddo
                if(cnt_intersect > 5) call img_col%write(trim(int2str(n_cc))//'ImageColNonRedundant'//trim(int2str(j))//'.mrc')
            enddo
        enddo
        call img_col%kill
    end subroutine identify_atom_columns

    subroutine identify_atom_planes(self)
        class(nanoparticle), intent(inout) :: self
        real, pointer      :: rmat_cc(:,:,:), rmat_plane(:,:,:)
        type(image)        :: img_plane
        integer, parameter :: N_DISCRET = 1000
        integer, parameter :: N_PLANES  = 1    ! number of planes to consider for each atom
        integer :: n_cc, loc(1), i, j
        integer :: t, s
        integer :: cnt_intersect, cnt
        logical :: mask(self%n_cc),mask_plane(N_DISCRET)
        logical :: flag_matrix(N_PLANES,self%n_cc,self%n_cc)
        real    :: dir_1(3), dir_2(3), vec(3), m(3), dist
        real    :: line(3,N_DISCRET)
        real, allocatable :: plane(:,:,:), points(:,:), distances_totheplane(:), radii(:)
        real    :: t_vec(N_DISCRET), s_vec(N_DISCRET), denominator
        allocate(plane(3, N_DISCRET, N_DISCRET), source = 0.)
        do i = 1, N_DISCRET/2
            t_vec(i) = -real(i)/10.
        enddo
        t_vec(N_DISCRET/2+1:N_DISCRET) = -t_vec(1:N_DISCRET/2)
        s_vec(:) = t_vec(:)
        call self%img_cc%get_rmat_ptr(rmat_cc)
        call img_plane%new(self%ldim, self%smpd)
        call img_plane%get_rmat_ptr(rmat_plane)
        flag_matrix(:,:,:) = .false.  ! initialization
        rmat_plane(:,:,:) = 0.
        write(logfhandle, *) 'planes identification'
        do n_cc = 205, 205!1, self%n_cc ! for each atom
            mask   = .true.
            do j = 1, N_PLANES
                cnt_intersect = 0  ! how many atoms does the line intersect
                ! find nearest neighbour
                dist   = pixels_dist(self%centers(:,n_cc), self%centers,'min', mask, loc)
                mask(loc(1)) = .false.                                 ! next time pick other neigh
                if(dist*self%smpd >  4.*self%theoretical_radius) cycle ! disregard if too far
                ! direction vectors of the plane
                dir_1    = self%centers(:,n_cc) - self%centers(:,loc(1))
                dist   = pixels_dist(self%centers(:,n_cc), self%centers,'min', mask, loc)
                mask(loc(1)) = .false.                                 ! next time pick other neigh
                if(dist*self%smpd >  4.*self%theoretical_radius) cycle ! disregard if too far
                dir_2    = self%centers(:,n_cc) - self%centers(:,loc(1))
                do t = 1, N_DISCRET
                    do s = 1, N_DISCRET
                        plane(1,t,s) = self%centers(1,n_cc) + t_vec(t)* dir_1(1) + s_vec(s)* dir_2(1)
                        plane(2,t,s) = self%centers(2,n_cc) + t_vec(t)* dir_1(2) + s_vec(s)* dir_2(2)
                        plane(3,t,s) = self%centers(3,n_cc) + t_vec(t)* dir_1(3) + s_vec(s)* dir_2(3)
                    enddo
                enddo
                ! calculate how many atoms does the plane intersect and flag them
                do i = 1, self%n_cc
                    mask_plane = .true.
                    do t = 1, N_DISCRET
                        do s = 1, N_DISCRET
                            dist   = euclid(self%centers(:3,i),plane(:3,t,s))
                            if(dist*self%smpd <= 0.9*self%theoretical_radius) then ! it intersects atoms i
                                 flag_matrix(j,n_cc,i) = .true. !flags also itself
                             endif
                         enddo
                    enddo
                enddo
            enddo
        enddo
        write(logfhandle, *) 'removing redundancies'
        ! remove redundancies
        do j = 1, N_PLANES
            do t = 1, N_PLANES
                if(j /= t) then
                    do n_cc = 1, self%n_cc     ! fix the first row
                        do i = 1, self%n_cc   ! fix the second row
                            if(all(flag_matrix(j,n_cc,:) .eqv. flag_matrix(t,i,:))) then
                                flag_matrix(t,i,:)= .false. ! unflag them
                            endif
                        enddo
                    enddo
                else
                    do n_cc = 1, self%n_cc-1       ! fix the first row
                        do i = n_cc+1, self%n_cc   ! fix the second row
                            if(all(flag_matrix(j,n_cc,:) .eqv. flag_matrix(t,i,:))) then
                                flag_matrix(t,i,:)= .false. ! unflag them
                            endif
                        enddo
                    enddo
                endif
            enddo
        enddo
        write(logfhandle, *) 'generating volumes for visualization'
        ! generate volume for visualisation
        do n_cc = 205, 205!1, self%n_cc
            call progress(n_cc, self%n_cc)
            ! reset
            rmat_plane(:,:,:) = 0.
            cnt_intersect   = 0
            do i = 1, self%n_cc
                if(flag_matrix(1,n_cc,i)) then
                    cnt_intersect = cnt_intersect + 1
                    where(abs(rmat_cc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) - real(i)) < TINY) rmat_plane(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 1.
                endif
            enddo
            if(cnt_intersect > 20) call img_plane%write(trim(int2str(n_cc))//'ImagePlane'//trim(int2str(j))//'.mrc')
        enddo
        ! TO MOVE
        allocate(points(3, count(flag_matrix)), source = 0.)
        m = self%nanopart_masscen()
        cnt = 0
        do n_cc = 205,205 ! 1, self%n_cc
            do j = 1, N_PLANES
                do i = 1, self%n_cc
                    if(flag_matrix(j,n_cc,i)) then
                        cnt = cnt + 1
                        points(:3,cnt) = self%centers(:3,i)-m(:)
                    endif
                enddo
            enddo
        enddo
        vec = plane_from_points(points)
        allocate(distances_totheplane(cnt), source = 0.)
        allocate(radii(cnt), source = 0.) ! which radius is the atom center belonging to
        cnt = 0
        denominator = sqrt(vec(1)*2+vec(2)*2+1.)
        ! TO PUT TOGETHER WITH THE PREVIOUS LOOP
        do n_cc = 205,205 ! 1, self%n_cc
            do j = 1, N_PLANES
                do i = 1, self%n_cc
                    if(flag_matrix(j,n_cc,i)) then
                        cnt = cnt + 1
                        points(:3,cnt) = self%centers(:3,i)-m(:)
                        ! formula for distance of a point to a plane
                        distances_totheplane(cnt) = abs(vec(1)*points(1,cnt)+vec(2)*points(2,cnt)-points(3,cnt)+vec(3))/denominator
                        radii(cnt) = euclid(self%centers(:,i), m)*self%smpd
                    endif
                enddo
            enddo
        enddo
        distances_totheplane = (distances_totheplane)*self%smpd
        open(119, file='Radii')
        write (119,*) 'r=[...'
        do i = 1, cnt
            write (119,'(A)', advance='no') trim(real2str(radii(i)))
            if(i < cnt) write (119,'(A)', advance='no') ', '
        end do
        write (119,*) '];'
        close(119)
        open(119, file='DistancesToThePlane')
        write (119,*) 'd=[...'
        do i = 1, cnt
            write (119,'(A)', advance='no') trim(real2str(distances_totheplane(i)))
            if(i < cnt) write (119,'(A)', advance='no') ', '
        end do
        write (119,*) '];'
        close(119)
        call img_plane%kill
        deallocate(plane)
    contains
        ! Find the plane that minimises the distance between
        ! a given set of points.
        ! It consists in a solution of a overdetermined system with
        ! the left pseudo inverse.
        ! SOURCE :
        ! https://stackoverflow.com/questions/1400213/3d-least-squares-plane
        ! The output plane will have cartesian equation
        ! vec(1)x + vec(2)y - z = -vec(3).
        ! FORMULA
        ! sol = inv(transpose(M)*M)*transpose(M)*b
        function plane_from_points(points) result(sol)
            real,                intent(inout) :: points(:,:) !input
            real    :: sol(3)  !vec(1)x + vec(2)y - z = -vec(3).
            real    :: M(size(points, dim = 2),3), b(size(points, dim = 2)), invM(3,size(points, dim = 2))
            real    :: prod(3,3), prod_inv(3,3), prod1(3,size(points, dim = 2))
            integer :: errflg ! if manages to find inverse matrix
            integer :: p
            integer :: N ! number of points
            if(size(points, dim=1) /=3) then
                write(logfhandle,*) 'Need to input points in 3D!; plane_from_points'
                return
            endif
            if(size(points, dim=2) < 3) then
                write(logfhandle,*) 'Not enough input points to fit a plane!; plane_from_points'
                return
            endif
            N = size(points, dim=2)
            do p = 1, N
                M(p,1) =  points(1,p)
                M(p,2) =  points(2,p)
                M(p,3) =  1.
                b(p)   =  points(3,p)
            enddo
            prod  = matmul(transpose(M),M)
            call matinv(prod,prod_inv,3,errflg)
            if( errflg /= 0 ) THROW_HARD('Couldn t find inverse matrix! ;plane_from_points')
            prod1 = matmul(prod_inv,transpose(M))
            sol   = matmul(prod1,b)
            print *, 'Solution'
            print *, sol
        end function plane_from_points
    end subroutine identify_atom_planes

    subroutine kill_nanoparticle(self)
        class(nanoparticle), intent(inout) :: self
        self%ldim(3)           = 0
        self%smpd              = 0.
        self%nanop_mass_cen(3) = 0.
        self%avg_dist_atoms    = 0.
        self%n_cc              = 0
        call self%img%kill()
        call self%img_raw%kill
        call self%img_bin%kill()
        call self%img_cc%kill()
        call self%centers_pdb%kill
        if(allocated(self%centers))          deallocate(self%centers)
        if(allocated(self%ratios))           deallocate(self%ratios)
        if(allocated(self%dists))            deallocate(self%dists)
        if(allocated(self%loc_longest_dist)) deallocate(self%loc_longest_dist)
        if(allocated(self%ang_var))          deallocate(self%ang_var)
    end subroutine kill_nanoparticle
end module simple_nanoparticles_mod

    ! !Output on file, Matlab Compatible
    ! open(119, file='RMSDhist')
    ! write (119,*) 'r=[...'
    ! do i = 1, size(dist_no_zero)
    !     write (119,'(A)', advance='no') trim(real2str(dist_no_zero(i)))
    !     if(i < size(dist_no_zero)) write (119,'(A)', advance='no') ', '
    ! end do
    ! write (119,*) '];'
    ! close(119)

    ! This subroutine performs size filtering on the connected
    ! components image of the nanoparticle. It calculates the
    ! size of all the connected components, the average and the standard
    ! deviation. It hypothesises gaussian distribution, so it discards
    ! the connected components which are outside the range [-2sigma,2sigma].
    ! subroutine size_filtering(self)
    !     class(nanoparticle), intent(inout) :: self
    !     integer, allocatable :: sz(:)
    !     real, pointer :: rmat(:,:,:)
    !     real, pointer :: rmat_cc(:,:,:)
    !     integer :: cc
    !     real    :: avg_sz
    !     real    :: stdev_sz
    !     write(logfhandle, *) '****size filtering, init'
    !     call self%img_bin%get_rmat_ptr(rmat)
    !     call self%img_cc%get_rmat_ptr(rmat_cc)
    !     sz = self%img_cc%size_connected_comps()
    !     avg_sz = real(sum(sz))/real(size(sz))
    !     stdev_sz = 0.
    !     do cc = 1, size(sz)
    !         stdev_sz = stdev_sz + (sz(cc) - avg_sz)**2
    !     enddo
    !     stdev_sz = sqrt(stdev_sz/(real(size(sz)-1)))
    !     !assuming Gaussian distrib, 95% is in [-2sigma,2sigma] and 68% is in [-sigma,sigma]
    !     !big ccs have already been removed by erosion. Now we need to remove too small
    !     !ccs, they usually represent background noise.
    !     do cc = 1, size(sz)
    !         if(sz(cc)<avg_sz-2.*stdev_sz ) then
    !             where(abs(rmat_cc-real(cc)) < TINY)
    !                 rmat_cc = 0.
    !                 rmat    = 0.
    !             endwhere
    !          endif
    !     enddo
    !     !update img_cc: re-order ccs
    !     call self%img_cc%order_cc()
    !     ! UPDATE IMG_BIN???
    !     call self%img_cc%write('AfterSizeFiltering.mrc')
    !     !update n_cc
    !     self%n_cc = nint(maxval(rmat_cc))
    !     !update centers
    !     call self%find_centers()
    !     write(logfhandle, *) '****size filtering, completed'
    !     deallocate(sz)
    ! end subroutine size_filtering
