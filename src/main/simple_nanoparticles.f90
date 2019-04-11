!USAGE 1: simple_exec prg=detect_atoms smpd=0.358 vol1='particle1.mrc'
!USAGE 2: simple_exec prg=compare_nano vol1=vol_c1_aligned2_c3axis.mrc vol2=particle1_c3.mrc
module simple_nanoparticles_mod
include 'simple_lib.f08'
use simple_image,         only : image
use simple_picker_chiara, only : pixels_dist, get_pixel_pos, polish_cc
use simple_atoms,         only : atoms

implicit none

private
public :: nanoparticle

#include "simple_local_flags.inc"

! module global constants
integer, parameter :: N_THRESH   = 20      !number of thresholds for binarization
logical, parameter :: DEBUG_HERE = .false. !for debugging purposes
real,    parameter :: VDW_RAD_PT = 1.77   !Platinum Van der Waals radius in Amstrongs

type :: nanoparticle
    private
    ! these image objects are part of the instance to avoid excessive memory re-allocations
    type(atoms) :: centers_pdb
    type(image) :: img, img_bin, img_cc, img_over_smp
    integer     :: ldim(3)           = 0
    real        :: smpd              = 0.
    real        :: nanop_mass_cen(3) = 0.!coordinates of the center of mass of the nanoparticle
    real        :: avg_dist_atoms    = 0.
    real        :: SCALE_FACTOR      = 0.!for oversampling
    integer     :: n_cc              = 0 !number of atoms (connected components)
    integer     :: dim_over(3)       = 0 !logical dimensions of the oversampled
    real,    allocatable  :: centers(:,:)
    real,    allocatable  :: ratios(:)
    real,    allocatable  :: circularity(:)
    real,    allocatable  :: ang_var(:)
    real,    allocatable  :: dists(:)
    integer, allocatable  :: loc_longest_dist(:,:)   !for indentific of the vxl that determins the longest dim of the atom
    character(len=STDLEN) :: partname = ''   !fname
    character(len=STDLEN) :: fbody    = ''   !fbody
  contains
    ! constructor
    procedure          :: new => new_nanoparticle
    ! getters/setters
    procedure          :: get_img
    procedure          :: set_img
    procedure          :: set_partname
    ! segmentation and statistics
    procedure          :: binarize => nanopart_binarization
    procedure          :: find_centers
    procedure, private :: nanopart_masscen
    procedure, private :: calc_aspect_ratio_private
    procedure          :: calc_aspect_ratio
    procedure          :: calc_circularity
    procedure          :: discard_outliers
    procedure, private :: distances_distribution
    ! clustering
    procedure          :: cluster => nanopart_cluster
    procedure, private :: affprop_cluster_ar
    procedure, private :: affprop_cluster_ang
    procedure, private :: affprop_cluster_circ
    procedure, private :: affprop_cluster_dist_distr
    procedure          :: search_polarization
    ! comparison
    procedure          :: compare_atomic_models
    ! visualization purpose
    procedure          :: print_asym_unit
    ! others
    procedure          :: make_soft_mask
    procedure, private :: over_sample
    ! kill
    procedure          :: kill => kill_nanoparticle
end type nanoparticle

contains

    !constructor
    subroutine new_nanoparticle(self, fname, sc_fac)
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        real, optional,      intent(in)    :: sc_fac
        integer :: nptcls
        real    :: ssc_fac
        call self%kill
        ssc_fac = 1.
        if(present(sc_fac)) ssc_fac = sc_fac
        self%SCALE_FACTOR = ssc_fac
        call self%set_partname(fname)
        self%fbody = get_fbody(trim(fname), trim(fname2ext(fname)))
        call find_ldim_nptcls(self%partname,  self%ldim, nptcls, self%smpd)
        call self%img%new         (self%ldim, self%smpd)
        call self%img_bin%new     (int(real(self%ldim)*self%SCALE_FACTOR), self%smpd/self%SCALE_FACTOR)
        call self%img_over_smp%new(int(real(self%ldim)*self%SCALE_FACTOR), self%smpd/self%SCALE_FACTOR)
        self%dim_over(1) = int(real(self%ldim(1))*self%SCALE_FACTOR)
        self%dim_over(2) = int(real(self%ldim(2))*self%SCALE_FACTOR)
        self%dim_over(3) = int(real(self%ldim(3))*self%SCALE_FACTOR)
        call self%img%read(fname)
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

     !set one of the images of the nanoparticle type
     subroutine set_img( self, img, which )
         class(nanoparticle), intent(inout)    :: self
         type(image),         intent(inout) :: img
         character(len=*),    intent(in)    :: which
         integer :: ldim(3)
         ldim = img%get_ldim()
         if(self%ldim(1) .ne. ldim(1) .or. self%ldim(2) .ne. ldim(2) .or. self%ldim(3) .ne. ldim(3)) THROW_HARD('Wrong dimension in input img;   set_img')
         select case(which)
         case('img')
             self%img = img
         case('img_bin')
             self%img_bin = img
         case('img_cc')
             self%img_cc = img
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

    ! This subrotuine takes in input a nanoparticle and
    ! binarizes it by thresholding. The gray level histogram is split
    ! in 20 parts, which corrispond to 20 possible threshold.
    ! Among those threshold, the selected one is the for which
    ! the size of the connected component of the correspondent
    ! thresholded nanoparticle has the maximum median value.
    ! The idea is that if the threshold is wrong, than the binarization
    ! produces a few huge ccs (connected components) and a lot of very small
    ! ones (dust). So the threshold with the maximum median value would
    ! correspond to the most consistent one, meaning that the size of the ccs
    ! would have a gaussian distribution.
    subroutine nanopart_binarization( self )
        class(nanoparticle), intent(inout) :: self
        real, allocatable :: rmat(:,:,:), rmat_t(:,:,:)
        real    :: step          !histogram disretization step
        real    :: thresh        !binarization threshold
        real    :: seleted_t(1)  !selected threshold for nanoparticle binarization
        real    :: x_thresh(N_THRESH-1), y_med(N_THRESH-1)
        integer, allocatable :: sz(:) !size of the ccs and correspondent label
        integer ::  i
        write(logfhandle,*) '****binarization, init'
        call self%over_sample()
        rmat = self%img_over_smp%get_rmat()
        allocate(rmat_t(self%dim_over(1), self%dim_over(2), self%dim_over(3)), source = 0.)
        step = maxval(rmat)/real(N_THRESH)
        do i = 1, N_THRESH-1
            call progress(i, N_THRESH-1)
            thresh = real(i)*step
            where(rmat > thresh)
                rmat_t = 1.
            elsewhere
                rmat_t = 0.
            endwhere
            call self%img_bin%set_rmat(rmat_t)
            call self%img_bin%find_connected_comps(self%img_cc)
            sz          = self%img_cc%size_connected_comps()
            x_thresh(i) = thresh
            y_med(i)    = median(real(sz))
        enddo
        seleted_t(:) = x_thresh(maxloc(y_med))
        if(DEBUG_HERE) write(logfhandle,*)  'SELECTED THRESHOLD = ', seleted_t(1)
        where(rmat > seleted_t(1))
            rmat_t = 1.
        elsewhere
            rmat_t = 0.
        endwhere
        call self%img_bin%set_rmat(rmat_t)
        deallocate(rmat, rmat_t, sz)
        write(logfhandle,*) '****binarization, completed'
    end subroutine nanopart_binarization

    ! This subroutine takes in input 2 nanoparticle and their
    ! correspondent filenames.
    ! The output is the rmsd of their atomic models.
    subroutine compare_atomic_models(nano1,nano2)
        class(nanoparticle), intent(inout) :: nano1, nano2 !nanoparticles to compare
        real, allocatable :: centers1(:,:), centers2(:,:)
        real :: ada  !minimum average dist atoms between nano1 and nano2
        real :: rmsd
        integer :: i
        call nano1%binarize()
        call nano2%binarize()
        ! Atoms centers identification
        call nano1%find_centers(centers1)
        call nano2%find_centers(centers2)
        ! Outliers discarding
        call nano1%discard_outliers()
        call nano2%discard_outliers()
        ! Ouput file
        open(121, file='CompareAtomicModels') !, position='append')
        write(unit = 121, fmt = '(a)') ''
        write(unit = 121, fmt = '(a)') 'Comparing atomic models of particles'
        write(unit = 121, fmt = '(a)') trim(nano1%fbody)
        write(unit = 121, fmt = '(a)') 'and'
        write(unit = 121, fmt = '(a)') trim(nano2%fbody)
        ! RMSD calculation
        write(logfhandle,*) 'avg dist between atoms in vol1 = ',nano1%avg_dist_atoms*nano1%smpd,'A'
        write(logfhandle,*) 'avg dist between atoms in vol2 = ',nano2%avg_dist_atoms*nano2%smpd,'A'
        ada = min(nano1%avg_dist_atoms, nano2%avg_dist_atoms)
        call atomic_position_rmsd(nano1,nano2, ada, rmsd)
        write(unit = 121, fmt = '(a,f6.3,a)') 'RMSD = ', rmsd*nano1%smpd, "A"
        write(unit = 121, fmt = '(a)') '***comparison completed'
        close(121)
        ! Circularities calculation
        call nano1%calc_circularity()
        call nano2%calc_circularity()
        ! kill
        call nano1%kill
        call nano2%kill
        deallocate(centers1,centers2)
    contains

        ! This subroutine takes in input coords of the atomic positions
        ! in the nanoparticles to compare and calculates the rmsd
        ! between them.
        ! See formula
        ! https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
        subroutine atomic_position_rmsd(nano1,nano2,ada,r)
            use simple_atoms, only : atoms
            use gnufor2
            class(nanoparticle), intent(inout) :: nano1, nano2 !nanoparticles to compare
            real,           intent(in)  :: ada !minimum average dist atoms between nano1 and nano2
            real, optional, intent(out) :: r   !rmsd calculated
            logical, allocatable :: mask(:)
            real,    allocatable :: dist(:), dist_sq(:), dist_no_zero(:)
            integer :: location(1)
            integer :: i, j
            integer :: N_min !min{nb atoms in nano1, nb atoms in nano2}
            integer :: N_max !max{nb atoms in nano1, nb atoms in nano2}
            integer :: cnt, cnt_zeros
            real    :: sum, rmsd
            type(atoms) :: centers_coupled1, centers_coupled2 !visualization purposes
            ! If they don't have the same nb of atoms
            if(size(nano1%centers, dim = 2) <= size(nano2%centers, dim = 2)) then
                N_min = size(nano1%centers, dim = 2)
                N_max = size(nano2%centers, dim = 2)
                write(logfhandle,*) 'NB atoms in vol1   ', N_min
                write(logfhandle,*) 'NB atoms in vol2   ', N_max
                write(logfhandle,*) abs(N_max-N_min), 'atoms do NOT correspond'
                call centers_coupled1%new(N_min, dummy=.true.)
                call centers_coupled2%new(N_max, dummy=.true.)
                allocate(dist(N_min), dist_sq(N_min), source = 0.)
                allocate(mask(size(nano2%centers, dim = 2)), source = .true.)
                cnt = 0
                do i = 1, N_min !compare based on centers1
                    dist(i) = pixels_dist(nano1%centers(:,i),nano2%centers(:,:),'min',mask,location)
                    if(dist(i)*nano1%smpd > VDW_RAD_PT+0.2*VDW_RAD_PT) then
                        dist(i) = 0. !it means there is no correspondent atom in the other nano
                        cnt = cnt+1  !to discard them in the rmsd calculation
                    else
                         dist_sq(i) = dist(i)**2 !formula wants them square, could improve performance here
                         if(DEBUG_HERE) then
                             write(logfhandle,*) 'ATOM', i,'coords: ', nano1%centers(:,i), 'coupled with '
                             write(logfhandle,*) '    ',location, 'coordinates: ', nano2%centers(:,location(1)), 'DIST^2= ', dist_sq(i), 'DIST = ', dist(i)
                         endif
                         call centers_coupled1%set_coord(i,(nano1%centers(:,i)-1.)*0.358/nano1%SCALE_FACTOR)
                         call centers_coupled2%set_coord(i,(nano2%centers(:,location(1))-1.)*0.358/nano2%SCALE_FACTOR)
                         mask(location(1)) = .false. ! not to consider the same atom more than once
                     endif
                enddo
            else
                N_min = size(nano2%centers, dim = 2)
                N_max = size(nano1%centers, dim = 2)
                write(logfhandle,*) 'NB atoms in vol1   ', N_max
                write(logfhandle,*) 'NB atoms in vol2   ', N_min
                write(logfhandle,*) abs(N_max-N_min), 'atoms do NOT correspond'
                allocate(dist(N_min), dist_sq(N_min), source = 0.)
                call centers_coupled1%new(N_max, dummy=.true.)
                call centers_coupled2%new(N_min, dummy=.true.)
                allocate(mask(size(nano1%centers, dim = 2)), source = .true.)
                cnt = 0
                do i = 1, N_min !compare based on centers2
                    dist(i) = pixels_dist(nano2%centers(:,i),nano1%centers(:,:),'min',mask,location)
                    if(dist(i)*nano1%smpd > VDW_RAD_PT+0.2*VDW_RAD_PT) then
                        dist(i) = 0.
                        cnt = cnt+1
                    else
                        dist_sq(i) = dist(i)**2        !formula wants them square
                         if(DEBUG_HERE) then
                            write(logfhandle,*) 'ATOM', i,'coordinates: ', nano2%centers(:,i), 'coupled with '
                            write(logfhandle,*) '    ',location, 'coordinates: ', nano1%centers(:,location(1)), 'DIST^2= ', dist_sq(i), 'DIST = ', dist(i)
                         endif
                        call centers_coupled2%set_coord(i,(nano2%centers(:,i)-1.)*0.358/nano2%SCALE_FACTOR)
                        call centers_coupled1%set_coord(i,(nano1%centers(:,location(1))-1.)*0.358/nano1%SCALE_FACTOR)
                        mask(location(1)) = .false. ! not to consider the same atom more than once
                    endif
                enddo
            endif
            call centers_coupled1%writepdb('atom_couples_vol1')
            call centers_coupled2%writepdb('atom_couples_vol2')
            write(logfhandle,*) cnt, 'atoms correspond with error bigger tahn 5pm correspond. They are ', real(cnt)/real(N_min)*100., '% of the # of atoms'
            rmsd = sqrt(sum(dist_sq)/real(N_min-cnt))
            allocate(dist_no_zero(N_min-cnt), source = 0.) ! not to consider in the hist the uncoupled atoms
            dist_no_zero = pack(dist, dist>TINY)
            dist_no_zero = dist_no_zero*nano1%smpd ! report distances in Amstrongs
            call hist(dist_no_zero, 50)
            write(logfhandle,*)count(dist_no_zero <= 0.5),'atoms correspond withing 5pm. They are about', count(dist_no_zero <= 0.5)*100./(N_min-cnt), '% of the coupled atoms'
            write(logfhandle,*) 'RMSD =', rmsd*nano1%smpd, 'A'
            if(present(r)) r=rmsd
            deallocate(dist, dist_sq, dist_no_zero, mask)
        end subroutine atomic_position_rmsd
    end subroutine compare_atomic_models

    ! This subroutine plots the histogram of the within-atoms
    ! distances distribution within the nanoparticle nano.
    ! To each atom the distance assigned is the mean distance
    ! between the atom and its 6 nearest neighbours.
    subroutine distances_distribution(self)
        use simple_atoms, only : atoms
        class(nanoparticle), intent(inout) :: self
        integer :: i, nsz
        real    :: neigh_coords(3,6)
        real    :: avg_dist(self%n_cc)
        logical :: mask(6)
        mask = .true.
        do i = 1, self%n_cc
            call find_nearest_neigh(self,i, neigh_coords)
            avg_dist(i) =  pixels_dist(self%centers(:,i), neigh_coords, 'sum', mask)
            avg_dist(i) = avg_dist(i)/6. ! average
        enddo
        do i = 1,self%n_cc
            call self%centers_pdb%set_beta(i, avg_dist(i))
        enddo
        call self%centers_pdb%writepdb(trim(self%fbody)//'DistancesDistr')
    contains     ! This subroutine identifies the coordinates of the centers of the 6
        ! nearest neighbour of the atom identified by label in nano.
        subroutine find_nearest_neigh(nano, label, neigh_coords)
            type(nanoparticle),  intent(inout) :: nano
            integer,             intent(in)    :: label
            real,                intent(out)   :: neigh_coords(3,6)
            integer :: i, location(1)
            real    :: dist
            logical :: mask(nano%n_cc)
            mask = .true. !initialise
            mask(label) = .false. !do not consider the selected atom
            do i = 1, 6
                dist =  pixels_dist(nano%centers(:,label), nano%centers(:,:), 'min', mask, location)
                neigh_coords(1:3,i) = nano%centers(1:3,location(1))
                mask(location(1))  = .false.
            enddo
        end subroutine find_nearest_neigh
    end subroutine distances_distribution

    ! This subroutine has visualization purpose only.
    ! It prints out the 3 asymmetric units in a nanoparticle
    ! with c3 symmetry.
    subroutine print_asym_unit(self)
        class(nanoparticle), intent(inout) :: self
        real, allocatable :: rmat1(:,:,:), rmat2(:,:,:), rmat3(:,:,:)
        type(image) :: img_asym
        integer     :: i, j
        integer     :: x, y
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
                vec(:) = real([i-x,j-y])  !translate the vector in the origin
                ang1 = ang2D_vecs(vec1,vec)
                if(ang1 > 120. .or. i <= x) then
                    rmat1(i,j,:) = 0.    !set to 0 the all line
                endif
                if(ang1 > 120. .or. i > x) then
                    rmat2(i,j,:) = 0.    !set to 0 the all line
                endif
                if(ang1 < 120.) rmat3(i,j,:) = 0.
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
    subroutine find_centers(self, coords)
       use simple_atoms, only: atoms
       use gnufor2
       class(nanoparticle),          intent(inout) :: self
       real, optional, allocatable,  intent(out)   :: coords(:,:)
       logical, allocatable :: mask(:)
       integer, allocatable :: sz(:) !cc size to discard bottom 10/15%
       real,    allocatable :: rmat_cc(:,:,:)
       real,    allocatable :: rmat(:,:,:)
       integer :: i, n_discard, label(1)
       ! next lines are too expensive, work on that
       call self%img_bin%find_connected_comps(self%img_cc)
       !Connected components size polishing, calculate the size and
       !discard bottom 10% connected components
       !It takes a lot of time and anyway I discard
       !outliers wrt the contact score, but I don't think it's enough.
       rmat_cc = self%img_cc%get_rmat()
       rmat    = self%img_bin%get_rmat()
       self%n_cc = int(maxval(rmat_cc))
       n_discard = self%n_cc*15/100
       sz = self%img_cc%size_connected_comps()
       allocate(mask(size(sz)), source = .true.)
       do i = 1, n_discard
           label(:) = minloc(sz,mask)
           where(abs(rmat_cc-label(1)) < TINY)
               rmat_cc = 0.
               rmat    = 0.
           endwhere
       mask(label(1)) = .false.
       enddo
       call self%img_cc%set_rmat(rmat_cc) !connected components clean up
       call self%img_bin%set_rmat(rmat)
       ! re-order ccs
       call self%img_cc%order_cc()
       if(allocated(rmat_cc)) deallocate(rmat_cc)
       rmat_cc   = self%img_cc%get_rmat()
       call self%img_bin%set_rmat(rmat)
       self%n_cc = nint(maxval(rmat_cc))
       ! global variables allocation
       allocate(self%centers(3, self%n_cc),source = 0.)     !global variable
       do i=1,self%n_cc
           self%centers(:,i) = atom_masscen(self,i)
       enddo
       ! saving centers coordinates, optional
       if(present(coords)) allocate(coords(3, self%n_cc), source = self%centers)
       !sum_dist_closest_atom = 0. !initialise
       if(allocated(mask)) deallocate(mask)
       allocate(mask (self%n_cc), source = .true.)
       allocate(self%dists(self%n_cc), source = 0.)
       do i = 1, self%n_cc
           mask(:) = .true.  ! restore
           mask(i) = .false. ! not to consider the pixel itself, but I think its unnecessary
           self%dists(i) = pixels_dist(self%centers(:,i), self%centers,'min', mask)
       enddo
       call hist(self%dists*self%smpd, 20)
       self%avg_dist_atoms = sum(self%dists(:))/real(self%n_cc)
       deallocate(rmat,rmat_cc,mask)
    contains

       !This function calculates the centers of mass of an
       !atom. It takes in input the image, its connected
       !component (cc) image and the label of the cc that
       !identifies the atom whose center of mass is to be
       !calculated. The result is stored in m.
       function atom_masscen(self, label) result(m)
           type(nanoparticle), intent(inout) :: self
           integer,            intent(in)    :: label
           real(sp)             :: m(3)  !mass center coords
           real,    allocatable :: rmat_in(:,:,:), rmat_cc_in(:,:,:)
           integer :: i, j, k
           integer :: sz !sz of the cc identified by label
           rmat_in    = self%img_bin%get_rmat()
           rmat_cc_in = self%img_cc%get_rmat()
           where(     abs(rmat_cc_in-real(label)) > TINY) rmat_in = 0.
           sz = count(abs(rmat_cc_in-real(label)) < TINY)
           m = 0.
           !ASK HOW TO PARALLELISE IT
           !omp do collapse(3) reduction(+:m) private(i,j,k) reduce schedule(static)
           do i = 1, int(real(self%ldim(1))*self%SCALE_FACTOR)
               do j = 1, int(real(self%ldim(2))*self%SCALE_FACTOR)
                   do k = 1, int(real(self%ldim(3))*self%SCALE_FACTOR)
                       if(abs(rmat_in(i,j,k))> TINY) m = m + rmat_in(i,j,k)*[i,j,k]
                   enddo
               enddo
           enddo
           !omp end do
           m = m/real(sz)
           deallocate(rmat_in, rmat_cc_in)
       end function atom_masscen
    end subroutine find_centers

     ! This subroutine discard outliers that resisted binarization.
     ! It calculates the contact score of each atom and discards the bottom
     ! 5% of the atoms according to the contact score.
     ! It modifies the img_bin and img_cc instances deleting the
     ! identified outliers.
     ! Contact score: fix a radius and an atom A. Count the number N of atoms
     ! in the sphere centered in A of that radius. Define contact score of A
     ! equal to N. Do it for each atom.
     subroutine discard_outliers(self, radius)
         class(nanoparticle), intent(inout) :: self
         real, optional,      intent(in)    :: radius !radius of the sphere to consider
         real    :: dist
         integer, allocatable :: contact_scores(:)
         real,    allocatable :: rmat(:,:,:), rmat_cc(:,:,:)
         real,    allocatable :: centers_polished(:,:)
         logical :: mask(self%n_cc)
         integer :: cc
         integer :: cnt      !contact_score
         integer :: loc(1)
         integer :: n_discard
         integer :: label(1)
         real    :: rradius  !default
         rradius = 2.*self%avg_dist_atoms
         if(present(radius)) rradius = radius
         allocate(contact_scores(self%n_cc), source = 0)
         do cc = 1, self%n_cc !fix the atom
             dist = 0.
             mask(1:self%n_cc) = .true.
             cnt = 0
             mask(cc)  = .false.
             do while(dist < rradius)
                 dist = pixels_dist(self%centers(:,cc), self%centers,'min', mask, loc)
                 mask(loc) = .false.
                 cnt = cnt + 1
             enddo
             contact_scores(cc) = cnt - 1 !-1 because while loop counts one extra before exiting
             !call self%centers_pdb%set_beta(cc, real(contact_scores(cc)))
         enddo
         ! call self%centers_pdb%writepdb('contact_scores')
         n_discard = self%n_cc*10/100 ! discard bottom 10%??
         mask(1:self%n_cc) = .true.
         rmat    = self%img_bin%get_rmat()
         rmat_cc = self%img_cc%get_rmat()
         ! Removing outliers from the binary image and the connected components image
         do cc = 1, n_discard
             label = minloc(contact_scores, mask)
             where(abs(rmat_cc-real(label(1))) < TINY)
                rmat = 0.   ! discard atom corresponding to label_discard(i)
                rmat_cc = 0.
             endwhere
             mask(label(1)) = .false.
             self%centers(:,label) = -1. !identify to discard them
         enddo
         call self%img_bin%set_rmat(rmat)
         call self%img_cc%set_rmat(rmat_cc)
         call self%img_cc%order_cc()
         rmat_cc = self%img_cc%get_rmat()
         self%n_cc = nint(maxval(rmat_cc)) !update
         ! fix centers, removing the discarded ones
         allocate(centers_polished(3, self%n_cc), source = 0.)
         cnt = 0
         do cc = 1, size(self%centers, dim = 2)
             if(abs(self%centers(1,cc)+1.) > TINY) then
                 cnt = cnt + 1
                 centers_polished(:,cnt) = self%centers(:,cc)
             endif
         enddo
         deallocate(self%centers)
         allocate(self%centers(3,self%n_cc),source = centers_polished)
         ! centers plot
         call self%centers_pdb%new(self%n_cc, dummy=.true.)
         do cc=1,self%n_cc
             call self%centers_pdb%set_coord(cc,(self%centers(:,cc)-1.)*self%smpd/self%SCALE_FACTOR)
         enddo
         call self%centers_pdb%writepdb(trim(self%fbody)//'_atom_centers')
         call self%img_bin%write(trim(self%fbody)//'BIN.mrc')
         call self%distances_distribution() !pdb file colored wrt atom distances distribuition
         deallocate(rmat, rmat_cc, centers_polished)
     end subroutine discard_outliers

     ! calc the avg of the centers coords
     function nanopart_masscen(self) result(m)
         class(nanoparticle), intent(inout) :: self
         real    :: m(3)  !mass center coords
         integer :: i, j, k
         m = 0.
         do i = 1, self%n_cc
              m = m + 1.*self%centers(:,i)
         enddo
         m = m/real(self%n_cc)
         if(self%ldim(3) == 1) m(3) = 0. !for 2D imgs
     end function nanopart_masscen

    ! This function estimates the circularity of each blob
    ! (atom) in the following way:
    !    c := 6*pi*V/S.
    ! c = 1 <=> the blob is a perfect sphere, otherwise it's < 1.
    ! See circularity definition (in shape factor page) and adapt
    ! it to 3D volumes.
    subroutine calc_circularity(self)
        class(nanoparticle), intent(inout) :: self
        real,    allocatable :: rmat_cc(:,:,:)
        logical, allocatable :: border(:,:,:)
        integer, allocatable :: imat(:,:,:)
        integer              :: label
         integer, allocatable :: v(:), s(:) !volumes and surfaces of each atom
        real, parameter :: pi = 3.14159265359
        real            :: avg_circularity
        allocate(self%circularity(self%n_cc), source = 0.)
        rmat_cc = self%img_cc%get_rmat()
        allocate(imat(self%dim_over(1),self%dim_over(2),self%dim_over(3)),source = nint(rmat_cc)) !for function get_pixel_pos
        allocate(v(self%n_cc), s(self%n_cc), source = 0)
        avg_circularity = 0.
        do label = 1, self%n_cc
             v(label) = count(imat == label) !just count the number of vxls in the cc
            call self%img_cc%border_mask(border, label) !the subroutine border_mask allocated/deallocates itself matrix border
            s(label) = count(border .eqv. .true.)
            self%circularity(label) = (6.*sqrt(pi)*real(v(label)))/sqrt(real(s(label))**3)
            if(DEBUG_HERE) write(logfhandle,*) 'atom = ', label, 'vol = ', v(label), 'surf = ', s(label), 'circ = ', self%circularity(label)*self%smpd
            avg_circularity = avg_circularity + self%circularity(label)
        enddo
        avg_circularity = avg_circularity / self%n_cc
        write(logfhandle,*) 'AVG CIRCULARITY particle = ', self%fbody,  avg_circularity
        deallocate(v, s, imat, border, rmat_cc)
    end subroutine calc_circularity

    ! This subroutine takes in input a connected component (cc) image
    ! and the label of one of its ccs and calculates its aspect ratio, which
    ! is defined as the ratio of the width and the height.
    ! The idea behind this is that the center of the cc is calculated,
    ! than everything is deleted except the borders of the cc. Finally,
    ! in order to calculate the width and the height, the min/max
    ! distances between the center and the borders are calculated. The
    ! aspect ratio is the ratio of those 2 distances.
    subroutine calc_aspect_ratio_private(self, label, ratio)
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: label
        real,                intent(out)   :: ratio
        integer, allocatable :: imat(:,:,:)
        integer, allocatable :: pos(:,:)
        real,    allocatable :: rmat_cc(:,:,:)
        logical, allocatable :: border(:,:,:)
        logical, allocatable :: mask_dist(:) !for min and max dist calculation
        integer :: location(1) !location of the farest vxls of the atom from its center
        integer :: i
        real    :: shortest_dist, longest_dist
        rmat_cc = self%img_cc%get_rmat()
        allocate(imat( self%dim_over(1),self%dim_over(2),self%dim_over(3) ),source = nint(rmat_cc)) !for function get_pixel_pos
        call self%img_cc%border_mask(border, label, .true.) !use 4neigh instead of 8neigh
        where(border .eqv. .true.)
            imat = 1
        elsewhere
            imat = 0
        endwhere
        call get_pixel_pos(imat,pos)   !pxls positions of the shell
        if(allocated(mask_dist)) deallocate(mask_dist)
        allocate(mask_dist(size(pos, dim = 2)), source = .true. )
        shortest_dist = pixels_dist(self%centers(:,label), real(pos),'min', mask_dist, location)
        longest_dist  = pixels_dist(self%centers(:,label), real(pos),'max', mask_dist, location)
        self%loc_longest_dist(:3, label) =  pos(:3,location(1))
        if(abs(longest_dist) > TINY) then
            ratio = shortest_dist/longest_dist
        else
            ratio = 0.
            if(DEBUG_HERE) write(logfhandle,*) 'cc ', label, 'LONGEST DIST = 0'
        endif
        !if(DEBUG_HERE) then
            write(logfhandle,*) 'ATOM #          ', label
            write(logfhandle,*) 'shortest dist = ', shortest_dist*self%smpd
            write(logfhandle,*) 'longest  dist = ', longest_dist*self%smpd
            write(logfhandle,*) 'RATIO         = ', ratio
        !endif
        deallocate(rmat_cc, border, imat, pos, mask_dist)
    end subroutine calc_aspect_ratio_private

    subroutine calc_aspect_ratio(self)
        use gnufor2
        class(nanoparticle), intent(inout) :: self
        integer :: label
        allocate(self%ratios(self%n_cc),             source = 0.)
        allocate(self%loc_longest_dist(3,self%n_cc), source = 0 )
        do label = 1, self%n_cc
            call calc_aspect_ratio_private(self, label, self%ratios(label))
        enddo
        call hist(self%ratios, 20)
        ! To dump some of the analysis on aspect ratios on file compatible
        ! with Matlab.
         open(123, file='AspectRatios')
         write (123,*) 'AR=[...'
         do label = 1, size(self%ratios)
             write (123,'(A)', advance='no') trim(real2str(self%ratios(label)))
             if(label < size(self%ratios)) write (123,'(A)', advance='no') ', '
         end do
         write (123,*) '];'
         close(123)
        write(logfhandle,*)'**aspect ratios calculations completed'
    end subroutine calc_aspect_ratio

    ! This subrotuine indentifies the 'rows' of atoms in the z direction.
    ! The idea is to look at the projection of the nanoparticle on the
    ! xy plane and identify all the atoms whose center has (almost, see
    ! MAX_DIST_CENTERS) the same x and y coords.
    ! The inputs are:
    ! -) centers, coordinates of the centers of mass of the atoms;
    ! For how it is built, self%centers(:3,i) contains the coords of the
    ! center of mass of the i-th cc.
    subroutine make_soft_mask(self) !change the name
        use gnufor2
        class(nanoparticle), intent(inout) :: self
        call self%img_bin%grow_bins(nint(0.5*self%avg_dist_atoms)+1)
        call self%img_bin%cos_edge(6)
        call self%img_bin%write(trim(self%fbody)//'SoftMask.mrc')
    end subroutine make_soft_mask

    ! This subroutine takes in input 2 2D vectors, centered in the origin
    ! and it gives as an output the angle between them, IN DEGREES.
    function ang2D_vecs(vec1, vec2) result(ang)
        real, intent(inout) :: vec1(2), vec2(2)
        real :: ang        !output angle
        real :: ang_rad    !angle in radians
        real :: mod1, mod2
        real :: dot_prod
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

    ! Polarization search via angle variance. The considered angle is
    ! the angle between the vector [0,0,1] and the direction of the
    ! longest dim of the atom.
    ! This is done for each of the atoms and then there is a search
    ! of the similarity between the calculated angles.
    subroutine search_polarization(self)
        class(nanoparticle), intent(inout) :: self
        real, allocatable :: loc_ld_real(:,:)
        integer :: k, i
        real    :: vec_fixed(3)     !vector indentifying z direction
        real    :: theta, mod_1, mod_2, dot_prod
        logical :: mask(self%n_cc) !printing purposes
        if(allocated(self%ang_var)) deallocate(self%ang_var)
           allocate (self%ang_var(self%n_cc), source = 0.)
        if(allocated(loc_ld_real)) deallocate(loc_ld_real)
           allocate (loc_ld_real(3,size(self%loc_longest_dist, dim=2)), source = 0.)
        loc_ld_real(:3,:) = real(self%loc_longest_dist(:3,:))
        !consider fixed vector [0,0,1] (z direction)
        vec_fixed(1) = 0.
        vec_fixed(2) = 0.
        vec_fixed(3) = 1.
        write(logfhandle,*)'>>>>>>>>>>>>>>>> calculating angles wrt the vector [0,0,1]'
        do k = 1, self%n_cc
            self%ang_var(k) = ang3D_vecs(vec_fixed(:),loc_ld_real(:,k))
            if(DEBUG_HERE) write(logfhandle,*) 'ATOM ', k, 'angle between direction longest dim and vec [0,0,1] ', self%ang_var(k)
        enddo
        open(129, file='AnglesLongestDims')
        write (129,*) 'ang=[...'
        do k = 1, self%n_cc
                write (129,'(A)', advance='no') trim(int2str(k))
                write (129,'(A)', advance='no') ', '
                write (129,'(A)', advance='no') trim(real2str(self%ang_var(k)))
                write (129,'(A)')'; ...'
        end do
        write (129,*) '];'
        close(129)
        mask = .true. !not to overprint
        ! Searching coincident angles
        do i = 1,  self%n_cc-1
            do k = i+1, self%n_cc
                if(abs(self%ang_var(i)-self%ang_var(k))< 1.5 .and. mask(i)) then !1.5 degrees allowed
                    if(DEBUG_HERE) write(logfhandle,*)'Atoms', i, 'and', k,'have the same ang_var wrt vector [0,0,1]'
                    mask(k) = .false.
                endif
            enddo
        enddo
        deallocate(loc_ld_real)
    end subroutine search_polarization

    ! Affinity propagation clustering based on teh distribuition of
    ! the atom distances within the nanoparticle.
    ! I want to study the dustribuition of the atom distances within the nanoparticle.
    ! For example, I would expect to have atoms close to eachother in the center
    ! and far apart in the surface.
    subroutine affprop_cluster_dist_distr(self)
        use simple_aff_prop
        use simple_atoms, only : atoms
        class(nanoparticle), intent(inout) :: self
        type(aff_prop)       :: ap_nano
        real, allocatable    :: simmat(:,:)
        real                 :: simsum
        integer, allocatable :: centers_ap(:), labels_ap(:)
        integer, allocatable :: centers_merged(:), labels_merged(:)
        integer              :: i, j, ncls, nerr
        integer :: dim
        integer :: cnt
        logical, allocatable :: mask(:)
        dim = self%n_cc
        allocate(simmat(dim, dim), source = 0.)
        do i=1,dim-1
            do j=i+1,dim
                simmat(i,j) = -sqrt((self%dists(i)-self%dists(j))**2)
                simmat(j,i) = simmat(i,j)
            end do
        end do
        call ap_nano%new(dim, simmat)
        call ap_nano%propagate(centers_ap, labels_ap, simsum)
        ncls = size(centers_ap)
        write(logfhandle,*) 'NR OF CLUSTERS FOUND DISTS DISTR:', ncls
        do i = 1, self%n_cc
            call self%centers_pdb%set_beta(i, real(labels_ap(i)))
        enddo
        write(logfhandle,*) 'CENTERS DISTS DISTR'
        do i=1,ncls
            write(logfhandle,*) self%dists(centers_ap(i))
        end do
        call self%centers_pdb%writepdb(trim(self%fbody)//'ClusterWRTDistDistr')
        deallocate(simmat, labels_ap, centers_ap)
    end subroutine affprop_cluster_dist_distr

    !Affinity propagation clustering based on aspect ratios
    subroutine affprop_cluster_ar(self)
        use simple_aff_prop
        class(nanoparticle), intent(inout) :: self
        type(aff_prop)       :: ap_nano
        real, allocatable    :: simmat(:,:)
        real                 :: simsum
        integer, allocatable :: centers_ap(:), labels_ap(:)
        integer              :: i, j, ncls, nerr
        integer :: dim
        real, allocatable :: rmat(:,:,:)
        type(image)       :: img_clusters
        dim = size(self%ratios)
        allocate(simmat(dim, dim), source = 0.)
        write(logfhandle,'(a)') '**info(simple_aff_prop_unit_test): testing all functionality'
        do i=1,dim-1
            do j=i+1,dim
                simmat(i,j) = -sqrt((self%ratios(i)-self%ratios(j))**2)
                simmat(j,i) = simmat(i,j)
            end do
        end do
        call ap_nano%new(dim, simmat)
        call ap_nano%propagate(centers_ap, labels_ap, simsum)
        ncls = size(centers_ap)
        write(logfhandle,*) 'NR OF CLUSTERS FOUND AR:', ncls
        do i = 1, self%n_cc
            call self%centers_pdb%set_beta(i, real(labels_ap(i)))
        enddo
        write(logfhandle,*) 'CENTERS AR'
        do i=1,size(centers_ap)
            write(logfhandle,*) self%ratios(centers_ap(i))
        end do
        call self%centers_pdb%writepdb(trim(self%fbody)//'ClustersAspectRatio')
        deallocate(simmat, labels_ap, centers_ap)
    end subroutine affprop_cluster_ar

    ! Affinity propagation clustering based on ahgles wrt vec [0,0,1].
    subroutine affprop_cluster_ang(self)
        use simple_aff_prop
        use simple_atoms, only : atoms
        class(nanoparticle), intent(inout) :: self
        type(aff_prop)       :: ap_nano
        real, allocatable    :: simmat(:,:)
        real                 :: simsum
        integer, allocatable :: centers_ap(:), labels_ap(:)
        integer              :: i, j, ncls, nerr
        integer :: dim
        dim = self%n_cc
        allocate(simmat(dim, dim), source = 0.)
        do i=1,dim-1
            do j=i+1,dim
                simmat(i,j) = -sqrt((self%ang_var(i)-self%ang_var(j))**2)
                simmat(j,i) = simmat(i,j)
            end do
        end do
        call ap_nano%new(dim, simmat)
        call ap_nano%propagate(centers_ap, labels_ap, simsum)
        ncls = size(centers_ap)
        write(logfhandle,*) 'NR OF CLUSTERS FOUND ANG:', ncls
        do i = 1, self%n_cc
            call self%centers_pdb%set_beta(i, real(labels_ap(i)))
        enddo
        write(logfhandle,*) 'CENTERS ANG'
        do i=1,ncls
            write(logfhandle,*) self%ang_var(centers_ap(i))
        end do
        call self%centers_pdb%writepdb(trim(self%fbody)//'ClusteWRTdirLongestDim')
        deallocate(simmat, labels_ap, centers_ap)
    end subroutine affprop_cluster_ang

    !Affinity propagation clustering based on circularity.
    subroutine affprop_cluster_circ(self)
      use simple_aff_prop
      use simple_atoms, only : atoms
      class(nanoparticle), intent(inout) :: self
      type(aff_prop)       :: ap_nano
      real, allocatable    :: simmat(:,:)
      real                 :: simsum
      integer, allocatable :: centers_ap(:), labels_ap(:)
      integer              :: i, j, ncls, nerr
      integer :: dim
      dim = self%n_cc
      allocate(simmat(dim, dim), source = 0.)
      do i=1,dim-1
      do j=i+1,dim
              simmat(i,j) = -sqrt((self%circularity(i)-self%circularity(j))**2)
              simmat(j,i) = simmat(i,j)
          end do
      end do
      call ap_nano%new(dim, simmat)
      call ap_nano%propagate(centers_ap, labels_ap, simsum)
      ncls = size(centers_ap)
      write(logfhandle,*) 'NR OF CLUSTERS FOUND CIRC:', ncls
      do i = 1, self%n_cc
          call self%centers_pdb%set_beta(i, real(labels_ap(i)))
      enddo
      write(logfhandle,*) 'CENTERS CIRC'
      do i=1,ncls
          write(logfhandle,*) self%ang_var(centers_ap(i))
      end do
      call self%centers_pdb%writepdb(trim(self%fbody)//'ClustersWRTcircularity')
      deallocate(simmat, labels_ap, centers_ap)
  end subroutine affprop_cluster_circ

    subroutine nanopart_cluster(self)
        class(nanoparticle), intent(inout) :: self
        ! clustering wrt to angle longest_dim-vector z=[0,0,1]
        call self%affprop_cluster_ang()
        ! clustering wrt aspect ratio
        call self%affprop_cluster_ar()
        ! clustering wrt circularity
        call self%affprop_cluster_circ()
    end subroutine nanopart_cluster

    subroutine over_sample(self)
        class(nanoparticle), intent(inout) :: self
        call self%img%fft()
        call self%img%pad(self%img_over_smp)
        call self%img_over_smp%ifft()
        if(DEBUG_HERE) call self%img_over_smp%write(trim(self%fbody)//'OverSampled.mrc')
    end subroutine over_sample

    subroutine kill_nanoparticle(self)
        class(nanoparticle), intent(inout) :: self
        self%ldim(3)           = 0
        self%smpd              = 0.
        self%nanop_mass_cen(3) = 0.
        self%avg_dist_atoms    = 0.
        self%SCALE_FACTOR      = 0.
        self%n_cc              = 0
        self%dim_over(3)       = 0
        call self%img%kill()
        call self%img_bin%kill()
        call self%img_cc%kill()
        call self%img_over_smp%kill()
        call self%centers_pdb%kill
        if(allocated(self%centers))          deallocate(self%centers)
        if(allocated(self%ratios))           deallocate(self%ratios)
        if(allocated(self%dists))            deallocate(self%dists)
        if(allocated(self%loc_longest_dist)) deallocate(self%loc_longest_dist)
        if(allocated(self%circularity))      deallocate(self%circularity)
        if(allocated(self%ang_var))          deallocate(self%ang_var)
    end subroutine kill_nanoparticle
end module simple_nanoparticles_mod

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !STANDARD DEVIATION CALCULATIONS WITHIN CENTERS OF THE CLUSTERS AND
! !AMONG DIFFERENT CENTERS
! ! Affinity propagation clustering based on ahgles wrt vec [0,0,1].
! subroutine affprop_cluster_ang(self)
!   use simple_aff_prop
!   use simple_atoms, only : atoms
!   class(nanoparticle), intent(inout) :: self
!   type(aff_prop)       :: ap_nano
!   real, allocatable    :: simmat(:,:)
!   real                 :: simsum
!   integer, allocatable :: centers_ap(:), labels_ap(:)
!   integer              :: i, j, ncls, nerr
!   integer :: dim
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   real    :: avg, stdev
!   real,    allocatable :: avg_within(:), stdev_within(:)
!   integer, allocatable :: cnt(:)
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   dim = self%n_cc
!   allocate(simmat(dim, dim), source = 0.)
!   do i=1,dim-1
!       do j=i+1,dim
!           simmat(i,j) = -sqrt((self%ang_var(i)-self%ang_var(j))**2)
!           simmat(j,i) = simmat(i,j)
!       end do
!   end do
!   call ap_nano%new(dim, simmat)
!   call ap_nano%propagate(centers_ap, labels_ap, simsum)
!   ncls = size(centers_ap)
!   write(logfhandle,*) 'NR OF CLUSTERS FOUND ANG:', ncls
!   do i = 1, self%n_cc
!       call self%centers_pdb%set_beta(i, real(labels_ap(i)))
!   enddo
!   write(logfhandle,*) 'CENTERS'
!   do i=1,ncls
!       write(logfhandle,*) self%ang_var(centers_ap(i))
!   end do
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !standard deviation within each center
!   allocate(  avg_within(ncls), source = 0.)
!   allocate(stdev_within(ncls), source = 0.)
!   allocate(cnt(ncls), source = 0)
!   do i = 1, ncls
!       do j = 1, self%n_cc
!           if(labels_ap(j) == i) then
!               cnt(i) = cnt(i) + 1 !cnt is how many atoms I have in each class
!               avg_within(i) = avg_within(i) + self%ang_var(j)
!           endif
!       enddo
!   enddo
!   avg_within = avg_within/cnt
!   write(logfhandle,*) 'number of atoms in each class = ', cnt, 'total atoms = ', sum(cnt)
!   write(logfhandle,*) 'avg_within = ', avg_within
!   do i = 1, ncls
!       do j = 1, self%n_cc
!           if(labels_ap(j) == i) stdev_within(i) = stdev_within(i) + (self%ang_var(j)-avg_within(i))**2
!       enddo
!   enddo
!   stdev_within = stdev_within/(cnt-1.)
!   stdev_within = sqrt(stdev_within)
!   write(logfhandle,*) 'stdev_within = ', stdev_within
!   !standard deviation among different centers
!   avg = 0.
!   do i = 1, ncls
!       avg = avg + self%ang_var(centers_ap(i))
!   enddo
!   avg = avg/real(ncls)
!   write(logfhandle,*) 'avg centers ', avg
!   stdev = 0.
!   do i = 1, ncls
!       stdev = stdev + (self%ang_var(centers_ap(i)) - avg)**2
!   enddo
!   stdev = stdev/(real(ncls)-1.)
!   stdev = sqrt(stdev)
!   write(logfhandle,*) 'standard deviation among centers: ', stdev
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   call self%centers_pdb%writepdb('ClustersAnglesWRTzvec')
!   deallocate(simmat, labels_ap, centers_ap)
! end subroutine affprop_cluster_ang

! ! TRYING TO MERGE CLASSES, STILL DOESN'T WORK BUT CLOSE
! ! Affinity propagation clustering based on teh distribuition of
! ! the atom distances within the nanoparticle.
! ! I want to study the dustribuition of the atom distances within the nanoparticle.
! ! For example, I would expect to have atoms close to eachother in the center
! ! and far apart in the surface.
! subroutine affprop_cluster_dist_distr(self)
!     use simple_aff_prop
!     use simple_atoms, only : atoms
!     class(nanoparticle), intent(inout) :: self
!     type(aff_prop)       :: ap_nano
!     real, allocatable    :: simmat(:,:)
!     real                 :: simsum
!     integer, allocatable :: centers_ap(:), labels_ap(:)
!     integer, allocatable :: centers_merged(:), labels_merged(:)
!     integer              :: i, j, ncls, nerr
!     integer :: dim
!     integer :: cnt
!     logical, allocatable :: mask(:)
!     dim = self%n_cc
!     allocate(simmat(dim, dim), source = 0.)
!     do i=1,dim-1
!         do j=i+1,dim
!             simmat(i,j) = -sqrt((self%dists(i)-self%dists(j))**2)
!             simmat(j,i) = simmat(i,j)
!         end do
!     end do
!     call ap_nano%new(dim, simmat)
!     call ap_nano%propagate(centers_ap, labels_ap, simsum)
!     ncls = size(centers_ap)
!     write(logfhandle,*) 'NR OF CLUSTERS FOUND DISTS DISTR:', ncls
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write(logfhandle,*)'>>>>>>>Merging similar classes>>>>>>'
!     allocate(mask(ncls), source = .true.) ! not to compare more than once the same couple
!     do i = 1, ncls-1
!         do j = i+1, ncls
!             if(mask(i) .and. mask(j)) then
!                 mask(i) = .false.
!                 if (abs(self%dists(centers_ap(i))-self%dists(centers_ap(j))) < 1.) then
!                     labels_ap(j)  = labels_ap(i) ! merge the classes
!                     centers_ap(j) = centers_ap(i)
!                     mask(j) = .false.
!                 endif
!             endif
!         enddo
!     enddo
!     write(logfhandle,*)'I have ', size(centers_ap), 'centers. And they are', centers_ap
!     write(logfhandle,*) 'CENTERS DISTS DISTR'
!     do i=1,ncls
!         write(logfhandle,*) self%dists(centers_ap(i))
!     end do
!     mask = .true. ! restore
!     cnt = 0
!     do i = 1, ncls-1
!         do j = i+1, ncls
!             if(mask(i) .and. mask(j)) then
!                 mask(i) = .false.
!                 if(centers_ap(i)==centers_ap(j)) then
!                     cnt = cnt + 1
!                     mask(j) = .false.
!                 endif
!             endif
!         enddo
!     enddo
!     write(logfhandle,*)'cnt = ', cnt
!     allocate(centers_merged(ncls-cnt), source = 0)
!     cnt = 0
!     do i = 1, ncls
!         if(.not. any(centers_merged(:)==centers_ap(i))) then
!             cnt = cnt + 1
!             centers_merged(cnt) = centers_ap(i)
!         endif
!     enddo
!     write(logfhandle,*)'Now I have ', size(centers_merged), 'centers. And they are', centers_merged
!     deallocate(mask)
!     write(logfhandle,*)'I have ', size(centers_ap), 'centers. And they are', centers_ap
!     write(logfhandle,*) 'CENTERS DISTS DISTR'
!     do i=1,ncls
!         write(logfhandle,*) self%dists(centers_ap(i))
!     end do
!     !!!!!!!!!!!!!!!!!!
!     if(allocated(mask)) deallocate(mask)
!     deallocate(centers_merged)
!     allocate(mask(self%n_cc), source = .true. )
!     cnt = 0
!     do i = 1, self%n_cc-1
!         do j = i+1, self%n_cc
!             if(mask(i) .and. mask(j)) then
!                 mask(i) = .false.
!                 if(labels_ap(i)==labels_ap(j)) then
!                     cnt = cnt + 1
!                     mask(j) = .false.
!                 endif
!             endif
!         enddo
!     enddo
!     write(logfhandle,*)'cnt = ', cnt
!     allocate(labels_merged(self%n_cc-cnt), source = 0)
!     cnt = 0
!     do i = 1, self%n_cc
!         if(.not. any(labels_merged(:)==labels_ap(i))) then
!             cnt = cnt + 1
!             labels_merged(cnt) = labels_ap(i)
!         endif
!     enddo
!     write(logfhandle,*)'I have ', size(labels_ap), 'labels'
!     write(logfhandle,*)'Now I have ', size(labels_merged), 'labels. And they are', labels_merged
!     deallocate(mask)
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     do i = 1, self%n_cc-cnt
!         call self%centers_pdb%set_beta(i, real(labels_merged(i)))
!         ! call self%centers_pdb%set_beta(i, real(labels_ap(i)))
!     enddo
!     call self%centers_pdb%writepdb(trim(self%fbody)//'ClusterMergedWRTDistDistr')
!     deallocate(labels_merged)
!     deallocate(simmat, labels_ap, centers_ap)
! end subroutine affprop_cluster_dist_distr
