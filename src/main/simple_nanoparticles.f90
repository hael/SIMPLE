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
integer, parameter :: N_THRESH   = 20      !number of thresholds for binarization
logical, parameter :: DEBUG_HERE = .false. !for debugging purposes

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
    !integer     :: dim_over(3)       = 0 !logical dimensions of the oversampled
    real,    allocatable  :: centers(:,:)
    real,    allocatable  :: ratios(:)
    real,    allocatable  :: ang_var(:)
    real,    allocatable  :: dists(:)
    integer, allocatable  :: loc_longest_dist(:,:)   !for indentific of the vxl that determins the longest dim of the atom
    character(len=STDLEN) :: partname   = ''   !fname
    character(len=STDLEN) :: fbody      = ''   !fbody
    character(len=STDLEN) :: output_dir = ''

  contains
    ! constructor
    procedure          :: new => new_nanoparticle
    ! getters/setters
    procedure          :: get_img
    procedure          :: set_img
    procedure          :: set_partname
    ! segmentation and statistics
    procedure          :: binarize => nanopart_binarization
    procedure          :: size_filtering
    procedure          :: find_centers
    procedure, private :: nanopart_masscen
    procedure, private :: calc_aspect_ratio_private
    procedure          :: calc_aspect_ratio
    procedure          :: discard_outliers
    procedure, private :: distances_distribution
    ! clustering
    procedure          :: cluster => nanopart_cluster
    procedure, private :: affprop_cluster_ar
    procedure, private :: affprop_cluster_ang
    procedure, private :: affprop_cluster_dist_distr
    procedure          :: search_polarization
    ! execution
    procedure          :: detect_atoms
    ! comparison
    procedure          :: compare_atomic_models
    ! visualization and output
    procedure          :: write_centers
    procedure          :: print_asym_unit
    ! others
    procedure          :: make_soft_mask
    procedure, private :: over_sample
    ! kill
    procedure          :: kill => kill_nanoparticle
end type nanoparticle

contains

    !constructor
    subroutine new_nanoparticle(self, fname, cline_smpd, sc_fac)
        use simple_syslib
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        real,                intent(in)    :: cline_smpd
        real, optional,      intent(in)    :: sc_fac
        integer :: nptcls
        real    :: ssc_fac
        real    :: smpd
        call self%kill
        call simple_getcwd(self%output_dir)
        ssc_fac = 1.
        if(present(sc_fac)) ssc_fac = sc_fac
        self%SCALE_FACTOR = ssc_fac
        call self%set_partname(fname)
        self%fbody = get_fbody(trim(fname), trim(fname2ext(fname)))
        self%smpd = cline_smpd
        call find_ldim_nptcls(self%partname,  self%ldim, nptcls, smpd)
        call self%img%new         (self%ldim, self%smpd)
        call self%img_bin%new     (int(real(self%ldim)*self%SCALE_FACTOR), self%smpd/self%SCALE_FACTOR)
        !call self%img_over_smp%new(int(real(self%ldim)*self%SCALE_FACTOR), self%smpd/self%SCALE_FACTOR)
        ! self%dim_over(1) = int(real(self%ldim(1))*self%SCALE_FACTOR)
        ! self%dim_over(2) = int(real(self%ldim(2))*self%SCALE_FACTOR)
        ! self%dim_over(3) = int(real(self%ldim(3))*self%SCALE_FACTOR)
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
        type(image)       :: img_bin_thresh(N_THRESH-1)
        type(image)       :: img_ccs_thresh(N_THRESH-1)
        real,    allocatable :: rmat(:,:,:), rmat_t(:,:,:)
        integer, allocatable :: rmat_int(:,:,:)
        integer, allocatable :: sz(:) !size of the ccs and correspondent label
        real    :: step               !histogram disretization step
        real    :: thresh             !binarization threshold
        integer :: ind(3)             !selected indexes corresponding to threshold for nanoparticle binarization
        integer :: position(1)
        integer :: cent_outside_atom(3)
        real    :: x_thresh(N_THRESH-1), y_med(N_THRESH-1)
        integer ::  i, cc
        logical, allocatable :: yes_no(:)
        real :: avg_sz, stdev_sz
        real :: t
        write(logfhandle,*) '****binarization, init'
        !call self%over_sample()
        rmat = self%img%get_rmat()
        allocate(rmat_t(self%ldim(1), self%ldim(2), self%ldim(3)), source = 0.)
        step = maxval(rmat)/real(N_THRESH)
        !$omp do collapse(1) schedule(static) private(i)
        do i = 1, N_THRESH-1
            call progress(i, N_THRESH-1)
            call img_bin_thresh(i)%new(self%ldim, self%smpd)
            thresh = real(i)*step
            where(rmat > thresh)
                rmat_t = 1.
            elsewhere
                rmat_t = 0.
            endwhere
            call img_bin_thresh(i)%set_rmat(rmat_t)
            call img_bin_thresh(i)%find_connected_comps(img_ccs_thresh(i))
            sz          = img_ccs_thresh(i)%size_connected_comps()
            x_thresh(i) = thresh
            y_med(i)    = median(real(sz))
        enddo
        !$omp end do
        write(logfhandle,*)'Intial threshold selected, starting refinement'
        !two consecutive thresholds after the ideal one
        ind(1:1) = maxloc(y_med)
        ind(2:2) = ind(1:1)+1
        ind(3:3) = ind(2:2)+1
        !$omp do collapse(1) schedule(static) private(i)
        do i = 1, 3
            rmat_t = img_ccs_thresh(ind(i))%get_rmat()
            self%n_cc = int(maxval(rmat_t))
            call self%find_centers(img_bin_thresh(ind(i)),img_ccs_thresh(ind(i)))
            allocate(yes_no(self%n_cc))
            yes_no = is_center_inside_atom(self,img_bin_thresh(ind(i)))
            cent_outside_atom(i) = count(yes_no .eqv. .false.)
            write(logfhandle,*)'For thresh ',x_thresh(ind(i)),', ', trim(int2str(cent_outside_atom(i))), ' centers are not inside the atom'
            deallocate(yes_no)
        enddo
        !$omp end do
        position(:) =  minloc(cent_outside_atom)
        write(logfhandle,*) 'Final threshold for binarization', x_thresh(ind(position(1)))
        call self%img_bin%copy(img_bin_thresh(ind(position(1))))
        call  self%img_cc%copy(img_ccs_thresh(ind(position(1))))
        !kill images
        !$omp do collapse(1) schedule(static) private(i)
        do i = 1, N_THRESH-1
            call img_bin_thresh(i)%kill
            call img_ccs_thresh(i)%kill
        enddo
        !$omp end do
        write(logfhandle,*) 'Final threshold selected, starting connected atoms erosion'
        sz = self%img_cc%size_connected_comps()
        avg_sz = real(sum(sz))/real(size(sz))
        stdev_sz = 0.
        do cc = 1, size(sz)
            stdev_sz = stdev_sz + (sz(cc) - avg_sz)**2
        enddo
        stdev_sz = sqrt(stdev_sz/(real(size(sz)-1)))
        t = avg_sz + 2.*stdev_sz !assuming Gaussian distrib, 95% is in [-2sigma,2sigma]
        write(logfhandle,*) 'Starting erosion of big atoms'
        if(DEBUG_HERE) then
            write(logfhandle,*)'avg   atom size   = ', avg_sz
            write(logfhandle,*)"stdev atom size   = ", stdev_sz
            write(logfhandle,*)'Erosion threshold = ', t
        endif
        rmat_int = int(self%img_cc%get_rmat())
        !$omp do collapse(1) schedule(static) private(cc)
        do cc = 1, size(sz)
            if(sz(cc) > t) then
                 call self%img_cc%erosion(cc)
            endif
        enddo
        !$omp end do
        ! update
        rmat = 0.
        rmat_int = int(self%img_cc%get_rmat())
        where(rmat_int > 0) rmat = 1.
        !update img bin
        call self%img_bin%set_rmat(rmat)
        !update image cc
        call self%img_bin%find_connected_comps(self%img_cc)
        !update centers
        call self%find_centers()
        ! deallocate
        if(allocated(sz))       deallocate(sz)
        if(allocated(rmat))     deallocate(rmat)
        if(allocated(rmat_t))   deallocate(rmat_t)
        if(allocated(rmat_int)) deallocate(rmat_int)
        if(allocated(yes_no))   deallocate(yes_no)
        write(logfhandle,*) '****binarization, completed'
    contains

        ! This function checks for each atom is the identified
        ! centers lies inside the atom or not. The result is
        ! saved in yes_no(:). The entry of yes_no corresponds
        ! to the connected component label of the atom.
        function is_center_inside_atom(self, img_bin) result(yes_no)
            class(nanoparticle), intent(inout) :: self
            type(image), intent(in) :: img_bin
            logical :: yes_no(self%n_cc)
            integer :: cc
            integer, allocatable :: imat(:,:,:)
            imat = int(img_bin%get_rmat())
            !$omp do collapse(1) schedule(static) private(cc)
            do cc = 1, self%n_cc
                if(imat(nint(self%centers(1,cc)),nint(self%centers(2,cc)),nint(self%centers(3,cc))) > 0.5) then
                    yes_no(cc) = .true.
                else
                    yes_no(cc) = .false.
                endif
            enddo
            !$omp end do
            deallocate(imat)
        end function is_center_inside_atom
    end subroutine nanopart_binarization

    ! This subroutine performs size filtering on the connected
    ! components image of the nanoparticle. It calculates the
    ! size of all the connected components, the average and the standard
    ! deviation. It hypothesises gaussian distribution, so it discards
    ! the connected components which are outside the range [-2sigma,2sigma].
    subroutine size_filtering(self)
        class(nanoparticle), intent(inout) :: self
        real,    allocatable :: rmat(:,:,:)
        integer, allocatable :: imat_cc(:,:,:)
        integer, allocatable :: sz(:)
        integer :: cc
        real    :: avg_sz, stdev_sz
        write(logfhandle, *) '****size filtering, init'
        imat_cc = int(self%img_cc%get_rmat())
        rmat    = self%img_bin%get_rmat()
        sz = self%img_cc%size_connected_comps()
        avg_sz = real(sum(sz))/real(size(sz))
        stdev_sz = 0.
        do cc = 1, size(sz)
            stdev_sz = stdev_sz + (sz(cc) - avg_sz)**2
        enddo
        stdev_sz = sqrt(stdev_sz/(real(size(sz)-1)))
        !assuming Gaussian distrib, 95% is in [-2sigma,2sigma] and 68% is in [-sigma,sigma]
        !big ccs have already been removed by erosion. Now we need to remove too small
        !ccs, they usually represent background noise.
        !$omp do collapse(1) schedule(static) private(cc)
        do cc = 1, size(sz)
            if(sz(cc)<avg_sz-1.*stdev_sz ) then
                where(imat_cc==cc)
                    imat_cc = 0
                    rmat    = 0.
                endwhere
             endif
        enddo
        !$omp end do
        call self%img_cc%set_rmat(real(imat_cc)) !connected components clean up
        call self%img_bin%set_rmat(rmat)
        !update img_cc: re-order ccs
        call self%img_cc%order_cc()
        if(allocated(imat_cc)) deallocate(imat_cc)
        !update img_bin
        call self%img_bin%set_rmat(rmat)
        imat_cc   = int(self%img_cc%get_rmat())
        !update n_cc
        self%n_cc = maxval(imat_cc)
        !update centers
        call self%find_centers()
        write(logfhandle, *) '****size filtering, completed'
        deallocate(rmat,imat_cc,sz)
    end subroutine size_filtering

   ! This is the subroutine that executes all the steps
   ! for the structural analysis of one nananoparticle
   ! 3D reconstruction.
    subroutine detect_atoms(self)
        class(nanoparticle), intent(inout) :: self
        ! Nanoparticle binarization
        call self%binarize() !scale factor is for over sampling purposes
        ! Outliers discarding
        call self%discard_outliers()
        ! Atom-to-atom distances distribution estimation
        call self%distances_distribution()
        ! Aspect ratios calculations
        call self%calc_aspect_ratio()
        ! Make soft mask
        call self%make_soft_mask()
        ! Polarization search
        call self%search_polarization()
        ! Clustering
        call self%cluster
    end subroutine detect_atoms

    ! This subroutine takes in input 2 nanoparticle and their
    ! correspondent filenames.
    ! The output is the rmsd of their atomic models.
    subroutine compare_atomic_models(nano1,nano2)
        class(nanoparticle), intent(inout) :: nano1, nano2 !nanoparticles to compare
        real, allocatable :: centers1(:,:), centers2(:,:)
        real    :: ada       ! minimum average dist atoms between nano1 and nano2
        real    :: rmsd
        integer :: i
        logical :: print_ar ! for printing aspect ratios statistics
        logical :: print_as ! for printing asymmetric units in c3-sym nanoparticles
        print_ar = .false.
        print_as = .false.
        if(print_as) then
            call nano1%print_asym_unit()
            call nano2%print_asym_unit()
        endif
        open(121, file='CompareAtomicModels')
        write(unit = 121, fmt = '(a)') '>>>>>>>>>   COMPARE NANO   >>>>>>>>>'
        write(unit = 121, fmt = '(a)') ''
        write(unit = 121, fmt = '(a)') 'Comparing atomic models of particles'
        write(unit = 121, fmt = '(a,a)') trim(nano1%fbody), ' ---> vol1'
        write(unit = 121, fmt = '(a)') 'and'
        write(unit = 121, fmt = '(a,a)') trim(nano2%fbody), ' ---> vol2'
        write(unit = 121, fmt = '(a)')  '>>>>>>>>>VOLUME COMPARISION>>>>>>>>'
        call nano1%binarize()
        call nano2%binarize()
        ! Outliers discarding (and size filtering)
        call nano1%discard_outliers(1)
        call nano2%discard_outliers(2)
        ! Aspect ratios calculations
        call nano1%calc_aspect_ratio(print_ar)
        call nano2%calc_aspect_ratio(print_ar)
        ! Ouput file, come back to initial folder
        call simple_chdir(trim(nano1%output_dir), errmsg="simple_nanoparticles :: compare_atomic_models, simple_chdir; ")
        write(unit = 121, fmt = '(a,f6.3,a)') 'avg dist between atoms in vol1:', nano1%avg_dist_atoms*nano1%smpd,' A'
        write(unit = 121, fmt = '(a,f6.3,a)') 'avg dist between atoms in vol2:', nano2%avg_dist_atoms*nano2%smpd,' A'
        ! RMSD calculation
        ada = min(nano1%avg_dist_atoms, nano2%avg_dist_atoms)
        call atomic_position_rmsd(nano1,nano2, ada, rmsd)
        write(logfhandle,*) '***comparison completed'
        close(121)
        call nano1%kill
        call nano2%kill
        if(allocated(centers1)) deallocate(centers1)
        if(allocated(centers2)) deallocate(centers2)
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
            real,    allocatable :: dist(:), dist_sq(:), dist_no_zero(:), dist_close(:)
            integer :: location(1)
            integer :: i, j
            integer :: N_min !min{nb atoms in nano1, nb atoms in nano2}
            integer :: N_max !max{nb atoms in nano1, nb atoms in nano2}
            integer :: cnt,cnt2,cnt3
            real    :: sum, rmsd
            real :: coord(3)
            real :: avg, stdev, m(3), tmp_max, d !for statistics calculation
            type(atoms) :: centers_coupled1, centers_coupled2 !visualization purposes
            type(atoms) :: centers_close1, centers_close2
            type(atoms) :: couples1
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
                cnt = 0
                cnt2= 0
                !$omp do collapse(1) schedule(static) private(i)
                do i = 1, N_max !compare based on centers2
                    dist(i) = pixels_dist(nano2%centers(:,i),nano1%centers(:,:),'min',mask,location)
                    if(dist(i)*nano2%smpd > 2.2) then
                        dist(i) = 0. !it means there is no correspondent atom in the other nano
                        cnt = cnt+1  !to discard them in the rmsd calculation
                        call couples1%set_coord(i,(nano2%centers(:,location(1))-1.)*nano2%smpd/nano2%SCALE_FACTOR)
                        ! remove the atoms from the pdb file
                        call centers_coupled1%set_occupancy(i,0.)
                        call centers_coupled2%set_occupancy(i,0.)
                        call centers_close1%set_occupancy(i,0.)
                        call centers_close2%set_occupancy(i,0.)
                        call couples1%set_occupancy(i,0.)
                    elseif(dist(i)*nano2%smpd < 0.5) then
                        dist_close(i) = dist(i)**2
                        call centers_close2%set_coord(i,(nano2%centers(:,i)-1.)*nano2%smpd/nano2%SCALE_FACTOR)
                        call centers_close1%set_coord(i,(nano1%centers(:,location(1))-1.)*nano1%smpd/nano1%SCALE_FACTOR)
                        ! remove the atoms from the pdb file
                        call centers_coupled1%set_occupancy(i,0.)
                        call centers_coupled2%set_occupancy(i,0.)
                        call couples1%set_occupancy(i,0.)
                    elseif(dist(i)*nano2%smpd > 0.5 .and. dist(i)*nano2%smpd<=2.2 ) then  !to save the atoms which correspond with a precision in the range [0,220] pm
                        cnt2 = cnt2 + 1
                        call centers_coupled2%set_coord(i,(nano2%centers(:,i)-1.)*nano2%smpd/nano2%SCALE_FACTOR)
                        call centers_coupled1%set_coord(i,(nano1%centers(:,location(1))-1.)*nano1%smpd/nano1%SCALE_FACTOR)
                        call couples1%set_coord(i,(nano2%centers(:,location(1))-1.)*nano2%smpd/nano2%SCALE_FACTOR)
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
                !$omp end do
            else
                N_min = size(nano2%centers, dim = 2)
                N_max = size(nano1%centers, dim = 2)
                write(unit = 121, fmt = '(a,i3)') 'NB atoms in vol1   ', N_max
                write(unit = 121, fmt = '(a,i3)') 'NB atoms in vol2   ', N_min
                write(unit = 121, fmt = '(i3,a)') abs(N_max-N_min), 'atoms do NOT correspond'
                allocate(dist(N_max), dist_sq(N_max), source = 0.)
                allocate(dist_close(N_max), source = 0.) ! there are going to be unused entry of the vector
                call centers_coupled1%new(N_max, dummy=.true.)
                call centers_coupled2%new(N_max, dummy=.true.)
                call centers_close1%new  (N_max, dummy=.true.)
                call centers_close2%new  (N_max, dummy=.true.)
                allocate(mask(N_min), source = .true.)
                cnt = 0
                cnt2= 0
                !$omp do collapse(1) schedule(static) private(i)
                do i = 1, N_max !compare based on centers1
                    dist(i) = pixels_dist(nano1%centers(:,i),nano2%centers(:,:),'min',mask,location)
                    if(dist(i)*nano2%smpd > 2.2) then ! 2.2 is the biggest lattice spacing they found in the paper
                        dist(i) = 0.
                        cnt = cnt+1
                        call couples1%set_coord(i,(nano1%centers(:,location(1))-1.)*nano1%smpd/nano2%SCALE_FACTOR)
                        ! remove the atoms from the pdb file
                        call centers_coupled1%set_occupancy(i,0.)
                        call centers_coupled2%set_occupancy(i,0.)
                        call centers_close1%set_occupancy(i,0.)
                        call centers_close2%set_occupancy(i,0.)
                        call couples1%set_occupancy(i,0.)
                    elseif(dist(i)*nano2%smpd <= 0.5) then
                            dist_close(i) = dist(i)**2
                            call centers_close1%set_coord(i,(nano1%centers(:,i)-1.)*nano1%smpd/nano1%SCALE_FACTOR)
                            call centers_close2%set_coord(i,(nano2%centers(:,location(1))-1.)*nano2%smpd/nano2%SCALE_FACTOR)
                            ! remove the atoms from the pdb file
                            call centers_coupled1%set_occupancy(i,0.)
                            call centers_coupled2%set_occupancy(i,0.)
                            call couples1%set_occupancy(i,0.)
                    elseif(dist(i)*nano2%smpd > 0.5 .and. dist(i)*nano2%smpd<=2.2 ) then  !to save the atoms which correspond with a precision in the range [0,220] pm
                        cnt2 = cnt2 + 1
                            call centers_coupled1%set_coord(i,(nano1%centers(:,i)-1.)*nano1%smpd/nano1%SCALE_FACTOR)
                            call centers_coupled2%set_coord(i,(nano2%centers(:,location(1))-1.)*nano2%smpd/nano2%SCALE_FACTOR)
                            call couples1%set_coord(i,(nano1%centers(:,location(1))-1.)*nano1%smpd/nano2%SCALE_FACTOR)
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
                !$omp end do
            endif
            ! remove unused atoms from the pdb file
            !$omp do collapse(1) schedule(static) private(i)
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
            !$omp end do
            call centers_close1%writepdb  (trim(nano1%fbody)//'_atom_close_couples')
            call centers_close2%writepdb  (trim(nano2%fbody)//'_atom_close_couples')
            call centers_coupled1%writepdb(trim(nano1%fbody)//'_atom_couples')
            call centers_coupled2%writepdb(trim(nano2%fbody)//'_atom_couples')
            call couples1%writepdb('extra_atoms')
            !Avg dist and stdev symmetry breaking atoms from the center
            !Max dist atoms from the center
            avg = 0.
            cnt3 = 0
            m(:) = nano1%nanopart_masscen()
            m(:) = (m(:)-1.)*nano1%smpd
            tmp_max = 0.
            !$omp do collapse(1) schedule(static) private(i)
            do i = 1, N_max
                coord(:) = centers_coupled1%get_coord(i)
                if(coord(1)>TINY .and. coord(2)>TINY .and. coord(3)>TINY) then
                    cnt3 = cnt3 + 1
                    avg = avg + euclid(coord,m)
                    if(euclid(coord,m) > tmp_max) tmp_max = euclid(coord,m)
                endif
            enddo
            !$omp end do
            avg = avg/real(cnt3)
            cnt3 = 0
            stdev = 0.
            !$omp do collapse(1) schedule(static) private(i)
            do i = 1, N_max
                coord(:) = centers_coupled1%get_coord(i)
                if(coord(1)>TINY .and. coord(2)>TINY .and. coord(3)>TINY) then
                    cnt3 = cnt3 + 1
                    stdev = stdev + (euclid(coord,m)-avg)**2
                endif
            enddo
            !$omp end do
            stdev = sqrt(stdev/(real(cnt3)-1.))
            write(unit = 121, fmt = '(a,f6.3,a)')'AVG     DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', avg, ' A'
            write(unit = 121, fmt = '(a,f6.3,a)')'STDEV   DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', stdev, ' A'
            write(unit = 121, fmt = '(a,f6.3,a)')'MAX     DIST ATOMS THAT BREAK THE SYMMETRY TO THE CENTER: ', tmp_max, ' A'
            tmp_max = 0. ! reset
            !$omp do collapse(1) schedule(static) private(i)
            do i = 1, size(nano1%centers, dim = 2)
                coord(:) = (nano1%centers(:,i)-1.)*nano1%smpd
                    d =  euclid(coord,m)
                    if(d > tmp_max) tmp_max = d
            enddo
            !$omp end do
            write(unit = 121, fmt = '(a,f6.3,a)')'MAX     DIST ATOMS IN VOL1 3D RECONSTRUCT. TO THE CENTER: ', tmp_max, ' A'
            ! kill atoms instances
            call centers_close1%kill
            call centers_close2%kill
            call centers_coupled1%kill
            call centers_coupled2%kill
            call couples1%kill
            !RMSD
            rmsd = sqrt(sum(dist_sq)/real(count(dist_sq > TINY)))
            dist_no_zero = pack(dist, dist>TINY)
            dist_no_zero = dist_no_zero*nano1%smpd ! report distances in Amstrongs
            call hist(dist_no_zero, 50)
            ! !Output on file, Matlab Compatible
            ! open(119, file='RMSDhist')
            ! write (119,*) 'r=[...'
            ! do i = 1, size(dist_no_zero)
            !     write (119,'(A)', advance='no') trim(real2str(dist_no_zero(i)))
            !     if(i < size(dist_no_zero)) write (119,'(A)', advance='no') ', '
            ! end do
            ! write (119,*) '];'
            ! close(119)
            !For SCV files
            open(117, file='RMSDhist')
            write (117,*) 'r'
            do i = 1, size(dist_no_zero)
                write (117,'(A)', advance='yes') trim(real2str(dist_no_zero(i)))
            end do
            close(117)
            write(unit = 121, fmt = '(i3,a,i2,a)') count(dist_no_zero <= 0.5),' atoms correspond withing      50 pm. (',count(dist_no_zero <= 0.5)*100/(N_min), '% of the atoms )'
            write(unit = 121, fmt = '(i3,a,i2,a)')  cnt-(N_max-N_min),        ' atoms have error bigger than 220 pm. (', (cnt-N_max+N_min)*100/N_min, '% of the atoms )'
            write(unit = 121, fmt = '(a,f6.3,a)') 'RMSD CALCULATED CONSIDERING ALL ATOMS   = ', rmsd*nano1%smpd, ' A'
            write(unit = 121, fmt = '(a,f6.3,a)') 'RMSD ATOMS THAT CORRESPOND WITHIN 50 PM = ', (sqrt(sum(dist_close)/real(count(dist_close > TINY))))*nano1%smpd, ' A'
            if(present(r)) r=rmsd
            deallocate(dist, dist_sq, dist_no_zero, mask)
        end subroutine atomic_position_rmsd
    end subroutine compare_atomic_models

    ! This subroutine plots the histogram of the within-atoms
    ! distances distribution within the nanoparticle nano.
    ! To each atom the distance assigned is the min distance
    ! to the other atoms. There is a threshols (3A) for
    ! outliers discarding.
    ! If coords in input, then it considers just the atom-to-atom
    ! distance between the atoms with centers in coords.
    subroutine distances_distribution(self,coords,volume)
        use simple_atoms, only : atoms
        use simple_stat,  only : mad_gau
        class(nanoparticle), intent(inout) :: self
        real,    optional,      intent(in) :: coords(:,:)
        integer, optional,      intent(in) :: volume
        real, allocatable :: dist(:)
        integer :: i, j, n_discard
        real    :: stdev
        real    :: med, robust_stdev
        logical :: mask(self%n_cc)
        ! Initialisations
        mask = .true.
        stdev      = 0.
        n_discard  = 0
        if(present(coords)) then
            allocate(dist(size(coords,dim=2)), source = 0.)
            !$omp do collapse(1) schedule(static) private(i)
            do i = 1, size(coords,dim=2)
                dist(i) =  pixels_dist(coords(:,i), self%centers(:,:), 'min', mask) !I have to use all the atoms when
                mask(:) = .true. ! restore
                !Discard outliers
                if(dist(i)*self%smpd > 3.0 ) then !2.5 is the average interatomic distance
                    dist(i) = 0.
                    n_discard = n_discard + 1
                endif
            enddo
            !$omp end do
            self%avg_dist_atoms = sum(dist)/real(size(coords,dim=2)-n_discard)
            !$omp do collapse(1) schedule(static) private(i)
            do i = 1, size(coords,dim=2)
                if(dist(i)*self%smpd <=3.0) stdev = stdev + (dist(i)-self%avg_dist_atoms)**2
            enddo
            !$omp end do
            stdev = sqrt(stdev/real(size(coords,dim=2)-1-n_discard))
            med = median(dist)
            robust_stdev = mad_gau(dist,med)
            if(present(volume)) then
                write(unit = 11, fmt = '(a,a,f6.3,a)') 'Average dist atoms vol ', trim(int2str(volume)), self%avg_dist_atoms*self%smpd, ' A'
                write(unit = 11, fmt = '(a,a,a,f6.3,a)') 'StDev   dist atoms vol ', trim(int2str(volume)),':', stdev*self%smpd, ' A'
                write(unit = 11, fmt = '(a,a,a,f6.3,a)') 'Median  dist atoms vol', trim(int2str(volume)),':', med*self%smpd, ' A'
                write(unit = 11, fmt = '(a,a,a,f6.3,a)') 'MadGau  dist atoms vol ', trim(int2str(volume)),':', robust_stdev*self%smpd, ' A'
            else
                write(unit = 11, fmt = '(a,f6.3,a)') 'Average dist atoms: ', self%avg_dist_atoms*self%smpd, ' A'
                write(unit = 11, fmt = '(a,f6.3,a)') 'StDev   dist atoms: ', stdev*self%smpd, ' A'
                write(unit = 11, fmt = '(a,f6.3,a)') 'Median  dist atoms: ', med*self%smpd, ' A'
                write(unit = 11, fmt = '(a,f6.3,a)') 'MadGau  dist atoms: ', robust_stdev*self%smpd,  ' A'
            endif
            deallocate(dist)
        else
            ! Report statistics and images in dedicated directory
            !coma back to root folder
            call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: distances_distribution, simple_chdir; ")
            call simple_mkdir(trim(self%output_dir)//'/ClusterDistDistr',errmsg="simple_nanoparticles :: distances_distribution, simple_mkdir; ")
            call simple_chdir(trim(self%output_dir)//'/ClusterDistDistr',errmsg="simple_nanoparticles :: distances_distribution, simple_chdir; ")
            open(15, file='DistancesDistr')
            allocate(self%dists(size(self%centers, dim = 2)), source = 0.)
            !$omp do collapse(1) schedule(static) private(i)
            do i = 1, size(self%centers, dim = 2)
                self%dists(i) =  pixels_dist(self%centers(:,i), self%centers(:,:), 'min', mask) !I have to use all the atoms when
                mask(:) = .true. ! restore
                !Discard outliers
                if(self%dists(i)*self%smpd > 3.0 ) then !2.5 is the average interatomic distance
                    self%dists(i) = 0.
                    n_discard = n_discard + 1
                endif
            enddo
            !$omp end do
            self%avg_dist_atoms = sum(self%dists)/real(size(self%centers, dim = 2)-n_discard)
            !$omp do collapse(1) schedule(static) private(i)
            do i = 1, size(self%centers, dim = 2)
                if(self%dists(i)*self%smpd <=3.0) stdev = stdev + (self%dists(i)-self%avg_dist_atoms)**2
            enddo
            !$omp end do
            stdev = sqrt(stdev/real(size(self%centers, dim = 2)-1-n_discard))
            med = median(self%dists)
            robust_stdev = mad_gau(self%dists,med)
            write(unit = 15, fmt = '(a,f6.3,a)') 'Average dist atoms: ', self%avg_dist_atoms*self%smpd, ' A'
            write(unit = 15, fmt = '(a,f6.3,a)') 'StDev   dist atoms: ', stdev*self%smpd, ' A'
            write(unit = 15, fmt = '(a,f6.3,a)') 'Median  dist atoms: ', med*self%smpd, ' A'
            write(unit = 15, fmt = '(a,f6.3,a)') 'MadGau  dist atoms: ', robust_stdev*self%smpd, ' A'
            do i = 1,self%n_cc
                call self%centers_pdb%set_beta(i, self%dists(i)*self%smpd)
            enddo
            call self%centers_pdb%writepdb('DistancesDistr')
        endif
        close(15)
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
        !$omp do collapse(2) schedule(static) private(i,j)
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
        !$omp end do
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
    ! It also plot the hist of the interatomic distances
    ! distribution.
    ! SLOWWWW!!
    subroutine find_centers(self, img_bin, img_cc, coords)
       use simple_atoms, only: atoms
       class(nanoparticle),          intent(inout) :: self
       type(image), optional,        intent(inout) :: img_bin, img_cc
       real, optional, allocatable,  intent(out)   :: coords(:,:)
       logical, allocatable :: mask(:)
       integer, allocatable :: sz(:) !cc size to discard bottom 10/15%
       real,    allocatable :: rmat_cc(:,:,:)
       real,    allocatable :: rmat(:,:,:)
       integer :: i, n_discard, label(1)
       ! sanity check
       if(present(img_bin) .and. .not. present(img_cc)) THROW_HARD('img_bin and img_cc have to be both present in input')
       ! global variables allocation
       if(allocated(self%centers)) deallocate(self%centers)
       allocate(self%centers(3,self%n_cc),source = 0.)     !global variable
       !$omp do collapse(1) schedule(static) private(i)
       do i=1,self%n_cc
           if(present(img_bin)) then
               self%centers(:,i) = atom_masscen(self,i, img_bin, img_cc)
           else
               self%centers(:,i) = atom_masscen(self,i)
           endif
       enddo
       !$omp end do
       ! saving centers coordinates, optional
       if(present(coords)) allocate(coords(3, self%n_cc), source = self%centers)
       !sum_dist_closest_atom = 0. !initialise
       if(allocated(mask)) deallocate(mask)
       ! allocate(mask (self%n_cc), source = .true.)
       ! if(allocated(self%dists)) deallocate(self%dists)
       ! allocate(self%dists(self%n_cc), source = 0.)
       ! do i = 1, self%n_cc
       !     mask(:) = .true.  ! restore
       !     mask(i) = .false. ! not to consider the pixel itself, but I think its unnecessary
       !     self%dists(i) = pixels_dist(self%centers(:,i), self%centers,'min', mask)
       ! enddo
       !Output on file, Matlab Compatible
       ! open(113, file=trim(self%fbody)//'DistancesDistr')
       ! write (113,*) 'd=[...'
       ! do i = 1, size(self%dists)
       !     write (113,'(A)', advance='no') trim(real2str(self%dists(i)))
       !     if(i < size(self%dists)) write (113,'(A)', advance='no') ', '
       ! end do
       ! write (113,*) '];'
       ! close(113)
       ! self%avg_dist_atoms = sum(self%dists(:))/real(self%n_cc)
       ! deallocate(mask)
    contains

       !This function calculates the centers of mass of an
       !atom. It takes in input the image, its connected
       !component (cc) image and the label of the cc that
       !identifies the atom whose center of mass is to be
       !calculated. The result is stored in m.
       function atom_masscen(self, label, img_bin, img_cc) result(m)
           type(nanoparticle),    intent(inout) :: self
           integer,               intent(in)    :: label
           type(image), optional, intent(inout) :: img_bin, img_cc
           real(sp)             :: m(3)  !mass center coords
           real,    allocatable :: rmat_in(:,:,:), rmat_cc_in(:,:,:)
           integer :: i, j, k
           integer :: sz !sz of the cc identified by label
           if(present(img_bin)) then
               rmat_in    = img_bin%get_rmat()
               rmat_cc_in = img_cc%get_rmat()
           else
               rmat_in    = self%img_bin%get_rmat()
               rmat_cc_in = self%img_cc%get_rmat()
           endif
           where(     abs(rmat_cc_in-real(label)) > TINY) rmat_in = 0.
           sz = count(abs(rmat_cc_in-real(label)) < TINY)
           m = 0.
           !$omp do collapse(3) schedule(static) reduction(+:m) private(i,j,k)
           do i = 1, int(real(self%ldim(1))*self%SCALE_FACTOR)
               do j = 1, int(real(self%ldim(2))*self%SCALE_FACTOR)
                   do k = 1, int(real(self%ldim(3))*self%SCALE_FACTOR)
                       if(abs(rmat_in(i,j,k))> TINY) m = m + rmat_in(i,j,k)*[i,j,k]
                   enddo
               enddo
           enddo
           !$omp end do
           m = m/real(sz)
           deallocate(rmat_in, rmat_cc_in)
       end function atom_masscen
    end subroutine find_centers

    subroutine write_centers(self, fname)
        class(nanoparticle),        intent(inout) :: self
        character(len=*), optional, intent(in)    :: fname
        integer :: cc
        call self%centers_pdb%kill
        call self%centers_pdb%new(self%n_cc, dummy=.true.)
        !$omp do collapse(1) schedule(static) private(cc)
        do cc=1,self%n_cc
            call self%centers_pdb%set_coord(cc,(self%centers(:,cc)-1.)*self%smpd/self%SCALE_FACTOR)
        enddo
        !$omp end do
        if(present(fname)) then
            call self%centers_pdb%writepdb(fname)
        else
            call self%centers_pdb%writepdb(trim(self%fbody)//'_atom_centers')
        endif
    end subroutine write_centers

     ! This subroutine discard outliers that resisted binarization.
     ! It calculates the contact score of each atom and discards the bottom
     ! 5% of the atoms according to the contact score.
     ! It modifies the img_bin and img_cc instances deleting the
     ! identified outliers.
     ! Contact score: fix a radius and an atom A. Count the number N of atoms
     ! in the sphere centered in A of that radius. Define contact score of A
     ! equal to N. Do it for each atom.
     subroutine discard_outliers(self, volume)
         class(nanoparticle), intent(inout) :: self
         integer, optional,      intent(in) :: volume !volume identifier
         real    :: dist
         integer, allocatable :: contact_scores(:)
         real,    allocatable :: rmat(:,:,:), rmat_cc(:,:,:)
         real,    allocatable :: centers_polished(:,:)
         logical, allocatable :: mask(:)
         integer :: cc
         integer :: cnt      !contact_score
         integer :: loc(1)
         integer :: n_discard
         integer :: label(1)
         real    :: radius  !radius of the sphere to consider
         type(atoms) :: radial_atoms5A, radial_atoms7A, radial_atoms9A,radial_atoms12A
         integer :: cnt5, cnt7, cnt9, cnt12 ! nu,ber of atoms in radial shells
         real, allocatable :: coords(:,:), coord_no_zero(:,:) !coordinates of the centers of the atoms according to radial distances
         integer :: N    !number of atoms c3 map
         real    :: m(3) !mass center of c3 map
         real    :: d    !distance atoms from the center
         ! Size filter.
         call self%size_filtering()
         if(DEBUG_HERE) call self%img_bin%write('SizeFiltering.mrc')
         write(logfhandle, *) '****outliers discarding, init'
         ! Outliers removal using contact score
         radius = 2.*self%avg_dist_atoms
         allocate(contact_scores(self%n_cc), source = 0)
         if(allocated(mask)) deallocate(mask)
         allocate(mask(self%n_cc), source = .true.)
         !$omp do collapse(1) schedule(static) private(cc)
         do cc = 1, self%n_cc !fix the atom
             dist = 0.
             mask(1:self%n_cc) = .true.
             cnt = 0
             mask(cc)  = .false.
             do while(dist < radius)
                 dist = pixels_dist(self%centers(:,cc), self%centers,'min', mask, loc)
                 mask(loc) = .false.
                 cnt = cnt + 1
             enddo
             contact_scores(cc) = cnt - 1 !-1 because while loop counts one extra before exiting
         enddo
         !$omp end do
         n_discard = self%n_cc*5/100 ! discard bottom 10 instead%??
         mask(1:self%n_cc) = .true.
         rmat    = self%img_bin%get_rmat()
         rmat_cc = self%img_cc%get_rmat()
         ! Removing outliers from the binary image and the connected components image
         !$omp do collapse(1) schedule(static) private(cc)
         do cc = 1, n_discard
             label = minloc(contact_scores, mask)
             where(abs(rmat_cc-real(label(1))) < TINY)
                rmat = 0.   ! discard atom corresponding to label_discard(i)
             endwhere
             mask(label(1)) = .false.
             self%centers(:,label) = -1. !identify to discard them
         enddo
         !$omp end do
         call self%img_bin%set_rmat(rmat)
         call self%img_bin%find_connected_comps(self%img_cc)
         rmat_cc = self%img_cc%get_rmat()
         self%n_cc = nint(maxval(rmat_cc)) !update
         call self%find_centers()
         call self%write_centers()
         write(logfhandle, *) '****outliers discarding, completed'
         write(logfhandle, *) '****atom-to-atom distances estimation, init'
         ! Figure showing radial dependent symmetry
         N = self%n_cc
         call radial_atoms5A%new (N, dummy=.true.)
         call radial_atoms7A%new (N, dummy=.true.)
         call radial_atoms9A%new (N, dummy=.true.)
         call radial_atoms12A%new(N, dummy=.true.)
         cnt5  = 0
         cnt7  = 0
         cnt9  = 0
         cnt12 = 0
         m = self%nanopart_masscen()
         !$omp do collapse(1) schedule(static) private(cc)
         do cc = 1, self%n_cc
             d = euclid(self%centers(:,cc), m)*self%smpd
            if(d<=5.) then
                call radial_atoms5A%set_coord(cc,(self%centers(:,cc)-1.)*self%smpd/self%SCALE_FACTOR)
                cnt5 = cnt5+1
                cnt7 = cnt7+1
                cnt9 = cnt9+1
                cnt12 = cnt12+1
            elseif(d>5. .and. d<=7.) then
                cnt7 = cnt7+1
                cnt9 = cnt9+1
                cnt12 = cnt12+1
                call radial_atoms7A%set_coord(cc,(self%centers(:,cc)-1.)*self%smpd/self%SCALE_FACTOR)
            elseif(d>7. .and. d<=9.) then
                cnt9 = cnt9+1
                cnt12 = cnt12+1
                call radial_atoms9A%set_coord(cc,(self%centers(:,cc)-1.)*self%smpd/self%SCALE_FACTOR)
            elseif(d>9. .and. d<=12.) then
                cnt12 = cnt12+1
                call radial_atoms12A%set_coord(cc,(self%centers(:,cc)-1.)*self%smpd/self%SCALE_FACTOR)
            endif
         enddo
         !$omp end do
         ! Report statistics in dedicated directory
         !coma back to root folder
         call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: discard_outliers, simple_chdir; ")
         call simple_mkdir(trim(self%output_dir)//'/RadialDependentStat',errmsg="simple_nanoparticles :: discard_outliers, simple_mkdir; ")
         call simple_chdir(trim(self%output_dir)//'/RadialDependentStat',errmsg="simple_nanoparticles :: discard_outliers, simple_chdir; ")
         open(11, file='RadialDependentStat', position= 'append')
         if(present(volume)) then
             call radial_atoms5A%writepdb('vol'//int2str(volume)//'_radial_atoms5A')
             call radial_atoms7A%writepdb('vol'//int2str(volume)//'_radial_atoms7A')
             call radial_atoms9A%writepdb('vol'//int2str(volume)//'_radial_atoms9A')
             call radial_atoms12A%writepdb('vol'//int2str(volume)//'_radial_atoms12A')
         else
             call radial_atoms5A%writepdb('radial_atoms5A')
             call radial_atoms7A%writepdb('radial_atoms7A')
             call radial_atoms9A%writepdb('radial_atoms9A')
             call radial_atoms12A%writepdb('radial_atoms12A')
         endif
         ! Estimation of avg distance and stdev among atoms in radial dependent shells
         if(present(volume)) then
             write(unit = 11, fmt = '(a,a)') 'Estimation of atom-to-atom statistics in 5A radius shell vol', trim(int2str(volume))
         else
             write(unit = 11, fmt = '(a)') 'Estimation of atom-to-atom statistics in 5A radius shell'
         endif
         allocate(coords(3,cnt5), source = 0.)
         cnt = 0
         !$omp do collapse(1) schedule(static) private(cc)
         do cc = 1, self%n_cc
             d = euclid(self%centers(:,cc), m)*self%smpd
             if(d<=5.) then
                 cnt = cnt + 1
                 coords(:3,cnt) = self%centers(:,cc)
             endif
         enddo
         !$omp end do
         if(present(volume)) then
             call self%distances_distribution(coords, volume)
         else
             call self%distances_distribution(coords)
         endif
         deallocate(coords)
         if(present(volume)) then
             write(unit = 11, fmt = '(a,a)') 'Estimation of atom-to-atom statistics in 7A radius shell vol', trim(int2str(volume))
         else
             write(unit = 11, fmt = '(a)') 'Estimation of atom-to-atom statistics in 7A radius shell'
         endif
         allocate(coords(3,cnt7), source = 0.)
         cnt = 0
         !$omp do collapse(1) schedule(static) private(cc)
         do cc = 1, self%n_cc
             d = euclid(self%centers(:,cc), m)*self%smpd
             if(d<=7.) then
                 cnt = cnt + 1
                 coords(:3,cnt) = self%centers(:,cc)
             endif
         enddo
         !$omp end do
         if(present(volume)) then
             call self%distances_distribution(coords, volume)
         else
             call self%distances_distribution(coords)
         endif
         deallocate(coords)
         if(present(volume)) then
             write(unit = 11, fmt = '(a,a)') 'Estimation of atom-to-atom statistics in 9A radius shell vol', trim(int2str(volume))
         else
             write(unit = 11, fmt = '(a)') 'Estimation of atom-to-atom statistics in 9A radius shell'
         endif
         allocate(coords(3,cnt9), source = 0.)
         cnt = 0
         !$omp do collapse(1) schedule(static) private(cc)
         do cc = 1, self%n_cc
             d = euclid(self%centers(:,cc), m)*self%smpd
             if(d<=9.) then
                 cnt = cnt + 1
                 coords(:3,cnt) = self%centers(:,cc)
             endif
         enddo
         !$omp end do
         if(present(volume)) then
             call self%distances_distribution(coords, volume)
         else
             call self%distances_distribution(coords)
         endif
         deallocate(coords)
         if(present(volume)) then
             write(unit = 11, fmt = '(a,a)') 'Estimation of atom-to-atom statistics in 12A radius shell vol', trim(int2str(volume))
         else
             write(unit = 11, fmt = '(a)') 'Estimation of atom-to-atom statistics in 12A radius shell'
         endif
         allocate(coords(3,cnt12), source = 0.)
         cnt = 0
         !$omp do collapse(1) schedule(static) private(cc)
         do cc = 1, self%n_cc
             d = euclid(self%centers(:,cc), m)*self%smpd
             if(d<=12.) then
                 cnt = cnt + 1
                 coords(:3,cnt) = self%centers(:,cc)
             endif
         enddo
         !$omp end do
         if(present(volume)) then
             call self%distances_distribution(coords, volume)
         else
             call self%distances_distribution(coords)
         endif
         deallocate(coords)
         close(11)
         call radial_atoms5A%kill
         call radial_atoms7A%kill
         call radial_atoms9A%kill
         call radial_atoms12A%kill
         ! Come back to root directory
         !call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: discard_outliers, simple_chdir; ")
         call self%img_bin%write(trim(self%fbody)//'BIN.mrc')
         write(logfhandle, *) '****atom-to-atom distances estimation, completed'
         if(allocated(rmat))    deallocate(rmat)
         if(allocated(rmat_cc)) deallocate(rmat_cc)
     end subroutine discard_outliers

     ! calc the avg of the centers coords
     function nanopart_masscen(self) result(m)
         class(nanoparticle), intent(inout) :: self
         real    :: m(3)  !mass center coords
         integer :: i, j, k
         m = 0.
         !$omp do collapse(1) schedule(static) reduction(+:m) private(i)
         do i = 1, self%n_cc
              m = m + 1.*self%centers(:,i)
         enddo
         !$omp end do
         m = m/real(self%n_cc)
         if(self%ldim(3) == 1) m(3) = 0. !for 2D imgs
     end function nanopart_masscen

    ! This subroutine takes in input a connected component (cc) image
    ! and the label of one of its ccs and calculates its aspect ratio, which
    ! is defined as the ratio of the width and the height.
    ! The idea behind this is that the center of the cc is calculated,
    ! than everything is deleted except the borders of the cc. Finally,
    ! in order to calculate the width and the height, the min/max
    ! distances between the center and the borders are calculated. The
    ! aspect ratio is the ratio of those 2 distances.
    subroutine calc_aspect_ratio_private(self, label, ratio, print_ar)
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: label
        real,                intent(out)   :: ratio
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
            longest_dist = shortest_dist
            self%loc_longest_dist(:3, label) = self%centers(:,label)
        else
            longest_dist  = pixels_dist(self%centers(:,label), real(pos),'max', mask_dist, location)
            self%loc_longest_dist(:3, label) =  pos(:3,location(1))
        endif
        if(abs(longest_dist) > TINY) then
            ratio = shortest_dist/longest_dist
        else
            ratio = 0.
            if(DEBUG_HERE) write(logfhandle,*) 'cc ', label, 'LONGEST DIST = 0'
        endif
        if(present(print_ar) .and. (print_ar .eqv. .true.)) then
            write(logfhandle,*) 'ATOM #          ', label
            write(logfhandle,*) 'shortest dist = ', shortest_dist*self%smpd
            write(logfhandle,*) 'longest  dist = ', longest_dist*self%smpd
            write(logfhandle,*) 'RATIO         = ', ratio
        endif
        deallocate(imat_cc, border, pos, mask_dist)
    end subroutine calc_aspect_ratio_private

    subroutine calc_aspect_ratio(self, print_ar)
        use gnufor2
        class(nanoparticle), intent(inout) :: self
        logical, optional,   intent(in)    :: print_ar !print longest/shortest dim and ratio
        integer :: label
        integer, allocatable :: imat_cc(:,:,:)
        call self%img_bin%find_connected_comps(self%img_cc)
        imat_cc = int(self%img_cc%get_rmat())
        self%n_cc = int(maxval(imat_cc))
        deallocate(imat_cc)
        allocate(self%ratios(self%n_cc),             source = 0.)
        allocate(self%loc_longest_dist(3,self%n_cc), source = 0 )
        !$omp do collapse(1) schedule(static) private(label)
        do label = 1, self%n_cc
            call calc_aspect_ratio_private(self, label, self%ratios(label), print_ar)
        enddo
        !$omp end do
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
        !$omp do collapse(1) schedule(static) private(k)
        do k = 1, self%n_cc
            self%ang_var(k) = ang3D_vecs(vec_fixed(:),loc_ld_real(:,k))
            if(DEBUG_HERE) write(logfhandle,*) 'ATOM ', k, 'angle between direction longest dim and vec [0,0,1] ', self%ang_var(k)
        enddo
        !$omp end do
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
        integer,          allocatable :: imat(:,:,:), imat_onecls(:,:,:,:)
        type(image),      allocatable :: img_clusters(:)
        dim = self%n_cc
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
        ! Report clusters on images in dedicated directory
        call simple_chdir(trim(self%output_dir),errmsg="simple_nanoparticles :: affprop_cluster_dist_distr, simple_chdir1; ")
        !call simple_mkdir(trim(self%output_dir)//'/ClusterDistDistr',errmsg="simple_nanoparticles :: affprop_cluster_dist_distr, simple_mkdir; ")
        call simple_chdir(trim(self%output_dir)//'/ClusterDistDistr',errmsg="simple_nanoparticles :: affprop_cluster_dist_distr, simple_chdir; ")
        allocate(imat      (self%ldim(1),self%ldim(2),self%ldim(3)),        source = int(self%img_cc%get_rmat()))
        allocate(imat_onecls(self%ldim(1),self%ldim(2),self%ldim(3), ncls), source = 0)
        allocate(img_clusters(ncls))
        !$omp do collapse(2) schedule(static) private(i,j)
        do i = 1, self%n_cc
            do j = 1, ncls
                call img_clusters(j)%new(self%ldim,self%smpd)
                if(labels_ap(i) == j) then
                    where(imat==i) imat_onecls(:,:,:,j) = 1
                endif
            enddo
        enddo
        !$omp end do
        !$omp do collapse(1) schedule(static) private(j)
        do j = 1, ncls
            call img_clusters(j)%set_rmat(real(imat_onecls(:,:,:,j)))
            call img_clusters(j)%write(int2str(j)//'ClusterDistDistr.mrc')
            call img_clusters(j)%kill
        enddo
        !$omp end do
        if(allocated(img_clusters)) deallocate(img_clusters)
        if(allocated(imat_onecls))  deallocate(imat_onecls)
        if(allocated(imat))         deallocate(imat)
        open(125, file = 'ClusterDistDistr')
        write(unit = 125, fmt = '(a,i2.1)') 'NR OF CLUSTERS FOUND DISTS DISTR:', ncls
        do i = 1, self%n_cc
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

    !Affinity propagation clustering based on aspect ratios
    subroutine affprop_cluster_ar(self)
        use simple_aff_prop
        class(nanoparticle), intent(inout) :: self
        type(aff_prop)       :: ap_nano
        real, allocatable    :: simmat(:,:)
        real                 :: simsum
        integer, allocatable :: centers_ap(:), labels_ap(:)
        integer              :: i, j, ncls, nerr
        integer              :: dim
        integer,          allocatable :: imat(:,:,:), imat_onecls(:,:,:,:)
        type(image),      allocatable :: img_clusters(:)
        dim = size(self%ratios)
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
        allocate(imat      (self%ldim(1),self%ldim(2),self%ldim(3)),        source = int(self%img_cc%get_rmat()))
        allocate(imat_onecls(self%ldim(1),self%ldim(2),self%ldim(3), ncls), source = 0)
        allocate(img_clusters(ncls))
        !$omp do collapse(2) schedule(static) private(i,j)
        do i = 1, self%n_cc
            do j = 1, ncls
                call img_clusters(j)%new(self%ldim,self%smpd)
                if(labels_ap(i) == j) then
                    where(imat==i) imat_onecls(:,:,:,j) = 1
                endif
            enddo
        enddo
        !$omp end do
        !$omp do collapse(1) schedule(static) private(j)
        do j = 1, ncls
            call img_clusters(j)%set_rmat(real(imat_onecls(:,:,:,j)))
            call img_clusters(j)%write(int2str(j)//'ClusterAr.mrc')
            call img_clusters(j)%kill
        enddo
        !$omp end do
        if(allocated(img_clusters)) deallocate(img_clusters)
        if(allocated(imat_onecls))  deallocate(imat_onecls)
        if(allocated(imat))         deallocate(imat)
        open(119, file='ClusterAspectRatio')
        write(unit = 119,fmt ='(a,i2.1)') 'NR OF CLUSTERS FOUND AR:', ncls
        do i = 1, self%n_cc
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

    ! Affinity propagation clustering based on agles wrt vec [0,0,1].
    subroutine affprop_cluster_ang(self)
        use simple_aff_prop
        use simple_atoms, only : atoms
        class(nanoparticle), intent(inout) :: self
        type(aff_prop)       :: ap_nano
        real, allocatable    :: simmat(:,:)
        real                 :: simsum
        integer, allocatable :: centers_ap(:), labels_ap(:)
        integer              :: i, j, ncls, nerr
        integer              :: dim
        integer,          allocatable :: imat(:,:,:), imat_onecls(:,:,:,:)
        type(image),      allocatable :: img_clusters(:)
        real                 :: avg, stdev
        real,    allocatable :: avg_within(:), stdev_within(:)
        integer, allocatable :: cnt(:)
        dim = self%n_cc
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
        allocate(imat      (self%ldim(1),self%ldim(2),self%ldim(3)),        source = int(self%img_cc%get_rmat()))
        allocate(imat_onecls(self%ldim(1),self%ldim(2),self%ldim(3), ncls), source = 0)
        allocate(img_clusters(ncls))
        !$omp do collapse(2) schedule(static) private(i,j)
        do i = 1, self%n_cc
            do j = 1, ncls
                call img_clusters(j)%new(self%ldim,self%smpd)
                if(labels_ap(i) == j) then
                    where(imat==i) imat_onecls(:,:,:,j) = 1
                endif
            enddo
        enddo
        !$omp end do
        !$omp do collapse(1) schedule(static) private(j)
        do j = 1, ncls
            call img_clusters(j)%set_rmat(real(imat_onecls(:,:,:,j)))
            call img_clusters(j)%write(int2str(j)//'ClusterAng.mrc')
            call img_clusters(j)%kill
        enddo
        !$omp end do
        if(allocated(img_clusters)) deallocate(img_clusters)
        if(allocated(imat_onecls))  deallocate(imat_onecls)
        if(allocated(imat))         deallocate(imat)
        open(111, file='ClusterAngLongestDim')
        write(unit = 111,fmt ='(a,i2.1)') 'NR OF CLUSTERS FOUND ANG:', ncls
        do i = 1, self%n_cc
            call self%centers_pdb%set_beta(i, real(labels_ap(i)))
        enddo
        write(unit = 111,fmt ='(a)') 'CENTERS ANG'
        do i=1,ncls
            write(unit = 111,fmt ='(f6.3)') self%ang_var(centers_ap(i))
        end do
        !standard deviation within each center
        allocate(  avg_within(ncls), source = 0.)
        allocate(stdev_within(ncls), source = 0.)
        allocate(cnt(ncls), source = 0)
        !$omp do collapse(2) schedule(static) private(i,j)
        do i = 1, ncls
            do j = 1, self%n_cc
                if(labels_ap(j) == i) then
                    cnt(i) = cnt(i) + 1 !cnt is how many atoms I have in each class
                    avg_within(i) = avg_within(i) + self%ang_var(j)
            endif
            enddo
        enddo
        !$omp end do
        avg_within = avg_within/cnt
        write(unit = 111,fmt ='(a,i3)')     'number of atoms in each class = ', cnt(1)
        do i =2, ncls
            write(unit = 111,fmt ='(a,i3)') '                                ', cnt(i)
        enddo
        write(unit = 111,fmt ='(a,i3)')     'total atoms =                   ', sum(cnt)
        !$omp do collapse(2) schedule(static) private(i,j)
        do i = 1, ncls
            do j = 1, self%n_cc
                if(labels_ap(j) == i) stdev_within(i) = stdev_within(i) + (self%ang_var(j)-avg_within(i))**2
            enddo
        enddo
        !$omp end do
        stdev_within = stdev_within/(cnt-1.)
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
        write(unit = 111,fmt ='(a,f5.3)') 'standard deviation among centers: ', stdev
        call self%centers_pdb%writepdb('ClusterAngLongestDim')
        close(111)
        write(logfhandle,*) '****clustering wrt direction longest dim, completed'
        deallocate(simmat, labels_ap, centers_ap)
        if(allocated(cnt))          deallocate(cnt)
        if(allocated(avg_within))   deallocate(avg_within)
        if(allocated(stdev_within)) deallocate(stdev_within)
    end subroutine affprop_cluster_ang

    subroutine nanopart_cluster(self)
        class(nanoparticle), intent(inout) :: self
        ! clustering wrt to angle longest_dim-vector z=[0,0,1]
        call self%affprop_cluster_ang()
        ! clustering wrt aspect ratio
        call self%affprop_cluster_ar()
        ! clustering wrt interatomic distances distribution
        call self%affprop_cluster_dist_distr()
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
        !self%dim_over(3)       = 0
        call self%img%kill()
        call self%img_bin%kill()
        call self%img_cc%kill()
        call self%img_over_smp%kill()
        call self%centers_pdb%kill
        if(allocated(self%centers))          deallocate(self%centers)
        if(allocated(self%ratios))           deallocate(self%ratios)
        if(allocated(self%dists))            deallocate(self%dists)
        if(allocated(self%loc_longest_dist)) deallocate(self%loc_longest_dist)
        if(allocated(self%ang_var))          deallocate(self%ang_var)
    end subroutine kill_nanoparticle
end module simple_nanoparticles_mod
