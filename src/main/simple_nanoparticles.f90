!USAGE : simple_exec prg=detect_atoms smpd=0.358 vol1='particle1.mrc' (optional nnn=40)
module simple_nanoparticles_mod
include 'simple_lib.f08'
use simple_image, only : image
use simple_picker_chiara, only: pixels_dist, get_pixel_pos
use simple_cmdline,       only: cmdline
use simple_parameters,    only: parameters

implicit none

private
public :: nanoparticle

#include "simple_local_flags.inc"

! module global constants
integer, parameter :: N_THRESH         = 20      !number of thresholds for binarization
real, parameter    :: MAX_DIST_CENTERS = 3.1     !to determin atoms belonging to the same row
logical, parameter :: DEBUGG           = .false. !for debuggin purposes

type :: nanoparticle
    private
    ! these image objects are part of the instance to avoid excessive memory re-allocations
    type(image) :: img, img_bin, img_cc
    type(image) :: img_rows
    integer     :: ldim(3)
    real        :: smpd
    real        :: nanop_mass_cen(3)                    !coordinates of the center of mass of the nanoparticle
    integer     :: n_cc   !number of atoms (connected components)
    real,    allocatable  :: centers(:,:)
    real,    allocatable  :: ratios(:)
    integer, allocatable  :: loc_longest_dist(:,:)   !for indentific of the vxl that determins the longest dim of the atom
    type(parameters)      :: p
    type(cmdline)         :: cline
    character(len=STDLEN) :: partname
  contains
    ! constructor
    procedure          :: new => new_nanoparticle
    ! getters/setters
    procedure          :: get_img
    procedure          :: set_img
    procedure          :: set_partname
    ! segmentation
    procedure          :: binarize => nanopart_binarization
    procedure          :: find_centers
    procedure, private :: nanopart_masscen
    procedure, private :: calc_aspect_ratio_private
    procedure          :: calc_aspect_ratio
    procedure, private :: atom_rows
    ! polarization search
    procedure          :: polar_masscenter
    procedure, private :: find_nearest_neigh
    procedure, private :: search_neigh_polarization
    ! clustering
    procedure          :: affprop_cluster
    ! running
    procedure, private :: run_nanoparticle_job
    procedure          :: run
    ! kill
    procedure, private ::  kill => kill_nanoparticle
end type nanoparticle

contains

    !constructor
     subroutine new_nanoparticle(self)
         class(nanoparticle), intent(inout) :: self
         integer :: nptcls
         call find_ldim_nptcls(self%partname,  self%ldim, nptcls, self%smpd)
         call self%img%new         (self%ldim, self%smpd)
         call self%img_bin%new     (self%ldim, self%smpd)
         call self%img_cc%new      (self%ldim, self%smpd) !don't think I need to initializate it
         call self%img_rows%new    (self%ldim, self%smpd) !don't know if need it
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
         rmat = self%img%get_rmat()
         allocate(rmat_t(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0.)
         step = maxval(rmat)/real(N_THRESH)
         do i = 1, N_THRESH-1
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
         if(DEBUGG) write(logfhandle,*)  'SELECTED THRESHOLD = ', seleted_t(1)
         where(rmat > seleted_t(1))
             rmat_t = 1.
         elsewhere
             rmat_t = 0.
         endwhere
         call self%img_bin%set_rmat(rmat_t)
         deallocate(rmat, rmat_t, sz)
         write(logfhandle,*) '****binarization, completed'
     end subroutine nanopart_binarization

     ! Find the centers coordinates of the atoms in the particle
     ! and save it in the global variable centers.
     subroutine find_centers(self)
         use simple_atoms, only: atoms
         class(nanoparticle),         intent(inout) :: self
         type(atoms) :: centers_pdb
         real, allocatable :: rmat_cc(:,:,:)
         integer :: i,label
         ! next lines are too expensive, work on that
         call self%img_bin%find_connected_comps(self%img_cc)
         call self%img_cc%elim_cc([5, 800]) !connected components clean up
         call self%img_cc%bin(0.5) !binarize to re-calculate the cc, to have them in order again
         self%img_bin = self%img_cc
         call self%img_cc%kill
         call self%img_bin%find_connected_comps(self%img_cc) !cc re-calculation
         if(DEBUGG) call self%img_cc%write('PolishedOrderedCC.mrc')
         rmat_cc   = self%img_cc%get_rmat()
         self%n_cc = nint(maxval(rmat_cc))
         ! global variables allocation
         allocate(self%centers         (3, self%n_cc) , source = 0.)     !global variable
         allocate(self%loc_longest_dist(3, self%n_cc) , source = 0 )     !global variable
         ! centers plot
         call centers_pdb%new(self%n_cc, dummy=.true.)
         do i=1,self%n_cc
             self%centers(:,i) = atom_masscen(self,i)
             call centers_pdb%set_coord(i,(self%centers(:,i)-1.)*self%smpd) !-1 because of the different indexing: starting from 0 or 1
         enddo
         call centers_pdb%writepdb('atom_centers')

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
             !ASK HANS HOW TO PARALLELISE IT
             !omp do collapse(3) reduction(+:m) private(i,j,k) reduce schedule(static)
             do i = 1, self%ldim(1)
                 do j = 1, self%ldim(2)
                     do k = 1, self%ldim(3)
                         if(abs(rmat_in(i,j,k))> TINY) m = m + rmat_in(i,j,k)*[i,j,k]
                     enddo
                 enddo
             enddo
             !omp end do
             m = m/real(sz)
             if(self%ldim(3) == 1) m(3) = 0. !for 2D imgs
         end function atom_masscen
     end subroutine find_centers

    ! calc the avg of the centers coords
    function nanopart_masscen(self) result(m)
        class(nanoparticle), intent(inout) :: self
        real             :: m(3)  !mass center coords
        integer :: i, j, k
        m = 0.
        do i = 1, self%n_cc
             m = m + 1.*self%centers(:,i)
        enddo
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
    subroutine calc_aspect_ratio_private(self, label, ratio)
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: label
        real,                intent(out)   :: ratio
        integer, allocatable :: imat(:,:,:)
        integer, allocatable :: pos(:,:)
        real,    allocatable :: rmat(:,:,:), rmat_cc(:,:,:)
        logical, allocatable :: border(:,:,:)
        logical, allocatable :: mask_dist(:) !for min and max dist calculation
        integer :: location(1) !location of the farest vxls of the atom from its center
        integer :: i
        real    :: shortest_dist, longest_dist
        rmat    = self%img%get_rmat()
        rmat_cc = self%img_cc%get_rmat()
        allocate(imat(self%ldim(1),self%ldim(2),self%ldim(3)),source = nint(rmat_cc)) !for function get_pixel_pos
        call self%img_cc%border_mask(border, label)
        where(border .eqv. .true.)
            imat = 1
        elsewhere
            imat = 0
        endwhere
        call get_pixel_pos(imat,pos)                       !pxls positions of the shell
        if(allocated(mask_dist)) deallocate(mask_dist)
        allocate(mask_dist(size(pos, dim = 2)), source = .true. )
        shortest_dist = pixels_dist(self%centers(:,label), real(pos),'min', mask_dist)
        longest_dist  = pixels_dist(self%centers(:,label), real(pos),'max', mask_dist, location)
        self%loc_longest_dist(:3, label) =  pos(:3,location(1))
        if(abs(longest_dist) > TINY) then
            ratio = shortest_dist/longest_dist
        else
            ratio = 0.
            if(DEBUGG) write(logfhandle,*) 'cc ', label, 'LONGEST DIST = 0'
        endif
        write(logfhandle,*) 'CC # ', label
        write(logfhandle,*) 'shortest dist = ', shortest_dist
        write(logfhandle,*) 'longest dist = ',   longest_dist
        write(logfhandle,*) 'RATIO = ', self%ratios(label)
        deallocate(rmat, rmat_cc, border, imat, pos, mask_dist)
    end subroutine calc_aspect_ratio_private

    subroutine calc_aspect_ratio(self)
        use gnufor2
        class(nanoparticle), intent(inout) :: self
        integer :: label
        allocate(self%ratios(self%n_cc), source = 0.)
        do label = 1, self%n_cc
            call calc_aspect_ratio_private(self, label, self%ratios(label))
        enddo
        call hist(self%ratios, 20)
        ! To dump some of the analysis on aspect ratios on file compatible
        ! with Matlab.
         open(123, file='Ratios', position='append')
         write (123,*) 'ratios=[...'
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
    subroutine atom_rows(self) !change the name
        use gnufor2
        class(nanoparticle), intent(inout) :: self
          integer, allocatable :: imat(:,:,:)
          real, allocatable    :: rmat_rows(:,:,:), rmat_cc(:,:,:), rmat_aux_cc(:,:,:)
          real, allocatable    :: ang_var(:), cc_in_the_row(:)
          real, allocatable    :: loc_ld_real(:,:)
          logical, allocatable :: flag(:)    !not to check more than once the same center
          logical, allocatable :: do_it(:)   ! not to do it more than once, cnt sometimes is not increased
          type(image) :: img_aux, img_aux_cc !to store the cc of just one row of atoms
          integer :: n_row                   !number of atoms rows
          integer :: i, j, h, k
          integer :: cnt, cnt_additional
          integer :: label_max, label1, label2
          integer :: m(3)
          real    :: label_i, label_j
          real    :: theta, mod_1, mod_2, dot_prod
          write(logfhandle,*)'atom rows analysis'
          call img_aux%new   (self%ldim, self%smpd)
          call img_aux_cc%new(self%ldim, self%smpd)
          allocate(flag(self%n_cc),  source = .true. )
          allocate(imat     (self%ldim(1),self%ldim(2),self%ldim(3)), source = 0 )
          allocate(rmat_rows(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0.)
          allocate(loc_ld_real(3,size(self%loc_longest_dist, dim=2)), source = 0. )
          rmat_cc = self%img_cc%get_rmat()
          cnt = 0
          cnt_additional = 0
          allocate(do_it(self%n_cc), source = .true.)  !initialise
          do i = 1, self%n_cc -1
              if(flag(i)) cnt = cnt + 1
              do j = 2, self%n_cc
                  if(j > i .and. flag(i) .and. flag(j)) then
                      if(sqrt((self%centers(1,i)-self%centers(1,j))**2+(self%centers(2,i)-self%centers(2,j))**2) < MAX_DIST_CENTERS) then
                         if(DEBUGG ) then
                             write(logfhandle,*)'i = ', i, '& j = ', j, 'PARALLEL,&
                                        &  tot dist = ', &
                                        & sqrt((self%centers(1,i)-self%centers(1,j))**2+(self%centers(2,i)-self%centers(2,j))**2+(self%centers(3,i)-self%centers(3,j))**2)
                         endif
                         imat(nint(self%centers(1,i)),&
                              & nint(self%centers(2,i)),nint(self%centers(3,i))) = cnt
                         imat(nint(self%centers(1,j)),&
                              & nint(self%centers(2,j)),nint(self%centers(3,j))) = cnt
                         flag(j) = .false.
                         ! Set rmat_rows=cnt for the whole atom (not just the center)
                         label_j = rmat_cc(nint(self%centers(1,j)),nint(self%centers(2,j)),nint(self%centers(3,j)))
                         label_i = rmat_cc(nint(self%centers(1,i)),nint(self%centers(2,i)),nint(self%centers(3,i)))
                         if(abs(label_j) > TINY) then
                             where(abs(rmat_cc-label_j)<TINY)
                                  rmat_rows = cnt
                             endwhere
                         endif
                      endif
                  endif
              enddo
              if(abs(label_i)>TINY) then
                  where(abs(rmat_cc-label_i)<TINY) rmat_rows = real(cnt)
              endif
              ! surf each row
              where(abs(rmat_rows - real(cnt)) > TINY) rmat_rows = 0.
              if(flag(i) .and. any(rmat_rows > TINY)) then
                  cnt_additional = cnt_additional + 1
                  call self%img_rows%set_rmat(rmat_rows)
                  call self%img_rows%write('FULLRow'//int2str(cnt_additional)//'xy.mrc')
              endif
          enddo
          call self%img_rows%set_rmat(rmat_rows)
          call self%img_rows%write('ImgROWS.mrc')
          !take each atom row and analyse it for polarization
          do n_row = 1, cnt
              if(do_it(n_row)) then
                  do_it(n_row) = .false.
                  rmat_rows = self%img_rows%get_rmat()   ! restore rmat_rows
                  rmat_cc   = self%img_cc%get_rmat()     ! restore rmat_cc
                  where(abs(rmat_rows - real(n_row)) > TINY)
                      rmat_rows = 0.
                      rmat_cc = 0.
                  endwhere
                  call img_aux%set_rmat(rmat_cc)
                  call img_aux%find_connected_comps(img_aux_cc)
                  rmat_aux_cc = img_aux_cc%get_rmat()
                  cnt_additional = 0
                  label_max = nint(maxval(rmat_aux_cc)) ! cc which corresponds to the highest label in the row
                  !label_max has to be at least 4 to have significance (rows with enough number of atoms)
                  if(label_max > 3) then
                      allocate(ang_var      (label_max-1), source = 0.)
                      allocate(cc_in_the_row(label_max-1), source = 0.)
                      do k = 1, label_max-1
                          m = minloc(abs(rmat_aux_cc(:,:,:)-real(k))) !select one pxl labeled k
                          label1 = nint(rmat_cc(m(1),m(2),m(3)))      !identify the corresponding label in rmat_cc
                          ! calculate for k just the first time
                          if (k == 1) self%loc_longest_dist(:3,k) = self%loc_longest_dist(:3,k)-nint(self%centers(:3,label1))+1 !vxl identifying the longest dim of the atom translated into the origin
                          m = minloc(abs(rmat_aux_cc(:,:,:)-real(k+1))) !select one pxl labeled k+1
                          label2 = nint(rmat_cc(m(1),m(2),m(3)))
                          !translating initial point(self%centers) of the vector center-loc_longest_dist in the origin(1,1,1)
                          self%loc_longest_dist(:3,k+1) = self%loc_longest_dist(:3,k+1)-nint(self%centers(:3,label2))+1 !+1 because the origin is not (0,0,0)
                          loc_ld_real = real(self%loc_longest_dist) !translation in reals
                          dot_prod = loc_ld_real(1,k)*loc_ld_real(1,k+1)+ &
                          &          loc_ld_real(2,k)*loc_ld_real(2,k+1)+ &
                          &          loc_ld_real(3,k)*loc_ld_real(3,k+1)
                          mod_1 = sqrt(loc_ld_real(1,k  )**2+loc_ld_real(2,k  )**2+loc_ld_real(3,k  )**2)
                          mod_2 = sqrt(loc_ld_real(1,k+1)**2+loc_ld_real(2,k+1)**2+loc_ld_real(3,k+1)**2)
                          if(dot_prod/(mod_1*mod_2) > 1. .or. dot_prod/(mod_1*mod_2)< -1.) then
                              THROW_WARN('Out of the domain of definition of arccos; atom_rows')
                              theta = 0.
                          else
                              theta = acos(dot_prod/(mod_1*mod_2))
                          endif
                          write(logfhandle,*) 'neighbour' , k, 'angles in degrees = ', rad2deg(theta)
                          if(DEBUGG) then
                              write(logfhandle,*)'>>>>>>>>>>>>>>>>>>>>>>'
                              write(logfhandle,*)'dot_prod = ', dot_prod
                              write(logfhandle,*)'mod_1 = ', mod_1
                              write(logfhandle,*)'mod_2 = ', mod_2
                              write(logfhandle,*)'mod1*mod2=',(mod_1*mod_2), 'frac = ', dot_prod/(mod_1*mod_2)
                          endif
                          ang_var(k) = rad2deg(theta) !radians
                          cc_in_the_row(k) = real(k)
                      enddo
                      if(DEBUGG) then
                          write(logfhandle,*)'>>>>>>>>>>>>>>>>>>>>>>'
                          write(logfhandle,*)'N_ROW = ', n_row
                          write(logfhandle,*)'cc_in_the_row=', cc_in_the_row
                          write(logfhandle,*)'ang_var=', ang_var
                      endif
                      open(129, file='Polarization_xy',position = 'append')
                      write (129,*) 'Ang_vars'//int2str(n_row)//'=[...'
                      do k = 1, label_max-1
                              write (129,'(A)', advance='no') trim(real2str(cc_in_the_row(k)))
                              write (129,'(A)', advance='no') ', '
                              write (129,'(A)', advance='no') trim(real2str(ang_var(k)))
                              write (129,'(A)')'; ...'
                      end do
                      write (129,*) '];'
                      close(129)
                      deallocate(ang_var,cc_in_the_row)
                  endif
              endif
          enddo
          write(logfhandle,*)'atom rows analysis completed'
    end subroutine atom_rows

   ! This subroutine performs polarization search with respect
   ! to the center of mass of the nanoparticle.
    subroutine polar_masscenter(self, n)
        use simple_atoms, only: atoms
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: n !number of neighbours to consider
        integer            :: labels_out(n)
        integer            :: i, pdb_sz
        type(atoms)        :: masscen_pdb, neigh_pdb
        self%nanop_mass_cen = self%nanopart_masscen()
        pdb_sz = 1
        call masscen_pdb%new(pdb_sz, dummy=.true.)
        call masscen_pdb%set_coord(1,(self%nanop_mass_cen(:)-1.)*self%smpd) !-1 because of the different indexing: starting from 0 or 1
        call masscen_pdb%writepdb('MassCenter')
        call self%find_nearest_neigh(labels_out)
        pdb_sz = n
        call neigh_pdb%new(pdb_sz, dummy=.true.)
        do i=1,n
            call neigh_pdb%set_coord(i,(self%centers(:,labels_out(i))-1.)*self%smpd) !-1 because of the different indexing: starting from 0 or 1
        enddo
        call neigh_pdb%writepdb(int2str(n)//'NeighMassCenter')
        call self%search_neigh_polarization(labels_out)
    end subroutine polar_masscenter

    ! This subroutine finds the labels of the connected components
    ! whose centers are the m nearest neighbours of the atom
    ! identified by label_in and stores them in labels_out.
    ! The size of labels_out determines the numbers of neighbours
    ! to consider. Remember that for how it is built, self%centers(:3,i)
    ! contains the coords of the
    ! center of mass of the i-th cc.
    subroutine find_nearest_neigh(self, labels_out, label_in)
        class(nanoparticle),  intent(inout) :: self
        integer,              intent(inout) :: labels_out(:)
        integer, optional,    intent(in)    :: label_in      !If label in is NOT present, then it considers the neighbors of the center of mass
        integer :: neigh_nb
        integer :: cnt
        real    :: dist
        integer :: location(1)
        logical :: mask(self%n_cc)
        if(.not. allocated(self%centers)) THROW_HARD('Atom centers have to be calculated at first; find_nearest_neigh')
        if(size(labels_out) < 2) THROW_HARD('Too small number of neighbours; find_nearest_neigh')
        if(size(labels_out) > self%n_cc) THROW_HARD('Cannot look for too many neighbours; find_nearest_neigh')
        if(size(labels_out) == self%n_cc) then
            labels_out(:) = nint(self%centers(:,2))
            return
        endif
        !initialise
        mask = .true.                     ! every atom center could be neighbour
        neigh_nb = size(labels_out)       ! number of neighs to find
        if(.not. present(label_in)) then  ! find neighbous of center of mass
            do cnt = 1, neigh_nb
                dist =  pixels_dist(self%nanop_mass_cen(:), self%centers, 'min', mask, location)
                labels_out(cnt) = location(1)
                mask(location)  = .false. ! do not consider this atom again
            enddo
        else
            mask(label_in) = .false.      !do not consider the center of the selected atom
            do cnt = 1, neigh_nb
                dist =  pixels_dist(self%centers(:,label_in), self%centers, 'min', mask, location)
                labels_out(cnt) = location(1)
                mask(location)  = .false.
            enddo
        endif
    end subroutine find_nearest_neigh

    ! polarization search via angle variance. The considered angle is:
    ! -)  IF label_atom is present: the angle between the direction of
    !    the longest dim in the atom identified by label_atom and the
    !    direction of the longest dim of its neighbours, identified
    !    by labels_neigh
    ! -) IF label_atom is NOT present: the angle between the direction
    !    [1,1,1] and the direction of the longest dim of the atoms identified
    !    by labels_neigh, which are the neighbours of the masscenter.
    subroutine search_neigh_polarization(self, labels_neigh, label_atom)
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: labels_neigh(:)
        integer, optional,   intent(in)    :: label_atom   !if it is not present it calculates the angles wrt a fixed direction
        real, allocatable :: ang_var(:)
        real, allocatable :: loc_ld_real(:,:)
        integer :: nb_neigh
        integer :: k, i
        real    :: theta, mod_1, mod_2, dot_prod
        logical :: mask(size(labels_neigh)) !printing purposes
        nb_neigh = size(labels_neigh)
        if(allocated(ang_var)) deallocate(ang_var)
        allocate(ang_var(nb_neigh), source = 0.)
        if(allocated(loc_ld_real)) deallocate(loc_ld_real)
        allocate(loc_ld_real(3,size(self%loc_longest_dist, dim=2)), source = 0.)
        if(present(label_atom)) then
            self%loc_longest_dist(:3,label_atom) = self%loc_longest_dist(:3,label_atom)-nint(self%centers(:3,label_atom))+1 !vxl identifying the longest dim of the atom translated into the origin
            write(logfhandle,*)'>>>>> Variance angles wrt the direction of the longest dim in the aom identified by', label_atom
            do k = 1, nb_neigh
                self%loc_longest_dist(:3,labels_neigh(k)) = self%loc_longest_dist(:3,labels_neigh(k))-nint(self%centers(:3,labels_neigh(k)))+1 !vxl identifying the longest dim of the atom translated into the origin
                ! https://www.intmath.com/vectors/7-vectors-in-3d-space.php FORMULAE HERE
                loc_ld_real = real(self%loc_longest_dist) !translation in reals
                dot_prod = loc_ld_real(1,label_atom)*loc_ld_real(1,labels_neigh(k))+ &
                         & loc_ld_real(2,label_atom)*loc_ld_real(2,labels_neigh(k))+ &
                         & loc_ld_real(3,label_atom)*loc_ld_real(3,labels_neigh(k))
                mod_1 = sqrt(loc_ld_real(1,label_atom     )**2+loc_ld_real(2,label_atom     )**2+loc_ld_real(3,label_atom     )**2)
                mod_2 = sqrt(loc_ld_real(1,labels_neigh(k))**2+loc_ld_real(2,labels_neigh(k))**2+loc_ld_real(3,labels_neigh(k))**2)
                if(dot_prod/(mod_1*mod_2) > 1. .or. dot_prod/(mod_1*mod_2)< -1.) then
                    THROW_WARN('Out of the domain of definition of arccos; search_neigh_polarization')
                    theta = 0.
                else
                    theta = acos(dot_prod/(mod_1*mod_2)) !the output of acos in RADIANS
                endif
                write(logfhandle,*) 'neighbour' , k, 'angles in degrees = ', rad2deg(theta)
                if(DEBUGG) then
                    write(logfhandle,*)'>>>>>>>>>>>>>>>>>>>>>>'
                    write(logfhandle,*)'dot_prod = ', dot_prod
                    write(logfhandle,*)'mod_1 = ', mod_1
                    write(logfhandle,*)'mod_2 = ', mod_2
                    write(logfhandle,*)'mod1*mod2=',(mod_1*mod_2), 'frac = ', dot_prod/(mod_1*mod_2)
                endif
                ang_var(k) = rad2deg(theta) !radians
            enddo
            open(129, file='Polarization_NEIGH',position = 'append')
            write (129,*) 'ang_var'//int2str(label_atom)//'=[...'
            do k = 1, nb_neigh
                    write (129,'(A)', advance='no') trim(int2str(k))
                    write (129,'(A)', advance='no') ', '
                    write (129,'(A)', advance='no') trim(real2str(ang_var(k)))
                    write (129,'(A)')'; ...'
            end do
            write (129,*) '];'
            close(129)
        else !consider fixed direction [1,1,1] (it should be the direction of maximum atoms density)
            write(logfhandle,*)'>>>>>>>>>>>>>>>>>>>>>>variance angles wrt the direction [1,1,1]'
            do k = 1, nb_neigh
                self%loc_longest_dist(:3,labels_neigh(k)) = self%loc_longest_dist(:3,labels_neigh(k))-nint(self%centers(:3,labels_neigh(k)))+1
                loc_ld_real = real(self%loc_longest_dist) !casting
                dot_prod = 1.*loc_ld_real(1,labels_neigh(k))+ &
                         & 1.*loc_ld_real(2,labels_neigh(k))+ &
                         & 1.*loc_ld_real(3,labels_neigh(k))
                mod_1 = sqrt(3.) !modulus of vec [1,1,1]
                mod_2 = sqrt(loc_ld_real(1,labels_neigh(k))**2+loc_ld_real(2,labels_neigh(k))**2+loc_ld_real(3,labels_neigh(k))**2)
                if(dot_prod/(mod_1*mod_2) > 1. .or. dot_prod/(mod_1*mod_2)< -1.) then
                    THROW_WARN('Out of the domain of definition of arccos; search_neigh_polarization')
                    theta = 0.
                else
                    theta = acos(dot_prod/(mod_1*mod_2))
                endif
                write(logfhandle,*) 'neighbour' , k, 'angles in degrees = ', rad2deg(theta)
                if(DEBUGG) then
                    write(logfhandle,*)'>>>>>>>>>>>>>>>>>>>>>>'
                    write(logfhandle,*)'dot_prod = ', dot_prod
                    write(logfhandle,*)'mod_1 = ', mod_1
                    write(logfhandle,*)'mod_2 = ', mod_2
                    write(logfhandle,*)'mod1*mod2=',(mod_1*mod_2), 'frac = ', dot_prod/(mod_1*mod_2)
                endif
                ang_var(k) = rad2deg(theta) !radians
            enddo
            open(129, file='Polarization_MASSCEN',position = 'append')
            write (129,*) 'ang_var=[...'
            do k = 1, nb_neigh
                    write (129,'(A)', advance='no') trim(int2str(k))
                    write (129,'(A)', advance='no') ', '
                    write (129,'(A)', advance='no') trim(real2str(ang_var(k)))
                    write (129,'(A)')'; ...'
            end do
            write (129,*) '];'
            close(129)
        endif
        mask = .true. !not to overprint
        do i = 1,  nb_neigh-1
            do k = i+1,  nb_neigh
                if(abs(ang_var(i)-ang_var(k))< 0.05 .and. mask(i)) then
                    if(present(label_atom)) then
                        write(logfhandle,*)'Neighbours', i, 'and', k, 'of cc ',label_atom, 'have the same ang_var'
                        mask(k) = .false.
                    else
                        write(logfhandle,*)'Neighbours', i, 'and', k,'have the same ang_var wrt direction [1,1,1]'
                        mask(k) = .false.
                    endif
                endif
            enddo
        enddo
        ! TO SAVE THE COUPLES OF ANGLES
        ! coords(:) = nint(self%centers(:,labels_neigh(1)))
        ! rmat(coords(1),coords(2),coords(3)) = 1
        ! coords(:) = nint(self%centers(:,labels_neigh(9)))
        ! rmat(coords(1),coords(2),coords(3)) = 1
        ! coords(:) = nint(self%centers(:,labels_neigh(36)))
        ! rmat(coords(1),coords(2),coords(3)) = 1
        ! call img_couple%set_rmat(rmat)
        ! call img_couple%write('2Couples1_'//int2str(label_atom)//'.mrc')
        ! rmat = 0
        deallocate(ang_var, loc_ld_real)
    end subroutine search_neigh_polarization

    !Affinity propagation clustering based on aspect ratios
    subroutine affprop_cluster(self)
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
      write(logfhandle,*) 'NR OF CLUSTERS FOUND:', ncls
      write(logfhandle,*) 'CENTERS'
      do i=1,size(centers_ap)
          write(logfhandle,*) self%ratios(centers_ap(i))
      end do
      call img_clusters%new(self%ldim, self%smpd)
      rmat = img_clusters%get_rmat() !allocation
      do i = 1, self%n_cc
          rmat(nint(self%centers(1,i)),nint(self%centers(2,i)),nint(self%centers(3,i))) = labels_ap(i)
      enddo
      call img_clusters%set_rmat(rmat)
      call img_clusters%write('AspectRatioClusters.mrc')
      deallocate(simmat, labels_ap, centers_ap)
    end subroutine affprop_cluster

   ! To test it.
    subroutine run_nanoparticle_job(self)
      class(nanoparticle), intent(inout) :: self
      !binarization
       call self%img%read(self%partname)
       call self%binarize()
      !centers identifications
      call self%find_centers()
      !aspect ratios
      call self%calc_aspect_ratio()
      !polarization search wrt to the center of mass of the nanoparticle
      call self%polar_masscenter(10)
    end subroutine run_nanoparticle_job

    subroutine run(self)
        class(nanoparticle), intent(inout) :: self
        call self%cline%parse_oldschool()
        call self%cline%checkvar('smpd', 1)
        call self%cline%checkvar('vol1', 2)
        call self%cline%check()
        call self%p%new(self%cline)
        call self%new()
        call self%run_nanoparticle_job() !just one nanoparticle
        call self%kill()
    end subroutine run

    subroutine kill_nanoparticle(self)
        class(nanoparticle), intent(inout) :: self
        call self%img%kill()
        call self%img_bin%kill()
        call self%img_cc%kill()
        call self%img_rows%kill()
        deallocate(self%centers,self%ratios,self%loc_longest_dist)
    end subroutine kill_nanoparticle
end module simple_nanoparticles_mod
