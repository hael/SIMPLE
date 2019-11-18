module simple_bin_cc_image
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'

use simple_image,      only: image
implicit none

public :: binimage
private
#include "simple_local_flags.inc"


type, extends(image) :: binimage
    private
    real(kind=c_float), pointer :: brmat(:,:,:)=>null()  ! Pointer to rmat
    integer,        allocatable :: bimat(:,:,:)          ! Integer version, no Fourier extra indeces
    integer                     :: bdim(3)=0             ! Logical dimensions
    real                        :: bsmpd=0.              ! Sampling distance
    integer                     :: nccs=0                ! Number of connected components

  contains
    ! CONSTRUCTORS
    procedure          :: new_bimage
    procedure          :: copy_bimage
    ! SETTERS/GETTERS
    procedure          :: set_imat
    procedure          :: get_imat
    procedure          :: get_nccs

    ! procedure          :: get_imat
    ! I/O
    procedure          :: write_bimage
    procedure          :: read_bimage
    ! CONNECTED COMPONENTS
    procedure          :: find_ccs
    procedure          :: size_ccs
    procedure          :: elim_ccs
    procedure          :: order_ccs
    procedure          :: polish_ccs
    procedure          :: filling_holes
    procedure          :: diameter_ccs

    ! BINARY IMAGE METHODS
    procedure          :: grow_bins_bimage
    procedure          :: cos_edge_bimage
    procedure          :: border_msk
    procedure          :: ccs2bin

    ! MORPHOLOGICAL OPERATIONS
    procedure          :: erode
    procedure          :: dilatate

    ! NEIGHBORS
    procedure, private :: neigh_8_1
    procedure, private :: neigh_8_2
    generic            :: neigh_8 => neigh_8_1, neigh_8_2
    procedure          :: neigh_8_3D
    procedure, private :: neigh_4_3D_1
    procedure, private :: neigh_4_3D_2
    generic            :: neigh_4_3D => neigh_4_3D_1, neigh_4_3D_2

    ! KILL
    procedure          :: kill_bimage
end type

contains

    !>  \brief  Constructor for simple_image class
    subroutine new_bimage( self, ldim, smpd, wthreads )
        class(binimage),   intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd
        logical, optional, intent(in)    :: wthreads
        call self%kill_bimage
        call self%new(ldim, smpd, wthreads)
        call self%get_rmat_ptr(self%brmat)
        self%bdim  = self%get_ldim()
        self%bsmpd = self%get_smpd()
        allocate(self%bimat(1:self%bdim(1),1:self%bdim(2),1:self%bdim(3)))
        self%bimat = nint(self%brmat(1:self%bdim(1),1:self%bdim(2),1:self%bdim(3)))
    end subroutine new_bimage

    !>  \brief copy is a constructor that copies the input object
    subroutine copy_bimage( self, self_in )
        class(binimage), intent(inout) :: self
        class(binimage), intent(in)    :: self_in
        call self%new_bimage(self_in%bdim, self_in%bsmpd)
        self%bimat = self_in%bimat
    end subroutine copy_bimage

    ! SETTERS/GETTERS
    subroutine set_imat(self, imat)
      class(binimage),   intent(inout) :: self
      integer, optional, intent(in)    :: imat(:,:,:)
      if(present(imat)) then
        self%bimat = imat
        ! sanity check
        if(any(self%bdim /= shape(imat))) then
          write(logfhandle, *) 'self%bdim: ', self%bdim, 'shape(imat):', shape(imat)
          THROW_HARD('Uncompatible dimensions of the input matrix; set_imat')
        endif
      else   ! set imat to be equal to rmat
        self%bimat = nint(self%brmat(1:self%bdim(1),1:self%bdim(2),1:self%bdim(3)))
      endif
    end subroutine set_imat

    subroutine get_imat(self, imat)
      class(binimage),      intent(inout) :: self
      integer, allocatable, intent(inout) :: imat(:,:,:)
      if(allocated(imat)) deallocate(imat)
      allocate(imat(self%bdim(1),self%bdim(2),self%bdim(3)), source = self%bimat)
    end subroutine get_imat

    subroutine get_nccs(self, n)
      class(binimage), intent(inout)  :: self
      integer,         intent(out)    :: n
      n = maxval(self%bimat)
      ! update the instance of the class
      self%nccs = n
    end subroutine get_nccs

    ! Finds connected components in the binary input image and saves
    ! them in the output connected component image.
    ! Black = yes finds the black connected components instead of the white ones
    subroutine find_ccs(bimage, ccimage, black)
        class(binimage),   intent(inout) :: bimage
        class(binimage),   intent(inout) :: ccimage
        logical, optional, intent(in)    :: black
        type(binimage)       :: ccimage_unordered   !in this img will be stored the cc with no specific order
        integer, allocatable :: mat4compare(:,:,:), neigh_8(:)
        integer              :: i, j, k, n_it, n_maxit, nsz, cnt, diff, tmp
        logical              :: finished_job
        call ccimage_unordered%new_bimage(bimage%bdim,bimage%bsmpd)
        call ccimage%new_bimage          (bimage%bdim,bimage%bsmpd)
        if(bimage%bdim(3) > 1) then
            allocate(neigh_8(27), source = 0)
        else
            allocate(neigh_8(9),  source = 0)
        endif
        if(present(black) .and. black .eqv. .true.) then
            ! enumerate black pixels
            cnt     = 0 ! # labels
            n_maxit = 0
            do i = 1, bimage%bdim(1)
                do j = 1, bimage%bdim(2)
                    do k = 1, bimage%bdim(3)
                        if( bimage%bimat(i,j,k) < 1 )then
                            cnt = cnt + 1
                            ccimage_unordered%bimat(i,j,k) = cnt
                            n_maxit = max(cnt,n_maxit)
                        endif
                    enddo
                enddo
            enddo
            ! find connected components in parallel
            finished_job = .false.
            allocate(mat4compare(bimage%bdim(1),bimage%bdim(2),bimage%bdim(3)), source = 0)
            !$omp parallel default(shared) private(i,j,k,neigh_8,nsz) proc_bind(close)
            do n_it = 1, n_maxit
                if( .not. finished_job )then
                    !$omp workshare
                    mat4compare = ccimage_unordered%bimat
                    !$omp end workshare nowait
                    !$omp single
                    diff = 0
                    !$omp end single nowait
                    do i = 1, bimage%bdim(1)
                        do j = 1, bimage%bdim(2)
                            do k = 1, bimage%bdim(3)
                                if( bimage%bimat(i,j,k) < 1) then ! background
                                    if(bimage%bdim(3) > 1) then
                                        call ccimage_unordered%neigh_8_3D([i,j,k], neigh_8, nsz)
                                    else
                                        call ccimage_unordered%neigh_8  ([i,j,1], neigh_8, nsz)
                                    endif
                                    ccimage_unordered%bimat(i,j,k) = minval(neigh_8(:nsz), neigh_8(:nsz) > 0)
                                    diff = diff + abs(mat4compare(i,j,k) - ccimage_unordered%bimat(i,j,k))
                                endif
                            enddo
                        enddo
                    enddo
                    !$omp single
                    if( diff <= 0 ) finished_job = .true.
                    !$omp end single nowait
                endif
            enddo
            !$omp end parallel
            ! enumerate connected components
            cnt = 0
            do i = 1, bimage%bdim(1)
                do j = 1, bimage%bdim(2)
                    do k = 1, bimage%bdim(3)
                        if( ccimage_unordered%bimat(i,j,k) > 0 ) then  !rmat == 0  --> background
                            cnt = cnt + 1
                            tmp = ccimage_unordered%bimat(i,j,k)
                            where(ccimage_unordered%bimat == tmp)
                                ccimage%bimat = cnt
                                ccimage_unordered%bimat  = 0           !Not to consider this cc again
                            endwhere
                        endif
                    enddo
                enddo
            enddo
        else
            ! enumerate white pixels
            cnt     = 0 ! # labels
            n_maxit = 0
            do i = 1, bimage%bdim(1)
                do j = 1, bimage%bdim(2)
                    do k = 1, bimage%bdim(3)
                        if( bimage%bimat(i,j,k) > 0 )then
                            cnt = cnt + 1
                            ccimage_unordered%bimat(i,j,k) = cnt
                            n_maxit = max(cnt,n_maxit)
                        endif
                    enddo
                enddo
            enddo
            ! ALL GOOD UNITL HERE
            ! find connected components in parallel
            finished_job = .false.
            allocate(mat4compare(bimage%bdim(1),bimage%bdim(2),bimage%bdim(3)), source = 0)
            !$omp parallel default(shared) private(i,j,k,neigh_8,nsz) proc_bind(close)
            do n_it = 1, n_maxit
                if( .not. finished_job )then
                    !$omp workshare
                    mat4compare = ccimage_unordered%bimat
                    !$omp end workshare nowait
                    !$omp single
                    diff = 0
                    !$omp end single nowait
                    do i = 1, bimage%bdim(1)
                        do j = 1, bimage%bdim(2)
                            do k = 1, bimage%bdim(3)
                                if( bimage%bimat(i,j,k) > 0) then ! not background
                                    if(bimage%bdim(3) > 1) then
                                        call ccimage_unordered%neigh_8_3D([i,j,k], neigh_8, nsz)
                                    else
                                        call ccimage_unordered%neigh_8  ([i,j,1], neigh_8, nsz)
                                    endif
                                    ccimage_unordered%bimat(i,j,k) = minval(neigh_8(:nsz), neigh_8(:nsz) > 0)
                                    diff = diff + abs(mat4compare(i,j,k) - ccimage_unordered%bimat(i,j,k))
                                endif
                            enddo
                        enddo
                    enddo
                    !$omp single
                    if( diff <= 0 ) finished_job = .true.
                    !$omp end single nowait
                endif
            enddo
            !$omp end parallel
            ! enumerate connected components
            cnt = 0
            do i = 1, bimage%bdim(1)
                do j = 1, bimage%bdim(2)
                    do k = 1, bimage%bdim(3)
                        if( ccimage_unordered%bimat(i,j,k) > 0) then  !rmat == 0  --> background
                            cnt = cnt + 1
                            tmp = ccimage_unordered%bimat(i,j,k)
                            where(ccimage_unordered%bimat == tmp)
                                ccimage%bimat = cnt
                                ccimage_unordered%bimat  = 0          !Not to consider this cc again
                            endwhere
                        endif
                    enddo
                enddo
            enddo
        endif
        ccimage%nccs = maxval(ccimage%bimat)
        deallocate(mat4compare)
        call ccimage_unordered%kill_bimage
    end subroutine find_ccs

    ! The result of the function is the size(# of pixels) of each cc.
    !  (cc = connected component)
    function size_ccs(self) result(sz)
        class(binimage), intent(in) :: self
        integer, allocatable :: sz(:)
        integer :: n_cc,imax
        if(allocated(sz)) deallocate(sz)
        imax = maxval(self%bimat)
        allocate(sz(imax), source = 0)
        do n_cc = 1,imax
            sz(n_cc) = count(self%bimat == n_cc)
        enddo
    end function size_ccs

    ! This subroutine takes in input a connected component (cc) image
    ! and sets to 0 the cc which has size (# pixels) smaller than min_sz = range(1)
    ! or bigger than max_sz = range(2).
    ! It is created in order to prepare a micrograph for picking particles.
    subroutine elim_ccs(self, range)
        class(binimage), intent(inout) :: self
        integer,         intent(in)    :: range(2)
        integer, allocatable :: sz(:)
        integer              :: n_cc
        sz = self%size_ccs()
        do n_cc = 1, size(sz)  !for each cc
            if(sz(n_cc) < range(1) .or. sz(n_cc) > range(2)) then  !if the cc has size < min_sz or > max_sz
                where(self%bimat == n_cc)  !rmat == label
                    self%bimat = 0         !set to 0
                endwhere
            endif
        enddo
        ! re-oder cc
        call self%order_ccs()
        ! update number of connected components
        self%nccs = maxval(self%bimat)
        deallocate(sz)
    end subroutine elim_ccs

    !This function is meant to re-order the connected components (cc)
    !after some of them have been eliminated, so that they have an
    !increasing labelling (1,2,3..) with NO holes (NOT 1,2,5..).
    !Self is meant to be a connected component image.
    subroutine order_ccs(self)
        class(binimage), intent(inout) :: self
        integer, allocatable :: imat_aux(:,:,:)
        integer :: cnt, i
        cnt = 0
        allocate(imat_aux(self%bdim(1),self%bdim(2),self%bdim(3)), source = self%bimat)
        do i = 1, maxval(self%bimat)    !for each cc
            if(any(imat_aux == i)) then !there is cc labelled i
                cnt = cnt + 1           !increasing order cc
                where(imat_aux == i) self%bimat = cnt
            endif
        enddo
        deallocate(imat_aux)
    end subroutine order_ccs

    ! This subroutine takes in input a connected components (cc)
    ! image and eliminates some of the ccs according to their size.
    ! The decision method consists in calculate the avg size of the ccs
    ! and their standar deviation.
    ! Elimin ccs which have size: > ave + 2.*stdev
    !                             < ave - 2.*stdev
    ! If present min_rad and max_rad, it cleans up ccs more.
    subroutine polish_ccs(self, min_rad, max_rad)
        class(binimage),  intent(inout) :: self
        real, optional, intent(in)    :: min_rad, max_rad
        integer, allocatable :: sz(:)
        real    :: lt, ht !low and high thresh for ccs polising
        real    :: ave, stdev ! avg and stdev of the size od the ccs
        integer :: n_cc
        if(present(min_rad) .and. .not. present(max_rad)) THROW_HARD('Both min_rad and max_rad need to be defined! polish_ccs')
        if(present(max_rad) .and. .not. present(min_rad)) THROW_HARD('Both min_rad and max_rad need to be defined! polish_ccs')
        ! Assuming gaussian distribution 95% of the particles
        ! are in [-2sigma, 2sigma]
        sz = self%size_ccs()
        ave = sum(sz)/size(sz)
        stdev = 0.
        do n_cc = 1, size(sz)
            stdev = stdev + (real(sz(n_cc))-ave)**2
         enddo
        stdev = sqrt(stdev/real(size(sz)-1))
        call self%elim_ccs([ floor(ave-2.*stdev) , ceiling(ave+2.*stdev) ])
        !call img_cc%order_cc() it is already done in the elim_cc subroutine
        ! Use particle radius. Hypothesise the biggest possible area is when
        ! particle is circular (with rad 3*max_rad), the smallest one is when the particle is rectangular
        ! with lengths min_rad/2 and max_rad/2
        if(present(min_rad)) call self%elim_ccs([ int(min_rad*max_rad/4.) , int(2*3.14*(3*max_rad)**2) ])
        ! update number of connected components
        self%nccs = maxval(self%bimat)
    end subroutine polish_ccs

    ! This subroutine calculates the diamenter of the
    ! connected component labelled n_cc in the connected
    ! component image img_cc
    subroutine diameter_ccs(self, n_cc, diam)
        class(binimage), intent(inout) :: self
        integer,         intent(in)    :: n_cc
        real,            intent(out)   :: diam
        integer, allocatable :: pos(:,:)         ! position of the pixels of a fixed cc
        integer, allocatable :: imat_cc(:,:,:)   ! auxiliary
        logical, allocatable :: msk(:) ! For using function pixels_dist
        real  :: center_of_mass(3)     ! geometrical center of mass
        real  :: radius
        allocate(imat_cc(self%bdim(1),self%bdim(2),self%bdim(3)), source = self%bimat)
        where(imat_cc .ne. n_cc) imat_cc = 0 ! keep just the considered cc
        if(.not. any(imat_cc > 0)) THROW_HARD('Inputted non-existent cc')
        ! Find center of mass of the cc
        call get_pixel_pos(imat_cc,pos)
        center_of_mass(1) = sum(pos(1,:))/real(size(pos,dim = 2))
        center_of_mass(2) = sum(pos(2,:))/real(size(pos,dim = 2))
        center_of_mass(3) = 1.
        allocate(msk(size(pos, dim =2)), source = .true.)
        ! Calculate maximim radius
        radius = pixels_dist(center_of_mass,real(pos),'max',msk)
        ! Return diameter
        diam = 2.*radius
        deallocate(msk, pos, imat_cc)
    end subroutine diameter_ccs

    subroutine diameter_global( self, diam )
        class(binimage), intent(inout) :: self
        real,            intent(out)   :: diam
        integer :: i, ii, j, jj, k, kk, maxdistsq, distsq
        maxdistsq = 0
        !$omp parallel do collapse(3) default(shared) private(i,j,k,ii,jj,kk,distsq) schedule(static)&
        !$omp proc_bind(close) reduction(max:maxdistsq)
        do i=1,self%bdim(1)
            do j=1,self%bdim(2)
                do k=1,self%bdim(3)
                    if( self%bimat(i,j,k) < 1 ) cycle
                    do ii=1,self%bdim(1)
                        do jj=1,self%bdim(2)
                            do kk=1,self%bdim(3)
                                if( self%bimat(ii,jj,kk) < 1 ) cycle
                                distsq = sum(([i,j,k] - [ii,jj,kk])**2)
                                if( distsq > maxdistsq ) maxdistsq = distsq
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        diam = sqrt(real(maxdistsq))
    end subroutine diameter_global

    ! This subroutine modifies the input cc image by making it bin
    ! where just the cc n_cc is kept.
    subroutine ccs2bin(self, n_cc)
        class(binimage), intent(inout) :: self  !cc image
        integer,         intent(in)    :: n_cc  !label of the cc to keep
        where(self%bimat .ne. n_cc) self%bimat = 0
        where(self%bimat > 0)       self%bimat = 1
    end subroutine ccs2bin

    ! This subroutine is ment for 2D binary images. It implements
    ! the morphological operation dilatation.
    subroutine dilatate(self)
        class(binimage), intent(inout) :: self
        type(binimage) :: self_copy
        integer        :: neigh_8(3,8), i, j, k, nsz
        call self_copy%copy_bimage(self)
        if(self%bdim(3) /= 1) THROW_HARD('Non implemented for volumes!; dilatate')
        do i = 1, self%bdim(1)
            do j = 1, self%bdim(2)
                if(self_copy%bimat(i,j,1) == 1) then  !just for white pixels
                    call self%neigh_8([i,j,1], neigh_8, nsz)
                    do k = 1, nsz
                       self%bimat(neigh_8(1,k),neigh_8(2,k),1) = 1 !self%rmat(neigh_8) = 1
                    enddo
                endif
            enddo
        enddo
        call self_copy%kill_bimage
    end subroutine dilatate

    ! This subroutine is ment for 2D binary images. It implements
    ! the morphological operation erosion.
    subroutine erode(self, label)
        class(binimage),   intent(inout) :: self
        integer, optional, intent(in)    :: label
        logical, allocatable :: border(:,:,:)
        if(present(label)) then
            call self%border_mask(border,label, .true.)
        else
            call self%border_mask(border)
        endif
        where(border) self%bimat = 0
    end subroutine erode

    ! This subroutine builds the logical array 'border' of the same dims of
    ! the input image self. Border is true in corrispondence of
    ! the border pixels in self. Self is meant to be binary.
    ! It is necessary for erosion operation.
    ! Optional parameter neigh_4 decides whether to perform calculations
    ! using 4neighbours of 8neighbours. It's just for 3D volumes. It is
    ! useful in the case of very small objects, to make so that
    ! surface does not coincide con the volume itself.
    subroutine border_msk(self, border, label, four)
        class(binimage),      intent(in)    :: self
        logical, allocatable, intent(inout) :: border(:,:,:)
        integer, optional,    intent(in)    :: label
        logical, optional,    intent(in)    :: four
        integer, allocatable :: neigh_8(:)
        integer              :: i, j, k, nsz, llabel
        logical              :: ffour
        llabel = 1
        if(present(label)) llabel = label
        ffour = .false.
        if(present(four)) ffour = four
        if(present(four) .and. self%bdim(3)==1) THROW_HARD('4-neighbours identification hasn t been implemented for 2D images; border_msk')
        if(allocated(border)) deallocate(border)
        allocate(border(1:self%bdim(1), 1:self%bdim(2), 1:self%bdim(3)), source = .false.)
        if(self%bdim(3) == 1) then
            allocate(neigh_8(9), source = 0)
        elseif(.not. ffour) then
            allocate(neigh_8(27), source = 0)
        elseif(four) then
            allocate(neigh_8(6), source = 0) !actually this is neigh_4
        endif
        do i = 1,self%bdim(1)
            do j = 1, self%bdim(2)
                do k = 1, self%bdim(3)
                    if(self%bimat(i,j,k) == llabel ) then !white pixels
                        if(self%bdim(3) == 1) then
                            call self%neigh_8  ([i,j,k], neigh_8, nsz)
                        elseif(.not. ffour) then
                            call self%neigh_8_3D([i,j,k], neigh_8, nsz)
                        elseif(ffour) then
                            call self%neigh_4_3D([i,j,k], neigh_8, nsz)
                        endif
                        if(any(neigh_8(:nsz)<1)) border(i,j,k) = .true. !the image is supposed to be either binary or connected components image
                    endif
                enddo
            enddo
        enddo
        deallocate(neigh_8)
   end subroutine border_msk

   ! This subroutine fills the holes in a binary image
   ! Idea: A hole is a set of background pixels that cannot be
   ! reached by filling in the background from the edge of the image.
   ! Find the cc to which the background belongs. Whatever is not background,
   ! set it to foreground.
   subroutine filling_holes(self)
       class(binimage), intent(inout) :: self
       type (binimage)       :: img_cc ! connected component image
       real,    pointer     :: rmat(:,:,:)
       integer, allocatable :: imat_cc(:,:,:)
       integer              :: seed
       if(self%bdim(3) > 1) THROW_HARD('Not implemented for volumes! filling_holes')
       call self%find_ccs(img_cc, black=.true.) !detect the connnected components in the background as well
       imat_cc = nint(img_cc%get_rmat())
       seed = img_cc%bimat(1,1,1) !the pxl (1,1,1) should belong to the background
       ! Set to foreground whatever is not background
       where(img_cc%bimat == seed)
            self%bimat = 0
       elsewhere
            self%bimat = 1
       endwhere
       if(allocated(imat_cc)) deallocate(imat_cc)
       call img_cc%kill_bimage
     end subroutine filling_holes

     !> \brief grow_bins adds nlayers of pixels bordering the background in a binary image
     subroutine grow_bins_bimage( self, nlayers )
         class(binimage), intent(inout) :: self
         integer,         intent(in)    :: nlayers
         integer                        :: i,j,k, tsz(3,2), win(3,2), pdsz(3,2)
         logical, allocatable           :: add_pixels(:,:,:), template(:,:,:)
         tsz(:,1) = -nlayers
         tsz(:,2) = nlayers
         if(self%is_2d())tsz(3,:) = 1
         allocate( template(tsz(1,1):tsz(1,2), tsz(2,1):tsz(2,2), tsz(3,1):tsz(3,2)), stat=alloc_stat )
         if(alloc_stat/=0)call allocchk('grow_bins_bimage; 2')
         pdsz(:,1) = 1 - nlayers
         pdsz(:,2) = self%bdim + nlayers
         if(self%is_2d())pdsz(3,:) = 1
         allocate( add_pixels(pdsz(1,1):pdsz(1,2), pdsz(2,1):pdsz(2,2),&
             &pdsz(3,1):pdsz(3,2)), stat=alloc_stat )
         if(alloc_stat/=0)call allocchk('grow_bins_bimage; 1')
         ! template matrix
         template = .true.
         do i = tsz(1,1), tsz(1,2)
             do j = tsz(2,1), tsz(2,2)
                 if(self%is_2d())then
                     if(dot_product([i,j], [i,j]) > nlayers**2) template(i,j,1) = .false.
                 else
                     do k = tsz(3,1), tsz(3,2)
                         if(dot_product([i,j,k],[i,j,k]) > nlayers**2) template(i,j,k) = .false.
                     enddo
                 endif
             enddo
         enddo
         ! init paddedd logical array
         add_pixels = .false.
         forall( i=1:self%bdim(1), j=1:self%bdim(2), k=1:self%bdim(3), self%bimat(i,j,k)==1 )&
                 & add_pixels(i,j,k) = .true.
         ! cycle
         if( self%is_3d() )then
             do i = 1, self%bdim(1)
                 if( .not.any(self%bimat(i,:,:) > 0) )cycle
                 do j = 1, self%bdim(2)
                     if( .not.any(self%bimat(i,j,:) > 0) )cycle
                     win(1:2,1) = [i, j] - nlayers
                     win(1:2,2) = [i, j] + nlayers
                     do k = 1, self%bdim(3)
                         if (self%bimat(i,j,k) < 1)cycle
                         win(3,1) = k - nlayers
                         win(3,2) = k + nlayers
                         add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
                             &add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2))&
                             &.or.template
                     enddo
                 enddo
             enddo
         else
             do i=1,self%bdim(1)
                 if( .not.any(self%bimat(i,:,1) > 0) )cycle
                 do j=1,self%bdim(2)
                     win(1:2,1) = [i, j] - nlayers
                     win(1:2,2) = [i, j] + nlayers
                     if (self%bimat(i,j,1) < 1)cycle
                     add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), 1) =&
                         &add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), 1).or.template(:,:,1)
                 enddo
             enddo
         endif
         ! finalize
         self%bimat = 0
         forall( i=1:self%bdim(1), j=1:self%bdim(2), k=1:self%bdim(3), add_pixels(i,j,k) ) &
             & self%bimat(i,j,k) = 1
         deallocate( template, add_pixels )
     end subroutine grow_bins_bimage

     !>  \brief cos_edge applies cosine squared edge to a binary image
     subroutine cos_edge_bimage( self, cos_img, falloff )
         class(binimage), intent(inout) :: self    ! input
         type(image),     intent(inout) :: cos_img ! output
         integer,         intent(in)    :: falloff
         integer, allocatable :: imat(:,:,:)
         real,    pointer     :: rmat(:,:,:)
         real                 :: rfalloff, scalefactor
         integer              :: i, j, k, is, js, ks, ie, je, ke
         integer              :: il, ir, jl, jr, kl, kr, falloff_sq
         if( falloff<=0 ) THROW_HARD('stictly positive values for edge fall-off allowed; cos_edge')
         rfalloff    = real( falloff )
         falloff_sq  = falloff**2
         scalefactor = PI / rfalloff
         allocate( imat(self%bdim(1),self%bdim(2),self%bdim(3)),stat=alloc_stat )
         if(alloc_stat /= 0)call allocchk("In simple_image::cos_edge")
         imat = self%bimat
         call cos_img%new(self%bdim,self%bsmpd)
         call cos_img%get_rmat_ptr(rmat)
         rmat(1:self%bdim(1),1:self%bdim(2),1:self%bdim(3)) = real(imat)
         rmat = rmat/maxval(rmat(:self%bdim(1),:self%bdim(2),:self%bdim(3)))
         if( self%is_2d() )then
             ! 2d
             do i=1,self%bdim(1)
                 is = max(1,i-1)                  ! left neighbour
                 ie = min(i+1,self%bdim(1))       ! right neighbour
                 il = max(1,i-falloff)            ! left bounding box limit
                 ir = min(i+falloff,self%bdim(1)) ! right bounding box limit
                 if( .not. any(imat == 1) )cycle
                 do j=1,self%bdim(2)
                     if( imat(i,j,1) /= 1 ) cycle
                     js = max(1,j-1)
                     je = min(j+1,self%bdim(2))
                     if( any( imat(is:ie,js:je,1) < 1) )then
                         jl = max(1,j-falloff)
                         jr = min(j+falloff,self%bdim(2))
                         call update_mask_2d
                     endif
                 end do
             end do
         else
             ! 3d
             do i=1,self%bdim(1)
                 is = max(1,i-1)                  ! left neighbour
                 ie = min(i+1,self%bdim(1))       ! right neighbour
                 il = max(1,i-falloff)            ! left bounding box limit
                 ir = min(i+falloff,self%bdim(1)) ! right bounding box limit
                 if( .not. any(imat == 1) )cycle ! no values equal to one
                 do j=1,self%bdim(2)
                     js = max(1,j-1)
                     je = min(j+1,self%bdim(2))
                     jl = max(1,j-falloff)
                     jr = min(j+falloff,self%bdim(2))
                     if(.not. any(imat == 1) )cycle ! cycle if equal to one
                     do k=1,self%bdim(3)
                         if(imat(i,j,k) /= 1)cycle
                         ! within mask region
                         ks = max(1,k-1)
                         ke = min(k+1,self%bdim(3))
                         if( any( imat(is:ie,js:je,ks:ke) < 1) )then
                             ! update since has a masked neighbour
                             kl = max(1,k-falloff)
                             kr = min(k+falloff,self%bdim(3))
                             call update_mask_3d
                         endif
                     end do
                 end do
             end do
         endif
         deallocate(imat)
     contains

         !> updates neighbours with cosine weight
         subroutine update_mask_2d
             integer :: ii, jj, di_sq, dist_sq
             do ii=il,ir
                 di_sq = (ii-i)**2                 ! 1D squared distance in x dim
                 do jj=jl,jr
                     dist_sq = di_sq + (jj-j)**2   ! 2D squared distance in x & y dim
                     if(dist_sq > falloff_sq)cycle
                     ! masked neighbour
                     if( rmat(ii,jj,1)<1 )&
                         &rmat(ii,jj,1) = max(local_versine(real(dist_sq)), rmat(ii,jj,1))
                 enddo
             enddo
         end subroutine update_mask_2d

         !> updates neighbours with cosine weight
         subroutine update_mask_3d
             integer :: ii, jj, kk, di_sq, dij_sq, dist_sq
             do ii=il,ir
                 di_sq = (ii-i)**2
                 do jj=jl,jr
                     dij_sq = di_sq+(jj-j)**2
                     do kk=kl,kr
                         dist_sq = dij_sq + (kk-k)**2
                         if(dist_sq > falloff_sq)cycle
                         if( rmat(ii,jj,kk)<1 )&
                             &rmat(ii,jj,kk) = max(local_versine(real(dist_sq)), rmat(ii,jj,kk))
                     enddo
                 enddo
             enddo
         end subroutine update_mask_3d

         !> Local elemental cosine edge function
         elemental real function local_versine( r_sq )result( c )
             real, intent(in) :: r_sq
             c = 0.5 * (1. - cos(scalefactor*(sqrt(r_sq)-rfalloff)) )
         end function local_versine
     end subroutine cos_edge_bimage

    ! NEIGHBORS
    ! Returns 8-neighborhoods of the pixel position px in self
    ! it returns the INTENSITY values of the 8-neigh in a CLOCKWISE order, starting from any 4-neigh
    ! and the value of the pixel itself (the last one).
    subroutine neigh_8_1(self, px, neigh_8, nsz )
        class(binimage), intent(in)    :: self
        integer,         intent(in)    :: px(3)
        integer,         intent(inout) :: neigh_8(9)
        integer,         intent(out)   :: nsz
        integer :: i, j
        i = px(1)
        j = px(2) ! Assumes to have a 2-dim matrix
        ! identify neighborhood
        if( i-1 < 1 .and. j-1 < 1 )then                            ! NW corner
            neigh_8(1) = self%bimat(i+1,j,1)
            neigh_8(2) = self%bimat(i+1,j+1,1)
            neigh_8(3) = self%bimat(i,j+1,1)
            neigh_8(4) = self%bimat(i,j,1)
            nsz = 4
        else if (j+1 > self%bdim(2) .and. i+1 > self%bdim(1)) then ! SE corner
            neigh_8(1) = self%bimat(i-1,j,1)
            neigh_8(2) = self%bimat(i-1,j-1,1)
            neigh_8(3) = self%bimat(i,j-1,1)
            neigh_8(4) = self%bimat(i,j,1)
            nsz = 4
        else if (j-1 < 1  .and. i+1 >self%bdim(1)) then            ! SW corner
            neigh_8(3) = self%bimat(i-1,j,1)
            neigh_8(2) = self%bimat(i-1,j+1,1)
            neigh_8(1) = self%bimat(i,j+1,1)
            neigh_8(4) = self%bimat(i,j,1)
            nsz = 4
        else if (j+1 > self%bdim(2) .and. i-1 < 1) then            ! NE corner
            neigh_8(1) = self%bimat(i,j-1,1)
            neigh_8(2) = self%bimat(i+1,j-1,1)
            neigh_8(3) = self%bimat(i+1,j,1)
            neigh_8(4) = self%bimat(i,j,1)
            nsz = 4
        else if( j-1 < 1 ) then                                    ! N border
            neigh_8(5) = self%bimat(i-1,j,1)
            neigh_8(4) = self%bimat(i-1,j+1,1)
            neigh_8(3) = self%bimat(i,j+1,1)
            neigh_8(2) = self%bimat(i+1,j+1,1)
            neigh_8(1) = self%bimat(i+1,j,1)
            neigh_8(6) = self%bimat(i,j,1)
            nsz = 6
        else if ( j+1 > self%bdim(2) ) then                        ! S border
            neigh_8(1) = self%bimat(i-1,j,1)
            neigh_8(2) = self%bimat(i-1,j-1,1)
            neigh_8(3) = self%bimat(i,j-1,1)
            neigh_8(4) = self%bimat(i+1,j-1,1)
            neigh_8(5) = self%bimat(i+1,j,1)
            neigh_8(6) = self%bimat(i,j,1)
            nsz = 6
        else if ( i-1 < 1 ) then                                   ! W border
            neigh_8(1) = self%bimat(i,j-1,1)
            neigh_8(2) = self%bimat(i+1,j-1,1)
            neigh_8(3) = self%bimat(i+1,j,1)
            neigh_8(4) = self%bimat(i+1,j+1,1)
            neigh_8(5) = self%bimat(i,j+1,1)
            neigh_8(6) = self%bimat(i,j,1)
            nsz = 6
        else if ( i+1 > self%bdim(1) ) then                       ! E border
            neigh_8(1) = self%bimat(i,j+1,1)
            neigh_8(2) = self%bimat(i-1,j+1,1)
            neigh_8(3) = self%bimat(i-1,j,1)
            neigh_8(4) = self%bimat(i-1,j-1,1)
            neigh_8(5) = self%bimat(i,j-1,1)
            neigh_8(6) = self%bimat(i,j,1)
            nsz = 6
        else                                                     ! DEFAULT
            neigh_8(1) = self%bimat(i-1,j-1,1)
            neigh_8(2) = self%bimat(i,j-1,1)
            neigh_8(3) = self%bimat(i+1,j-1,1)
            neigh_8(4) = self%bimat(i+1,j,1)
            neigh_8(5) = self%bimat(i+1,j+1,1)
            neigh_8(6) = self%bimat(i,j+1,1)
            neigh_8(7) = self%bimat(i-1,j+1,1)
            neigh_8(8) = self%bimat(i-1,j,1)
            neigh_8(9) = self%bimat(i,j,1)
            nsz = 9
        endif
    end subroutine neigh_8_1

    ! Returns 8-neighborhoods of the pixel position px in self
    ! it returns the pixel INDECES of the 8-neigh in a CLOCKWISE order,
    ! starting from any 8-neigh. It doesn't consider the pixel itself.
    subroutine neigh_8_2(self, px, neigh_8, nsz)
        class(binimage), intent(in)   :: self
        integer,         intent(in)   :: px(3)
        integer,         intent(inout):: neigh_8(3,8)
        integer,         intent(out)  :: nsz
        integer :: i, j
        i = px(1)
        j = px(2)            !Assumes to have a 2-dim matrix
        if ( i-1 < 1 .and. j-1 < 1 ) then
            neigh_8(1:3,1) = [i+1,j,1]
            neigh_8(1:3,2) = [i+1,j+1,1]
            neigh_8(1:3,3) = [i,j+1,1]
            nsz = 3
        else if (j+1 > self%bdim(2) .and. i+1 > self%bdim(1)) then
            neigh_8(1:3,1) = [i-1,j,1]
            neigh_8(1:3,2) = [i-1,j-1,1]
            neigh_8(1:3,3) = [i,j-1,1]
            nsz = 3
        else if (j-1 < 1  .and. i+1 >self%bdim(1)) then
            neigh_8(1:3,3) = [i-1,j,1]
            neigh_8(1:3,2) = [i-1,j+1,1]
            neigh_8(1:3,1) = [i,j+1,1]
            nsz = 3
        else if (j+1 > self%bdim(2) .and. i-1 < 1) then
            neigh_8(1:3,1) = [i,j-1,1]
            neigh_8(1:3,2) = [i+1,j-1,1]
            neigh_8(1:3,3) = [i+1,j,1]
            nsz = 3
        else if( j-1 < 1 ) then
            neigh_8(1:3,5) = [i-1,j,1]
            neigh_8(1:3,4) = [i-1,j+1,1]
            neigh_8(1:3,3) = [i,j+1,1]
            neigh_8(1:3,2) = [i+1,j+1,1]
            neigh_8(1:3,1) = [i+1,j,1]
            nsz = 5
        else if ( j+1 > self%bdim(2) ) then
            neigh_8(1:3,1) = [i-1,j,1]
            neigh_8(1:3,2) = [i-1,j-1,1]
            neigh_8(1:3,3) = [i,j-1,1]
            neigh_8(1:3,4) = [i+1,j-1,1]
            neigh_8(1:3,5) = [i+1,j,1]
            nsz = 5
        else if ( i-1 < 1 ) then
            neigh_8(1:3,1) = [i,j-1,1]
            neigh_8(1:3,2) = [i+1,j-1,1]
            neigh_8(1:3,3) = [i+1,j,1]
            neigh_8(1:3,4) = [i+1,j+1,1]
            neigh_8(1:3,5) = [i,j+1,1]
            nsz = 5
        else if ( i+1 > self%bdim(1) ) then
            neigh_8(1:3,1) = [i,j+1,1]
            neigh_8(1:3,2) = [i-1,j+1,1]
            neigh_8(1:3,3) = [i-1,j,1]
            neigh_8(1:3,4) = [i-1,j-1,1]
            neigh_8(1:3,5) = [i,j-1,1]
            nsz = 5
        else
            neigh_8(1:3,1) = [i-1,j-1,1]
            neigh_8(1:3,2) = [i,j-1,1]
            neigh_8(1:3,3) = [i+1,j-1,1]
            neigh_8(1:3,4) = [i+1,j,1]
            neigh_8(1:3,5) = [i+1,j+1,1]
            neigh_8(1:3,6) = [i,j+1,1]
            neigh_8(1:3,7) = [i-1,j+1,1]
            neigh_8(1:3,8) = [i-1,j,1]
            nsz = 8
        endif
    end subroutine neigh_8_2


    ! Returns 8-neighborhoods (in 3D they are 27) of the pixel position px in self
    ! it returns the INTENSITY values of the 8-neigh in a CLOCKWISE order, starting from any 4-neigh
    ! of the first slice, then central slice and finally third slice.
    ! The value of the pixel itself is saved as the last one.
    ! This function is for volumes.
    subroutine neigh_8_3D(self, px, neigh_8, nsz )
        class(binimage), intent(in)    :: self
        integer,         intent(in)    :: px(3)
        integer,         intent(inout) :: neigh_8(27)
        integer,         intent(out)   :: nsz
        integer :: i, j, k
        logical :: i_okay, j_okay, k_okay
        i = px(1)
        j = px(2)
        k = px(3)
        i_okay = (i-1 > 0 .and. i+1 <= self%bdim(1))
        j_okay = (j-1 > 0 .and. j+1 <= self%bdim(2))
        k_okay = (k-1 > 0 .and. k+1 <= self%bdim(3))
        if( i-1 < 1 .and. j-1 < 1 .and. k-1 < 1 )then
            neigh_8(1) = self%bimat(i+1,j,k)
            neigh_8(2) = self%bimat(i+1,j+1,k)
            neigh_8(3) = self%bimat(i,j+1,k)
            neigh_8(4) = self%bimat(i+1,j,k+1)
            neigh_8(5) = self%bimat(i+1,j+1,k+1)
            neigh_8(6) = self%bimat(i,j+1,k+1)
            neigh_8(7) = self%bimat(i,j,k+1)
            neigh_8(8) = self%bimat(i,j,k)
            nsz = 8
            return
        elseif(i+1 > self%bdim(1) .and. j+1 > self%bdim(2) .and. k+1 > self%bdim(3) )then
            neigh_8(1) = self%bimat(i-1,j,k)
            neigh_8(2) = self%bimat(i-1,j-1,k)
            neigh_8(3) = self%bimat(i,j-1,k)
            neigh_8(4) = self%bimat(i-1,j,k-1)
            neigh_8(5) = self%bimat(i-1,j-1,k-1)
            neigh_8(6) = self%bimat(i,j-1,k-1)
            neigh_8(7) = self%bimat(i,j,k-1)
            neigh_8(8) = self%bimat(i,j,k)
            nsz = 8
            return
        elseif(i-1 < 1 .and. j-1 < 1 .and. k+1 > self%bdim(3) )then
            neigh_8(1) = self%bimat(i+1,j,k-1)
            neigh_8(2) = self%bimat(i+1,j+1,k-1)
            neigh_8(3) = self%bimat(i,j+1,k-1)
            neigh_8(4) = self%bimat(i,j,k-1)
            neigh_8(5) = self%bimat(i+1,j,k)
            neigh_8(6) = self%bimat(i+1,j+1,k)
            neigh_8(7) = self%bimat(i,j+1,k)
            neigh_8(8) = self%bimat(i,j,k)
            nsz = 8
            return
        elseif(i+1 > self%bdim(1) .and. j-1 < 1 .and. k-1 < 1 ) then
            neigh_8(1) = self%bimat(i-1,j,k)
            neigh_8(2) = self%bimat(i-1,j+1,k)
            neigh_8(3) = self%bimat(i,j+1,k)
            neigh_8(4) = self%bimat(i-1,j,k+1)
            neigh_8(5) = self%bimat(i-1,j+1,k+1)
            neigh_8(6) = self%bimat(i,j+1,k+1)
            neigh_8(7) = self%bimat(i,j,k+1)
            neigh_8(8) = self%bimat(i,j,k)
            nsz = 8
            return
        elseif(i-1 < 1 .and. j+1 > self%bdim(2) .and. k-1 < 1 ) then
            neigh_8(1) = self%bimat(i+1,j,k)
            neigh_8(2) = self%bimat(i+1,j-1,k)
            neigh_8(3) = self%bimat(i+1,j,k+1)
            neigh_8(4) = self%bimat(i,j,k+1)
            neigh_8(5) = self%bimat(i,j,k)
            neigh_8(6) = self%bimat(i+1,j-1,k+1)
            neigh_8(7) = self%bimat(i,j-1,k+1)
            neigh_8(8) = self%bimat(i,j-1,k)
            nsz = 8
            return
        elseif(i+1 > self%bdim(1) .and. j-1 < 1 .and. k+1 > self%bdim(3) ) then
            neigh_8(1) = self%bimat(i-1,j,k)
            neigh_8(2) = self%bimat(i-1,j+1,k)
            neigh_8(3) = self%bimat(i,j+1,k)
            neigh_8(4) = self%bimat(i-1,j,k-1)
            neigh_8(5) = self%bimat(i-1,j+1,k-1)
            neigh_8(6) = self%bimat(i,j+1,k-1)
            neigh_8(7) = self%bimat(i,j,k-1)
            neigh_8(8) = self%bimat(i,j,k)
            nsz = 8
            return
        elseif(i-1 < 1 .and. j+1 > self%bdim(2) .and. k+1 > self%bdim(3) ) then
            neigh_8(1) = self%bimat(i+1,j,k)
            neigh_8(2) = self%bimat(i+1,j-1,k)
            neigh_8(3) = self%bimat(i+1,j,k-1)
            neigh_8(4) = self%bimat(i,j,k-1)
            neigh_8(5) = self%bimat(i,j,k)
            neigh_8(6) = self%bimat(i+1,j-1,k-1)
            neigh_8(7) = self%bimat(i,j-1,k-1)
            neigh_8(8) = self%bimat(i,j-1,k)
            nsz = 8
            return
        elseif(i+1 > self%bdim(1) .and. j+1 > self%bdim(2) .and. k-1 < 1 ) then
            neigh_8(1) = self%bimat(i-1,j,k)
            neigh_8(2) = self%bimat(i-1,j-1,k)
            neigh_8(3) = self%bimat(i-1,j,k+1)
            neigh_8(4) = self%bimat(i,j,k+1)
            neigh_8(5) = self%bimat(i,j,k)
            neigh_8(6) = self%bimat(i-1,j-1,k+1)
            neigh_8(7) = self%bimat(i,j-1,k+1)
            neigh_8(8) = self%bimat(i,j-1,k)
            nsz = 8
            return
        elseif( i-1 < 1 .and. k-1 < 1 .and. j_okay )then
            neigh_8(1) = self%bimat(i,j-1,k)
            neigh_8(2) = self%bimat(i,j,k)
            neigh_8(3) = self%bimat(i,j+1,k)
            neigh_8(4) = self%bimat(i+1,j-1,k)
            neigh_8(5) = self%bimat(i+1,j,k)
            neigh_8(6) = self%bimat(i+1,j+1,k)
            neigh_8(7) = self%bimat(i,j-1,k+1)
            neigh_8(8) = self%bimat(i,j,k+1)
            neigh_8(9) = self%bimat(i,j+1,k+1)
            neigh_8(10) = self%bimat(i+1,j-1,k+1)
            neigh_8(11) = self%bimat(i+1,j,k+1)
            neigh_8(12) = self%bimat(i+1,j+1,k+1)
            nsz = 12
            return
        elseif( i-1 < 1 .and. k+1 > self%bdim(3) .and. j_okay )then
            neigh_8(1) = self%bimat(i,j-1,k)
            neigh_8(2) = self%bimat(i,j,k)
            neigh_8(3) = self%bimat(i,j+1,k)
            neigh_8(4) = self%bimat(i+1,j-1,k)
            neigh_8(5) = self%bimat(i+1,j,k)
            neigh_8(6) = self%bimat(i+1,j+1,k)
            neigh_8(7) = self%bimat(i,j-1,k-1)
            neigh_8(8) = self%bimat(i,j,k-1)
            neigh_8(9) = self%bimat(i,j+1,k-1)
            neigh_8(10) = self%bimat(i+1,j-1,k-1)
            neigh_8(11) = self%bimat(i+1,j,k-1)
            neigh_8(12) = self%bimat(i+1,j+1,k-1)
            nsz = 12
            return
        elseif( i+1 > self%bdim(1) .and. k+1 > self%bdim(3) .and. j_okay )then
            neigh_8(1) = self%bimat(i,j-1,k)
            neigh_8(2) = self%bimat(i,j,k)
            neigh_8(3) = self%bimat(i,j+1,k)
            neigh_8(4) = self%bimat(i-1,j-1,k)
            neigh_8(5) = self%bimat(i-1,j,k)
            neigh_8(6) = self%bimat(i-1,j+1,k)
            neigh_8(7) = self%bimat(i,j-1,k-1)
            neigh_8(8) = self%bimat(i,j,k-1)
            neigh_8(9) = self%bimat(i,j+1,k-1)
            neigh_8(10) = self%bimat(i-1,j-1,k-1)
            neigh_8(11) = self%bimat(i-1,j,k-1)
            neigh_8(12) = self%bimat(i-1,j+1,k-1)
            nsz = 12
            return
        elseif( i+1 > self%bdim(1) .and. k-1 < 1 .and. j_okay )then
            neigh_8(1) = self%bimat(i,j-1,k)
            neigh_8(2) = self%bimat(i,j,k)
            neigh_8(3) = self%bimat(i,j+1,k)
            neigh_8(4) = self%bimat(i-1,j-1,k)
            neigh_8(5) = self%bimat(i-1,j,k)
            neigh_8(6) = self%bimat(i-1,j+1,k)
            neigh_8(7) = self%bimat(i,j-1,k+1)
            neigh_8(8) = self%bimat(i,j,k+1)
            neigh_8(9) = self%bimat(i,j+1,k+1)
            neigh_8(10) = self%bimat(i-1,j-1,k+1)
            neigh_8(11) = self%bimat(i-1,j,k+1)
            neigh_8(12) = self%bimat(i-1,j+1,k+1)
            nsz = 12
            return
        elseif( j-1 < 1 .and. k-1 < 1 .and. i_okay )then
            neigh_8(1) = self%bimat(i-1,j,k)
            neigh_8(2) = self%bimat(i,j,k)
            neigh_8(3) = self%bimat(i+1,j,k)
            neigh_8(4) = self%bimat(i-1,j+1,k)
            neigh_8(5) = self%bimat(i,j+1,k)
            neigh_8(6) = self%bimat(i+1,j+1,k)
            neigh_8(7) = self%bimat(i-1,j,k+1)
            neigh_8(8) = self%bimat(i,j,k+1)
            neigh_8(9) = self%bimat(i+1,j,k+1)
            neigh_8(10) = self%bimat(i-1,j+1,k+1)
            neigh_8(11) = self%bimat(i,j+1,k+1)
            neigh_8(12) = self%bimat(i+1,j+1,k+1)
            nsz = 12
            return
        elseif( j+1 > self%bdim(2) .and. k-1 < 1 .and. i_okay )then
            neigh_8(1) = self%bimat(i-1,j,k)
            neigh_8(2) = self%bimat(i,j,k)
            neigh_8(3) = self%bimat(i+1,j,k)
            neigh_8(4) = self%bimat(i-1,j-1,k)
            neigh_8(5) = self%bimat(i,j-1,k)
            neigh_8(6) = self%bimat(i+1,j-1,k)
            neigh_8(7) = self%bimat(i-1,j,k+1)
            neigh_8(8) = self%bimat(i,j,k+1)
            neigh_8(9) = self%bimat(i+1,j,k+1)
            neigh_8(10) = self%bimat(i-1,j-1,k+1)
            neigh_8(11) = self%bimat(i,j-1,k+1)
            neigh_8(12) = self%bimat(i+1,j-1,k+1)
            nsz = 12
            return
        elseif( j+1 > self%bdim(2) .and. k+1 > self%bdim(3) .and. i_okay )then
            neigh_8(1) = self%bimat(i-1,j,k)
            neigh_8(2) = self%bimat(i,j,k)
            neigh_8(3) = self%bimat(i+1,j,k)
            neigh_8(4) = self%bimat(i-1,j-1,k)
            neigh_8(5) = self%bimat(i,j-1,k)
            neigh_8(6) = self%bimat(i+1,j-1,k)
            neigh_8(7) = self%bimat(i-1,j,k-1)
            neigh_8(8) = self%bimat(i,j,k-1)
            neigh_8(9) = self%bimat(i+1,j,k-1)
            neigh_8(10) = self%bimat(i-1,j-1,k-1)
            neigh_8(11) = self%bimat(i,j-1,k-1)
            neigh_8(12) = self%bimat(i+1,j-1,k-1)
            nsz = 12
            return
        elseif( j-1 < 1 .and. k+1 > self%bdim(3) .and. i_okay )then
            neigh_8(1) = self%bimat(i-1,j,k)
            neigh_8(2) = self%bimat(i,j,k)
            neigh_8(3) = self%bimat(i+1,j,k)
            neigh_8(4) = self%bimat(i-1,j+1,k)
            neigh_8(5) = self%bimat(i,j+1,k)
            neigh_8(6) = self%bimat(i+1,j+1,k)
            neigh_8(7) = self%bimat(i-1,j,k-1)
            neigh_8(8) = self%bimat(i,j,k-1)
            neigh_8(9) = self%bimat(i+1,j,k-1)
            neigh_8(10) = self%bimat(i-1,j+1,k-1)
            neigh_8(11) = self%bimat(i,j+1,k-1)
            neigh_8(12) = self%bimat(i+1,j+1,k-1)
            nsz = 12
            return
        elseif( i-1 < 1 .and. j-1 < 1 .and. k_okay )then                            ! NW corner
            neigh_8(1) = self%bimat(i+1,j,k-1)
            neigh_8(2) = self%bimat(i+1,j+1,k-1)
            neigh_8(3) = self%bimat(i,j+1,k-1)
            neigh_8(4) = self%bimat(i,j,k-1)
            neigh_8(5) = self%bimat(i+1,j,k)
            neigh_8(6) = self%bimat(i+1,j+1,k)
            neigh_8(7) = self%bimat(i,j+1,k)
            neigh_8(8) = self%bimat(i+1,j,k+1)
            neigh_8(9) = self%bimat(i+1,j+1,k+1)
            neigh_8(10) = self%bimat(i,j+1,k+1)
            neigh_8(11) = self%bimat(i,j,k+1)
            neigh_8(12) = self%bimat(i,j,k)
            nsz = 12
            return
        else if ( j+1 > self%bdim(2) .and. i+1 > self%bdim(1) .and. k_okay ) then ! SE corner
            neigh_8(1) = self%bimat(i-1,j,k-1)
            neigh_8(2) = self%bimat(i-1,j-1,k-1)
            neigh_8(3) = self%bimat(i,j-1,k-1)
            neigh_8(4) = self%bimat(i,j,k-1)
            neigh_8(5) = self%bimat(i-1,j,k)
            neigh_8(6) = self%bimat(i-1,j-1,k)
            neigh_8(7) = self%bimat(i,j-1,k)
            neigh_8(8) = self%bimat(i-1,j,k+1)
            neigh_8(9) = self%bimat(i-1,j-1,k+1)
            neigh_8(10) = self%bimat(i,j-1,k+1)
            neigh_8(11) = self%bimat(i,j,k+1)
            neigh_8(12) = self%bimat(i,j,k)
            nsz = 12
            return
        else if ( j-1 < 1  .and. i+1 >self%bdim(1) .and. k_okay ) then            ! SW corner
            neigh_8(1) = self%bimat(i,j+1,k-1)
            neigh_8(2) = self%bimat(i-1,j+1,k-1)
            neigh_8(3) = self%bimat(i-1,j,k-1)
            neigh_8(4) = self%bimat(i,j,k-1)
            neigh_8(5) = self%bimat(i,j+1,k)
            neigh_8(6) = self%bimat(i-1,j+1,k)
            neigh_8(7) = self%bimat(i-1,j,k)
            neigh_8(8) = self%bimat(i,j+1,k+1)
            neigh_8(9) = self%bimat(i-1,j+1,k+1)
            neigh_8(10) = self%bimat(i-1,j,k+1)
            neigh_8(11) = self%bimat(i,j,k+1)
            neigh_8(12) = self%bimat(i,j,k)
            nsz = 12
            return
        else if ( j+1 > self%bdim(2) .and. i-1 < 1 .and. k_okay ) then            ! NE corner
            neigh_8(1) = self%bimat(i,j-1,k-1)
            neigh_8(2) = self%bimat(i+1,j-1,k-1)
            neigh_8(3) = self%bimat(i+1,j,k-1)
            neigh_8(4) = self%bimat(i,j,k-1)
            neigh_8(5) = self%bimat(i,j-1,k)
            neigh_8(6) = self%bimat(i+1,j-1,k)
            neigh_8(7) = self%bimat(i+1,j,k)
            neigh_8(8) = self%bimat(i,j-1,k+1)
            neigh_8(9) = self%bimat(i+1,j-1,k+1)
            neigh_8(10) = self%bimat(i+1,j,k+1)
            neigh_8(11) = self%bimat(i,j,k+1)
            neigh_8(12) = self%bimat(i,j,k)
            nsz = 12
            return
        else if( j-1 < 1 .and. i_okay .and. k_okay ) then                                    ! N border
            neigh_8(1) = self%bimat(i+1,j,k-1)
            neigh_8(2) = self%bimat(i+1,j+1,k-1)
            neigh_8(3) = self%bimat(i,j+1,k-1)
            neigh_8(4) = self%bimat(i-1,j+1,k-1)
            neigh_8(5) = self%bimat(i-1,j,k-1)
            neigh_8(6) = self%bimat(i,j,k-1)
            neigh_8(7) = self%bimat(i+1,j,k)
            neigh_8(8) = self%bimat(i+1,j+1,k)
            neigh_8(9) = self%bimat(i,j+1,k)
            neigh_8(10) = self%bimat(i-1,j+1,k)
            neigh_8(11) = self%bimat(i-1,j,k)
            neigh_8(12) = self%bimat(i+1,j,k+1)
            neigh_8(13) = self%bimat(i+1,j+1,k+1)
            neigh_8(14) = self%bimat(i,j+1,k+1)
            neigh_8(15) = self%bimat(i-1,j+1,k+1)
            neigh_8(16) = self%bimat(i-1,j,k+1)
            neigh_8(17) = self%bimat(i,j,k+1)
            neigh_8(18) = self%bimat(i,j,k)
            nsz = 18
            return
        else if ( j+1 > self%bdim(2) .and. i_okay .and. k_okay ) then                        ! S border
            neigh_8(1) = self%bimat(i-1,j,k-1)
            neigh_8(2) = self%bimat(i-1,j-1,k-1)
            neigh_8(3) = self%bimat(i,j-1,k-1)
            neigh_8(4) = self%bimat(i+1,j-1,k-1)
            neigh_8(5) = self%bimat(i+1,j,k-1)
            neigh_8(6) = self%bimat(i,j,k-1)
            neigh_8(7) = self%bimat(i-1,j,k)
            neigh_8(8) = self%bimat(i-1,j-1,k)
            neigh_8(9) = self%bimat(i,j-1,k)
            neigh_8(10) = self%bimat(i+1,j-1,k)
            neigh_8(11) = self%bimat(i+1,j,k)
            neigh_8(12) = self%bimat(i-1,j,k+1)
            neigh_8(13) = self%bimat(i-1,j-1,k+1)
            neigh_8(14) = self%bimat(i,j-1,k+1)
            neigh_8(15) = self%bimat(i+1,j-1,k+1)
            neigh_8(16) = self%bimat(i+1,j,k+1)
            neigh_8(17) = self%bimat(i,j,k+1)
            neigh_8(18) = self%bimat(i,j,k)
            nsz = 18
            return
        else if ( i-1 < 1 .and. j_okay .and. k_okay  ) then                                   ! W border
            neigh_8(1) = self%bimat(i,j-1,k-1)
            neigh_8(2) = self%bimat(i+1,j-1,k-1)
            neigh_8(3) = self%bimat(i+1,j,k-1)
            neigh_8(4) = self%bimat(i+1,j+1,k-1)
            neigh_8(5) = self%bimat(i,j+1,k-1)
            neigh_8(6) = self%bimat(i,j,k-1)
            neigh_8(7) = self%bimat(i,j-1,k)
            neigh_8(8) = self%bimat(i+1,j-1,k)
            neigh_8(9) = self%bimat(i+1,j,k)
            neigh_8(10) = self%bimat(i+1,j+1,k)
            neigh_8(11) = self%bimat(i,j+1,k)
            neigh_8(12) = self%bimat(i,j-1,k+1)
            neigh_8(13) = self%bimat(i+1,j-1,k+1)
            neigh_8(14) = self%bimat(i+1,j,k+1)
            neigh_8(15) = self%bimat(i+1,j+1,k+1)
            neigh_8(16) = self%bimat(i,j+1,k+1)
            neigh_8(17) = self%bimat(i,j,k+1)
            neigh_8(18) = self%bimat(i,j,k)
            nsz = 18
            return
        else if ( i+1 > self%bdim(1) .and. j_okay .and. k_okay  ) then                       ! E border
            neigh_8(1) = self%bimat(i,j+1,k-1)
            neigh_8(2) = self%bimat(i-1,j+1,k-1)
            neigh_8(3) = self%bimat(i-1,j,k-1)
            neigh_8(4) = self%bimat(i-1,j-1,k-1)
            neigh_8(5) = self%bimat(i,j-1,k-1)
            neigh_8(6) = self%bimat(i,j,k-1)
            neigh_8(7) = self%bimat(i,j+1,k)
            neigh_8(8) = self%bimat(i-1,j+1,k)
            neigh_8(9) = self%bimat(i-1,j,k)
            neigh_8(10) = self%bimat(i-1,j-1,k)
            neigh_8(11) = self%bimat(i,j-1,k)
            neigh_8(12) = self%bimat(i,j+1,k+1)
            neigh_8(13) = self%bimat(i-1,j+1,k+1)
            neigh_8(14) = self%bimat(i-1,j,k+1)
            neigh_8(15) = self%bimat(i-1,j-1,k+1)
            neigh_8(16) = self%bimat(i,j-1,k+1)
            neigh_8(17) = self%bimat(i,j,k+1)
            neigh_8(18) = self%bimat(i,j,k)
            nsz = 18
            return
        else if ( k-1 < 1 .and. i_okay .and. j_okay ) then
            neigh_8(1) = self%bimat(i-1,j-1,k)
            neigh_8(2) = self%bimat(i-1,j-1,k+1)
            neigh_8(3) = self%bimat(i-1,j,k)
            neigh_8(4) = self%bimat(i-1,j+1,k)
            neigh_8(5) = self%bimat(i-1,j+1,k+1)
            neigh_8(6) = self%bimat(i-1,j,k+1)
            neigh_8(7) = self%bimat(i,j-1,k)
            neigh_8(8) = self%bimat(i+1,j-1,k)
            neigh_8(9) = self%bimat(i+1,j,k)
            neigh_8(10) = self%bimat(i+1,j+1,k)
            neigh_8(11) = self%bimat(i,j+1,k)
            neigh_8(12) = self%bimat(i,j-1,k+1)
            neigh_8(13) = self%bimat(i+1,j-1,k+1)
            neigh_8(14) = self%bimat(i+1,j,k+1)
            neigh_8(15) = self%bimat(i+1,j+1,k+1)
            neigh_8(16) = self%bimat(i,j+1,k+1)
            neigh_8(17) = self%bimat(i,j,k+1)
            neigh_8(18) = self%bimat(i,j,k)
            nsz = 18
            return
        else if ( k+1 > self%bdim(3) .and. i_okay .and. j_okay) then
            neigh_8(1) = self%bimat(i-1,j-1,k)
            neigh_8(2) = self%bimat(i-1,j-1,k-1)
            neigh_8(3) = self%bimat(i-1,j,k)
            neigh_8(4) = self%bimat(i-1,j+1,k)
            neigh_8(5) = self%bimat(i-1,j+1,k-1)
            neigh_8(6) = self%bimat(i-1,j,k-1)
            neigh_8(7) = self%bimat(i,j-1,k)
            neigh_8(8) = self%bimat(i+1,j-1,k)
            neigh_8(9) = self%bimat(i+1,j,k)
            neigh_8(10) = self%bimat(i+1,j+1,k)
            neigh_8(11) = self%bimat(i,j+1,k)
            neigh_8(12) = self%bimat(i,j-1,k-1)
            neigh_8(13) = self%bimat(i+1,j-1,k-1)
            neigh_8(14) = self%bimat(i+1,j,k-1)
            neigh_8(15) = self%bimat(i+1,j+1,k-1)
            neigh_8(16) = self%bimat(i,j+1,k-1)
            neigh_8(17) = self%bimat(i,j,k-1)
            neigh_8(18) = self%bimat(i,j,k)
            nsz = 18
            return
        else if(i_okay .and. j_okay .and. k_okay) then                                                   ! DEFAULT
            neigh_8(1) = self%bimat(i-1,j-1,k-1)
            neigh_8(2) = self%bimat(i,j-1,k-1)
            neigh_8(3) = self%bimat(i+1,j-1,k-1)
            neigh_8(4) = self%bimat(i+1,j,k-1)
            neigh_8(5) = self%bimat(i+1,j+1,k-1)
            neigh_8(6) = self%bimat(i,j+1,k-1)
            neigh_8(7) = self%bimat(i-1,j+1,k-1)
            neigh_8(8) = self%bimat(i-1,j,k-1)
            neigh_8(9) = self%bimat(i,j,k-1)
            neigh_8(10) = self%bimat(i-1,j-1,k)
            neigh_8(11) = self%bimat(i,j-1,k)
            neigh_8(12) = self%bimat(i+1,j-1,k)
            neigh_8(13) = self%bimat(i+1,j,k)
            neigh_8(14) = self%bimat(i+1,j+1,k)
            neigh_8(15) = self%bimat(i,j+1,k)
            neigh_8(16) = self%bimat(i-1,j+1,k)
            neigh_8(17) = self%bimat(i-1,j,k)
            neigh_8(18) = self%bimat(i-1,j-1,k+1)
            neigh_8(19) = self%bimat(i,j-1,k+1)
            neigh_8(20) = self%bimat(i+1,j-1,k+1)
            neigh_8(21) = self%bimat(i+1,j,k+1)
            neigh_8(22) = self%bimat(i+1,j+1,k+1)
            neigh_8(23) = self%bimat(i,j+1,k+1)
            neigh_8(24) = self%bimat(i-1,j+1,k+1)
            neigh_8(25) = self%bimat(i-1,j,k+1)
            neigh_8(26) = self%bimat(i,j,k+1)
            neigh_8(27) = self%bimat(i,j,k)
            nsz = 27
            return
        else
            write(logfhandle, *) 'i, j, k =  ', i, j, k
            THROW_HARD('Case not covered!; calc3D_neigh_8')
        endif
    end subroutine neigh_8_3D


    ! Returns 4-neighborhoods (in 3D they are 6) of the pixel position px in self
    ! it returns the INTENSITY values of the 8-neigh in a fixed order.
    ! The value of the pixel itself is NOT saved.
    ! This function is for volumes.
    subroutine neigh_4_3D_1(self, px, neigh_4, nsz )
        class(binimage), intent(in)    :: self
        integer,         intent(in)    :: px(3)
        integer,         intent(inout) :: neigh_4(6)
        integer,         intent(out)   :: nsz
        integer :: i, j, k
        i = px(1)
        j = px(2)
        k = px(3)
        if(i+1<self%bdim(1) .and. i-1>0 .and. j+1<self%bdim(2) .and. j-1>0 .and. k+1<self%bdim(3) .and. k-1>0) then
            neigh_4(1) = self%bimat(i,j,k+1)
            neigh_4(2) = self%bimat(i,j,k-1)
            neigh_4(3) = self%bimat(i,j+1,k)
            neigh_4(4) = self%bimat(i,j-1,k)
            neigh_4(5) = self%bimat(i+1,j,k)
            neigh_4(6) = self%bimat(i-1,j,k)
            nsz = 6
            return
        endif
        if( i == 1 .and. j == 1 .and. k == 1) then
            neigh_4(1) = self%bimat(i,j,k+1)
            neigh_4(2) = self%bimat(i,j+1,k)
            neigh_4(3) = self%bimat(i+1,j,k)
            nsz = 3
            return
        elseif( i == 1 .and. j == 1) then
            neigh_4(1) = self%bimat(i,j,k+1)
            neigh_4(2) = self%bimat(i,j,k-1)
            neigh_4(3) = self%bimat(i,j+1,k)
            neigh_4(4) = self%bimat(i+1,j,k)
            nsz = 4
            return
        elseif( i == 1 .and. k == 1) then
            neigh_4(1) = self%bimat(i,j,k+1)
            neigh_4(3) = self%bimat(i,j+1,k)
            neigh_4(3) = self%bimat(i,j-1,k)
            neigh_4(4) = self%bimat(i+1,j,k)
            nsz = 4
            return
        elseif( j == 1 .and. k == 1) then
            neigh_4(1) = self%bimat(i,j,k+1)
            neigh_4(2) = self%bimat(i,j+1,k)
            neigh_4(3) = self%bimat(i+1,j,k)
            neigh_4(4) = self%bimat(i-1,j,k)
            nsz = 4
            return
        endif
        if( i+1 == self%bdim(1) .and. j+1 == self%bdim(2) .and. k+1 == self%bdim(3)) then
            neigh_4(1) = self%bimat(i,j,k-1)
            neigh_4(2) = self%bimat(i,j-1,k)
            neigh_4(3) = self%bimat(i-1,j,k)
            nsz = 3
            return
        elseif( i+1 == self%bdim(1) .and. j+1 == self%bdim(2)) then
            neigh_4(1) = self%bimat(i,j,k+1)
            neigh_4(2) = self%bimat(i,j,k-1)
            neigh_4(3) = self%bimat(i,j-1,k)
            neigh_4(4) = self%bimat(i-1,j,k)
            nsz = 4
            return
        elseif( i+1 == self%bdim(1) .and. k+1 == self%bdim(3)) then
            neigh_4(1) = self%bimat(i,j,k-1)
            neigh_4(3) = self%bimat(i,j+1,k)
            neigh_4(3) = self%bimat(i,j-1,k)
            neigh_4(4) = self%bimat(i-1,j,k)
            nsz = 4
            return
        elseif( j+1 == self%bdim(2) .and. k+1 == self%bdim(3)) then
            neigh_4(1) = self%bimat(i,j,k-1)
            neigh_4(2) = self%bimat(i,j-1,k)
            neigh_4(3) = self%bimat(i+1,j,k)
            neigh_4(4) = self%bimat(i-1,j,k)
            nsz = 4
            return
        endif
    end subroutine neigh_4_3D_1

    ! Returns 4-neighborhoods (in 3D they are 6) of the pixel position px in self
    ! it returns the COORDINATES of the 8-neigh in a fixed order.
    ! The value of the pixel itself is NOT saved.
    ! This function is for volumes.
    subroutine neigh_4_3D_2(self, px, neigh_4, nsz )
        class(binimage), intent(in)    :: self
        integer,         intent(in)    :: px(3)
        integer,         intent(inout) :: neigh_4(3,6)
        integer,         intent(out)   :: nsz
        integer :: i, j, k
        i = px(1)
        j = px(2)
        k = px(3)
         if( i == 1 .and. j == 1 .and. k == 1) then
            neigh_4(1:3,1) = [i,j,k+1]
            neigh_4(1:3,2) = [i,j+1,k]
            neigh_4(1:3,3) = [i+1,j,k]
            nsz = 3
            return
         elseif( i == 1 .and. j == 1) then
            neigh_4(1:3,1) = [i,j,k+1]
            neigh_4(1:3,2) = [i,j,k-1]
            neigh_4(1:3,3) = [i,j+1,k]
            neigh_4(1:3,4) = [i+1,j,k]
            nsz = 4
            return
         elseif( i == 1 .and. k == 1) then
            neigh_4(1:3,1) = [i,j,k+1]
            neigh_4(1:3,3) = [i,j+1,k]
            neigh_4(1:3,3) = [i,j-1,k]
            neigh_4(1:3,4) = [i+1,j,k]
            nsz = 4
            return
         elseif( j == 1 .and. k == 1) then
            neigh_4(1:3,1) = [i,j,k+1]
            neigh_4(1:3,2) = [i,j+1,k]
            neigh_4(1:3,3) = [i+1,j,k]
            neigh_4(1:3,4) = [i-1,j,k]
            nsz = 4
            return
         endif
         if( i+1 == self%bdim(1) .and. j+1 == self%bdim(2) .and. k+1 == self%bdim(3)) then
            neigh_4(1:3,1) = [i,j,k-1]
            neigh_4(1:3,2) = [i,j-1,k]
            neigh_4(1:3,3) = [i-1,j,k]
            nsz = 3
            return
         elseif( i+1 == self%bdim(1) .and. j+1 == self%bdim(2)) then
            neigh_4(1:3,1) = [i,j,k+1]
            neigh_4(1:3,2) = [i,j,k-1]
            neigh_4(1:3,3) = [i,j-1,k]
            neigh_4(1:3,4) = [i-1,j,k]
            nsz = 4
            return
         elseif( i+1 == self%bdim(1) .and. k+1 == self%bdim(3)) then
            neigh_4(1:3,1) = [i,j,k-1]
            neigh_4(1:3,3) = [i,j+1,k]
            neigh_4(1:3,3) = [i,j-1,k]
            neigh_4(1:3,4) = [i-1,j,k]
            nsz = 4
            return
         elseif( j+1 == self%bdim(2) .and. k+1 == self%bdim(3)) then
            neigh_4(1:3,1) = [i,j,k-1]
            neigh_4(1:3,2) = [i,j-1,k]
            neigh_4(1:3,3) = [i+1,j,k]
            neigh_4(1:3,4) = [i-1,j,k]
            nsz = 4
            return
         endif
        neigh_4(1:3,1) = [i,j,k+1]
        neigh_4(1:3,2) = [i,j,k-1]
        neigh_4(1:3,3) = [i,j+1,k]
        neigh_4(1:3,4) = [i,j-1,k]
        neigh_4(1:3,5) = [i+1,j,k]
        neigh_4(1:3,6) = [i-1,j,k]
        nsz = 6
    end subroutine neigh_4_3D_2

    ! I/O subroutines

    subroutine write_bimage(self, fname)
      class(binimage),  intent(inout) :: self
      character(len=*), intent(in)    :: fname
      call self%set_rmat(real(self%bimat))
      call self%write(fname)
    end subroutine write_bimage

    subroutine read_bimage(self, fname)
      class(binimage),  intent(inout) :: self
      character(len=*), intent(in)    :: fname
      call self%read(fname)
      call self%set_imat()
    end subroutine read_bimage

    !>  \brief  is a destructor
    subroutine kill_bimage( self )
        class(binimage), intent(inout) :: self
        self%brmat=>null()
        if(allocated(self%bimat)) deallocate(self%bimat)
        call self%kill
    end subroutine kill_bimage
end module simple_bin_cc_image
