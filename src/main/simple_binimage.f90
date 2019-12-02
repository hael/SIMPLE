module simple_binimage
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image, only: image
use simple_neighs
implicit none

public :: binimage
private
#include "simple_local_flags.inc"

type, extends(image) :: binimage
    private
    integer,  allocatable :: bimat(:,:,:)           ! Integer version, no extra Fourier pixels
    integer               :: bldim(3)     = 0       ! Logical dimension
    real                  :: bsmpd        = 0.      ! Sampling distance
    integer               :: nccs         = 0       ! # connected components
    logical               :: bimat_is_set = .false. ! bimat set flag
  contains
    ! CONSTRUCTORS
    procedure          :: new_bimg
    procedure          :: copy_bimg
    ! SETTERS/GETTERS
    procedure          :: set_imat
    procedure          :: get_imat
    procedure          :: get_nccs
    procedure          :: update_img_rmat ! for updating the rmat in the extended image class
    ! I/O
    procedure          :: write_bimg
    procedure          :: read_bimg
    ! CONNECTED COMPONENTS
    procedure          :: find_ccs
    procedure          :: size_ccs
    procedure          :: elim_ccs
    procedure          :: order_ccs
    procedure          :: polish_ccs
    procedure          :: fill_holes
    procedure          :: diameter_cc
    ! BINARY IMAGE METHODS
    procedure          :: diameter_bin
    procedure          :: grow_bins
    procedure          :: cos_edge
    procedure          :: border_mask
    procedure          :: cc2bin
    ! MORPHOLOGICAL OPERATIONS
    procedure          :: dilatate
    procedure          :: erode
    ! DESTRUCT
    procedure          :: kill_bimg
end type

contains

    ! CONSTRUCTORS

    subroutine new_bimg( self, ldim, smpd, wthreads )
        class(binimage),   intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd
        logical, optional, intent(in)    :: wthreads
        call self%kill_bimg
        call self%new(ldim, smpd, wthreads)
        self%bldim = self%get_ldim()
        self%bsmpd = self%get_smpd()
        allocate(self%bimat(self%bldim(1),self%bldim(2),self%bldim(3)), source=0)
        self%bimat_is_set = .false.
    end subroutine new_bimg

    subroutine copy_bimg( self, self_in )
        class(binimage),         intent(inout) :: self
        class(binimage), target, intent(in)    :: self_in
        real(kind=c_float), pointer :: brmat(:,:,:)=>null()
        call self%new_bimg(self_in%bldim, self_in%bsmpd)
        self%bimat = self_in%bimat
        call self%update_img_rmat
    end subroutine copy_bimg

    ! SETTERS/GETTERS

    subroutine set_imat(self, imat)
        class(binimage),   intent(inout) :: self
        integer, optional, intent(in)    :: imat(self%bldim(1),self%bldim(2),self%bldim(3))
        real(kind=c_float), pointer      :: brmat(:,:,:)=>null()
        if( present(imat) ) then
            self%bimat = imat
        else ! set imat to be equal to rmat
            call self%get_rmat_ptr(brmat)
            self%bimat = nint(brmat(1:self%bldim(1),1:self%bldim(2),1:self%bldim(3)))
            self%bimat_is_set = .true.
        endif
        call self%update_img_rmat
    end subroutine set_imat

    subroutine get_imat(self, imat)
        class(binimage),      intent(inout) :: self
        integer, allocatable, intent(inout) :: imat(:,:,:)
        if( allocated(imat) ) deallocate(imat)
        allocate(imat(self%bldim(1),self%bldim(2),self%bldim(3)), source=self%bimat)
    end subroutine get_imat

    subroutine get_nccs(self, n)
        class(binimage), intent(inout)  :: self
        integer,         intent(out)    :: n
        n = maxval(self%bimat)
        ! update the instance of the class
        self%nccs = n
    end subroutine get_nccs

    subroutine update_img_rmat( self )
        class(binimage), intent(inout) :: self
        real(kind=c_float), pointer :: brmat(:,:,:)=>null()
        call self%get_rmat_ptr(brmat)
        brmat(1:self%bldim(1),1:self%bldim(2),1:self%bldim(3)) = real(self%bimat)
    end subroutine update_img_rmat

    ! I/O

    subroutine write_bimg(self, fname)
        class(binimage),  intent(inout) :: self
        character(len=*), intent(in)    :: fname
        if( .not. self%bimat_is_set ) call self%set_imat
        call self%set_rmat(real(self%bimat))
        call self%write(fname)
    end subroutine write_bimg

    subroutine read_bimg(self, fname)
        class(binimage),  intent(inout) :: self
        character(len=*), intent(in)    :: fname
        call self%read(fname)
        call self%set_imat()
    end subroutine read_bimg

    ! CONNECTED COMPONENTS

    ! Finds connected components in the binary input image and saves
    ! them in the output connected component image.
    ! Black = .true. finds the black connected components instead of the white ones
    subroutine find_ccs(self, ccimage, black)
        class(binimage),   intent(inout) :: self
        class(binimage),   intent(inout) :: ccimage
        logical, optional, intent(in)    :: black
        type(binimage)       :: ccimage_unordered ! in this img will be stored the cc with no specific order
        integer, allocatable :: mat4compare(:,:,:), neigh_8_pixs(:)
        integer :: i, j, k, n_it, n_maxit, nsz, cnt, diff, tmp
        logical :: finished_job, black_present
        if( .not. self%bimat_is_set ) call self%set_imat
        black_present = present(black)
        call ccimage_unordered%new_bimg(self%bldim,self%bsmpd)
        call ccimage%new_bimg          (self%bldim,self%bsmpd)
        if(self%bldim(3) > 1) then
            allocate(neigh_8_pixs(27), source = 0)
        else
            allocate(neigh_8_pixs(9),  source = 0)
        endif
        if( black_present .and. black .eqv. .true.) then
            ! flip foreground to background and vice versa
            self%bimat = -1 * (self%bimat - 1)
        endif
        ! enumerate white pixels
        cnt     = 0 ! # labels
        n_maxit = 0
        do i = 1, self%bldim(1)
            do j = 1, self%bldim(2)
                do k = 1, self%bldim(3)
                    if( self%bimat(i,j,k) > 0 )then
                        cnt = cnt + 1
                        ccimage_unordered%bimat(i,j,k) = cnt
                        n_maxit = max(cnt,n_maxit)
                    endif
                enddo
            enddo
        enddo
        ! find connected components in parallel
        finished_job = .false.
        allocate(mat4compare(self%bldim(1),self%bldim(2),self%bldim(3)), source = 0)
        !$omp parallel default(shared) private(i,j,k,neigh_8_pixs,nsz) proc_bind(close)
        do n_it = 1, n_maxit
            if( .not. finished_job )then
                !$omp workshare
                mat4compare = ccimage_unordered%bimat
                !$omp end workshare nowait
                !$omp single
                diff = 0
                !$omp end single nowait
                do i = 1, self%bldim(1)
                    do j = 1, self%bldim(2)
                        do k = 1, self%bldim(3)
                            if( self%bimat(i,j,k) > 0) then ! not background
                                if(self%bldim(3) > 1) then
                                    call neigh_8_3D(ccimage_unordered%bldim, ccimage_unordered%bimat, [i,j,k], neigh_8_pixs, nsz)
                                else
                                    call neigh_8   (ccimage_unordered%bldim, ccimage_unordered%bimat, [i,j,1], neigh_8_pixs, nsz)
                                endif
                                ccimage_unordered%bimat(i,j,k) = minval(neigh_8_pixs(:nsz), neigh_8_pixs(:nsz) > 0)
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
        do i = 1, self%bldim(1)
            do j = 1, self%bldim(2)
                do k = 1, self%bldim(3)
                    if( ccimage_unordered%bimat(i,j,k) > 0) then ! rmat == 0  --> background
                        cnt = cnt + 1
                        tmp = ccimage_unordered%bimat(i,j,k)
                        where(ccimage_unordered%bimat == tmp)
                            ccimage%bimat = cnt
                            ccimage_unordered%bimat  = 0         ! Not to consider this cc again
                        endwhere
                    endif
                enddo
            enddo
        enddo
        ! update instance
        ccimage%bimat_is_set = .true.
        ccimage%nccs = maxval(ccimage%bimat)
        ! kill
        deallocate(mat4compare)
        call ccimage_unordered%kill_bimg
        if( black_present .and. black .eqv. .true.) then
            ! flip foreground to background and vice versa
            self%bimat = -1 * (self%bimat - 1)
        endif
        call self%update_img_rmat
    end subroutine find_ccs

    ! The result is the size(# of pixels) of each cc.
    !  (cc = connected component)
    function size_ccs( self ) result( sz )
        class(binimage), intent(in) :: self
        integer, allocatable :: sz(:)
        integer :: n_cc, imax
        if(.not. any(self%bimat > 0)) THROW_HARD('Inputted non-existent cc; size_ccs')
        if( allocated(sz) ) deallocate( sz )
        imax = maxval(self%bimat)
        allocate(sz(imax), source = 0)
        do n_cc = 1, imax
            sz(n_cc) = count(self%bimat == n_cc)
        enddo
    end function size_ccs

    ! This subroutine takes in input a connected component (cc) image
    ! and sets to 0 the cc which has size (# pixels) smaller than minmaxsz(1)
    ! or bigger than minmaxsz(2).
    ! It is created in order to prepare a micrograph for picking particles.
    subroutine elim_ccs( self, minmaxsz )
        class(binimage), intent(inout) :: self
        integer,         intent(in)    :: minmaxsz(2)
        integer, allocatable :: sz(:)
        integer :: n_cc
        if(.not. any(self%bimat > 0)) THROW_HARD('Inputted non-existent cc; elim_ccs')
        sz = self%size_ccs()
        do n_cc = 1, size(sz) ! for each cc
            if(sz(n_cc) < minmaxsz(1) .or. sz(n_cc) > minmaxsz(2)) then ! if the cc has size < min_sz or > max_sz
                where(self%bimat == n_cc)  ! rmat == label
                    self%bimat = 0         ! set to 0
                endwhere
            endif
        enddo
        ! re-oder cc
        call self%order_ccs()
        ! update number of connected components
        self%nccs = maxval(self%bimat)
        deallocate(sz)
        call self%update_img_rmat
    end subroutine elim_ccs

    ! This function is reorders the connected components (cc:s)
    ! after some of them have been eliminated, so that they have contiguous
    ! labelling (1,2,3..) withouty absenses (NOT 1,2,5..).
    ! Self is should be a connected component image.
    subroutine order_ccs(self)
        class(binimage), intent(inout) :: self
        integer, allocatable :: imat_aux(:,:,:)
        integer :: cnt, i
        if(.not. any(self%bimat > 0)) THROW_HARD('Inputted non-existent cc; order_ccs')
        cnt = 0
        allocate(imat_aux(self%bldim(1),self%bldim(2),self%bldim(3)), source = self%bimat)
        do i = 1, maxval(self%bimat)    ! for each cc
            if(any(imat_aux == i)) then ! there is cc labelled i
                cnt = cnt + 1           ! increasing order cc
                where(imat_aux == i) self%bimat = cnt
            endif
        enddo
        deallocate(imat_aux)
        call self%update_img_rmat
    end subroutine order_ccs

    ! This subroutine takes in input a connected components (cc)
    ! image and eliminates some of the ccs according to their size.
    ! The decision method calculates the avg size of the ccs
    ! and their standar deviation.
    ! Elimin ccs which have size: > ave + 2.*stdev
    !                             < ave - 2.*stdev
    ! If present(minmax_rad(2)), it cleans up ccs more.
    subroutine polish_ccs( self, minmax_rad )
        class(binimage), intent(inout) :: self
        real, optional,  intent(in)    :: minmax_rad(2)
        integer, allocatable :: sz(:)
        real    :: lt, ht     ! low and high thresh for ccs polising
        real    :: ave, stdev ! avg and stdev of the size od the ccs
        integer :: n_cc
        if(.not. any(self%bimat > 0)) THROW_HARD('Inputted non-existent cc; polish_ccs')
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
        ! call img_cc%order_cc() is already done in the elim_cc subroutine
        ! Use particle radius. Hypothesis is the biggest possible area is when
        ! particle is circular (with rad 3*minmax_rad(2)), the smallest one is when the particle is rectangular
        ! with lengths minmax_rad(1)/2 and minmax_rad(2)/2
        if( present(minmax_rad) )then
            call self%elim_ccs([int(minmax_rad(1)*minmax_rad(2)/4.), int(2.*PI*(3.*minmax_rad(2))**2.)])
        endif
        ! update number of connected components
        self%nccs = maxval(self%bimat)
        call self%update_img_rmat
    end subroutine polish_ccs

    ! This subroutine calculates the diamenter of the
    ! connected component labelled n_cc in the connected
    ! component image img_cc
    subroutine diameter_cc( self, n_cc, diam )
        class(binimage), intent(inout) :: self
        integer,         intent(in)    :: n_cc
        real,            intent(out)   :: diam
        integer, allocatable :: pos(:,:)       ! position of the pixels of a fixed cc
        integer, allocatable :: imat_cc(:,:,:) ! auxiliary
        logical, allocatable :: msk(:)         ! For using function pixels_dist
        real  :: center_of_mass(3)             ! geometrical center of mass
        real  :: radius
        if(.not. any(self%bimat > 0)) THROW_HARD('Inputted non-existent cc; diameter_cc')
        allocate(imat_cc(self%bldim(1),self%bldim(2),self%bldim(3)), source=self%bimat)
        where(imat_cc .ne. n_cc) imat_cc = 0   ! keep just the considered cc
        ! Find center of mass of the cc
        call get_pixel_pos(imat_cc,pos)
        center_of_mass(1) = sum(pos(1,:))/real(size(pos,dim = 2))
        center_of_mass(2) = sum(pos(2,:))/real(size(pos,dim = 2))
        center_of_mass(3) = 1.
        allocate(msk(size(pos, dim =2)), source=.true.)
        ! Calculate maximim radius
        radius = pixels_dist(center_of_mass,real(pos),'max',msk)
        ! Return diameter
        diam = 2.*radius
        deallocate(msk, pos, imat_cc)
    end subroutine diameter_cc

    subroutine diameter_bin( self, diam )
        class(binimage), intent(inout) :: self
        real,            intent(out)   :: diam
        integer :: i, ii, j, jj, k, kk, maxdistsq, distsq
        if( .not. self%bimat_is_set ) call self%set_imat
        maxdistsq = 0
        !$omp parallel do collapse(3) default(shared) private(i,j,k,ii,jj,kk,distsq) schedule(static)&
        !$omp proc_bind(close) reduction(max:maxdistsq)
        do i=1,self%bldim(1)
            do j=1,self%bldim(2)
                do k=1,self%bldim(3)
                    if( self%bimat(i,j,k) < 1 ) cycle
                    do ii=1,self%bldim(1)
                        do jj=1,self%bldim(2)
                            do kk=1,self%bldim(3)
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
    end subroutine diameter_bin

    ! This subroutine modifies the input cc image by making it bin
    ! where just the cc n_cc is kept.
    subroutine cc2bin( self, n_cc )
        class(binimage), intent(inout) :: self  ! cc image
        integer,         intent(in)    :: n_cc  ! label of the cc to keep
        if(.not. any(self%bimat > 0)) THROW_HARD('Inputted non-existent cc; cc2bin')
        where(self%bimat .ne. n_cc) self%bimat = 0
        where(self%bimat > 0)       self%bimat = 1
        self%bimat_is_set = .true.
        call self%update_img_rmat
    end subroutine cc2bin

    ! MORPHOLOGICAL OPERATIONS

    ! This subroutine implements the morphological operation dilatation for 2D binary images
    subroutine dilatate( self )
        class(binimage), intent(inout) :: self
        type(binimage) :: self_copy
        integer        :: neigh_8_inds(3,8), i, j, k, nsz
        if( .not. self%bimat_is_set ) call self%set_imat
        call self_copy%copy_bimg(self)
        if(self%bldim(3) /= 1) THROW_HARD('Non implemented for volumes!; dilatate')
        do i = 1, self%bldim(1)
            do j = 1, self%bldim(2)
                if(self_copy%bimat(i,j,1) == 1) then  ! just for white pixels
                    call neigh_8(self%bldim, [i,j,1], neigh_8_inds, nsz)
                    do k = 1, nsz
                        self%bimat(neigh_8_inds(1,k),neigh_8_inds(2,k),1) = 1
                    enddo
                endif
            enddo
        enddo
        call self_copy%kill_bimg
        call self%update_img_rmat
    end subroutine dilatate

    ! This subroutine implements the morphological operation erosion for 2D binary images
    subroutine erode(self, label)
        class(binimage),   intent(inout) :: self
        integer, optional, intent(in)    :: label
        logical, allocatable :: border(:,:,:)
        if( .not. self%bimat_is_set ) call self%set_imat
        if( present(label) ) then
            call self%border_mask(border,label, .true.)
        else
            call self%border_mask(border)
        endif
        where(border) self%bimat = 0
        call self%update_img_rmat
    end subroutine erode

    ! This subroutine builds the logical array 'border' of the same dims of
    ! the input image self. Border is true in the border pixels in the binary image self.
    ! It is necessary for the morphological operation of erosion.
    ! Optional parameter neigh_4 decides whether to perform calculations
    ! using 4neighbours of 8neighbours. It's just for 3D volumes. It is
    ! useful in the case of very small objects, to make sure that the
    ! surface does not coincide with the volume itself.
    subroutine border_mask( self, border, label, four )
        class(binimage),      intent(inout) :: self
        logical, allocatable, intent(inout) :: border(:,:,:)
        integer, optional,    intent(in)    :: label
        logical, optional,    intent(in)    :: four
        integer, allocatable :: neigh_8_pixs(:)
        integer :: i, j, k, nsz, llabel
        logical :: ffour
        if( .not. self%bimat_is_set ) call self%set_imat
        llabel = 1
        if(present(label)) llabel = label
        ffour  = .false.
        if(present(four)) ffour = four
        if(present(four) .and. self%bldim(3)==1) THROW_HARD('4-neighbours identification is not implemented for 2D images; border_mask')
        if( allocated(border) ) deallocate(border)
        allocate(border(self%bldim(1),self%bldim(2),self%bldim(3)), source=.false.)
        if(self%bldim(3) == 1) then
            allocate(neigh_8_pixs(9),  source=0)
            do i = 1,self%bldim(1)
                do j = 1, self%bldim(2)
                    do k = 1, self%bldim(3)
                        if(self%bimat(i,j,k) == llabel ) then
                            call neigh_8(self%bldim, self%bimat, [i,j,k], neigh_8_pixs, nsz)
                            if(any(neigh_8_pixs(:nsz)<1)) border(i,j,k) = .true.
                        endif
                    enddo
                enddo
            enddo
        elseif(.not. ffour) then
            allocate(neigh_8_pixs(27), source=0)
            do i = 1,self%bldim(1)
                do j = 1, self%bldim(2)
                    do k = 1, self%bldim(3)
                        if(self%bimat(i,j,k) == llabel ) then
                            call neigh_8_3D(self%bldim, self%bimat, [i,j,k], neigh_8_pixs, nsz)
                            if(any(neigh_8_pixs(:nsz)<1)) border(i,j,k) = .true.
                        endif
                    enddo
                enddo
            enddo
        elseif(ffour) then
            allocate(neigh_8_pixs(6),  source=0) ! actually this is neigh_4
            do i = 1,self%bldim(1)
                do j = 1, self%bldim(2)
                    do k = 1, self%bldim(3)
                        if(self%bimat(i,j,k) == llabel ) then
                            call neigh_4_3D(self%bldim, self%bimat, [i,j,k], neigh_8_pixs, nsz)
                            if(any(neigh_8_pixs(:nsz)<1)) border(i,j,k) = .true.
                        endif
                    enddo
                enddo
            enddo
        endif
        deallocate(neigh_8_pixs)
    end subroutine border_mask

    ! This subroutine fills the holes in a binary image
    ! Idea: A hole is a set of background pixels that cannot be
    ! reached by filling in the background from the edge of the image.
    ! Find the cc to which the background belongs. Whatever is not background,
    ! set it to foreground.
    subroutine fill_holes( self )
        class(binimage), intent(inout) :: self
        type (binimage)      :: img_cc ! connected component image
        real,    pointer     :: rmat(:,:,:)
        integer, allocatable :: imat_cc(:,:,:)
        integer :: seed
        if(self%bldim(3) > 1) THROW_HARD('Not implemented for volumes! fill_holes')
        call self%find_ccs(img_cc, black=.true.) ! detect the connnected components in the background as well
        imat_cc = nint(img_cc%get_rmat())
        seed    = img_cc%bimat(1,1,1) ! the pxl (1,1,1) should belong to the background
        ! Set to foreground whatever is not background
        where(img_cc%bimat == seed)
            self%bimat = 0
        elsewhere
            self%bimat = 1
        endwhere
        if( allocated(imat_cc) ) deallocate(imat_cc)
        call img_cc%kill_bimg
        call self%update_img_rmat
    end subroutine fill_holes

    !> \brief grow_bins adds nlayers of pixels bordering the background in a binary image
    subroutine grow_bins( self, nlayers )
        class(binimage), intent(inout) :: self
        integer,         intent(in)    :: nlayers
        logical, allocatable :: add_pixels(:,:,:), template(:,:,:)
        integer :: i,j,k, tsz(3,2), win(3,2), pdsz(3,2)
        if( .not. self%bimat_is_set ) call self%set_imat
        tsz(:,1) = -nlayers
        tsz(:,2) = nlayers
        if( self%bldim(3) == 1 ) tsz(3,:) = 1
        allocate(template(tsz(1,1):tsz(1,2), tsz(2,1):tsz(2,2), tsz(3,1):tsz(3,2)))
        pdsz(:,1) = 1 - nlayers
        pdsz(:,2) = self%bldim + nlayers
        if( self%bldim(3) == 1 ) pdsz(3,:) = 1
        allocate(add_pixels(pdsz(1,1):pdsz(1,2), pdsz(2,1):pdsz(2,2), pdsz(3,1):pdsz(3,2)))
        ! template matrix
        template = .true.
        do i = tsz(1,1), tsz(1,2)
            do j = tsz(2,1), tsz(2,2)
                if( self%bldim(3) == 1 )then
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
        forall( i=1:self%bldim(1), j=1:self%bldim(2), k=1:self%bldim(3), self%bimat(i,j,k)==1 ) add_pixels(i,j,k) = .true.
            ! cycle
            if( self%bldim(3) > 1 )then
                do i = 1, self%bldim(1)
                    if( .not.any(self%bimat(i,:,:) > 0) ) cycle
                    do j = 1, self%bldim(2)
                        if( .not.any(self%bimat(i,j,:) > 0) ) cycle
                        win(1:2,1) = [i, j] - nlayers
                        win(1:2,2) = [i, j] + nlayers
                        do k = 1, self%bldim(3)
                            if (self%bimat(i,j,k) < 1) cycle
                            win(3,1) = k - nlayers
                            win(3,2) = k + nlayers
                            add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
                            &add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)).or.template
                        enddo
                    enddo
                enddo
            else
                do i=1,self%bldim(1)
                    if( .not.any(self%bimat(i,:,1) > 0) ) cycle
                    do j=1,self%bldim(2)
                        win(1:2,1) = [i, j] - nlayers
                        win(1:2,2) = [i, j] + nlayers
                        if (self%bimat(i,j,1) < 1) cycle
                        add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), 1) =&
                        &add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), 1).or.template(:,:,1)
                    enddo
                enddo
            endif
            ! finalize
            self%bimat = 0
            forall( i=1:self%bldim(1), j=1:self%bldim(2), k=1:self%bldim(3), add_pixels(i,j,k) ) self%bimat(i,j,k) = 1
                deallocate(template, add_pixels)
            call self%update_img_rmat
        end subroutine grow_bins

        !>  \brief cos_edge applies cosine squared edge to a binary image
        subroutine cos_edge( self, cos_img, falloff )
            class(binimage), intent(inout)  :: self    ! input
            class(image),     intent(inout) :: cos_img ! output
            integer,         intent(in)     :: falloff
            integer, allocatable :: imat(:,:,:)
            real,    pointer     :: rmat(:,:,:)
            real                 :: rfalloff, scalefactor
            integer              :: i, j, k, is, js, ks, ie, je, ke
            integer              :: il, ir, jl, jr, kl, kr, falloff_sq
            if( falloff<=0 ) THROW_HARD('only stictly positive values for edge fall-off allowed; cos_edge')
            if( .not. self%bimat_is_set ) call self%set_imat
            rfalloff    = real( falloff )
            falloff_sq  = falloff**2
            scalefactor = PI / rfalloff
            allocate(imat(self%bldim(1),self%bldim(2),self%bldim(3)), source=self%bimat)
            call cos_img%new(self%bldim,self%bsmpd)
            call cos_img%get_rmat_ptr(rmat)
            rmat(1:self%bldim(1),1:self%bldim(2),1:self%bldim(3)) = real(imat)
            rmat = rmat/maxval(rmat(:self%bldim(1),:self%bldim(2),:self%bldim(3)))
            if( self%bldim(3) == 1 )then
                ! 2d
                do i=1,self%bldim(1)
                    is = max(1,i-1)                   ! left neighbour
                    ie = min(i+1,self%bldim(1))       ! right neighbour
                    il = max(1,i-falloff)             ! left bounding box limit
                    ir = min(i+falloff,self%bldim(1)) ! right bounding box limit
                    if( .not. any(imat == 1) ) cycle
                    do j=1,self%bldim(2)
                        if( imat(i,j,1) /= 1 ) cycle
                        js = max(1,j-1)
                        je = min(j+1,self%bldim(2))
                        if( any( imat(is:ie,js:je,1) < 1) )then
                            jl = max(1,j-falloff)
                            jr = min(j+falloff,self%bldim(2))
                            call update_mask_2d
                        endif
                    end do
                end do
            else
                ! 3d
                do i=1,self%bldim(1)
                    is = max(1,i-1)                   ! left neighbour
                    ie = min(i+1,self%bldim(1))       ! right neighbour
                    il = max(1,i-falloff)             ! left bounding box limit
                    ir = min(i+falloff,self%bldim(1)) ! right bounding box limit
                    if( .not. any(imat == 1) ) cycle  ! no values equal to one
                    do j=1,self%bldim(2)
                        js = max(1,j-1)
                        je = min(j+1,self%bldim(2))
                        jl = max(1,j-falloff)
                        jr = min(j+falloff,self%bldim(2))
                        if(.not. any(imat == 1) ) cycle ! cycle if equal to one
                        do k=1,self%bldim(3)
                            if(imat(i,j,k) /= 1) cycle
                            ! within mask region
                            ks = max(1,k-1)
                            ke = min(k+1,self%bldim(3))
                            if( any( imat(is:ie,js:je,ks:ke) < 1) )then
                                ! update since has a masked neighbour
                                kl = max(1,k-falloff)
                                kr = min(k+falloff,self%bldim(3))
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
                        if( rmat(ii,jj,1)<1 ) rmat(ii,jj,1) = max(local_versine(real(dist_sq)), rmat(ii,jj,1))
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
                            if(dist_sq > falloff_sq) cycle
                            if( rmat(ii,jj,kk)<1 ) rmat(ii,jj,kk) = max(local_versine(real(dist_sq)), rmat(ii,jj,kk))
                        enddo
                    enddo
                enddo
            end subroutine update_mask_3d

            !> Local elemental cosine edge function
            elemental real function local_versine( r_sq )result( c )
                real, intent(in) :: r_sq
                c = 0.5 * (1. - cos(scalefactor*(sqrt(r_sq)-rfalloff)) )
            end function local_versine

    end subroutine cos_edge

    ! DESTRUCT

    subroutine kill_bimg( self )
        class(binimage), intent(inout) :: self
        if(allocated(self%bimat)) deallocate(self%bimat)
        call self%kill
    end subroutine kill_bimg

end module simple_binimage
