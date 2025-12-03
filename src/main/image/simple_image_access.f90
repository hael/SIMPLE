submodule (simple_image) simple_image_access
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none
contains

    !===========================
    ! Basic shape / metadata
    !===========================

    module pure function get_array_shape( self ) result( shape)
        class(image), intent(in) :: self
        integer :: shape(3)
        shape = self%array_shape
    end function get_array_shape

    module pure function get_ldim( self ) result( ldim )
        class(image), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function get_ldim

    module pure integer function get_box( self )
        class(image), intent(in) :: self
        get_box = self%ldim(1)
    end function get_box

    module pure function get_smpd( self ) result( smpd )
        class(image), intent(in) :: self
        real :: smpd
        smpd = self%smpd
    end function get_smpd

    module pure function get_nyq( self ) result( nyq )
        class(image), intent(in) :: self
        integer :: nyq
        nyq = fdim(self%ldim(1)) - 1
    end function get_nyq

    module pure function get_filtsz( self ) result( n )
        class(image), intent(in) :: self
        integer :: n
        n = fdim(self%ldim(1)) - 1
    end function get_filtsz

    module pure function get_shconst( self ) result( shconst )
        class(image), intent(in) :: self
        real :: shconst(3)
        shconst = self%shconst
    end function get_shconst

    module subroutine set_smpd( self, smpd )
        class(image), intent(inout) :: self
        real,         intent(in)    :: smpd
        self%smpd = smpd
        call self%fit%new(self%ldim, self%smpd)
    end subroutine set_smpd

    !===========================
    ! Getters
    !===========================

    module function get( self, logi ) result( val )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real :: val
        if( logi(1) > self%ldim(1) .or. logi(1) < 1 )then
            val = 0.
            return
        endif
        if( logi(2) > self%ldim(2) .or. logi(2) < 1 )then
            val = 0.
            return
        endif
        if( logi(3) > self%ldim(3) .or. logi(3) < 1 )then
            val = 0.
            return
        endif
        val = self%rmat(logi(1),logi(2),logi(3))
    end function get

    module pure function get_rmat( self ) result( rmat )
        class(image), intent(in) :: self
        real, allocatable :: rmat(:,:,:)
        integer :: ldim(3)
        ldim = self%ldim
        allocate(rmat(ldim(1),ldim(2),ldim(3)), source=self%rmat(:ldim(1),:ldim(2),:ldim(3)))
    end function get_rmat

    module subroutine get_mat_ptrs( self, mat_ptrs )
        class(image),      target, intent(in)  :: self
        class(image_ptr),          intent(out) :: mat_ptrs
        mat_ptrs%cmat => self%cmat
        mat_ptrs%rmat => self%rmat
    end subroutine get_mat_ptrs

    module subroutine get_rmat_ptr( self, rmat_ptr )
        class(image), target,        intent(in)  :: self
        real(kind=c_float), pointer, intent(out) :: rmat_ptr(:,:,:)
        rmat_ptr => self%rmat
    end subroutine get_rmat_ptr

    module pure subroutine get_rmat_sub( self, rmat )
        class(image), intent(in)  :: self
        real,         intent(out) :: rmat(self%ldim(1),self%ldim(2),self%ldim(3))
        rmat = self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))
    end subroutine get_rmat_sub

    module pure function get_rmat_at_1( self, logi ) result( val )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: logi(3)
        real :: val
        val = self%rmat(logi(1),logi(2),logi(3))
    end function get_rmat_at_1

    module pure function get_rmat_at_2( self, i,j,k ) result( val )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: i,j,k
        real :: val
        val = self%rmat(i,j,k)
    end function get_rmat_at_2

    module pure function get_cmat( self ) result( cmat )
        class(image), intent(in) :: self
        complex, allocatable :: cmat(:,:,:)
        allocate(cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3)), source=self%cmat)
    end function get_cmat

    module subroutine get_cmat_ptr( self, cmat_ptr )
        class(image), target,                   intent(in)  :: self
        complex(kind=c_float_complex), pointer, intent(out) :: cmat_ptr(:,:,:)
        cmat_ptr => self%cmat
    end subroutine get_cmat_ptr

    module pure subroutine get_cmat_sub( self, cmat )
        class(image), intent(in)  :: self
        complex,      intent(out) :: cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        cmat=self%cmat
    end subroutine get_cmat_sub

    module pure function get_cmat_at_1( self, phys ) result( comp )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: phys(3)
        complex :: comp
        comp = self%cmat(phys(1),phys(2),phys(3))
    end function get_cmat_at_1

    module pure function get_cmat_at_2( self, h,k,l ) result( comp )
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k,l
        complex :: comp
        comp = self%cmat(h,k,l)
    end function get_cmat_at_2

    module pure function get_fcomp( self, logi, phys ) result( comp )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: logi(3), phys(3)
        complex :: comp
        comp = merge(conjg(self%cmat(phys(1),phys(2),phys(3))),&
            &self%cmat(phys(1),phys(2),phys(3)), logi(1)<0)
    end function get_fcomp

    module elemental complex function get_fcomp2D(self, h, k)
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k
        integer :: phys1, phys2
        if (h .ge. 0) then
            phys1 = h + 1
            phys2 = k + 1 + merge(self%ldim(2),0, k<0)
            get_fcomp2D = self%cmat(phys1,phys2,1)
        else
            phys1 = -h + 1
            phys2 = -k + 1 + merge(self%ldim(2),0, -k<0)
            get_fcomp2D = conjg(self%cmat(phys1,phys2,1))
        endif
    end function get_fcomp2D

    !===========================
    ! Setters
    !===========================

    module subroutine set_1( self, logi, val )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: val
        if( logi(1) <= self%ldim(1) .and. logi(1) >= 1 .and. logi(2) <= self%ldim(2)&
            .and. logi(2) >= 1 .and. logi(3) <= self%ldim(3) .and. logi(3) >= 1 )then
            self%rmat(logi(1),logi(2),logi(3)) = val
        endif
    end subroutine set_1

    module subroutine set_2( self, self2set )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2set
        self%ft       = self2set%ft
        self%rmat     = self2set%rmat
        self%wthreads = self2set%wthreads
    end subroutine set_2

    module subroutine set_rmat( self, rmat, ft )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rmat(:,:,:)
        logical,      intent(in)    :: ft
        self%rmat = 0.
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))
        self%ft   = ft
    end subroutine set_rmat

    module pure subroutine set_rmat_at( self, i,j,k, val )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i,j,k
        real,         intent(in)    :: val
        self%rmat(i,j,k) = val
    end subroutine set_rmat_at

    module subroutine set_cmat( self, cmat )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cmat(:,:,:)
        integer :: cdim(3)
        cdim(1) = size(cmat,1)
        cdim(2) = size(cmat,2)
        cdim(3) = size(cmat,3)
        if( all(self%array_shape .eq. cdim) )then
            self%ft   = .true.
            self%cmat = cmplx(0.,0.)
            self%cmat(:cdim(1),:cdim(2),:cdim(3)) = cmat
        else
            write(logfhandle,*) 'dim(cmat): ', cdim
            write(logfhandle,*) 'dim(self%cmat): ', self%array_shape
            THROW_HARD('nonconforming dims; set_cmat')
        endif
    end subroutine set_cmat

    module pure subroutine set_cmat_1( self, cmat )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        self%ft   = .true.
        self%cmat = cmat
    end subroutine set_cmat_1

    module pure subroutine set_cmat_2( self, cval )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cval
        self%ft   = .true.
        self%cmat = cval
    end subroutine set_cmat_2

    module pure subroutine set_cmat_3( self, self2copy )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2copy
        self%ft   = .true.
        self%cmat = self2copy%cmat
    end subroutine set_cmat_3

    module pure subroutine set_cmat_4( self, cmat )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cmat(self%array_shape(1),self%array_shape(2))
        self%ft          = .true.
        self%cmat(:,:,1) = cmat
    end subroutine set_cmat_4

    module pure subroutine set_cmat_at_1( self, phys ,comp)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        complex,      intent(in)    :: comp
        self%cmat(phys(1),phys(2),phys(3)) = comp
    end subroutine set_cmat_at_1

    module pure subroutine set_cmat_at_2( self, h, k, l, comp)
        class(image), intent(inout) :: self
        integer, intent(in) :: h,k,l
        complex, intent(in) :: comp
        self%cmat(h,k,l) =  comp
    end subroutine set_cmat_at_2

    module subroutine set_fcomp( self, logi, phys, comp )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3), phys(3)
        complex,      intent(in)    :: comp
        complex :: comp_here
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        self%cmat(phys(1),phys(2),phys(3)) = comp_here
    end subroutine set_fcomp

    !>  \brief  set pixels to value within a sphere
    module subroutine set_within( self, xyz, radius, val )
        class(image), intent(inout) :: self
        real,         intent(in)    :: xyz(3), radius, val
        real    :: rpos(3), vec(3), dist_sq, radius_sq
        integer :: i,j,k, win(3,2)
        radius_sq = radius**2.
        rpos      = xyz / self%smpd
        win(:,1)  = 1 + floor(rpos - radius)
        where( win(:,1) < 1 )win(:,1) = 1
        win(:,2)  = 1 + ceiling(rpos + radius)
        if(win(1,2) > self%ldim(1)) win(1,2) = self%ldim(1)
        if(win(2,2) > self%ldim(2)) win(2,2) = self%ldim(2)
        if(win(3,2) > self%ldim(3)) win(3,2) = self%ldim(3)
        do i = win(1,1),win(1,2)
            do j = win(2,1),win(2,2)
                do k = win(3,1),win(3,2)
                    vec     = real([i,j,k] - 1) * self%smpd - xyz
                    dist_sq = dot_product(vec,vec)
                    if(dist_sq <= radius_sq)self%rmat(i,j,k) = val
                end do
            end do
      end do
    end subroutine set_within

    module subroutine set_cmats_from_cmats( self1 , self2 , self3, self4, self2set1, self2set2, lims, expcmat3, expcmat4)
        class(image), intent(in)    :: self1, self2,self3,self4
        class(image), intent(inout) :: self2set1, self2set2
        integer,      intent(in)    :: lims(3,2)
        real,         intent(inout) :: expcmat3(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
        real,         intent(inout) :: expcmat4(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
        integer :: h, k, logi(3), phys(3)
        !$omp parallel default(shared) private(h,k,logi,phys) proc_bind(close)
        !$omp workshare
        self1%cmat = self2set1%cmat
        self2%cmat = self2set2%cmat
        !$omp end workshare
        !$omp do collapse(2) schedule(static)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                logi = [h,k,0]
                phys = self1%comp_addr_phys(logi)
                self3%cmat(phys(1),phys(2),phys(3)) = cmplx(expcmat3(h,k),0.)
                self4%cmat(phys(1),phys(2),phys(3)) = cmplx(expcmat4(h,k),0.)
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine set_cmats_from_cmats

    !===========================
    ! Slices, sub-images, freq info
    !===========================

    module subroutine get_slice( self3D, slice, self2D )
        class(image), intent(in)    :: self3D
        integer,      intent(in)    :: slice
        class(image), intent(inout) :: self2D
        self2D%rmat(:,:,1) = self3D%rmat(:,:,slice)
    end subroutine get_slice

    module subroutine set_slice( self3D, slice, self2D )
        class(image), intent(in)    :: self2D
        integer,      intent(in)    :: slice
        class(image), intent(inout) :: self3D
        self3D%rmat(:,:,slice) = self2D%rmat(:,:,1)
    end subroutine set_slice

    module subroutine get_subimg(self, binning, xoffset, yoffset, img)
        class(image), intent(in)  :: self
        integer,      intent(in)  :: binning, xoffset, yoffset
        class(image), intent(out) :: img
        integer :: i,j,ii,jj,ldim(3)
        if( .not.is_even(binning) .and. binning > 0 )then
            THROW_HARD('Binning must be even; get_sub_img')
        endif
        ldim(1:2) = self%ldim(1:2) / binning
        ldim(3)   = 1
        call img%new(ldim, real(binning)*self%smpd)
        do j=1,ldim(2)
            jj = binning*(j-1) + yoffset + 1
            do i=1,ldim(1)
                ii = binning*(i-1) + xoffset + 1
                img%rmat(i,j,1) = self%rmat(ii,jj,1)
            enddo
        enddo
    end subroutine get_subimg

    module pure function get_lfny( self, which ) result( fnyl )
        class(image), intent(in) :: self
        integer,      intent(in) :: which
        integer :: fnyl
        fnyl = self%fit%get_lfny(which)
    end function get_lfny

    module pure function get_lhp( self, which ) result( hpl )
        class(image), intent(in) :: self
        integer,      intent(in) :: which
        integer :: hpl
        hpl = self%fit%get_lhp(which)
    end function get_lhp

    module pure function get_lp( self, ind ) result( lp )
        class(image), intent(in) :: self
        integer,      intent(in) :: ind
        real                     :: lp
        lp = self%fit%get_lp(1, ind)
    end function get_lp

    module pure function get_spat_freq( self, ind ) result( spat_freq )
        class(image), intent(in) :: self
        integer,      intent(in) :: ind
        real                     :: spat_freq
        spat_freq = self%fit%get_spat_freq(1, ind)
    end function get_spat_freq

    module pure function get_find( self, res ) result( ind )
        class(image), intent(in) :: self
        real,         intent(in) :: res
        integer :: ind
        ind = self%fit%get_find(1, res)
    end function get_find

    !===========================
    ! Flags
    !===========================

    module function rmat_associated( self ) result( assoc )
        class(image), intent(in) :: self
        logical :: assoc
        assoc = associated(self%rmat)
    end function rmat_associated

    module function cmat_associated( self ) result( assoc )
        class(image), intent(in) :: self
        logical :: assoc
        assoc = associated(self%cmat)
    end function cmat_associated

    module function is_wthreads( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
        is = self%wthreads
    end function is_wthreads

    module subroutine set_ft( self, is )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: is
        self%ft = is
    end subroutine set_ft

    !===========================
    ! Serialization / misc
    !===========================

    module function serialize_1( self ) result( vec )
        class(image), intent(in) :: self
        real,    allocatable :: vec(:)
        complex, allocatable :: cvec(:)
        if( self%is_ft() )then
            cvec = pack(self%cmat(:self%array_shape(1),:self%array_shape(2),:self%array_shape(3)), mask=.true.)
            vec  = sqrt(real(cvec * conjg(cvec)))
        else
            vec = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=.true.)
        endif
    end function serialize_1

    module function serialize_2( self, l_msk )result( pcavec )
        class(image), intent(in) :: self
        logical,      intent(in) :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3))
        real, allocatable :: pcavec(:)
        integer :: sz, cnt, i, j, k
        sz = count(l_msk)
        allocate(pcavec(sz))
        cnt = 0
        do k=1,self%ldim(3)
            do j=1,self%ldim(2)
                do i=1,self%ldim(1)
                    if( l_msk(i,j,k) )then
                        cnt         = cnt + 1
                        pcavec(cnt) = self%rmat(i,j,k)
                    endif
                end do
            end do
        end do
    end function serialize_2

    module function serialize_3( self, thres ) result( vec )
        class(image), intent(in) :: self
        real,         intent(in) :: thres
        real, allocatable :: vec(:)
        vec = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > thres)
    end function serialize_3

    module subroutine unserialize( self, pcavec, l_msk )
        class(image),      intent(inout) :: self
        real,              intent(in)    :: pcavec(:)
        logical, optional, intent(in)    :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3))
        integer :: sz, sz_msk, i, j, k, cnt
        if( present(l_msk) )then
            sz     = size(pcavec)
            sz_msk = count(l_msk)
            if( sz /= sz_msk )then
                write(logfhandle,*) 'ERROR! Nonconforming sizes'
                write(logfhandle,*) 'sizeof(pcavec): ', sz
                write(logfhandle,*) 'sizeof(l_msk) : ', sz_msk
                THROW_HARD('unserialize')
            endif
        endif
        if( self%ft ) self%ft = .false.
        self%rmat = 0.
        cnt = 0
        if( present(l_msk) )then
            do k=1,self%ldim(3)
                do j=1,self%ldim(2)
                    do i=1,self%ldim(1)
                        if( l_msk(i,j,k) )then
                            cnt = cnt + 1
                            self%rmat(i,j,k) =  pcavec(cnt)
                        endif
                    end do
                end do
            end do
        else
            do k=1,self%ldim(3)
                do j=1,self%ldim(2)
                    do i=1,self%ldim(1)
                        cnt = cnt + 1
                        self%rmat(i,j,k) =  pcavec(cnt)
                    end do
                end do
            end do
        endif
    end subroutine unserialize

    !>  \brief winserialize is for packing/unpacking a serialized image vector for convolutional pca analysis
    module subroutine winserialize( self, coord, winsz, pcavec )
        class(image),      intent(inout) :: self
        real, allocatable, intent(inout) :: pcavec(:)
        integer,           intent(in)    :: coord(:), winsz
        integer :: i, j, k, cnt, npix
        logical :: pack
        if( self%ft ) THROW_HARD('winserialization not yet implemented for Fourier')
        if( self%is_2d() )then
            npix = winsz**2
            call set_action
            cnt = 0
            do i=coord(1),coord(1)+winsz-1
                do j=coord(2),coord(2)+winsz-1
                    cnt = cnt+1
                    if( pack )then
                        if( i > self%ldim(1) .or. j > self%ldim(2) )then
                            pcavec(cnt) = 0.
                        else
                            pcavec(cnt) = self%rmat(i,j,1)
                        endif
                    else
                        if( i > self%ldim(1) .or. j > self%ldim(2) )then
                        else
                            self%rmat(i,j,1) = self%rmat(i,j,1)+pcavec(cnt)
                        endif
                    endif
                end do
            end do
        else
            if( size(coord) < 3 ) THROW_HARD('need a 3D coordinate for a 3D image; winserialize')
            npix = winsz**3
            call set_action
            cnt = 0
            do i=coord(1),coord(1)+winsz-1
                do j=coord(2),coord(2)+winsz-1
                    do k=coord(3),coord(3)+winsz-1
                        cnt = cnt+1
                        if( pack )then
                            if( i > self%ldim(1) .or. j > self%ldim(2) .or. k > self%ldim(3) )then
                                pcavec(cnt) = 0.
                            else
                                pcavec(cnt) = self%rmat(i,j,k)
                            endif
                        else
                            if( i > self%ldim(1) .or. j > self%ldim(2) .or. k > self%ldim(3) )then
                            else
                                self%rmat(i,j,k) = self%rmat(i,j,k)+pcavec(cnt)
                            endif
                        endif
                    end do
                end do
            end do
        endif

        contains

            subroutine set_action
                if( allocated(pcavec) )then
                    if( size(pcavec) /= npix ) THROW_HARD('size mismatch mask/npix; winserialize')
                    pack = .false.
                else
                    pack = .true.
                    allocate( pcavec(npix) )
                    pcavec = 0.
                endif
            end subroutine set_action

    end subroutine winserialize

end submodule simple_image_access
