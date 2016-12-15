!>  \brief  SIMPLE Spider volume/image stack class
module simple_spider
use simple_defs    ! singleton
use simple_jiffys, only: progress, fopen_err, alloc_err
implicit none

public :: spider, test_spider
private

type spider
    private
    real, allocatable     :: spistkhed(:,:)               !< stack or volume header
    real, allocatable     :: spiimghed(:,:)               !< image header 
    real, allocatable     :: pixrow(:)                    !< row of pixels
    integer, allocatable  :: indextab(:,:)                !< table of record indices in a volume file
    character(len=STDLEN) :: endconv=''                   !< endian conversion string
    integer               :: filnum=0                     !< file handle
    integer               :: recsz=0                      !< record size
    integer               :: lenbyt=0, labrec=0, labbyt=0 !< spider constants
    integer               :: ldim(3)=[1,1,1]              !< logical image dimensions
    logical               :: existence=.false.            !< indicates existence
  contains
    ! CONSTRUCTORS
    procedure          :: new
    ! I/O
    procedure          :: open
    procedure          :: close
    procedure          :: is_open
    procedure, private :: read_1
    procedure, private :: read_2
    generic :: read => read_1, read_2
    procedure, private :: write_1
    procedure, private :: write_2
    generic :: write => write_1, write_2
    ! PRINTERS/GETTERS
    procedure          :: print_img_hed
    procedure          :: print_hed
    procedure          :: get_nptcls
    procedure          :: get_ldim
    ! CALCULATORS
    procedure          :: copy
    procedure          :: make_pattern_stack
    procedure          :: pad
    procedure          :: clip
    procedure          :: norm
    procedure          :: stats
    procedure          :: neg
    procedure          :: rankify
    procedure          :: masscen
    procedure          :: cure
    procedure          :: resize
    procedure          :: acf
    procedure          :: add_noise
    procedure          :: frameavg
    procedure          :: fwd_ft
    procedure          :: ft2img
    procedure          :: bwd_ft
    procedure          :: shift
    procedure          :: bp
    procedure          :: phase_rand
    procedure          :: apply_ctf
    procedure          :: ctf2img
    procedure          :: shrot
    procedure          :: mask
    procedure          :: bin
    procedure          :: make_avg
    procedure          :: make_cavgs
    procedure          :: make_cavg
    procedure, private :: print_internal_img_hed
    procedure, private :: make_headers
    procedure, private :: imgnr2recpos
    procedure :: find_endconv
    ! DESTRUCTOR
    procedure :: kill
end type

interface spider
    module procedure constructor 
end interface

contains

    ! CONSTRUCTORS
    
    !>  \brief is is for constructing a new spider object
    function constructor( ldim ) result( self )
        integer, intent(in) :: ldim(3)
        type(spider)        :: self   
        call self%new( ldim ) 
    end function

    !>  \brief is for constructing a new spider object
    subroutine new( self, ldim )
        class(spider), intent(inout) :: self  
        integer, intent(in)          :: ldim(3)
        integer                      :: alloc_stat, cnt, j
        call self%kill
        ! set const:
        self%ldim = ldim
        ! from spider:
        self%lenbyt = self%ldim(1)*4
        self%labrec = 1024/self%lenbyt
        if( mod(1024,self%lenbyt) /= 0 ) self%labrec = self%labrec+1
        self%labbyt = self%labrec*self%lenbyt
        ! allocate:
        allocate( self%pixrow(self%ldim(1)), self%indextab(self%ldim(3),2),&
        self%spistkhed(self%labrec,self%ldim(1)), self%spiimghed(self%labrec,self%ldim(1)),&
        stat=alloc_stat )
        call alloc_err('new; simple_spider', alloc_stat)
        self%pixrow     = 0.
        self%spistkhed    = 0.
        self%spiimghed = 0.
        inquire( iolength=self%recsz ) self%pixrow
        ! make index array for extracting the slices of the volume
        cnt = self%labrec
        do j=1,self%ldim(3)
            self%indextab(j,1) = cnt+1 ! from
            cnt = cnt+self%ldim(2) 
            self%indextab(j,2) = cnt   ! to
        end do
        self%existence = .true.
    end subroutine

    ! I/O
    
    !>  \brief  is for opening a spider file
    subroutine open( self, name, replace )
        use simple_jiffys, only: get_fileunit, del_binfile
        class(spider), intent(inout)  :: self 
        character(len=*), intent(in)  :: name
        character(len=*), intent(in), optional :: replace
        real    :: hed(self%labrec,self%ldim(1))
        logical :: here, test(6), foundit
        integer :: j, ier
        if( self%is_open() ) call self%close
        inquire(FILE=name, EXIST=here)
        if( present(replace) .and. here )then
            call del_binfile(name)
            here = .false.
        endif
        if( here )then
            foundit = self%find_endconv( name )
            self%filnum = get_fileunit( )
            open(unit=self%filnum, convert=self%endconv, status='old', action='readwrite',& !***CONVERT
            file=name, access='direct', form='unformatted', recl=self%recsz, iostat=ier)
            call fopen_err('open; simple_spider', ier)
            do j=1,self%labrec
                read(self%filnum, rec=j) hed(j,:)
            end do
            test(1) = self%lenbyt  == nint(hed(1,23))
            test(2) = self%labrec  == nint(hed(1,13))
            test(3) = self%labbyt  == nint(hed(1,22))
            test(4) = self%ldim(1) == nint(hed(1,12)) 
            test(5) = self%ldim(2) == nint(hed(1,2))
            test(6) = self%ldim(3) == nint(hed(1,1))
            if( .not. all(test) ) stop 'stack object does not conform with file; open; simple_spider'
            self%spistkhed = hed
        else        
            ! associate a file with this object:
            self%filnum = get_fileunit( )
            open(unit=self%filnum, convert='LITTLE_ENDIAN', status='new', action='readwrite',& !***CONVERT
            file=name, access='direct', form='unformatted', recl=self%recsz, iostat=ier)
            ! write stack header:
            call self%make_headers( self%ldim, .false. )
            do j=1,self%labrec
                write(self%filnum, rec=j) self%spistkhed(j,:)
            end do
        endif
    end subroutine
        
    !>  \brief  is for closing the spider file opened above
    subroutine close( self )
        class(spider), intent(in) :: self !< object
        if( self%is_open() ) close(unit=self%filnum)
    end subroutine
    
    !>  \brief  is for checking if spider file is open
    function is_open( self ) result( is )
        class(spider), intent(in) :: self
        logical :: is
        inquire( unit=self%filnum, opened=is )
    end function

    !>  \brief  is for reading a single image, with auto normalization
    subroutine read_1( self, img, nr )
        use simple_image, only: image
        class(spider), intent(inout) :: self
        class(image), intent(inout)  :: img  
        integer, intent(in)          :: nr
        integer                      :: hedinds(2), iminds(2), i, j, cnt, nans
        real                         :: hed(self%labrec,self%ldim(1))
        logical                      :: fted
        if( .not. self%is_open() )stop 'need to open stack for read; read_1; simple_spider'
        if( img%is_3d() ) stop 'cannot read 3d images with this method; read_1; simple_spider'
        ! make sure that dimensions are same
        if( .not. img%same_dims(self%ldim) )then
            stop 'incompatible dims; read_1; simple_spider'
        endif
        call self%imgnr2recpos( nr, hedinds, iminds )
        ! read header
        cnt = 0
        do i=hedinds(1),hedinds(2)
            cnt = cnt +1
            read(self%filnum, rec=i) hed(cnt,:)
        end do
        ! find out if image is fted
        fted = .false.
        if( nint(hed(1,5)) < 0 ) fted = .true.
        call img%set_ft(fted)
        ! read data
        cnt = 0
        nans = 0
        do j=iminds(1),iminds(2)
            cnt = cnt+1
            read(self%filnum, rec=j) self%pixrow
            do i=1,self%ldim(1)
                if( isnan(self%pixrow(i)) )then
                    nans = nans+1
                    self%pixrow(i) = 0.
                endif 
                call img%set([i,cnt,1], self%pixrow(i))
            end do
        end do
        if( nans > 0 )then
            write(*,'(a,1x,I5,a)') 'WARNING!', nans, ' NAN:s during read'
        endif
        if( .not. img%is_ft() ) call img%norm
    end subroutine
    
    !>  \brief  is for reading a spider volume, with auto normalization
    subroutine read_2( self, img )
        use simple_image, only: image
        class(spider), intent(inout) :: self
        class(image), intent(inout)  :: img
        integer                      :: i, j, k, cnt, nans
        real                         :: hed(self%labrec,self%ldim(1))
        logical                      :: fted
        if( .not. self%is_open() ) stop 'need to open volume for read; read_2; simple_spider'
        if( img%is_2d() ) stop 'cannot read 2d images with this method; read_2; simple_spider'
        ! make sure that dimensions are same
        if( .not. img%same_dims(self%ldim) )then
            stop 'incompatible dims; read_2; simple_spider'
        endif
        ! read header
        cnt = 0
        do i=1,self%labrec
            cnt = cnt+1
            read(self%filnum, rec=i) hed(cnt,:)
        end do
        ! set ft state
        fted = .false.
        if( self%spistkhed(1,5) < 0 ) fted = .true.
        call img%set_ft(fted)
        ! read volume
        nans = 0
        do k=1,self%ldim(3)
            cnt = 0
            do j=self%indextab(k,1),self%indextab(k,2) 
                cnt = cnt+1
                read(self%filnum, rec=j) self%pixrow
                do i=1,self%ldim(1)
                    if( isnan(self%pixrow(i)) )then
                        nans = nans+1
                        self%pixrow(i) = 0.
                    endif
                    call img%set([i,cnt,k], self%pixrow(i))
                end do
            end do
        end do
        if( nans > 0 )then
            write(*,'(a,1x,I5,a)') 'WARNING!', nans, ' NAN:s during read'
        endif
        if( .not. img%is_ft() ) call img%norm
    end subroutine
    
    !>  \brief  is for writing an image to spider stack
    subroutine write_1( self, img, nr )
        use simple_image, only: image
        class(spider), intent(inout) :: self 
        class(image), intent(inout)  :: img  
        integer, intent(in)          :: nr
        real                         :: hed(self%labrec,self%ldim(1))
        integer                      :: hedinds(2), iminds(2), i, j, cnt
        if( .not. self%is_open() ) stop 'need to open stack for write; write_1; simple_spider'
        if( img%is_3d() ) stop 'cannot write 3d images with this method; write_1; simple_spider'
        ! make sure that dimensions are same
        if( .not. img%same_dims(self%ldim) )then
            stop 'incompatible dims; write_2; simple_spider'
        endif
        call self%make_headers( self%ldim, img%is_ft() )
        call self%imgnr2recpos(nr, hedinds, iminds)
        ! if nr > nptcls in stack, update header
        read(self%filnum, rec=1) hed(1,:)
        if( nr > nint(hed(1,26)) ) hed(1,26) = real(nr)
        self%spiimghed(1,27) = hed(1,26)
        write(self%filnum, rec=1) hed(1,:)
        read(self%filnum, rec=1) hed(1,:)
        ! write image header
        cnt = 0
        do j=hedinds(1),hedinds(2)
            cnt = cnt+1
            write(self%filnum, rec=j) self%spiimghed(cnt,:)
        end do
        ! write image data
        cnt = 0
        do j=iminds(1),iminds(2)
            cnt = cnt+1
            do i=1,self%ldim(1)
                self%pixrow(i) = img%get([i,cnt,1])
            end do
            write(self%filnum, rec=j) self%pixrow
        end do
    end subroutine
    
    !>  \brief  is for writing a spider volume to file 
    subroutine write_2( self, img )
        use simple_image, only: image
        class(spider), intent(inout) :: self
        class(image), intent(inout)  :: img 
        integer :: i, j, k, cnt
        if( .not. self%is_open() ) stop 'need to open volume for write; write_2; simple_spider'
        if( img%is_2d() ) stop 'cannot write 2d images with this method; write_2; simple_spider'
        ! make sure that dimensions are same
        if( .not. img%same_dims(self%ldim) )then
            write(*,*) 'img dim:', img%get_ldim()
            write(*,*) 'stack dim:', self%ldim
            stop 'incompatible dims; write_2; simple_spider'
        endif
        call self%make_headers( self%ldim, img%is_ft() )
        ! write volume header
        do i=1,self%labrec
            write(self%filnum, rec=i) self%spistkhed(i,:)
        end do
        ! write volume data 
        do k=1,self%ldim(3)
            cnt = 0
            do j=self%indextab(k,1),self%indextab(k,2)
                cnt = cnt+1
                do i=1,self%ldim(1)
                    self%pixrow(i) = img%get([i,cnt,k])
                end do
                write(self%filnum, rec=j) self%pixrow
            end do
        end do
    end subroutine
    
    ! CALCULATORS
    
    !>  \brief  is for stack/volume copying
    function copy( self, name, nptcls, fromto ) result( self_copy )
        use simple_image, only: image
        class(spider), intent(inout)  :: self
        character(len=*), intent(in)  :: name
        integer, intent(in), optional :: nptcls
        integer, intent(in), optional :: fromto(2)
        type(spider) :: self_copy
        type(image) :: img
        integer :: n, i, cnt
        if( .not. self%is_open() ) stop 'need to open file to be copied; copy; simple_spider'
        call self_copy%new(self%ldim)
        self_copy%spistkhed  = self%spistkhed
        self_copy%spiimghed  = self%spiimghed
        self_copy%indextab   = self%indextab
        self_copy%endconv    = self%endconv
        self_copy%recsz      = self%recsz
        self_copy%lenbyt     = self%lenbyt
        self_copy%labrec     = self%labrec
        self_copy%labbyt     = self%labbyt
        self_copy%ldim       = self%ldim
        if( present(nptcls) ) then
            n = nptcls
        else
            n = self%get_nptcls()
        endif
        call self_copy%open(name, 'replace')
        call img%new(self%ldim,1.)
        if( self%ldim(3) == 1 )then
            if( n >= 1 )then
                write(*,'(a)') '>>> COPYING IMAGES'
                if( present(fromto) )then
                    cnt = 0
                    do i=fromto(1),fromto(2)
                        cnt = cnt+1
                        call progress(cnt, fromto(2)-fromto(1)+1)
                        call self%read(img, i)
                        call self_copy%write(img, cnt)
                    end do
                else                
                    do i=1,n
                        call progress(i, n)
                        call self%read(img, i)
                        call self_copy%write(img, i)
                    end do
                endif
            endif
        else
            call self%read(img)
            call self_copy%write(img)
        endif
        call img%kill
    end function
    
    !>  \brief  is for making a stack of normalized vectors for pca/DBN analysis    
    subroutine make_pattern_stack( self, mskrad, name, D, recsz, avg, hfun )
        use simple_image,  only: image
        use simple_jiffys, only: get_fileunit
        use simple_stat,   only: normalize, normalize_net
        class(spider), intent(inout)             :: self
        real, intent(in)                         :: mskrad
        character(len=*), intent(in)             :: name
        integer, intent(out)                     :: D, recsz
        real, allocatable, intent(out), optional :: avg(:)
        character(len=*), intent(in), optional   :: hfun
        type(image)                              :: img
        real, allocatable                        :: pcavec(:)
        integer                                  :: n, fnum, ier, i, alloc_stat
        logical                                  :: err
        if( .not. self%is_open() ) stop 'need to open stack first; make_pattern_stack; simple_spider'
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        D = img%get_npix(mskrad)
        allocate(pcavec(D), stat=alloc_stat)
        call alloc_err('make_pattern_stack; simple_spider, 1', alloc_stat)
        pcavec = 0.
        inquire(iolength=recsz) pcavec
        deallocate(pcavec)
        if( present(avg) )then
            allocate(avg(D), stat=alloc_stat)
            call alloc_err('make_pattern_stack; simple_spider, 2', alloc_stat)
            avg = 0.
        endif
        fnum = get_fileunit()
        open(unit=fnum, status='replace', action='readwrite', file=name,&
        access='direct', form='unformatted', recl=recsz, iostat=ier)
        call fopen_err('make_pattern_stack; simple_spider', ier)
        write(*,'(a)') '>>> MAKING PATTERN STACK'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            if( img%is_ft() ) call img%bwd_ft
            call img%serialize(mskrad, pcavec)
            err = .false.
            if( present(hfun) )then
                call normalize_net(pcavec, hfun)
            else
                call normalize(pcavec, err)
            endif
            if( err ) write(*,'(a,i7)') 'WARNING: variance zero! image nr: ', i 
            if( present(avg) ) avg = avg+pcavec
            write(fnum,rec=i) pcavec
            deallocate(pcavec)
        end do
        if( present(avg) )then
            avg = avg/real(n)
            allocate(pcavec(D), stat=alloc_stat)
            call alloc_err('make_pattern_stack; simple_spider, 3', alloc_stat)
            do i=1,n
                read(fnum,rec=i) pcavec
                pcavec = pcavec-avg
                write(fnum,rec=i) pcavec
            end do
        endif
        close(unit=fnum)
        call img%kill
    end subroutine
    
    !>  \brief  is for stack padding
    subroutine pad( self, ldim, outstk )
        use simple_image, only: image
        class(spider), intent(inout) :: self, outstk
        integer, intent(in)          :: ldim(3)
        type(image)  :: img, img_pad
        integer      :: n, i
        if( .not. self%is_open() ) stop 'need to open stack to be padded; pad; simple_spider'
        if( .not. outstk%is_open() ) stop 'need to open ouput stack; pad; simple_spider'
        if( ldim(1) >= self%ldim(1) .and. ldim(2) >= self%ldim(2)&
        .and. ldim(3) >= self%ldim(3) )then
            call img%new(self%ldim,1.)
            call img_pad%new(ldim,1.)
            n = self%get_nptcls()
            write(*,'(a)') '>>> PADDING IMAGES'
            do i=1,n
                call progress(i,n)
                call self%read(img, i)
                call img%pad(img_pad) ! FT state preserved
                call outstk%write(img_pad, i)
            end do
        end if
        call img%kill
        call img_pad%kill
    end subroutine
    
    !>  \brief  is for stack clipping
    subroutine clip( self, ldim, outstk )
        use simple_image, only: image
        class(spider), intent(inout) :: self, outstk
        integer, intent(in)          :: ldim(3)
        type(image) :: img, img_clip
        integer     :: n, i
        if( .not. self%is_open() ) stop 'need to open stack to be clipped; clip; simple_spider'
        if( .not. outstk%is_open() ) stop 'need to open output stack; clip; simple_spider'
        if( ldim(1) <= self%ldim(1) .and. ldim(2) <= self%ldim(2)&
        .and. ldim(3) <= self%ldim(3) )then
            call img%new(self%ldim,1.)
            call img_clip%new(ldim,1.)
            n = self%get_nptcls()
            write(*,'(a)') '>>> CLIPPING IMAGES'
            do i=1,n
                call progress(i,n)
                call self%read(img, i)
                call img%clip(img_clip) ! FT state preserved
                call outstk%write(img_clip, i)
            end do
        end if
        call img%kill
        call img_clip%kill
    end subroutine
    
    !>  \brief  is for stack resizing
    subroutine resize( self, ldim, outstk )
        use simple_image, only: image
        class(spider), intent(inout) :: self, outstk
        integer, intent(in)          :: ldim(3)
        type(image) :: img, img_resized
        integer     :: n, i
        if( .not. self%is_open() ) stop 'need to open stack to be resized; resize; simple_spider'
        if( .not. outstk%is_open() ) stop 'need to open output stack; resize; simple_spider'
        call img%new(self%ldim,1.)
        call img_resized%new(ldim,1.)
        n = self%get_nptcls()
        write(*,'(a)') '>>> RESIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%fwd_ft
            if( ldim(1) <= self%ldim(1) .and. ldim(2) <= self%ldim(2)&
            .and. ldim(3) <= self%ldim(3) )then
                call img%clip(img_resized)
            else
                call img%pad(img_resized)
            endif
            call img_resized%bwd_ft
            call outstk%write(img_resized, i)
        end do
        call img%kill
        call img_resized%kill
    end subroutine
    
    !>  \brief  is for stack normalization
    subroutine norm( self, outstk, hfun )
        use simple_image, only: image
        class(spider), intent(inout)           :: self
        class(spider), intent(inout), optional :: outstk
        character(len=*), intent(in), optional :: hfun
        type(image) :: img
        integer     :: i, n
        if( .not. self%is_open() ) stop 'need to open stack to be normalized; norm; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; norm; simple_spider'
        endif
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        write(*,'(a)') '>>> NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%norm(hfun)
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        call img%kill
    end subroutine
    
    !>  \brief  is for stack normalization
    subroutine stats( self, which, ave, sdev, var, med, msk )
        use simple_image, only: image
        class(spider), intent(inout) :: self
        character(len=*), intent(in) :: which
        real, intent(out)            :: ave, sdev, var, med
        real, intent(in), optional   :: msk
        real        :: ave_loc, sdev_loc, var_loc, med_loc
        type(image) :: img
        integer     :: i, n
        if( .not. self%is_open() ) stop 'need to open stack to be analyzed; stats; simple_spider'
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        write(*,'(a)') '>>> CALCULATING STACK STATISTICS'
        ave  = 0.
        sdev = 0.
        var  = 0.
        med  = 0.
        do i=1,n
            call progress(i,n)
            call self%read(img,i)
            call img%stats(which, ave_loc, sdev_loc, var_loc, med_loc, msk=msk) 
            ave  = ave  + ave_loc
            sdev = sdev + sdev_loc
            var  = var  + var_loc
            med  = med  + med_loc
        end do
        ave  = ave/real(n)
        sdev = sdev/real(n)
        var  = var/real(n)
        med  = med/real(n)
        call img%kill
    end subroutine
    
    !>  \brief  is for stack normalization
    subroutine neg( self, outstk )
        use simple_image, only: image
        class(spider), intent(inout)           :: self
        class(spider), intent(inout), optional :: outstk
        type(image) :: img
        integer     :: i, n
        if( .not. self%is_open() ) stop 'need to open stack to be contrast inverted; neg; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; neg; simple_spider'
        endif
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        write(*,'(a)') '>>> INVERTING IMAGE CONTRAST'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%neg
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        call img%kill
    end subroutine
    
    !>  \brief  is for stack rank transformation
    subroutine rankify( self, outstk )
        use simple_image, only: image
        class(spider), intent(inout)           :: self
        class(spider), intent(inout), optional :: outstk
        type(image) :: img
        integer     :: i, n
        if( .not. self%is_open() ) stop 'need to open stack to be rankified; norm; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; rankify; simple_spider'
        endif
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        write(*,'(a)') '>>> RANKIFYING IMAGES'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%fwd_ft
            call img%rankify
            call img%bwd_ft
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        call img%kill
    end subroutine
    
    !>  \brief  is for centering based on center of mass
    subroutine masscen( self, smpd, lp, msk, tres, outstk )
        use simple_image, only: image
        class(spider), intent(inout)           :: self
        real, intent(in)                       :: smpd, lp
        real, intent(in), optional             :: msk, tres
        class(spider), intent(inout), optional :: outstk
        type(image) :: img
        integer     :: i, n
        if( .not. self%is_open() ) stop 'need to open stack to be processed; masscen; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; masscen; simple_spider'
        endif
        n = self%get_nptcls()
        call img%new(self%ldim,smpd)
        write(*,'(a)') '>>> CENTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%center(lp, msk, tres)
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        call img%kill
    end subroutine
    
    !>  \brief  is for stack curing
    subroutine cure( self, outstk )
        use simple_image,  only: image
        use simple_jiffys, only: get_fileunit
        class(spider), intent(inout) :: self
        class(spider), intent(inout), optional :: outstk
        type(image) :: img
        integer     :: i, n, n_nans, filnum
        real        :: maxv, minv, ave, sdev
        if( .not. self%is_open() ) stop 'need to open stack to be cured; cure; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; cure; simple_spider'
        endif
        filnum = get_fileunit()
        open(unit=filnum, status='replace', file='cure_stats.txt')
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        write(*,'(a)') '>>> CURING IMAGES'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%cure(maxv, minv, ave, sdev, n_nans)
            write(filnum,'(A,1X,f12.2,1X,A,1X,f12.2,1X,A,1X,f12.2,1X,A,1X,f12.2,1X,A,1X,I9)') 'MAX:', maxv,&
            'MIN:', minv, 'AVE:', ave, 'SDEV:', sdev, 'NANS:', n_nans
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        close(filnum)
        call img%kill
    end subroutine
    
     !>  \brief  is for calculating the acf of a stack
    subroutine acf( self, outstk )
        use simple_image, only: image
        class(spider), intent(inout)           :: self
        class(spider), intent(inout), optional :: outstk
        type(image) :: img
        integer     :: i, n
        if( .not. self%is_open() ) stop 'need to open stack to be acf:ed; acf; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; acf; simple_spider'
        endif
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        write(*,'(a)') '>>> CALCULATING ACF:S OF THE IMAGES'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%acf
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        call img%kill
    end subroutine
    
    !>  \brief  is for adding noise
    subroutine add_noise( self, snr, outstk )
        use simple_image, only: image
        class(spider), intent(inout)           :: self
        real, intent(in)                       :: snr
        class(spider), intent(inout), optional :: outstk
        type(image) :: img
        integer :: i, n
        if( .not. self%is_open() ) stop 'need to open stack to add noise; add_noise; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; add_noise; simple_spider'
        endif
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        write(*,'(a)') '>>> ADDING NOISE TO THE IMAGES'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%add_gauran(snr)
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        call img%kill
    end subroutine
    
    !>  \brief  is for making frame averages of dose-fractionated image series
    subroutine frameavg( self_ptcls, outstk, navg )
        use simple_image, only: image
        class(spider), intent(inout) :: self_ptcls, outstk
        integer, intent(in) :: navg
        type(image) :: img, avg
        integer :: i, n, cnt, cnt2
        if( .not. self_ptcls%is_open() ) stop 'need to open ptcl stack 4 frame-averaging; frameavg; simple_spider'
        if( .not. outstk%is_open()  ) stop 'need to open avg  stack 4 frame-averaging; frameavg; simple_spider'
        n = self_ptcls%get_nptcls()
        call img%new(self_ptcls%ldim,1.)
        call avg%new(self_ptcls%ldim,1.)
        cnt = 0
        cnt2 = 0
        write(*,'(a)') '>>> AVERAGING FRAMES'
        do i=1,n
            call progress(i,n)
            cnt = cnt+1
            call self_ptcls%read(img, i)
            if( cnt <= navg )then
                call avg%add(img)
            endif
            if( cnt == navg )then
                cnt2 = cnt2+1
                call outstk%write(avg, cnt2)
                cnt = 0
                avg = 0.
            endif
        end do
        call img%kill
        call avg%kill
    end subroutine
    
    !>  \brief  is for stack Fourier transformation
    subroutine fwd_ft( self )
        use simple_image, only: image
        class(spider), intent(inout) :: self
        type(image) :: img
        integer     :: n, i
        if( self%ldim(3) > 1 ) stop 'not for volumes; fwd_ft; simple_spider'
        if( .not. self%is_open() ) stop 'need to open stack to be transformed; fwd_ft; simple_spider'
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        do i=1,n
            call self%read(img, i)
            if( .not. img%is_ft() )then
                call img%fwd_ft
                call self%write(img, i)
            endif 
        end do
        call img%kill
        ! write stack header:
        call self%make_headers( self%ldim, .true. )
        self%spistkhed(1,26) = self%get_nptcls()
        do i=1,self%labrec
            write(self%filnum, rec=i) self%spistkhed(i,:)
        end do
    end subroutine
    
    !>  \brief  is for stack Fourier transformation visualization
    subroutine ft2img( self, which, outstk )
        use simple_image, only: image
        class(spider), intent(inout)           :: self
        character(len=*), intent(in)           :: which
        class(spider), intent(inout), optional :: outstk
        type(image) :: img, img2
        integer     :: n, i
        if( self%ldim(3) > 1 ) stop 'not for volumes; ft2img; simple_spider'
        if( .not. self%is_open() ) stop 'need to open stack to be transformed; ft2img; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; ft2img; simple_spider'
        endif
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        do i=1,n
            call self%read(img, i)
            img2 = img%ft2img(which)
            if( present(outstk) )then
                call outstk%write(img2, i)
            else
                call self%write(img2, i)
            endif
        end do
        call img%kill
        call img2%kill
    end subroutine
    
    !>  \brief  is for stack back Fourier transformation
    subroutine bwd_ft( self )
        use simple_image, only: image
        class(spider), intent(inout) :: self
        type(image) :: img
        integer     :: n, i
        if( self%ldim(3) > 1 ) stop 'not for volumes; bwd_ft; simple_spider'
        if( .not. self%is_open() ) stop 'need to open stack to be transformed; bwd_ft; simple_spider'
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        do i=1,n
            call self%read(img, i)
            if( img%is_ft() )then
                call img%bwd_ft
                call self%write(img, i)
            endif 
        end do
        call img%kill
        ! write stack header:
        call self%make_headers( self%ldim, .false. )
        self%spistkhed(1,26) = self%get_nptcls()
        do i=1,self%labrec
            write(self%filnum, rec=i) self%spistkhed(i,:)
        end do
    end subroutine
    
    !>  \brief  is for origin shifting a stack according to info in o
    subroutine shift( self, o, smpd, mul, round, outstk )
        use simple_image, only: image
        use simple_oris,  only: oris
        class(spider), intent(inout)           :: self
        class(oris), intent(inout)             :: o
        real, intent(in)                       :: smpd
        real, intent(in), optional             :: mul
        logical, intent(in), optional          :: round
        class(spider), intent(inout), optional :: outstk
        type(image) :: img
        integer     :: n, i
        real        :: x, y, xhere, yhere
        logical     :: rround
        if( .not. self%is_open() ) stop 'need to open stack to be shifted; shift; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; shift; simple_spider'
        endif
        rround = .false.
        if( present(round) ) rround = round
        n = self%get_nptcls()
        if( n /= o%get_noris() ) stop 'inconsistent nr entries; shift; simple_spider'
        call img%new(self%ldim,smpd)
        write(*,'(a)') '>>> SHIFTING IMAGES'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            x = o%get(i, 'x')
            y = o%get(i, 'y')
            if( present(mul) )then
                xhere = x*mul
                yhere = y*mul
            else
                xhere = x
                yhere = y
            endif
            if( rround )then
                xhere = real(nint(xhere))
                yhere = real(nint(yhere))
            endif
            call img%shift(-xhere, -yhere)
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        call img%kill
    end subroutine
    
    !>  \brief  is for filtering of a spider stack
    subroutine bp( self, smpd, hp, lp, outstk )
        use simple_image, only: image
        class(spider), intent(inout)           :: self
        real, intent(in)                       :: smpd, hp, lp
        class(spider), intent(inout), optional :: outstk
        type(image) :: img
        integer     :: n, i
        if( .not. self%is_open() ) stop 'need to open stack to be filtered; bp; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; bp; simple_spider'
        endif
        n = self%get_nptcls()
        call img%new(self%ldim,smpd)
        write(*,'(a)') '>>> BAND-PASS FILTERING IMAGES'   
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%bp(hp,lp,width=10.)
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        call img%kill
    end subroutine
    
    !>  \brief  is for phase randomization of a spider stack
    subroutine phase_rand( self, smpd, lp, outstk )
        use simple_image, only: image
        class(spider), intent(inout)           :: self
        real, intent(in)                       :: smpd, lp
        class(spider), intent(inout), optional :: outstk
        type(image)                  :: img
        integer                      :: n, i
        if( .not. self%is_open() ) stop 'need to open stack to be filtered; phase_rand; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; phase_rand; simple_spider'
        endif
        n = self%get_nptcls()
        call img%new(self%ldim,smpd)
        write(*,'(a)') '>>> PHASE RANDOMIZING IMAGES'     
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%phase_rand(lp)
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        call img%kill
    end subroutine
    
    !>  \brief  is for phase flipping a stack according to info in o
    subroutine apply_ctf( self, o, smpd, tfun, mode, bfac, outstk )
        use simple_math,  only: deg2rad
        use simple_image, only: image
        use simple_oris,  only: oris
        use simple_ctf,   only: ctf
        class(spider), intent(inout)           :: self   !< instance
        class(oris), intent(inout)             :: o      !< orientations object
        real, intent(in)                       :: smpd   !< sampling distance, a/pix
        class(ctf), intent(inout)              :: tfun   !< transfer function 
        character(len=*), intent(in)           :: mode   !< abs, ctf, flip, flipneg, neg, square
        real, intent(in), optional             :: bfac   !< bfactor
        class(spider), intent(inout), optional :: outstk !< output stack
        type(image) :: img
        integer     :: i, n
        if( .not. self%is_open() ) stop 'need to open stack to be modified; apply_ctf; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; apply_ctf; simple_spider'
        endif
        n = self%get_nptcls()
        if( n /= o%get_noris() ) stop 'inconsistent nr entries; apply_ctf; simple_spider'
        call img%new(self%ldim,smpd)
        write(*,'(a)') '>>> APPLYING CTF TO IMAGES' 
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%fwd_ft
            if( o%isthere('dfy') )then ! astigmatic CTF   
                call tfun%apply(img, o%get(i,'dfx'), mode, dfy=o%get(i,'dfy'), angast=o%get(i,'angast'), bfac=bfac)
            else ! non-astigmatic CTF
                call tfun%apply(img, o%get(i,'dfx'), mode, bfac=bfac)
            endif
            call img%bwd_ft
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
    end subroutine
    
    !>  \brief  is for visualization of the CTF power spectrum
    subroutine ctf2img( self, o, smpd, tfun, bfac )
        use simple_math,  only: deg2rad
        use simple_image, only: image
        use simple_oris,  only: oris
        use simple_ctf,   only: ctf
        class(spider), intent(inout) :: self !< instance
        class(oris), intent(inout)   :: o    !< orientations object
        real, intent(in)             :: smpd !< sampling distance, a/pix
        class(ctf), intent(inout)    :: tfun !< transfer function 
        real, intent(in), optional   :: bfac !< bfactor
        type(image) :: img, img2
        integer     :: i, n
        n = o%get_noris()
        if( .not. self%is_open() ) stop 'need to open stack to write to it; apply_ctf; simple_spider'
        call img%new(self%ldim,smpd)
        call img2%new(self%ldim,smpd)
        write(*,'(a)') '>>> MAKING CTF POWSPEC IMAGES' 
        do i=1,n
            call progress(i,n)
            if( o%isthere('dfy') )then ! astigmatic CTF   
                call tfun%ctf2img(img, o%get(i,'dfx'), 'square', dfy=o%get(i,'dfy'), angast=o%get(i,'angast'), bfac=bfac)
            else
                call tfun%ctf2img(img, o%get(i,'dfx'), 'square', bfac=bfac)
            endif
            img2 = img%ft2img('real')
            call self%write(img2, i)
        end do
        call img%kill
        call img2%kill
    end subroutine    
    
    !>  \brief  is for origin shifting and rotating a stack
    !!          according to info in o 
    subroutine shrot( self, o, smpd, mul, outstk )
        use simple_image, only: image
        use simple_oris,  only: oris
        class(spider), intent(inout)           :: self
        class(oris), intent(inout)             :: o
        real, intent(in)                       :: smpd
        real, intent(in), optional             :: mul
        class(spider), intent(inout), optional :: outstk
        type(image) :: img, img_rot
        integer     :: n, i
        real        :: x, y
        if( .not. self%is_open() ) stop 'need to open stack to be rotated; shrot; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; shrot; simple_spider'
        endif
        n = self%get_nptcls()       
        if( n /= o%get_noris() ) stop 'inconsistent nr entries; shrot; simple_spider'
        call img%new(self%ldim,smpd)
        call img_rot%new(self%ldim,smpd)
        write(*,'(a)') '>>> SHIFTING AND ROTATING IMAGES' 
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%fwd_ft
            x = o%get(i, 'x')
            y = o%get(i, 'y')
            if( present(mul) )then
                call img%shift(-x*mul, -y*mul)
            else
                call img%shift(-x, -y)
            endif
            call img%bwd_ft
            call img%rtsq(-o%e3get(i), 0., 0., img_rot)
            if( present(outstk) )then
                call outstk%write(img_rot, i)
            else
                call self%write(img_rot, i)
            endif
        end do
        call img%kill
        call img_rot%kill
    end subroutine
    
    !>  \brief  is for applying a soft circular mask to all images in stack
    subroutine mask( self, mskrad, inner, width, outstk )
        use simple_image, only: image
        class(spider), intent(inout)           :: self
        real, intent(in)                       :: mskrad
        real, intent(in), optional             :: inner, width
        class(spider), intent(inout), optional :: outstk
        type(image) :: img
        integer     :: n, i
        if( .not. self%is_open() ) stop 'need to open stack to be masked; mask; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; mask; simple_spider'
        endif
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        write(*,'(a)') '>>> MASKING IMAGES' 
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call img%mask(mskrad, 'soft', inner, width) ! FT state preserved 
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        call img%kill
    end subroutine
    
    !>  \brief  is for binarizing all images in the stack
    subroutine bin( self, tres, outstk )
        use simple_image, only: image
        class(spider), intent(inout)           :: self
        real, intent(in), optional             :: tres
        class(spider), intent(inout), optional :: outstk
        type(image) :: img
        integer     :: n, i
        logical     :: didft
        if( .not. self%is_open() ) stop 'need to open stack to be binarized; bin; simple_spider'
        if( present(outstk) )then
            if( .not. outstk%is_open() ) stop 'need to open output stack; bin; simple_spider'
        endif
        n = self%get_nptcls()
        call img%new(self%ldim,1.)
        write(*,'(a)') '>>> BINARIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            didft = .false.
            if( img%is_ft() )then
                call img%bwd_ft
                didft = .true.
            endif
            if( present(tres) )then
                call img%bin(tres)
            else
                call img%bin
            endif
            if( didft ) call img%fwd_ft
            if( present(outstk) )then
                call outstk%write(img, i)
            else
                call self%write(img, i)
            endif
        end do
        call img%kill
    end subroutine
    
    !>  \brief  is for making an average of the entire stack
    subroutine make_avg( self, name )
        use simple_image, only: image
        class(spider), intent(inout) :: self
        character(len=*), intent(in) :: name
        type(spider)                 :: avgstk
        type(image)                  :: avg, img
        integer                      :: i, n
        call avg%new(self%ldim, 1.)
        call img%new(self%ldim, 1.)
        avg = 0.
        n = self%get_nptcls()
        write(*,'(a)') '>>> MAKING GLOBAL STACK AVERAGE'
        do i=1,n
            call progress(i,n)
            call self%read(img, i)
            call avg%add(img)       
        end do
        call avg%div(real(n))
        call avgstk%new(self%ldim)
        call avgstk%open(name, 'replace')
        call avgstk%write(avg,1)
        call avgstk%kill
    end subroutine
    
    !>  \brief  is for making class averages according to the 
    !!          classification info in o
    subroutine make_cavgs( self, o, minp, name, self_cls, list )
        use simple_image, only: image
        use simple_oris,  only: oris
        use simple_sll,   only: sll
        class(spider), intent(inout)     :: self
        class(oris), intent(inout)       :: o
        integer, intent(in)              :: minp
        character(len=*), intent(in)     :: name
        type(spider), intent(out)        :: self_cls
        type(sll), intent(out), optional :: list
        type(image)                      :: avg
        logical                          :: empty
        integer :: i, ncls, ncls_here, pop
        ncls = o%get_ncls()
        call self_cls%new(self%ldim)
        call self_cls%open(name, 'replace')
        if( present(list) ) call list%new
        ncls_here = 0
        write(*,'(a)') '>>> MAKING CLASS AVERAGES'
        do i=1,ncls
            call progress(i,ncls)
            pop = o%get_clspop(i)
            if( pop >= minp )then
                if( present(list) ) call list%add(iarr=[i])
                ncls_here = ncls_here+1
                call self%make_cavg(o, i, avg, empty )
                call self_cls%write(avg, ncls_here)
            endif        
        end do
    end subroutine
    
    !>  \brief  is for making one class average according to the 
    !!          classification info in o and the inputted class
    subroutine make_cavg( self, o, class, avg, empty, smpd )
        use simple_image, only: image
        use simple_oris,  only: oris
        class(spider), intent(inout) :: self
        class(oris), intent(inout)   :: o
        integer, intent(in)          :: class
        type(image), intent(inout)   :: avg
        logical, intent(out)         :: empty
        real, intent(in), optional   :: smpd
        type(image)                  :: img
        real                         :: ssmpd
        integer :: j, ncls, pop, n, mycls
        if( .not. self%is_open() ) stop 'need to open particle stack; make_cavg; simple_spider'
        n = self%get_nptcls()
        if( n /= o%get_noris() ) stop 'inconsistent nr entries; make_cavg; simple_spider'
        ssmpd = 1.
        if( present(smpd) ) ssmpd = smpd
        ncls = o%get_ncls()
        pop = o%get_clspop(class) 
        empty = .false.
        if( pop == 0 )then
            empty = .true.
        else
            call avg%new(self%ldim, ssmpd)
            call img%new(self%ldim, ssmpd)
            avg = 0.         
            do j=1,n
                mycls = nint(o%get(j, 'class'))
                if( mycls == class )then
                    call self%read(img, j)
                    call avg%add(img)
                endif
            end do
            call avg%div(real(pop))
            call img%kill
        endif
    end subroutine
    
    ! PRINTERS/GETTERS
    
    !>  \brief  is for printing the image header
    subroutine print_img_hed( self, nr )
        class(spider), intent(in) :: self
        integer, intent(in)       :: nr
        real                      :: hed(self%labrec,self%ldim(1))
        integer                   :: i, hedinds(2), iminds(2), cnt
        if( .not. self%is_open() ) stop 'need to open file for printing; print_img_hed; simple_spider'
        write(*,'(a)') ">>> IMAGE HEADER"
        call self%imgnr2recpos( nr, hedinds, iminds )
        cnt = 0
        do i=hedinds(1),hedinds(2)
            cnt = cnt+1
            read(self%filnum, rec=i) hed(cnt,:)
        end do
        do i=1,self%ldim(1)
            write(*,*) 'header pos:', i, 'val:', self%spiimghed(1,i)
        end do
    end subroutine

    !>  \brief  is for printing the stack/volume header
    subroutine print_hed( self )
        class(spider), intent(in) :: self
        real                      :: hed(self%labrec,self%ldim(1))
        integer                   :: i
        if( .not. self%is_open() ) stop 'need to open file for printing; print_hed; simple_spider'
        write(*,'(a)') ">>> STACK/VOLUME HEADER"
        do i=1,self%labrec
            read(self%filnum, rec=i) hed(i,:)
        end do
        do i=1,self%ldim(1)
            write(*,*) 'header pos:', i, 'val:', self%spistkhed(1,i)
        end do
    end subroutine
    
    !>  \brief  is for getting the number of ptcls in the stack
    function get_nptcls( self ) result( n )
        class(spider), intent(in) :: self
        real                      :: hed(self%labrec,self%ldim(1))
        integer                   :: n
        if( .not. self%is_open() ) stop 'need to open file for getting nr of ptcls; get_nptcls; simple_spider'
        read(self%filnum, rec=1) hed(1,:)
        if( hed(1,1) > 1 )then
            n = 1
        else
            n = nint(hed(1,26))
        endif
    end function
    
    !>  \brief  is for getting the logical dimension of the stack
    function get_ldim( self ) result( ldim )
        class(spider), intent(in) :: self
        integer                   :: ldim(3)
        ldim = self%ldim
    end function

    ! PRIVATE STUFF
    
    !>  \brief  is for printing the internal image header
    subroutine print_internal_img_hed( self )
        class(spider) :: self
        integer       :: i
        write(*,'(a)') ">>> IMAGE HEADER"
        do i=1,self%ldim(1)
            write(*,*) 'header pos:', i, 'val:', self%spiimghed(1,i)
        end do
    end subroutine
    
    !>  \brief  is for making the stack/volume header
    subroutine make_headers( self, ldim, is_ft )
    !*(1) nz, number of slices (planes) in volume (=1 for an image) (ldim(3))
    !*(2) ny, number of rows per slice (ldim(2))
    !*(3) irec, total number of records (including header records) in each image of a simple image or stacked image file (ny+labrec=ldim(2)+labrec)
    !*(5) 3=3D vol, -11=2D FT odd, -12=2D FT even, -21=3D FT odd, -22=3D FT even
    ! (6) 1/0 indicates if max,min,avg,stdev have been computed
    ! (7) maximum data val
    ! (8) minimum data val
    ! (9) average data val
    !*(10) standard dev -1 indicates not computed 
    !*(12) nx, number of pixels per line (ldim(1))
    !*(13) labrec (nr of records in file header)
    ! (14) 1/0 indicates existence or not of Euler angles in header
    ! (15) e1
    ! (16) e2
    ! (17) e3
    ! (18) xsh
    ! (19) ysh
    ! (20) zsh
    !*(21) scale (default=1)
    !*(22) labbyt (total nr of bytes in header)
    !*(23) lenbyt (record lenght in bytes)
!    *(24) A value > 0 indicates a stack of images
    !*(26) Highest number in stack (1 for volumes)
    ! (27) Position of image in stack, or 0 if not used
    ! THERE ARE OTHER PARAMS AS WELL BUT THIS IS ENOUGH 4 NOW
        use simple_math,  only: is_even
        class(spider), intent(inout) :: self
        integer, intent(in)          :: ldim(3)
        logical, intent(in)          :: is_ft
        logical                      :: even_dims, is_3d
        if( ldim(3) == 1 )then
            even_dims = is_even(ldim(:2))
        else
            even_dims = is_even(ldim)
        endif
        is_3d = .false.
        if( ldim(3) > 1 ) is_3d = .true.
        self%spistkhed = 0.
        self%spistkhed(1,1) = real(self%ldim(3)) 
        self%spistkhed(1,2) = real(self%ldim(2)) 
        self%spistkhed(1,3) = real(self%ldim(2)+self%labrec)
        if( .not. is_3d .and. .not. is_ft )then
            self%spistkhed(1,5) = 1.
        else if( is_3d .and. .not. is_ft )then
            self%spistkhed(1,5) = 3.
        else if( .not. is_3d .and. is_ft .and. .not. even_dims )then
            self%spistkhed(1,5) = -11.
        else if( .not. is_3d .and. is_ft .and. even_dims )then
            self%spistkhed(1,5) = -12.
        else if( is_3d .and. is_ft .and. .not. even_dims )then
            self%spistkhed(1,5) = -21.
        else if( is_3d .and. is_ft .and. even_dims )then
            self%spistkhed(1,5) = -22.
        else
            stop 'undefined file type specifier, 1; make_header; simple_spider'
        endif
        self%spistkhed(1,10) = -1.
        self%spistkhed(1,12) = real(ldim(1))
        self%spistkhed(1,13) = real(self%labrec)
        self%spistkhed(1,21) = 1.
        self%spistkhed(1,22) = real(self%labbyt)
        self%spistkhed(1,23) = real(self%lenbyt)
        if( ldim(3) == 1 ) self%spistkhed(1,24) = 2. ! stack header
        if( ldim(3) > 1 )then
            self%spistkhed(1,26) = 1.
        else
            self%spistkhed(1,26) = 0. ! maximum nr of images in stack
        !******! (26) need to be updated on write (stack)
        endif
        ! make an image header
        self%spiimghed = 0.
        if( ldim(3) == 1 )then
            self%spiimghed(1,1) = 1.                   
            self%spiimghed(1,2) = real(self%ldim(2))
            self%spiimghed(1,3) = real(self%ldim(2)+self%labrec)
            if( .not. is_ft )then
                self%spiimghed(1,5) = 1.
            else if( is_ft .and. .not. even_dims )then
                self%spiimghed(1,5) = -11.
            else if( is_ft .and. even_dims )then
                self%spiimghed(1,5) = -12.
            else
                stop 'undefined file type specifier, 1; make_header; simple_spider'
            endif
            self%spiimghed(1,10) = -1.
            self%spiimghed(1,12) = real(ldim(1))
            self%spiimghed(1,13) = real(self%labrec)
            self%spiimghed(1,22) = real(self%labbyt)
            self%spiimghed(1,23) = real(self%lenbyt)
            self%spiimghed(1,27) = 0.  ! indicates position of image in stack or 0 if image is not used
            !******! (27) need to be updated on write (image)
        endif
    end subroutine
    
    !>  \brief  for translating an image index to record indices in the stack
    subroutine imgnr2recpos( self, nr, hedinds, iminds )
        class(spider), intent(in) :: self 
        integer, intent(in)       :: nr
        integer, intent(out)      :: hedinds(2), iminds(2)
        integer                   :: cnt, j
        cnt = self%labrec
        do j=1,nr
            hedinds(1) = cnt+1   ! hed from
            cnt = cnt+self%labrec 
            hedinds(2) = cnt     ! hed to
            iminds(1) = cnt+1    ! im from
            cnt = cnt+self%ldim(2)
            iminds(2) = cnt      ! im to
        end do
    end subroutine
    
    !>  \brief for finding endianness conversion of spider file
    function find_endconv( self, fname, box ) result( foundbox )
        use simple_jiffys, only: get_fileunit
        class(spider), intent(inout)  :: self
        character(len=*), intent(in)  :: fname
        integer, intent(in), optional :: box
        integer                       :: i, ier, filnum
        logical                       :: here, foundbox
        inquire(FILE=fname, EXIST=here)
        self%endconv = ''
        call check( 'LITTLE_ENDIAN' )
        call check( 'BIG_ENDIAN' )
        foundbox = .false.
        if( present(box) )then
            if( self%endconv == 'LITTLE_ENDIAN' .or. self%endconv == 'BIG_ENDIAN' )then
                foundbox = .true.
            endif
        else
            if( self%endconv == '' )then
                write(*,'(a)') 'ERROR, could not determine endianess of file!'
                write(*,'(a)') 'In: find_endconv, module: simple_spider'
                write(*,'(a)') 'File corrupt or logical dimension incorrect?'
                stop
            else
                foundbox = .true.
            endif
        endif
        
        contains
            
            subroutine check( str )
                character(len=*), intent(in) :: str
                filnum = get_fileunit( )
                open(unit=filnum, convert=str, status='old', action='read', file=fname,&
                access='direct', form='unformatted', recl=self%recsz, iostat=ier)
!                call fopen_err('find_endconv; simple_spider', ier)
                ! read header
                do i=1,self%labrec
                    read(filnum, rec=i) self%spistkhed(i,:)
                end do
                close(filnum)
!                print *, 'endianness:', str
!                print *, 'hed(1,1):', self%spistkhed(1,1), 'ldim(3):', self%ldim(3)
!                print *, 'hed(1,2):', self%spistkhed(1,2), 'ldim(2):', self%ldim(2)
                if( present(box) )then
                    if( int(self%spistkhed(1,2)) == box ) self%endconv = str
                else
                    if( int(self%spistkhed(1,1)) == self%ldim(3) .and.&
                    int(self%spistkhed(1,2)) == self%ldim(2) ) self%endconv = str
                endif
                close(unit=filnum)
            end subroutine
            
    end function
    
    
    ! DOES NOT WORK WITH GFORTRAN
    !>  \brief for determining endianness conversion of spider file
!    subroutine determine_endconv( self, fname )
!        class(spider), intent(inout)  :: self
!        character(len=*), intent(in)  :: fname
!        inquire(FILE=fname, CONVERT=self%endconv) !***CONVERT
!        print *, 'determined endconv:', self%endconv
!    end subroutine
    
    ! UNIT TEST
    
    subroutine test_spider
        write(*,'(a)') '**info(simple_spider_unit_test): testing square dimensions'
        call test_spider_local(100,100,100)
        write(*,'(a)') '**info(simple_spider_unit_test): testing non-square dimensions'
        call test_spider_local(120,90,80)
        write(*,'(a)') 'SIMPLE_SPIDER_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine
    
    subroutine test_spider_local( ld1, ld2, ld3 )
        use simple_image, only: image
        use simple_oris,  only: oris
        integer, intent(in) :: ld1, ld2, ld3
        type(spider) :: spistk, spistk_pad, spistk_sh, spistk_clip
        type(spider) :: spistk_cavgs, spivol, spistk_copy, spistk_rot, spistk_cmeds
        type(image)  :: img, img_read, vol, vol_read
        type(oris)   :: o
        real         :: corr1, corr2, corr
        logical      :: passed
        integer      :: i, ld_here(3)
        write(*,'(a)') '**info(simple_spider_unit_test, part 1): testing writing/reading spider volumes'
        passed = .false.
        call vol%new([ld1,ld2,ld3], 2.)
        call vol%square(20)
        call vol_read%copy(vol)
        call spivol%new(vol%get_ldim())
        call spivol%open('test_spider_cube.spi', 'replace')
        call spivol%write(vol)
        call spivol%read(vol)
        call vol%fwd_ft
        call vol_read%fwd_ft
        corr = vol%corr(vol_read)
        if( corr > 0.99 ) passed = .true.
        if( .not. passed ) stop 'writing/reading spider volumes test failed'
        call spivol%open('test_spider_ftcube.spi', 'replace')
        ! EMAN cannot read the ftcube file, but Spider produces a (98x100,100) cube upon 
        ! reverse transformation if we exclude the shifting
        call spivol%write(vol)
        call spivol%read(vol_read)
        corr = vol%corr(vol_read)
        if( corr > 0.99 ) passed = .true.
        if( .not. passed ) stop 'writing/reading FTed volumes with the spider class failed'
        call vol_read%bwd_ft
        call spivol%open('test_spider_ft_rev_ftcube.spi', 'replace')
        call spivol%write(vol_read)
        call spivol%kill
        write(*,'(a)') '**info(simple_spider_unit_test, part 2): testing writing/reading images to/from spider stacks'
        call img%new([ld1,ld2,1], 2.)
        call img_read%copy(img)
        call img%square(20)
        call spistk%new(img%get_ldim())
        call spistk%open('test_spider_squares.spi', 'replace')
        do i=1,10
            call spistk%write( img, i )
        end do
        call img%fwd_ft
        do i=1,10
            passed = .false.
            call spistk%read( img_read, i )
            call img_read%fwd_ft
            if( img%corr(img_read) > 0.99 ) passed = .true.
            if( .not. passed ) stop 'writing/reading images to spider stack test failed'
        end do
        write(*,'(a)') '**info(simple_spider_unit_test, part 3): testing writing/reading 2D FTs to/from spider stacks'
        call spistk%open('test_spider_ftsquares.spi', 'replace')
        do i=1,10
            call spistk%write( img, i )
        end do
        do i=1,10
            passed = .false.
            call spistk%read( img_read, i )
            corr = img%corr(img_read)
            if( corr > 0.99 ) passed = .true.
            if( .not. passed ) stop 'writing/reading 2D FTs to spider stack test failed'
        end do
        passed = .false.
        if( spistk%get_nptcls() == 10 ) passed = .true.     
        write(*,'(a)') '**info(simple_spider_unit_test, part 4): testing stack calculators'
        call img%square(20)
        call spistk%open('test_spider_squares.spi')
        spistk_copy = spistk%copy('test_spider_squares_copy.spi')
        do i=1,10
            passed = .false.
            call spistk_copy%read( img_read, i )
            if( img%corr(img_read) > 0.99 ) passed = .true.
            if( .not. passed ) stop 'stack copy test failed'
        end do       
        ld_here = img%get_ldim()
        ld_here(1) = 2*ld_here(1)
        ld_here(2) = 2*ld_here(2)
        call spistk_pad%new(ld_here)
        call spistk_clip%new(img%get_ldim())
        call spistk_pad%open('test_spider_padded.spi','replace')
        call spistk_clip%open('test_spider_clipped.spi','replace')
        call spistk%pad(ld_here,outstk=spistk_pad)
        call spistk_pad%clip(img%get_ldim(),outstk=spistk_clip)
        do i=1,10
            passed = .false.
            call spistk_clip%read(img_read, i)
            if( img%corr(img_read) > 0.99 ) passed = .true.
            if( .not. passed ) stop 'stack pad/clip test failed'
        end do
        call spistk_clip%fwd_ft
        call img%fwd_ft
        do i=1,10
            passed = .false.
            call spistk_clip%read(img_read, i)
            if( .not. img_read%is_ft() ) stop 'fwd transform test failed, 1'
            if( img%corr(img_read) > 0.99 ) passed = .true.
            if( .not. passed ) stop 'fwd transform test failed, 2'
        end do
        call img%bwd_ft
        call spistk_clip%bwd_ft
        do i=1,10
            passed = .false.
            call spistk_clip%read(img_read, i)
            if( img_read%is_ft() ) stop 'bwd transform test failed, 1'
            if( img%corr(img_read) > 0.99 ) passed = .true.
            if( .not. passed ) stop 'bwd transform test failed, 2'
        end do
        call o%new(10)
        call o%rnd_oris([[0.,360.],[0.,180.],[0.,360.]], 5.)
        call img%gauimg(20)
        call spistk%open('test_spider_gaussians_shifted.spi', 'replace')
        do i=1,10
            call spistk%write(img,i)
        end do       
        spistk_sh = spistk%copy( 'test_spider_gaussians.spi' )
        call spistk%shift(o, 2.)
        call o%revshsgn
        call spistk%shift(o, 2.)
        call spistk%open('test_spider_squares.spi')
        if( ld1 == ld2 )then
            spistk_rot = spistk%copy('test_spider_squares_rot.spi')
            call o%zero_shifts
            call spistk_rot%shrot(o, 2.)       
            call o%revorisgn
            call spistk_rot%shrot(o, 2.)
            call img%square(20)
            do i=1,10
                passed = .false.
                call spistk%read(img_read,i)
                if( img%corr(img_read) > 0.99 ) passed = .true.
                if( .not. passed ) stop 'rot stack test failed'
            end do
        endif
        call img%ran
        call spistk%open('test_spider_masked.spi', 'replace')
        do i=1,10
            call spistk%write(img,i)
        end do
        call spistk%mask(40.)
        call spistk%open('test_spider_ptcls.spi', 'replace')
        do i=1,100
            call img%square(20)
            call img%add_gauran(0.1)
            call spistk%write(img,i)
        end do
        do i=101,200
            call img%gauimg(10)
            call img%add_gauran(0.1)
            call spistk%write(img,i)
        end do
        call o%new(200)
        do i=1,100
            call o%set(i,'class',1.)
        end do
        do i=101,200
            call o%set(i,'class',2.)
        end do
        call spistk%make_cavgs(o, 1, 'test_spider_classavgs.spi', spistk_cavgs)
        passed = .false.
        call spistk_cavgs%read(img_read,1)
        call img%square(20)
        corr1 = img%corr(img_read,10.)
        call spistk_cavgs%read(img_read,2)
        call img%gauimg(10)
        corr2 = img%corr(img_read,10.)
        if( corr1+corr2 > 1. ) passed = .true.
        if( .not. passed ) stop 'making classavgs test failed'
        call spistk_cavgs%close
        call spistk%kill
        call spistk_pad%kill
        call spistk_sh%kill
        call spistk_clip%kill
        call spistk_cavgs%kill
        call spistk_cmeds%kill
        call spivol%kill
        call spistk_copy%kill
        call spistk_rot%kill
        call img%kill
        call img_read%kill
        call vol%kill
        call vol_read%kill
    end subroutine
    
    ! DESTRUCTOR
    
    !>  \brief is a destructor
    subroutine kill( self )
        class(spider), intent(inout) :: self
        if( self%existence )then
            deallocate( self%spistkhed, self%spiimghed, self%pixrow, self%indextab )
            if( self%is_open() ) call self%close
            self%existence = .false.
        endif
    end subroutine
    
end module simple_spider
