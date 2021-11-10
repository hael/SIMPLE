module simple_eer_factory
use, intrinsic :: iso_c_binding
use simple_defs
use simple_strings
use simple_syslib
use simple_fileio
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image,      only: image
use simple_map_reduce, only: split_nobjs_even
#ifdef USING_TIFF
use simple_tifflib
#endif
implicit none

public :: eer_decoder
private
#include "simple_local_flags.inc"

integer, parameter :: EER_IMAGE_WIDTH  = 4096
integer, parameter :: EER_IMAGE_HEIGHT = 4096
integer, parameter :: EER_IMAGE_PIXELS = EER_IMAGE_WIDTH * EER_IMAGE_HEIGHT
integer, parameter :: TIFF_COMPRESSION_EER8bit = 65000
integer, parameter :: TIFF_COMPRESSION_EER7bit = 65001

type eer_decoder
    private
    type(c_ptr)                   :: fhandle = c_null_ptr       !< File handle
    character(len=:), allocatable :: fname                      !< File name
    integer(kind=1),  allocatable :: raw_bytes(:)               !< parsed raw bytes
    integer(kind=8),  allocatable :: frames_start(:), frames_sz(:) !< temporal dimensions
    real                          :: smpd          = 0.         !< physical pixel size
    real                          :: osmpd         = 0.         !< output pixel size
    integer(kind=8)               :: filesz        = 0          !< file size, upper bound to # of electrons
    integer                       :: nx, ny        = 0          !< base dimensions
    integer                       :: onx, ony      = 0          !< physical output dimensions
    integer                       :: upsampling    = 1          !< desired sampling: 1=4K, 2=8K
    integer                       :: nframes       = 0          !< # of frames
    integer(kind=8)               :: nmax_el       = 0          !< max # of electrons
    logical                       :: l_hasbeenread = .false.    !< whether file has been parsed
    logical                       :: l_7bit        = .false.    !< compression
    logical                       :: l_exists      = .false.    !< whether object exists
    contains
        ! Constructor
        procedure          :: new
        procedure, private :: read
        ! Doers
        procedure          :: decode
        procedure, private :: decode_frames
        procedure          :: prep_gainref
        ! Setter
        procedure          :: set_dims
        ! Getter
        procedure          :: get_nframes
        procedure          :: get_ldim
        procedure          :: get_smpd_out
        procedure          :: get_smpd
        ! Destructor
        procedure          :: kill
end type eer_decoder

contains

    !>  Constructor & parsing
    subroutine new( self, fname, smpd, upsampling)
        class(eer_decoder), intent(inout) :: self
        character(len=*),   intent(in)    :: fname
        real,               intent(in)    :: smpd       ! always refers to physical pixel size
        integer,            intent(in)    :: upsampling
        integer :: compression
        call self%kill
#ifdef USING_TIFF
        if( .not.file_exists(fname) ) THROW_HARD('File could not be found: '//trim(fname))
        self%fname = trim(fname)
        inquire(file=self%fname,size=self%filesz)
        ! open
        call TIFFMuteWarnings
        self%fhandle = TIFFOpen(toCstring(self%fname),toCstring('r'))
        call TIFFUnMuteWarnings
        ! header
        self%nx      = TIFFGetWidth(self%fhandle)
        self%ny      = TIFFGetLength(self%fhandle)
        compression  = TIFFGetCompression(self%fhandle)
        self%nframes = TIFFNumDirectories(self%fhandle)
        ! dimensions & sampling
        if( self%nx /= EER_IMAGE_WIDTH ) THROW_HARD('Invalid width: '//int2str(self%nx))
        if( self%ny /= EER_IMAGE_HEIGHT )THROW_HARD('Invalid height: '//int2str(self%nx))
        if( self%nframes < 1 ) THROW_HARD('Invalid # of frames: '//int2str(self%nframes))
        self%smpd = smpd
        call self%set_dims(upsampling)
        ! Compression
        select case(compression)
            case(TIFF_COMPRESSION_EER8bit)
                self%l_7bit = .false.
            case(TIFF_COMPRESSION_EER7bit)
                self%l_7bit = .true.
            case DEFAULT
                THROW_HARD('Invalid compression')
        end select
        ! final  init
        allocate(self%frames_start(0:self%nframes-1),self%frames_sz(0:self%nframes-1), source=int(0,kind=8))
        self%l_hasbeenread = .false.
        self%l_exists = .true.
        call self%read
#else
        THROW_HARD('EER support requires TIFF support!')
#endif
    end subroutine new

    !>  Parse & store raw bytes
    subroutine read( self )
        class(eer_decoder), intent(inout) :: self
#ifdef USING_TIFF
        type(c_ptr)               :: stripbuffer_ptr
        integer(kind=1),  pointer :: byte_array8(:)
        integer(kind=8)           :: pos, stripsz8
        integer                   :: iframe, iostat, nstrips, istrip, stripsz, nbytes
        stripbuffer_ptr = TIFFAllocateStripBuffer(self%fhandle)
        call TIFFMuteWarnings
        allocate(self%raw_bytes(0:self%filesz-1), source=int(0,kind=1))
        pos = int(0,kind=8)
        do iframe = 0,self%nframes-1
            iostat  = TIFFSetDirectory(self%fhandle, int(iframe,kind=c_int16_t))
            nstrips = TIFFNumberOfStrips(self%fhandle)
            self%frames_start(iframe) = pos
            self%frames_sz(iframe)    = int(0,kind=8)
            do istrip = 0,nstrips-1
                stripsz  = TIFFRawStripSizer(self%fhandle, istrip)
                stripsz8 = int(stripsz,kind=8)
                nbytes   = TIFFReadRawStrip(self%fhandle, istrip, stripbuffer_ptr, stripsz)
                call c_f_pointer(stripbuffer_ptr, byte_array8, [nbytes])
                self%raw_bytes(pos:pos+nbytes-1) = byte_array8(:)
                pos = pos + stripsz8
                self%frames_sz(iframe) = self%frames_sz(iframe) + stripsz8
            enddo
        enddo
        call TIFFUnMuteWarnings
        self%nmax_el = maxval(self%frames_sz)
        ! iostat = TIFFfree(stripbuffer_ptr)
        call TIFFClose(self%fhandle)
        nullify(byte_array8)
        self%l_hasbeenread = .true.
#endif
    end subroutine

    !>  reset sampling, output dimensions & pixel size
    subroutine set_dims( self, upsampling )
        class(eer_decoder), intent(inout) :: self
        integer,            intent(in) :: upsampling
        select case(upsampling)
        case(1)
            self%onx   = self%nx
            self%ony   = self%ny
            self%osmpd = self%smpd
        case(2)
            self%onx   = 2 * self%nx
            self%ony   = 2 * self%ny
            self%osmpd = self%smpd / 2.
        case(3)
            self%onx   = 4 * self%nx
            self%ony   = 4 * self%ny
            self%osmpd = self%smpd / 4.
        case DEFAULT
            THROW_HARD('Invalid upsampling: '//int2str(upsampling))
        end select
        self%upsampling = upsampling
    end subroutine set_dims

    ! DOERS

    subroutine decode( self, imgs, fraction )
        class(eer_decoder), intent(in)    :: self
        class(image),       intent(inout) :: imgs(:)
        integer,            intent(in)    :: fraction
        integer :: i,nimgs
        if( .not.self%l_exists ) THROW_HARD('EER decoder object has not been instantiated! decode')
        nimgs  = size(imgs)
        if( nimgs*fraction > self%nframes )then
            THROW_HARD('EER fraction is too large! decode')
        endif
        ! here adopting relion's convention for fractionating where the last overhanging
        ! raw frames (and most dose-damaged) are simply abandoned
        !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
        do i = 1,nimgs
            call imgs(i)%kill
            call self%decode_frames(imgs(i), (i-1)*fraction+1, i*fraction)
        enddo
        !$omp end parallel do
    end subroutine decode

    subroutine decode_frames( self, img, istart, iend )
        class(eer_decoder), intent(in)    :: self
        class(image),       intent(inout) :: img
        integer,            intent(in)    :: istart, iend ! base 1
#ifdef USING_TIFF
        integer(kind=2), parameter :: one2 = int(1,kind=2)
        real,              pointer :: prmat(:,:,:)
        integer(kind=1),   pointer :: byte_array(:)
        type(c_ptr)     :: cptr
        integer         :: symbols(self%nmax_el), pos_els(self%nmax_el)
        integer(kind=2) :: imat(self%onx,self%ony)
        integer(kind=8) :: pos, first_byte
        integer         :: iframe, n_el, p1,s1,p2,s2, x,y, i, iostat, bit_pos, pos_el
        if( .not.self%l_exists ) THROW_HARD('EER decoder object has not been instantiated!')
        if( .not.self%l_hasbeenread ) THROW_HARD('File needs be read prior to decoding!')
        cptr = TIFFmalloc(4)
        call c_f_pointer(cptr, byte_array, [4])
        imat = int(0,kind=2)
        do iframe = istart-1,iend-1
            pos     = self%frames_start(iframe)
            pos_els = 0
            pos_el  = 0
            n_el    = 0
            ! works out indexed positions
            if( self%l_7bit )then
                ! 7-bit: 2 * (4+7) in 2*(2*8)
                bit_pos = 0
                do while( .true. )
                    first_byte      = pos + shiftr(bit_pos,3)
                    byte_array(1:4) = self%raw_bytes(first_byte:first_byte+3)
                    call EERDecode_7bit(cptr, bit_pos, p1, s1, p2, s2)
                    ! First electron
                    bit_pos = bit_pos + 7
                    pos_el  = pos_el  + p1
                    if( pos_el >= EER_IMAGE_PIXELS ) exit
                    if( p1 == 127 ) cycle
                    bit_pos       = bit_pos + 11
                    n_el          = n_el    + 1
                    pos_els(n_el) = pos_el
                    symbols(n_el) = s1
                    ! Second electron
                    pos_el        = pos_el  + p2 + 1
                    if( pos_el >= EER_IMAGE_PIXELS ) exit
                    if( p2  == 127 ) cycle
                    bit_pos       = bit_pos + 4
                    n_el          = n_el    + 1
                    pos_els(n_el) = pos_el
                    symbols(n_el) = s2
                    pos_el        = pos_el  + 1
                enddo
            else
                ! 8-bit: 2 * (4+8) in 3*8; untested
                do while( pos < self%frames_start(iframe) + self%frames_sz(iframe))
                    call EERDecode_8bit(self%raw_bytes(pos),self%raw_bytes(pos+1),self%raw_bytes(pos+2), p1, s1, p2, s2)
                    pos_el = pos_el + p1
                    if( pos_el >= EER_IMAGE_PIXELS ) exit
                    if( p1 < 255 )then
                        n_el          = n_el + 1
                        pos_els(n_el) = pos_el
                        symbols(n_el) = s1
                        pos_el        = pos_el + 1
                    endif
                    pos_el = pos_el + p2
                    if( pos_el >= EER_IMAGE_PIXELS ) exit
                    if( p2 < 255 )then
                        n_el          = n_el + 1
                        pos_els(n_el) = pos_el
                        symbols(n_el) = s2
                        pos_el        = pos_el + 1
                    endif
                    pos = pos + 3
                enddo
            endif
            ! indexed positions to coordinates & counts
            select case(self%upsampling)
            case(1)
                do i = 1,n_el
                    call EERdecodePos4K(pos_els(i), x,y)
                    imat(x,y) = imat(x,y) + one2
                enddo
            case(2)
                do i = 1,n_el
                    call EERdecodePos8K(pos_els(i), symbols(i), x,y)
                    imat(x,y) = imat(x,y) + one2
                enddo
            case(3)
                do i = 1,n_el
                    call EERdecodePos16K(pos_els(i), symbols(i), x,y)
                    imat(x,y) = imat(x,y) + one2
                enddo
            end select
        enddo
        ! generates image
        call img%new([self%onx,self%ony,1], self%osmpd)
        call img%zero_and_unflag_ft
        call img%get_rmat_ptr(prmat)
        prmat(:self%onx,:self%ony,1) = real(imat)
        ! cleanup
        nullify(prmat,byte_array)
        iostat = TIFFfree(cptr)
#endif
    end subroutine decode_frames

    subroutine prep_gainref( self, fname, gain )
        use simple_math,    only: is_zero
        use simple_imghead, only: find_ldim_nptcls
        class(eer_decoder), intent(in)    :: self
        character(len=*),   intent(in)    :: fname
        class(image),       intent(inout) :: gain
        type(image)   :: tmp
        real, pointer :: prmat(:,:,:)
        real          :: avg, val
        integer       :: ldim_gain(3), ifoo, i,j, ii,jj
        logical       :: dotgain
        dotgain = (fname2format(fname) == 'J') .or. (fname2format(fname) == 'L')
        call find_ldim_nptcls(fname,ldim_gain,ifoo)
        ! gain always at desired sampling of frames
        call gain%new([self%onx,self%ony,1],self%osmpd)
        if( ldim_gain(1)==self%onx .and. ldim_gain(2)==self%ony )then
            ! gain dimensions = desired frames dimensions
            call gain%read(fname)
            if( dotgain ) call flipY(gain)
        else
            call tmp%new(ldim_gain, 1.)
            call tmp%read(fname)
            if( dotgain ) call flipY(tmp)
            call gain%zero
            if( ldim_gain(1)==EER_IMAGE_WIDTH .and. self%upsampling==2 )then
                ! gain 4K, images 8K
                !$omp parallel do default(shared) private(i,j,ii,jj) proc_bind(close) schedule(static) collapse(2)
                do j = 1,2*EER_IMAGE_HEIGHT
                    do i = 1,2*EER_IMAGE_WIDTH
                        jj = floor(real(j-1)/2.)+1
                        ii = floor(real(i-1)/2.)+1
                        call gain%set([i,j,1], tmp%get([ii,jj,1]))
                    enddo
                enddo
                !$omp end parallel do
            else if( ldim_gain(1)==2*EER_IMAGE_WIDTH .and. self%upsampling==1 )then
                ! gain 8K, images 4K
                call tmp%get_rmat_ptr(prmat)
                !$omp parallel do default(shared) private(i,j,ii,jj) proc_bind(close) schedule(static) collapse(2)
                do j = 1,self%ony
                    do i = 1,self%onx
                        jj  = 2*(j-1) + 1
                        ii  = 2*(i-1) + 1
                        call gain%set([i,j,1], sum(prmat(ii:ii+1,jj:jj+1,1)))
                    enddo
                enddo
                !$omp end parallel do
                nullify(prmat)
            else
                THROW_HARD('Unsupported combination of eer sampling & gain reference dimensions; prep_gainref')
            endif
            call tmp%kill
        endif
        if( .not.dotgain )then
            ! .mrc, taking inverse times average
            call gain%get_rmat_ptr(prmat)
            !$omp workshare
            avg = real(sum(real(prmat(:self%onx,:self%ony,1),dp)) / real(self%onx*self%ony,dp))
            where( is_zero(prmat(:self%onx,:self%ony,1)) )
                ! zeroes are preserved and dealt with at outliers curation level
            else where
                prmat(:self%onx,:self%ony,1) = avg / prmat(:self%onx,:self%ony,1)
            end where
            !$omp end workshare
        endif

        contains

            subroutine flipY( img )
                class(image), intent(inout) :: img
                call img%get_rmat_ptr(prmat)
                do j = 1,ldim_gain(2)/2
                    jj = ldim_gain(2) - j - 1
                    do i = 1,ldim_gain(1)
                        val           = prmat(i,j,1)
                        prmat(i,j,1)  = prmat(i,jj,1)
                        prmat(i,jj,1) = val
                    enddo
                enddo
                nullify(prmat)
            end subroutine flipY

    end subroutine prep_gainref


    ! SETTERS

    integer function get_nframes( self )
        class(eer_decoder), intent(in) :: self
        get_nframes = self%nframes
    end function get_nframes

    function get_ldim( self )result( ldim )
        class(eer_decoder), intent(in)  :: self
        integer :: ldim(3)
        ldim(1) = self%onx
        ldim(2) = self%ony
        ldim(3) = 1
    end function get_ldim

    real function get_smpd_out( self )
        class(eer_decoder), intent(in) :: self
        get_smpd_out = self%osmpd
    end function get_smpd_out

    real function get_smpd( self )
        class(eer_decoder), intent(in) :: self
        get_smpd = self%smpd
    end function get_smpd

    !>  Destructor
    subroutine kill( self )
        class(eer_decoder), intent(inout) :: self
        self%fhandle = c_null_ptr
        self%fname   = ''
        if( allocated(self%raw_bytes) )    deallocate(self%raw_bytes)
        if( allocated(self%frames_start) ) deallocate(self%frames_start)
        if( allocated(self%frames_sz) )    deallocate(self%frames_sz)
        self%smpd          = 0.
        self%osmpd         = 0.
        self%nx            = 0
        self%ny            = 0
        self%onx           = 0
        self%ony           = 0
        self%upsampling    = 1
        self%nframes       = 0
        self%nmax_el       = 0
        self%l_hasbeenread = .false.
        self%l_7bit        = .false.
        self%l_exists      = .false.
    end subroutine kill

end module simple_eer_factory
