!@descr: versioned dose-series manifest and validation of its unmasked real-space half maps
module simple_dose_series_manifest
use simple_defs,   only: dp, logfhandle
use simple_string, only: string
use simple_syslib, only: simple_abspath, file_exists, simple_atomic_replace, simple_sync_file
use simple_imghead, only: MrcImgHead, MRCHEADSZ, MRC_MODE_FLOAT32, MRC_MODE_FLOAT16
use, intrinsic :: iso_fortran_env, only: iostat_end, iostat_eor, int8, int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
private
public :: dose_series_manifest

character(len=*), parameter :: VERSION_LINE = '# SIMPLE dose series 1'
integer, parameter :: LINE_CHUNK = 1024, MAX_LINE = 65536, MAX_PATH = 4095
real(dp), parameter :: SMPD_ATOL = 1.e-4_dp, VALUE_RTOL = 1.e-5_dp, KV_ATOL = 1.e-3_dp
real(dp), parameter :: ANGLE_ATOL = 1.e-3_dp

type :: dose_window
    type(string) :: even, odd
    integer :: first = 0, last = 0
    real(dp) :: lo = 0._dp, hi = 0._dp
end type dose_window

type :: dose_series_manifest
    private
    type(dose_window), allocatable :: windows(:)
    type(string) :: digest
    integer :: box = 0, state = 1, nptcls = 0, neven = 0, nodd = 0
    real(dp) :: smpd = 0._dp, kv = 0._dp, dose_per_fraction = 0._dp
    logical :: paths_resolved = .false.
  contains
    procedure :: new
    procedure :: add_window
    procedure :: read
    procedure :: write
    procedure :: select_windows
    procedure :: get_geometry
    procedure :: get_nwindows
    procedure :: get_window
    procedure :: kill
    procedure, private :: validate
    procedure, private :: validate_files
end type dose_series_manifest

contains

    !> Initialize metadata; validation is deferred until publication so raw outputs can survive a bad dose.
    subroutine new(self, box, smpd, kv, dpf, nptcls, neven, nodd, digest)
        class(dose_series_manifest), intent(inout) :: self
        integer, intent(in) :: box, nptcls, neven, nodd
        real(dp), intent(in) :: smpd, kv, dpf
        class(string), intent(in) :: digest
        call self%kill
        self%box = box
        self%smpd = smpd
        self%kv = kv
        self%dose_per_fraction = dpf
        self%nptcls = nptcls
        self%neven = neven
        self%nodd = nodd
        self%digest = digest
    end subroutine new

    !> Append one actual interval, retaining a possible final overlap for the consumer to select.
    subroutine add_window(self, a, b, dlo, dhi, even, odd)
        class(dose_series_manifest), intent(inout) :: self
        integer, intent(in) :: a, b
        real(dp), intent(in) :: dlo, dhi
        class(string), intent(in) :: even, odd
        type(dose_window), allocatable :: grown(:)
        integer :: n
        n = self%get_nwindows()
        allocate(grown(n+1))
        if( n > 0 ) grown(:n) = self%windows
        grown(n+1)%first = a
        grown(n+1)%last = b
        grown(n+1)%lo = dlo
        grown(n+1)%hi = dhi
        grown(n+1)%even = even
        grown(n+1)%odd = odd
        call move_alloc(grown, self%windows)
    end subroutine add_window

    !> Read and validate metadata and the actual MRC headers without terminating the caller.
    logical function read(self, fname, message) result(ok)
        class(dose_series_manifest), intent(inout) :: self
        class(string), intent(in) :: fname
        type(string), intent(out) :: message
        type(string) :: directory, line, tokens(6), even, odd
        character(len=:), allocatable :: txt, key, val
        logical :: seen(10), first, parsed
        integer :: unit, ios, ntok, slot, a, b
        real(dp) :: dlo, dhi
        call self%kill
        ok = .false.
        message = 'Cannot open dose series manifest: '//fname%to_char()
        open(newunit=unit, file=fname%to_char(), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        directory = parent_directory(fname)
        seen = .false.
        first = .true.
        do
            call read_line(unit, line, ios)
            if( ios == iostat_end ) exit
            if( ios /= 0 )then
                message = 'Unreadable or overlong manifest line'
                exit
            endif
            txt = trim(adjustl(line%to_char()))
            if( len(txt) == 0 ) cycle
            if( first )then
                first = .false.
                if( txt /= VERSION_LINE )then
                    message = 'Unknown dose series manifest version'
                    exit
                endif
                cycle
            endif
            if( txt(1:1) == '#' )then
                call split_tokens(txt(2:), tokens, ntok)
                if( ntok < 1 ) cycle
                key = tokens(1)%to_char()
                slot = header_slot(key)
                if( slot == 0 ) cycle
                if( ntok /= 2 )then
                    message = 'Manifest header needs one value: '//key
                    exit
                endif
                if( seen(slot) )then
                    message = 'Duplicate manifest header: '//key
                    exit
                endif
                seen(slot) = .true.
                val = tokens(2)%to_char()
                parsed = .true.
                select case(slot)
                    case(1)
                        parsed = val == 'refine_motion_model'
                    case(2)
                        parsed = parse_int(val, self%state)
                    case(3)
                        parsed = parse_int(val, self%box)
                    case(4)
                        parsed = parse_real(val, self%smpd)
                    case(5)
                        parsed = parse_real(val, self%kv)
                    case(6)
                        parsed = parse_int(val, self%nptcls)
                    case(7)
                        parsed = parse_int(val, self%neven)
                    case(8)
                        parsed = parse_int(val, self%nodd)
                    case(9)
                        self%digest = val
                    case(10)
                        parsed = parse_real(val, self%dose_per_fraction)
                end select
                if( .not. parsed )then
                    message = 'Invalid manifest header: '//key
                    exit
                endif
            else
                call split_tokens(txt, tokens, ntok)
                if( ntok /= 6 )then
                    message = 'Manifest row needs exactly six whitespace-separated columns'
                    exit
                endif
                parsed = parse_int(tokens(1)%to_char(), a)
                if( parsed ) parsed = parse_int(tokens(2)%to_char(), b)
                if( parsed ) parsed = parse_real(tokens(3)%to_char(), dlo)
                if( parsed ) parsed = parse_real(tokens(4)%to_char(), dhi)
                if( .not. parsed )then
                    message = 'Invalid fraction or dose in manifest row'
                    exit
                endif
                even = tokens(5)
                odd = tokens(6)
                call self%add_window(a, b, dlo, dhi, even, odd)
            endif
        enddo
        close(unit)
        if( ios /= iostat_end )then
            call self%kill
            return
        endif
        if( first )then
            message = 'Dose series manifest has no version line'
        else if( .not. all(seen(:10)) )then
            message = 'Dose series manifest is missing required metadata'
        else
            ok = self%validate(message)
            if( ok ) ok = self%validate_files(directory, message)
        endif
        if( .not. ok )then
            call self%kill
            return
        endif
        ! Resolve paths before params%new can change the working directory.
        do a = 1, self%get_nwindows()
            self%windows(a)%even = resolve_path(directory, self%windows(a)%even)
            self%windows(a)%odd  = resolve_path(directory, self%windows(a)%odd)
        enddo
        self%paths_resolved = .true.
        message = ''
    end function read

    !> Atomically publish all rows only after the complete metadata and referenced pair set validate.
    logical function write(self, fname, message) result(ok)
        class(dose_series_manifest), intent(in) :: self
        class(string), intent(in) :: fname
        type(string), intent(out) :: message
        type(string) :: directory, tmp, even, odd
        integer :: unit, ios, close_status, i
        ok = self%validate(message)
        if( .not. ok ) return
        directory = parent_directory(fname)
        ok = self%validate_files(directory, message)
        if( .not. ok ) return
        ok = .false.
        tmp = fname//'.tmp'
        message = 'Cannot create manifest temporary file: '//tmp%to_char()
        open(newunit=unit, file=tmp%to_char(), status='new', action='write', iostat=ios)
        if( ios /= 0 ) return
        write(unit,'(A)',iostat=ios) VERSION_LINE
        if( ios == 0 ) write(unit,'(A)',iostat=ios) '# program refine_motion_model'
        if( ios == 0 ) write(unit,'(A,I0)',iostat=ios) '# state ', self%state
        if( ios == 0 ) write(unit,'(A,I0)',iostat=ios) '# box ', self%box
        if( ios == 0 ) write(unit,'(A,ES24.16E3)',iostat=ios) '# smpd ', self%smpd
        if( ios == 0 ) write(unit,'(A,ES24.16E3)',iostat=ios) '# kv ', self%kv
        if( ios == 0 ) write(unit,'(A,I0)',iostat=ios) '# nptcls ', self%nptcls
        if( ios == 0 ) write(unit,'(A,I0)',iostat=ios) '# neven ', self%neven
        if( ios == 0 ) write(unit,'(A,I0)',iostat=ios) '# nodd ', self%nodd
        if( ios == 0 ) write(unit,'(A)',iostat=ios) '# eo_digest '//self%digest%to_char()
        if( ios == 0 ) write(unit,'(A,ES24.16E3)',iostat=ios) '# dose_per_fraction ', self%dose_per_fraction
        if( ios == 0 ) write(unit,'(A)',iostat=ios) '# columns from to dose_lo dose_hi even odd'
        do i = 1, self%get_nwindows()
            if( ios /= 0 ) exit
            even = portable_path(directory, self%windows(i)%even)
            odd  = portable_path(directory, self%windows(i)%odd)
            if( .not. valid_path(even) .or. .not. valid_path(odd) )then
                ios = 1
                exit
            endif
            write(unit,'(2(I0,1X),2(ES24.16E3,1X),A,1X,A)',iostat=ios) self%windows(i)%first, &
                &self%windows(i)%last, self%windows(i)%lo, self%windows(i)%hi, even%to_char(), odd%to_char()
        enddo
        if( ios /= 0 )then
            close(unit, status='delete', iostat=close_status)
            message = 'Cannot write dose series manifest'
            return
        endif
        close(unit, iostat=ios)
        if( ios == 0 ) call simple_sync_file(tmp, ios)
        if( ios == 0 ) call simple_atomic_replace(tmp, fname, ios)
        if( ios /= 0 )then
            call remove_temporary(tmp)
            message = 'Cannot publish dose series manifest: '//fname%to_char()
            return
        endif
        message = ''
        ok = .true.
    end function write

    !> Select disjoint windows, keeping the full-width last interval when it overlaps its predecessor.
    logical function select_windows(self, message) result(ok)
        class(dose_series_manifest), intent(inout) :: self
        type(string), intent(out) :: message
        type(dose_window), allocatable :: selected(:)
        integer :: n, drop
        ok = self%validate(message)
        if( .not. ok ) return
        n = self%get_nwindows()
        drop = overlap_predecessor(self)
        if( drop > 0 )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> DOSE SERIES: WINDOW ', self%windows(drop)%first, '-', &
                &self%windows(drop)%last, ' DROPPED, IT OVERLAPS THE LAST WINDOW'
            allocate(selected(n-1))
            if( drop > 1 ) selected(:drop-1) = self%windows(:drop-1)
            selected(drop) = self%windows(n)
            call move_alloc(selected, self%windows)
        endif
    end function select_windows

    !> Return the native geometry and physical voltage recorded by the producer.
    subroutine get_geometry(self, box, smpd, kv)
        class(dose_series_manifest), intent(in) :: self
        integer, intent(out) :: box
        real(dp), intent(out) :: smpd, kv
        box = self%box
        smpd = self%smpd
        kv = self%kv
    end subroutine get_geometry

    !> Return the number of recorded or selected windows.
    integer function get_nwindows(self) result(n)
        class(dose_series_manifest), intent(in) :: self
        n = 0
        if( allocated(self%windows) ) n = size(self%windows)
    end function get_nwindows

    !> Return a window; callers iterate from one through get_nwindows.
    subroutine get_window(self, i, a, b, dlo, dhi, even, odd)
        class(dose_series_manifest), intent(in) :: self
        integer, intent(in) :: i
        integer, intent(out) :: a, b
        real(dp), intent(out) :: dlo, dhi
        type(string), intent(out) :: even, odd
        a = self%windows(i)%first
        b = self%windows(i)%last
        dlo = self%windows(i)%lo
        dhi = self%windows(i)%hi
        even = self%windows(i)%even
        odd = self%windows(i)%odd
    end subroutine get_window

    !> Release interval storage and reset the metadata.
    subroutine kill(self)
        class(dose_series_manifest), intent(inout) :: self
        if( allocated(self%windows) ) deallocate(self%windows)
        call self%digest%kill
        self%box = 0
        self%state = 1
        self%nptcls = 0
        self%neven = 0
        self%nodd = 0
        self%smpd = 0._dp
        self%kv = 0._dp
        self%dose_per_fraction = 0._dp
        self%paths_resolved = .false.
    end subroutine kill

    !> Check numeric metadata, interval doses and the disjoint subset without changing any rows.
    logical function validate(self, message) result(ok)
        class(dose_series_manifest), intent(in) :: self
        type(string), intent(out) :: message
        character(len=:), allocatable :: hash
        real(dp) :: lo, hi
        integer :: i, n, drop, previous
        ok = .false.
        message = 'Dose series needs a positive even box and state 1'
        if( self%box < 2 .or. mod(self%box,2) /= 0 .or. self%state /= 1 ) return
        message = 'Dose series sampling or dose per fraction is not finite and positive'
        if( .not. ieee_is_finite(self%smpd) .or. .not. ieee_is_finite(self%dose_per_fraction) ) return
        if( self%smpd <= 0._dp .or. self%dose_per_fraction <= 0._dp ) return
        message = 'Dose series voltage must be 300, 200 or 100 kV'
        if( .not. ieee_is_finite(self%kv) ) return
        if( minval(abs(self%kv-[100._dp,200._dp,300._dp])) > KV_ATOL ) return
        message = 'Dose series particle/half counts are inconsistent or empty'
        if( self%neven < 1 .or. self%nodd < 1 .or. self%nptcls < 2 ) return
        if( int(self%neven,int64)+int(self%nodd,int64) /= int(self%nptcls,int64) ) return
        message = 'Dose series partition fingerprint must contain 16 hexadecimal digits'
        hash = self%digest%to_char()
        if( len(hash) /= 16 ) return
        if( verify(hash,'0123456789abcdefABCDEF') /= 0 ) return
        n = self%get_nwindows()
        message = 'A dose series needs at least two windows'
        if( n < 2 ) return
        do i = 1, n
            message = 'Invalid dose series fraction range'
            if( self%windows(i)%first < 1 .or. self%windows(i)%last < self%windows(i)%first ) return
            message = 'Dose series has nonfinite or nonpositive dose span'
            lo = self%windows(i)%lo
            hi = self%windows(i)%hi
            if( .not. ieee_is_finite(lo) .or. .not. ieee_is_finite(hi) ) return
            if( lo < 0._dp .or. hi <= lo ) return
            message = 'Dose per fraction overflows the accumulated exposure'
            if( self%dose_per_fraction > huge(1._dp)/real(self%windows(i)%last,dp) ) return
            message = 'Dose bounds disagree with original fraction indices and dose per fraction'
            if( .not. close_value(lo,real(self%windows(i)%first-1,dp)*self%dose_per_fraction) ) return
            if( .not. close_value(hi,real(self%windows(i)%last,dp)*self%dose_per_fraction) ) return
            message = 'Manifest half-map path is empty, too long or contains whitespace'
            if( self%paths_resolved )then
                if( self%windows(i)%even%strlen() > MAX_PATH .or. self%windows(i)%odd%strlen() > MAX_PATH ) return
            else
                if( .not. valid_path(self%windows(i)%even) .or. .not. valid_path(self%windows(i)%odd) ) return
            endif
            if( i > 1 )then
                message = 'Dose series intervals or doses are not strictly increasing'
                if( self%windows(i)%first <= self%windows(i-1)%first ) return
                if( self%windows(i)%last <= self%windows(i-1)%last ) return
                if( lo <= self%windows(i-1)%lo ) return
                if( i < n )then
                    message = 'Only the final dose interval may overlap its predecessor'
                    if( self%windows(i)%first <= self%windows(i-1)%last ) return
                endif
            endif
        enddo
        drop = overlap_predecessor(self)
        message = 'Pruning the final overlap leaves fewer than two windows'
        if( drop > 0 .and. n < 3 ) return
        previous = 0
        do i = 1, n
            if( i == drop ) cycle
            if( previous > 0 )then
                message = 'Dose series still overlaps after final-window pruning'
                if( self%windows(i)%first <= self%windows(previous)%last ) return
            endif
            previous = i
        enddo
        message = ''
        ok = .true.
    end function validate

    !> Verify distinct canonical MRC paths against the manifest and a common coordinate system.
    logical function validate_files(self, directory, message) result(ok)
        class(dose_series_manifest), intent(in) :: self
        class(string), intent(in) :: directory
        type(string), intent(out) :: message
        type(string), allocatable :: paths(:)
        integer :: i, j, coords(6), reference(6)
        ok = .false.
        allocate(paths(2*self%get_nwindows()))
        reference = 0
        do i = 1, self%get_nwindows()
            paths(2*i-1) = resolve_path(directory, self%windows(i)%even)
            paths(2*i)   = resolve_path(directory, self%windows(i)%odd)
        enddo
        do i = 1, size(paths)
            do j = 1, i-1
                message = 'A half-map file is used more than once in the dose series: '//paths(i)%to_char()
                if( paths(i) == paths(j) ) return
            enddo
            if( .not. check_header(paths(i), self%box, self%smpd, coords, message) ) return
            if( i == 1 )then
                reference = coords
            else
                message = 'Dose-series MRC origins/start indices disagree: '//paths(i)%to_char()
                if( any(coords /= reference) ) return
            endif
            ! Available producer provenance must describe an unmasked, unregularized gridding map.
            if( .not. check_support(paths(i), message) ) return
        enddo
        message = ''
        ok = .true.
    end function validate_files

    !> Read header bytes with recoverable I/O errors before using the common MRC decoder.
    logical function check_header(fname, box, smpd, coords, message) result(ok)
        class(string), intent(in) :: fname
        integer, intent(in) :: box
        real(dp), intent(in) :: smpd
        integer, intent(out) :: coords(6)
        type(string), intent(out) :: message
        type(MrcImgHead) :: header
        integer(int8) :: bytes(MRCHEADSZ)
        integer :: unit, ios, bpp
        integer(int64) :: length, available
        real(dp) :: sampling(3), angles(3)
        real :: origins(3)
        ok = .false.
        coords = 0
        message = 'Cannot read MRC header: '//fname%to_char()
        open(newunit=unit, file=fname%to_char(), status='old', action='read', access='stream', &
            &form='unformatted', iostat=ios)
        if( ios /= 0 ) return
        inquire(unit=unit, size=length, iostat=ios)
        if( ios == 0 ) read(unit, pos=1, iostat=ios) bytes
        close(unit)
        if( ios /= 0 ) return
        call header%new
        call header%decode(bytes)
        message = 'MRC header lacks MAP signature: '//fname%to_char()
        if( header%map /= 'MAP ' )then
            call header%kill
            return
        endif
        message = 'MRC dimensions disagree with dose series box: '//fname%to_char()
        if( any([header%nx,header%ny,header%nz] /= box) )then
            call header%kill
            return
        endif
        message = 'Dose series requires real-space float32/float16 MRC maps: '//fname%to_char()
        bpp = 0
        if( header%mode == MRC_MODE_FLOAT32 ) bpp = 4
        if( header%mode == MRC_MODE_FLOAT16 ) bpp = 2
        if( bpp == 0 )then
            call header%kill
            return
        endif
        message = 'MRC grid/cell does not define positive isotropic sampling: '//fname%to_char()
        if( any([header%mx,header%my,header%mz] <= 0) )then
            call header%kill
            return
        endif
        sampling = real([header%cella1,header%cella2,header%cella3],dp) / &
            &real([header%mx,header%my,header%mz],dp)
        angles = real([header%cellb1,header%cellb2,header%cellb3],dp)
        if( .not. all(ieee_is_finite(sampling)) .or. .not. all(ieee_is_finite(angles)) )then
            call header%kill
            return
        endif
        if( any(sampling <= 0._dp) .or. any(abs(sampling-smpd) > max(SMPD_ATOL,VALUE_RTOL*smpd)) )then
            call header%kill
            return
        endif
        message = 'Dose series requires orthogonal XYZ MRC axes: '//fname%to_char()
        if( any([header%mapc,header%mapr,header%maps] /= [1,2,3]) .or. any(abs(angles-90._dp) > ANGLE_ATOL) )then
            call header%kill
            return
        endif
        origins = transfer([header%originx,header%originy,header%originz],origins)
        message = 'MRC origin is nonfinite: '//fname%to_char()
        if( .not. all(ieee_is_finite(origins)) )then
            call header%kill
            return
        endif
        coords = [header%nxstart,header%nystart,header%nzstart,header%originx,header%originy,header%originz]
        message = 'MRC volume data are truncated or its extended-header size is invalid: '//fname%to_char()
        if( header%nsymbt < 0 )then
            call header%kill
            return
        endif
        available = length-int(MRCHEADSZ,int64)-int(header%nsymbt,int64)
        call header%kill
        if( available < 0_int64 ) return
        if( int(box,int64) > available/int(bpp,int64)/int(box,int64)/int(box,int64) ) return
        message = ''
        ok = .true.
    end function check_header

    !> Reject explicitly masked/regularized producer provenance while allowing external maps without a sidecar.
    logical function check_support(fname, message) result(ok)
        use simple_halfmap_diagnostics, only: support_provenance_fname
        class(string), intent(in) :: fname
        type(string), intent(out) :: message
        type(string) :: sidecar, line
        character(len=:), allocatable :: txt
        logical :: support, kind
        integer :: unit, ios
        ok = .true.
        sidecar = support_provenance_fname(fname)
        if( .not. file_exists(sidecar) ) return
        ok = .false.
        message = 'Invalid unmasked-half support provenance: '//sidecar%to_char()
        open(newunit=unit, file=sidecar%to_char(), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        support = .false.
        kind = .false.
        do
            call read_line(unit, line, ios)
            if( ios /= 0 ) exit
            txt = trim(line%to_char())
            if( txt == 'solve_support=none' ) support = .true.
            if( txt == 'solve_kind=gridding' ) kind = .true.
            if( index(txt,'solve_support=') == 1 .and. txt /= 'solve_support=none' ) exit
            if( index(txt,'solve_kind=') == 1 .and. txt /= 'solve_kind=gridding' ) exit
        enddo
        close(unit)
        if( ios /= iostat_end ) return
        ok = support .and. kind
        if( ok ) message = ''
    end function check_support

    !> Identify only the final predecessor; validate checks all other overlap cases separately.
    integer function overlap_predecessor(self) result(drop)
        class(dose_series_manifest), intent(in) :: self
        integer :: n
        drop = 0
        n = self%get_nwindows()
        if( n < 2 ) return
        if( self%windows(n)%first <= self%windows(n-1)%last ) drop = n-1
    end function overlap_predecessor

    !> Match physical metadata to a relative tolerance that respects the saved single-precision dose.
    logical function close_value(a,b) result(close)
        real(dp), intent(in) :: a,b
        close = .false.
        if( .not. ieee_is_finite(a) .or. .not. ieee_is_finite(b) ) return
        close = abs(a-b) <= VALUE_RTOL*max(1._dp,abs(a),abs(b))
    end function close_value

    !> Resolve existing paths canonically so aliases cannot reuse one file as several independent halves.
    function resolve_path(directory, path) result(resolved)
        class(string), intent(in) :: directory, path
        type(string) :: resolved
        character(len=:), allocatable :: text
        integer :: status
        text = path%to_char()
        resolved = path
        if( len(text) == 0 ) return
        if( text(1:1) /= '/' ) resolved = directory//path
        resolved = simple_abspath(resolved, status=status, check_exists=.false.)
        if( status /= 0 ) resolved = ''
    end function resolve_path

    !> Make paths beneath the manifest directory portable; preserve external absolute paths.
    function portable_path(directory, path) result(portable)
        class(string), intent(in) :: directory, path
        type(string) :: portable, absolute
        character(len=:), allocatable :: root, text
        absolute = resolve_path(directory, path)
        root = directory%to_char()
        text = absolute%to_char()
        portable = absolute
        if( index(text,root) == 1 ) portable = text(len(root)+1:)
    end function portable_path

    !> Obtain an absolute directory with a trailing separator before the consumer changes directory.
    function parent_directory(fname) result(directory)
        class(string), intent(in) :: fname
        type(string) :: directory, absolute
        character(len=:), allocatable :: text
        integer :: pos, status
        absolute = simple_abspath(fname, status=status, check_exists=.false.)
        text = absolute%to_char()
        pos = index(text,'/',back=.true.)
        directory = text(:pos)
    end function parent_directory

    !> Enforce version-1 path tokens without silently truncating or interpreting quoted whitespace.
    logical function valid_path(path) result(valid)
        class(string), intent(in) :: path
        character(len=:), allocatable :: text
        text = path%to_char()
        valid = len(text) > 0 .and. len(text) <= MAX_PATH
        if( scan(text,' '//achar(9)//achar(10)//achar(13)//achar(0)) /= 0 ) valid = .false.
    end function valid_path

    !> Read a complete logical line in chunks; exceeding the explicit format bound is an error.
    subroutine read_line(unit, line, ios)
        integer, intent(in) :: unit
        type(string), intent(out) :: line
        integer, intent(out) :: ios
        character(len=LINE_CHUNK) :: chunk
        character(len=:), allocatable :: text
        integer :: got
        text = ''
        do
            read(unit,'(A)',advance='no',size=got,iostat=ios) chunk
            if( got > 0 ) text = text//chunk(:got)
            if( len(text) > MAX_LINE )then
                ios = 1
                exit
            endif
            if( ios == iostat_eor )then
                ios = 0
                exit
            endif
            if( ios /= 0 ) exit
        enddo
        line = text
    end subroutine read_line

    !> Split all whitespace tokens, reporting excess columns instead of discarding them.
    subroutine split_tokens(text, tokens, n)
        character(len=*), intent(in) :: text
        type(string), intent(out) :: tokens(:)
        integer, intent(out) :: n
        integer :: i, first, length
        n = 0
        i = 1
        length = len_trim(text)
        do while( i <= length )
            if( text(i:i) == ' ' .or. text(i:i) == achar(9) )then
                i = i+1
                cycle
            endif
            first = i
            do while( i <= length )
                if( text(i:i) == ' ' .or. text(i:i) == achar(9) ) exit
                i = i+1
            enddo
            n = n+1
            if( n > size(tokens) ) return
            tokens(n) = text(first:i-1)
        enddo
    end subroutine split_tokens

    !> Return the index of required metadata; other comment lines are descriptive only.
    integer function header_slot(key) result(slot)
        character(len=*), intent(in) :: key
        slot = 0
        select case(key)
            case('program'); slot = 1
            case('state'); slot = 2
            case('box'); slot = 3
            case('smpd'); slot = 4
            case('kv'); slot = 5
            case('nptcls'); slot = 6
            case('neven'); slot = 7
            case('nodd'); slot = 8
            case('eo_digest'); slot = 9
            case('dose_per_fraction'); slot = 10
        end select
    end function header_slot

    !> Parse exactly one integer token, rejecting list-directed separators and repetitions.
    logical function parse_int(text,value) result(ok)
        character(len=*), intent(in) :: text
        integer, intent(out) :: value
        integer :: ios
        ok = .false.
        if( scan(text,',/*') /= 0 ) return
        read(text,*,iostat=ios) value
        ok = ios == 0
    end function parse_int

    !> Parse exactly one finite real token, rejecting list-directed separators and repetitions.
    logical function parse_real(text,value) result(ok)
        character(len=*), intent(in) :: text
        real(dp), intent(out) :: value
        integer :: ios
        ok = .false.
        if( scan(text,',/*') /= 0 ) return
        read(text,*,iostat=ios) value
        if( ios /= 0 ) return
        ok = ieee_is_finite(value)
    end function parse_real

    !> Remove only a temporary manifest owned by this attempted write.
    subroutine remove_temporary(fname)
        class(string), intent(in) :: fname
        integer :: unit, ios
        open(newunit=unit,file=fname%to_char(),status='old',iostat=ios)
        if( ios == 0 ) close(unit,status='delete',iostat=ios)
    end subroutine remove_temporary

end module simple_dose_series_manifest
