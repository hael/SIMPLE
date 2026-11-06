!@descr: the measured dose curve dose_filter3D writes: per-window contrast and B-factor knots in dose, read per movie fraction
! Model and file: doc/implementation_notes/planned/dose_filter3D.md, section 12
module simple_dose_curve
use simple_defs
use simple_error,        only: simple_exception
use simple_fileio,       only: fopen, fclose, fileiochk
use simple_string,       only: string
use simple_string_utils, only: int2str, realdp2str
use, intrinsic :: iso_fortran_env, only: iostat_end
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
#include "simple_local_flags.inc"

public :: dose_curve, scaled_dose
private

character(len=*), parameter :: VERSION_LINE = '# SIMPLE dose curve 1'
integer,  parameter :: NHEADER   = 7         !< program kv dose_per_fraction band line contrast nknots
real(dp), parameter :: KV_ATOL   = 1.e-3_dp
real(dp), parameter :: DOSE_RTOL = 1.e-5_dp  !< the saved dose per fraction is single precision

type :: dose_curve
    private
    integer  :: nknots = 0
    real(dp) :: kv = 0._dp, dose_per_fraction = 0._dp, hp = 0._dp, lp = 0._dp
    real(dp) :: b0 = 0._dp, beta = 0._dp, lnc0 = 0._dp, gamma = 0._dp
    integer,  allocatable :: ffrom(:), fto(:)  !< (nknots) movie fractions of each window
    real(dp), allocatable :: dose(:)           !< (nknots) scaled dose at the window centre
    real(dp), allocatable :: c(:), b(:)        !< (nknots) contrast and B-factor relative to window 1
    logical,  allocatable :: used(:)           !< (nknots) the window's attenuation was fitted
    logical  :: l_exists = .false.
  contains
    procedure :: new
    procedure :: set_knot
    procedure :: write
    procedure :: read
    procedure :: check_movie
    procedure :: fraction_weights
    procedure :: kill
    procedure, private :: validate
    procedure, private :: evaluate
end type dose_curve

contains

    !> the dose coordinate of calc_dose_weights: physical dose divided by 0.8 at 200 kV and by 0.64 at 100 kV
    pure real(dp) function scaled_dose( dose, kv ) result( x )
        real(dp), intent(in) :: dose, kv
        x = dose
        if( abs(kv - 200._dp) < 1._dp ) x = dose / 0.8_dp
        if( abs(kv - 100._dp) < 1._dp ) x = dose / 0.64_dp
    end function scaled_dose

    !> nknots windows recorded at kv with a physical dose per fraction; band in A; in the scaled dose D the
    !! line is B = b0 + beta D and the contrast ln c = lnc0 - gamma D
    subroutine new( self, nknots, kv, dose_per_fraction, hp, lp, b0, beta, lnc0, gamma )
        class(dose_curve), intent(inout) :: self
        integer,           intent(in)    :: nknots
        real(dp),          intent(in)    :: kv, dose_per_fraction, hp, lp, b0, beta, lnc0, gamma
        call self%kill
        if( nknots < 2 ) THROW_HARD('A dose curve needs at least two windows')
        self%nknots            = nknots
        self%kv                = kv
        self%dose_per_fraction = dose_per_fraction
        self%hp                = hp
        self%lp                = lp
        self%b0                = b0
        self%beta              = beta
        self%lnc0              = lnc0
        self%gamma             = gamma
        allocate(self%ffrom(nknots), self%fto(nknots), source=0)
        allocate(self%dose(nknots), self%c(nknots), self%b(nknots), source=0._dp)
        allocate(self%used(nknots), source=.false.)
        self%l_exists = .true.
    end subroutine new

    !> one window: its fractions, scaled centre dose, contrast and B-factor, and whether its fit is used
    subroutine set_knot( self, i, ffrom, fto, dose, c, b, used )
        class(dose_curve), intent(inout) :: self
        integer,           intent(in)    :: i, ffrom, fto
        real(dp),          intent(in)    :: dose, c, b
        logical,           intent(in)    :: used
        if( .not. self%l_exists )        THROW_HARD('Dose curve not constructed')
        if( i < 1 .or. i > self%nknots ) THROW_HARD('Dose curve window index out of range')
        self%ffrom(i) = ffrom
        self%fto(i)   = fto
        self%dose(i)  = dose
        self%c(i)     = c
        self%b(i)     = b
        self%used(i)  = used
    end subroutine set_knot

    subroutine write( self, fname )
        class(dose_curve), intent(in) :: self
        class(string),     intent(in) :: fname
        type(string) :: message
        integer      :: funit, i, io_stat
        if( .not. self%validate(message) ) THROW_HARD(message%to_char())
        call fopen(funit, file=fname, status='replace', action='write', form='formatted')
        write(funit,'(A)',iostat=io_stat) VERSION_LINE
        if( io_stat == 0 ) write(funit,'(A)',iostat=io_stat) '# program dose_filter3D'
        if( io_stat == 0 ) write(funit,'(A)',iostat=io_stat) '# units: dose in e/A^2 scaled as calc_dose_weights; B in A^2; band in A'
        if( io_stat == 0 ) write(funit,'(A,ES24.16E3)',iostat=io_stat) '# kv ', self%kv
        if( io_stat == 0 ) write(funit,'(A,ES24.16E3)',iostat=io_stat) '# dose_per_fraction ', self%dose_per_fraction
        if( io_stat == 0 ) write(funit,'(A,2(1X,ES24.16E3))',iostat=io_stat) '# band', self%hp, self%lp
        if( io_stat == 0 ) write(funit,'(A,2(1X,ES24.16E3))',iostat=io_stat) '# line', self%b0, self%beta
        if( io_stat == 0 ) write(funit,'(A,2(1X,ES24.16E3))',iostat=io_stat) '# contrast', self%lnc0, self%gamma
        if( io_stat == 0 ) write(funit,'(A,I0)',iostat=io_stat) '# nknots ', self%nknots
        if( io_stat == 0 ) write(funit,'(A)',iostat=io_stat) '# columns window from to dose c b used'
        do i = 1, self%nknots
            if( io_stat /= 0 ) exit
            write(funit,'(3(I0,1X),3(ES24.16E3,1X),I0)',iostat=io_stat) i, self%ffrom(i), self%fto(i), &
                &self%dose(i), self%c(i), self%b(i), merge(1, 0, self%used(i))
        enddo
        call fileiochk('Writing the dose curve '//fname%to_char(), io_stat)
        call fclose(funit)
    end subroutine write

    !> parse and validate a curve file; on failure the curve is left empty and message says why
    logical function read( self, fname, message ) result( ok )
        class(dose_curve), intent(inout) :: self
        class(string),     intent(in)    :: fname
        type(string),      intent(out)   :: message
        character(len=LONGSTRLEN)     :: line
        character(len=:), allocatable :: txt
        character(len=32)             :: key, word
        integer,  allocatable :: ffrom(:), fto(:), iused(:)
        real(dp), allocatable :: dose(:), c(:), b(:)
        real(dp) :: kv, dpf, band(2), line_par(2), contrast(2)
        integer  :: unit, ios, nk, nrows, irow, i, slot
        logical  :: seen(NHEADER)
        call self%kill
        ok = .false.
        open(newunit=unit, file=fname%to_char(), status='old', action='read', iostat=ios)
        if( ios /= 0 )then
            message = 'Cannot open the dose curve '//fname%to_char()
            return
        endif
        message = ''
        seen    = .false.
        nk      = 0
        nrows   = 0
        read(unit,'(A)',iostat=ios) line
        if( ios /= 0 )then
            message = 'Empty dose curve '//fname%to_char()
        else if( trim(adjustl(line)) /= VERSION_LINE )then
            message = 'Not a version 1 dose curve: '//fname%to_char()
        endif
        do while( message%strlen_trim() == 0 )
            read(unit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            txt = trim(adjustl(line))
            if( len(txt) == 0 ) cycle
            if( txt(1:1) == '#' )then
                read(txt(2:),*,iostat=ios) key
                if( ios /= 0 ) cycle
                slot = header_slot(key)
                if( slot == 0 ) cycle
                if( seen(slot) )then
                    message = 'Duplicate dose curve header: '//trim(key)
                    exit
                endif
                seen(slot) = .true.
                if( scan(txt(2:), ',/*') /= 0 .or. ntokens(txt(2:)) /= header_ntokens(slot) )then
                    message = 'Malformed dose curve header: '//trim(key)
                    exit
                endif
                select case(slot)
                    case(1)
                        read(txt(2:),*,iostat=ios) key, word
                        if( ios == 0 .and. trim(word) /= 'dose_filter3D' ) ios = 1
                    case(2)
                        read(txt(2:),*,iostat=ios) key, kv
                    case(3)
                        read(txt(2:),*,iostat=ios) key, dpf
                    case(4)
                        read(txt(2:),*,iostat=ios) key, band
                    case(5)
                        read(txt(2:),*,iostat=ios) key, line_par
                    case(6)
                        read(txt(2:),*,iostat=ios) key, contrast
                    case(7)
                        read(txt(2:),*,iostat=ios) key, nk
                        if( ios == 0 .and. nk < 2 ) ios = 1
                        if( ios == 0 ) allocate(ffrom(nk), fto(nk), iused(nk), dose(nk), c(nk), b(nk))
                end select
                if( ios /= 0 )then
                    message = 'Invalid dose curve header: '//trim(key)
                    exit
                endif
            else
                if( .not. seen(7) )then
                    message = 'Dose curve window rows precede nknots'
                    exit
                endif
                if( nrows == nk )then
                    message = 'Dose curve has more windows than nknots'
                    exit
                endif
                if( scan(txt, ',/*') /= 0 .or. ntokens(txt) /= 7 )then
                    message = 'Dose curve window rows need seven columns'
                    exit
                endif
                nrows = nrows + 1
                read(txt,*,iostat=ios) irow, ffrom(nrows), fto(nrows), dose(nrows), c(nrows), b(nrows), iused(nrows)
                if( ios /= 0 .or. irow /= nrows .or. iused(nrows) < 0 .or. iused(nrows) > 1 )then
                    message = 'Invalid dose curve window row '//int2str(nrows)
                    exit
                endif
            endif
        enddo
        if( message%strlen_trim() == 0 .and. ios /= iostat_end ) message = 'Unreadable dose curve line'
        close(unit)
        if( message%strlen_trim() == 0 )then
            if( .not. all(seen) )then
                message = 'Dose curve is missing required headers'
            else if( nrows /= nk )then
                message = 'Dose curve has fewer windows than nknots'
            endif
        endif
        if( message%strlen_trim() > 0 ) return
        call self%new(nk, kv, dpf, band(1), band(2), line_par(1), line_par(2), contrast(1), contrast(2))
        do i = 1, nk
            call self%set_knot(i, ffrom(i), fto(i), dose(i), c(i), b(i), iused(i) == 1)
        enddo
        ok = self%validate(message)
        if( .not. ok ) call self%kill

      contains

        integer function header_slot( k ) result( s )
            character(len=*), intent(in) :: k
            select case(trim(k))
                case('program');           s = 1
                case('kv');                s = 2
                case('dose_per_fraction'); s = 3
                case('band');              s = 4
                case('line');              s = 5
                case('contrast');          s = 6
                case('nknots');            s = 7
                case default;              s = 0
            end select
        end function header_slot

        !> the key and its values
        integer function header_ntokens( s ) result( n )
            integer, intent(in) :: s
            select case(s)
                case(4,5,6)
                    n = 3
                case default
                    n = 2
            end select
        end function header_ntokens

    end function read

    !> the curve holds for a movie at its voltage and dose rate, the single rate of the producer
    logical function check_movie( self, kv, dose_per_frame, message ) result( ok )
        class(dose_curve), intent(in)  :: self
        real(dp),          intent(in)  :: kv, dose_per_frame
        type(string),      intent(out) :: message
        ok = .false.
        if( .not. self%l_exists )then
            message = 'Dose curve not read'
        else if( .not. ieee_is_finite(kv) )then
            message = 'Movie voltage must be finite'
        else if( .not. ieee_is_finite(dose_per_frame) .or. dose_per_frame <= 0._dp )then
            message = 'Movie dose per fraction must be finite and positive'
        else if( abs(kv - self%kv) > KV_ATOL )then
            message = 'Movie voltage '//trim(realdp2str(kv))//' kV differs from the dose curve '// &
                &trim(realdp2str(self%kv))//' kV'
        else if( abs(dose_per_frame - self%dose_per_fraction) > DOSE_RTOL * max(dose_per_frame, self%dose_per_fraction) )then
            message = 'Movie dose per fraction '//trim(realdp2str(dose_per_frame))//' differs from the dose curve '// &
                &trim(realdp2str(self%dose_per_fraction))
        else
            message = ''
            ok      = .true.
        endif
    end function check_movie

    !> contrast and B-factor of each movie fraction, at the fraction's centre dose
    subroutine fraction_weights( self, nframes, dose_per_frame, kv, cfrac, bfrac )
        class(dose_curve),     intent(in)    :: self
        integer,               intent(in)    :: nframes
        real(dp),              intent(in)    :: dose_per_frame, kv
        real(dp), allocatable, intent(inout) :: cfrac(:), bfrac(:)
        type(string) :: message
        integer      :: f
        if( .not. self%validate(message) )                        THROW_HARD(message%to_char())
        if( .not. self%check_movie(kv, dose_per_frame, message) ) THROW_HARD(message%to_char())
        if( nframes < 1 ) THROW_HARD('A movie needs at least one fraction')
        if( allocated(cfrac) ) deallocate(cfrac)
        if( allocated(bfrac) ) deallocate(bfrac)
        allocate(cfrac(nframes), bfrac(nframes))
        do f = 1, nframes
            call self%evaluate(scaled_dose((real(f,dp) - 0.5_dp) * dose_per_frame, kv), cfrac(f), bfrac(f))
        enddo
    end subroutine fraction_weights

    !> contrast and B-factor at scaled dose x: linear in ln c and B between fitted windows, the first segment
    !! continued before the first window, and after the last the line's slopes, clamped so that an
    !! unmeasured fraction never gains on the last measured one
    pure subroutine evaluate( self, x, c, b )
        class(dose_curve), intent(in)  :: self
        real(dp),          intent(in)  :: x
        real(dp),          intent(out) :: c, b
        integer, allocatable :: k(:)
        real(dp) :: t
        integer  :: i, j, m
        k = pack([(i, i = 1, self%nknots)], self%used)
        m = size(k)
        if( x >= self%dose(k(m)) )then
            c = self%c(k(m)) * exp(-max(0._dp, self%gamma) * (x - self%dose(k(m))))
            b = self%b(k(m)) + max(0._dp, self%beta) * (x - self%dose(k(m)))
            return
        endif
        j = 1
        do i = 1, m - 1
            if( x >= self%dose(k(i)) ) j = i
        enddo
        t = (x - self%dose(k(j))) / (self%dose(k(j+1)) - self%dose(k(j)))
        c = exp((1._dp - t) * log(self%c(k(j))) + t * log(self%c(k(j+1))))
        b = (1._dp - t) * self%b(k(j)) + t * self%b(k(j+1))
    end subroutine evaluate

    !> the invariants the extraction relies on
    logical function validate( self, message ) result( ok )
        class(dose_curve), intent(in)  :: self
        type(string),      intent(out) :: message
        integer :: i
        ok      = .false.
        message = ''
        if( .not. self%l_exists )then
            message = 'Dose curve not constructed'
        else if( .not. ieee_is_finite(self%kv) )then
            message = 'Dose curve voltage must be finite'
        else if( minval(abs(self%kv - [100._dp, 200._dp, 300._dp])) > KV_ATOL )then
            message = 'Dose curve voltage must be 100, 200 or 300 kV'
        else if( .not. ieee_is_finite(self%dose_per_fraction) .or. self%dose_per_fraction <= 0._dp )then
            message = 'Dose curve dose per fraction must be positive'
        else if( .not. all(ieee_is_finite([self%hp, self%lp, self%b0, self%beta, self%lnc0, self%gamma])) )then
            message = 'Dose curve band, line and contrast must be finite'
        else if( .not. all(ieee_is_finite(self%dose)) )then
            message = 'Dose curve window doses must be finite'
        else if( any(self%ffrom < 1) .or. any(self%fto < self%ffrom) )then
            message = 'Dose curve fraction ranges are invalid'
        else if( .not. self%used(1) .or. abs(self%c(1) - 1._dp) > 1.e-9_dp .or. abs(self%b(1)) > 1.e-9_dp )then
            message = 'Dose curve window 1 must be the gauge: fitted, c = 1 and B = 0'
        else if( count(self%used) < 2 )then
            message = 'Dose curve needs at least two fitted windows'
        endif
        if( message%strlen_trim() > 0 ) return
        do i = 2, self%nknots
            if( self%dose(i) <= self%dose(i-1) )then
                message = 'Dose curve windows must increase in dose'
                return
            endif
        enddo
        do i = 1, self%nknots
            if( .not. self%used(i) ) cycle
            if( .not. ieee_is_finite(self%c(i)) .or. .not. ieee_is_finite(self%b(i)) .or. self%c(i) <= 0._dp )then
                message = 'Dose curve window '//int2str(i)//' has no valid contrast and B-factor'
                return
            endif
        enddo
        ok = .true.
    end function validate

    !> number of blank-separated tokens of a line
    pure integer function ntokens( text ) result( n )
        character(len=*), intent(in) :: text
        logical :: inside
        integer :: i
        n      = 0
        inside = .false.
        do i = 1, len_trim(text)
            if( text(i:i) == ' ' .or. text(i:i) == achar(9) )then
                inside = .false.
            else if( .not. inside )then
                inside = .true.
                n      = n + 1
            endif
        enddo
    end function ntokens

    subroutine kill( self )
        class(dose_curve), intent(inout) :: self
        if( .not. self%l_exists ) return
        if( allocated(self%ffrom) ) deallocate(self%ffrom)
        if( allocated(self%fto)   ) deallocate(self%fto)
        if( allocated(self%dose)  ) deallocate(self%dose)
        if( allocated(self%c)     ) deallocate(self%c)
        if( allocated(self%b)     ) deallocate(self%b)
        if( allocated(self%used)  ) deallocate(self%used)
        self%nknots = 0
        self%kv = 0._dp; self%dose_per_fraction = 0._dp; self%hp = 0._dp; self%lp = 0._dp
        self%b0 = 0._dp; self%beta = 0._dp; self%lnc0 = 0._dp; self%gamma = 0._dp
        self%l_exists = .false.
    end subroutine kill

end module simple_dose_curve
