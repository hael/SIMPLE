module simple_string_tester
use simple_string
use simple_defs
use simple_test_utils
implicit none
private
public :: run_all_string_tests

contains

    subroutine run_all_string_tests
        write(*,'(A)') '**** running all string tests ****'
        call test_default_and_kill()
        call test_assign_and_ctor()
        call test_append_operator()
        call test_strlen()
        call test_to_char()
        call test_numeric_conversion()
        call test_has_substr()
        call test_substr_remove()
        call test_ends_with()
        call test_substr_ind()
        call test_substr_replace()
        call test_is_blank_and_is_allocated()
        call test_equality_operators()
        call test_readline_writeline()
        call test_readline_multiple_lines()
        call test_readline_empty_line()
        call test_readline_long_line()
        call test_readline_eof_behavior()
        call test_writeline_unallocated()
        ! call report_summary()
    end subroutine run_all_string_tests

    !---------------- basic lifecycle ----------------

    subroutine test_default_and_kill()
        type(string) :: s
        write(*,'(A)') 'test_default_and_kill'
        ! default: unallocated, blank, length 0
        call assert_true(.not. s%is_allocated(), 'default string not allocated')
        call assert_true(      s%is_blank(),     'default string is_blank() true')
        ! strlen on unallocated must be 0 (per semantics)
        call assert_int(0, s%strlen(),           'strlen on unallocated = 0')
        call assert_int(0, s%strlen_trim(),      'strlen_trim on unallocated = 0')
        s = 'abc'
        call assert_true(s%is_allocated(),       'allocated after assignment')
        call assert_true(.not. s%is_blank(),     'non-empty after assignment not blank')
        call s%kill()
        call assert_true(.not. s%is_allocated(), 'kill deallocates')
        call assert_true(      s%is_blank(),     'killed string is_blank() true')
    end subroutine test_default_and_kill

    !---------------- assignment & ctor (no trimming in ctor) ----------------

    subroutine test_assign_and_ctor()
        type(string) :: s1, s2
        integer :: len1
        write(*,'(A)') 'test_assign_and_ctor'
        ! ctor must preserve whitespace in the buffer
        s1 = string('  hello  ')
        ! buffer should contain both leading and trailing blanks
        len1 = s1%strlen()
        call assert_int(len('  hello  '), len1,        'ctor preserves full length with spaces')
        ! raw() returns underlying buffer including spaces
        call assert_char('  hello  ', s1%raw(),        'raw() preserves leading & trailing spaces')
        ! substring index should reflect leading spaces
        call assert_int(3, s1%substr_ind('hello'),     'ctor preserves leading spaces')
        ! to_char() drops trailing but not leading blanks (due to trim())
        call assert_char('  hello', s1%to_char(),      'to_char() keeps leading, drops trailing blanks')
        ! assignment from string should COPY, not steal, and preserve buffer contents
        s2 = s1
        call assert_true( s2%is_allocated(),           'assigned string is allocated')
        call assert_true( s1%is_allocated(),           'source string still allocated after assignment (copy semantics)')
        call assert_int(s1%strlen(), s2%strlen(),      'assignment preserves length')
        call assert_int(3, s2%substr_ind('hello'),     'assignment preserves leading spaces in copy')
        ! assign character(*), must also preserve whitespace
        s1 = '  world  '
        call assert_int(len('  world  '), s1%strlen(), 'assignment from char preserves spaces')
        call assert_int(3, s1%substr_ind('world'),     'assignment from char preserves leading spaces')
    end subroutine test_assign_and_ctor

    !---------------- concatenation (no trimming) ----------------

    subroutine test_append_operator()
        type(string) :: s1, s2, s3
        write(*,'(A)') 'test_append_operator'
        s1 = string('foo')
        s2 = string('  bar  ')
        s3 = s1 // s2
        ! result buffer should be exactly 'foo  bar  '
        call assert_int(len('foo  bar  '), s3%strlen(), 'append preserves spaces in result length')
        call assert_int(1, s3%substr_ind('foo'),        'append preserves foo at start')
        call assert_int(4, s3%substr_ind('  bar'),      'append preserves spaces before bar')
        call assert_int(6, s3%substr_ind('bar'),        'bar starts at position 6')
        ! string // char
        s3 = s1 // '  baz  '
        call assert_int(len('foo  baz  '), s3%strlen(), 'append string//char preserves spaces')
        call assert_int(4, s3%substr_ind('  baz'),      'append string//char leading spaces before baz')
    end subroutine test_append_operator

    !---------------- length ----------------

    subroutine test_strlen()
        type(string) :: s
        write(*,'(A)') 'test_strlen'
        s = string('abc   ')
        call assert_int(len('abc   '), s%strlen(),           'strlen counts full buffer length including trailing spaces')
        call assert_int(len_trim('abc   '), s%strlen_trim(), 'strlen_trim matches len_trim of buffer text')
    end subroutine test_strlen

    !---------------- to_char overloads ----------------

    subroutine test_to_char()
        type(string) :: s
        character(:), allocatable :: sub
        write(*,'(A)') 'test_to_char'
        s = string('abcdef  ')
        ! to_char() should drop trailing spaces but keep all characters before them
        call assert_char('abcdef', s%to_char(), 'to_char() trims trailing blanks only')
        ! to_char([i,j]) returns exact slice of buffer (including spaces if any)
        s = string('  abcdef  ')
        sub = s%to_char([3,5])  ! positions 3:5 are 'abc'
        call assert_char('abc', sub,            'to_char([3,5]) substring from buffer')
        sub = s%to_char([1,2])
        call assert_char('  ', sub,             'to_char([1,2]) returns leading spaces')
    end subroutine test_to_char

    !---------------- numeric conversions ----------------

    subroutine test_numeric_conversion()
        type(string) :: s
        integer :: ival, ios
        real    :: rval
        real(dp):: dval
        integer :: ios_r, ios_d
        write(*,'(A)') 'test_numeric_conversion'
        s = string('123')
        ival = s%to_int(io_stat=ios)
        call assert_int(0, ios, 'to_int success ios=0')
        call assert_int(123, ival, 'to_int 123')
        s = string('  1.5  ')
        rval = s%to_real(io_stat=ios_r)
        dval = s%to_dble(io_stat=ios_d)
        call assert_int(0, ios_r, 'to_real success ios=0')
        call assert_int(0, ios_d, 'to_dble success ios=0')
        call assert_true(abs(rval - 1.5) < 1.0e-6, 'to_real value')
        call assert_true(abs(real(dval,kind=dp) - 1.5_dp) < 1.0e-12_dp, 'to_dble value')
        ! invalid numeric
        s = string('abc')
        ival = s%to_int(io_stat=ios)
        call assert_true(ios /= 0, 'to_int invalid sets nonzero ios')
        call assert_int(0, ival, 'to_int invalid returns 0')
        ! constructors from numeric types
        s = string(42)
        call assert_char('42', s%to_char(), 'string(int) constructor')
        s = string(1.5)
        call assert_true(abs(s%to_real() - 1.5) < 1.0e-6, 'string(real) constructor')
        s = string(1.5_dp)
        call assert_true(abs(s%to_dble() - 1.5_dp) < 1.0e-12_dp, 'string(dble) constructor')
    end subroutine test_numeric_conversion

    !---------------- substring presence ----------------

    subroutine test_has_substr()
        type(string) :: s, sub
        write(*,'(A)') 'test_has_substr'
        s   = string(' hello  world ')
        sub = string('world')
        call assert_true(s%has_substr(sub),         'has_substr(string) true')
        call assert_true(s%has_substr('world'),     'has_substr(char) true')
        call assert_true(.not. s%has_substr('xxx'), 'has_substr(char) false')
    end subroutine test_has_substr

    !---------------- substr_remove ----------------

    subroutine test_substr_remove()
        type(string) :: s, sub, res
        write(*,'(A)') 'test_substr_remove'
        s   = string('abc123abc123')
        sub = string('123')
        res = s%substr_remove(sub)
        call assert_char('abcabc', res%to_char(),   'substr_remove removes all occurrences')
        s   = string('xxxx')
        sub = string('xxxx')
        res = s%substr_remove(sub)
        call assert_true(res%is_blank(),            'substr_remove whole string -> blank/empty')
        s   = string('nochange')
        sub = string('zzz')
        res = s%substr_remove(sub)
        call assert_char('nochange', res%to_char(), 'substr_remove with no match returns same')
    end subroutine test_substr_remove

    !---------------- ends_with_substr ----------------

    subroutine test_ends_with()
        type(string) :: s, sub
         write(*,'(A)') 'test_ends_with'
        s   = string('foo_bar  ')
        sub = string('bar  ')
        call assert_true(s%ends_with_substr(sub),         'ends_with_substr(string) true including spaces')
        call assert_true(.not. s%ends_with_substr('bar'), 'ends_with_substr(char) must match exact tail')
    end subroutine test_ends_with

    !---------------- substr_ind (index) ----------------

    subroutine test_substr_ind()
        type(string) :: s, sub
        integer :: pos
        write(*,'(A)') 'test_substr_ind'
        s   = string('  abcabc  ')
        sub = string('bc')
        pos = s%substr_ind('bc')
        call assert_int(4, pos, 'substr_ind(char) first occurrence (with leading spaces)')
        pos = s%substr_ind('bc', back=.true.)
        call assert_int(7, pos, 'substr_ind(char,back) last occurrence')
        pos = s%substr_ind(sub)
        call assert_int(4, pos, 'substr_ind(string) first occurrence')
        pos = s%substr_ind(sub, back=.true.)
        call assert_int(7, pos, 'substr_ind(string,back) last occurrence')
        pos = s%substr_ind('zzz')
        call assert_int(0, pos, 'substr_ind not found -> 0')
    end subroutine test_substr_ind

    !---------------- substr_replace ----------------

    subroutine test_substr_replace()
        type(string) :: s
        write(*,'(A)') 'test_substr_replace'
        s = string('1-2-3-2-1')
        call s%substr_replace('2', 'X')
        call assert_char('1-X-3-X-1', s%to_char(), 'substr_replace all occurrences')
        s = string('1-2-3-2-1')
        call s%substr_replace('2', 'Y', one=.true.)
        call assert_char('1-Y-3-2-1', s%to_char(), 'substr_replace one=.true.')
        s = string('1-2-3-2-1')
        call s%substr_replace('2', 'Z', one=.true., back=.true.)
        call assert_char('1-2-3-Z-1', s%to_char(), 'substr_replace one=.true., back=.true.')
        s = string('aaaa')
        call s%substr_replace('aa', 'b')
        call assert_char('bb', s%to_char(),        'substr_replace overlapping semantics')
    end subroutine test_substr_replace

    !---------------- blank / allocated ----------------

    subroutine test_is_blank_and_is_allocated()
        type(string) :: s
        write(*,'(A)') 'test_is_blank_and_is_allocated'
        ! default: unallocated, but considered blank
        call assert_true(.not. s%is_allocated(), 'uninitialized not allocated')
        call assert_true(      s%is_blank(),     'uninitialized is_blank() true')
        s = string('   ')
        call assert_true(s%is_allocated(),       'spaces only is allocated')
        call assert_true(s%is_blank(),           'spaces only is_blank() true (len_trim==0)')
        s = string('  a ')
        call assert_true(.not. s%is_blank(),     'nonblank not is_blank()')
    end subroutine test_is_blank_and_is_allocated

    !---------------- equality / inequality ----------------

    subroutine test_equality_operators()
        type(string) :: s1, s2
        write(*,'(A)') 'test_equality_operators'
        s1 = string('abc')
        s2 = string('abc')
        call assert_true(       s1 .eq. s2,    'string .eq. same contents')
        call assert_true(.not. (s1 .ne. s2),   'string .ne. same is false')
        s2 = string('abc  ')
        call assert_true(       s1 .eq. s2,    'string .eq. same contents, despite trailing spaces')
        call assert_true(.not. (s1 .ne. s2),   'string .ne. same is false, despite trailing spaces')
        s2 = string('  abc')
        call assert_true(       s1 .eq. s2,    'string .eq. same contents, despite leading spaces')
        call assert_true(.not. (s1 .ne. s2),   'string .ne. same is false, despite leading spaces')
        s2 = string('  abc  ')
        call assert_true(       s1 .eq. s2,    'string .eq. same contents, despite leading & trailing spaces')
        call assert_true(.not. (s1 .ne. s2),   'string .ne. same is false, despite leading & trailing spaces')
        ! unallocated vs blank must not be equal
        call s1%kill()
        s2 = string('')
        call assert_true(.not. s1 .eq. s2,     'unallocated == explicit empty string is false')
        call assert_true(      s1 .ne. s2,     'unallocated /= explicit empty string is true')
        ! string vs char: equality uses exact buffer contents
        s1 = string('foo')
        call assert_true( s1 .eq. 'foo',       'string .eq. char exact match')
        call assert_true( s1  ==  'foo',       'string == char exact match')
        s2 = string('foo ')
        call assert_true( .not. s2 .ne. 'foo', 'string with extra space == char without space')
        call assert_true( .not. s2 /= 'foo',   'string with extra space == char without space')
    end subroutine test_equality_operators

    !---------------- I/O ----------------

    subroutine test_readline_writeline()
        type(string) :: s_in, s_out
        integer :: unit, ios
        write(*,'(A)') 'test_readline_writeline'
        ! create scratch file
        open(newunit=unit, status='scratch', action='readwrite', iostat=ios)
        call assert_int(0, ios, 'open scratch file')
        ! writeline: should write buffer exactly as stored
        s_out = string('  hello world  ')
        call s_out%writeline(unit, ios)
        call assert_int(0, ios, 'writeline ios=0')
        rewind(unit)
        call s_in%readline(unit, ios)
        call assert_int(0, ios, 'readline ios=0')
        ! After readline, buffer should preserve leading and trailing spaces
        call assert_int(len('hello world'), s_in%strlen(), 'readline removes leading and trailing spaces')
        call assert_int(4, s_in%substr_ind('lo'),          'readline removes leading and trailing spaces')
        close(unit)
    end subroutine test_readline_writeline

    subroutine test_readline_multiple_lines()
        type(string) :: s_in, s_out
        integer :: unit, ios, i
        character(:), allocatable :: line
        write(*,'(A)') 'test_readline_multiple_lines'
        open(newunit=unit, status='scratch', action='readwrite', iostat=ios)
        call assert_int(0, ios, 'open scratch file (multi)')
        ! write 3 lines with various spaces
        do i = 1, 3
            select case(i)
            case(1)
                s_out = string('  line1  ')
            case(2)
                s_out = string('   line2   ')
            case(3)
                s_out = string('line3')
            end select
            call s_out%writeline(unit, ios)
            call assert_int(0, ios, 'writeline ios=0 (multi)')
        end do
        rewind(unit)
        ! read #1
        call s_in%readline(unit, ios)
        call assert_int(0, ios, 'readline #1 ios=0')
        line = s_in%to_char()
        call assert_char('line1', line, 'readline #1 trim(adjustl)')
        ! read #2
        call s_in%readline(unit, ios)
        call assert_int(0, ios, 'readline #2 ios=0')
        line = s_in%to_char()
        call assert_char('line2', line, 'readline #2 trim(adjustl)')
        ! read #3
        call s_in%readline(unit, ios)
        call assert_int(0, ios, 'readline #3 ios=0')
        line = s_in%to_char()
        call assert_char('line3', line, 'readline #3 unchanged')
        close(unit)
    end subroutine test_readline_multiple_lines

    subroutine test_readline_empty_line()
        type(string) :: s_in, s_out
        integer :: unit, ios
        character(:), allocatable :: line
        write(*,'(A)') 'test_readline_empty_line'
        open(newunit=unit, status='scratch', action='readwrite', iostat=ios)
        call assert_int(0, ios, 'open scratch file (empty line)')
        ! write: nonempty, empty, nonempty
        s_out = string('first')
        call s_out%writeline(unit, ios)
        call assert_int(0, ios, 'writeline first')
        s_out = string('')   ! empty string -> writes empty line
        call s_out%writeline(unit, ios)
        call assert_int(0, ios, 'writeline empty')
        s_out = string('third')
        call s_out%writeline(unit, ios)
        call assert_int(0, ios, 'writeline third')
        rewind(unit)
        ! read first
        call s_in%readline(unit, ios)
        call assert_int(0, ios, 'readline first ios=0')
        line = s_in%to_char()
        call assert_char('first', line, 'first line text')
        ! read empty line -> buffer_static all blanks -> trim(adjustl) -> ''
        call s_in%readline(unit, ios)
        call assert_int(0, ios, 'readline empty ios=0')
        call assert_true(s_in%is_allocated(), 'empty line allocates buffer')
        call assert_true(s_in%is_blank(), 'empty line gives blank string')
        ! read third
        call s_in%readline(unit, ios)
        call assert_int(0, ios, 'readline third ios=0')
        line = s_in%to_char()
        call assert_char('third', line, 'third line text')
        close(unit)
    end subroutine test_readline_empty_line

    subroutine test_readline_long_line()
        type(string) :: s_in, s_out
        integer :: unit, ios, i, long_len
        character(:), allocatable :: expected, actual
        write(*,'(A)') 'test_readline_long_line'
        ! pick a length safely below XLONGSTRLEN
        long_len = min(XLONGSTRLEN - 1, 4096)
        allocate(character(len=long_len) :: expected)
        do i = 1, long_len
            expected(i:i) = 'a'
        end do
        open(newunit=unit, status='scratch', action='readwrite', iostat=ios)
        call assert_int(0, ios, 'open scratch file (long line)')
        s_out = expected
        call s_out%writeline(unit, ios)
        call assert_int(0, ios, 'writeline long line')
        rewind(unit)
        call s_in%readline(unit, ios)
        call assert_int(0, ios, 'readline long line ios=0')
        actual = s_in%to_char()
        call assert_int(len_trim(expected), len_trim(actual), 'long line length preserved')
        call assert_char(expected, actual, 'long line content preserved')
        close(unit)
        deallocate(expected)
    end subroutine test_readline_long_line

    subroutine test_readline_eof_behavior()
        type(string) :: s_in, s_out
        integer :: unit, ios
        character(:), allocatable :: line
        write(*,'(A)') 'test_readline_eof_behavior'
        open(newunit=unit, status='scratch', action='readwrite', iostat=ios)
        call assert_int(0, ios, 'open scratch file (eof)')
        s_out = string('only')
        call s_out%writeline(unit, ios)
        call assert_int(0, ios, 'writeline only')
        rewind(unit)
        ! first read: OK
        call s_in%readline(unit, ios)
        call assert_int(0, ios, 'readline only ios=0')
        line = s_in%to_char()
        call assert_char('only', line, 'only line content')
        ! second read: EOF -> ios /= 0, and kill leaves it unallocated
        call s_in%readline(unit, ios)
        call assert_true(ios /= 0, 'readline at EOF sets nonzero ios')
        call assert_true(.not. s_in%is_allocated(), 'after EOF, buffer is unallocated')
        call assert_true(s_in%is_blank(), 'after EOF, string is_blank() true')
        close(unit)
    end subroutine test_readline_eof_behavior

    subroutine test_writeline_unallocated()
        type(string) :: s
        integer :: unit, ios
        write(*,'(A)') 'test_writeline_unallocated'
        ! s is default: unallocated
        call assert_true(.not. s%is_allocated(), 'default string not allocated')
        open(newunit=unit, status='scratch', action='readwrite', iostat=ios)
        call assert_int(0, ios, 'open scratch file (writeline unallocated)')
        call s%writeline(unit, ios)
        call assert_int(-1, ios, 'writeline on unallocated sets ios=-1')
        close(unit)
    end subroutine test_writeline_unallocated

end module simple_string_tester
