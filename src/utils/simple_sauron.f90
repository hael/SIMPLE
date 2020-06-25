! SAURON: SIMPLE Attempt to a Unified Resources and Orientations Notebook
! provides support for one-line per one particle input/output
module simple_sauron
use simple_defs
use simple_strings, only: str2format, str2real, str2int, parsestr, split_str, compact
use simple_error,   only: allocchk
use simple_hash,    only: hash
use simple_chash,   only: chash
use simple_syslib,  only: sscanf
implicit none

contains

    subroutine sauron_line_parser( line, htab, chtab )
        character(len=*), intent(inout)    :: line
        class(hash),      intent(inout)    :: htab
        class(chash),     intent(inout)    :: chtab
        character(len=32),     allocatable :: keys(:)
        character(len=STDLEN), allocatable :: vals(:)
        character(len=:),      allocatable :: line_trimmed, form
        character(len=STDLEN) :: args(128), args_pair(2)
        real    :: rval
        integer :: nargs, iarg, ival
        allocate(line_trimmed, source=trim(line))
        call parsestr(line_trimmed,' ', args, nargs)
        allocate(keys(nargs), vals(nargs), stat=alloc_stat)
        if(alloc_stat /= 0)call allocchk("simple_sauron::sauron_line_parser ",alloc_stat)
        do iarg=1,nargs
            args_pair(2) = args(iarg)
            call split_str(args_pair(2), '=', args_pair(1))
            keys(iarg) = trim(args_pair(1))
            vals(iarg) = args_pair(2)
            call str2format(vals(iarg),form,rval,ival)
            select case(form)
                case( 'real' )
                    call htab%set(keys(iarg), rval)
                case( 'file', 'dir', 'char' )
                    call chtab%set(trim(keys(iarg)), trim(vals(iarg)))
                case( 'int'  )
                    call htab%set(keys(iarg), real(ival))
            end select
        end do
        deallocate(keys, vals, line_trimmed)
    end subroutine sauron_line_parser

    !> This parser is optimised for real/character output types.
    !  Refer to sauron_ori_parser_old for more generic cases
    subroutine sauron_ori_parser( line, htab, chtab )
        character(len=*), intent(inout) :: line
        class(hash),      intent(inout) :: htab
        class(chash),     intent(inout) :: chtab
        logical               :: is_real(96)
        character(len=32)     :: keys(96)
        character(len=STDLEN) :: args(96), vals(96), arg
        real                  :: rvals(96)
        integer               :: eq_pos(96), nargs, iarg, lenstr, ipos, ipos_prev, nreals, i
        call compact(line)
        lenstr = len_trim(line)
        if(lenstr == 0) return
        ! split entries
        nargs = 0
        ipos_prev = 1
        do ipos=2,lenstr-1
            if( line(ipos:ipos) == ' ' )then
                i = scan(line(ipos_prev:ipos-1),'=')
                if( i == 0 ) cycle
                nargs         = nargs+1
                eq_pos(nargs) = i
                args(nargs)   = line(ipos_prev:ipos-1)
                ipos_prev     = ipos+1
            endif
        enddo
        i = scan(line(ipos_prev:lenstr),'=')
        if( i > 0 )then
            nargs         = nargs+1
            eq_pos(nargs) = i
            args(nargs)   = line(ipos_prev:lenstr)
        endif
        ! split key/values pairs
        do iarg=1,nargs
            arg  = trim(args(iarg))
            ipos = eq_pos(iarg)
            keys(iarg) = adjustl(trim(arg(1:ipos-1)))
            vals(iarg) = adjustl(arg(ipos+1:))
            call str2real_local(vals(iarg), is_real(iarg), rvals(iarg))
        end do
        nreals = count(is_real(:nargs))
        ! alloc
        if(nreals > 0) call htab%new(nreals)
        if(nreals /= nargs ) call chtab%new(nargs-nreals)
        ! fill
        do iarg=1,nargs
            if(is_real(iarg))then
                call htab%push(keys(iarg), rvals(iarg))
            else
                call chtab%push(keys(iarg),vals(iarg))
            endif
        end do
    contains

        !> whether a string represents a real (and its value) or characters
        subroutine str2real_local( str, isreal, rvar)
            use simple_strings, only: char_is_a_letter
            character(len=*), intent(in)  :: str
            logical         , intent(out) :: isreal
            real,             intent(out) :: rvar
            integer :: i, l
            isreal = .false.
            i = index(str, '.')
            l = len_trim(str)
            ! file name
            if( i /= 0 .and. i < l )then
                if( char_is_a_letter(str(i+1:i+1)) ) return
            endif
            ! directory
            if( index(str, PATH_SEPARATOR) /= 0 ) return
            ! real
            ! read(str,*,iostat=iostat) rvar
            ! isreal = (iostat == 0)
            isreal = sscanf(str(1:l)//C_NULL_CHAR,"%f"//C_NULL_CHAR,rvar) /= 0
        end subroutine str2real_local

    end subroutine sauron_ori_parser

    subroutine sauron_ori_parser_old( line, htab, chtab )
        character(len=*), intent(inout) :: line
        class(hash),      intent(inout) :: htab
        class(chash),     intent(inout) :: chtab
        character(len=:), allocatable :: form
        character(len=32)     :: keys(128)
        character(len=STDLEN) :: args(128), vals(128), val
        real                  :: rvals(128), rval
        logical               :: real_mask(128)
        integer               :: nargs, iarg, ival, io_stat, lenstr, ipos, ipos_prev, nreals
        call compact(line)
        lenstr = len_trim(line)
        if(lenstr == 0) return
        ! split entries
        nargs = 0
        ipos_prev = 1
        do ipos=2,lenstr-1
            if( line(ipos:ipos) == ' ' )then
                nargs       = nargs+1
                args(nargs) = line(ipos_prev:ipos-1)
                ipos_prev   = ipos+1
            endif
        enddo
        nargs       = nargs+1
        args(nargs) = line(ipos_prev:lenstr)
        ! split key/values pairs
        real_mask = .true.
        nreals    = 0
        do iarg=1,nargs
            val  = trim(args(iarg))
            ipos = scan(val,'=')
            if(ipos==0)cycle
            keys(iarg) = adjustl(trim(val(1:ipos-1)))
            val  = adjustl(val(ipos+1:))
            call str2format(val,form,rval,ival)
            select case(form)
                case('real')
                    rvals(iarg) = rval
                    nreals      = nreals + 1
                case('file', 'dir', 'char')
                    real_mask(iarg) = .false.
                    vals(iarg)      = val
                case('int')
                    rvals(iarg) = real(ival)
                    nreals      = nreals + 1
            end select
        end do
        ! alloc
        if(nreals > 0) call htab%new(nreals)
        if(nreals /= nargs ) call chtab%new(nargs-nreals)
        ! fill
        do iarg=1,nargs
            if(real_mask(iarg))then
                call htab%push(keys(iarg), rvals(iarg))
            else
                call chtab%push(keys(iarg),vals(iarg))
            endif
        end do
    end subroutine sauron_ori_parser_old

end module simple_sauron
