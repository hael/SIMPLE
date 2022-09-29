! SAURON: SIMPLE Attempt to a Unified Resources and Orientations Notebook
! provides support for one-line per one particle input/output
module simple_sauron
use simple_defs
use simple_defs_ori
use simple_strings, only: str2format, str2real, str2int, parsestr, split_str, compact
use simple_hash,    only: hash
use simple_chash,   only: chash
use simple_syslib,  only: sscanf
implicit none

public :: sauron_ori_parser
private

interface sauron_ori_parser
    module procedure sauron_ori_parser_1
    module procedure sauron_ori_parser_2
end interface

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
        ! On Apple M1 some combinations of XCode & GNU compilers make sscanf fail systematically,
        ! Fortran instrinsic is a bit slower but safer it seems
        read(str, *, iostat=i) rvar
        isreal = i==0
        ! isreal = .true.
        ! rvar = str2real(str)
        ! isreal = sscanf(str(1:l)//C_NULL_CHAR,"%f"//C_NULL_CHAR,rvar) /= 0
    end subroutine str2real_local

    !> This parser is optimised for real/character output types.
    !  Refer to sauron_ori_parser_old for more generic cases
    subroutine sauron_ori_parser_1( line, htab, chtab )
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
    end subroutine sauron_ori_parser_1

    subroutine sauron_ori_parser_2( line, pparms, htab, chtab )
        character(len=*), intent(inout) :: line
        real,             intent(inout) :: pparms(N_PTCL_ORIPARAMS)
        class(hash),      intent(inout) :: htab
        class(chash),     intent(inout) :: chtab
        logical               :: is_real(96)
        character(len=32)     :: keys(96)
        character(len=STDLEN) :: args(96), vals(96), arg
        real                  :: rvals(96)
        integer               :: eq_pos(96), nargs, iarg, lenstr, ipos, ipos_prev, nreals, i, ind
        call compact(line)
        lenstr = len_trim(line)
        if( lenstr == 0 ) return
        ! split entries
        nargs     = 0
        ipos_prev = 1
        do ipos = 2, lenstr - 1
            if( line(ipos:ipos) == ' ' )then
                i = scan(line(ipos_prev:ipos-1), '=')
                if( i == 0 ) cycle
                nargs         = nargs + 1
                eq_pos(nargs) = i
                args(nargs)   = line(ipos_prev:ipos - 1)
                ipos_prev     = ipos + 1
            endif
        enddo
        i = scan(line(ipos_prev:lenstr), '=')
        if( i > 0 )then
            nargs         = nargs + 1
            eq_pos(nargs) = i
            args(nargs)   = line(ipos_prev:lenstr)
        endif
        ! split key/values pairs
        do iarg=1,nargs
            arg  = trim(args(iarg))
            ipos = eq_pos(iarg)
            keys(iarg) = adjustl(trim(arg(1:ipos - 1)))
            vals(iarg) = adjustl(arg(ipos + 1:))
            call str2real_local(vals(iarg), is_real(iarg), rvals(iarg))
        end do
        nreals = 0
        do iarg = 1, nargs
            if( is_real(iarg) )then
                ind = get_oriparam_ind(trim(keys(iarg)))
                if( ind == 0 ) nreals = nreals + 1
            endif
        end do
        ! alloc
        if( nreals > 0 )      call htab%new(nreals)
        if( nreals /= nargs ) call chtab%new(nargs - nreals)
        ! fill
        do iarg = 1, nargs
            if( is_real(iarg) )then
                ind = get_oriparam_ind(trim(keys(iarg)))
                if( ind == 0 )then
                    call htab%push(keys(iarg), rvals(iarg))
                else
                    pparms(ind) = rvals(iarg)
                endif
            else
                call chtab%push(keys(iarg), vals(iarg))
            endif
        end do
    end subroutine sauron_ori_parser_2

end module simple_sauron
