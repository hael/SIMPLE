!@descr: byte identity of two test run directories over declared outputs, with the declared run-specific fields ignored
! Run-specific fields excepted: in solve3D_manifest.txt the lines run_id, layout_digest,
! stack_digest and checksum and the path ending each artifact line; in sigma2_state.bin the layout
! digest (header word 14, bytes 121-128). Project files and text files hold absolute paths, so the
! records of every project segment, and the lines of a text file, are compared with the run
! directory's own prefix removed. A declared pattern that selects no file fails.
module simple_test_identity
use, intrinsic :: iso_fortran_env, only: int8, int64
use simple_string, only: string
implicit none

public :: compare_run_dirs, same_bytes, same_manifest, manifest_line_kept, SIGMA2_DIGEST_BYTES
private
#include "simple_local_flags.inc"

character(len=*), parameter :: MANIFEST_FNAME = 'solve3D_manifest.txt'
character(len=*), parameter :: SIGMA2_FNAME   = 'sigma2_state.bin'
character(len=*), parameter :: SIGMA2_MAGIC   = 'SIMPLE_SIGMA2_V2'
integer,          parameter :: SIGMA2_DIGEST_BYTES(2) = [121, 128] !< after the 16-byte magic, word 14 of 8 bytes

contains

    !> Compares the files the patterns (blank-separated, relative to the run directories) select in
    !! dir1 and in dir2; one SAME or DIFF line per file goes to report. ndiff counts the files that
    !! differ or exist in one directory only.
    subroutine compare_run_dirs( dir1, dir2, patterns, report, nsame, ndiff )
        use simple_syslib,       only: simple_list_files
        use simple_string_utils, only: parsestr
        character(len=*), intent(in)  :: dir1, dir2, patterns, report
        integer,          intent(out) :: nsame, ndiff
        type(string), allocatable :: list1(:), list2(:)
        character(len=512) :: pats(64)
        character(len=len(patterns)) :: patcopy
        character(len=:), allocatable :: rel
        integer :: funit, npat, ipat, i
        logical :: l_same
        nsame = 0
        ndiff = 0
        open(newunit=funit, file=report, action='write', status='replace')
        write(funit,'(a)') 'identity of '//trim(dir1)//' and '//trim(dir2)
        patcopy = patterns
        call parsestr(patcopy, ' ', pats, npat)
        do ipat = 1, npat
            call simple_list_files(trim(dir1)//'/'//trim(pats(ipat)), list1)
            call simple_list_files(trim(dir2)//'/'//trim(pats(ipat)), list2)
            ! a declaration that matches nothing compares nothing: it fails
            if( size(list1) == 0 .and. size(list2) == 0 ) call record('NO MATCH: '//trim(pats(ipat)), .false.)
            do i = 1, size(list1)
                rel = relative(list1(i)%to_char(), dir1)
                if( has(list2, dir2, rel) )then
                    l_same = same_file(trim(dir1)//'/'//rel, trim(dir2)//'/'//rel, dir1, dir2)
                    call record(merge('SAME ', 'DIFF ', l_same)//rel, l_same)
                else
                    call record('MISSING in dir2: '//rel, .false.)
                endif
            enddo
            do i = 1, size(list2)
                rel = relative(list2(i)%to_char(), dir2)
                if( .not. has(list1, dir1, rel) ) call record('MISSING in dir1: '//rel, .false.)
            enddo
        enddo
        write(funit,'(a,i0,a,i0)') 'files the same: ', nsame, '; differing or missing: ', ndiff
        close(funit)

      contains

        subroutine record( line, ok )
            character(len=*), intent(in) :: line
            logical,          intent(in) :: ok
            write(funit,'(a)') line
            if( ok )then
                nsame = nsame + 1
            else
                ndiff = ndiff + 1
            endif
        end subroutine record

        function relative( path, dir ) result( rel )
            character(len=*), intent(in)  :: path, dir
            character(len=:), allocatable :: rel
            rel = path(len_trim(dir)+2:)
        end function relative

        logical function has( list, dir, rel )
            type(string),     intent(in) :: list(:)
            character(len=*), intent(in) :: dir, rel
            integer :: j
            has = .false.
            do j = 1, size(list)
                if( list(j)%to_char() == trim(dir)//'/'//rel )then
                    has = .true.
                    return
                endif
            enddo
        end function has

    end subroutine compare_run_dirs

    !> the comparison rule of the file's name; dir1 and dir2 are the run directories the two files
    !! lie in, whose prefixes the records and lines of projects and text files are compared without
    logical function same_file( f1, f2, dir1, dir2 )
        use simple_fileio, only: basename
        character(len=*), intent(in) :: f1, f2, dir1, dir2
        type(string) :: bname
        integer :: n
        bname = basename(string(f1))
        n     = len_trim(f1)
        if( bname%to_char() == MANIFEST_FNAME )then
            same_file = same_manifest(f1, f2)
        else if( bname%to_char() == SIGMA2_FNAME )then
            if( has_sigma2_magic(f1) )then
                same_file = same_bytes(f1, f2, SIGMA2_DIGEST_BYTES)
            else
                same_file = same_bytes(f1, f2)
            endif
        else if( n > 7 .and. f1(max(1,n-6):n) == '.simple' )then
            same_file = same_project_records(f1, f2, dir1, dir2)
        else if( n > 4 .and. f1(max(1,n-3):n) == '.txt' )then
            same_file = same_text_lines(f1, f2, dir1, dir2)
        else
            same_file = same_bytes(f1, f2)
        endif
    end function same_file

    !> byte identity, bytes skip(1) to skip(2) (1-based) excepted
    logical function same_bytes( f1, f2, skip )
        character(len=*),  intent(in) :: f1, f2
        integer, optional, intent(in) :: skip(2)
        integer, parameter :: CHUNK = 1048576
        integer(int8), allocatable :: b1(:), b2(:)
        integer(int64) :: sz1, sz2, pos, nread, k
        integer :: u1, u2, ios1, ios2
        same_bytes = .false.
        inquire(file=f1, size=sz1)
        inquire(file=f2, size=sz2)
        if( sz1 < 0 .or. sz1 /= sz2 ) return
        open(newunit=u1, file=f1, access='stream', form='unformatted', action='read', status='old', iostat=ios1)
        open(newunit=u2, file=f2, access='stream', form='unformatted', action='read', status='old', iostat=ios2)
        if( ios1 == 0 .and. ios2 == 0 )then
            allocate(b1(CHUNK), b2(CHUNK))
            same_bytes = .true.
            pos = 0_int64
            do while( pos < sz1 )
                nread = min(int(CHUNK, int64), sz1 - pos)
                read(u1, iostat=ios1) b1(1:nread)
                read(u2, iostat=ios2) b2(1:nread)
                if( ios1 /= 0 .or. ios2 /= 0 )then
                    same_bytes = .false.
                    exit
                endif
                if( present(skip) )then
                    do k = max(1_int64, skip(1) - pos), min(nread, skip(2) - pos)
                        b2(k) = b1(k)
                    enddo
                endif
                if( any(b1(1:nread) /= b2(1:nread)) )then
                    same_bytes = .false.
                    exit
                endif
                pos = pos + nread
            enddo
        endif
        if( ios1 == 0 ) close(u1)
        if( ios2 == 0 ) close(u2)
    end function same_bytes

    !> the manifests agree line by line, the run-specific lines and the artifact paths left out
    logical function same_manifest( f1, f2 )
        character(len=*), intent(in) :: f1, f2
        character(len=4096) :: l1, l2
        integer :: u1, u2, ios1, ios2
        same_manifest = .false.
        open(newunit=u1, file=f1, action='read', status='old', iostat=ios1)
        if( ios1 /= 0 ) return
        open(newunit=u2, file=f2, action='read', status='old', iostat=ios2)
        if( ios2 /= 0 )then
            close(u1)
            return
        endif
        do
            call next_kept(u1, l1, ios1)
            call next_kept(u2, l2, ios2)
            if( ios1 /= 0 .or. ios2 /= 0 ) exit
            if( l1 /= l2 ) exit
        enddo
        ! both ended together, with every kept line equal
        same_manifest = ios1 /= 0 .and. ios2 /= 0
        close(u1)
        close(u2)

      contains

        subroutine next_kept( u, line, ios )
            integer,          intent(in)  :: u
            character(len=*), intent(out) :: line
            integer,          intent(out) :: ios
            do
                read(u,'(a)',iostat=ios) line
                if( ios /= 0 ) return
                if( manifest_line_kept(line) ) exit
            enddo
        end subroutine next_kept

    end function same_manifest

    !> false for a run-specific manifest line; an artifact line loses its path (its last field)
    logical function manifest_line_kept( line )
        character(len=*), intent(inout) :: line
        character(len=32) :: key
        integer :: pos
        key = adjustl(line)
        pos = index(key, ' ')
        if( pos > 0 ) key = key(:pos-1)
        select case(trim(key))
            case('run_id', 'layout_digest', 'stack_digest', 'checksum')
                manifest_line_kept = .false.
            case('artifact')
                manifest_line_kept = .true.
                pos = index(trim(line), ' ', back=.true.)
                if( pos > 0 ) line = line(:pos-1)
            case default
                manifest_line_kept = .true.
        end select
    end function manifest_line_kept

    logical function has_sigma2_magic( f )
        character(len=*), intent(in) :: f
        character(len=len(SIGMA2_MAGIC)) :: magic
        integer :: u, ios
        has_sigma2_magic = .false.
        open(newunit=u, file=f, access='stream', form='unformatted', action='read', status='old', iostat=ios)
        if( ios /= 0 ) return
        read(u, iostat=ios) magic
        close(u)
        has_sigma2_magic = ios == 0 .and. magic == SIGMA2_MAGIC
    end function has_sigma2_magic

    !> every segment of two projects agrees, record by record, with each run directory's prefix
    !! removed from the records
    logical function same_project_records( f1, f2, dir1, dir2 )
        use simple_sp_project, only: sp_project
        use simple_oris,       only: oris
        character(len=*), intent(in) :: f1, f2, dir1, dir2
        integer, parameter :: NSEGS = 11
        class(sp_project), allocatable :: p1, p2
        type(string) :: pre1, pre2
        integer :: iseg, i
        allocate(p1, p2)
        call p1%read(string(f1))
        call p2%read(string(f2))
        pre1 = string(trim(dir1))
        pre2 = string(trim(dir2))
        same_project_records = .true.
        do iseg = 1, NSEGS
            if( .not. same_oris(segment(p1, iseg), segment(p2, iseg)) )then
                same_project_records = .false.
                exit
            endif
        enddo
        call p1%kill
        call p2%kill
        deallocate(p1, p2)

      contains

        function segment( p, iseg ) result( os )
            class(sp_project), intent(in) :: p
            integer,           intent(in) :: iseg
            type(oris) :: os
            select case(iseg)
                case(1);  os = p%os_mic
                case(2);  os = p%os_stk
                case(3);  os = p%os_ptcl2D
                case(4);  os = p%os_cls2D
                case(5);  os = p%os_cls3D
                case(6);  os = p%os_ptcl3D
                case(7);  os = p%os_out
                case(8);  os = p%os_optics
                case(9);  os = p%projinfo
                case(10); os = p%jobproc
                case default; os = p%compenv
            end select
        end function segment

        logical function same_oris( a, b )
            type(oris), intent(in) :: a, b
            type(string) :: sa, sb
            same_oris = a%get_noris() == b%get_noris()
            if( .not. same_oris ) return
            do i = 1, a%get_noris()
                sa = a%ori2str(i)
                sb = b%ori2str(i)
                sa = sa%substr_remove(pre1)
                sb = sb%substr_remove(pre2)
                if( .not. (sa == sb) )then
                    same_oris = .false.
                    return
                endif
            enddo
        end function same_oris

    end function same_project_records

    !> two text files agree line by line, each run directory's prefix removed from its lines
    logical function same_text_lines( f1, f2, dir1, dir2 )
        character(len=*), intent(in) :: f1, f2, dir1, dir2
        character(len=4096) :: l1, l2
        type(string) :: s1, s2, pre1, pre2
        integer :: u1, u2, ios1, ios2
        same_text_lines = .false.
        open(newunit=u1, file=f1, action='read', status='old', iostat=ios1)
        if( ios1 /= 0 ) return
        open(newunit=u2, file=f2, action='read', status='old', iostat=ios2)
        if( ios2 /= 0 )then
            close(u1)
            return
        endif
        pre1 = string(trim(dir1))
        pre2 = string(trim(dir2))
        do
            read(u1,'(a)',iostat=ios1) l1
            read(u2,'(a)',iostat=ios2) l2
            if( ios1 /= 0 .or. ios2 /= 0 ) exit
            s1 = string(trim(l1))
            s2 = string(trim(l2))
            s1 = s1%substr_remove(pre1)
            s2 = s2%substr_remove(pre2)
            if( .not. (s1 == s2) )then
                close(u1)
                close(u2)
                return
            endif
        enddo
        close(u1)
        close(u2)
        ! both ended together
        same_text_lines = ios1 /= 0 .and. ios2 /= 0
    end function same_text_lines

end module simple_test_identity
