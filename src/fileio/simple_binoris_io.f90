module simple_binoris_io
use simple_defs
use simple_strings
use simple_fileio
use simple_oris,       only: oris
use simple_sp_project, only: sp_project
implicit none

contains

    subroutine binread_oritab( fname, spproj, a, fromto )
        character(len=*),  intent(in)    :: fname
        class(sp_project), intent(inout) :: spproj
        class(oris),       intent(inout) :: a
        integer,           intent(in)    :: fromto(2)
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_oritab_1'
        endif
        select case(fname2format(fname))
            case('O')
                call spproj%read(fname)
            case('T')
                call a%read(fname, fromto=fromto)
            case DEFAULT
                write(*,*) 'file: ', trim(fname)
                stop 'format unsupported; simple_binoris_io :: binread_oritab'
        end select
    end subroutine binread_oritab

    subroutine binread_ctfparams_state_eo( fname, spproj, a, fromto )
        character(len=*),  intent(in)    :: fname
        class(sp_project), intent(inout) :: spproj
        class(oris),       intent(inout) :: a
        integer,           intent(in), optional    :: fromto(2)
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_ctfparams_and_state'
        endif
        select case(fname2format(fname))
            case('O')
                call spproj%read_ctfparams_state_eo(fname)
            case('T')
                call a%read_ctfparams_state_eo(fname)
            case DEFAULT
                write(*,*) 'file: ', trim(fname)
                stop 'format unsupported; simple_binoris_io :: binread_ctfparams_state_eo'
        end select
    end subroutine binread_ctfparams_state_eo

    function binread_nlines( fname ) result( nl )
        use simple_binoris,    only: binoris
        use simple_parameters, only: params_glob
        character(len=*), intent(in) :: fname
        integer       :: nl
        type(binoris) :: bos
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_nlines'
        endif
        select case(fname2format(fname))
            case('O')
                call bos%open(fname)
                nl = bos%get_n_records(params_glob%spproj_iseg)
                call bos%close
            case('T')
                nl = nlines(fname)
            case DEFAULT
                write(*,*) 'file: ', trim(fname)
                stop 'format unsupported; simple_binoris_io :: binread_nlines'
        end select
    end function binread_nlines

    subroutine binwrite_oritab( fname, spproj, a, fromto, isegment )
        character(len=*),  intent(in)    :: fname
        class(sp_project), intent(inout) :: spproj
        class(oris),       intent(inout) :: a
        integer,           intent(in)    :: fromto(2)
        integer, optional, intent(in)    :: isegment
        select case(fname2format(fname))
            case('O')
                call spproj%write(fname, fromto, isegment)
            case('T')
                call a%write(fname, fromto)
            case DEFAULT
                write(*,*) 'file: ', trim(fname)
                stop 'format unsupported; simple_binoris_io :: binwrite_oritab'
        end select
    end subroutine binwrite_oritab

    subroutine binwrite_o_peaks( fname, fromto, o_peaks, isegment )
        use simple_binoris, only: binoris
        character(len=*), intent(in)    :: fname
        integer,          intent(in)    :: fromto(2), isegment
        class(oris),      intent(inout) :: o_peaks(fromto(1):fromto(2))
        type(str4arr),    allocatable   :: o_peaks_str(:,:)
        type(oris)      :: o_peaks_on_disc(fromto(1):fromto(2))
        type(binoris)   :: bos
        integer(kind=8) :: ibytes
        integer         :: strlen, strlen_max, i, ipeak, npeaks, npeaks_on_disc, nrec
        real            :: ow
        logical         :: l_did_read
        ! read peak data if file exists
        l_did_read = .false.
        if( file_exists(fname) )then
            call binread_o_peaks(fname, fromto, o_peaks_on_disc, isegment)
            l_did_read = .true.
        endif
        ! create merged string representation (o_peaks + o_peaks_on_disc)
        allocate(o_peaks_str(fromto(1):fromto(2),NPEAKS2REFINE))
        strlen_max = 0
        nrec       = 0
        do i=fromto(1),fromto(2)
            npeaks = o_peaks(i)%get_noris()
            npeaks_on_disc = 0
            if( l_did_read ) npeaks_on_disc = o_peaks_on_disc(i)%get_noris()
            if( npeaks > 1 )then
                do ipeak=1,npeaks
                    ow = o_peaks(i)%get(ipeak, 'ow')
                    if( ow > TINY )then
                        o_peaks_str(i,ipeak)%str = o_peaks(i)%ori2str(ipeak)
                        strlen     = len_trim(o_peaks_str(i,ipeak)%str)
                        strlen_max = max(strlen,strlen_max)
                        nrec       = nrec + 1
                    endif
                end do
            else if( npeaks_on_disc > 1 )then
                do ipeak=1,npeaks_on_disc
                    ow = o_peaks_on_disc(i)%get(ipeak, 'ow')
                    if( ow > TINY )then
                        o_peaks_str(i,ipeak)%str = o_peaks_on_disc(i)%ori2str(ipeak)
                        strlen     = len_trim(o_peaks_str(i,ipeak)%str)
                        strlen_max = max(strlen,strlen_max)
                        nrec       = nrec + 1
                    endif
                end do
            endif
        enddo
        ! be mindful of memory use
        if( l_did_read )then
            do i=fromto(1),fromto(2)
                call o_peaks_on_disc(i)%kill
            enddo
        endif
        ! write data to binary file
        call bos%open(fname, del_if_exists=.true.)
        ! add segment to stack
        call bos%add_segment(isegment, [1,nrec], strlen_max)
        ! update byte ranges in header
        call bos%update_byte_ranges
        ! init byte counter
        ibytes = bos%get_first_data_byte(isegment)
        ! write records
        do i=fromto(1),fromto(2)
            do ipeak=1,NPEAKS2REFINE
                if( allocated(o_peaks_str(i,ipeak)%str) )then
                    call bos%write_record(isegment, ibytes, o_peaks_str(i,ipeak)%str)
                    deallocate(o_peaks_str(i,ipeak)%str)
                endif
            enddo
        enddo
        deallocate(o_peaks_str)
        call bos%write_header
        call bos%close
    end subroutine binwrite_o_peaks

    subroutine binread_o_peaks( fname, fromto, o_peaks, isegment )
        use simple_binoris, only: binoris
        use simple_ori,     only: ori
        character(len=*), intent(in)    :: fname
        integer,          intent(in)    :: fromto(2), isegment
        class(oris),      intent(inout) :: o_peaks(fromto(1):fromto(2))
        type(str4arr),    allocatable   :: os_strings(:)
        type(binoris) :: bos
        type(oris)    :: os_tmp
        type(ori)     :: o
        integer       :: fromto_in_file(2), cnt, i, j, iptcl, iptcl_prev
        ! temporary container
        call os_tmp%new(NPEAKS2REFINE)
        ! read segment into string representation
        call bos%open(fname)
        fromto_in_file = bos%get_fromto(isegment)
        allocate(os_strings(fromto_in_file(1):fromto_in_file(2)))
        call bos%read_segment(isegment, os_strings)
        call bos%close
        ! transfer orientation info to o_peaks
        iptcl = 0
        cnt   = 0
        do i=fromto_in_file(1),fromto_in_file(2)
            iptcl_prev = iptcl
            call o%str2ori(os_strings(i)%str)
            deallocate(os_strings(i)%str)
            ! it is assumed that all entries have an iptcl index
            iptcl = nint(o%get('iptcl'))
            if( iptcl == 0 )then
                write(*,*) 'ERROR! iptcl index == 0 not allowed, index probably not set'
                stop 'simple_binoris_io :: binread_o_peaks'
            endif
            if( iptcl_prev == 0 .or. iptcl == iptcl_prev )then
                cnt = cnt + 1
                call os_tmp%set_ori(cnt, o)
            else
                ! transfer to o_peaks
                call o_peaks(iptcl_prev)%new(cnt)
                do j=1,cnt
                    call o_peaks(iptcl_prev)%set_ori(j, os_tmp%get_ori(j))
                end do
                ! reset counter and update os_tmp
                cnt = 1
                call os_tmp%set_ori(cnt, o)
            endif
        end do
        ! transfer last set to o_peaks
        call o_peaks(iptcl_prev)%new(cnt)
        do j=1,cnt
            call o_peaks(iptcl_prev)%set_ori(j, os_tmp%get_ori(j))
        end do
        ! destruct
        call os_tmp%kill
        deallocate(os_strings)
    end subroutine binread_o_peaks

end module simple_binoris_io
