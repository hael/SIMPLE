#include "simple_lib.f08"
module simple_binoris_io
use simple_defs
use simple_binoris,      only: binoris
use simple_ori,          only: ori
use simple_oris,         only: oris
use simple_prime3D_srch, only: prime3D_srch
use simple_fileio,       only: file_exists, nlines, fileio_errmsg
use simple_syslib,       only: alloc_errchk
use simple_strings,      only: str_has_substr
implicit none

interface binread_oritab
    module procedure binread_oritab_1
    module procedure binread_oritab_2
end interface

interface binwrite_oritab
    module procedure binwrite_oritab_1
!    module procedure binwrite_oritab_2
end interface

contains

    subroutine binread_oritab_1( fname, a, fromto, nst )
        character(len=*),      intent(in)    :: fname
        class(oris),           intent(inout) :: a
        integer,               intent(in)    :: fromto(2)
        integer,     optional, intent(out)   :: nst
        type(binoris) :: bos
        integer       :: irec
        ! if( .not. file_exists(fname) )&
!            call fileio_errmsg('binoris_io :: binread_oritab_1  '//  trim(fname) // 'does not exist in cwd') 
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_oritab_1'
        endif
        if( str_has_substr(fname,'.txt') )then
            call a%read(fname, nst)
        else
            call bos%open(fname)
            do irec=fromto(1),fromto(2)
                call bos%read_record(irec, a, nst=nst)
            end do
            call bos%close
        endif
    end subroutine binread_oritab_1

    subroutine binread_oritab_2( fname, a, fromto, primesrch3D, nst )
        character(len=*),    intent(in)    :: fname
        class(oris),         intent(inout) :: a
        integer,             intent(in)    :: fromto(2)
        class(prime3D_srch), intent(inout) :: primesrch3D(fromto(1):fromto(2))
        integer, optional,   intent(out)   :: nst
        type(binoris) :: bos
        integer       :: irec
        type(oris)    :: os_peak
        !if( .not. file_exists(fname))&
         !   call fileio_errmsg('binoris_io :: binread_oritab_2  '//  trim(fname) // 'does not exist in cwd')
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_oritab_2'
        endif 
        call bos%open(fname)
        do irec=fromto(1),fromto(2)
            call bos%read_record(irec, a, os_peak, nst)
            call primesrch3D(irec)%set_o_peaks(os_peak)
        end do
        call bos%close
    end subroutine binread_oritab_2
 
    subroutine binread_ctfparams_and_state( fname, a, fromto )
        character(len=*), intent(in)    :: fname
        class(oris),      intent(inout) :: a
        integer,          intent(in)    :: fromto(2)
        type(binoris) :: bos
        integer       :: irec
     !   if( .not. file_exists(fname) )&
     !       call fileio_errmsg('binoris_io :: binread_ctfparams_and_state  '//  trim(fname) // 'does not exist in cwd')
                if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_ctfparams_and_state'
        endif
        if( str_has_substr(fname,'.txt') )then
            call a%read_ctfparams_and_state(fname)
        else
            call bos%open(fname)
            do irec=fromto(1),fromto(2)
                call bos%read_ctfparams_and_state(irec, a)
            end do
            call bos%close
        endif
    end subroutine binread_ctfparams_and_state

    function binread_nlines( fname ) result( nl )
        character(len=*), intent(in) :: fname
        type(binoris) :: bos
        integer       :: nl
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_nlines'
        endif
        if( str_has_substr(fname,'.txt') )then
            nl = nlines(fname)
        else
            call bos%open(fname)
            nl = bos%get_n_records()
            call bos%close
        endif
    end function binread_nlines
 
    subroutine binwrite_oritab_1( fname, a, fromto, mask, fname_fill_in )
        character(len=*),           intent(in)    :: fname
        class(oris),                intent(inout) :: a
        integer,                    intent(in)    :: fromto(2)
        logical,          optional, intent(in)    :: mask(fromto(1):fromto(2))
        character(len=*), optional, intent(in)    :: fname_fill_in
        type(binoris) :: bos
        type(binoris) :: bos_fill_in
        integer       :: irec
        logical       :: mmask(fromto(1):fromto(2)), l_fill_in
        if( str_has_substr(fname,'.txt') )then
            call a%write(fname, fromto)
            return
        endif
        mmask     = .true.
        l_fill_in = .false.
        if( present(mask) )then
            mmask = mask
            if( .not. present(fname_fill_in) )&
            &stop 'need fill in file in conjunction with mask; binoris_io :: binwrite_oritab_2'
           ! if( .not. file_exists(fname_fill_in) )&
            !    call fileio_errmsg('binoris_io :: binwrite_oritab_1  '//trim(fname_fill_in)// 'does not exist in cwd')
            if( .not. file_exists(fname_fill_in) )then
                write(*,*) 'file: ', trim(fname_fill_in)
                stop 'does not exist in cwd; binoris_io :: binwrite_oritab_2'
            endif
            l_fill_in = .true.
        endif
        ! establish outfile header
        call bos%new(a, fromto)
        call bos%open(fname, del_if_exists=.true.)
        call bos%write_header()
        if( l_fill_in )then
            ! establish fill-in filehandler
            call bos_fill_in%open(fname_fill_in)
        endif
        do irec=fromto(1),fromto(2)
            if( mmask(irec) )then
                ! particle has been subjected to update
                call bos%write_record(irec, a)
            else
                ! info in fill-in file is needed
                call bos_fill_in%read_record(irec)
                call bos%write_record(irec, bos_fill_in)
            endif
        end do
        call bos%close
    end subroutine binwrite_oritab_1
 
 !    subroutine binwrite_oritab_2( fname, a, fromto, primesrch3D, mask, fname_fill_in )
!         character(len=*),           intent(in)    :: fname
!         class(oris),                intent(inout) :: a
!         integer,                    intent(in)    :: fromto(2)
!         class(prime3D_srch),        intent(inout) :: primesrch3D(fromto(1):fromto(2))
!         logical,          optional, intent(in)    :: mask(fromto(1):fromto(2))
!         character(len=*), optional, intent(in)    :: fname_fill_in
!         type(binoris) :: bos
!         type(binoris) :: bos_fill_in
!         type(oris)    :: os_peak, os_peak_fill_in, os_peak_conforming
!         type(ori)     :: o_single
!         integer       :: npeaks, npeaks_fill_in, irec
!         logical       :: reshape_fill_in, mmask(fromto(1):fromto(2)), l_fill_in
!         if( str_has_substr(fname,'.txt') )then
!             stop 'not intended for text file output; binoris_io :: binwrite_oritab_2'
!         endif
!         mmask     = .true.
!         l_fill_in = .false.
!         if( present(mask) )then
!             mmask = mask
!             if( .not. present(fname_fill_in) )&
!             &stop 'need fill in file in conjunction with mask; binoris_io :: binwrite_oritab_2'
!           !  if( .not. file_exists(fname_fill_in) )&
!            !     call fileio_errmsg('binoris_io :: binwrite_oritab_2  '//trim(fname_fill_in)//'does not exist in cwd')
!                         if( .not. file_exists(fname_fill_in) )then
!                 write(*,*) 'file: ', trim(fname_fill_in)
!                 stop 'does not exist in cwd; binoris_io :: binwrite_oritab_2'
!             endif
!             l_fill_in = .true.
!         endif
!         ! establish outfile header
!         o_single = a%get_ori(fromto(1))
!         call primesrch3D(fromto(1))%get_oris(os_peak, o_single)
!         npeaks = os_peak%get_noris()
!         call bos%new(a, fromto, os_peak)
!         call bos%open(fname, del_if_exists=.true.)
!         call bos%write_header()
!         if( l_fill_in )then
!             ! establish fill-in filehandler
!             call bos_fill_in%open(fname_fill_in)
!             npeaks_fill_in  = bos_fill_in%get_n_peaks()
!             reshape_fill_in = npeaks_fill_in /= npeaks
!         endif
!         do irec=fromto(1),fromto(2)
!             if( mmask(irec) )then
!                 ! particle has been subjected to update
!                 o_single = a%get_ori(irec)
!                 call primesrch3D(irec)%get_oris(os_peak, o_single)
!                 call bos%write_record(irec, a, os_peak)
!             else
!                 ! info in fill-in file is needed
!                 if( reshape_fill_in )then
!                     call bos_fill_in%read_record(irec, a, os_peak_fill_in)
!                     os_peak_conforming =&
!                     &os_peak_fill_in%create_conforming_npeaks_set(npeaks)
!                     call bos%write_record(irec, a, os_peak_conforming)
!                 else
!                     call bos_fill_in%read_record(irec)
!                     call bos%write_record(irec, bos_fill_in)
!                 endif
!             endif
!         end do
!         call bos%kill
!         call bos_fill_in%kill
!         call os_peak%kill
!         call os_peak_fill_in%kill
!         call os_peak_conforming%kill
!     end subroutine binwrite_oritab_2
 
    !
end module simple_binoris_io
