module simple_binoris_io
use simple_binoris,      only: binoris
use simple_ori,          only: ori
use simple_oris,         only: oris
use simple_prime3D_srch, only: prime3D_srch
use simple_fileio,       only: file_exists
implicit none
 
interface binread_oritab
    module procedure binread_oritab_1
    module procedure binread_oritab_2
end interface
 
interface binwrite_oritab
    module procedure binwrite_oritab_1
    module procedure binwrite_oritab_2
end interface
 
contains
 
    subroutine binread_oritab_1( fname, a, fromto, nst )
        character(len=*),      intent(in)    :: fname
        class(oris),           intent(inout) :: a
        integer,               intent(in)    :: fromto(2)
        integer,     optional, intent(out)   :: nst
        type(binoris) :: bos
        integer       :: irec
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_oritab_1'
        endif
        call bos%open(fname)
        do irec=fromto(1),fromto(2)
            call bos%read_record(irec, a, nst=nst)
        end do
        call bos%close
    end subroutine binread_oritab_1
 
     subroutine binread_oritab_2( fname, a, fromto, primesrch3D, mask )
         character(len=*),    intent(in)    :: fname
         class(oris),         intent(inout) :: a
         integer,             intent(in)    :: fromto(2)
         class(prime3D_srch), intent(inout) :: primesrch3D(fromto(1):fromto(2))
         logical,             intent(in)    :: mask(fromto(1):fromto(2))
!         type(binoris) :: bos
!         type(binoris) :: bos_fill_in
!         if( .not. file_exists(fname) )then
!             write(*,*) 'file: ', trim(fname)
!             stop 'does not exist in cwd; binoris_io :: binread_oritab_2'
!         endif
!         ! establish input file handler
!         call bos%open(fname)
 
! !!!!!!!!!!!!!!! 2 BE CONTINUED
 
     end subroutine binread_oritab_2
 
    subroutine binread_ctfparams_and_state( fname, a, fromto )
        character(len=*), intent(in)    :: fname
        class(oris),      intent(inout) :: a
        integer,          intent(in)    :: fromto(2)
        type(binoris) :: bos
        integer       :: irec
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_ctfparams_and_state'
        endif
        call bos%open(fname)
        do irec=fromto(1),fromto(2)
            call bos%read_ctfparams_and_state(irec, a)
        end do
        call bos%close
    end subroutine binread_ctfparams_and_state
 
    subroutine binwrite_oritab_1( fname, a, fromto )
        character(len=*), intent(in)    :: fname
        class(oris),      intent(inout) :: a
        integer,          intent(in)    :: fromto(2)
        type(binoris) :: bos
        integer       :: irec
        call bos%new(a, fromto)
        call bos%open(fname, del_if_exists=.true.)
        call bos%write_header()
        do irec=fromto(1),fromto(2)
            call bos%write_record(irec, a)
        end do
        call bos%close
    end subroutine binwrite_oritab_1
 
    subroutine binwrite_oritab_2( fname, fname_fill_in, a, fromto, primesrch3D, mask )
        character(len=*),    intent(in)    :: fname, fname_fill_in
        class(oris),         intent(inout) :: a
        integer,             intent(in)    :: fromto(2)
        class(prime3D_srch), intent(inout) :: primesrch3D(fromto(1):fromto(2))
        logical,             intent(in)    :: mask(fromto(1):fromto(2))
        type(binoris) :: bos
        type(binoris) :: bos_fill_in
        type(oris)    :: os_peak, os_peak_fill_in, os_peak_conforming
        type(ori)     :: o_single
        integer       :: npeaks, npeaks_fill_in, i
        logical       :: reshape_fill_in
        if( .not. file_exists(fname_fill_in) )then
            write(*,*) 'file: ', trim(fname_fill_in)
            stop 'does not exist in cwd; binoris_io :: binwrite_oritab_2'
        endif
        ! establish outfile header
        o_single = a%get_ori(fromto(1))
        call primesrch3D(fromto(1))%get_oris(os_peak, o_single)
        npeaks = os_peak%get_noris()
        call bos%new(a, fromto, os_peak)
        call bos%open(fname, del_if_exists=.true.)
        call bos%write_header()
        ! establish fill-in filehandler
        call bos_fill_in%open(fname_fill_in)
        npeaks_fill_in  = bos_fill_in%get_n_peaks()
        reshape_fill_in = npeaks_fill_in /= npeaks
        do i=fromto(1),fromto(2)
            if( mask(i) )then
                ! particle has been subjected to update
                o_single = a%get_ori(i)
                call primesrch3D(i)%get_oris(os_peak, o_single)
                call bos%write_record(i, a, os_peak)
            else
                ! info in fill-in file is needed
                if( reshape_fill_in )then
                    call bos_fill_in%read_record(i, a, os_peak_fill_in)
                    os_peak_conforming =&
                    &os_peak_fill_in%create_conforming_npeaks_set(npeaks)
                    call bos%write_record(i, a, os_peak_conforming)
                else
                    call bos_fill_in%read_record(i)
                    call bos%write_record(i, bos_fill_in)
                endif
            endif
        end do
        call bos%kill
        call bos_fill_in%kill
        call os_peak%kill
        call os_peak_fill_in%kill
        call os_peak_conforming%kill
    end subroutine binwrite_oritab_2
 
end module simple_binoris_io
