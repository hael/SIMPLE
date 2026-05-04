!@descr: particle polar Fourier line streaming I/O
submodule (simple_polarft_calc) simple_polarft_ops_io
implicit none
#include "simple_local_flags.inc"
contains

    module subroutine write_ptcl_pft_range( self, fname, nptcls_total, iwrite_from, iwrite_to )
        class(polarft_calc), intent(in) :: self
        class(string),       intent(in) :: fname
        integer,             intent(in) :: nptcls_total
        integer,             intent(in) :: iwrite_from, iwrite_to
        integer :: funit, io_stat, nwrite
        integer(kind=8) :: header_bytes, pft_bytes, write_pos
        if( .not. self%existence ) THROW_HARD('polarft_calc does not exist; write_ptcl_pft_range')
        if( nptcls_total < 1 ) THROW_HARD('nptcls_total must be > 0; write_ptcl_pft_range')
        if( iwrite_from < 1 .or. iwrite_to < iwrite_from .or. iwrite_to > nptcls_total )then
            write(logfhandle,*) 'iwrite_from, iwrite_to, nptcls_total: ', iwrite_from, iwrite_to, nptcls_total
            THROW_HARD('invalid output range; write_ptcl_pft_range')
        endif
        nwrite = iwrite_to - iwrite_from + 1
        if( nwrite /= self%nptcls )then
            write(logfhandle,*) 'nwrite, self%nptcls: ', nwrite, self%nptcls
            THROW_HARD('current particle batch does not match requested output range; write_ptcl_pft_range')
        endif
        if( iwrite_from == 1 )then
            call open_pft_or_ctf2_array_for_write(fname, funit)
            write(unit=funit,pos=1) [self%pftsz, self%kfromto(1), self%interpklim, nptcls_total]
        else
            call fopen(funit, fname, access='STREAM', action='WRITE', status='OLD', iostat=io_stat)
            call fileiochk('write_ptcl_pft_range: '//fname%to_char(), io_stat)
        endif
        header_bytes = 4_int8 * int(sizeof(funit), kind=8)
        pft_bytes    = int(self%pftsz, kind=8) * int(self%interpklim - self%kfromto(1) + 1, kind=8) * &
                       2_int8 * int(storage_size(real(1._sp))/8, kind=8)
        write_pos    = header_bytes + 1_int8 + int(iwrite_from - 1, kind=8) * pft_bytes
        write(unit=funit,pos=write_pos) cmplx(self%pfts_ptcls(:,:,1:self%nptcls), kind=sp)
        call fclose(funit)
    end subroutine write_ptcl_pft_range

    module subroutine write_ref_pfts( self, fname, iseven )
        class(polarft_calc), intent(in) :: self
        class(string),       intent(in) :: fname
        logical,             intent(in) :: iseven
        integer :: funit
        integer :: header(4)
        integer(kind=8) :: payload_pos
        if( .not. self%existence ) THROW_HARD('polarft_calc does not exist; write_ref_pfts')
        header      = [self%pftsz, self%kfromto(1), self%kfromto(2), self%nrefs]
        payload_pos = int(sizeof(header), kind=8) + 1_int8
        call open_pft_or_ctf2_array_for_write(fname, funit)
        write(unit=funit,pos=1) header
        if( iseven )then
            write(unit=funit,pos=payload_pos) self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),:)
        else
            write(unit=funit,pos=payload_pos) self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),:)
        endif
        call fclose(funit)
    end subroutine write_ref_pfts

    module subroutine read_ref_pfts( self, fname, iseven )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: fname
        logical,             intent(in)    :: iseven
        integer :: funit, io_stat
        integer :: header(4)
        integer(kind=8) :: payload_pos
        if( .not. self%existence ) THROW_HARD('polarft_calc does not exist; read_ref_pfts')
        if( .not. file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist; read_ref_pfts')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_ref_pfts: '//fname%to_char(), io_stat)
        read(unit=funit,pos=1) header
        if( header(1) /= self%pftsz .or. header(2) /= self%kfromto(1) .or. &
            &header(3) /= self%kfromto(2) .or. header(4) /= self%nrefs )then
            write(logfhandle,*) 'reprojection model header mismatch in: ', fname%to_char()
            write(logfhandle,*) 'expected: ', self%pftsz, self%kfromto(1), self%kfromto(2), self%nrefs
            write(logfhandle,*) 'found:    ', header
            THROW_HARD('incompatible reprojection model header; read_ref_pfts')
        endif
        payload_pos = int(sizeof(header), kind=8) + 1_int8
        if( iseven )then
            read(unit=funit,pos=payload_pos) self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),:)
        else
            read(unit=funit,pos=payload_pos) self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),:)
        endif
        call fclose(funit)
    end subroutine read_ref_pfts

    subroutine open_pft_or_ctf2_array_for_write( fname, funit )
        class(string), intent(in)  :: fname
        integer,       intent(out) :: funit
        integer :: io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("open_pft_or_ctf2_array_for_write: "//fname%to_char(),io_stat)
    end subroutine open_pft_or_ctf2_array_for_write

end submodule simple_polarft_ops_io
