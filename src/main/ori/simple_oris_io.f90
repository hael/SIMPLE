submodule (simple_oris) simple_oris_io
use simple_ori_api
use simple_ori, only: ori
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine print( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        call self%o(i)%print_ori()
    end subroutine print

    module subroutine print_matrices( self )
        class(oris), intent(inout) :: self
        integer :: i
        write(logfhandle,*) 'ORDER OF ROTATION MATRIX ELEMENTS: (1,1) (1,2) (1,3) (2,1) (2,2) (2,3) (3,1) (3,2) (3,3)'
        do i=1,self%n
            call self%o(i)%print_mat()
        end do
    end subroutine print_matrices

    module subroutine read( self, orifile, fromto, nst )
        class(oris),       intent(inout) :: self
        class(string),     intent(in)    :: orifile
        integer, optional, intent(in)    :: fromto(2)
        integer, optional, intent(out)   :: nst
        character(len=100) :: io_message
        integer :: file_stat, i, fnr, state, istart, iend
        if( .not. file_exists(orifile) )then
            THROW_HARD("the file you are trying to read: "//orifile%to_char()//' does not exist in cwd' )
        endif
        if( fname2ext(orifile) == 'bin' )then
            THROW_HARD('this method does not support binary files; read')
        endif
        io_message='No error'
        call fopen(fnr, FILE=orifile, STATUS='OLD', action='READ', iostat=file_stat,iomsg=io_message)
        call fileiochk("oris ; read ,Error when opening file for reading: "//orifile%to_char()//':'//trim(io_message), file_stat)
        if( present(nst) ) nst = 0
        if( present(fromto) )then
            istart = fromto(1)
            iend   = fromto(2)
            if(istart <      1) THROW_HARD('Invalid index; read')
            if(iend   > self%n) THROW_HARD('Invalid index; read')
        else
            istart = 1
            iend   = self%n
        endif
        do i = istart, iend
            call self%o(i)%read(fnr)
            if( present(nst) )then
                state = self%o(i)%get_state()
                nst   = max(1,max(state,nst))
            endif
        end do
        call fclose(fnr)
    end subroutine read

    module subroutine read_ctfparams_state_eo( self, ctfparamfile )
        class(oris),    intent(inout) :: self
        class(string),  intent(in)    :: ctfparamfile
        logical    :: params_are_there(10)
        integer    :: i
        type(oris) :: os_tmp
        if( .not. file_exists(ctfparamfile) )then
            THROW_HARD ("read_ctfparams_state_eo; The file you are trying to read: "//ctfparamfile%to_char()//' does not exist')
        endif
        if( ctfparamfile%has_substr('.bin') )then
            THROW_HARD('this method does not support binary files; read_ctfparams_state_eo')
        endif
        call os_tmp%new(self%n, self%o(1)%is_particle())
        call os_tmp%read(ctfparamfile)
        params_are_there(1)  = os_tmp%isthere('smpd')
        params_are_there(2)  = os_tmp%isthere('kv')
        params_are_there(3)  = os_tmp%isthere('cs')
        params_are_there(4)  = os_tmp%isthere('fraca')
        params_are_there(5)  = os_tmp%isthere('phshift')
        params_are_there(6)  = os_tmp%isthere('dfx')
        params_are_there(7)  = os_tmp%isthere('dfy')
        params_are_there(8)  = os_tmp%isthere('angast')
        params_are_there(9)  = os_tmp%isthere('state')
        params_are_there(10) = os_tmp%isthere('eo')
        do i=1,self%n
            if( params_are_there(1)  ) call self%set(i, 'smpd',    os_tmp%get(i, 'smpd')   )
            if( params_are_there(2)  ) call self%set(i, 'kv',      os_tmp%get(i, 'kv')     )
            if( params_are_there(3)  ) call self%set(i, 'cs',      os_tmp%get(i, 'cs')     )
            if( params_are_there(4)  ) call self%set(i, 'fraca',   os_tmp%get(i, 'fraca')  )
            if( params_are_there(5)  ) call self%set(i, 'phshift', os_tmp%get(i, 'phshift'))
            if( params_are_there(6)  ) call self%set_dfx(i,        os_tmp%get_dfx(i)       )
            if( params_are_there(7)  ) call self%set_dfy(i,        os_tmp%get_dfy(i)       )
            if( params_are_there(8)  ) call self%set(i, 'angast',  os_tmp%get(i, 'angast') )
            if( params_are_there(9)  ) call self%set(i, 'state',   os_tmp%get(i, 'state')  )
            if( params_are_there(10) ) call self%set(i, 'eo',      os_tmp%get(i, 'eo')     )
        end do
        call os_tmp%kill
    end subroutine read_ctfparams_state_eo

    module subroutine write_1( self, orifile, fromto )
        class(oris),       intent(in) :: self
        class(string),     intent(in) :: orifile
        integer, optional, intent(in) :: fromto(2)
        character(len=100) :: io_message
        integer            :: file_stat, fnr, i, ffromto(2), cnt
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        call fopen(fnr, orifile, status='REPLACE', action='WRITE', iostat=file_stat, iomsg=io_message)
        call fileiochk(' Error opening file for writing: '//orifile%to_char()//' ; '//trim(io_message), file_stat)
        cnt = 0
        do i=ffromto(1),ffromto(2)
            cnt = cnt + 1
            call self%o(i)%write(fnr)
        end do
        call fclose(fnr)
    end subroutine write_1

    module subroutine write_2( self, i, orifile )
        class(oris),   intent(inout) :: self
        class(string), intent(in)    :: orifile
        integer,       intent(in)    :: i
        integer :: fnr, file_stat
        call fopen(fnr, orifile, status='UNKNOWN', action='WRITE', position='APPEND', iostat=file_stat)
        call fileiochk( 'In: write_2, module: simple_oris.f90  opening '//orifile%to_char(), file_stat )
        call self%o(i)%write(fnr)
        call fclose(fnr)
    end subroutine write_2

    module subroutine write2bild( self, file )
        class(oris),   intent(inout) :: self
        class(string), intent(in)    :: file
        integer :: i,funit, file_stat
        call fopen(funit, file, status='REPLACE', action='WRITE',iostat=file_stat)
        call fileiochk( 'In: write2bild, module: simple_oris.f90  opening '//file%to_char(), file_stat )
        ! header
        write(funit,'(A)')".translate 0.0 0.0 0.0"
        write(funit,'(A)')".scale 10"
        write(funit,'(A)')".comment -- unit sphere --"
        write(funit,'(A)')".color 0.8 0.8 0.8"
        write(funit,'(A)')".sphere 0 0 0 1.0"
        write(funit,'(A)')".comment -- planes --"
        write(funit,'(A)')".color 0.3 0.3 0.3"
        write(funit,'(A)')".cylinder -0.02 0 0 0.02 0 0 1.02"
        write(funit,'(A)')".cylinder 0 -0.02 0 0 0.02 0 1.02"
        write(funit,'(A)')".cylinder 0 0 -0.02 0 0 0.02 1.02"
        write(funit,'(A)')".comment -- x-axis --"
        write(funit,'(A)')".color 1 0 0"
        write(funit,'(A)')".cylinder -1.5 0 0 1.5 0 0 0.02"
        write(funit,'(A)')".comment -- y-axis --"
        write(funit,'(A)')".color 0 1 0"
        write(funit,'(A)')".cylinder 0 -1.5 0 0 1.5 0 0.02"
        write(funit,'(A)')".comment -- z-axis --"
        write(funit,'(A)')".color 0 0 1"
        write(funit,'(A)')".cylinder 0 0 -1.5 0 0 1.5 0.02"
        write(funit,'(A)')".comment -- north pole --"
        write(funit,'(A)')".color 0 0 1"
        write(funit,'(A)')".sphere 0 0 1.5 0.1"
        ! body
        write(funit,'(A)')".color 0.4 0.4 0.4"
        do i=1,self%n
            call self%o(i)%write2bild(funit)
        enddo
        call fclose(funit)
    end subroutine write2bild

end submodule simple_oris_io
