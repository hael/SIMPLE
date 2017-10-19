module simple_stktab_handler
use simple_defs
use simple_fileio
use simple_imghead, only: find_ldim_nptcls

implicit none

public :: stktab_handler, test_stktab_handler
private

type :: mic_stk
    character(len=:), allocatable :: stkname
    integer :: fromp=0, top=0
end type mic_stk

type :: ptcl_entry
    type(mic_stk), pointer :: msp => null()
end type ptcl_entry

type :: stktab_handler
    private
    integer                       :: nmics   = 0
    integer                       :: nptcls  = 0
    integer                       :: ldim(3) = [0,0,0]
    type(ptcl_entry), allocatable :: pens(:)
    type(mic_stk),    allocatable :: mics(:)
    logical                       :: exists=.false.
  contains
    ! constructor
    procedure :: new
    ! getters
    procedure :: get_nmics
    procedure :: get_nptcls
    procedure :: get_ldim
    procedure :: get_stkname
    procedure :: get_stkname_and_ind
    ! setters
    procedure :: add_scale_tag
    procedure :: del_scale_tag
    ! I/O
    procedure :: write_stktab
    procedure :: del_stktab_files
    ! destructor
    procedure :: kill
end type stktab_handler

contains

    subroutine new( self, filetabname, part )
        class(stktab_handler), target, intent(inout) :: self
        character(len=*),              intent(in)    :: filetabname
        integer, optional,             intent(in)    :: part
        character(len=STDLEN), allocatable :: micnames(:)
        integer :: imic, ldim(3), pind_cnt, istart, istop, nptcls, iptcl, fnr, ppart
        ppart = 1
        if( present(part) ) ppart = part
        ! take care of micrograph stacks
        call read_filetable(filetabname, micnames)
        self%nmics = size(micnames)
        allocate(self%mics(self%nmics))
        istart = 1
        istop  = 0
        if( ppart == 1 ) call fopen(fnr, FILE='stktab_info.txt', STATUS='REPLACE', action='WRITE')
        if( ppart == 1 ) write(fnr,'(a)') '#MIC STKNAME #NPTCLS IN STK #FROMP #TOP'
        do imic=1,self%nmics
            allocate(self%mics(imic)%stkname, source=trim(micnames(imic)))
            call find_ldim_nptcls(trim(micnames(imic)), ldim, nptcls)
            ldim(3) = 1
            if( imic == 1 )then
                self%ldim = ldim
            else
                if( .not. all(self%ldim == ldim) )then
                    print *, 'micrograph stack #:   ', imic
                    print *, 'ldim in object:       ', self%ldim
                    print *, 'ldim read from stack: ', ldim
                    stop 'inconsistent logical dimensions; stktab_handler :: new'
                endif
            endif
            istop  = istop + nptcls
            self%mics(imic)%fromp = istart
            self%mics(imic)%top   = istop
            istart = istart + nptcls
            if( ppart == 1 ) write(fnr,'(a,1x,i6,1x,i10,1x,i10)') self%mics(imic)%stkname, nptcls, self%mics(imic)%fromp, self%mics(imic)%top
        end do
        if( ppart == 1 ) call fclose(fnr)
        self%nptcls = self%mics(self%nmics)%top
        ! take care of particle entries
        allocate(self%pens(self%nptcls))
        pind_cnt = 0
        do imic=1,self%nmics
            do iptcl=self%mics(imic)%fromp,self%mics(imic)%top
                pind_cnt = pind_cnt + 1
                self%pens(pind_cnt)%msp => self%mics(imic)
            end do
        end do
        self%exists = .true.
    end subroutine new

    integer function get_nmics( self )
        class(stktab_handler), intent(in) :: self
        get_nmics = self%nmics
    end function get_nmics

    integer function get_nptcls( self )
        class(stktab_handler), intent(in) :: self
        get_nptcls = self%nptcls
    end function get_nptcls

    function get_ldim( self ) result( ldim )
        class(stktab_handler), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function get_ldim

    function get_stkname( self, imic ) result( stkname )
        class(stktab_handler), intent(in) :: self
        integer,               intent(in) :: imic
        character(len=:), allocatable :: stkname
        ! sanity check
        if( imic < 1 .or. imic > self%nmics )then
            print *, 'imic:       ', imic
            print *, 'self%nmics: ', self%nmics
            stop 'imic index out of range; simple_stktab_handler :: get_stkname'
        endif
        allocate(stkname, source=trim(self%mics(imic)%stkname))
    end function get_stkname

    subroutine get_stkname_and_ind( self, iptcl, stkname, ind )
        class(stktab_handler),         intent(in)  :: self
        integer,                       intent(in)  :: iptcl
        character(len=:), allocatable, intent(out) :: stkname
        integer,                       intent(out) :: ind
        integer :: ldim(3), mic_stk_ind
        real    :: smpd
        ! first sanity check
        if( iptcl < 1 .or. iptcl > self%nptcls )then
            print *, 'iptcl:       ', iptcl
            print *, 'self%nptcls: ', self%nptcls
            stop 'iptcl index out of overall range; simple_stktab_handler :: get_stkname_and_ind'
        endif
        ! second sanity check
        if( iptcl < self%pens(iptcl)%msp%fromp .or. iptcl > self%pens(iptcl)%msp%top )then
            print *, 'iptcl:             ', iptcl
            print *, 'prange for micstk: ', self%pens(iptcl)%msp%fromp, self%pens(iptcl)%msp%top
            stop 'iptcl index out of micstk range; simple_stktab_handler :: get_stkname_and_ind'
        endif
        ! calculate index in mic stk
        ind = iptcl - self%pens(iptcl)%msp%fromp + 1
        ! read image
        if( allocated(stkname) ) deallocate(stkname)
        allocate(stkname, source=trim(self%pens(iptcl)%msp%stkname))
    end subroutine get_stkname_and_ind

    subroutine add_scale_tag( self )
        class(stktab_handler), intent(inout) :: self
        character(len=:), allocatable :: ext, newname
        integer :: imic
        do imic=1,self%nmics
            ext     = fname2ext(trim(self%mics(imic)%stkname))
            newname = add2fbody(self%mics(imic)%stkname, '.'//ext, '_sc')
            deallocate(self%mics(imic)%stkname)
            allocate(self%mics(imic)%stkname, source=newname)
        end do
    end subroutine add_scale_tag

    subroutine del_scale_tag( self )
        class(stktab_handler), intent(inout) :: self
        character(len=:), allocatable :: ext, newname
        integer :: imic
        do imic=1,self%nmics
            ext     = fname2ext(trim(self%mics(imic)%stkname))
            newname = del_from_fbody(self%mics(imic)%stkname, '.'//ext, '_sc')
            deallocate(self%mics(imic)%stkname)
            allocate(self%mics(imic)%stkname, source=newname)
        end do
    end subroutine del_scale_tag

    subroutine write_stktab( self, tabname )
        class(stktab_handler), intent(in) :: self
        character(len=*),      intent(in) :: tabname
        integer :: fnr, imic
        call fopen(fnr, file=trim(tabname), status='replace', action='write')
        do imic=1,self%nmics
            write(fnr,'(a)') self%mics(imic)%stkname
        end do
        call fclose(fnr)
    end subroutine write_stktab

    subroutine del_stktab_files( self )
        class(stktab_handler), intent(in) :: self
        integer :: imic
        do imic=1,self%nmics
            call del_file(self%mics(imic)%stkname)
        end do
    end subroutine del_stktab_files

    subroutine kill( self )
        class(stktab_handler), intent(inout) :: self
        integer :: imic, iptcl
        if( self%exists )then
            do imic=1,self%nmics
                deallocate(self%mics(imic)%stkname)
            end do
            do iptcl=1,self%nptcls
                self%pens(iptcl)%msp => null()
            end do
            self%nmics  = 0
            self%nptcls = 0
            self%ldim   = 0
            deallocate(self%pens, self%mics)
            self%exists = .false.
        endif
    end subroutine kill

    subroutine test_stktab_handler
        character(len=STDLEN), allocatable :: names(:)
        type(stktab_handler) :: stkhandle
        integer, parameter   :: NMICS=10
        integer :: imic
        names = make_filenames('stack_mic', NMICS, '.mrc')
        allocate(stkhandle%mics(NMICS))
        stkhandle%nmics = NMICS
        do imic=1,NMICS
            allocate(stkhandle%mics(imic)%stkname, source=trim(names(imic)))
            print *, stkhandle%mics(imic)%stkname
        end do
        print *, '*******************'
        call stkhandle%add_scale_tag
        do imic=1,NMICS
            print *, stkhandle%mics(imic)%stkname
        end do
        print *, '*******************'
        call stkhandle%del_scale_tag
        do imic=1,NMICS
            print *, stkhandle%mics(imic)%stkname
        end do
    end subroutine test_stktab_handler

end module simple_stktab_handler
