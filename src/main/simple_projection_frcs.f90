module simple_projection_frcs
use simple_filterer, only: resample_filter
use simple_oris,     only: oris
implicit none

type projection_frcs
    private
    integer           :: nprojs     = 0
    integer           :: filtsz     = 0 
    integer           :: box_target = 0
    integer           :: nstates    = 0
    real              :: smpd       = 0.0
    real, allocatable :: res_target(:)
    real, allocatable :: frcs(:,:,:)
    type(oris)        :: e_space
    logical           :: exists = .false.
contains
    ! constructor
    procedure :: new
    ! exception
    procedure :: raise_exception
    ! setters/getters
    procedure :: set_frc
    procedure :: get_frc
    procedure :: get_ssnr
    procedure :: estimate_res
    ! I/O
    procedure :: read
    procedure :: write
    ! destructor
    procedure :: kill
end type projection_frcs

contains

    ! constructor

    subroutine new( self, e_space, box_target, smpd, nstates )
        use simple_syslib, only: alloc_errchk
        use simple_math,   only: fdim, get_resarr
        class(projection_frcs), intent(inout) :: self
        class(oris),            intent(inout) :: e_space 
        integer,                intent(in)    :: box_target
        real,                   intent(in)    :: smpd
        integer, optional,      intent(in)    :: nstates
        integer :: alloc_stat
        call self%kill
        self%e_space    = e_space
        self%box_target = box_target
        self%smpd       = smpd
        self%nprojs     = self%e_space%get_noris()
        self%filtsz     = fdim(box_target)-1
        self%nstates    = 1
        self%res_target = get_resarr(self%box_target, self%smpd)
        if( present(nstates) ) self%nstates = nstates
        allocate( self%frcs(self%nstates,self%nprojs,self%filtsz), stat=alloc_stat)
        call alloc_errchk('new; simple_projection_frcs', alloc_stat)
        self%frcs   = 1.0
        self%exists = .true.
    end subroutine new

    ! exception

    subroutine raise_exception( self, state, proj, msg )
        class(projection_frcs), intent(in) :: self
        integer,                intent(in) :: state, proj
        character(len=*),       intent(in) :: msg
        if( state < 1 .or. state > self%nstates )  write(*,*) 'state: ', state
        if( proj  < 1 .or. proj  > self%nprojs  )  write(*,*) 'proj: ', proj
        write(*,'(a)') msg
        stop 'simple_projection_frcs :: raise_exception'
    end subroutine raise_exception

    ! setters/getters

    subroutine set_frc( self, box, state, proj, frc )
        use simple_math, only: get_resarr
        class(projection_frcs), intent(inout) :: self
        integer,                intent(in)    :: box, state, proj
        real,                   intent(in)    :: frc(:)
        real, allocatable :: res(:)
        call self%raise_exception( state, proj, 'ERROR, out of bounds in set_frc')
        if( box /= self%box_target )then
            res = get_resarr(box, self%smpd)
            self%frcs(state,proj,:) = resample_filter(frc, res, self%res_target) 
        else
            self%frcs(state,proj,:) = frc
        endif
    end subroutine set_frc

    function get_frc( self, state, proj ) result( frc )
        class(projection_frcs), intent(in) :: self
        integer,                intent(in) :: state, proj
        real, allocatable :: frc(:)
        call self%raise_exception( state, proj, 'ERROR, out of bounds in get_frc')
        allocate(frc(self%filtsz), source=self%frcs(state,proj,:))
    end function get_frc

    function get_ssnr( self, state, proj ) result( ssnr )
        use simple_estimate_ssnr, only: fsc2ssnr
        class(projection_frcs), intent(in) :: self
        integer,                intent(in) :: state, proj
        real, allocatable :: ssnr(:)
        call self%raise_exception( state, proj, 'ERROR, out of bounds in get_ssnr')
        ssnr = self%frcs(state,proj,:)
    end function get_ssnr

    subroutine estimate_res( self, state, proj, frc05, frc0143 )
        use simple_math, only: get_resolution
        class(projection_frcs), intent(in)  :: self
        integer,                intent(in)  :: state, proj
        real,                   intent(out) :: frc05, frc0143
        call self%raise_exception( state, proj, 'ERROR, out of bounds in estimate_res')
        call get_resolution(self%frcs(state,proj,:), self%res_target, frc05, frc0143 )
    end subroutine estimate_res

    ! I/O

    subroutine read( self, fname )
        use simple_fileio, only: fopen, fclose
        class(projection_frcs), intent(inout) :: self
        character(len=*),       intent(in)    :: fname
        integer          :: funit, io_stat
        character(len=7) :: stat_str
        if(.not.fopen(funit,fname,access='STREAM',action='READ',status='OLD', iostat=io_stat))&
        &call fileio_errmsg('projection_frcs; read; open for read '// trim(fname), io_stat)
        read(unit=funit,pos=1,iostat=io_stat) self%frcs
        call fileio_errmsg('projection_frcs; read; actual read', io_stat)
        if(.not.fclose(funit, iostat=io_stat))&
        &call fileio_errmsg('projection_frcs; read; fhandle cose', io_stat)
    end subroutine read

    subroutine write( self, fname )
        use simple_fileio, only: fopen, fclose
        class(projection_frcs), intent(in) :: self
        character(len=*),       intent(in) :: fname
        integer          :: funit, io_stat
        character(len=7) :: stat_str
        if(.not.fopen(funit,fname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat))&
        &call fileio_errmsg('projection_frcs; write; open for write '// trim(fname), io_stat)
        write(unit=funit,pos=1,iostat=io_stat) self%frcs
        call fileio_errmsg('projection_frcs; write; actual write', io_stat)
        if(.not.fclose(funit, iostat=io_stat))&
        &call fileio_errmsg('projection_frcs; write; fhandle cose', io_stat)
    end subroutine write

    ! destructor

    subroutine kill( self )
        class(projection_frcs), intent(inout) :: self
        if( self%exists )then
            deallocate(self%res_target, self%frcs)
            call self%e_space%kill
            self%nprojs     = 0 
            self%filtsz     = 0
            self%box_target = 0
            self%exists     = .false.
        endif
    end subroutine kill

end module simple_projection_frcs
