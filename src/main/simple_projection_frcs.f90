module simple_projection_frcs
#include "simple_lib.f08"
use simple_estimate_ssnr, only: resample_filter
use simple_oris,          only: oris

implicit none

type projection_frcs
    private
    integer           :: nprojs       = 0
    integer           :: filtsz       = 0
    integer           :: box4frc_calc = 0
    integer           :: nstates      = 1
    integer           :: headsz       = 0
    real              :: file_header(4)
    real              :: smpd         = 0.0
    real              :: dstep        = 0.0
    real, allocatable :: res4frc_calc(:)
    real, allocatable :: frcs(:,:,:)
    logical           :: exists       = .false.
contains
    ! constructor
    procedure          :: new
    ! exception
    procedure, private :: raise_exception
    procedure, private :: bound_res
    ! setters/getters
    procedure          :: get_nprojs
    procedure          :: set_frc
    procedure          :: get_frc
    procedure, private :: estimate_res_1
    procedure, private :: estimate_res_2
    procedure, private :: estimate_res_3
    generic            :: estimate_res => estimate_res_1, estimate_res_2, estimate_res_3
    ! I/O
    procedure          :: read
    procedure          :: write
    ! destructor
    procedure          :: kill
end type projection_frcs

contains

    ! constructor

    subroutine new( self, nprojs, box4frc_calc, smpd, nstates )
        class(projection_frcs), intent(inout) :: self
        integer,                intent(in)    :: nprojs
        integer,                intent(in)    :: box4frc_calc
        real,                   intent(in)    :: smpd
        integer, optional,      intent(in)    :: nstates
        integer :: alloc_stat
        ! init
        call self%kill
        self%nprojs       = nprojs
        self%box4frc_calc = box4frc_calc
        self%smpd         = smpd
        self%filtsz       = fdim(box4frc_calc) - 1
        self%dstep        = real(box4frc_calc - 1) * smpd ! first wavelength of FT
        self%res4frc_calc = get_resarr(box4frc_calc, self%smpd)
        self%nstates      = 1
        if( present(nstates) ) self%nstates = nstates
        ! prep file header
        self%file_header(1) = real(self%nprojs)
        self%file_header(2) = real(self%box4frc_calc)
        self%file_header(3) = self%smpd
        self%file_header(4) = real(self%nstates)
        self%headsz         = sizeof(self%file_header)
        ! alloc
        allocate( self%frcs(self%nstates,self%nprojs,self%filtsz), stat=alloc_stat)
        call alloc_errchk('new; simple_projection_frcs', alloc_stat)
        self%frcs   = 1.0
        self%exists = .true.
    end subroutine new

    ! exception

    subroutine raise_exception( self, proj, state, msg )
        class(projection_frcs), intent(in) :: self
        integer,                intent(in) :: proj, state
        character(len=*),       intent(in) :: msg
        logical :: l_outside
        l_outside = .false.
        if( proj  < 1 .or. proj  > self%nprojs  )then
            write(*,*) 'proj: ', proj
            l_outside = .true.
        endif
        if( state < 1 .or. state > self%nstates ) then
            write(*,*) 'state: ', state
            l_outside = .true.
        endif
        if( l_outside )then
            write(*,'(a)') msg
            stop 'simple_projection_frcs :: raise_exception'
        endif
    end subroutine raise_exception

    ! bound res

    subroutine bound_res( self, frc, frc05, frc0143 )
        class(projection_frcs), intent(in)    :: self
        real,                   intent(in)    :: frc(:)
        real,                   intent(inout) :: frc05, frc0143  
        if( all(frc >  0.5) )then
            frc05   = 2. * self%smpd
            frc0143 = frc05
            return
        endif
        if( all(frc >  0.143) )then
            frc0143 = 2. * self%smpd
            return
        endif
        if( all(frc <=  0.143) )then
            frc05   = self%dstep
            frc0143 = frc05
            return
        endif
        if( all(frc <=  0.5) )then
            frc0143 = self%dstep
            return
        endif
    end subroutine bound_res

    ! setters/getters

    integer function get_nprojs( self )
        class(projection_frcs), intent(in) :: self
        get_nprojs = self%nprojs
    end function get_nprojs

    subroutine set_frc( self, proj, frc, state )
        class(projection_frcs), intent(inout) :: self
        integer,                intent(in)    :: proj
        real,                   intent(in)    :: frc(:)
        integer, optional,      intent(in)    :: state
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( proj, sstate, 'ERROR, out of bounds in set_frc')
        if( size(frc) /= self%filtsz )then
            stop 'size of input frc not conforming; simple_projection_frcs :: set_frc'
        else
            self%frcs(sstate,proj,:) = frc
        endif
    end subroutine set_frc

    function get_frc( self, proj, box, state ) result( frc )
        class(projection_frcs), intent(in) :: self
        integer,                intent(in) :: proj, box
        integer, optional,      intent(in) :: state
        real, allocatable :: frc(:)
        real, allocatable :: res(:)
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( proj, sstate, 'ERROR, out of bounds in get_frc')
        if( box /= self%box4frc_calc )then
            res = get_resarr(box, self%smpd)
            frc = resample_filter(self%frcs(sstate,proj,:), self%res4frc_calc, res)
        else
            allocate(frc(self%filtsz), source=self%frcs(sstate,proj,:))
        endif
    end function get_frc

    subroutine estimate_res_1( self, proj, frc05, frc0143, state )
        class(projection_frcs), intent(in)  :: self
        integer,                intent(in)  :: proj
        real,                   intent(out) :: frc05, frc0143
        integer, optional,      intent(in)  :: state
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( proj, sstate, 'ERROR, out of bounds in estimate_res')
        call get_resolution(self%frcs(sstate,proj,:), self%res4frc_calc, frc05, frc0143)
        call self%bound_res(self%frcs(sstate,proj,:), frc05, frc0143)
    end subroutine estimate_res_1

    subroutine estimate_res_2( self, state )
        class(projection_frcs), intent(in)  :: self
        integer, optional,      intent(in)  :: state
        integer :: alloc_stat, sstate, j
        real    :: frc05, frc0143
        real, allocatable :: frc(:)
        allocate( frc(self%filtsz), stat=alloc_stat )
        call alloc_errchk( 'estimate_res_2; simple_projection_frcs', alloc_stat )
        sstate = 1
        if( present(state) ) sstate = state
        frc = sum(self%frcs(sstate,:,:),dim=1) / real(self%nprojs)
        call get_resolution(frc, self%res4frc_calc, frc05, frc0143)
        do j=1,size(self%res4frc_calc)
            write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', self%res4frc_calc(j), '>>> CORRELATION:', frc(j)
        end do
        call self%bound_res(frc, frc05, frc0143)
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FRC=0.500 DETERMINED TO:', frc05
        write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FRC=0.143 DETERMINED TO:', frc0143 
    end subroutine estimate_res_2

    subroutine estimate_res_3( self, nbest, frc, res, frc05, frc0143, state )
        use simple_math, only: get_resolution, hpsort, get_lplim
        class(projection_frcs), intent(in)  :: self
        integer,                intent(in)  :: nbest
        real, allocatable,      intent(out) :: frc(:), res(:)
        real,                   intent(out) :: frc05, frc0143
        integer, optional,      intent(in)  :: state
        integer :: alloc_stat, sstate, iproj
        real,    allocatable :: lplims(:)
        integer, allocatable :: order(:)
        if( allocated(frc) ) deallocate(frc)
        if( allocated(res) ) deallocate(res)
        allocate( frc(self%filtsz), res(self%filtsz),&
            &lplims(self%nprojs), order(self%nprojs), stat=alloc_stat )
        allocchk( 'estimate_res_2; simple_projection_frcs allocating frc, res, order and lplims' )
        sstate = 1
        if( present(state) ) sstate = state
        ! order FRCs according to low-pass limit (best first)
        do iproj=1,self%nprojs
            lplims(iproj) = get_lplim(self%frcs(sstate,iproj,:))
        end do
        order = (/(iproj,iproj=1,self%nprojs)/)
        call hpsort(self%nprojs, lplims, order)
        ! average the nbest ones
        frc = 0.0
        do iproj=1,nbest
            frc = frc + self%frcs(sstate,order(iproj),:)
        end do
        ! prep output
        frc = frc / real(nbest)
        call get_resolution(frc, self%res4frc_calc, frc05, frc0143)
        call self%bound_res(frc, frc05, frc0143)
        res = self%res4frc_calc
    end subroutine estimate_res_3

    ! I/O

    subroutine read( self, fname )
        class(projection_frcs), intent(inout) :: self
        character(len=*),       intent(in)    :: fname
        integer          :: funit, io_stat
        character(len=7) :: stat_str
        call fopen(funit,fname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
        call fileio_errmsg('projection_frcs; read; open for read '//trim(fname), io_stat)
        read(unit=funit,pos=1) self%file_header
        ! re-create the object according to file_header info
        call self%new(nint(self%file_header(1)), nint(self%file_header(2)), self%file_header(3), nint(self%file_header(4)))
        read(unit=funit,pos=self%headsz + 1) self%frcs
        call fclose(funit, errmsg='projection_frcs; read; fhandle cose')
    end subroutine read

    subroutine write( self, fname )
        class(projection_frcs), intent(in) :: self
        character(len=*),       intent(in) :: fname
        integer          :: funit, io_stat
        character(len=7) :: stat_str
        call fopen(funit,fname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        call fileio_errmsg('projection_frcs; write; open for write '//trim(fname), io_stat)
        write(unit=funit,pos=1) self%file_header
        write(unit=funit,pos=self%headsz + 1) self%frcs
        call fclose(funit, errmsg='projection_frcs; write; fhandle cose')
    end subroutine write

    ! destructor

    subroutine kill( self )
        class(projection_frcs), intent(inout) :: self
        if( self%exists )then
            deallocate(self%res4frc_calc, self%frcs)
            self%nprojs       = 0 
            self%filtsz       = 0
            self%box4frc_calc = 0
            self%exists       = .false.
        endif
    end subroutine kill

end module simple_projection_frcs
