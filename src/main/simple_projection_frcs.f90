module simple_projection_frcs
#include "simple_lib.f08"
use simple_estimate_ssnr, only: resample_filter
use simple_oris,          only: oris
implicit none

public :: projection_frcs
private

type projection_frcs
    private
    integer           :: nprojs       = 0
    integer           :: filtsz       = 0
    integer           :: box4frc_calc = 0
    integer           :: nstates      = 1
    integer           :: headsz       = 0
    integer           :: hpind_frc    = 0
    real              :: file_header(4)
    real              :: smpd         = 0.0
    real              :: dstep        = 0.0
    real, allocatable :: res4frc_calc(:)
    real, allocatable :: frcs(:,:,:)
    logical           :: phaseplate   = .false.
    logical           :: exists       = .false.
contains
    ! constructor
    procedure          :: new
    ! exception
    procedure, private :: raise_exception
    procedure, private :: bound_res
    ! setters/getters
    procedure          :: get_nprojs
    procedure          :: get_filtsz
    procedure          :: set_frc
    procedure          :: get_frc
    procedure          :: frc_getter
    procedure          :: estimate_res
    procedure          :: estimate_find_for_eoavg
    procedure          :: estimate_lp_for_align
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

    subroutine bound_res( self, frc, res_frc05, res_frc0143 )
        class(projection_frcs), intent(in)    :: self
        real,                   intent(in)    :: frc(:)
        real,                   intent(inout) :: res_frc05, res_frc0143
        if( all(frc >  0.5) )then
            res_frc05   = 2. * self%smpd
            res_frc0143 = res_frc05
            return
        endif
        if( all(frc >  0.143) )then
            res_frc0143 = 2. * self%smpd
            return
        endif
        if( all(frc <=  0.143) )then
            res_frc05   = self%dstep
            res_frc0143 = res_frc05
            return
        endif
        if( all(frc <=  0.5) )then
            res_frc0143 = self%dstep
            return
        endif
    end subroutine bound_res

    ! setters/getters

    pure integer function get_nprojs( self )
        class(projection_frcs), intent(in) :: self
        get_nprojs = self%nprojs
    end function get_nprojs

    pure integer function get_filtsz( self )
        class(projection_frcs), intent(in) :: self
        get_filtsz = self%filtsz
    end function get_filtsz

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

    function get_frc( self, proj, box, hpind_fsc, phaseplate, state ) result( frc )
        class(projection_frcs), intent(in) :: self
        integer,                intent(in) :: proj, box, hpind_fsc
        logical,                intent(in) :: phaseplate
        integer, optional,      intent(in) :: state
        real, allocatable :: frc(:)
        real, allocatable :: res(:)
        real    :: frcmax
        integer :: sstate, find_plate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( proj, sstate, 'ERROR, out of bounds in get_frc')
        if( box /= self%box4frc_calc )then
            res = get_resarr(box, self%smpd)
            frc = resample_filter(self%frcs(sstate,proj,:), self%res4frc_calc, res)
        else
            allocate(frc(self%filtsz), source=self%frcs(sstate,proj,:))
        endif
        if( phaseplate ) call phaseplate_correct_fsc(frc, find_plate)
        if( hpind_fsc > 0 )then
            frcmax = maxval(frc)
            frc(:hpind_fsc) = frcmax
        endif
    end function get_frc

    subroutine frc_getter( self, proj, hpind_fsc, phaseplate, frc, state )
        class(projection_frcs), intent(in)  :: self
        integer,                intent(in)  :: proj, hpind_fsc
        logical,                intent(in)  :: phaseplate
        real,                   intent(out) :: frc(self%filtsz)
        integer, optional,      intent(in)  :: state
        real    :: frcmax
        integer :: sstate, find_plate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( proj, sstate, 'ERROR, out of bounds in get_frc')
        frc = self%frcs(sstate,proj,:)
        if( phaseplate ) call phaseplate_correct_fsc(frc, find_plate)
        if( hpind_fsc > 0 )then
            frcmax = maxval(frc)
            frc(:hpind_fsc) = frcmax
        endif
    end subroutine frc_getter

    subroutine estimate_res( self, proj, res_frc05, res_frc0143, state )
        class(projection_frcs), intent(in)  :: self
        integer,                intent(in)  :: proj
        real,                   intent(out) :: res_frc05, res_frc0143
        integer, optional,      intent(in)  :: state
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( proj, sstate, 'ERROR, out of bounds in estimate_res')
        call get_resolution(self%frcs(sstate,proj,:), self%res4frc_calc, res_frc05, res_frc0143)
        call self%bound_res(self%frcs(sstate,proj,:), res_frc05, res_frc0143)
    end subroutine estimate_res

    function estimate_find_for_eoavg( self, proj, state ) result( find )
        class(projection_frcs), intent(in)  :: self
        integer,                intent(in)  :: proj
        integer, optional,      intent(in)  :: state
        integer :: sstate, find
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( proj, sstate, 'ERROR, out of bounds in estimate_find_for_eoavg')
        find = max(K4EOAVGLB,get_lplim_at_corr(self%frcs(sstate,proj,:), FSC4EOAVG2D))
    end function estimate_find_for_eoavg

    function estimate_lp_for_align( self, state ) result( lp )
        class(projection_frcs), intent(in)  :: self
        integer, optional,      intent(in)  :: state
        real    :: lplims(self%nprojs)
        integer :: alloc_stat, sstate, iproj
        real    :: lp, res_frc05, res_frc0143
        sstate = 1
        if( present(state) ) sstate = state
        ! order FRCs according to low-pass limit (best first)
        do iproj=1,self%nprojs
            call get_resolution(self%frcs(sstate,iproj,:), self%res4frc_calc, res_frc05, res_frc0143)
            call self%bound_res(self%frcs(sstate,iproj,:), res_frc05, res_frc0143)
            lplims(iproj) = res_frc0143
        end do
        call hpsort(lplims)
        ! return median of top three clusters
        lp = median(lplims(1:3))
    end function estimate_lp_for_align

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
