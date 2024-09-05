module simple_class_frcs
include 'simple_lib.f08'
use simple_fsc
implicit none

public :: class_frcs
private
#include "simple_local_flags.inc"

type class_frcs
    private 
    integer           :: ncls         = 0
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
    procedure          :: get_ncls
    procedure          :: get_filtsz
    procedure          :: set_frc
    procedure          :: get_frc
    procedure          :: frc_getter
    procedure          :: avg_frc_getter
    procedure          :: getter
    procedure          :: estimate_res
    procedure          :: estimate_find_for_eoavg
    procedure          :: estimate_lp_for_align
    procedure          :: pad
    ! I/O
    procedure          :: read
    procedure          :: write
    procedure          :: print_frcs
    ! destructor
    procedure          :: kill
end type class_frcs

contains

    ! constructor

    subroutine new( self, ncls, box4frc_calc, smpd, nstates )
        class(class_frcs), intent(inout) :: self
        integer,           intent(in)    :: ncls
        integer,           intent(in)    :: box4frc_calc
        real,              intent(in)    :: smpd
        integer, optional, intent(in)    :: nstates
        ! init
        call self%kill
        self%ncls         = ncls
        self%box4frc_calc = box4frc_calc
        self%smpd         = smpd
        self%filtsz       = fdim(box4frc_calc) - 1
        self%dstep        = real(box4frc_calc - 1) * smpd ! first wavelength of FT
        self%res4frc_calc = get_resarr(box4frc_calc, self%smpd)
        self%nstates      = 1
        if( present(nstates) ) self%nstates = nstates
        ! prep file header
        self%file_header(1) = real(self%ncls)
        self%file_header(2) = real(self%box4frc_calc)
        self%file_header(3) = self%smpd
        self%file_header(4) = real(self%nstates)
        self%headsz         = sizeof(self%file_header)
        ! alloc
        allocate(self%frcs(self%nstates,self%ncls,self%filtsz), source=0.0)
        self%exists = .true.
    end subroutine new

    ! exception

    subroutine raise_exception( self, cls, state, msg )
        class(class_frcs), intent(in) :: self
        integer,           intent(in) :: cls, state
        character(len=*),  intent(in) :: msg
        logical :: l_outside
        l_outside = .false.
        if( cls  < 1 .or. cls  > self%ncls  )then
            write(logfhandle,*) self%ncls
            write(logfhandle,*) 'exists: ', self%exists
            write(logfhandle,*) 'cls: ', cls
            l_outside = .true.
        endif
        if( state < 1 .or. state > self%nstates ) then
            write(logfhandle,*) 'state: ', state
            l_outside = .true.
        endif
        if( l_outside ) THROW_HARD(trim(msg)//'; raise_exception')
    end subroutine raise_exception

    ! bound res

    subroutine bound_res( self, frc, res_frc05, res_frc0143 )
        class(class_frcs), intent(in)    :: self
        real,              intent(in)    :: frc(:)
        real,              intent(inout) :: res_frc05, res_frc0143
        if( res_frc05 < 0.0001 .and. res_frc0143 < 0.0001 )then
            res_frc0143 = self%res4frc_calc(1)
            res_frc05   = res_frc0143
            return
        endif
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

    pure integer function get_ncls( self )
        class(class_frcs), intent(in) :: self
        get_ncls = self%ncls
    end function get_ncls

    pure integer function get_filtsz( self )
        class(class_frcs), intent(in) :: self
        get_filtsz = self%filtsz
    end function get_filtsz

    subroutine set_frc( self, cls, frc, state )
        class(class_frcs), intent(inout) :: self
        integer,           intent(in)    :: cls
        real,              intent(in)    :: frc(:)
        integer, optional, intent(in)    :: state
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( cls, sstate, 'ERROR, out of bounds in set_frc')
        if( size(frc) /= self%filtsz )then
            THROW_HARD('size of input frc not conforming; set_frc')
        else
            self%frcs(sstate,cls,:) = frc
        endif
    end subroutine set_frc

    !> \brief  re-samples a filter array
    function resample_filter( filt_orig, res_orig, res_new ) result( filt_resamp )
        real, intent(in)  :: filt_orig(:), res_orig(:), res_new(:)
        real, allocatable :: filt_resamp(:) !< output filter array
        integer :: filtsz_orig, filtsz_resamp, k, ind
        real    :: dist
        filtsz_orig   = size(filt_orig)
        filtsz_resamp = size(res_new)
        allocate(filt_resamp(filtsz_resamp))
        do k=1,filtsz_resamp
            call find(res_orig, filtsz_orig, res_new(k), ind, dist)
            filt_resamp(k) = filt_orig(ind)
        end do
    end function resample_filter

    !>  getter for values interepreted as FRCs
    function get_frc( self, cls, box, state ) result( frc )
        class(class_frcs), intent(in) :: self
        integer,           intent(in) :: cls, box
        integer, optional, intent(in) :: state
        real, allocatable :: frc(:)
        real, allocatable :: res(:)
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( cls, sstate, 'ERROR, out of bounds in get_frc')
        if( box /= self%box4frc_calc )then
            res = get_resarr(box, self%smpd)
            frc = resample_filter(self%frcs(sstate,cls,:), self%res4frc_calc, res)
        else
            allocate(frc(self%filtsz), source=self%frcs(sstate,cls,:))
        endif
    end function get_frc

    !>  getter for values interepreted as FRCs
    subroutine frc_getter( self, cls, hpind_fsc, phaseplate, frc, state )
        class(class_frcs), intent(in)  :: self
        integer,           intent(in)  :: cls, hpind_fsc
        logical,           intent(in)  :: phaseplate
        real,              intent(out) :: frc(self%filtsz)
        integer, optional, intent(in)  :: state
        integer :: sstate, find_plate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( cls, sstate, 'ERROR, out of bounds in frc_getter')
        frc = self%frcs(sstate,cls,:)
        if( phaseplate )then
            if( any(frc > 0.5) )call phaseplate_correct_fsc(frc, find_plate)
        endif
        if( hpind_fsc > 0 ) frc(:hpind_fsc) = frc(hpind_fsc + 1)
    end subroutine frc_getter

    subroutine avg_frc_getter( self, frcs_avg, states, state )
        class(class_frcs), intent(in)  :: self
        real,              intent(out) :: frcs_avg(self%filtsz)
        integer,           intent(in)  :: states(self%ncls)
        integer, optional, intent(in)  :: state
        integer :: sstate, icls, k
        sstate = 1
        if( present(state) ) sstate = state
        do k = 1,self%filtsz
            frcs_avg(k) = sum(self%frcs(sstate,:,k), mask=states > 0 .and. self%frcs(sstate,:,k) > 0.)
        enddo
        frcs_avg = frcs_avg / real(count(states > 0))
    end subroutine avg_frc_getter

    !>  getter for raw values
    subroutine getter( self, cls, frc, state )
        class(class_frcs), intent(in)  :: self
        integer,           intent(in)  :: cls
        real,              intent(out) :: frc(self%filtsz)
        integer, optional, intent(in)  :: state
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( cls, sstate, 'ERROR, out of bounds in frc_getter')
        frc = self%frcs(sstate,cls,:)
    end subroutine getter

    subroutine estimate_res( self, cls, res_frc05, res_frc0143, state )
        class(class_frcs), intent(in)  :: self
        integer,           intent(in)  :: cls
        real,              intent(out) :: res_frc05, res_frc0143
        integer, optional, intent(in)  :: state
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( cls, sstate, 'ERROR, out of bounds in estimate_res')
        call get_resolution(self%frcs(sstate,cls,:), self%res4frc_calc, res_frc05, res_frc0143)
        call self%bound_res(self%frcs(sstate,cls,:), res_frc05, res_frc0143)
    end subroutine estimate_res

    function estimate_find_for_eoavg( self, cls, state ) result( find )
        class(class_frcs), intent(in)  :: self
        integer,           intent(in)  :: cls
        integer, optional, intent(in)  :: state
        integer :: sstate, find
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( cls, sstate, 'ERROR, out of bounds in estimate_find_for_eoavg')
        find = max(K4EOAVGLB,get_find_at_corr(self%frcs(sstate,cls,:), FSC4EOAVG2D))
    end function estimate_find_for_eoavg

    function estimate_lp_for_align( self, state ) result( lp )
        class(class_frcs), intent(in)  :: self
        integer, optional, intent(in)  :: state
        real    :: lplims(self%ncls),lp3(3)
        integer :: sstate, icls
        real    :: lp, res_frc05, res_frc0143
        sstate = 1
        if( present(state) ) sstate = state
        ! order FRCs according to low-pass limit (best first)
        do icls=1,self%ncls
            call get_resolution(self%frcs(sstate,icls,:), self%res4frc_calc, res_frc05, res_frc0143)
            call self%bound_res(self%frcs(sstate,icls,:), res_frc05, res_frc0143)
            lplims(icls) = res_frc0143
        end do
        lp3 = min3(lplims)
        ! return median of top three clusters
        lp = median(lp3)
    end function estimate_lp_for_align

    subroutine pad( self, newsmpd, newbox, self_out )
        class(class_frcs), intent(in)  :: self
        real,              intent(in)  :: newsmpd
        integer,           intent(in)  :: newbox
        type(class_frcs),  intent(out) :: self_out
        if( newbox < self%box4frc_calc )then
            THROW_HARD('New <= old filter size; downsample')
        else if( newbox == self%box4frc_calc )then
            ! copy
            call self_out%new(self%ncls, newbox, newsmpd, self%nstates)
            self_out%frcs(:,:,:) = self%frcs(:,:,:)
        else
            ! zero padding
            call self_out%new(self%ncls, newbox, newsmpd, self%nstates)
            self_out%frcs(:,:,:self%filtsz) = self%frcs(:,:,:)
            if( self_out%filtsz > self%filtsz ) self_out%frcs(:,:,self%filtsz+1:) = 0.0
        endif
    end subroutine pad

    ! I/O

    subroutine read( self, fname )
        class(class_frcs), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        integer          :: funit, io_stat
        call fopen(funit,fname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
        call fileiochk('class_frcs; read; open for read '//trim(fname), io_stat)
        read(unit=funit,pos=1) self%file_header
        ! re-create the object according to file_header info
        call self%new(nint(self%file_header(1)), nint(self%file_header(2)), self%file_header(3), nint(self%file_header(4)))
        read(unit=funit,pos=self%headsz + 1) self%frcs
        call fclose(funit)
    end subroutine read

    subroutine write( self, fname )
        class(class_frcs), intent(in) :: self
        character(len=*),  intent(in) :: fname
        integer          :: funit, io_stat
        call fopen(funit,fname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        call fileiochk('class_frcs; write; open for write '//trim(fname), io_stat)
        write(unit=funit,pos=1) self%file_header
        write(unit=funit,pos=self%headsz + 1) self%frcs
        call fclose(funit)
    end subroutine write

    subroutine print_frcs( self, fname, state )
        class(class_frcs), intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        integer, optional, intent(in)    :: state
        real, allocatable :: res(:)
        integer :: j, sstate, icls
        sstate = 1
        if( present(state) ) sstate = state
        call self%read(fname)
        res = get_resarr(self%box4frc_calc, self%smpd)
        do icls=1,self%ncls
            write(logfhandle,'(A,1X,I4)') '>>> FRC FOR clsECTION INDEX:', icls
            write(logfhandle,*) ''
            do j=1,size(res)
               write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', self%frcs(sstate,icls,j)
            end do
            write(logfhandle,*) ''
        end do
    end subroutine print_frcs

    ! destructor

    subroutine kill( self )
        class(class_frcs), intent(inout) :: self
        if( self%exists )then
            deallocate(self%res4frc_calc, self%frcs)
            self%ncls       = 0
            self%filtsz       = 0
            self%box4frc_calc = 0
            self%exists       = .false.
        endif
    end subroutine kill

end module simple_class_frcs
