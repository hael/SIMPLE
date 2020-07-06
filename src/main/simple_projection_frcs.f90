module simple_projection_frcs
include 'simple_lib.f08'
implicit none

public :: projection_frcs
private
#include "simple_local_flags.inc"

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
    procedure          :: getter
    procedure          :: estimate_res
    procedure          :: estimate_find_for_eoavg
    procedure          :: estimate_lp_for_align
    procedure          :: downsample
    procedure          :: upsample
    ! I/O
    procedure          :: read
    procedure          :: write
    procedure          :: print_frcs
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
        if(alloc_stat .ne. 0)call allocchk('new; simple_projection_frcs', alloc_stat)
        self%frcs   = 0.0
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
            write(logfhandle,*) self%nprojs
            write(logfhandle,*) 'exists: ', self%exists
            write(logfhandle,*) 'proj: ', proj
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
        class(projection_frcs), intent(in)    :: self
        real,                   intent(in)    :: frc(:)
        real,                   intent(inout) :: res_frc05, res_frc0143
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
            THROW_HARD('size of input frc not conforming; set_frc')
        else
            self%frcs(sstate,proj,:) = frc
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
        allocate(filt_resamp(filtsz_resamp),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("simple_estimate_ssnr::resample_filter ",alloc_stat)
        do k=1,filtsz_resamp
            call find(res_orig, filtsz_orig, res_new(k), ind, dist)
            filt_resamp(k) = filt_orig(ind)
        end do
    end function resample_filter

    !>  getter for values interepreted as FRCs
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

    !>  getter for values interepreted as FRCs
    subroutine frc_getter( self, proj, hpind_fsc, phaseplate, frc, state )
        class(projection_frcs), intent(in)  :: self
        integer,                intent(in)  :: proj, hpind_fsc
        logical,                intent(in)  :: phaseplate
        real,                   intent(out) :: frc(self%filtsz)
        integer, optional,      intent(in)  :: state
        integer :: sstate, find_plate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( proj, sstate, 'ERROR, out of bounds in frc_getter')
        frc = self%frcs(sstate,proj,:)
        if( phaseplate )then
            if( any(frc > 0.5) )call phaseplate_correct_fsc(frc, find_plate)
        endif
        if( hpind_fsc > 0 ) frc(:hpind_fsc) = frc(hpind_fsc + 1)
    end subroutine frc_getter

    !>  getter for raw values
    subroutine getter( self, proj, frc, state )
        class(projection_frcs), intent(in)  :: self
        integer,                intent(in)  :: proj
        real,                   intent(out) :: frc(self%filtsz)
        integer, optional,      intent(in)  :: state
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        call self%raise_exception( proj, sstate, 'ERROR, out of bounds in frc_getter')
        frc = self%frcs(sstate,proj,:)
    end subroutine getter

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
        real    :: lplims(self%nprojs),lp3(3)
        integer :: sstate, iproj
        real    :: lp, res_frc05, res_frc0143
        sstate = 1
        if( present(state) ) sstate = state
        ! order FRCs according to low-pass limit (best first)
        do iproj=1,self%nprojs
            call get_resolution(self%frcs(sstate,iproj,:), self%res4frc_calc, res_frc05, res_frc0143)
            call self%bound_res(self%frcs(sstate,iproj,:), res_frc05, res_frc0143)
            lplims(iproj) = res_frc0143
        end do
        lp3 = min3(lplims)
        ! return median of top three clusters
        lp = median(lp3)
    end function estimate_lp_for_align

    subroutine downsample( self, newbox, self_out )
        use simple_estimate_ssnr, only: subsample_optlp
        class(projection_frcs), intent(in)  :: self
        integer,                intent(in)  :: newbox
        type(projection_frcs),  intent(out) :: self_out
        integer :: istate, iproj
        real    :: new_smpd
        if( newbox > self%box4frc_calc ) THROW_HARD('New > old filter size; downsample')
        new_smpd = self%smpd * real(self%box4frc_calc) / real(newbox)
        call self_out%new(self%nprojs, newbox, new_smpd, self%nstates)
        do istate = 1,self%nstates
            do iproj = 1,self%nprojs
                call subsample_optlp(self%filtsz, self_out%filtsz, self%frcs(istate, iproj,:), self_out%frcs(istate, iproj,:))
            enddo
        enddo
    end subroutine downsample

    subroutine upsample( self, newsmpd, newbox, self_out )
        class(projection_frcs), intent(in)  :: self
        real,                   intent(in)  :: newsmpd
        integer,                intent(in)  :: newbox
        type(projection_frcs),  intent(out) :: self_out
        real    :: x, d, scale
        integer :: istate, iproj, sh, l,r
        if( newbox <= self%box4frc_calc ) THROW_HARD('New <= old filter size; downsample')
        call self_out%new(self%nprojs, newbox, newsmpd, self%nstates)
        scale = real(self%filtsz) / real(self_out%filtsz)
        do istate = 1,self%nstates
            do iproj = 1,self%nprojs
                do sh = 1,self_out%filtsz
                    x = scale*real(sh)
                    l = max(1, floor(x))
                    r = min(self%filtsz, ceiling(x))
                    if( l==r )then
                        self_out%frcs(istate, iproj,sh) = self%frcs(istate, iproj,l)
                    else
                        d = x-real(l)
                        self_out%frcs(istate, iproj,sh) = (1.-d) * self%frcs(istate, iproj,l) +&
                                                        &     d  * self%frcs(istate, iproj,r)
                    endif
                enddo
            enddo
        enddo
    end subroutine upsample

    ! I/O

    subroutine read( self, fname )
        class(projection_frcs), intent(inout) :: self
        character(len=*),       intent(in)    :: fname
        integer          :: funit, io_stat
        call fopen(funit,fname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
        call fileiochk('projection_frcs; read; open for read '//trim(fname), io_stat)
        read(unit=funit,pos=1) self%file_header
        ! re-create the object according to file_header info
        call self%new(nint(self%file_header(1)), nint(self%file_header(2)), self%file_header(3), nint(self%file_header(4)))
        read(unit=funit,pos=self%headsz + 1) self%frcs
        call fclose(funit, errmsg='projection_frcs; read; fhandle close')
    end subroutine read

    subroutine write( self, fname )
        class(projection_frcs), intent(in) :: self
        character(len=*),       intent(in) :: fname
        integer          :: funit, io_stat
        call fopen(funit,fname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        call fileiochk('projection_frcs; write; open for write '//trim(fname), io_stat)
        write(unit=funit,pos=1) self%file_header
        write(unit=funit,pos=self%headsz + 1) self%frcs
        call fclose(funit, errmsg='projection_frcs; write; fhandle close')
    end subroutine write

    subroutine print_frcs( self, fname, state )
        class(projection_frcs), intent(inout) :: self
        character(len=*),       intent(in)    :: fname
        integer,      optional, intent(in)    :: state
        real, allocatable :: res(:)
        integer :: j, sstate, iproj
        sstate = 1
        if( present(state) ) sstate = state
        call self%read(fname)
        res = get_resarr(self%box4frc_calc, self%smpd)
        do iproj=1,self%nprojs
            write(logfhandle,'(A,1X,I4)') '>>> FRC FOR PROJECTION INDEX:', iproj
            write(logfhandle,*) ''
            do j=1,size(res)
               write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', self%frcs(sstate,iproj,j)
            end do
            write(logfhandle,*) ''
        end do
    end subroutine print_frcs

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
