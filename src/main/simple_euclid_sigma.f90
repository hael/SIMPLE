module simple_euclid_sigma
include 'simple_lib.f08'
use simple_parameters,       only: params_glob
use simple_polarft_corrcalc, only: polarft_corrcalc, pftcc_glob
use simple_oris,             only: oris
implicit none

public :: euclid_sigma, eucl_sigma_glob
private
#include "simple_local_flags.inc"

type euclid_sigma
    private
    integer               :: file_header(4)
    integer               :: nptcls      = 0
    logical               :: exists      = .false.
    integer               :: headsz      = 0
    integer               :: sigmassz    = 0
    character(len=STDLEN) :: fname
    logical               :: do_divide   = .false.
    real,    allocatable  :: divide_by(:)
    real,    allocatable  :: sigma2_noise(:,:)
    logical, allocatable  :: sigma2_exists_msk(:)
    integer, allocatable  :: pinds(:)
contains
    ! constructor
    procedure          :: new
    ! I/O
    procedure          :: read
    procedure          :: create_empty
    procedure          :: calc_and_write_sigmas
    procedure, private :: open_and_check_header
    ! getters / setters
    procedure          :: sigma2_exists
    procedure          :: set_do_divide
    procedure          :: get_do_divide
    procedure          :: set_divide_by
    procedure          :: get_divide_by
    ! destructor
    procedure          :: kill
end type euclid_sigma

class(euclid_sigma), pointer :: eucl_sigma_glob => null()

contains

    subroutine new( self, fname )
        class(euclid_sigma), target, intent(inout) :: self
        character(len=*),            intent(in)    :: fname
        real(sp)                                   :: r
        call self%kill
        allocate( self%sigma2_noise(params_glob%kfromto(1):params_glob%kfromto(2),1:pftcc_glob%get_nptcls()),&
            self%sigma2_exists_msk(1:pftcc_glob%get_nptcls()),&
            self%pinds(params_glob%fromp:params_glob%top),&
            self%divide_by(params_glob%kfromto(1):params_glob%kfromto(2)) )
        call pftcc_glob%assign_sigma2_noise(self%sigma2_noise, self%sigma2_exists_msk)
        call pftcc_glob%assign_pinds(self%pinds)
        self%fname = trim(fname)
        self%file_header(1)    = params_glob%fromp
        self%file_header(2)    = params_glob%top
        self%file_header(3)    = params_glob%kfromto(1)
        self%file_header(4)    = params_glob%kfromto(2)
        self%headsz            = sizeof(self%file_header)
        self%sigmassz          = sizeof(r)*(params_glob%kfromto(2)-params_glob%kfromto(1)+1)
        self%do_divide         = .false.
        self%exists            = .true.
        self%sigma2_noise      = 0.
        self%sigma2_exists_msk = .false.
        eucl_sigma_glob        => self
    end subroutine new

    ! I/O

    subroutine read( self, ptcl_mask )
        class(euclid_sigma),    intent(inout) :: self
        logical, optional,      intent(in)    :: ptcl_mask(params_glob%fromp:params_glob%top)
        integer                               :: funit
        integer                               :: iptcl
        real(sp)                              :: sigma2_noise_n(params_glob%kfromto(1):params_glob%kfromto(2))
        integer                               :: addr
        logical                               :: success
        call pftcc_glob%assign_pinds(self%pinds)
        if (.not. file_exists(trim(self%fname))) call self%create_empty()
        success = self%open_and_check_header(funit)
        if (success) then
            do iptcl = params_glob%fromp, params_glob%top
                if (.not. ptcl_mask(iptcl)) cycle
                addr = self%headsz + (iptcl - params_glob%fromp) * self%sigmassz + 1
                read(unit=funit,pos=addr) sigma2_noise_n
                self%sigma2_noise(:,self%pinds(iptcl))    = sigma2_noise_n
                self%sigma2_exists_msk(self%pinds(iptcl)) = (.not. any(sigma2_noise_n == 0.))
            end do
            call fclose(funit, errmsg='euclid_sigma; read; fhandle cose')
        else
            sigma2_noise_n = 0.
            do iptcl = params_glob%fromp, params_glob%top
                if (.not. ptcl_mask(iptcl)) cycle
                self%sigma2_noise(:,self%pinds(iptcl))    = sigma2_noise_n
                self%sigma2_exists_msk(self%pinds(iptcl)) = .false.
            end do
        end if
    end subroutine read

    subroutine calc_and_write_sigmas( self, os, o_peaks, ptcl_mask )
        class(euclid_sigma),              intent(inout) :: self
        class(oris),                      intent(inout) :: os
        type(oris),          allocatable, intent(inout) :: o_peaks(:)
        logical,             allocatable, intent(in)    :: ptcl_mask(:)
        integer :: iptcl, ipeak, iref, npeaks
        real :: sigmas_tmp(params_glob%kfromto(1):params_glob%kfromto(2))
        real :: sigma_contrib(params_glob%kfromto(1):params_glob%kfromto(2))
        real :: shvec(2)
        integer :: irot
        real    :: weight
        integer :: funit, io_stat
        integer :: addr
        logical :: success
        integer :: phys_i !< physical particle index
        if (.not. file_exists(trim(self%fname))) then
            call self%create_empty()
        else
            success = self%open_and_check_header(funit)
            call fclose(funit, errmsg='euclid_sigma; write ')
            if (.not. success) then
                THROW_HARD('sigmas file has wrong dimensions')
            end if
        end if
        self%sigma2_noise = 0.
!omp parallel do default(shared) schedule(static) private(iref,iptcl,ipeak,weight,shvec,irot,npeaks,sigmas_tmp,sigma_contrib) proc_bind(close)
        do iptcl = params_glob%fromp,params_glob%top
            if (.not. ptcl_mask(iptcl)) cycle
            if ( os%get_state(iptcl)==0 ) cycle
            sigmas_tmp = 0.
            npeaks = o_peaks(iptcl)%get_noris()
            do ipeak = 1,npeaks
                weight = o_peaks(iptcl)%get(ipeak, 'ow')
                if (weight < TINY) cycle
                iref   = int(o_peaks(iptcl)%get(ipeak, 'proj'))
                shvec  = o_peaks(iptcl)%get_2Dshift(ipeak)
                irot   = pftcc_glob%get_roind(360. - o_peaks(iptcl)%e3get(ipeak))
                call pftcc_glob%gencorr_sigma_contrib(iref, iptcl, shvec, irot, sigma_contrib)
                sigmas_tmp = sigmas_tmp + 0.5 * weight * sigma_contrib
            end do
            self%sigma2_noise(:,self%pinds(iptcl))    = sigmas_tmp
            self%sigma2_exists_msk(self%pinds(iptcl)) = (.not. any(sigmas_tmp == 0.))
        end do
!omp end parallel do
        call fopen(funit,trim(self%fname),access='STREAM',iostat=io_stat)
        call fileiochk('euclid_sigma; write; open for write '//trim(self%fname), io_stat)
        do iptcl = params_glob%fromp,params_glob%top
            if (.not. ptcl_mask(iptcl)) cycle
            if( os%get_state(iptcl)==0 ) cycle
            addr = self%headsz + (iptcl - params_glob%fromp) * self%sigmassz + 1
            write (funit,pos=addr) self%sigma2_noise(:,self%pinds(iptcl))
        end do
        call fclose(funit, errmsg='euclid_sigma; write; fhandle cose')
    end subroutine calc_and_write_sigmas

    function open_and_check_header( self, funit ) result ( success )
        class(euclid_sigma), intent(inout) :: self
        integer,             intent(out)   :: funit
        logical                            :: success
        integer                            :: io_stat
        integer                            :: fromp_here, top_here, kfromto_here(2)
        if (.not. file_exists(trim(self%fname))) then
            success = .false.
            return
        end if
        call fopen(funit,trim(self%fname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        call fileiochk('euclid_sigma; read_sigmas; open for read '//trim(self%fname), io_stat)
        read(unit=funit,pos=1) self%file_header
        fromp_here      = self%file_header(1)
        top_here        = self%file_header(2)
        kfromto_here(:) = self%file_header(3:4)
        if ((fromp_here.ne.params_glob%fromp) .or. (top_here.ne.params_glob%top) .or. &
            (kfromto_here(1).ne.params_glob%kfromto(1)) .or. (kfromto_here(2).ne.params_glob%kfromto(2))) then
            THROW_WARN('dimensions in sigmas file do not match')
            write (*,*) 'params_glob%fromp: ',   params_glob%fromp,   ' ; in sigmas file: ', fromp_here
            write (*,*) 'params_glob%top: ',     params_glob%top,     ' ; in sigmas file: ', top_here
            write (*,*) 'params_glob%kfromto: ', params_glob%kfromto, ' ; in sigmas file: ', kfromto_here
            write (*,*) 'resorting to cross-correlations.'
            call fclose(funit)
            success = .false.
        else
            success = .true.
        end if
    end function open_and_check_header

    subroutine create_empty( self )
        class(euclid_sigma), intent(in) :: self
        real(sp), allocatable           :: sigmas_empty(:,:)
        integer                         :: funit, io_stat
        allocate(sigmas_empty(params_glob%kfromto(1):params_glob%kfromto(2), params_glob%fromp:params_glob%top))
        sigmas_empty = 0.
        call fopen(funit,trim(self%fname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) self%file_header
        write(unit=funit,pos=self%headsz + 1) sigmas_empty
        call fclose(funit)
        deallocate(sigmas_empty)
    end subroutine create_empty

    ! getters / setters

    function sigma2_exists( self, iptcl ) result( l_flag )
        class(euclid_sigma), intent(in) :: self
        integer,             intent(in) :: iptcl
        logical                         :: l_flag
        l_flag = .false.
        if( .not. self%exists ) return
        if( self%pinds(iptcl) .eq. 0 )then
            write (*,*) 'iptcl = ', iptcl
            THROW_HARD('sigma2_exists. iptcl index wrong!')
        else
            l_flag = self%sigma2_exists_msk(self%pinds(iptcl))
        end if
    end function sigma2_exists

    subroutine set_do_divide( self, do_divide )
        class(euclid_sigma), intent(inout) :: self
        logical,             intent(in)    :: do_divide
        self%do_divide = do_divide
    end subroutine set_do_divide

    function get_do_divide( self ) result( res )
        class(euclid_sigma), intent(in) :: self
        logical                         :: res
        if( self%exists )then
            res = self%do_divide
        else
            res = .false.
        endif
    end function get_do_divide

    subroutine set_divide_by( self, iptcl )
        class(euclid_sigma), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        if (self%pinds(iptcl) .eq. 0) then
            write (*,*) ' iptcl = ', iptcl
            THROW_HARD('set_divide_by. iptcl index wrong!')
        else
            self%divide_by(:) = self%sigma2_noise(:,self%pinds(iptcl))
        end if
    end subroutine set_divide_by

    function get_divide_by( self ) result( res )
        class(euclid_sigma), intent(in) :: self
        real, allocatable               :: res(:)
        res = self%divide_by
    end function get_divide_by

    ! destructor

    subroutine kill( self )
        class(euclid_sigma), intent(inout) :: self
        if( self%exists )then
            if ( allocated(self%sigma2_noise) )      deallocate(self%sigma2_noise)
            if ( allocated(self%sigma2_exists_msk) ) deallocate(self%sigma2_exists_msk)
            if ( allocated(self%pinds) )             deallocate(self%pinds)
            if ( allocated(self%divide_by) )         deallocate(self%divide_by)
            self%file_header = 0
            self%headsz      = 0
            self%sigmassz    = 0
            self%nptcls      = 0
            self%do_divide   = .false.
            self%exists      = .false.
            eucl_sigma_glob  => null()
        endif
    end subroutine kill

end module simple_euclid_sigma
