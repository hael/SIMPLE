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
    real,    allocatable  :: divide_by(:)
    real,    allocatable  :: sigma2_noise(:,:)
    logical, allocatable  :: sigma2_exists_msk(:)
    integer, allocatable  :: pinds(:)
    integer               :: file_header(4) = 0
    integer               :: kfromto(2)     = 0
    integer               :: headsz         = 0
    integer               :: sigmassz       = 0
    character(len=STDLEN) :: fname
    logical               :: do_divide = .false.
    logical               :: exists    = .false.

contains
    ! constructor
    procedure          :: new
    ! I/O
    procedure          :: read
    procedure, private :: create_empty
    procedure          :: calc_and_write_sigmas
    procedure          :: calc_and_write_sigmas2D
    procedure, private :: open_and_check_header
    procedure, private :: write
    ! getters / setters
    procedure          :: sigma2_exists
    procedure          :: set_do_divide
    procedure          :: get_do_divide
    procedure          :: set_sigma2
    procedure          :: get_sigma2
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
        self%kfromto = params_glob%kfromto
        allocate( self%sigma2_noise(self%kfromto(1):self%kfromto(2),1:pftcc_glob%get_nptcls()),&
            self%sigma2_exists_msk(1:pftcc_glob%get_nptcls()),&
            self%pinds(params_glob%fromp:params_glob%top),&
            self%divide_by(self%kfromto(1):self%kfromto(2)) )
        call pftcc_glob%assign_sigma2_noise(self%sigma2_noise, self%sigma2_exists_msk)
        call pftcc_glob%assign_pinds(self%pinds)
        self%fname = trim(fname)
        self%file_header(1)    = params_glob%fromp
        self%file_header(2)    = params_glob%top
        self%file_header(3:4)  = self%kfromto
        self%headsz            = sizeof(self%file_header)
        self%sigmassz          = sizeof(r)*(self%kfromto(2)-self%kfromto(1)+1)
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
        real(sp)                              :: sigma2_noise_n(self%kfromto(1):self%kfromto(2))
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

    !>  For soft assignment
    subroutine calc_and_write_sigmas( self, os, o_peaks, ptcl_mask )
        class(euclid_sigma),              intent(inout) :: self
        class(oris),                      intent(inout) :: os
        type(oris),          allocatable, intent(inout) :: o_peaks(:)
        logical,             allocatable, intent(in)    :: ptcl_mask(:)
        integer :: iptcl, ipeak, iref, npeaks, irot, funit
        real    :: sigmas_tmp(self%kfromto(1):self%kfromto(2))
        real    :: sigma_contrib(self%kfromto(1):self%kfromto(2))
        real    :: shvec(2), weight
        logical :: success
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
        !$omp parallel do default(shared) schedule(static) proc_bind(close)&
        !$omp private(iref,iptcl,ipeak,weight,shvec,irot,npeaks,sigmas_tmp,sigma_contrib)
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
        !$omp end parallel do
        call self%write(funit, os, ptcl_mask)
    end subroutine calc_and_write_sigmas

    !>  For 2D (hard assignment)
    subroutine calc_and_write_sigmas2D( self, os, ptcl_mask )
        class(euclid_sigma),              intent(inout) :: self
        class(oris),                      intent(inout) :: os
        logical,             allocatable, intent(in)    :: ptcl_mask(:)
        real    :: sigma_contrib(self%kfromto(1):self%kfromto(2))
        real    :: shvec(2), weight
        integer :: i,iptcl, iref, irot, funit
        logical :: success
        if (.not. file_exists(trim(self%fname))) then
            call self%create_empty()
        else
            success = self%open_and_check_header(funit)
            call fclose(funit, errmsg='euclid_sigma; write ')
            if (.not. success) THROW_HARD('sigmas file has wrong dimensions')
        endif
        self%sigma2_noise = 0.
        !$omp parallel do default(shared) schedule(static) private(iref,iptcl,weight,shvec,irot,i,sigma_contrib)&
        !$omp proc_bind(close)
        do iptcl = params_glob%fromp,params_glob%top
            if (.not. ptcl_mask(iptcl)  ) cycle
            if ( os%get_state(iptcl)==0 ) cycle
            weight = os%get(iptcl, 'w')
            if (weight < TINY) cycle
            i      = self%pinds(iptcl)
            iref   = nint(os%get(iptcl, 'class'))
            shvec  = os%get_2Dshift(iptcl)
            irot   = pftcc_glob%get_roind(360. - os%e3get(iptcl))
            call pftcc_glob%gencorr_sigma_contrib(iref, iptcl, shvec, irot, sigma_contrib)
            sigma_contrib = 0.5 * weight * sigma_contrib
            self%sigma2_noise(:,i)    = sigma_contrib
            self%sigma2_exists_msk(i) = (.not. any(sigma_contrib == 0.))
        end do
        !$omp end parallel do
        call self%write(funit, os, ptcl_mask)
    end subroutine calc_and_write_sigmas2D

    function open_and_check_header( self, funit ) result ( success )
        class(euclid_sigma), intent(inout) :: self
        integer,             intent(out)   :: funit
        logical                            :: success
        integer                            :: io_stat
        integer                            :: fromp_here, top_here
        if (.not. file_exists(trim(self%fname))) then
            success = .false.
            return
        end if
        call fopen(funit,trim(self%fname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        call fileiochk('euclid_sigma; read_sigmas; open for read '//trim(self%fname), io_stat)
        read(unit=funit,pos=1) self%file_header
        fromp_here      = self%file_header(1)
        top_here        = self%file_header(2)
        self%kfromto    = self%file_header(3:4)
        if ((fromp_here.ne.params_glob%fromp) .or. (top_here.ne.params_glob%top) .or. &
            (self%kfromto(1).ne.self%kfromto(1)) .or. (self%kfromto(2).ne.self%kfromto(2))) then
            THROW_WARN('dimensions in sigmas file do not match')
            write (*,*) 'params_glob%fromp: ',   params_glob%fromp,   ' ; in sigmas file: ', fromp_here
            write (*,*) 'params_glob%top: ',     params_glob%top,     ' ; in sigmas file: ', top_here
            write (*,*) 'self%kfromto: ', self%kfromto, ' ; in sigmas file: ', self%kfromto
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
        allocate(sigmas_empty(self%kfromto(1):self%kfromto(2), params_glob%fromp:params_glob%top))
        sigmas_empty = 0.
        call fopen(funit,trim(self%fname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) self%file_header
        write(unit=funit,pos=self%headsz + 1) sigmas_empty
        call fclose(funit)
        deallocate(sigmas_empty)
    end subroutine create_empty

    subroutine write( self, funit, os, ptcl_mask )
        class(euclid_sigma),  intent(inout) :: self
        integer,              intent(inout) :: funit
        class(oris),          intent(in)    :: os
        logical, allocatable, intent(in)    :: ptcl_mask(:)
        integer :: io_stat, addr, iptcl
        call fopen(funit,trim(self%fname),access='STREAM',iostat=io_stat)
        call fileiochk('euclid_sigma; write; open for write '//trim(self%fname), io_stat)
        do iptcl = params_glob%fromp,params_glob%top
            if (.not. ptcl_mask(iptcl) ) cycle
            if( os%get_state(iptcl)==0 ) cycle
            addr = self%headsz + (iptcl - params_glob%fromp) * self%sigmassz + 1
            write(funit,pos=addr) self%sigma2_noise(:,self%pinds(iptcl))
        end do
        call fclose(funit, errmsg='euclid_sigma; write; fhandle close')
    end subroutine

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

    logical function get_do_divide( self )
        class(euclid_sigma), intent(in) :: self
        get_do_divide = .false.
        if( self%exists ) get_do_divide = self%do_divide
    end function get_do_divide

    subroutine set_sigma2( self, iptcl )
        class(euclid_sigma), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        if (self%pinds(iptcl) .eq. 0) then
            write (*,*) ' iptcl = ', iptcl
            THROW_HARD('set_divide_by. iptcl index wrong!')
        else
            self%divide_by(:) = self%sigma2_noise(:,self%pinds(iptcl))
        end if
    end subroutine set_sigma2

    subroutine get_sigma2( self, nyq, sigma2 )
        class(euclid_sigma), intent(in)  :: self
        integer,             intent(in)  :: nyq
        real,                intent(out) :: sigma2(:)
        real    :: scale, loc, ld
        integer :: k, lk, kstart
        sigma2 = 1.
        if( .not.self%exists )return
        if( nyq == self%kfromto(2) )then
            sigma2(self%kfromto(1):self%kfromto(2)) = self%divide_by(:)
        else if( nyq > self%kfromto(2) )then
            ! requires interpolation, assumed that kfromto(2) is nyquist index
            scale       = real(self%kfromto(2)) / real(nyq)
            kstart      = ceiling(real(self%kfromto(1))/scale)
            sigma2(nyq) = self%divide_by(self%kfromto(2))
            do k = kstart,nyq-1
                loc = real(k)*scale
                lk  = floor(loc)
                ld  = loc-real(lk)
                sigma2(k) = ld*self%divide_by(lk+1) + (1.-ld)*self%divide_by(lk)
            enddo
        else
            THROW_HARD('Incompatible requested size; get_sigma2')
        endif
    end subroutine get_sigma2

    ! destructor

    subroutine kill( self )
        class(euclid_sigma), intent(inout) :: self
        if( self%exists )then
            if ( allocated(self%sigma2_noise) )      deallocate(self%sigma2_noise)
            if ( allocated(self%sigma2_exists_msk) ) deallocate(self%sigma2_exists_msk)
            if ( allocated(self%pinds) )             deallocate(self%pinds)
            if ( allocated(self%divide_by) )         deallocate(self%divide_by)
            self%file_header  = 0
            self%headsz       = 0
            self%sigmassz     = 0
            self%kfromto      = 0
            self%do_divide    = .false.
            self%exists       = .false.
            eucl_sigma_glob  => null()
        endif
    end subroutine kill

end module simple_euclid_sigma
