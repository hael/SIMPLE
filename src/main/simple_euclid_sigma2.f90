module simple_euclid_sigma2
include 'simple_lib.f08'
use simple_parameters,       only: params_glob
use simple_polarft_corrcalc, only: polarft_corrcalc, pftcc_glob
use simple_oris,             only: oris
implicit none

public :: euclid_sigma2, eucl_sigma2_glob
private
#include "simple_local_flags.inc"

type euclid_sigma2
    private
    real,    allocatable  :: divide_by(:)
    real,    allocatable  :: sigma2_noise(:,:)      !< for reading & calculating individual contributions
    real,    allocatable  :: mic_sigma2_noise(:,:)  !< weighted average for euclidian distance calculation and reconstruction
    logical, allocatable  :: sigma2_exists_msk(:)   !< particle &
    integer, allocatable  :: pinds(:)
    integer, allocatable  :: micinds(:)
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
    procedure, private :: open_and_check_header
    procedure, private :: write
    ! getters / setters
    procedure          :: sigma2_exists
    procedure, private :: set_do_divide
    procedure          :: get_do_divide
    procedure          :: set_sigma2
    procedure          :: get_sigma2
    ! destructor
    procedure          :: kill_ptclsigma2
    procedure          :: kill
end type euclid_sigma2

class(euclid_sigma2), pointer :: eucl_sigma2_glob => null()

contains

    subroutine new( self, fname )
        class(euclid_sigma2), target, intent(inout) :: self
        character(len=*),            intent(in)    :: fname
        real(sp)                                   :: r
        call self%kill
        self%kfromto = params_glob%kfromto
        allocate( self%sigma2_noise(self%kfromto(1):self%kfromto(2),1:pftcc_glob%get_nptcls()),&
            self%sigma2_exists_msk(1:pftcc_glob%get_nptcls()),&
            self%pinds(params_glob%fromp:params_glob%top),&
            self%divide_by(self%kfromto(1):self%kfromto(2)) )
        call pftcc_glob%assign_sigma2_noise(self%sigma2_noise )
        call pftcc_glob%assign_pinds(self%pinds)
        self%fname            = trim(fname)
        self%file_header(1)   = params_glob%fromp
        self%file_header(2)   = params_glob%top
        self%file_header(3:4) = self%kfromto
        self%headsz           = sizeof(self%file_header)
        self%sigmassz         = sizeof(r)*(self%kfromto(2)-self%kfromto(1)+1)
        self%sigma2_noise     = 0.
        self%do_divide          = .false.
        self%sigma2_exists_msk  = .false.
        self%exists             = .true.
        eucl_sigma2_glob       => self
    end subroutine new

    ! I/O

    subroutine read( self, os, ptcl_mask )
        class(euclid_sigma2),    intent(inout) :: self
        class(oris),            intent(inout) :: os
        logical,                intent(in)    :: ptcl_mask(params_glob%fromp:params_glob%top)
        real, allocatable :: wsums(:)
        real(sp) :: sigma2_noise_n(self%kfromto(1):self%kfromto(2)),w
        integer  :: fromm, tom, funit, iptcl, addr, pind, imic
        logical  :: success, defined
        call pftcc_glob%assign_pinds(self%pinds)
        self%sigma2_noise      = 0.
        self%sigma2_exists_msk = .false.
        if (.not. file_exists(trim(self%fname))) call self%create_empty()
        if( params_glob%cc_objfun /= OBJFUN_EUCLID ) return
        success = self%open_and_check_header(funit)
        if (success) then
            ! mic related init
            allocate(self%micinds(params_glob%fromp:params_glob%top),source=0)
            fromm = huge(fromm)
            tom   = -1
            do iptcl = params_glob%fromp, params_glob%top
                imic  = nint(os%get(iptcl, 'stkind'))
                fromm = min(fromm, imic)
                tom   = max(tom, imic)
            enddo
            allocate(self%mic_sigma2_noise(self%kfromto(1):self%kfromto(2),fromm:tom),&
            &wsums(fromm:tom),source=0.)
            ! read & sum
            do iptcl = params_glob%fromp, params_glob%top
                if( os%get_state(iptcl) == 0 ) cycle
                addr = self%headsz + (iptcl - params_glob%fromp) * self%sigmassz + 1
                read(unit=funit,pos=addr) sigma2_noise_n
                defined = any(sigma2_noise_n > TINY)
                if( .not.defined .and. ptcl_mask(iptcl) )then
                    print *,'UNDEFINED SIGMA2:', iptcl
                endif
                ! sets mask
                if( ptcl_mask(iptcl) ) self%sigma2_exists_msk(self%pinds(iptcl)) = .true.
                ! accumulates weighted sum
                w    = os%get(iptcl, 'w')
                imic = nint(os%get(iptcl, 'stkind'))
                self%micinds(iptcl) = imic
                wsums(imic) = wsums(imic) + w
                self%mic_sigma2_noise(:,imic) = self%mic_sigma2_noise(:,imic) + w*sigma2_noise_n
            end do
            call fclose(funit, errmsg='euclid_sigma2; read; fhandle close')
            ! micrograph weighted average
            do imic = fromm,tom
                w = wsums(imic)
                if( w <= TINY ) cycle
                self%mic_sigma2_noise(:,imic) = self%mic_sigma2_noise(:,imic) / w
            enddo
            ! transfer micrograph sigma2 to masked ptcls only for euclidian distance calculation
            do iptcl = params_glob%fromp, params_glob%top
                if( .not.ptcl_mask(iptcl) )cycle
                imic = self%micinds(iptcl)
                pind = self%pinds(iptcl)
                self%sigma2_noise(:,pind) = self%mic_sigma2_noise(:,imic)
            end do
            deallocate(wsums)
        end if
    end subroutine read

    !>  For soft assignment
    subroutine calc_and_write_sigmas( self, os, o_peaks, ptcl_mask )
        class(euclid_sigma2),              intent(inout) :: self
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
            call fclose(funit, errmsg='euclid_sigma2; write ')
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
                iref   = nint(o_peaks(iptcl)%get(ipeak, 'proj'))
                shvec  = o_peaks(iptcl)%get_2Dshift(ipeak)
                irot   = pftcc_glob%get_roind(360. - o_peaks(iptcl)%e3get(ipeak))
                call pftcc_glob%gencorr_sigma_contrib(iref, iptcl, shvec, irot, sigma_contrib)
                sigmas_tmp = sigmas_tmp + 0.5 * weight * sigma_contrib
            end do
            self%sigma2_noise(:,self%pinds(iptcl)) = sigmas_tmp
        end do
        !$omp end parallel do
        call self%write(funit, os, ptcl_mask)
    end subroutine calc_and_write_sigmas

    function open_and_check_header( self, funit ) result ( success )
        class(euclid_sigma2), intent(inout) :: self
        integer,             intent(out)   :: funit
        logical                            :: success
        integer                            :: io_stat
        integer                            :: fromp_here, top_here
        if (.not. file_exists(trim(self%fname))) then
            success = .false.
            return
        end if
        call fopen(funit,trim(self%fname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        call fileiochk('euclid_sigma2; read_sigmas; open for read '//trim(self%fname), io_stat)
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
        class(euclid_sigma2), intent(in) :: self
        integer  :: funit, io_stat
        real(sp) :: sigmas_empty(self%kfromto(1):self%kfromto(2), params_glob%fromp:params_glob%top)
        sigmas_empty = 0.
        call fopen(funit,trim(self%fname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) self%file_header
        write(unit=funit,pos=self%headsz + 1) sigmas_empty
        call fclose(funit)
    end subroutine create_empty

    subroutine write( self, funit, os, ptcl_mask )
        class(euclid_sigma2),  intent(inout) :: self
        integer,              intent(inout) :: funit
        class(oris),          intent(in)    :: os
        logical, allocatable, intent(in)    :: ptcl_mask(:)
        integer :: io_stat, addr, iptcl
        call fopen(funit,trim(self%fname),access='STREAM',iostat=io_stat)
        call fileiochk('euclid_sigma2; write; open for write '//trim(self%fname), io_stat)
        do iptcl = params_glob%fromp,params_glob%top
            if (.not. ptcl_mask(iptcl) ) cycle
            if( os%get_state(iptcl)==0 ) cycle
            addr = self%headsz + (iptcl - params_glob%fromp) * self%sigmassz + 1
            write(funit,pos=addr) self%sigma2_noise(:,self%pinds(iptcl))
        end do
        call fclose(funit, errmsg='euclid_sigma2; write; fhandle close')
    end subroutine

    ! getters / setters

    logical function sigma2_exists( self, iptcl )
        class(euclid_sigma2), intent(in) :: self
        integer,             intent(in) :: iptcl
        sigma2_exists = .false.
        if( .not. self%exists ) return
        if( self%pinds(iptcl) .eq. 0 )then
            write (*,*) 'iptcl = ', iptcl
            THROW_HARD('sigma2_exists. iptcl index wrong!')
        else
            sigma2_exists = self%sigma2_exists_msk(self%pinds(iptcl))
        end if
    end function sigma2_exists

    subroutine set_do_divide( self, do_divide )
        class(euclid_sigma2), intent(inout) :: self
        logical,             intent(in)    :: do_divide
        self%do_divide = do_divide
    end subroutine set_do_divide

    logical function get_do_divide( self )
        class(euclid_sigma2), intent(in) :: self
        get_do_divide = .false.
        if( self%exists ) get_do_divide = self%do_divide
    end function get_do_divide

    !>  push sigma2 average to buffer
    subroutine set_sigma2( self, iptcl )
        class(euclid_sigma2), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        integer :: micind
        self%do_divide = .false.
        if( .not. self%exists ) return
        if( self%sigma2_exists(iptcl) )then
            micind = self%micinds(iptcl)
            if( micind <= 0 ) THROW_HARD('set_sigma2; mic index wrong! micind='//int2str(micind))
            self%divide_by(:) = self%mic_sigma2_noise(:,micind)
            self%do_divide    = .true.
        else
            self%divide_by(:) = 1.
        end if
    end subroutine set_sigma2

    !>  fetch sigma2 average from buffer
    subroutine get_sigma2( self, nyq, sigma2 )
        class(euclid_sigma2), intent(in)  :: self
        integer,              intent(in)  :: nyq ! oversampling nyquist limit
        real,                 intent(out) :: sigma2(1:2*nyq)
        real    :: scale, loc, ld
        integer :: k, kstart, lk
        sigma2 = 1.
        if( .not.self%exists )return
        if( nyq == self%kfromto(2) )then
            sigma2(self%kfromto(1):self%kfromto(2)) = self%divide_by
        else if( nyq > self%kfromto(2) )then
            ! resampling
            scale  = real(self%kfromto(2)) / real(nyq)
            kstart = ceiling(real(self%kfromto(1))/scale)
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

    subroutine kill_ptclsigma2( self )
        class(euclid_sigma2), intent(inout) :: self
        if(allocated(self%sigma2_noise)) deallocate(self%sigma2_noise)
    end subroutine kill_ptclsigma2

    subroutine kill( self )
        class(euclid_sigma2), intent(inout) :: self
        if( self%exists )then
            call self%kill_ptclsigma2
            if(allocated(self%sigma2_exists_msk)) deallocate(self%sigma2_exists_msk)
            if(allocated(self%pinds))             deallocate(self%pinds)
            if(allocated(self%divide_by))         deallocate(self%divide_by)
            if(allocated(self%mic_sigma2_noise))  deallocate(self%mic_sigma2_noise)
            if(allocated(self%micinds))           deallocate(self%micinds)
            self%file_header  = 0
            self%headsz       = 0
            self%sigmassz     = 0
            self%kfromto      = 0
            self%do_divide    = .false.
            self%exists       = .false.
            eucl_sigma2_glob  => null()
        endif
    end subroutine kill

end module simple_euclid_sigma2
