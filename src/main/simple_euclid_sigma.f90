module simple_euclid_sigma
include 'simple_lib.f08'
use simple_parameters,               only: params_glob
use simple_polarft_corrcalc,         only: polarft_corrcalc
!use simple_strategy3D_alloc,         only: s3D              !cannot use here bc of cyclical dependencies
use simple_oris, only: oris
implicit none

public :: euclid_sigma
private

type euclid_sigma
    private
    integer                         :: file_header(4)
    integer                         :: nptcls      = 0
    logical                         :: exists      = .false.
    integer                         :: headsz      = 0
    integer                         :: sigmassz    = 0
    character(len=STDLEN)           :: fname
    real,    allocatable            :: sigma2_noise(:,:)
    logical, allocatable            :: sigma2_exists_msk(:)
    integer, pointer                :: pinds(:) => null()
contains
    ! constructor
    procedure                       :: new
    ! I/O
    procedure                       :: read
    procedure                       :: create_empty
    procedure                       :: calc_and_write_sigmas
    procedure, private              :: open_and_check_header
    ! destructor
    procedure                       :: kill
end type euclid_sigma

contains

    ! constructor

    subroutine new( self, fname, pftcc )
        class(euclid_sigma),    intent(inout) :: self
        character(len=*),       intent(in)    :: fname
        type(polarft_corrcalc), intent(inout) :: pftcc
        real(sp)                              :: r
        call self%kill
        allocate( self%sigma2_noise(params_glob%kfromto(1):params_glob%kfromto(2),1:pftcc%get_nptcls()),&
            self%sigma2_exists_msk(1:pftcc%get_nptcls()) )
        call pftcc%assign_sigma2_noise(self%sigma2_noise, self%sigma2_exists_msk)
        call pftcc%assign_pinds(self%pinds)
        self%fname = trim(fname)
        self%file_header(1) = params_glob%fromp
        self%file_header(2) = params_glob%top
        self%file_header(3) = params_glob%kfromto(1)
        self%file_header(4) = params_glob%kfromto(2)
        self%headsz         = sizeof(self%file_header)
        self%sigmassz       = sizeof(r)*(params_glob%kfromto(2)-params_glob%kfromto(1)+1)
        self%exists         = .true.
    end subroutine new

    ! I/O

    subroutine read( self, pftcc, ptcl_mask )
        class(euclid_sigma),    intent(inout) :: self
        type(polarft_corrcalc), intent(inout) :: pftcc
        logical, optional,      intent(in)    :: ptcl_mask(params_glob%fromp:params_glob%top)
        integer                               :: funit
        integer                               :: iptcl
        real(sp)                              :: sigma2_noise_n(params_glob%kfromto(1):params_glob%kfromto(2))
        integer                               :: addr
        logical                               :: success
        call pftcc%assign_pinds(self%pinds)
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
                self%sigma2_noise(:,self%pinds(iptcl)) = sigma2_noise_n
                self%sigma2_exists_msk(self%pinds(iptcl)) = .false.
            end do
        end if
    end subroutine read

    subroutine calc_and_write_sigmas( self, pftcc, os, o_peaks, ptcl_mask )
        class(euclid_sigma),              intent(inout) :: self
        class(polarft_corrcalc), pointer, intent(in)    :: pftcc
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
                write (*,*) 'ERROR: sigmas file has wrong dimensions. stop'
                stop 1
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
                irot   = pftcc%get_roind(360. - o_peaks(iptcl)%e3get(ipeak))
                call pftcc%gencorr_sigma_contrib(iref, iptcl, shvec, irot, sigma_contrib)
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
            write (*,*) 'WARNING: dimensions in sigmas file do not match.'
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

    ! destructor

    subroutine kill( self )
        class(euclid_sigma), intent(inout) :: self
        if( self%exists )then
            if ( allocated(self%sigma2_noise) )      deallocate(self%sigma2_noise)
            if ( allocated(self%sigma2_exists_msk) ) deallocate(self%sigma2_exists_msk)
            self%pinds       => null()
            self%file_header = 0
            self%headsz      = 0
            self%sigmassz    = 0
            self%nptcls      = 0
            self%exists      = .false.
        endif
    end subroutine kill

end module simple_euclid_sigma
