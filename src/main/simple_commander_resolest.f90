! concrete commander: resolest for resolution estimation
module simple_commander_resolest
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,     only: parameters, params_glob
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_masker,         only: masker
implicit none

public :: fsc_commander
public :: opt_3D_filter_commander
public :: opt_2D_filter_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: fsc_commander
  contains
    procedure :: execute      => exec_fsc
end type fsc_commander

type, extends(commander_base) :: opt_2D_filter_commander
  contains
    procedure :: execute      => exec_opt_2D_filter
end type opt_2D_filter_commander

type, extends(commander_base) :: opt_3D_filter_commander
  contains
    procedure :: execute      => exec_opt_3D_filter
end type opt_3D_filter_commander

contains

    !> calculates Fourier shell correlation from Even/Odd Volume pairs
    subroutine exec_fsc( self, cline )
        use simple_estimate_ssnr, only: fsc2optlp_sub, fsc2TVfilt
        class(fsc_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd
        type(masker)      :: mskvol
        integer           :: j, find_plate, k_hp, k_lp, nyq
        real              :: res_fsc05, res_fsc0143
        real, allocatable :: res(:), corrs(:)
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! read even/odd pair
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        if( cline%defined('mskfile') )then
            if( file_exists(params%mskfile) )then
                call mskvol%new([params%box,params%box,params%box], params%smpd)
                call mskvol%read(params%mskfile)
                call even%zero_background
                call odd%zero_background
                call even%mul(mskvol)
                call odd%mul(mskvol)
            else
                THROW_HARD('mskfile: '//trim(params%mskfile)//' does not exist in cwd; exec_fsc')
            endif
        else
            ! spherical masking
            call even%mask(params%msk, 'soft')
            call odd%mask(params%msk, 'soft')
        endif
        ! forward FT
        call even%fft()
        call odd%fft()
        ! calculate FSC
        res = even%get_res()
        nyq = even%get_filtsz()
        allocate(corrs(nyq))
        call even%fsc(odd, corrs)
        if( params%l_phaseplate ) call phaseplate_correct_fsc(corrs, find_plate)
        do j=1,nyq
           write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
        end do
        call get_resolution(corrs, res, res_fsc05, res_fsc0143)
        k_hp = calc_fourier_index(params%hp, params%box, params%smpd)
        k_lp = calc_fourier_index(params%lp, params%box, params%smpd)
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
        write(logfhandle,'(A,1X,F8.4)') '>>> MEDIAN FSC (SPECSCORE):', median_nocopy(corrs(k_hp:k_lp))
        call even%kill
        call odd%kill
        call arr2file(corrs, 'FSC.bin')
        ! end gracefully
        call simple_end('**** SIMPLE_FSC NORMAL STOP ****')
    end subroutine exec_fsc

    subroutine exec_opt_3D_filter(self, cline)
        use simple_opt_filter, only: opt_filter_3D
        class(opt_3D_filter_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: even, odd, mskvol
        logical           :: have_mask_file
        character(len=90) :: file_tag
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call odd %new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd %read(params%vols(1))
        call even%read(params%vols(2))
        if( params%l_nonuniform )then
            file_tag = 'nonuniform_3D_filter_ext_'//int2str(params%smooth_ext)
        else
            file_tag = 'uniform_3D_filter_ext_'//int2str(params%smooth_ext)
        endif
        have_mask_file = .false.
        if( cline%defined('mskfile') )then
            if( file_exists(params%mskfile) )then
                call mskvol%new([params%box,params%box,params%box], params%smpd)
                call mskvol%read(params%mskfile)
                call even%zero_background
                call odd%zero_background
                call even%mul(mskvol)
                call odd%mul(mskvol)
                call mskvol%one_at_edge ! to expand before masking of reference internally
                have_mask_file = .true.
            else
                THROW_HARD('mskfile: '//trim(params%mskfile)//' does not exist in cwd; exec_opt_3D_filter')
            endif
        else
            ! spherical masking
            call even%mask(params%msk, 'soft')
            call odd%mask(params%msk, 'soft')
            call mskvol%disc([params%box,params%box,params%box], params%smpd,&
            &real(min(params%box/2, int(params%msk + COSMSKHALFWIDTH))))
        endif        
        call opt_filter_3D(odd, even, mskvol)
        if( have_mask_file )then
            call mskvol%read(params%mskfile) ! restore the soft edge
            call even%mul(mskvol)
            call odd%mul(mskvol)
        else
            call even%mask(params%msk, 'soft')
            call odd%mask(params%msk, 'soft')
        endif
        call odd%write(trim(file_tag)//'_odd.mrc')
        call even%write(trim(file_tag)//'_even.mrc')
        call odd%add(even)
        call odd%mul(0.5)
        call odd%write(trim(file_tag)//'_avg.mrc')
        ! destruct
        call odd%kill
        call even%kill
        call mskvol%kill
        ! end gracefully
        call simple_end('**** SIMPLE_OPT_3D_FILTER NORMAL STOP ****')
    end subroutine exec_opt_3D_filter
    
    subroutine exec_opt_2D_filter( self, cline )
        use simple_opt_filter, only: opt_2D_filter_sub
        use simple_tvfilter,   only: tvfilter
        use simple_class_frcs, only: class_frcs
        class(opt_2D_filter_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        character(len=:), allocatable :: file_tag
        type(image),      allocatable :: even(:), odd(:)
        type(parameters) :: params
        integer          :: iptcl
        ! init
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',    'yes')
        if( .not. cline%defined('smooth_ext') ) call cline%set('smooth_ext', 20.)
        call params%new(cline) 
        call find_ldim_nptcls(params%stk, params%ldim, params%nptcls)
        params%ldim(3) = 1 ! because we operate on stacks
        if( params%l_nonuniform )then
            file_tag = 'nonuniform_2D_filter_ext_'//int2str(params%smooth_ext)
        else
            file_tag = 'uniform_2D_filter_ext_'//int2str(params%smooth_ext)
        endif
        ! allocate
        allocate(odd(params%nptcls), even(params%nptcls))
        ! construct & read
        do iptcl = 1, params%nptcls
            call odd( iptcl)%new(params%ldim, params%smpd, .false.)
            call even(iptcl)%new(params%ldim, params%smpd, .false.)
            call odd( iptcl)%read(params%stk2, iptcl)
            call even(iptcl)%read(params%stk,  iptcl)
        enddo
        ! filter
        call opt_2D_filter_sub( even, odd )
        ! destruct
        do iptcl = 1, params%nptcls
            call odd( iptcl)%write(trim(file_tag)//'_odd.mrc',  iptcl)
            call even(iptcl)%write(trim(file_tag)//'_even.mrc', iptcl)
            call odd( iptcl)%kill()
            call even(iptcl)%kill()
        enddo
        call cline%set('odd_stk',  trim(file_tag)//'_odd.mrc')
        call cline%set('even_stk', trim(file_tag)//'_even.mrc')
        ! end gracefully
        call simple_end('**** SIMPLE_OPT_2D_FILTER NORMAL STOP ****')
    end subroutine exec_opt_2D_filter

end module simple_commander_resolest
