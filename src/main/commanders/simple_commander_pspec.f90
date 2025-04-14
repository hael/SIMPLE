module simple_commander_pspec
include 'simple_lib.f08'
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_pspecs,         only: pspecs
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters, params_glob
use simple_sp_project,     only: sp_project
use simple_builder,        only: builder, build_glob
use simple_stack_io,       only: stack_io
use simple_strategy2D3D_common
implicit none

public :: analyze_pspecs_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: analyze_pspecs_commander
    contains
    procedure :: execute => exec_analyze_pspecs
end type analyze_pspecs_commander

contains

    subroutine exec_analyze_pspecs( self, cline )
        class(analyze_pspecs_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)   :: params
        type(builder)      :: build
        type(pspecs)       :: pows
        integer            :: nptcls, iptcl_batch, iptcl, ind_in_stk, batch_end, ithr
        integer            :: batchsz, nbatches, batchsz_max, ibatch, batch_start, state
        logical, parameter :: DEBUG = .true.
        character(len=:), allocatable :: stkname
        integer,          allocatable :: batches(:,:)
        if( .not. cline%defined('hp')      ) call cline%set('hp',       20.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',        6.)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'yes')
        ! build and params
        call build%init_params_and_build_general_tbox(cline, params)
        call build%spproj%update_projinfo(cline)
        call build%spproj%write_segment_inside('projinfo')
        ! create object for power spectra
        nptcls = build%spproj%get_nptcls()
        ! call pows%new(nptcls, build%img, params%hp, params%lp )
        ! prep batch alignment
        batchsz     = params%nthr * BATCHTHRSZ
        nbatches    = ceiling(real(nptcls)/real(batchsz))
        batches     = split_nobjs_even(nptcls, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
        call prepimgbatch(batchsz_max)

        if( DEBUG )then
            print *, 'nthr        ', params%nthr
            print *, 'nptcls      ', nptcls
            print *, 'batchsz     ', batchsz
            print *, 'nbatches    ', nbatches
            print *, 'batchsz_max ', batchsz_max
        endif

        ! do ibatch=1,nbatches
        !     batch_start = batches(ibatch,1)
        !     batch_end   = batches(ibatch,2)
        !     batchsz     = batch_end - batch_start + 1
        !     call read_imgbatch([batch_start,batch_end])
        !     !$omp parallel do default(shared) private(iptcl_batch,iptcl,ithr,state)&
        !     !$omp schedule(static) proc_bind(close)
        !     do iptcl_batch = 1, batchsz               ! particle batch index
        !         ithr  = omp_get_thread_num() + 1      ! thread index
        !         iptcl = batch_start + iptcl_batch - 1 ! particle index
        !         call pows%set_pspec(iptcl, build%imgbatch(iptcl_batch), params%msk)
        !         state = build%spproj%os_ptcl2D%get_state(iptcl)
        !         if( state > 0 )then
        !             call pows%set_class_good(iptcl)
        !         else
        !             call pows%set_class_bad(iptcl)
        !         endif
        !     end do
        !     !$omp end parallel do
        ! end do
        ! deallocate(batches)
        ! call killimgbatch

        ! print *, '# empty ', pows%count_empty()

        ! call pows%master

        call simple_end('**** SIMPLE_ANALYZE_PSPECS NORMAL STOP ****')
    end subroutine exec_analyze_pspecs

end module simple_commander_pspec