module simple_reconstructor_tester
#include "simple_lib.f08"
use simple_build,          only: build
use simple_params,         only: params
use simple_cmdline,        only: cmdline
use simple_kbinterpol,     only: kbinterpol
use simple_prep4cgrid,     only: prep4cgrid
use simple_ori,            only: ori
use simple_oris,           only: oris
use simple_image,          only: image
use simple_hadamard_common ! use all in there
use simple_timer           ! use all in there
use simple_defs            ! use all in there
implicit none

public :: exec_old_school_rec, exec_rec_batch_gridprep
private

integer(timer_int_kind) :: t_init, t_rec, t_norm, t_tot, t_read, t_gridprep, t_gridding
real(timer_int_kind)    :: rt_init, rt_rec, rt_norm, rt_tot, rt_read, rt_gridprep, rt_gridding
character(len=STDLEN)   :: benchfname
integer                 :: fnr

contains

    subroutine exec_old_school_rec( b, p, cline )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        type(ori)        :: orientation
        type(kbinterpol) :: kbwin
        type(prep4cgrid) :: gridprep
        integer          :: iptcl
        real             :: reslim
        t_init = tic()
        t_tot  = t_init
        ! make the gridding prepper
        kbwin = b%eorecvols(1)%get_kbwin()
        call gridprep%new(b%img, kbwin, [p%boxpd,p%boxpd,1])
        ! init volumes
        call preprecvols(b, p)
        rt_init = toc(t_init)
        ! reconstruction
        t_rec = tic()
        do iptcl=1,p%nptcls
            orientation = b%a%get_ori(iptcl)
            call read_img_and_norm(b, p, iptcl)
            call gridprep%prep(b%img, b%img_pad)
            call grid_ptcl(b, p, b%se, orientation)
        end do
        rt_rec = toc(t_rec)
        ! normalise structure factors
        t_norm = tic()
        call eonorm_struct_facts(b, p, cline, reslim)
        rt_norm = toc(t_norm)
        ! recvols % gridprep not needed anymore
        call killrecvols(b, p)
        call gridprep%kill
        rt_tot = toc(t_tot)
        call gridprep%kill
        benchfname = 'REC_TEST_OLD_SCHOOL_BENCH.txt'
        call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
        write(fnr,'(a)') '*** TIMINGS (s) ***'
        write(fnr,'(a,1x,f9.2)') 'initialisation  : ', rt_init
        write(fnr,'(a,1x,f9.2)') 'gridding        : ', rt_rec
        write(fnr,'(a,1x,f9.2)') 'eonorm          : ', rt_norm
        write(fnr,'(a,1x,f9.2)') 'total time      : ', rt_tot
        write(fnr,'(a)') ''
        write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
        write(fnr,'(a,1x,f9.2)') 'initialisation  : ', (rt_init/rt_tot) * 100.
        write(fnr,'(a,1x,f9.2)') 'gridding        : ', (rt_rec/rt_tot)  * 100.
        write(fnr,'(a,1x,f9.2)') 'eonorm          : ', (rt_norm/rt_tot) * 100.
        write(fnr,'(a,1x,f9.2)') '% accounted for : ', ((rt_init+rt_rec+rt_norm)/rt_tot) * 100.
        call fclose(fnr)
    end subroutine exec_old_school_rec

    subroutine exec_rec_batch_gridprep( b, p, cline )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        type(ori)                :: orientation
        type(kbinterpol)         :: kbwin
        type(prep4cgrid)         :: gridprep
        integer                  :: iptcl, i, iptcl_batch, batchlims(2), ibatch
        real                     :: reslim
        type(image), allocatable :: rec_imgs(:)
        t_init = tic()
        t_tot  = t_init
        ! make the gridding prepper
        kbwin = b%eorecvols(1)%get_kbwin()
        call gridprep%new(b%img, kbwin, [p%boxpd,p%boxpd,1])
        ! init volumes
        call preprecvols(b, p)
        ! prep rec imgs
        allocate(rec_imgs(MAXIMGBATCHSZ))
        do i=1,MAXIMGBATCHSZ
            call rec_imgs(i)%new([p%boxpd, p%boxpd, 1], p%smpd)
        end do
        ! prep batch imgs
        call prepimgbatch(b, p, MAXIMGBATCHSZ)
        rt_init = toc(t_init)
        ! reconstruction
        t_rec       = tic()
        rt_read     = 0.
        rt_gridprep = 0.
        rt_gridding = 0.
        do iptcl_batch=1,p%nptcls,MAXIMGBATCHSZ
            batchlims = [iptcl_batch,min(p%nptcls,iptcl_batch + MAXIMGBATCHSZ - 1)]
            t_read = tic()
            call read_imgbatch(b, p, batchlims)
            rt_read    = rt_read + toc(t_read)
            t_gridprep = tic()
            !$omp parallel do default(shared) private(iptcl)&
            !$omp schedule(static) proc_bind(close)
            do iptcl=batchlims(1),batchlims(2)
                ibatch = iptcl - batchlims(1) + 1
                call gridprep%prep_serial_no_fft(b%imgbatch(ibatch), rec_imgs(ibatch))
            end do
            !$omp end parallel do
            rt_gridprep = rt_gridprep + toc(t_gridprep)
            t_gridding  = tic()
            do iptcl=batchlims(1),batchlims(2)
                ibatch = iptcl - batchlims(1) + 1
                orientation = b%a%get_ori(iptcl)
                call grid_ptcl_tst(b, p, rec_imgs(ibatch), orientation )
            end do
            rt_gridding = rt_gridding + toc(t_gridding)
        end do
        rt_rec = toc(t_rec)
        ! normalise structure factors
        t_norm = tic()
        call eonorm_struct_facts(b, p, cline, reslim)
        rt_norm = toc(t_norm)
        ! destruct
        call killrecvols(b, p)
        call gridprep%kill
        do ibatch=1,MAXIMGBATCHSZ
            call rec_imgs(ibatch)%kill
            call b%imgbatch(ibatch)%kill
        end do
        deallocate(rec_imgs, b%imgbatch)
        rt_tot = toc(t_tot)
        call gridprep%kill
        benchfname = 'REC_TEST_BATCH_GRIDPREP.txt'
        call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
        write(fnr,'(a)') '*** TIMINGS (s) ***'
        write(fnr,'(a,1x,f9.2)') 'initialisation  : ', rt_init
        write(fnr,'(a,1x,f9.2)') 'reconstruction  : ', rt_rec
        write(fnr,'(a,1x,f9.2)') '* read          : ', rt_read
        write(fnr,'(a,1x,f9.2)') '* gridprep      : ', rt_gridprep
        write(fnr,'(a,1x,f9.2)') '* gridding      : ', rt_gridding
        write(fnr,'(a,1x,f9.2)') 'eonorm          : ', rt_norm
        write(fnr,'(a,1x,f9.2)') 'total time      : ', rt_tot
        write(fnr,'(a)') ''
        write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
        write(fnr,'(a,1x,f9.2)') 'initialisation  : ', (rt_init/rt_tot) * 100.
        write(fnr,'(a,1x,f9.2)') 'reconstruction  : ', (rt_rec/rt_tot)  * 100.
        write(fnr,'(a,1x,f9.2)') '* read          : ', (rt_read/rt_tot)     * 100.
        write(fnr,'(a,1x,f9.2)') '* gridprep      : ', (rt_gridprep/rt_tot) * 100.
        write(fnr,'(a,1x,f9.2)') '* gridding      : ', (rt_gridding/rt_tot) * 100.
        write(fnr,'(a,1x,f9.2)') 'eonorm          : ', (rt_norm/rt_tot) * 100.
        write(fnr,'(a,1x,f9.2)') '% accounted for : ', ((rt_init+rt_rec+rt_norm)/rt_tot) * 100.
        call fclose(fnr)
    end subroutine exec_rec_batch_gridprep

end module simple_reconstructor_tester
