!@descr: standard EM image processing
module simple_commanders_imgproc
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_ctfops
  contains
    procedure :: execute      => exec_ctfops
end type commander_ctfops

type, extends(commander_base) :: commander_ctf_correct
  contains
    procedure :: execute      => exec_ctf_correct
end type commander_ctf_correct

type, extends(commander_base) :: commander_estimate_diam
  contains
    procedure :: execute      => exec_estimate_diam
end type commander_estimate_diam

contains

    !> for applying CTF to stacked images
    subroutine exec_ctfops( self, cline )
        use simple_procimgstk, only: apply_ctf_imgfile
        class(commander_ctfops), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        if( cline%defined('oritab') .or. cline%defined('deftab') )then
        else
            THROW_HARD('oritab/deftab with CTF info needed for phase flipping/multiplication/CTF image generation')
        endif
        if( .not. cline%defined('stk') ) call cline%set('box', 256.)
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        if( params%ctf .ne. 'no' )then
            select case( params%ctf )
                case( 'flip' )
                    if( params%neg .eq. 'yes' )then
                        call apply_ctf_imgfile(params%stk, params%outstk, build%spproj_field, params%smpd, 'flipneg')
                    else
                        call apply_ctf_imgfile(params%stk, params%outstk, build%spproj_field, params%smpd, 'flip')
                    endif
                case( 'yes' )
                    if( params%neg .eq. 'yes' )then
                        call apply_ctf_imgfile(params%stk, params%outstk, build%spproj_field, params%smpd, 'neg')
                    else
                        call apply_ctf_imgfile(params%stk, params%outstk, build%spproj_field, params%smpd, 'ctf')
                    endif
                case DEFAULT
                    THROW_HARD('Unknown ctf argument')
            end select
        else
            THROW_HARD('Nothing to do!')
        endif
        ! cleanup
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_CTFOPS NORMAL STOP ****')
    end subroutine exec_ctfops

    subroutine exec_ctf_correct( self, cline )
        use simple_matcher_ptcl_io,     only: prepimgbatch, killimgbatch, read_imgbatch
        use simple_ori,                 only: ori
        use simple_ctf,                 only: ctf
        use simple_estimate_ssnr,       only: fsc2wiener_regularizer
        use simple_fileio,              only: file2rarr
        use simple_memoize_ft_maps,     only: memoize_ft_maps, forget_ft_maps
        class(commander_ctf_correct),   intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(image), allocatable :: imgs(:)
        type(stack_io)   :: stkio_r, stkio_w
        type(ctf)        :: tfun
        type(ctfparams)  :: ctfparms
        type(parameters) :: params
        type(builder)    :: build
        real, allocatable :: fsc(:), noise_to_signal(:), noise_to_signal_pad(:)
        integer, allocatable :: pinds(:)
        integer          :: nptcls, ldim(3), ldim_wpad(3), iptcl, ndone, report_freq, wiener_nobs
        integer          :: batchsz, batch_from, batch_to, n_batch, i, ibatch, nbatches
        logical          :: use_stk_input, heuristic_logged
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'ctf_corrected'//STK_EXT)
        write(logfhandle,'(A)') '>>> CTF-CORRECT INITIALIZING'
        call flush(logfhandle)
        call build%init_params_and_build_general_tbox(cline,params)
        select case(trim(params%ctf_correct_mode))
            case('phaseflip', 'wiener')
                ! all good
            case DEFAULT
                THROW_HARD('ctf_correct_mode must be phaseflip or wiener')
        end select
        use_stk_input = cline%defined('stk')
        if( use_stk_input )then
            call find_ldim_nptcls(params%stk, ldim, nptcls)
        else
            nptcls = build%spproj%get_nptcls()
            ldim   = build%img%get_ldim()
        endif
        if( nptcls < 1 ) THROW_HARD('no particles found; ctf_correct')
        if( trim(params%ctf_correct_mode) == 'wiener' )then
            ldim_wpad    = 2 * ldim
            ldim_wpad(3) = 1
            call memoize_ft_maps(ldim_wpad, params%smpd)
        else
            call memoize_ft_maps(ldim, params%smpd)
        endif
        batchsz          = max(1, min(nptcls, max(1, params%nthr) * BATCHTHRSZ))
        nbatches         = ceiling(real(nptcls) / real(batchsz))
        if( .not. use_stk_input )then
            allocate(pinds(nptcls), source=(/(i,i=1,nptcls)/))
            call prepimgbatch(params, build, batchsz)
        else
            call stkio_r%open(params%stk, params%smpd, 'read')
        endif
        ndone            = 0
        report_freq      = max(1, nptcls / 10)
        heuristic_logged = .false.
        write(logfhandle,'(A,I0,A,I0,A)') '>>> CTF-CORRECTING ', nptcls, ' PARTICLES IN ', nbatches, ' BATCHES'
        call flush(logfhandle)
        call stkio_w%open(params%outstk, params%smpd, 'write', box=ldim(1))
        ibatch = 0
        do batch_from = 1, nptcls, batchsz
            ibatch   = ibatch + 1
            batch_to = min(nptcls, batch_from + batchsz - 1)
            n_batch  = batch_to - batch_from + 1
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> BATCH ', ibatch, '/', nbatches, ': rows ', batch_from, '-', batch_to
            call flush(logfhandle)
            allocate(imgs(n_batch))
            if( .not. use_stk_input )then
                call read_imgbatch(params, build, nptcls, pinds, [batch_from,batch_to])
            endif
            do i = 1, n_batch
                iptcl = batch_from + i - 1
                call imgs(i)%new([ldim(1),ldim(2),1], params%smpd)
                if( use_stk_input )then
                    call stkio_r%read(iptcl, imgs(i))
                else
                    call imgs(i)%copy(build%imgbatch(i))
                endif
            enddo
            call correct_particles_batch(batch_from, n_batch, heuristic_logged)
            do i = 1, n_batch
                iptcl = batch_from + i - 1
                call stkio_w%write(iptcl, imgs(i))
                call imgs(i)%kill
            enddo
            deallocate(imgs)
        enddo
        if( use_stk_input )then
            call stkio_r%close
        else
            call killimgbatch(build)
            deallocate(pinds)
        endif
        call stkio_w%close
        call forget_ft_maps
        call update_local_project_copy()
        ! cleanup
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_CTF_CORRECT NORMAL STOP ****')

    contains

        subroutine correct_particles_batch( batch_from, n_batch, heuristic_logged )
            integer, intent(in) :: batch_from, n_batch
            logical, intent(inout) :: heuristic_logged
            integer :: k, i
            if( trim(params%ctf_correct_mode) == 'wiener' .and. cline%defined('fsc') )then
                if( .not. allocated(noise_to_signal_pad) )then
                    wiener_nobs = build%spproj%count_state_gt_zero()
                    if( wiener_nobs < 1 ) THROW_HARD('no active particles in project; ctf_correct')
                    fsc = file2rarr(params%fsc)
                    noise_to_signal = fsc2wiener_regularizer(fsc, wiener_nobs)
                    allocate(noise_to_signal_pad(2 * imgs(1)%get_filtsz()), source=0.)
                    do k = 1, size(noise_to_signal_pad)
                        noise_to_signal_pad(k) = noise_to_signal(min(size(noise_to_signal), max(1,nint(real(k) / 2.))))
                    end do
                    write(logfhandle,'(A,I0,A)') '>>> USING FSC-DERIVED WIENER REGULARIZATION FROM ', wiener_nobs, ' ACTIVE PARTICLES'
                    call flush(logfhandle)
                endif
            else if( trim(params%ctf_correct_mode) == 'wiener' )then
                if( .not. heuristic_logged )then
                    write(logfhandle,'(A)') '>>> USING GRIGORIEFF HEURISTIC WIENER REGULARIZATION (10 PERCENT OF MEAN CTF SQUARED)'
                    call flush(logfhandle)
                    heuristic_logged = .true.
                endif
            endif
            !$omp parallel do private(i,iptcl,ctfparms,tfun) default(shared) proc_bind(close) schedule(static)
            do i = 1, n_batch
                iptcl = batch_from + i - 1
                ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
                tfun     = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                select case(trim(params%ctf_correct_mode))
                    case('phaseflip')
                        call imgs(i)%fft
                        call imgs(i)%apply_ctf(tfun, 'flip', ctfparms)
                        call imgs(i)%ifft
                    case('wiener')
                        if( allocated(noise_to_signal_pad) )then
                            call imgs(i)%apply_ctf_wiener_wpad(tfun, ctfparms, 0., noise_to_signal_pad)
                        else
                            call imgs(i)%apply_ctf_wiener_wpad(tfun, ctfparms, 0.)
                        endif
                end select
                !$omp critical(ctf_correct_progress)
                ndone = ndone + 1
                if( mod(ndone, report_freq) == 0 .or. ndone == nptcls )then
                    write(logfhandle,'(A,I3,A)') '>>> ', nint(100.0 * real(ndone) / real(nptcls)), '% done'
                    call flush(logfhandle)
                endif
                !$omp end critical(ctf_correct_progress)
            end do
            !$omp end parallel do
        end subroutine correct_particles_batch

        subroutine update_local_project_copy()
            type(sp_project) :: spproj_out
            type(oris)       :: os_ptcl2D_tmp, os_ptcl3D_tmp
            type(ctfparams)  :: ctfparms_out
            integer :: nptcls_out, nptcls_proj, nstks, iptcl
            integer :: ldim_out(3)
            call spproj_out%read(params%projfile)
            nptcls_proj = spproj_out%os_ptcl2D%get_noris()
            if( nptcls_proj < 1 ) THROW_HARD('project has no ptcl2D entries; ctf_correct')
            call find_ldim_nptcls(params%outstk, ldim_out, nptcls_out)
            if( nptcls_out /= nptcls_proj )then
                write(logfhandle,*) 'nptcls in project ptcl2D : ', nptcls_proj
                write(logfhandle,*) 'nptcls in output stack   : ', nptcls_out
                THROW_HARD('output stack does not match project particle count; ctf_correct')
            endif
            nstks = spproj_out%os_stk%get_noris()
            if( nstks < 1 ) THROW_HARD('project has no stack metadata; ctf_correct')
            ctfparms_out = spproj_out%get_ctfparams('stk', 1)
            select case(trim(params%ctf_correct_mode))
                case('phaseflip')
                    ctfparms_out%ctfflag = CTFFLAG_FLIP
                case('wiener')
                    ctfparms_out%ctfflag = CTFFLAG_NO
                case DEFAULT
                    THROW_HARD('unsupported ctf_correct_mode for project update; ctf_correct')
            end select
            os_ptcl2D_tmp = spproj_out%os_ptcl2D
            os_ptcl3D_tmp = spproj_out%os_ptcl3D
            call spproj_out%os_stk%kill
            call spproj_out%os_ptcl2D%kill
            call spproj_out%os_ptcl3D%kill
            call spproj_out%add_stk(params%outstk, ctfparms_out)
            spproj_out%os_ptcl2D = os_ptcl2D_tmp
            spproj_out%os_ptcl3D = os_ptcl3D_tmp
            call os_ptcl2D_tmp%kill
            call os_ptcl3D_tmp%kill
            do iptcl = 1,spproj_out%os_ptcl2D%get_noris()
                call spproj_out%os_ptcl2D%set(iptcl, 'stkind', 1)
                call spproj_out%os_ptcl2D%set(iptcl, 'indstk', iptcl)
            enddo
            do iptcl = 1,spproj_out%os_ptcl3D%get_noris()
                call spproj_out%os_ptcl3D%set(iptcl, 'stkind', 1)
                call spproj_out%os_ptcl3D%set(iptcl, 'indstk', iptcl)
            enddo
            call spproj_out%write(params%projfile)
            write(logfhandle,'(A,A)') '>>> UPDATED LOCAL PROJECT COPY: ', params%projfile%to_char()
            call spproj_out%kill
        end subroutine update_local_project_copy
    end subroutine exec_ctf_correct

    subroutine exec_estimate_diam( self, cline )
        use simple_segmentation
        class(commander_estimate_diam), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        ! constants
        character(len=*), parameter :: FILT   = 'filtered.mrc'
        character(len=*), parameter :: BIN    = 'binarized.mrc'
        character(len=*), parameter :: MASKED = 'masked.mrc'
        ! varables
        type(parameters)            :: params
        type(image),    allocatable :: imgs(:) 
        type(stats_struct)          :: diamstats    ! stats struct
        real,           allocatable :: diams(:), diams_nonzero(:), shifts(:,:)
        integer :: i, funit
        real    :: med_diam
        if( .not. cline%defined('lp')      ) call cline%set('lp',        7.0)
        if( .not. cline%defined('automsk') ) call cline%set('automsk', 'yes')
        if( .not. cline%defined('amsklp') )then
            call cline%set('amsklp', cline%get_rarg('lp'))
        else if( cline%get_rarg('amsklp') <= 0. )then
            call cline%set('amsklp', cline%get_rarg('lp'))
        endif
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        call params%new(cline)
        ! allocate & read cavgs
        allocate(imgs(params%nptcls))
        do i=1,params%nptcls
            call imgs(i)%new([params%box,params%box,1],params%smpd)
            call imgs(i)%read(params%stk, i)
        end do
        call automask2D(params, imgs, 0, 0, params%edge, diams, shifts)
        diams_nonzero = pack(diams, mask=diams > TINY)
        call calc_stats(diams_nonzero, diamstats)
        ! output
        med_diam = median(diams)
        write(logfhandle,'(A,2F6.1)') '>>> AVG    DIAMETER (IN A & pix): ', diamstats%avg,  diamstats%avg/params%smpd
        write(logfhandle,'(A,2F6.1)') '>>> SDEV   DIAMETER (IN A & pix): ', diamstats%sdev, diamstats%sdev/params%smpd
        write(logfhandle,'(A,2F6.1)') '>>> MEDIAN DIAMETER (IN A & pix): ', med_diam,       med_diam/params%smpd
        write(logfhandle,'(A,2F6.1)') '>>> MAX    DIAMETER (IN A & pix): ', diamstats%maxv, diamstats%maxv/params%smpd
        write(logfhandle,'(A,2F6.1)') '>>> MIN    DIAMETER (IN A & pix): ', diamstats%minv, diamstats%minv/params%smpd
        call fopen(funit, file=string('diameter_stats.txt'), status='replace')
        write(funit,     '(A,2F6.1)') '>>> AVG    DIAMETER (IN A & pix): ', diamstats%avg,  diamstats%avg/params%smpd
        write(funit,     '(A,2F6.1)') '>>> SDEV   DIAMETER (IN A & pix): ', diamstats%sdev, diamstats%sdev/params%smpd
        write(funit,     '(A,2F6.1)') '>>> MEDIAN DIAMETER (IN A & pix): ', med_diam,       med_diam/params%smpd
        write(funit,     '(A,2F6.1)') '>>> MAX    DIAMETER (IN A & pix): ', diamstats%maxv, diamstats%maxv/params%smpd
        write(funit,     '(A,2F6.1)') '>>> MIN    DIAMETER (IN A & pix): ', diamstats%minv, diamstats%minv/params%smpd
        call fclose(funit)
        ! output the minimum and maximum diameter value in the command line object
        call cline%set('min_diam', diamstats%minv)
        call cline%set('max_diam', diamstats%maxv)
        ! destruct
        do i=1,size(imgs)
            call imgs(i)%kill
        end do
        deallocate(imgs)
        ! end gracefully
        call simple_end('**** SIMPLE_ESTIMATE_DIAM NORMAL STOP ****')
    end subroutine exec_estimate_diam

end module simple_commanders_imgproc
