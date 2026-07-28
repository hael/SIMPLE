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
        use simple_matcher_ptcl_io,     only: read_imgbatch
        use simple_ori,                 only: ori
        use simple_ctf,                 only: ctf
        use simple_estimate_ssnr,       only: fsc2wiener_regularizer
        use simple_fileio,              only: file2rarr
        use simple_memoize_ft_maps,     only: memoize_ft_maps, forget_ft_maps
        class(commander_ctf_correct),   intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(image), allocatable :: imgs(:)
        type(stack_io)   :: stkio_w
        type(ctf)        :: tfun
        type(ctfparams)  :: ctfparms
        type(parameters) :: params
        type(builder)    :: build
        real, allocatable :: fsc(:), noise_to_signal(:), noise_to_signal_pad(:)
        integer          :: nptcls, ldim(3), ldim_wpad(3), iptcl, ndone, report_freq, wiener_nobs
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'ctf_corrected'//STK_EXT)
        call build%init_params_and_build_general_tbox(cline,params)
        select case(trim(params%ctf_correct_mode))
            case('phaseflip', 'wiener')
                continue
            case DEFAULT
                THROW_HARD('ctf_correct_mode must be phaseflip or wiener')
        end select
        if( cline%defined('stk') )then
            call find_ldim_nptcls(params%stk, ldim, nptcls)
            allocate(imgs(nptcls))
            do iptcl = 1, nptcls
                call imgs(iptcl)%new([ldim(1),ldim(2),1],params%smpd)
                call imgs(iptcl)%read(params%stk, iptcl)
            enddo
            if( trim(params%ctf_correct_mode) == 'wiener' )then
                ldim_wpad    = 2 * imgs(1)%get_ldim()
                ldim_wpad(3) = 1
                call memoize_ft_maps(ldim_wpad, imgs(1)%get_smpd())
            else
                call memoize_ft_maps(imgs(1)%get_ldim(), imgs(1)%get_smpd())
            endif
            call correct_particles()
            do iptcl = 1, nptcls
                call imgs(iptcl)%write(params%outstk, iptcl)
                call imgs(iptcl)%kill
            enddo
            call forget_ft_maps
        else
            nptcls = build%spproj%get_nptcls()
            ldim   = build%img%get_ldim()
            allocate(imgs(nptcls))
            do iptcl = 1, nptcls
                call imgs(iptcl)%new([ldim(1),ldim(2),1], params%smpd)
                call read_imgbatch(params, build, iptcl, imgs(iptcl))
            enddo
            if( trim(params%ctf_correct_mode) == 'wiener' )then
                ldim_wpad    = 2 * imgs(1)%get_ldim()
                ldim_wpad(3) = 1
                call memoize_ft_maps(ldim_wpad, imgs(1)%get_smpd())
            else
                call memoize_ft_maps(ldim, params%smpd)
            endif
            call correct_particles()
            call stkio_w%open(params%outstk, params%smpd, 'write', box=ldim(1))
            do iptcl = 1, nptcls
                call stkio_w%write(iptcl, imgs(iptcl))
                call imgs(iptcl)%kill
            enddo
            call stkio_w%close
            call forget_ft_maps
        endif
        ! cleanup
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_CTF_CORRECT NORMAL STOP ****')

    contains

        subroutine correct_particles()
            integer :: k
            if( trim(params%ctf_correct_mode) == 'wiener' .and. cline%defined('fsc') )then
                wiener_nobs = build%spproj%count_state_gt_zero()
                if( wiener_nobs < 1 ) THROW_HARD('no active particles in project; ctf_correct')
                fsc = file2rarr(params%fsc)
                noise_to_signal = fsc2wiener_regularizer(fsc, wiener_nobs)
                allocate(noise_to_signal_pad(2 * imgs(1)%get_filtsz()), source=0.)
                do k = 1, size(noise_to_signal_pad)
                    noise_to_signal_pad(k) = noise_to_signal(min(size(noise_to_signal), max(1,nint(real(k) / 2.))))
                end do
                write(logfhandle,'(A,I0,A)') '>>> USING FSC-DERIVED WIENER REGULARIZATION FROM ', wiener_nobs, ' ACTIVE PARTICLES'
            endif
            write(logfhandle,'(A,I0,A)') '>>> CTF-CORRECTING ', nptcls, ' PARTICLES'
            call flush(logfhandle)
            ndone       = 0
            report_freq = max(1, nptcls / 10)
            !$omp parallel do private(iptcl,ctfparms,tfun) default(shared) proc_bind(close) schedule(static)
            do iptcl = 1, nptcls
                ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
                tfun     = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                select case(trim(params%ctf_correct_mode))
                    case('phaseflip')
                        call imgs(iptcl)%fft
                        call imgs(iptcl)%apply_ctf(tfun, 'flip', ctfparms)
                        call imgs(iptcl)%ifft
                    case('wiener')
                        if( allocated(noise_to_signal_pad) )then
                            call imgs(iptcl)%apply_ctf_wiener_wpad(tfun, ctfparms, params%wiener_const, noise_to_signal_pad)
                        else
                            call imgs(iptcl)%apply_ctf_wiener_wpad(tfun, ctfparms, params%wiener_const)
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
        end subroutine correct_particles
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
        if( .not. cline%defined('amsklp')  ) call cline%set('amsklp', cline%get_rarg('lp'))
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
