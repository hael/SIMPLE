!@descr: standard EM image processing
module simple_commanders_imgproc
use simple_commander_module_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_ctfops
  contains
    procedure :: execute      => exec_ctfops
end type commander_ctfops

type, extends(commander_base) :: commander_ctf_phaseflip
  contains
    procedure :: execute      => exec_ctf_phaseflip
end type commander_ctf_phaseflip

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

    subroutine exec_ctf_phaseflip( self, cline )
        use simple_strategy2D3D_common, only: read_imgbatch
        use simple_ori,                 only: ori
        use simple_ctf,                 only: ctf
        use simple_memoize_ft_maps,     only: memoize_ft_maps, forget_ft_maps
        class(commander_ctf_phaseflip), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(image), allocatable :: imgs(:)
        type(stack_io)   :: stkio_w
        type(ctf)        :: tfun
        type(ctfparams)  :: ctfparms
        type(parameters) :: params
        type(builder)    :: build
        integer          :: nptcls, ldim(3), iptcl
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'phaseflipped'//STK_EXT)
        call build%init_params_and_build_general_tbox(cline,params)
        if( cline%defined('stk') )then
            call find_ldim_nptcls(params%stk, ldim, nptcls)
            allocate(imgs(nptcls))
            do iptcl = 1, nptcls
                call imgs(iptcl)%new([ldim(1),ldim(2),1],params%smpd)
                call imgs(iptcl)%read(params%stk, iptcl)
            enddo
            call memoize_ft_maps(imgs(1)%get_ldim())
            print *, 'FINISHED READING...'
            !$omp parallel do private(iptcl,ctfparms,tfun) default(shared) proc_bind(close) schedule(static)
            do iptcl = 1, nptcls
                call imgs(iptcl)%fft
                ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
                tfun     = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call imgs(iptcl)%apply_ctf(tfun, 'flip', ctfparms)
                call imgs(iptcl)%ifft
            end do
            !$omp end parallel do
            print *, 'START WRITING...'
            do iptcl = 1, nptcls
                call imgs(iptcl)%write(params%outstk, iptcl)
                call imgs(iptcl)%kill
            enddo
            call forget_ft_maps
        else
            nptcls = build%spproj%get_nptcls()
            ldim   = build%img%get_ldim()
            call memoize_ft_maps(build%img%get_ldim())
            call stkio_w%open(params%outstk, params%smpd, 'write', box=ldim(1))
            do iptcl = 1, nptcls
                call read_imgbatch(iptcl, build%img)
                call build%img%fft
                ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
                tfun     = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call imgs(iptcl)%apply_ctf(tfun, 'flip', ctfparms)
                call build%img%ifft
                call stkio_w%write(iptcl, build%img)
            end do
            call stkio_w%close
            call forget_ft_maps
        endif
        ! cleanup
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_CTF_PHASEFLIP NORMAL STOP ****')
    end subroutine exec_ctf_phaseflip

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
        call automask2D(imgs, 0, 0, params%edge, diams, shifts)
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
