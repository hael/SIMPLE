!@descr: stuff that didn't fit elsewhere
module simple_commanders_misc
use simple_commanders_api
use simple_simple_volinterp, only: rotvol
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_print_fsc
  contains
    procedure :: execute       => exec_print_fsc
end type commander_print_fsc

type, extends(commander_base) :: commander_print_magic_boxes
  contains
    procedure :: execute       => exec_print_magic_boxes
end type commander_print_magic_boxes

type, extends(commander_base) :: commander_print_dose_weights
  contains
    procedure :: execute       => exec_print_dose_weights
end type commander_print_dose_weights

type, extends(commander_base) :: commander_kstest
  contains
    procedure :: execute       => exec_kstest
end type commander_kstest

type, extends(commander_base) :: commander_pearsn
  contains
    procedure :: execute       => exec_pearsn
end type commander_pearsn

type, extends(commander_base) :: commander_mkdir
  contains
    procedure :: execute       => exec_mkdir
end type commander_mkdir

contains

    !>  for printing the binary FSC files produced by PRIME3D
    subroutine exec_print_fsc( self, cline )
        use simple_fsc,        only: plot_fsc
        use simple_class_frcs, only: class_frcs
        class(commander_print_fsc), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        real,  allocatable :: res(:), fsc(:)
        type(parameters) :: params
        type(image)      :: img
        type(class_frcs) :: frcs
        type(string)     :: tmpl_fname
        integer          :: k,n
        real             :: res0143, res05
        logical :: l_fsc, l_frcs
        call params%new(cline)
        l_fsc  = cline%defined('fsc')
        l_frcs = cline%defined('frcs')
        if( .not.l_fsc .and. .not.l_frcs ) THROW_HARD('FSC or FRCS must be defined!')
        if( l_fsc )then
            if( .not.cline%defined('smpd') .or. .not.cline%defined('box')) then
                THROW_HARD('SMPD and BOX must be defined!')
            endif
            call img%new([params%box,params%box,1], params%smpd)
            res = img%get_res()
            fsc = file2rarr(params%fsc)
            n = size(fsc)
            do k=1,n
                write(logfhandle,'(A,1X,F6.2,1X,A,1X,F15.3)') '>>> RESOLUTION:', res(k), '>>> FSC:', fsc(k)
            end do
            ! get & print resolution
            call get_resolution(fsc, res, res05, res0143)
            write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res0143
            write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res05
            ! plot
            tmpl_fname = get_fbody(params%fsc,BIN_EXT,separator=.false.)
            call plot_fsc(n, fsc, res, params%smpd, tmpl_fname%to_char())
            call img%kill
        endif
        if( l_frcs )then
            call frcs%read(params%frcs)
            call frcs%plot_frcs(string('frcs'))
            call frcs%kill
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PRINT_FSC NORMAL STOP ****')
    end subroutine exec_print_fsc

    !> for printing magic box sizes (fast FFT)
    subroutine exec_print_magic_boxes( self, cline )
        class(commander_print_magic_boxes), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters) :: params
        call params%new(cline)
        call print_magic_box_range(params%smpd, params%moldiam )
        ! end gracefully
        call simple_end('**** SIMPLE_PRINT_MAGIC_BOXES NORMAL STOP ****')
    end subroutine exec_print_magic_boxes

    subroutine exec_print_dose_weights( self, cline )
        class(commander_print_dose_weights), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters)  :: params
        real, allocatable :: weights(:,:), res(:)
        integer           :: iframe, k, filtsz
        call params%new(cline)
        call calc_dose_weights(params%nframes, params%box, params%smpd, params%kV, params%total_dose, weights)
        filtsz = size(weights, dim=2)
        res = get_resarr(params%box, params%smpd)
        write(logfhandle,'(A)') 'RESOLUTION, DOSE_WEIGHTS'
        do k = 1,filtsz
            write(logfhandle, '(F7.1,A)', advance='no') res(k), ', '
            do iframe = 1,params%nframes - 1
                write(logfhandle, '(f3.1,A)', advance='no') weights(iframe,k), ', '
            end do
            write(logfhandle, '(f3.1,1X)') weights(iframe,k)
        end do
        write(logfhandle,*)
        ! end gracefully
        call simple_end('**** SIMPLE_PRINT_DOSE_WEIGHTS_NORMAL STOP ****')
    end subroutine exec_print_dose_weights

    subroutine exec_kstest( self, cline )
        class(commander_kstest), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)   :: params
        integer            :: ndat1, ndat2
        real               :: ksstat, prob, ave1, sdev1, var, ave2, sdev2
        real, allocatable  :: dat1(:), dat2(:)
        logical            :: err
        call params%new(cline)
        call read_nrs_dat(params%infile,  dat1, ndat1)
        call read_nrs_dat(params%infile2, dat2, ndat2)
        write(logfhandle,'(a)') '>>> STATISTICS OF THE TWO DISTRIBUTIONS'
        call moment(dat1, ave1, sdev1, var, err)
        call moment(dat2, ave2, sdev2, var, err)
        write(logfhandle,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile : ', ave1, sdev1
        write(logfhandle,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile2: ', ave2, sdev2
        write(logfhandle,'(a)') '>>> KOLMOGOROV-SMIRNOV TEST TO DEDUCE EQUIVALENCE OR NON-EQUIVALENCE BETWEEN TWO DISTRIBUTIONS'
        call kstwo(dat1, ndat1, dat2, ndat2, ksstat, prob)
        write(logfhandle,'(a,1x,f4.2)') 'K-S statistic = ', ksstat
        write(logfhandle,'(a,1x,f4.2)') 'P             = ', prob
        write(logfhandle,'(a)') 'P represents the significance level for the null hypothesis that the two data sets are drawn from the same distribution'
        write(logfhandle,'(a)') 'Small P values show that the cumulative distribution functions of the two data sets differ significantly'
        ! end gracefully
        call simple_end('**** SIMPLE_KSTEST NORMAL STOP ****')
    end subroutine exec_kstest

    subroutine exec_pearsn( self, cline )
        class(commander_pearsn), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)   :: params
        integer            :: ndat1, ndat2
        real               :: corr, ave1, sdev1, var, ave2, sdev2
        real, allocatable  :: dat1(:), dat2(:)
        logical            :: err
        call params%new(cline)
        call read_nrs_dat(params%infile,  dat1, ndat1)
        call read_nrs_dat(params%infile2, dat2, ndat2)
        if( ndat1 /= ndat2 ) THROW_HARD('Input distributions not identical')
        call moment(dat1, ave1, sdev1, var, err)
        call moment(dat2, ave2, sdev2, var, err)
        write(logfhandle,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile : ', ave1, sdev1
        write(logfhandle,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile2: ', ave2, sdev2
        write(logfhandle,'(a)') '>>> PEARSON CORRELATION OF THE TWO DISTRIBUTIONS'
        corr = pearsn(dat1, dat2)
        write(logfhandle,'(a,1x,f4.2)') 'CORR = ', corr
        write(logfhandle,'(a)') 'P represents the significance level for the null hypothesis that the two data sets are drawn from the same distribution'
        write(logfhandle,'(a)') 'Small P values show that the cumulative distribution functions of the two data sets differ significantly'
        ! end gracefully
        call simple_end('**** SIMPLE_PEARSN NORMAL STOP ****')
    end subroutine exec_pearsn

    subroutine exec_mkdir( self, cline )
        class(commander_mkdir), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        call cline%set('mkdir', 'yes')
        call params%new(cline)
    end subroutine exec_mkdir

    ! utility routines

    subroutine read_nrs_dat( filename, arr, ndat )
        class(string),     intent(in)  :: filename
        real, allocatable, intent(out) :: arr(:)
        integer,           intent(out) :: ndat
        integer            :: ndatlines, nrecs, i, j, cnt
        real, allocatable  :: line(:)
        type(nrtxtfile)    :: nrsfile
        call nrsfile%new(filename, 1)
        ndatlines = nrsfile%get_ndatalines()
        nrecs     = nrsfile%get_nrecs_per_line()
        ndat = ndatlines * nrecs
        allocate( line(nrecs), arr(ndat) )
        cnt = 0
        do i=1,ndatlines
            call nrsfile%readNextDataLine(line)
            do j=1,nrecs
                cnt = cnt + 1
                arr(cnt) = line(j)
            end do
        end do
        deallocate (line)
        call nrsfile%kill
    end subroutine read_nrs_dat

end module simple_commanders_misc
