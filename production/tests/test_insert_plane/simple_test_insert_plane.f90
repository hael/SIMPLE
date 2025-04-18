program simple_test_insert_plane
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_cmdline,           only: cmdline
use simple_builder,           only: builder
use simple_image,             only: image
use simple_parameters,        only: parameters
use simple_fplane
use simple_ori
implicit none
type(cmdline)             :: cline
type(builder)             :: b
type(parameters)          :: p
type(ori)                 :: o
type(image)               :: img
type(ctfparams)           :: ctfparms
type(fplane), allocatable :: fpls(:)
real,         allocatable :: rhoexp_singlethread(:,:,:), rhoexp_multithread(:,:,:)
real,         pointer     :: ptr(:,:,:)
real                      :: error
integer                   :: i,j, nvoxels, nerrors, stride
integer, parameter :: BOX      = 64
integer, parameter :: NREPEATS = 50

if( command_argument_count() < 2 )then
    stop 'ERROR! execution is: simple_test_insert_plane prg=dummy nthr=NTHR'
endif
call cline%parse_oldschool
if( .not.cline%defined('nthr') )then
    stop 'ERROR! execution is: simple_test_insert_plane prg=dummy nthr=NTHR'
endif
call cline%set('ml_reg',  'no')
call cline%set('objfun',  'cc')
call cline%set('ctf',     'no')
call cline%set('pgrp',    'c1')
call cline%set('oritype', 'ptcl3D')
call cline%set('smpd',     1.0)
call cline%set('smpd_crop',1.0)
if( .not.cline%defined('box')    )call cline%set('box',    BOX)
if( .not.cline%defined('nptcls') )call cline%set('nptcls', 10)
call b%init_params_and_build_strategy3D_tbox(cline, p)

! dummy inits
call b%spproj%os_stk%new(1,.false.)
call b%spproj%os_stk%set(1,'fromp',1)
call b%spproj%os_stk%set(1,'top',p%nptcls)
call b%spproj%os_ptcl3D%new(p%nptcls,.true.)
call b%spproj%os_ptcl3D%set_all2single('stkind', 1.0)
ctfparms%ctfflag = CTFFLAG_NO
ctfparms%smpd    = p%smpd

! prep reconstruction objects
call b%eorecvols(1)%new(b%spproj, expand=.true.)
allocate(fpls(p%nptcls))
call img%new([p%box,p%box,1], p%smpd)
do i = 1,p%nptcls
    call img%zero_and_unflag_ft
    call img%gauran(0.,1.)
    call img%fft
    call fpls(i)%new(img)
    call fpls(i)%gen_planes(img, ctfparms, [ran3(), ran3()], i)
enddo

call o%new_ori(.true.)
do stride = 1,16

    error   = 0.0
    nerrors = 0
    do j = 1,NREPEATS

        ! Single threaded reconstruction, guaranteed no race condition
        call b%eorecvols(1)%reset_all
        call omp_set_num_threads(1)         !!!
        do i = 1,p%nptcls
            call o%rnd_euler
            call b%spproj%os_ptcl3D%set_ori(i, o)
            call b%eorecvols(1)%test_grid_plane(b%pgrpsyms, o, fpls(i), 0, 1.0, 1)
        enddo
        call omp_set_num_threads(p%nthr)    !!!
        call b%eorecvols(1)%get_rhoexp_ptr('even',ptr)
        rhoexp_singlethread = ptr(:,:,:)
        ! number of voxels involved in interpolation
        nvoxels = count(rhoexp_singlethread>0.5)

        ! Multi threaded
        call b%eorecvols(1)%reset_all
        do i = 1,p%nptcls
            call b%spproj%os_ptcl3D%get_ori(i, o)
            call b%eorecvols(1)%test_grid_plane(b%pgrpsyms, o, fpls(i), 0, 1.0, stride)
        enddo
        call b%eorecvols(1)%get_rhoexp_ptr('even',ptr)
        rhoexp_multithread = ptr(:,:,:)

        ! number of voxels differently updated between single and multi-threaded
        nerrors = count(abs(rhoexp_multithread-rhoexp_singlethread)>0.5)
        error   = error + real(nerrors) / real(nvoxels)
        deallocate(rhoexp_singlethread,rhoexp_multithread)
    enddo
    error = 100.0*error/real(p%nptcls)/real(NREPEATS)
    print *,'box stride %error ',p%box, stride, error
enddo

end program simple_test_insert_plane