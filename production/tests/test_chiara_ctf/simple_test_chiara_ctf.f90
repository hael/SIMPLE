module construct_ctf
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_ctf,        only: ctf
use simple_image,      only: image
use simple_parameters, only: parameters
implicit none
contains
    function build_ctf(smpd, dfx,dfy,angast, cs,kv, fraca) result(img_ctf)
        integer, parameter     :: box=256
        real                   :: dfx, dfy, angast, smpd, cs, kv, fraca
        type(image)            :: ctf_image,img4viz,img_ctf
        type(ctf)              :: tfun
        type(parameters)       :: params
        integer, dimension(3)  :: ldim
        call ctf_image%new([box,box,1], smpd)
        call img4viz%new([box,box,1],   smpd)
        tfun = ctf(smpd, kv, cs, fraca)
        call tfun%ctf2img(ctf_image, dfx, dfy, angast)
        call ctf_image%ft2img('real', img4viz)
        img_ctf = img4viz
    end function build_ctf
end module construct_ctf

program simple_chiara_ctf
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_ctf,        only: ctf
use simple_image,      only: image
use simple_parameters, only: parameters
use construct_ctf
implicit none
type(cmdline)          :: cline
integer, parameter     :: box=256
real,    parameter     :: dfx=1.54, dfy =1.72, angast=30, smpd1=1.0
!logical, parameter     :: is = .true.
type(image)            :: ctf_image, img4viz,img
type(ctf)              :: tfun
type(parameters)       :: params
integer, dimension(3)  :: ldim

if( command_argument_count() < 1 )then
    write(logfhandle,'(a)',advance='no') 'simple_test_chiara_ctf2 smpd=<sampling distance(in A)>'
    write(logfhandle,'(a)')              ' [kv=<acceleration voltage(in kV){300.}>] [fraca=<fraction &
        & of amplitude{0.1}>] [cs=<spherical aberration constant(in mm){2.7}>]'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('smpd', 1)       !<Set the required variable
call cline%check                     !<check if the required variable is is_present
!Set defaults
if( .not. cline%defined('kv') )    call cline%set('kv',     300.)
if( .not. cline%defined('fraca') ) call cline%set('fraca',   0.1)
if( .not. cline%defined('cs') )    call cline%set('cs',      2.7)
call params%new(cline)              !<read cline parameters
call ctf_image%new([box,box,1], params%smpd)
call img4viz%new  ([box,box,1], params%smpd)
tfun = ctf(params%smpd, params%kv, params%cs, params%fraca)
call tfun%ctf2img(ctf_image, dfx, dfy, angast)
call ctf_image%ft2img('real', img4viz)

call img4viz%write('ctfimg.mrc')
call img4viz%write_jpg('ctfimg.jpg')
img4viz = 1.
call img4viz%write('ones.mrc')

call img4viz%read('ctfimg.mrc')
call img4viz%write('ctfimg2.mrc')

ldim = img4viz%get_ldim()
call img%new([box+10,box+10,1],1.)
img=1.
call img4viz%pad(img)
call img%vis
end program simple_chiara_ctf

!program another_way
!include 'simple_lib.f08'
!use simple_cmdline,    only: cmdline
!use simple_ctf,        only: ctf
!use simple_image,      only: image
!use simple_parameters, only: parameters
!implicit none
!type(cmdline)          :: cline
!integer, parameter     :: box=256
!real,    parameter     :: dfx=1.54, dfy =1.72, angast=30, smpd1=1.0
!logical, parameter     :: is = .true.
!type(image)            :: ctf_image, img4viz
!type(ctf)              :: tfun
!type(parameters)       :: params
!if( command_argument_count() < 1 )then
!    write(logfhandle,'(a)',advance='no') 'simple_test_chiara_ctf2 smpd=<sampling distance(in A)>'
!    write(logfhandle,'(a)')              ' [kv=<acceleration voltage(in kV){300.}>] [fraca=<fraction &
!                                & of amplitude{0.1}>] [cs=<spherical aberration constant(in mm){2.7}>]'
!    stop
!endif
!call cline%parse_oldschool
!call cline%checkvar('smpd', 1)       !<Set the required variable
!call cline%check                     !<check if the required variable is is_present
!Set defaults
!if( .not. cline%defined('kv') )    call cline%set('kv',     300.)
!if( .not. cline%defined('fraca') ) call cline%set('fraca',   0.1)
!if( .not. cline%defined('cs') )    call cline%set('cs',      2.7)
!call params%new(cline)              !<read cline parameters
!call ctf_image%new([box,box,1], params%smpd)
!call img4viz%new([box,box,1],   params%smpd)
!tfun = ctf(params%smpd, params%kv, params%cs, params%fraca)
!call tfun%ctf2img(ctf_image, dfx, dfy, angast)
!if(ctf_image .eqdims. img4viz) then  !useless, just to try
!call ctf_image%ft2img('real', img4viz)
!call img4viz%vis
!end if
!end program another_way
