program simple_map3dshift22d
include 'simple_lib.f08'
use simple_image,          only: image
use simple_oris,           only: oris
use simple_projector_hlev, only: reproject
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline, cmdline_err
use simple_parameters,     only: parameters, params_glob
use simple_rec_master,     only: exec_rec

implicit none

character(len=KEYLEN) :: keys_required(MAXNKEYS)='', keys_optional(MAXNKEYS)=''
character(len=STDLEN) :: xarg, prg, entire_line
type(cmdline)         :: cline
type(parameters)      :: params
type(builder)         :: build
type(ctfparams)       :: ctfparms
integer               :: cmdstat, cmdlen, pos
type(oris)  :: ref_os,os
type(image), allocatable :: imgs(:)
real        :: shvec(3)
integer     :: i

! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
keys_required(1) = 'vol1'
keys_required(2) = 'smpd'
call cline%set('nptcls',100.)
call cline%set('eo','no')
call cline%set('lp','4.')
call cline%parse_private
call cline%set('oritype','ptcl3D')
call build%init_params_and_build_general_tbox(cline, params)

params_glob%eo = 'no'

call build%vol%read(params%vols(1))
call ref_os%new(100)
call ref_os%spiral
call ref_os%set_all2single('state',1.)
call ref_os%set_all2single('w',1.)
call ref_os%write('reforis.txt')

ctfparms%smpd = params%smpd
call build%spproj%os_ptcl3D%kill
call build%spproj%add_single_stk('imgs.mrc', ctfparms, ref_os)


imgs = reproject(build%vol, ref_os)
do i=1,100
    call imgs(i)%write('imgs.mrc',i)
enddo


shvec = [-5.,10.,-15.]
call build%vol%shift(shvec)
call build%vol%write('shifted.mrc')
call build%spproj_field%map3dshift22d(-shvec)
call build%spproj_field%write('shifted.txt')
call build%build_rec_eo_tbox(params)    ! reconstruction objects built
call exec_rec



! shift back
call build%spproj_field%map3dshift22d(shvec)
call build%spproj_field%write('backshifted.txt')


end program simple_map3dshift22d
