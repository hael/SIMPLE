program simple_test_simnano
include 'simple_lib.f08'
use simple_oris,       only: oris
use simple_ori,        only: ori
use simple_image,      only: image
use simple_kbinterpol, only: kbinterpol
use simple_projector,  only: projector
use simple_ctf,        only: ctf
use simple_simulator,  only: simimg
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_commander_sim, only: simulate_atoms_commander
implicit none
type(simulate_atoms_commander) :: xsim_atoms
type(parameters) :: params
type(cmdline)    :: cline, cline_graphene, cline_particle
type(image)      :: graphene, graphene_vol, particle_vol, particle, img
type(ori)        :: orientation
type(oris)       :: spiral
type(ctf)        :: tfun
type(projector)  :: vol_pad
character(len=:), allocatable :: path
real             :: snr_pink, snr_detector, x,y, crosssection_C, crosssection, scalefactor
integer          :: i,envstat,pad
character(len=LONGSTRLEN), parameter :: graphene_fname = 'sheet.mrc'
character(len=LONGSTRLEN), parameter :: particle_fname = 'ptcl.mrc'
call cline%checkvar('smpd',      1)
call cline%checkvar('nptcls',    2)
call cline%checkvar('snr',       3)
call cline%checkvar('nthr',      4)
call cline%checkvar('bfac',      5)
call cline%checkvar('moldiam',   6)
call cline%checkvar('element',   7)
call cline%parse_oldschool
call cline%check
call cline%set('prg','simnano')
call cline%set('ctf','yes')
call cline%set('box', 512.)
call params%new(cline)
! init
snr_pink     = params%snr/0.2
snr_detector = params%snr/0.8
call spiral%new(params%nptcls)
call spiral%spiral
! simulate graphene
path = simple_getenv('SIMPLE_PATH',envstat)
path = trim(path)//'/../production/tests/test_simnano/graphene_trans.pdb'
call cline_graphene%set('mkdir','no')
call cline_graphene%set('smpd',params%smpd)
call cline_graphene%set('pdbfile', trim(path))
call cline_graphene%set('outvol',  graphene_fname)
call cline_graphene%set('box',real(params%box))
call cline_graphene%set('nthr',real(params%nthr))
call xsim_atoms%execute(cline_graphene)
! simulate nano-particle
call cline_particle%set('mkdir','no')
call cline_particle%set('smpd',params%smpd)
call cline_particle%set('element', trim(params%element))
call cline_particle%set('outvol',  trim(particle_fname))
call cline_particle%set('moldiam', real(params%moldiam))
call cline_particle%set('box',real(params%box))
call cline_particle%set('nthr',real(params%nthr))
call xsim_atoms%execute(cline_particle)
! cross sections
crosssection_C = 1.79 !(x 2.8e-21 m2)
crosssection   = crosssection_C
select case(uppercase(trim(params%element)))
case('FE')
    crosssection = 1.43
case('PT')
    crosssection = 5.12
case('PD')
    crosssection = 2.88
case('AU')
    crosssection = 5.14
end select
scalefactor = crosssection_C / crosssection
pad = len(int2str(params%nptcls))
! graphene slice
call graphene_vol%new([params%box,params%box,params%box],params%smpd)
call graphene%new([params%boxpd,params%boxpd,1],params%smpd)
call orientation%new
call orientation%set_euler([0.,0.,0.])
call graphene_vol%read(graphene_fname)
call graphene_vol%neg
call vol_pad%new([params%boxpd, params%boxpd, params%boxpd], params%smpd)
call graphene_vol%pad(vol_pad)
if( params%griddev.eq.'yes' ) call vol_pad%div_w_instrfun(params%alpha)
call graphene_vol%kill
call del_file(graphene_fname)
call vol_pad%fft
call vol_pad%expand_cmat(params%alpha)
call vol_pad%fproject(orientation, graphene)
! prep particle
tfun = ctf(params%smpd, 300., 0.0, 0.4)
call vol_pad%new([params%boxpd, params%boxpd, params%boxpd], params%smpd)
call particle_vol%new([params%box,params%box,params%box],params%smpd)
call particle%new([params%boxpd,params%boxpd,1],params%smpd)
print *,params%box,params%boxpd
call img%new([params%box,params%box,1],params%smpd)
call particle_vol%read(particle_fname)
call particle_vol%neg
call particle_vol%pad(vol_pad)
if( params%griddev.eq.'yes' ) call vol_pad%div_w_instrfun(params%alpha)
call vol_pad%mul(scalefactor)
call particle_vol%kill
! call del_file(particle_fname)
call vol_pad%fft
call vol_pad%expand_cmat(params%alpha)
do i=1,params%nptcls
    call progress(i,params%nptcls)
    ! zero images
    particle = cmplx(0.,0.)
    img = 0.
    ! extract ori
    call spiral%set(i,'dfx',0.03+(ran3()-0.5)/5000.)
    call spiral%get_ori(i, orientation)
    ! project volume
    call vol_pad%fproject(orientation, particle)
    ! shift
    x = 20.*cos(real(i-1)/real(params%nptcls)*PI)        +(ran3()-0.5)/2.
    y = 15.*cos(real(i-1)/real(params%nptcls)*PI+PI*0.25)+(ran3()-0.5)/2.
    call spiral%set(i,'x',x)
    call spiral%set(i,'y',y)
    call particle%shift([x,y,0.])
    ! add
    call particle%add(graphene)
    ! simulate
    call simimg(particle, orientation, tfun, params%ctf, params%snr, snr_pink, snr_detector, params%bfac)
    call particle%ifft
    ! clip & write
    call particle%clip(img)
    call img%write('frame_'//int2str_pad(i,pad)//'.mrc')
end do
call spiral%write('trajectory.txt')
end program simple_test_simnano
