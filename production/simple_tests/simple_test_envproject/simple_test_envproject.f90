program simple_test_envproject
use simple_image,  only: image
use simple_masker, only: masker
use simple_oris,   only: oris
use simple_ori,    only: ori
use simple_defs
implicit none
character(len=STDLEN) :: fname = 'automask.mrc'
type(masker)          :: mskobj
integer, parameter    :: BOX=256, NIMGS=200
real,    parameter    :: SMPD=1.2156, MSK=72.
type(image)           :: vol_amask_soft, img
type(oris)            :: os
type(ori)             :: o
integer               :: i
call vol_amask_soft%new([BOX,BOX,BOX], smpd)
call vol_amask_soft%read(trim(fname))
mskobj = vol_amask_soft
call mskobj%init_envmask2D(MSK)
call img%new([BOX,BOX,1], smpd)
call os%new(NIMGS)
call os%spiral
do i=1,NIMGS
	img = 1.0
	o = os%get_ori(i)
	call mskobj%apply_envmask2D(o, img)
	call img%write('envprojs.mrc', i)
end do
end program simple_test_envproject
