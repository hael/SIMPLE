program simple_test_polarops
use simple_pftc_srch_api
use simple_builder, only: builder
implicit none
complex,   allocatable :: pft(:,:)
integer,     parameter :: N=128
integer,     parameter :: NIMGS=200
integer,     parameter :: NCLS=5
type(image)            :: tmpl_img, img, cavgs(NCLS)
type(cmdline)          :: cline
type(polarft_calc)     :: pftc
type(parameters)       :: p
type(builder)          :: b
real    :: ang, shift(2), shifts(2,NIMGS)
integer :: pinds(NIMGS), i, eo, icls
! dummy structure
call tmpl_img%soft_ring([N,N,1], 1., 8.)
call tmpl_img%fft
call tmpl_img%shift2Dserial([ 8.,-16.])
call img%soft_ring([N,N,1], 1., 12.)
call img%fft
call img%shift2Dserial([ 32., 0.])
call tmpl_img%add(img)
call img%soft_ring([N,N,1], 1., 16.)
call img%fft
call img%shift2Dserial([ -16., 8.])
call tmpl_img%add(img)
call img%soft_ring([N,N,1], 1., 32.)
call img%fft
call tmpl_img%add(img)
call tmpl_img%ifft
call tmpl_img%write(string('template.mrc'))
! init of options & parameters
call cline%set('prg',    'xxx')
call cline%set('objfun', 'cc')
call cline%set('smpd',   1.0)
call cline%set('box',    N)
call cline%set('ctf',    'no')
call cline%set('oritype','ptcl2D')
call cline%set('ncls',    NCLS)
call cline%set('nptcls',  NIMGs)
call cline%set('lp',      3.)
call cline%set('nthr',    8)
call cline%set('mskdiam', real(N)/2-10.)
call cline%set('ref_type', 'polar_cavg')
! Calculators
call b%init_params_and_build_strategy2D_tbox(cline, p)
call pftc%new(NCLS, [1,NIMGS], p%kfromto)
pinds = (/(i,i=1,NIMGS)/)
call b%img_crop%memoize4polarize(pftc%get_pdim())
pft = pftc%allocate_pft()
do i = 1,NIMGS
    shift = 10.*[ran3(), ran3()] - 5.
    ! ang   = 360. * ran3()
    ang   = 0.
    eo    = 0
    if( .not.is_even(i) ) eo = 1
    icls  = ceiling(ran3()*4.)
    call img%copy_fast(tmpl_img)
    call img%fft
    call img%shift2Dserial(-shift)
    call img%ifft
    call img%rtsq(ang, 0.,0.)
    call img%add_gauran(2.)
    call img%write(string('rotimgs.mrc'), i)
    call img%fft
    call b%spproj_field%set_euler(i, [0.,0.,ang])
    call b%spproj_field%set_shift(i, shift)
    call b%spproj_field%set(i,'w',1.0)
    call b%spproj_field%set(i,'state',1)
    call b%spproj_field%set(i,'class', icls)
    call b%spproj_field%set(i,'eo',eo)
    shifts(:,i) = -shift
    call img%polarize(pft, mask=b%l_resmsk)
    call pftc%set_ptcl_pft(i, pft)
enddo
call pftc%polar_cavger_new(.false.)
call pftc%polar_cavger_update_sums(NIMGS, pinds, b%spproj, shifts)
call pftc%polar_cavger_merge_eos_and_norm2D
call pftc%polar_cavger_calc_and_write_frcs_and_eoavg(b%clsfrcs, b%spproj_field%get_update_frac(), string(FRCS_FILE), cline)
! write
call pftc%polar_cavger_write(string('cavgs_even.bin'), 'even')
call pftc%polar_cavger_write(string('cavgs_odd.bin'),  'odd')
call pftc%polar_cavger_write(string('cavgs.bin'),      'merged')
call pftc%polar_cavger_refs2cartesian(cavgs, 'even')
call write_imgarr(cavgs, string('cavgs_even.mrc'))
call pftc%polar_cavger_refs2cartesian(cavgs, 'odd')
call write_imgarr(cavgs, string('cavgs_odd.mrc'))
call pftc%polar_cavger_refs2cartesian(cavgs, 'merged')
call write_imgarr(cavgs, string('cavgs_merged.mrc'))
call pftc%polar_cavger_kill
! read & write again
call pftc%polar_cavger_new(.false.)
call pftc%polar_cavger_read(string('cavgs_even.bin'), 'even')
call pftc%polar_cavger_read(string('cavgs_odd.bin'),  'odd')
call pftc%polar_cavger_read(string('cavgs.bin'),      'merged')
call pftc%polar_cavger_refs2cartesian(cavgs, 'even')
call write_imgarr(cavgs, string('cavgs2_even.mrc'))
call pftc%polar_cavger_refs2cartesian(cavgs, 'odd')
call write_imgarr(cavgs, string('cavgs2_odd.mrc'))
call pftc%polar_cavger_refs2cartesian(cavgs, 'merged')
call write_imgarr(cavgs, string('cavgs2_merged.mrc'))
call pftc%polar_cavger_kill
if( allocated(pft) ) deallocate(pft)
end program simple_test_polarops
