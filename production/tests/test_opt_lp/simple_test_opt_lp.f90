program simple_test_opt_lp
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_commander_volops,   only: reproject_commander
use simple_image,              only: image
use simple_opt_filter,         only: butterworth_filter
implicit none
type(parameters)              :: p
type(cmdline)                 :: cline, cline_projection
type(image)                   :: img, noise, res_map, img_noisy, img_filt
type(reproject_commander)     :: xreproject
integer                       :: k, nyq, nptcls, rc, iptcl, find_stop, find_start, n_bin, n_vec, find_cur
real                          :: ave, sdev, maxv, minv, med, vec_mean, vec_std
character(len=:), allocatable :: cmd
logical                       :: mrc_exists
real,             allocatable :: cur_fil(:), vec_noise(:), xhist(:), yest(:)
integer,          allocatable :: yhist(:)
real,             pointer     :: rmat_img_noisy(:,:,:), rmat_img_filt(:,:,:)
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'Usage: simple_test_2D_opt_filt smpd=xx nthr=yy stk=stk.mrc, mskdiam=zz'
    write(logfhandle,'(a)') 'Example: projections of https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
    write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
    inquire(file="1JYX.mrc", exist=mrc_exists)
    if( .not. mrc_exists )then
        write(*, *) 'Downloading the example dataset...'
        cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
        write(*, *) 'Converting .pdb to .mrc...'
        cmd = 'e2pdb2mrc.py 1JYX.pdb 1JYX.mrc'
        call execute_command_line(cmd, exitstat=rc)
        cmd = 'rm 1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
        write(*, *) 'Projecting 1JYX.mrc...'
        call cline_projection%set('vol1'      , '1JYX.mrc')
        call cline_projection%set('smpd'      , 1.)
        call cline_projection%set('pgrp'      , 'c1')
        call cline_projection%set('mskdiam'   , 180.)
        call cline_projection%set('nspace'    , 6.)
        call cline_projection%set('nthr'      , 16.)
        call xreproject%execute(cline_projection)
    endif
    call cline%set('smpd'   , 1.)
    call cline%set('nthr'   , 16.)
    call cline%set('stk'    , 'reprojs.mrcs')
    call cline%set('mskdiam', 180.)
else
    call cline%parse_oldschool
endif
call cline%checkvar('smpd',    1)
call cline%checkvar('nthr',    2)
call cline%checkvar('stk' ,    3)
call cline%checkvar('mskdiam', 4)
call cline%check
call p%new(cline)
call find_ldim_nptcls(p%stk, p%ldim, nptcls)
p%ldim(3) = 1 ! because we operate on stacks
n_vec     = p%ldim(1)*p%ldim(2)
n_bin     = int(n_vec/1000.)
print *, n_vec, n_bin
call img%new(  p%ldim, p%smpd)
call noise%new(p%ldim, p%smpd)
allocate(cur_fil(p%ldim(1)), source=0.)
allocate(vec_noise(n_vec),   source=0.)
allocate(xhist(n_bin),       source=0.)
allocate(yest( n_bin),       source=0.)
allocate(yhist(n_bin),       source=0)
find_stop  = calc_fourier_index(p%lpstart,   p%ldim(1), p%smpd)
find_start = calc_fourier_index(p%lp_lowres, p%ldim(1), p%smpd)
do iptcl = 1, 1
    write(*, *) 'Particle # ', iptcl
    cur_fil = 0.
    call img%read(p%stk, iptcl)
    call img_noisy%copy(img)
    ! spherical masking
    call img%mask(p%msk, 'soft')
    ! img stats
    call img%stats('foreground', ave, sdev, maxv, minv)
    ! add noise in a small center region of the even
    call noise%gauran(0., .2 * sdev)
    call noise%mask(1.5 * p%msk, 'soft')
    call img_noisy%add(noise)
    call img_noisy%write('stk_noisy.mrc', iptcl)
    call img_noisy%get_rmat_ptr(rmat_img_noisy)
    ! find_cur = int((find_start + find_stop)/2.)
    ! find_cur = int(find_start + 50)
    do find_cur = find_start + 10, find_stop - 10, 5
        call img_filt%copy(img_noisy)
        call img_noisy%fft
        call img_filt%fft
        call butterworth_filter(find_cur - 5, find_cur + 5, cur_fil)
        call img_filt%apply_filter(cur_fil)
        call img_filt%add(img_noisy)
        call img_filt%ifft
        call img_noisy%ifft
        call img_filt%write('stk_filt.mrc', iptcl)
        call img_filt%get_rmat_ptr( rmat_img_filt )
        vec_noise = [rmat_img_noisy] - [rmat_img_filt]
        call create_hist_vector( vec_noise, n_bin, xhist, yhist )
        yest = real(yhist)
        call avg_sdev(yest, vec_mean, vec_std)
        yest = exp(-(yest-vec_mean)**2/2./vec_std**2)/2./pi/vec_std
        yest = yest*sum(xhist*yhist)/sum(xhist*yest)
        print *, sum( (yest - real(yhist))**2 )
    enddo
    ! print *, yhist
    call img%zero_and_unflag_ft
    call img_noisy%zero_and_unflag_ft
    call noise%zero_and_unflag_ft
enddo
call img%kill()
call img_noisy%kill()
call noise%kill()
end program simple_test_opt_lp
