program simple_test_fftw_plan_many
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_commander_volops,   only: reproject_commander
use simple_image,              only: image
use simple_fftw3
implicit none
type(parameters)              :: p
type(cmdline)                 :: cline, cline_projection
type(reproject_commander)     :: xreproject
integer                       :: k, l, nptcls, rc, iptcl, cnt, c_shape(3)
type(image)     , allocatable :: imgs(:)
character(len=:), allocatable :: cmd
complex,          pointer     :: cmat_img(:,:,:)
real,             pointer     :: rmat_img(:,:,:)
logical                       :: mrc_exists
integer(timer_int_kind)       ::  t_tot
real(timer_int_kind)          :: rt_tot
type(c_ptr)                   :: plan_fwd, plan_bwd, ptr             !< fftw plan for the image (fwd)
real(   kind=c_float),         pointer ::  in(:,:,:)
complex(kind=c_float_complex), pointer :: out(:,:,:)
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'Usage: simple_test_fftw_plan_many smpd=xx nthr=yy stk=stk.mrc, mskdiam=zz'
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
        call cline_projection%set('nspace'    , 2.)
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
nptcls = 1
p%ldim(3) = 1 ! because we operate on stacks
allocate(imgs(nptcls))
do iptcl = 1, nptcls
    call imgs(iptcl)%new(p%ldim, p%smpd)
    call imgs(iptcl)%read(p%stk, iptcl)
enddo
t_tot = tic()
do iptcl = 1, nptcls
    call imgs(iptcl)%fft()
    call imgs(iptcl)%get_cmat_ptr(cmat_img)
    if( iptcl == 1 ) print *, cmat_img(1,162:166,1)
    cmat_img = cmat_img/sum(cmat_img)
    call imgs(iptcl)%ifft()
enddo
rt_tot = toc(t_tot)
print *, 'total time = ', rt_tot
do iptcl = 1, nptcls
    call imgs(iptcl)%write('test1.mrc', iptcl)
enddo
! with fft_plan_many
c_shape = [p%ldim(2), p%ldim(1)/2+1, nptcls]
ptr = fftwf_alloc_complex(int(product(c_shape),c_size_t))
! Set up the complex array which will point at the allocated memory
call c_f_pointer(ptr,out,c_shape)
! Work out the shape of the real array
c_shape(2) = c_shape(2)*2
! Set up the real array
call c_f_pointer(ptr,in,c_shape)
c_shape(2) = c_shape(2)/2
!$omp critical
call fftwf_plan_with_nthreads(nthr_glob)
plan_fwd = fftwf_plan_many_dft_r2c(2,  [p%ldim(2), p%ldim(1)],nptcls,&
                                  &in ,[p%ldim(2), p%ldim(1)],1,product([p%ldim(1), p%ldim(2)]),&
                                  &out,[p%ldim(2), p%ldim(1)],1,product([p%ldim(1), p%ldim(2)]),FFTW_ESTIMATE)
plan_bwd = fftwf_plan_many_dft_c2r(2,  [p%ldim(2), p%ldim(1)],nptcls,&
                                  &out,[p%ldim(2), p%ldim(1)],1,product([p%ldim(1), p%ldim(2)]),&
                                  &in ,[p%ldim(2), p%ldim(1)],1,product([p%ldim(1), p%ldim(2)]),FFTW_ESTIMATE)
!$omp end critical
do iptcl = 1, nptcls
    call imgs(iptcl)%new(p%ldim, p%smpd)
    call imgs(iptcl)%read(p%stk, iptcl)
enddo
cnt = 1
t_tot = tic()
do iptcl = 1, nptcls
    call imgs(iptcl)%get_rmat_ptr(rmat_img)
    do k = 1, p%ldim(1)+2
        do l = 1, p%ldim(2)
            in(l,k,iptcl) = rmat_img(k,l,1)
            cnt = cnt + 1
        enddo
    enddo
enddo
call fftwf_execute_dft_r2c(plan_fwd,in,out)
rt_tot = toc(t_tot)
cnt = 1
do iptcl = 1, nptcls
    call imgs(iptcl)%fft()
    call imgs(iptcl)%get_cmat_ptr(cmat_img)
    do k = 1, p%ldim(1)/2+1
        do l = 1, p%ldim(2)
            cmat_img(k,l,1) = out(l,k,iptcl)/product([p%ldim(1), p%ldim(2)])
            if( mod(k+l,2) == 1 ) cmat_img(k,l,1) = -cmat_img(k,l,1)
            cnt = cnt + 1
        enddo
    enddo
    if( iptcl == 1 ) print *, cmat_img(1,162:166,1)
    call imgs(iptcl)%ifft()
enddo
!print *, out(1:10)
do iptcl = 1, nptcls
    call imgs(iptcl)%write('test2.mrc', iptcl)
enddo
print *, 'total time = ', rt_tot
deallocate(imgs)
end program simple_test_fftw_plan_many