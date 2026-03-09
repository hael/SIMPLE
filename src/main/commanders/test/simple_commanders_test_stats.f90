!@descr: for all stats tests
module simple_commanders_test_stats
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_class_sample_test
  contains
    procedure :: execute      => exec_test_class_sample_test
end type commander_test_class_sample_test

type, extends(commander_base) :: commander_test_clustering
  contains
    procedure :: execute      => exec_test_clustering
end type commander_test_clustering

type, extends(commander_base) :: commander_test_ctf_test
  contains
    procedure :: execute      => exec_test_ctf_test
end type commander_test_ctf_test

type, extends(commander_base) :: commander_test_eo_diff
  contains
    procedure :: execute      => exec_test_eo_diff
end type commander_test_eo_diff

type, extends(commander_base) :: commander_test_extr_frac
  contains
    procedure :: execute      => exec_test_extr_frac
end type commander_test_extr_frac

type, extends(commander_base) :: commander_test_multinomal_test
  contains
    procedure :: execute      => exec_test_multinomal_test
end type commander_test_multinomal_test

type, extends(commander_base) :: commander_test_pca_all
  contains
    procedure :: execute      => exec_test_pca_all
end type commander_test_pca_all

type, extends(commander_base) :: commander_test_pca_imgvar
  contains
    procedure :: execute      => exec_test_pca_imgvar
end type commander_test_pca_imgvar

type, extends(commander_base) :: commander_test_sp_project
  contains
    procedure :: execute      => exec_test_sp_project
end type commander_test_sp_project

contains

subroutine exec_test_class_sample_test( self, cline )
    use simple_class_sample_io
    class(commander_test_class_sample_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    type(class_sample), allocatable :: cs(:), cs2(:)
    call make_cs(cs)
    call make_cs(cs2)
    call print_class_sample(cs(1))
    call print_class_sample(cs(2))
    call print_class_sample(cs(3))
    print *, '*********************'
    call write_class_samples(cs, string('clssmp.bin'))
    call read_class_samples(cs2, string('clssmp.bin'))
    print *, class_samples_same(cs(1), cs2(1))
    print *, class_samples_same(cs(2), cs2(2))
    print *, class_samples_same(cs(3), cs2(3))
    call simple_end('**** SIMPLE_TEST_CLASS_SAMPLE_TEST_WORKFLOW NORMAL STOP ****')
    contains

        subroutine make_cs( cs )
            type(class_sample), allocatable, intent(inout) :: cs(:)
            allocate(cs(3))
            cs(1)%clsind = 1
            cs(2)%clsind = 2
            cs(3)%clsind = 3
            cs(1)%pop    = 1
            cs(2)%pop    = 2
            cs(3)%pop    = 3
            allocate(cs(1)%pinds(1), source=[1])
            allocate(cs(2)%pinds(2), source=[1,2])
            allocate(cs(3)%pinds(3), source=[1,2,3])
            allocate(cs(1)%ccs(1),   source=[1.])
            allocate(cs(2)%ccs(2),   source=[1.,0.])
            allocate(cs(3)%ccs(3),   source=[1.,0.,-1.])
        end subroutine make_cs


end subroutine exec_test_class_sample_test

subroutine exec_test_clustering( self, cline )
    use simple_aff_prop
    use simple_linalg
    class(commander_test_clustering),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call test_aff_prop
    call simple_end('**** SIMPLE_TEST_CLUSTERING_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_clustering

subroutine exec_test_ctf_test( self, cline )
    use simple_image,           only: image
    use simple_ctf,             only: ctf
    use simple_memoize_ft_maps, only: memoize_ft_maps
    class(commander_test_ctf_test),    intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    integer, parameter :: LDIM(3) = [256,256,1]
    real,    parameter :: SMPD = 1.0, DFX = 2.0, DFY = 2.0, ANGAST = 0., KV = 300., CS = 2.0, AC = 0.1
    type(image)        :: img, img_spec
    type(ctf)          :: tfun
    call img%new(LDIM, SMPD)
    call img_spec%new(LDIM, SMPD)
    tfun = ctf(SMPD, KV, CS, AC)
    call memoize_ft_maps(LDIM, SMPD)
    call img%ctf2img(tfun, DFX, DFY, ANGAST)
    call img%ft2img('real', img_spec)
    call img_spec%write(string('ctfimg.mrc'))
    call simple_end('**** SIMPLE_TEST_CTF_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_ctf_test

subroutine exec_test_eo_diff( self, cline )
    use simple_image, only: image
    class(commander_test_eo_diff),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    character(len=*), parameter :: name_vol      = 'recvol_state01.mrc'
    character(len=*), parameter :: name_vol_even = 'recvol_state01_even.mrc'
    character(len=*), parameter :: name_vol_odd  = 'recvol_state01_odd.mrc'
    integer,          parameter :: box           = 300
    integer,          parameter :: ldim(3)       = [box,box,box]
    real,             parameter :: smpd          = 1.2156
    type(image)       :: vol, vol_even, vol_odd, vol_noise
    integer           :: filtsz
    real, allocatable :: res(:), corrs(:)
    call vol%new(ldim, smpd)
    call vol_even%new(ldim, smpd)
    call vol_odd%new(ldim, smpd)
    call vol%read(string(name_vol))
    call vol_even%read(string(name_vol_even))
    call vol_odd%read(string(name_vol_odd))
    call vol_noise%copy(vol_even)
    call vol_noise%subtr(vol_odd)
    call vol_noise%write(string('noisevol_state01.mrc'))
    res = vol_even%get_res()
    call vol_even%fft
    call vol_odd%fft
    call vol_noise%fft
    filtsz = vol_even%get_filtsz()
    allocate(corrs(filtsz))
    call vol_even%ran_phases_below_noise_power(vol_odd)
    call vol_even%ifft
    call vol_even%write(string('new_impl.mrc'))
    call simple_end('**** SIMPLE_TEST_EO_DIFF_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_eo_diff

subroutine exec_test_extr_frac( self, cline )
    use simple_decay_funs
    class(commander_test_extr_frac),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    integer, parameter :: MAXITS = 100 ! upper iteration bound
    integer            :: i, nptcls
    ! nptcls = 5000
    ! print *, calc_nsampl_fromto( nptcls )
    ! nptcls = 10000
    ! print *, calc_nsampl_fromto( nptcls )
    nptcls = 36000
    ! print *, calc_nsampl_fromto( nptcls )
    ! nptcls = 80000
    ! print *, calc_nsampl_fromto( nptcls )
    ! nptcls = 200000
    ! print *, calc_nsampl_fromto( nptcls )
    ! nptcls = 1000000
    ! print *, calc_nsampl_fromto( nptcls )
    ! do i = 1,MAXITS
    !     nsampl      = nsampl_decay( i, MAXITS, NPTCLS)
    !     update_frac = real(nsampl) / real(NPTCLS)
    !     print *, i, nsampl, update_frac
    ! end do
    do i = 1,MAXITS
        print *, i, cos_decay( i, maxits, [0.05,1.0] )
    end do
    ! do i = 1,MAXITS
    !     nsampl      = inv_nsampl_decay( i, MAXITS, NPTCLS, NSAMPLE_MINMAX_DEFAULT)
    !     update_frac = real(nsampl) / real(NPTCLS)
    !     print *, i, nsampl, update_frac
    ! end do
    call simple_end('**** SIMPLE_TEST_EXTR_FRAC_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_extr_frac

subroutine exec_test_multinomal_test( self, cline )
    use simple_rnd
    class(commander_test_multinomal_test), intent(inout) :: self
    class(cmdline),                        intent(inout) :: cline
    real    :: pvec(7), cnts(7)
    integer :: which, i
    pvec(7) = 10.
    pvec(6) = 20.
    pvec(5) = 30.
    pvec(4) = 100.
    pvec(3) = 80.
    pvec(2) = 20.
    pvec(1) = 10.
    pvec    = pvec / sum(pvec)
    print *, pvec
    ! sample the distribution and calculate frequencies
    call seed_rnd
    cnts = 0.
    do i=1,10000
        which = multinomal( pvec )
        cnts(which) = cnts(which) + 1.0
    end do
    pvec = cnts / sum(cnts)
    print *, pvec
    call simple_end('**** SIMPLE_TEST_MULTINOMAL_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_multinomal_test

subroutine exec_test_pca_all( self, cline )
    use simple_ppca_inmem, only: ppca_inmem
    use simple_pca_svd,    only: pca_svd
    use simple_kpca_svd,   only: kpca_svd
    use simple_cmdline,    only: cmdline
    use simple_parameters, only: parameters
    class(commander_test_pca_all),    intent(inout) :: self
    class(cmdline),                   intent(inout) :: cline
    integer, parameter :: NP = 3, NS = 4, NC = 3, MAXPCAITS = 15
    type(ppca_inmem)   :: prob_pca
    type(pca_svd)      :: pca_obj
    type(kpca_svd)     :: kpca_obj
    type(parameters)   :: params
    integer :: j
    real    :: data_ori(NP, NS), avg(NP), tmpvec(NP), data_pca(NP, NS), E_zn(NC, NS), data_cen(NP, NS)
    call params%new(cline)
    data_ori(1,:) = [ 1, 2, 3, 4]
    data_ori(2,:) = [ 3, 1, 5, 8]
    data_ori(3,:) = [-1, 0, 4, 10]
    ! data_ori(4,:) = [ 0, 0, 7, 10]
    ! data_ori(5,:) = [-2, 0, 1, 10]
    print *, 'Original data:'
    print *, data_ori(1,:)
    print *, data_ori(2,:)
    print *, data_ori(3,:)
    avg = sum(data_ori, dim=2) / real(NS)
    do j = 1, NP
         data_cen(j,:) = data_ori(j,:) - avg(j)
    enddo
    call prob_pca%new(NS, NP, NC)
    call prob_pca%master(data_cen, MAXPCAITS)
    print *, 'PPCA eigenvalues:'
    !$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
    do j = 1, NS
        call prob_pca%generate(j, avg, tmpvec)
        data_pca(:,j) = tmpvec
    end do
    !$omp end parallel do
    print *, 'Pre-imaged data using PPCA:'
    do j = 1, NP
        print *, data_pca(j,:)
    enddo
    print *, 'Feature vecs using PPCA:'
    do j = 1, NS
        E_zn(:,j) = prob_pca%get_feat(j)
    enddo
    do j = 1, NC
        print *, E_zn(j, :)
    enddo
    print *, '---------------------------------------------------'
    ! PCA test
    print *, 'PCA eigenvalues/eigenvectors:'
    call pca_obj%new(NS, NP, NC)
    call pca_obj%master(data_cen)
    print *, 'Feature vecs using PCA:'
    do j = 1, NS
        E_zn(:,j) = pca_obj%get_feat(j)
    enddo
    do j = 1, NC
        print *, E_zn(j, :)
    enddo
    !$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
    do j = 1, NS
        call pca_obj%generate(j, avg, tmpvec)
        data_pca(:,j) = tmpvec
    end do
    !$omp end parallel do
    print *, 'Pre-imaged data using PCA:'
    do j = 1, NP
        print *, data_pca(j,:)
    enddo
    print *, '---------------------------------------------------'
    ! kPCA test
    call kpca_obj%new(NS, NP, NC)
    call kpca_obj%set_params(params%nthr, params%kpca_ker, params%kpca_target)
    call kpca_obj%master(data_cen)
    !$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
    do j = 1, NS
        call kpca_obj%generate(j, avg, tmpvec)
        data_pca(:,j) = tmpvec
    end do
    !$omp end parallel do
    print *, 'Pre-imaged data using kPCA:'
    do j = 1, NP
        print *, data_pca(j,:)
    enddo
    call simple_end('**** SIMPLE_TEST_PCA_ALL_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_pca_all

subroutine exec_test_pca_imgvar( self, cline )
    use simple_ppca_inmem, only: ppca_inmem
    use simple_pca_svd,    only: pca_svd
    use simple_kpca_svd,   only: kpca_svd
    use simple_cmdline,    only: cmdline
    use simple_parameters, only: parameters
    class(commander_test_pca_imgvar),   intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    integer, parameter :: NX = 5, NY = 5, NP = NX*NY, NC = 2, MAXPCAITS = 15
    type(ppca_inmem)   :: prob_pca
    integer :: i, j, cnt
    real    :: imgs(NX, NY, NC), flat_img(NP), dist_x(NP), dist_y(NP)
    real    :: data_ori(NP, NP), avg(NP), E_zn(NC, NP), data_cen(NP, NP), tmpvec(NP)
    real    :: lists(4), result

    lists  = [2., 1., 3., 4.]
    result = 0.
    do i = 1, size(lists)
        do j = 1, size(lists)
            result = result + lists(i) * lists(j)
        enddo
    enddo
    print *, result

    imgs(1,:,1) = [ 0., 2., 3., 4., 0.]
    imgs(2,:,1) = [ 0., 1., 5., 0., 0.]
    imgs(3,:,1) = [ 0., 0., 0., 0., 0.]
    imgs(4,:,1) = [ 0., 0., 0., 0., 0.]
    imgs(5,:,1) = [ 0., 0., 0., 0., 0.]

    imgs(1,:,2) = [ 0.,  0.,  0.,  0., 0.]
    imgs(2,:,2) = [ 0.,  0.,  0.,  0., 0.]
    imgs(3,:,2) = [ 0.,  0.,  0.,  0., 0.]
    imgs(4,:,2) = [ 0.,  0., -5., -9., 0.]
    imgs(5,:,2) = [ 0.,-11., -1., -2., 0.]

    cnt = 1
    do i = 1, NY
        do j = 1, NX
            flat_img(cnt) = imgs(i,j,1) + imgs(i,j,2)
            dist_x(cnt)   = real(i)
            dist_y(cnt)   = real(j)
            cnt           = cnt + 1
        enddo
    enddo

    data_ori = 0.
    do i = 1, NP
        do j = 1, NP
            if( i == j ) cycle
            data_ori(i,j) = abs(flat_img(i) - flat_img(j)) / abs(real(i - j))
        enddo
    enddo
    avg = sum(data_ori, dim=2) / real(NP)
    do j = 1, NP
        data_cen(j,:) = data_ori(j,:) - avg(j)
    enddo
    call prob_pca%new(NP, NP, NC)
    call prob_pca%master(data_cen, MAXPCAITS)
    !$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
    do j = 1, NP
        call prob_pca%generate(j, avg, tmpvec)
    end do
    !$omp end parallel do
    print *, 'Feature vecs using PPCA:'
    do j = 1, NP
        E_zn(:,j) = prob_pca%get_feat(j)
    enddo
    do j = 1, 1
        print *, E_zn(j, 1:5)
        print *, E_zn(j, 6:10)
        print *, E_zn(j,11:15)
        print *, E_zn(j,16:20)
        print *, E_zn(j,21:25)
    enddo
    do j = 1, 1
        print *, imgs(1,:,1)
        print *, imgs(2,:,1)
        print *, imgs(3,:,1)
        print *, imgs(4,:,1)
        print *, imgs(5,:,1)
    enddo
    call simple_end('**** SIMPLE_TEST_PCA_IMGVAR_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_pca_imgvar

subroutine exec_test_sp_project( self, cline )
    use simple_sp_project, only: sp_project
    use simple_binoris_io, only: binwrite_oritab
    class(commander_test_sp_project),   intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    integer, parameter :: NMICS = 87
    integer, parameter :: NPTCLS = 5646
    type(sp_project)   :: project1, project2, project3
    type(string)       :: template, str1, str2
    real    :: eullims(3,2)
    integer :: i
    call seed_rnd
    eullims(1,:) = [0.,360.]
    eullims(2,:) = [0.,180.]
    eullims(3,:) = [0.,360.]
    ! prepare dummy project
    call project1%os_mic%new(NMICS, is_ptcl=.false.)
    do i = 1,NMICS
        template = 'FoilHole_'//int2str_pad(i,16)
        call project1%os_mic%set(i, 'movie',    template%to_char()//'_fractions.tiff')
        call project1%os_mic%set(i, 'intg',     template%to_char()//'_intg.mrc')
        call project1%os_mic%set(i, 'forctf',   template%to_char()//'_forctf.mrc')
        call project1%os_mic%set(i, 'thumb',    template%to_char()//'_thumb.jpg')
        call project1%os_mic%set(i, 'smpd',     1.34)
        call project1%os_mic%set(i, 'kv',       300.)
        call project1%os_mic%set(i, 'cs',       2.7 )
        call project1%os_mic%set(i, 'fraca',    0.1 )
        call project1%os_mic%set(i, 'dfx',      1.+ran3())
        call project1%os_mic%set(i, 'dfy',      1.+ran3())
        call project1%os_mic%set(i, 'angast',   360.*ran3())
        call project1%os_mic%set(i, 'phshift',  0.)
    enddo
    ! prepare dummy project
    call project1%os_ptcl2D%new(NPTCLS, .true.)
    do i = 1,NPTCLS
        call project1%os_ptcl2D%set_dfx(i, 1.+ran3())
        call project1%os_ptcl2D%set_dfy(i, 1.+ran3())
        call project1%os_ptcl2D%set(i, 'corr',    ran3())
        call project1%os_ptcl2D%set_class(i, nint(ran3()*100.))
        call project1%os_ptcl2D%set_state(i, nint(ran3()))
    enddo
    call project1%os_ptcl2D%set_all2single('w', 1.0)
    project1%os_ptcl3D = project1%os_ptcl2D
    call project1%os_ptcl2D%rnd_oris( trs=5.0, eullims=eullims)
    call project1%os_ptcl2D%delete_entry('z')
    call project1%os_ptcl3D%rnd_oris( trs=5.0, eullims=eullims)
    call project1%os_ptcl3D%delete_entry('class')
    call project1%os_ptcl3D%delete_entry('z')
    do i = 1,NPTCLS
        call project1%os_ptcl3D%set(i, 'proj', nint(1000.*ran3()))
    enddo
    call project1%update_projinfo(string('myproject.simple'))
    ! write/read
    call project1%write(string('myproject.simple'))
    call project2%read(string('myproject.simple'))
    ! compare
    do i = 1,NMICS
        str1 = project1%os_mic%ori2str(i)
        str2 = project2%os_mic%ori2str(i)
        if( project1%os_mic%ori2str(i) /= project2%os_mic%ori2str(i) )then
            write(*,*)'1 TEST FAILED COMPARING ', str1%to_char(),' AND ', str2%to_char()
            stop
        endif
        call str1%kill
        call str2%kill
    enddo
    write(*,*)'TEST SUCCES WRITE/READ/COMPARE os_mic'
    do i = 1,NPTCLS
        str1 = project1%os_ptcl2D%ori2str(i)
        str2 = project2%os_ptcl2D%ori2str(i)
        if( project1%os_ptcl2D%ori2str(i) /= project2%os_ptcl2D%ori2str(i) )then
            write(*,*)'2 TEST FAILED COMPARING ', str1%to_char(),' AND ', str2%to_char()
            stop
        endif
        call str1%kill
        call str2%kill
    enddo
    write(*,*)'TEST SUCCES WRITE/READ/COMPARE os_ptcl2D'
    do i = 1,NPTCLS
        str1 = project1%os_ptcl3D%ori2str(i)
        str2 = project2%os_ptcl3D%ori2str(i)
        if( project1%os_ptcl3D%ori2str(i) /= project2%os_ptcl3D%ori2str(i) )then
            write(*,*)'3 TEST FAILED COMPARING ', str1%to_char(),' AND ', str2%to_char()
            stop
        endif
        call str1%kill
        call str2%kill
    enddo
    write(*,*)'TEST SUCCES WRITE/READ/COMPARE os_ptcl3D'
    call project2%kill
    ! merging
    call project3%os_ptcl2D%new(NPTCLS, .true.)
    call project3%update_projinfo(string('project3.simple'))
    call binwrite_oritab(string('doc_1.simple'), project1, project1%os_ptcl2D, [  1,   1882],  isegment=PTCL2D_SEG)
    call binwrite_oritab(string('doc_2.simple'), project1, project1%os_ptcl2D, [1883,  3764],  isegment=PTCL2D_SEG)
    call binwrite_oritab(string('doc_3.simple'), project1, project1%os_ptcl2D, [3765, NPTCLS], isegment=PTCL2D_SEG)
    call project3%merge_algndocs(NPTCLS, 3, 'ptcl2D', 'doc_', 1 )
    call del_file('doc_1.simple')
    call del_file('doc_2.simple')
    call del_file('doc_3.simple')
    do i = 1,NPTCLS
        str1 = project1%os_ptcl2D%ori2str(i)
        str2 = project3%os_ptcl2D%ori2str(i)
        if( project1%os_ptcl2D%ori2str(i) /= project3%os_ptcl2D%ori2str(i) )then
            write(*,*)'TEST FAILED COMPARING ',str1%to_char(),' AND ', str2%to_char()
            stop
        endif
    enddo
    write(*,*)'TEST SUCCES WRITE PARTS/MERGE/COMPARE os_ptcl2D'
    call project3%kill
    ! some i/o functionalities
    call project1%print_segment_json( 'ptcl3D', string('myproject.simple'), [NPTCLS-1,NPTCLS])
    call project1%write_segment2txt('mic', string('some_mics.txt'), [2,5])
    write(*,*)'TEST SUCCES WRITE/PRINT'
    call del_file('myproject.simple')
    call del_file('project3.simple')
    call del_file('some_mics.txt')
    call simple_end('**** SIMPLE_TEST_SP_PROJECT_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_sp_project

end module simple_commanders_test_stats
