!==Program simple_cluster
!
! <cluster/begin> is a program for image clustering based on reference-free in-plane alignment \citep{Penczek:1992aa} 
! and probabilistic principal component analysis (PCA) for generation of feature vectors \citep{Bishop:2006}. 
! Agglomerative hierarchical clustering (HAC) is used for grouping of feature vectors. Refinement of the clustering 
! solution is done using center-based k-means clustering. \prgname{simple\_cluster} in-plane aligns the input image 
! stack. Bicubic interpolation is used for shifting and rotating the stack before extraction of the pixels within the 
! circular mask defined by mask radius \texttt{msk}. Next, the probabilistic PCA method generates feature vectors from 
! the vectors of extracted pixels. The minimum cluster population (\texttt{minp}) prevents clusters below population 
! \texttt{minp} to be represented by an output average. <comment/begin> The setup allows for quick testing of the number 
! of clusters. One pass produces the file \texttt{pdfile.bin} containing the matrix of all pair-wise feature vector 
! distances. Using the optional parameter \texttt{dopca=no} in a second round of execution from the same directory will 
! make the program read the previously generated distances and re-do the clustering using whatever settings inputted for 
! parameters \texttt{ncls} \& \texttt{minp}. The optional parameter \texttt{oritab} is used to provide in-plane parameters 
! for the clustering (provided by program \prgname{simple\_prime}). This option can be used for generating class averages 
! that are going to be subjected to heterogeneity analysis. The default setting uses 30 eigenvectors, if you are 
! \textit{not} inputting in-plane parameters via optional parameter \texttt{oritab}, and 60 eigenvectors if you do input 
! in-plane parameters. Note that the distance matrix is kept in RAM, so for large data sets you need LOTS of internal memory. 
! This quirk can be addressed by using a random sample of the data for initial clustering by HAC. This is done by setting 
! \texttt{nran} to some number < \texttt{nptcls}. In this setting, the HAC centres generated from the random sample are used 
! to extend the clustering to the entire data set with k-means. This overcomes the well-known initialisation problem of k-means, 
! and enables clustering of many hundreds of thousands of particle images. SIMPLE has been used to cluster 300,000 images with
! a box size of 100 using a random subset of 60,000 images on a machine with 96 GB RAM. <comment/end> <cluster/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2011
!
program simple_cluster
use simple_defs         ! singleton
use simple_cmdline      ! singleton
use simple_procimgfile  ! singleton
use simple_jiffys,      only: alloc_err, simple_end, add2fbody
use simple_math,        only: hpsort, euclid
use simple_build,       only: build
use simple_params,      only: params
use simple_image,       only: image
use simple_procimgfile, only: make_pattern_stack, norm_imgfile
implicit none
type(params)          :: p
type(build)           :: b
type(image)           :: img
character(len=STDLEN) :: stknams(4)
character(len=STDLEN) :: eulnams(4)
integer               :: pcarecsz, i, utst
real, allocatable     :: AVG(:)
if( command_argument_count() < 3 )then
    write(*,'(a)', advance='no') 'SIMPLE_CLUSTER stk=<stack.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)', advance='no') ' ncls=<nr of clusters{500})> [msk=<mask radius(in pixels){box/2}>]'
    write(*,'(a)', advance='no') ' [minp=<minimum ptcls in cluster{10}>] [nran=<size of random sample{nptcls}>]'
    write(*,'(a)') ' [oritab=<SIMPLE alignment doc>] [nthr=<nr of openMP threads{1}>]'
    write(*,'(a)') '**less commonly used**' 
    write(*,'(a)', advance='no') '[nvars=<nr of eigen vectors{30/60}>] [clsdoc=<SPIDER clustering doc>]'
    write(*,'(a)', advance='no') ' [kmeans=<yes|no{yes}>] [dopca=<yes|no{yes}>] [doalign=<yes|no{yes}>]'
    write(*,'(a)') ' [rnd=<yes|no{no}>] [utst=<unit test nr(1-4){0}>]'
    stop
endif
call parse_cmdline
if( defined_cmd_arg('utst') )then
    utst = nint(get_cmd_rarg('utst'))
    if( utst < 1 .or. utst > 5 )then
        stop 'unrecognized unit test number; simple_cluster'
    else
        call utst_set_cmdlin(utst)
    endif
endif
call cmdcheckvar('stk',  1)
call cmdcheckvar('smpd', 2)
call cmdcheckvar('ncls', 3)
! set those that lacks standardized values in params and were not set by the user
if( .not. defined_cmd_arg('nvars') )then
    if( defined_cmd_arg('oritab') .or. p%doalign .eq. 'no' ) call set_cmdline('nvars', 60.)
endif
call cmdcheck                   ! checks args and prints cmdline to cmdline.dat
p = params()                    ! constants & derived constants produced
call b%build_general_tbox(p)    ! general objects built
call b%build_cluster_tbox(p)    ! objects for clustering built
! create image
call img%new([p%box,p%box,1],p%smpd)
if( p%utst > 0 )then
    write(*,'(A,1X,I1)') '>>> RUNNING UNIT TEST NR:', p%utst
    if( p%utst == 1 ) write(*,'(A,1X,A)') '>>> UNIT TEST INFO:', 'HOMOGENEOUS POLYMERASE WITH 0 INPLANE ANGLE (15K PTCLS)'
    if( p%utst == 2 ) write(*,'(A,1X,A)') '>>> UNIT TEST INFO:', 'HOMOGENEOUS POLYMERASE WITH RANDOM INPLANE ANGLE (15K PTCLS)'
    if( p%utst == 3 ) write(*,'(A,1X,A)') '>>> UNIT TEST INFO:', 'HETEROGENEOUS POLYMERASE WITH 0 INPLANE ANGLE (21K PTCLS)'
    if( p%utst == 4 ) write(*,'(A,1X,A)') '>>> UNIT TEST INFO:', 'HETEROGENEOUS RIBOSOME WITH 0 INPLANE ANGLE (10K PTCLS)'
endif
if( p%rnd .eq. 'yes' )then
    call b%a%rnd_classes(p%ncls)
    goto 199
endif
write(*,'(A)') '**** CLUSTERING PROCEDURE BEGINS ****'
if( p%clsdoc == '' )then
    if( p%doalign == 'yes' )then
        ! DO ROTATIONAL ALIGNMENT
        if( p%oritab == '' )then
            call b%ralgn%align
            call b%ralgn%kill
            write(*,'(A)') '>>> WROTE ALIGNMENT REFS TO SPIDER STACK: roalgn_refs.ext'
            call b%a%write('cluster_algndoc.dat')
            write(*,'(A)') '>>> WROTE ALIGNMENT DOC: cluster_algndoc.dat'
        endif
        ! shift and rotate stack
        call shrot_imgfile(p%stk, 'inplalgnstk'//p%ext, b%a, p%smpd)
        p%stk = 'inplalgnstk'//p%ext       
        write(*,'(A)') '>>> WROTE ALIGNED IMAGES TO SPIDER STACK: inplalgnstk'//p%ext
        write(*,'(A)') '**** FINISHED ROTATIONAL ALIGNMENT ****'
    endif
    if( p%dopca == 'yes' )then
        ! DO PROBABILISTIC PCA
        ! create pca vec stack
        call make_pattern_stack(p%stk, p%pcastk, p%msk, p%ncomps, pcarecsz, AVG )
        ! do the choka choka  
        call b%pca%master(p%pcastk, pcarecsz, p%kmstk, 100)
        call b%build_read_features(p)
        write(*,'(A)') '**** FINISHED GENERATIVE ITERATIVE PCA ****'
        write(*,'(A)') '>>> WROTE MANIFOLD SAMPLES TO FILE: generative'//p%ext
        ! stop here if nr of clusters is <= 1
        if( p%ncls <= 1 ) stop
        ! build distance matrix
        call b%pd%build(distfun)
        write(*,'(A)') '**** FINISHED BUILDING DISTANCE PD TABLE ****'
    else
        ! stop here if nr of clusters is <= 1
        if( p%ncls <= 1 ) stop
        call b%pd%read(p%pdfile)
        call b%build_read_features(p)
    endif
    ! DO HIERARCHICAL CLUSTERING
    call b%hacls%cluster(b%pd, p%pdfile, p%ncls, p%minp)
    call b%a%write_clsdoc('hcl.txt')
    write(*,'(A)') '>>> WROTE FIRST CLUSTERING SOLUTION TO FILE: hcl.txt'
    if( p%kmeans == 'yes' )then
        call b%kmcls%refine(100)
        call b%a%write_clsdoc('kmeans.txt')
        write(*,'(A)') '>>> WROTE REFINED CLUSTERING SOLUTION TO FILE: kmeans.txt'
    end if
    write(*,'(A)') '**** FINISHED CLUSTERING ****'
else
    ! READ IN A CLUSTERING SOLUTION
    call b%a%read_clsdoc(p%clsdoc)
endif
! MAKE CLASS AVERAGES
call make_cavgs_imgfile(p%stk, 'cavgstk'//p%ext, b%a, p%minp)
write(*,'(A)') 'WROTE CLASS AVERAGES TO SPIDER STACK: cavgstk'//p%ext
199 if( p%utst /= 0 ) call utst_calc_stat ! calculates unit test statistics if needed
! END GRACEFULLY
call simple_end('**** SIMPLE_CLUSTER NORMAL STOP ****')

contains

    function distfun( ia, ib ) result( dist )
        integer, intent(in) :: ia, ib
        real :: dist
        dist = euclid(b%features(ia,:), b%features(ib,:))
    end function
    
    subroutine utst_set_cmdlin( nr )
        integer, intent(in) :: nr
        character(len=STDLEN) :: simpledir
        ! set the command line
        call get_environment_variable('SIMPLEPATH', simpledir)
        stknams(1) = trim(simpledir)//'/unit_test_data/STKtest1'//p%ext     ! homogeneous polymerase with 0 inplane angle      (15k)
        stknams(2) = trim(simpledir)//'/unit_test_data/STKtest2'//p%ext     ! homogeneous polymerase with random inplane angle (15k)
        stknams(3) = trim(simpledir)//'/unit_test_data/STKtest3'//p%ext     ! heterogeneous polymerase with 0 inplane angle    (21k)
        stknams(4) = trim(simpledir)//'/unit_test_data/STKtest4'//p%ext     ! heterogeneous ribosome with 0 inplane angle      (10k)
        eulnams(1) = trim(simpledir)//'/unit_test_data/eulers_test1'//p%ext ! homogeneous polymerase with 0 inplane angle      (15k)
        eulnams(2) = trim(simpledir)//'/unit_test_data/eulers_test2'//p%ext ! homogeneous polymerase with random inplane angle (15k)
        eulnams(3) = trim(simpledir)//'/unit_test_data/eulers_test3'//p%ext ! heterogeneous polymerase with 0 inplane angle    (21k)
        eulnams(4) = trim(simpledir)//'/unit_test_data/eulers_test4'//p%ext ! heterogeneous ribosome with 0 inplane angle      (10k)     
                      call set_cmdline('stk', stknams(nr))
        if( nr <= 3 ) call set_cmdline('box',        100.)
        if( nr == 4 ) call set_cmdline('box',        128.)
        if( nr <= 3 ) call set_cmdline('smpd',       2.33)
        if( nr == 4 ) call set_cmdline('smpd',       2.81)
        if( nr <= 3 ) call set_cmdline('ring2',       42.)
        if( nr == 4 ) call set_cmdline('ring2',       55.)
                      call set_cmdline('nvars',       30.)
        if( nr == 2 )then
                      call set_cmdline('doalign',   'yes')
        else
                      call set_cmdline('doalign',    'no')
        endif
    end subroutine
    
    subroutine utst_calc_stat
        ! STATE PARTITIONINGS:
        !data set 3: 1-7000, 7001-14000, 14001-21000
        !data set 4: 1-5000, 5001-10000
        integer              :: ncls_here, cnt, alloc_stat
        real                 :: avg, sdev, med, sdevmed, avg10best, sdev10best, avg10worst, sdev10worst, homo
        real, allocatable    :: avgs(:), sdevs(:)
        integer, allocatable :: order(:)
        ncls_here = b%a%get_ncls(p%minp)
        if( p%utst <= 3 )then ! we are interested in angular variations
            call b%a%read_spidoc(eulnams(p%utst))
            allocate(avgs(ncls_here), sdevs(ncls_here), order(ncls_here), stat=alloc_stat)
            call alloc_err('utst_calc_stat; simple_cluster', alloc_stat)
            ! collect statistics
            cnt = 0
            do i=1,p%ncls
                if( b%a%get_clspop(i) >= p%minp )then
                    cnt = cnt+1
                    order(cnt) = cnt
                    call b%a%class_dist_stat(i, avgs(cnt), sdevs(cnt))
                endif
            end do
            call hpsort(ncls_here, avgs, order)
            ! median
            med     = avgs(ncls_here/2)
            sdevmed = sdevs(order(ncls_here/2))
            ! average
            avg  = sum(avgs)/real(ncls_here)
            sdev = sum(sdevs)/real(ncls_here)
            ! 10 best
            avg10best  = 0.
            sdev10best = 0.
            do i=1,10
                avg10best  = avg10best+avgs(i)
                sdev10best = sdev10best+sdevs(order(i))
            end do
            avg10best  = avg10best/10.
            sdev10best = sdev10best/10.
            ! 10 worst
            avg10worst  = 0.
            sdev10worst = 0.
            do i=ncls_here,ncls_here-9,-1
                avg10worst  = avg10worst+avgs(i)
                sdev10worst = sdev10worst+sdevs(order(i))
            end do
            avg10worst  = avg10worst/10.
            sdev10worst = sdev10worst/10.
            write(*,'(a)') '**** UNIT TEST RESULTS, ANGULAR SPREAD WITHIN CLASSES ****'
            write(*,'(a)') '|  median   |  average  |  10 best  |  10 worst |'
            write(*,'(a)') '|avg    sdev|avg    sdev|avg    sdev|avg    sdev|'
            write(*,'(1x,f4.1,3x,f4.1,1x,f4.1,3x,f4.1,1x,f4.1,3x,f4.1,1x,f4.1,3x,f4.1)') &
            med, sdevmed, avg, sdev, avg10best, sdev10best, avg10worst, sdev10worst
            deallocate( avgs, sdevs, order )
        endif
        if( p%utst >= 3 )then ! we are interested in state homogeneity
            homo = 0.
            cnt  = 0
            do i=1,p%ncls
                if( b%a%get_clspop(i) >= p%minp )then
                    if( p%utst == 3 )then
                        homo = homo+max(b%a%class_calc_frac(i, [1,7000]),max(b%a%class_calc_frac(i,&
                        [7001,14000]),b%a%class_calc_frac(i, [14001,21000])))
                        cnt = cnt+1
                    else
                        homo = homo+max(b%a%class_calc_frac(i, [1,5000]),b%a%class_calc_frac(i, [5001,10000]))
                        cnt = cnt+1
                    endif
                endif
            end do
            homo = homo/real(cnt)
            write(*,'(A)') '**** UNIT TEST RESULTS, STATE HOMOGENEITY ****'
            write(*,'(a,1x,f4.3)') '|  homogeneity =', homo
        endif
    end subroutine
    
end program simple_cluster