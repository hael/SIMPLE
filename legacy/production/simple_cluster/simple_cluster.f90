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
use simple_jiffys,      ! singleton
use simple_math,        only: hpsort, euclid
use simple_build,       only: build
use simple_params,      only: params
use simple_image,       only: image
use simple_procimgfile, only: make_pattern_stack, norm_imgfile
implicit none
type(params)       :: p
type(build)        :: b
integer            :: pcarecsz, iptcl
real, allocatable  :: AVG(:)
logical, parameter :: debug=.false.
if( command_argument_count() < 3 )then
    write(*,'(a)', advance='no') 'SIMPLE_CLUSTER stk=<stack.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)', advance='no') ' ncls=<nr of clusters> msk=<mask radius(in pixels)>'
    write(*,'(a)', advance='no') ' oritab=<PRIME3D doc> [nthr=<nr of openMP threads{1}>]'
    write(*,'(a)') '**less commonly used**' 
    write(*,'(a)') '[nvars=<nr of eigen vectors{30}>]'
    stop
endif
call parse_cmdline
call cmdcheckvar('stk',    1)
call cmdcheckvar('smpd',   2)
call cmdcheckvar('ncls',   3)
call cmdcheckvar('msk',    4)
call cmdcheckvar('oritab', 5)
call cmdcheck                ! checks args and prints cmdline to cmdline.dat
p = params()                 ! constants & derived constants produced
call b%build_general_tbox(p) ! general objects built
call b%build_cluster_tbox(p) ! objects for clustering built
write(*,'(A)') '**** CLUSTERING PROCEDURE BEGINS ****'
! create pattern stack (pcastk) for probabilistic PCA (pPCA)
call make_pattern_stack(p%stk, p%pcastk, p%msk, p%ncomps, pcarecsz, AVG, b%a)
! generate a pPCA solution using Maximum Likelihood estimation 
call b%pca%master(p%pcastk, pcarecsz, p%featstk, p%maxits)
call b%build_read_features(p)
if( debug )then
    do iptcl=1,p%nptcls
        print *, '>>>> PARTICLE: ', iptcl
        print *, b%features(iptcl,:)
    end do
endif
write(*,'(A)') '**** FINISHED PROBABILISTIC PCA ****'
! make a random clustering start
call b%a%rnd_cls(p%ncls, srch_inpl=.false.) ! .false. means no in-plane randomisation
call b%cenclust%srch_shc(p%maxits) ! maximum number of iterations is 500
! write output
call b%a%write('shc_clustering_ncls'//int2str(p%ncls)//'.txt')
write(*,'(A,1X,A)') '>>> WROTE CONVERGED CLUSTERING SOLUTION TO FILE:',&
&'shc_clustering_ncls'//int2str(p%ncls)//'.txt'
! end gracefully
call simple_end('**** SIMPLE_CLUSTER NORMAL STOP ****')
end program simple_cluster