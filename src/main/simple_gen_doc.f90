! automatic documentation generation
module simple_gen_doc
implicit none

contains

    subroutine print_doc_automask2D
        write(*,'(A)', advance='no') 'is a program for solvent flattening of class averages. The algorithm for backgro'
        write(*,'(A)', advance='no') 'und removal is based on low-pass filtering and binarization. First, the class av'
        write(*,'(A)', advance='no') 'erages are low-pass filtered to amsklp. Binary representatives are then generate'
        write(*,'(A)', advance='no') 'd by assigning foreground pixels using sortmeans. A cosine function softens the'
        write(*,'(A)', advance='no') 'edge of the binary mask before it is  multiplied with the unmasked input average'
        write(*,'(A)') 's to accomplish flattening'
        stop
    end subroutine print_doc_automask2D

    subroutine print_doc_binarise
        write(*,'(A)') 'is a program for binarisation of stacks and volumes'
        stop
    end subroutine print_doc_binarise

    subroutine print_doc_boxconvs
        write(*,'(A)', advance='no') 'is a program for averaging overlapping boxes across a micrograph in order to che'
        write(*,'(A)') 'ck if gain correction was appropriately done'
        stop
    end subroutine print_doc_boxconvs

    subroutine print_doc_cavgassemble
        write(*,'(A)', advance='no') 'is a program that assembles class averages when the clustering program (cluster2'
        write(*,'(A)') 'D) has been executed in distributed mode'
        stop
    end subroutine print_doc_cavgassemble

    subroutine print_doc_cenvol
        write(*,'(A)', advance='no') 'is a program for centering a volume and mapping the shift parameters back to the'
        write(*,'(A)') ' particle images'
        stop
    end subroutine print_doc_cenvol

    subroutine print_doc_check2D_conv
        write(*,'(A)', advance='no') 'is a program for checking if a cluster2D run has converged. The statistics outpu'
        write(*,'(A)', advance='no') 'tted include (1) the overlap between the distribution of parameters for succesiv'
        write(*,'(A)', advance='no') 'e runs. (2) The percentage of search space scanned, i.e. how many reference imag'
        write(*,'(A)', advance='no') 'es are evaluated on average. (3) The average correlation between the images and'
        write(*,'(A)', advance='no') 'their corresponding best matching reference section. If convergence to a local o'
        write(*,'(A)', advance='no') 'ptimum is achieved, the fraction increases. Convergence is achieved if the param'
        write(*,'(A)', advance='no') 'eter distribution overlap is larger than 0.95 and more than 99% of the reference'
        write(*,'(A)') ' sections need to be searched to find an improving solution'
        stop
    end subroutine print_doc_check2D_conv

    subroutine print_doc_check3D_conv
        write(*,'(A)', advance='no') 'is a program for checking if a PRIME3D run has converged. The statistics outputt'
        write(*,'(A)', advance='no') 'ed include (1) angle of feasible region, which is proportional to the angular re'
        write(*,'(A)', advance='no') 'solution of the set of discrete projection directions being searched. (2) The av'
        write(*,'(A)', advance='no') 'erage angular distance between orientations in the present and previous iteratio'
        write(*,'(A)', advance='no') 'n. In the early iterations, the distance is large because a diverse set of orien'
        write(*,'(A)', advance='no') 'tations is explored. If convergence to a local optimum is achieved, the distance'
        write(*,'(A)', advance='no') ' decreases. (3) The percentage of search space scanned, i.e. how many reference'
        write(*,'(A)', advance='no') 'images are evaluated on average. (4) The average correlation between the images'
        write(*,'(A)', advance='no') 'and their corresponding best matching reference sections. (5) The average standa'
        write(*,'(A)', advance='no') 'rd deviation of the Euler angles. Convergence is achieved if the angular distanc'
        write(*,'(A)', advance='no') 'e between the orientations in successive iterations falls significantly below th'
        write(*,'(A)', advance='no') 'e angular resolution of the search space and more than 99% of the reference sect'
        write(*,'(A)') 'ions need to be matched on average'
        stop
    end subroutine print_doc_check3D_conv

    subroutine print_doc_check_box
        write(*,'(A)', advance='no') 'is a program for checking the image dimensions of MRC and SPIDER stacks and volu'
        write(*,'(A)') 'mes'
        stop
    end subroutine print_doc_check_box

    subroutine print_doc_check_nptcls
        write(*,'(A)') 'is a program for checking the number of images in MRC and SPIDER stacks'
        stop
    end subroutine print_doc_check_nptcls

    subroutine print_doc_cluster2D
        write(*,'(A)', advance='no') 'is a reference-free 2D alignment/clustering algorithm adopted from the prime3D p'
        write(*,'(A)') 'robabilistic ab initio 3D reconstruction algorithm'
        stop
    end subroutine print_doc_cluster2D

    subroutine print_doc_cluster_cavgs
        write(*,'(A)', advance='no') 'is a program for analyzing class averages with affinity propagation, in order to'
        write(*,'(A)') ' get a better understanding of the view distribution'
        stop
    end subroutine print_doc_cluster_cavgs

    subroutine print_doc_cluster_oris
        write(*,'(A)') 'is a program for clustering orientations based on geodesic distance'
        stop
    end subroutine print_doc_cluster_oris

    subroutine print_doc_cluster_smat
        write(*,'(A)', advance='no') 'is a program for clustering a similarity matrix and use an combined cluster vali'
        write(*,'(A)', advance='no') 'dation index to assess the quality of the clustering based on the number of clus'
        write(*,'(A)') 'ters'
        stop
    end subroutine print_doc_cluster_smat

    subroutine print_doc_comlin_smat
        write(*,'(A)', advance='no') 'is a program for creating a similarity matrix based on common line correlation.'
        write(*,'(A)', advance='no') 'The idea being that it should be possible to cluster images based on their 3D si'
        write(*,'(A)', advance='no') 'milarity witout having a 3D model by only operating on class averages and find a'
        write(*,'(A)') 'verages that fit well together in 3D'
        stop
    end subroutine print_doc_comlin_smat

    subroutine print_doc_convert
        write(*,'(A)') 'is a program for converting between SPIDER and MRC formats'
        stop
    end subroutine print_doc_convert

    subroutine print_doc_corrcompare
        write(*,'(A)', advance='no') 'is a program for comparing stacked images using real-space and Fourier-based app'
        write(*,'(A)') 'roaches'
        stop
    end subroutine print_doc_corrcompare

    subroutine print_doc_ctf_estimate
        write(*,'(A)') 'is a SIMPLE program for fitting the CTF'
        stop
    end subroutine print_doc_ctf_estimate

    subroutine print_doc_ctffind
        write(*,'(A)') 'is a wrapper program for CTFFIND4 (Grigorieff lab)'
        stop
    end subroutine print_doc_ctffind

    subroutine print_doc_ctfops
        write(*,'(A)') 'is a program for applying CTF to stacked images'
        stop
    end subroutine print_doc_ctfops

    subroutine print_doc_dock_volpair
        write(*,'(A)') 'is a program for docking a pair of volumes. vol1 is reference and vol2 target.'
        stop
    end subroutine print_doc_dock_volpair

    subroutine print_doc_dsymsrch
        write(*,'(A)', advance='no') 'is a program for identifying rotational symmetries in class averages of D-symmet'
        write(*,'(A)') 'ric molecules and generating a cylinder that matches the shape.'
        stop
    end subroutine print_doc_dsymsrch

    subroutine print_doc_eo_volassemble
        write(*,'(A)', advance='no') 'is a program that assembles volume(s) when the reconstruction program (recvol wi'
        write(*,'(A)', advance='no') 'th eo=yes) has been executed in distributed mode. inner applies a soft-edged inn'
        write(*,'(A)', advance='no') 'er mask. An inner mask is used for icosahedral virus reconstruction, because the'
        write(*,'(A)', advance='no') ' DNA or RNA core is often unordered and  if not removed it may negatively impact'
        write(*,'(A)', advance='no') ' the alignment. The width parameter controls the fall-off of the edge of the inn'
        write(*,'(A)') 'er mask'
        stop
    end subroutine print_doc_eo_volassemble

    subroutine print_doc_extract
        write(*,'(A)', advance='no') 'is a program that extracts particle images from DDD movies or integrated movies.'
        write(*,'(A)', advance='no') ' Boxfiles are assumed to be in EMAN format but we provide a conversion script (r'
        write(*,'(A)', advance='no') 'elion2emanbox.pl) for *.star files containing particle coordinates obtained with'
        write(*,'(A)', advance='no') ' Relion. The program creates one stack per movie frame as well as a stack of cor'
        write(*,'(A)', advance='no') 'rected framesums. In addition to single-particle image stacks, the program produ'
        write(*,'(A)', advance='no') 'ces a parameter file extract_params.ext that can be used in conjunction with oth'
        write(*,'(A)') 'er SIMPLE programs. We obtain CTF parameters with CTFFIND4'
        stop
    end subroutine print_doc_extract

    subroutine print_doc_filter
        write(*,'(A)') 'is a program for filtering stacks/volumes'
        stop
    end subroutine print_doc_filter

    subroutine print_doc_fsc
        write(*,'(A)', advance='no') 'is a program for calculating the FSC between the two input volumes. No modificat'
        write(*,'(A)', advance='no') 'ions are done to the volumes, which allow you to test different masking options'
        write(*,'(A)') 'and see how they affect the FSCs'
        stop
    end subroutine print_doc_fsc

    subroutine print_doc_image_diff
        write(*,'(A)') 'is a program for comparing stacked images using differences'
        stop
    end subroutine print_doc_image_diff

    subroutine print_doc_image_smat
        write(*,'(A)', advance='no') 'is a program for creating a similarity matrix based on common line correlation.'
        write(*,'(A)', advance='no') 'The idea being that it should be possible to cluster images based on their 3D si'
        write(*,'(A)', advance='no') 'milarity witout having a 3D model by only operating on class averages and find a'
        write(*,'(A)') 'verages that fit well together in 3D'
        stop
    end subroutine print_doc_image_smat

    subroutine print_doc_iminfo
        write(*,'(A)', advance='no') 'is a program for printing header information in MRC and SPIDER stacks and volume'
        write(*,'(A)') 's'
        stop
    end subroutine print_doc_iminfo

    subroutine print_doc_intgpeaks
        write(*,'(A)') 'is a program for'
        stop
    end subroutine print_doc_intgpeaks

    subroutine print_doc_makecavgs
        write(*,'(A)', advance='no') 'is used  to produce class averages or initial random references for cluster2D ex'
        write(*,'(A)') 'ecution.'
        stop
    end subroutine print_doc_makecavgs

    subroutine print_doc_makedeftab
        write(*,'(A)', advance='no') 'is a program for creating a SIMPLE conformant file of CTF parameter values (deft'
        write(*,'(A)', advance='no') 'ab). Input is either an earlier SIMPLE deftab/oritab. The purpose is to get the'
        write(*,'(A)', advance='no') 'kv, cs, and fraca parameters as part of the CTF input doc as that is the new con'
        write(*,'(A)', advance='no') 'vention. The other alternative is to input a plain text file with CTF parameters'
        write(*,'(A)', advance='no') ' dfx, dfy, angast, phshift according to the Frealign convention. Unit conversion'
        write(*,'(A)', advance='no') 's are dealt with using optional variables. The units refer to the units in the i'
        write(*,'(A)') 'nputted document'
        stop
    end subroutine print_doc_makedeftab

    subroutine print_doc_makeoris
        write(*,'(A)', advance='no') 'is a program for making SIMPLE orientation/parameter files (text files containin'
        write(*,'(A)', advance='no') 'g input parameters and/or parameters estimated by cluster2D or refine3D). The pr'
        write(*,'(A)', advance='no') 'ogram generates random Euler angles e1.in.[0,360], e2.in.[0,180], and e3.in.[0,3'
        write(*,'(A)', advance='no') '60] and random origin shifts x.in.[-trs,yrs] and y.in.[-trs,yrs]. If ndiscrete i'
        write(*,'(A)', advance='no') 's set to an integer number > 0, the orientations produced are randomly sampled f'
        write(*,'(A)', advance='no') 'rom the set of ndiscrete quasi-even projection directions, and the in-plane para'
        write(*,'(A)', advance='no') 'meters are assigned randomly. If even=yes, then all nptcls orientations are assi'
        write(*,'(A)', advance='no') 'gned quasi-even projection directions,and random in-plane parameters. If nstates'
        write(*,'(A)', advance='no') ' is set to some integer number > 0, then states are assigned randomly .in.[1,nst'
        write(*,'(A)', advance='no') 'ates]. If zero=yes in this mode of execution, the projection directions are zero'
        write(*,'(A)', advance='no') 'ed and only the in-plane parameters are kept intact. If errify=yes and astigerr'
        write(*,'(A)', advance='no') 'is defined, then uniform random astigmatism errors are introduced .in.[-astigerr'
        write(*,'(A)') ',astigerr]'
        stop
    end subroutine print_doc_makeoris

    subroutine print_doc_makepickrefs
        write(*,'(A)') 'is a program for generating references for template-based particle picking'
        stop
    end subroutine print_doc_makepickrefs

    subroutine print_doc_map2ptcls
        write(*,'(A)', advance='no') 'is a program for mapping parameters that have been obtained using class averages'
        write(*,'(A)') ' to the individual particle images'
        stop
    end subroutine print_doc_map2ptcls

    subroutine print_doc_mask
        write(*,'(A)', advance='no') 'is a program for masking images and volumes. If you want to mask your images wit'
        write(*,'(A)') 'h a spherical mask with a soft falloff, set msk to the radius in pixels'
        stop
    end subroutine print_doc_mask

    subroutine print_doc_masscen
        write(*,'(A)') 'is a program for centering images acccording to their centre of mass'
        stop
    end subroutine print_doc_masscen

    subroutine print_doc_merge_algndocs
        write(*,'(A)', advance='no') 'is a program for merging alignment documents from SIMPLE runs in distributed mod'
        write(*,'(A)') 'e'
        stop
    end subroutine print_doc_merge_algndocs

    subroutine print_doc_merge_nnmat
        write(*,'(A)', advance='no') 'is a program for merging partial nearest neighbour matrices calculated in distri'
        write(*,'(A)') 'buted mode'
        stop
    end subroutine print_doc_merge_nnmat

    subroutine print_doc_merge_similarities
        write(*,'(A)', advance='no') 'is a program for merging similarities calculated between pairs of objects into a'
        write(*,'(A)') ' similarity matrix that can be inputted to cluster_smat'
        stop
    end subroutine print_doc_merge_similarities

    subroutine print_doc_motion_correct
        write(*,'(A)', advance='no') 'is a program for movie alignment or motion_correctring based the same principal'
        write(*,'(A)', advance='no') 'strategy as Grigorieffs program (hence the name). There are two important differ'
        write(*,'(A)', advance='no') 'ences: automatic weighting of the frames using a correlation-based M-estimator a'
        write(*,'(A)', advance='no') 'nd continuous optimisation of the shift parameters. Input is a textfile with abs'
        write(*,'(A)', advance='no') 'olute paths to movie files in addition to a few input parameters, some of which'
        write(*,'(A)', advance='no') 'deserve a comment. If dose_rate and exp_time are given the individual frames wil'
        write(*,'(A)', advance='no') 'l be low-pass filtered accordingly (dose-weighting strategy). If scale is given,'
        write(*,'(A)', advance='no') ' the movie will be Fourier cropped according to the down-scaling factor (for sup'
        write(*,'(A)', advance='no') 'er-resolution movies). If nframesgrp is given the frames will be pre-averaged in'
        write(*,'(A)', advance='no') ' the given chunk size (Falcon 3 movies). If fromf/tof are given, a contiguous su'
        write(*,'(A)') 'bset of frames will be averaged without any dose-weighting applied.'
        stop
    end subroutine print_doc_motion_correct

    subroutine print_doc_motion_correct_ctffind
        write(*,'(A)') 'is a pipelined motion_correct + ctffind program'
        stop
    end subroutine print_doc_motion_correct_ctffind

    subroutine print_doc_multiptcl_init
        write(*,'(A)', advance='no') 'is a program for generating random initial models for initialisation of PRIME3D'
        write(*,'(A)') 'when run in multiparticle mode'
        stop
    end subroutine print_doc_multiptcl_init

    subroutine print_doc_noiseimgs
        write(*,'(A)') 'is a program for generating pure noise images'
        stop
    end subroutine print_doc_noiseimgs

    subroutine print_doc_norm
        write(*,'(A)', advance='no') 'is a program for normalization of MRC or SPIDER stacks and volumes. If you want'
        write(*,'(A)', advance='no') 'to normalise your images inputted with stk, set norm=yes. If you want to perform'
        write(*,'(A)', advance='no') ' noise normalisation of the images set noise_norm=yes given a mask radius msk (p'
        write(*,'(A)', advance='no') 'ixels). If you want to normalise your images or volume (vol1) with respect to th'
        write(*,'(A)') 'eir power spectrum set shell_norm=yes'
        stop
    end subroutine print_doc_norm

    subroutine print_doc_npeaks
        write(*,'(A)', advance='no') 'is a program for checking the number of nonzero orientation weights (number of c'
        write(*,'(A)') 'orrelation peaks included in the weighted reconstruction)'
        stop
    end subroutine print_doc_npeaks

    subroutine print_doc_nspace
        write(*,'(A)', advance='no') 'is a program for calculating the expected resolution obtainable with different v'
        write(*,'(A)', advance='no') 'alues of nspace (number of discrete projection directions used for discrete sear'
        write(*,'(A)') 'ch)'
        stop
    end subroutine print_doc_nspace

    subroutine print_doc_orisops
        write(*,'(A)', advance='no') 'is a program for analyzing SIMPLE orientation/parameter files (text files contai'
        write(*,'(A)', advance='no') 'ning input parameters and/or parameters estimated by cluster2D or refine3D). If'
        write(*,'(A)', advance='no') 'only oritab is inputted, there are a few options available. If errify=yes, then'
        write(*,'(A)', advance='no') 'the program introduces uniform random angular errors .in.[-angerr,angerr], and u'
        write(*,'(A)', advance='no') 'niform origin shift errors .in.[-sherr,sherr], and uniform random defocus errors'
        write(*,'(A)', advance='no') ' .in.[-dferr,dferr]. If nstates > 1 then random states are assigned .in.[1,nstat'
        write(*,'(A)', advance='no') 'es]. If mirr=2d, then the Euler angles in oritab are mirrored according to the r'
        write(*,'(A)', advance='no') 'elation e1=e1, e2=180.+e2, e3=-e3. If mirr=3d, then the Euler angles in oritab a'
        write(*,'(A)', advance='no') 're mirrored according to the relation R=M(M*R), where R is the rotation matrix c'
        write(*,'(A)', advance='no') 'alculated from the Euler angle triplet and M is a 3D reflection matrix (like a u'
        write(*,'(A)', advance='no') 'nit matrix but with the 3,3-component sign swapped). If e1, e2, or e3 is inputte'
        write(*,'(A)', advance='no') 'd, the orientations in oritab are rotated correspondingly. If you input state as'
        write(*,'(A)', advance='no') ' well, you rotate only the orientations assigned to state state. If mul is defin'
        write(*,'(A)', advance='no') 'ed, you multiply the origin shifts with mul. If zero=yes, then the shifts are ze'
        write(*,'(A)', advance='no') 'roed. If none of the above described parameter are defined, and oritab is still'
        write(*,'(A)', advance='no') 'defined, the program projects the 3D orientation into the xy-plane and plots the'
        write(*,'(A)') ' resulting vector (this is useful for checking orientation coverage)'
        stop
    end subroutine print_doc_orisops

    subroutine print_doc_oristats
        write(*,'(A)', advance='no') 'is a program for analyzing SIMPLE orientation/parameter files (text files contai'
        write(*,'(A)', advance='no') 'ning input parameters and/or parameters estimated by cluster2D or refine3D). If'
        write(*,'(A)', advance='no') 'two orientation tables (oritab and oritab2) are inputted, the program provides s'
        write(*,'(A)', advance='no') 'tatistics of the distances between the orientations in the two documents. These'
        write(*,'(A)', advance='no') 'statistics include the sum of angular distances between the orientations, the av'
        write(*,'(A)', advance='no') 'erage angular distance between the orientations, the standard deviation of angul'
        write(*,'(A)') 'ar distances, the minimum angular distance, and the maximum angular distance'
        stop
    end subroutine print_doc_oristats

    subroutine print_doc_pick
        write(*,'(A)') 'is a template-based picker program'
        stop
    end subroutine print_doc_pick

    subroutine print_doc_postprocess
        write(*,'(A)') 'is a program for post-processing of volumes'
        stop
    end subroutine print_doc_postprocess

    subroutine print_doc_powerspecs
        write(*,'(A)') 'is a program for generating powerspectra from a stack or filetable'
        stop
    end subroutine print_doc_powerspecs

    subroutine print_doc_preproc
        write(*,'(A)') 'is a program that executes motion_correct, ctffind and pick in sequence'
        stop
    end subroutine print_doc_preproc

    subroutine print_doc_prime3D_init
        write(*,'(A)', advance='no') 'is a program for generating a random initial model for initialisation of PRIME3D'
        write(*,'(A)', advance='no') '. If the data set is large (>5000 images), generating a random model can be slow'
        write(*,'(A)', advance='no') '. To speedup, set nran to some smaller number, resulting in nran images selected'
        write(*,'(A)') ' randomly for reconstruction'
        stop
    end subroutine print_doc_prime3D_init

    subroutine print_doc_print_cmd_dict
        write(*,'(A)') 'is a program for printing the command line key dictonary'
        stop
    end subroutine print_doc_print_cmd_dict

    subroutine print_doc_print_dose_weights
        write(*,'(A)') 'is a program for printing the dose weights applied to individual frames'
        stop
    end subroutine print_doc_print_dose_weights

    subroutine print_doc_print_fsc
        write(*,'(A)') 'is a program for printing the binary FSC files produced by PRIME3D'
        stop
    end subroutine print_doc_print_fsc

    subroutine print_doc_print_magic_boxes
        write(*,'(A)') 'is a program for printing magic box sizes (fast FFT)'
        stop
    end subroutine print_doc_print_magic_boxes

    subroutine print_doc_print_project_info
        write(*,'(A)') 'prints information abourt a *.simple project file'
        stop
    end subroutine print_doc_print_project_info

    subroutine print_doc_project2txt
        write(*,'(A)') 'converts a binary *.simple project file to a text oritab'
        stop
    end subroutine print_doc_project2txt

    subroutine print_doc_projvol
        write(*,'(A)', advance='no') 'is a program for projecting a volume using interpolation in Fourier space. Input'
        write(*,'(A)', advance='no') ' is a SPIDER or MRC volume. Output is a stack of projection images of the same f'
        write(*,'(A)', advance='no') 'ormat as the inputted volume. Projections are generated by extraction of central'
        write(*,'(A)', advance='no') ' sections from the Fourier volume and back transformation of the 2D FTs. nspace'
        write(*,'(A)', advance='no') 'controls the number of projection images generated with quasi-even projection di'
        write(*,'(A)', advance='no') 'rections. The oritab parameter allows you to input the orientations that you wis'
        write(*,'(A)', advance='no') 'h to have your volume projected in. If rnd=yes, random rather than quasi-even pr'
        write(*,'(A)', advance='no') 'ojections are generated, trs then controls the halfwidth of the random origin sh'
        write(*,'(A)', advance='no') 'ift. Less commonly used parameters are pgrp, which controls the point-group symm'
        write(*,'(A)', advance='no') 'etry c (rotational), d (dihedral), t (tetrahedral), o (octahedral) or i (icosahe'
        write(*,'(A)', advance='no') 'dral). The point-group symmetry is used to restrict the set of projections to wi'
        write(*,'(A)') 'thin the asymmetric unit. neg inverts the contrast of the projections.'
        stop
    end subroutine print_doc_projvol

    subroutine print_doc_rank_cavgs
        write(*,'(A)', advance='no') 'is a program for ranking class averages by decreasing population, given the stac'
        write(*,'(A)', advance='no') 'k of class averages (stk argument) and the 2D orientations document (oritab) gen'
        write(*,'(A)') 'erated by cluster2D'
        stop
    end subroutine print_doc_rank_cavgs

    subroutine print_doc_recvol
        write(*,'(A)', advance='no') 'is a program for reconstructing volumes from MRC and SPIDER stacks, given input'
        write(*,'(A)', advance='no') 'orientations and state assignments. The algorithm is based on direct Fourier inv'
        write(*,'(A)', advance='no') 'ersion with a Kaiser-Bessel (KB) interpolation kernel. This window function redu'
        write(*,'(A)', advance='no') 'ces the real-space ripple artifacts associated with direct moving windowed-sinc'
        write(*,'(A)', advance='no') 'interpolation. The feature sought when implementing this algorithm was to enable'
        write(*,'(A)', advance='no') ' quick, reliable reconstruction from aligned individual particle images. mul is'
        write(*,'(A)', advance='no') 'used to scale the origin shifts if down-sampled were used for alignment and the'
        write(*,'(A)', advance='no') 'original images are used for reconstruction. ctf=yes or ctf=flip turns on the Wi'
        write(*,'(A)', advance='no') 'ener restoration. If the images were phase-flipped set ctf=flip. The inner param'
        write(*,'(A)', advance='no') 'eter controls the radius of the soft-edged mask used to remove the unordered DNA'
        write(*,'(A)') '/RNA core of spherical icosahedral viruses'
        stop
    end subroutine print_doc_recvol

    subroutine print_doc_refine3D
        write(*,'(A)', advance='no') 'is an ab inito reconstruction/refinement program based on probabilistic projecti'
        write(*,'(A)', advance='no') 'on matching. There are a daunting number of options in refine3D. If you are proc'
        write(*,'(A)', advance='no') 'essing class averages we recommend that you instead use the simple_distr_exec pr'
        write(*,'(A)') 'g= initial_3Dmodel route.'
        stop
    end subroutine print_doc_refine3D

    subroutine print_doc_res
        write(*,'(A)', advance='no') 'is a program for checking the low-pass resolution limit for a given Fourier inde'
        write(*,'(A)') 'x'
        stop
    end subroutine print_doc_res

    subroutine print_doc_resmask
        write(*,'(A)') 'is a program for 3D envelope masking for resolution estimation'
        stop
    end subroutine print_doc_resmask

    subroutine print_doc_resrange
        write(*,'(A)', advance='no') 'is a program for estimating the resolution range used in the heuristic resolutio'
        write(*,'(A)', advance='no') 'n-stepping scheme in the PRIME3D initial model production procedure. The initial'
        write(*,'(A)', advance='no') ' low-pass limit is set so that each image receives ten nonzero orientation weigh'
        write(*,'(A)', advance='no') 'ts. When quasi-convergence has been reached, the limit is updated one Fourier in'
        write(*,'(A)', advance='no') 'dex at the time, until PRIME reaches the condition where six nonzero orientation'
        write(*,'(A)', advance='no') ' weights are assigned to each image. FSC-based filtering is unfortunately not po'
        write(*,'(A)', advance='no') 'ssible to do in the ab initio 3D reconstruction step, because when the orientati'
        write(*,'(A)', advance='no') 'ons are mostly random, the FSC overestimates the resolution. This program is use'
        write(*,'(A)', advance='no') 'd internally when executing PRIME in distributed mode. We advise you to check th'
        write(*,'(A)', advance='no') 'e starting and stopping low-pass limits before executing PRIME3D using this prog'
        write(*,'(A)', advance='no') 'ram. The resolution range estimate depends on the molecular diameter, which is e'
        write(*,'(A)', advance='no') 'stimated based on the box size. If you want to override this estimate, set moldi'
        write(*,'(A)', advance='no') 'am to the desired value (in A). This may be necessary if your images have a lot'
        write(*,'(A)', advance='no') 'of background padding. However, for starting model generation it is probably bet'
        write(*,'(A)', advance='no') 'ter to clip the images snugly around the particle, because smaller images equal'
        write(*,'(A)') 'less computation'
        stop
    end subroutine print_doc_resrange

    subroutine print_doc_rotmats2oris
        write(*,'(A)', advance='no') 'converts a text file (9 records per line) describing rotation matrices into a SI'
        write(*,'(A)') 'MPLE oritab'
        stop
    end subroutine print_doc_rotmats2oris

    subroutine print_doc_scale
        write(*,'(A)', advance='no') 'is a program that provides re-scaling and clipping routines for MRC or SPIDER st'
        write(*,'(A)') 'acks and volumes'
        stop
    end subroutine print_doc_scale

    subroutine print_doc_select
        write(*,'(A)') 'is a program for selecting files based on image correlation matching'
        stop
    end subroutine print_doc_select

    subroutine print_doc_select_frames
        write(*,'(A)') 'is a program for selecting contiguous segments of frames from DDD movies'
        stop
    end subroutine print_doc_select_frames

    subroutine print_doc_shift
        write(*,'(A)') 'is a program for shifting a stack according to shifts in oritab'
        stop
    end subroutine print_doc_shift

    subroutine print_doc_simimgs
        write(*,'(A)', advance='no') 'is a program for simulating cryo-EM images. It is not a very sophisticated simul'
        write(*,'(A)', advance='no') 'ator, but it is nevertheless useful for testing purposes. It does not do any mul'
        write(*,'(A)', advance='no') 'ti-slice simulation and it cannot be used for simulating molecules containing he'
        write(*,'(A)', advance='no') 'avy atoms. It does not even accept a PDB file as an input. Input is a cryo-EM ma'
        write(*,'(A)', advance='no') 'p, which we usually generate from a PDB file using EMANs program pdb2mrc. simimg'
        write(*,'(A)', advance='no') 's then projects the volume using Fourier interpolation, adds 20% of the total no'
        write(*,'(A)', advance='no') 'ise to the images (pink noise), Fourier transforms them, and multiplies them wit'
        write(*,'(A)', advance='no') 'h astigmatic CTF and B-factor. The images are inverse FTed before the remaining'
        write(*,'(A)') '80% of the noise (white noise) is added'
        stop
    end subroutine print_doc_simimgs

    subroutine print_doc_simmovie
        write(*,'(A)', advance='no') 'is a program for crude simulation of a DDD movie. Input is a set of projection i'
        write(*,'(A)', advance='no') 'mages to place. Movie frames are then generated related by randomly shifting the'
        write(*,'(A)') ' base image and applying noise'
        stop
    end subroutine print_doc_simmovie

    subroutine print_doc_simsubtomo
        write(*,'(A)') 'is a program for crude simulation of a subtomograms'
        stop
    end subroutine print_doc_simsubtomo

    subroutine print_doc_split
        write(*,'(A)', advance='no') 'is a program for splitting of image stacks into partitions for parallel executio'
        write(*,'(A)') 'n. This is done to reduce I/O latency'
        stop
    end subroutine print_doc_split

    subroutine print_doc_split_pairs
        write(*,'(A)', advance='no') 'is a program for splitting calculations between pairs of objects into balanced p'
        write(*,'(A)') 'artitions'
        stop
    end subroutine print_doc_split_pairs

    subroutine print_doc_stack
        write(*,'(A)') 'is a program for stacking individual images or multiple stacks into one'
        stop
    end subroutine print_doc_stack

    subroutine print_doc_stackops
        write(*,'(A)', advance='no') 'is a program that provides standard single-particle image processing routines th'
        write(*,'(A)', advance='no') 'at are applied to MRC or SPIDER stacks. If you want to extract a particular stat'
        write(*,'(A)', advance='no') 'e, give an alignment document (oritab) and set state to the state that you want'
        write(*,'(A)', advance='no') 'to extract. If you want to select the fraction of best particles (according to t'
        write(*,'(A)', advance='no') 'he goal function), input an alignment doc (oritab) and set frac. You can combine'
        write(*,'(A)', advance='no') ' the state and frac options. If you want to apply noise to images, give the desi'
        write(*,'(A)', advance='no') 'red signal-to-noise ratio via snr. If you want to calculate the autocorrelation'
        write(*,'(A)', advance='no') 'function of your images set acf=yes. If you want to extract a contiguous subset'
        write(*,'(A)', advance='no') 'of particle images from the stack, set fromp and top. If you want to fish out a'
        write(*,'(A)', advance='no') 'number of particle images from your stack at random, set nran to some nonzero in'
        write(*,'(A)', advance='no') 'teger number less than nptcls. With avg=yes the global average of the inputted s'
        write(*,'(A)', advance='no') 'tack is calculated. If you define nframesgrp to some integer number larger than'
        write(*,'(A)', advance='no') 'one averages with chunk sizes of nframesgrp are produced, which may be useful fo'
        write(*,'(A)', advance='no') 'r analysis of dose-fractionated image series. neg inverts the contrast of the im'
        write(*,'(A)') 'ages'
        stop
    end subroutine print_doc_stackops

    subroutine print_doc_stktabinfo
        write(*,'(A)') 'is a program for for printing information about stktab'
        stop
    end subroutine print_doc_stktabinfo

    subroutine print_doc_sym_aggregate
        write(*,'(A)', advance='no') 'is a program for robust identifiaction of the symmetry axis of a map using image'
        write(*,'(A)') '-to-volume simiarity validation of the axis'
        stop
    end subroutine print_doc_sym_aggregate

    subroutine print_doc_symsrch
        write(*,'(A)', advance='no') 'is a program for searching for the principal symmetry axis of a volume reconstru'
        write(*,'(A)', advance='no') 'cted without assuming any point-group symmetry. The program takes as input an as'
        write(*,'(A)', advance='no') 'ymmetrical 3D reconstruction. The alignment document for all the particle images'
        write(*,'(A)', advance='no') ' that have gone into the 3D reconstruction and the desired point-group symmetry'
        write(*,'(A)', advance='no') 'needs to be inputted. The 3D reconstruction is then projected in 50 (default opt'
        write(*,'(A)', advance='no') 'ion) even directions, common lines-based optimisation is used to identify the pr'
        write(*,'(A)', advance='no') 'incipal symmetry axis, the rotational transformation is applied to the inputted'
        write(*,'(A)', advance='no') 'orientations, and a new alignment document is produced. Input this document to r'
        write(*,'(A)', advance='no') 'ecvol together with the images and the point-group symmetry to generate a symmet'
        write(*,'(A)') 'rised map.'
        stop
    end subroutine print_doc_symsrch

    subroutine print_doc_tseries_backgr_subtr
        write(*,'(A)', advance='no') 'is a program for background subtraction in time-series data. The goal is to subt'
        write(*,'(A)', advance='no') 'ract the two graphene peaks @ 2.14 A and @ 1.23 A. This is done by band-pass fil'
        write(*,'(A)', advance='no') 'tering the background image, recommended (and default settings) are hp=5.0 lp=1.'
        write(*,'(A)') '1 and width=5.0.'
        stop
    end subroutine print_doc_tseries_backgr_subtr

    subroutine print_doc_tseries_extract
        write(*,'(A)', advance='no') 'is a program for creating overlapping chunks of nframesgrp frames from time-seri'
        write(*,'(A)') 'es data'
        stop
    end subroutine print_doc_tseries_extract

    subroutine print_doc_tseries_split
        write(*,'(A)') 'is a program for splitting a time-series stack and its associated orientations'
        stop
    end subroutine print_doc_tseries_split

    subroutine print_doc_tseries_track
        write(*,'(A)') 'is a program for particle tracking in time-series data'
        stop
    end subroutine print_doc_tseries_track

    subroutine print_doc_txt2project
        write(*,'(A)') 'adds or replaces a text oritab in a binary *.simple project file'
        stop
    end subroutine print_doc_txt2project

    subroutine print_doc_vizoris
        write(*,'(A)', advance='no') 'extract projection direction from an orientation direction for visualization in'
        write(*,'(A)') 'UCSF Chimera'
        stop
    end subroutine print_doc_vizoris

    subroutine print_doc_volassemble
        write(*,'(A)', advance='no') 'is a program that assembles volume(s) when the reconstruction program (recvol) h'
        write(*,'(A)', advance='no') 'as been executed in distributed mode. odd is used to assemble the odd reconstruc'
        write(*,'(A)', advance='no') 'tion, even is used to assemble the even reconstruction, eo is used to assemble b'
        write(*,'(A)', advance='no') 'oth the even and the odd reconstruction and state is used to assemble the inputt'
        write(*,'(A)', advance='no') 'ed state. Normally, you do not fiddle with these parameters. They are used inter'
        write(*,'(A)') 'nally'
        stop
    end subroutine print_doc_volassemble

    subroutine print_doc_volaverager
        write(*,'(A)') 'is a program for averaging volumes according to state label in oritab'
        stop
    end subroutine print_doc_volaverager

    subroutine print_doc_volops
        write(*,'(A)', advance='no') 'provides standard single-particle image processing routines that are applied to'
        write(*,'(A)') 'MRC or SPIDER volumes'
        stop
    end subroutine print_doc_volops

    subroutine print_doc_volume_smat
        write(*,'(A)') 'is a program for creating a similarity matrix based on volume2volume correlation'
        stop
    end subroutine print_doc_volume_smat

    subroutine print_doc_cluster2D_stream
        write(*,'(A)', advance='no') 'is a distributed workflow implementing reference-free 2D alignment/clustering al'
        write(*,'(A)', advance='no') 'gorithm adopted from the prime3D probabilistic ab initio 3D reconstruction algor'
        write(*,'(A)') 'ithm'
        stop
    end subroutine print_doc_cluster2D_stream

    subroutine print_doc_het
        write(*,'(A)') 'is a distributed workflow for heterogeneity analysis'
        stop
    end subroutine print_doc_het

    subroutine print_doc_het_refine
        write(*,'(A)') 'is a distributed workflow for heterogeneity analysis refinement'
        stop
    end subroutine print_doc_het_refine

    subroutine print_doc_initial_3Dmodel
        write(*,'(A)', advance='no') 'is a distributed workflow for generating an initial 3D model from class averages'
        write(*,'(A)') ' obtained with cluster2D'
        stop
    end subroutine print_doc_initial_3Dmodel

    subroutine print_doc_motion_correct_tomo
        write(*,'(A)', advance='no') 'is a distributed workflow for movie alignment or motion_correctring of tomograph'
        write(*,'(A)', advance='no') 'ic movies. Input is a textfile with absolute paths to movie files in addition to'
        write(*,'(A)', advance='no') ' a few input parameters, some of which deserve a comment. The exp_doc document s'
        write(*,'(A)', advance='no') 'hould contain per line exp_time=X and dose_rate=Y. It is asssumed that the input'
        write(*,'(A)', advance='no') ' list of movies (one per tilt) are ordered temporally. This is necessary for cor'
        write(*,'(A)', advance='no') 'rect dose-weighting of tomographic tilt series. If scale is given, the movie wil'
        write(*,'(A)', advance='no') 'l be Fourier cropped according to the down-scaling factor (for super-resolution'
        write(*,'(A)', advance='no') 'movies). If nframesgrp is given the frames will be pre-averaged in the given chu'
        write(*,'(A)') 'nk size (Falcon 3 movies).'
        stop
    end subroutine print_doc_motion_correct_tomo

    subroutine print_doc_scale_stk_parts
        write(*,'(A)') 'is a distributed workflow for scaling partial stacks'
        stop
    end subroutine print_doc_scale_stk_parts

    subroutine list_all_simple_programs
        write(*,'(A)') 'automask2D'
        write(*,'(A)') 'binarise'
        write(*,'(A)') 'boxconvs'
        write(*,'(A)') 'cavgassemble'
        write(*,'(A)') 'cenvol'
        write(*,'(A)') 'check2D_conv'
        write(*,'(A)') 'check3D_conv'
        write(*,'(A)') 'check_box'
        write(*,'(A)') 'check_nptcls'
        write(*,'(A)') 'cluster2D'
        write(*,'(A)') 'cluster_cavgs'
        write(*,'(A)') 'cluster_oris'
        write(*,'(A)') 'cluster_smat'
        write(*,'(A)') 'comlin_smat'
        write(*,'(A)') 'convert'
        write(*,'(A)') 'corrcompare'
        write(*,'(A)') 'ctf_estimate'
        write(*,'(A)') 'ctffind'
        write(*,'(A)') 'ctfops'
        write(*,'(A)') 'dock_volpair'
        write(*,'(A)') 'dsymsrch'
        write(*,'(A)') 'eo_volassemble'
        write(*,'(A)') 'extract'
        write(*,'(A)') 'filter'
        write(*,'(A)') 'fsc'
        write(*,'(A)') 'image_diff'
        write(*,'(A)') 'image_smat'
        write(*,'(A)') 'iminfo'
        write(*,'(A)') 'intgpeaks'
        write(*,'(A)') 'makecavgs'
        write(*,'(A)') 'makedeftab'
        write(*,'(A)') 'makeoris'
        write(*,'(A)') 'makepickrefs'
        write(*,'(A)') 'map2ptcls'
        write(*,'(A)') 'map2ptcls_doc'
        write(*,'(A)') 'mask'
        write(*,'(A)') 'masscen'
        write(*,'(A)') 'merge_algndocs'
        write(*,'(A)') 'merge_nnmat'
        write(*,'(A)') 'merge_similarities'
        write(*,'(A)') 'motion_correct'
        write(*,'(A)') 'motion_correct_ctffind'
        write(*,'(A)') 'multiptcl_init'
        write(*,'(A)') 'noiseimgs'
        write(*,'(A)') 'norm'
        write(*,'(A)') 'npeaks'
        write(*,'(A)') 'nspace'
        write(*,'(A)') 'orisops'
        write(*,'(A)') 'oristats'
        write(*,'(A)') 'pick'
        write(*,'(A)') 'postprocess'
        write(*,'(A)') 'powerspecs'
        write(*,'(A)') 'preproc'
        write(*,'(A)') 'prime3D_init'
        write(*,'(A)') 'print_cmd_dict'
        write(*,'(A)') 'print_dose_weights'
        write(*,'(A)') 'print_fsc'
        write(*,'(A)') 'print_magic_boxes'
        write(*,'(A)') 'print_project_info'
        write(*,'(A)') 'project2txt'
        write(*,'(A)') 'projvol'
        write(*,'(A)') 'rank_cavgs'
        write(*,'(A)') 'rec_test'
        write(*,'(A)') 'recvol'
        write(*,'(A)') 'refine3D'
        write(*,'(A)') 'res'
        write(*,'(A)') 'resmask'
        write(*,'(A)') 'resrange'
        write(*,'(A)') 'rotmats2oris'
        write(*,'(A)') 'scale'
        write(*,'(A)') 'select'
        write(*,'(A)') 'select_frames'
        write(*,'(A)') 'shift'
        write(*,'(A)') 'simimgs'
        write(*,'(A)') 'simmovie'
        write(*,'(A)') 'simsubtomo'
        write(*,'(A)') 'split'
        write(*,'(A)') 'split_pairs'
        write(*,'(A)') 'stack'
        write(*,'(A)') 'stackops'
        write(*,'(A)') 'stktabinfo'
        write(*,'(A)') 'sym_aggregate'
        write(*,'(A)') 'symsrch'
        write(*,'(A)') 'tseries_backgr_subtr'
        write(*,'(A)') 'tseries_extract'
        write(*,'(A)') 'tseries_split'
        write(*,'(A)') 'tseries_track'
        write(*,'(A)') 'txt2project'
        write(*,'(A)') 'vizoris'
        write(*,'(A)') 'volassemble'
        write(*,'(A)') 'volaverager'
        write(*,'(A)') 'volops'
        write(*,'(A)') 'volume_smat'
        stop
    end subroutine list_all_simple_programs

    subroutine list_all_simple_distr_programs
        write(*,'(A)') 'cluster2D'
        write(*,'(A)') 'cluster2D_stream'
        write(*,'(A)') 'comlin_smat'
        write(*,'(A)') 'ctf_estimate'
        write(*,'(A)') 'ctffind'
        write(*,'(A)') 'het'
        write(*,'(A)') 'het_refine'
        write(*,'(A)') 'initial_3Dmodel'
        write(*,'(A)') 'makecavgs'
        write(*,'(A)') 'motion_correct'
        write(*,'(A)') 'motion_correct_ctffind'
        write(*,'(A)') 'motion_correct_tomo'
        write(*,'(A)') 'pick'
        write(*,'(A)') 'powerspecs'
        write(*,'(A)') 'preproc'
        write(*,'(A)') 'prime3D_init'
        write(*,'(A)') 'recvol'
        write(*,'(A)') 'refine3D'
        write(*,'(A)') 'scale_stk_parts'
        write(*,'(A)') 'symsrch'
        write(*,'(A)') 'tseries_track'
        stop
    end subroutine list_all_simple_distr_programs

end module simple_gen_doc
