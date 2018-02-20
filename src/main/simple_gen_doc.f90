! automatic documentation generation
module simple_gen_doc
implicit none

contains

    subroutine print_doc_center
        write(*,'(A)', advance='no') 'is a program for centering a volume and mapping the shift parameters back to the'
        write(*,'(A)') ' particle images'
        stop
    end subroutine print_doc_center

    subroutine print_doc_cluster_cavgs
        write(*,'(A)', advance='no') 'is a program for analyzing class averages with affinity propagation, in order to'
        write(*,'(A)') ' get a better understanding of the view distribution'
        stop
    end subroutine print_doc_cluster_cavgs

    subroutine print_doc_convert
        write(*,'(A)') 'is a program for converting between SPIDER and MRC formats'
        stop
    end subroutine print_doc_convert

    subroutine print_doc_ctfops
        write(*,'(A)') 'is a program for applying CTF to stacked images'
        stop
    end subroutine print_doc_ctfops

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

    subroutine print_doc_info_image
        write(*,'(A)', advance='no') 'is a program for printing header information in MRC and SPIDER stacks and volume'
        write(*,'(A)') 's'
        stop
    end subroutine print_doc_info_image

    subroutine print_doc_info_stktab
        write(*,'(A)') 'is a program for for printing information about stktab'
        stop
    end subroutine print_doc_info_stktab

    subroutine print_doc_make_deftab
        write(*,'(A)', advance='no') 'is a program for creating a SIMPLE conformant file of CTF parameter values (deft'
        write(*,'(A)', advance='no') 'ab). Input is either an earlier SIMPLE deftab/oritab. The purpose is to get the'
        write(*,'(A)', advance='no') 'kv, cs, and fraca parameters as part of the CTF input doc as that is the new con'
        write(*,'(A)', advance='no') 'vention. The other alternative is to input a plain text file with CTF parameters'
        write(*,'(A)', advance='no') ' dfx, dfy, angast, phshift according to the Frealign convention. Unit conversion'
        write(*,'(A)', advance='no') 's are dealt with using optional variables. The units refer to the units in the i'
        write(*,'(A)') 'nputted document'
        stop
    end subroutine print_doc_make_deftab

    subroutine print_doc_make_oris
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
    end subroutine print_doc_make_oris

    subroutine print_doc_make_pickrefs
        write(*,'(A)') 'is a program for generating references for template-based particle picking'
        stop
    end subroutine print_doc_make_pickrefs

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

    subroutine print_doc_normalize
        write(*,'(A)', advance='no') 'is a program for normalization of MRC or SPIDER stacks and volumes. If you want'
        write(*,'(A)', advance='no') 'to normalize your images inputted with stk, set norm=yes. If you want to perform'
        write(*,'(A)', advance='no') ' noise normalisation of the images set noise_norm=yes given a mask radius msk (p'
        write(*,'(A)', advance='no') 'ixels). If you want to normalize your images or volume (vol1) with respect to th'
        write(*,'(A)') 'eir power spectrum set shell_norm=yes'
        stop
    end subroutine print_doc_normalize

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

    subroutine print_doc_postprocess
        write(*,'(A)') 'is a program for post-processing of volumes'
        stop
    end subroutine print_doc_postprocess

    subroutine print_doc_print_cmd_dict
        write(*,'(A)') 'is a program for printing the command line key dictonary'
        stop
    end subroutine print_doc_print_cmd_dict

    subroutine print_doc_print_fsc
        write(*,'(A)') 'is a program for printing the binary FSC files produced by PRIME3D'
        stop
    end subroutine print_doc_print_fsc

    subroutine print_doc_print_magic_boxes
        write(*,'(A)') 'is a program for printing magic box sizes (fast FFT)'
        stop
    end subroutine print_doc_print_magic_boxes

    subroutine print_doc_project
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
    end subroutine print_doc_project

    subroutine print_doc_scale
        write(*,'(A)', advance='no') 'is a program that provides re-scaling and clipping routines for MRC or SPIDER st'
        write(*,'(A)') 'acks and volumes'
        stop
    end subroutine print_doc_scale

    subroutine print_doc_select
        write(*,'(A)') 'is a program for selecting files based on image correlation matching'
        stop
    end subroutine print_doc_select

    subroutine print_doc_shift
        write(*,'(A)') 'is a program for shifting a stack according to shifts in oritab'
        stop
    end subroutine print_doc_shift

    subroutine print_doc_simulate_movie
        write(*,'(A)', advance='no') 'is a program for crude simulation of a DDD movie. Input is a set of projection i'
        write(*,'(A)', advance='no') 'mages to place. Movie frames are then generated related by randomly shifting the'
        write(*,'(A)') ' base image and applying noise'
        stop
    end subroutine print_doc_simulate_movie

    subroutine print_doc_simulate_noise
        write(*,'(A)') 'is a program for generating pure noise images'
        stop
    end subroutine print_doc_simulate_noise

    subroutine print_doc_simulate_particles
        write(*,'(A)', advance='no') 'is a program for simulating cryo-EM images. It is not a very sophisticated simul'
        write(*,'(A)', advance='no') 'ator, but it is nevertheless useful for testing purposes. It does not do any mul'
        write(*,'(A)', advance='no') 'ti-slice simulation and it cannot be used for simulating molecules containing he'
        write(*,'(A)', advance='no') 'avy atoms. It does not even accept a PDB file as an input. Input is a cryo-EM ma'
        write(*,'(A)', advance='no') 'p, which we usually generate from a PDB file using EMANs program pdb2mrc. simula'
        write(*,'(A)', advance='no') 'te_particles then projects the volume using Fourier interpolation, adds 20% of t'
        write(*,'(A)', advance='no') 'he total noise to the images (pink noise), Fourier transforms them, and multipli'
        write(*,'(A)', advance='no') 'es them with astigmatic CTF and B-factor. The images are inverse FTed before the'
        write(*,'(A)') ' remaining 80% of the noise (white noise) is added'
        stop
    end subroutine print_doc_simulate_particles

    subroutine print_doc_simulate_subtomogram
        write(*,'(A)') 'is a program for crude simulation of a subtomograms'
        stop
    end subroutine print_doc_simulate_subtomogram

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

    subroutine print_doc_vizoris
        write(*,'(A)', advance='no') 'extract projection direction from an orientation direction for visualization in'
        write(*,'(A)') 'UCSF Chimera'
        stop
    end subroutine print_doc_vizoris

    subroutine print_doc_volops
        write(*,'(A)', advance='no') 'provides standard single-particle image processing routines that are applied to'
        write(*,'(A)') 'MRC or SPIDER volumes'
        stop
    end subroutine print_doc_volops

    subroutine print_doc_cluster2D
        write(*,'(A)', advance='no') 'is a distributed workflow implementing a reference-free 2D alignment/clustering'
        write(*,'(A)', advance='no') 'algorithm adopted from the prime3D probabilistic ab initio 3D reconstruction alg'
        write(*,'(A)') 'orithm'
        stop
    end subroutine print_doc_cluster2D

    subroutine print_doc_cluster2D_stream
        write(*,'(A)', advance='no') 'is a distributed workflow implementing reference-free 2D alignment/clustering al'
        write(*,'(A)', advance='no') 'gorithm adopted from the prime3D probabilistic ab initio 3D reconstruction algor'
        write(*,'(A)') 'ithm'
        stop
    end subroutine print_doc_cluster2D_stream

    subroutine print_doc_cluster3D
        write(*,'(A)') 'is a distributed workflow for heterogeneity analysis by 3D clustering'
        stop
    end subroutine print_doc_cluster3D

    subroutine print_doc_cluster3D_refine
        write(*,'(A)') 'is a distributed workflow for refinement of heterogeneity analysis by cluster3D'
        stop
    end subroutine print_doc_cluster3D_refine

    subroutine print_doc_comlin_smat
        write(*,'(A)', advance='no') 'is a distributed workflow for creating a similarity matrix based on common line'
        write(*,'(A)', advance='no') 'correlation. The idea being that it should be possible to cluster images based o'
        write(*,'(A)', advance='no') 'n their 3D similarity witout having a 3D model by only operating on class averag'
        write(*,'(A)') 'es and find averages that fit well together in 3D'
        stop
    end subroutine print_doc_comlin_smat

    subroutine print_doc_ctf_estimate
        write(*,'(A)') 'is a distributed SIMPLE workflow for fitting the CTF'
        stop
    end subroutine print_doc_ctf_estimate

    subroutine print_doc_ctffind
        write(*,'(A)') 'is a distributed workflow that wraps CTFFIND4 (Grigorieff lab)'
        stop
    end subroutine print_doc_ctffind

    subroutine print_doc_initial_3Dmodel
        write(*,'(A)', advance='no') 'is a distributed workflow for generating an initial 3D model from class averages'
        write(*,'(A)') ' obtained with cluster2D'
        stop
    end subroutine print_doc_initial_3Dmodel

    subroutine print_doc_make_cavgs
        write(*,'(A)', advance='no') 'is a distributed workflow used for producing class averages or initial random re'
        write(*,'(A)') 'ferences for cluster2D execution.'
        stop
    end subroutine print_doc_make_cavgs

    subroutine print_doc_motion_correct
        write(*,'(A)', advance='no') 'is a distributed workflow for movie alignment or motion_correctring based the sa'
        write(*,'(A)', advance='no') 'me principal strategy as Grigorieffs program (hence the name). There are two imp'
        write(*,'(A)', advance='no') 'ortant differences: automatic weighting of the frames using a correlation-based'
        write(*,'(A)', advance='no') 'M-estimator and continuous optimisation of the shift parameters. Input is a text'
        write(*,'(A)', advance='no') 'file with absolute paths to movie files in addition to a few input parameters, s'
        write(*,'(A)', advance='no') 'ome of which deserve a comment. If dose_rate and exp_time are given the individu'
        write(*,'(A)', advance='no') 'al frames will be low-pass filtered accordingly (dose-weighting strategy). If sc'
        write(*,'(A)', advance='no') 'ale is given, the movie will be Fourier cropped according to the down-scaling fa'
        write(*,'(A)', advance='no') 'ctor (for super-resolution movies). If nframesgrp is given the frames will be pr'
        write(*,'(A)', advance='no') 'e-averaged in the given chunk size (Falcon 3 movies). If fromf/tof are given, a'
        write(*,'(A)', advance='no') 'contiguous subset of frames will be averaged without any dose-weighting applied.'
        write(*,'(A)') ''
        stop
    end subroutine print_doc_motion_correct

    subroutine print_doc_motion_correct_ctffind
        write(*,'(A)') 'is a pipelined distributed workflow: motion_correct + ctffind program'
        stop
    end subroutine print_doc_motion_correct_ctffind

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

    subroutine print_doc_pick
        write(*,'(A)') 'is a distributed workflow for template-based particle picking'
        stop
    end subroutine print_doc_pick

    subroutine print_doc_powerspecs
        write(*,'(A)') 'is a program for generating powerspectra from a stack or filetable'
        stop
    end subroutine print_doc_powerspecs

    subroutine print_doc_preprocess
        write(*,'(A)', advance='no') 'is a distributed workflow that executes motion_correct, ctffind and pick in sequ'
        write(*,'(A)') 'ence and in streaming mode as the microscope collects the data'
        stop
    end subroutine print_doc_preprocess

    subroutine print_doc_prime3D_init
        write(*,'(A)', advance='no') 'is a distributed workflow for generating a random initial model for initialisati'
        write(*,'(A)', advance='no') 'on of PRIME3D. If the data set is large (>5000 images), generating a random mode'
        write(*,'(A)', advance='no') 'l can be slow. To speedup, set nran to some smaller number, resulting in nran im'
        write(*,'(A)') 'ages selected randomly for reconstruction'
        stop
    end subroutine print_doc_prime3D_init

    subroutine print_doc_reconstruct3D
        write(*,'(A)', advance='no') 'is a distributed workflow for reconstructing volumes from MRC and SPIDER stacks,'
        write(*,'(A)', advance='no') ' given input orientations and state assignments. The algorithm is based on direc'
        write(*,'(A)', advance='no') 't Fourier inversion with a Kaiser-Bessel (KB) interpolation kernel. This window'
        write(*,'(A)', advance='no') 'function reduces the real-space ripple artifacts associated with direct moving w'
        write(*,'(A)', advance='no') 'indowed-sinc interpolation. The feature sought when implementing this algorithm'
        write(*,'(A)', advance='no') 'was to enable quick, reliable reconstruction from aligned individual particle im'
        write(*,'(A)', advance='no') 'ages. mul is used to scale the origin shifts if down-sampled were used for align'
        write(*,'(A)', advance='no') 'ment and the original images are used for reconstruction. ctf=yes or ctf=flip tu'
        write(*,'(A)', advance='no') 'rns on the Wiener restoration. If the images were phase-flipped set ctf=flip. am'
        write(*,'(A)', advance='no') 'sklp, mw, and edge control the solvent mask: the low-pass limit used to generate'
        write(*,'(A)', advance='no') ' the envelope; the molecular weight of the molecule (protein assumed but it work'
        write(*,'(A)', advance='no') 's reasonably well also for RNA; slight modification of mw might be needed). The'
        write(*,'(A)', advance='no') 'inner parameter controls the radius of the soft-edged mask used to remove the un'
        write(*,'(A)') 'ordered DNA/RNA core of spherical icosahedral viruses'
        stop
    end subroutine print_doc_reconstruct3D

    subroutine print_doc_refine3D
        write(*,'(A)', advance='no') 'is a distributed workflow for ab inito reconstruction/refinement based on probab'
        write(*,'(A)', advance='no') 'ilistic projection matching. PRIME is short for PRobabilistic Initial 3D Model g'
        write(*,'(A)', advance='no') 'eneration for Single-particle cryo-Electron microscopy. There are a daunting num'
        write(*,'(A)', advance='no') 'ber of options in refine3D. If you are processing class averages we recommend th'
        write(*,'(A)', advance='no') 'at you instead use the simple_distr_exec prg=initial_3Dmodel route for executing'
        write(*,'(A)', advance='no') ' refine3D. Automated workflows for single- and multi-particle refinement using r'
        write(*,'(A)') 'efine3D are planned for the next release (3.0)'
        stop
    end subroutine print_doc_refine3D

    subroutine print_doc_scale_stk_parts
        write(*,'(A)') 'is a distributed workflow for scaling partial stacks'
        stop
    end subroutine print_doc_scale_stk_parts

    subroutine print_doc_symsrch
        write(*,'(A)', advance='no') 'is a distributed workflow for searching for the principal symmetry axis of a vol'
        write(*,'(A)', advance='no') 'ume reconstructed without assuming any point-group symmetry. The program takes a'
        write(*,'(A)', advance='no') 's input an asymmetrical 3D reconstruction. The alignment document for all the pa'
        write(*,'(A)', advance='no') 'rticle images that have gone into the 3D reconstruction and the desired point-gr'
        write(*,'(A)', advance='no') 'oup symmetry needs to be inputted. The 3D reconstruction is then projected in 50'
        write(*,'(A)', advance='no') ' (default option) even directions, common lines-based optimisation is used to id'
        write(*,'(A)', advance='no') 'entify the principal symmetry axis, the rotational transformation is applied to'
        write(*,'(A)', advance='no') 'the inputted orientations, and a new alignment document is produced. Input this'
        write(*,'(A)', advance='no') 'document to reconstruct3D together with the images and the point-group symmetry'
        write(*,'(A)', advance='no') 'to generate a symmetrised map. If you are unsure about the point-group, you shou'
        write(*,'(A)', advance='no') 'ld use the compare=yes mode and input the highest conceviable point-group. The p'
        write(*,'(A)') 'rogram then calculates probabilities for all lower groups inclusive.'
        stop
    end subroutine print_doc_symsrch

    subroutine print_doc_tseries_track
        write(*,'(A)') 'is a distributed workflow for particle tracking in time-series data'
        stop
    end subroutine print_doc_tseries_track

    subroutine list_all_simple_programs
        write(*,'(A)') 'center'
        write(*,'(A)') 'cluster_cavgs'
        write(*,'(A)') 'convert'
        write(*,'(A)') 'ctfops'
        write(*,'(A)') 'extract'
        write(*,'(A)') 'filter'
        write(*,'(A)') 'fsc'
        write(*,'(A)') 'info_image'
        write(*,'(A)') 'info_stktab'
        write(*,'(A)') 'make_deftab'
        write(*,'(A)') 'make_oris'
        write(*,'(A)') 'make_pickrefs'
        write(*,'(A)') 'map2ptcls'
        write(*,'(A)') 'mask'
        write(*,'(A)') 'normalize'
        write(*,'(A)') 'orisops'
        write(*,'(A)') 'oristats'
        write(*,'(A)') 'postprocess'
        write(*,'(A)') 'print_cmd_dict'
        write(*,'(A)') 'print_fsc'
        write(*,'(A)') 'print_magic_boxes'
        write(*,'(A)') 'project'
        write(*,'(A)') 'scale'
        write(*,'(A)') 'select'
        write(*,'(A)') 'shift'
        write(*,'(A)') 'simulate_movie'
        write(*,'(A)') 'simulate_noise'
        write(*,'(A)') 'simulate_particles'
        write(*,'(A)') 'simulate_subtomogram'
        write(*,'(A)') 'stack'
        write(*,'(A)') 'stackops'
        write(*,'(A)') 'vizoris'
        write(*,'(A)') 'volops'
        stop
    end subroutine list_all_simple_programs

    subroutine list_all_simple_distr_programs
        write(*,'(A)') 'cluster2D'
        write(*,'(A)') 'cluster2D_stream'
        write(*,'(A)') 'cluster3D'
        write(*,'(A)') 'cluster3D_refine'
        write(*,'(A)') 'comlin_smat'
        write(*,'(A)') 'ctf_estimate'
        write(*,'(A)') 'ctffind'
        write(*,'(A)') 'initial_3Dmodel'
        write(*,'(A)') 'make_cavgs'
        write(*,'(A)') 'motion_correct'
        write(*,'(A)') 'motion_correct_ctffind'
        write(*,'(A)') 'motion_correct_tomo'
        write(*,'(A)') 'pick'
        write(*,'(A)') 'powerspecs'
        write(*,'(A)') 'preprocess'
        write(*,'(A)') 'prime3D_init'
        write(*,'(A)') 'reconstruct3D'
        write(*,'(A)') 'refine3D'
        write(*,'(A)') 'scale_stk_parts'
        write(*,'(A)') 'symsrch'
        write(*,'(A)') 'tseries_track'
        stop
    end subroutine list_all_simple_distr_programs

end module simple_gen_doc
