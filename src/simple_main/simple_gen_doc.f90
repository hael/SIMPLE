module simple_gen_doc
implicit none

contains

    subroutine print_doc_noiseimgs
        write(*,'(A)') ' is a program for generating noise images'
        stop
    end subroutine print_doc_noiseimgs

    subroutine print_doc_simimgs
        write(*,'(A)', advance='no') ' is a program for simulating cryo-EM images. It is not a very sophisticated simu'
        write(*,'(A)', advance='no') 'lator, but it is nevertheless useful for testing purposes. It does not do any mu'
        write(*,'(A)', advance='no') 'lti-slice simulation and it cannot be used for simulating molecules containing h'
        write(*,'(A)', advance='no') 'eavy atoms. It does not even accept a PDB file as an input. Input is a cryo-EM m'
        write(*,'(A)', advance='no') 'ap, which we usually generate from a PDB file using EMANs program pdb2mrc. simim'
        write(*,'(A)', advance='no') 'gs then projects the volume using Fourier interpolation, applies 20% of the tota'
        write(*,'(A)', advance='no') 'l noise to the images (pink noise), Fourier transforms them, and multiplies them'
        write(*,'(A)', advance='no') ' with astigmatic CTF and B-factor. The images are inverse FTed before the remain'
        write(*,'(A)') 'ing 80% of the noise (white noise) is added.'
        stop
    end subroutine print_doc_simimgs

    subroutine print_doc_simmovie
        write(*,'(A)', advance='no') ' is a program for crude simulation of a DDD movie. Input is a set of projection'
        write(*,'(A)', advance='no') 'images to place. Movie frames are then generated related by randomly shifting th'
        write(*,'(A)') 'e base image and applying two different noise sources: shot and detector noise.'
        stop
    end subroutine print_doc_simmovie

    subroutine print_doc_simsubtomo
        write(*,'(A)') ' is a program for crude simulation of a subtomograms.'
        stop
    end subroutine print_doc_simsubtomo

    subroutine print_doc_select_frames
        write(*,'(A)') ' is a program for selecting contiguous segments of frames from DDD movies.'
        stop
    end subroutine print_doc_select_frames

    subroutine print_doc_boxconvs
        write(*,'(A)', advance='no') 'is a program for averaging overlapping boxes across a micrograph in order to che'
        write(*,'(A)') 'ck if gain correction was appropriately done.'
        stop
    end subroutine print_doc_boxconvs

    subroutine print_doc_integrate_movies
        write(*,'(A)') ' is a program for integrating DDD movies.'
        stop
    end subroutine print_doc_integrate_movies

    subroutine print_doc_powerspecs
        write(*,'(A)') ' is a program for generating powerspectra from a stack or filetable.'
        stop
    end subroutine print_doc_powerspecs

    subroutine print_doc_unblur_movies
        write(*,'(A)', advance='no') ' is a program for movie alignment or unblurring. Input is a textfile with absolu'
        write(*,'(A)', advance='no') 'te paths to movie files in addition to a few obvious input parameters. Output is'
        write(*,'(A)') ' (x,y) shift parameters for every frame of the movie.'
        stop
    end subroutine print_doc_unblur_movies

    subroutine print_doc_select
        write(*,'(A)') ' is a program for selecting files based on image correlation matching.'
        stop
    end subroutine print_doc_select

    subroutine print_doc_pick
        write(*,'(A)', advance='no') ' is a program that implements robust template-based particle picking from integr'
        write(*,'(A)', advance='no') 'ated movies (micrographs) generated with unblur. Use make_pick_refs to generate'
        write(*,'(A)') 'the refs templates (from a volume or from class averages)'
        stop
    end subroutine print_doc_pick

    subroutine print_doc_extract
        write(*,'(A)', advance='no') ' is a program that extracts particle images from DDD movies or integrated movies'
        write(*,'(A)', advance='no') '. Boxfiles are assumed to be in EMAN format but we provide a conversion script ('
        write(*,'(A)', advance='no') 'relion2emanbox.pl) for *.star files containing particle coordinates obtained wit'
        write(*,'(A)', advance='no') 'h Relion. The program creates one stack per movie frame as well as a stack of co'
        write(*,'(A)', advance='no') 'rrected framesums. In addition to single-particle image stacks, the program prod'
        write(*,'(A)', advance='no') 'uces a parameter file extract_params.txt that can be used in conjunction with ot'
        write(*,'(A)', advance='no') 'her SIMPLE programs. We obtain CTF parameters with CTFFIND4 using the script (ex'
        write(*,'(A)', advance='no') 'ec_ctffind.pl) but if you have already obtained CTF parameters from CTFFIND4, pl'
        write(*,'(A)') 'ease see section "CTF parameters convention" in the manual.'
        stop
    end subroutine print_doc_extract

    subroutine print_doc_prime2D_init
        write(*,'(A)', advance='no') ' is used  to produce the initial random references for prime2D execution. The ra'
        write(*,'(A)', advance='no') 'ndom clustering and in-plane alignment is printed in the file prime2D_startdoc.t'
        write(*,'(A)', advance='no') 'xt produced by the program. This file is used together with the initial referenc'
        write(*,'(A)') 'es (startcavgs.ext) to execute prime2D.'
        stop
    end subroutine print_doc_prime2D_init

    subroutine print_doc_prime2D
        write(*,'(A)', advance='no') ' is a reference-free 2D alignment/clustering algorithm adopted from the prime3D'
        write(*,'(A)', advance='no') 'probabilistic  ab initio 3D reconstruction algorithm. Do not search the origin s'
        write(*,'(A)', advance='no') 'hifts initially, when the cluster centers are of low quality. If your images are'
        write(*,'(A)', advance='no') ' far off centre, use stackops with option shalgn=yes instead to shiftalign the i'
        write(*,'(A)', advance='no') 'mages beforehand (the algorithm implemented is the same as EMANs cenalignint pro'
        write(*,'(A)') 'gram).'
        stop
    end subroutine print_doc_prime2D

    subroutine print_doc_classrefine
        write(*,'(A)', advance='no') ' is a program for multi-resolution within class refinement with validation. Outp'
        write(*,'(A)') 'ut is the refined class average and a doc with updated parameters.'
        stop
    end subroutine print_doc_classrefine

    subroutine print_doc_cavgassemble
        write(*,'(A)', advance='no') ' is a program that assembles class averages when the clustering program (prime2D'
        write(*,'(A)') ') has been executed in distributed mode.'
        stop
    end subroutine print_doc_cavgassemble

    subroutine print_doc_check2D_conv
        write(*,'(A)', advance='no') ' is a program for checking if a PRIME2D run has converged. The statistics output'
        write(*,'(A)', advance='no') 'ted include (1) the overlap between the distribution of parameters for succesive'
        write(*,'(A)', advance='no') ' runs. (2) The percentage of search space scanned, i.e. how many reference image'
        write(*,'(A)', advance='no') 's are evaluated on average. (3) The average correlation between the images and t'
        write(*,'(A)', advance='no') 'heir corresponding best matching reference section. If convergence to a local op'
        write(*,'(A)', advance='no') 'timum is achieved, the fraction increases. Convergence is achieved if the parame'
        write(*,'(A)', advance='no') 'ter distribution overlap is larger than 0.95 and more than 99% of the reference'
        write(*,'(A)') 'sections need to be searched to find an improving solution.'
        stop
    end subroutine print_doc_check2D_conv

    subroutine print_doc_rank_cavgs
        write(*,'(A)', advance='no') ' is a program for ranking class averages by decreasing population, given the sta'
        write(*,'(A)', advance='no') 'ck of class averages (stk argument) and the 2D orientations document (oritab) ge'
        write(*,'(A)') 'nerated by prime2D.'
        stop
    end subroutine print_doc_rank_cavgs

    subroutine print_doc_resrange
        write(*,'(A)', advance='no') ' is a program for estimating the resolution range used in the heuristic resoluti'
        write(*,'(A)', advance='no') 'on-stepping scheme in the PRIME3D initial model production procedure. The initia'
        write(*,'(A)', advance='no') 'l low-pass limit is set so that each image receives ten nonzero orientation weig'
        write(*,'(A)', advance='no') 'hts. When quasi-convergence has been reached, the limit is updated one Fourier i'
        write(*,'(A)', advance='no') 'ndex at the time, until PRIME reaches the condition where six nonzero orientatio'
        write(*,'(A)', advance='no') 'n weights are assigned to each image. FSC-based filtering is unfortunately not p'
        write(*,'(A)', advance='no') 'ossible to do in the ab initio 3D reconstruction step, because when the orientat'
        write(*,'(A)', advance='no') 'ions are mostly random, the FSC overestimates the resolution. This program is us'
        write(*,'(A)', advance='no') 'ed internally when executing PRIME in distributed mode. We advise you to check t'
        write(*,'(A)', advance='no') 'he starting and stopping low-pass limits before executing PRIME3D using this pro'
        write(*,'(A)', advance='no') 'gram. The resolution range estimate depends on the molecular diameter, which is'
        write(*,'(A)', advance='no') 'estimated based on the box size. If you want to override this estimate, set mold'
        write(*,'(A)', advance='no') 'iam to the desired value (in \AA{). This may be necessary if your images have a'
        write(*,'(A)', advance='no') 'lot of background padding. However, for starting model generation it is probably'
        write(*,'(A)', advance='no') ' better to clip the images snugly around the particle, because smaller images eq'
        write(*,'(A)') 'ual less computation.'
        stop
    end subroutine print_doc_resrange

    subroutine print_doc_npeaks
        write(*,'(A)', advance='no') ' is a program for checking the number of nonzero orientation weights (number of'
        write(*,'(A)') 'correlation peaks included in the weighted reconstruction).'
        stop
    end subroutine print_doc_npeaks

    subroutine print_doc_nspace
        write(*,'(A)', advance='no') ' is a program for calculating the expected resolution obtainable with different'
        write(*,'(A)', advance='no') 'values of nspace (number of discrete projection directions used for discrete sea'
        write(*,'(A)') 'rch).'
        stop
    end subroutine print_doc_nspace

    subroutine print_doc_shellweight3D
        write(*,'(A)', advance='no') ' is a program for calculating the shell-by-shell resolution weights in a global'
        write(*,'(A)', advance='no') 'sense, so that particles that do contribute with higher resolution information ('
        write(*,'(A)') 'as measure by the FRC) are given the appropriate weight.'
        stop
    end subroutine print_doc_shellweight3D

    subroutine print_doc_prime3D_init
        write(*,'(A)', advance='no') ' is a program for generating a random initial model for initialisation of PRIME3'
        write(*,'(A)', advance='no') 'D. If the data set is large (>5000 images), generating a random model can be slo'
        write(*,'(A)', advance='no') 'w. To speedup, set nran to some smaller number, resulting in nran images selecte'
        write(*,'(A)') 'd randomly for reconstruction.'
        stop
    end subroutine print_doc_prime3D_init

    subroutine print_doc_multiptcl_init
        write(*,'(A)', advance='no') ' is a program for generating random initial models for initialisation of PRIME3D'
        write(*,'(A)') ' when run in multiparticle mode.'
        stop
    end subroutine print_doc_multiptcl_init

    subroutine print_doc_prime3D
        write(*,'(A)', advance='no') ' is an ab inito reconstruction/refinement program based on probabilistic project'
        write(*,'(A)', advance='no') 'ion matching. PRIME is short for PRobabilistic Initial 3D Model generation for S'
        write(*,'(A)', advance='no') 'ingle- particle cryo-Electron microscopy. Do not search the origin shifts initia'
        write(*,'(A)', advance='no') 'lly, when the model is of very low quality. If your images are far off centre, u'
        write(*,'(A)', advance='no') 'se stackops with option shalgn=yes instead to shiftalign the images beforehand ('
        write(*,'(A)', advance='no') 'the algorithm implemented is the same as EMANs cenalignint program). We recommen'
        write(*,'(A)', advance='no') 'd running the first round of PRIME with the default dynamic resolution stepping'
        write(*,'(A)', advance='no') 'dynlp=yes. The dynlp option implements a heuristic resolution weighting/update s'
        write(*,'(A)', advance='no') 'cheme. The initial low-pass limit is set so that each image receives ten nonzero'
        write(*,'(A)', advance='no') ' orientation weights. When quasi-convergence has been reached, the limit is upda'
        write(*,'(A)', advance='no') 'ted one Fourier index at the time until PRIME reaches the condition where six no'
        write(*,'(A)', advance='no') 'nzero orientation weights are assigned to each image. FSC-based filtering is unf'
        write(*,'(A)', advance='no') 'ortunately not possible to do in the ab initio reconstruction step, because when'
        write(*,'(A)', advance='no') ' the orientations are mostly random, the FSC overestimates the resolution. Once'
        write(*,'(A)', advance='no') 'the initial model has converged, we recommend start searching the shifts (by set'
        write(*,'(A)', advance='no') 'ting trs to some nonzero value) and applying the FSC for resolution- weighting ('
        write(*,'(A)', advance='no') 'by setting eo=yes). In order to be able to use Wiener restoration, give the ctf'
        write(*,'(A)', advance='no') 'flag on the command line to indicate what has been done to the images. You then'
        write(*,'(A)', advance='no') 'also need to input CTF parameters, for example via deftab=defocus_values.txt. Re'
        write(*,'(A)', advance='no') 'member that the defocus values should be given in microns and the astigmatism an'
        write(*,'(A)', advance='no') 'gle in degrees (one row of the file defocus_values.txt may look like: dfx=3.5 df'
        write(*,'(A)', advance='no') 'y=3.3 angast=20.0). Note that we do not assume any point-group symmetry in the i'
        write(*,'(A)', advance='no') 'nitial runs. However, the symsrch program can be used to align the 3D reconstruc'
        write(*,'(A)', advance='no') 'tion to its symmetry axis so that future searches can be restricted to the asymm'
        write(*,'(A)', advance='no') 'etric unit. Less commonly used and less obvious input parameters are nspace, whi'
        write(*,'(A)', advance='no') 'ch  controls the number of reference projections, amsklp, which controls the low'
        write(*,'(A)', advance='no') '-pass limit used in the automask routine, maxits, which controls the maximum num'
        write(*,'(A)', advance='no') 'ber of iterations executed, pgrp, which controls the point- group symmetry, assu'
        write(*,'(A)', advance='no') 'ming that the starting volume is aligned to its principal symmetry axis, edge, w'
        write(*,'(A)') 'hich controls the size of the softening edge in the automask routine.'
        stop
    end subroutine print_doc_prime3D

    subroutine print_doc_cont3D
        write(*,'(A)') ' is a continuous refinement code under development.'
        stop
    end subroutine print_doc_cont3D

    subroutine print_doc_check3D_conv
        write(*,'(A)', advance='no') ' is a program for checking if a PRIME3D run has converged. The statistics output'
        write(*,'(A)', advance='no') 'ted include (1) angle of feasible region, which is proportional to the angular r'
        write(*,'(A)', advance='no') 'esolution of the set of discrete projection directions being searched. (2) The a'
        write(*,'(A)', advance='no') 'verage angular distance between orientations in the present and previous iterati'
        write(*,'(A)', advance='no') 'on. In the early iterations, the distance is large because a diverse set of orie'
        write(*,'(A)', advance='no') 'ntations is explored. If convergence to a local optimum is achieved, the distanc'
        write(*,'(A)', advance='no') 'e decreases. (3) The percentage of search space scanned, i.e. how many reference'
        write(*,'(A)', advance='no') ' images are evaluated on average. (4) The average correlation between the images'
        write(*,'(A)', advance='no') ' and their corresponding best matching reference sections. (5) The average stand'
        write(*,'(A)', advance='no') 'ard deviation of the Euler angles. Convergence is achieved if the angular distan'
        write(*,'(A)', advance='no') 'ce between the orientations in successive iterations falls significantly below t'
        write(*,'(A)', advance='no') 'he angular resolution of the search space and more than 99% of the reference sec'
        write(*,'(A)') 'tions need to be matched on average.'
        stop
    end subroutine print_doc_check3D_conv

    subroutine print_doc_comlin_smat
        write(*,'(A)', advance='no') ' is a program for creating a similarity matrix based on common line correlation.'
        write(*,'(A)', advance='no') ' The idea being that it should be possible to cluster images based on their 3D s'
        write(*,'(A)', advance='no') 'imilarity witout having a 3D model by only operating on class averages and find'
        write(*,'(A)') 'averages that fit well together in 3D.'
        stop
    end subroutine print_doc_comlin_smat

    subroutine print_doc_symsrch
        write(*,'(A)', advance='no') ' is a program for searching for the principal symmetry axis of a volume reconstr'
        write(*,'(A)', advance='no') 'ucted without assuming any point-group symmetry or assessing the degree of symme'
        write(*,'(A)', advance='no') 'try of class averages or individual particles of higher pointgroups (dn,t,o,i).'
        write(*,'(A)', advance='no') 'The program takes as input an asymmetrical reconstruction or stack of class aver'
        write(*,'(A)', advance='no') 'ages/individual particles. For volumes, the alignment document for all the parti'
        write(*,'(A)', advance='no') 'cle images that have gone into the 3D reconstruction and the desired point-group'
        write(*,'(A)', advance='no') ' symmetry needs to be inputted. The 3D reconstruction is then projected in 20 (d'
        write(*,'(A)', advance='no') 'efault option) even directions, common lines-based optimisation is used to ident'
        write(*,'(A)', advance='no') 'ify the principal symmetry axis, the rotational transformation is applied to the'
        write(*,'(A)', advance='no') ' inputted orientations, and a new alignment document is produced. Input this doc'
        write(*,'(A)', advance='no') 'ument to recvol or eo_recvol together with the images and the point-group symmet'
        write(*,'(A)', advance='no') 'ry to generate a symmetrised map. If you are unsure about the point-group, you s'
        write(*,'(A)', advance='no') 'hould use the compare=yes mode and input the highest conceviable point-group. Th'
        write(*,'(A)', advance='no') 'e program then calculates probabilities for all lower groups inclusive. The clas'
        write(*,'(A)', advance='no') 's average/particle option operates in an equivalent fashion but with individual'
        write(*,'(A)', advance='no') 'images. The output is then a per-image correlation value that informs about how'
        write(*,'(A)', advance='no') 'well the image conforms to to inputted point-group. The state parameter allows y'
        write(*,'(A)') 'ou to apply symmetry for the given state.'
        stop
    end subroutine print_doc_symsrch

    subroutine print_doc_mask
        write(*,'(A)', advance='no') ' is a program for masking images and volumes. If you want to mask your images wi'
        write(*,'(A)') 'th a spherical mask with a soft falloff, set msk to the radius in pixels.'
        stop
    end subroutine print_doc_mask

    subroutine print_doc_automask2D
        write(*,'(A)', advance='no') ' is a program for solvent flattening of class averages (MRC or SPIDER). The algo'
        write(*,'(A)', advance='no') 'rithm for background removal is based on low-pass filtering and binarization. Fi'
        write(*,'(A)', advance='no') 'rst, the class averages are low-pass filtered to amsklp. Binary representatives'
        write(*,'(A)', advance='no') 'are then generated by assigning foreground pixels using sortmeans. A cosine func'
        write(*,'(A)', advance='no') 'tion softens the edge of the binary mask before it is  multiplied with the unmas'
        write(*,'(A)') 'ked input averages to accomplish flattening.'
        stop
    end subroutine print_doc_automask2D

    subroutine print_doc_automask3D
        write(*,'(A)', advance='no') ' is a program for solvent flattening of a volume (MRC or SPIDER). The algorithm'
        write(*,'(A)', advance='no') 'for background removal is based on low-pass filtering and binarization. First, t'
        write(*,'(A)', advance='no') 'he volume is low-pass filtered to amsklp. A binary volume is then generated by a'
        write(*,'(A)', advance='no') 'ssigning foreground pixels (=1) based on the volume calculated from the molecula'
        write(*,'(A)', advance='no') 'r weight. A cosine function softens the edge of the binary mask before it is  mu'
        write(*,'(A)') 'ltiplied with the unmasked input to generate the flattened map.'
        stop
    end subroutine print_doc_automask3D

    subroutine print_doc_eo_recvol
        write(*,'(A)', advance='no') ' is a program for reconstructing volumes from MRC or SPIDER stacks, given input'
        write(*,'(A)', advance='no') 'orientations and state assignments (obtained by program prime3D). The algorithm'
        write(*,'(A)', advance='no') 'is based on direct Fourier inversion with a Kaiser-Bessel (KB) interpolation ker'
        write(*,'(A)', advance='no') 'nel. This window function reduces the real-space ripple artefacts associated wit'
        write(*,'(A)', advance='no') 'h direct moving windowed-sinc interpolation. The feature sought when implementin'
        write(*,'(A)', advance='no') 'g this algorithm was to enable quick, reliable reconstruction from aligned indiv'
        write(*,'(A)', advance='no') 'idual particle images. The even and odd pairs are automatically reconstructed, t'
        write(*,'(A)', advance='no') 'he FSC calculated, and the Wiener filter formalism used for image restoration (C'
        write(*,'(A)', advance='no') 'TF correction). mul is used to scale the origin shifts if down-sampled images we'
        write(*,'(A)', advance='no') 're used for alignment and the original images are used for reconstruction. ctf,'
        write(*,'(A)', advance='no') 'kv, fraca, cs and deftab are used to communicate CTF information to the program.'
        write(*,'(A)', advance='no') ' ctf=yes, ctf=flip or ctf=mul turns on the Wiener restoration. If you input CTF'
        write(*,'(A)', advance='no') 'info to the program, please ensure that the correct kV, Cs and fraca (fraction o'
        write(*,'(A)', advance='no') 'f amplitude contrast) parameters are inputted as well. If the images were pre-mu'
        write(*,'(A)', advance='no') 'ltiplied with the CTF, set ctf=mul or if the images were phase-flipped set ctf=f'
        write(*,'(A)', advance='no') 'lip. amsklp and mw parameters control the solvent mask: the low-pass limit used'
        write(*,'(A)', advance='no') 'to generate the envelope; the molecular weight of the molecule (protein assumed'
        write(*,'(A)', advance='no') 'but it works reasonably well also for RNA; slight modification of mw might be ne'
        write(*,'(A)', advance='no') 'eded). The inner parameter controls the radius of the soft-edged mask used to re'
        write(*,'(A)') 'move the unordered DNA/RNA core of spherical icosahedral viruses.'
        stop
    end subroutine print_doc_eo_recvol

    subroutine print_doc_eo_volassemble
        write(*,'(A)', advance='no') ' is a program that assembles volume(s) when the reconstruction program (eo_recvo'
        write(*,'(A)', advance='no') 'l) has been executed in distributed mode using distr_simple.pl. inner applies a'
        write(*,'(A)', advance='no') 'soft-edged inner mask. An inner mask is used for icosahedral virus reconstructio'
        write(*,'(A)', advance='no') 'n, because the DNA or RNA core is often unordered and  if not removed it may neg'
        write(*,'(A)', advance='no') 'atively impact the alignment. The width parameter controls the fall-off of the e'
        write(*,'(A)') 'dge of the inner mask.'
        stop
    end subroutine print_doc_eo_volassemble

    subroutine print_doc_recvol
        write(*,'(A)', advance='no') ' is a program for reconstructing volumes from MRC and SPIDER stacks, given input'
        write(*,'(A)', advance='no') ' orientations and state assignments. The algorithm is based on direct Fourier in'
        write(*,'(A)', advance='no') 'version with a Kaiser-Bessel (KB) interpolation kernel. This window function red'
        write(*,'(A)', advance='no') 'uces the real-space ripple artifacts associated with direct moving windowed-sinc'
        write(*,'(A)', advance='no') ' interpolation. The feature sought when implementing this algorithm was to enabl'
        write(*,'(A)', advance='no') 'e quick, reliable reconstruction from aligned individual particle images. mul is'
        write(*,'(A)', advance='no') ' used to scale the origin shifts if down-sampled were used for alignment and the'
        write(*,'(A)', advance='no') ' original images are used for reconstruction. ctf, kv, fraca, cs and deftab are'
        write(*,'(A)', advance='no') 'used to communicate CTF information to the program. ctf=yes or ctf=flip turns on'
        write(*,'(A)', advance='no') ' the Wiener restoration. If the images were phase-flipped set ctf=flip. amsklp,'
        write(*,'(A)', advance='no') 'mw, and edge control the solvent mask: the low-pass limit used to generate the e'
        write(*,'(A)', advance='no') 'nvelope; the molecular weight of the molecule (protein assumed but it works reas'
        write(*,'(A)', advance='no') 'onably well also for RNA; slight modification of mw might be needed). The inner'
        write(*,'(A)', advance='no') 'parameter controls the radius of the soft-edged mask used to remove the unordere'
        write(*,'(A)', advance='no') 'd DNA/RNA core of spherical icosahedral viruses. The even and odd parameters all'
        write(*,'(A)') 'ow you to reconstruct either the even or the odd pair.'
        stop
    end subroutine print_doc_recvol

    subroutine print_doc_volassemble
        write(*,'(A)', advance='no') ' is a program that assembles volume(s) when the reconstruction program (recvol)'
        write(*,'(A)', advance='no') 'has been executed in distributed mode. odd is used to assemble the odd reconstru'
        write(*,'(A)', advance='no') 'ction, even is used to assemble the even reconstruction, eo is used to assemble'
        write(*,'(A)', advance='no') 'both the even and the odd reconstruction and state is used to assemble the input'
        write(*,'(A)', advance='no') 'ted state. Normally, you do not fiddle with these parameters. They are used inte'
        write(*,'(A)') 'rnally.'
        stop
    end subroutine print_doc_volassemble

    subroutine print_doc_check_box
        write(*,'(A)', advance='no') ' is a program for checking the image dimensions of MRC and SPIDER stacks and vol'
        write(*,'(A)') 'umes.'
        stop
    end subroutine print_doc_check_box

    subroutine print_doc_check_nptcls
        write(*,'(A)') ' is a program for checking the number of images in MRC and SPIDER stacks.'
        stop
    end subroutine print_doc_check_nptcls

    subroutine print_doc_iminfo
        write(*,'(A)', advance='no') ' is a program for printing header information in MRC and SPIDER stacks and volum'
        write(*,'(A)') 'es.'
        stop
    end subroutine print_doc_iminfo

    subroutine print_doc_cenvol
        write(*,'(A)', advance='no') ' is a program for centering a volume and mapping the shift parameters back to th'
        write(*,'(A)', advance='no') 'e particle images, Often, when class averages are used for 3D processing and the'
        write(*,'(A)', advance='no') ' paramaters are mapped back to the particles, the reconstructed volume is off-ce'
        write(*,'(A)', advance='no') 'ntre. This program is useful for making sure that the mask does not cut off-cent'
        write(*,'(A)') 're volumes.'
        stop
    end subroutine print_doc_cenvol

    subroutine print_doc_projvol
        write(*,'(A)', advance='no') ' is a program for projecting a volume using interpolation in Fourier space. Inpu'
        write(*,'(A)', advance='no') 't is a SPIDER or MRC volume. Output is a stack of projection images of the same'
        write(*,'(A)', advance='no') 'format as the inputted volume. Projections are generated by extraction of centra'
        write(*,'(A)', advance='no') 'l sections from the Fourier volume and back transformation of the 2D FTs. nspace'
        write(*,'(A)', advance='no') ' controls the number of projection images generated with quasi-even projection d'
        write(*,'(A)', advance='no') 'irections.The oritab parameter allows you to input the orientations that you wis'
        write(*,'(A)', advance='no') 'h to have your volume projected in. If rnd=yes, random rather than quasi-even pr'
        write(*,'(A)', advance='no') 'ojections are generated, trs then controls the halfwidth of the random origin sh'
        write(*,'(A)', advance='no') 'ift. Less commonly used parameters are pgrp, which controls the point-group symm'
        write(*,'(A)', advance='no') 'etry c (rotational), d (dihedral), t (tetrahedral), o (octahedral) or i (icosahe'
        write(*,'(A)', advance='no') 'dral). The point-group symmetry is used to restrict the set of projections to wi'
        write(*,'(A)', advance='no') 'thin the asymmetric unit. ctf=yes allows you to apply CTF to the images, using c'
        write(*,'(A)', advance='no') 'onstant defocus and no astigmatism. If you want to do this you need to define th'
        write(*,'(A)', advance='no') 'e parameters kv, fraca, cs, defocus and bfac. neg inverts the contrast of the pr'
        write(*,'(A)', advance='no') 'ojections. mirr=yes mirrors the projection by modifying the Euler angles. If mir'
        write(*,'(A)', advance='no') 'r=x or mirr=y the projection is physically mirrored after it has been generated.'
        write(*,'(A)') ''
        stop
    end subroutine print_doc_projvol

    subroutine print_doc_volaverager
        write(*,'(A)') ' is a program for averaging volumes according to state label in oritab.'
        stop
    end subroutine print_doc_volaverager

    subroutine print_doc_volops
        write(*,'(A)', advance='no') ' provides standard single-particle image processing routines that are applied to'
        write(*,'(A)') ' MRC or SPIDER volumes.'
        stop
    end subroutine print_doc_volops

    subroutine print_doc_volume_smat
        write(*,'(A)', advance='no') ' is a program for creating a similarity matrix based on volume2volume correlatio'
        write(*,'(A)') 'n.'
        stop
    end subroutine print_doc_volume_smat

    subroutine print_doc_binarise
        write(*,'(A)') ' is a program for binarisation of stacks and volumes.'
        stop
    end subroutine print_doc_binarise

    subroutine print_doc_convert
        write(*,'(A)') ' is a program for converting between SPIDER and MRC formats.'
        stop
    end subroutine print_doc_convert

    subroutine print_doc_corrcompare
        write(*,'(A)', advance='no') ' is a program for comparing stacked images using real-space and Fourier-based ap'
        write(*,'(A)') 'proaches'
        stop
    end subroutine print_doc_corrcompare

    subroutine print_doc_ctfops
        write(*,'(A)') ' is a program for applying CTF to stacked images.'
        stop
    end subroutine print_doc_ctfops

    subroutine print_doc_filter
        write(*,'(A)') ' is a program for filtering stacks/volumes.'
        stop
    end subroutine print_doc_filter

    subroutine print_doc_image_smat
        write(*,'(A)', advance='no') ' is a program for creating a similarity matrix based on common line correlation.'
        write(*,'(A)', advance='no') ' The idea being that it should be possible to cluster images based on their 3D s'
        write(*,'(A)', advance='no') 'imilarity witout having a 3D model by only operating on class averages and find'
        write(*,'(A)') 'averages that fit well together in 3D.'
        stop
    end subroutine print_doc_image_smat

    subroutine print_doc_norm
        write(*,'(A)', advance='no') ' is a program for normalization of MRC or SPIDER stacks and volumes. If you want'
        write(*,'(A)', advance='no') ' to normalise your images inputted with stk, set norm=yes. hfun (e.g. hfun=sigm)'
        write(*,'(A)', advance='no') ' controls the normalisation function. If you want to perform noise normalisation'
        write(*,'(A)', advance='no') ' of the images set noise_norm=yes given a mask radius msk (pixels). If you want'
        write(*,'(A)', advance='no') 'to normalise your images or volume (vol1) with respect to their power spectrum s'
        write(*,'(A)') 'et shell_norm=yes.'
        stop
    end subroutine print_doc_norm

    subroutine print_doc_scale
        write(*,'(A)', advance='no') ' is a program that provides re-scaling and clipping routines for MRC or SPIDER s'
        write(*,'(A)') 'tacks and volumes.'
        stop
    end subroutine print_doc_scale

    subroutine print_doc_stack
        write(*,'(A)') 'is a program for stacking individual images or multiple stacks into one.'
        stop
    end subroutine print_doc_stack

    subroutine print_doc_stackops
        write(*,'(A)', advance='no') ' is a program that provides standard single-particle image processing routines t'
        write(*,'(A)', advance='no') 'hat are applied to MRC or SPIDER stacks. If you want to extract a particular sta'
        write(*,'(A)', advance='no') 'te, give an alignment document (oritab) and set state to the state that you want'
        write(*,'(A)', advance='no') ' to extract. If you want to select the fraction of best particles (according to'
        write(*,'(A)', advance='no') 'the goal function), input an alignment doc (oritab) and set frac. You can combin'
        write(*,'(A)', advance='no') 'e the state and frac options. If you want to apply noise to images, give the des'
        write(*,'(A)', advance='no') 'ired signal-to-noise ratio via snr. If you want to calculate the autocorrelation'
        write(*,'(A)', advance='no') ' function of your images set acf=yes. If you want to extract a contiguous subset'
        write(*,'(A)', advance='no') ' of particle images from the stack, set fromp and top. If you want to fish out a'
        write(*,'(A)', advance='no') ' number of particle images from your stack at random, set nran to some nonzero i'
        write(*,'(A)', advance='no') 'nteger number less than nptcls. With avg=yes the global average of the inputted'
        write(*,'(A)', advance='no') 'stack is calculated. If you define frameavg to some integer number larger than o'
        write(*,'(A)', advance='no') 'ne averages with chunk sizes of frameavg are produced, which may be useful for a'
        write(*,'(A)', advance='no') 'nalysis of dose-fractionated image series. neg inverts the contrast of the image'
        write(*,'(A)') 's.'
        stop
    end subroutine print_doc_stackops

    subroutine print_doc_tseries_split
        write(*,'(A)') ' is a program for splitting time series.'
        stop
    end subroutine print_doc_tseries_split

    subroutine print_doc_cluster_smat
        write(*,'(A)', advance='no') ' is a program for clustering a similarity matrix and use an combined cluster val'
        write(*,'(A)', advance='no') 'idation index to assess the quality of the clustering based on the number of clu'
        write(*,'(A)') 'sters.'
        stop
    end subroutine print_doc_cluster_smat

    subroutine print_doc_find_nnimgs
        write(*,'(A)', advance='no') ' is a program for cidentifying the nnn nearest neighbor images for each image in'
        write(*,'(A)') ' the inputted stack.'
        stop
    end subroutine print_doc_find_nnimgs

    subroutine print_doc_masscen
        write(*,'(A)') ' is a program for centering images acccording to their centre of mass'
        stop
    end subroutine print_doc_masscen

    subroutine print_doc_print_cmd_dict
        write(*,'(A)') 'is a program for printing the command line key dictonary.'
        stop
    end subroutine print_doc_print_cmd_dict

    subroutine print_doc_print_dose_weights
        write(*,'(A)') ' is a program for printing the dose weights applied to individual frames.'
        stop
    end subroutine print_doc_print_dose_weights

    subroutine print_doc_print_fsc
        write(*,'(A)') ' is a program for printing the binary FSC files produced by PRIME3D.'
        stop
    end subroutine print_doc_print_fsc

    subroutine print_doc_res
        write(*,'(A)', advance='no') ' is a program for checking the low-pass resolution limit for a given Fourier ind'
        write(*,'(A)') 'ex.'
        stop
    end subroutine print_doc_res

    subroutine print_doc_shift
        write(*,'(A)') ' is a program for shifting a stack according to shifts in oritab'
        stop
    end subroutine print_doc_shift

    subroutine print_doc_cluster_oris
        write(*,'(A)') ' is a program for clustering orientations based on geodesic distance'
        stop
    end subroutine print_doc_cluster_oris

    subroutine print_doc_makeoris
        write(*,'(A)', advance='no') ' is a program for analyzing SIMPLE orientation/parameter files (text files conta'
        write(*,'(A)', advance='no') 'ining input parameters and/or parameters estimated by prime2D or prime3D). The p'
        write(*,'(A)', advance='no') 'rogram generates random Euler angles e1.in.[0,360], e2.in.[0,180], and e3.in.[0,'
        write(*,'(A)', advance='no') '360] and random origin shifts x.in.[-trs,yrs] and y.in.[-trs,yrs]. If ndiscrete'
        write(*,'(A)', advance='no') 'is set to an integer number > 0, the orientations produced are randomly sampled'
        write(*,'(A)', advance='no') 'from the set of ndiscrete quasi-even projection directions, and the in-plane par'
        write(*,'(A)', advance='no') 'ameters are assigned randomly. If even=yes, then all nptcls orientations are ass'
        write(*,'(A)', advance='no') 'igned quasi-even projection directions,and random in-plane parameters. If nstate'
        write(*,'(A)', advance='no') 's is set to some integer number > 0, then states are assigned randomly .in.[1,ns'
        write(*,'(A)', advance='no') 'tates]. If zero=yes in this mode of execution, the projection directions are zer'
        write(*,'(A)', advance='no') 'oed and only the in-plane parameters are kept intact. If errify=yes and astigerr'
        write(*,'(A)', advance='no') ' is defined, then uniform random astigmatism errors are introduced .in.[-astiger'
        write(*,'(A)') 'r,astigerr].'
        stop
    end subroutine print_doc_makeoris

    subroutine print_doc_map2ptcls
        write(*,'(A)', advance='no') ' is a program for mapping parameters that have been obtained using class average'
        write(*,'(A)') 's to the individual particle images.'
        stop
    end subroutine print_doc_map2ptcls

    subroutine print_doc_orisops
        write(*,'(A)', advance='no') ' is a program for analyzing SIMPLE orientation/parameter files (text files conta'
        write(*,'(A)', advance='no') 'ining input parameters and/or parameters estimated by prime2D or prime3D). If on'
        write(*,'(A)', advance='no') 'ly oritab is inputted, there are a few options available. If errify=yes, then th'
        write(*,'(A)', advance='no') 'e program introduces uniform random angular errors .in.[-angerr,angerr], and uni'
        write(*,'(A)', advance='no') 'form origin shift errors .in.[-sherr,sherr], and uniform random defocus errors .'
        write(*,'(A)', advance='no') 'in.[-dferr,dferr]. If nstates > 1 then random states are assigned .in.[1,nstates'
        write(*,'(A)', advance='no') ']. If mirr=2d, then the Euler angles in oritab are mirrored according to the rel'
        write(*,'(A)', advance='no') 'ation e1=e1, e2=180.+e2, e3=-e3. If mirr=3d, then the Euler angles in oritab are'
        write(*,'(A)', advance='no') ' mirrored according to the relation R=M(M*R), where R is the rotation matrix cal'
        write(*,'(A)', advance='no') 'culated from the Euler angle triplet and M is a 3D reflection matrix (like a uni'
        write(*,'(A)', advance='no') 't matrix but with the 3,3-component sign swapped). If e1, e2, or e3 is inputted,'
        write(*,'(A)', advance='no') ' the orientations in oritab are rotated correspondingly. If you input state as w'
        write(*,'(A)', advance='no') 'ell, you rotate only the orientations assigned to state state. If mul is defined'
        write(*,'(A)', advance='no') ', you multiply the origin shifts with mul. If zero=yes, then the shifts are zero'
        write(*,'(A)', advance='no') 'ed. If none of the above described parameter are defined, and oritab is still de'
        write(*,'(A)', advance='no') 'fined, the program projects the 3D orientation into the xy-plane and plots the r'
        write(*,'(A)') 'esulting vector (this is useful for checking orientation coverage).'
        stop
    end subroutine print_doc_orisops

    subroutine print_doc_oristats
        write(*,'(A)', advance='no') ' is a program for analyzing SIMPLE orientation/parameter files (text files conta'
        write(*,'(A)', advance='no') 'ining input parameters and/or parameters estimated by prime2D or prime3D). If tw'
        write(*,'(A)', advance='no') 'o orientation tables (oritab and oritab2) are inputted, the program provides sta'
        write(*,'(A)', advance='no') 'tistics of the distances between the orientations in the two documents. These st'
        write(*,'(A)', advance='no') 'atistics include the sum of angular distances between the orientations, the aver'
        write(*,'(A)', advance='no') 'age angular distance between the orientations, the standard deviation of angular'
        write(*,'(A)') ' distances, the minimum angular distance, and the maximum angular distance.'
        stop
    end subroutine print_doc_oristats

    subroutine print_doc_rotmats2oris
        write(*,'(A)', advance='no') ' convertes a text file (9 records per line) describing rotation matrices into a'
        write(*,'(A)') 'SIMPLE oritab'
        stop
    end subroutine print_doc_rotmats2oris

    subroutine print_doc_merge_algndocs
        write(*,'(A)', advance='no') ' is a program for merging alignment documents SIMPLE was run when run in distrib'
        write(*,'(A)') 'uted mode.'
        stop
    end subroutine print_doc_merge_algndocs

    subroutine print_doc_merge_nnmat
        write(*,'(A)', advance='no') ' is a program for merging partial nearest neighbour matrices calculated in distr'
        write(*,'(A)') 'ibuted mode'
        stop
    end subroutine print_doc_merge_nnmat

    subroutine print_doc_merge_shellweights
        write(*,'(A)', advance='no') ' is a program for merging partial shellweight matrices calculated in distributed'
        write(*,'(A)') ' mode'
        stop
    end subroutine print_doc_merge_shellweights

    subroutine print_doc_merge_similarities
        write(*,'(A)', advance='no') ' is a program for merging similarities calculated between pairs of objects into'
        write(*,'(A)') 'a similarity matrix that can be inoutted to cluster_smat.'
        stop
    end subroutine print_doc_merge_similarities

    subroutine print_doc_split_pairs
        write(*,'(A)', advance='no') ' is a program for splitting calculations between pairs of objects into balanced'
        write(*,'(A)') 'partitions.'
        stop
    end subroutine print_doc_split_pairs

    subroutine print_doc_split
        write(*,'(A)', advance='no') ' is a program for splitting of image stacks into partitions for parallel executi'
        write(*,'(A)') 'on. This is done to reduce I/O latency.'
        stop
    end subroutine print_doc_split

    subroutine list_all_simple_programs
        write(*,'(A)') 'simple_automask2D'
        write(*,'(A)') 'simple_automask3D'
        write(*,'(A)') 'simple_binarise'
        write(*,'(A)') 'simple_boxconvs'
        write(*,'(A)') 'simple_cavgassemble'
        write(*,'(A)') 'simple_cenvol'
        write(*,'(A)') 'simple_check2D_conv'
        write(*,'(A)') 'simple_check3D_conv'
        write(*,'(A)') 'simple_check_box'
        write(*,'(A)') 'simple_check_nptcls'
        write(*,'(A)') 'simple_classrefine'
        write(*,'(A)') 'simple_cluster_oris'
        write(*,'(A)') 'simple_cluster_smat'
        write(*,'(A)') 'simple_comlin_smat'
        write(*,'(A)') 'simple_cont3D'
        write(*,'(A)') 'simple_convert'
        write(*,'(A)') 'simple_corrcompare'
        write(*,'(A)') 'simple_ctfops'
        write(*,'(A)') 'simple_eo_recvol'
        write(*,'(A)') 'simple_eo_volassemble'
        write(*,'(A)') 'simple_extract'
        write(*,'(A)') 'simple_filter'
        write(*,'(A)') 'simple_find_nnimgs'
        write(*,'(A)') 'simple_image_smat'
        write(*,'(A)') 'simple_iminfo'
        write(*,'(A)') 'simple_integrate_movies'
        write(*,'(A)') 'simple_makeoris'
        write(*,'(A)') 'simple_map2ptcls'
        write(*,'(A)') 'simple_mask'
        write(*,'(A)') 'simple_masscen'
        write(*,'(A)') 'simple_merge_algndocs'
        write(*,'(A)') 'simple_merge_nnmat'
        write(*,'(A)') 'simple_merge_shellweights'
        write(*,'(A)') 'simple_merge_similarities'
        write(*,'(A)') 'simple_multiptcl_init'
        write(*,'(A)') 'simple_noiseimgs'
        write(*,'(A)') 'simple_norm'
        write(*,'(A)') 'simple_npeaks'
        write(*,'(A)') 'simple_nspace'
        write(*,'(A)') 'simple_orisops'
        write(*,'(A)') 'simple_oristats'
        write(*,'(A)') 'simple_pick'
        write(*,'(A)') 'simple_postproc_vol'
        write(*,'(A)') 'simple_powerspecs'
        write(*,'(A)') 'simple_prime2D'
        write(*,'(A)') 'simple_prime2D_init'
        write(*,'(A)') 'simple_prime3D'
        write(*,'(A)') 'simple_prime3D_init'
        write(*,'(A)') 'simple_print_cmd_dict'
        write(*,'(A)') 'simple_print_dose_weights'
        write(*,'(A)') 'simple_print_fsc'
        write(*,'(A)') 'simple_projvol'
        write(*,'(A)') 'simple_rank_cavgs'
        write(*,'(A)') 'simple_recvol'
        write(*,'(A)') 'simple_res'
        write(*,'(A)') 'simple_resrange'
        write(*,'(A)') 'simple_rotmats2oris'
        write(*,'(A)') 'simple_scale'
        write(*,'(A)') 'simple_select'
        write(*,'(A)') 'simple_select_frames'
        write(*,'(A)') 'simple_shellweight3D'
        write(*,'(A)') 'simple_shift'
        write(*,'(A)') 'simple_simimgs'
        write(*,'(A)') 'simple_simmovie'
        write(*,'(A)') 'simple_simsubtomo'
        write(*,'(A)') 'simple_split'
        write(*,'(A)') 'simple_split_pairs'
        write(*,'(A)') 'simple_stack'
        write(*,'(A)') 'simple_stackops'
        write(*,'(A)') 'simple_symsrch'
        write(*,'(A)') 'simple_tseries_split'
        write(*,'(A)') 'simple_unblur_movies'
        write(*,'(A)') 'simple_volassemble'
        write(*,'(A)') 'simple_volaverager'
        write(*,'(A)') 'simple_volops'
        write(*,'(A)') 'simple_volume_smat'
        stop
    end subroutine list_all_simple_programs

end module simple_gen_doc
