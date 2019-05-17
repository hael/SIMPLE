import * as fs from 'fs'

class Module {
	
	constructor(){
		console.log("Loaded module Relion")
	}
	
	metadata = {
			"moduletitle" :"Relion",
			"modulename" :"relion",
			"tasks" : {
				"Refine_MPI" : {
					"name":"Refine_MPI",
					"descr_short":"Refine MPI",
					"descr_long":"Refine MPI",
					"executable":"relion_refine_mpi",
					"pages" : [
						{
							"title":"General options",
							"keys" : [
								{
									"key":"i",
									"keytype":"file",
									"descr_short":"Input images (in a star-file or a stack)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"o",
									"keytype":"text",
									"descr_short":"Output rootname",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"iter",
									"keytype":"num",
									"descr_short":"Maximum number of iterations to perform",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"angpix",
									"keytype":"num",
									"descr_short":"Pixel size (in Angstroms)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"tau2_fudge",
									"keytype":"num",
									"descr_short":"Regularisation parameter (values higher than 1 give more weight to the data)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"K",
									"keytype":"num",
									"descr_short":"Number of references to be refined",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"particle_diameter",
									"keytype":"num",
									"descr_short":"Diameter of the circular mask that will be applied to the experimental images (in Angstroms)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"zero_mask",
									"keytype":"binary",
									"descr_short":"Mask surrounding background in particles to zero (by default the solvent area is filled with random noise)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"tau",
									"keytype":"file",
									"descr_short":"STAR file with input tau2-spectrum (to be kept constant)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"flatten_solvent",
									"keytype":"binary",
									"descr_short":"Perform masking on the references as well?",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"solvent_mask",
									"keytype":"file",
									"descr_short":"User-provided mask for the references (default is to use spherical mask with particle_diameter)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"solvent_mask2",
									"keytype":"file",
									"descr_short":"User-provided secondary mask (with its own average density)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"local_symmetry",
									"keytype":"file",
									"descr_short":"Local symmetry description file containing list of masks and their operators",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"split_random_halves",
									"keytype":"binary",
									"descr_short":"Refine two random halves of the data completely separately",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"low_resol_join_halves",
									"keytype":"binary",
									"descr_short":"Resolution (in Angstrom) up to which the two random half-reconstructions will not be independent to prevent diverging orientations",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								}
							]
						},
						{
							"title":"Initialisation",
							"keys" : [
								{
									"key":"ref",
									"keytype":"file",
									"descr_short":"Image, stack or star-file with the reference(s). (Compulsory for 3D refinement!)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"denovo_3dref",
									"keytype":"binary",
									"descr_short":"Make an initial 3D model from randomly oriented 2D particles",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"offset",
									"keytype":"num",
									"descr_short":"Initial estimated stddev for the origin offsets",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"firstiter_cc",
									"keytype":"binary",
									"descr_short":"Perform CC-calculation in the first iteration (use this if references are not on the absolute intensity scale)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"ini_high",
									"keytype":"file",
									"descr_short":"Resolution (in Angstroms) to which to limit refinement in the first iteration",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								}
							]
						},
						{
							"title":"Orientations",
							"keys" : [
								{
									"key":"oversampling",
									"keytype":"num",
									"descr_short":"Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"healpix_order",
									"keytype":"num",
									"descr_short":"Healpix order for the angular sampling (before oversampling) on the (3D) sphere: hp2=15deg, hp3=7.5deg, etc",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"psi_step",
									"keytype":"num",
									"descr_short":"Sampling rate (before oversampling) for the in-plane angle (default=10deg for 2D, hp sampling for 3D)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"limit_tilt",
									"keytype":"num",
									"descr_short":"Limited tilt angle: positive for keeping side views, negative for keeping top views",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"sym",
									"keytype":"text",
									"descr_short":"Symmetry group",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"offset_range",
									"keytype":"file",
									"descr_short":"Search range for origin offsets (in pixels)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"offset_step",
									"keytype":"num",
									"descr_short":"Sampling rate (before oversampling) for origin offsets (in pixels)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"helical_offset_step",
									"keytype":"num",
									"descr_short":"Sampling rate (before oversampling) for offsets along helical axis (in pixels)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"perturb",
									"keytype":"num",
									"descr_short":"Perturbation factor for the angular sampling (0=no perturb; 0.5=perturb)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"auto_refine",
									"keytype":"binary",
									"descr_short":"Perform 3D auto-refine procedure?",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"auto_local_healpix_order",
									"keytype":"num",
									"descr_short":"Minimum healpix order (before oversampling) from which autosampling procedure will use local searches",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"sigma_ang",
									"keytype":"num",
									"descr_short":"Stddev on all three Euler angles for local angular searches (of +/- 3 stddev)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"sigma_rot",
									"keytype":"num",
									"descr_short":"Stddev on the first Euler angle for local angular searches (of +/- 3 stddev)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"sigma_tilt",
									"keytype":"num",
									"descr_short":"Stddev on the second Euler angle for local angular searches (of +/- 3 stddev)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"sigma_psi",
									"keytype":"num",
									"descr_short":"Stddev on the in-plane angle for local angular searches (of +/- 3 stddev)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"skip_align",
									"keytype":"binary",
									"descr_short":"Skip orientational assignment (only classify)?",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"skip_rotate",
									"keytype":"binary",
									"descr_short":"Skip rotational assignment (only translate and classify)?",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								},
								{
									"key":"bimodal_psi",
									"keytype":"binary",
									"descr_short":" Do bimodal searches of psi angle?",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								}
							]
						},
						{
							"title":"Helical reconstruction",
							"keys" : [
								{
									"key":"oversampling",
									"keytype":"num",
									"descr_short":"Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)",
									"descr_long":"",
									"descr_placeholder":"",
									"required":""
								}
							]
						}
					]
				},
				"Refine" : {
					"name":"Refine",
					"descr_short":"Refine MPI",
					"descr_long":"Refine MPI",
					"executable":"relion_refine_mpi",
					"pages" : [
					]
				}
			}
		}
		
	relion_refine_mpi(modules, arg) {
		var command = "relion_refine_mpi"
		var keys = Object.keys(arg);
		for(var key of keys){
			if(arg[key] != ""){
				command += " --" + key.replace('key', '') + " " + arg[key]
			}
		}
		return (modules['available']['core']['createTask'](modules, arg))
	}
}

module.exports = new Module()

/*

====== Helical reconstruction (in development...) ===== 
                    helix (false) : Perform 3D classification or refinement for helices?
  --ignore_helical_symmetry (false) : Ignore helical symmetry?
               --helical_nr_asu (1) : Number of new helical asymmetric units (asu) per box (1 means no helical symmetry is present)
       --helical_twist_initial (0.) : Helical twist (in degrees, positive values for right-handedness)
           --helical_twist_min (0.) : Minimum helical twist (in degrees, positive values for right-handedness)
           --helical_twist_max (0.) : Maximum helical twist (in degrees, positive values for right-handedness)
       --helical_twist_inistep (0.) : Initial step of helical twist search (in degrees)
        --helical_rise_initial (0.) : Helical rise (in Angstroms)
            --helical_rise_min (0.) : Minimum helical rise (in Angstroms)
            --helical_rise_max (0.) : Maximum helical rise (in Angstroms)
        --helical_rise_inistep (0.) : Initial step of helical rise search (in Angstroms)
       --helical_z_percentage (0.3) : This box length along the center of Z axis contains good information of the helix. Important in imposing and refining symmetry
     --helical_inner_diameter (-1.) : Inner diameter of helical tubes in Angstroms (for masks of helical references and particles)
     --helical_outer_diameter (-1.) : Outer diameter of helical tubes in Angstroms (for masks of helical references and particles)
  --helical_symmetry_search (false) : Perform local refinement of helical symmetry?
     --helical_sigma_distance (-1.) : Sigma of distance along the helical tracks
  --helical_keep_tilt_prior_fixed (false) : Keep helical tilt priors fixed (at 90 degrees) in global angular searches?
====== Corrections ===== 
                      --ctf (false) : Perform CTF correction?
    --ctf_intact_first_peak (false) : Ignore CTFs until their first peak?
        --ctf_corrected_ref (false) : Have the input references been CTF-amplitude corrected?
        --ctf_phase_flipped (false) : Have the data been CTF phase-flipped?
           --ctf_multiplied (false) : Have the data been premultiplied with their CTF?
         --only_flip_phases (false) : Only perform CTF phase-flipping? (default is full amplitude-correction)
                     --norm (false) : Perform normalisation-error correction?
                    --scale (false) : Perform intensity-scale corrections on image groups?
                  --no_norm (false) : Switch off normalisation-error correction?
                 --no_scale (false) : Switch off intensity-scale corrections on image groups?
====== Stochastic Gradient Descent ===== 
                      --sgd (false) : Perform stochastic gradient descent instead of default expectation-maximization
                         --mu (0.9) : Momentum parameter for SGD updates
                 --subset_size (-1) : Size of the subsets for SGD
               --sgd_stepsize (0.5) : Step size parameter for SGD updates
                 --max_subsets (-1) : Stop SGD after processing this many subsets (possibly more than 1 iteration)
      --sgd_sigma2fudge_initial (8) : Initial factor by which the noise variance will be multiplied for SGD (not used if halftime is negative)
    --sgd_sigma2fudge_halflife (-1) : Initialise SGD with 8x higher noise-variance, and reduce with this half-life in # of particles (default is keep normal variance)
               --write_subsets (-1) : Write out model every so many subsets (default is not writing any)
          --strict_highres_sgd (20) : Resolution limit (in Angstrom) to restrict probability calculations in SGD
====== Computation ===== 
                         --pool (1) : Number of images to pool for each thread task
                            --j (1) : Number of threads to run in parallel (only useful on multi-core machines)
  --dont_combine_weights_via_disc (false) : Send the large arrays of summed weights through the MPI network, instead of writing large files to disc
          --onthefly_shifts (false) : Calculate shifted images on-the-fly, do not store precalculated ones in memory
      --no_parallel_disc_io (false) : Do NOT let parallel (MPI) processes access the disc simultaneously (use this option with NFS)
           --preread_images (false) : Use this to let the master process read all particles into memory. Be careful you have enough RAM for large data sets!
                   --scratch_dir () : If provided, particle stacks will be copied to this local scratch disk prior to refinement.
           --keep_free_scratch (10) : Space available for copying particle stacks (in Gb)
            --reuse_scratch (false) : Re-use data on scratchdir, instead of wiping it and re-copying all data.
                      --gpu (false) : Use available gpu resources for some calculations
              --free_gpu_memory (0) : GPU device memory (in Mb) to leave free after allocation.
====== Expert options ===== 
                          --pad (2) : Oversampling factor for the Fourier transforms of the references
                       --NN (false) : Perform nearest-neighbour instead of linear Fourier-space interpolation?
                    --r_min_nn (10) : Minimum number of Fourier shells to perform linear Fourier-space interpolation
                         --verb (1) : Verbosity (1=normal, 0=silent)
                 --random_seed (-1) : Number for the random seed generator
                 --coarse_size (-1) : Maximum image size for the first pass of the adaptive sampling approach
        --adaptive_fraction (0.999) : Fraction of the weights to be considered in the first pass of adaptive oversampling 
                     --maskedge (5) : Width of the soft edge of the spherical mask (in pixels)
          --fix_sigma_noise (false) : Fix the experimental noise spectra?
         --fix_sigma_offset (false) : Fix the stddev in the origin offsets?
                   --incr_size (10) : Number of Fourier shells beyond the current resolution to be included in refinement
    --print_metadata_labels (false) : Print a table with definitions of all metadata labels, and exit
       --print_symmetry_ops (false) : Print all symmetry transformation matrices, and exit
          --strict_highres_exp (-1) : Resolution limit (in Angstrom) to restrict probability calculations in the expectation step
          --dont_check_norm (false) : Skip the check whether the images are normalised correctly
                --always_cc (false) : Perform CC-calculation in all iterations (useful for faster denovo model generation?)
      --solvent_correct_fsc (false) : Correct FSC curve for the effects of the solvent mask?
            --skip_maximize (false) : Skip maximization step (only write out data.star file)?
          --failsafe_threshold (40) : Maximum number of particles permitted to be handled by fail-safe mode, due to zero sum of weights, before exiting with an error (GPU only).
====== MPI options ===== 
  --only_do_unfinished_movies (false) : When processing movies on a per-micrograph basis, ignore those movies for which the output STAR file already exists.
*/
