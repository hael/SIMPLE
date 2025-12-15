module simple_core_module_api
use simple_binoris,         only: binoris, binoris_seginfo
use simple_chash,           only: chash
use simple_class_sample_io, only: print_class_sample, class_samples_same, write_class_samples, read_class_samples, deallocate_class_samples
use simple_error,           only: simple_exception
use simple_estimate_ssnr,   only: fsc2optlp, fsc2optlp_sub, gaussian_filter, fsc2boostfilter, get_resolution, mskdiam2lplimits,&
                                 &get_resolution_at_fsc, lpstages, lpstages_fast, edit_lpstages4polar, mskdiam2streamresthreshold,&
                                 &calc_dose_weights
use simple_fileio,          only: add2fbody, append2basename, arr2file, arr2txtfile, basename, del_file, del_files, fclose, file2drarr, file2rarr,&
                                 &file2rmat, file_exists, fileiochk, filepath, fname2ext, fname2format, fname_new_ext, fopen, get_fbody, get_fpath,&
                                 &move_files2dir, nlines, read_filetable, rmat2file, simple_abspath, simple_chdir, simple_chmod, simple_copy_file,&
                                 &simple_getenv, simple_list_dirs, simple_list_files, simple_list_files_regexp, simple_rename, stemname, swap_suffix,&
                                 &wait_for_closure, simple_touch, simple_rmdir, simple_getcwd, write_filetable, write_singlelineoftext, move_files_in_cwd,&
                                 &read_exit_code
use simple_hash,            only: hash
use simple_is_check_assert, only: is_a_number, is_zero, is_gt_zero, is_equal, is_even, check4nans3D, check4nans2D, check4nans, assert_eq,&
                                 &is_odd, is_even
use simple_imghead,         only: ImgHead, MrcImgHead, SpiImgHead, TiffImgHead, find_ldim_nptcls, update_stack_nimgs
use simple_jiffys,          only: progress, progress_gfortran, simple_end, swap
use simple_kbinterpol,      only: kbinterpol
use simple_linalg,          only: eigsrt, jacobi, matinv, norm_2, svdcmp, svdfit, svd_multifit, euclid, hyp, myacos, deg2rad, rad2deg, pythag,&
                                 &eigh, arg, fit_lsq_plane, fit_straight_line, plane_from_points, projz, trace, ang2vox, vox2ang
use simple_magic_boxes,     only: magic_pftsz, find_larger_magic_box, find_magic_box, print_magic_box_range, autoscale
use simple_map_reduce,      only: split_nobjs_even
use simple_math,            only: otsu, pixels_dist, equispaced_vals, put_last, bounds_from_mask3D, elim_dup, mode, sortmeans,&
                                 &quantize_vec_serial, quantize_vec, quadri, create_hist_vector, round2even, round2odd, rotmat2d, gauwfun,&
                                 &gaussian1D, gaussian2D, gaussian3D, shft, cross, get_pixel_pos
use simple_math_ft,         only: csq_fast, csq, calc_fourier_index, calc_graphene_mask, calc_lowpass_lim, cyci_1d, cyci_1d_static,&
                                 &fdim, get_find_at_res, get_find_at_crit, get_resarr, mycabs, phase_angle, resang
use simple_nrtxtfile,       only: nrtxtfile
use simple_ori,             only: ori, euler2m, m2euler, m2euler_fast, euler_dist, euler_inplrotdist, euler_compose, euler_mirror, geodesic_frobdev
use simple_oris,            only: oris
use simple_ran_tabu,        only: ran_tabu
use simple_rnd,             only: ran3arr, gasdev, greedy_sampling, seed_rnd, ran3, irnd_uni, irnd_uni_pair, shcloc, rnd_inds, mnorm_smp,&
                                 &multinomal
use simple_edges_sqwins,    only: cosedge, cosedge_inner, hardedge, hardedge_inner, sqwin_1d, sqwin_2d, sqwin_3d
use simple_srch_sort_loc,   only: find, hpsort, locate, reverse, maxnloc, min3, minnloc, scores2order, dists2order, mask2inds, reorder, selec,&
                                 &reverse_f
use simple_stat,            only: avg_sdev, moment, skewness, kurtosis, pearsn, normalize, normalize_minmax, merge_dmats,&
                                 &avg_frac_smallest, pearsn_serial, std_mean_diff, calc_stats, corrs2weights, kstwo, analyze_smat,&
                                 &dmat2smat, smat2dmat, scores2scores_percen, dists2scores_percen, merge_smats, medoid_from_smat,&
                                 &medoid_from_dmat, median, median_nocopy, conv2rank_weights, rank_sum_weights, rank_centroid_weights,&
                                 &rank_exponent_weights, rank_inverse_weights, z_scores, mad_gau, robust_scaling
use simple_string,          only: string
use simple_string_utils,    only: str2format, str2int, str2real, real2str, findloc_str, spaces, char_is_a_letter, char_is_a_number,&
                               &str_has_substr, list_of_ints2arr, int2str, int2str_pad, map_str_nrs, to_cstring, lex_sort, upperCase,&
                               &str_pad, lowercase, parsestr, split, split_str
use simple_sym,             only: sym
use simple_syslib,          only: is_open, syslib_symlink, exec_cmdline, simple_mkdir, get_process_id, find_next_int_dir_prefix, dir_exists,&
                                 &simple_file_stat
use simple_timer,           only: timer_int_kind, tic, toc, simple_gettime, cast_time_char
end module simple_core_module_api