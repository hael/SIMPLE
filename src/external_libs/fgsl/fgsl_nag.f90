module fgsl
!-------------------------------------------------------------------------------
! Interface module for use of GSL from Fortran
! Author: R. Bader, Leibniz Supercomputing Centre, Garching, Germany
! Notes:
! (1) Conventions used:
! api fgsl.f90 fgsl_nag.f90 fgsl_utils.c interface preprocess_fgsl.sh GSL_* must be replaced by FGSL_* for each API call, abstract
! data type and parameters (with exception of the M_* mathematical
! constants)
! api fgsl.f90 fgsl_nag.f90 fgsl_utils.c interface preprocess_fgsl.sh Some names were changed due to UC/LC aliasing. See the documentation
! chapter on special functions for details
! api fgsl.f90 fgsl_nag.f90 fgsl_utils.c interface preprocess_fgsl.sh Type matching:
! - real(fgsl_double) is used for double
! - real(fgsl_float) is used for float
! - integer(fgsl_int) for integer
! - integer(fgsl_long) for long integer
! - integer(fgsl_size_t) for size_t integer
! - complex(fgsl_double_complex) for gsl_complex
! - character(fgsl_char) for characters
! - no value attributes and mostly no pointers in Fortran calls
! - unsigned int must be converted to integer(fgsl_long).
! - char api fgsl.f90 fgsl_nag.f90 fgsl_utils.c interface preprocess_fgsl.sh results are converted to fixed length strings. Use TRIM.
! (2) Additional routines
! api fgsl.f90 fgsl_nag.f90 fgsl_utils.c interface preprocess_fgsl.sh Generic interface fgsl_well_defined for checking status of FGSL
! objects (which are typically opaque).
! api fgsl.f90 fgsl_nag.f90 fgsl_utils.c interface preprocess_fgsl.sh See api/array.finc for array alignment routines.
! api fgsl.f90 fgsl_nag.f90 fgsl_utils.c interface preprocess_fgsl.sh See api/math.finc for function object constructors.
! api fgsl.f90 fgsl_nag.f90 fgsl_utils.c interface preprocess_fgsl.sh See api/io.finc for I/O related add-ons.
! (3) Only the manual-documented interface are implemented. The
! C include files may contain more stuff which may only be meant
! for internal use. Or is not documented yet. Tough.
! (4) Inlining of GSL routines is of course not possible
! (5) Macros
! api fgsl.f90 fgsl_nag.f90 fgsl_utils.c interface preprocess_fgsl.sh are not supported
! api fgsl.f90 fgsl_nag.f90 fgsl_utils.c interface preprocess_fgsl.sh macro values are replicated as parameters
! api fgsl.f90 fgsl_nag.f90 fgsl_utils.c interface preprocess_fgsl.sh Inf/Nan need to use IEEE_VALUE (if available)
!-------------------------------------------------------------------------------
use, intrinsic :: iso_c_binding
implicit none
private
!
! Exported entities
!
public :: fgsl_set_error_handler_off, fgsl_set_error_handler, fgsl_strerror, &
fgsl_error_handler_init, fgsl_name, fgsl_error
public :: fgsl_well_defined, fgsl_obj_c_ptr, fgsl_sizeof
public :: fgsl_open, fgsl_close, fgsl_stdin, fgsl_stdout, fgsl_stderr, fgsl_flush
public :: assignment(=)
public :: fgsl_isnan, fgsl_isinf, fgsl_finite, fgsl_log1p, fgsl_expm1, &
fgsl_hypot, fgsl_acosh, fgsl_asinh, fgsl_atanh, fgsl_ldexp, fgsl_frexp, &
fgsl_fcmp, fgsl_function_init, fgsl_function_fdf_init, fgsl_function_free, &
fgsl_function_fdf_free, fgsl_fn_eval, fgsl_fn_fdf_eval_f, fgsl_fn_fdf_eval_df, &
fgsl_fn_fdf_eval_f_df
! complex functions
public :: fgsl_complex_arg, fgsl_complex_logabs, fgsl_complex_log10, &
fgsl_complex_log_b, fgsl_complex_arcsin, fgsl_complex_arcsin_real, &
fgsl_complex_arccos, fgsl_complex_arccos_real, fgsl_complex_arctan, &
fgsl_complex_arcsec, fgsl_complex_arcsec_real, fgsl_complex_arccsc, &
fgsl_complex_arccsc_real, fgsl_complex_arccot, fgsl_complex_arcsinh, &
fgsl_complex_arccosh, fgsl_complex_arccosh_real, fgsl_complex_arctanh, &
fgsl_complex_arctanh_real, fgsl_complex_arcsech, fgsl_complex_arccsch, &
fgsl_complex_arccoth
! polynomials
public :: fgsl_poly_eval, fgsl_poly_complex_eval, fgsl_complex_poly_complex_eval, &
fgsl_poly_eval_derivs, fgsl_poly_dd_init, fgsl_poly_dd_eval, fgsl_poly_dd_taylor, &
fgsl_poly_solve_quadratic, fgsl_poly_complex_solve_quadratic, &
fgsl_poly_solve_cubic, fgsl_poly_complex_solve_cubic, &
fgsl_poly_complex_workspace_alloc, fgsl_poly_complex_workspace_free, &
fgsl_poly_complex_solve
! special functions
public :: fgsl_sf_airy_ai, fgsl_sf_airy_ai_e, fgsl_sf_airy_bi, fgsl_sf_airy_bi_e, &
fgsl_sf_airy_ai_scaled, fgsl_sf_airy_ai_scaled_e, fgsl_sf_airy_bi_scaled, &
fgsl_sf_airy_bi_scaled_e, fgsl_sf_airy_ai_deriv, fgsl_sf_airy_bi_deriv, &
fgsl_sf_airy_ai_deriv_scaled, fgsl_sf_airy_bi_deriv_scaled, &
fgsl_sf_airy_ai_deriv_e, fgsl_sf_airy_bi_deriv_e, &
fgsl_sf_airy_ai_deriv_scaled_e, fgsl_sf_airy_bi_deriv_scaled_e, &
fgsl_sf_airy_zero_ai, fgsl_sf_airy_zero_ai_e, &
fgsl_sf_airy_zero_bi, fgsl_sf_airy_zero_bi_e, &
fgsl_sf_airy_zero_ai_deriv, fgsl_sf_airy_zero_ai_deriv_e, &
fgsl_sf_airy_zero_bi_deriv, fgsl_sf_airy_zero_bi_deriv_e
public :: fgsl_sf_bessel_jc0, fgsl_sf_bessel_jc0_e, &
fgsl_sf_bessel_jc1, fgsl_sf_bessel_jc1_e, fgsl_sf_bessel_jcn, &
fgsl_sf_bessel_jcn_e, fgsl_sf_bessel_jcn_array, &
fgsl_sf_bessel_yc0, fgsl_sf_bessel_yc0_e, &
fgsl_sf_bessel_yc1, fgsl_sf_bessel_yc1_e, fgsl_sf_bessel_ycn, &
fgsl_sf_bessel_ycn_e, fgsl_sf_bessel_ycn_array, &
fgsl_sf_bessel_ic0, fgsl_sf_bessel_ic0_e, &
fgsl_sf_bessel_ic1, fgsl_sf_bessel_ic1_e, fgsl_sf_bessel_icn, &
fgsl_sf_bessel_icn_e, fgsl_sf_bessel_icn_array, &
fgsl_sf_bessel_ic0_scaled, fgsl_sf_bessel_ic0_scaled_e, &
fgsl_sf_bessel_ic1_scaled, fgsl_sf_bessel_ic1_scaled_e, fgsl_sf_bessel_icn_scaled, &
fgsl_sf_bessel_icn_scaled_e, fgsl_sf_bessel_icn_scaled_array, &
fgsl_sf_bessel_kc0, fgsl_sf_bessel_kc0_e, &
fgsl_sf_bessel_kc1, fgsl_sf_bessel_kc1_e, fgsl_sf_bessel_kcn, &
fgsl_sf_bessel_kcn_e, fgsl_sf_bessel_kcn_array, &
fgsl_sf_bessel_kc0_scaled, fgsl_sf_bessel_kc0_scaled_e, &
fgsl_sf_bessel_kc1_scaled, fgsl_sf_bessel_kc1_scaled_e, fgsl_sf_bessel_kcn_scaled, &
fgsl_sf_bessel_kcn_scaled_e, fgsl_sf_bessel_kcn_scaled_array
public :: fgsl_sf_bessel_js0, fgsl_sf_bessel_js0_e, fgsl_sf_bessel_js1, fgsl_sf_bessel_js1_e, &
fgsl_sf_bessel_js2, fgsl_sf_bessel_js2_e, fgsl_sf_bessel_jsl, fgsl_sf_bessel_jsl_e, &
fgsl_sf_bessel_jsl_array, fgsl_sf_bessel_jsl_steed_array, &
fgsl_sf_bessel_ys0, fgsl_sf_bessel_ys0_e, fgsl_sf_bessel_ys1, fgsl_sf_bessel_ys1_e, &
fgsl_sf_bessel_ys2, fgsl_sf_bessel_ys2_e, fgsl_sf_bessel_ysl, fgsl_sf_bessel_ysl_e, &
fgsl_sf_bessel_ysl_array, fgsl_sf_bessel_is0_scaled, fgsl_sf_bessel_is0_scaled_e, &
fgsl_sf_bessel_is1_scaled, fgsl_sf_bessel_is1_scaled_e, &
fgsl_sf_bessel_is2_scaled, fgsl_sf_bessel_is2_scaled_e, fgsl_sf_bessel_isl_scaled, &
fgsl_sf_bessel_isl_scaled_e, fgsl_sf_bessel_isl_scaled_array, &
fgsl_sf_bessel_ks0_scaled, fgsl_sf_bessel_ks0_scaled_e, &
fgsl_sf_bessel_ks1_scaled, fgsl_sf_bessel_ks1_scaled_e, &
fgsl_sf_bessel_ks2_scaled, fgsl_sf_bessel_ks2_scaled_e, fgsl_sf_bessel_ksl_scaled, &
fgsl_sf_bessel_ksl_scaled_e, fgsl_sf_bessel_ksl_scaled_array
public :: fgsl_sf_bessel_jnu, fgsl_sf_bessel_jnu_e, fgsl_sf_bessel_sequence_jnu_e, &
fgsl_sf_bessel_ynu, fgsl_sf_bessel_ynu_e, fgsl_sf_bessel_inu, fgsl_sf_bessel_inu_e, &
fgsl_sf_bessel_inu_scaled, fgsl_sf_bessel_inu_scaled_e, &
fgsl_sf_bessel_knu, fgsl_sf_bessel_knu_e, &
fgsl_sf_bessel_lnknu, fgsl_sf_bessel_lnknu_e, &
fgsl_sf_bessel_knu_scaled, fgsl_sf_bessel_knu_scaled_e, &
fgsl_sf_bessel_zero_jc0, fgsl_sf_bessel_zero_jc0_e, &
fgsl_sf_bessel_zero_jc1, fgsl_sf_bessel_zero_jc1_e, &
fgsl_sf_bessel_zero_jnu, fgsl_sf_bessel_zero_jnu_e
public :: fgsl_sf_clausen, fgsl_sf_clausen_e, fgsl_sf_hydrogenicr_1, &
fgsl_sf_hydrogenicr_1_e, fgsl_sf_hydrogenicr, fgsl_sf_hydrogenicr_e, &
fgsl_sf_coulomb_wave_fg_e, fgsl_sf_coulomb_wave_f_array, fgsl_sf_coulomb_wave_fg_array, &
fgsl_sf_coulomb_wave_fgp_array, fgsl_sf_coulomb_wave_sphf_array, &
fgsl_sf_coulomb_cl_e, fgsl_sf_coulomb_cl_array
public :: fgsl_sf_coupling_3j, fgsl_sf_coupling_3j_e, fgsl_sf_coupling_6j, &
fgsl_sf_coupling_6j_e, fgsl_sf_coupling_9j, fgsl_sf_coupling_9j_e, &
fgsl_sf_dawson, fgsl_sf_dawson_e, fgsl_sf_debye_1, fgsl_sf_debye_1_e, &
fgsl_sf_debye_2, fgsl_sf_debye_2_e, fgsl_sf_debye_3, fgsl_sf_debye_3_e, &
fgsl_sf_debye_4, fgsl_sf_debye_4_e, fgsl_sf_debye_5, fgsl_sf_debye_5_e, &
fgsl_sf_debye_6, fgsl_sf_debye_6_e, fgsl_sf_dilog, fgsl_sf_dilog_e, &
fgsl_sf_complex_dilog_e, fgsl_sf_multiply_e, fgsl_sf_multiply_err_e
public :: fgsl_sf_ellint_kcomp, fgsl_sf_ellint_kcomp_e, fgsl_sf_ellint_ecomp, &
fgsl_sf_ellint_ecomp_e, fgsl_sf_ellint_f, fgsl_sf_ellint_pcomp, &
fgsl_sf_ellint_pcomp_e, fgsl_sf_ellint_f_e, &
fgsl_sf_ellint_e, fgsl_sf_ellint_e_e, fgsl_sf_ellint_p, fgsl_sf_ellint_p_e, &
fgsl_sf_ellint_d, fgsl_sf_ellint_d_e, fgsl_sf_ellint_rc, fgsl_sf_ellint_rc_e, &
fgsl_sf_ellint_rd, fgsl_sf_ellint_rd_e, fgsl_sf_ellint_rf, fgsl_sf_ellint_rf_e, &
fgsl_sf_ellint_rj, fgsl_sf_ellint_rj_e, fgsl_sf_elljac_e
public :: fgsl_sf_erf, fgsl_sf_erf_e, fgsl_sf_erfc, fgsl_sf_erfc_e, fgsl_sf_log_erfc, &
fgsl_sf_log_erfc_e, fgsl_sf_erf_z, fgsl_sf_erf_z_e, fgsl_sf_erf_q, fgsl_sf_erf_q_e, &
fgsl_sf_hazard, fgsl_sf_hazard_e, fgsl_sf_exp, fgsl_sf_exp_e, fgsl_sf_exp_e10_e, &
fgsl_sf_exp_mult, fgsl_sf_exp_mult_e, fgsl_sf_exp_mult_e10_e, fgsl_sf_expm1, &
fgsl_sf_expm1_e, fgsl_sf_exprel, fgsl_sf_exprel_e, fgsl_sf_exprel_2, fgsl_sf_exprel_2_e, &
fgsl_sf_exprel_n, fgsl_sf_exprel_n_e, fgsl_sf_exp_err_e, fgsl_sf_exp_err_e10_e, &
fgsl_sf_exp_mult_err_e, fgsl_sf_exp_mult_err_e10_e
public :: fgsl_sf_expint_e1, fgsl_sf_expint_e1_e, fgsl_sf_expint_e2, fgsl_sf_expint_e2_e, &
fgsl_sf_expint_en, fgsl_sf_expint_en_e, &
fgsl_sf_expint_ei, fgsl_sf_expint_ei_e, fgsl_sf_shi, fgsl_sf_shi_e, fgsl_sf_chi, &
fgsl_sf_chi_e, fgsl_sf_expint_3, fgsl_sf_expint_3_e, fgsl_sf_si, fgsl_sf_si_e, &
fgsl_sf_ci, fgsl_sf_ci_e, fgsl_sf_atanint, fgsl_sf_atanint_e
public :: fgsl_sf_fermi_dirac_m1, fgsl_sf_fermi_dirac_m1_e, fgsl_sf_fermi_dirac_0, &
fgsl_sf_fermi_dirac_0_e, fgsl_sf_fermi_dirac_1, fgsl_sf_fermi_dirac_1_e, &
fgsl_sf_fermi_dirac_2, fgsl_sf_fermi_dirac_2_e, fgsl_sf_fermi_dirac_int, &
fgsl_sf_fermi_dirac_int_e, fgsl_sf_fermi_dirac_mhalf, fgsl_sf_fermi_dirac_mhalf_e, &
fgsl_sf_fermi_dirac_half, fgsl_sf_fermi_dirac_half_e, fgsl_sf_fermi_dirac_3half, &
fgsl_sf_fermi_dirac_3half_e, fgsl_sf_fermi_dirac_inc_0, &
fgsl_sf_fermi_dirac_inc_0_e
public :: fgsl_sf_gamma, fgsl_sf_gamma_e, fgsl_sf_lngamma, fgsl_sf_lngamma_e, &
fgsl_sf_lngamma_sgn_e, fgsl_sf_gammastar, fgsl_sf_gammastar_e, fgsl_sf_gammainv, &
fgsl_sf_gammainv_e, fgsl_sf_lngamma_complex_e, fgsl_sf_fact, fgsl_sf_fact_e, &
fgsl_sf_doublefact, fgsl_sf_doublefact_e, fgsl_sf_lnfact, fgsl_sf_lnfact_e, &
fgsl_sf_lndoublefact, fgsl_sf_lndoublefact_e, fgsl_sf_choose, fgsl_sf_choose_e, &
fgsl_sf_lnchoose, fgsl_sf_lnchoose_e, fgsl_sf_taylorcoeff, fgsl_sf_taylorcoeff_e, &
fgsl_sf_poch, fgsl_sf_poch_e, fgsl_sf_lnpoch, fgsl_sf_lnpoch_e, &
fgsl_sf_lnpoch_sgn_e, fgsl_sf_pochrel, fgsl_sf_pochrel_e, &
fgsl_sf_gamma_inc, fgsl_sf_gamma_inc_e, &
fgsl_sf_gamma_inc_q, fgsl_sf_gamma_inc_q_e, fgsl_sf_gamma_inc_p, &
fgsl_sf_gamma_inc_p_e, fgsl_sf_beta, fgsl_sf_beta_e, fgsl_sf_lnbeta, &
fgsl_sf_lnbeta_e, fgsl_sf_beta_inc, fgsl_sf_beta_inc_e
public :: fgsl_sf_gegenpoly_1, fgsl_sf_gegenpoly_1_e, fgsl_sf_gegenpoly_2, &
fgsl_sf_gegenpoly_2_e, fgsl_sf_gegenpoly_3, fgsl_sf_gegenpoly_3_e, &
fgsl_sf_gegenpoly_n, fgsl_sf_gegenpoly_n_e, fgsl_sf_gegenpoly_array, &
fgsl_sf_hyperg_0f1, fgsl_sf_hyperg_0f1_e, fgsl_sf_hyperg_1f1_int, &
fgsl_sf_hyperg_1f1_int_e, fgsl_sf_hyperg_1f1, fgsl_sf_hyperg_1f1_e, &
fgsl_sf_hyperg_u_int, fgsl_sf_hyperg_u_int_e, fgsl_sf_hyperg_u_int_e10_e, &
fgsl_sf_hyperg_u, fgsl_sf_hyperg_u_e, fgsl_sf_hyperg_u_e10_e, &
fgsl_sf_hyperg_2f1, fgsl_sf_hyperg_2f1_e, fgsl_sf_hyperg_2f1_conj, &
fgsl_sf_hyperg_2f1_conj_e, fgsl_sf_hyperg_2f1_renorm, fgsl_sf_hyperg_2f1_renorm_e, &
fgsl_sf_hyperg_2f1_conj_renorm, fgsl_sf_hyperg_2f1_conj_renorm_e, &
fgsl_sf_hyperg_2f0, fgsl_sf_hyperg_2f0_e
public :: fgsl_sf_laguerre_1, fgsl_sf_laguerre_1_e, fgsl_sf_laguerre_2, &
fgsl_sf_laguerre_2_e, fgsl_sf_laguerre_3, fgsl_sf_laguerre_3_e, &
fgsl_sf_laguerre_n, fgsl_sf_laguerre_n_e
public :: fgsl_sf_lambert_w0, fgsl_sf_lambert_w0_e, fgsl_sf_lambert_wm1, &
fgsl_sf_lambert_wm1_e, fgsl_sf_legendre_p1, fgsl_sf_legendre_p1_e, &
fgsl_sf_legendre_p2, fgsl_sf_legendre_p2_e, fgsl_sf_legendre_p3, &
fgsl_sf_legendre_p3_e, fgsl_sf_legendre_pl, fgsl_sf_legendre_pl_e, &
fgsl_sf_legendre_pl_array, fgsl_sf_legendre_pl_deriv_array, &
fgsl_sf_legendre_q0, fgsl_sf_legendre_q0_e, &
fgsl_sf_legendre_q1, fgsl_sf_legendre_q1_e, fgsl_sf_legendre_ql, &
fgsl_sf_legendre_ql_e, fgsl_sf_legendre_plm, fgsl_sf_legendre_plm_e, &
fgsl_sf_legendre_plm_array, fgsl_sf_legendre_plm_deriv_array, &
fgsl_sf_legendre_sphplm, fgsl_sf_legendre_sphplm_e, &
fgsl_sf_legendre_sphplm_array, fgsl_sf_legendre_sphplm_deriv_array, &
fgsl_sf_legendre_array_size
public :: fgsl_sf_conicalp_half, fgsl_sf_conicalp_half_e, fgsl_sf_conicalp_mhalf, &
fgsl_sf_conicalp_mhalf_e, fgsl_sf_conicalp_0, fgsl_sf_conicalp_0_e, &
fgsl_sf_conicalp_1, fgsl_sf_conicalp_1_e, fgsl_sf_conicalp_sph_reg, &
fgsl_sf_conicalp_sph_reg_e, fgsl_sf_conicalp_cyl_reg, fgsl_sf_conicalp_cyl_reg_e, &
fgsl_sf_legendre_h3d_0, fgsl_sf_legendre_h3d_0_e, fgsl_sf_legendre_h3d_1, &
fgsl_sf_legendre_h3d_1_e, fgsl_sf_legendre_h3d, fgsl_sf_legendre_h3d_e, &
fgsl_sf_legendre_h3d_array
public :: fgsl_sf_log, fgsl_sf_log_e, fgsl_sf_log_abs, fgsl_sf_log_abs_e, &
fgsl_sf_complex_log_e, fgsl_sf_log_1plusx, fgsl_sf_log_1plusx_e, &
fgsl_sf_log_1plusx_mx, fgsl_sf_log_1plusx_mx_e, fgsl_sf_psi_int, &
fgsl_sf_psi_1_int, fgsl_sf_psi_1_int_e, fgsl_sf_psi_1, fgsl_sf_psi_1_e,&
fgsl_sf_psi_n, fgsl_sf_psi_n_e,&
fgsl_sf_psi_int_e, fgsl_sf_psi, fgsl_sf_psi_e, fgsl_sf_psi_1piy, &
fgsl_sf_psi_1piy_e, fgsl_sf_synchrotron_1, fgsl_sf_synchrotron_1_e, &
fgsl_sf_synchrotron_2, fgsl_sf_synchrotron_2_e, fgsl_sf_transport_2, &
fgsl_sf_transport_2_e, fgsl_sf_transport_3, fgsl_sf_transport_3_e, &
fgsl_sf_transport_4, fgsl_sf_transport_4_e, fgsl_sf_transport_5, &
fgsl_sf_transport_5_e, fgsl_sf_hypot, fgsl_sf_hypot_e, fgsl_sf_sinc, &
fgsl_sf_sinc_e, fgsl_sf_complex_sin_e, fgsl_sf_complex_cos_e, &
fgsl_sf_complex_logsin_e, fgsl_sf_lnsinh, fgsl_sf_lnsinh_e, &
fgsl_sf_lncosh, fgsl_sf_lncosh_e, fgsl_sf_polar_to_rect, &
fgsl_sf_rect_to_polar, fgsl_sf_angle_restrict_symm, fgsl_sf_angle_restrict_symm_e, &
fgsl_sf_angle_restrict_pos, fgsl_sf_angle_restrict_pos_e, fgsl_sf_sin_err_e, &
fgsl_sf_cos_err_e, fgsl_sf_zeta_int, fgsl_sf_zeta_int_e, &
fgsl_sf_eta_int, fgsl_sf_eta_int_e, fgsl_sf_zeta, &
fgsl_sf_zeta_e, fgsl_sf_zetam1_int, fgsl_sf_zetam1_int_e, fgsl_sf_zetam1, &
fgsl_sf_zetam1_e, fgsl_sf_hzeta, fgsl_sf_hzeta_e, fgsl_sf_eta, fgsl_sf_eta_e
! array processing
public :: fgsl_vector_init, fgsl_vector_align, fgsl_vector_free
public :: fgsl_matrix_init, fgsl_matrix_align, fgsl_matrix_free
! interpolation
public :: fgsl_interp_alloc, fgsl_interp_init, &
fgsl_interp_free, fgsl_interp_eval, fgsl_interp_eval_e, &
fgsl_interp_eval_deriv, fgsl_interp_eval_deriv_e, &
fgsl_interp_eval_deriv2, fgsl_interp_eval_deriv2_e, &
fgsl_interp_eval_integ, fgsl_interp_eval_integ_e, &
fgsl_interp_min_size,fgsl_interp_name
public :: fgsl_spline_alloc, fgsl_spline_init, fgsl_spline_free, &
fgsl_spline_name, fgsl_spline_min_size, fgsl_spline_eval, &
fgsl_spline_eval_e, fgsl_spline_eval_deriv, fgsl_spline_eval_deriv_e, &
fgsl_spline_eval_deriv2, fgsl_spline_eval_deriv2_e, &
fgsl_spline_eval_integ, fgsl_spline_eval_integ_e
public :: fgsl_interp_accel_alloc, fgsl_interp_accel_free, &
fgsl_interp_accel_find, fgsl_interp_bsearch
! permutations and combinations
public :: fgsl_permutation_alloc, fgsl_permutation_calloc, fgsl_permutation_init, &
fgsl_permutation_free, fgsl_permutation_memcpy, fgsl_permutation_get, &
fgsl_permutation_swap, fgsl_permutation_size, fgsl_permutation_data, &
fgsl_permutation_valid, fgsl_permutation_reverse, fgsl_permutation_inverse, &
fgsl_permutation_next, fgsl_permutation_prev, &
fgsl_permute, fgsl_permute_inverse, fgsl_permute_vector, &
fgsl_permute_vector_inverse, fgsl_permutation_mul, &
fgsl_permutation_linear_to_canonical, fgsl_permutation_canonical_to_linear, &
fgsl_permutation_inversions, fgsl_permutation_linear_cycles, &
fgsl_permutation_canonical_cycles, fgsl_permutation_fwrite, &
fgsl_permutation_fread, fgsl_permutation_fprintf, &
fgsl_permutation_fscanf, fgsl_combination_alloc, &
fgsl_combination_calloc, fgsl_combination_init_first, fgsl_combination_init_last, &
fgsl_combination_free, fgsl_combination_memcpy, fgsl_combination_get, &
fgsl_combination_n, fgsl_combination_k, fgsl_combination_data, &
fgsl_combination_valid, fgsl_combination_next, &
fgsl_combination_prev, fgsl_combination_fwrite, fgsl_combination_fread, &
fgsl_combination_fprintf, fgsl_combination_fscanf
! sorting
public :: fgsl_heapsort, fgsl_heapsort_index, fgsl_sort, fgsl_sort_index, &
fgsl_sort_smallest, fgsl_sort_smallest_index, &
fgsl_sort_largest, fgsl_sort_largest_index
! linear algebra
public :: fgsl_linalg_lu_decomp, fgsl_linalg_complex_lu_decomp, &
fgsl_linalg_lu_solve, fgsl_linalg_complex_lu_solve, &
fgsl_linalg_lu_svx, fgsl_linalg_complex_lu_svx, &
fgsl_linalg_lu_refine, fgsl_linalg_complex_lu_refine, &
fgsl_linalg_lu_invert, fgsl_linalg_complex_lu_invert, &
fgsl_linalg_lu_det, fgsl_linalg_complex_lu_det, &
fgsl_linalg_lu_lndet, fgsl_linalg_complex_lu_lndet, &
fgsl_linalg_lu_sgndet, fgsl_linalg_complex_lu_sgndet, &
fgsl_linalg_qr_decomp, fgsl_linalg_qr_solve, fgsl_linalg_qr_svx, &
fgsl_linalg_qr_lssolve, fgsl_linalg_qr_qtvec, fgsl_linalg_qr_qvec, &
fgsl_linalg_qr_qtmat, fgsl_linalg_qr_rsolve, fgsl_linalg_qr_rsvx, &
fgsl_linalg_qr_unpack, fgsl_linalg_qr_qrsolve, fgsl_linalg_qr_update, &
fgsl_linalg_r_solve, fgsl_linalg_r_svx, fgsl_linalg_qrpt_decomp, &
fgsl_linalg_qrpt_decomp2, fgsl_linalg_qrpt_solve, fgsl_linalg_qrpt_svx, &
fgsl_linalg_qrpt_qrsolve, fgsl_linalg_qrpt_update, &
fgsl_linalg_qrpt_rsolve, fgsl_linalg_qrpt_rsvx, &
fgsl_linalg_sv_decomp, fgsl_linalg_sv_decomp_mod, &
fgsl_linalg_sv_decomp_jacobi, fgsl_linalg_sv_solve, &
fgsl_linalg_cholesky_decomp, fgsl_linalg_complex_cholesky_decomp, &
fgsl_linalg_cholesky_solve, fgsl_linalg_complex_cholesky_solve, &
fgsl_linalg_cholesky_svx, fgsl_linalg_complex_cholesky_svx, &
fgsl_linalg_cholesky_invert, &
fgsl_linalg_symmtd_decomp, fgsl_linalg_symmtd_unpack, &
fgsl_linalg_symmtd_unpack_t, fgsl_linalg_hermtd_decomp, &
fgsl_linalg_hermtd_unpack, fgsl_linalg_hermtd_unpack_t, &
fgsl_linalg_hessenberg_decomp, fgsl_linalg_hessenberg_unpack, &
fgsl_linalg_hessenberg_unpack_accum, fgsl_linalg_hessenberg_set_zero, &
fgsl_linalg_hesstri_decomp, &
fgsl_linalg_bidiag_decomp, fgsl_linalg_bidiag_unpack, &
fgsl_linalg_bidiag_unpack2, fgsl_linalg_bidiag_unpack_b, &
fgsl_linalg_householder_transform, &
fgsl_linalg_complex_householder_transform, &
fgsl_linalg_householder_hm, fgsl_linalg_complex_householder_hm, &
fgsl_linalg_householder_mh, fgsl_linalg_complex_householder_mh, &
fgsl_linalg_householder_hv, fgsl_linalg_complex_householder_hv, &
fgsl_linalg_hh_solve, fgsl_linalg_hh_svx, fgsl_linalg_solve_tridiag, &
fgsl_linalg_solve_symm_tridiag, fgsl_linalg_solve_cyc_tridiag, &
fgsl_linalg_solve_symm_cyc_tridiag, fgsl_linalg_balance_matrix
! eigensystems
public :: fgsl_eigen_symm_alloc, fgsl_eigen_symm_free, fgsl_eigen_symm, &
fgsl_eigen_symmv_alloc, fgsl_eigen_symmv_free, fgsl_eigen_symmv, &
fgsl_eigen_herm_alloc, fgsl_eigen_herm_free, fgsl_eigen_herm, &
fgsl_eigen_hermv_alloc, fgsl_eigen_hermv_free, fgsl_eigen_hermv, &
fgsl_eigen_nonsymm_alloc, fgsl_eigen_nonsymm_free, fgsl_eigen_nonsymm, &
fgsl_eigen_nonsymmv_alloc, fgsl_eigen_nonsymmv_free, fgsl_eigen_nonsymmv, &
fgsl_eigen_nonsymm_params, fgsl_eigen_nonsymm_z, fgsl_eigen_nonsymmv_z, &
fgsl_eigen_gensymm_alloc, fgsl_eigen_gensymm_free, fgsl_eigen_gensymm, &
fgsl_eigen_gensymmv_alloc, fgsl_eigen_gensymmv_free, fgsl_eigen_gensymmv, &
fgsl_eigen_genherm_alloc, fgsl_eigen_genherm_free, fgsl_eigen_genherm, &
fgsl_eigen_genhermv_alloc, fgsl_eigen_genhermv_free, fgsl_eigen_genhermv, &
fgsl_eigen_gen_alloc, fgsl_eigen_gen_free, fgsl_eigen_gen, &
fgsl_eigen_genv_alloc, fgsl_eigen_genv_free, fgsl_eigen_genv, &
fgsl_eigen_gen_params, fgsl_eigen_gen_qz, fgsl_eigen_genv_qz, &
fgsl_eigen_symmv_sort, fgsl_eigen_hermv_sort, &
fgsl_eigen_nonsymmv_sort, fgsl_eigen_gensymmv_sort, &
fgsl_eigen_genhermv_sort, fgsl_eigen_genv_sort
! FFT
public :: fgsl_fft_complex_radix2_forward, fgsl_fft_complex_radix2_transform, &
fgsl_fft_complex_radix2_backward, fgsl_fft_complex_radix2_inverse, &
fgsl_fft_complex_radix2_dif_forward, fgsl_fft_complex_radix2_dif_transform, &
fgsl_fft_complex_radix2_dif_backward, fgsl_fft_complex_radix2_dif_inverse, &
fgsl_fft_complex_wavetable_alloc, fgsl_fft_complex_wavetable_free, &
fgsl_fft_complex_workspace_alloc, fgsl_fft_complex_workspace_free, &
fgsl_fft_complex_forward, fgsl_fft_complex_transform, &
fgsl_fft_complex_backward, fgsl_fft_complex_inverse, &
fgsl_fft_real_radix2_transform, fgsl_fft_halfcomplex_radix2_inverse, &
fgsl_fft_halfcomplex_radix2_backward, fgsl_fft_real_wavetable_alloc, &
fgsl_fft_real_wavetable_free, fgsl_fft_halfcomplex_wavetable_alloc, &
fgsl_fft_halfcomplex_wavetable_free, fgsl_fft_real_transform, &
fgsl_fft_halfcomplex_transform, fgsl_fft_real_unpack, &
fgsl_fft_halfcomplex_unpack
! numerical integration
public :: fgsl_integration_qng, fgsl_integration_workspace_alloc, &
fgsl_integration_workspace_free, &
fgsl_integration_qag, fgsl_integration_qagi, fgsl_integration_qags, &
fgsl_integration_qagp, fgsl_integration_qagiu, fgsl_integration_qagil, &
fgsl_integration_qawc, fgsl_integration_qaws_table_alloc, &
fgsl_integration_qaws_table_set, fgsl_integration_qaws_table_free, &
fgsl_integration_qaws, &
fgsl_integration_qawo_table_alloc, fgsl_integration_qawo_table_set, &
fgsl_integration_qawo_table_set_length, fgsl_integration_qawo, &
fgsl_integration_qawo_table_free, fgsl_integration_qawf
! random numbers, quasi-random numbers, distribution functions
public :: fgsl_rng_alloc, fgsl_rng_set, fgsl_rng_free, fgsl_rng_get, fgsl_rng_uniform, &
fgsl_rng_uniform_pos, fgsl_rng_uniform_int, fgsl_rng_name, fgsl_rng_max, &
fgsl_rng_min, fgsl_rng_env_setup, fgsl_rng_memcpy, fgsl_rng_clone, &
fgsl_rng_fwrite, fgsl_rng_fread, fgsl_qrng_alloc, &
fgsl_qrng_free, fgsl_qrng_init, fgsl_qrng_get, fgsl_qrng_name, fgsl_qrng_memcpy, &
fgsl_qrng_clone
public :: fgsl_ran_gaussian, fgsl_ran_gaussian_pdf, fgsl_ran_gaussian_ziggurat, &
fgsl_ran_gaussian_ratio_method, fgsl_ran_ugaussian, fgsl_ran_ugaussian_pdf, &
fgsl_ran_ugaussian_ratio_method, fgsl_cdf_gaussian_p, fgsl_cdf_gaussian_q, &
fgsl_cdf_gaussian_pinv, fgsl_cdf_gaussian_qinv, fgsl_cdf_ugaussian_p, &
fgsl_cdf_ugaussian_q, fgsl_cdf_ugaussian_pinv, fgsl_cdf_ugaussian_qinv, &
fgsl_ran_gaussian_tail, fgsl_ran_gaussian_tail_pdf, fgsl_ran_ugaussian_tail, &
fgsl_ran_ugaussian_tail_pdf, fgsl_ran_bivariate_gaussian, fgsl_ran_exponential, &
fgsl_ran_exponential_pdf, fgsl_cdf_exponential_p, fgsl_cdf_exponential_q, &
fgsl_cdf_exponential_pinv, fgsl_cdf_exponential_qinv, fgsl_ran_laplace, &
fgsl_ran_laplace_pdf, fgsl_cdf_laplace_p, fgsl_cdf_laplace_q, &
fgsl_cdf_laplace_pinv, fgsl_cdf_laplace_qinv, fgsl_ran_exppow, &
fgsl_ran_exppow_pdf, fgsl_cdf_exppow_p, fgsl_cdf_exppow_q, fgsl_ran_cauchy, &
fgsl_ran_cauchy_pdf, fgsl_cdf_cauchy_p, fgsl_cdf_cauchy_q, fgsl_cdf_cauchy_pinv, &
fgsl_cdf_cauchy_qinv, fgsl_ran_rayleigh, fgsl_ran_rayleigh_pdf, &
fgsl_cdf_rayleigh_p, fgsl_cdf_rayleigh_q, fgsl_cdf_rayleigh_pinv, &
fgsl_cdf_rayleigh_qinv, fgsl_ran_rayleigh_tail, fgsl_ran_rayleigh_tail_pdf, &
fgsl_ran_landau, fgsl_ran_landau_pdf, fgsl_ran_levy, fgsl_ran_levy_skew
public :: fgsl_ran_gamma, fgsl_ran_gamma_mt, fgsl_ran_gamma_pdf, fgsl_cdf_gamma_p, &
fgsl_cdf_gamma_q, fgsl_cdf_gamma_pinv, fgsl_cdf_gamma_qinv, fgsl_ran_flat, &
fgsl_ran_flat_pdf, fgsl_cdf_flat_p, fgsl_cdf_flat_q, fgsl_cdf_flat_pinv, &
fgsl_cdf_flat_qinv, fgsl_ran_lognormal, fgsl_ran_lognormal_pdf, &
fgsl_cdf_lognormal_p, fgsl_cdf_lognormal_q, fgsl_cdf_lognormal_pinv, &
fgsl_cdf_lognormal_qinv, fgsl_ran_chisq, fgsl_ran_chisq_pdf, fgsl_cdf_chisq_p, &
fgsl_cdf_chisq_q, fgsl_cdf_chisq_pinv, fgsl_cdf_chisq_qinv, fgsl_ran_fdist, &
fgsl_ran_fdist_pdf, fgsl_cdf_fdist_p, fgsl_cdf_fdist_q, fgsl_cdf_fdist_pinv, &
fgsl_cdf_fdist_qinv, fgsl_ran_tdist, fgsl_ran_tdist_pdf, fgsl_cdf_tdist_p, &
fgsl_cdf_tdist_q, fgsl_cdf_tdist_pinv, fgsl_cdf_tdist_qinv, fgsl_ran_beta, &
fgsl_ran_beta_pdf, fgsl_cdf_beta_p, fgsl_cdf_beta_q, fgsl_cdf_beta_pinv, &
fgsl_cdf_beta_qinv, fgsl_ran_logistic, fgsl_ran_logistic_pdf, fgsl_cdf_logistic_p, &
fgsl_cdf_logistic_q, fgsl_cdf_logistic_pinv, fgsl_cdf_logistic_qinv, &
fgsl_ran_pareto, fgsl_ran_pareto_pdf, fgsl_cdf_pareto_p, fgsl_cdf_pareto_q, &
fgsl_cdf_pareto_pinv, fgsl_cdf_pareto_qinv, fgsl_ran_dir_2d, &
fgsl_ran_dir_2d_trig_method, fgsl_ran_dir_3d, fgsl_ran_dir_nd, &
fgsl_ran_weibull, fgsl_ran_weibull_pdf, fgsl_cdf_weibull_p, fgsl_cdf_weibull_q, &
fgsl_cdf_weibull_pinv, fgsl_cdf_weibull_qinv, fgsl_ran_gumbel1, &
fgsl_ran_gumbel1_pdf, fgsl_cdf_gumbel1_p, fgsl_cdf_gumbel1_q, &
fgsl_cdf_gumbel1_pinv, fgsl_cdf_gumbel1_qinv, fgsl_ran_gumbel2, &
fgsl_ran_gumbel2_pdf, fgsl_cdf_gumbel2_p, fgsl_cdf_gumbel2_q, &
fgsl_cdf_gumbel2_pinv, fgsl_cdf_gumbel2_qinv, fgsl_ran_dirichlet, &
fgsl_ran_dirichlet_pdf, fgsl_ran_dirichlet_lnpdf, fgsl_ran_discrete_preproc, &
fgsl_ran_discrete, fgsl_ran_discrete_pdf, fgsl_ran_poisson, fgsl_ran_poisson_pdf, &
fgsl_cdf_poisson_p, fgsl_cdf_poisson_q, fgsl_ran_bernoulli, &
fgsl_ran_bernoulli_pdf, fgsl_ran_binomial, fgsl_ran_binomial_pdf, &
fgsl_cdf_binomial_p, fgsl_cdf_binomial_q, fgsl_ran_multinomial, &
fgsl_ran_multinomial_pdf, fgsl_ran_multinomial_lnpdf, fgsl_ran_negative_binomial, &
fgsl_ran_negative_binomial_pdf, fgsl_cdf_negative_binomial_p, &
fgsl_cdf_negative_binomial_q, fgsl_ran_pascal, fgsl_ran_pascal_pdf, &
fgsl_cdf_pascal_p, fgsl_cdf_pascal_q, fgsl_ran_geometric, fgsl_ran_geometric_pdf, &
fgsl_cdf_geometric_p, fgsl_cdf_geometric_q, fgsl_ran_hypergeometric, &
fgsl_ran_hypergeometric_pdf, fgsl_cdf_hypergeometric_p, fgsl_cdf_hypergeometric_q, &
fgsl_ran_logarithmic, fgsl_ran_logarithmic_pdf, &
fgsl_ran_shuffle, fgsl_ran_choose, fgsl_ran_sample
! simulated annealing
public :: fgsl_siman_params_init, fgsl_siman_params_free, fgsl_siman_solve
! ordinary differential equations
public :: fgsl_odeiv_system_init, fgsl_odeiv_system_free, &
fgsl_odeiv_step_alloc, fgsl_odeiv_step_status, fgsl_odeiv_system_status, &
fgsl_odeiv_step_reset, fgsl_odeiv_step_free, fgsl_odeiv_step_name, &
fgsl_odeiv_step_order, fgsl_odeiv_step_apply, fgsl_odeiv_control_standard_new, &
fgsl_odeiv_control_y_new, fgsl_odeiv_control_yp_new, &
fgsl_odeiv_control_scaled_new, fgsl_odeiv_control_init, &
fgsl_odeiv_control_free, fgsl_odeiv_control_hadjust, &
fgsl_odeiv_control_name, fgsl_odeiv_evolve_alloc, fgsl_odeiv_evolve_apply, &
fgsl_odeiv_evolve_reset, fgsl_odeiv_evolve_free
! Monte Carlo
public :: fgsl_monte_function_init, fgsl_monte_function_free, &
fgsl_monte_plain_alloc, fgsl_monte_plain_init, &
fgsl_monte_plain_integrate, fgsl_monte_plain_free, &
fgsl_monte_miser_alloc, fgsl_monte_miser_init, &
fgsl_monte_miser_integrate, fgsl_monte_miser_free, &
fgsl_monte_vegas_alloc, fgsl_monte_vegas_init, &
fgsl_monte_vegas_integrate, fgsl_monte_vegas_free, &
fgsl_monte_vegas_chisq, fgsl_monte_vegas_runval, &
fgsl_monte_miser_setparams, fgsl_monte_miser_getparams, &
fgsl_monte_vegas_setparams, fgsl_monte_vegas_getparams
! Histograms
public :: fgsl_histogram_alloc, fgsl_histogram_set_ranges, &
fgsl_histogram_set_ranges_uniform, fgsl_histogram_free, &
fgsl_histogram_memcpy, fgsl_histogram_clone, fgsl_histogram_increment, &
fgsl_histogram_accumulate, fgsl_histogram_get, fgsl_histogram_get_range, &
fgsl_histogram_max, fgsl_histogram_min, fgsl_histogram_bins, &
fgsl_histogram_reset, fgsl_histogram_find, fgsl_histogram_max_val, &
fgsl_histogram_min_val, fgsl_histogram_max_bin, &
fgsl_histogram_min_bin, fgsl_histogram_mean, &
fgsl_histogram_sigma, fgsl_histogram_sum, fgsl_histogram_equal_bins_p, &
fgsl_histogram_add, fgsl_histogram_sub, fgsl_histogram_mul, &
fgsl_histogram_div, fgsl_histogram_scale, fgsl_histogram_shift, &
fgsl_histogram_fwrite, fgsl_histogram_fread, fgsl_histogram_fprintf, &
fgsl_histogram_fscanf, fgsl_histogram_pdf_alloc, fgsl_histogram_pdf_init, &
fgsl_histogram_pdf_free, fgsl_histogram_pdf_sample, &
fgsl_histogram2d_alloc, fgsl_histogram2d_set_ranges, &
fgsl_histogram2d_set_ranges_uniform, fgsl_histogram2d_free, &
fgsl_histogram2d_memcpy, fgsl_histogram2d_clone, fgsl_histogram2d_increment, &
fgsl_histogram2d_accumulate, fgsl_histogram2d_get, fgsl_histogram2d_get_xrange, &
fgsl_histogram2d_get_yrange, fgsl_histogram2d_xmax, &
fgsl_histogram2d_xmin, fgsl_histogram2d_ymax, &
fgsl_histogram2d_ymin, fgsl_histogram2d_nx, fgsl_histogram2d_ny, &
fgsl_histogram2d_reset, fgsl_histogram2d_find, fgsl_histogram2d_max_val, &
fgsl_histogram2d_min_val, fgsl_histogram2d_max_bin, &
fgsl_histogram2d_min_bin, fgsl_histogram2d_xmean, fgsl_histogram2d_ymean, &
fgsl_histogram2d_xsigma, fgsl_histogram2d_ysigma, fgsl_histogram2d_cov, &
fgsl_histogram2d_sum, fgsl_histogram2d_equal_bins_p, &
fgsl_histogram2d_add, fgsl_histogram2d_sub, fgsl_histogram2d_mul, &
fgsl_histogram2d_div, fgsl_histogram2d_scale, fgsl_histogram2d_shift, &
fgsl_histogram2d_fwrite, fgsl_histogram2d_fread, fgsl_histogram2d_fprintf, &
fgsl_histogram2d_fscanf, fgsl_histogram2d_pdf_alloc, fgsl_histogram2d_pdf_init, &
fgsl_histogram2d_pdf_free, fgsl_histogram2d_pdf_sample
! Ntuples
public :: fgsl_ntuple_create, fgsl_ntuple_open, fgsl_ntuple_write, &
fgsl_ntuple_bookdata, fgsl_ntuple_read, fgsl_ntuple_close, &
fgsl_ntuple_select_fn_init, fgsl_ntuple_value_fn_init, &
fgsl_ntuple_select_fn_free, fgsl_ntuple_value_fn_free, &
fgsl_ntuple_project, fgsl_ntuple_data, fgsl_ntuple_size
! Numerical derivatives
public :: fgsl_deriv_central, fgsl_deriv_forward, fgsl_deriv_backward
! Chebyshev approximations
public :: fgsl_cheb_alloc, fgsl_cheb_free, fgsl_cheb_init, fgsl_cheb_order, &
fgsl_cheb_size, fgsl_cheb_coeffs, fgsl_cheb_eval, &
fgsl_cheb_eval_err, fgsl_cheb_eval_n, fgsl_cheb_eval_n_err, fgsl_cheb_calc_deriv, &
fgsl_cheb_calc_integ
! Series acceleration
public :: fgsl_sum_levin_u_alloc, fgsl_sum_levin_u_free, fgsl_sum_levin_u_accel, &
fgsl_sum_levin_utrunc_alloc, fgsl_sum_levin_utrunc_free, fgsl_sum_levin_utrunc_accel
! Wavelet transforms
public :: fgsl_wavelet_alloc, fgsl_wavelet_name, fgsl_wavelet_free, &
fgsl_wavelet_workspace_alloc, fgsl_wavelet_workspace_free, fgsl_wavelet_transform, &
fgsl_wavelet_transform_forward, fgsl_wavelet_transform_inverse, &
fgsl_wavelet2d_transform, fgsl_wavelet2d_transform_forward, &
fgsl_wavelet2d_transform_inverse, fgsl_wavelet2d_nstransform, &
fgsl_wavelet2d_nstransform_forward, fgsl_wavelet2d_nstransform_inverse
! Hankel Transform
public :: fgsl_dht_alloc, fgsl_dht_init, fgsl_dht_new, fgsl_dht_free, &
fgsl_dht_apply, fgsl_dht_x_sample, fgsl_dht_k_sample
! One-dimensional root finding
public :: fgsl_root_fsolver_alloc, fgsl_root_fdfsolver_alloc, fgsl_root_fsolver_set, &
fgsl_root_fdfsolver_set, fgsl_root_fsolver_free, fgsl_root_fdfsolver_free, &
fgsl_root_fsolver_name, fgsl_root_fdfsolver_name, fgsl_root_fsolver_iterate, &
fgsl_root_fdfsolver_iterate, fgsl_root_fsolver_x_lower, fgsl_root_fsolver_x_upper, &
fgsl_root_test_interval, fgsl_root_test_delta, fgsl_root_test_residual, &
fgsl_root_fsolver_root, fgsl_root_fdfsolver_root
! One-dimensional minimization
public :: fgsl_min_fminimizer_alloc, fgsl_min_fminimizer_free, fgsl_min_fminimizer_set, &
fgsl_min_fminimizer_set_with_values, fgsl_min_fminimizer_iterate, fgsl_min_fminimizer_name, &
fgsl_min_fminimizer_x_minimum, fgsl_min_fminimizer_x_lower, fgsl_min_fminimizer_x_upper, &
fgsl_min_fminimizer_f_minimum, fgsl_min_fminimizer_f_lower, fgsl_min_fminimizer_f_upper, &
fgsl_min_test_interval
! Multi-root
public :: fgsl_multiroot_function_init, fgsl_multiroot_function_free, &
fgsl_multiroot_fsolver_alloc, fgsl_multiroot_fsolver_free, fgsl_multiroot_fsolver_name, &
fgsl_multiroot_fsolver_iterate, fgsl_multiroot_fsolver_root, fgsl_multiroot_fsolver_f, &
fgsl_multiroot_fsolver_dx, fgsl_multiroot_test_delta, fgsl_multiroot_test_residual, &
fgsl_multiroot_fsolver_set, fgsl_multiroot_fdfsolver_alloc, fgsl_multiroot_fdfsolver_free, &
fgsl_multiroot_fdfsolver_name, fgsl_multiroot_function_fdf_init, fgsl_multiroot_function_fdf_free, &
fgsl_multiroot_fdfsolver_iterate, fgsl_multiroot_fdfsolver_root, fgsl_multiroot_fdfsolver_f, &
fgsl_multiroot_fdfsolver_dx, fgsl_multiroot_fdfsolver_set
! Multi-minimization
public :: fgsl_multimin_function_init, fgsl_multimin_function_fdf_init, &
fgsl_multimin_function_free, fgsl_multimin_function_fdf_free, fgsl_multimin_fminimizer_alloc, &
fgsl_multimin_fdfminimizer_alloc, fgsl_multimin_fminimizer_free, fgsl_multimin_fdfminimizer_free, &
fgsl_multimin_fminimizer_set, fgsl_multimin_fdfminimizer_set, fgsl_multimin_fminimizer_name, &
fgsl_multimin_fdfminimizer_name, fgsl_multimin_fminimizer_iterate, &
fgsl_multimin_fdfminimizer_iterate, &
fgsl_multimin_fminimizer_minimum, fgsl_multimin_fdfminimizer_minimum, &
fgsl_multimin_fdfminimizer_gradient, fgsl_multimin_fminimizer_size, &
fgsl_multimin_fdfminimizer_restart, fgsl_multimin_fminimizer_x, &
fgsl_multimin_test_size, fgsl_multimin_fdfminimizer_x, fgsl_multimin_test_gradient
! Linear and nonlinear fitting
public :: fgsl_fit_linear, fgsl_fit_wlinear, fgsl_fit_linear_est, fgsl_fit_mul, &
fgsl_fit_wmul, fgsl_fit_mul_est, fgsl_multifit_linear_alloc, fgsl_multifit_linear_free, &
fgsl_multifit_linear, fgsl_multifit_linear_svd, fgsl_multifit_wlinear, &
fgsl_multifit_wlinear_svd, fgsl_multifit_linear_est, fgsl_multifit_linear_residuals
public :: fgsl_multifit_function_init, fgsl_multifit_function_fdf_init, &
fgsl_multifit_function_free, fgsl_multifit_function_fdf_free, fgsl_multifit_fsolver_alloc, &
fgsl_multifit_fdfsolver_alloc, fgsl_multifit_fsolver_free, fgsl_multifit_fdfsolver_free, &
fgsl_multifit_fsolver_set, fgsl_multifit_fdfsolver_set, fgsl_multifit_fsolver_name, &
fgsl_multifit_fdfsolver_name, fgsl_multifit_fsolver_iterate, fgsl_multifit_fdfsolver_iterate, &
fgsl_multifit_fsolver_position, fgsl_multifit_fdfsolver_position, &
fgsl_multifit_fdfsolver_dx, fgsl_multifit_fdfsolver_f, fgsl_multifit_fdfsolver_jac, &
fgsl_multifit_test_delta, fgsl_multifit_test_gradient, fgsl_multifit_gradient, &
fgsl_multifit_covar
! statistics
public :: fgsl_stats_mean, fgsl_stats_variance, fgsl_stats_variance_m, &
fgsl_stats_sd, fgsl_stats_sd_m, fgsl_stats_variance_with_fixed_mean, &
fgsl_stats_sd_with_fixed_mean, fgsl_stats_absdev, fgsl_stats_absdev_m, &
fgsl_stats_skew, fgsl_stats_skew_m_sd, fgsl_stats_kurtosis, &
fgsl_stats_kurtosis_m_sd, fgsl_stats_lag1_autocorrelation, fgsl_stats_lag1_autocorrelation_m, &
fgsl_stats_covariance, fgsl_stats_correlation, fgsl_stats_covariance_m, &
fgsl_stats_wmean, fgsl_stats_wvariance, fgsl_stats_wvariance_m, &
fgsl_stats_wsd, fgsl_stats_wsd_m, &
fgsl_stats_wvariance_with_fixed_mean, fgsl_stats_wsd_with_fixed_mean, &
fgsl_stats_wabsdev, fgsl_stats_wabsdev_m, fgsl_stats_wskew, &
fgsl_stats_wskew_m_sd, fgsl_stats_wkurtosis, fgsl_stats_wkurtosis_m_sd, fgsl_stats_max, &
fgsl_stats_min, fgsl_stats_minmax, fgsl_stats_max_index, &
fgsl_stats_min_index, fgsl_stats_minmax_index, fgsl_stats_median_from_sorted_data, &
fgsl_stats_quantile_from_sorted_data
! B-splines
public :: fgsl_bspline_alloc, fgsl_bspline_free, fgsl_bspline_knots, &
fgsl_bspline_knots_uniform, fgsl_bspline_eval, &
fgsl_bspline_eval_nonzero, fgsl_bspline_ncoeffs, &
fgsl_bspline_deriv_alloc, fgsl_bspline_deriv_free, &
fgsl_bspline_deriv_eval, fgsl_bspline_deriv_eval_nonzero, &
fgsl_bspline_greville_abscissa

! IEEE
public :: fgsl_ieee_fprintf, fgsl_ieee_printf, fgsl_ieee_env_setup
!
!
! Kind and length parameters are default integer
!
integer, parameter, public :: fgsl_double = c_double
integer, parameter, public :: fgsl_double_complex = c_double_complex
! integer, parameter, public :: fgsl_extended = selected_real_kind(18)
integer, parameter, public :: fgsl_extended = selected_real_kind(13)
! FIXME - c_long_double unsupported, selected_real_kind(30) unsupported in g95
integer, parameter, public :: fgsl_float = c_float
integer, parameter, public :: fgsl_int = c_int
integer, parameter, public :: fgsl_long = c_long
integer, parameter, public :: fgsl_size_t = c_size_t
integer, parameter, public :: fgsl_char = c_char
integer, parameter, public :: fgsl_strmax = 128
integer, parameter, public :: fgsl_pathmax = 2048
!
! Version strings
!
character(kind=fgsl_char, len=*), public, parameter :: fgsl_version='0.9.2'
character(kind=fgsl_char, len=*), public, parameter :: fgsl_gslbase='1.13'
!
! Error codes
!
integer(fgsl_int), parameter, public :: fgsl_success = 0
integer(fgsl_int), parameter, public :: fgsl_failure = -1
integer(fgsl_int), parameter, public :: fgsl_continue = -2 ! iteration has not converged
integer(fgsl_int), parameter, public :: fgsl_edom = 1 ! input domain error, e.g. sqrt(-1)
integer(fgsl_int), parameter, public :: fgsl_erange = 2 ! output range error, e.g. exp(1e100)
integer(fgsl_int), parameter, public :: fgsl_efault = 3 ! invalid pointer
integer(fgsl_int), parameter, public :: fgsl_einval = 4 ! invalid argument supplied by user
integer(fgsl_int), parameter, public :: fgsl_efactor = 6 ! generic failure
integer(fgsl_int), parameter, public :: fgsl_esanity = 7 ! sanity check failed - shouldn't happen
integer(fgsl_int), parameter, public :: fgsl_enomem = 8 ! malloc failed
integer(fgsl_int), parameter, public :: fgsl_ebadfunc = 9 ! problem with user-supplied function
integer(fgsl_int), parameter, public :: fgsl_erunaway = 10 ! iterative process is out of control
integer(fgsl_int), parameter, public :: fgsl_emaxiter = 11 ! exceeded max number of iterations
integer(fgsl_int), parameter, public :: fgsl_ezerodiv = 12 ! tried to divide by zero
integer(fgsl_int), parameter, public :: fgsl_ebadtol = 13 ! user specified an invalid tolerance
integer(fgsl_int), parameter, public :: fgsl_etol = 14 ! failed to reach the specified tolerance
integer(fgsl_int), parameter, public :: fgsl_eundrflw = 15 ! underflow
integer(fgsl_int), parameter, public :: fgsl_eovrflw = 16 ! overflow
integer(fgsl_int), parameter, public :: fgsl_eloss = 17 ! loss of accuracy
integer(fgsl_int), parameter, public :: fgsl_eround = 18 ! failed because of roundoff error
integer(fgsl_int), parameter, public :: fgsl_ebadlen = 19 ! matrix, vector lengths are not conformant
integer(fgsl_int), parameter, public :: fgsl_enotsqr = 20 ! matrix not square
integer(fgsl_int), parameter, public :: fgsl_esing = 21 ! apparent singularity detected
integer(fgsl_int), parameter, public :: fgsl_ediverge = 22 ! integral or series is divergent
integer(fgsl_int), parameter, public :: fgsl_eunsup = 23 ! no hw support for requested feature
integer(fgsl_int), parameter, public :: fgsl_eunimpl = 24 ! requested feature not (yet) implemented
integer(fgsl_int), parameter, public :: fgsl_ecache = 25 ! cache limit exceeded
integer(fgsl_int), parameter, public :: fgsl_etable = 26 ! table limit exceeded
integer(fgsl_int), parameter, public :: fgsl_enoprog = 27 ! iteration: no progress towards solution
integer(fgsl_int), parameter, public :: fgsl_enoprogj = 28 ! jacobian evals not improving the solution
integer(fgsl_int), parameter, public :: fgsl_etolf = 29 ! can't reach specified tolerance in F
integer(fgsl_int), parameter, public :: fgsl_etolx = 30 ! can't reach specified tolerance in X
integer(fgsl_int), parameter, public :: fgsl_etolg = 31 ! can't reach specified tolerance in gradient
integer(fgsl_int), parameter, public :: fgsl_eof = 32 ! end of file
!
! mathematical constants from gsl_math.h
!
real(fgsl_extended), parameter, public :: m_e = 2.71828182845904523536028747135_fgsl_extended
real(fgsl_extended), parameter, public :: m_log2e = 1.44269504088896340735992468100_fgsl_extended
real(fgsl_extended), parameter, public :: m_log10e = 0.43429448190325182765112891892_fgsl_extended
real(fgsl_extended), parameter, public :: m_sqrt2 = 1.41421356237309504880168872421_fgsl_extended
real(fgsl_extended), parameter, public :: m_sqrt1_2 = 0.70710678118654752440084436210_fgsl_extended
real(fgsl_extended), parameter, public :: m_sqrt3 = 1.73205080756887729352744634151_fgsl_extended
real(fgsl_extended), parameter, public :: m_pi = 3.14159265358979323846264338328_fgsl_extended
real(fgsl_extended), parameter, public :: m_pi_2 = 1.57079632679489661923132169164_fgsl_extended
real(fgsl_extended), parameter, public :: m_pi_4 = 0.78539816339744830961566084582_fgsl_extended
real(fgsl_extended), parameter, public :: m_sqrtpi = 1.77245385090551602729816748334_fgsl_extended
real(fgsl_extended), parameter, public :: m_2_sqrtpi = 1.12837916709551257389615890312_fgsl_extended
real(fgsl_extended), parameter, public :: m_1_pi = 0.31830988618379067153776752675_fgsl_extended
real(fgsl_extended), parameter, public :: m_2_pi = 0.63661977236758134307553505349_fgsl_extended
real(fgsl_extended), parameter, public :: m_ln10 = 2.30258509299404568401799145468_fgsl_extended
real(fgsl_extended), parameter, public :: m_ln2 = 0.69314718055994530941723212146_fgsl_extended
real(fgsl_extended), parameter, public :: m_lnpi = 1.14472988584940017414342735135_fgsl_extended
real(fgsl_extended), parameter, public :: m_euler = 0.57721566490153286060651209008_fgsl_extended
! the following provokes warnings from g95 ... may need to change if refused by other compilers
! real(fgsl_double), parameter, public :: fgsl_posinf = 1.0_fgsl_double / 0.0_fgsl_double
! real(fgsl_double), parameter, public :: fgsl_neginf = -1.0_fgsl_double / 0.0_fgsl_double
! real(fgsl_double), parameter, public :: fgsl_nan = 0.0_fgsl_double / 0.0_fgsl_double
! probably should throw this out - use IEEE_VALUE intrinsic if these are needed.
!
! Numerical constants
!
real(fgsl_double), parameter, public :: fgsl_const_num_fine_structure = 7.297352533E-3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_avogadro = 6.02214199E23_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_yotta = 1e24_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_zetta = 1e21_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_exa = 1e18_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_peta = 1e15_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_tera = 1e12_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_giga = 1e9_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_mega = 1e6_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_kilo = 1e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_milli = 1e-3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_micro = 1e-6_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_nano = 1e-9_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_pico = 1e-12_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_femto = 1e-15_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_atto = 1e-18_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_zepto = 1e-21_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_num_yocto = 1e-24_fgsl_double
!
! MKSA physical units
!
real(fgsl_double), parameter, public :: fgsl_const_mksa_speed_of_light = 2.99792458e8_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_gravitational_constant = 6.673e-11_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_plancks_constant_h = 6.62606896e-34_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_plancks_constant_hbar = 1.05457162825e-34_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_astronomical_unit = 1.49597870691e11_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_light_year = 9.46053620707e15_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_parsec = 3.08567758135e16_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_grav_accel = 9.80665e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_electron_volt = 1.602176487e-19_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_mass_electron = 9.10938188e-31_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_mass_muon = 1.88353109e-28_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_mass_proton = 1.67262158e-27_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_mass_neutron = 1.67492716e-27_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_rydberg = 2.17987196968e-18_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_boltzmann = 1.3806504e-23_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_bohr_magneton = 9.27400899e-24_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_nuclear_magneton = 5.05078317e-27_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_electron_magnetic_moment = 9.28476362e-24_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_proton_magnetic_moment = 1.410606633e-26_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_molar_gas = 8.314472e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_standard_gas_volume = 2.2710981e-2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_minute = 6e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_hour = 3.6e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_day = 8.64e4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_week = 6.048e5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_inch = 2.54e-2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_foot = 3.048e-1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_yard = 9.144e-1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_mile = 1.609344e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_nautical_mile = 1.852e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_fathom = 1.8288e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_mil = 2.54e-5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_point = 3.52777777778e-4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_texpoint = 3.51459803515e-4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_micron = 1e-6_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_angstrom = 1e-10_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_hectare = 1e4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_acre = 4.04685642241e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_barn = 1e-28_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_liter = 1e-3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_us_gallon = 3.78541178402e-3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_quart = 9.46352946004e-4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_pint = 4.73176473002e-4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_cup = 2.36588236501e-4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_fluid_ounce = 2.95735295626e-5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_tablespoon = 1.47867647813e-5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_teaspoon = 4.92892159375e-6_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_canadian_gallon = 4.54609e-3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_uk_gallon = 4.546092e-3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_miles_per_hour = 4.4704e-1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_kilometers_per_hour = 2.77777777778e-1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_knot = 5.14444444444e-1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_pound_mass = 4.5359237e-1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_ounce_mass = 2.8349523125e-2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_ton = 9.0718474e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_metric_ton = 1e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_uk_ton = 1.0160469088e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_troy_ounce = 3.1103475e-2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_carat = 2e-4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_unified_atomic_mass = 1.660538782e-27_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_gram_force = 9.80665e-3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_pound_force = 4.44822161526e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_kilopound_force = 4.44822161526e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_poundal = 1.38255e-1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_calorie = 4.1868e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_btu = 1.05505585262e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_therm = 1.05506e8_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_horsepower = 7.457e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_bar = 1e5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_std_atmosphere = 1.01325e5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_torr = 1.33322368421e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_meter_of_mercury = 1.33322368421e5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_inch_of_mercury = 3.38638815789e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_inch_of_water = 2.490889e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_psi = 6.89475729317e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_poise = 1e-1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_stokes = 1e-4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_faraday = 9.64853429775e4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_electron_charge = 1.602176487e-19_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_gauss = 1e-4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_stilb = 1e4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_lumen = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_lux = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_phot = 1e4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_footcandle = 1.076e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_lambert = 1e4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_footlambert = 1.07639104e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_curie = 3.7e10_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_roentgen = 2.58e-4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_rad = 1e-2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_solar_mass = 1.98892e30_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_bohr_radius = 5.291772083e-11_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_newton = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_dyne = 1e-5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_joule = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_erg = 1e-7_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_stefan_boltzmann_constant = 5.67040047374e-8_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_thomson_cross_section = 6.65245893699e-29_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_vacuum_permittivity = 8.854187817e-12_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_vacuum_permeability = 1.25663706144e-6_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_mksa_debye = 3.33564095198e-30_fgsl_double
!
! CGSM physical constants
!
real(fgsl_double), parameter, public :: fgsl_const_cgsm_speed_of_light = 2.99792458e10_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_gravitational_constant = 6.673e-8_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_plancks_constant_h = 6.62606896e-27_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_plancks_constant_hbar = 1.05457162825e-27_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_astronomical_unit = 1.49597870691e13_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_light_year = 9.46053620707e17_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_parsec = 3.08567758135e18_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_grav_accel = 9.80665e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_electron_volt = 1.602176487e-12_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_mass_electron = 9.10938188e-28_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_mass_muon = 1.88353109e-25_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_mass_proton = 1.67262158e-24_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_mass_neutron = 1.67492716e-24_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_rydberg = 2.17987196968e-11_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_boltzmann = 1.3806504e-16_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_bohr_magneton = 9.27400899e-21_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_nuclear_magneton = 5.05078317e-24_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_electron_magnetic_moment = 9.28476362e-21_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_proton_magnetic_moment = 1.410606633e-23_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_molar_gas = 8.314472e7_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_standard_gas_volume = 2.2710981e4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_minute = 6e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_hour = 3.6e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_day = 8.64e4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_week = 6.048e5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_inch = 2.54e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_foot = 3.048e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_yard = 9.144e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_mile = 1.609344e5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_nautical_mile = 1.852e5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_fathom = 1.8288e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_mil = 2.54e-3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_point = 3.52777777778e-2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_texpoint = 3.51459803515e-2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_micron = 1e-4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_angstrom = 1e-8_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_hectare = 1e8_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_acre = 4.04685642241e7_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_barn = 1e-24_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_liter = 1e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_us_gallon = 3.78541178402e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_quart = 9.46352946004e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_pint = 4.73176473002e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_cup = 2.36588236501e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_fluid_ounce = 2.95735295626e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_tablespoon = 1.47867647813e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_teaspoon = 4.92892159375e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_canadian_gallon = 4.54609e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_uk_gallon = 4.546092e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_miles_per_hour = 4.4704e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_kilometers_per_hour = 2.77777777778e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_knot = 5.14444444444e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_pound_mass = 4.5359237e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_ounce_mass = 2.8349523125e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_ton = 9.0718474e5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_metric_ton = 1e6_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_uk_ton = 1.0160469088e6_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_troy_ounce = 3.1103475e1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_carat = 2e-1_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_unified_atomic_mass = 1.660538782e-24_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_gram_force = 9.80665e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_pound_force = 4.44822161526e5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_kilopound_force = 4.44822161526e8_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_poundal = 1.38255e4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_calorie = 4.1868e7_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_btu = 1.05505585262e10_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_therm = 1.05506e15_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_horsepower = 7.457e9_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_bar = 1e6_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_std_atmosphere = 1.01325e6_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_torr = 1.33322368421e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_meter_of_mercury = 1.33322368421e6_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_inch_of_mercury = 3.38638815789e4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_inch_of_water = 2.490889e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_psi = 6.89475729317e4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_poise = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_stokes = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_faraday = 9.64853429775e3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_electron_charge = 1.602176487e-20_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_gauss = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_stilb = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_lumen = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_lux = 1e-4_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_phot = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_footcandle = 1.076e-3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_lambert = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_footlambert = 1.07639104e-3_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_curie = 3.7e10_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_roentgen = 2.58e-8_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_rad = 1e2_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_solar_mass = 1.98892e33_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_bohr_radius = 5.291772083e-9_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_newton = 1e5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_dyne = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_joule = 1e7_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_erg = 1e0_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_stefan_boltzmann_constant = 5.67040047374e-5_fgsl_double
real(fgsl_double), parameter, public :: fgsl_const_cgsm_thomson_cross_section = 6.65245893699e-25_fgsl_double
!
! Types : Error treatment
!
type, public :: fgsl_error_handler_t
private
type(c_funptr) :: gsl_error_handler_t = c_null_funptr
end type fgsl_error_handler_t
!
! Types: I/O Add-ons
!
type, public :: fgsl_file
private
type(c_ptr) :: gsl_file = c_null_ptr
end type fgsl_file
!
! Types: Mathematical Functions
!
type, public :: fgsl_function
private
type(c_ptr) :: gsl_function = c_null_ptr
end type fgsl_function
type, public :: fgsl_function_fdf
private
type(c_ptr) :: gsl_function_fdf = c_null_ptr
end type fgsl_function_fdf
!
! Types: Polynomial
!
! FIXME ifort refuses = overload if not public
type, public, bind(c) :: gsl_complex
real(c_double) :: dat(2)
end type
type, public :: fgsl_poly_complex_workspace
private
type(c_ptr) :: gsl_poly_complex_workspace
end type fgsl_poly_complex_workspace
!
! Types: Special Functions
!
type, public :: fgsl_sf_result
real(fgsl_double) :: val, err
end type fgsl_sf_result
! FIXME ifort refuses = overload if not public
type, public, bind(c) :: gsl_sf_result
real(c_double) :: val, err
end type
type, public :: fgsl_sf_result_e10
real(fgsl_double) :: val, err
integer(fgsl_int) :: e10
end type fgsl_sf_result_e10
! FIXME ifort refuses = overload if not public
type, public, bind(c) :: gsl_sf_result_e10
real(c_double) :: val, err
integer(c_int) :: e10
end type
type, public :: fgsl_mode_t
private
integer(c_int) :: gsl_mode = 0
end type fgsl_mode_t
type(fgsl_mode_t), parameter, public :: &
fgsl_prec_double = fgsl_mode_t(0), &
fgsl_prec_single = fgsl_mode_t(1), &
fgsl_prec_approx = fgsl_mode_t(2)
!
! Types : Array support
!
type, public :: fgsl_vector
private
type(c_ptr) :: gsl_vector = c_null_ptr
end type fgsl_vector
type, public :: fgsl_matrix
private
type(c_ptr) :: gsl_matrix = c_null_ptr
end type fgsl_matrix
type, public :: fgsl_vector_complex
private
type(c_ptr) :: gsl_vector_complex = c_null_ptr
end type fgsl_vector_complex
type, public :: fgsl_matrix_complex
private
type(c_ptr) :: gsl_matrix_complex = c_null_ptr
end type fgsl_matrix_complex
!
! Types : Interpolation
!
! integer, parameter :: interp_maxnum = 6
type, public :: fgsl_interp_type
private
integer(fgsl_int) :: which = 0
end type fgsl_interp_type
type(fgsl_interp_type), parameter, public :: &
fgsl_interp_linear = fgsl_interp_type(1), &
fgsl_interp_polynomial = fgsl_interp_type(2), &
fgsl_interp_cspline = fgsl_interp_type(3), &
fgsl_interp_cspline_periodic = fgsl_interp_type(4), &
fgsl_interp_akima = fgsl_interp_type(5), &
fgsl_interp_akima_periodic = fgsl_interp_type(6)
type, public :: fgsl_interp
private
type(c_ptr) :: gsl_interp = c_null_ptr
end type fgsl_interp
type, public :: fgsl_interp_accel
private
type(c_ptr) :: gsl_interp_accel = c_null_ptr
end type fgsl_interp_accel
type, public :: fgsl_spline
private
type(c_ptr) :: gsl_spline = c_null_ptr
end type fgsl_spline
!
! Types: Permutations and Combinations
!
type, public :: fgsl_permutation
private
type(c_ptr) :: gsl_permutation = c_null_ptr
end type fgsl_permutation
type, public :: fgsl_combination
private
type(c_ptr) :: gsl_combination = c_null_ptr
end type fgsl_combination
!
! Types: Eigensystems
!
type, public :: fgsl_eigen_symm_workspace
private
type(c_ptr) :: gsl_eigen_symm_workspace = c_null_ptr
end type fgsl_eigen_symm_workspace
type, public :: fgsl_eigen_symmv_workspace
private
type(c_ptr) :: gsl_eigen_symmv_workspace = c_null_ptr
end type fgsl_eigen_symmv_workspace
type, public :: fgsl_eigen_herm_workspace
private
type(c_ptr) :: gsl_eigen_herm_workspace = c_null_ptr
end type fgsl_eigen_herm_workspace
type, public :: fgsl_eigen_hermv_workspace
private
type(c_ptr) :: gsl_eigen_hermv_workspace = c_null_ptr
end type fgsl_eigen_hermv_workspace
type, public :: fgsl_eigen_nonsymm_workspace
private
type(c_ptr) :: gsl_eigen_nonsymm_workspace = c_null_ptr
end type fgsl_eigen_nonsymm_workspace
type, public :: fgsl_eigen_nonsymmv_workspace
private
type(c_ptr) :: gsl_eigen_nonsymmv_workspace = c_null_ptr
end type fgsl_eigen_nonsymmv_workspace
type, public :: fgsl_eigen_gensymm_workspace
private
type(c_ptr) :: gsl_eigen_gensymm_workspace = c_null_ptr
end type fgsl_eigen_gensymm_workspace
type, public :: fgsl_eigen_gensymmv_workspace
private
type(c_ptr) :: gsl_eigen_gensymmv_workspace = c_null_ptr
end type fgsl_eigen_gensymmv_workspace
type, public :: fgsl_eigen_genherm_workspace
private
type(c_ptr) :: gsl_eigen_genherm_workspace = c_null_ptr
end type fgsl_eigen_genherm_workspace
type, public :: fgsl_eigen_genhermv_workspace
private
type(c_ptr) :: gsl_eigen_genhermv_workspace = c_null_ptr
end type fgsl_eigen_genhermv_workspace
type, public :: fgsl_eigen_gen_workspace
private
type(c_ptr) :: gsl_eigen_gen_workspace = c_null_ptr
end type fgsl_eigen_gen_workspace
type, public :: fgsl_eigen_genv_workspace
private
type(c_ptr) :: gsl_eigen_genv_workspace = c_null_ptr
end type fgsl_eigen_genv_workspace
integer(c_int), parameter, public :: fgsl_eigen_sort_val_asc = 0
integer(c_int), parameter, public :: fgsl_eigen_sort_val_desc = 1
integer(c_int), parameter, public :: fgsl_eigen_sort_abs_asc = 2
integer(c_int), parameter, public :: fgsl_eigen_sort_abs_desc = 3
!
! Types: FFT
!
type, public :: fgsl_fft_complex_wavetable
private
type(c_ptr) :: gsl_fft_complex_wavetable = c_null_ptr
end type fgsl_fft_complex_wavetable
type, public :: fgsl_fft_real_wavetable
private
type(c_ptr) :: gsl_fft_real_wavetable = c_null_ptr
end type fgsl_fft_real_wavetable
type, public :: fgsl_fft_halfcomplex_wavetable
private
type(c_ptr) :: gsl_fft_halfcomplex_wavetable = c_null_ptr
end type fgsl_fft_halfcomplex_wavetable
type, public :: fgsl_fft_complex_workspace
private
type(c_ptr) :: gsl_fft_complex_workspace = c_null_ptr
end type fgsl_fft_complex_workspace
type, public :: fgsl_fft_real_workspace
private
type(c_ptr) :: gsl_fft_real_workspace = c_null_ptr
end type fgsl_fft_real_workspace
!
! Types: Numerical Integration
!
type, public :: fgsl_integration_workspace
private
type(c_ptr) :: gsl_integration_workspace = c_null_ptr
logical :: status = .false.
end type fgsl_integration_workspace
integer(fgsl_int), parameter, public :: fgsl_integ_gauss15 = 1
integer(fgsl_int), parameter, public :: fgsl_integ_gauss21 = 2
integer(fgsl_int), parameter, public :: fgsl_integ_gauss31 = 3
integer(fgsl_int), parameter, public :: fgsl_integ_gauss41 = 4
integer(fgsl_int), parameter, public :: fgsl_integ_gauss51 = 5
integer(fgsl_int), parameter, public :: fgsl_integ_gauss61 = 6
type, public :: fgsl_integration_qaws_table
private
type(c_ptr) :: gsl_integration_qaws_table
logical :: status = .false.
end type fgsl_integration_qaws_table
type, public :: fgsl_integration_qawo_table
private
type(c_ptr) :: gsl_integration_qawo_table
logical :: status = .false.
end type fgsl_integration_qawo_table
integer(fgsl_int), parameter, public :: fgsl_integ_cosine = 0
integer(fgsl_int), parameter, public :: fgsl_integ_sine = 1
!
! Types: Random and Quasi-random numbers
!
type, public :: fgsl_rng
private
type(c_ptr) :: gsl_rng
end type fgsl_rng
type, public :: fgsl_rng_type
private
type(c_ptr) :: gsl_rng_type
integer(fgsl_int) :: type = 0
end type fgsl_rng_type
! Note: we need a dynamic component here, since
! fgsl_rng_default needs to change at run time.
! fgsl_rng_default will be set by fgsl_rng_env_setup,
! and the static objects will be all set at the
! first call of fgsl_rng_alloc
integer, parameter :: rngmax = 61
! check new g95 type(fgsl_rng_type) :: fgsl_rng_allgen(rngmax) = (/(fgsl_rng_type(c_null_ptr, i)), i=1, rngmax/)
! cannot have protected attribute since modified by fgsl_rng_alloc
type(fgsl_rng_type), public :: &
fgsl_rng_default = fgsl_rng_type(c_null_ptr, -1), &
fgsl_rng_borosh13 = fgsl_rng_type(c_null_ptr, 1), &
fgsl_rng_coveyou = fgsl_rng_type(c_null_ptr, 2), &
fgsl_rng_cmrg = fgsl_rng_type(c_null_ptr, 3), &
fgsl_rng_fishman18 = fgsl_rng_type(c_null_ptr, 4), &
fgsl_rng_fishman20 = fgsl_rng_type(c_null_ptr, 5), &
fgsl_rng_fishman2x = fgsl_rng_type(c_null_ptr, 6), &
fgsl_rng_gfsr4 = fgsl_rng_type(c_null_ptr, 7), &
fgsl_rng_knuthran = fgsl_rng_type(c_null_ptr, 8), &
fgsl_rng_knuthran2 = fgsl_rng_type(c_null_ptr, 9), &
fgsl_rng_lecuyer21 = fgsl_rng_type(c_null_ptr, 10), &
fgsl_rng_minstd = fgsl_rng_type(c_null_ptr, 11), &
fgsl_rng_mrg = fgsl_rng_type(c_null_ptr, 12), &
fgsl_rng_mt19937 = fgsl_rng_type(c_null_ptr, 13), &
fgsl_rng_mt19937_1999 = fgsl_rng_type(c_null_ptr, 14), &
fgsl_rng_mt19937_1998 = fgsl_rng_type(c_null_ptr, 15), &
fgsl_rng_r250 = fgsl_rng_type(c_null_ptr, 16), &
fgsl_rng_ran0 = fgsl_rng_type(c_null_ptr, 17), &
fgsl_rng_ran1 = fgsl_rng_type(c_null_ptr, 18), &
fgsl_rng_ran2 = fgsl_rng_type(c_null_ptr, 19), &
fgsl_rng_ran3 = fgsl_rng_type(c_null_ptr, 20), &
fgsl_rng_rand = fgsl_rng_type(c_null_ptr, 21), &
fgsl_rng_rand48 = fgsl_rng_type(c_null_ptr, 22), &
fgsl_rng_random128_bsd = fgsl_rng_type(c_null_ptr, 23), &
fgsl_rng_random128_glibc2 = fgsl_rng_type(c_null_ptr, 24), &
fgsl_rng_random128_libc5 = fgsl_rng_type(c_null_ptr, 25), &
fgsl_rng_random256_bsd = fgsl_rng_type(c_null_ptr, 26), &
fgsl_rng_random256_glibc2 = fgsl_rng_type(c_null_ptr, 27), &
fgsl_rng_random256_libc5 = fgsl_rng_type(c_null_ptr, 28), &
fgsl_rng_random32_bsd = fgsl_rng_type(c_null_ptr, 29), &
fgsl_rng_random32_glibc2 = fgsl_rng_type(c_null_ptr, 30), &
fgsl_rng_random32_libc5 = fgsl_rng_type(c_null_ptr, 31), &
fgsl_rng_random64_bsd = fgsl_rng_type(c_null_ptr, 32), &
fgsl_rng_random64_glibc2 = fgsl_rng_type(c_null_ptr, 33), &
fgsl_rng_random64_libc5 = fgsl_rng_type(c_null_ptr, 34), &
fgsl_rng_random8_bsd = fgsl_rng_type(c_null_ptr, 35)
type(fgsl_rng_type), public :: &
fgsl_rng_random8_glibc2 = fgsl_rng_type(c_null_ptr, 36), &
fgsl_rng_random8_libc5 = fgsl_rng_type(c_null_ptr, 37), &
fgsl_rng_random_bsd = fgsl_rng_type(c_null_ptr, 38), &
fgsl_rng_random_glibc2 = fgsl_rng_type(c_null_ptr, 39), &
fgsl_rng_random_libc5 = fgsl_rng_type(c_null_ptr, 40), &
fgsl_rng_randu = fgsl_rng_type(c_null_ptr, 41), &
fgsl_rng_ranf = fgsl_rng_type(c_null_ptr, 42), &
fgsl_rng_ranlux = fgsl_rng_type(c_null_ptr, 43), &
fgsl_rng_ranlux389 = fgsl_rng_type(c_null_ptr, 44), &
fgsl_rng_ranlxd1 = fgsl_rng_type(c_null_ptr, 45), &
fgsl_rng_ranlxd2 = fgsl_rng_type(c_null_ptr, 46), &
fgsl_rng_ranlxs0 = fgsl_rng_type(c_null_ptr, 47), &
fgsl_rng_ranlxs1 = fgsl_rng_type(c_null_ptr, 48), &
fgsl_rng_ranlxs2 = fgsl_rng_type(c_null_ptr, 49), &
fgsl_rng_ranmar = fgsl_rng_type(c_null_ptr, 50), &
fgsl_rng_slatec = fgsl_rng_type(c_null_ptr, 51), &
fgsl_rng_taus = fgsl_rng_type(c_null_ptr, 52), &
fgsl_rng_taus2 = fgsl_rng_type(c_null_ptr, 53), &
fgsl_rng_taus113 = fgsl_rng_type(c_null_ptr, 54), &
fgsl_rng_transputer = fgsl_rng_type(c_null_ptr, 55), &
fgsl_rng_tt800 = fgsl_rng_type(c_null_ptr, 56), &
fgsl_rng_uni = fgsl_rng_type(c_null_ptr, 57), &
fgsl_rng_uni32 = fgsl_rng_type(c_null_ptr, 58), &
fgsl_rng_vax = fgsl_rng_type(c_null_ptr, 59), &
fgsl_rng_waterman14 = fgsl_rng_type(c_null_ptr, 60), &
fgsl_rng_zuf = fgsl_rng_type(c_null_ptr, 61), &
fgsl_rng_knuthran2002 = fgsl_rng_type(c_null_ptr, 62)
integer(fgsl_long), public, bind(c, name='gsl_rng_default_seed') :: fgsl_rng_default_seed
type, public :: fgsl_qrng
private
type(c_ptr) :: gsl_qrng
end type fgsl_qrng
type, public :: fgsl_qrng_type
private
integer(fgsl_int) :: type = 0
end type fgsl_qrng_type
type(fgsl_qrng_type), parameter, public :: &
fgsl_qrng_niederreiter_2 = fgsl_qrng_type(1), &
fgsl_qrng_sobol = fgsl_qrng_type(2), &
fgsl_qrng_halton = fgsl_qrng_type(3), &
fgsl_qrng_reversehalton = fgsl_qrng_type(4)
type, public :: fgsl_ran_discrete_t
private
type(c_ptr) :: gsl_ran_discrete_t
end type fgsl_ran_discrete_t
!
! Types: Histograms
!
type, public :: fgsl_histogram
private
type(c_ptr) :: gsl_histogram = c_null_ptr
end type fgsl_histogram
type, public :: fgsl_histogram_pdf
private
type(c_ptr) :: gsl_histogram_pdf = c_null_ptr
end type fgsl_histogram_pdf
type, public :: fgsl_histogram2d
private
type(c_ptr) :: gsl_histogram2d = c_null_ptr
end type fgsl_histogram2d
type, public :: fgsl_histogram2d_pdf
private
type(c_ptr) :: gsl_histogram2d_pdf = c_null_ptr
end type fgsl_histogram2d_pdf
!
! Types: Ntuples
!
type, public :: fgsl_ntuple
private
type(c_ptr) :: gsl_ntuple = c_null_ptr
end type fgsl_ntuple
type, public :: fgsl_ntuple_select_fn
private
type(c_ptr) :: gsl_ntuple_select_fn = c_null_ptr
end type fgsl_ntuple_select_fn
type, public :: fgsl_ntuple_value_fn
private
type(c_ptr) :: gsl_ntuple_value_fn = c_null_ptr
end type fgsl_ntuple_value_fn
!
! Types: Monte Carlo integration
!
type, public :: fgsl_monte_function
private
type(c_ptr) :: gsl_monte_function = c_null_ptr
end type fgsl_monte_function
type, public :: fgsl_monte_plain_state
private
type(c_ptr) :: gsl_monte_plain_state = c_null_ptr
end type fgsl_monte_plain_state
type, public :: fgsl_monte_miser_state
private
type(c_ptr) :: gsl_monte_miser_state = c_null_ptr
end type fgsl_monte_miser_state
type, public :: fgsl_monte_vegas_state
private
type(c_ptr) :: gsl_monte_vegas_state = c_null_ptr
end type fgsl_monte_vegas_state
! NOTE: not all compilers support enum yet
integer(c_int), parameter, public :: fgsl_vegas_mode_importance = 1
integer(c_int), parameter, public :: fgsl_vegas_mode_importance_only = 0
integer(c_int), parameter, public :: fgsl_vegas_mode_stratified = -1
!
! Types: Simulated Annealing
!
type, bind(c) :: gsl_siman_params_t
integer(c_int) :: n_tries, iters_fixed_t
real(c_double) :: step_size, k, t_initial, mu_t, t_min
end type gsl_siman_params_t
type, public :: fgsl_siman_params_t
private
type(gsl_siman_params_t), pointer :: gsl_siman_params_t => null()
end type fgsl_siman_params_t
!
! Types: Ordinary Differential Equations
!
type, public :: fgsl_odeiv_system
private
type(c_ptr) :: gsl_odeiv_system = c_null_ptr
end type fgsl_odeiv_system
type, public :: fgsl_odeiv_step_type
private
integer(c_int) :: which = 0
end type fgsl_odeiv_step_type
type(fgsl_odeiv_step_type), parameter, public :: &
fgsl_odeiv_step_rk2 = fgsl_odeiv_step_type(1), &
fgsl_odeiv_step_rk4 = fgsl_odeiv_step_type(2), &
fgsl_odeiv_step_rkf45 = fgsl_odeiv_step_type(3), &
fgsl_odeiv_step_rkck = fgsl_odeiv_step_type(4), &
fgsl_odeiv_step_rk8pd = fgsl_odeiv_step_type(5), &
fgsl_odeiv_step_rk2imp = fgsl_odeiv_step_type(6), &
fgsl_odeiv_step_rk2simp = fgsl_odeiv_step_type(7), &
fgsl_odeiv_step_rk4imp = fgsl_odeiv_step_type(8), &
fgsl_odeiv_step_bsimp = fgsl_odeiv_step_type(9), &
fgsl_odeiv_step_gear1 = fgsl_odeiv_step_type(10), &
fgsl_odeiv_step_gear2 = fgsl_odeiv_step_type(11)
type, public :: fgsl_odeiv_step
type(c_ptr) :: gsl_odeiv_step = c_null_ptr
end type fgsl_odeiv_step
type, public :: fgsl_odeiv_control
type(c_ptr) :: gsl_odeiv_control = c_null_ptr
end type fgsl_odeiv_control
integer(fgsl_int), parameter, public :: fgsl_odeiv_hadj_inc = 1
integer(fgsl_int), parameter, public :: fgsl_odeiv_hadj_nil = 0
integer(fgsl_int), parameter, public :: fgsl_odeiv_hadj_dec = -1
type, public :: fgsl_odeiv_evolve
type(c_ptr) :: gsl_odeiv_evolve
end type fgsl_odeiv_evolve
!
! Types: Chebyshev approximation
!
type, public :: fgsl_cheb_series
private
type(c_ptr) :: gsl_cheb_series = c_null_ptr
end type fgsl_cheb_series
!
! Types: Series acceleration
!
type, public :: fgsl_sum_levin_u_workspace
private
type(c_ptr) :: gsl_sum_levin_u_workspace = c_null_ptr
end type fgsl_sum_levin_u_workspace
type, public :: fgsl_sum_levin_utrunc_workspace
private
type(c_ptr) :: gsl_sum_levin_utrunc_workspace = c_null_ptr
end type fgsl_sum_levin_utrunc_workspace
!
! Types: Wavelet transforms
!
type, public :: fgsl_wavelet
private
type(c_ptr) :: gsl_wavelet = c_null_ptr
end type fgsl_wavelet
type, public :: fgsl_wavelet_type
private
integer(c_int) :: which = 0
end type fgsl_wavelet_type
type(fgsl_wavelet_type), public, parameter :: &
fgsl_wavelet_daubechies = fgsl_wavelet_type(1), &
fgsl_wavelet_daubechies_centered = fgsl_wavelet_type(2), &
fgsl_wavelet_haar = fgsl_wavelet_type(3), &
fgsl_wavelet_haar_centered = fgsl_wavelet_type(4), &
fgsl_wavelet_bspline = fgsl_wavelet_type(5), &
fgsl_wavelet_bspline_centered = fgsl_wavelet_type(6)
type, public :: fgsl_wavelet_workspace
private
type(c_ptr) :: gsl_wavelet_workspace
end type fgsl_wavelet_workspace
!
! Types: Hankel transforms
!
type, public :: fgsl_dht
private
type(c_ptr) :: gsl_dht = c_null_ptr
end type fgsl_dht
!
! Types: Root finding
!
type, public :: fgsl_root_fsolver_type
private
integer(c_int) :: which = 0
end type fgsl_root_fsolver_type
type(fgsl_root_fsolver_type), public, parameter :: &
fgsl_root_fsolver_bisection = fgsl_root_fsolver_type(1), &
fgsl_root_fsolver_brent = fgsl_root_fsolver_type(2), &
fgsl_root_fsolver_falsepos = fgsl_root_fsolver_type(3)
type, public :: fgsl_root_fdfsolver_type
private
integer(c_int) :: which = 0
end type fgsl_root_fdfsolver_type
type(fgsl_root_fdfsolver_type), public, parameter :: &
fgsl_root_fdfsolver_newton = fgsl_root_fdfsolver_type(1), &
fgsl_root_fdfsolver_secant = fgsl_root_fdfsolver_type(2), &
fgsl_root_fdfsolver_steffenson = fgsl_root_fdfsolver_type(3)
type, public :: fgsl_root_fsolver
private
type(c_ptr) :: gsl_root_fsolver = c_null_ptr
end type fgsl_root_fsolver
type, public :: fgsl_root_fdfsolver
private
type(c_ptr) :: gsl_root_fdfsolver = c_null_ptr
end type fgsl_root_fdfsolver
!
! Types: Minimization
!
type, public :: fgsl_min_fminimizer_type
private
integer(c_int) :: which = 0
end type fgsl_min_fminimizer_type
type(fgsl_min_fminimizer_type), public, parameter :: &
fgsl_min_fminimizer_goldensection = fgsl_min_fminimizer_type(1), &
fgsl_min_fminimizer_brent = fgsl_min_fminimizer_type(2), &
fgsl_min_fminimizer_quad_golden = fgsl_min_fminimizer_type(3)
type, public :: fgsl_min_fminimizer
private
type(c_ptr) :: gsl_min_fminimizer = c_null_ptr
end type fgsl_min_fminimizer
!
! Types: Multi-Root
!
type, public :: fgsl_multiroot_function
private
type(c_ptr) :: gsl_multiroot_function = c_null_ptr
end type fgsl_multiroot_function
type, public :: fgsl_multiroot_function_fdf
private
type(c_ptr) :: gsl_multiroot_function_fdf = c_null_ptr
end type fgsl_multiroot_function_fdf
type, public :: fgsl_multiroot_fsolver
private
type(c_ptr) :: gsl_multiroot_fsolver = c_null_ptr
end type fgsl_multiroot_fsolver
type, public :: fgsl_multiroot_fsolver_type
private
integer(c_int) :: which = 0
end type fgsl_multiroot_fsolver_type
type(fgsl_multiroot_fsolver_type), public, parameter :: &
fgsl_multiroot_fsolver_dnewton = fgsl_multiroot_fsolver_type(1), &
fgsl_multiroot_fsolver_broyden = fgsl_multiroot_fsolver_type(2), &
fgsl_multiroot_fsolver_hybrid = fgsl_multiroot_fsolver_type(3), &
fgsl_multiroot_fsolver_hybrids = fgsl_multiroot_fsolver_type(4)
type, public :: fgsl_multiroot_fdfsolver
private
type(c_ptr) :: gsl_multiroot_fdfsolver = c_null_ptr
end type fgsl_multiroot_fdfsolver
type, public :: fgsl_multiroot_fdfsolver_type
private
integer(c_int) :: which = 0
end type fgsl_multiroot_fdfsolver_type
type(fgsl_multiroot_fdfsolver_type), public, parameter :: &
fgsl_multiroot_fdfsolver_newton = fgsl_multiroot_fdfsolver_type(1), &
fgsl_multiroot_fdfsolver_gnewton = fgsl_multiroot_fdfsolver_type(2), &
fgsl_multiroot_fdfsolver_hybridj = fgsl_multiroot_fdfsolver_type(3), &
fgsl_multiroot_fdfsolver_hybridsj = fgsl_multiroot_fdfsolver_type(4)
!
! Types: Multi-Min
!
type, public :: fgsl_multimin_function
private
type(c_ptr) :: gsl_multimin_function = c_null_ptr
end type fgsl_multimin_function
type, public :: fgsl_multimin_function_fdf
private
type(c_ptr) :: gsl_multimin_function_fdf = c_null_ptr
end type fgsl_multimin_function_fdf
type, public :: fgsl_multimin_fminimizer
private
type(c_ptr) :: gsl_multimin_fminimizer = c_null_ptr
end type fgsl_multimin_fminimizer
type, public :: fgsl_multimin_fminimizer_type
private
integer(c_int) :: which = 0
end type fgsl_multimin_fminimizer_type
type(fgsl_multimin_fminimizer_type), public, parameter :: &
fgsl_multimin_fminimizer_nmsimplex = fgsl_multimin_fminimizer_type(1), &
fgsl_multimin_fminimizer_nmsimplex2 = fgsl_multimin_fminimizer_type(2), &
fgsl_multimin_fminimizer_nmsimplex2rand = fgsl_multimin_fminimizer_type(3)
type, public :: fgsl_multimin_fdfminimizer
private
type(c_ptr) :: gsl_multimin_fdfminimizer = c_null_ptr
end type fgsl_multimin_fdfminimizer
type, public :: fgsl_multimin_fdfminimizer_type
private
integer(c_int) :: which = 0
end type fgsl_multimin_fdfminimizer_type
type(fgsl_multimin_fdfminimizer_type), public, parameter :: &
fgsl_multimin_fdfminimizer_steepest_descent = fgsl_multimin_fdfminimizer_type(1), &
fgsl_multimin_fdfminimizer_conjugate_pr = fgsl_multimin_fdfminimizer_type(2), &
fgsl_multimin_fdfminimizer_conjugate_fr = fgsl_multimin_fdfminimizer_type(3), &
fgsl_multimin_fdfminimizer_vector_bfgs = fgsl_multimin_fdfminimizer_type(4), &
fgsl_multimin_fdfminimizer_vector_bfgs2 = fgsl_multimin_fdfminimizer_type(5)
!
! Types: Fitting
!
type, public :: fgsl_multifit_linear_workspace
private
type(c_ptr) :: gsl_multifit_linear_workspace = c_null_ptr
end type fgsl_multifit_linear_workspace
type, public :: fgsl_multifit_function
private
type(c_ptr) :: gsl_multifit_function = c_null_ptr
end type fgsl_multifit_function
type, public :: fgsl_multifit_function_fdf
private
type(c_ptr) :: gsl_multifit_function_fdf = c_null_ptr
end type fgsl_multifit_function_fdf
type, public :: fgsl_multifit_fsolver
private
type(c_ptr) :: gsl_multifit_fsolver = c_null_ptr
end type fgsl_multifit_fsolver
type, public :: fgsl_multifit_fsolver_type
private
integer(c_int) :: which = 0
end type fgsl_multifit_fsolver_type
type, public :: fgsl_multifit_fdfsolver
private
type(c_ptr) :: gsl_multifit_fdfsolver = c_null_ptr
end type fgsl_multifit_fdfsolver
type, public :: fgsl_multifit_fdfsolver_type
private
integer(c_int) :: which = 0
end type fgsl_multifit_fdfsolver_type
type(fgsl_multifit_fdfsolver_type), public, parameter :: &
fgsl_multifit_fdfsolver_lmder = fgsl_multifit_fdfsolver_type(1), &
fgsl_multifit_fdfsolver_lmsder = fgsl_multifit_fdfsolver_type(2)
!
! Types: B-Splines
!
type, public :: fgsl_bspline_workspace
private
type(c_ptr) :: gsl_bspline_workspace
end type fgsl_bspline_workspace
type, public :: fgsl_bspline_deriv_workspace
private
type(c_ptr) :: gsl_bspline_deriv_workspace
end type fgsl_bspline_deriv_workspace
!
! required C interfaces
! FGSL names occurring here are auxiliary routines
! needed to transfer static C information to the Fortran subsystem
interface
!-*-f90-*-
!
!  Interfaces: Error treatment
!
     function gsl_set_error_handler(new_handler) bind(c)
       import
       type(c_funptr), value :: new_handler
       type(c_funptr) :: gsl_set_error_handler
     end function gsl_set_error_handler
     function gsl_set_error_handler_off() bind(c)
       import
       type(c_funptr) :: gsl_set_error_handler_off
     end function gsl_set_error_handler_off
     function gsl_strerror(errno) bind(c)
       import
       integer(c_int), value :: errno
       type(c_ptr) :: gsl_strerror
     end function gsl_strerror
     subroutine gsl_error(reason, file, line, gsl_errno) bind(c)
       import
       character(c_char), dimension(*) :: reason, file
       integer(c_int), value :: line
       integer(c_int), value :: gsl_errno
     end subroutine gsl_error
!-*-f90-*-
!
!  Interfaces: miscellaneous additions
!
  function gsl_aux_sizeof_double() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_double
  end function gsl_aux_sizeof_double
  function gsl_aux_sizeof_float() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_float
  end function gsl_aux_sizeof_float
  function gsl_aux_sizeof_int() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_int
  end function gsl_aux_sizeof_int
  function gsl_aux_sizeof_long() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_long
  end function gsl_aux_sizeof_long
  function gsl_aux_sizeof_size_t() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_size_t
  end function gsl_aux_sizeof_size_t
  function gsl_aux_sizeof_char() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_char
  end function gsl_aux_sizeof_char
!-*-f90-*-
!
!  Interfaces: I/O Add-ons
!
  function fopen(path, mode) bind(c)
    import :: c_char, c_ptr
    character(kind=c_char), dimension(*) :: path, mode
    type(c_ptr) :: fopen
  end function fopen
  function fclose(fd) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: fd
    integer(c_int) :: fclose
  end function fclose
  function fgsl_cstdin() bind(c)
    import :: c_ptr
    type(c_ptr) :: fgsl_cstdin
  end function fgsl_cstdin
  function fgsl_cstdout() bind(c)
    import :: c_ptr
    type(c_ptr) :: fgsl_cstdout
  end function fgsl_cstdout
  function fgsl_cstderr() bind(c)
    import :: c_ptr
    type(c_ptr) :: fgsl_cstderr
  end function fgsl_cstderr
  function fflush(stream) bind(c)
    import :: c_int, c_ptr
    type(c_ptr), value :: stream
    integer(c_int) :: fflush
  end function fflush
!-*-f90-*-
!
!  Interfaces: Mathematical Functions
!
     function gsl_isnan(x) bind(c)
       import
       real(c_double), value :: x
       integer(c_int) :: gsl_isnan
     end function gsl_isnan
     function gsl_isinf(x) bind(c)
       import
       real(c_double), value :: x
       integer(c_int) :: gsl_isinf
     end function gsl_isinf
     function gsl_finite(x) bind(c)
       import
       real(c_double), value :: x
       integer(c_int) :: gsl_finite
     end function gsl_finite
     function gsl_log1p(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_log1p
     end function gsl_log1p
     function gsl_expm1(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_expm1
     end function gsl_expm1
     function gsl_hypot(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_hypot
     end function gsl_hypot
     function gsl_acosh(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_acosh
     end function gsl_acosh
     function gsl_asinh(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_asinh
     end function gsl_asinh
     function gsl_atanh(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_atanh
     end function gsl_atanh
     function gsl_ldexp(x,e) bind(c)
       import
       real(c_double), value :: x
       integer(c_int), value :: e
       real(c_double) :: gsl_ldexp
     end function gsl_ldexp
     function gsl_frexp(x,e) bind(c)
       import
       real(c_double), value :: x
       integer(c_int), intent(out) :: e
       real(c_double) :: gsl_frexp
     end function gsl_frexp
     function gsl_fcmp(x,y,eps) bind(c)
       import
       real(c_double), value :: x, y, eps
       integer(c_int) :: gsl_fcmp
     end function gsl_fcmp
! constructors for abstract types
     function fgsl_function_cinit(func, params) bind(c)
       import
       type(c_funptr), value :: func
       type(c_ptr), value :: params
       type(c_ptr) :: fgsl_function_cinit
     end function fgsl_function_cinit
     function fgsl_function_fdf_cinit(f, df, fdf, params) bind(c)
       import
       type(c_funptr), value :: f, df, fdf
       type(c_ptr), value :: params
       type(c_ptr) :: fgsl_function_fdf_cinit
     end function fgsl_function_fdf_cinit
     subroutine fgsl_function_cfree(sfunc) bind(c)
       import
       type(c_ptr), value :: sfunc
     end subroutine fgsl_function_cfree
     subroutine fgsl_function_fdf_cfree(sfunc) bind(c)
       import
       type(c_ptr), value :: sfunc
     end subroutine fgsl_function_fdf_cfree
! auxiliary routines
     function fgsl_fn_eval_aux(f, x) bind(c)
       import
       type(c_ptr), value :: f
       real(c_double), value :: x
       real(c_double) :: fgsl_fn_eval_aux
     end function fgsl_fn_eval_aux
     function fgsl_fn_fdf_eval_f_aux(f, x) bind(c)
       import
       type(c_ptr), value :: f
       real(c_double), value :: x
       real(c_double) :: fgsl_fn_fdf_eval_f_aux
     end function fgsl_fn_fdf_eval_f_aux
     function fgsl_fn_fdf_eval_df_aux(f, x) bind(c)
       import
       type(c_ptr), value :: f
       real(c_double), value :: x
       real(c_double) :: fgsl_fn_fdf_eval_df_aux
     end function fgsl_fn_fdf_eval_df_aux
     subroutine fgsl_fn_fdf_eval_f_df_aux(f, x, y, dy) bind(c)
       import
       type(c_ptr), value :: f
       real(c_double), value :: x
       real(c_double), intent(out) :: y, dy
     end subroutine fgsl_fn_fdf_eval_f_df_aux
!-*-f90-*-
!
!  Interfaces: Complex numbers
!
     function gsl_complex_arg(z) bind(c)
       import
       type(gsl_complex), value :: z
       real(c_double) :: gsl_complex_arg
     end function gsl_complex_arg
     function gsl_complex_logabs(z) bind(c)
       import
       type(gsl_complex), value :: z
       real(c_double) :: gsl_complex_logabs
     end function gsl_complex_logabs
     function gsl_complex_log10(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_log10
     end function gsl_complex_log10
     function gsl_complex_log_b(z,b) bind(c)
       import
       type(gsl_complex), value :: z, b
       type(gsl_complex) :: gsl_complex_log_b
     end function gsl_complex_log_b
     function gsl_complex_arcsin(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arcsin
     end function gsl_complex_arcsin
     function gsl_complex_arcsin_real(r) bind(c)
       import
       real(c_double), value :: r
       type(gsl_complex) :: gsl_complex_arcsin_real
     end function gsl_complex_arcsin_real
     function gsl_complex_arccos(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arccos
     end function gsl_complex_arccos
     function gsl_complex_arccos_real(r) bind(c)
       import
       real(c_double), value :: r
       type(gsl_complex) :: gsl_complex_arccos_real
     end function gsl_complex_arccos_real
     function gsl_complex_arctan(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arctan
     end function gsl_complex_arctan
     function gsl_complex_arcsec(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arcsec
     end function gsl_complex_arcsec
     function gsl_complex_arcsec_real(r) bind(c)
       import
       real(c_double), value :: r
       type(gsl_complex) :: gsl_complex_arcsec_real
     end function gsl_complex_arcsec_real
     function gsl_complex_arccsc(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arccsc
     end function gsl_complex_arccsc
     function gsl_complex_arccsc_real(r) bind(c)
       import
       real(c_double), value :: r
       type(gsl_complex) :: gsl_complex_arccsc_real
     end function gsl_complex_arccsc_real
     function gsl_complex_arccot(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arccot
     end function gsl_complex_arccot
     function gsl_complex_arcsinh(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arcsinh
     end function gsl_complex_arcsinh
     function gsl_complex_arccosh(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arccosh
     end function gsl_complex_arccosh
     function gsl_complex_arccosh_real(r) bind(c)
       import
       real(c_double), value :: r
       type(gsl_complex) :: gsl_complex_arccosh_real
     end function gsl_complex_arccosh_real
     function gsl_complex_arctanh(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arctanh
     end function gsl_complex_arctanh
     function gsl_complex_arctanh_real(r) bind(c)
       import
       real(c_double), value :: r
       type(gsl_complex) :: gsl_complex_arctanh_real
     end function gsl_complex_arctanh_real
     function gsl_complex_arcsech(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arcsech
     end function gsl_complex_arcsech
     function gsl_complex_arccsch(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arccsch
     end function gsl_complex_arccsch
     function gsl_complex_arccoth(z) bind(c)
       import
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_arccoth
     end function gsl_complex_arccoth
!-*-f90-*-
!
!  Interfaces: Polynomials
!
     function gsl_poly_eval(c, len, x) bind(c)
       import
       real(c_double), dimension(*), intent(in) :: c
       integer(c_int), value :: len
       real(c_double), value :: x
       real(c_double) :: gsl_poly_eval
     end function gsl_poly_eval
     function gsl_poly_complex_eval(c, len, z) bind(c)
       import
       real(c_double), dimension(*), intent(in) :: c
       integer(c_int), value :: len
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_poly_complex_eval
     end function gsl_poly_complex_eval
     function gsl_complex_poly_complex_eval(c, len, z) bind(c)
       import
       type(gsl_complex), dimension(*), intent(in) :: c
       integer(c_int), value :: len
       type(gsl_complex), value :: z
       type(gsl_complex) :: gsl_complex_poly_complex_eval
     end function gsl_complex_poly_complex_eval
     function gsl_poly_eval_derivs(c, lenc, x, res, lenres) bind(c)
       import :: c_double, c_size_t, c_int
       integer(c_int) :: gsl_poly_eval_derivs
       integer(c_size_t), value :: lenc, lenres
       real(c_double), dimension(*), intent(in) :: c
       real(c_double), dimension(*) :: res
       real(c_double), value :: x
     end function gsl_poly_eval_derivs
     function gsl_poly_dd_init(dd, x, y, size) bind(c)
       import
       real(c_double), dimension(*), intent(inout) :: dd
       real(c_double), dimension(*), intent(in) :: x, y
       integer(c_size_t), value :: size
       integer(c_int) :: gsl_poly_dd_init
     end function gsl_poly_dd_init
     function gsl_poly_dd_eval(dd, xa, size, x) bind(c)
       import
       real(c_double), dimension(*), intent(in) :: dd, xa
       real(c_double), value :: x
       integer(c_size_t), value :: size
       real(c_double) :: gsl_poly_dd_eval
     end function gsl_poly_dd_eval
     function gsl_poly_dd_taylor(c, xp, dd, x, size, w) bind(c)
       import
       real(c_double), dimension(*), intent(inout) :: c
       real(c_double), value :: xp
       real(c_double), dimension(*), intent(in) :: dd, x
       real(c_double), dimension(*), intent(out) :: w
       integer(c_size_t), value :: size
       integer(c_int) :: gsl_poly_dd_taylor
     end function gsl_poly_dd_taylor
     function gsl_poly_solve_quadratic(a, b, c, x0, x1) bind(c)
       import
       real(c_double), value :: a, b, c
       real(c_double), intent(out) :: x0, x1
       integer(c_int) :: gsl_poly_solve_quadratic
     end function gsl_poly_solve_quadratic
     function gsl_poly_complex_solve_quadratic(a, b, c, x0, x1) bind(c, name='gsl_poly_complex_solve_quadratic')
       import
       real(c_double), value :: a, b, c
       type(gsl_complex), intent(out) :: x0, x1
       integer(c_int) :: gsl_poly_complex_solve_quadratic
     end function gsl_poly_complex_solve_quadratic
     function gsl_poly_solve_cubic(a, b, c, x0, x1, x2) bind(c)
       import
       real(c_double), value :: a, b, c
       real(c_double), intent(out) :: x0, x1, x2
       integer(c_int) :: gsl_poly_solve_cubic
     end function gsl_poly_solve_cubic
     function gsl_poly_complex_solve_cubic(a, b, c, x0, x1, x2) bind(c, name='gsl_poly_complex_solve_cubic')
       import
       real(c_double), value :: a, b, c
       type(gsl_complex), intent(out) :: x0, x1, x2
       integer(c_int) :: gsl_poly_complex_solve_cubic
     end function gsl_poly_complex_solve_cubic
     function gsl_poly_complex_workspace_alloc(n) bind(c, name='gsl_poly_complex_workspace_alloc')
       import
       integer(c_size_t), value :: n
       type(c_ptr) :: gsl_poly_complex_workspace_alloc
     end function gsl_poly_complex_workspace_alloc
     subroutine gsl_poly_complex_workspace_free(w) bind(c, name='gsl_poly_complex_workspace_free')
       import
       type(c_ptr), value :: w
     end subroutine gsl_poly_complex_workspace_free
     function gsl_poly_complex_solve(a, n, w, z) bind(c, name='gsl_poly_complex_solve')
       import
       real(c_double), dimension(*), intent(in) :: a
       integer(c_size_t), value :: n
       type(c_ptr), value :: w
       real(c_double), dimension(*), intent(out) :: z
       integer(c_int) :: gsl_poly_complex_solve
     end function gsl_poly_complex_solve
!-*-f90-*-
!
!  Interfaces: Special Functions
!
     function gsl_sf_airy_ai(x, mode) bind(c, name='gsl_sf_airy_Ai')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_airy_ai
     end function gsl_sf_airy_ai
     function gsl_sf_airy_ai_e(x, mode, result) bind(c, name='gsl_sf_airy_Ai_e')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_airy_ai_e
     end function gsl_sf_airy_ai_e
     function gsl_sf_airy_bi(x, mode) bind(c, name='gsl_sf_airy_Bi')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_airy_bi
     end function gsl_sf_airy_bi
     function gsl_sf_airy_bi_e(x, mode, result) bind(c, name='gsl_sf_airy_Bi_e')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_airy_bi_e
     end function gsl_sf_airy_bi_e
     function gsl_sf_airy_ai_scaled(x, mode) bind(c, name='gsl_sf_airy_Ai_scaled')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_airy_ai_scaled
     end function gsl_sf_airy_ai_scaled
     function gsl_sf_airy_ai_scaled_e(x, mode, result) bind(c, name='gsl_sf_airy_Ai_scaled_e')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_airy_ai_scaled_e
     end function gsl_sf_airy_ai_scaled_e
     function gsl_sf_airy_bi_scaled(x, mode) bind(c, name='gsl_sf_airy_Bi_scaled')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_airy_bi_scaled
     end function gsl_sf_airy_bi_scaled
     function gsl_sf_airy_bi_scaled_e(x, mode, result) bind(c, name='gsl_sf_airy_Bi_scaled_e')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_airy_bi_scaled_e
     end function gsl_sf_airy_bi_scaled_e
!!!!!
     function gsl_sf_airy_ai_deriv(x, mode) bind(c, name='gsl_sf_airy_Ai_deriv')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_airy_ai_deriv
     end function gsl_sf_airy_ai_deriv
     function gsl_sf_airy_ai_deriv_e(x, mode, result) bind(c, name='gsl_sf_airy_Ai_deriv_e')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_airy_ai_deriv_e
     end function gsl_sf_airy_ai_deriv_e
     function gsl_sf_airy_bi_deriv(x, mode) bind(c, name='gsl_sf_airy_Bi_deriv')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_airy_bi_deriv
     end function gsl_sf_airy_bi_deriv
     function gsl_sf_airy_bi_deriv_e(x, mode, result) bind(c, name='gsl_sf_airy_Bi_deriv_e')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_airy_bi_deriv_e
     end function gsl_sf_airy_bi_deriv_e
     function gsl_sf_airy_ai_deriv_scaled(x, mode) bind(c, name='gsl_sf_airy_Ai_deriv_scaled')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_airy_ai_deriv_scaled
     end function gsl_sf_airy_ai_deriv_scaled
     function gsl_sf_airy_ai_deriv_scaled_e(x, mode, result) bind(c, name='gsl_sf_airy_Ai_deriv_scaled_e')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_airy_ai_deriv_scaled_e
     end function gsl_sf_airy_ai_deriv_scaled_e
     function gsl_sf_airy_bi_deriv_scaled(x, mode) bind(c, name='gsl_sf_airy_Bi_deriv_scaled')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_airy_bi_deriv_scaled
     end function gsl_sf_airy_bi_deriv_scaled
     function gsl_sf_airy_bi_deriv_scaled_e(x, mode, result) bind(c, name='gsl_sf_airy_Bi_deriv_scaled_e')
       import
       real(c_double), value :: x
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_airy_bi_deriv_scaled_e
     end function gsl_sf_airy_bi_deriv_scaled_e
     function gsl_sf_airy_zero_ai(s)  bind(c, name='gsl_sf_airy_zero_Ai')
       import
       integer(c_int), value :: s
       real(c_double) :: gsl_sf_airy_zero_ai
     end function gsl_sf_airy_zero_ai
     function gsl_sf_airy_zero_ai_e(s, result) bind(c, name='gsl_sf_airy_zero_Ai_e')
       import
       integer(c_int), value :: s
       integer(c_int) :: gsl_sf_airy_zero_ai_e
       type(gsl_sf_result), intent(out) :: result
     end function gsl_sf_airy_zero_ai_e
     function gsl_sf_airy_zero_bi(s)  bind(c, name='gsl_sf_airy_zero_Bi')
       import
       integer(c_int), value :: s
       real(c_double) :: gsl_sf_airy_zero_bi
     end function gsl_sf_airy_zero_bi
     function gsl_sf_airy_zero_bi_e (s, result) bind(c, name='gsl_sf_airy_zero_Bi_e')
       import
       integer(c_int), value :: s
       integer(c_int) :: gsl_sf_airy_zero_bi_e
       type(gsl_sf_result), intent(out) :: result
     end function gsl_sf_airy_zero_bi_e
     function gsl_sf_airy_zero_ai_deriv(s)  bind(c, name='gsl_sf_airy_zero_Ai_deriv')
       import
       integer(c_int), value :: s
       real(c_double) :: gsl_sf_airy_zero_ai_deriv
     end function gsl_sf_airy_zero_ai_deriv
     function gsl_sf_airy_zero_ai_deriv_e (s, result) bind(c, name='gsl_sf_airy_zero_Ai_deriv_e')
       import
       integer(c_int), value :: s
       integer(c_int) :: gsl_sf_airy_zero_ai_deriv_e
       type(gsl_sf_result), intent(out) :: result
     end function gsl_sf_airy_zero_ai_deriv_e
     function gsl_sf_airy_zero_bi_deriv(s)  bind(c, name='gsl_sf_airy_zero_Bi_deriv')
       import
       integer(c_int), value :: s
       real(c_double) :: gsl_sf_airy_zero_bi_deriv
     end function gsl_sf_airy_zero_bi_deriv
     function gsl_sf_airy_zero_bi_deriv_e (s, result) bind(c, name='gsl_sf_airy_zero_Bi_deriv_e')
       import
       integer(c_int), value :: s
       integer(c_int) :: gsl_sf_airy_zero_bi_deriv_e
       type(gsl_sf_result), intent(out) :: result
     end function gsl_sf_airy_zero_bi_deriv_e
!     function gsl_sf_bessel_jc0(x) bind(c, name='gsl_sf_bessel_J0')
     function gsl_sf_bessel_jc0(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_jc0
     end function gsl_sf_bessel_jc0
!     function gsl_sf_bessel_jc0_e(x, result) bind(c, name='gsl_sf_bessel_J0_e')
     function gsl_sf_bessel_jc0_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_jc0_e
     end function gsl_sf_bessel_jc0_e
!     function gsl_sf_bessel_jc1(x) bind(c, name='gsl_sf_bessel_J1')
     function gsl_sf_bessel_jc1(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_jc1
     end function gsl_sf_bessel_jc1
!     function gsl_sf_bessel_jc1_e(x, result) bind(c, name='gsl_sf_bessel_J1_e')
     function gsl_sf_bessel_jc1_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_jc1_e
     end function gsl_sf_bessel_jc1_e
!     function gsl_sf_bessel_jcn(n,x) bind(c, name='gsl_sf_bessel_Jn')
     function gsl_sf_bessel_jcn(n,x) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_jcn
     end function gsl_sf_bessel_jcn
!     function gsl_sf_bessel_jcn_e(n, x, result) bind(c, name='gsl_sf_bessel_Jn_e')
     function gsl_sf_bessel_jcn_e(n, x, result) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_jcn_e
     end function gsl_sf_bessel_jcn_e
!     function gsl_sf_bessel_jcn_array(nmin, nmax, x, result) bind(c, name='gsl_sf_bessel_Jn_array')
     function gsl_sf_bessel_jcn_array(nmin, nmax, x, result) bind(c)
       import
       integer(c_int), value :: nmin, nmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_jcn_array
     end function gsl_sf_bessel_jcn_array
!     function gsl_sf_bessel_yc0(x) bind(c, name='gsl_sf_bessel_Y0')
     function gsl_sf_bessel_yc0(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_yc0
     end function gsl_sf_bessel_yc0
!     function gsl_sf_bessel_yc0_e(x, result) bind(c, name='gsl_sf_bessel_Y0_e')
     function gsl_sf_bessel_yc0_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_yc0_e
     end function gsl_sf_bessel_yc0_e
!     function gsl_sf_bessel_yc1(x) bind(c, name='gsl_sf_bessel_Y1')
     function gsl_sf_bessel_yc1(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_yc1
     end function gsl_sf_bessel_yc1
!     function gsl_sf_bessel_yc1_e(x, result) bind(c, name='gsl_sf_bessel_Y1_e')
     function gsl_sf_bessel_yc1_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_yc1_e
     end function gsl_sf_bessel_yc1_e
!     function gsl_sf_bessel_ycn(n,x) bind(c, name='gsl_sf_bessel_Yn')
     function gsl_sf_bessel_ycn(n,x) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ycn
     end function gsl_sf_bessel_ycn
!     function gsl_sf_bessel_ycn_e(n, x, result) bind(c, name='gsl_sf_bessel_Yn_e')
     function gsl_sf_bessel_ycn_e(n, x, result) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ycn_e
     end function gsl_sf_bessel_ycn_e
!     function gsl_sf_bessel_ycn_array(nmin, nmax, x, result) bind(c, name='gsl_sf_bessel_Yn_array')
     function gsl_sf_bessel_ycn_array(nmin, nmax, x, result) bind(c)
       import
       integer(c_int), value :: nmin, nmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ycn_array
     end function gsl_sf_bessel_ycn_array
!     function gsl_sf_bessel_ic0(x) bind(c, name='gsl_sf_bessel_I0')
     function gsl_sf_bessel_ic0(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ic0
     end function gsl_sf_bessel_ic0
!     function gsl_sf_bessel_ic0_e(x, result) bind(c, name='gsl_sf_bessel_I0_e')
     function gsl_sf_bessel_ic0_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ic0_e
     end function gsl_sf_bessel_ic0_e
!     function gsl_sf_bessel_ic1(x) bind(c, name='gsl_sf_bessel_I1')
     function gsl_sf_bessel_ic1(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ic1
     end function gsl_sf_bessel_ic1
!     function gsl_sf_bessel_ic1_e(x, result) bind(c, name='gsl_sf_bessel_I1_e')
     function gsl_sf_bessel_ic1_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ic1_e
     end function gsl_sf_bessel_ic1_e
!     function gsl_sf_bessel_icn(n,x) bind(c, name='gsl_sf_bessel_In')
     function gsl_sf_bessel_icn(n,x) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_icn
     end function gsl_sf_bessel_icn
!     function gsl_sf_bessel_icn_e(n, x, result) bind(c, name='gsl_sf_bessel_In_e')
     function gsl_sf_bessel_icn_e(n, x, result) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_icn_e
     end function gsl_sf_bessel_icn_e
!     function gsl_sf_bessel_icn_array(nmin, nmax, x, result) bind(c, name='gsl_sf_bessel_In_array')
     function gsl_sf_bessel_icn_array(nmin, nmax, x, result) bind(c)
       import
       integer(c_int), value :: nmin, nmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_icn_array
     end function gsl_sf_bessel_icn_array
!     function gsl_sf_bessel_ic0_scaled(x) bind(c, name='gsl_sf_bessel_I0_scaled')
     function gsl_sf_bessel_ic0_scaled(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ic0_scaled
     end function gsl_sf_bessel_ic0_scaled
!     function gsl_sf_bessel_ic0_scaled_e(x, result) bind(c, name='gsl_sf_bessel_I0_scaled_e')
     function gsl_sf_bessel_ic0_scaled_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ic0_scaled_e
     end function gsl_sf_bessel_ic0_scaled_e
!     function gsl_sf_bessel_ic1_scaled(x) bind(c, name='gsl_sf_bessel_I1_scaled')
     function gsl_sf_bessel_ic1_scaled(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ic1_scaled
     end function gsl_sf_bessel_ic1_scaled
!     function gsl_sf_bessel_ic1_scaled_e(x, result) bind(c, name='gsl_sf_bessel_I1_scaled_e')
     function gsl_sf_bessel_ic1_scaled_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ic1_scaled_e
     end function gsl_sf_bessel_ic1_scaled_e
!     function gsl_sf_bessel_icn_scaled(n,x) bind(c, name='gsl_sf_bessel_In_scaled')
     function gsl_sf_bessel_icn_scaled(n,x) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_icn_scaled
     end function gsl_sf_bessel_icn_scaled
!     function gsl_sf_bessel_icn_scaled_e(n, x, result) bind(c, name='gsl_sf_bessel_In_scaled_e')
     function gsl_sf_bessel_icn_scaled_e(n, x, result) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_icn_scaled_e
     end function gsl_sf_bessel_icn_scaled_e
!     function gsl_sf_bessel_icn_scaled_array(nmin, nmax, x, result) bind(c, name='gsl_sf_bessel_In_scaled_array')
     function gsl_sf_bessel_icn_scaled_array(nmin, nmax, x, result) bind(c)
       import
       integer(c_int), value :: nmin, nmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_icn_scaled_array
     end function gsl_sf_bessel_icn_scaled_array
!     function gsl_sf_bessel_kc0(x) bind(c, name='gsl_sf_bessel_K0')
     function gsl_sf_bessel_kc0(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_kc0
     end function gsl_sf_bessel_kc0
!     function gsl_sf_bessel_kc0_e(x, result) bind(c, name='gsl_sf_bessel_K0_e')
     function gsl_sf_bessel_kc0_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_kc0_e
     end function gsl_sf_bessel_kc0_e
!     function gsl_sf_bessel_kc1(x) bind(c, name='gsl_sf_bessel_K1')
     function gsl_sf_bessel_kc1(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_kc1
     end function gsl_sf_bessel_kc1
!     function gsl_sf_bessel_kc1_e(x, result) bind(c, name='gsl_sf_bessel_K1_e')
     function gsl_sf_bessel_kc1_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_kc1_e
     end function gsl_sf_bessel_kc1_e
!     function gsl_sf_bessel_kcn(n,x) bind(c, name='gsl_sf_bessel_Kn')
     function gsl_sf_bessel_kcn(n,x) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_kcn
     end function gsl_sf_bessel_kcn
!     function gsl_sf_bessel_kcn_e(n, x, result) bind(c, name='gsl_sf_bessel_Kn_e')
     function gsl_sf_bessel_kcn_e(n, x, result) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_kcn_e
     end function gsl_sf_bessel_kcn_e
!     function gsl_sf_bessel_kcn_array(nmin, nmax, x, result) bind(c, name='gsl_sf_bessel_Kn_array')
     function gsl_sf_bessel_kcn_array(nmin, nmax, x, result) bind(c)
       import
       integer(c_int), value :: nmin, nmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_kcn_array
     end function gsl_sf_bessel_kcn_array
!     function gsl_sf_bessel_kc0_scaled(x) bind(c, name='gsl_sf_bessel_K0_scaled')
     function gsl_sf_bessel_kc0_scaled(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_kc0_scaled
     end function gsl_sf_bessel_kc0_scaled
!     function gsl_sf_bessel_kc0_scaled_e(x, result) bind(c, name='gsl_sf_bessel_K0_scaled_e')
     function gsl_sf_bessel_kc0_scaled_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_kc0_scaled_e
     end function gsl_sf_bessel_kc0_scaled_e
!     function gsl_sf_bessel_kc1_scaled(x) bind(c, name='gsl_sf_bessel_K1_scaled')
     function gsl_sf_bessel_kc1_scaled(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_kc1_scaled
     end function gsl_sf_bessel_kc1_scaled
!     function gsl_sf_bessel_kc1_scaled_e(x, result) bind(c, name='gsl_sf_bessel_K1_scaled_e')
     function gsl_sf_bessel_kc1_scaled_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_kc1_scaled_e
     end function gsl_sf_bessel_kc1_scaled_e
!     function gsl_sf_bessel_kcn_scaled(n,x) bind(c, name='gsl_sf_bessel_Kn_scaled')
     function gsl_sf_bessel_kcn_scaled(n,x) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_kcn_scaled
     end function gsl_sf_bessel_kcn_scaled
!     function gsl_sf_bessel_kcn_scaled_e(n, x, result) bind(c, name='gsl_sf_bessel_Kn_scaled_e')
     function gsl_sf_bessel_kcn_scaled_e(n, x, result) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_kcn_scaled_e
     end function gsl_sf_bessel_kcn_scaled_e
!     function gsl_sf_bessel_kcn_scaled_array(nmin, nmax, x, result) bind(c, name='gsl_sf_bessel_Kn_scaled_array')
     function gsl_sf_bessel_kcn_scaled_array(nmin, nmax, x, result) bind(c)
       import
       integer(c_int), value :: nmin, nmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_kcn_scaled_array
     end function gsl_sf_bessel_kcn_scaled_array
!    spherical bessel functions
     function gsl_sf_bessel_js0(x) bind(c, name='gsl_sf_bessel_j0')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_js0
     end function gsl_sf_bessel_js0
     function gsl_sf_bessel_js0_e(x, result) bind(c, name='gsl_sf_bessel_j0_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_js0_e
     end function gsl_sf_bessel_js0_e
     function gsl_sf_bessel_js1(x) bind(c, name='gsl_sf_bessel_j1')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_js1
     end function gsl_sf_bessel_js1
     function gsl_sf_bessel_js1_e(x, result) bind(c, name='gsl_sf_bessel_j1_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_js1_e
     end function gsl_sf_bessel_js1_e
     function gsl_sf_bessel_js2(x) bind(c, name='gsl_sf_bessel_j2')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_js2
     end function gsl_sf_bessel_js2
     function gsl_sf_bessel_js2_e(x, result) bind(c, name='gsl_sf_bessel_j2_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_js2_e
     end function gsl_sf_bessel_js2_e
     function gsl_sf_bessel_jsl(n,x) bind(c, name='gsl_sf_bessel_jl')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_jsl
     end function gsl_sf_bessel_jsl
     function gsl_sf_bessel_jsl_e(n, x, result) bind(c, name='gsl_sf_bessel_jl_e')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_jsl_e
     end function gsl_sf_bessel_jsl_e
     function gsl_sf_bessel_jsl_array(lmax, x, result) bind(c, name='gsl_sf_bessel_jl_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_jsl_array
     end function gsl_sf_bessel_jsl_array
     function gsl_sf_bessel_jsl_steed_array(lmax, x, result) bind(c, name='gsl_sf_bessel_jl_steed_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_jsl_steed_array
     end function gsl_sf_bessel_jsl_steed_array
     function gsl_sf_bessel_ys0(x) bind(c, name='gsl_sf_bessel_y0')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ys0
     end function gsl_sf_bessel_ys0
     function gsl_sf_bessel_ys0_e(x, result) bind(c, name='gsl_sf_bessel_y0_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ys0_e
     end function gsl_sf_bessel_ys0_e
     function gsl_sf_bessel_ys1(x) bind(c, name='gsl_sf_bessel_y1')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ys1
     end function gsl_sf_bessel_ys1
     function gsl_sf_bessel_ys1_e(x, result) bind(c, name='gsl_sf_bessel_y1_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ys1_e
     end function gsl_sf_bessel_ys1_e
     function gsl_sf_bessel_ys2(x) bind(c, name='gsl_sf_bessel_y2')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ys2
     end function gsl_sf_bessel_ys2
     function gsl_sf_bessel_ys2_e(x, result) bind(c, name='gsl_sf_bessel_y2_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ys2_e
     end function gsl_sf_bessel_ys2_e
     function gsl_sf_bessel_ysl(n,x) bind(c, name='gsl_sf_bessel_yl')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ysl
     end function gsl_sf_bessel_ysl
     function gsl_sf_bessel_ysl_e(n, x, result) bind(c, name='gsl_sf_bessel_yl_e')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ysl_e
     end function gsl_sf_bessel_ysl_e
     function gsl_sf_bessel_ysl_array(lmax, x, result) bind(c, name='gsl_sf_bessel_yl_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ysl_array
     end function gsl_sf_bessel_ysl_array
     function gsl_sf_bessel_is0_scaled(x) bind(c, name='gsl_sf_bessel_i0_scaled')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_is0_scaled
     end function gsl_sf_bessel_is0_scaled
     function gsl_sf_bessel_is0_scaled_e(x, result) bind(c, name='gsl_sf_bessel_i0_scaled_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_is0_scaled_e
     end function gsl_sf_bessel_is0_scaled_e
     function gsl_sf_bessel_is1_scaled(x) bind(c, name='gsl_sf_bessel_i1_scaled')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_is1_scaled
     end function gsl_sf_bessel_is1_scaled
     function gsl_sf_bessel_is1_scaled_e(x, result) bind(c, name='gsl_sf_bessel_i1_scaled_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_is1_scaled_e
     end function gsl_sf_bessel_is1_scaled_e
     function gsl_sf_bessel_is2_scaled(x) bind(c, name='gsl_sf_bessel_i2_scaled')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_is2_scaled
     end function gsl_sf_bessel_is2_scaled
     function gsl_sf_bessel_is2_scaled_e(x, result) bind(c, name='gsl_sf_bessel_i2_scaled_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_is2_scaled_e
     end function gsl_sf_bessel_is2_scaled_e
     function gsl_sf_bessel_isl_scaled(n,x) bind(c, name='gsl_sf_bessel_il_scaled')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_isl_scaled
     end function gsl_sf_bessel_isl_scaled
     function gsl_sf_bessel_isl_scaled_e(n, x, result) bind(c, name='gsl_sf_bessel_il_scaled_e')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_isl_scaled_e
     end function gsl_sf_bessel_isl_scaled_e
     function gsl_sf_bessel_isl_scaled_array(lmax, x, result) bind(c, name='gsl_sf_bessel_il_scaled_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_isl_scaled_array
     end function gsl_sf_bessel_isl_scaled_array
     function gsl_sf_bessel_ks0_scaled(x) bind(c, name='gsl_sf_bessel_k0_scaled')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ks0_scaled
     end function gsl_sf_bessel_ks0_scaled
     function gsl_sf_bessel_ks0_scaled_e(x, result) bind(c, name='gsl_sf_bessel_k0_scaled_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ks0_scaled_e
     end function gsl_sf_bessel_ks0_scaled_e
     function gsl_sf_bessel_ks1_scaled(x) bind(c, name='gsl_sf_bessel_k1_scaled')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ks1_scaled
     end function gsl_sf_bessel_ks1_scaled
     function gsl_sf_bessel_ks1_scaled_e(x, result) bind(c, name='gsl_sf_bessel_k1_scaled_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ks1_scaled_e
     end function gsl_sf_bessel_ks1_scaled_e
     function gsl_sf_bessel_ks2_scaled(x) bind(c, name='gsl_sf_bessel_k2_scaled')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ks2_scaled
     end function gsl_sf_bessel_ks2_scaled
     function gsl_sf_bessel_ks2_scaled_e(x, result) bind(c, name='gsl_sf_bessel_k2_scaled_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ks2_scaled_e
     end function gsl_sf_bessel_ks2_scaled_e
     function gsl_sf_bessel_ksl_scaled(n,x) bind(c, name='gsl_sf_bessel_kl_scaled')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ksl_scaled
     end function gsl_sf_bessel_ksl_scaled
     function gsl_sf_bessel_ksl_scaled_e(n, x, result) bind(c, name='gsl_sf_bessel_kl_scaled_e')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ksl_scaled_e
     end function gsl_sf_bessel_ksl_scaled_e
     function gsl_sf_bessel_ksl_scaled_array(lmax, x, result) bind(c, name='gsl_sf_bessel_kl_scaled_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ksl_scaled_array
     end function gsl_sf_bessel_ksl_scaled_array
!  fractional order bessel functions
     function gsl_sf_bessel_jnu(n,x) bind(c, name='gsl_sf_bessel_Jnu')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_jnu
     end function gsl_sf_bessel_jnu
     function gsl_sf_bessel_jnu_e(n, x, result) bind(c, name='gsl_sf_bessel_Jnu_e')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_jnu_e
     end function gsl_sf_bessel_jnu_e
     function gsl_sf_bessel_sequence_jnu_e(nu, mode, size, v) bind(c, name='gsl_sf_bessel_sequence_Jnu_e')
       import
       real(c_double), value :: nu
       integer(c_int), value :: mode
       integer(c_size_t), value :: size
       real(c_double), dimension(*), intent(inout) :: v
       integer(c_int) :: gsl_sf_bessel_sequence_jnu_e
     end function gsl_sf_bessel_sequence_jnu_e
     function gsl_sf_bessel_ynu(n,x) bind(c, name='gsl_sf_bessel_Ynu')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_ynu
     end function gsl_sf_bessel_ynu
     function gsl_sf_bessel_ynu_e(n, x, result) bind(c, name='gsl_sf_bessel_Ynu_e')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_ynu_e
     end function gsl_sf_bessel_ynu_e
     function gsl_sf_bessel_inu(n,x) bind(c, name='gsl_sf_bessel_Inu')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_inu
     end function gsl_sf_bessel_inu
     function gsl_sf_bessel_inu_e(n, x, result) bind(c, name='gsl_sf_bessel_Inu_e')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_inu_e
     end function gsl_sf_bessel_inu_e
     function gsl_sf_bessel_inu_scaled(n,x) bind(c, name='gsl_sf_bessel_Inu_scaled')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_inu_scaled
     end function gsl_sf_bessel_inu_scaled
     function gsl_sf_bessel_inu_scaled_e(n, x, result) bind(c, name='gsl_sf_bessel_Inu_scaled_e')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_inu_scaled_e
     end function gsl_sf_bessel_inu_scaled_e
     function gsl_sf_bessel_knu(n,x) bind(c, name='gsl_sf_bessel_Knu')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_knu
     end function gsl_sf_bessel_knu
     function gsl_sf_bessel_knu_e(n, x, result) bind(c, name='gsl_sf_bessel_Knu_e')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_knu_e
     end function gsl_sf_bessel_knu_e
     function gsl_sf_bessel_lnknu(n,x) bind(c, name='gsl_sf_bessel_lnKnu')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_lnknu
     end function gsl_sf_bessel_lnknu
     function gsl_sf_bessel_lnknu_e(n, x, result) bind(c, name='gsl_sf_bessel_lnKnu_e')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_lnknu_e
     end function gsl_sf_bessel_lnknu_e
     function gsl_sf_bessel_knu_scaled(n,x) bind(c, name='gsl_sf_bessel_Knu_scaled')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_bessel_knu_scaled
     end function gsl_sf_bessel_knu_scaled
     function gsl_sf_bessel_knu_scaled_e(n, x, result) bind(c, name='gsl_sf_bessel_Knu_scaled_e')
       import
       real(c_double), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_knu_scaled_e
     end function gsl_sf_bessel_knu_scaled_e
     function gsl_sf_bessel_zero_jc0 (s) bind(c, name='gsl_sf_bessel_zero_J0')
       import
       integer(c_int), value :: s
       real(c_double) :: gsl_sf_bessel_zero_jc0
     end function gsl_sf_bessel_zero_jc0
     function gsl_sf_bessel_zero_jc0_e (s, result) bind(c, name='gsl_sf_bessel_zero_J0_e')
       import
       integer(c_int), value :: s
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_zero_jc0_e
     end function gsl_sf_bessel_zero_jc0_e
     function gsl_sf_bessel_zero_jc1 (s) bind(c, name='gsl_sf_bessel_zero_J1')
       import
       integer(c_int), value :: s
       real(c_double) :: gsl_sf_bessel_zero_jc1
     end function gsl_sf_bessel_zero_jc1
     function gsl_sf_bessel_zero_jc1_e (s, result) bind(c, name='gsl_sf_bessel_zero_J1_e')
       import
       integer(c_int), value :: s
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_zero_jc1_e
     end function gsl_sf_bessel_zero_jc1_e
     function gsl_sf_bessel_zero_jnu(nu, s) bind(c, name='gsl_sf_bessel_zero_Jnu')
       import
       real(c_double), value :: nu
       integer(c_int), value :: s
       real(c_double) :: gsl_sf_bessel_zero_jnu
     end function gsl_sf_bessel_zero_jnu
     function gsl_sf_bessel_zero_jnu_e (nu, s, result) bind(c, name='gsl_sf_bessel_zero_Jnu_e')
       import
       real(c_double), value :: nu
       integer(c_int), value :: s
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_bessel_zero_jnu_e
     end function gsl_sf_bessel_zero_jnu_e
!
     function gsl_sf_clausen(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_clausen
     end function gsl_sf_clausen
     function gsl_sf_clausen_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_clausen_e
     end function gsl_sf_clausen_e
     function gsl_sf_hydrogenicr_1(z, r) bind(c, name='gsl_sf_hydrogenicR_1')
       import
       real(c_double), value :: z, r
       real(c_double) :: gsl_sf_hydrogenicr_1
     end function gsl_sf_hydrogenicr_1
     function gsl_sf_hydrogenicr_1_e(z, r, result) bind(c, name='gsl_sf_hydrogenicR_1_e')
       import
       real(c_double), value :: z, r
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hydrogenicr_1_e
     end function gsl_sf_hydrogenicr_1_e
     function gsl_sf_hydrogenicr(n, l, z, r) bind(c, name='gsl_sf_hydrogenicR')
       import
       integer(c_int), value :: n, l
       real(c_double), value :: z, r
       real(c_double) :: gsl_sf_hydrogenicr
     end function gsl_sf_hydrogenicr
     function gsl_sf_hydrogenicr_e(n, l, z, r, result) &
          bind(c, name='gsl_sf_hydrogenicR_e')
       import
       integer(c_int), value :: n, l
       real(c_double), value :: z, r
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hydrogenicr_e
     end function gsl_sf_hydrogenicr_e
     function gsl_sf_coulomb_wave_fg_e(eta, x, l_f, k, f, fp, g, gp, exp_f, exp_g) &
          bind(c, name='gsl_sf_coulomb_wave_FG_e')
       import
       real(c_double), value :: eta, x, l_f
       integer(c_int), value :: k
       type(gsl_sf_result), intent(out) :: f, fp, g, gp
       real(c_double), intent(out) :: exp_f, exp_g
       integer(c_int) :: gsl_sf_coulomb_wave_fg_e
     end function gsl_sf_coulomb_wave_fg_e
     function gsl_sf_coulomb_wave_f_array (l_min, kmax, eta, x, fc_array, &
          f_exponent) bind(c, name='gsl_sf_coulomb_wave_F_array') 
       import
       real(c_double), value :: l_min, eta, x
       integer(c_int), value :: kmax
       real(c_double), dimension(*), intent(out) :: fc_array
       real(c_double), intent(out) :: f_exponent
       integer(c_int) :: gsl_sf_coulomb_wave_f_array
     end function gsl_sf_coulomb_wave_f_array
     function gsl_sf_coulomb_wave_fg_array (l_min, kmax, eta, x, fc_array, &
          gc_array, f_exponent, g_exponent) bind(c, name='gsl_sf_coulomb_wave_FG_array') 
       import
       real(c_double), value :: l_min, eta, x
       integer(c_int), value :: kmax
       real(c_double), dimension(*), intent(out) :: fc_array, gc_array
       real(c_double), intent(out) :: f_exponent, g_exponent
       integer(c_int) :: gsl_sf_coulomb_wave_fg_array
     end function gsl_sf_coulomb_wave_fg_array
     function gsl_sf_coulomb_wave_fgp_array (l_min, kmax, eta, x, fc_array, fcp_array, &
          gc_array, gcp_array, f_exponent, g_exponent) &
          bind(c, name='gsl_sf_coulomb_wave_FGp_array') 
       import
       real(c_double), value :: l_min, eta, x
       integer(c_int), value :: kmax
       real(c_double), dimension(*), intent(out) :: fc_array, gc_array, &
            fcp_array, gcp_array
       real(c_double), intent(out) :: f_exponent, g_exponent
       integer(c_int) :: gsl_sf_coulomb_wave_fgp_array
     end function gsl_sf_coulomb_wave_fgp_array
     function gsl_sf_coulomb_wave_sphf_array (l_min, kmax, eta, x, fc_array, &
          f_exponent) bind(c, name='gsl_sf_coulomb_wave_sphF_array') 
       import
       real(c_double), value :: l_min, eta, x
       integer(c_int), value :: kmax
       real(c_double), dimension(*), intent(out) :: fc_array
       real(c_double), intent(out) :: f_exponent
       integer(c_int) :: gsl_sf_coulomb_wave_sphf_array
     end function gsl_sf_coulomb_wave_sphf_array
     function gsl_sf_coulomb_cl_e(l, eta, result) bind(c, name='gsl_sf_coulomb_CL_e')
       import
       real(c_double), value :: eta, l
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_coulomb_cl_e
     end function gsl_sf_coulomb_cl_e
     function gsl_sf_coulomb_cl_array (l_min, kmax, eta, cl) &
          bind(c, name='gsl_sf_coulomb_CL_array')
       import
       real(c_double), value :: l_min, eta
       integer(c_int), value :: kmax
       real(c_double), dimension(*), intent(out) :: cl
       integer(c_int) :: gsl_sf_coulomb_cl_array
     end function gsl_sf_coulomb_cl_array
     function gsl_sf_coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc) bind(c)
       import
       integer(c_int), value :: two_ja, two_jb, two_jc, two_ma, two_mb, two_mc
       real(c_double) :: gsl_sf_coupling_3j
     end function gsl_sf_coupling_3j
     function gsl_sf_coupling_3j_e(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc, result) &
          bind(c)
       import
       integer(c_int), value :: two_ja, two_jb, two_jc, two_ma, two_mb, two_mc
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_coupling_3j_e
     end function gsl_sf_coupling_3j_e
     function gsl_sf_coupling_6j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf) bind(c)
       import
       integer(c_int), value :: two_ja, two_jb, two_jc, two_jd, two_je, two_jf
       real(c_double) :: gsl_sf_coupling_6j
     end function gsl_sf_coupling_6j
     function gsl_sf_coupling_6j_e(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, result) &
          bind(c)
       import
       integer(c_int), value :: two_ja, two_jb, two_jc, two_jd, two_je, two_jf
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_coupling_6j_e
     end function gsl_sf_coupling_6j_e
     function gsl_sf_coupling_9j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, &
          two_jg, two_jh, two_ji) bind(c)
       import
       integer(c_int), value :: two_ja, two_jb, two_jc, two_jd, two_je, two_jf, &
            two_jg, two_jh, two_ji
       real(c_double) :: gsl_sf_coupling_9j
     end function gsl_sf_coupling_9j
     function gsl_sf_coupling_9j_e(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, &
          two_jg, two_jh, two_ji, result) bind(c)
       import
       integer(c_int), value :: two_ja, two_jb, two_jc, two_jd, two_je, two_jf, &
            two_jg, two_jh, two_ji
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_coupling_9j_e
     end function gsl_sf_coupling_9j_e
     function gsl_sf_dawson(x) bind(c, name='gsl_sf_dawson')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_dawson
     end function gsl_sf_dawson
     function gsl_sf_dawson_e(x, result) bind(c, name='gsl_sf_dawson_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_dawson_e
     end function gsl_sf_dawson_e
     function gsl_sf_debye_1(x) bind(c, name='gsl_sf_debye_1')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_debye_1
     end function gsl_sf_debye_1
     function gsl_sf_debye_1_e(x, result) bind(c, name='gsl_sf_debye_1_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_debye_1_e
     end function gsl_sf_debye_1_e
     function gsl_sf_debye_2(x) bind(c, name='gsl_sf_debye_2')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_debye_2
     end function gsl_sf_debye_2
     function gsl_sf_debye_2_e(x, result) bind(c, name='gsl_sf_debye_2_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_debye_2_e
     end function gsl_sf_debye_2_e
     function gsl_sf_debye_3(x) bind(c, name='gsl_sf_debye_3')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_debye_3
     end function gsl_sf_debye_3
     function gsl_sf_debye_3_e(x, result) bind(c, name='gsl_sf_debye_3_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_debye_3_e
     end function gsl_sf_debye_3_e
     function gsl_sf_debye_4(x) bind(c, name='gsl_sf_debye_4')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_debye_4
     end function gsl_sf_debye_4
     function gsl_sf_debye_4_e(x, result) bind(c, name='gsl_sf_debye_4_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_debye_4_e
     end function gsl_sf_debye_4_e
     function gsl_sf_debye_5(x) bind(c, name='gsl_sf_debye_5')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_debye_5
     end function gsl_sf_debye_5
     function gsl_sf_debye_5_e(x, result) bind(c, name='gsl_sf_debye_5_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_debye_5_e
     end function gsl_sf_debye_5_e
     function gsl_sf_debye_6(x) bind(c, name='gsl_sf_debye_6')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_debye_6
     end function gsl_sf_debye_6
     function gsl_sf_debye_6_e(x, result) bind(c, name='gsl_sf_debye_6_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_debye_6_e
     end function gsl_sf_debye_6_e
     function gsl_sf_dilog(x) bind(c, name='gsl_sf_dilog')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_dilog
     end function gsl_sf_dilog
     function gsl_sf_dilog_e(x, result) bind(c, name='gsl_sf_dilog_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_dilog_e
     end function gsl_sf_dilog_e
     function gsl_sf_complex_dilog_e(r, theta, result_re, result_im) bind(c) 
       import
       real(c_double), value :: r, theta
       type(gsl_sf_result), intent(out) :: result_re, result_im
       integer(c_int) :: gsl_sf_complex_dilog_e
     end function gsl_sf_complex_dilog_e
     function gsl_sf_multiply_e(x, y, result) bind(c)
       import
       real(c_double), value :: x, y
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_multiply_e       
     end function gsl_sf_multiply_e
     function gsl_sf_multiply_err_e(x, dx, y, dy, result) bind(c)
       import
       real(c_double), value :: x, y, dx, dy
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_multiply_err_e       
     end function gsl_sf_multiply_err_e
     function gsl_sf_ellint_kcomp(k, mode) bind(c, name='gsl_sf_ellint_Kcomp')
       import
       real(c_double), value :: k
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_ellint_kcomp
     end function gsl_sf_ellint_kcomp
     function gsl_sf_ellint_kcomp_e(k, mode, result) bind(c, name='gsl_sf_ellint_Kcomp_e')
       import
       real(c_double), value :: k
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ellint_kcomp_e
     end function gsl_sf_ellint_kcomp_e
     function gsl_sf_ellint_ecomp(k, mode) bind(c, name='gsl_sf_ellint_Ecomp')
       import
       real(c_double), value :: k
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_ellint_ecomp
     end function gsl_sf_ellint_ecomp
     function gsl_sf_ellint_ecomp_e(k, mode, result) bind(c, name='gsl_sf_ellint_Ecomp_e')
       import
       real(c_double), value :: k
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ellint_ecomp_e
     end function gsl_sf_ellint_ecomp_e
     function gsl_sf_ellint_pcomp(k, n, mode) bind(c, name='gsl_sf_ellint_Pcomp')
       import
       real(c_double), value :: k, n
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_ellint_pcomp
     end function gsl_sf_ellint_pcomp
     function gsl_sf_ellint_pcomp_e(k, n, mode, result) bind(c, name='gsl_sf_ellint_Pcomp_e')
       import
       real(c_double), value :: k, n
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ellint_pcomp_e
     end function gsl_sf_ellint_pcomp_e
     function gsl_sf_ellint_f(phi, k, mode) bind(c, name='gsl_sf_ellint_F')
       import
       real(c_double), value :: phi, k
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_ellint_f
     end function gsl_sf_ellint_f
     function gsl_sf_ellint_f_e(phi, k, mode, result) bind(c, name='gsl_sf_ellint_F_e')
       import
       real(c_double), value :: phi, k
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ellint_f_e
     end function gsl_sf_ellint_f_e
     function gsl_sf_ellint_e(phi, k, mode) bind(c, name='gsl_sf_ellint_E')
       import
       real(c_double), value :: phi, k
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_ellint_e
     end function gsl_sf_ellint_e
     function gsl_sf_ellint_e_e(phi, k, mode, result) bind(c, name='gsl_sf_ellint_E_e')
       import
       real(c_double), value :: phi, k
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ellint_e_e
     end function gsl_sf_ellint_e_e
     function gsl_sf_ellint_p(phi, k, n, mode) bind(c, name='gsl_sf_ellint_P')
       import
       real(c_double), value :: phi, k, n
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_ellint_p
     end function gsl_sf_ellint_p
     function gsl_sf_ellint_p_e(phi, k, n, mode, result) bind(c, name='gsl_sf_ellint_P_e')
       import
       real(c_double), value :: phi, k, n
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ellint_p_e
     end function gsl_sf_ellint_p_e
     function gsl_sf_ellint_d(phi, k, n, mode) bind(c, name='gsl_sf_ellint_D')
       import
       real(c_double), value :: phi, k, n
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_ellint_d
     end function gsl_sf_ellint_d
     function gsl_sf_ellint_d_e(phi, k, n, mode, result) bind(c, name='gsl_sf_ellint_D_e')
       import
       real(c_double), value :: phi, k, n
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ellint_d_e
     end function gsl_sf_ellint_d_e
     function gsl_sf_ellint_rc(x, y, mode) bind(c, name='gsl_sf_ellint_RC')
       import
       real(c_double), value :: x, y
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_ellint_rc
     end function gsl_sf_ellint_rc
     function gsl_sf_ellint_rc_e(x, y, mode, result) bind(c, name='gsl_sf_ellint_RC_e')
       import
       real(c_double), value :: x, y
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ellint_rc_e
     end function gsl_sf_ellint_rc_e
     function gsl_sf_ellint_rd(x, y, z, mode) bind(c, name='gsl_sf_ellint_RD')
       import
       real(c_double), value :: x, y, z
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_ellint_rd
     end function gsl_sf_ellint_rd
     function gsl_sf_ellint_rd_e(x, y, z, mode, result) bind(c, name='gsl_sf_ellint_RD_e')
       import
       real(c_double), value :: x, y, z
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ellint_rd_e
     end function gsl_sf_ellint_rd_e
     function gsl_sf_ellint_rf(x, y, z, mode) bind(c, name='gsl_sf_ellint_RF')
       import
       real(c_double), value :: x, y, z
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_ellint_rf
     end function gsl_sf_ellint_rf
     function gsl_sf_ellint_rf_e(x, y, z, mode, result) bind(c, name='gsl_sf_ellint_RF_e')
       import
       real(c_double), value :: x, y, z
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ellint_rf_e
     end function gsl_sf_ellint_rf_e
     function gsl_sf_ellint_rj(x, y, z, p, mode) bind(c, name='gsl_sf_ellint_RJ')
       import
       real(c_double), value :: x, y, z, p
       integer(c_int), value :: mode
       real(c_double) :: gsl_sf_ellint_rj
     end function gsl_sf_ellint_rj
     function gsl_sf_ellint_rj_e(x, y, z, p, mode, result) bind(c, name='gsl_sf_ellint_RJ_e')
       import
       real(c_double), value :: x, y, z, p
       integer(c_int), value :: mode
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ellint_rj_e
     end function gsl_sf_ellint_rj_e
     function gsl_sf_elljac_e(u, m, sn, cn, dn) bind(c)
       import
       real(c_double), value :: u, m
       real(c_double), intent(out) :: sn, cn, dn
       integer(c_int) :: gsl_sf_elljac_e
     end function gsl_sf_elljac_e
     function gsl_sf_erf(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_erf
     end function gsl_sf_erf
     function gsl_sf_erf_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_erf_e
     end function gsl_sf_erf_e
     function gsl_sf_erfc(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_erfc
     end function gsl_sf_erfc
     function gsl_sf_erfc_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_erfc_e
     end function gsl_sf_erfc_e
     function gsl_sf_log_erfc(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_log_erfc
     end function gsl_sf_log_erfc
     function gsl_sf_log_erfc_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_log_erfc_e
     end function gsl_sf_log_erfc_e
     function gsl_sf_erf_z(x) bind(c, name='gsl_sf_erf_Z')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_erf_Z
     end function gsl_sf_erf_z
     function gsl_sf_erf_z_e(x, result) bind(c, name='gsl_sf_erf_Z_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_erf_z_e
     end function gsl_sf_erf_z_e
     function gsl_sf_erf_q(x) bind(c, name='gsl_sf_erf_Q')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_erf_q
     end function gsl_sf_erf_q
     function gsl_sf_erf_q_e(x, result) bind(c, name='gsl_sf_erf_Q_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_erf_q_e
     end function gsl_sf_erf_q_e
     function gsl_sf_hazard(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_hazard
     end function gsl_sf_hazard
     function gsl_sf_hazard_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hazard_e
     end function gsl_sf_hazard_e
     function gsl_sf_exp(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_exp
     end function gsl_sf_exp
     function gsl_sf_exp_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_exp_e
     end function gsl_sf_exp_e
     function gsl_sf_exp_e10_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result_e10), intent(out) :: result
       integer(c_int) :: gsl_sf_exp_e10_e
     end function gsl_sf_exp_e10_e
     function gsl_sf_exp_mult(x, y) bind(c)
       import
       real(c_double), value :: x, y
       real(c_double) :: gsl_sf_exp_mult
     end function gsl_sf_exp_mult
     function gsl_sf_exp_mult_e(x, y, result) bind(c)
       import
       real(c_double), value :: x, y
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_exp_mult_e
     end function gsl_sf_exp_mult_e
     function gsl_sf_exp_mult_e10_e(x, y, result) bind(c)
       import
       real(c_double), value :: x, y
       type(gsl_sf_result_e10), intent(out) :: result
       integer(c_int) :: gsl_sf_exp_mult_e10_e
     end function gsl_sf_exp_mult_e10_e
     function gsl_sf_expm1(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_expm1
     end function gsl_sf_expm1
     function gsl_sf_expm1_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_expm1_e
     end function gsl_sf_expm1_e
     function gsl_sf_exprel(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_exprel
     end function gsl_sf_exprel
     function gsl_sf_exprel_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_exprel_e
     end function gsl_sf_exprel_e
     function gsl_sf_exprel_2(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_exprel_2
     end function gsl_sf_exprel_2
     function gsl_sf_exprel_2_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_exprel_2_e
     end function gsl_sf_exprel_2_e
     function gsl_sf_exprel_n(n, x) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_exprel_n
     end function gsl_sf_exprel_n
     function gsl_sf_exprel_n_e(n, x, result) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_exprel_n_e
     end function gsl_sf_exprel_n_e
     function gsl_sf_exp_err_e(x, dx, result) bind(c)
       import
       real(c_double), value :: x, dx
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_exp_err_e
     end function gsl_sf_exp_err_e
     function gsl_sf_exp_err_e10_e(x, dx, result) bind(c)
       import
       real(c_double), value :: x, dx
       type(gsl_sf_result_e10), intent(out) :: result
       integer(c_int) :: gsl_sf_exp_err_e10_e
     end function gsl_sf_exp_err_e10_e
     function gsl_sf_exp_mult_err_e(x, dx, y, dy, result) bind(c)
       import
       real(c_double), value :: x, dx, y, dy
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_exp_mult_err_e
     end function gsl_sf_exp_mult_err_e
     function gsl_sf_exp_mult_err_e10_e(x, dx, y, dy, result) bind(c)
       import
       real(c_double), value :: x, dx, y, dy
       type(gsl_sf_result_e10), intent(out) :: result
       integer(c_int) :: gsl_sf_exp_mult_err_e10_e
     end function gsl_sf_exp_mult_err_e10_e
     function gsl_sf_expint_e1(x) bind(c, name='gsl_sf_expint_E1')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_expint_e1
     end function gsl_sf_expint_e1
     function gsl_sf_expint_e1_e(x, result) bind(c, name='gsl_sf_expint_E1_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_expint_e1_e
     end function gsl_sf_expint_e1_e
     function gsl_sf_expint_e2(x) bind(c, name='gsl_sf_expint_E2')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_expint_e2
     end function gsl_sf_expint_e2
     function gsl_sf_expint_e2_e(x, result) bind(c, name='gsl_sf_expint_E2_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_expint_e2_e
     end function gsl_sf_expint_e2_e
     function gsl_sf_expint_en(n, x) bind(c, name='gsl_sf_expint_En')
       import 
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_expint_en
     end function gsl_sf_expint_en
     function gsl_sf_expint_en_e(n, x, result) bind(c, name='gsl_sf_expint_En_e')
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_expint_en_e
     end function gsl_sf_expint_en_e
      function gsl_sf_expint_ei(x) bind(c, name='gsl_sf_expint_Ei')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_expint_ei
     end function gsl_sf_expint_ei
     function gsl_sf_expint_ei_e(x, result) bind(c, name='gsl_sf_expint_Ei_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_expint_ei_e
     end function gsl_sf_expint_ei_e
     function gsl_sf_shi(x) bind(c, name='gsl_sf_Shi')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_shi
     end function gsl_sf_shi
     function gsl_sf_shi_e(x, result) bind(c, name='gsl_sf_Shi_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_shi_e
     end function gsl_sf_shi_e
     function gsl_sf_chi(x) bind(c, name='gsl_sf_Chi')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_chi
     end function gsl_sf_chi
     function gsl_sf_chi_e(x, result) bind(c, name='gsl_sf_Chi_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_chi_e
     end function gsl_sf_chi_e
     function gsl_sf_expint_3(x) bind(c, name='gsl_sf_expint_3')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_expint_3
     end function gsl_sf_expint_3
     function gsl_sf_expint_3_e(x, result) bind(c, name='gsl_sf_expint_3_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_expint_3_e
     end function gsl_sf_expint_3_e
     function gsl_sf_si(x) bind(c, name='gsl_sf_Si')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_si
     end function gsl_sf_si
     function gsl_sf_si_e(x, result) bind(c, name='gsl_sf_Si_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_si_e
     end function gsl_sf_si_e
     function gsl_sf_ci(x) bind(c, name='gsl_sf_Ci')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_ci
     end function gsl_sf_ci
     function gsl_sf_ci_e(x, result) bind(c, name='gsl_sf_Ci_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_ci_e
     end function gsl_sf_ci_e
     function gsl_sf_atanint(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_atanint
     end function gsl_sf_atanint
     function gsl_sf_atanint_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_atanint_e
     end function gsl_sf_atanint_e
     function gsl_sf_fermi_dirac_m1(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_fermi_dirac_m1
     end function gsl_sf_fermi_dirac_m1
     function gsl_sf_fermi_dirac_m1_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_fermi_dirac_m1_e
     end function gsl_sf_fermi_dirac_m1_e
     function gsl_sf_fermi_dirac_0(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_fermi_dirac_0
     end function gsl_sf_fermi_dirac_0
     function gsl_sf_fermi_dirac_0_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_fermi_dirac_0_e
     end function gsl_sf_fermi_dirac_0_e
     function gsl_sf_fermi_dirac_1(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_fermi_dirac_1
     end function gsl_sf_fermi_dirac_1
     function gsl_sf_fermi_dirac_1_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_fermi_dirac_1_e
     end function gsl_sf_fermi_dirac_1_e
     function gsl_sf_fermi_dirac_2(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_fermi_dirac_2
     end function gsl_sf_fermi_dirac_2
     function gsl_sf_fermi_dirac_2_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_fermi_dirac_2_e
     end function gsl_sf_fermi_dirac_2_e
     function gsl_sf_fermi_dirac_int(i, x) bind(c)
       import
       integer(c_int), value :: i
       real(c_double), value :: x
       real(c_double) :: gsl_sf_fermi_dirac_int
     end function gsl_sf_fermi_dirac_int
     function gsl_sf_fermi_dirac_int_e(i, x, result) bind(c)
       import
       integer(c_int), value :: i
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_fermi_dirac_int_e
     end function gsl_sf_fermi_dirac_int_e
     function gsl_sf_fermi_dirac_mhalf(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_fermi_dirac_mhalf
     end function gsl_sf_fermi_dirac_mhalf
     function gsl_sf_fermi_dirac_mhalf_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_fermi_dirac_mhalf_e
     end function gsl_sf_fermi_dirac_mhalf_e
     function gsl_sf_fermi_dirac_half(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_fermi_dirac_half
     end function gsl_sf_fermi_dirac_half
     function gsl_sf_fermi_dirac_half_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_fermi_dirac_half_e
     end function gsl_sf_fermi_dirac_half_e
     function gsl_sf_fermi_dirac_3half(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_fermi_dirac_3half
     end function gsl_sf_fermi_dirac_3half
     function gsl_sf_fermi_dirac_3half_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_fermi_dirac_3half_e
     end function gsl_sf_fermi_dirac_3half_e
     function gsl_sf_fermi_dirac_inc_0(x, b) bind(c)
       import
       real(c_double), value :: x, b
       real(c_double) :: gsl_sf_fermi_dirac_inc_0
     end function gsl_sf_fermi_dirac_inc_0
     function gsl_sf_fermi_dirac_inc_0_e(x, b, result) bind(c)
       import
       real(c_double), value :: x, b
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_fermi_dirac_inc_0_e
     end function gsl_sf_fermi_dirac_inc_0_e
     function gsl_sf_gamma(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_gamma
     end function gsl_sf_gamma
     function gsl_sf_gamma_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_gamma_e
     end function gsl_sf_gamma_e
     function gsl_sf_lngamma(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_lngamma
     end function gsl_sf_lngamma
     function gsl_sf_lngamma_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_lngamma_e
     end function gsl_sf_lngamma_e
     function gsl_sf_lngamma_sgn_e(x, result_lg, sgn) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result_lg
       real(c_double), intent(out) :: sgn
       integer(c_int) :: gsl_sf_lngamma_sgn_e
     end function gsl_sf_lngamma_sgn_e
     function gsl_sf_gammastar(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_gammastar
     end function gsl_sf_gammastar
     function gsl_sf_gammastar_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_gammastar_e
     end function gsl_sf_gammastar_e
     function gsl_sf_gammainv(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_gammainv
     end function gsl_sf_gammainv
     function gsl_sf_gammainv_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_gammainv_e
     end function gsl_sf_gammainv_e
     function gsl_sf_lngamma_complex_e(zr, zi, lnr, arg) bind(c)
       import
       real(c_double), value :: zr, zi
       type(gsl_sf_result), intent(out) :: lnr, arg
       integer(c_int) :: gsl_sf_lngamma_complex_e
     end function gsl_sf_lngamma_complex_e
     function gsl_sf_fact(n) bind(c)
       import
       integer(c_int), value :: n
       real(c_double) :: gsl_sf_fact
     end function gsl_sf_fact
     function gsl_sf_fact_e(n, result) bind(c)
       import
       integer(c_int), value :: n
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_fact_e
     end function gsl_sf_fact_e
     function gsl_sf_doublefact(n) bind(c)
       import
       integer(c_int), value :: n
       real(c_double) :: gsl_sf_doublefact
     end function gsl_sf_doublefact
     function gsl_sf_doublefact_e(n, result) bind(c)
       import
       integer(c_int), value :: n
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_doublefact_e
     end function gsl_sf_doublefact_e
     function gsl_sf_lnfact(n) bind(c)
       import
       integer(c_int), value :: n
       real(c_double) :: gsl_sf_lnfact
     end function gsl_sf_lnfact
     function gsl_sf_lnfact_e(n, result) bind(c)
       import
       integer(c_int), value :: n
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_lnfact_e
     end function gsl_sf_lnfact_e
     function gsl_sf_lndoublefact(n) bind(c)
       import
       integer(c_int), value :: n
       real(c_double) :: gsl_sf_lndoublefact
     end function gsl_sf_lndoublefact
     function gsl_sf_lndoublefact_e(n, result) bind(c)
       import
       integer(c_int), value :: n
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_lndoublefact_e
     end function gsl_sf_lndoublefact_e
     function gsl_sf_choose(n, m) bind(c)
       import
       integer(c_int), value :: n, m
       real(c_double) :: gsl_sf_choose
     end function gsl_sf_choose
     function gsl_sf_choose_e(n, m, result) bind(c)
       import
       integer(c_int), value :: n, m
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_choose_e
     end function gsl_sf_choose_e
     function gsl_sf_lnchoose(n, m) bind(c)
       import
       integer(c_int), value :: n, m
       real(c_double) :: gsl_sf_lnchoose
     end function gsl_sf_lnchoose
     function gsl_sf_lnchoose_e(n, m, result) bind(c)
       import
       integer(c_int), value :: n, m
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_lnchoose_e
     end function gsl_sf_lnchoose_e
     function gsl_sf_taylorcoeff(n, x) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_taylorcoeff
     end function gsl_sf_taylorcoeff
     function gsl_sf_taylorcoeff_e(n, x, result) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_taylorcoeff_e
     end function gsl_sf_taylorcoeff_e
     function gsl_sf_poch(a, x) bind(c)
       import
       real(c_double), value :: a, x
       real(c_double) :: gsl_sf_poch
     end function gsl_sf_poch
     function gsl_sf_poch_e(a, x, result) bind(c)
       import
       real(c_double), value :: a, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_poch_e
     end function gsl_sf_poch_e
     function gsl_sf_lnpoch(a, x) bind(c)
       import
       real(c_double), value :: a, x
       real(c_double) :: gsl_sf_lnpoch
     end function gsl_sf_lnpoch
     function gsl_sf_lnpoch_e(a, x, result) bind(c)
       import
       real(c_double), value :: a, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_lnpoch_e
     end function gsl_sf_lnpoch_e
     function gsl_sf_lnpoch_sgn_e(a, x, result_lg, sgn) bind(c)
       import
       real(c_double), value :: a, x
       type(gsl_sf_result), intent(out) :: result_lg
       real(c_double), intent(out) :: sgn
       integer(c_int) :: gsl_sf_lnpoch_sgn_e
     end function gsl_sf_lnpoch_sgn_e
     function gsl_sf_pochrel(a, x) bind(c)
       import
       real(c_double), value :: a, x
       real(c_double) :: gsl_sf_pochrel
     end function gsl_sf_pochrel
     function gsl_sf_pochrel_e(a, x, result) bind(c)
       import
       real(c_double), value :: a, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_pochrel_e
     end function gsl_sf_pochrel_e
     function gsl_sf_gamma_inc(a, x) bind(c)
       import
       real(c_double), value :: a, x
       real(c_double) :: gsl_sf_gamma_inc
     end function gsl_sf_gamma_inc
     function gsl_sf_gamma_inc_e(a, x, result) bind(c)
       import
       real(c_double), value :: a, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_gamma_inc_e
     end function gsl_sf_gamma_inc_e
     function gsl_sf_gamma_inc_q(a, x) bind(c, name='gsl_sf_gamma_inc_Q')
       import
       real(c_double), value :: a, x
       real(c_double) :: gsl_sf_gamma_inc_q
     end function gsl_sf_gamma_inc_q
     function gsl_sf_gamma_inc_q_e(a, x, result) bind(c, name='gsl_sf_gamma_inc_Q_e')
       import
       real(c_double), value :: a, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_gamma_inc_q_e
     end function gsl_sf_gamma_inc_q_e
     function gsl_sf_gamma_inc_p(a, x) bind(c, name='gsl_sf_gamma_inc_P')
       import
       real(c_double), value :: a, x
       real(c_double) :: gsl_sf_gamma_inc_p
     end function gsl_sf_gamma_inc_p
     function gsl_sf_gamma_inc_p_e(a, x, result) bind(c, name='gsl_sf_gamma_inc_P_e')
       import
       real(c_double), value :: a, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_gamma_inc_p_e
     end function gsl_sf_gamma_inc_p_e
     function gsl_sf_beta(a, b) bind(c)
       import
       real(c_double), value :: a, b
       real(c_double) :: gsl_sf_beta
     end function gsl_sf_beta
     function gsl_sf_beta_e(a, b, result) bind(c)
       import
       real(c_double), value :: a, b
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_beta_e
     end function gsl_sf_beta_e
     function gsl_sf_lnbeta(a, b) bind(c)
       import
       real(c_double), value :: a, b
       real(c_double) :: gsl_sf_lnbeta
     end function gsl_sf_lnbeta
     function gsl_sf_lnbeta_e(a, b, result) bind(c)
       import
       real(c_double), value :: a, b
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_lnbeta_e
     end function gsl_sf_lnbeta_e
     function gsl_sf_beta_inc(a, b, x) bind(c)
       import
       real(c_double), value :: a, b, x
       real(c_double) :: gsl_sf_beta_inc
     end function gsl_sf_beta_inc
     function gsl_sf_beta_inc_e(a, b, x, result) bind(c)
       import
       real(c_double), value :: a, b, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_beta_inc_e
     end function gsl_sf_beta_inc_e
     function gsl_sf_gegenpoly_1(lambda, x) bind(c)
       import
       real(c_double), value :: lambda, x
       real(c_double) :: gsl_sf_gegenpoly_1
     end function gsl_sf_gegenpoly_1
     function gsl_sf_gegenpoly_1_e(lambda, x, result) bind(c)
       import
       real(c_double), value :: lambda, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_gegenpoly_1_e
     end function gsl_sf_gegenpoly_1_e
     function gsl_sf_gegenpoly_2(lambda, x) bind(c)
       import
       real(c_double), value :: lambda, x
       real(c_double) :: gsl_sf_gegenpoly_2
     end function gsl_sf_gegenpoly_2
     function gsl_sf_gegenpoly_2_e(lambda, x, result) bind(c)
       import
       real(c_double), value :: lambda, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_gegenpoly_2_e
     end function gsl_sf_gegenpoly_2_e
     function gsl_sf_gegenpoly_3(lambda, x) bind(c)
       import
       real(c_double), value :: lambda, x
       real(c_double) :: gsl_sf_gegenpoly_3
     end function gsl_sf_gegenpoly_3
     function gsl_sf_gegenpoly_3_e(lambda, x, result) bind(c)
       import
       real(c_double), value :: lambda, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_gegenpoly_3_e
     end function gsl_sf_gegenpoly_3_e
     function gsl_sf_gegenpoly_n(n, lambda, x) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: lambda, x
       real(c_double) :: gsl_sf_gegenpoly_n
     end function gsl_sf_gegenpoly_n
     function gsl_sf_gegenpoly_n_e(n, lambda, x, result) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: lambda, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_gegenpoly_n_e
     end function gsl_sf_gegenpoly_n_e
     function gsl_sf_gegenpoly_array(nmax, lambda, x, result_array) bind(c)
       import
       integer(c_int), value :: nmax
       real(c_double), value :: lambda, x
       real(c_double), dimension(*), intent(out) :: result_array
       integer(c_int) :: gsl_sf_gegenpoly_array
     end function gsl_sf_gegenpoly_array
     function gsl_sf_hyperg_0f1(c, x) bind(c, name='gsl_sf_hyperg_0F1')
       import
       real(c_double), value :: c, x
       real(c_double) :: gsl_sf_hyperg_0F1
     end function gsl_sf_hyperg_0f1
     function gsl_sf_hyperg_0f1_e(c, x, result) bind(c, name='gsl_sf_hyperg_0F1_e')
       import
       real(c_double), value :: c, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_0F1_e
     end function gsl_sf_hyperg_0f1_e
     function gsl_sf_hyperg_1f1_int(m, n, x) bind(c, name='gsl_sf_hyperg_1F1_int')
       import
       integer(c_int), value :: m, n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_hyperg_1F1_int
     end function gsl_sf_hyperg_1f1_int
     function gsl_sf_hyperg_1f1_int_e(m, n, x, result) bind(c, name='gsl_sf_hyperg_1F1_int_e')
       import
       integer(c_int), value :: m, n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_1F1_int_e
     end function gsl_sf_hyperg_1f1_int_e
     function gsl_sf_hyperg_1f1(a, b, x) bind(c, name='gsl_sf_hyperg_1F1')
       import
       real(c_double), value :: a, b, x
       real(c_double) :: gsl_sf_hyperg_1f1
     end function gsl_sf_hyperg_1f1
     function gsl_sf_hyperg_1f1_e(a, b, x, result) bind(c, name='gsl_sf_hyperg_1F1_e')
       import
       real(c_double), value :: a, b, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_1f1_e
     end function gsl_sf_hyperg_1f1_e
     function gsl_sf_hyperg_u_int(m, n, x) bind(c, name='gsl_sf_hyperg_U_int')
       import
       integer(c_int), value :: m, n
       real(c_double), value :: x
       real(c_double) :: gsl_sf_hyperg_u_int
     end function gsl_sf_hyperg_u_int
     function gsl_sf_hyperg_u_int_e(m, n, x, result) bind(c, name='gsl_sf_hyperg_U_int_e')
       import
       integer(c_int), value :: m, n
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_u_int_e
     end function gsl_sf_hyperg_u_int_e
     function gsl_sf_hyperg_u_int_e10_e(m, n, x, result) bind(c, name='gsl_sf_hyperg_U_int_e10_e')
       import
       integer(c_int), value :: m, n
       real(c_double), value :: x
       type(gsl_sf_result_e10), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_u_int_e10_e
     end function gsl_sf_hyperg_u_int_e10_e
     function gsl_sf_hyperg_u(a, b, x) bind(c, name='gsl_sf_hyperg_U')
       import
       real(c_double), value :: a, b, x
       real(c_double) :: gsl_sf_hyperg_u
     end function gsl_sf_hyperg_u
     function gsl_sf_hyperg_u_e(a, b, x, result) bind(c, name='gsl_sf_hyperg_U_e')
       import
       real(c_double), value :: a, b, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_u_e
     end function gsl_sf_hyperg_u_e
     function gsl_sf_hyperg_u_e10_e(a, b, x, result) bind(c, name='gsl_sf_hyperg_U_e10_e')
       import
       real(c_double), value :: a, b, x
       type(gsl_sf_result_e10), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_u_e10_e
     end function gsl_sf_hyperg_u_e10_e
     function gsl_sf_hyperg_2f1(a, b, c, x) bind(c, name='gsl_sf_hyperg_2F1')
       import
       real(c_double), value :: a, b, c, x
       real(c_double) :: gsl_sf_hyperg_2f1
     end function gsl_sf_hyperg_2f1
     function gsl_sf_hyperg_2f1_e(a, b, c, x, result) bind(c, name='gsl_sf_hyperg_2F1_e')
       import
       real(c_double), value :: a, b, c, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_2f1_e
     end function gsl_sf_hyperg_2f1_e
     function gsl_sf_hyperg_2f1_conj(ar, ai, c, x) bind(c, name='gsl_sf_hyperg_2F1_conj')
       import
       real(c_double), value :: ar, ai, c, x
       real(c_double) :: gsl_sf_hyperg_2f1_conj
     end function gsl_sf_hyperg_2f1_conj
     function gsl_sf_hyperg_2f1_conj_e(ar, ai, c, x, result) &
          bind(c, name='gsl_sf_hyperg_2F1_conj_e')
       import
       real(c_double), value :: ar, ai, c, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_2f1_conj_e
     end function gsl_sf_hyperg_2f1_conj_e
     function gsl_sf_hyperg_2f1_renorm(a, b, c, x) bind(c, name='gsl_sf_hyperg_2F1_renorm')
       import
       real(c_double), value :: a, b, c, x
       real(c_double) :: gsl_sf_hyperg_2f1_renorm
     end function gsl_sf_hyperg_2f1_renorm
     function gsl_sf_hyperg_2f1_renorm_e(a, b, c, x, result) &
          bind(c, name='gsl_sf_hyperg_2F1_renorm_e')
       import
       real(c_double), value :: a, b, c, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_2f1_renorm_e
     end function gsl_sf_hyperg_2f1_renorm_e
     function gsl_sf_hyperg_2f1_conj_renorm(ar, ai, c, x) &
          bind(c, name='gsl_sf_hyperg_2F1_conj_renorm')
       import
       real(c_double), value :: ar, ai, c, x
       real(c_double) :: gsl_sf_hyperg_2f1_conj_renorm
     end function gsl_sf_hyperg_2f1_conj_renorm
     function gsl_sf_hyperg_2f1_conj_renorm_e(ar, ai, c, x, result) &
          bind(c, name='gsl_sf_hyperg_2F1_conj_renorm_e')
       import
       real(c_double), value :: ar, ai, c, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_2f1_conj_renorm_e
     end function gsl_sf_hyperg_2f1_conj_renorm_e
     function gsl_sf_hyperg_2f0(a, b, x) bind(c, name='gsl_sf_hyperg_2F0')
       import
       real(c_double), value :: a, b, x
       real(c_double) :: gsl_sf_hyperg_2f0
     end function gsl_sf_hyperg_2f0
     function gsl_sf_hyperg_2f0_e(a, b, x, result) bind(c, name='gsl_sf_hyperg_2F0_e')
       import
       real(c_double), value :: a, b, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hyperg_2f0_e
     end function gsl_sf_hyperg_2f0_e
     function gsl_sf_laguerre_1(a, x) bind(c)
       import
       real(c_double), value :: a, x
       real(c_double) :: gsl_sf_laguerre_1
     end function gsl_sf_laguerre_1
     function gsl_sf_laguerre_1_e(a, x, result) bind(c)
       import
       real(c_double), value :: a, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_laguerre_1_e
     end function gsl_sf_laguerre_1_e
     function gsl_sf_laguerre_2(a, x) bind(c)
       import
       real(c_double), value :: a, x
       real(c_double) :: gsl_sf_laguerre_2
     end function gsl_sf_laguerre_2
     function gsl_sf_laguerre_2_e(a, x, result) bind(c)
       import
       real(c_double), value :: a, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_laguerre_2_e
     end function gsl_sf_laguerre_2_e
     function gsl_sf_laguerre_3(a, x) bind(c)
       import
       real(c_double), value :: a, x
       real(c_double) :: gsl_sf_laguerre_3
     end function gsl_sf_laguerre_3
     function gsl_sf_laguerre_3_e(a, x, result) bind(c)
       import
       real(c_double), value :: a, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_laguerre_3_e
     end function gsl_sf_laguerre_3_e
     function gsl_sf_laguerre_n(n, a, x) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: a, x
       real(c_double) :: gsl_sf_laguerre_n
     end function gsl_sf_laguerre_n
     function gsl_sf_laguerre_n_e(n, a, x, result) bind(c)
       import
       integer(c_int), value :: n
       real(c_double), value :: a, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_laguerre_n_e
     end function gsl_sf_laguerre_n_e
     function gsl_sf_lambert_w0(x) bind(c, name='gsl_sf_lambert_W0')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_lambert_w0
     end function gsl_sf_lambert_w0
     function gsl_sf_lambert_w0_e(x, result) bind(c, name='gsl_sf_lambert_W0_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_lambert_w0_e
     end function gsl_sf_lambert_w0_e
     function gsl_sf_lambert_wm1(x) bind(c, name='gsl_sf_lambert_Wm1')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_lambert_wm1
     end function gsl_sf_lambert_wm1
     function gsl_sf_lambert_wm1_e(x, result) bind(c, name='gsl_sf_lambert_Wm1_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_lambert_wm1_e
     end function gsl_sf_lambert_wm1_e
     function gsl_sf_legendre_p1(x) bind(c, name='gsl_sf_legendre_P1')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_legendre_p1
     end function gsl_sf_legendre_p1
     function gsl_sf_legendre_p1_e(x, result) bind(c, name='gsl_sf_legendre_P1_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_p1_e
     end function gsl_sf_legendre_p1_e
     function gsl_sf_legendre_p2(x) bind(c, name='gsl_sf_legendre_P2')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_legendre_p2
     end function gsl_sf_legendre_p2
     function gsl_sf_legendre_p2_e(x, result) bind(c, name='gsl_sf_legendre_P2_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_p2_e
     end function gsl_sf_legendre_p2_e
     function gsl_sf_legendre_p3(x) bind(c, name='gsl_sf_legendre_P3')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_legendre_p3
     end function gsl_sf_legendre_p3
     function gsl_sf_legendre_p3_e(x, result) bind(c, name='gsl_sf_legendre_P3_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_p3_e
     end function gsl_sf_legendre_p3_e
     function gsl_sf_legendre_pl(l, x) bind(c, name='gsl_sf_legendre_Pl')
       import
       integer(c_int), value :: l
       real(c_double), value :: x
       real(c_double) :: gsl_sf_legendre_pl
     end function gsl_sf_legendre_pl
     function gsl_sf_legendre_pl_e(l, x, result) bind(c, name='gsl_sf_legendre_Pl_e')
       import
       integer(c_int), value :: l
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_pl_e
     end function gsl_sf_legendre_pl_e
     function gsl_sf_legendre_pl_array(lmax, x, res_arr) bind(c, name='gsl_sf_legendre_Pl_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: res_arr
       integer(c_int) :: gsl_sf_legendre_pl_array
     end function gsl_sf_legendre_pl_array
     function gsl_sf_legendre_pl_deriv_array(lmax, x, res_arr, der_arr) &
          bind(c, name='gsl_sf_legendre_Pl_deriv_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: res_arr, der_arr
       integer(c_int) :: gsl_sf_legendre_pl_deriv_array
     end function gsl_sf_legendre_pl_deriv_array
     function gsl_sf_legendre_q0(x) bind(c, name='gsl_sf_legendre_Q0')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_legendre_q0
     end function gsl_sf_legendre_q0
     function gsl_sf_legendre_q0_e(x, result) bind(c, name='gsl_sf_legendre_Q0_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_q0_e
     end function gsl_sf_legendre_q0_e
     function gsl_sf_legendre_q1(x) bind(c, name='gsl_sf_legendre_Q1')
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_legendre_q1
     end function gsl_sf_legendre_q1
     function gsl_sf_legendre_q1_e(x, result) bind(c, name='gsl_sf_legendre_Q1_e')
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_q1_e
     end function gsl_sf_legendre_q1_e
     function gsl_sf_legendre_ql(l, x) bind(c, name='gsl_sf_legendre_Ql')
       import
       integer(c_int), value :: l
       real(c_double), value :: x
       real(c_double) :: gsl_sf_legendre_ql
     end function gsl_sf_legendre_ql
     function gsl_sf_legendre_ql_e(l, x, result) bind(c, name='gsl_sf_legendre_Ql_e')
       import
       integer(c_int), value :: l
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_ql_e
     end function gsl_sf_legendre_ql_e
     function gsl_sf_legendre_plm(l, m, x) bind(c, name='gsl_sf_legendre_Plm')
       import
       integer(c_int), value :: l, m
       real(c_double), value :: x
       real(c_double) :: gsl_sf_legendre_plm
     end function gsl_sf_legendre_plm
     function gsl_sf_legendre_plm_e(l, m, x, result) bind(c, name='gsl_sf_legendre_Plm_e')
       import
       integer(c_int), value :: l, m
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_plm_e
     end function gsl_sf_legendre_plm_e
     function gsl_sf_legendre_plm_array(lmax, m, x, res_arr) &
          bind(c, name='gsl_sf_legendre_Plm_array')
       import
       integer(c_int), value :: lmax, m
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: res_arr
       integer(c_int) :: gsl_sf_legendre_plm_array
     end function gsl_sf_legendre_plm_array
     function gsl_sf_legendre_plm_deriv_array(lmax, m, x, res_arr, der_arr) &
          bind(c, name='gsl_sf_legendre_Plm_deriv_array')
       import
       integer(c_int), value :: lmax, m
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: res_arr, der_arr
       integer(c_int) :: gsl_sf_legendre_plm_deriv_array
     end function gsl_sf_legendre_plm_deriv_array
     function gsl_sf_legendre_sphplm(l, m, x) bind(c, name='gsl_sf_legendre_sphPlm')
       import
       integer(c_int), value :: l, m
       real(c_double), value :: x
       real(c_double) :: gsl_sf_legendre_sphplm
     end function gsl_sf_legendre_sphplm
     function gsl_sf_legendre_sphplm_e(l, m, x, result) bind(c, name='gsl_sf_legendre_sphPlm_e')
       import
       integer(c_int), value :: l, m
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_sphplm_e
     end function gsl_sf_legendre_sphplm_e
     function gsl_sf_legendre_sphplm_array(lmax, m, x, res_arr) &
          bind(c, name='gsl_sf_legendre_sphPlm_array')
       import
       integer(c_int), value :: lmax, m
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: res_arr
       integer(c_int) :: gsl_sf_legendre_sphplm_array
     end function gsl_sf_legendre_sphplm_array
     function gsl_sf_legendre_sphplm_deriv_array(lmax, m, x, res_arr, der_arr) &
          bind(c, name='gsl_sf_legendre_sphPlm_deriv_array')
       import
       integer(c_int), value :: lmax, m
       real(c_double), value :: x
       real(c_double), dimension(*), intent(out) :: res_arr, der_arr
       integer(c_int) :: gsl_sf_legendre_sphplm_deriv_array
     end function gsl_sf_legendre_sphplm_deriv_array
     function gsl_sf_legendre_array_size(lmax, m) bind(c)
       import
       integer(c_int), value :: lmax, m
       integer(c_int) :: gsl_sf_legendre_array_size
     end function gsl_sf_legendre_array_size
     function gsl_sf_conicalp_half(lambda, x) bind(c, name='gsl_sf_conicalP_half')
       import
       real(c_double), value :: lambda, x
       real(c_double) :: gsl_sf_conicalp_half
     end function gsl_sf_conicalp_half
     function gsl_sf_conicalp_half_e(lambda, x, result) bind(c, name='gsl_sf_conicalP_half_e')
       import
       real(c_double), value :: lambda, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_conicalp_half_e
     end function gsl_sf_conicalp_half_e
     function gsl_sf_conicalp_mhalf(lambda, x) bind(c, name='gsl_sf_conicalP_mhalf')
       import
       real(c_double), value :: lambda, x
       real(c_double) :: gsl_sf_conicalp_mhalf
     end function gsl_sf_conicalp_mhalf
     function gsl_sf_conicalp_mhalf_e(lambda, x, result) bind(c, name='gsl_sf_conicalP_mhalf_e')
       import
       real(c_double), value :: lambda, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_conicalp_mhalf_e
     end function gsl_sf_conicalp_mhalf_e
     function gsl_sf_conicalp_0(lambda, x) bind(c, name='gsl_sf_conicalP_0')
       import
       real(c_double), value :: lambda, x
       real(c_double) :: gsl_sf_conicalp_0
     end function gsl_sf_conicalp_0
     function gsl_sf_conicalp_0_e(lambda, x, result) bind(c, name='gsl_sf_conicalP_0_e')
       import
       real(c_double), value :: lambda, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_conicalp_0_e
     end function gsl_sf_conicalp_0_e
     function gsl_sf_conicalp_1(lambda, x) bind(c, name='gsl_sf_conicalP_1')
       import
       real(c_double), value :: lambda, x
       real(c_double) :: gsl_sf_conicalp_1
     end function gsl_sf_conicalp_1
     function gsl_sf_conicalp_1_e(lambda, x, result) bind(c, name='gsl_sf_conicalP_1_e')
       import
       real(c_double), value :: lambda, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_conicalp_1_e
     end function gsl_sf_conicalp_1_e
     function gsl_sf_conicalp_sph_reg(l, lambda, x) bind(c, name='gsl_sf_conicalP_sph_reg')
       import
       integer(c_int), value :: l
       real(c_double), value :: lambda, x
       real(c_double) :: gsl_sf_conicalp_sph_reg
     end function gsl_sf_conicalp_sph_reg
     function gsl_sf_conicalp_sph_reg_e(l, lambda, x, result) &
          bind(c, name='gsl_sf_conicalP_sph_reg_e')
       import
       integer(c_int), value :: l
       real(c_double), value :: lambda, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_conicalp_sph_reg_e
     end function gsl_sf_conicalp_sph_reg_e
     function gsl_sf_conicalp_cyl_reg(l, lambda, x) bind(c, name='gsl_sf_conicalP_cyl_reg')
       import
       integer(c_int), value :: l
       real(c_double), value :: lambda, x
       real(c_double) :: gsl_sf_conicalp_cyl_reg
     end function gsl_sf_conicalp_cyl_reg
     function gsl_sf_conicalp_cyl_reg_e(l, lambda, x, result) &
          bind(c, name='gsl_sf_conicalP_cyl_reg_e')
       import
       integer(c_int), value :: l
       real(c_double), value :: lambda, x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_conicalp_cyl_reg_e
     end function gsl_sf_conicalp_cyl_reg_e
     function gsl_sf_legendre_h3d_0(lambda, eta) bind(c, name='gsl_sf_legendre_H3d_0')
       import
       real(c_double), value :: lambda, eta
       real(c_double) :: gsl_sf_legendre_h3d_0
     end function gsl_sf_legendre_h3d_0
     function gsl_sf_legendre_h3d_0_e(lambda, eta, result) bind(c, name='gsl_sf_legendre_H3d_0_e')
       import
       real(c_double), value :: lambda, eta
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_h3d_0_e
     end function gsl_sf_legendre_h3d_0_e
     function gsl_sf_legendre_h3d_1(lambda, eta) bind(c, name='gsl_sf_legendre_H3d_1')
       import
       real(c_double), value :: lambda, eta
       real(c_double) :: gsl_sf_legendre_h3d_1
     end function gsl_sf_legendre_h3d_1
     function gsl_sf_legendre_h3d_1_e(lambda, eta, result) bind(c, name='gsl_sf_legendre_H3d_1_e')
       import
       real(c_double), value :: lambda, eta
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_h3d_1_e
     end function gsl_sf_legendre_h3d_1_e
     function gsl_sf_legendre_h3d(l, lambda, eta) bind(c, name='gsl_sf_legendre_H3d')
       import
       integer(c_int), value :: l
       real(c_double), value :: lambda, eta
       real(c_double) :: gsl_sf_legendre_h3d
     end function gsl_sf_legendre_h3d
     function gsl_sf_legendre_h3d_e(l, lambda, eta, result) bind(c, name='gsl_sf_legendre_H3d_e')
       import
       integer(c_int), value :: l
       real(c_double), value :: lambda, eta
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_legendre_h3d_e
     end function gsl_sf_legendre_h3d_e
     function gsl_sf_legendre_h3d_array(lmax, lambda, eta, res_arr) &
          bind(c, name='gsl_sf_legendre_H3d_array')
       import
       integer(c_int), value :: lmax
       real(c_double), value :: lambda, eta
       real(c_double), dimension(*), intent(out) :: res_arr
       integer(c_int) :: gsl_sf_legendre_h3d_array
     end function gsl_sf_legendre_h3d_array
     function gsl_sf_log(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_log
     end function gsl_sf_log
     function gsl_sf_log_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_log_e
     end function gsl_sf_log_e
     function gsl_sf_log_abs(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_log_abs
     end function gsl_sf_log_abs
     function gsl_sf_log_abs_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_log_abs_e
     end function gsl_sf_log_abs_e
     function gsl_sf_complex_log_e(zr, zi, lnr, theta) bind(c)
       import
       real(c_double), value :: zr, zi
       type(gsl_sf_result), intent(out) :: lnr, theta
       integer(c_int) :: gsl_sf_complex_log_e
     end function gsl_sf_complex_log_e
     function gsl_sf_log_1plusx(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_log_1plusx
     end function gsl_sf_log_1plusx
     function gsl_sf_log_1plusx_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_log_1plusx_e
     end function gsl_sf_log_1plusx_e
     function gsl_sf_log_1plusx_mx(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_log_1plusx_mx
     end function gsl_sf_log_1plusx_mx
     function gsl_sf_log_1plusx_mx_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_log_1plusx_mx_e
     end function gsl_sf_log_1plusx_mx_e
     function gsl_sf_psi_int(n) bind(c)
       import
       integer(c_int), value :: n
       real(c_double) :: gsl_sf_psi_int
     end function gsl_sf_psi_int
     function gsl_sf_psi_int_e(n, result) bind(c)
       import
       integer(c_int), value :: n
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_psi_int_e
     end function gsl_sf_psi_int_e
     function gsl_sf_psi(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_psi
     end function gsl_sf_psi
     function gsl_sf_psi_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_psi_e
     end function gsl_sf_psi_e
     function gsl_sf_psi_1_int(n) bind(c)
       import
       integer(c_int), value :: n
       real(c_double) :: gsl_sf_psi_1_int
     end function gsl_sf_psi_1_int
     function gsl_sf_psi_1_int_e(n, result) bind(c)
       import
       integer(c_int), value :: n
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_psi_1_int_e
     end function gsl_sf_psi_1_int_e
     function gsl_sf_psi_1(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_psi_1
     end function gsl_sf_psi_1
     function gsl_sf_psi_1_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_psi_1_e
     end function gsl_sf_psi_1_e
     function gsl_sf_psi_n(m, x) bind(c)
       import
       integer(c_int), value :: m
       real(c_double), value :: x
       real(c_double) :: gsl_sf_psi_n
     end function gsl_sf_psi_n
     function gsl_sf_psi_n_e(m, x, result) bind(c)
       import
       integer(c_int), value :: m
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_psi_n_e
     end function gsl_sf_psi_n_e
     function gsl_sf_psi_1piy(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_psi_1piy
     end function gsl_sf_psi_1piy
     function gsl_sf_psi_1piy_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_psi_1piy_e
     end function gsl_sf_psi_1piy_e
     function gsl_sf_synchrotron_1(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_synchrotron_1
     end function gsl_sf_synchrotron_1
     function gsl_sf_synchrotron_1_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_synchrotron_1_e
     end function gsl_sf_synchrotron_1_e
     function gsl_sf_synchrotron_2(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_synchrotron_2
     end function gsl_sf_synchrotron_2
     function gsl_sf_synchrotron_2_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_synchrotron_2_e
     end function gsl_sf_synchrotron_2_e
     function gsl_sf_transport_2(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_transport_2
     end function gsl_sf_transport_2
     function gsl_sf_transport_2_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_transport_2_e
     end function gsl_sf_transport_2_e
     function gsl_sf_transport_3(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_transport_3
     end function gsl_sf_transport_3
     function gsl_sf_transport_3_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_transport_3_e
     end function gsl_sf_transport_3_e
     function gsl_sf_transport_4(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_transport_4
     end function gsl_sf_transport_4
     function gsl_sf_transport_4_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_transport_4_e
     end function gsl_sf_transport_4_e
     function gsl_sf_transport_5(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_transport_5
     end function gsl_sf_transport_5
     function gsl_sf_transport_5_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_transport_5_e
     end function gsl_sf_transport_5_e
     function gsl_sf_hypot(x, y) bind(c)
       import
       real(c_double), value :: x, y
       real(c_double) :: gsl_sf_hypot
     end function gsl_sf_hypot
     function gsl_sf_hypot_e(x, y, result) bind(c)
       import
       real(c_double), value :: x, y
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hypot_e
     end function gsl_sf_hypot_e
     function gsl_sf_sinc(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_sinc
     end function gsl_sf_sinc
     function gsl_sf_sinc_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_sinc_e
     end function gsl_sf_sinc_e
     function gsl_sf_complex_sin_e(zr, zi, szr, szi) bind(c)
       import
       real(c_double), value :: zr, zi
       type(gsl_sf_result), intent(out) :: szr, szi
       integer(c_int) :: gsl_sf_complex_sin_e
     end function gsl_sf_complex_sin_e
     function gsl_sf_complex_cos_e(zr, zi, czr, czi) bind(c)
       import
       real(c_double), value :: zr, zi
       type(gsl_sf_result), intent(out) :: czr, czi
       integer(c_int) :: gsl_sf_complex_cos_e
     end function gsl_sf_complex_cos_e
     function gsl_sf_complex_logsin_e(zr, zi, lszr, lszi) bind(c)
       import
       real(c_double), value :: zr, zi
       type(gsl_sf_result), intent(out) :: lszr, lszi
       integer(c_int) :: gsl_sf_complex_logsin_e
     end function gsl_sf_complex_logsin_e
     function gsl_sf_lnsinh(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_lnsinh
     end function gsl_sf_lnsinh
     function gsl_sf_lnsinh_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_lnsinh_e
     end function gsl_sf_lnsinh_e
     function gsl_sf_lncosh(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_lncosh
     end function gsl_sf_lncosh
     function gsl_sf_lncosh_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_lncosh_e
     end function gsl_sf_lncosh_e
     function gsl_sf_polar_to_rect(r, theta, x, y) bind(c)
       import
       real(c_double), value :: r, theta
       type(gsl_sf_result), intent(out) :: x, y
       integer(c_int) :: gsl_sf_polar_to_rect
     end function gsl_sf_polar_to_rect
     function gsl_sf_rect_to_polar(x, y, r, theta) bind(c)
       import
       real(c_double), value :: x, y
       type(gsl_sf_result), intent(out) :: r, theta
       integer(c_int) :: gsl_sf_rect_to_polar
     end function gsl_sf_rect_to_polar
     function gsl_sf_angle_restrict_symm(theta) bind(c)
       import
       real(c_double), value :: theta
       real(c_double) :: gsl_sf_angle_restrict_symm
     end function gsl_sf_angle_restrict_symm
     function gsl_sf_angle_restrict_symm_e(theta) bind(c)
       import
       real(c_double), intent(inout) :: theta
       integer(c_int) :: gsl_sf_angle_restrict_symm_e
     end function gsl_sf_angle_restrict_symm_e
     function gsl_sf_angle_restrict_pos(theta) bind(c)
       import
       real(c_double), value :: theta
       real(c_double) :: gsl_sf_angle_restrict_pos
     end function gsl_sf_angle_restrict_pos
     function gsl_sf_angle_restrict_pos_e(theta) bind(c)
       import
       real(c_double), intent(inout) :: theta
       integer(c_int) :: gsl_sf_angle_restrict_pos_e
     end function gsl_sf_angle_restrict_pos_e
     function gsl_sf_sin_err_e(x, dx, result) bind(c)
       import
       real(c_double), value :: x, dx
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_sin_err_e
     end function gsl_sf_sin_err_e
     function gsl_sf_cos_err_e(x, dx, result) bind(c)
       import
       real(c_double), value :: x, dx
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_cos_err_e
     end function gsl_sf_cos_err_e
     function gsl_sf_zeta_int(n) bind(c)
       import
       integer(c_int), value :: n
       real(c_double) :: gsl_sf_zeta_int
     end function gsl_sf_zeta_int
     function gsl_sf_zeta_int_e(n, result) bind(c)
       import
       integer(c_int), value :: n
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_zeta_int_e
     end function gsl_sf_zeta_int_e
     function gsl_sf_zeta(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_zeta
     end function gsl_sf_zeta
     function gsl_sf_zeta_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_zeta_e
     end function gsl_sf_zeta_e
     function gsl_sf_zetam1_int(n) bind(c)
       import
       integer(c_int), value :: n
       real(c_double) :: gsl_sf_zetam1_int
     end function gsl_sf_zetam1_int
     function gsl_sf_zetam1_int_e(n, result) bind(c)
       import
       integer(c_int), value :: n
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_zetam1_int_e
     end function gsl_sf_zetam1_int_e
     function gsl_sf_zetam1(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_zetam1
     end function gsl_sf_zetam1
     function gsl_sf_zetam1_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_zetam1_e
     end function gsl_sf_zetam1_e
     function gsl_sf_hzeta(s, q) bind(c)
       import
       real(c_double), value :: s, q
       real(c_double) :: gsl_sf_hzeta
     end function gsl_sf_hzeta
     function gsl_sf_hzeta_e(s, q, result) bind(c)
       import
       real(c_double), value :: s, q
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_hzeta_e
     end function gsl_sf_hzeta_e
     function gsl_sf_eta_int(n) bind(c)
       import
       integer(c_int), value :: n
       real(c_double) :: gsl_sf_eta_int
     end function gsl_sf_eta_int
     function gsl_sf_eta_int_e(n, result) bind(c)
       import
       integer(c_int), value :: n
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_eta_int_e
     end function gsl_sf_eta_int_e
     function gsl_sf_eta(x) bind(c)
       import
       real(c_double), value :: x
       real(c_double) :: gsl_sf_eta
     end function gsl_sf_eta
     function gsl_sf_eta_e(x, result) bind(c)
       import
       real(c_double), value :: x
       type(gsl_sf_result), intent(out) :: result
       integer(c_int) :: gsl_sf_eta_e
     end function gsl_sf_eta_e
!-*-f90-*-
!
!  Interfaces: Array support
!
  function gsl_vector_get(v, i) bind(c)
    import
    type(c_ptr), value :: v
    integer(c_size_t), value :: i
    real(c_double) :: gsl_vector_get
  end function gsl_vector_get
  function gsl_vector_ptr(v, i) bind(c)
    import
    type(c_ptr), value :: v
    integer(c_size_t), value :: i
    type(c_ptr) :: gsl_vector_ptr
  end function gsl_vector_ptr
  function gsl_vector_complex_get(v, i) bind(c)
    import
    type(c_ptr), value :: v
    integer(c_size_t), value :: i
    type(gsl_complex) :: gsl_vector_complex_get
  end function gsl_vector_complex_get
  function gsl_vector_complex_ptr(v, i) bind(c)
    import
    type(c_ptr), value :: v
    integer(c_size_t), value :: i
    type(c_ptr) :: gsl_vector_complex_ptr
  end function gsl_vector_complex_ptr
  function gsl_matrix_get(v, j, i) bind(c)
    import
    type(c_ptr), value :: v
    integer(c_size_t), value :: j, i
    real(c_double) :: gsl_matrix_get
  end function gsl_matrix_get
  function gsl_matrix_complex_get(v, j, i) bind(c)
    import
    type(c_ptr), value :: v
    integer(c_size_t), value :: j, i
    type(gsl_complex) :: gsl_matrix_complex_get
  end function gsl_matrix_complex_get
!
! auxiliary functions within FGSL only
!
  function fgsl_aux_vector_double_init() bind(c)
    import
    type(c_ptr) :: fgsl_aux_vector_double_init
  end function fgsl_aux_vector_double_init
  subroutine fgsl_aux_vector_double_free(v) bind(c)
    import
    type(c_ptr), value :: v
  end subroutine fgsl_aux_vector_double_free
  function fgsl_aux_vector_double_align(a, len, fvec, size, offset, stride) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: a
    type(c_ptr), value :: fvec
    integer(c_size_t), value :: len, size, offset, stride
    integer(c_int) :: fgsl_aux_vector_double_align
  end function fgsl_aux_vector_double_align
  function fgsl_aux_vector_double_size(fvec) bind(c)
    import
    type(c_ptr), value :: fvec
    integer(c_size_t) fgsl_aux_vector_double_size
  end function fgsl_aux_vector_double_size
  function fgsl_aux_vector_double_stride(fvec) bind(c)
    import
    type(c_ptr), value :: fvec
    integer(c_size_t) fgsl_aux_vector_double_stride
  end function fgsl_aux_vector_double_stride
  function fgsl_aux_matrix_double_init() bind(c)
    import
    type(c_ptr) :: fgsl_aux_matrix_double_init
  end function fgsl_aux_matrix_double_init
  subroutine fgsl_aux_matrix_double_free(v) bind(c)
    import
    type(c_ptr), value :: v
  end subroutine fgsl_aux_matrix_double_free
  function fgsl_aux_matrix_double_align(a, lda, n, m, fvec) bind(c)
    import
    integer(c_size_t), value :: lda, n, m
    type(c_ptr), value :: a
!    real(c_double), dimension(lda, *), intent(in) :: a
    type(c_ptr), value :: fvec
    integer(c_int) :: fgsl_aux_matrix_double_align
  end function fgsl_aux_matrix_double_align
  subroutine fgsl_aux_matrix_double_size(fmat, lda, m, n) bind(c)
    import
    type(c_ptr), value :: fmat
    integer(c_size_t), intent(out) :: lda, m, n
  end subroutine fgsl_aux_matrix_double_size
  function gsl_matrix_ptr(m, i, j) bind(c)
    import
    type(c_ptr), value :: m
    integer(c_size_t), value :: i, j
    type(c_ptr) :: gsl_matrix_ptr
  end function gsl_matrix_ptr
  function gsl_aux_sizeof_vector() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_vector
  end function gsl_aux_sizeof_vector
  function gsl_aux_sizeof_matrix() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_matrix
  end function gsl_aux_sizeof_matrix
!
! complex variants
!
  function fgsl_aux_vector_complex_init() bind(c)
    import
    type(c_ptr) :: fgsl_aux_vector_complex_init
  end function fgsl_aux_vector_complex_init
  subroutine fgsl_aux_vector_complex_free(v) bind(c)
    import
    type(c_ptr), value :: v
  end subroutine fgsl_aux_vector_complex_free
  function fgsl_aux_vector_complex_align(a, len, fvec, size, offset, stride) bind(c)
    import
!    complex(c_double), dimension(*), intent(in) :: a
    type(c_ptr), value :: a
    type(c_ptr), value :: fvec
    integer(c_size_t), value :: len, size, offset, stride
    integer(c_int) :: fgsl_aux_vector_complex_align
  end function fgsl_aux_vector_complex_align
  function fgsl_aux_vector_complex_size(fvec) bind(c)
    import
    type(c_ptr), value :: fvec
    integer(c_size_t) fgsl_aux_vector_complex_size
  end function fgsl_aux_vector_complex_size
  function fgsl_aux_vector_complex_stride(fvec) bind(c)
    import
    type(c_ptr), value :: fvec
    integer(c_size_t) fgsl_aux_vector_complex_stride
  end function fgsl_aux_vector_complex_stride
  function fgsl_aux_matrix_complex_init() bind(c)
    import
    type(c_ptr) :: fgsl_aux_matrix_complex_init
  end function fgsl_aux_matrix_complex_init
  subroutine fgsl_aux_matrix_complex_free(v) bind(c)
    import
    type(c_ptr), value :: v
  end subroutine fgsl_aux_matrix_complex_free
  function fgsl_aux_matrix_complex_align(a, lda, n, m, fvec) bind(c)
    import
    integer(c_size_t), value :: lda, n, m
!    complex(c_double), dimension(lda, *), intent(in) :: a
    type(c_ptr), value :: a
    type(c_ptr), value :: fvec
    integer(c_int) :: fgsl_aux_matrix_complex_align
  end function fgsl_aux_matrix_complex_align
  subroutine fgsl_aux_matrix_complex_size(fmat, lda, m, n) bind(c)
    import
    type(c_ptr), value :: fmat
    integer(c_size_t), intent(out) :: lda, m, n
  end subroutine fgsl_aux_matrix_complex_size
  function gsl_matrix_complex_ptr(m, i, j) bind(c)
    import
    type(c_ptr), value :: m
    integer(c_size_t), value :: i, j
    type(c_ptr) :: gsl_matrix_complex_ptr
  end function gsl_matrix_complex_ptr
  function gsl_aux_sizeof_vector_complex() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_vector_complex
  end function gsl_aux_sizeof_vector_complex
  function gsl_aux_sizeof_matrix_complex() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_matrix_complex
  end function gsl_aux_sizeof_matrix_complex
!-*-f90-*-
!
!  Interfaces: Interpolation 
!
     function gsl_interp_eval(interp, xa, ya, x, acc) bind(c)
       import
       type(c_ptr), value :: interp
       real(c_double), dimension(*), intent(in) :: xa, ya
       real(c_double), value :: x
       type(c_ptr), value :: acc
       real(c_double) :: gsl_interp_eval
     end function gsl_interp_eval
     function gsl_interp_eval_e(interp, xa, ya, x, acc, y) bind(c)
       import
       type(c_ptr), value :: interp
       real(c_double), dimension(*), intent(in) :: xa, ya
       real(c_double), value :: x
       real(c_double), intent(out) :: y
       type(c_ptr), value :: acc
       integer(c_int) :: gsl_interp_eval_e
     end function gsl_interp_eval_e
     function gsl_interp_eval_integ(interp, xa, ya, a, b, acc) bind(c)
       import
       type(c_ptr), value :: interp
       real(c_double), dimension(*), intent(in) :: xa, ya
       real(c_double), value :: a, b
       type(c_ptr), value :: acc
       real(c_double) :: gsl_interp_eval_integ
     end function gsl_interp_eval_integ
     function gsl_interp_eval_integ_e(interp, xa, ya, a, b, acc, result) bind(c)
       import
       type(c_ptr), value :: interp
       real(c_double), dimension(*), intent(in) :: xa, ya
       real(c_double), value :: a, b
       real(c_double), intent(out) :: result
       type(c_ptr), value :: acc
       integer(c_int) :: gsl_interp_eval_integ_e
     end function gsl_interp_eval_integ_e
     function gsl_interp_eval_deriv(interp, xa, ya, x, acc) bind(c)
       import
       type(c_ptr), value :: interp
       real(c_double), dimension(*), intent(in) :: xa, ya
       real(c_double), value :: x
       type(c_ptr), value :: acc
       real(c_double) :: gsl_interp_eval_deriv
     end function gsl_interp_eval_deriv
     function gsl_interp_eval_deriv_e(interp, xa, ya, x, acc, y) bind(c)
       import
       type(c_ptr), value :: interp
       real(c_double), dimension(*), intent(in) :: xa, ya
       real(c_double), value :: x
       real(c_double), intent(out) :: y
       type(c_ptr), value :: acc
       integer(c_int) :: gsl_interp_eval_deriv_e
     end function gsl_interp_eval_deriv_e
     function gsl_interp_eval_deriv2(interp, xa, ya, x, acc) bind(c)
       import
       type(c_ptr), value :: interp
       real(c_double), dimension(*), intent(in) :: xa, ya
       real(c_double), value :: x
       type(c_ptr), value :: acc
       real(c_double) :: gsl_interp_eval_deriv2
     end function gsl_interp_eval_deriv2
     function gsl_interp_eval_deriv2_e(interp, xa, ya, x, acc, y) bind(c)
       import
       type(c_ptr), value :: interp
       real(c_double), dimension(*), intent(in) :: xa, ya
       real(c_double), value :: x
       real(c_double), intent(out) :: y
       type(c_ptr), value :: acc
       integer(c_int) :: gsl_interp_eval_deriv2_e
     end function gsl_interp_eval_deriv2_e
     function fgsl_aux_interp_alloc(int_interp) bind(c)
       import
       integer(fgsl_int), value :: int_interp
       type(c_ptr) :: fgsl_aux_interp_alloc
     end function fgsl_aux_interp_alloc
     function gsl_interp_alloc(interp_type, size) bind(c)
       import
       type(c_ptr), value :: interp_type
       integer(c_size_t), value :: size
       type(c_ptr) :: gsl_interp_alloc
     end function gsl_interp_alloc
     subroutine gsl_interp_free(interp) bind(c)
       import
       type(c_ptr), value :: interp
     end subroutine gsl_interp_free
     function gsl_interp_init(interp, xa, ya, size) bind(c)
       import
       type(c_ptr), value :: interp
       real(c_double), dimension(*), intent(in) :: xa, ya
       integer(c_size_t), value :: size
       integer(c_int) gsl_interp_init
     end function gsl_interp_init
     function gsl_interp_accel_alloc() bind(c)
       import
       type(c_ptr) :: gsl_interp_accel_alloc
     end function gsl_interp_accel_alloc
     subroutine gsl_interp_accel_free(acc) bind(c)
       import
       type(c_ptr), value :: acc
     end subroutine gsl_interp_accel_free
     function gsl_interp_name(interp) bind(c)
       import
       type(c_ptr), value :: interp
       type(c_ptr) :: gsl_interp_name
     end function gsl_interp_name
     function fgsl_aux_interp_min_size(interp) bind(c)
       import
       type(c_ptr), value :: interp
       integer(c_long_long) :: fgsl_aux_interp_min_size
     end function fgsl_aux_interp_min_size
     function gsl_interp_bsearch(xa, x, index_lo, index_hi) bind(c)
       import
       real(c_double), dimension(*), intent(in) :: xa
       real(c_double), value :: x
       integer(c_size_t), value :: index_lo, index_hi
       integer(c_size_t) ::  gsl_interp_bsearch
     end function gsl_interp_bsearch
     function gsl_interp_accel_find(acc, xa, size, x) bind(c)
       import
       type(c_ptr), value :: acc
       real(c_double), dimension(*), intent(in) :: xa
       integer(c_size_t), value :: size
       real(c_double), value :: x
       integer(c_size_t) ::  gsl_interp_accel_find
     end function gsl_interp_accel_find
     function gsl_spline_alloc(interp_type, size) bind(c)
       import
       type(c_ptr), value :: interp_type
       integer(c_size_t), value :: size
       type(c_ptr) :: gsl_spline_alloc
     end function gsl_spline_alloc
     function gsl_spline_init(spline, xa, ya, size) bind(c)
       import
       type(c_ptr), value :: spline
       real(c_double), dimension(*), intent(in) :: xa, ya
       integer(c_size_t), value :: size
       integer(c_int) gsl_spline_init
     end function gsl_spline_init
     function gsl_spline_name(spline) bind(c)
       import
       type(c_ptr), value :: spline
       type(c_ptr) :: gsl_spline_name
     end function gsl_spline_name
     function fgsl_aux_spline_min_size(spline) bind(c)
       import
       type(c_ptr), value :: spline
       integer(c_long_long) :: fgsl_aux_spline_min_size
     end function fgsl_aux_spline_min_size
     function gsl_spline_eval(spline, x, acc) bind(c)
       import
       type(c_ptr), value :: spline
       real(c_double), value :: x
       type(c_ptr), value :: acc
       real(c_double) :: gsl_spline_eval
     end function gsl_spline_eval
     function gsl_spline_eval_e(spline, x, acc, y) bind(c)
       import
       type(c_ptr), value :: spline
       real(c_double), value :: x
       real(c_double), intent(out) :: y
       type(c_ptr), value :: acc
       integer(c_int) :: gsl_spline_eval_e
     end function gsl_spline_eval_e
     function gsl_spline_eval_deriv(spline, x, acc) bind(c)
       import
       type(c_ptr), value :: spline
       real(c_double), value :: x
       type(c_ptr), value :: acc
       real(c_double) :: gsl_spline_eval_deriv
     end function gsl_spline_eval_deriv
     function gsl_spline_eval_deriv_e(spline, x, acc, y) bind(c)
       import
       type(c_ptr), value :: spline
       real(c_double), value :: x
       real(c_double), intent(out) :: y
       type(c_ptr), value :: acc
       integer(c_int) :: gsl_spline_eval_deriv_e
     end function gsl_spline_eval_deriv_e
     function gsl_spline_eval_deriv2(spline, x, acc) bind(c)
       import
       type(c_ptr), value :: spline
       real(c_double), value :: x
       type(c_ptr), value :: acc
       real(c_double) :: gsl_spline_eval_deriv2
     end function gsl_spline_eval_deriv2
     function gsl_spline_eval_deriv2_e(spline, x, acc, y) bind(c)
       import
       type(c_ptr), value :: spline
       real(c_double), value :: x
       real(c_double), intent(out) :: y
       type(c_ptr), value :: acc
       integer(c_int) :: gsl_spline_eval_deriv2_e
     end function gsl_spline_eval_deriv2_e
     function gsl_spline_eval_integ(spline, a, b, acc) bind(c)
       import
       type(c_ptr), value :: spline
       real(c_double), value :: a, b
       type(c_ptr), value :: acc
       real(c_double) :: gsl_spline_eval_integ
     end function gsl_spline_eval_integ
     function gsl_spline_eval_integ_e(spline, a, b, acc, y) bind(c)
       import
       type(c_ptr), value :: spline
       real(c_double), value :: a, b
       real(c_double), intent(out) :: y
       type(c_ptr), value :: acc
       integer(c_int) :: gsl_spline_eval_integ_e
     end function gsl_spline_eval_integ_e
     subroutine gsl_spline_free(spline) bind(c)
       import
       type(c_ptr), value :: spline
     end subroutine gsl_spline_free
     function gsl_aux_sizeof_interp() bind(c)
       import :: c_size_t
       integer(c_size_t) :: gsl_aux_sizeof_interp
     end function gsl_aux_sizeof_interp
!-*-f90-*-
!
! Interfaces: Permutations and Combinations
!
  function gsl_permutation_alloc(n) bind(c)
    import
    integer(c_size_t), value :: n
    type(c_ptr) :: gsl_permutation_alloc
  end function gsl_permutation_alloc
  function gsl_permutation_calloc(n) bind(c)
    import
    integer(c_size_t), value :: n
    type(c_ptr) :: gsl_permutation_calloc
  end function gsl_permutation_calloc
  subroutine gsl_permutation_init(p) bind(c)
    import
    type(c_ptr), value :: p
  end subroutine gsl_permutation_init
  subroutine gsl_permutation_free(p) bind(c)
    import
    type(c_ptr), value :: p
  end subroutine gsl_permutation_free
  function gsl_permutation_memcpy(dest, src) bind(c)
    import
    type(c_ptr), value :: dest
    type(c_ptr), value :: src
    integer(c_int) :: gsl_permutation_memcpy
  end function gsl_permutation_memcpy
  function gsl_permutation_get(p, i) bind(c)
    import
    type(c_ptr), value :: p
    integer(c_size_t), value :: i
    integer(c_size_t) :: gsl_permutation_get
  end function gsl_permutation_get
  function gsl_permutation_swap(p, i, j) bind(c)
    import
    type(c_ptr), value :: p
    integer(c_size_t), value :: i, j
    integer(c_int) :: gsl_permutation_swap
  end function gsl_permutation_swap
  function gsl_permutation_size(p) bind(c)   
    import
    type(c_ptr), value :: p
    integer(c_size_t) :: gsl_permutation_size
  end function gsl_permutation_size
  function gsl_permutation_data(p) bind(c)   
    import
    type(c_ptr), value :: p
    type(c_ptr) :: gsl_permutation_data
  end function gsl_permutation_data
  function gsl_permutation_valid(p) bind(c)
    import
    type(c_ptr), value :: p
    integer(c_int) :: gsl_permutation_valid
  end function gsl_permutation_valid
  subroutine gsl_permutation_reverse(p) bind(c)
    import
    type(c_ptr), value :: p
  end subroutine gsl_permutation_reverse
  function gsl_permutation_inverse(inv, p) bind(c)
    import
    type(c_ptr), value :: inv
    type(c_ptr), value :: p
    integer(c_int) :: gsl_permutation_inverse
  end function gsl_permutation_inverse
  function gsl_permutation_next(p) bind(c)   
    import
    type(c_ptr), value :: p
    integer(c_int) :: gsl_permutation_next
  end function gsl_permutation_next
  function gsl_permutation_prev(p) bind(c)   
    import
    type(c_ptr), value :: p
    integer(c_int) :: gsl_permutation_prev
  end function gsl_permutation_prev
  function gsl_permute(p, data, stride, n) bind(c)
    import
    integer(c_size_t), dimension(*), intent(in) :: p
    integer(c_size_t), value :: stride, n
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_int) :: gsl_permute
  end function gsl_permute
  function gsl_permute_long(p, data, stride, n) bind(c)
    import
    integer(c_size_t), dimension(*), intent(in) :: p
    integer(c_size_t), value :: stride, n
    integer(c_long), dimension(*), intent(inout) :: data
    integer(c_int) :: gsl_permute_long
  end function gsl_permute_long
  function gsl_permute_inverse(p, data, stride, n) bind(c)
    import
    integer(c_size_t), dimension(*), intent(in) :: p
    integer(c_size_t), value :: stride, n
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_int) :: gsl_permute_inverse
  end function gsl_permute_inverse
  function gsl_permute_long_inverse(p, data, stride, n) bind(c)
    import
    integer(c_size_t), dimension(*), intent(in) :: p
    integer(c_size_t), value :: stride, n
    integer(c_long), dimension(*), intent(inout) :: data
    integer(c_int) :: gsl_permute_long_inverse
  end function gsl_permute_long_inverse
  function gsl_permute_vector(p,v) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: p,v
    integer(c_int) :: gsl_permute_vector
  end function gsl_permute_vector
  function gsl_permute_vector_inverse(p,v) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: p,v
    integer(c_int) :: gsl_permute_vector_inverse
  end function gsl_permute_vector_inverse
  function gsl_permutation_mul(p, pa, pb) bind(c)
    import
    type(c_ptr), value :: p
    type(c_ptr), value :: pa, pb
    integer(c_int) :: gsl_permutation_mul
  end function gsl_permutation_mul
  function gsl_permutation_fwrite(stream, p) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, p
    integer(c_int) :: gsl_permutation_fwrite
  end function gsl_permutation_fwrite
  function gsl_permutation_fread(stream, p) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, p
    integer(c_int) :: gsl_permutation_fread
  end function gsl_permutation_fread
  function gsl_permutation_fprintf(stream, p, format) bind(c)
    import :: c_ptr, c_int, c_char
    type(c_ptr), value :: stream, p
    character(kind=c_char), dimension(*) :: format
    integer(c_int) :: gsl_permutation_fprintf
  end function gsl_permutation_fprintf
  function gsl_permutation_fscanf(stream, p) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, p
    integer(c_int) :: gsl_permutation_fscanf
  end function gsl_permutation_fscanf
  function gsl_permutation_linear_to_canonical(q, p) bind(c)
    import
    type(c_ptr), value :: q
    type(c_ptr), value :: p
    integer(c_int) :: gsl_permutation_linear_to_canonical
  end function gsl_permutation_linear_to_canonical
  function gsl_permutation_canonical_to_linear(p, q) bind(c)
    import
    type(c_ptr), value :: p
    type(c_ptr), value :: q
    integer(c_int) :: gsl_permutation_canonical_to_linear
  end function gsl_permutation_canonical_to_linear
  function gsl_permutation_inversions(p) bind(c)   
    import
    type(c_ptr), value :: p
    integer(c_size_t) :: gsl_permutation_inversions
  end function gsl_permutation_inversions
  function gsl_permutation_linear_cycles(p) bind(c)   
    import
    type(c_ptr), value :: p
    integer(c_size_t) :: gsl_permutation_linear_cycles
  end function gsl_permutation_linear_cycles
  function gsl_permutation_canonical_cycles(p) bind(c)   
    import
    type(c_ptr), value :: p
    integer(c_size_t) :: gsl_permutation_canonical_cycles
  end function gsl_permutation_canonical_cycles
!
  function gsl_combination_alloc(n, k) bind(c)
    import
    integer(c_size_t), value :: n, k
    type(c_ptr) :: gsl_combination_alloc
  end function gsl_combination_alloc
  function gsl_combination_calloc(n, k) bind(c)
    import
    integer(c_size_t), value :: n, k
    type(c_ptr) :: gsl_combination_calloc
  end function gsl_combination_calloc
  subroutine gsl_combination_init_first(c) bind(c)
    import
    type(c_ptr), value :: c
  end subroutine gsl_combination_init_first
  subroutine gsl_combination_init_last(c) bind(c)
    import
    type(c_ptr), value :: c
  end subroutine gsl_combination_init_last
  subroutine gsl_combination_free(c) bind(c)
    import
    type(c_ptr), value :: c
  end subroutine gsl_combination_free
  function gsl_combination_memcpy(dest, src) bind(c)
    import
    type(c_ptr), value :: dest
    type(c_ptr), value :: src
    integer(c_int) :: gsl_combination_memcpy
  end function gsl_combination_memcpy
  function gsl_combination_get(c, i) bind(c)
    import
    type(c_ptr), value :: c
    integer(c_size_t), value :: i
    integer(c_size_t) :: gsl_combination_get
  end function gsl_combination_get
  function gsl_combination_n(c) bind(c)   
    import
    type(c_ptr), value :: c
    integer(c_size_t) :: gsl_combination_n
  end function gsl_combination_n
  function gsl_combination_k(c) bind(c)   
    import
    type(c_ptr), value :: c
    integer(c_size_t) :: gsl_combination_k
  end function gsl_combination_k
  function gsl_combination_data(c) bind(c)   
    import
    type(c_ptr), value :: c
    type(c_ptr) :: gsl_combination_data
  end function gsl_combination_data
  function gsl_combination_valid(c) bind(c)
    import
    type(c_ptr), value :: c
    integer(c_int) :: gsl_combination_valid
  end function gsl_combination_valid
  function gsl_combination_next(c) bind(c)   
    import
    type(c_ptr), value :: c
    integer(c_int) :: gsl_combination_next
  end function gsl_combination_next
  function gsl_combination_prev(c) bind(c)   
    import
    type(c_ptr), value :: c
    integer(c_int) :: gsl_combination_prev
  end function gsl_combination_prev
  function gsl_combination_fwrite(stream, c) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, c
    integer(c_int) :: gsl_combination_fwrite
  end function gsl_combination_fwrite
  function gsl_combination_fread(stream, c) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, c
    integer(c_int) :: gsl_combination_fread
  end function gsl_combination_fread
  function gsl_combination_fprintf(stream, c, format) bind(c)
    import :: c_ptr, c_int, c_char
    type(c_ptr), value :: stream, c
    character(kind=c_char), dimension(*) :: format
    integer(c_int) :: gsl_combination_fprintf
  end function gsl_combination_fprintf
  function gsl_combination_fscanf(stream, c) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, c
    integer(c_int) :: gsl_combination_fscanf
  end function gsl_combination_fscanf
!
  function gsl_aux_sizeof_permutation() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_permutation
  end function gsl_aux_sizeof_permutation
  function gsl_aux_sizeof_combination() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_combination
  end function gsl_aux_sizeof_combination
    
!-*-f90-*-
!
!  Interfaces: Sorting
!
  subroutine gsl_heapsort(array, count, size, compare) bind(c)
    import
    type(c_ptr), value :: array
    integer(c_size_t), value :: count, size
    type(c_funptr), value :: compare
  end subroutine gsl_heapsort
  function gsl_heapsort_index(p, array, count, size, compare) bind(c)
    import
    type(c_ptr), value :: array
    integer(c_size_t), value :: count, size
    integer(c_size_t), intent(out) :: p(count)
    type(c_funptr), value :: compare
    integer(c_int) :: gsl_heapsort_index
  end function gsl_heapsort_index
  subroutine gsl_sort(data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_size_t), value :: stride, n
  end subroutine gsl_sort
  subroutine gsl_sort_index(p, data, stride, n) bind(c)
    import
    integer(c_size_t), dimension(*), intent(out) :: p
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
  end subroutine gsl_sort_index
  function gsl_sort_smallest(dest, k, src, stride, n) bind(c)
    import
    integer(c_size_t), value :: k, stride, n
    real(c_double), intent(out) :: dest(k)
    real(c_double), dimension(*), intent(in) :: src
    integer(c_int) :: gsl_sort_smallest
  end function gsl_sort_smallest
  function gsl_sort_smallest_index(p, k, src, stride, n) bind(c)
    import
    integer(c_size_t), value :: k, stride, n
    integer(c_size_t), intent(out) :: p(k)
    real(c_double), dimension(*), intent(in) :: src
    integer(c_int) :: gsl_sort_smallest_index
  end function gsl_sort_smallest_index
  function gsl_sort_largest(dest, k, src, stride, n) bind(c)
    import
    integer(c_size_t), value :: k, stride, n
    real(c_double), intent(out) :: dest(k)
    real(c_double), dimension(*), intent(in) :: src
    integer(c_int) :: gsl_sort_largest
  end function gsl_sort_largest
  function gsl_sort_largest_index(p, k, src, stride, n) bind(c)
    import
    integer(c_size_t), value :: k, stride, n
    integer(c_size_t), intent(out) :: p(k)
    real(c_double), dimension(*), intent(in) :: src
    integer(c_int) :: gsl_sort_largest_index
  end function gsl_sort_largest_index
  subroutine gsl_sort_long(data, stride, n) bind(c)
    import
    integer(c_long), dimension(*), intent(inout) :: data
    integer(c_size_t), value :: stride, n
  end subroutine gsl_sort_long
  subroutine gsl_sort_long_index(p, data, stride, n) bind(c)
    import
    integer(c_size_t), dimension(*), intent(out) :: p
    integer(c_long), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
  end subroutine gsl_sort_long_index
  function gsl_sort_long_smallest(dest, k, src, stride, n) bind(c)
    import
    integer(c_size_t), value :: k, stride, n
    integer(c_long), intent(out) :: dest(k)
    integer(c_long), dimension(*), intent(in) :: src
    integer(c_int) :: gsl_sort_long_smallest
  end function gsl_sort_long_smallest
  function gsl_sort_long_smallest_index(p, k, src, stride, n) bind(c)
    import
    integer(c_size_t), value :: k, stride, n
    integer(c_size_t), intent(out) :: p(k)
    integer(c_long), dimension(*), intent(in) :: src
    integer(c_int) :: gsl_sort_long_smallest_index
  end function gsl_sort_long_smallest_index
  function gsl_sort_long_largest(dest, k, src, stride, n) bind(c)
    import
    integer(c_size_t), value :: k, stride, n
    integer(c_long), intent(out) :: dest(k)
    integer(c_long), dimension(*), intent(in) :: src
    integer(c_int) :: gsl_sort_long_largest
  end function gsl_sort_long_largest
  function gsl_sort_long_largest_index(p, k, src, stride, n) bind(c)
    import
    integer(c_size_t), value :: k, stride, n
    integer(c_size_t), intent(out) :: p(k)
    integer(c_long), dimension(*), intent(in) :: src
    integer(c_int) :: gsl_sort_long_largest_index
  end function gsl_sort_long_largest_index
    
    
!-*-f90-*-
!
!  Interfaces: Linear Algebra support
!
! LU
!
  function gsl_linalg_lu_decomp (a, p, signum) &
       bind(c, name='gsl_linalg_LU_decomp')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, p
    integer(c_int) :: signum
    integer(c_int) :: gsl_linalg_lu_decomp
  end function gsl_linalg_lu_decomp
  function gsl_linalg_complex_lu_decomp (a, p, signum) &
       bind(c, name='gsl_linalg_complex_LU_decomp')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, p
    integer(c_int) :: signum
    integer(c_int) :: gsl_linalg_complex_lu_decomp
  end function gsl_linalg_complex_lu_decomp
  function gsl_linalg_lu_solve (lu, p, b, x) &
       bind(c, name='gsl_linalg_LU_solve')
    import :: c_ptr, c_int
    type(c_ptr), value :: lu, p, b, x
    integer(c_int) :: gsl_linalg_lu_solve
  end function gsl_linalg_lu_solve
  function gsl_linalg_complex_lu_solve (lu, p, b, x) &
       bind(c, name='gsl_linalg_complex_LU_solve')
    import :: c_ptr, c_int
    type(c_ptr), value :: lu, p, b, x
    integer(c_int) :: gsl_linalg_complex_lu_solve
  end function gsl_linalg_complex_lu_solve
  function gsl_linalg_lu_svx (lu, p, x) &
       bind(c, name='gsl_linalg_LU_svx')
    import :: c_ptr, c_int
    type(c_ptr), value :: lu, p, x
    integer(c_int) :: gsl_linalg_lu_svx
  end function gsl_linalg_lu_svx
  function gsl_linalg_complex_lu_svx (lu, p, x) &
       bind(c, name='gsl_linalg_complex_LU_svx')
    import :: c_ptr, c_int
    type(c_ptr), value :: lu, p, x
    integer(c_int) :: gsl_linalg_complex_lu_svx
  end function gsl_linalg_complex_lu_svx
  function gsl_linalg_lu_refine (a, lu, p, b, x, residual) &
       bind(c, name='gsl_linalg_LU_refine')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, lu, p, b, x, residual
    integer(c_int) :: gsl_linalg_lu_refine
  end function gsl_linalg_lu_refine
  function gsl_linalg_complex_lu_refine (a, lu, p, b, x, residual) &
       bind(c, name='gsl_linalg_complex_LU_refine')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, lu, p, b, x, residual
    integer(c_int) :: gsl_linalg_complex_lu_refine
  end function gsl_linalg_complex_lu_refine
  function gsl_linalg_lu_invert (lu, p, inv) &
       bind(c, name='gsl_linalg_LU_invert')
    import :: c_ptr, c_int
    type(c_ptr), value :: lu, p, inv
    integer(c_int) :: gsl_linalg_lu_invert
  end function gsl_linalg_lu_invert
  function gsl_linalg_complex_lu_invert (lu, p, inv) &
       bind(c, name='gsl_linalg_complex_LU_invert')
    import :: c_ptr, c_int
    type(c_ptr), value :: lu, p, inv
    integer(c_int) :: gsl_linalg_complex_lu_invert
  end function gsl_linalg_complex_lu_invert
  function gsl_linalg_lu_det(lu, signum) bind(c, name='gsl_linalg_LU_det')
    import :: c_ptr, c_double, c_int
    type(c_ptr), value :: lu
    integer(c_int), value :: signum
    real(c_double) :: gsl_linalg_lu_det
  end function gsl_linalg_lu_det
  function gsl_linalg_complex_lu_det(lu, signum) &
       bind(c, name='gsl_linalg_complex_LU_det')
    import :: c_ptr, c_int, gsl_complex
    type(c_ptr), value :: lu
    integer(c_int), value :: signum
    type(gsl_complex) :: gsl_linalg_complex_lu_det
  end function gsl_linalg_complex_lu_det
  function gsl_linalg_lu_lndet(lu) bind(c, name='gsl_linalg_LU_lndet')
    import :: c_ptr, c_double
    type(c_ptr), value :: lu
    real(c_double) :: gsl_linalg_lu_lndet
  end function gsl_linalg_lu_lndet
  function gsl_linalg_complex_lu_lndet(lu) &
       bind(c, name='gsl_linalg_complex_LU_lndet')
    import :: c_ptr, c_double
    type(c_ptr), value :: lu
    real(c_double) :: gsl_linalg_complex_lu_lndet
  end function gsl_linalg_complex_lu_lndet
  function gsl_linalg_lu_sgndet(lu, signum) bind(c, name='gsl_linalg_LU_sgndet')
    import :: c_ptr, c_double, c_int
    type(c_ptr), value :: lu
    integer(c_int), value :: signum
    integer(c_int) :: gsl_linalg_lu_sgndet
  end function gsl_linalg_lu_sgndet
  function gsl_linalg_complex_lu_sgndet(lu, signum) &
       bind(c, name='gsl_linalg_complex_LU_sgndet')
    import :: c_ptr, c_int, gsl_complex
    type(c_ptr), value :: lu
    integer(c_int), value :: signum
    type(gsl_complex) :: gsl_linalg_complex_lu_sgndet
  end function gsl_linalg_complex_lu_sgndet
!
! QR
!
  function gsl_linalg_qr_decomp (a, tau) bind(c, name='gsl_linalg_QR_decomp')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, tau
    integer(c_int) :: gsl_linalg_qr_decomp
  end function gsl_linalg_qr_decomp
  function gsl_linalg_qr_solve (qr, tau, b, x) &
       bind(c, name='gsl_linalg_QR_solve')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, tau, b, x
    integer(c_int) :: gsl_linalg_qr_solve
  end function gsl_linalg_qr_solve
  function gsl_linalg_qr_svx (qr, tau, x) &
       bind(c, name='gsl_linalg_QR_svx')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, tau, x
    integer(c_int) :: gsl_linalg_qr_svx
  end function gsl_linalg_qr_svx
  function gsl_linalg_qr_lssolve (qr, tau, b, x, residual) &
       bind(c, name='gsl_linalg_QR_lssolve')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, tau, b, x, residual
    integer(c_int) :: gsl_linalg_qr_lssolve
  end function gsl_linalg_qr_lssolve
  function gsl_linalg_qr_qtvec (qr, tau, v) &
       bind(c, name='gsl_linalg_QR_QTvec')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, tau, v
    integer(c_int) :: gsl_linalg_qr_qtvec
  end function gsl_linalg_qr_qtvec
  function gsl_linalg_qr_qvec (qr, tau, v) &
       bind(c, name='gsl_linalg_QR_Qvec')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, tau, v
    integer(c_int) :: gsl_linalg_qr_qvec
  end function gsl_linalg_qr_qvec
  function gsl_linalg_qr_qtmat (qr, tau, a) &
       bind(c, name='gsl_linalg_QR_QTmat')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, tau, a
    integer(c_int) :: gsl_linalg_qr_qtmat
  end function gsl_linalg_qr_qtmat
  function gsl_linalg_qr_rsolve (qr, b, x) &
       bind(c, name='gsl_linalg_QR_Rsolve')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, b, x
    integer(c_int) :: gsl_linalg_qr_rsolve
  end function gsl_linalg_qr_rsolve
  function gsl_linalg_qr_rsvx (qr, x) &
       bind(c, name='gsl_linalg_QR_Rsvx')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, x
    integer(c_int) :: gsl_linalg_qr_rsvx
  end function gsl_linalg_qr_rsvx
  function gsl_linalg_qr_unpack (qr, tau, q, r) &
       bind(c, name='gsl_linalg_QR_unpack')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, tau, q, r
    integer(c_int) :: gsl_linalg_qr_unpack
  end function gsl_linalg_qr_unpack
  function gsl_linalg_qr_qrsolve (q, r, b, x) &
       bind(c, name='gsl_linalg_QR_QRsolve')
    import :: c_ptr, c_int
    type(c_ptr), value :: q, r, b, x
    integer(c_int) :: gsl_linalg_qr_qrsolve
  end function gsl_linalg_qr_qrsolve
  function gsl_linalg_qr_update (q, r, w, v) &
       bind(c, name='gsl_linalg_QR_update')
    import :: c_ptr, c_int
    type(c_ptr), value :: q, r, w, v
    integer(c_int) :: gsl_linalg_qr_update
  end function gsl_linalg_qr_update
  function gsl_linalg_r_solve (r, b, x) &
       bind(c, name='gsl_linalg_R_solve')
    import :: c_ptr, c_int
    type(c_ptr), value :: r, b, x
    integer(c_int) :: gsl_linalg_r_solve
  end function gsl_linalg_r_solve
  function gsl_linalg_r_svx (r, x) &
       bind(c, name='gsl_linalg_R_svx')
    import :: c_ptr, c_int
    type(c_ptr), value :: r, x
    integer(c_int) :: gsl_linalg_r_svx
  end function gsl_linalg_r_svx
  function gsl_linalg_qrpt_decomp (a, tau, p, signum, norm) &
       bind(c, name='gsl_linalg_QRPT_decomp')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, tau, p, norm
    integer(c_int), intent(out) :: signum
    integer(c_int) :: gsl_linalg_qrpt_decomp
  end function gsl_linalg_qrpt_decomp
  function gsl_linalg_qrpt_decomp2 (a, q, r, tau, p, signum, norm) &
       bind(c, name='gsl_linalg_QRPT_decomp2')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, q, r, tau, p, norm
    integer(c_int), intent(out) :: signum
    integer(c_int) :: gsl_linalg_qrpt_decomp2
  end function gsl_linalg_qrpt_decomp2
  function gsl_linalg_qrpt_solve (qr, tau, p, b, x) &
       bind(c, name='gsl_linalg_QRPT_solve')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, tau, p, b, x
    integer(c_int) :: gsl_linalg_qrpt_solve
  end function gsl_linalg_qrpt_solve
  function gsl_linalg_qrpt_svx (qr, tau, p, x) &
       bind(c, name='gsl_linalg_QRPT_svx')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, tau, p, x
    integer(c_int) :: gsl_linalg_qrpt_svx
  end function gsl_linalg_qrpt_svx
  function gsl_linalg_qrpt_qrsolve (q, r, p, b, x) &
       bind(c, name='gsl_linalg_QRPT_QRsolve')
    import :: c_ptr, c_int
    type(c_ptr), value :: q, r, p, b, x
    integer(c_int) :: gsl_linalg_qrpt_qrsolve
  end function gsl_linalg_qrpt_qrsolve
  function gsl_linalg_qrpt_update (q, r, p, w, v) &
       bind(c, name='gsl_linalg_QRPT_update')
    import :: c_ptr, c_int
    type(c_ptr), value :: q, r, p, w, v
    integer(c_int) :: gsl_linalg_qrpt_update
  end function gsl_linalg_qrpt_update
  function gsl_linalg_qrpt_rsolve (qr, p, b, x) &
       bind(c, name='gsl_linalg_QRPT_Rsolve')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, p, b, x
    integer(c_int) :: gsl_linalg_qrpt_rsolve
  end function gsl_linalg_qrpt_rsolve
  function gsl_linalg_qrpt_rsvx (qr, p, x) &
       bind(c, name='gsl_linalg_QRPT_Rsvx')
    import :: c_ptr, c_int
    type(c_ptr), value :: qr, p, x
    integer(c_int) :: gsl_linalg_qrpt_rsvx
  end function gsl_linalg_qrpt_rsvx
!
! SVD
!
  function gsl_linalg_sv_decomp (a, v, s, work) &
       bind(c, name='gsl_linalg_SV_decomp')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, v, s, work
    integer(c_int) :: signum
    integer(c_int) :: gsl_linalg_sv_decomp
  end function gsl_linalg_sv_decomp
  function gsl_linalg_sv_decomp_mod (a, x, v, s, work) &
       bind(c, name='gsl_linalg_SV_decomp_mod')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, x, v, s, work
    integer(c_int) :: signum
    integer(c_int) :: gsl_linalg_sv_decomp_mod
  end function gsl_linalg_sv_decomp_mod
  function gsl_linalg_sv_decomp_jacobi (a, v, s) &
       bind(c, name='gsl_linalg_SV_decomp_jacobi')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, v, s
    integer(c_int) :: signum
    integer(c_int) :: gsl_linalg_sv_decomp_jacobi
  end function gsl_linalg_sv_decomp_jacobi
  function gsl_linalg_sv_solve (u, v, s, b, x) &
       bind(c, name='gsl_linalg_SV_solve')
    import :: c_ptr, c_int
    type(c_ptr), value :: u, v, s, b, x
    integer(c_int) :: signum
    integer(c_int) :: gsl_linalg_sv_solve
  end function gsl_linalg_sv_solve
!
! Cholesky
!
  function gsl_linalg_cholesky_decomp (a) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a
    integer(c_int) :: gsl_linalg_cholesky_decomp
  end function gsl_linalg_cholesky_decomp
  function gsl_linalg_complex_cholesky_decomp (a) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a
    integer(c_int) :: gsl_linalg_complex_cholesky_decomp
  end function gsl_linalg_complex_cholesky_decomp
  function gsl_linalg_cholesky_solve (chol, b, x) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: chol, b, x
    integer(c_int) :: gsl_linalg_cholesky_solve
  end function gsl_linalg_cholesky_solve
  function gsl_linalg_complex_cholesky_solve (chol, b, x) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: chol, b, x
    integer(c_int) :: gsl_linalg_complex_cholesky_solve
  end function gsl_linalg_complex_cholesky_solve
  function gsl_linalg_cholesky_svx (chol, x) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: chol, x
    integer(c_int) :: gsl_linalg_cholesky_svx
  end function gsl_linalg_cholesky_svx
  function gsl_linalg_complex_cholesky_svx (chol, x) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: chol, x
    integer(c_int) :: gsl_linalg_complex_cholesky_svx
  end function gsl_linalg_complex_cholesky_svx
  function gsl_linalg_cholesky_invert (chol) bind(c)
    import :: c_ptr, c_int
    integer(c_int) :: gsl_linalg_cholesky_invert
    type(c_ptr), value :: chol
  end function gsl_linalg_cholesky_invert
!
! Tridiag
!
  function gsl_linalg_symmtd_decomp (a, tau) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a, tau
    integer(c_int) :: gsl_linalg_symmtd_decomp
  end function gsl_linalg_symmtd_decomp
  function gsl_linalg_symmtd_unpack (a, tau, q, diag, subdiag) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a, tau, q, diag, subdiag
    integer(c_int) :: gsl_linalg_symmtd_unpack
  end function gsl_linalg_symmtd_unpack
  function gsl_linalg_symmtd_unpack_t (a, diag, subdiag) &
       bind(c, name='gsl_linalg_symmtd_unpack_T')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, diag, subdiag
    integer(c_int) :: gsl_linalg_symmtd_unpack_t
  end function gsl_linalg_symmtd_unpack_t
  function gsl_linalg_hermtd_decomp (a, tau) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a, tau
    integer(c_int) :: gsl_linalg_hermtd_decomp
  end function gsl_linalg_hermtd_decomp
  function gsl_linalg_hermtd_unpack (a, tau, q, diag, subdiag) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a, tau, q, diag, subdiag
    integer(c_int) :: gsl_linalg_hermtd_unpack
  end function gsl_linalg_hermtd_unpack
  function gsl_linalg_hermtd_unpack_t (a, diag, subdiag) &
       bind(c, name='gsl_linalg_hermtd_unpack_T')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, diag, subdiag
    integer(c_int) :: gsl_linalg_hermtd_unpack_t
  end function gsl_linalg_hermtd_unpack_t
!
! Hessenberg
!
  function gsl_linalg_hessenberg_decomp (a, tau) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a, tau
    integer(c_int) :: gsl_linalg_hessenberg_decomp
  end function gsl_linalg_hessenberg_decomp
  function gsl_linalg_hessenberg_unpack (h, tau, u) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: h, tau, u
    integer(c_int) :: gsl_linalg_hessenberg_unpack
  end function gsl_linalg_hessenberg_unpack
  function gsl_linalg_hessenberg_unpack_accum (h, tau, v) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: h, tau, v
    integer(c_int) :: gsl_linalg_hessenberg_unpack_accum
  end function gsl_linalg_hessenberg_unpack_accum
  function gsl_linalg_hessenberg_set_zero (h) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: h
    integer(c_int) :: gsl_linalg_hessenberg_set_zero
  end function gsl_linalg_hessenberg_set_zero
  function gsl_linalg_hesstri_decomp (a, b, u, v, work) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a, b, u, v, work
    integer(c_int) :: gsl_linalg_hesstri_decomp
  end function gsl_linalg_hesstri_decomp
!
! Bidiag
!
  function gsl_linalg_bidiag_decomp (a, tau_u, tau_v) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a, tau_u, tau_v
    integer(c_int) :: gsl_linalg_bidiag_decomp
  end function gsl_linalg_bidiag_decomp
  function gsl_linalg_bidiag_unpack (a, tau_u, u, tau_v, v, diag, & 
       superdiag) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a, tau_u, tau_v, u, v, diag, superdiag
    integer(c_int) :: gsl_linalg_bidiag_unpack
  end function gsl_linalg_bidiag_unpack
  function gsl_linalg_bidiag_unpack2 (a, tau_u, tau_v, v) &
       bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a, tau_u, tau_v, v
    integer(c_int) :: gsl_linalg_bidiag_unpack2
  end function gsl_linalg_bidiag_unpack2
  function gsl_linalg_bidiag_unpack_b (a, diag, superdiag) &
       bind(c, name='gsl_linalg_bidiag_unpack_B')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, diag, superdiag
    integer(c_int) :: gsl_linalg_bidiag_unpack_b
  end function gsl_linalg_bidiag_unpack_b
!
! Householder
!
  function gsl_linalg_householder_transform (v) bind(c)
    import :: c_ptr, c_double
    type(c_ptr), value :: v
    real(c_double) :: gsl_linalg_householder_transform
  end function gsl_linalg_householder_transform
  function gsl_linalg_complex_householder_transform (v) bind(c)
    import :: c_ptr, gsl_complex
    type(c_ptr), value :: v
    type(gsl_complex) :: gsl_linalg_complex_householder_transform
  end function gsl_linalg_complex_householder_transform
  function gsl_linalg_householder_hm (tau, v, a) bind(c)
    import :: c_ptr, c_double, c_int
    real(c_double), value :: tau
    type(c_ptr), value :: v, a
    integer(c_int) :: gsl_linalg_householder_hm
  end function gsl_linalg_householder_hm
  function gsl_linalg_complex_householder_hm (tau, v, a) bind(c)
    import :: c_ptr, c_int, gsl_complex
    type(gsl_complex), value :: tau
    type(c_ptr), value :: v, a
    integer(c_int) :: gsl_linalg_complex_householder_hm
  end function gsl_linalg_complex_householder_hm
  function gsl_linalg_householder_mh (tau, v, a) bind(c)
    import :: c_ptr, c_double, c_int
    real(c_double), value :: tau
    type(c_ptr), value :: v, a
    integer(c_int) :: gsl_linalg_householder_mh
  end function gsl_linalg_householder_mh
  function gsl_linalg_complex_householder_mh (tau, v, a) bind(c)
    import :: c_ptr, c_int, gsl_complex
    type(gsl_complex), value :: tau
    type(c_ptr), value :: v, a
    integer(c_int) :: gsl_linalg_complex_householder_mh
  end function gsl_linalg_complex_householder_mh
  function gsl_linalg_householder_hv (tau, v, w) bind(c)
    import :: c_ptr, c_double, c_int
    real(c_double), value :: tau
    type(c_ptr), value :: v, w
    integer(c_int) :: gsl_linalg_householder_hv
  end function gsl_linalg_householder_hv
  function gsl_linalg_complex_householder_hv (tau, v, w) bind(c)
    import :: c_ptr, c_int, gsl_complex
    type(gsl_complex), value :: tau
    type(c_ptr), value :: v, w
    integer(c_int) :: gsl_linalg_complex_householder_hv
  end function gsl_linalg_complex_householder_hv
  function gsl_linalg_hh_solve (a, b, x) &
       bind(c, name='gsl_linalg_HH_solve')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, b, x
    integer(c_int) :: gsl_linalg_hh_solve
  end function gsl_linalg_hh_solve
  function gsl_linalg_hh_svx (a, x) &
       bind(c, name='gsl_linalg_HH_svx')
    import :: c_ptr, c_int
    type(c_ptr), value :: a, x
    integer(c_int) :: gsl_linalg_hh_svx
  end function gsl_linalg_hh_svx
!
! Tridiagonal
!
  function gsl_linalg_solve_tridiag(diag, e, f, b, x) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: diag, e, f, b, x
    integer(c_int) :: gsl_linalg_solve_tridiag
  end function gsl_linalg_solve_tridiag
  function gsl_linalg_solve_symm_tridiag(diag, e, b, x) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: diag, e, b, x
    integer(c_int) :: gsl_linalg_solve_symm_tridiag
  end function gsl_linalg_solve_symm_tridiag
  function gsl_linalg_solve_cyc_tridiag(diag, e, f, b, x) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: diag, e, f, b, x
    integer(c_int) :: gsl_linalg_solve_cyc_tridiag
  end function gsl_linalg_solve_cyc_tridiag
  function gsl_linalg_solve_symm_cyc_tridiag(diag, e, b, x) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: diag, e, b, x
    integer(c_int) :: gsl_linalg_solve_symm_cyc_tridiag
  end function gsl_linalg_solve_symm_cyc_tridiag
  function gsl_linalg_balance_matrix (a, d) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: a, d
    integer(c_int) :: gsl_linalg_balance_matrix
  end function gsl_linalg_balance_matrix


!-*-f90-*-
!
!  Interfaces: Eigensystem support
!
function gsl_eigen_symm_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_symm_alloc
end function gsl_eigen_symm_alloc
subroutine gsl_eigen_symm_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_symm_free
function gsl_eigen_symm(a, eval, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, eval, w
  integer(c_int) :: gsl_eigen_symm
end function gsl_eigen_symm
function gsl_eigen_symmv_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_symmv_alloc
end function gsl_eigen_symmv_alloc
subroutine gsl_eigen_symmv_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_symmv_free
function gsl_eigen_symmv(a, eval, evec, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, eval, evec, w
  integer(c_int) :: gsl_eigen_symmv
end function gsl_eigen_symmv
function gsl_eigen_herm_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_herm_alloc
end function gsl_eigen_herm_alloc
subroutine gsl_eigen_herm_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_herm_free
function gsl_eigen_herm(a, eval, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, eval, w
  integer(c_int) :: gsl_eigen_herm
end function gsl_eigen_herm
function gsl_eigen_hermv_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_hermv_alloc
end function gsl_eigen_hermv_alloc
subroutine gsl_eigen_hermv_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_hermv_free
function gsl_eigen_hermv(a, eval, evec, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, eval, evec, w
  integer(c_int) :: gsl_eigen_hermv
end function gsl_eigen_hermv
function gsl_eigen_nonsymm_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_nonsymm_alloc
end function gsl_eigen_nonsymm_alloc
subroutine gsl_eigen_nonsymm_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_nonsymm_free
subroutine  gsl_eigen_nonsymm_params (compute_t, balance, w) bind(c)
  import :: c_int, c_ptr
  integer(c_int), value :: compute_t, balance
  type(c_ptr), value :: w
end subroutine gsl_eigen_nonsymm_params
function gsl_eigen_nonsymm(a, eval, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, eval, w
  integer(c_int) :: gsl_eigen_nonsymm
end function gsl_eigen_nonsymm
function gsl_eigen_nonsymm_z(a, eval, z, w) bind(c, name='gsl_eigen_nonsymm_Z')
  import :: c_int, c_ptr
  type(c_ptr), value :: a, eval, z, w
  integer(c_int) :: gsl_eigen_nonsymm_z
end function gsl_eigen_nonsymm_z
function gsl_eigen_nonsymmv_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_nonsymmv_alloc
end function gsl_eigen_nonsymmv_alloc
subroutine gsl_eigen_nonsymmv_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_nonsymmv_free
function gsl_eigen_nonsymmv(a, eval, evec, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, eval, evec, w
  integer(c_int) :: gsl_eigen_nonsymmv
end function gsl_eigen_nonsymmv
function gsl_eigen_nonsymmv_z(a, eval, evec, z, w) &
     bind(c, name='gsl_eigen_nonsymmv_Z')
  import :: c_int, c_ptr
  type(c_ptr), value :: a, eval, evec, z, w
  integer(c_int) :: gsl_eigen_nonsymmv_z
end function gsl_eigen_nonsymmv_z
function gsl_eigen_gensymm_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_gensymm_alloc
end function gsl_eigen_gensymm_alloc
subroutine gsl_eigen_gensymm_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_gensymm_free
function gsl_eigen_gensymm(a, b, eval, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, b, eval, w
  integer(c_int) :: gsl_eigen_gensymm
end function gsl_eigen_gensymm
function gsl_eigen_gensymmv_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_gensymmv_alloc
end function gsl_eigen_gensymmv_alloc
subroutine gsl_eigen_gensymmv_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_gensymmv_free
function gsl_eigen_gensymmv(a, b, eval, evec, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, b, eval, evec, w
  integer(c_int) :: gsl_eigen_gensymmv
end function gsl_eigen_gensymmv
function gsl_eigen_genherm_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_genherm_alloc
end function gsl_eigen_genherm_alloc
subroutine gsl_eigen_genherm_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_genherm_free
function gsl_eigen_genherm(a, b, eval, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, b, eval, w
  integer(c_int) :: gsl_eigen_genherm
end function gsl_eigen_genherm
function gsl_eigen_genhermv_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_genhermv_alloc
end function gsl_eigen_genhermv_alloc
subroutine gsl_eigen_genhermv_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_genhermv_free
function gsl_eigen_genhermv(a, b, eval, evec, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, b, eval, evec, w
  integer(c_int) :: gsl_eigen_genhermv
end function gsl_eigen_genhermv
function gsl_eigen_gen_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_gen_alloc
end function gsl_eigen_gen_alloc
subroutine gsl_eigen_gen_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_gen_free
subroutine  gsl_eigen_gen_params (compute_s, compute_t, balance, w) bind(c)
  import :: c_int, c_ptr
  integer(c_int), value :: compute_s, compute_t, balance
  type(c_ptr), value :: w
end subroutine gsl_eigen_gen_params
function gsl_eigen_gen(a, b, alpha, beta, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, b, alpha, beta, w
  integer(c_int) :: gsl_eigen_gen
end function gsl_eigen_gen
function gsl_eigen_gen_qz(a, b, alpha, beta, q, z, w) &
     bind(c, name='gsl_eigen_gen_QZ')
  import :: c_int, c_ptr
  type(c_ptr), value :: a, b, alpha, beta, q, z, w
  integer(c_int) :: gsl_eigen_gen_qz
end function gsl_eigen_gen_qz
function gsl_eigen_genv_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_eigen_genv_alloc
end function gsl_eigen_genv_alloc
subroutine gsl_eigen_genv_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_eigen_genv_free
function gsl_eigen_genv(a, b, alpha, beta, evec, w) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: a, b, alpha, beta, evec, w
  integer(c_int) :: gsl_eigen_genv
end function gsl_eigen_genv
function gsl_eigen_genv_qz(a, b, alpha, beta, evec, q, z, w) &
     bind(c, name='gsl_eigen_genv_QZ')
  import :: c_int, c_ptr
  type(c_ptr), value :: a, b, alpha, beta, evec, q, z, w
  integer(c_int) :: gsl_eigen_genv_qz
end function gsl_eigen_genv_qz
function gsl_eigen_symmv_sort (eval, evec, sort_type) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: eval, evec
  integer(c_int), value :: sort_type
  integer(c_int) :: gsl_eigen_symmv_sort
end function gsl_eigen_symmv_sort
function gsl_eigen_hermv_sort (eval, evec, sort_type) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: eval, evec
  integer(c_int), value :: sort_type
  integer(c_int) :: gsl_eigen_hermv_sort
end function gsl_eigen_hermv_sort
function gsl_eigen_nonsymmv_sort (eval, evec, sort_type) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: eval, evec
  integer(c_int), value :: sort_type
  integer(c_int) :: gsl_eigen_nonsymmv_sort
end function gsl_eigen_nonsymmv_sort
function gsl_eigen_gensymmv_sort (eval, evec, sort_type) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: eval, evec
  integer(c_int), value :: sort_type
  integer(c_int) :: gsl_eigen_gensymmv_sort
end function gsl_eigen_gensymmv_sort
function gsl_eigen_genhermv_sort (eval, evec, sort_type) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: eval, evec
  integer(c_int), value :: sort_type
  integer(c_int) :: gsl_eigen_genhermv_sort
end function gsl_eigen_genhermv_sort
function gsl_eigen_genv_sort (alpha, beta, evec, sort_type) bind(c)
  import :: c_int, c_ptr
  type(c_ptr), value :: alpha, beta, evec
  integer(c_int), value :: sort_type
  integer(c_int) :: gsl_eigen_genv_sort
end function gsl_eigen_genv_sort


!-*-f90-*-
!
!  Interface: Fast FT support
!
function gsl_fft_complex_radix2_forward(data, stride, n) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data
  integer(c_int) :: gsl_fft_complex_radix2_forward
end function gsl_fft_complex_radix2_forward
function gsl_fft_complex_radix2_transform(data, stride, n, sign) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data
  integer(c_int), value :: sign
  integer(c_int) :: gsl_fft_complex_radix2_transform
end function gsl_fft_complex_radix2_transform
function gsl_fft_complex_radix2_backward(data, stride, n) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data
  integer(c_int) :: gsl_fft_complex_radix2_backward
end function gsl_fft_complex_radix2_backward
function gsl_fft_complex_radix2_inverse(data, stride, n) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data
  integer(c_int) :: gsl_fft_complex_radix2_inverse
end function gsl_fft_complex_radix2_inverse
function gsl_fft_complex_radix2_dif_forward(data, stride, n) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data
  integer(c_int) :: gsl_fft_complex_radix2_dif_forward
end function gsl_fft_complex_radix2_dif_forward
function gsl_fft_complex_radix2_dif_transform(data, stride, n, sign) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data
  integer(c_int), value :: sign
  integer(c_int) :: gsl_fft_complex_radix2_dif_transform
end function gsl_fft_complex_radix2_dif_transform
function gsl_fft_complex_radix2_dif_backward(data, stride, n) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data
  integer(c_int) :: gsl_fft_complex_radix2_dif_backward
end function gsl_fft_complex_radix2_dif_backward
function gsl_fft_complex_radix2_dif_inverse(data, stride, n) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data
  integer(c_int) :: gsl_fft_complex_radix2_dif_inverse
end function gsl_fft_complex_radix2_dif_inverse
!
function gsl_fft_complex_wavetable_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_fft_complex_wavetable_alloc
end function gsl_fft_complex_wavetable_alloc
subroutine gsl_fft_complex_wavetable_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_fft_complex_wavetable_free
function gsl_fft_complex_workspace_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_fft_complex_workspace_alloc
end function gsl_fft_complex_workspace_alloc
subroutine gsl_fft_complex_workspace_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_fft_complex_workspace_free
function gsl_fft_complex_forward(data, stride, n, wavetable, work) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data, wavetable, work
  integer(c_int) :: gsl_fft_complex_forward
end function gsl_fft_complex_forward
function gsl_fft_complex_transform(data, stride, n, wavetable, work, sign) &
     bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data, wavetable, work
  integer(c_int), value :: sign
  integer(c_int) :: gsl_fft_complex_transform
end function gsl_fft_complex_transform
function gsl_fft_complex_backward(data, stride, n, wavetable, work) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data, wavetable, work
  integer(c_int) :: gsl_fft_complex_backward
end function gsl_fft_complex_backward
function gsl_fft_complex_inverse(data, stride, n, wavetable, work) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data, wavetable, work
  integer(c_int) :: gsl_fft_complex_inverse
end function gsl_fft_complex_inverse
function gsl_fft_real_radix2_transform(data, stride, n) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data
  integer(c_int) :: gsl_fft_real_radix2_transform
end function gsl_fft_real_radix2_transform
function gsl_fft_halfcomplex_radix2_inverse(data, stride, n) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data
  integer(c_int) :: gsl_fft_halfcomplex_radix2_inverse
end function gsl_fft_halfcomplex_radix2_inverse
function gsl_fft_halfcomplex_radix2_backward(data, stride, n) bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data
  integer(c_int) :: gsl_fft_halfcomplex_radix2_backward
end function gsl_fft_halfcomplex_radix2_backward
function gsl_fft_real_wavetable_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_fft_real_wavetable_alloc
end function gsl_fft_real_wavetable_alloc
subroutine gsl_fft_real_wavetable_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_fft_real_wavetable_free
function gsl_fft_halfcomplex_wavetable_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_fft_halfcomplex_wavetable_alloc
end function gsl_fft_halfcomplex_wavetable_alloc
subroutine gsl_fft_halfcomplex_wavetable_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_fft_halfcomplex_wavetable_free
function gsl_fft_real_workspace_alloc(n) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: n
  type(c_ptr) :: gsl_fft_real_workspace_alloc
end function gsl_fft_real_workspace_alloc
subroutine gsl_fft_real_workspace_free(w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_fft_real_workspace_free
function gsl_fft_real_transform(data, stride, n, wavetable, work) &
     bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data, wavetable, work
  integer(c_int) :: gsl_fft_real_transform
end function gsl_fft_real_transform
function gsl_fft_halfcomplex_transform(data, stride, n, wavetable, work) &
     bind(c)
  import :: c_size_t, c_ptr, c_int
  integer(c_size_t), value :: stride, n
  type(c_ptr), value :: data, wavetable, work
  integer(c_int) :: gsl_fft_halfcomplex_transform
end function gsl_fft_halfcomplex_transform
function gsl_fft_real_unpack (real_coefficient, complex_coefficient, &
     stride, n) bind(c)
  import :: c_size_t, c_ptr, c_int
  type(c_ptr), value :: real_coefficient, complex_coefficient
  integer(c_size_t), value :: stride, n
  integer(c_int) :: gsl_fft_real_unpack
end function gsl_fft_real_unpack
function gsl_fft_halfcomplex_unpack (halfcomplex_coefficient, complex_coefficient, &
     stride, n) bind(c)
  import :: c_size_t, c_ptr, c_int
  type(c_ptr), value :: halfcomplex_coefficient, complex_coefficient
  integer(c_size_t), value :: stride, n
  integer(c_int) :: gsl_fft_halfcomplex_unpack
end function gsl_fft_halfcomplex_unpack
!-*-f90-*-
!
! Interfaces: Numerical Integration
!
    function gsl_integration_qng(f, a, b, epsabs, epsrel, result, abserr, neval) bind(c)
      import
      type(c_ptr), value :: f
      real(c_double), value :: a, b, epsabs, epsrel
      real(c_double), intent(out) :: result, abserr
      integer(c_size_t), intent(inout) :: neval
      integer(c_int) :: gsl_integration_qng
    end function gsl_integration_qng
    function gsl_integration_workspace_alloc (n) bind(c)
      import
      integer(c_size_t), value :: n
      type(c_ptr) :: gsl_integration_workspace_alloc
    end function gsl_integration_workspace_alloc
    subroutine gsl_integration_workspace_free (w) bind(c)
      import
      type(c_ptr), value :: w
    end subroutine gsl_integration_workspace_free
    function gsl_integration_qag(f, a, b, epsabs, epsrel, limit, key, &
         workspace, result, abserr) bind(c)
      import
      type(c_ptr), value :: f
      real(c_double), value :: a, b, epsabs, epsrel
      integer(c_size_t), value :: limit
      integer(c_int), value :: key
      type(c_ptr), value :: workspace
      real(c_double), intent(out) :: result, abserr
      integer(c_int) :: gsl_integration_qag
    end function gsl_integration_qag
    function gsl_integration_qags(f, a, b, epsabs, epsrel, limit, &
         workspace, result, abserr) bind(c)
      import
      type(c_ptr), value :: f
      real(c_double), value :: a, b, epsabs, epsrel
      integer(c_size_t), value :: limit
      type(c_ptr), value :: workspace
      real(c_double), intent(out) :: result, abserr
      integer(c_int) :: gsl_integration_qags
    end function gsl_integration_qags
    function gsl_integration_qagp(f, pts, npts, epsabs, epsrel, limit, &
         workspace, result, abserr) bind(c)
      import
      type(c_ptr), value :: f
      real(c_double), dimension(*), intent(in) :: pts
      integer(c_size_t), value :: npts
      real(c_double), value :: epsabs, epsrel
      integer(c_size_t), value :: limit
      type(c_ptr), value :: workspace
      real(c_double), intent(out) :: result, abserr
      integer(c_int) :: gsl_integration_qagp
    end function gsl_integration_qagp
    function gsl_integration_qagi(f, epsabs, epsrel, limit, &
         workspace, result, abserr) bind(c)
      import
      type(c_ptr), value :: f
      real(c_double), value :: epsabs, epsrel
      integer(c_size_t), value :: limit
      type(c_ptr), value :: workspace
      real(c_double), intent(out) :: result, abserr
      integer(c_int) :: gsl_integration_qagi
    end function gsl_integration_qagi
    function gsl_integration_qagiu(f, a, epsabs, epsrel, limit, &
         workspace, result, abserr) bind(c)
      import
      type(c_ptr), value :: f
      real(c_double), value :: a, epsabs, epsrel
      integer(c_size_t), value :: limit
      type(c_ptr), value :: workspace
      real(c_double), intent(out) :: result, abserr
      integer(c_int) :: gsl_integration_qagiu
    end function gsl_integration_qagiu
    function gsl_integration_qagil(f, b, epsabs, epsrel, limit, &
         workspace, result, abserr) bind(c)
      import
      type(c_ptr), value :: f
      real(c_double), value :: b, epsabs, epsrel
      integer(c_size_t), value :: limit
      type(c_ptr), value :: workspace
      real(c_double), intent(out) :: result, abserr
      integer(c_int) :: gsl_integration_qagil
    end function gsl_integration_qagil
    function gsl_integration_qawc(f, a, b, c, epsabs, epsrel, limit, &
         workspace, result, abserr) bind(c)
      import
      type(c_ptr), value :: f
      real(c_double), value :: a, b, c, epsabs, epsrel
      integer(c_size_t), value :: limit
      type(c_ptr), value :: workspace
      real(c_double), intent(out) :: result, abserr
      integer(c_int) :: gsl_integration_qawc
    end function gsl_integration_qawc
    function gsl_integration_qaws_table_alloc (alpha, beta, mu, nu) bind(c)
      import
      real(c_double), value :: alpha, beta
      integer(c_int), value :: mu, nu
      type(c_ptr) :: gsl_integration_qaws_table_alloc
    end function gsl_integration_qaws_table_alloc
    function gsl_integration_qaws_table_set(t, alpha, beta, mu, nu) bind(c)
      import
      type(c_ptr), value :: t
      real(c_double), value :: alpha, beta
      integer(c_int), value :: mu, nu
      integer(c_int) :: gsl_integration_qaws_table_set
    end function gsl_integration_qaws_table_set
    subroutine gsl_integration_qaws_table_free (w) bind(c)
      import
      type(c_ptr), value :: w
    end subroutine gsl_integration_qaws_table_free
    function gsl_integration_qaws(f, a, b, t, epsabs, epsrel, limit, workspace, &
         result, abserr) bind(c)
      import
      type(c_ptr), value :: f
      real(c_double), value :: a, b, epsabs, epsrel
      type(c_ptr), value :: t
      integer(c_size_t), value :: limit
      type(c_ptr), value :: workspace
      real(c_double), intent(out) :: result, abserr
      integer(c_int) :: gsl_integration_qaws
    end function gsl_integration_qaws
    function gsl_integration_qawo_table_alloc(omega, l, sine, n) bind(c)
      import
      real(c_double), value :: omega, l
      integer(c_int), value :: sine
      integer(c_size_t), value :: n
      type(c_ptr) :: gsl_integration_qawo_table_alloc
    end function gsl_integration_qawo_table_alloc
    function gsl_integration_qawo_table_set(t, omega, l, sine) bind(c)
      import
      type(c_ptr), value :: t
      real(c_double), value :: omega, l
      integer(c_int), value :: sine
      integer(c_int) :: gsl_integration_qawo_table_set
    end function gsl_integration_qawo_table_set
    function gsl_integration_qawo_table_set_length(t, l) bind(c)
      import
      type(c_ptr), value :: t
      real(c_double), value :: l
      integer(c_int) :: gsl_integration_qawo_table_set_length
    end function gsl_integration_qawo_table_set_length
    subroutine gsl_integration_qawo_table_free (w) bind(c)
      import
      type(c_ptr), value :: w
    end subroutine gsl_integration_qawo_table_free
    function gsl_integration_qawo(f, a, epsabs, epsrel, limit, workspace, &
         wf, result, abserr) bind(c)
      import
      type(c_ptr), value :: f
      real(c_double), value :: a,  epsabs, epsrel
      integer(c_size_t), value :: limit
      type(c_ptr), value :: workspace, wf
      real(c_double), intent(out) :: result, abserr
      integer(c_int) :: gsl_integration_qawo
    end function gsl_integration_qawo
    function gsl_integration_qawf(f, a, epsabs, limit, workspace, cyc_workspace, &
         wf, result, abserr) bind(c)
      import
      type(c_ptr), value :: f
      real(c_double), value :: a,  epsabs
      integer(c_size_t), value :: limit
      type(c_ptr), value :: workspace, cyc_workspace, wf
      real(c_double), intent(out) :: result, abserr
      integer(c_int) :: gsl_integration_qawf
    end function gsl_integration_qawf
    function gsl_aux_sizeof_integration_workspace() bind(c)
      import :: c_size_t
      integer(c_size_t) :: gsl_aux_sizeof_integration_workspace
    end function gsl_aux_sizeof_integration_workspace
    function gsl_aux_sizeof_integration_qaws_table() bind(c)
      import :: c_size_t
      integer(c_size_t) :: gsl_aux_sizeof_integration_qaws_table
    end function gsl_aux_sizeof_integration_qaws_table
    function gsl_aux_sizeof_integration_qawo_table() bind(c)
      import :: c_size_t
      integer(c_size_t) :: gsl_aux_sizeof_integration_qawo_table
    end function gsl_aux_sizeof_integration_qawo_table
!-*-f90-*-
!
! Interfaces: Random and Quasi-random numbers, distributions
!
  function gsl_rng_alloc(t) bind(c)
    import
    type(c_ptr), value :: t
    type(c_ptr) :: gsl_rng_alloc
  end function gsl_rng_alloc
  subroutine gsl_rng_set(r, s) bind(c)
    import
    type(c_ptr), value :: r
    integer(c_long), value :: s
  end subroutine gsl_rng_set
  subroutine gsl_rng_free(r) bind(c)
    import
    type(c_ptr), value :: r
  end subroutine gsl_rng_free
  function gsl_rng_get(r) bind(c)
    import
    type(c_ptr), value :: r
    integer(c_long) :: gsl_rng_get
  end function gsl_rng_get
  function gsl_rng_uniform(r) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double) :: gsl_rng_uniform
  end function gsl_rng_uniform
  function gsl_rng_uniform_pos(r) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double) :: gsl_rng_uniform_pos
  end function gsl_rng_uniform_pos
  function gsl_rng_uniform_int(r, n) bind(c)
    import
    type(c_ptr), value :: r
    integer(c_long), value :: n
    integer(c_long) :: gsl_rng_uniform_int
  end function gsl_rng_uniform_int
  function gsl_rng_name(r) bind(c)
    import
    type(c_ptr), value :: r
    type(c_ptr) :: gsl_rng_name
  end function gsl_rng_name
  function gsl_rng_max(r) bind(c)
    import
    type(c_ptr), value :: r
    integer(c_long) :: gsl_rng_max
  end function gsl_rng_max
  function gsl_rng_min(r) bind(c)
    import
    type(c_ptr), value :: r
    integer(c_long) :: gsl_rng_min
  end function gsl_rng_min
  function gsl_rng_env_setup() bind(c)
    import
    type(c_ptr) :: gsl_rng_env_setup
  end function gsl_rng_env_setup
  function gsl_rng_memcpy(cpy, src) bind(c)
    import
    type(c_ptr), value :: cpy
    type(c_ptr), value :: src
    integer(c_int) :: gsl_rng_memcpy
  end function gsl_rng_memcpy
  function gsl_rng_clone(r) bind(c)
    import
    type(c_ptr), value :: r
    type(c_ptr) :: gsl_rng_clone
  end function gsl_rng_clone
  function gsl_rng_fwrite(stream, r) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, r
    integer(c_int) :: gsl_rng_fwrite
  end function gsl_rng_fwrite
  function gsl_rng_fread(stream, r) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, r
    integer(c_int) :: gsl_rng_fread
  end function gsl_rng_fread
  function fgsl_aux_rng_assign(i) bind(c)
    import
    integer(c_int), value :: i
    type(c_ptr) :: fgsl_aux_rng_assign
  end function fgsl_aux_rng_assign
!
  function gsl_qrng_alloc(t, d) bind(c)
    import
    type(c_ptr), value :: t
    integer(c_int), value :: d
    type(c_ptr) :: gsl_qrng_alloc
  end function gsl_qrng_alloc
!
  subroutine gsl_qrng_free(q) bind(c)
    import
    type(c_ptr), value :: q
  end subroutine gsl_qrng_free
  subroutine gsl_qrng_init(q) bind(c)
    import
    type(c_ptr), value :: q
  end subroutine gsl_qrng_init
  function gsl_qrng_get(q, x) bind(c)
    import
    type(c_ptr), value :: q
    real(c_double), dimension(*), intent(out) :: x
    integer(c_int) :: gsl_qrng_get
  end function gsl_qrng_get
  function gsl_qrng_name(q) bind(c)
    import
    type(c_ptr), value :: q
    type(c_ptr) :: gsl_qrng_name
  end function gsl_qrng_name
  function gsl_qrng_memcpy(cpy, src) bind(c)
    import
    type(c_ptr), value :: cpy
    type(c_ptr), value :: src
    integer(c_int) :: gsl_qrng_memcpy
  end function gsl_qrng_memcpy
  function gsl_qrng_clone(q) bind(c)
    import
    type(c_ptr), value :: q
    type(c_ptr) :: gsl_qrng_clone
  end function gsl_qrng_clone

!  
  function fgsl_aux_qrng_assign(i) bind(c)
    import
    integer(c_int), value :: i
    type(c_ptr) :: fgsl_aux_qrng_assign
  end function fgsl_aux_qrng_assign
!
  function gsl_ran_gaussian(r, sigma) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: sigma
    real(c_double) :: gsl_ran_gaussian
  end function gsl_ran_gaussian
  function gsl_ran_gaussian_pdf(x, sigma) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: sigma
    real(c_double) :: gsl_ran_gaussian_pdf
  end function gsl_ran_gaussian_pdf
  function gsl_ran_gaussian_ziggurat(r, sigma) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: sigma
    real(c_double) :: gsl_ran_gaussian_ziggurat
  end function gsl_ran_gaussian_ziggurat
  function gsl_ran_gaussian_ratio_method(r, sigma) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: sigma
    real(c_double) :: gsl_ran_gaussian_ratio_method
  end function gsl_ran_gaussian_ratio_method
  function gsl_ran_ugaussian(r) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double) :: gsl_ran_ugaussian
  end function gsl_ran_ugaussian
  function gsl_ran_ugaussian_pdf(x) bind(c)
    import
    real(c_double), value :: x
    real(c_double) :: gsl_ran_ugaussian_pdf
  end function gsl_ran_ugaussian_pdf
  function gsl_ran_ugaussian_ratio_method(r) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double) :: gsl_ran_ugaussian_ratio_method
  end function gsl_ran_ugaussian_ratio_method
  function gsl_cdf_gaussian_p(x, sigma) bind(c, name='gsl_cdf_gaussian_P')
    import
    real(c_double), value :: x
    real(c_double), value :: sigma
    real(c_double) :: gsl_cdf_gaussian_p
  end function gsl_cdf_gaussian_p
  function gsl_cdf_gaussian_q(x, sigma) bind(c, name='gsl_cdf_gaussian_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: sigma
    real(c_double) :: gsl_cdf_gaussian_q
  end function gsl_cdf_gaussian_q
  function gsl_cdf_gaussian_pinv(p, sigma) bind(c, name='gsl_cdf_gaussian_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: sigma
    real(c_double) :: gsl_cdf_gaussian_pinv
  end function gsl_cdf_gaussian_pinv
  function gsl_cdf_gaussian_qinv(q, sigma) bind(c, name='gsl_cdf_gaussian_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: sigma
    real(c_double) :: gsl_cdf_gaussian_qinv
  end function gsl_cdf_gaussian_qinv
  function gsl_cdf_ugaussian_p(x) bind(c, name='gsl_cdf_ugaussian_P')
    import
    real(c_double), value :: x
    real(c_double) :: gsl_cdf_ugaussian_p
  end function gsl_cdf_ugaussian_p
  function gsl_cdf_ugaussian_q(x) bind(c, name='gsl_cdf_ugaussian_Q')
    import
    real(c_double), value :: x
    real(c_double) :: gsl_cdf_ugaussian_q
  end function gsl_cdf_ugaussian_q
  function gsl_cdf_ugaussian_pinv(p) bind(c, name='gsl_cdf_ugaussian_Pinv')
    import
    real(c_double), value :: p
    real(c_double) :: gsl_cdf_ugaussian_pinv
  end function gsl_cdf_ugaussian_pinv
  function gsl_cdf_ugaussian_qinv(q) bind(c, name='gsl_cdf_ugaussian_Qinv')
    import
    real(c_double), value :: q
    real(c_double) :: gsl_cdf_ugaussian_qinv
  end function gsl_cdf_ugaussian_qinv
  function gsl_ran_gaussian_tail(r, a, sigma) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a, sigma
    real(c_double) :: gsl_ran_gaussian_tail
  end function gsl_ran_gaussian_tail
  function gsl_ran_gaussian_tail_pdf(x, a, sigma) bind(c)
    import
    real(c_double), value :: x, a
    real(c_double), value :: sigma
    real(c_double) :: gsl_ran_gaussian_tail_pdf
  end function gsl_ran_gaussian_tail_pdf
  function gsl_ran_ugaussian_tail(r, a) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a
    real(c_double) :: gsl_ran_ugaussian_tail
  end function gsl_ran_ugaussian_tail
  function gsl_ran_ugaussian_tail_pdf(x, a) bind(c)
    import
    real(c_double), value :: x, a
    real(c_double) :: gsl_ran_ugaussian_tail_pdf
  end function gsl_ran_ugaussian_tail_pdf
  subroutine gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y, rho, x, y) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: sigma_x, sigma_y, rho
    real(c_double), intent(out) :: x, y
  end subroutine gsl_ran_bivariate_gaussian
  function gsl_ran_bivariate_gaussian_pdf(x, y, sigma_x, sigma_y, rho) bind(c)
    import
    real(c_double), value :: x, y, sigma_x, sigma_y, rho
    real(c_double) :: gsl_ran_bivariate_gaussian_pdf
  end function gsl_ran_bivariate_gaussian_pdf
  function gsl_ran_exponential(r, mu) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: mu
    real(c_double) :: gsl_ran_exponential
  end function gsl_ran_exponential
  function gsl_ran_exponential_pdf(x, mu) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: mu
    real(c_double) :: gsl_ran_exponential_pdf
  end function gsl_ran_exponential_pdf
  function gsl_cdf_exponential_p(x, mu) bind(c, name='gsl_cdf_exponential_P')
    import
    real(c_double), value :: x
    real(c_double), value :: mu
    real(c_double) :: gsl_cdf_exponential_p
  end function gsl_cdf_exponential_p
  function gsl_cdf_exponential_q(x, mu) bind(c, name='gsl_cdf_exponential_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: mu
    real(c_double) :: gsl_cdf_exponential_q
  end function gsl_cdf_exponential_q
  function gsl_cdf_exponential_pinv(p, mu) bind(c, name='gsl_cdf_exponential_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: mu
    real(c_double) :: gsl_cdf_exponential_pinv
  end function gsl_cdf_exponential_pinv
  function gsl_cdf_exponential_qinv(q, mu) bind(c, name='gsl_cdf_exponential_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: mu
    real(c_double) :: gsl_cdf_exponential_qinv
  end function gsl_cdf_exponential_qinv
  function gsl_ran_laplace(r, a) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a
    real(c_double) :: gsl_ran_laplace
  end function gsl_ran_laplace
  function gsl_ran_laplace_pdf(x, a) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: a
    real(c_double) :: gsl_ran_laplace_pdf
  end function gsl_ran_laplace_pdf
  function gsl_cdf_laplace_p(x, a) bind(c, name='gsl_cdf_laplace_P')
    import
    real(c_double), value :: x
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_laplace_p
  end function gsl_cdf_laplace_p
  function gsl_cdf_laplace_q(x, a) bind(c, name='gsl_cdf_laplace_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_laplace_q
  end function gsl_cdf_laplace_q
  function gsl_cdf_laplace_pinv(p, a) bind(c, name='gsl_cdf_laplace_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_laplace_pinv
  end function gsl_cdf_laplace_pinv
  function gsl_cdf_laplace_qinv(q, a) bind(c, name='gsl_cdf_laplace_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_laplace_qinv
  end function gsl_cdf_laplace_qinv
  function gsl_ran_exppow(r, a, b) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_exppow
  end function gsl_ran_exppow
  function gsl_ran_exppow_pdf(x, a, b) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_exppow_pdf
  end function gsl_ran_exppow_pdf
  function gsl_cdf_exppow_p(x, a, b) bind(c, name='gsl_cdf_exppow_P')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_exppow_p
  end function gsl_cdf_exppow_p
  function gsl_cdf_exppow_q(x, a, b) bind(c, name='gsl_cdf_exppow_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_exppow_q
  end function gsl_cdf_exppow_q
  function gsl_ran_cauchy(r, a) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a
    real(c_double) :: gsl_ran_cauchy
  end function gsl_ran_cauchy
  function gsl_ran_cauchy_pdf(x, a) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: a
    real(c_double) :: gsl_ran_cauchy_pdf
  end function gsl_ran_cauchy_pdf
  function gsl_cdf_cauchy_p(x, a) bind(c, name='gsl_cdf_cauchy_P')
    import
    real(c_double), value :: x
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_cauchy_p
  end function gsl_cdf_cauchy_p
  function gsl_cdf_cauchy_q(x, a) bind(c, name='gsl_cdf_cauchy_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_cauchy_q
  end function gsl_cdf_cauchy_q
  function gsl_cdf_cauchy_pinv(p, a) bind(c, name='gsl_cdf_cauchy_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_cauchy_pinv
  end function gsl_cdf_cauchy_pinv
  function gsl_cdf_cauchy_qinv(q, a) bind(c, name='gsl_cdf_cauchy_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_cauchy_qinv
  end function gsl_cdf_cauchy_qinv
  function gsl_ran_rayleigh(r, sigma) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: sigma
    real(c_double) :: gsl_ran_rayleigh
  end function gsl_ran_rayleigh
  function gsl_ran_rayleigh_pdf(x, sigma) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: sigma
    real(c_double) :: gsl_ran_rayleigh_pdf
  end function gsl_ran_rayleigh_pdf
  function gsl_cdf_rayleigh_p(x, sigma) bind(c, name='gsl_cdf_rayleigh_P')
    import
    real(c_double), value :: x
    real(c_double), value :: sigma
    real(c_double) :: gsl_cdf_rayleigh_p
  end function gsl_cdf_rayleigh_p
  function gsl_cdf_rayleigh_q(x, sigma) bind(c, name='gsl_cdf_rayleigh_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: sigma
    real(c_double) :: gsl_cdf_rayleigh_q
  end function gsl_cdf_rayleigh_q
  function gsl_cdf_rayleigh_pinv(p, sigma) bind(c, name='gsl_cdf_rayleigh_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: sigma
    real(c_double) :: gsl_cdf_rayleigh_pinv
  end function gsl_cdf_rayleigh_pinv
  function gsl_cdf_rayleigh_qinv(q, sigma) bind(c, name='gsl_cdf_rayleigh_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: sigma
    real(c_double) :: gsl_cdf_rayleigh_qinv
  end function gsl_cdf_rayleigh_qinv
  function gsl_ran_rayleigh_tail(r, a, sigma) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a, sigma
    real(c_double) :: gsl_ran_rayleigh_tail
  end function gsl_ran_rayleigh_tail
  function gsl_ran_rayleigh_tail_pdf(x, a, sigma) bind(c)
    import
    real(c_double), value :: x, a
    real(c_double), value :: sigma
    real(c_double) :: gsl_ran_rayleigh_tail_pdf
  end function gsl_ran_rayleigh_tail_pdf
  function gsl_ran_landau(r) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double) :: gsl_ran_landau
  end function gsl_ran_landau
  function gsl_ran_landau_pdf(x) bind(c)
    import
    real(c_double), value :: x
    real(c_double) :: gsl_ran_landau_pdf
  end function gsl_ran_landau_pdf
  function gsl_ran_levy(r, c, alpha) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: c, alpha
    real(c_double) :: gsl_ran_levy
  end function gsl_ran_levy
  function gsl_ran_levy_skew(r, c, alpha, beta) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: c, alpha, beta
    real(c_double) :: gsl_ran_levy_skew
  end function gsl_ran_levy_skew
  function gsl_ran_gamma(r, a, b) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_gamma
  end function gsl_ran_gamma
  function gsl_ran_gamma_mt(r, a, b) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_gamma_mt
  end function gsl_ran_gamma_mt
  function gsl_ran_gamma_pdf(x, a, b) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_gamma_pdf
  end function gsl_ran_gamma_pdf
  function gsl_cdf_gamma_p(x, a, b) bind(c, name='gsl_cdf_gamma_P')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gamma_p
  end function gsl_cdf_gamma_p
  function gsl_cdf_gamma_q(x, a, b) bind(c, name='gsl_cdf_gamma_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gamma_q
  end function gsl_cdf_gamma_q
  function gsl_cdf_gamma_pinv(p, a, b) bind(c, name='gsl_cdf_gamma_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gamma_pinv
  end function gsl_cdf_gamma_pinv
  function gsl_cdf_gamma_qinv(q, a, b) bind(c, name='gsl_cdf_gamma_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gamma_qinv
  end function gsl_cdf_gamma_qinv
  function gsl_ran_flat(r, a, b) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_flat
  end function gsl_ran_flat
  function gsl_ran_flat_pdf(x, a, b) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_flat_pdf
  end function gsl_ran_flat_pdf
  function gsl_cdf_flat_p(x, a, b) bind(c, name='gsl_cdf_flat_P')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_flat_p
  end function gsl_cdf_flat_p
  function gsl_cdf_flat_q(x, a, b) bind(c, name='gsl_cdf_flat_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_flat_q
  end function gsl_cdf_flat_q
  function gsl_cdf_flat_pinv(p, a, b) bind(c, name='gsl_cdf_flat_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_flat_pinv
  end function gsl_cdf_flat_pinv
  function gsl_cdf_flat_qinv(q, a, b) bind(c, name='gsl_cdf_flat_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_flat_qinv
  end function gsl_cdf_flat_qinv
  function gsl_ran_lognormal(r, zeta, sigma) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: zeta, sigma
    real(c_double) :: gsl_ran_lognormal
  end function gsl_ran_lognormal
  function gsl_ran_lognormal_pdf(x, zeta, sigma) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: zeta, sigma
    real(c_double) :: gsl_ran_lognormal_pdf
  end function gsl_ran_lognormal_pdf
  function gsl_cdf_lognormal_p(x, zeta, sigma) bind(c, name='gsl_cdf_lognormal_P')
    import
    real(c_double), value :: x
    real(c_double), value :: zeta, sigma
    real(c_double) :: gsl_cdf_lognormal_p
  end function gsl_cdf_lognormal_p
  function gsl_cdf_lognormal_q(x, zeta, sigma) bind(c, name='gsl_cdf_lognormal_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: zeta, sigma
    real(c_double) :: gsl_cdf_lognormal_q
  end function gsl_cdf_lognormal_q
  function gsl_cdf_lognormal_pinv(p, zeta, sigma) bind(c, name='gsl_cdf_lognormal_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: zeta, sigma
    real(c_double) :: gsl_cdf_lognormal_pinv
  end function gsl_cdf_lognormal_pinv
  function gsl_cdf_lognormal_qinv(q, zeta, sigma) bind(c, name='gsl_cdf_lognormal_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: zeta, sigma
    real(c_double) :: gsl_cdf_lognormal_qinv
  end function gsl_cdf_lognormal_qinv
  function gsl_ran_chisq(r, nu) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: nu
    real(c_double) :: gsl_ran_chisq
  end function gsl_ran_chisq
  function gsl_ran_chisq_pdf(x, nu) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: nu
    real(c_double) :: gsl_ran_chisq_pdf
  end function gsl_ran_chisq_pdf
  function gsl_cdf_chisq_p(x, nu) bind(c, name='gsl_cdf_chisq_P')
    import
    real(c_double), value :: x
    real(c_double), value :: nu
    real(c_double) :: gsl_cdf_chisq_p
  end function gsl_cdf_chisq_p
  function gsl_cdf_chisq_q(x, nu) bind(c, name='gsl_cdf_chisq_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: nu
    real(c_double) :: gsl_cdf_chisq_q
  end function gsl_cdf_chisq_q
  function gsl_cdf_chisq_pinv(p, nu) bind(c, name='gsl_cdf_chisq_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: nu
    real(c_double) :: gsl_cdf_chisq_pinv
  end function gsl_cdf_chisq_pinv
  function gsl_cdf_chisq_qinv(q, nu) bind(c, name='gsl_cdf_chisq_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: nu
    real(c_double) :: gsl_cdf_chisq_qinv
  end function gsl_cdf_chisq_qinv
  function gsl_ran_fdist(r, nu1, nu2) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: nu1, nu2
    real(c_double) :: gsl_ran_fdist
  end function gsl_ran_fdist
  function gsl_ran_fdist_pdf(x, nu1, nu2) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: nu1, nu2
    real(c_double) :: gsl_ran_fdist_pdf
  end function gsl_ran_fdist_pdf
  function gsl_cdf_fdist_p(x, nu1, nu2) bind(c, name='gsl_cdf_fdist_P')
    import
    real(c_double), value :: x
    real(c_double), value :: nu1, nu2
    real(c_double) :: gsl_cdf_fdist_p
  end function gsl_cdf_fdist_p
  function gsl_cdf_fdist_q(x, nu1, nu2) bind(c, name='gsl_cdf_fdist_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: nu1, nu2
    real(c_double) :: gsl_cdf_fdist_q
  end function gsl_cdf_fdist_q
  function gsl_cdf_fdist_pinv(p, nu1, nu2) bind(c, name='gsl_cdf_fdist_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: nu1, nu2
    real(c_double) :: gsl_cdf_fdist_pinv
  end function gsl_cdf_fdist_pinv
  function gsl_cdf_fdist_qinv(q, nu1, nu2) bind(c, name='gsl_cdf_fdist_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: nu1, nu2
    real(c_double) :: gsl_cdf_fdist_qinv
  end function gsl_cdf_fdist_qinv
  function gsl_ran_tdist(r, nu) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: nu
    real(c_double) :: gsl_ran_tdist
  end function gsl_ran_tdist
  function gsl_ran_tdist_pdf(x, nu) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: nu
    real(c_double) :: gsl_ran_tdist_pdf
  end function gsl_ran_tdist_pdf
  function gsl_cdf_tdist_p(x, nu) bind(c, name='gsl_cdf_tdist_P')
    import
    real(c_double), value :: x
    real(c_double), value :: nu
    real(c_double) :: gsl_cdf_tdist_p
  end function gsl_cdf_tdist_p
  function gsl_cdf_tdist_q(x, nu) bind(c, name='gsl_cdf_tdist_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: nu
    real(c_double) :: gsl_cdf_tdist_q
  end function gsl_cdf_tdist_q
  function gsl_cdf_tdist_pinv(p, nu) bind(c, name='gsl_cdf_tdist_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: nu
    real(c_double) :: gsl_cdf_tdist_pinv
  end function gsl_cdf_tdist_pinv
  function gsl_cdf_tdist_qinv(q, nu) bind(c, name='gsl_cdf_tdist_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: nu
    real(c_double) :: gsl_cdf_tdist_qinv
  end function gsl_cdf_tdist_qinv
  function gsl_ran_beta(r, a, b) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_beta
  end function gsl_ran_beta
  function gsl_ran_beta_pdf(x, a, b) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_beta_pdf
  end function gsl_ran_beta_pdf
  function gsl_cdf_beta_p(x, a, b) bind(c, name='gsl_cdf_beta_P')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_beta_p
  end function gsl_cdf_beta_p
  function gsl_cdf_beta_q(x, a, b) bind(c, name='gsl_cdf_beta_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_beta_q
  end function gsl_cdf_beta_q
  function gsl_cdf_beta_pinv(p, a, b) bind(c, name='gsl_cdf_beta_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_beta_pinv
  end function gsl_cdf_beta_pinv
  function gsl_cdf_beta_qinv(q, a, b) bind(c, name='gsl_cdf_beta_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_beta_qinv
  end function gsl_cdf_beta_qinv
  function gsl_ran_logistic(r, a) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a
    real(c_double) :: gsl_ran_logistic
  end function gsl_ran_logistic
  function gsl_ran_logistic_pdf(x, a) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: a
    real(c_double) :: gsl_ran_logistic_pdf
  end function gsl_ran_logistic_pdf
  function gsl_cdf_logistic_p(x, a) bind(c, name='gsl_cdf_logistic_P')
    import
    real(c_double), value :: x
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_logistic_p
  end function gsl_cdf_logistic_p
  function gsl_cdf_logistic_q(x, a) bind(c, name='gsl_cdf_logistic_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_logistic_q
  end function gsl_cdf_logistic_q
  function gsl_cdf_logistic_pinv(p, a) bind(c, name='gsl_cdf_logistic_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_logistic_pinv
  end function gsl_cdf_logistic_pinv
  function gsl_cdf_logistic_qinv(q, a) bind(c, name='gsl_cdf_logistic_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: a
    real(c_double) :: gsl_cdf_logistic_qinv
  end function gsl_cdf_logistic_qinv
  function gsl_ran_pareto(r, a, b) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_pareto
  end function gsl_ran_pareto
  function gsl_ran_pareto_pdf(x, a, b) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_pareto_pdf
  end function gsl_ran_pareto_pdf
  function gsl_cdf_pareto_p(x, a, b) bind(c, name='gsl_cdf_pareto_P')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_pareto_p
  end function gsl_cdf_pareto_p
  function gsl_cdf_pareto_q(x, a, b) bind(c, name='gsl_cdf_pareto_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_pareto_q
  end function gsl_cdf_pareto_q
  function gsl_cdf_pareto_pinv(p, a, b) bind(c, name='gsl_cdf_pareto_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_pareto_pinv
  end function gsl_cdf_pareto_pinv
  function gsl_cdf_pareto_qinv(q, a, b) bind(c, name='gsl_cdf_pareto_Qinv')
    import
    real(c_double), value :: q
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_pareto_qinv
  end function gsl_cdf_pareto_qinv
  subroutine gsl_ran_dir_2d(r, x, y) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), intent(out) :: x, y
  end subroutine gsl_ran_dir_2d
  subroutine gsl_ran_dir_2d_trig_method(r, x, y) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), intent(out) :: x, y
  end subroutine gsl_ran_dir_2d_trig_method
  subroutine gsl_ran_dir_3d(r, x, y, z) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), intent(out) :: x, y, z
  end subroutine gsl_ran_dir_3d
  subroutine gsl_ran_dir_nd(r, n, x) bind(c)
    import
    type(c_ptr), value :: r
    integer(c_size_t), value :: n
    real(c_double), intent(out) :: x
  end subroutine gsl_ran_dir_nd
  function gsl_ran_weibull(r, a, b) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_weibull
  end function gsl_ran_weibull
  function gsl_ran_weibull_pdf(x, a, b) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_weibull_pdf
  end function gsl_ran_weibull_pdf
  function gsl_cdf_weibull_p(x, a, b) bind(c, name='gsl_cdf_weibull_P')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_weibull_p
  end function gsl_cdf_weibull_p
  function gsl_cdf_weibull_q(x, a, b) bind(c, name='gsl_cdf_weibull_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_weibull_q
  end function gsl_cdf_weibull_q
  function gsl_cdf_weibull_pinv(p, a, b) bind(c, name='gsl_cdf_weibull_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_weibull_pinv
  end function gsl_cdf_weibull_pinv
  function gsl_cdf_weibull_qinv(p, a, b) bind(c, name='gsl_cdf_weibull_Qinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_weibull_qinv
  end function gsl_cdf_weibull_qinv
  function gsl_ran_gumbel1(r, a, b) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_gumbel1
  end function gsl_ran_gumbel1
  function gsl_ran_gumbel1_pdf(x, a, b) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_gumbel1_pdf
  end function gsl_ran_gumbel1_pdf
  function gsl_cdf_gumbel1_p(x, a, b) bind(c, name='gsl_cdf_gumbel1_P')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gumbel1_p
  end function gsl_cdf_gumbel1_p
  function gsl_cdf_gumbel1_q(x, a, b) bind(c, name='gsl_cdf_gumbel1_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gumbel1_q
  end function gsl_cdf_gumbel1_q
  function gsl_cdf_gumbel1_pinv(p, a, b) bind(c, name='gsl_cdf_gumbel1_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gumbel1_pinv
  end function gsl_cdf_gumbel1_pinv
  function gsl_cdf_gumbel1_qinv(p, a, b) bind(c, name='gsl_cdf_gumbel1_Qinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gumbel1_qinv
  end function gsl_cdf_gumbel1_qinv
  function gsl_ran_gumbel2(r, a, b) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_gumbel2
  end function gsl_ran_gumbel2
  function gsl_ran_gumbel2_pdf(x, a, b) bind(c)
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_ran_gumbel2_pdf
  end function gsl_ran_gumbel2_pdf
  function gsl_cdf_gumbel2_p(x, a, b) bind(c, name='gsl_cdf_gumbel2_P')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gumbel2_p
  end function gsl_cdf_gumbel2_p
  function gsl_cdf_gumbel2_q(x, a, b) bind(c, name='gsl_cdf_gumbel2_Q')
    import
    real(c_double), value :: x
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gumbel2_q
  end function gsl_cdf_gumbel2_q
  function gsl_cdf_gumbel2_pinv(p, a, b) bind(c, name='gsl_cdf_gumbel2_Pinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gumbel2_pinv
  end function gsl_cdf_gumbel2_pinv
  function gsl_cdf_gumbel2_qinv(p, a, b) bind(c, name='gsl_cdf_gumbel2_Qinv')
    import
    real(c_double), value :: p
    real(c_double), value :: a, b
    real(c_double) :: gsl_cdf_gumbel2_qinv
  end function gsl_cdf_gumbel2_qinv
  subroutine gsl_ran_dirichlet(r, k, alpha, theta) bind(c)
    import
    type(c_ptr), value :: r
    integer(c_size_t), value :: k
    real(c_double), dimension(*), intent(in) :: alpha
    real(c_double), dimension(*), intent(out) :: theta
  end subroutine gsl_ran_dirichlet
  function gsl_ran_dirichlet_pdf(k, alpha, theta) bind(c)
    import
    integer(c_size_t), value :: k
    real(c_double), dimension(*), intent(in) :: alpha
    real(c_double), dimension(*), intent(in) :: theta
    real(c_double) :: gsl_ran_dirichlet_pdf
  end function gsl_ran_dirichlet_pdf
  function gsl_ran_dirichlet_lnpdf(k, alpha, theta) bind(c)
    import
    integer(c_size_t), value :: k
    real(c_double), dimension(*), intent(in) :: alpha
    real(c_double), dimension(*), intent(in) :: theta
    real(c_double) :: gsl_ran_dirichlet_lnpdf
  end function gsl_ran_dirichlet_lnpdf
  function gsl_ran_discrete_preproc(k, p) bind(c)
    import
    integer(c_size_t), value :: k
    real(c_double), dimension(*), intent(in) :: p
    type(c_ptr) :: gsl_ran_discrete_preproc
  end function gsl_ran_discrete_preproc
  function gsl_ran_discrete(r, g) bind(c)
    import
    type(c_ptr), value :: r, g
    integer(c_size_t) :: gsl_ran_discrete
  end function gsl_ran_discrete
  function gsl_ran_discrete_pdf(k, g) bind(c)
    import
    integer(c_size_t), value :: k
    type(c_ptr), value :: g
    real(c_double) :: gsl_ran_discrete_pdf
  end function gsl_ran_discrete_pdf
  subroutine gsl_ran_discrete_free(g) bind(c)
    import
    type(c_ptr), value :: g
  end subroutine gsl_ran_discrete_free
  function gsl_ran_poisson(r, mu) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: mu
    integer(c_int) :: gsl_ran_poisson
  end function gsl_ran_poisson
  function gsl_ran_poisson_pdf(k, mu) bind(c)
    import
    integer(c_int), value :: k
    real(c_double), value :: mu
    real(c_double) :: gsl_ran_poisson_pdf
  end function gsl_ran_poisson_pdf
  function gsl_cdf_poisson_p(k, mu) bind(c, name='gsl_cdf_poisson_P')
    import
    integer(c_int), value :: k
    real(c_double), value :: mu
    real(c_double) :: gsl_cdf_poisson_p
  end function gsl_cdf_poisson_p
  function gsl_cdf_poisson_q(k, mu) bind(c, name='gsl_cdf_poisson_Q')
    import
    integer(c_int), value :: k
    real(c_double), value :: mu
    real(c_double) :: gsl_cdf_poisson_q
  end function gsl_cdf_poisson_q
  function gsl_ran_bernoulli(r, p) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: p
    integer(c_int) :: gsl_ran_bernoulli
  end function gsl_ran_bernoulli
  function gsl_ran_bernoulli_pdf(k, p) bind(c)
    import
    integer(c_int), value :: k
    real(c_double), value :: p
    real(c_double) :: gsl_ran_bernoulli_pdf
  end function gsl_ran_bernoulli_pdf
  function gsl_ran_binomial(r, p, n) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: p
    integer(c_int), value :: n
    real(c_double) :: gsl_ran_binomial
  end function gsl_ran_binomial
  function gsl_ran_binomial_pdf(k, p, n) bind(c)
    import
    integer(c_int), value :: k, n
    real(c_double), value :: p
    real(c_double) :: gsl_ran_binomial_pdf
  end function gsl_ran_binomial_pdf
  function gsl_cdf_binomial_p(k, p, n) bind(c, name='gsl_cdf_binomial_P')
    import
    integer(c_int), value :: k, n
    real(c_double), value :: p
    real(c_double) :: gsl_cdf_binomial_p
  end function gsl_cdf_binomial_p
  function gsl_cdf_binomial_q(k, p, n) bind(c, name='gsl_cdf_binomial_Q')
    import
    integer(c_int), value :: k, n
    real(c_double), value :: p
    real(c_double) :: gsl_cdf_binomial_q
  end function gsl_cdf_binomial_q
  subroutine gsl_ran_multinomial(r, k, nn, p, n) bind(c)
    import
    type(c_ptr), value :: r
    integer(c_size_t), value :: k
    integer(c_int), value :: nn
    real(c_double), dimension(*), intent(in) :: p
    integer(c_int), dimension(*), intent(out) :: n
  end subroutine gsl_ran_multinomial
  function gsl_ran_multinomial_pdf(k, p, n) bind(c)
    import
    integer(c_size_t), value :: k
    real(c_double), dimension(*), intent(in) :: p
    integer(c_int), dimension(*), intent(in) :: n
    real(c_double) :: gsl_ran_multinomial_pdf
  end function gsl_ran_multinomial_pdf
  function gsl_ran_multinomial_lnpdf(k, p, n) bind(c)
    import
    integer(c_size_t), value :: k
    real(c_double), dimension(*), intent(in) :: p
    integer(c_int), dimension(*), intent(in) :: n
    real(c_double) :: gsl_ran_multinomial_lnpdf
  end function gsl_ran_multinomial_lnpdf
  function gsl_ran_negative_binomial(r, p, n) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: p, n
    integer(c_int) :: gsl_ran_negative_binomial
  end function gsl_ran_negative_binomial
  function gsl_ran_negative_binomial_pdf(k, p, n) bind(c)
    import
    integer(c_int), value :: k
    real(c_double), value :: p, n
    real(c_double) :: gsl_ran_negative_binomial_pdf
  end function gsl_ran_negative_binomial_pdf
  function gsl_cdf_negative_binomial_p(k, p, n) bind(c, name='gsl_cdf_negative_binomial_P')
    import
    integer(c_int), value :: k
    real(c_double), value :: p, n
    real(c_double) :: gsl_cdf_negative_binomial_p
  end function gsl_cdf_negative_binomial_p
  function gsl_cdf_negative_binomial_q(k, p, n) bind(c, name='gsl_cdf_negative_binomial_Q')
    import
    integer(c_int), value :: k
    real(c_double), value :: p, n
    real(c_double) :: gsl_cdf_negative_binomial_q
  end function gsl_cdf_negative_binomial_q
  function gsl_ran_pascal(r, p, n) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: p, n
    integer(c_int) :: gsl_ran_pascal
  end function gsl_ran_pascal
  function gsl_ran_pascal_pdf(k, p, n) bind(c)
    import
    integer(c_int), value :: k
    real(c_double), value :: p, n
    real(c_double) :: gsl_ran_pascal_pdf
  end function gsl_ran_pascal_pdf
  function gsl_cdf_pascal_p(k, p, n) bind(c, name='gsl_cdf_pascal_P')
    import
    integer(c_int), value :: k
    real(c_double), value :: p, n
    real(c_double) :: gsl_cdf_pascal_p
  end function gsl_cdf_pascal_p
  function gsl_cdf_pascal_q(k, p, n) bind(c, name='gsl_cdf_pascal_Q')
    import
    integer(c_int), value :: k
    real(c_double), value :: p, n
    real(c_double) :: gsl_cdf_pascal_q
  end function gsl_cdf_pascal_q
  function gsl_ran_geometric(r, p) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: p
    integer(c_int) :: gsl_ran_geometric
  end function gsl_ran_geometric
  function gsl_ran_geometric_pdf(k, p) bind(c)
    import
    integer(c_int), value :: k
    real(c_double), value :: p
    real(c_double) :: gsl_ran_geometric_pdf
  end function gsl_ran_geometric_pdf
  function gsl_cdf_geometric_p(k, p) bind(c, name='gsl_cdf_geometric_P')
    import
    integer(c_int), value :: k
    real(c_double), value :: p
    real(c_double) :: gsl_cdf_geometric_p
  end function gsl_cdf_geometric_p
  function gsl_cdf_geometric_q(k, p) bind(c, name='gsl_cdf_geometric_Q')
    import
    integer(c_int), value :: k
    real(c_double), value :: p
    real(c_double) :: gsl_cdf_geometric_q
  end function gsl_cdf_geometric_q
  function gsl_ran_hypergeometric(r, n1, n2, t) bind(c)
    import
    type(c_ptr), value :: r
    integer(c_int), value :: n1, n2, t
    integer(c_int) :: gsl_ran_hypergeometric
  end function gsl_ran_hypergeometric
  function gsl_ran_hypergeometric_pdf(k, n1, n2, t) bind(c)
    import
    integer(c_int), value :: k
    integer(c_int), value :: n1, n2, t
    real(c_double) :: gsl_ran_hypergeometric_pdf
  end function gsl_ran_hypergeometric_pdf
  function gsl_cdf_hypergeometric_p(k, n1, n2, t) bind(c, name='gsl_cdf_hypergeometric_P')
    import
    integer(c_int), value :: k
    integer(c_int), value :: n1, n2, t
    real(c_double) :: gsl_cdf_hypergeometric_p
  end function gsl_cdf_hypergeometric_p
  function gsl_cdf_hypergeometric_q(k, n1, n2, t) bind(c, name='gsl_cdf_hypergeometric_Q')
    import
    integer(c_int), value :: k
    integer(c_int), value :: n1, n2, t
    real(c_double) :: gsl_cdf_hypergeometric_q
  end function gsl_cdf_hypergeometric_q
  function gsl_ran_logarithmic(r, p) bind(c)
    import
    type(c_ptr), value :: r
    real(c_double), value :: p
    integer(c_int) :: gsl_ran_logarithmic
  end function gsl_ran_logarithmic
  function gsl_ran_logarithmic_pdf(k, p) bind(c)
    import
    integer(c_int), value :: k
    real(c_double), value :: p
    real(c_double) :: gsl_ran_logarithmic_pdf
  end function gsl_ran_logarithmic_pdf
  subroutine gsl_ran_shuffle(r, base, n, size) bind(c)
    import :: c_ptr, c_size_t
    type(c_ptr), value :: r, base
    integer(c_size_t), value :: n, size
  end subroutine gsl_ran_shuffle
  function gsl_ran_choose(r, dest, k, src, n, size) bind(c)
    import :: c_ptr, c_size_t, c_int
    type(c_ptr), value :: r, dest, src
    integer(c_size_t), value :: k, n, size
    integer(c_int) :: gsl_ran_choose
  end function gsl_ran_choose
  subroutine gsl_ran_sample(r, dest, k, src, n, size) bind(c)
    import :: c_ptr, c_size_t
    type(c_ptr), value :: r, dest, src
    integer(c_size_t), value :: k, n, size
  end subroutine gsl_ran_sample
!-*-f90-*-
!
! Interfaces: Statistics
!
  function gsl_stats_mean (data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double) :: gsl_stats_mean
  end function gsl_stats_mean
  function gsl_stats_variance (data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double) :: gsl_stats_variance
  end function gsl_stats_variance
  function gsl_stats_variance_m (data, stride, n, mean) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double), value :: mean
    real(c_double) :: gsl_stats_variance_m
  end function gsl_stats_variance_m
  function gsl_stats_sd (data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double) :: gsl_stats_sd
  end function gsl_stats_sd
  function gsl_stats_sd_m (data, stride, n, mean) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double), value :: mean
    real(c_double) :: gsl_stats_sd_m
  end function gsl_stats_sd_m
  function gsl_stats_variance_with_fixed_mean (data, stride, n, mean) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double), value :: mean
    real(c_double) :: gsl_stats_variance_with_fixed_mean
  end function gsl_stats_variance_with_fixed_mean
  function gsl_stats_sd_with_fixed_mean (data, stride, n, mean) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double), value :: mean
    real(c_double) :: gsl_stats_sd_with_fixed_mean
  end function gsl_stats_sd_with_fixed_mean
  function gsl_stats_absdev (data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double) :: gsl_stats_absdev
  end function gsl_stats_absdev
  function gsl_stats_absdev_m (data, stride, n, mean) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double), value :: mean
    real(c_double) :: gsl_stats_absdev_m
  end function gsl_stats_absdev_m
  function gsl_stats_skew (data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double) :: gsl_stats_skew
  end function gsl_stats_skew
  function gsl_stats_skew_m_sd (data, stride, n, mean, sd) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double), value :: mean, sd
    real(c_double) :: gsl_stats_skew_m_sd
  end function gsl_stats_skew_m_sd
  function gsl_stats_kurtosis (data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double) :: gsl_stats_kurtosis
  end function gsl_stats_kurtosis
  function gsl_stats_kurtosis_m_sd (data, stride, n, mean, sd) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double), value :: mean, sd
    real(c_double) :: gsl_stats_kurtosis_m_sd
  end function gsl_stats_kurtosis_m_sd
  function gsl_stats_lag1_autocorrelation (data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double) :: gsl_stats_lag1_autocorrelation
  end function gsl_stats_lag1_autocorrelation
  function gsl_stats_lag1_autocorrelation_m (data, stride, n, mean) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double), value :: mean
    real(c_double) :: gsl_stats_lag1_autocorrelation_m
  end function gsl_stats_lag1_autocorrelation_m
  function gsl_stats_covariance(data1, stride1, data2, stride2, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data1, data2
    integer(c_size_t), value :: stride1, stride2, n
    real(c_double) :: gsl_stats_covariance
  end function gsl_stats_covariance
  function gsl_stats_covariance_m(data1, stride1, data2, stride2, n, &
       mean1, mean2) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data1, data2
    integer(c_size_t), value :: stride1, stride2, n
    real(c_double), value :: mean1, mean2
    real(c_double) :: gsl_stats_covariance_m
  end function gsl_stats_covariance_m
  function gsl_stats_correlation(data1, stride1, data2, stride2, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data1, data2
    integer(c_size_t), value :: stride1, stride2, n
    real(c_double) :: gsl_stats_correlation
  end function gsl_stats_correlation
  function gsl_stats_wmean(w, wstride, data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double) :: gsl_stats_wmean
  end function gsl_stats_wmean
  function gsl_stats_wvariance(w, wstride, data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double) :: gsl_stats_wvariance
  end function gsl_stats_wvariance
  function gsl_stats_wvariance_m(w, wstride, data, stride, n, mean) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double), value :: mean
    real(c_double) :: gsl_stats_wvariance_m
  end function gsl_stats_wvariance_m
  function gsl_stats_wsd(w, wstride, data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double) :: gsl_stats_wsd
  end function gsl_stats_wsd
  function gsl_stats_wsd_m(w, wstride, data, stride, n, mean) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double), value :: mean
    real(c_double) :: gsl_stats_wsd_m
  end function gsl_stats_wsd_m
  function gsl_stats_wvariance_with_fixed_mean(w, wstride, data, stride, n, mean) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double), value :: mean
    real(c_double) :: gsl_stats_wvariance_with_fixed_mean
  end function gsl_stats_wvariance_with_fixed_mean
  function gsl_stats_wsd_with_fixed_mean(w, wstride, data, stride, n, mean) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double), value :: mean
    real(c_double) :: gsl_stats_wsd_with_fixed_mean
  end function gsl_stats_wsd_with_fixed_mean
  function gsl_stats_wabsdev(w, wstride, data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double) :: gsl_stats_wabsdev
  end function gsl_stats_wabsdev
  function gsl_stats_wabsdev_m(w, wstride, data, stride, n, mean) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double), value :: mean
    real(c_double) :: gsl_stats_wabsdev_m
  end function gsl_stats_wabsdev_m
  function gsl_stats_wskew(w, wstride, data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double) :: gsl_stats_wskew
  end function gsl_stats_wskew
  function gsl_stats_wskew_m_sd(w, wstride, data, stride, n, mean, sd) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double), value :: mean, sd
    real(c_double) :: gsl_stats_wskew_m_sd
  end function gsl_stats_wskew_m_sd
  function gsl_stats_wkurtosis(w, wstride, data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double) :: gsl_stats_wkurtosis
  end function gsl_stats_wkurtosis
  function gsl_stats_wkurtosis_m_sd(w, wstride, data, stride, n, mean, sd) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: w, data
    integer(c_size_t), value :: wstride, stride, n
    real(c_double), value :: mean, sd
    real(c_double) :: gsl_stats_wkurtosis_m_sd
  end function gsl_stats_wkurtosis_m_sd
  function gsl_stats_max(data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double) :: gsl_stats_max
  end function gsl_stats_max
  function gsl_stats_min(data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double) :: gsl_stats_min
  end function gsl_stats_min
  subroutine gsl_stats_minmax(min, max, data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double), intent(out) :: min, max
  end subroutine gsl_stats_minmax
  function gsl_stats_max_index(data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    integer(c_size_t) :: gsl_stats_max_index
  end function gsl_stats_max_index
  function gsl_stats_min_index(data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    integer(c_size_t) :: gsl_stats_min_index
  end function gsl_stats_min_index
  subroutine gsl_stats_minmax_index(min_index, max_index, data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    integer(c_size_t), intent(out) :: min_index, max_index
  end subroutine gsl_stats_minmax_index
  function gsl_stats_median_from_sorted_data(data, stride, n) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    integer(c_size_t), value :: stride, n
    real(c_double) :: gsl_stats_median_from_sorted_data
  end function gsl_stats_median_from_sorted_data
  function gsl_stats_quantile_from_sorted_data(data, stride, n, f) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: data
    real(c_double), value :: f
    integer(c_size_t), value :: stride, n
    real(c_double) :: gsl_stats_quantile_from_sorted_data
  end function gsl_stats_quantile_from_sorted_data
!-*-f90-*-
!
! Interfaces: Histograms
!
  function gsl_histogram_alloc(n) bind(c)
    import
    integer(c_size_t), value :: n
    type(c_ptr) :: gsl_histogram_alloc
  end function gsl_histogram_alloc
  function gsl_histogram_set_ranges(h, range, size) bind(c)
    import 
    type(c_ptr), value :: h
    integer(c_size_t), value :: size
    real(c_double), dimension(*), intent(in) :: range
    integer(c_int) :: gsl_histogram_set_ranges
  end function gsl_histogram_set_ranges
  function gsl_histogram_set_ranges_uniform(h, xmin, xmax) bind(c)
    import
    type(c_ptr), value :: h
    real(c_double), value :: xmin, xmax
    integer(c_int) :: gsl_histogram_set_ranges_uniform
  end function gsl_histogram_set_ranges_uniform
  subroutine gsl_histogram_free(h) bind(c)
    import
    type(c_ptr), value :: h
  end subroutine gsl_histogram_free
  function gsl_histogram_memcpy(dest, src) bind(c)
    import
    type(c_ptr), value :: dest, src
    integer(c_int) :: gsl_histogram_memcpy
  end function gsl_histogram_memcpy
  function gsl_histogram_clone(src) bind(c)
    import
    type(c_ptr), value :: src
    type(c_ptr) :: gsl_histogram_clone
  end function gsl_histogram_clone
  function gsl_histogram_increment(h, x) bind(c)
    import
    type(c_ptr), value :: h
    real(c_double), value :: x
    integer(c_int) :: gsl_histogram_increment
  end function gsl_histogram_increment
  function gsl_histogram_accumulate(h, x, weight) bind(c)
    import
    type(c_ptr), value :: h
    real(c_double), value :: x, weight
    integer(c_int) :: gsl_histogram_accumulate
  end function gsl_histogram_accumulate
  function gsl_histogram_get(h, i) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t), value :: i
    real(c_double) :: gsl_histogram_get
  end function gsl_histogram_get
  function gsl_histogram_get_range(h, i, lower, upper) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t), value :: i
    real(c_double), intent(out) :: lower, upper
    integer(c_int) :: gsl_histogram_get_range
  end function gsl_histogram_get_range
  function gsl_histogram_max(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram_max
  end function gsl_histogram_max
  function gsl_histogram_min(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram_min
  end function gsl_histogram_min
  function gsl_histogram_bins(h) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t) :: gsl_histogram_bins
  end function gsl_histogram_bins
  subroutine gsl_histogram_reset(h) bind(c)
    import
    type(c_ptr), value :: h
  end subroutine gsl_histogram_reset
  function gsl_histogram_find(h, x, i) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double), value :: x    
    integer(c_size_t), intent(out) :: i
    integer(c_int) :: gsl_histogram_find
  end function gsl_histogram_find
  function gsl_histogram_max_val(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram_max_val
  end function gsl_histogram_max_val
  function gsl_histogram_max_bin(h) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t) :: gsl_histogram_max_bin
  end function gsl_histogram_max_bin
  function gsl_histogram_min_val(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram_min_val
  end function gsl_histogram_min_val
  function gsl_histogram_min_bin(h) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t) :: gsl_histogram_min_bin
  end function gsl_histogram_min_bin
  function gsl_histogram_mean(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram_mean
  end function gsl_histogram_mean
  function gsl_histogram_sigma(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram_sigma
  end function gsl_histogram_sigma
  function gsl_histogram_sum(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram_sum
  end function gsl_histogram_sum
  function gsl_histogram_equal_bins_p(h1, h2) bind(c) 
    import
    type(c_ptr), value :: h1, h2
    integer(c_int) :: gsl_histogram_equal_bins_p
  end function gsl_histogram_equal_bins_p
  function gsl_histogram_add(h1, h2) bind(c) 
    import
    type(c_ptr), value :: h1, h2
    integer(c_int) :: gsl_histogram_add
  end function gsl_histogram_add
  function gsl_histogram_sub(h1, h2) bind(c) 
    import
    type(c_ptr), value :: h1, h2
    integer(c_int) :: gsl_histogram_sub
  end function gsl_histogram_sub
  function gsl_histogram_mul(h1, h2) bind(c) 
    import
    type(c_ptr), value :: h1, h2
    integer(c_int) :: gsl_histogram_mul
  end function gsl_histogram_mul
  function gsl_histogram_div(h1, h2) bind(c) 
    import
    type(c_ptr), value :: h1, h2
    integer(c_int) :: gsl_histogram_div
  end function gsl_histogram_div
  function gsl_histogram_scale(h, scale) bind(c)
    import
    type(c_ptr), value :: h
    real(c_double), value :: scale
    integer(c_int) :: gsl_histogram_scale
  end function gsl_histogram_scale
  function gsl_histogram_shift(h, offset) bind(c)
    import
    type(c_ptr), value :: h
    real(c_double), value :: offset
    integer(c_int) :: gsl_histogram_shift
  end function gsl_histogram_shift
  function gsl_histogram_fwrite(stream, h) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, h
    integer(c_int) :: gsl_histogram_fwrite
  end function gsl_histogram_fwrite
  function gsl_histogram_fread(stream, h) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, h
    integer(c_int) :: gsl_histogram_fread
  end function gsl_histogram_fread
  function gsl_histogram_fprintf(stream, h, range_format, bin_format) bind(c)
    import :: c_ptr, c_int, c_char
    type(c_ptr), value :: stream, h
    character(kind=c_char), dimension(*) :: range_format, bin_format
    integer(c_int) :: gsl_histogram_fprintf
  end function gsl_histogram_fprintf
  function gsl_histogram_fscanf(stream, h) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, h
    integer(c_int) :: gsl_histogram_fscanf
  end function gsl_histogram_fscanf
  function gsl_histogram_pdf_alloc(n) bind(c)
    import
    integer(c_size_t), value :: n
    type(c_ptr) :: gsl_histogram_pdf_alloc
  end function gsl_histogram_pdf_alloc
  function gsl_histogram_pdf_init(p, h) bind(c)
    import
    type(c_ptr), value :: p, h
    integer(c_int) :: gsl_histogram_pdf_init
  end function gsl_histogram_pdf_init
  subroutine gsl_histogram_pdf_free(p) bind(c)
    import
    type(c_ptr), value :: p
  end subroutine gsl_histogram_pdf_free
  function gsl_histogram_pdf_sample(p, r) bind(c)
    import
    type(c_ptr), value :: p
    real(c_double), value :: r
    real(c_double) :: gsl_histogram_pdf_sample
  end function gsl_histogram_pdf_sample
  function gsl_histogram2d_alloc(nx, ny) bind(c)
    import
    integer(c_size_t), value :: nx, ny
    type(c_ptr) :: gsl_histogram2d_alloc
  end function gsl_histogram2d_alloc
  function gsl_histogram2d_set_ranges(h, xrange, xsize, yrange, ysize) bind(c)
    import 
    type(c_ptr), value :: h
    integer(c_size_t), value :: xsize, ysize
    real(c_double), dimension(*), intent(in) :: xrange, yrange
    integer(c_int) :: gsl_histogram2d_set_ranges
  end function gsl_histogram2d_set_ranges
  function gsl_histogram2d_set_ranges_uniform(h, xmin, xmax, ymin, ymax) bind(c)
    import
    type(c_ptr), value :: h
    real(c_double), value :: xmin, xmax, ymin, ymax
    integer(c_int) :: gsl_histogram2d_set_ranges_uniform
  end function gsl_histogram2d_set_ranges_uniform
  subroutine gsl_histogram2d_free(h) bind(c)
    import
    type(c_ptr), value :: h
  end subroutine gsl_histogram2d_free
  function gsl_histogram2d_memcpy(dest, src) bind(c)
    import
    type(c_ptr), value :: dest, src
    integer(c_int) :: gsl_histogram2d_memcpy
  end function gsl_histogram2d_memcpy
  function gsl_histogram2d_clone(src) bind(c)
    import
    type(c_ptr), value :: src
    type(c_ptr) :: gsl_histogram2d_clone
  end function gsl_histogram2d_clone
  function gsl_histogram2d_increment(h, x, y) bind(c)
    import
    type(c_ptr), value :: h
    real(c_double), value :: x, y
    integer(c_int) :: gsl_histogram2d_increment
  end function gsl_histogram2d_increment
  function gsl_histogram2d_accumulate(h, x, y, weight) bind(c)
    import
    type(c_ptr), value :: h
    real(c_double), value :: x, y, weight
    integer(c_int) :: gsl_histogram2d_accumulate
  end function gsl_histogram2d_accumulate
  function gsl_histogram2d_get(h, i, j) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t), value :: i, j
    real(c_double) :: gsl_histogram2d_get
  end function gsl_histogram2d_get
  function gsl_histogram2d_get_xrange(h, i, xlower, xupper) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t), value :: i
    real(c_double), intent(out) :: xlower, xupper
    integer(c_int) :: gsl_histogram2d_get_xrange
  end function gsl_histogram2d_get_xrange
  function gsl_histogram2d_get_yrange(h, i, ylower, yupper) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t), value :: i
    real(c_double), intent(out) :: ylower, yupper
    integer(c_int) :: gsl_histogram2d_get_yrange
  end function gsl_histogram2d_get_yrange
  function gsl_histogram2d_xmax(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_xmax
  end function gsl_histogram2d_xmax
  function gsl_histogram2d_xmin(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_xmin
  end function gsl_histogram2d_xmin
  function gsl_histogram2d_nx(h) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t) :: gsl_histogram2d_nx
  end function gsl_histogram2d_nx
  function gsl_histogram2d_ymax(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_ymax
  end function gsl_histogram2d_ymax
  function gsl_histogram2d_ymin(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_ymin
  end function gsl_histogram2d_ymin
  function gsl_histogram2d_ny(h) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t) :: gsl_histogram2d_ny
  end function gsl_histogram2d_ny
  subroutine gsl_histogram2d_reset(h) bind(c)
    import
    type(c_ptr), value :: h
  end subroutine gsl_histogram2d_reset
  function gsl_histogram2d_find(h, x, y, i, j) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double), value :: x, y    
    integer(c_size_t), intent(out) :: i, j
    integer(c_int) :: gsl_histogram2d_find
  end function gsl_histogram2d_find
  function gsl_histogram2d_max_val(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_max_val
  end function gsl_histogram2d_max_val
  subroutine gsl_histogram2d_max_bin(h, i, j) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t), intent(out) :: i, j
  end subroutine gsl_histogram2d_max_bin
  function gsl_histogram2d_min_val(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_min_val
  end function gsl_histogram2d_min_val
  subroutine gsl_histogram2d_min_bin(h, i, j) bind(c) 
    import
    type(c_ptr), value :: h
    integer(c_size_t), intent(out) :: i, j
  end subroutine gsl_histogram2d_min_bin
  function gsl_histogram2d_xmean(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_xmean
  end function gsl_histogram2d_xmean
  function gsl_histogram2d_ymean(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_ymean
  end function gsl_histogram2d_ymean
  function gsl_histogram2d_xsigma(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_xsigma
  end function gsl_histogram2d_xsigma
  function gsl_histogram2d_ysigma(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_ysigma
  end function gsl_histogram2d_ysigma
  function gsl_histogram2d_cov(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_cov
  end function gsl_histogram2d_cov
  function gsl_histogram2d_sum(h) bind(c) 
    import
    type(c_ptr), value :: h
    real(c_double) :: gsl_histogram2d_sum
  end function gsl_histogram2d_sum
  function gsl_histogram2d_equal_bins_p(h1, h2) bind(c) 
    import
    type(c_ptr), value :: h1, h2
    integer(c_int) :: gsl_histogram2d_equal_bins_p
  end function gsl_histogram2d_equal_bins_p
  function gsl_histogram2d_add(h1, h2) bind(c) 
    import
    type(c_ptr), value :: h1, h2
    integer(c_int) :: gsl_histogram2d_add
  end function gsl_histogram2d_add
  function gsl_histogram2d_sub(h1, h2) bind(c) 
    import
    type(c_ptr), value :: h1, h2
    integer(c_int) :: gsl_histogram2d_sub
  end function gsl_histogram2d_sub
  function gsl_histogram2d_mul(h1, h2) bind(c) 
    import
    type(c_ptr), value :: h1, h2
    integer(c_int) :: gsl_histogram2d_mul
  end function gsl_histogram2d_mul
  function gsl_histogram2d_div(h1, h2) bind(c) 
    import
    type(c_ptr), value :: h1, h2
    integer(c_int) :: gsl_histogram2d_div
  end function gsl_histogram2d_div
  function gsl_histogram2d_scale(h, scale) bind(c)
    import
    type(c_ptr), value :: h
    real(c_double), value :: scale
    integer(c_int) :: gsl_histogram2d_scale
  end function gsl_histogram2d_scale
  function gsl_histogram2d_shift(h, offset) bind(c)
    import
    type(c_ptr), value :: h
    real(c_double), value :: offset
    integer(c_int) :: gsl_histogram2d_shift
  end function gsl_histogram2d_shift
  function gsl_histogram2d_fwrite(stream, h) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, h
    integer(c_int) :: gsl_histogram2d_fwrite
  end function gsl_histogram2d_fwrite
  function gsl_histogram2d_fread(stream, h) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, h
    integer(c_int) :: gsl_histogram2d_fread
  end function gsl_histogram2d_fread
  function gsl_histogram2d_fprintf(stream, h, range_format, bin_format) bind(c)
    import :: c_ptr, c_int, c_char
    type(c_ptr), value :: stream, h
    character(kind=c_char), dimension(*) :: range_format, bin_format
    integer(c_int) :: gsl_histogram2d_fprintf
  end function gsl_histogram2d_fprintf
  function gsl_histogram2d_fscanf(stream, h) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: stream, h
    integer(c_int) :: gsl_histogram2d_fscanf
  end function gsl_histogram2d_fscanf
  function gsl_histogram2d_pdf_alloc(nx, ny) bind(c)
    import
    integer(c_size_t), value :: nx, ny
    type(c_ptr) :: gsl_histogram2d_pdf_alloc
  end function gsl_histogram2d_pdf_alloc
  function gsl_histogram2d_pdf_init(p, h) bind(c)
    import
    type(c_ptr), value :: p, h
    integer(c_int) :: gsl_histogram2d_pdf_init
  end function gsl_histogram2d_pdf_init
  subroutine gsl_histogram2d_pdf_free(p) bind(c)
    import
    type(c_ptr), value :: p
  end subroutine gsl_histogram2d_pdf_free
  function gsl_histogram2d_pdf_sample(p, r1, r2, x, y) bind(c)
    import
    type(c_ptr), value :: p
    real(c_double), value :: r1, r2
    real(c_double), intent(out) :: x, y
    integer(c_int) :: gsl_histogram2d_pdf_sample
  end function gsl_histogram2d_pdf_sample
!-*-f90-*-
!
! Interfaces: Ntuples
! 
  function gsl_ntuple_create(fname, data, size) bind(c)
    import
    character(kind=c_char), dimension(*) :: fname
    type(c_ptr), value :: data
    integer(c_size_t), value :: size
    type(c_ptr) :: gsl_ntuple_create
  end function gsl_ntuple_create
  function gsl_ntuple_open(fname, data, size) bind(c)
    import
    character(kind=c_char), dimension(*) :: fname
    type(c_ptr), value :: data
    integer(c_size_t), value :: size
    type(c_ptr) :: gsl_ntuple_open
  end function gsl_ntuple_open
  function gsl_ntuple_write(ntuple) bind(c)
    import
    type(c_ptr), value :: ntuple
    integer(c_int) :: gsl_ntuple_write
  end function gsl_ntuple_write
  function gsl_ntuple_read(ntuple) bind(c)
    import
    type(c_ptr), value :: ntuple
    integer(c_int) :: gsl_ntuple_read
  end function gsl_ntuple_read
  function gsl_ntuple_close(ntuple) bind(c)
    import
    type(c_ptr), value :: ntuple
    integer(c_int) :: gsl_ntuple_close
  end function gsl_ntuple_close
  function fgsl_ntuple_select_fn_cinit(func, params) bind(c)
    import
    type(c_funptr), value :: func
    type(c_ptr), value :: params
    type(c_ptr) :: fgsl_ntuple_select_fn_cinit
  end function fgsl_ntuple_select_fn_cinit
  function fgsl_ntuple_value_fn_cinit(func, params) bind(c)
    import
    type(c_funptr), value :: func
    type(c_ptr), value :: params
    type(c_ptr) :: fgsl_ntuple_value_fn_cinit
  end function fgsl_ntuple_value_fn_cinit
  subroutine fgsl_ntuple_select_fn_cfree(sfunc) bind(c)
    import
    type(c_ptr), value :: sfunc
  end subroutine fgsl_ntuple_select_fn_cfree
  subroutine fgsl_ntuple_value_fn_cfree(sfunc) bind(c)
    import
    type(c_ptr), value :: sfunc
  end subroutine fgsl_ntuple_value_fn_cfree
  function gsl_ntuple_project(h, ntuple, value_func, select_func) bind(c)
    import
    type(c_ptr), value :: h, ntuple, value_func, select_func
    integer(c_int) :: gsl_ntuple_project
  end function gsl_ntuple_project
!
  function fgsl_aux_ntuple_data(ntuple) bind(c)
    import
    type(c_ptr), value :: ntuple
    type(c_ptr) :: fgsl_aux_ntuple_data
  end function fgsl_aux_ntuple_data
  function fgsl_aux_ntuple_size(ntuple) bind(c)
    import
    type(c_ptr), value :: ntuple
    integer(c_size_t) :: fgsl_aux_ntuple_size
  end function fgsl_aux_ntuple_size
!-*-f90-*-
!
! Interfaces: Monte Carlo integration
!
  function fgsl_monte_function_cinit(func, dim, params) bind(c)
    import :: c_funptr, c_ptr, c_size_t
    type(c_funptr), value :: func
    integer(c_size_t), value :: dim
    type(c_ptr), value :: params
    type(c_ptr) :: fgsl_monte_function_cinit
  end function fgsl_monte_function_cinit
  subroutine fgsl_monte_function_cfree(func) bind(c)
    import
    type(c_ptr), value :: func
  end subroutine fgsl_monte_function_cfree
  function gsl_monte_plain_alloc(dim) bind(c)
    import
    integer(c_size_t), value :: dim
    type(c_ptr) :: gsl_monte_plain_alloc
  end function gsl_monte_plain_alloc
  function gsl_monte_plain_init(s) bind(c)
    import
    type(c_ptr), value :: s
    integer(c_int) :: gsl_monte_plain_init
  end function gsl_monte_plain_init
  function gsl_monte_plain_integrate(f, xl, xu, dim, calls, r, s, result, abserr) bind(c)
    import
    type(c_ptr), value :: f
    integer(c_size_t), value :: dim, calls
    real(c_double), intent(in) :: xl(dim), xu(dim)
    type(c_ptr), value :: r, s
    real(c_double), intent(out) :: result, abserr
    integer(c_int) :: gsl_monte_plain_integrate
  end function gsl_monte_plain_integrate
  subroutine gsl_monte_plain_free(s) bind(c)
    import
    type(c_ptr), value :: s
  end subroutine gsl_monte_plain_free
  function gsl_monte_miser_alloc(dim) bind(c)
    import
    integer(c_size_t), value :: dim
    type(c_ptr) :: gsl_monte_miser_alloc 
  end function gsl_monte_miser_alloc
  function gsl_monte_miser_init(s) bind(c)
    import
    type(c_ptr), value :: s
    integer(c_int) :: gsl_monte_miser_init
  end function gsl_monte_miser_init
  function gsl_monte_miser_integrate(f, xl, xu, dim, calls, r, s, result, abserr) bind(c)
    import
    type(c_ptr), value :: f
    integer(c_size_t), value :: dim, calls
    real(c_double), intent(in) :: xl(dim), xu(dim)
    type(c_ptr), value :: r, s
    real(c_double), intent(out) :: result, abserr
    integer(c_int) :: gsl_monte_miser_integrate
  end function gsl_monte_miser_integrate
  subroutine gsl_monte_miser_free(s) bind(c)
    import
    type(c_ptr), value :: s
  end subroutine gsl_monte_miser_free
  function gsl_monte_vegas_alloc(dim) bind(c)
    import
    integer(c_size_t), value :: dim
    type(c_ptr) :: gsl_monte_vegas_alloc 
  end function gsl_monte_vegas_alloc
  function gsl_monte_vegas_init(s) bind(c)
    import
    type(c_ptr), value :: s
    integer(c_int) :: gsl_monte_vegas_init
  end function gsl_monte_vegas_init
  function gsl_monte_vegas_integrate(f, xl, xu, dim, calls, r, s, result, abserr) bind(c)
    import
    type(c_ptr), value :: f
    integer(c_size_t), value :: dim, calls
    real(c_double), intent(in) :: xl(dim), xu(dim)
    type(c_ptr), value :: r, s
    real(c_double), intent(out) :: result, abserr
    integer(c_int) :: gsl_monte_vegas_integrate
  end function gsl_monte_vegas_integrate
  subroutine gsl_monte_vegas_free(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
  end subroutine gsl_monte_vegas_free
  function gsl_monte_vegas_chisq(s) bind(c)
    import :: c_double, c_ptr 
    real(c_double) :: gsl_monte_vegas_chisq
    type(c_ptr), value :: s
  end function gsl_monte_vegas_chisq
  subroutine gsl_monte_vegas_runval(s, result, sigma) bind(c)
    import :: c_double, c_ptr 
    type(c_ptr), value :: s
    real(c_double) :: result, sigma
  end subroutine gsl_monte_vegas_runval
! add-on routines
  subroutine fgsl_monte_miser_csetparams(s, estimate_frac, min_calls, &
       min_calls_per_bisection, alpha, dither) bind(c)
    import :: c_size_t, c_double, c_ptr
    type(c_ptr), value :: s
    real(c_double), value :: estimate_frac, alpha, dither
    integer(c_size_t), value ::  min_calls, min_calls_per_bisection
  end subroutine fgsl_monte_miser_csetparams
  subroutine fgsl_monte_miser_cgetparams(s, estimate_frac, min_calls, &
       min_calls_per_bisection, alpha, dither) bind(c)
    import :: c_size_t, c_double, c_ptr
    type(c_ptr), value :: s
    real(c_double), intent(out) :: estimate_frac, alpha, dither
    integer(c_size_t), intent(out) ::  min_calls, min_calls_per_bisection
  end subroutine fgsl_monte_miser_cgetparams
  subroutine fgsl_monte_vegas_csetparams(s, result, sigma, chisq, alpha, &
       iterations, stage, mode, verbose, ostream) bind(c)
    import :: c_size_t, c_double, c_ptr, c_int
    type(c_ptr), value :: s, ostream
    real(c_double), value :: result, sigma, chisq, alpha
    integer(c_size_t), value ::  iterations
    integer(c_int), value :: stage, mode, verbose
  end subroutine fgsl_monte_vegas_csetparams
  subroutine fgsl_monte_vegas_cgetparams(s, result, sigma, chisq, alpha, &
       iterations, stage, mode, verbose, ostream) bind(c)
    import :: c_size_t, c_double, c_ptr, c_int
    type(c_ptr), value :: s, ostream
    real(c_double), intent(out) :: result, sigma, chisq, alpha
    integer(c_size_t), intent(out) ::  iterations
    integer(c_int), intent(out) :: stage, mode, verbose
  end subroutine fgsl_monte_vegas_cgetparams
!-*-f90-*-
!
!  Interfaces:  Simulated annealing
!
  subroutine gsl_siman_solve(rng, x0_p, ef, take_step, distance, &
       print_position, copy_func, copy_constructor, destructor, &
       element_size, params) bind(c)
    import :: c_size_t, c_ptr, c_funptr, gsl_siman_params_t
    type(c_funptr), value :: ef, take_step, distance, print_position, &
         copy_func, copy_constructor, destructor
    type(c_ptr), value :: rng, x0_p
    integer(c_size_t), value :: element_size
    type(gsl_siman_params_t), value :: params
  end subroutine gsl_siman_solve
!-*-f90-*-
!
!  Interfaces: Ordinary differential equations
!
     function gsl_odeiv_step_alloc(t, dim) bind(c)
       import
       type(c_ptr), value :: t
       integer(c_size_t), value :: dim
       type(c_ptr) :: gsl_odeiv_step_alloc
     end function gsl_odeiv_step_alloc
     function fgsl_aux_odeiv_step_alloc(step_type) bind(c)
       import
       integer(c_int), value :: step_type
       type(c_ptr) :: fgsl_aux_odeiv_step_alloc
     end function fgsl_aux_odeiv_step_alloc
     function fgsl_odeiv_system_cinit(func, dimension, params, jacobian) bind(c)
       import
       type(c_funptr), value :: func
       integer(c_size_t), value :: dimension
       type(c_ptr), value :: params
       type(c_funptr), value :: jacobian
       type(c_ptr) :: fgsl_odeiv_system_cinit
     end function fgsl_odeiv_system_cinit
     subroutine fgsl_odeiv_system_cfree(system) bind(c)
       import
       type(c_ptr), value :: system
     end subroutine fgsl_odeiv_system_cfree
     function gsl_odeiv_step_reset(s) bind(c)
       import 
       type(c_ptr), value :: s
       integer(c_int) :: gsl_odeiv_step_reset
     end function gsl_odeiv_step_reset
     subroutine gsl_odeiv_step_free(s) bind(c)
       import 
       type(c_ptr), value :: s
     end subroutine gsl_odeiv_step_free
     function gsl_odeiv_step_name (s) bind(c)
       import
       type(c_ptr), value :: s
       type(c_ptr) :: gsl_odeiv_step_name
     end function gsl_odeiv_step_name
     function gsl_odeiv_step_order(s) bind(c)
       import
       type(c_ptr), value :: s
       integer(c_int) :: gsl_odeiv_step_order
     end function gsl_odeiv_step_order
     function gsl_odeiv_step_apply(s, t, h, y, yerr, dydt_in, dydt_out, dydt) bind(c)
       import
       type(c_ptr), value :: s
       real(c_double), value :: t, h
       real(c_double), dimension(*), intent(inout) :: y, yerr, dydt_in, dydt_out
       type(c_ptr), value :: dydt
       integer(c_int) :: gsl_odeiv_step_apply
     end function gsl_odeiv_step_apply
     function gsl_odeiv_control_standard_new(eps_abs, eps_rel, a_y, a_dydt) bind(c)
       import
       real(c_double), value :: eps_abs, eps_rel, a_y, a_dydt
       type(c_ptr) :: gsl_odeiv_control_standard_new
     end function gsl_odeiv_control_standard_new
     function gsl_odeiv_control_y_new(eps_abs, eps_rel) bind(c)
       import
       real(c_double), value :: eps_abs, eps_rel
       type(c_ptr) :: gsl_odeiv_control_y_new
     end function gsl_odeiv_control_y_new
     function gsl_odeiv_control_yp_new(eps_abs, eps_rel) bind(c)
       import
       real(c_double), value :: eps_abs, eps_rel
       type(c_ptr) :: gsl_odeiv_control_yp_new
     end function gsl_odeiv_control_yp_new
     function gsl_odeiv_control_scaled_new(eps_abs, eps_rel, a_y, a_dydt, &
          scale_abs, dim) bind(c)
       import
       real(c_double), value :: eps_abs, eps_rel, a_y, a_dydt
       real(c_double), dimension(*), intent(in) :: scale_abs
       integer(c_size_t), value :: dim
       type(c_ptr) :: gsl_odeiv_control_scaled_new
     end function gsl_odeiv_control_scaled_new
! gsl_odeiv_control_alloc presently not attached
     function gsl_odeiv_control_init(c, eps_abs, eps_rel, a_y, a_dydt) bind(c)
       import
       type(c_ptr), value :: c
       real(c_double), value :: eps_abs, eps_rel, a_y, a_dydt
       integer(c_int) :: gsl_odeiv_control_init
     end function gsl_odeiv_control_init
     subroutine gsl_odeiv_control_free(c) bind(c)
       import
       type(c_ptr), value :: c
     end subroutine gsl_odeiv_control_free
     function gsl_odeiv_control_hadjust(c, s, y0, yerr, dydt, h) bind(c)
       import
       type(c_ptr), value :: c, s
       real(c_double), dimension(*), intent(in) :: y0, yerr, dydt
       real(c_double), dimension(*), intent(inout) :: h
       integer(c_int) :: gsl_odeiv_control_hadjust
     end function gsl_odeiv_control_hadjust
     function gsl_odeiv_control_name (s) bind(c)
       import
       type(c_ptr), value :: s
       type(c_ptr) :: gsl_odeiv_control_name
     end function gsl_odeiv_control_name
     function gsl_odeiv_evolve_alloc(dim) bind(c)
       import
       integer(c_size_t), value :: dim
       type(c_ptr) gsl_odeiv_evolve_alloc
     end function gsl_odeiv_evolve_alloc
     function gsl_odeiv_evolve_apply(e, con, step, dydt, t, t1, h, y) bind(c)
       import
       type(c_ptr), value :: e, con, step, dydt
       real(c_double), dimension(*), intent(inout) :: y
       real(c_double), intent(inout) :: t, h
       real(c_double), value :: t1
       integer(c_int) :: gsl_odeiv_evolve_apply
     end function gsl_odeiv_evolve_apply
     function gsl_odeiv_evolve_reset(s) bind(c)
       import 
       type(c_ptr), value :: s
       integer(c_int) :: gsl_odeiv_evolve_reset
     end function gsl_odeiv_evolve_reset
     subroutine gsl_odeiv_evolve_free(s) bind(c)
       import 
       type(c_ptr), value :: s
     end subroutine gsl_odeiv_evolve_free
!-*-f90-*-
!
! Interfaces: Numerical Differentiation
!
  function gsl_deriv_central(f, x, h, result, abserr) bind(c)
    import
    type(c_ptr), value :: f
    real(c_double), value :: x, h
    real(c_double), intent(out) :: result, abserr
    integer(c_int) :: gsl_deriv_central
  end function gsl_deriv_central
  function gsl_deriv_forward(f, x, h, result, abserr) bind(c)
    import
    type(c_ptr), value :: f
    real(c_double), value :: x, h
    real(c_double), intent(out) :: result, abserr
    integer(c_int) :: gsl_deriv_forward
  end function gsl_deriv_forward
  function gsl_deriv_backward(f, x, h, result, abserr) bind(c)
    import
    type(c_ptr), value :: f
    real(c_double), value :: x, h
    real(c_double), intent(out) :: result, abserr
    integer(c_int) :: gsl_deriv_backward
  end function gsl_deriv_backward
!-*-f90-*-
!
! Interfaces: Chebyshev Approximations
! 
  function gsl_cheb_alloc(n) bind(c)
    import
    integer(c_int), value :: n
    type(c_ptr) :: gsl_cheb_alloc
  end function gsl_cheb_alloc
  subroutine gsl_cheb_free(cs) bind(c)
    import
    type(c_ptr), value :: cs
  end subroutine gsl_cheb_free
  function gsl_cheb_init(cs, f, a, b) bind(c)
    import
    type(c_ptr), value :: cs, f
    real(c_double), value :: a, b
    integer(c_int) :: gsl_cheb_init
  end function gsl_cheb_init
  function gsl_cheb_order(cs) bind(c)
    import :: c_ptr, c_size_t
    integer(c_size_t) :: gsl_cheb_order
    type(c_ptr), value :: cs
  end function gsl_cheb_order
  function gsl_cheb_size(cs) bind(c)
    import :: c_ptr, c_size_t
    integer(c_size_t) :: gsl_cheb_size
    type(c_ptr), value :: cs
  end function gsl_cheb_size
  function gsl_cheb_coeffs(cs) bind(c)
    import :: c_ptr
    type(c_ptr) :: gsl_cheb_coeffs
    type(c_ptr), value :: cs
  end function gsl_cheb_coeffs
  function gsl_cheb_eval(cs, x) bind(c)
    import
    type(c_ptr), value :: cs
    real(c_double), value :: x
    real(c_double) :: gsl_cheb_eval
  end function gsl_cheb_eval
  function gsl_cheb_eval_err(cs, x, result, abserr) bind(c)
    import
    type(c_ptr), value :: cs
    real(c_double), value :: x
    real(c_double), intent(out) :: result, abserr
    integer(c_int) :: gsl_cheb_eval_err
  end function gsl_cheb_eval_err
  function gsl_cheb_eval_n(cs, order, x) bind(c)
    import
    type(c_ptr), value :: cs
    integer(c_size_t), value :: order
    real(c_double), value :: x
    real(c_double) :: gsl_cheb_eval_n
  end function gsl_cheb_eval_n
  function gsl_cheb_eval_n_err(cs, order, x, result, abserr) bind(c)
    import
    type(c_ptr), value :: cs
    integer(c_size_t), value :: order
    real(c_double), value :: x
    real(c_double), intent(out) :: result, abserr
    integer(c_int) :: gsl_cheb_eval_n_err
  end function gsl_cheb_eval_n_err
  function gsl_cheb_calc_deriv(deriv, cs) bind(c)
    import
    type(c_ptr), value :: deriv, cs
    integer(c_int) :: gsl_cheb_calc_deriv
  end function gsl_cheb_calc_deriv
  function gsl_cheb_calc_integ(integ, cs) bind(c)
    import
    type(c_ptr), value :: integ, cs
    integer(c_int) :: gsl_cheb_calc_integ
  end function gsl_cheb_calc_integ
!-*-f90-*-
!
! Interfaces: series acceleration
! 
  function gsl_sum_levin_u_alloc(n) bind(c)
    import
    integer(c_size_t), value :: n
    type(c_ptr) :: gsl_sum_levin_u_alloc
  end function gsl_sum_levin_u_alloc
  function gsl_sum_levin_u_free(w) bind(c)
    import
    type(c_ptr), value :: w
    integer(c_int) :: gsl_sum_levin_u_free
  end function gsl_sum_levin_u_free
  function gsl_sum_levin_u_accel(array, array_size, w, sum_accel, abserr) bind(c)
    import
    integer(c_size_t), value :: array_size
    real(c_double), intent(in) :: array(array_size)
    type(c_ptr), value :: w
    real(c_double), intent(out) :: sum_accel, abserr
    integer(c_int) :: gsl_sum_levin_u_accel
  end function gsl_sum_levin_u_accel
  function gsl_sum_levin_utrunc_alloc(n) bind(c)
    import
    integer(c_size_t), value :: n
    type(c_ptr) :: gsl_sum_levin_utrunc_alloc
  end function gsl_sum_levin_utrunc_alloc
  function gsl_sum_levin_utrunc_free(w) bind(c)
    import
    type(c_ptr), value :: w
    integer(c_int) :: gsl_sum_levin_utrunc_free
  end function gsl_sum_levin_utrunc_free
  function gsl_sum_levin_utrunc_accel(array, array_size, w, sum_accel, abserr) bind(c)
    import
    integer(c_size_t), value :: array_size
    real(c_double), intent(in) :: array(array_size)
    type(c_ptr), value :: w
    real(c_double), intent(out) :: sum_accel, abserr
    integer(c_int) :: gsl_sum_levin_utrunc_accel
  end function gsl_sum_levin_utrunc_accel
  
!-*-f90-*-
!
! Interfaces: wavelet transforms
! 
  function gsl_wavelet_alloc(t, k) bind(c)
    import
    type(c_ptr), value :: t
    integer(c_size_t), value :: k
    type(c_ptr) :: gsl_wavelet_alloc
  end function gsl_wavelet_alloc
  function gsl_wavelet_name(wavelet) bind(c)
    import
    type(c_ptr), value :: wavelet
    type(c_ptr) :: gsl_wavelet_name
  end function gsl_wavelet_name
  subroutine gsl_wavelet_free(w) bind(c)
    import
    type(c_ptr), value :: w
  end subroutine gsl_wavelet_free
  function gsl_wavelet_workspace_alloc(n) bind(c)
    import
    integer(c_size_t), value :: n
    type(c_ptr) :: gsl_wavelet_workspace_alloc
  end function gsl_wavelet_workspace_alloc
  subroutine gsl_wavelet_workspace_free(w) bind(c)
    import
    type(c_ptr), value :: w
  end subroutine gsl_wavelet_workspace_free
  function gsl_wavelet_transform(w, data, stride, n, dir, work) bind(c)
    import
    type(c_ptr), value :: w
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_size_t), value :: stride, n
    integer(c_int), value :: dir
    type(c_ptr), value :: work
    integer(c_int) :: gsl_wavelet_transform
  end function gsl_wavelet_transform
  function gsl_wavelet_transform_forward(w, data, stride, n, work) bind(c)
    import
    type(c_ptr), value :: w
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_size_t), value :: stride, n
    type(c_ptr), value :: work
    integer(c_int) :: gsl_wavelet_transform_forward
  end function gsl_wavelet_transform_forward
  function gsl_wavelet_transform_inverse(w, data, stride, n, work) bind(c)
    import
    type(c_ptr), value :: w
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_size_t), value :: stride, n
    type(c_ptr), value :: work
    integer(c_int) :: gsl_wavelet_transform_inverse
  end function gsl_wavelet_transform_inverse
  function gsl_wavelet2d_transform(w, data, tda, size1, size2, dir, work) bind(c)
    import
    type(c_ptr), value :: w
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_size_t), value :: tda, size1, size2
    integer(c_int), value :: dir
    type(c_ptr), value :: work
    integer(c_int) :: gsl_wavelet2d_transform
  end function gsl_wavelet2d_transform
  function gsl_wavelet2d_transform_forward(w, data, tda, size1, size2, work) bind(c)
    import
    type(c_ptr), value :: w
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_size_t), value :: tda, size1, size2
    type(c_ptr), value :: work
    integer(c_int) :: gsl_wavelet2d_transform_forward
  end function gsl_wavelet2d_transform_forward
  function gsl_wavelet2d_transform_inverse(w, data, tda, size1, size2, work) bind(c)
    import
    type(c_ptr), value :: w
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_size_t), value :: tda, size1, size2
    type(c_ptr), value :: work
    integer(c_int) :: gsl_wavelet2d_transform_inverse
  end function gsl_wavelet2d_transform_inverse
  function gsl_wavelet2d_transform_matrix(w, m, dir, work) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: w, m, work
    integer(c_int), value :: dir
    integer(c_int) :: gsl_wavelet2d_transform_matrix
  end function gsl_wavelet2d_transform_matrix
  function gsl_wavelet2d_transform_matrix_forward(w, m, work) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: w, m, work
    integer(c_int) :: gsl_wavelet2d_transform_matrix_forward
  end function gsl_wavelet2d_transform_matrix_forward
  function gsl_wavelet2d_transform_matrix_inverse(w, m, work) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: w, m, work
    integer(c_int) :: gsl_wavelet2d_transform_matrix_inverse
  end function gsl_wavelet2d_transform_matrix_inverse
  function gsl_wavelet2d_nstransform(w, data, tda, size1, size2, dir, work) bind(c)
    import
    type(c_ptr), value :: w
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_size_t), value :: tda, size1, size2
    integer(c_int), value :: dir
    type(c_ptr), value :: work
    integer(c_int) :: gsl_wavelet2d_nstransform
  end function gsl_wavelet2d_nstransform
  function gsl_wavelet2d_nstransform_forward(w, data, tda, size1, size2, work) bind(c)
    import
    type(c_ptr), value :: w
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_size_t), value :: tda, size1, size2
    type(c_ptr), value :: work
    integer(c_int) :: gsl_wavelet2d_nstransform_forward
  end function gsl_wavelet2d_nstransform_forward
  function gsl_wavelet2d_nstransform_inverse(w, data, tda, size1, size2, work) bind(c)
    import
    type(c_ptr), value :: w
    real(c_double), dimension(*), intent(inout) :: data
    integer(c_size_t), value :: tda, size1, size2
    type(c_ptr), value :: work
    integer(c_int) :: gsl_wavelet2d_nstransform_inverse
  end function gsl_wavelet2d_nstransform_inverse
  function gsl_wavelet2d_nstransform_matrix(w, m, dir, work) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: w, m, work
    integer(c_int), value :: dir
    integer(c_int) :: gsl_wavelet2d_nstransform_matrix
  end function gsl_wavelet2d_nstransform_matrix
  function gsl_wavelet2d_nstransform_matrix_forward(w, m, work) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: w, m, work
    integer(c_int) :: gsl_wavelet2d_nstransform_matrix_forward
  end function gsl_wavelet2d_nstransform_matrix_forward
  function gsl_wavelet2d_nstransform_matrix_inverse(w, m, work) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: w, m, work
    integer(c_int) :: gsl_wavelet2d_nstransform_matrix_inverse
  end function gsl_wavelet2d_nstransform_matrix_inverse
!
  function fgsl_aux_wavelet_alloc(i) bind(c)
    import
    integer(c_int), value :: i
    type(c_ptr) :: fgsl_aux_wavelet_alloc
  end function fgsl_aux_wavelet_alloc
  function gsl_aux_sizeof_wavelet_workspace() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_wavelet_workspace
  end function gsl_aux_sizeof_wavelet_workspace
  function gsl_aux_sizeof_wavelet() bind(c)
    import :: c_size_t
    integer(c_size_t) :: gsl_aux_sizeof_wavelet
  end function gsl_aux_sizeof_wavelet
    
!-*-f90-*-
!
! Interfaces: Hankel transforms
! 
  function gsl_dht_alloc(size) bind(c)
    import
    integer(c_size_t), value :: size
    type(c_ptr) :: gsl_dht_alloc
  end function gsl_dht_alloc
  function gsl_dht_init(t, nu, xmax) bind(c)
    import
    type(c_ptr), value :: t
    real(c_double), value :: nu, xmax
    integer(c_int) :: gsl_dht_init
  end function gsl_dht_init
  function gsl_dht_new(size, nu, xmax) bind(c)
    import
    integer(c_size_t), value :: size
    real(c_double), value :: nu, xmax
    type(c_ptr) :: gsl_dht_new
  end function gsl_dht_new
  subroutine gsl_dht_free(t) bind(c)
    import
    type(c_ptr), value :: t
  end subroutine gsl_dht_free
  function gsl_dht_apply(t, f_in, f_out) bind(c)
    import
    type(c_ptr), value :: t
    real(c_double), dimension(*), intent(in) :: f_in
    real(c_double), dimension(*), intent(out) :: f_out 
    integer(c_int) :: gsl_dht_apply
  end function gsl_dht_apply
  function gsl_dht_x_sample(t, n) bind(c)
    import
    type(c_ptr), value :: t
    integer(c_int), value :: n
    real(c_double) :: gsl_dht_x_sample
  end function gsl_dht_x_sample
  function gsl_dht_k_sample(t, n) bind(c)
    import
    type(c_ptr), value :: t
    integer(c_int), value :: n
    real(c_double) :: gsl_dht_k_sample
  end function gsl_dht_k_sample
    
    
!-*-f90-*-
!
!  Interfaces: Root finding
!
  function gsl_root_fsolver_alloc(t) bind(c)
    import
    type(c_ptr), value :: t
    type(c_ptr) :: gsl_root_fsolver_alloc
  end function gsl_root_fsolver_alloc
  function fgsl_aux_fsolver_alloc(it) bind(c)
    import
    integer(c_int), value :: it
    type(c_ptr) :: fgsl_aux_fsolver_alloc
  end function fgsl_aux_fsolver_alloc
  function gsl_root_fdfsolver_alloc(t) bind(c)
    import
    type(c_ptr), value :: t
    type(c_ptr) :: gsl_root_fdfsolver_alloc
  end function gsl_root_fdfsolver_alloc
  function fgsl_aux_fdfsolver_alloc(it) bind(c)
    import
    integer(c_int), value :: it
    type(c_ptr) :: fgsl_aux_fdfsolver_alloc
  end function fgsl_aux_fdfsolver_alloc
  function gsl_root_fsolver_set(s, f, x_lower, x_upper) bind(c)
    import
    type(c_ptr), value :: s, f
    real(c_double), value :: x_lower, x_upper
    integer(c_int) :: gsl_root_fsolver_set
  end function gsl_root_fsolver_set
  function gsl_root_fdfsolver_set(s, f, x) bind(c)
    import
    type(c_ptr), value :: s, f
    real(c_double), value :: x
    integer(c_int) :: gsl_root_fdfsolver_set
  end function gsl_root_fdfsolver_set
  subroutine gsl_root_fsolver_free(s) bind(c)
    import
    type(c_ptr), value :: s
  end subroutine gsl_root_fsolver_free
  subroutine gsl_root_fdfsolver_free(s) bind(c)
    import
    type(c_ptr), value :: s
  end subroutine gsl_root_fdfsolver_free
  function gsl_root_fsolver_name(s) bind(c)
    import
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_root_fsolver_name
  end function gsl_root_fsolver_name
  function gsl_root_fdfsolver_name(s) bind(c)
    import
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_root_fdfsolver_name
  end function gsl_root_fdfsolver_name
  function gsl_root_fsolver_iterate(s) bind(c)
    import
    type(c_ptr), value :: s
    integer(c_int) :: gsl_root_fsolver_iterate
  end function gsl_root_fsolver_iterate
  function gsl_root_fdfsolver_iterate(s) bind(c)
    import
    type(c_ptr), value :: s
    integer(c_int) :: gsl_root_fdfsolver_iterate
  end function gsl_root_fdfsolver_iterate
  function gsl_root_fsolver_root(s) bind(c)
    import
    type(c_ptr), value :: s
    real(c_double) :: gsl_root_fsolver_root
  end function gsl_root_fsolver_root
  function gsl_root_fdfsolver_root(s) bind(c)
    import
    type(c_ptr), value :: s
    real(c_double) :: gsl_root_fdfsolver_root
  end function gsl_root_fdfsolver_root
  function gsl_root_fsolver_x_lower(s) bind(c)
    import
    type(c_ptr), value :: s
    real(c_double) :: gsl_root_fsolver_x_lower
  end function gsl_root_fsolver_x_lower
  function gsl_root_fsolver_x_upper(s) bind(c)
    import
    type(c_ptr), value :: s
    real(c_double) :: gsl_root_fsolver_x_upper
  end function gsl_root_fsolver_x_upper
  function gsl_root_test_interval(x_lower, x_upper, epsabs, epsrel) bind(c)
    import
    real(c_double), value :: x_lower, x_upper, epsabs, epsrel
    integer(c_int) :: gsl_root_test_interval
  end function gsl_root_test_interval
  function gsl_root_test_delta(x1, x0, epsabs, epsrel) bind(c)
    import
    real(c_double), value :: x1, x0, epsabs, epsrel
    integer(c_int) :: gsl_root_test_delta
  end function gsl_root_test_delta
  function gsl_root_test_residual(f, epsabs) bind(c)
    import
    real(c_double), value :: f, epsabs
    integer(c_int) :: gsl_root_test_residual
  end function gsl_root_test_residual
!-*-f90-*-
!
! Interfaces: Minimization
!
  function gsl_min_fminimizer_alloc(t) bind(c)
    import
    type(c_ptr), value :: t
    type(c_ptr) :: gsl_min_fminimizer_alloc
  end function gsl_min_fminimizer_alloc
  function fgsl_aux_fminimizer_alloc(it) bind(c)
    import
    integer(c_int), value :: it
    type(c_ptr) :: fgsl_aux_fminimizer_alloc
  end function fgsl_aux_fminimizer_alloc
  subroutine gsl_min_fminimizer_free(s) bind(c)
    import
    type(c_ptr), value :: s
  end subroutine gsl_min_fminimizer_free
  function gsl_min_fminimizer_set(s, f, x_minimum, x_lower, x_upper) bind(c)
    import
    type(c_ptr), value :: s, f
    real(c_double), value :: x_minimum, x_lower, x_upper
    integer(c_int) :: gsl_min_fminimizer_set
  end function gsl_min_fminimizer_set
  function gsl_min_fminimizer_set_with_values(s, f, x_minimum, f_minimum, &
           x_lower, f_lower, x_upper, f_upper) bind(c)
    import
    type(c_ptr), value :: s, f
    real(c_double), value :: x_minimum, f_minimum, &
           x_lower, f_lower, x_upper, f_upper
    integer(c_int) :: gsl_min_fminimizer_set_with_values
  end function gsl_min_fminimizer_set_with_values
  function  gsl_min_fminimizer_iterate(s) bind(c)
    import 
    type(c_ptr), value :: s
    integer(c_int) :: gsl_min_fminimizer_iterate
  end function gsl_min_fminimizer_iterate
  function gsl_min_fminimizer_name(s) bind(c)
    import
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_min_fminimizer_name
  end function gsl_min_fminimizer_name
  function gsl_min_fminimizer_x_minimum(s) bind(c)
    import
    type(c_ptr), value :: s
    real(c_double) ::  gsl_min_fminimizer_x_minimum
  end function gsl_min_fminimizer_x_minimum
  function gsl_min_fminimizer_x_lower(s) bind(c)
    import
    type(c_ptr), value :: s
    real(c_double) ::  gsl_min_fminimizer_x_lower
  end function gsl_min_fminimizer_x_lower
  function gsl_min_fminimizer_x_upper(s) bind(c)
    import
    type(c_ptr), value :: s
    real(c_double) ::  gsl_min_fminimizer_x_upper
  end function gsl_min_fminimizer_x_upper
  function gsl_min_fminimizer_f_minimum(s) bind(c)
    import
    type(c_ptr), value :: s
    real(c_double) ::  gsl_min_fminimizer_f_minimum
  end function gsl_min_fminimizer_f_minimum
  function gsl_min_fminimizer_f_lower(s) bind(c)
    import
    type(c_ptr), value :: s
    real(c_double) ::  gsl_min_fminimizer_f_lower
  end function gsl_min_fminimizer_f_lower
  function gsl_min_fminimizer_f_upper(s) bind(c)
    import
    type(c_ptr), value :: s
    real(c_double) ::  gsl_min_fminimizer_f_upper
  end function gsl_min_fminimizer_f_upper
  function gsl_min_test_interval(x_lower, x_upper, epsabs, epsrel) bind(c)
    import
    real(c_double), value :: x_lower, x_upper, epsabs, epsrel
    integer(c_int) :: gsl_min_test_interval
  end function gsl_min_test_interval
!-*-f90-*-
!
!  Interfaces: multi-dimensional root finding
!
  function gsl_multiroot_fsolver_alloc(t, n) bind(c)
    import
    type(c_ptr), value :: t
    integer(c_size_t), value :: n
    type(c_ptr) :: gsl_multiroot_fsolver_alloc
  end function gsl_multiroot_fsolver_alloc
  function gsl_multiroot_fdfsolver_alloc(t, n) bind(c)
    import
    type(c_ptr), value :: t
    integer(c_size_t), value :: n
    type(c_ptr) :: gsl_multiroot_fdfsolver_alloc
  end function gsl_multiroot_fdfsolver_alloc
  function gsl_multiroot_fsolver_set(s, f, x) bind(c)
    import
    type(c_ptr), value :: s, f, x
    integer(c_int) :: gsl_multiroot_fsolver_set
  end function gsl_multiroot_fsolver_set
  function gsl_multiroot_fdfsolver_set(s, f, x) bind(c)
    import
    type(c_ptr), value :: s, f, x
    integer(c_int) :: gsl_multiroot_fdfsolver_set
  end function gsl_multiroot_fdfsolver_set
  subroutine gsl_multiroot_fsolver_free(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
  end subroutine gsl_multiroot_fsolver_free
  subroutine gsl_multiroot_fdfsolver_free(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
  end subroutine gsl_multiroot_fdfsolver_free
  function gsl_multiroot_fsolver_name(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multiroot_fsolver_name
  end function gsl_multiroot_fsolver_name
  function gsl_multiroot_fdfsolver_name(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multiroot_fdfsolver_name
  end function gsl_multiroot_fdfsolver_name
  function gsl_multiroot_fsolver_iterate(s) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: s
    integer(c_int) :: gsl_multiroot_fsolver_iterate
  end function gsl_multiroot_fsolver_iterate
  function gsl_multiroot_fdfsolver_iterate(s) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: s
    integer(c_int) :: gsl_multiroot_fdfsolver_iterate
  end function gsl_multiroot_fdfsolver_iterate
  function gsl_multiroot_fsolver_root(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multiroot_fsolver_root
  end function gsl_multiroot_fsolver_root
  function gsl_multiroot_fsolver_f(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multiroot_fsolver_f
  end function gsl_multiroot_fsolver_f
  function gsl_multiroot_fsolver_dx(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multiroot_fsolver_dx
  end function gsl_multiroot_fsolver_dx
  function gsl_multiroot_fdfsolver_root(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multiroot_fdfsolver_root
  end function gsl_multiroot_fdfsolver_root
  function gsl_multiroot_fdfsolver_f(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multiroot_fdfsolver_f
  end function gsl_multiroot_fdfsolver_f
  function gsl_multiroot_fdfsolver_dx(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multiroot_fdfsolver_dx
  end function gsl_multiroot_fdfsolver_dx
  function gsl_multiroot_test_delta(dx, x, epsabs, epsrel) bind(c)
    import :: c_ptr, c_int, c_double
    type(c_ptr), value :: dx, x
    real(c_double), value :: epsabs, epsrel
    integer(c_int) :: gsl_multiroot_test_delta
  end function gsl_multiroot_test_delta
!    int gsl_multiroot_test_residual (const gsl_vector * f, double epsabs)
  function gsl_multiroot_test_residual(f, epsabs) bind(c)
    import :: c_ptr, c_int, c_double
    type(c_ptr), value :: f
    real(c_double), value :: epsabs
    integer(c_int) :: gsl_multiroot_test_residual
  end function gsl_multiroot_test_residual
!
  function fgsl_multiroot_function_cinit(fp, ndim, params) bind(c)
    import
    type(c_funptr), value :: fp
    integer(c_size_t), value :: ndim
    type(c_ptr), value :: params
    type(c_ptr) :: fgsl_multiroot_function_cinit
  end function fgsl_multiroot_function_cinit
  function fgsl_multiroot_function_fdf_cinit(fp, dfp, fdfp, ndim, params) bind(c)
    import
    type(c_funptr), value :: fp, dfp, fdfp
    integer(c_size_t), value :: ndim
    type(c_ptr), value :: params
    type(c_ptr) :: fgsl_multiroot_function_fdf_cinit
  end function fgsl_multiroot_function_fdf_cinit
  subroutine fgsl_multiroot_function_cfree(f) bind(c)
    import :: c_ptr
    type(c_ptr), value :: f
  end subroutine fgsl_multiroot_function_cfree
  subroutine fgsl_multiroot_function_fdf_cfree(f) bind(c)
    import :: c_ptr
    type(c_ptr), value :: f
  end subroutine fgsl_multiroot_function_fdf_cfree
  function fgsl_aux_multiroot_fsolver_alloc(it) bind(c)
    import
    integer(c_int), value :: it
    type(c_ptr) :: fgsl_aux_multiroot_fsolver_alloc
  end function fgsl_aux_multiroot_fsolver_alloc
  function fgsl_aux_multiroot_fdfsolver_alloc(it) bind(c)
    import
    integer(c_int), value :: it
    type(c_ptr) :: fgsl_aux_multiroot_fdfsolver_alloc
  end function fgsl_aux_multiroot_fdfsolver_alloc
!-*-f90-*-
!
!  Interfaces: multi-dimensional minimization
!
  function gsl_multimin_fminimizer_alloc(t, n) bind(c)
    import
    type(c_ptr), value :: t
    integer(c_size_t), value :: n
    type(c_ptr) :: gsl_multimin_fminimizer_alloc
  end function gsl_multimin_fminimizer_alloc
  function gsl_multimin_fdfminimizer_alloc(t, n) bind(c)
    import
    type(c_ptr), value :: t
    integer(c_size_t), value :: n
    type(c_ptr) :: gsl_multimin_fdfminimizer_alloc
  end function gsl_multimin_fdfminimizer_alloc
  function gsl_multimin_fminimizer_set(s, f, x, step) bind(c)
    import
    type(c_ptr), value :: s, f, x, step
    integer(c_int) :: gsl_multimin_fminimizer_set
  end function gsl_multimin_fminimizer_set
  function gsl_multimin_fdfminimizer_set(s, f, x, step, tol) bind(c)
    import
    type(c_ptr), value :: s, f, x
    real(c_double), value :: step, tol
    integer(c_int) :: gsl_multimin_fdfminimizer_set
  end function gsl_multimin_fdfminimizer_set
  subroutine gsl_multimin_fminimizer_free(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
  end subroutine gsl_multimin_fminimizer_free
  subroutine gsl_multimin_fdfminimizer_free(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
  end subroutine gsl_multimin_fdfminimizer_free
  function gsl_multimin_fminimizer_name(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multimin_fminimizer_name
  end function gsl_multimin_fminimizer_name
  function gsl_multimin_fdfminimizer_name(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multimin_fdfminimizer_name
  end function gsl_multimin_fdfminimizer_name
  function gsl_multimin_fminimizer_iterate(s) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: s
    integer(c_int) :: gsl_multimin_fminimizer_iterate
  end function gsl_multimin_fminimizer_iterate
  function gsl_multimin_fdfminimizer_iterate(s) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: s
    integer(c_int) :: gsl_multimin_fdfminimizer_iterate
  end function gsl_multimin_fdfminimizer_iterate
  function gsl_multimin_fminimizer_x(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multimin_fminimizer_x
  end function gsl_multimin_fminimizer_x
  function gsl_multimin_fdfminimizer_x(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multimin_fdfminimizer_x
  end function gsl_multimin_fdfminimizer_x
  function gsl_multimin_fminimizer_minimum(s) bind(c)
    import :: c_ptr, c_double
    type(c_ptr), value :: s
    real(c_double) :: gsl_multimin_fminimizer_minimum
  end function gsl_multimin_fminimizer_minimum
  function gsl_multimin_fdfminimizer_minimum(s) bind(c)
    import :: c_ptr, c_double
    type(c_ptr), value :: s
    real(c_double) :: gsl_multimin_fdfminimizer_minimum
  end function gsl_multimin_fdfminimizer_minimum
  function gsl_multimin_fdfminimizer_gradient(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multimin_fdfminimizer_gradient
  end function gsl_multimin_fdfminimizer_gradient
  function gsl_multimin_fminimizer_size(s) bind(c)
    import :: c_ptr, c_double
    type(c_ptr), value :: s
    real(c_double) :: gsl_multimin_fminimizer_size
  end function gsl_multimin_fminimizer_size
  function gsl_multimin_fdfminimizer_restart(s) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: s
    integer(c_int) :: gsl_multimin_fdfminimizer_restart
  end function gsl_multimin_fdfminimizer_restart
  function gsl_multimin_test_gradient(g, epsabs) bind(c)
    import :: c_ptr, c_int, c_double
    type(c_ptr), value :: g
    real(c_double), value :: epsabs
    integer(c_int) :: gsl_multimin_test_gradient
  end function gsl_multimin_test_gradient
  function gsl_multimin_test_size(size, epsabs) bind(c)
    import :: c_int, c_double
    real(c_double), value :: size, epsabs
    integer(c_int) :: gsl_multimin_test_size
  end function gsl_multimin_test_size
!
  function fgsl_multimin_function_cinit(fp, ndim, params) bind(c)
    import 
    type(c_funptr), value :: fp
    integer(c_size_t), value :: ndim
    type(c_ptr), value :: params
    type(c_ptr) :: fgsl_multimin_function_cinit
  end function fgsl_multimin_function_cinit
  function fgsl_multimin_function_fdf_cinit(fp, dfp, fdfp, ndim, params) bind(c)
    import
    type(c_funptr), value :: fp, dfp, fdfp
    integer(c_size_t), value :: ndim
    type(c_ptr), value :: params
    type(c_ptr) :: fgsl_multimin_function_fdf_cinit
  end function fgsl_multimin_function_fdf_cinit
  subroutine fgsl_multimin_function_cfree(f) bind(c)
    import :: c_ptr
    type(c_ptr), value :: f
  end subroutine fgsl_multimin_function_cfree
  subroutine fgsl_multimin_function_fdf_cfree(f) bind(c)
    import :: c_ptr
    type(c_ptr), value :: f
  end subroutine fgsl_multimin_function_fdf_cfree
  function fgsl_aux_multimin_fminimizer_alloc(it) bind(c)
    import
    integer(c_int), value :: it
    type(c_ptr) :: fgsl_aux_multimin_fminimizer_alloc
  end function fgsl_aux_multimin_fminimizer_alloc
  function fgsl_aux_multimin_fdfminimizer_alloc(it) bind(c)
    import
    integer(c_int), value :: it
    type(c_ptr) :: fgsl_aux_multimin_fdfminimizer_alloc
  end function fgsl_aux_multimin_fdfminimizer_alloc
!-*-f90-*-
!
!  Interfaces: Fitting
!
  function gsl_fit_linear(x, xstride, y, ystride, n, c0, c1, &
       cov00, cov01, cov11, sumsq) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: x, y
    integer(c_size_t), value :: xstride, ystride, n
    real(c_double), intent(out) :: c0, c1, cov00, cov01, cov11, sumsq
    integer(c_int) :: gsl_fit_linear
  end function gsl_fit_linear
  function gsl_fit_wlinear(x, xstride, w, wstride, y, ystride, n, c0, c1, &
       cov00, cov01, cov11, chisq) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: x, y, w
    integer(c_size_t), value :: xstride, ystride, wstride, n
    real(c_double), intent(out) :: c0, c1, cov00, cov01, cov11, chisq
    integer(c_int) :: gsl_fit_wlinear
  end function gsl_fit_wlinear
  function gsl_fit_linear_est(x, c0, c1, cov00, cov01, cov11, y, y_err) bind(c)
    import
    real(c_double), value :: x, c0, c1, cov00, cov01, cov11
    real(c_double), intent(out) ::  y, y_err
    integer(c_int) :: gsl_fit_linear_est
  end function gsl_fit_linear_est
  function gsl_fit_mul(x, xstride, y, ystride, n, c1, &
        cov11, sumsq) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: x, y
    integer(c_size_t), value :: xstride, ystride, n
    real(c_double), intent(out) :: c1, cov11, sumsq
    integer(c_int) :: gsl_fit_mul
  end function gsl_fit_mul
  function gsl_fit_wmul(x, xstride, w, wstride, y, ystride, n, c1, &
        cov11, chisq) bind(c)
    import
    real(c_double), dimension(*), intent(in) :: x, y, w
    integer(c_size_t), value :: xstride, ystride, wstride, n
    real(c_double), intent(out) :: c1, cov11, chisq
    integer(c_int) :: gsl_fit_wmul
  end function gsl_fit_wmul
  function gsl_fit_mul_est(x, c1, cov11, y, y_err) bind(c)
    import
    real(c_double), value :: x, c1, cov11
    real(c_double), intent(out) ::  y, y_err
    integer(c_int) :: gsl_fit_mul_est
  end function gsl_fit_mul_est
  function gsl_multifit_linear_alloc(n, p) bind(c)
    import
    integer(c_size_t), value :: n, p
    type(c_ptr) :: gsl_multifit_linear_alloc
  end function gsl_multifit_linear_alloc
  subroutine gsl_multifit_linear_free(w) bind(c)
    import
    type(c_ptr), value :: w
  end subroutine gsl_multifit_linear_free
  function gsl_multifit_linear(x, y, c, cov, chisq, work) bind(c)
    import :: c_ptr, c_int, c_double
    type(c_ptr), value :: x, y, c, cov, work
    real(c_double) :: chisq
    integer(c_int) :: gsl_multifit_linear
  end function gsl_multifit_linear
  function gsl_multifit_linear_svd(x, y, tol, rank, c, cov, chisq, work) bind(c)
    import :: c_ptr, c_int, c_double, c_size_t
    type(c_ptr), value :: x, y, c, cov, work
    real(c_double), value :: tol
    integer(c_size_t) :: rank
    real(c_double) :: chisq
    integer(c_int) :: gsl_multifit_linear_svd
  end function gsl_multifit_linear_svd
  function gsl_multifit_wlinear(x, w, y, c, cov, chisq, work) bind(c)
    import :: c_ptr, c_int, c_double
    type(c_ptr), value :: x, w, y, c, cov, work
    real(c_double) :: chisq
    integer(c_int) :: gsl_multifit_wlinear
  end function gsl_multifit_wlinear
  function gsl_multifit_wlinear_svd(x, w, y, tol, rank, c, cov, chisq, work) bind(c)
    import :: c_ptr, c_int, c_double, c_size_t
    type(c_ptr), value :: x, w, y, c, cov, work
    real(c_double), value :: tol
    integer(c_size_t) :: rank
    real(c_double) :: chisq
    integer(c_int) :: gsl_multifit_wlinear_svd
  end function gsl_multifit_wlinear_svd
  function gsl_multifit_linear_est(x, c, cov, y, y_err) bind(c)
    import :: c_ptr, c_double, c_int
    type(c_ptr), value :: x, c, cov
    real(c_double) :: y, y_err
    integer(c_int) :: gsl_multifit_linear_est
  end function gsl_multifit_linear_est
  function gsl_multifit_linear_residuals(x, y, c, r) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: x
    type(c_ptr), value :: y, c
    type(c_ptr), value :: r
    integer(c_int) :: gsl_multifit_linear_residuals
  end function gsl_multifit_linear_residuals
!-*-f90-*-
!
!  Interfaces: non-linear least-squares fitting
!
  function gsl_multifit_fsolver_alloc(t, n, p) bind(c)
    import
    type(c_ptr), value :: t
    integer(c_size_t), value :: n, p
    type(c_ptr) :: gsl_multifit_fsolver_alloc
  end function gsl_multifit_fsolver_alloc
  function gsl_multifit_fdfsolver_alloc(t, n, p) bind(c)
    import
    type(c_ptr), value :: t
    integer(c_size_t), value :: n, p
    type(c_ptr) :: gsl_multifit_fdfsolver_alloc
  end function gsl_multifit_fdfsolver_alloc
  function gsl_multifit_fsolver_set(s, f, x) bind(c)
    import
    type(c_ptr), value :: s, f, x
    integer(c_int) :: gsl_multifit_fsolver_set
  end function gsl_multifit_fsolver_set
  function gsl_multifit_fdfsolver_set(s, f, x) bind(c)
    import
    type(c_ptr), value :: s, f, x
    integer(c_int) :: gsl_multifit_fdfsolver_set
  end function gsl_multifit_fdfsolver_set
  subroutine gsl_multifit_fsolver_free(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
  end subroutine gsl_multifit_fsolver_free
  subroutine gsl_multifit_fdfsolver_free(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
  end subroutine gsl_multifit_fdfsolver_free
  function gsl_multifit_fsolver_name(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multifit_fsolver_name
  end function gsl_multifit_fsolver_name
  function gsl_multifit_fdfsolver_name(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multifit_fdfsolver_name
  end function gsl_multifit_fdfsolver_name
  function gsl_multifit_fsolver_iterate(s) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: s
    integer(c_int) :: gsl_multifit_fsolver_iterate
  end function gsl_multifit_fsolver_iterate
  function gsl_multifit_fdfsolver_iterate(s) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: s
    integer(c_int) :: gsl_multifit_fdfsolver_iterate
  end function gsl_multifit_fdfsolver_iterate
  function gsl_multifit_fsolver_position(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multifit_fsolver_position
  end function gsl_multifit_fsolver_position
  function gsl_multifit_fdfsolver_position(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multifit_fdfsolver_position
  end function gsl_multifit_fdfsolver_position
  function gsl_multifit_fdfsolver_dx(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multifit_fdfsolver_dx
  end function gsl_multifit_fdfsolver_dx
  function gsl_multifit_fdfsolver_f(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multifit_fdfsolver_f
  end function gsl_multifit_fdfsolver_f
  function gsl_multifit_fdfsolver_jac(s) bind(c)
    import :: c_ptr
    type(c_ptr), value :: s
    type(c_ptr) :: gsl_multifit_fdfsolver_jac
  end function gsl_multifit_fdfsolver_jac
  function gsl_multifit_test_delta(dx, x, epsabs, epsrel) bind(c)
    import :: c_ptr, c_int, c_double
    type(c_ptr), value :: dx, x
    real(c_double), value :: epsabs, epsrel
    integer(c_int) :: gsl_multifit_test_delta
  end function gsl_multifit_test_delta
  function gsl_multifit_test_gradient(g, epsabs) bind(c)
    import :: c_ptr, c_int, c_double
    type(c_ptr), value :: g
    real(c_double), value :: epsabs
    integer(c_int) :: gsl_multifit_test_gradient
  end function gsl_multifit_test_gradient
  function gsl_multifit_gradient(j, f, g) bind(c)
    import :: c_ptr, c_int
    type(c_ptr), value :: j, f, g
    integer(c_int) :: gsl_multifit_gradient
  end function gsl_multifit_gradient
  function gsl_multifit_covar(j, epsrel, cov) bind(c)
    import :: c_ptr, c_int, c_double
    type(c_ptr), value :: j, cov
    real(c_double), value :: epsrel
    integer(c_int) :: gsl_multifit_covar
  end function gsl_multifit_covar
!
  function fgsl_multifit_function_cinit(fp, ndim, p, params) bind(c)
    import 
    type(c_funptr), value :: fp
    integer(c_size_t), value :: ndim, p
    type(c_ptr), value :: params
    type(c_ptr) :: fgsl_multifit_function_cinit
  end function fgsl_multifit_function_cinit
  function fgsl_multifit_function_fdf_cinit(fp, dfp, fdfp, ndim, p, params) bind(c)
    import
    type(c_funptr), value :: fp, dfp, fdfp
    integer(c_size_t), value :: ndim, p
    type(c_ptr), value :: params
    type(c_ptr) :: fgsl_multifit_function_fdf_cinit
  end function fgsl_multifit_function_fdf_cinit
  subroutine fgsl_multifit_function_cfree(f) bind(c)
    import :: c_ptr
    type(c_ptr), value :: f
  end subroutine fgsl_multifit_function_cfree
  subroutine fgsl_multifit_function_fdf_cfree(f) bind(c)
    import :: c_ptr
    type(c_ptr), value :: f
  end subroutine fgsl_multifit_function_fdf_cfree
  function fgsl_aux_multifit_fsolver_alloc(it) bind(c)
    import
    integer(c_int), value :: it
    type(c_ptr) :: fgsl_aux_multifit_fsolver_alloc
  end function fgsl_aux_multifit_fsolver_alloc
  function fgsl_aux_multifit_fdfsolver_alloc(it) bind(c)
    import
    integer(c_int), value :: it
    type(c_ptr) :: fgsl_aux_multifit_fdfsolver_alloc
  end function fgsl_aux_multifit_fdfsolver_alloc
!-*-f90-*-
!
!  Interfaces: Basis splines
!
function gsl_bspline_alloc(k, nbreak) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: k, nbreak
  type(c_ptr) :: gsl_bspline_alloc
end function gsl_bspline_alloc
subroutine gsl_bspline_free (w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_bspline_free
function gsl_bspline_deriv_alloc(k) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t), value :: k
  type(c_ptr) :: gsl_bspline_deriv_alloc
end function gsl_bspline_deriv_alloc
subroutine gsl_bspline_deriv_free (w) bind(c)
  import :: c_ptr
  type(c_ptr), value :: w
end subroutine gsl_bspline_deriv_free
function gsl_bspline_knots(breakpts, w) bind(c)
  import :: c_int, c_ptr
  integer(c_int) :: gsl_bspline_knots
  type(c_ptr), value :: breakpts, w
end function gsl_bspline_knots
function gsl_bspline_knots_uniform(a, b, w) bind(c)
  import :: c_int, c_ptr, c_double
  integer(c_int) :: gsl_bspline_knots_uniform
  real(c_double), value :: a, b
  type(c_ptr), value :: w
end function gsl_bspline_knots_uniform
function gsl_bspline_eval(x, b, w) bind(c)
  import :: c_int, c_ptr, c_double
  integer(c_int) :: gsl_bspline_eval
  real(c_double), value :: x
  type(c_ptr), value :: b, w
end function gsl_bspline_eval
function gsl_bspline_eval_nonzero(x, b, istart, iend, w) bind(c)
  import :: c_int, c_ptr, c_double, c_size_t
  integer(c_int) :: gsl_bspline_eval_nonzero
  real(c_double), value :: x
  integer(c_size_t) :: istart, iend
  type(c_ptr), value :: b, w
end function gsl_bspline_eval_nonzero
function gsl_bspline_deriv_eval(x, nderiv, b, w, dw) bind(c)
  import :: c_int, c_ptr, c_double, c_size_t
  integer(c_int) :: gsl_bspline_deriv_eval
  integer(c_size_t), value :: nderiv
  real(c_double), value :: x
  type(c_ptr), value :: b, w, dw
end function gsl_bspline_deriv_eval
function gsl_bspline_deriv_eval_nonzero(x, nderiv, b, istart, &
     iend, w, dw) bind(c)
  import :: c_int, c_ptr, c_double, c_size_t
  integer(c_int) :: gsl_bspline_deriv_eval_nonzero
  real(c_double), value :: x
  integer(c_size_t), value :: nderiv
  integer(c_size_t) :: istart, iend
  type(c_ptr), value :: b, w, dw
end function gsl_bspline_deriv_eval_nonzero
function gsl_bspline_ncoeffs(w) bind(c)
  import :: c_size_t, c_ptr
  integer(c_size_t) :: gsl_bspline_ncoeffs
  type(c_ptr), value :: w
end function gsl_bspline_ncoeffs
function gsl_bspline_greville_abscissa(i, w) bind(c)
  import :: c_size_t, c_double, c_ptr
  real(c_double) :: gsl_bspline_greville_abscissa
  integer(c_size_t) :: i
  type(c_ptr), value :: w
end function gsl_bspline_greville_abscissa
!-*-f90-*-
!
!  Interfaces: IEEE
!
  subroutine gsl_ieee_fprintf_float(stream, x) bind(c)
    import :: c_ptr, c_float
    type(c_ptr), value :: stream
    real(c_float) :: x
  end subroutine gsl_ieee_fprintf_float
  subroutine gsl_ieee_fprintf_double(stream, x) bind(c)
    import :: c_ptr, c_double
    type(c_ptr), value :: stream
    real(c_double) :: x
  end subroutine gsl_ieee_fprintf_double
  subroutine gsl_ieee_printf_float(x) bind(c)
    import :: c_float
    real(c_float) :: x
  end subroutine gsl_ieee_printf_float
  subroutine gsl_ieee_printf_double(x) bind(c)
    import :: c_double
    real(c_double) :: x
  end subroutine gsl_ieee_printf_double
  subroutine gsl_ieee_env_setup() bind(c)
  end subroutine gsl_ieee_env_setup
end interface
!-*-f90-*-
!
!  Interfaces: Generics
!
!  FIXME: missing still
!  interface fgsl_error
!     module procedure fgsl_err_info
!     module procedure fgsl_err_noinfo
!  end interface
  interface fgsl_well_defined
     module procedure fgsl_vector_status
     module procedure fgsl_matrix_status
     module procedure fgsl_vector_complex_status
     module procedure fgsl_matrix_complex_status
     module procedure fgsl_cheb_series_status
     module procedure fgsl_interp_status
     module procedure fgsl_dht_status
     module procedure fgsl_error_handler_status
     module procedure fgsl_integration_workspace_status
     module procedure fgsl_integration_qawo_table_status
     module procedure fgsl_integration_qaws_table_status
     module procedure fgsl_interp_accel_status
     module procedure fgsl_spline_status
     module procedure fgsl_permutation_status
     module procedure fgsl_combination_status
     module procedure fgsl_odeiv_control_status
     module procedure fgsl_odeiv_evolve_status
     module procedure fgsl_odeiv_step_status
     module procedure fgsl_odeiv_system_status
     module procedure fgsl_poly_complex_workspace_stat
     module procedure fgsl_rng_status
     module procedure fgsl_qrng_status
     module procedure fgsl_ran_discrete_t_status
     module procedure fgsl_root_fsolver_status
     module procedure fgsl_root_fdfsolver_status
     module procedure fgsl_siman_params_t_status
     module procedure fgsl_min_fminimizer_status
     module procedure fgsl_histogram_status
     module procedure fgsl_ntuple_status
     module procedure fgsl_ntuple_value_fn_status
     module procedure fgsl_ntuple_select_fn_status
     module procedure fgsl_monte_function_status
     module procedure fgsl_monte_plain_status
     module procedure fgsl_monte_miser_status
     module procedure fgsl_monte_vegas_status
     module procedure fgsl_multiroot_fsolver_status
     module procedure fgsl_multiroot_fdfsolver_status
     module procedure fgsl_multimin_fminimizer_status
     module procedure fgsl_multimin_fdfminimizer_status
     module procedure fgsl_multifit_status
     module procedure fgsl_multifit_fsolver_status
     module procedure fgsl_multifit_fdfsolver_status
     module procedure fgsl_file_status
     module procedure fgsl_wavelet_status
     module procedure fgsl_wavelet_workspace_status
  end interface
  interface fgsl_sizeof
     module procedure fgsl_sizeof_double
     module procedure fgsl_sizeof_float

     !include "integer.finc"     
     ! The line below when compiling for 64 bit
     module procedure fgsl_sizeof_int
     module procedure fgsl_sizeof_size_t

     ! The line below when compiling for 32 bit
     !module procedure fgsl_sizeof_int
     
     module procedure fgsl_sizeof_char
     module procedure fgsl_sizeof_vector
     module procedure fgsl_sizeof_matrix
     module procedure fgsl_sizeof_vector_complex
     module procedure fgsl_sizeof_matrix_complex
     module procedure fgsl_sizeof_interp
     module procedure fgsl_sizeof_permutation
     module procedure fgsl_sizeof_combination
     module procedure fgsl_sizeof_integration_workspace
     module procedure fgsl_sizeof_integration_qaws_table
     module procedure fgsl_sizeof_integration_qawo_table
     module procedure fgsl_sizeof_wavelet
     module procedure fgsl_sizeof_wavelet_workspace
  end interface
  interface fgsl_obj_c_ptr
     module procedure fgsl_rng_c_ptr 
     module procedure fgsl_vector_c_ptr 
     module procedure fgsl_matrix_c_ptr 
  end interface
  interface assignment(=)
     module procedure fgsl_complex_to_complex
     module procedure complex_to_fgsl_complex
     module procedure gsl_sf_to_fgsl_sf
     module procedure gsl_sfe10_to_fgsl_sfe10
     module procedure fgsl_vector_to_array
     module procedure fgsl_vector_complex_to_array
     module procedure fgsl_matrix_to_array
     module procedure fgsl_matrix_complex_to_array
  end interface
!
! Array processing
!
  interface fgsl_vector_init
     module procedure fgsl_vector_init
     module procedure fgsl_vector_complex_init
  end interface
  interface fgsl_vector_free
     module procedure fgsl_vector_free
     module procedure fgsl_vector_complex_free
  end interface
  interface fgsl_matrix_init
     module procedure fgsl_matrix_init
     module procedure fgsl_matrix_complex_init
  end interface
  interface fgsl_matrix_free
     module procedure fgsl_matrix_free
     module procedure fgsl_matrix_complex_free
  end interface
  interface fgsl_vector_align
     module procedure fgsl_vector_align
     module procedure fgsl_vector_complex_align
     module procedure fgsl_vector_pointer_align
     module procedure fgsl_vector_complex_pointer_align
  end interface
  interface fgsl_matrix_align
     module procedure fgsl_matrix_align
     module procedure fgsl_matrix_pointer_align
     module procedure fgsl_matrix_complex_align
     module procedure fgsl_matrix_complex_pointer_align
  end interface
!
! Permutations and combinations
!
  interface fgsl_permute
     module procedure fgsl_permute
     module procedure fgsl_permute_long
  end interface
  interface fgsl_permute_inverse
     module procedure fgsl_permute_inverse
     module procedure fgsl_permute_long_inverse
  end interface
!
! Sorting
!
  interface fgsl_sort
     module procedure fgsl_sort_double
     module procedure fgsl_sort_long
  end interface
  interface fgsl_sort_index
     module procedure fgsl_sort_double_index
     module procedure fgsl_sort_long_index
  end interface
  interface fgsl_sort_smallest
     module procedure fgsl_sort_double_smallest
     module procedure fgsl_sort_long_smallest
  end interface
  interface fgsl_sort_smallest_index
     module procedure fgsl_sort_double_smallest_index
     module procedure fgsl_sort_long_smallest_index
  end interface
  interface fgsl_sort_largest
     module procedure fgsl_sort_double_largest
     module procedure fgsl_sort_long_largest
  end interface
  interface fgsl_sort_largest_index
     module procedure fgsl_sort_double_largest_index
     module procedure fgsl_sort_long_largest_index
  end interface
!
! Random number stuff
!
  interface fgsl_ran_shuffle
     module procedure fgsl_ran_shuffle
     module procedure fgsl_ran_shuffle_double
     module procedure fgsl_ran_shuffle_size_t
  end interface
!
! IEEE
!
  interface fgsl_ieee_fprintf
     module procedure fgsl_ieee_fprintf_float
     module procedure fgsl_ieee_fprintf_double
  end interface
  interface fgsl_ieee_printf
     module procedure fgsl_ieee_printf_float
     module procedure fgsl_ieee_printf_double
  end interface
     
     
contains
!-*-f90-*-
!
!  API: Error treatment
!  FIXME: introduction of own error handler requires
!
  function fgsl_set_error_handler(new_handler)
    type(fgsl_error_handler_t), intent(in) :: new_handler
    type(fgsl_error_handler_t) :: fgsl_set_error_handler
    fgsl_set_error_handler%gsl_error_handler_t = &
         gsl_set_error_handler(new_handler%gsl_error_handler_t)
  end function fgsl_set_error_handler
  function fgsl_set_error_handler_off()
    type(fgsl_error_handler_t) :: fgsl_set_error_handler_off
    fgsl_set_error_handler_off%gsl_error_handler_t = &
         gsl_set_error_handler_off()
  end function fgsl_set_error_handler_off
  function fgsl_strerror(errno) 
    integer(fgsl_int), intent(in) :: errno
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_strerror
!
    type(c_ptr) :: name
!
!    write(6, *) 'Error is ',errno
    name = gsl_strerror(errno)
    fgsl_strerror = fgsl_name(name)
  end function fgsl_strerror
  subroutine fgsl_error(reason, file, line, errno)
    character(kind=fgsl_char,len=*), intent(in) :: &
         reason, file
    integer(fgsl_int), intent(in) :: line, errno
    call gsl_error(reason//c_null_char, file//c_null_char, line, errno)
  end subroutine fgsl_error
  function fgsl_error_handler_status(error_handler_t)
    type(fgsl_error_handler_t), intent(in) :: error_handler_t
    logical :: fgsl_error_handler_status
    fgsl_error_handler_status = .true.
    if (.not. c_associated(error_handler_t%gsl_error_handler_t)) &
         fgsl_error_handler_status = .false.
  end function fgsl_error_handler_status
!
! initialize own error handler
!
  function fgsl_error_handler_init(handler_sr)
    interface
       subroutine handler_sr(reason, file, line, errno) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: reason, file
         integer(c_int), value :: line, errno
       end subroutine handler_sr
    end interface
    type(fgsl_error_handler_t) :: fgsl_error_handler_init
!
    type(c_funptr) :: fptr
    fptr = c_funloc(handler_sr)
    fgsl_error_handler_init%gsl_error_handler_t = fptr
  end function fgsl_error_handler_init
!-*-f90-*-
!
!  API: miscellaneous additions
!
! standardized C string to Fortran string conversion
  function fgsl_name(c_name)
    type(c_ptr), intent(in) :: c_name
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_name
!
    character(kind=fgsl_char,len=fgsl_strmax) :: result
    character(fgsl_char), pointer :: fc(:)
    integer :: len
!
    call c_f_pointer(c_name,fc,(/fgsl_strmax/))
    len = 1
    do 
       if (fc(len) == c_null_char .or. len > fgsl_strmax) exit
       result(len:len) = fc(len)
       len = len + 1
    end do
    if (fc(len) == c_null_char) len = len - 1
    fgsl_name = result(1:len)
  end function fgsl_name
! sizes of intrinsic types
  function fgsl_sizeof_double(x)
    real(fgsl_double), intent(in) :: x
    integer(fgsl_size_t) :: fgsl_sizeof_double
    fgsl_sizeof_double = gsl_aux_sizeof_double()
  end function fgsl_sizeof_double
  function fgsl_sizeof_float(x)
    real(fgsl_float), intent(in) :: x
    integer(fgsl_size_t) :: fgsl_sizeof_float
    fgsl_sizeof_float = gsl_aux_sizeof_float()
  end function fgsl_sizeof_float
  function fgsl_sizeof_int(x)
    integer(fgsl_int), intent(in) :: x
    integer(fgsl_size_t) :: fgsl_sizeof_int
    fgsl_sizeof_int = gsl_aux_sizeof_int()
  end function fgsl_sizeof_int
  function fgsl_sizeof_long(x)
    integer(fgsl_long), intent(in) :: x
    integer(fgsl_size_t) :: fgsl_sizeof_long
    fgsl_sizeof_long = gsl_aux_sizeof_long()
  end function fgsl_sizeof_long
  function fgsl_sizeof_size_t(x)
    integer(fgsl_size_t), intent(in) :: x
    integer(fgsl_size_t) :: fgsl_sizeof_size_t
    fgsl_sizeof_size_t = gsl_aux_sizeof_size_t()
  end function fgsl_sizeof_size_t
  function fgsl_sizeof_char(x)
    character(fgsl_char), intent(in) :: x
    integer(fgsl_size_t) :: fgsl_sizeof_char
    fgsl_sizeof_char = gsl_aux_sizeof_char()
  end function fgsl_sizeof_char
!-*-f90-*-
!
!  API: I/O Add-ons
!
  function fgsl_open(path, mode)
    character(kind=fgsl_char, len=*), intent(in) :: path, mode
!    character(kind=fgsl_char), dimension(*), intent(in) :: path, mode
    type(fgsl_file) :: fgsl_open
    character(kind=fgsl_char,len=fgsl_pathmax) :: lpath
    character(kind=fgsl_char,len=fgsl_strmax) :: lmode
    if (len(trim(path)) < fgsl_pathmax .and. len(trim(mode)) < fgsl_strmax) then
       lpath = trim(path) // c_null_char
       lmode = trim(mode) // c_null_char
       fgsl_open%gsl_file = fopen(lpath, lmode)
    else
       fgsl_open%gsl_file = c_null_ptr
    end if
  end function fgsl_open
  function fgsl_close(fd)
    type(fgsl_file) :: fd
    integer(fgsl_int) :: fgsl_close
!
    integer(c_int) :: status
    status = fclose(fd%gsl_file)
    fgsl_close = fgsl_success
    if (status /= 0) fgsl_close = fgsl_efault
  end function fgsl_close
  function fgsl_stdin()
    type(fgsl_file) :: fgsl_stdin
    fgsl_stdin%gsl_file = fgsl_cstdin()
  end function fgsl_stdin
  function fgsl_stdout()
    type(fgsl_file) :: fgsl_stdout
    fgsl_stdout%gsl_file = fgsl_cstdout()
  end function fgsl_stdout
  function fgsl_stderr()
    type(fgsl_file) :: fgsl_stderr
    fgsl_stderr%gsl_file = fgsl_cstderr()
  end function fgsl_stderr
  function fgsl_flush(file) 
    type(fgsl_file), intent(in) :: file
    integer(fgsl_int) :: fgsl_flush
    fgsl_flush = fflush(file%gsl_file)
  end function fgsl_flush
  function fgsl_file_status(file)
    type(fgsl_file), intent(in) :: file
    logical :: fgsl_file_status
    fgsl_file_status = .true. 
    if (.not. c_associated(file%gsl_file)) fgsl_file_status = .false.
  end function fgsl_file_status
!-*-f90-*-
!
!  API: Mathematical Functions
!
  function fgsl_isnan(x) 
    real(fgsl_double), intent(in) :: x
    integer(fgsl_int) :: fgsl_isnan
    fgsl_isnan = gsl_isnan(x)
  end function fgsl_isnan
  function fgsl_isinf(x)
    real(fgsl_double), intent(in) :: x
    integer(fgsl_int) :: fgsl_isinf
    fgsl_isinf = gsl_isinf(x)
  end function fgsl_isinf
  function fgsl_finite(x)
    real(fgsl_double), intent(in) :: x
    integer(fgsl_int) :: fgsl_finite
    fgsl_finite = gsl_finite(x)
  end function fgsl_finite
  function fgsl_log1p(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_log1p
    fgsl_log1p = gsl_log1p(x)
  end function fgsl_log1p
  function fgsl_expm1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_expm1
    fgsl_expm1 = gsl_expm1(x)
  end function fgsl_expm1
  function fgsl_hypot(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_hypot
    fgsl_hypot = gsl_hypot(x)
  end function fgsl_hypot
  function fgsl_acosh(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_acosh
    fgsl_acosh = gsl_acosh(x)
  end function fgsl_acosh
  function fgsl_asinh(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_asinh
    fgsl_asinh = gsl_asinh(x)
  end function fgsl_asinh
  function fgsl_atanh(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_atanh
    fgsl_atanh = gsl_atanh(x)
  end function fgsl_atanh
  function fgsl_ldexp(x,e) 
    real(fgsl_double), intent(in) :: x
    integer(fgsl_int), intent(in) :: e
    real(fgsl_double) :: fgsl_ldexp
    fgsl_ldexp = gsl_ldexp(x,e)
  end function fgsl_ldexp
  function fgsl_frexp(x,e) 
    real(fgsl_double), intent(in) :: x
    integer(fgsl_int), intent(out) :: e
    real(fgsl_double) :: fgsl_frexp
    fgsl_frexp = gsl_frexp(x,e)
  end function fgsl_frexp
  function fgsl_fcmp(x,y,eps)
    real(fgsl_double), intent(in) :: x, y, eps
    integer(fgsl_int) :: fgsl_fcmp
    fgsl_fcmp = gsl_fcmp(x,y,eps)
  end function fgsl_fcmp
! constructors for abstract types
  function fgsl_function_init(func, params)
    interface
       function func(x, params) bind(c)
         use, intrinsic :: iso_c_binding
         real(c_double), value :: x
         type(c_ptr), value :: params
         real(c_double) :: func
       end function func
    end interface
    type(c_ptr), intent(in) :: params
    type(fgsl_function) :: fgsl_function_init
!
    type(c_funptr) :: fp
    fp = c_funloc(func)
    fgsl_function_init%gsl_function = fgsl_function_cinit(fp, params)
  end function fgsl_function_init
  function fgsl_function_fdf_init(f, df, fdf, params)
    interface
       function f(x, params) bind(c)
         use, intrinsic :: iso_c_binding
         real(c_double), value :: x
         type(c_ptr), value :: params
         real(c_double) :: f
       end function f
       function df(x, params) bind(c)
         use, intrinsic :: iso_c_binding
         real(c_double), value :: x
         type(c_ptr), value :: params
         real(c_double) :: df
       end function df
       subroutine fdf(x, params, f, df) bind(c)
         use, intrinsic :: iso_c_binding
         real(c_double), value :: x
         type(c_ptr), value :: params
         real(c_double), intent(out) :: f, df
       end subroutine fdf
    end interface
    type(c_ptr), intent(in) :: params
    type(fgsl_function_fdf) :: fgsl_function_fdf_init
!
    type(c_funptr) :: fp, dfp, fdfp
    fp = c_funloc(f)
    dfp = c_funloc(df)
    fdfp = c_funloc(fdf)
    fgsl_function_fdf_init%gsl_function_fdf = fgsl_function_fdf_cinit(fp, dfp, fdfp, params)
  end function fgsl_function_fdf_init
  subroutine fgsl_function_free(sfunc)
    type(fgsl_function), intent(inout) :: sfunc
    call fgsl_function_cfree(sfunc%gsl_function)
  end subroutine fgsl_function_free
  subroutine fgsl_function_fdf_free(sfunc)
    type(fgsl_function_fdf), intent(inout) :: sfunc
    call fgsl_function_fdf_cfree(sfunc%gsl_function_fdf)
  end subroutine fgsl_function_fdf_free
  function fgsl_fn_eval(sfunc, x)
    type(fgsl_function), intent(inout) :: sfunc
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_fn_eval
    fgsl_fn_eval = fgsl_fn_eval_aux(sfunc%gsl_function, x)
  end function fgsl_fn_eval
  function fgsl_fn_fdf_eval_f(sfunc, x)
    type(fgsl_function_fdf), intent(inout) :: sfunc
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_fn_fdf_eval_f
    fgsl_fn_fdf_eval_f = fgsl_fn_fdf_eval_f_aux(sfunc%gsl_function_fdf, x)
  end function fgsl_fn_fdf_eval_f
  function fgsl_fn_fdf_eval_df(sfunc, x)
    type(fgsl_function_fdf), intent(inout) :: sfunc
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_fn_fdf_eval_df
    fgsl_fn_fdf_eval_df = fgsl_fn_fdf_eval_df_aux(sfunc%gsl_function_fdf, x)
  end function fgsl_fn_fdf_eval_df
  subroutine fgsl_fn_fdf_eval_f_df(sfunc, x, y, dy)
    type(fgsl_function_fdf), intent(inout) :: sfunc
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: y, dy
    call fgsl_fn_fdf_eval_f_df_aux(sfunc%gsl_function_fdf, x, y, dy)
  end subroutine fgsl_fn_fdf_eval_f_df
!-*-f90-*-
!
!  API: Complex numbers 
!       (internal conversion routines)
!
  function fgsl_complex_arg(z)
    complex(fgsl_double_complex), intent(in) :: z
    real(fgsl_double) :: fgsl_complex_arg
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arg = gsl_complex_arg(zz)
  end function fgsl_complex_arg
  function fgsl_complex_logabs(z)
    complex(fgsl_double_complex), intent(in) :: z
    real(fgsl_double) :: fgsl_complex_logabs
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_logabs = gsl_complex_logabs(zz)
  end function fgsl_complex_logabs
  function fgsl_complex_log10(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_log10
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_log10 = gsl_complex_log10(zz)
  end function fgsl_complex_log10
  function fgsl_complex_log_b(z,b)
    complex(fgsl_double_complex), intent(in) :: z, b
    complex(fgsl_double_complex) :: fgsl_complex_log_b
!
    type(gsl_complex) :: zz, bb
    zz = z
    bb = b
    fgsl_complex_log_b = gsl_complex_log_b(zz, bb)
  end function fgsl_complex_log_b
  function fgsl_complex_arcsin(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arcsin
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arcsin = gsl_complex_arcsin(zz)
  end function fgsl_complex_arcsin
  function fgsl_complex_arcsin_real(r)
    real(fgsl_double), intent(in) :: r
    complex(fgsl_double_complex) :: fgsl_complex_arcsin_real
    fgsl_complex_arcsin_real = gsl_complex_arcsin_real(r)
  end function fgsl_complex_arcsin_real
  function fgsl_complex_arccos(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arccos
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arccos = gsl_complex_arccos(zz)
  end function fgsl_complex_arccos
  function fgsl_complex_arccos_real(r)
    real(fgsl_double), intent(in) :: r
    complex(fgsl_double_complex) :: fgsl_complex_arccos_real
    fgsl_complex_arccos_real = gsl_complex_arccos_real(r)
  end function fgsl_complex_arccos_real
  function fgsl_complex_arctan(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arctan
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arctan = gsl_complex_arctan(zz)
  end function fgsl_complex_arctan
  function fgsl_complex_arcsec(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arcsec
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arcsec = gsl_complex_arcsec(zz)
  end function fgsl_complex_arcsec
  function fgsl_complex_arcsec_real(r)
    real(fgsl_double), intent(in) :: r
    complex(fgsl_double_complex) :: fgsl_complex_arcsec_real
    fgsl_complex_arcsec_real = gsl_complex_arcsec_real(r)
  end function fgsl_complex_arcsec_real
  function fgsl_complex_arccsc(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arccsc
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arccsc = gsl_complex_arccsc(zz)
  end function fgsl_complex_arccsc
  function fgsl_complex_arccsc_real(r)
    real(fgsl_double), intent(in) :: r
    complex(fgsl_double_complex) :: fgsl_complex_arccsc_real
    fgsl_complex_arccsc_real = gsl_complex_arccsc_real(r)
  end function fgsl_complex_arccsc_real
  function fgsl_complex_arccot(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arccot
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arccot = gsl_complex_arccot(zz)
  end function fgsl_complex_arccot
  function fgsl_complex_arcsinh(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arcsinh
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arcsinh = gsl_complex_arcsinh(zz)
  end function fgsl_complex_arcsinh
  function fgsl_complex_arccosh(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arccosh
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arccosh = gsl_complex_arccosh(zz)
  end function fgsl_complex_arccosh
  function fgsl_complex_arccosh_real(r)
    real(fgsl_double), intent(in) :: r
    complex(fgsl_double_complex) :: fgsl_complex_arccosh_real
    fgsl_complex_arccosh_real = gsl_complex_arccosh_real(r)
  end function fgsl_complex_arccosh_real
  function fgsl_complex_arctanh(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arctanh
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arctanh = gsl_complex_arctanh(zz)
  end function fgsl_complex_arctanh
  function fgsl_complex_arctanh_real(r)
    real(fgsl_double), intent(in) :: r
    complex(fgsl_double_complex) :: fgsl_complex_arctanh_real
    fgsl_complex_arctanh_real = gsl_complex_arctanh_real(r)
  end function fgsl_complex_arctanh_real
  function fgsl_complex_arcsech(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arcsech
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arcsech = gsl_complex_arcsech(zz)
  end function fgsl_complex_arcsech
  function fgsl_complex_arccsch(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arccsch
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arccsch = gsl_complex_arccsch(zz)
  end function fgsl_complex_arccsch
  function fgsl_complex_arccoth(z)
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_arccoth
!
    type(gsl_complex) :: zz
    zz = z
    fgsl_complex_arccoth = gsl_complex_arccoth(zz)
  end function fgsl_complex_arccoth
!
  elemental subroutine fgsl_complex_to_complex(result, source)
    complex(fgsl_double_complex), intent(out) :: result
    type(gsl_complex), intent(in) :: source
    result = source%dat(1) + (0.0_fgsl_double, 1.0_fgsl_double) * source%dat(2)
  end subroutine fgsl_complex_to_complex
  elemental subroutine complex_to_fgsl_complex(result, source)
    type(gsl_complex), intent(out) :: result
    complex(fgsl_double_complex), intent(in) :: source
    result%dat(1) = real(source)
    result%dat(2) = aimag(source)
  end subroutine complex_to_fgsl_complex  
!-*-f90-*-
!
! API: Polynomials
!
  function fgsl_poly_eval(c, len, x) 
    real(fgsl_double), intent(in) :: c(:)
    integer(fgsl_int), intent(in) :: len
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_poly_eval
    fgsl_poly_eval = gsl_poly_eval(c, len, x)
  end function fgsl_poly_eval
  function fgsl_poly_complex_eval(c, len, z) 
    real(fgsl_double), intent(in) :: c(:)
    integer(fgsl_int), intent(in) :: len
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_poly_complex_eval
    type(gsl_complex) :: zz
    zz = z
    fgsl_poly_complex_eval = gsl_poly_complex_eval(c, len, zz)
  end function fgsl_poly_complex_eval
  function fgsl_complex_poly_complex_eval(c, len, z) 
    complex(fgsl_double_complex), intent(in) :: c(:)
    integer(fgsl_int), intent(in) :: len
    complex(fgsl_double_complex), intent(in) :: z
    complex(fgsl_double_complex) :: fgsl_complex_poly_complex_eval
    type(gsl_complex) :: zz, cz(len)
    zz = z
    cz = c
    fgsl_complex_poly_complex_eval = gsl_complex_poly_complex_eval(cz, len, zz) 
  end function fgsl_complex_poly_complex_eval
  function fgsl_poly_eval_derivs(c, lenc, x, res, lenres) 
    integer(fgsl_int) :: fgsl_poly_eval_derivs
    integer(fgsl_size_t), intent(in) :: lenc, lenres
    real(fgsl_double), dimension(:), intent(in) :: c
    real(fgsl_double), dimension(:) :: res
    real(fgsl_double), intent(in) :: x
    fgsl_poly_eval_derivs = gsl_poly_eval_derivs(c, lenc, x, res, lenres) 
  end function fgsl_poly_eval_derivs
  function fgsl_poly_dd_init(dd, x, y, size)
    real(fgsl_double), intent(inout) :: dd(:)
    real(fgsl_double), intent(in) :: x(:), y(:)
    integer(fgsl_size_t), intent(in) :: size
    integer(fgsl_int) :: fgsl_poly_dd_init
    fgsl_poly_dd_init = gsl_poly_dd_init(dd, x, y, size)
  end function fgsl_poly_dd_init
  function fgsl_poly_dd_eval(dd, xa, size, x) 
    real(fgsl_double), intent(in) :: dd(:), xa(:)
    real(fgsl_double), intent(in) :: x
    integer(fgsl_size_t), intent(in) :: size
    real(fgsl_double) :: fgsl_poly_dd_eval
    fgsl_poly_dd_eval = gsl_poly_dd_eval(dd, xa, size, x)
  end function fgsl_poly_dd_eval
  function fgsl_poly_dd_taylor(c, xp, dd, x, size, w) 
    real(fgsl_double), intent(inout) :: c(:)
    real(fgsl_double), intent(in) :: xp
    real(fgsl_double), intent(in) :: dd(:), x(:)
    real(fgsl_double), intent(out) :: w(:)
    integer(fgsl_size_t), intent(in) :: size
    integer(fgsl_int) :: fgsl_poly_dd_taylor
    fgsl_poly_dd_taylor = gsl_poly_dd_taylor(c, xp, dd, x, size, w)
  end function fgsl_poly_dd_taylor
  function fgsl_poly_solve_quadratic(a, b, c, x0, x1) 
    real(fgsl_double), intent(in) :: a, b, c
    real(fgsl_double), intent(out) :: x0, x1
    integer(fgsl_int) :: fgsl_poly_solve_quadratic
    fgsl_poly_solve_quadratic = gsl_poly_solve_quadratic(a, b, c, x0, x1)
  end function fgsl_poly_solve_quadratic
  function fgsl_poly_complex_solve_quadratic(a, b, c, x0, x1) 
    real(fgsl_double), intent(in) :: a, b, c
    complex(fgsl_double_complex), intent(out) :: x0, x1
    integer(fgsl_int) :: fgsl_poly_complex_solve_quadratic
!
    type(gsl_complex) :: z0, z1
    fgsl_poly_complex_solve_quadratic = gsl_poly_complex_solve_quadratic(a, b, c, z0, z1)
    x0 = z0 
    x1 = z1
  end function fgsl_poly_complex_solve_quadratic
  function fgsl_poly_solve_cubic(a, b, c, x0, x1, x2) 
    real(fgsl_double), intent(in) :: a, b, c
    real(fgsl_double), intent(out) :: x0, x1, x2
    integer(fgsl_int) :: fgsl_poly_solve_cubic
    fgsl_poly_solve_cubic = gsl_poly_solve_cubic(a, b, c, x0, x1, x2)
  end function fgsl_poly_solve_cubic
  function fgsl_poly_complex_solve_cubic(a, b, c, x0, x1, x2) 
    real(fgsl_double), intent(in) :: a, b, c
    complex(fgsl_double_complex), intent(out) :: x0, x1, x2
    integer(fgsl_int) :: fgsl_poly_complex_solve_cubic
!
    type(gsl_complex) :: z0, z1, z2
    fgsl_poly_complex_solve_cubic = gsl_poly_complex_solve_cubic(a, b, c, z0, z1, z2)
    x0 = z0
    x1 = z1
    x2 = z2
  end function fgsl_poly_complex_solve_cubic
  function fgsl_poly_complex_workspace_alloc(n) 
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_poly_complex_workspace) :: fgsl_poly_complex_workspace_alloc
    fgsl_poly_complex_workspace_alloc%gsl_poly_complex_workspace = gsl_poly_complex_workspace_alloc(n)
  end function fgsl_poly_complex_workspace_alloc
  subroutine fgsl_poly_complex_workspace_free(w) 
    type(fgsl_poly_complex_workspace), intent(inout) :: w
    call gsl_poly_complex_workspace_free(w%gsl_poly_complex_workspace)
  end subroutine fgsl_poly_complex_workspace_free
  function fgsl_poly_complex_workspace_stat(w)
    type(fgsl_poly_complex_workspace), intent(in) :: w
    logical :: fgsl_poly_complex_workspace_stat
    fgsl_poly_complex_workspace_stat = .true.
    if (.not. c_associated(w%gsl_poly_complex_workspace)) then
       fgsl_poly_complex_workspace_stat = .false.
    end if
  end function fgsl_poly_complex_workspace_stat
  function fgsl_poly_complex_solve(a, n, w, z) 
    real(fgsl_double), intent(in) :: a(:)
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_poly_complex_workspace), intent(inout) :: w
    complex(fgsl_double_complex), intent(out) :: z(:)
    integer(fgsl_int) :: fgsl_poly_complex_solve
!
    real(fgsl_double), allocatable :: zz(:)
    integer :: istat, i
    allocate(zz(2*n-2), stat=istat)
    if (istat /= 0) then
       fgsl_poly_complex_solve = fgsl_enomem
    else
       fgsl_poly_complex_solve = gsl_poly_complex_solve(a, n, w%gsl_poly_complex_workspace, zz) 
       do i=1,n-1
          z(i) = zz(2*i-1) + (0.0_fgsl_double, 1.0_fgsl_double) * zz(2*i)
       end do
       deallocate(zz)
    end if
  end function fgsl_poly_complex_solve
!-*-f90-*-
!
! API: Special Functions
!
  function fgsl_sf_airy_ai(x, mode) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_airy_ai
    fgsl_sf_airy_ai = gsl_sf_airy_ai(x, mode%gsl_mode)
  end function fgsl_sf_airy_ai
  function fgsl_sf_airy_ai_e(x, mode, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_airy_ai_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_ai_e = gsl_sf_airy_ai_e(x, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_airy_ai_e
  function fgsl_sf_airy_bi(x, mode) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_airy_bi
    fgsl_sf_airy_bi = gsl_sf_airy_bi(x, mode%gsl_mode)
  end function fgsl_sf_airy_bi
  function fgsl_sf_airy_bi_e(x, mode, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_airy_bi_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_bi_e = gsl_sf_airy_bi_e(x, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_airy_bi_e
  function fgsl_sf_airy_ai_scaled(x, mode) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_airy_ai_scaled
    fgsl_sf_airy_ai_scaled = gsl_sf_airy_ai_scaled(x, mode%gsl_mode)
  end function fgsl_sf_airy_ai_scaled
  function fgsl_sf_airy_ai_scaled_e(x, mode, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_airy_ai_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_ai_scaled_e = gsl_sf_airy_ai_scaled_e(x, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_airy_ai_scaled_e
  function fgsl_sf_airy_bi_scaled(x, mode) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_airy_bi_scaled
    fgsl_sf_airy_bi_scaled = gsl_sf_airy_bi_scaled(x, mode%gsl_mode)
  end function fgsl_sf_airy_bi_scaled
  function fgsl_sf_airy_bi_scaled_e(x, mode, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_airy_bi_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_bi_scaled_e = gsl_sf_airy_bi_scaled_e(x, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_airy_bi_scaled_e
!!!!!
  function fgsl_sf_airy_ai_deriv(x, mode) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_airy_ai_deriv
    fgsl_sf_airy_ai_deriv = gsl_sf_airy_ai_deriv(x, mode%gsl_mode)
  end function fgsl_sf_airy_ai_deriv
  function fgsl_sf_airy_ai_deriv_e(x, mode, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_airy_ai_deriv_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_ai_deriv_e = gsl_sf_airy_ai_deriv_e(x, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_airy_ai_deriv_e
  function fgsl_sf_airy_bi_deriv(x, mode) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_airy_bi_deriv
    fgsl_sf_airy_bi_deriv = gsl_sf_airy_bi_deriv(x, mode%gsl_mode)
  end function fgsl_sf_airy_bi_deriv
  function fgsl_sf_airy_bi_deriv_e(x, mode, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_airy_bi_deriv_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_bi_deriv_e = gsl_sf_airy_bi_deriv_e(x, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_airy_bi_deriv_e
  function fgsl_sf_airy_ai_deriv_scaled(x, mode) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_airy_ai_deriv_scaled
    fgsl_sf_airy_ai_deriv_scaled = gsl_sf_airy_ai_deriv_scaled(x, mode%gsl_mode)
  end function fgsl_sf_airy_ai_deriv_scaled
  function fgsl_sf_airy_ai_deriv_scaled_e(x, mode, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_airy_ai_deriv_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_ai_deriv_scaled_e = gsl_sf_airy_ai_deriv_scaled_e(x, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_airy_ai_deriv_scaled_e
  function fgsl_sf_airy_bi_deriv_scaled(x, mode) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_airy_bi_deriv_scaled
    fgsl_sf_airy_bi_deriv_scaled = gsl_sf_airy_bi_deriv_scaled(x, mode%gsl_mode)
  end function fgsl_sf_airy_bi_deriv_scaled
  function fgsl_sf_airy_bi_deriv_scaled_e(x, mode, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_airy_bi_deriv_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_bi_deriv_scaled_e = gsl_sf_airy_bi_deriv_scaled_e(x, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_airy_bi_deriv_scaled_e
  function fgsl_sf_airy_zero_ai(s)
    integer(fgsl_int), intent(in) :: s
    real(fgsl_double) :: fgsl_sf_airy_zero_ai
    fgsl_sf_airy_zero_ai = gsl_sf_airy_zero_ai(s)
  end function fgsl_sf_airy_zero_ai
  function fgsl_sf_airy_zero_ai_e (s, result) 
    integer(fgsl_int), intent(in) :: s
    integer(fgsl_int) :: fgsl_sf_airy_zero_ai_e
    type(fgsl_sf_result), intent(out) :: result
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_zero_ai_e = gsl_sf_airy_zero_ai_e(s, res)
    result = res
  end function fgsl_sf_airy_zero_ai_e
  function fgsl_sf_airy_zero_bi(s)
    integer(fgsl_int), intent(in) :: s
    real(fgsl_double) :: fgsl_sf_airy_zero_bi
    fgsl_sf_airy_zero_bi = gsl_sf_airy_zero_bi(s)
  end function fgsl_sf_airy_zero_bi
  function fgsl_sf_airy_zero_bi_e (s, result)
    integer(fgsl_int), intent(in) :: s
    integer(fgsl_int) :: fgsl_sf_airy_zero_bi_e
    type(fgsl_sf_result), intent(out) :: result
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_zero_bi_e = gsl_sf_airy_zero_bi_e(s, res)
    result = res
  end function fgsl_sf_airy_zero_bi_e
  function fgsl_sf_airy_zero_ai_deriv(s)
    integer(fgsl_int), intent(in) :: s
    real(fgsl_double) :: fgsl_sf_airy_zero_ai_deriv
    fgsl_sf_airy_zero_ai_deriv = gsl_sf_airy_zero_ai_deriv(s)
  end function fgsl_sf_airy_zero_ai_deriv
  function fgsl_sf_airy_zero_ai_deriv_e (s, result) 
    integer(fgsl_int), intent(in) :: s
    integer(fgsl_int) :: fgsl_sf_airy_zero_ai_deriv_e
    type(fgsl_sf_result), intent(out) :: result
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_zero_ai_deriv_e = gsl_sf_airy_zero_ai_deriv_e(s, res)
    result = res
  end function fgsl_sf_airy_zero_ai_deriv_e
  function fgsl_sf_airy_zero_bi_deriv(s)
    integer(fgsl_int), intent(in) :: s
    real(fgsl_double) :: fgsl_sf_airy_zero_bi_deriv
    fgsl_sf_airy_zero_bi_deriv = gsl_sf_airy_zero_bi_deriv(s)
  end function fgsl_sf_airy_zero_bi_deriv
  function fgsl_sf_airy_zero_bi_deriv_e (s, result)
    integer(fgsl_int), intent(in) :: s
    integer(fgsl_int) :: fgsl_sf_airy_zero_bi_deriv_e
    type(fgsl_sf_result), intent(out) :: result
!
    type(gsl_sf_result) :: res
    fgsl_sf_airy_zero_bi_deriv_e = gsl_sf_airy_zero_bi_deriv_e(s, res)
    result = res
  end function fgsl_sf_airy_zero_bi_deriv_e
  function fgsl_sf_bessel_jc0(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_jc0
    fgsl_sf_bessel_jc0 = gsl_sf_bessel_jc0(x)
  end function fgsl_sf_bessel_jc0
  function fgsl_sf_bessel_jc0_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_jc0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_jc0_e = gsl_sf_bessel_jc0_e(x, res)
    result = res
  end function fgsl_sf_bessel_jc0_e
  function fgsl_sf_bessel_jc1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_jc1
    fgsl_sf_bessel_jc1 = gsl_sf_bessel_jc1(x)
  end function fgsl_sf_bessel_jc1
  function fgsl_sf_bessel_jc1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_jc1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_jc1_e = gsl_sf_bessel_jc1_e(x, res)
    result = res
  end function fgsl_sf_bessel_jc1_e
  function fgsl_sf_bessel_jcn(n, x)
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_jcn
    fgsl_sf_bessel_jcn = gsl_sf_bessel_jcn(n, x)
  end function fgsl_sf_bessel_jcn
  function fgsl_sf_bessel_jcn_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_jcn_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_jcn_e = gsl_sf_bessel_jcn_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_jcn_e
  function fgsl_sf_bessel_jcn_array(nmin, nmax, x, result) 
    integer(fgsl_int), intent(in) :: nmin, nmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result(:)
    integer(fgsl_int) :: fgsl_sf_bessel_jcn_array
    fgsl_sf_bessel_jcn_array = gsl_sf_bessel_jcn_array(nmin, nmax, x, result)
  end function fgsl_sf_bessel_jcn_array
  function fgsl_sf_bessel_yc0(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_yc0
    fgsl_sf_bessel_yc0 = gsl_sf_bessel_yc0(x)
  end function fgsl_sf_bessel_yc0
  function fgsl_sf_bessel_yc0_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_yc0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_yc0_e = gsl_sf_bessel_yc0_e(x, res)
    result = res
  end function fgsl_sf_bessel_yc0_e
  function fgsl_sf_bessel_yc1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_yc1
    fgsl_sf_bessel_yc1 = gsl_sf_bessel_yc1(x)
  end function fgsl_sf_bessel_yc1
  function fgsl_sf_bessel_yc1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_yc1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_yc1_e = gsl_sf_bessel_yc1_e(x, res)
    result = res
  end function fgsl_sf_bessel_yc1_e
  function fgsl_sf_bessel_ycn(n, x)
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ycn
    fgsl_sf_bessel_ycn = gsl_sf_bessel_ycn(n, x)
  end function fgsl_sf_bessel_ycn
  function fgsl_sf_bessel_ycn_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ycn_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ycn_e = gsl_sf_bessel_ycn_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_ycn_e
  function fgsl_sf_bessel_ycn_array(nmin, nmax, x, result) 
    integer(fgsl_int), intent(in) :: nmin, nmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result(:)
    integer(fgsl_int) :: fgsl_sf_bessel_ycn_array
    fgsl_sf_bessel_ycn_array = gsl_sf_bessel_ycn_array(nmin, nmax, x, result)
  end function fgsl_sf_bessel_ycn_array
  function fgsl_sf_bessel_ic0(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ic0
    fgsl_sf_bessel_ic0 = gsl_sf_bessel_ic0(x)
  end function fgsl_sf_bessel_ic0
  function fgsl_sf_bessel_ic0_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ic0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ic0_e = gsl_sf_bessel_ic0_e(x, res)
    result = res
  end function fgsl_sf_bessel_ic0_e
  function fgsl_sf_bessel_ic1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ic1
    fgsl_sf_bessel_ic1 = gsl_sf_bessel_ic1(x)
  end function fgsl_sf_bessel_ic1
  function fgsl_sf_bessel_ic1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ic1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ic1_e = gsl_sf_bessel_ic1_e(x, res)
    result = res
  end function fgsl_sf_bessel_ic1_e
  function fgsl_sf_bessel_icn(n, x)
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_icn
    fgsl_sf_bessel_icn = gsl_sf_bessel_icn(n, x)
  end function fgsl_sf_bessel_icn
  function fgsl_sf_bessel_icn_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_icn_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_icn_e = gsl_sf_bessel_icn_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_icn_e
  function fgsl_sf_bessel_icn_array(nmin, nmax, x, result) 
    integer(fgsl_int), intent(in) :: nmin, nmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result(:)
    integer(fgsl_int) :: fgsl_sf_bessel_icn_array
    fgsl_sf_bessel_icn_array = gsl_sf_bessel_icn_array(nmin, nmax, x, result)
  end function fgsl_sf_bessel_icn_array
  function fgsl_sf_bessel_ic0_scaled(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ic0_scaled
    fgsl_sf_bessel_ic0_scaled = gsl_sf_bessel_ic0_scaled(x)
  end function fgsl_sf_bessel_ic0_scaled
  function fgsl_sf_bessel_ic0_scaled_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ic0_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ic0_scaled_e = gsl_sf_bessel_ic0_scaled_e(x, res)
    result = res
  end function fgsl_sf_bessel_ic0_scaled_e
  function fgsl_sf_bessel_ic1_scaled(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ic1_scaled
    fgsl_sf_bessel_ic1_scaled = gsl_sf_bessel_ic1_scaled(x)
  end function fgsl_sf_bessel_ic1_scaled
  function fgsl_sf_bessel_ic1_scaled_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ic1_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ic1_scaled_e = gsl_sf_bessel_ic1_scaled_e(x, res)
    result = res
  end function fgsl_sf_bessel_ic1_scaled_e
  function fgsl_sf_bessel_icn_scaled(n, x)
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_icn_scaled
    fgsl_sf_bessel_icn_scaled = gsl_sf_bessel_icn_scaled(n, x)
  end function fgsl_sf_bessel_icn_scaled
  function fgsl_sf_bessel_icn_scaled_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_icn_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_icn_scaled_e = gsl_sf_bessel_icn_scaled_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_icn_scaled_e
  function fgsl_sf_bessel_icn_scaled_array(nmin, nmax, x, result) 
    integer(fgsl_int), intent(in) :: nmin, nmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result(:)
    integer(fgsl_int) :: fgsl_sf_bessel_icn_scaled_array
    fgsl_sf_bessel_icn_scaled_array = gsl_sf_bessel_icn_scaled_array(nmin, nmax, x, result)
  end function fgsl_sf_bessel_icn_scaled_array
  function fgsl_sf_bessel_kc0(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_kc0
    fgsl_sf_bessel_kc0 = gsl_sf_bessel_kc0(x)
  end function fgsl_sf_bessel_kc0
  function fgsl_sf_bessel_kc0_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_kc0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_kc0_e = gsl_sf_bessel_kc0_e(x, res)
    result = res
  end function fgsl_sf_bessel_kc0_e
  function fgsl_sf_bessel_kc1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_kc1
    fgsl_sf_bessel_kc1 = gsl_sf_bessel_kc1(x)
  end function fgsl_sf_bessel_kc1
  function fgsl_sf_bessel_kc1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_kc1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_kc1_e = gsl_sf_bessel_kc1_e(x, res)
    result = res
  end function fgsl_sf_bessel_kc1_e
  function fgsl_sf_bessel_kcn(n, x)
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_kcn
    fgsl_sf_bessel_kcn = gsl_sf_bessel_kcn(n, x)
  end function fgsl_sf_bessel_kcn
  function fgsl_sf_bessel_kcn_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_kcn_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_kcn_e = gsl_sf_bessel_kcn_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_kcn_e
  function fgsl_sf_bessel_kcn_array(nmin, nmax, x, result) 
    integer(fgsl_int), intent(in) :: nmin, nmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result(:)
    integer(fgsl_int) :: fgsl_sf_bessel_kcn_array
    fgsl_sf_bessel_kcn_array = gsl_sf_bessel_kcn_array(nmin, nmax, x, result)
  end function fgsl_sf_bessel_kcn_array
  function fgsl_sf_bessel_kc0_scaled(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_kc0_scaled
    fgsl_sf_bessel_kc0_scaled = gsl_sf_bessel_kc0_scaled(x)
  end function fgsl_sf_bessel_kc0_scaled
  function fgsl_sf_bessel_kc0_scaled_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_kc0_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_kc0_scaled_e = gsl_sf_bessel_kc0_scaled_e(x, res)
    result = res
  end function fgsl_sf_bessel_kc0_scaled_e
  function fgsl_sf_bessel_kc1_scaled(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_kc1_scaled
    fgsl_sf_bessel_kc1_scaled = gsl_sf_bessel_kc1_scaled(x)
  end function fgsl_sf_bessel_kc1_scaled
  function fgsl_sf_bessel_kc1_scaled_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_kc1_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_kc1_scaled_e = gsl_sf_bessel_kc1_scaled_e(x, res)
    result = res
  end function fgsl_sf_bessel_kc1_scaled_e
  function fgsl_sf_bessel_kcn_scaled(n, x)
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_kcn_scaled
    fgsl_sf_bessel_kcn_scaled = gsl_sf_bessel_kcn_scaled(n, x)
  end function fgsl_sf_bessel_kcn_scaled
  function fgsl_sf_bessel_kcn_scaled_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_kcn_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_kcn_scaled_e = gsl_sf_bessel_kcn_scaled_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_kcn_scaled_e
  function fgsl_sf_bessel_kcn_scaled_array(nmin, nmax, x, result) 
    integer(fgsl_int), intent(in) :: nmin, nmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result(:)
    integer(fgsl_int) :: fgsl_sf_bessel_kcn_scaled_array
    fgsl_sf_bessel_kcn_scaled_array = gsl_sf_bessel_kcn_scaled_array(nmin, nmax, x, result)
  end function fgsl_sf_bessel_kcn_scaled_array
! spherical bessel functions
  function fgsl_sf_bessel_js0(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_js0
    fgsl_sf_bessel_js0 = gsl_sf_bessel_js0(x)
  end function fgsl_sf_bessel_js0
  function fgsl_sf_bessel_js0_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_js0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_js0_e = gsl_sf_bessel_js0_e(x, res)
    result = res
  end function fgsl_sf_bessel_js0_e
  function fgsl_sf_bessel_js1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_js1
    fgsl_sf_bessel_js1 = gsl_sf_bessel_js1(x)
  end function fgsl_sf_bessel_js1
  function fgsl_sf_bessel_js1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_js1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_js1_e = gsl_sf_bessel_js1_e(x, res)
    result = res
  end function fgsl_sf_bessel_js1_e
  function fgsl_sf_bessel_js2(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_js2
    fgsl_sf_bessel_js2 = gsl_sf_bessel_js2(x)
  end function fgsl_sf_bessel_js2
  function fgsl_sf_bessel_js2_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_js2_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_js2_e = gsl_sf_bessel_js2_e(x, res)
    result = res
  end function fgsl_sf_bessel_js2_e
  function fgsl_sf_bessel_jsl(n, x)
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_jsl
    fgsl_sf_bessel_jsl = gsl_sf_bessel_jsl(n, x)
  end function fgsl_sf_bessel_jsl
  function fgsl_sf_bessel_jsl_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_jsl_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_jsl_e = gsl_sf_bessel_jsl_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_jsl_e
  function fgsl_sf_bessel_jsl_array(lmax, x, result) 
    integer(fgsl_int), intent(in) :: lmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result(:)
    integer(fgsl_int) :: fgsl_sf_bessel_jsl_array
    fgsl_sf_bessel_jsl_array = gsl_sf_bessel_jsl_array(lmax, x, result)
  end function fgsl_sf_bessel_jsl_array
  function fgsl_sf_bessel_jsl_steed_array(lmax, x, result) 
    integer(fgsl_int), intent(in) :: lmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result(:)
    integer(fgsl_int) :: fgsl_sf_bessel_jsl_steed_array
    fgsl_sf_bessel_jsl_steed_array = gsl_sf_bessel_jsl_steed_array(lmax, x, result)
  end function fgsl_sf_bessel_jsl_steed_array
  function fgsl_sf_bessel_ys0(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ys0
    fgsl_sf_bessel_ys0 = gsl_sf_bessel_ys0(x)
  end function fgsl_sf_bessel_ys0
  function fgsl_sf_bessel_ys0_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ys0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ys0_e = gsl_sf_bessel_ys0_e(x, res)
    result = res
  end function fgsl_sf_bessel_ys0_e
  function fgsl_sf_bessel_ys1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ys1
    fgsl_sf_bessel_ys1 = gsl_sf_bessel_ys1(x)
  end function fgsl_sf_bessel_ys1
  function fgsl_sf_bessel_ys1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ys1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ys1_e = gsl_sf_bessel_ys1_e(x, res)
    result = res
  end function fgsl_sf_bessel_ys1_e
  function fgsl_sf_bessel_ys2(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ys2
    fgsl_sf_bessel_ys2 = gsl_sf_bessel_ys2(x)
  end function fgsl_sf_bessel_ys2
  function fgsl_sf_bessel_ys2_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ys2_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ys2_e = gsl_sf_bessel_ys2_e(x, res)
    result = res
  end function fgsl_sf_bessel_ys2_e
  function fgsl_sf_bessel_ysl(n, x)
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ysl
    fgsl_sf_bessel_ysl = gsl_sf_bessel_ysl(n, x)
  end function fgsl_sf_bessel_ysl
  function fgsl_sf_bessel_ysl_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ysl_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ysl_e = gsl_sf_bessel_ysl_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_ysl_e
  function fgsl_sf_bessel_ysl_array(lmax, x, result) 
    integer(fgsl_int), intent(in) :: lmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result(:)
    integer(fgsl_int) :: fgsl_sf_bessel_ysl_array
    fgsl_sf_bessel_ysl_array = gsl_sf_bessel_ysl_array(lmax, x, result)
  end function fgsl_sf_bessel_ysl_array
  function fgsl_sf_bessel_is0_scaled(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_is0_scaled
    fgsl_sf_bessel_is0_scaled = gsl_sf_bessel_is0_scaled(x)
  end function fgsl_sf_bessel_is0_scaled
  function fgsl_sf_bessel_is0_scaled_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_is0_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_is0_scaled_e = gsl_sf_bessel_is0_scaled_e(x, res)
    result = res
  end function fgsl_sf_bessel_is0_scaled_e
  function fgsl_sf_bessel_is1_scaled(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_is1_scaled
    fgsl_sf_bessel_is1_scaled = gsl_sf_bessel_is1_scaled(x)
  end function fgsl_sf_bessel_is1_scaled
  function fgsl_sf_bessel_is1_scaled_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_is1_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_is1_scaled_e = gsl_sf_bessel_is1_scaled_e(x, res)
    result = res
  end function fgsl_sf_bessel_is1_scaled_e
  function fgsl_sf_bessel_is2_scaled(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_is2_scaled
    fgsl_sf_bessel_is2_scaled = gsl_sf_bessel_is2_scaled(x)
  end function fgsl_sf_bessel_is2_scaled
  function fgsl_sf_bessel_is2_scaled_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_is2_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_is2_scaled_e = gsl_sf_bessel_is2_scaled_e(x, res)
    result = res
  end function fgsl_sf_bessel_is2_scaled_e
  function fgsl_sf_bessel_isl_scaled(n, x)
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_isl_scaled
    fgsl_sf_bessel_isl_scaled = gsl_sf_bessel_isl_scaled(n, x)
  end function fgsl_sf_bessel_isl_scaled
  function fgsl_sf_bessel_isl_scaled_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_isl_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_isl_scaled_e = gsl_sf_bessel_isl_scaled_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_isl_scaled_e
  function fgsl_sf_bessel_isl_scaled_array(lmax, x, result) 
    integer(fgsl_int), intent(in) :: lmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result(:)
    integer(fgsl_int) :: fgsl_sf_bessel_isl_scaled_array
    fgsl_sf_bessel_isl_scaled_array = gsl_sf_bessel_isl_scaled_array(lmax, x, result)
  end function fgsl_sf_bessel_isl_scaled_array
  function fgsl_sf_bessel_ks0_scaled(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ks0_scaled
    fgsl_sf_bessel_ks0_scaled = gsl_sf_bessel_ks0_scaled(x)
  end function fgsl_sf_bessel_ks0_scaled
  function fgsl_sf_bessel_ks0_scaled_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ks0_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ks0_scaled_e = gsl_sf_bessel_ks0_scaled_e(x, res)
    result = res
  end function fgsl_sf_bessel_ks0_scaled_e
  function fgsl_sf_bessel_ks1_scaled(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ks1_scaled
    fgsl_sf_bessel_ks1_scaled = gsl_sf_bessel_ks1_scaled(x)
  end function fgsl_sf_bessel_ks1_scaled
  function fgsl_sf_bessel_ks1_scaled_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ks1_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ks1_scaled_e = gsl_sf_bessel_ks1_scaled_e(x, res)
    result = res
  end function fgsl_sf_bessel_ks1_scaled_e
  function fgsl_sf_bessel_ks2_scaled(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ks2_scaled
    fgsl_sf_bessel_ks2_scaled = gsl_sf_bessel_ks2_scaled(x)
  end function fgsl_sf_bessel_ks2_scaled
  function fgsl_sf_bessel_ks2_scaled_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ks2_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ks2_scaled_e = gsl_sf_bessel_ks2_scaled_e(x, res)
    result = res
  end function fgsl_sf_bessel_ks2_scaled_e
  function fgsl_sf_bessel_ksl_scaled(n, x)
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ksl_scaled
    fgsl_sf_bessel_ksl_scaled = gsl_sf_bessel_ksl_scaled(n, x)
  end function fgsl_sf_bessel_ksl_scaled
  function fgsl_sf_bessel_ksl_scaled_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ksl_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ksl_scaled_e = gsl_sf_bessel_ksl_scaled_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_ksl_scaled_e
  function fgsl_sf_bessel_ksl_scaled_array(lmax, x, result) 
    integer(fgsl_int), intent(in) :: lmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result(:)
    integer(fgsl_int) :: fgsl_sf_bessel_ksl_scaled_array
    fgsl_sf_bessel_ksl_scaled_array = gsl_sf_bessel_ksl_scaled_array(lmax, x, result)
  end function fgsl_sf_bessel_ksl_scaled_array
! fractional order bessel functions
  function fgsl_sf_bessel_jnu(n, x)
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_jnu
    fgsl_sf_bessel_jnu = gsl_sf_bessel_jnu(n, x)
  end function fgsl_sf_bessel_jnu
  function fgsl_sf_bessel_jnu_e(n, x, result) 
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_jnu_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_jnu_e = gsl_sf_bessel_jnu_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_jnu_e
  function fgsl_sf_bessel_sequence_jnu_e(nu, mode, size, v)
    real(fgsl_double), intent(in) :: nu
    type(fgsl_mode_t), intent(in) :: mode
    integer(fgsl_size_t), intent(in) :: size
    real(fgsl_double), intent(inout) :: v(:)
    integer(fgsl_int) :: fgsl_sf_bessel_sequence_jnu_e
    fgsl_sf_bessel_sequence_jnu_e = gsl_sf_bessel_sequence_jnu_e(nu, mode%gsl_mode, size, v)
  end function fgsl_sf_bessel_sequence_jnu_e
  function fgsl_sf_bessel_ynu(n, x)
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_ynu
    fgsl_sf_bessel_ynu = gsl_sf_bessel_ynu(n, x)
  end function fgsl_sf_bessel_ynu
  function fgsl_sf_bessel_ynu_e(n, x, result) 
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_ynu_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_ynu_e = gsl_sf_bessel_ynu_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_ynu_e
  function fgsl_sf_bessel_inu(n, x)
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_inu
    fgsl_sf_bessel_inu = gsl_sf_bessel_inu(n, x)
  end function fgsl_sf_bessel_inu
  function fgsl_sf_bessel_inu_e(n, x, result) 
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_inu_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_inu_e = gsl_sf_bessel_inu_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_inu_e
  function fgsl_sf_bessel_inu_scaled(n, x)
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_inu_scaled
    fgsl_sf_bessel_inu_scaled = gsl_sf_bessel_inu_scaled(n, x)
  end function fgsl_sf_bessel_inu_scaled
  function fgsl_sf_bessel_inu_scaled_e(n, x, result) 
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_inu_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_inu_scaled_e = gsl_sf_bessel_inu_scaled_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_inu_scaled_e
  function fgsl_sf_bessel_knu(n, x)
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_knu
    fgsl_sf_bessel_knu = gsl_sf_bessel_knu(n, x)
  end function fgsl_sf_bessel_knu
  function fgsl_sf_bessel_knu_e(n, x, result) 
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_knu_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_knu_e = gsl_sf_bessel_knu_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_knu_e
  function fgsl_sf_bessel_lnknu(n, x)
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_lnknu
    fgsl_sf_bessel_lnknu = gsl_sf_bessel_lnknu(n, x)
  end function fgsl_sf_bessel_lnknu
  function fgsl_sf_bessel_lnknu_e(n, x, result) 
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_lnknu_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_lnknu_e = gsl_sf_bessel_lnknu_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_lnknu_e
  function fgsl_sf_bessel_knu_scaled(n, x)
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_bessel_knu_scaled
    fgsl_sf_bessel_knu_scaled = gsl_sf_bessel_knu_scaled(n, x)
  end function fgsl_sf_bessel_knu_scaled
  function fgsl_sf_bessel_knu_scaled_e(n, x, result) 
    real(fgsl_double), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_knu_scaled_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_knu_scaled_e = gsl_sf_bessel_knu_scaled_e(n, x, res)
    result = res
  end function fgsl_sf_bessel_knu_scaled_e
  function fgsl_sf_bessel_zero_jc0(s) 
    integer(fgsl_int), intent(in) :: s
    real(fgsl_double) :: fgsl_sf_bessel_zero_jc0
    fgsl_sf_bessel_zero_jc0 = gsl_sf_bessel_zero_jc0(s) 
  end function fgsl_sf_bessel_zero_jc0
  function fgsl_sf_bessel_zero_jc0_e(s, result) 
    integer(fgsl_int), intent(in) :: s
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_zero_jc0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_zero_jc0_e = gsl_sf_bessel_zero_jc0_e(s, res) 
    result = res
  end function fgsl_sf_bessel_zero_jc0_e
  function fgsl_sf_bessel_zero_jc1(s) 
    integer(fgsl_int), intent(in) :: s
    real(fgsl_double) :: fgsl_sf_bessel_zero_jc1
    fgsl_sf_bessel_zero_jc1 = gsl_sf_bessel_zero_jc1(s) 
  end function fgsl_sf_bessel_zero_jc1
  function fgsl_sf_bessel_zero_jc1_e(s, result) 
    integer(fgsl_int), intent(in) :: s
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_zero_jc1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_zero_jc1_e = gsl_sf_bessel_zero_jc1_e(s, res) 
    result = res
  end function fgsl_sf_bessel_zero_jc1_e
  function fgsl_sf_bessel_zero_jnu(nu, s) 
    real(fgsl_double), intent(in) :: nu
    integer(fgsl_int), intent(in) :: s
    real(fgsl_double) :: fgsl_sf_bessel_zero_jnu
    fgsl_sf_bessel_zero_jnu = gsl_sf_bessel_zero_jnu(nu, s) 
  end function fgsl_sf_bessel_zero_jnu
  function fgsl_sf_bessel_zero_jnu_e(nu, s, result) 
    real(fgsl_double), intent(in) :: nu
    integer(fgsl_int), intent(in) :: s
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_bessel_zero_jnu_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_bessel_zero_jnu_e = gsl_sf_bessel_zero_jnu_e(nu, s, res) 
    result = res
  end function fgsl_sf_bessel_zero_jnu_e
!
  function fgsl_sf_clausen(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_clausen
    fgsl_sf_clausen = gsl_sf_clausen(x) 
  end function fgsl_sf_clausen
  function fgsl_sf_clausen_e(x, result)
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_clausen_e
!
    type(gsl_sf_result) :: res 
    fgsl_sf_clausen_e = gsl_sf_clausen_e(x, res)
    result = res
  end function fgsl_sf_clausen_e
  function fgsl_sf_hydrogenicr_1(z, r)
    real(fgsl_double), intent(in) :: z, r
    real(fgsl_double) :: fgsl_sf_hydrogenicr_1
    fgsl_sf_hydrogenicr_1 = gsl_sf_hydrogenicr_1(z, r)
  end function fgsl_sf_hydrogenicr_1
  function fgsl_sf_hydrogenicr_1_e(z, r, result)
    real(fgsl_double), intent(in) :: z, r
    integer(fgsl_int) :: fgsl_sf_hydrogenicr_1_e
    type(fgsl_sf_result), intent(out) :: result
!
    type(gsl_sf_result) :: res
    fgsl_sf_hydrogenicr_1_e = gsl_sf_hydrogenicr_1_e(z, r, res)
    result = res
  end function fgsl_sf_hydrogenicr_1_e
  function fgsl_sf_hydrogenicr(n, l, z, r)
    integer(fgsl_int), intent(in) :: n, l
    real(fgsl_double), intent(in) :: z, r
    real(fgsl_double) :: fgsl_sf_hydrogenicr
    fgsl_sf_hydrogenicr = gsl_sf_hydrogenicr(n, l, z, r)
  end function fgsl_sf_hydrogenicr
  function fgsl_sf_hydrogenicr_e(n, l, z, r, result)
    integer(fgsl_int), intent(in) :: n, l
    real(fgsl_double), intent(in) :: z, r
    integer(fgsl_int) :: fgsl_sf_hydrogenicr_e
    type(fgsl_sf_result), intent(out) :: result
!
    type(gsl_sf_result) :: res
    fgsl_sf_hydrogenicr_e = gsl_sf_hydrogenicr_e(n, l, z, r, res)
    result = res
  end function fgsl_sf_hydrogenicr_e
  function fgsl_sf_coulomb_wave_fg_e(eta, x, l_f, k, f, fp, g, gp, exp_f, exp_g) 
    real(fgsl_double), intent(in) :: eta, x, l_f
    integer(fgsl_int), intent(in) :: k
    type(fgsl_sf_result), intent(out) :: f, fp, g, gp
    real(fgsl_double), intent(out) :: exp_f, exp_g
    integer(fgsl_int) :: fgsl_sf_coulomb_wave_fg_e
!
    type(gsl_sf_result) :: fl, fpl, gl, gpl
    fgsl_sf_coulomb_wave_fg_e = gsl_sf_coulomb_wave_fg_e(eta, x, l_f, k, &
         fl, fpl, gl, gpl, exp_f, exp_g)
    f = fl
    fp = fpl
    g = gl
    gp = gpl
  end function fgsl_sf_coulomb_wave_fg_e
  function fgsl_sf_coulomb_wave_f_array (l_min, kmax, eta, x, fc_array, &
          f_exponent) 
    real(fgsl_double), intent(in) :: l_min, eta, x
    integer(fgsl_int), intent(in) :: kmax
    real(fgsl_double), intent(out) :: fc_array(:)
    real(fgsl_double), intent(out) :: f_exponent
    integer(fgsl_int) :: fgsl_sf_coulomb_wave_f_array
    fgsl_sf_coulomb_wave_f_array = gsl_sf_coulomb_wave_f_array (l_min, kmax, &
         eta, x, fc_array, f_exponent)
  end function fgsl_sf_coulomb_wave_f_array
  function fgsl_sf_coulomb_wave_fg_array(l_min, kmax, eta, x, fc_array, &
          gc_array, f_exponent, g_exponent) 
    real(fgsl_double), intent(in) :: l_min, eta, x
    integer(fgsl_int), intent(in) :: kmax
    real(fgsl_double), intent(out) :: fc_array(:), gc_array(:)
    real(fgsl_double), intent(out) :: f_exponent, g_exponent
    integer(fgsl_int) :: fgsl_sf_coulomb_wave_fg_array
    fgsl_sf_coulomb_wave_fg_array = gsl_sf_coulomb_wave_fg_array (l_min, kmax, &
         eta, x, fc_array, gc_array, f_exponent, g_exponent)
  end function fgsl_sf_coulomb_wave_fg_array
  function fgsl_sf_coulomb_wave_fgp_array(l_min, kmax, eta, x, fc_array, fcp_array, &
          gc_array, gcp_array, f_exponent, g_exponent) 
    real(fgsl_double), intent(in) :: l_min, eta, x
    integer(fgsl_int), intent(in) :: kmax
    real(fgsl_double), intent(out) :: fc_array(:), gc_array(:), fcp_array(:), gcp_array(:)
    real(fgsl_double), intent(out) :: f_exponent, g_exponent
    integer(fgsl_int) :: fgsl_sf_coulomb_wave_fgp_array
    fgsl_sf_coulomb_wave_fgp_array = gsl_sf_coulomb_wave_fgp_array (l_min, kmax, &
         eta, x, fc_array, fcp_array, gc_array, gcp_array, f_exponent, g_exponent)
  end function fgsl_sf_coulomb_wave_fgp_array
  function fgsl_sf_coulomb_wave_sphf_array(l_min, kmax, eta, x, fc_array, &
          f_exponent) 
    real(fgsl_double), intent(in) :: l_min, eta, x
    integer(fgsl_int), intent(in) :: kmax
    real(fgsl_double), intent(out) :: fc_array(:)
    real(fgsl_double), intent(out) :: f_exponent
    integer(fgsl_int) :: fgsl_sf_coulomb_wave_sphf_array
    fgsl_sf_coulomb_wave_sphf_array = gsl_sf_coulomb_wave_sphf_array (l_min, kmax, &
         eta, x, fc_array, f_exponent)
  end function fgsl_sf_coulomb_wave_sphf_array
  function fgsl_sf_coulomb_cl_e(l, eta, result) 
    real(fgsl_double), intent(in) :: eta, l
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_coulomb_cl_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_coulomb_cl_e = gsl_sf_coulomb_cl_e(l, eta, res)
    result = res
  end function fgsl_sf_coulomb_cl_e
  function fgsl_sf_coulomb_cl_array(l_min, kmax, eta, cl) 
    real(fgsl_double), intent(in) :: l_min, eta
    integer(fgsl_int), intent(in) :: kmax
    real(fgsl_double), intent(out) :: cl(:)
    integer(fgsl_int) :: fgsl_sf_coulomb_cl_array
    fgsl_sf_coulomb_cl_array = gsl_sf_coulomb_cl_array (l_min, kmax, eta, cl)
  end function fgsl_sf_coulomb_cl_array
  function fgsl_sf_coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc)
    integer(fgsl_int), intent(in) :: two_ja, two_jb, two_jc, two_ma, two_mb, two_mc
    real(fgsl_double) :: fgsl_sf_coupling_3j
    fgsl_sf_coupling_3j = gsl_sf_coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc)
  end function fgsl_sf_coupling_3j
  function fgsl_sf_coupling_3j_e(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc, result)
    integer(fgsl_int), intent(in) :: two_ja, two_jb, two_jc, two_ma, two_mb, two_mc
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_coupling_3j_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_coupling_3j_e = gsl_sf_coupling_3j_e(two_ja, two_jb, two_jc, two_ma, &
         two_mb, two_mc, res)
    result = res
  end function fgsl_sf_coupling_3j_e
  function fgsl_sf_coupling_6j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf)
    integer(fgsl_int), intent(in) :: two_ja, two_jb, two_jc, two_jd, two_je, two_jf
    real(fgsl_double) :: fgsl_sf_coupling_6j
    fgsl_sf_coupling_6j = gsl_sf_coupling_6j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf)
  end function fgsl_sf_coupling_6j
  function fgsl_sf_coupling_6j_e(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, result)
    integer(fgsl_int), intent(in) :: two_ja, two_jb, two_jc, two_jd, two_je, two_jf
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_coupling_6j_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_coupling_6j_e = gsl_sf_coupling_6j_e(two_ja, two_jb, two_jc, two_jd, &
         two_je, two_jf, res)
    result = res
  end function fgsl_sf_coupling_6j_e
  function fgsl_sf_coupling_9j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, &
          two_jg, two_jh, two_ji)
    integer(fgsl_int), intent(in) :: two_ja, two_jb, two_jc, two_jd, two_je, two_jf, &
          two_jg, two_jh, two_ji
    real(fgsl_double) :: fgsl_sf_coupling_9j
    fgsl_sf_coupling_9j = gsl_sf_coupling_9j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, &
          two_jg, two_jh, two_ji)
  end function fgsl_sf_coupling_9j
  function fgsl_sf_coupling_9j_e(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, &
          two_jg, two_jh, two_ji, result)
    integer(fgsl_int), intent(in) :: two_ja, two_jb, two_jc, two_jd, two_je, two_jf, &
          two_jg, two_jh, two_ji
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_coupling_9j_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_coupling_9j_e = gsl_sf_coupling_9j_e(two_ja, two_jb, two_jc, two_jd, &
         two_je, two_jf, two_jg, two_jh, two_ji, res)
    result = res
  end function fgsl_sf_coupling_9j_e
  function fgsl_sf_dawson(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_dawson
    fgsl_sf_dawson = gsl_sf_dawson(x)
  end function fgsl_sf_dawson
  function fgsl_sf_dawson_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_dawson_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_dawson_e = gsl_sf_dawson_e(x, res)
    result = res
  end function fgsl_sf_dawson_e
  function fgsl_sf_debye_1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_debye_1
    fgsl_sf_debye_1 = gsl_sf_debye_1(x)
  end function fgsl_sf_debye_1
  function fgsl_sf_debye_1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_debye_1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_debye_1_e = gsl_sf_debye_1_e(x, res)
    result = res
  end function fgsl_sf_debye_1_e
  function fgsl_sf_debye_2(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_debye_2
    fgsl_sf_debye_2 = gsl_sf_debye_2(x)
  end function fgsl_sf_debye_2
  function fgsl_sf_debye_2_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_debye_2_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_debye_2_e = gsl_sf_debye_2_e(x, res)
    result = res
  end function fgsl_sf_debye_2_e
  function fgsl_sf_debye_3(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_debye_3
    fgsl_sf_debye_3 = gsl_sf_debye_3(x)
  end function fgsl_sf_debye_3
  function fgsl_sf_debye_3_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_debye_3_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_debye_3_e = gsl_sf_debye_3_e(x, res)
    result = res
  end function fgsl_sf_debye_3_e
  function fgsl_sf_debye_4(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_debye_4
    fgsl_sf_debye_4 = gsl_sf_debye_4(x)
  end function fgsl_sf_debye_4
  function fgsl_sf_debye_4_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_debye_4_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_debye_4_e = gsl_sf_debye_4_e(x, res)
    result = res
  end function fgsl_sf_debye_4_e
  function fgsl_sf_debye_5(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_debye_5
    fgsl_sf_debye_5 = gsl_sf_debye_5(x)
  end function fgsl_sf_debye_5
  function fgsl_sf_debye_5_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_debye_5_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_debye_5_e = gsl_sf_debye_5_e(x, res)
    result = res
  end function fgsl_sf_debye_5_e
  function fgsl_sf_debye_6(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_debye_6
    fgsl_sf_debye_6 = gsl_sf_debye_6(x)
  end function fgsl_sf_debye_6
  function fgsl_sf_debye_6_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_debye_6_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_debye_6_e = gsl_sf_debye_6_e(x, res)
    result = res
  end function fgsl_sf_debye_6_e
  function fgsl_sf_dilog(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_dilog
    fgsl_sf_dilog = gsl_sf_dilog(x)
  end function fgsl_sf_dilog
  function fgsl_sf_dilog_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_dilog_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_dilog_e = gsl_sf_dilog_e(x, res)
    result = res
  end function fgsl_sf_dilog_e
  function fgsl_sf_complex_dilog_e(r, theta, result_re, result_im)
    real(fgsl_double), intent(in) :: r, theta
    type(fgsl_sf_result), intent(out) :: result_re, result_im
    integer(fgsl_int) :: fgsl_sf_complex_dilog_e
!
    type(gsl_sf_result) :: res_re, res_im
    fgsl_sf_complex_dilog_e = gsl_sf_complex_dilog_e(r, theta, res_re, res_im)
    result_re = res_re
    result_im = res_im
  end function fgsl_sf_complex_dilog_e
  function fgsl_sf_multiply_e(x, y, result) 
    real(fgsl_double), intent(in) :: x, y
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_multiply_e     
!
    type(gsl_sf_result) :: res
    fgsl_sf_multiply_e = gsl_sf_multiply_e(x, y, res)
    result = res
  end function fgsl_sf_multiply_e
  function fgsl_sf_multiply_err_e(x, dx, y, dy, result) 
    real(fgsl_double), intent(in) :: x, y, dx, dy
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_multiply_err_e     
!
    type(gsl_sf_result) :: res
    fgsl_sf_multiply_err_e = gsl_sf_multiply_err_e(x, dx, y, dy, res)
    result = res
  end function fgsl_sf_multiply_err_e
  function fgsl_sf_ellint_kcomp(k, mode)
    real(fgsl_double), intent(in) :: k
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_ellint_kcomp
    fgsl_sf_ellint_kcomp = gsl_sf_ellint_kcomp(k, mode%gsl_mode)
  end function fgsl_sf_ellint_kcomp
  function fgsl_sf_ellint_kcomp_e(k, mode, result)
    real(fgsl_double), intent(in) :: k
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ellint_kcomp_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ellint_kcomp_e = gsl_sf_ellint_kcomp_e(k, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_ellint_kcomp_e
  function fgsl_sf_ellint_ecomp(k, mode)
    real(fgsl_double), intent(in) :: k
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_ellint_ecomp
    fgsl_sf_ellint_ecomp = gsl_sf_ellint_ecomp(k, mode%gsl_mode)
  end function fgsl_sf_ellint_ecomp
  function fgsl_sf_ellint_ecomp_e(k, mode, result)
    real(fgsl_double), intent(in) :: k
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ellint_ecomp_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ellint_ecomp_e = gsl_sf_ellint_ecomp_e(k, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_ellint_ecomp_e
  function fgsl_sf_ellint_pcomp(k, n, mode)
    real(fgsl_double), intent(in) :: k, n
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_ellint_pcomp
    fgsl_sf_ellint_pcomp = gsl_sf_ellint_pcomp(k, n, mode%gsl_mode)
  end function fgsl_sf_ellint_pcomp
  function fgsl_sf_ellint_pcomp_e(k, n, mode, result)
    real(fgsl_double), intent(in) :: k, n
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ellint_pcomp_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ellint_pcomp_e = gsl_sf_ellint_pcomp_e(k, n, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_ellint_pcomp_e
  function fgsl_sf_ellint_f(phi, k, mode)
    real(fgsl_double), intent(in) :: phi, k
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_ellint_f
    fgsl_sf_ellint_f = gsl_sf_ellint_f(phi, k, mode%gsl_mode)
  end function fgsl_sf_ellint_f
  function fgsl_sf_ellint_f_e(phi, k, mode, result)
    real(fgsl_double), intent(in) :: phi, k
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ellint_f_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ellint_f_e = gsl_sf_ellint_f_e(phi, k, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_ellint_f_e
  function fgsl_sf_ellint_e(phi, k, mode)
    real(fgsl_double), intent(in) :: phi, k
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_ellint_e
    fgsl_sf_ellint_e = gsl_sf_ellint_e(phi, k, mode%gsl_mode)
  end function fgsl_sf_ellint_e
  function fgsl_sf_ellint_e_e(phi, k, mode, result)
    real(fgsl_double), intent(in) :: phi, k
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ellint_e_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ellint_e_e = gsl_sf_ellint_e_e(phi, k, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_ellint_e_e
  function fgsl_sf_ellint_p(phi, k, n, mode)
    real(fgsl_double), intent(in) :: phi, k, n
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_ellint_p
    fgsl_sf_ellint_p = gsl_sf_ellint_p(phi, k, n, mode%gsl_mode)
  end function fgsl_sf_ellint_p
  function fgsl_sf_ellint_p_e(phi, k, n, mode, result)
    real(fgsl_double), intent(in) :: phi, k, n
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ellint_p_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ellint_p_e = gsl_sf_ellint_p_e(phi, k, n, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_ellint_p_e
  function fgsl_sf_ellint_d(phi, k, n, mode)
    real(fgsl_double), intent(in) :: phi, k, n
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_ellint_d
    fgsl_sf_ellint_d = gsl_sf_ellint_d(phi, k, n, mode%gsl_mode)
  end function fgsl_sf_ellint_d
! NOTE: Parameter n appears to be superfluous in ellint_d
! it is not referenced in the implementation
  function fgsl_sf_ellint_d_e(phi, k, n, mode, result)
    real(fgsl_double), intent(in) :: phi, k, n
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ellint_d_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ellint_d_e = gsl_sf_ellint_d_e(phi, k, n, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_ellint_d_e
  function fgsl_sf_ellint_rc(x, y, mode)
    real(fgsl_double), intent(in) :: x, y
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_ellint_rc
    fgsl_sf_ellint_rc = gsl_sf_ellint_rc(x, y, mode%gsl_mode)
  end function fgsl_sf_ellint_rc
  function fgsl_sf_ellint_rc_e(x, y, mode, result)
    real(fgsl_double), intent(in) :: x, y
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ellint_rc_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ellint_rc_e = gsl_sf_ellint_rc_e(x, y, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_ellint_rc_e
  function fgsl_sf_ellint_rd(x, y, z, mode)
    real(fgsl_double), intent(in) :: x, y, z
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_ellint_rd
    fgsl_sf_ellint_rd = gsl_sf_ellint_rd(x, y, z, mode%gsl_mode)
  end function fgsl_sf_ellint_rd
  function fgsl_sf_ellint_rd_e(x, y, z, mode, result)
    real(fgsl_double), intent(in) :: x, y, z
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ellint_rd_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ellint_rd_e = gsl_sf_ellint_rd_e(x, y, z, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_ellint_rd_e
  function fgsl_sf_ellint_rf(x, y, z, mode)
    real(fgsl_double), intent(in) :: x, y, z
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_ellint_rf
    fgsl_sf_ellint_rf = gsl_sf_ellint_rf(x, y, z, mode%gsl_mode)
  end function fgsl_sf_ellint_rf
  function fgsl_sf_ellint_rf_e(x, y, z, mode, result)
    real(fgsl_double), intent(in) :: x, y, z
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ellint_rf_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ellint_rf_e = gsl_sf_ellint_rf_e(x, y, z, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_ellint_rf_e
  function fgsl_sf_ellint_rj(x, y, z, p, mode)
    real(fgsl_double), intent(in) :: x, y, z, p
    type(fgsl_mode_t), intent(in) :: mode
    real(fgsl_double) :: fgsl_sf_ellint_rj
    fgsl_sf_ellint_rj = gsl_sf_ellint_rj(x, y, z, p, mode%gsl_mode)
  end function fgsl_sf_ellint_rj
  function fgsl_sf_ellint_rj_e(x, y, z, p, mode, result)
    real(fgsl_double), intent(in) :: x, y, z, p
    type(fgsl_mode_t), intent(in) :: mode
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ellint_rj_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ellint_rj_e = gsl_sf_ellint_rj_e(x, y, z, p, mode%gsl_mode, res)
    result = res
  end function fgsl_sf_ellint_rj_e
  function fgsl_sf_elljac_e(u, m, sn, cn, dn)
    real(fgsl_double), intent(in) :: u, m
    real(fgsl_double), intent(out) :: sn, cn, dn
    integer(fgsl_int) :: fgsl_sf_elljac_e
    fgsl_sf_elljac_e = gsl_sf_elljac_e(u, m, sn, cn, dn)
  end function fgsl_sf_elljac_e
  function fgsl_sf_erf(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_erf
    fgsl_sf_erf = gsl_sf_erf(x)
  end function fgsl_sf_erf
  function fgsl_sf_erf_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_erf_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_erf_e = gsl_sf_erf_e(x, res)
    result = res
  end function fgsl_sf_erf_e
  function fgsl_sf_erfc(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_erfc
    fgsl_sf_erfc = gsl_sf_erfc(x)
  end function fgsl_sf_erfc
  function fgsl_sf_erfc_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_erfc_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_erfc_e = gsl_sf_erfc_e(x, res)
    result = res
  end function fgsl_sf_erfc_e
  function fgsl_sf_log_erfc(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_log_erfc
    fgsl_sf_log_erfc = gsl_sf_log_erfc(x)
  end function fgsl_sf_log_erfc
  function fgsl_sf_log_erfc_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_log_erfc_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_log_erfc_e = gsl_sf_log_erfc_e(x, res)
    result = res
  end function fgsl_sf_log_erfc_e
  function fgsl_sf_erf_z(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_erf_z
    fgsl_sf_erf_z = gsl_sf_erf_z(x)
  end function fgsl_sf_erf_z
  function fgsl_sf_erf_z_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_erf_z_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_erf_z_e = gsl_sf_erf_z_e(x, res)
    result = res
  end function fgsl_sf_erf_z_e
  function fgsl_sf_erf_q(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_erf_q
    fgsl_sf_erf_q = gsl_sf_erf_q(x)
  end function fgsl_sf_erf_q
  function fgsl_sf_erf_q_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_erf_q_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_erf_q_e = gsl_sf_erf_q_e(x, res)
    result = res
  end function fgsl_sf_erf_q_e
  function fgsl_sf_hazard(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_hazard
    fgsl_sf_hazard = gsl_sf_hazard(x)
  end function fgsl_sf_hazard
  function fgsl_sf_hazard_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hazard_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hazard_e = gsl_sf_hazard_e(x, res)
    result = res
  end function fgsl_sf_hazard_e
  function fgsl_sf_exp(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_exp
    fgsl_sf_exp = gsl_sf_exp(x)
  end function fgsl_sf_exp
  function fgsl_sf_exp_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_exp_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_exp_e = gsl_sf_exp_e(x, res)
    result = res
  end function fgsl_sf_exp_e
  function fgsl_sf_exp_e10_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result_e10), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_exp_e10_e
!
    type(gsl_sf_result_e10) :: res
    fgsl_sf_exp_e10_e = gsl_sf_exp_e10_e(x, res)
    result = res
  end function fgsl_sf_exp_e10_e
  function fgsl_sf_exp_mult(x, y) 
    real(fgsl_double), intent(in) :: x, y
    real(fgsl_double) :: fgsl_sf_exp_mult
    fgsl_sf_exp_mult = gsl_sf_exp_mult(x, y)
  end function fgsl_sf_exp_mult
  function fgsl_sf_exp_mult_e(x, y, result) 
    real(fgsl_double), intent(in) :: x, y
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_exp_mult_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_exp_mult_e = gsl_sf_exp_mult_e(x, y, res)
    result = res
  end function fgsl_sf_exp_mult_e
  function fgsl_sf_exp_mult_e10_e(x, y, result) 
    real(fgsl_double), intent(in) :: x, y
    type(fgsl_sf_result_e10), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_exp_mult_e10_e
!
    type(gsl_sf_result_e10) :: res
    fgsl_sf_exp_mult_e10_e = gsl_sf_exp_mult_e10_e(x, y, res)
    result = res
  end function fgsl_sf_exp_mult_e10_e
  function fgsl_sf_expm1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_expm1
    fgsl_sf_expm1 = gsl_sf_expm1(x)
  end function fgsl_sf_expm1
  function fgsl_sf_expm1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_expm1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_expm1_e = gsl_sf_expm1_e(x, res)
    result = res
  end function fgsl_sf_expm1_e
  function fgsl_sf_exprel(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_exprel
    fgsl_sf_exprel = gsl_sf_exprel(x)
  end function fgsl_sf_exprel
  function fgsl_sf_exprel_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_exprel_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_exprel_e = gsl_sf_exprel_e(x, res)
    result = res
  end function fgsl_sf_exprel_e
  function fgsl_sf_exprel_2(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_exprel_2
    fgsl_sf_exprel_2 = gsl_sf_exprel_2(x)
  end function fgsl_sf_exprel_2
  function fgsl_sf_exprel_2_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_exprel_2_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_exprel_2_e = gsl_sf_exprel_2_e(x, res)
    result = res
  end function fgsl_sf_exprel_2_e
  function fgsl_sf_exprel_n(n, x) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_exprel_n
    fgsl_sf_exprel_n = gsl_sf_exprel_n(n, x)
  end function fgsl_sf_exprel_n
  function fgsl_sf_exprel_n_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_exprel_n_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_exprel_n_e = gsl_sf_exprel_n_e(n, x, res)
    result = res
  end function fgsl_sf_exprel_n_e
  function fgsl_sf_exp_err_e(x, dx, result) 
    real(fgsl_double), intent(in) :: x, dx
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_exp_err_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_exp_err_e = gsl_sf_exp_err_e(x, dx, res)
    result = res
  end function fgsl_sf_exp_err_e
  function fgsl_sf_exp_err_e10_e(x, dx, result) 
    real(fgsl_double), intent(in) :: x, dx
    type(fgsl_sf_result_e10), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_exp_err_e10_e
!
    type(gsl_sf_result_e10) :: res
    fgsl_sf_exp_err_e10_e = gsl_sf_exp_err_e10_e(x, dx, res)
    result = res
  end function fgsl_sf_exp_err_e10_e
  function fgsl_sf_exp_mult_err_e(x, dx, y, dy, result) 
    real(fgsl_double), intent(in) :: x, dx, y, dy
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_exp_mult_err_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_exp_mult_err_e = gsl_sf_exp_mult_err_e(x, dx, y, dy, res)
    result = res
  end function fgsl_sf_exp_mult_err_e
  function fgsl_sf_exp_mult_err_e10_e(x, dx, y, dy, result) 
    real(fgsl_double), intent(in) :: x, dx, y, dy
    type(fgsl_sf_result_e10), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_exp_mult_err_e10_e
!
    type(gsl_sf_result_e10) :: res
    fgsl_sf_exp_mult_err_e10_e = gsl_sf_exp_mult_err_e10_e(x, dx, y, dy, res)
    result = res
  end function fgsl_sf_exp_mult_err_e10_e
  function fgsl_sf_expint_e1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_expint_e1
    fgsl_sf_expint_e1 = gsl_sf_expint_e1(x)
  end function fgsl_sf_expint_e1
  function fgsl_sf_expint_e1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_expint_e1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_expint_e1_e = gsl_sf_expint_e1_e(x, res)
    result = res
  end function fgsl_sf_expint_e1_e
  function fgsl_sf_expint_e2(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_expint_e2
    fgsl_sf_expint_e2 = gsl_sf_expint_e2(x)
  end function fgsl_sf_expint_e2
  function fgsl_sf_expint_e2_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_expint_e2_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_expint_e2_e = gsl_sf_expint_e2_e(x, res)
    result = res
  end function fgsl_sf_expint_e2_e
  function fgsl_sf_expint_en(n, x) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_expint_en
    fgsl_sf_expint_en = gsl_sf_expint_en(n, x)
  end function fgsl_sf_expint_en
  function fgsl_sf_expint_en_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_expint_en_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_expint_en_e = gsl_sf_expint_en_e(n, x, res)
    result = res
  end function fgsl_sf_expint_en_e
  function fgsl_sf_expint_ei(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_expint_ei
    fgsl_sf_expint_ei = gsl_sf_expint_ei(x)
  end function fgsl_sf_expint_ei
  function fgsl_sf_expint_ei_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_expint_ei_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_expint_ei_e = gsl_sf_expint_ei_e(x, res)
    result = res
  end function fgsl_sf_expint_ei_e
  function fgsl_sf_shi(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_shi
    fgsl_sf_shi = gsl_sf_shi(x)
  end function fgsl_sf_shi
  function fgsl_sf_shi_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_shi_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_shi_e = gsl_sf_shi_e(x, res)
    result = res
  end function fgsl_sf_shi_e
  function fgsl_sf_chi(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_chi
    fgsl_sf_chi = gsl_sf_chi(x)
  end function fgsl_sf_chi
  function fgsl_sf_chi_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_chi_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_chi_e = gsl_sf_chi_e(x, res)
    result = res
  end function fgsl_sf_chi_e
  function fgsl_sf_expint_3(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_expint_3
    fgsl_sf_expint_3 = gsl_sf_expint_3(x)
  end function fgsl_sf_expint_3
  function fgsl_sf_expint_3_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_expint_3_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_expint_3_e = gsl_sf_expint_3_e(x, res)
    result = res
  end function fgsl_sf_expint_3_e
  function fgsl_sf_si(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_si
    fgsl_sf_si = gsl_sf_si(x)
  end function fgsl_sf_si
  function fgsl_sf_si_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_si_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_si_e = gsl_sf_si_e(x, res)
    result = res
  end function fgsl_sf_si_e
  function fgsl_sf_ci(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_ci
    fgsl_sf_ci = gsl_sf_ci(x)
  end function fgsl_sf_ci
  function fgsl_sf_ci_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_ci_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_ci_e = gsl_sf_ci_e(x, res)
    result = res
  end function fgsl_sf_ci_e
  function fgsl_sf_atanint(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_atanint
    fgsl_sf_atanint = gsl_sf_atanint(x)
  end function fgsl_sf_atanint
  function fgsl_sf_atanint_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_atanint_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_atanint_e = gsl_sf_atanint_e(x, res)
    result = res
  end function fgsl_sf_atanint_e
  function fgsl_sf_fermi_dirac_m1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_fermi_dirac_m1
    fgsl_sf_fermi_dirac_m1 = gsl_sf_fermi_dirac_m1(x)
  end function fgsl_sf_fermi_dirac_m1
  function fgsl_sf_fermi_dirac_m1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_fermi_dirac_m1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_fermi_dirac_m1_e = gsl_sf_fermi_dirac_m1_e(x, res)
    result = res
  end function fgsl_sf_fermi_dirac_m1_e
  function fgsl_sf_fermi_dirac_0(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_fermi_dirac_0
    fgsl_sf_fermi_dirac_0 = gsl_sf_fermi_dirac_0(x)
  end function fgsl_sf_fermi_dirac_0
  function fgsl_sf_fermi_dirac_0_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_fermi_dirac_0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_fermi_dirac_0_e = gsl_sf_fermi_dirac_0_e(x, res)
    result = res
  end function fgsl_sf_fermi_dirac_0_e
  function fgsl_sf_fermi_dirac_1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_fermi_dirac_1
    fgsl_sf_fermi_dirac_1 = gsl_sf_fermi_dirac_1(x)
  end function fgsl_sf_fermi_dirac_1
  function fgsl_sf_fermi_dirac_1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_fermi_dirac_1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_fermi_dirac_1_e = gsl_sf_fermi_dirac_1_e(x, res)
    result = res
  end function fgsl_sf_fermi_dirac_1_e
  function fgsl_sf_fermi_dirac_2(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_fermi_dirac_2
    fgsl_sf_fermi_dirac_2 = gsl_sf_fermi_dirac_2(x)
  end function fgsl_sf_fermi_dirac_2
  function fgsl_sf_fermi_dirac_2_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_fermi_dirac_2_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_fermi_dirac_2_e = gsl_sf_fermi_dirac_2_e(x, res)
    result = res
  end function fgsl_sf_fermi_dirac_2_e
  function fgsl_sf_fermi_dirac_int(i, x) 
    integer(fgsl_int), intent(in) :: i
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_fermi_dirac_int
    fgsl_sf_fermi_dirac_int = gsl_sf_fermi_dirac_int(i, x)
  end function fgsl_sf_fermi_dirac_int
  function fgsl_sf_fermi_dirac_int_e(i, x, result) 
    integer(fgsl_int), intent(in) :: i
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_fermi_dirac_int_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_fermi_dirac_int_e = gsl_sf_fermi_dirac_int_e(i, x, res)
    result = res
  end function fgsl_sf_fermi_dirac_int_e
  function fgsl_sf_fermi_dirac_mhalf(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_fermi_dirac_mhalf
    fgsl_sf_fermi_dirac_mhalf = gsl_sf_fermi_dirac_mhalf(x)
  end function fgsl_sf_fermi_dirac_mhalf
  function fgsl_sf_fermi_dirac_mhalf_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_fermi_dirac_mhalf_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_fermi_dirac_mhalf_e = gsl_sf_fermi_dirac_mhalf_e(x, res)
    result = res
  end function fgsl_sf_fermi_dirac_mhalf_e
  function fgsl_sf_fermi_dirac_half(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_fermi_dirac_half
    fgsl_sf_fermi_dirac_half = gsl_sf_fermi_dirac_half(x)
  end function fgsl_sf_fermi_dirac_half
  function fgsl_sf_fermi_dirac_half_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_fermi_dirac_half_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_fermi_dirac_half_e = gsl_sf_fermi_dirac_half_e(x, res)
    result = res
  end function fgsl_sf_fermi_dirac_half_e
  function fgsl_sf_fermi_dirac_3half(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_fermi_dirac_3half
    fgsl_sf_fermi_dirac_3half = gsl_sf_fermi_dirac_3half(x)
  end function fgsl_sf_fermi_dirac_3half
  function fgsl_sf_fermi_dirac_3half_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_fermi_dirac_3half_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_fermi_dirac_3half_e = gsl_sf_fermi_dirac_3half_e(x, res)
    result = res
  end function fgsl_sf_fermi_dirac_3half_e
  function fgsl_sf_fermi_dirac_inc_0(x, b) 
    real(fgsl_double), intent(in) :: x, b
    real(fgsl_double) :: fgsl_sf_fermi_dirac_inc_0
    fgsl_sf_fermi_dirac_inc_0 = gsl_sf_fermi_dirac_inc_0(x, b)
  end function fgsl_sf_fermi_dirac_inc_0
  function fgsl_sf_fermi_dirac_inc_0_e(x, b, result) 
    real(fgsl_double), intent(in) :: x, b
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_fermi_dirac_inc_0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_fermi_dirac_inc_0_e = gsl_sf_fermi_dirac_inc_0_e(x, b, res)
    result = res
  end function fgsl_sf_fermi_dirac_inc_0_e
  function fgsl_sf_gamma(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_gamma
    fgsl_sf_gamma = gsl_sf_gamma(x)
  end function fgsl_sf_gamma
  function fgsl_sf_gamma_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_gamma_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_gamma_e = gsl_sf_gamma_e(x, res)
    result = res
  end function fgsl_sf_gamma_e
  function fgsl_sf_lngamma(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_lngamma
    fgsl_sf_lngamma = gsl_sf_lngamma(x)
  end function fgsl_sf_lngamma
  function fgsl_sf_lngamma_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_lngamma_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lngamma_e = gsl_sf_lngamma_e(x, res)
    result = res
  end function fgsl_sf_lngamma_e
  function fgsl_sf_lngamma_sgn_e(x, result_lg, sgn) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result_lg
    real(fgsl_double), intent(out) :: sgn
    integer(fgsl_int) :: fgsl_sf_lngamma_sgn_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lngamma_sgn_e = gsl_sf_lngamma_sgn_e(x, res, sgn)
    result_lg = res
  end function fgsl_sf_lngamma_sgn_e
  function fgsl_sf_gammastar(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_gammastar
    fgsl_sf_gammastar = gsl_sf_gammastar(x)
  end function fgsl_sf_gammastar
  function fgsl_sf_gammastar_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_gammastar_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_gammastar_e = gsl_sf_gammastar_e(x, res)
    result = res
  end function fgsl_sf_gammastar_e
  function fgsl_sf_gammainv(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_gammainv
    fgsl_sf_gammainv = gsl_sf_gammainv(x)
  end function fgsl_sf_gammainv
  function fgsl_sf_gammainv_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_gammainv_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_gammainv_e = gsl_sf_gammainv_e(x, res)
    result = res
  end function fgsl_sf_gammainv_e
  function fgsl_sf_lngamma_complex_e(zr, zi, lnr, arg)
    real(fgsl_double), intent(in) :: zr, zi
    type(fgsl_sf_result), intent(out) :: lnr, arg
    integer(fgsl_int) :: fgsl_sf_lngamma_complex_e
!
    type(gsl_sf_result) :: lnr_loc, arg_loc
    fgsl_sf_lngamma_complex_e = gsl_sf_lngamma_complex_e(zr, zi, lnr_loc, arg_loc)
    lnr = lnr_loc
    arg = arg_loc
  end function fgsl_sf_lngamma_complex_e
  function fgsl_sf_fact(n) 
    integer(c_int), intent(in) :: n
    real(fgsl_double) :: fgsl_sf_fact
    fgsl_sf_fact = gsl_sf_fact(n)
  end function fgsl_sf_fact
  function fgsl_sf_fact_e(n, result) 
    integer(c_int), intent(in) :: n
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_fact_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_fact_e = gsl_sf_fact_e(n, res)
    result = res
  end function fgsl_sf_fact_e
  function fgsl_sf_doublefact(n) 
    integer(c_int), intent(in) :: n
    real(fgsl_double) :: fgsl_sf_doublefact
    fgsl_sf_doublefact = gsl_sf_doublefact(n)
  end function fgsl_sf_doublefact
  function fgsl_sf_doublefact_e(n, result) 
    integer(c_int), intent(in) :: n
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_doublefact_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_doublefact_e = gsl_sf_doublefact_e(n, res)
    result = res
  end function fgsl_sf_doublefact_e
  function fgsl_sf_lnfact(n) 
    integer(c_int), intent(in) :: n
    real(fgsl_double) :: fgsl_sf_lnfact
    fgsl_sf_lnfact = gsl_sf_lnfact(n)
  end function fgsl_sf_lnfact
  function fgsl_sf_lnfact_e(n, result) 
    integer(c_int), intent(in) :: n
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_lnfact_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lnfact_e = gsl_sf_lnfact_e(n, res)
    result = res
  end function fgsl_sf_lnfact_e
  function fgsl_sf_lndoublefact(n) 
    integer(c_int), intent(in) :: n
    real(fgsl_double) :: fgsl_sf_lndoublefact
    fgsl_sf_lndoublefact = gsl_sf_lndoublefact(n)
  end function fgsl_sf_lndoublefact
  function fgsl_sf_lndoublefact_e(n, result) 
    integer(c_int), intent(in) :: n
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_lndoublefact_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lndoublefact_e = gsl_sf_lndoublefact_e(n, res)
    result = res
  end function fgsl_sf_lndoublefact_e
  function fgsl_sf_choose(n, m) 
    integer(c_int), intent(in) :: n, m
    real(fgsl_double) :: fgsl_sf_choose
    fgsl_sf_choose = gsl_sf_choose(n, m)
  end function fgsl_sf_choose
  function fgsl_sf_choose_e(n, m, result) 
    integer(c_int), intent(in) :: n, m
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_choose_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_choose_e = gsl_sf_choose_e(n, m, res)
    result = res
  end function fgsl_sf_choose_e
  function fgsl_sf_lnchoose(n, m) 
    integer(c_int), intent(in) :: n, m
    real(fgsl_double) :: fgsl_sf_lnchoose
    fgsl_sf_lnchoose = gsl_sf_lnchoose(n, m)
  end function fgsl_sf_lnchoose
  function fgsl_sf_lnchoose_e(n, m, result) 
    integer(c_int), intent(in) :: n, m
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_lnchoose_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lnchoose_e = gsl_sf_lnchoose_e(n, m, res)
    result = res
  end function fgsl_sf_lnchoose_e
  function fgsl_sf_taylorcoeff(n, x) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_taylorcoeff
    fgsl_sf_taylorcoeff = gsl_sf_taylorcoeff(n, x)
  end function fgsl_sf_taylorcoeff
  function fgsl_sf_taylorcoeff_e(n, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_taylorcoeff_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_taylorcoeff_e = gsl_sf_taylorcoeff_e(n, x, res)
    result = res
  end function fgsl_sf_taylorcoeff_e
  function fgsl_sf_poch(a, x) 
    real(fgsl_double), intent(in) :: a, x
    real(fgsl_double) :: fgsl_sf_poch
    fgsl_sf_poch = gsl_sf_poch(a, x)
  end function fgsl_sf_poch
  function fgsl_sf_poch_e(a, x, result) 
    real(fgsl_double), intent(in) :: a, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_poch_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_poch_e = gsl_sf_poch_e(a, x, res)
    result = res
  end function fgsl_sf_poch_e
  function fgsl_sf_lnpoch(a, x) 
    real(fgsl_double), intent(in) :: a, x
    real(fgsl_double) :: fgsl_sf_lnpoch
    fgsl_sf_lnpoch = gsl_sf_lnpoch(a, x)
  end function fgsl_sf_lnpoch
  function fgsl_sf_lnpoch_e(a, x, result) 
    real(fgsl_double), intent(in) :: a, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_lnpoch_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lnpoch_e = gsl_sf_lnpoch_e(a, x, res)
    result = res
  end function fgsl_sf_lnpoch_e
  function fgsl_sf_lnpoch_sgn_e(a, x, result_lg, sgn) 
    real(fgsl_double), intent(in) :: a, x
    type(fgsl_sf_result), intent(out) :: result_lg
    real(fgsl_double), intent(out) :: sgn
    integer(fgsl_int) :: fgsl_sf_lnpoch_sgn_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lnpoch_sgn_e = gsl_sf_lnpoch_sgn_e(a, x, res, sgn)
    result_lg = res
  end function fgsl_sf_lnpoch_sgn_e
  function fgsl_sf_pochrel(a, x) 
    real(fgsl_double), intent(in) :: a, x
    real(fgsl_double) :: fgsl_sf_pochrel
    fgsl_sf_pochrel = gsl_sf_pochrel(a, x)
  end function fgsl_sf_pochrel
  function fgsl_sf_pochrel_e(a, x, result) 
    real(fgsl_double), intent(in) :: a, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_pochrel_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_pochrel_e = gsl_sf_pochrel_e(a, x, res)
    result = res
  end function fgsl_sf_pochrel_e
  function fgsl_sf_gamma_inc(a, x) 
    real(fgsl_double), intent(in) :: a, x
    real(fgsl_double) :: fgsl_sf_gamma_inc
    fgsl_sf_gamma_inc = gsl_sf_gamma_inc(a, x)
  end function fgsl_sf_gamma_inc
  function fgsl_sf_gamma_inc_e(a, x, result) 
    real(fgsl_double), intent(in) :: a, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_gamma_inc_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_gamma_inc_e = gsl_sf_gamma_inc_e(a, x, res)
    result = res
  end function fgsl_sf_gamma_inc_e
  function fgsl_sf_gamma_inc_q(a, x) 
    real(fgsl_double), intent(in) :: a, x
    real(fgsl_double) :: fgsl_sf_gamma_inc_q
    fgsl_sf_gamma_inc_q = gsl_sf_gamma_inc_q(a, x)
  end function fgsl_sf_gamma_inc_q
  function fgsl_sf_gamma_inc_q_e(a, x, result) 
    real(fgsl_double), intent(in) :: a, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_gamma_inc_q_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_gamma_inc_q_e = gsl_sf_gamma_inc_q_e(a, x, res)
    result = res
  end function fgsl_sf_gamma_inc_q_e
  function fgsl_sf_gamma_inc_p(a, x) 
    real(fgsl_double), intent(in) :: a, x
    real(fgsl_double) :: fgsl_sf_gamma_inc_p
    fgsl_sf_gamma_inc_p = gsl_sf_gamma_inc_p(a, x)
  end function fgsl_sf_gamma_inc_p
  function fgsl_sf_gamma_inc_p_e(a, x, result) 
    real(fgsl_double), intent(in) :: a, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_gamma_inc_p_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_gamma_inc_p_e = gsl_sf_gamma_inc_p_e(a, x, res)
    result = res
  end function fgsl_sf_gamma_inc_p_e
  function fgsl_sf_beta(a, b) 
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_sf_beta
    fgsl_sf_beta = gsl_sf_beta(a, b)
  end function fgsl_sf_beta
  function fgsl_sf_beta_e(a, b, result) 
    real(fgsl_double), intent(in) :: a, b
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_beta_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_beta_e = gsl_sf_beta_e(a, b, res)
    result = res
  end function fgsl_sf_beta_e
  function fgsl_sf_lnbeta(a, b) 
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_sf_lnbeta
    fgsl_sf_lnbeta = gsl_sf_lnbeta(a, b)
  end function fgsl_sf_lnbeta
  function fgsl_sf_lnbeta_e(a, b, result) 
    real(fgsl_double), intent(in) :: a, b
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_lnbeta_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lnbeta_e = gsl_sf_lnbeta_e(a, b, res)
    result = res
  end function fgsl_sf_lnbeta_e
  function fgsl_sf_beta_inc(a, b, x) 
    real(fgsl_double), intent(in) :: a, b, x
    real(fgsl_double) :: fgsl_sf_beta_inc
    fgsl_sf_beta_inc = gsl_sf_beta_inc(a, b, x)
  end function fgsl_sf_beta_inc
  function fgsl_sf_beta_inc_e(a, b, x, result) 
    real(fgsl_double), intent(in) :: a, b, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_beta_inc_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_beta_inc_e = gsl_sf_beta_inc_e(a, b, x, res)
    result = res
  end function fgsl_sf_beta_inc_e
  function fgsl_sf_gegenpoly_1(lambda, x) 
    real(fgsl_double), intent(in) :: lambda, x
    real(fgsl_double) :: fgsl_sf_gegenpoly_1
    fgsl_sf_gegenpoly_1 = gsl_sf_gegenpoly_1(lambda, x)
  end function fgsl_sf_gegenpoly_1
  function fgsl_sf_gegenpoly_1_e(lambda, x, result) 
    real(fgsl_double), intent(in) :: lambda, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_gegenpoly_1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_gegenpoly_1_e = gsl_sf_gegenpoly_1_e(lambda, x, res)
    result = res
  end function fgsl_sf_gegenpoly_1_e
  function fgsl_sf_gegenpoly_2(lambda, x) 
    real(fgsl_double), intent(in) :: lambda, x
    real(fgsl_double) :: fgsl_sf_gegenpoly_2
    fgsl_sf_gegenpoly_2 = gsl_sf_gegenpoly_2(lambda, x)
  end function fgsl_sf_gegenpoly_2
  function fgsl_sf_gegenpoly_2_e(lambda, x, result) 
    real(fgsl_double), intent(in) :: lambda, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_gegenpoly_2_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_gegenpoly_2_e = gsl_sf_gegenpoly_2_e(lambda, x, res)
    result = res
  end function fgsl_sf_gegenpoly_2_e
  function fgsl_sf_gegenpoly_3(lambda, x) 
    real(fgsl_double), intent(in) :: lambda, x
    real(fgsl_double) :: fgsl_sf_gegenpoly_3
    fgsl_sf_gegenpoly_3 = gsl_sf_gegenpoly_3(lambda, x)
  end function fgsl_sf_gegenpoly_3
  function fgsl_sf_gegenpoly_3_e(lambda, x, result) 
    real(fgsl_double), intent(in) :: lambda, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_gegenpoly_3_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_gegenpoly_3_e = gsl_sf_gegenpoly_3_e(lambda, x, res)
    result = res
  end function fgsl_sf_gegenpoly_3_e
  function fgsl_sf_gegenpoly_n(n, lambda, x) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: lambda, x
    real(fgsl_double) :: fgsl_sf_gegenpoly_n
    fgsl_sf_gegenpoly_n = gsl_sf_gegenpoly_n(n, lambda, x)
  end function fgsl_sf_gegenpoly_n
  function fgsl_sf_gegenpoly_n_e(n, lambda, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: lambda, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_gegenpoly_n_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_gegenpoly_n_e = gsl_sf_gegenpoly_n_e(n, lambda, x, res)
    result = res
  end function fgsl_sf_gegenpoly_n_e
  function fgsl_sf_gegenpoly_array(nmax, lambda, x, result_array)
    integer(fgsl_int), intent(in) :: nmax
    real(fgsl_double), intent(in) :: lambda, x
    real(fgsl_double), intent(out) :: result_array(:)
    integer(fgsl_int) :: fgsl_sf_gegenpoly_array
    fgsl_sf_gegenpoly_array = gsl_sf_gegenpoly_array(nmax, lambda, x, result_array)
  end function fgsl_sf_gegenpoly_array
  function fgsl_sf_hyperg_0f1(c, x) 
    real(fgsl_double), intent(in) :: c, x
    real(fgsl_double) :: fgsl_sf_hyperg_0f1
    fgsl_sf_hyperg_0f1 = gsl_sf_hyperg_0f1(c, x)
  end function fgsl_sf_hyperg_0f1
  function fgsl_sf_hyperg_0f1_e(c, x, result) 
    real(fgsl_double), intent(in) :: c, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_0f1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hyperg_0f1_e = gsl_sf_hyperg_0f1_e(c, x, res)
    result = res
  end function fgsl_sf_hyperg_0f1_e
  function fgsl_sf_hyperg_1f1_int(m, n, x) 
    integer(fgsl_int), intent(in) :: m, n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_hyperg_1f1_int
    fgsl_sf_hyperg_1f1_int = gsl_sf_hyperg_1f1_int(m, n, x)
  end function fgsl_sf_hyperg_1f1_int
  function fgsl_sf_hyperg_1f1_int_e(m, n, x, result) 
    integer(fgsl_int), intent(in) :: m, n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_1f1_int_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hyperg_1f1_int_e = gsl_sf_hyperg_1f1_int_e(m, n, x, res)
    result = res
  end function fgsl_sf_hyperg_1f1_int_e
  function fgsl_sf_hyperg_1f1(a, b, x) 
    real(fgsl_double), intent(in) :: a, b, x
    real(fgsl_double) :: fgsl_sf_hyperg_1f1
    fgsl_sf_hyperg_1f1 = gsl_sf_hyperg_1f1(a, b, x)
  end function fgsl_sf_hyperg_1f1
  function fgsl_sf_hyperg_1f1_e(a, b, x, result) 
    real(fgsl_double), intent(in) :: a, b, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_1f1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hyperg_1f1_e = gsl_sf_hyperg_1f1_e(a, b, x, res)
    result = res
  end function fgsl_sf_hyperg_1f1_e
  function fgsl_sf_hyperg_u_int(m, n, x) 
    integer(fgsl_int), intent(in) :: m, n
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_hyperg_u_int
    fgsl_sf_hyperg_u_int = gsl_sf_hyperg_u_int(m, n, x)
  end function fgsl_sf_hyperg_u_int
  function fgsl_sf_hyperg_u_int_e(m, n, x, result) 
    integer(fgsl_int), intent(in) :: m, n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_u_int_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hyperg_u_int_e = gsl_sf_hyperg_u_int_e(m, n, x, res)
    result = res
  end function fgsl_sf_hyperg_u_int_e
  function fgsl_sf_hyperg_u_int_e10_e(m, n, x, result) 
    integer(fgsl_int), intent(in) :: m, n
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result_e10), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_u_int_e10_e
!
    type(gsl_sf_result_e10) :: res
    fgsl_sf_hyperg_u_int_e10_e = gsl_sf_hyperg_u_int_e10_e(m, n, x, res)
    result = res
  end function fgsl_sf_hyperg_u_int_e10_e
  function fgsl_sf_hyperg_u(a, b, x) 
    real(fgsl_double), intent(in) :: a, b, x
    real(fgsl_double) :: fgsl_sf_hyperg_u
    fgsl_sf_hyperg_u = gsl_sf_hyperg_u(a, b, x)
  end function fgsl_sf_hyperg_u
  function fgsl_sf_hyperg_u_e(a, b, x, result) 
    real(fgsl_double), intent(in) :: a, b, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_u_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hyperg_u_e = gsl_sf_hyperg_u_e(a, b, x, res)
    result = res
  end function fgsl_sf_hyperg_u_e
  function fgsl_sf_hyperg_u_e10_e(a, b, x, result) 
    real(fgsl_double), intent(in) :: a, b, x
    type(fgsl_sf_result_e10), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_u_e10_e
!
    type(gsl_sf_result_e10) :: res
    fgsl_sf_hyperg_u_e10_e = gsl_sf_hyperg_u_e10_e(a, b, x, res)
    result = res
  end function fgsl_sf_hyperg_u_e10_e
  function fgsl_sf_hyperg_2f1(a, b, c, x) 
    real(fgsl_double), intent(in) :: a, b, c, x
    real(fgsl_double) :: fgsl_sf_hyperg_2f1
    fgsl_sf_hyperg_2f1 = gsl_sf_hyperg_2f1(a, b, c, x)
  end function fgsl_sf_hyperg_2f1
  function fgsl_sf_hyperg_2f1_e(a, b, c, x, result) 
    real(fgsl_double), intent(in) :: a, b, c, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_2f1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hyperg_2f1_e = gsl_sf_hyperg_2f1_e(a, b, c, x, res)
    result = res
  end function fgsl_sf_hyperg_2f1_e
  function fgsl_sf_hyperg_2f1_conj(ar, ai, c, x) 
    real(fgsl_double), intent(in) :: ar, ai, c, x
    real(fgsl_double) :: fgsl_sf_hyperg_2f1_conj
    fgsl_sf_hyperg_2f1_conj = gsl_sf_hyperg_2f1_conj(ar, ai, c, x)
  end function fgsl_sf_hyperg_2f1_conj
  function fgsl_sf_hyperg_2f1_conj_e(ar, ai, c, x, result) 
    real(fgsl_double), intent(in) :: ar, ai, c, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_2f1_conj_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hyperg_2f1_conj_e = gsl_sf_hyperg_2f1_conj_e(ar, ai, c, x, res)
    result = res
  end function fgsl_sf_hyperg_2f1_conj_e
  function fgsl_sf_hyperg_2f1_renorm(a, b, c, x) 
    real(fgsl_double), intent(in) :: a, b, c, x
    real(fgsl_double) :: fgsl_sf_hyperg_2f1_renorm
    fgsl_sf_hyperg_2f1_renorm = gsl_sf_hyperg_2f1_renorm(a, b, c, x)
  end function fgsl_sf_hyperg_2f1_renorm
  function fgsl_sf_hyperg_2f1_renorm_e(a, b, c, x, result) 
    real(fgsl_double), intent(in) :: a, b, c, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_2f1_renorm_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hyperg_2f1_renorm_e = gsl_sf_hyperg_2f1_renorm_e(a, b, c, x, res)
    result = res
  end function fgsl_sf_hyperg_2f1_renorm_e
  function fgsl_sf_hyperg_2f1_conj_renorm(ar, ai, c, x) 
    real(fgsl_double), intent(in) :: ar, ai, c, x
    real(fgsl_double) :: fgsl_sf_hyperg_2f1_conj_renorm
    fgsl_sf_hyperg_2f1_conj_renorm = gsl_sf_hyperg_2f1_conj_renorm(ar, ai, c, x)
  end function fgsl_sf_hyperg_2f1_conj_renorm
  function fgsl_sf_hyperg_2f1_conj_renorm_e(ar, ai, c, x, result) 
    real(fgsl_double), intent(in) :: ar, ai, c, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_2f1_conj_renorm_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hyperg_2f1_conj_renorm_e = gsl_sf_hyperg_2f1_conj_renorm_e(ar, ai, c, x, res)
    result = res
  end function fgsl_sf_hyperg_2f1_conj_renorm_e
  function fgsl_sf_hyperg_2f0(a, b, x) 
    real(fgsl_double), intent(in) :: a, b, x
    real(fgsl_double) :: fgsl_sf_hyperg_2f0
    fgsl_sf_hyperg_2f0 = gsl_sf_hyperg_2f0(a, b, x)
  end function fgsl_sf_hyperg_2f0
  function fgsl_sf_hyperg_2f0_e(a, b, x, result) 
    real(fgsl_double), intent(in) :: a, b, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hyperg_2f0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hyperg_2f0_e = gsl_sf_hyperg_2f0_e(a, b, x, res)
    result = res
  end function fgsl_sf_hyperg_2f0_e
  function fgsl_sf_laguerre_1(a, x) 
    real(fgsl_double), intent(in) :: a, x
    real(fgsl_double) :: fgsl_sf_laguerre_1
    fgsl_sf_laguerre_1 = gsl_sf_laguerre_1(a, x)
  end function fgsl_sf_laguerre_1
  function fgsl_sf_laguerre_1_e(a, x, result) 
    real(fgsl_double), intent(in) :: a, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_laguerre_1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_laguerre_1_e = gsl_sf_laguerre_1_e(a, x, res)
    result = res
  end function fgsl_sf_laguerre_1_e
  function fgsl_sf_laguerre_2(a, x) 
    real(fgsl_double), intent(in) :: a, x
    real(fgsl_double) :: fgsl_sf_laguerre_2
    fgsl_sf_laguerre_2 = gsl_sf_laguerre_2(a, x)
  end function fgsl_sf_laguerre_2
  function fgsl_sf_laguerre_2_e(a, x, result) 
    real(fgsl_double), intent(in) :: a, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_laguerre_2_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_laguerre_2_e = gsl_sf_laguerre_2_e(a, x, res)
    result = res
  end function fgsl_sf_laguerre_2_e
  function fgsl_sf_laguerre_3(a, x) 
    real(fgsl_double), intent(in) :: a, x
    real(fgsl_double) :: fgsl_sf_laguerre_3
    fgsl_sf_laguerre_3 = gsl_sf_laguerre_3(a, x)
  end function fgsl_sf_laguerre_3
  function fgsl_sf_laguerre_3_e(a, x, result) 
    real(fgsl_double), intent(in) :: a, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_laguerre_3_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_laguerre_3_e = gsl_sf_laguerre_3_e(a, x, res)
    result = res
  end function fgsl_sf_laguerre_3_e
  function fgsl_sf_laguerre_n(n, a, x) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: a, x
    real(fgsl_double) :: fgsl_sf_laguerre_n
    fgsl_sf_laguerre_n = gsl_sf_laguerre_n(n, a, x)
  end function fgsl_sf_laguerre_n
  function fgsl_sf_laguerre_n_e(n, a, x, result) 
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: a, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_laguerre_n_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_laguerre_n_e = gsl_sf_laguerre_n_e(n, a, x, res)
    result = res
  end function fgsl_sf_laguerre_n_e
  function fgsl_sf_lambert_w0(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_lambert_w0
    fgsl_sf_lambert_w0 = gsl_sf_lambert_w0(x)
  end function fgsl_sf_lambert_w0
  function fgsl_sf_lambert_w0_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_lambert_w0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lambert_w0_e = gsl_sf_lambert_w0_e(x, res)
    result = res
  end function fgsl_sf_lambert_w0_e
  function fgsl_sf_lambert_wm1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_lambert_wm1
    fgsl_sf_lambert_wm1 = gsl_sf_lambert_wm1(x)
  end function fgsl_sf_lambert_wm1
  function fgsl_sf_lambert_wm1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_lambert_wm1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lambert_wm1_e = gsl_sf_lambert_wm1_e(x, res)
    result = res
  end function fgsl_sf_lambert_wm1_e
  function fgsl_sf_legendre_p1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_legendre_p1
    fgsl_sf_legendre_p1 = gsl_sf_legendre_p1(x)
  end function fgsl_sf_legendre_p1
  function fgsl_sf_legendre_p1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_p1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_p1_e = gsl_sf_legendre_p1_e(x, res)
    result = res
  end function fgsl_sf_legendre_p1_e
  function fgsl_sf_legendre_p2(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_legendre_p2
    fgsl_sf_legendre_p2 = gsl_sf_legendre_p2(x)
  end function fgsl_sf_legendre_p2
  function fgsl_sf_legendre_p2_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_p2_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_p2_e = gsl_sf_legendre_p2_e(x, res)
    result = res
  end function fgsl_sf_legendre_p2_e
  function fgsl_sf_legendre_p3(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_legendre_p3
    fgsl_sf_legendre_p3 = gsl_sf_legendre_p3(x)
  end function fgsl_sf_legendre_p3
  function fgsl_sf_legendre_p3_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_p3_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_p3_e = gsl_sf_legendre_p3_e(x, res)
    result = res
  end function fgsl_sf_legendre_p3_e
  function fgsl_sf_legendre_pl(l, x)
    integer(fgsl_int), intent(in) :: l
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_legendre_pl
    fgsl_sf_legendre_pl = gsl_sf_legendre_pl(l, x)
  end function fgsl_sf_legendre_pl
  function fgsl_sf_legendre_pl_e(l, x, result) 
    integer(fgsl_int), intent(in) :: l
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_pl_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_pl_e = gsl_sf_legendre_pl_e(l, x, res)
    result = res
  end function fgsl_sf_legendre_pl_e
  function fgsl_sf_legendre_pl_array(lmax, x, result_array)
    integer(fgsl_int), intent(in) :: lmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result_array(:)
    real(fgsl_double) :: fgsl_sf_legendre_pl_array
    fgsl_sf_legendre_pl_array = gsl_sf_legendre_pl_array(lmax, x, result_array)
  end function fgsl_sf_legendre_pl_array
  function fgsl_sf_legendre_pl_deriv_array(lmax, x, result_array, deriv_array)
    integer(fgsl_int), intent(in) :: lmax
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result_array(:), deriv_array(:)
    real(fgsl_double) :: fgsl_sf_legendre_pl_deriv_array
    fgsl_sf_legendre_pl_deriv_array = &
         gsl_sf_legendre_pl_deriv_array(lmax, x, result_array, deriv_array)
  end function fgsl_sf_legendre_pl_deriv_array
  function fgsl_sf_legendre_q0(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_legendre_q0
    fgsl_sf_legendre_q0 = gsl_sf_legendre_q0(x)
  end function fgsl_sf_legendre_q0
  function fgsl_sf_legendre_q0_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_q0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_q0_e = gsl_sf_legendre_q0_e(x, res)
    result = res
  end function fgsl_sf_legendre_q0_e
  function fgsl_sf_legendre_q1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_legendre_q1
    fgsl_sf_legendre_q1 = gsl_sf_legendre_q1(x)
  end function fgsl_sf_legendre_q1
  function fgsl_sf_legendre_q1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_q1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_q1_e = gsl_sf_legendre_q1_e(x, res)
    result = res
  end function fgsl_sf_legendre_q1_e
  function fgsl_sf_legendre_ql(l, x)
    integer(fgsl_int), intent(in) :: l
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_legendre_ql
    fgsl_sf_legendre_ql = gsl_sf_legendre_ql(l, x)
  end function fgsl_sf_legendre_ql
  function fgsl_sf_legendre_ql_e(l, x, result) 
    integer(fgsl_int), intent(in) :: l
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_ql_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_ql_e = gsl_sf_legendre_ql_e(l, x, res)
    result = res
  end function fgsl_sf_legendre_ql_e
  function fgsl_sf_legendre_plm(l, m, x)
    integer(fgsl_int), intent(in) :: l, m
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_legendre_plm
    fgsl_sf_legendre_plm = gsl_sf_legendre_plm(l, m, x)
  end function fgsl_sf_legendre_plm
  function fgsl_sf_legendre_plm_e(l, m, x, result) 
    integer(fgsl_int), intent(in) :: l, m
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_plm_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_plm_e = gsl_sf_legendre_plm_e(l, m, x, res)
    result = res
  end function fgsl_sf_legendre_plm_e
  function fgsl_sf_legendre_plm_array(lmax, m, x, result_array)
    integer(fgsl_int), intent(in) :: lmax, m
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result_array(:)
    real(fgsl_double) :: fgsl_sf_legendre_plm_array
    fgsl_sf_legendre_plm_array = gsl_sf_legendre_plm_array(lmax, m, x, result_array)
  end function fgsl_sf_legendre_plm_array
  function fgsl_sf_legendre_plm_deriv_array(lmax, m, x, result_array, deriv_array)
    integer(fgsl_int), intent(in) :: lmax, m
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result_array(:), deriv_array(:)
    real(fgsl_double) :: fgsl_sf_legendre_plm_deriv_array
    fgsl_sf_legendre_plm_deriv_array = &
         gsl_sf_legendre_plm_deriv_array(lmax, m, x, result_array, deriv_array)
  end function fgsl_sf_legendre_plm_deriv_array
  function fgsl_sf_legendre_sphplm(l, m, x)
    integer(fgsl_int), intent(in) :: l, m
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_legendre_sphplm
    fgsl_sf_legendre_sphplm = gsl_sf_legendre_sphplm(l, m, x)
  end function fgsl_sf_legendre_sphplm
  function fgsl_sf_legendre_sphplm_e(l, m, x, result) 
    integer(fgsl_int), intent(in) :: l, m
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_sphplm_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_sphplm_e = gsl_sf_legendre_sphplm_e(l, m, x, res)
    result = res
  end function fgsl_sf_legendre_sphplm_e
  function fgsl_sf_legendre_sphplm_array(lmax, m, x, result_array)
    integer(fgsl_int), intent(in) :: lmax, m
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result_array(:)
    real(fgsl_double) :: fgsl_sf_legendre_sphplm_array
    fgsl_sf_legendre_sphplm_array = gsl_sf_legendre_sphplm_array(lmax, m, x, result_array)
  end function fgsl_sf_legendre_sphplm_array
  function fgsl_sf_legendre_sphplm_deriv_array(lmax, m, x, result_array, deriv_array)
    integer(fgsl_int), intent(in) :: lmax, m
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result_array(:), deriv_array(:)
    real(fgsl_double) :: fgsl_sf_legendre_sphplm_deriv_array
    fgsl_sf_legendre_sphplm_deriv_array = &
         gsl_sf_legendre_sphplm_deriv_array(lmax, m, x, result_array, deriv_array)
  end function fgsl_sf_legendre_sphplm_deriv_array
  function fgsl_sf_legendre_array_size(lmax, m)
    integer(fgsl_int), intent(in) :: lmax, m
    integer(c_int) :: fgsl_sf_legendre_array_size
    fgsl_sf_legendre_array_size = gsl_sf_legendre_array_size(lmax, m)
  end function fgsl_sf_legendre_array_size
  function fgsl_sf_conicalp_half(lambda, x) 
    real(fgsl_double), intent(in) :: lambda, x
    real(fgsl_double) :: fgsl_sf_conicalp_half
    fgsl_sf_conicalp_half = gsl_sf_conicalp_half(lambda, x)
  end function fgsl_sf_conicalp_half
  function fgsl_sf_conicalp_half_e(lambda, x, result) 
    real(fgsl_double), intent(in) :: lambda, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_conicalp_half_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_conicalp_half_e = gsl_sf_conicalp_half_e(lambda, x, res)
    result = res
  end function fgsl_sf_conicalp_half_e
  function fgsl_sf_conicalp_mhalf(lambda, x) 
    real(fgsl_double), intent(in) :: lambda, x
    real(fgsl_double) :: fgsl_sf_conicalp_mhalf
    fgsl_sf_conicalp_mhalf = gsl_sf_conicalp_mhalf(lambda, x)
  end function fgsl_sf_conicalp_mhalf
  function fgsl_sf_conicalp_mhalf_e(lambda, x, result) 
    real(fgsl_double), intent(in) :: lambda, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_conicalp_mhalf_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_conicalp_mhalf_e = gsl_sf_conicalp_mhalf_e(lambda, x, res)
    result = res
  end function fgsl_sf_conicalp_mhalf_e
  function fgsl_sf_conicalp_0(lambda, x) 
    real(fgsl_double), intent(in) :: lambda, x
    real(fgsl_double) :: fgsl_sf_conicalp_0
    fgsl_sf_conicalp_0 = gsl_sf_conicalp_0(lambda, x)
  end function fgsl_sf_conicalp_0
  function fgsl_sf_conicalp_0_e(lambda, x, result) 
    real(fgsl_double), intent(in) :: lambda, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_conicalp_0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_conicalp_0_e = gsl_sf_conicalp_0_e(lambda, x, res)
    result = res
  end function fgsl_sf_conicalp_0_e
  function fgsl_sf_conicalp_1(lambda, x) 
    real(fgsl_double), intent(in) :: lambda, x
    real(fgsl_double) :: fgsl_sf_conicalp_1
    fgsl_sf_conicalp_1 = gsl_sf_conicalp_1(lambda, x)
  end function fgsl_sf_conicalp_1
  function fgsl_sf_conicalp_1_e(lambda, x, result) 
    real(fgsl_double), intent(in) :: lambda, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_conicalp_1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_conicalp_1_e = gsl_sf_conicalp_1_e(lambda, x, res)
    result = res
  end function fgsl_sf_conicalp_1_e
  function fgsl_sf_conicalp_sph_reg(l, lambda, x) 
    integer(fgsl_int), intent(in) :: l
    real(fgsl_double), intent(in) :: lambda, x
    real(fgsl_double) :: fgsl_sf_conicalp_sph_reg
    fgsl_sf_conicalp_sph_reg = gsl_sf_conicalp_sph_reg(l, lambda, x)
  end function fgsl_sf_conicalp_sph_reg
  function fgsl_sf_conicalp_sph_reg_e(l, lambda, x, result) 
    integer(fgsl_int), intent(in) :: l
    real(fgsl_double), intent(in) :: lambda, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_conicalp_sph_reg_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_conicalp_sph_reg_e = gsl_sf_conicalp_sph_reg_e(l, lambda, x, res)
    result = res
  end function fgsl_sf_conicalp_sph_reg_e
  function fgsl_sf_conicalp_cyl_reg(l, lambda, x) 
    integer(fgsl_int), intent(in) :: l
    real(fgsl_double), intent(in) :: lambda, x
    real(fgsl_double) :: fgsl_sf_conicalp_cyl_reg
    fgsl_sf_conicalp_cyl_reg = gsl_sf_conicalp_cyl_reg(l, lambda, x)
  end function fgsl_sf_conicalp_cyl_reg
  function fgsl_sf_conicalp_cyl_reg_e(l, lambda, x, result) 
    integer(fgsl_int), intent(in) :: l
    real(fgsl_double), intent(in) :: lambda, x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_conicalp_cyl_reg_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_conicalp_cyl_reg_e = gsl_sf_conicalp_cyl_reg_e(l, lambda, x, res)
    result = res
  end function fgsl_sf_conicalp_cyl_reg_e
  function fgsl_sf_legendre_h3d_0(lambda, eta) 
    real(fgsl_double), intent(in) :: lambda, eta
    real(fgsl_double) :: fgsl_sf_legendre_h3d_0
    fgsl_sf_legendre_h3d_0 = gsl_sf_legendre_h3d_0(lambda, eta)
  end function fgsl_sf_legendre_h3d_0
  function fgsl_sf_legendre_h3d_0_e(lambda, eta, result) 
    real(fgsl_double), intent(in) :: lambda, eta
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_h3d_0_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_h3d_0_e = gsl_sf_legendre_h3d_0_e(lambda, eta, res)
    result = res
  end function fgsl_sf_legendre_h3d_0_e
  function fgsl_sf_legendre_h3d_1(lambda, eta) 
    real(fgsl_double), intent(in) :: lambda, eta
    real(fgsl_double) :: fgsl_sf_legendre_h3d_1
    fgsl_sf_legendre_h3d_1 = gsl_sf_legendre_h3d_1(lambda, eta)
  end function fgsl_sf_legendre_h3d_1
  function fgsl_sf_legendre_h3d_1_e(lambda, eta, result) 
    real(fgsl_double), intent(in) :: lambda, eta
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_h3d_1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_h3d_1_e = gsl_sf_legendre_h3d_1_e(lambda, eta, res)
    result = res
  end function fgsl_sf_legendre_h3d_1_e
  function fgsl_sf_legendre_h3d(l, lambda, eta) 
    integer(fgsl_int), intent(in) :: l
    real(fgsl_double), intent(in) :: lambda, eta
    real(fgsl_double) :: fgsl_sf_legendre_h3d
    fgsl_sf_legendre_h3d = gsl_sf_legendre_h3d(l, lambda, eta)
  end function fgsl_sf_legendre_h3d
  function fgsl_sf_legendre_h3d_e(l, lambda, eta, result) 
    integer(fgsl_int), intent(in) :: l
    real(fgsl_double), intent(in) :: lambda, eta
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_legendre_h3d_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_legendre_h3d_e = gsl_sf_legendre_h3d_e(l, lambda, eta, res)
    result = res
  end function fgsl_sf_legendre_h3d_e
  function fgsl_sf_legendre_h3d_array(lmax, lambda, eta, result_array)
    integer(fgsl_int), intent(in) :: lmax
    real(fgsl_double), intent(in) :: lambda, eta
    real(fgsl_double), intent(out) :: result_array(:)
    integer(fgsl_int) :: fgsl_sf_legendre_h3d_array
    fgsl_sf_legendre_h3d_array = gsl_sf_legendre_h3d_array(lmax, lambda, eta, result_array)
  end function fgsl_sf_legendre_h3d_array
  function fgsl_sf_log(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_log
    fgsl_sf_log = gsl_sf_log(x)
  end function fgsl_sf_log
  function fgsl_sf_log_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_log_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_log_e = gsl_sf_log_e(x, res)
    result = res
  end function fgsl_sf_log_e
  function fgsl_sf_log_abs(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_log_abs
    fgsl_sf_log_abs = gsl_sf_log_abs(x)
  end function fgsl_sf_log_abs
  function fgsl_sf_log_abs_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_log_abs_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_log_abs_e = gsl_sf_log_abs_e(x, res)
    result = res
  end function fgsl_sf_log_abs_e
  function fgsl_sf_complex_log_e(zr, zi, lnr, theta) 
    real(fgsl_double), intent(in) :: zr, zi
    type(fgsl_sf_result), intent(out) :: lnr, theta
    integer(fgsl_int) :: fgsl_sf_complex_log_e
!
    type(gsl_sf_result) :: lnr_loc, theta_loc
    fgsl_sf_complex_log_e = gsl_sf_complex_log_e(zr, zi, lnr_loc, theta_loc)
    lnr = lnr_loc
    theta = theta_loc
  end function fgsl_sf_complex_log_e
  function fgsl_sf_log_1plusx(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_log_1plusx
    fgsl_sf_log_1plusx = gsl_sf_log_1plusx(x)
  end function fgsl_sf_log_1plusx
  function fgsl_sf_log_1plusx_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_log_1plusx_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_log_1plusx_e = gsl_sf_log_1plusx_e(x, res)
    result = res
  end function fgsl_sf_log_1plusx_e
  function fgsl_sf_log_1plusx_mx(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_log_1plusx_mx
    fgsl_sf_log_1plusx_mx = gsl_sf_log_1plusx_mx(x)
  end function fgsl_sf_log_1plusx_mx
  function fgsl_sf_log_1plusx_mx_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_log_1plusx_mx_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_log_1plusx_mx_e = gsl_sf_log_1plusx_mx_e(x, res)
    result = res
  end function fgsl_sf_log_1plusx_mx_e
  function fgsl_sf_psi_int(n) 
    integer(c_int), intent(in) :: n
    real(fgsl_double) :: fgsl_sf_psi_int
    fgsl_sf_psi_int = gsl_sf_psi_int(n)
  end function fgsl_sf_psi_int
  function fgsl_sf_psi_int_e(n, result) 
    integer(c_int), intent(in) :: n
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_psi_int_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_psi_int_e = gsl_sf_psi_int_e(n, res)
    result = res
  end function fgsl_sf_psi_int_e
  function fgsl_sf_psi(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_psi
    fgsl_sf_psi = gsl_sf_psi(x)
  end function fgsl_sf_psi
  function fgsl_sf_psi_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_psi_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_psi_e = gsl_sf_psi_e(x, res)
    result = res
  end function fgsl_sf_psi_e
  function fgsl_sf_psi_1_int(n) 
    integer(c_int), intent(in) :: n
    real(fgsl_double) :: fgsl_sf_psi_1_int
    fgsl_sf_psi_1_int = gsl_sf_psi_1_int(n)
  end function fgsl_sf_psi_1_int
  function fgsl_sf_psi_1_int_e(n, result) 
    integer(c_int), intent(in) :: n
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_psi_1_int_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_psi_1_int_e = gsl_sf_psi_1_int_e(n, res)
    result = res
  end function fgsl_sf_psi_1_int_e
  function fgsl_sf_psi_1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_psi_1
    fgsl_sf_psi_1 = gsl_sf_psi_1(x)
  end function fgsl_sf_psi_1
  function fgsl_sf_psi_1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_psi_1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_psi_1_e = gsl_sf_psi_1_e(x, res)
    result = res
  end function fgsl_sf_psi_1_e
  function fgsl_sf_psi_n(m, x)
    integer(fgsl_int), intent(in) :: m
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_psi_n
    fgsl_sf_psi_n = gsl_sf_psi_n(m, x)
  end function fgsl_sf_psi_n
  function fgsl_sf_psi_n_e(m, x, result) 
    integer(fgsl_int), intent(in) :: m
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_psi_n_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_psi_n_e = gsl_sf_psi_n_e(m, x, res)
    result = res
  end function fgsl_sf_psi_n_e
  function fgsl_sf_psi_1piy(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_psi_1piy
    fgsl_sf_psi_1piy = gsl_sf_psi_1piy(x)
  end function fgsl_sf_psi_1piy
  function fgsl_sf_psi_1piy_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_psi_1piy_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_psi_1piy_e = gsl_sf_psi_1piy_e(x, res)
    result = res
  end function fgsl_sf_psi_1piy_e
  function fgsl_sf_synchrotron_1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_synchrotron_1
    fgsl_sf_synchrotron_1 = gsl_sf_synchrotron_1(x)
  end function fgsl_sf_synchrotron_1
  function fgsl_sf_synchrotron_1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_synchrotron_1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_synchrotron_1_e = gsl_sf_synchrotron_1_e(x, res)
    result = res
  end function fgsl_sf_synchrotron_1_e
  function fgsl_sf_synchrotron_2(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_synchrotron_2
    fgsl_sf_synchrotron_2 = gsl_sf_synchrotron_2(x)
  end function fgsl_sf_synchrotron_2
  function fgsl_sf_synchrotron_2_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_synchrotron_2_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_synchrotron_2_e = gsl_sf_synchrotron_2_e(x, res)
    result = res
  end function fgsl_sf_synchrotron_2_e
  function fgsl_sf_transport_2(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_transport_2
    fgsl_sf_transport_2 = gsl_sf_transport_2(x)
  end function fgsl_sf_transport_2
  function fgsl_sf_transport_2_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_transport_2_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_transport_2_e = gsl_sf_transport_2_e(x, res)
    result = res
  end function fgsl_sf_transport_2_e
  function fgsl_sf_transport_3(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_transport_3
    fgsl_sf_transport_3 = gsl_sf_transport_3(x)
  end function fgsl_sf_transport_3
  function fgsl_sf_transport_3_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_transport_3_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_transport_3_e = gsl_sf_transport_3_e(x, res)
    result = res
  end function fgsl_sf_transport_3_e
  function fgsl_sf_transport_4(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_transport_4
    fgsl_sf_transport_4 = gsl_sf_transport_4(x)
  end function fgsl_sf_transport_4
  function fgsl_sf_transport_4_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_transport_4_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_transport_4_e = gsl_sf_transport_4_e(x, res)
    result = res
  end function fgsl_sf_transport_4_e
  function fgsl_sf_transport_5(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_transport_5
    fgsl_sf_transport_5 = gsl_sf_transport_5(x)
  end function fgsl_sf_transport_5
  function fgsl_sf_transport_5_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_transport_5_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_transport_5_e = gsl_sf_transport_5_e(x, res)
    result = res
  end function fgsl_sf_transport_5_e
  function fgsl_sf_hypot(x, y) 
    real(fgsl_double), intent(in) :: x, y
    real(fgsl_double) :: fgsl_sf_hypot
    fgsl_sf_hypot = gsl_sf_hypot(x, y)
  end function fgsl_sf_hypot
  function fgsl_sf_hypot_e(x, y, result) 
    real(fgsl_double), intent(in) :: x, y
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hypot_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hypot_e = gsl_sf_hypot_e(x, y, res)
    result = res
  end function fgsl_sf_hypot_e
  function fgsl_sf_sinc(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_sinc
    fgsl_sf_sinc = gsl_sf_sinc(x)
  end function fgsl_sf_sinc
  function fgsl_sf_sinc_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_sinc_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_sinc_e = gsl_sf_sinc_e(x, res)
    result = res
  end function fgsl_sf_sinc_e
  function fgsl_sf_complex_sin_e(zr, zi, szr, szi) 
    real(fgsl_double), intent(in) :: zr, zi
    type(fgsl_sf_result), intent(out) :: szr, szi
    integer(fgsl_int) :: fgsl_sf_complex_sin_e
!
    type(gsl_sf_result) :: r_loc, i_loc
    fgsl_sf_complex_sin_e = gsl_sf_complex_sin_e(zr, zi, r_loc, i_loc)
    szr = r_loc
    szi = i_loc
  end function fgsl_sf_complex_sin_e
  function fgsl_sf_complex_cos_e(zr, zi, czr, czi) 
    real(fgsl_double), intent(in) :: zr, zi
    type(fgsl_sf_result), intent(out) :: czr, czi
    integer(fgsl_int) :: fgsl_sf_complex_cos_e
!
    type(gsl_sf_result) :: r_loc, i_loc
    fgsl_sf_complex_cos_e = gsl_sf_complex_cos_e(zr, zi, r_loc, i_loc)
    czr = r_loc
    czi = i_loc
  end function fgsl_sf_complex_cos_e
  function fgsl_sf_complex_logsin_e(zr, zi, lszr, lszi) 
    real(fgsl_double), intent(in) :: zr, zi
    type(fgsl_sf_result), intent(out) :: lszr, lszi
    integer(fgsl_int) :: fgsl_sf_complex_logsin_e
!
    type(gsl_sf_result) :: r_loc, i_loc
    fgsl_sf_complex_logsin_e = gsl_sf_complex_logsin_e(zr, zi, r_loc, i_loc)
    lszr = r_loc
    lszi = i_loc
  end function fgsl_sf_complex_logsin_e
  function fgsl_sf_lnsinh(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_lnsinh
    fgsl_sf_lnsinh = gsl_sf_lnsinh(x)
  end function fgsl_sf_lnsinh
  function fgsl_sf_lnsinh_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_lnsinh_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lnsinh_e = gsl_sf_lnsinh_e(x, res)
    result = res
  end function fgsl_sf_lnsinh_e
  function fgsl_sf_lncosh(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_lncosh
    fgsl_sf_lncosh = gsl_sf_lncosh(x)
  end function fgsl_sf_lncosh
  function fgsl_sf_lncosh_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_lncosh_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_lncosh_e = gsl_sf_lncosh_e(x, res)
    result = res
  end function fgsl_sf_lncosh_e
  function fgsl_sf_polar_to_rect(r, theta, x, y) 
    real(fgsl_double), intent(in) :: r, theta
    type(fgsl_sf_result), intent(out) :: x, y
    integer(fgsl_int) :: fgsl_sf_polar_to_rect
!
    type(gsl_sf_result) :: x_loc, y_loc
    fgsl_sf_polar_to_rect = gsl_sf_polar_to_rect(r, theta, x_loc, y_loc)
    x = x_loc
    y = y_loc
  end function fgsl_sf_polar_to_rect
  function fgsl_sf_rect_to_polar(x, y, r, theta) 
    real(fgsl_double), intent(in) :: x, y
    type(fgsl_sf_result), intent(out) :: r, theta
    integer(fgsl_int) :: fgsl_sf_rect_to_polar
!
    type(gsl_sf_result) :: r_loc, th_loc
    fgsl_sf_rect_to_polar = gsl_sf_rect_to_polar(x, y, r_loc, th_loc)
    r = r_loc
    theta = th_loc
  end function fgsl_sf_rect_to_polar
  function fgsl_sf_angle_restrict_symm(theta) 
    real(fgsl_double), intent(in) :: theta
    real(fgsl_double) :: fgsl_sf_angle_restrict_symm
    fgsl_sf_angle_restrict_symm = gsl_sf_angle_restrict_symm(theta)
  end function fgsl_sf_angle_restrict_symm
  function fgsl_sf_angle_restrict_symm_e(theta) 
    real(fgsl_double), intent(inout) :: theta
    integer(fgsl_int) :: fgsl_sf_angle_restrict_symm_e
!
    fgsl_sf_angle_restrict_symm_e = gsl_sf_angle_restrict_symm_e(theta)
  end function fgsl_sf_angle_restrict_symm_e
  function fgsl_sf_angle_restrict_pos(theta) 
    real(fgsl_double), intent(in) :: theta
    real(fgsl_double) :: fgsl_sf_angle_restrict_pos
    fgsl_sf_angle_restrict_pos = gsl_sf_angle_restrict_pos(theta)
  end function fgsl_sf_angle_restrict_pos
  function fgsl_sf_angle_restrict_pos_e(theta) 
    real(fgsl_double), intent(inout) :: theta
    integer(fgsl_int) :: fgsl_sf_angle_restrict_pos_e
!
    fgsl_sf_angle_restrict_pos_e = gsl_sf_angle_restrict_pos_e(theta)
  end function fgsl_sf_angle_restrict_pos_e
  function fgsl_sf_sin_err_e(x, dx, result) 
    real(fgsl_double), intent(in) :: x, dx
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_sin_err_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_sin_err_e = gsl_sf_sin_err_e(x, dx, res)
    result = res
  end function fgsl_sf_sin_err_e
  function fgsl_sf_cos_err_e(x, dx, result) 
    real(fgsl_double), intent(in) :: x, dx
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_cos_err_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_cos_err_e = gsl_sf_cos_err_e(x, dx, res)
    result = res
  end function fgsl_sf_cos_err_e
  function fgsl_sf_zeta_int(n) 
    integer(c_int), intent(in) :: n
    real(fgsl_double) :: fgsl_sf_zeta_int
    fgsl_sf_zeta_int = gsl_sf_zeta_int(n)
  end function fgsl_sf_zeta_int
  function fgsl_sf_zeta_int_e(n, result) 
    integer(c_int), intent(in) :: n
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_zeta_int_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_zeta_int_e = gsl_sf_zeta_int_e(n, res)
    result = res
  end function fgsl_sf_zeta_int_e
  function fgsl_sf_zeta(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_zeta
    fgsl_sf_zeta = gsl_sf_zeta(x)
  end function fgsl_sf_zeta
  function fgsl_sf_zeta_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_zeta_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_zeta_e = gsl_sf_zeta_e(x, res)
    result = res
  end function fgsl_sf_zeta_e
  function fgsl_sf_zetam1_int(n) 
    integer(c_int), intent(in) :: n
    real(fgsl_double) :: fgsl_sf_zetam1_int
    fgsl_sf_zetam1_int = gsl_sf_zetam1_int(n)
  end function fgsl_sf_zetam1_int
  function fgsl_sf_zetam1_int_e(n, result) 
    integer(c_int), intent(in) :: n
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_zetam1_int_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_zetam1_int_e = gsl_sf_zetam1_int_e(n, res)
    result = res
  end function fgsl_sf_zetam1_int_e
  function fgsl_sf_zetam1(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_zetam1
    fgsl_sf_zetam1 = gsl_sf_zetam1(x)
  end function fgsl_sf_zetam1
  function fgsl_sf_zetam1_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_zetam1_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_zetam1_e = gsl_sf_zetam1_e(x, res)
    result = res
  end function fgsl_sf_zetam1_e
  function fgsl_sf_hzeta(s, q) 
    real(fgsl_double), intent(in) :: s, q
    real(fgsl_double) :: fgsl_sf_hzeta
    fgsl_sf_hzeta = gsl_sf_hzeta(s, q)
  end function fgsl_sf_hzeta
  function fgsl_sf_hzeta_e(s, q, result) 
    real(fgsl_double), intent(in) :: s, q
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_hzeta_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_hzeta_e = gsl_sf_hzeta_e(s, q, res)
    result = res
  end function fgsl_sf_hzeta_e
  function fgsl_sf_eta_int(n) 
    integer(c_int), intent(in) :: n
    real(fgsl_double) :: fgsl_sf_eta_int
    fgsl_sf_eta_int = gsl_sf_eta_int(n)
  end function fgsl_sf_eta_int
  function fgsl_sf_eta_int_e(n, result) 
    integer(c_int), intent(in) :: n
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_eta_int_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_eta_int_e = gsl_sf_eta_int_e(n, res)
    result = res
  end function fgsl_sf_eta_int_e
  function fgsl_sf_eta(x) 
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_sf_eta
    fgsl_sf_eta = gsl_sf_eta(x)
  end function fgsl_sf_eta
  function fgsl_sf_eta_e(x, result) 
    real(fgsl_double), intent(in) :: x
    type(fgsl_sf_result), intent(out) :: result
    integer(fgsl_int) :: fgsl_sf_eta_e
!
    type(gsl_sf_result) :: res
    fgsl_sf_eta_e = gsl_sf_eta_e(x, res)
    result = res
  end function fgsl_sf_eta_e
!
  elemental subroutine gsl_sf_to_fgsl_sf(result, source)
    type(fgsl_sf_result), intent(out) :: result
    type(gsl_sf_result), intent(in) :: source
    result%val = source%val
    result%err = source%err
  end subroutine gsl_sf_to_fgsl_sf
  elemental subroutine gsl_sfe10_to_fgsl_sfe10(result, source)
    type(fgsl_sf_result_e10), intent(out) :: result
    type(gsl_sf_result_e10), intent(in) :: source
    result%val = source%val
    result%err = source%err
    result%e10 = source%e10
  end subroutine gsl_sfe10_to_fgsl_sfe10
!-*-f90-*-
!
!  API: Array support
!
! vectors (real)
!
  function fgsl_vector_init(type)
    real(fgsl_double), intent(in) :: type
    type(fgsl_vector) :: fgsl_vector_init
    fgsl_vector_init%gsl_vector = fgsl_aux_vector_double_init()
  end function fgsl_vector_init
! Note: The following requires the actual argument
! array having the TARGET attribute. Otherwise
! being passed by reference is not guaranteed by the Fortran standard.
  function fgsl_vector_align(array, len, fvec, size, offset, stride)
    integer(fgsl_size_t), intent(in) :: len, size, offset, stride
    real(fgsl_double), dimension(len), target, intent(in) :: array
    type(fgsl_vector), intent(inout) :: fvec
    integer(fgsl_int) :: fgsl_vector_align
!
    fgsl_vector_align = fgsl_aux_vector_double_align(array, len, &
         fvec%gsl_vector, size, offset, stride)
  end function fgsl_vector_align
  function fgsl_vector_pointer_align(ptr, fvec)
    real(fgsl_double), pointer, intent(out) :: ptr(:)
    type(fgsl_vector), intent(in) :: fvec
    integer(fgsl_int) :: fgsl_vector_pointer_align
!
    real(fgsl_double), pointer :: fp_local(:)
    type(c_ptr) :: cp
    integer(fgsl_size_t) :: size, stride
! tests
!    real(fgsl_double) :: cc(3)
    size = fgsl_aux_vector_double_size(fvec%gsl_vector)
    stride = fgsl_aux_vector_double_stride(fvec%gsl_vector)
    if (stride == 0) then
       fgsl_vector_pointer_align = fgsl_einval
    else
       cp = gsl_vector_ptr(fvec%gsl_vector,0_fgsl_size_t)
!       cc(1) = gsl_vector_get(fvec%gsl_vector,0_c_size_t)
!       cc(2) = gsl_vector_get(fvec%gsl_vector,1_c_size_t)
!       cc(3) = gsl_vector_get(fvec%gsl_vector,2_c_size_t)
       call c_f_pointer(cp, fp_local, (/ size*stride /))
!       write(6, *) 'size, stride, fp_local: ',size,stride,fp_local(1:3),cc(1:3)
       ptr => fp_local(1:size*stride:stride)
       fgsl_vector_pointer_align = fgsl_success
    end if
  end function fgsl_vector_pointer_align
  subroutine fgsl_vector_to_array(result, source)
    real(fgsl_double), intent(inout) :: result(:)
    type(fgsl_vector), intent(in) :: source
!
    integer(fgsl_size_t) :: i, n, k
    k = size(result)
    n = min(k,fgsl_aux_vector_double_size(source%gsl_vector))
!    write(6,*) 'result length: ',size(result)
!    write(6,*) 'vector length: ', &
!    fgsl_aux_vector_double_size(source%gsl_vector)
    do i=1,n
       result(i) = gsl_vector_get(source%gsl_vector,i-1)
    end do
    do i=n+1,size(result)
       result(i) = 0.0_fgsl_double
    end do
  end subroutine fgsl_vector_to_array
  subroutine fgsl_vector_free(fvec)
    type(fgsl_vector), intent(inout) :: fvec
!    call gsl_vector_free(fvec%gsl_vector)
    call fgsl_aux_vector_double_free(fvec%gsl_vector)
  end subroutine fgsl_vector_free
  subroutine fgsl_vector_c_ptr(res, src) 
    type(c_ptr), intent(in) :: src
    type(fgsl_vector), intent(out) :: res
    res%gsl_vector = src
  end subroutine fgsl_vector_c_ptr
  function fgsl_vector_status(vector)
    type(fgsl_vector), intent(in) :: vector
    logical :: fgsl_vector_status
    fgsl_vector_status = .true.
    if (.not. c_associated(vector%gsl_vector)) fgsl_vector_status = .false.
  end function fgsl_vector_status
  function fgsl_sizeof_vector(w)
    type(fgsl_vector), intent(in) :: w
    integer(fgsl_size_t) :: fgsl_sizeof_vector
    fgsl_sizeof_vector = gsl_aux_sizeof_vector()
  end function fgsl_sizeof_vector
!
! vectors (complex)
!
  function fgsl_vector_complex_init(type)
    complex(fgsl_double_complex), intent(in) :: type
    type(fgsl_vector_complex) :: fgsl_vector_complex_init
    fgsl_vector_complex_init%gsl_vector_complex = fgsl_aux_vector_complex_init()
  end function fgsl_vector_complex_init
! Note: The following requires the actual argument
! array having the TARGET attribute. Otherwise
! being passed by reference is not guaranteed by the Fortran standard.
  function fgsl_vector_complex_align(array, len, fvec, size, offset, stride)
    integer(fgsl_size_t), intent(in) :: len, size, offset, stride
    complex(fgsl_double_complex), dimension(len), target, intent(in) :: array
    type(fgsl_vector_complex), intent(inout) :: fvec
    integer(fgsl_int) :: fgsl_vector_complex_align
!
    fgsl_vector_complex_align = &
         fgsl_aux_vector_complex_align(c_loc(array), len, &
         fvec%gsl_vector_complex, size, offset, stride)
  end function fgsl_vector_complex_align
  function fgsl_vector_complex_pointer_align(ptr, fvec)
    complex(fgsl_double_complex), pointer, intent(out) :: ptr(:)
    type(fgsl_vector_complex), intent(in) :: fvec
    integer(fgsl_int) :: fgsl_vector_complex_pointer_align
!
    complex(fgsl_double_complex), pointer :: fp_local(:)
    type(c_ptr) :: cp
    integer(fgsl_size_t) :: size, stride
! tests
!    real(fgsl_double) :: cc(3)
    size = fgsl_aux_vector_complex_size(fvec%gsl_vector_complex)
    stride = fgsl_aux_vector_complex_stride(fvec%gsl_vector_complex)
    if (stride == 0) then
       fgsl_vector_complex_pointer_align = fgsl_einval
    else
       cp = gsl_vector_complex_ptr(fvec%gsl_vector_complex,0_fgsl_size_t)
!       cc(1) = gsl_vector_complex_get(fvec%gsl_vector_complex,0_c_size_t)
!       cc(2) = gsl_vector_complex_get(fvec%gsl_vector_complex,1_c_size_t)
!       cc(3) = gsl_vector_complex_get(fvec%gsl_vector_complex,2_c_size_t)
       call c_f_pointer(cp, fp_local, (/ size*stride /))
!       write(6, *) 'size, stride, fp_local: ',size,stride,fp_local(1:3),cc(1:3)
       ptr => fp_local(1:size*stride:stride)
       fgsl_vector_complex_pointer_align = fgsl_success
    end if
  end function fgsl_vector_complex_pointer_align
  subroutine fgsl_vector_complex_to_array(result, source)
    complex(fgsl_double_complex), intent(inout) :: result(:)
    type(fgsl_vector_complex), intent(in) :: source
!    type(gsl_complex) :: aux
!
    integer(fgsl_size_t) :: i, n, k
    k = size(result)
    n = min(k,fgsl_aux_vector_complex_size(source%gsl_vector_complex))
!    write(6,*) 'result length: ',size(result)
!    write(6,*) 'vector_complex length: ', &
!         fgsl_aux_vector_complex_size(source%gsl_vector_complex)
    do i=1,n
       result(i) = gsl_vector_complex_get(source%gsl_vector_complex,i-1)
!       aux = gsl_vector_complex_get(source%gsl_vector_complex,i-1)
!       result(i) = aux
!       write(6, *) 'i=',i,' res = ',result(i)
    end do
    do i=n+1,size(result)
       result(i) = 0.0_fgsl_double
    end do
  end subroutine fgsl_vector_complex_to_array
  subroutine fgsl_vector_complex_free(fvec)
    type(fgsl_vector_complex), intent(inout) :: fvec
!    call gsl_vector_complex_free(fvec%gsl_vector_complex)
    call fgsl_aux_vector_complex_free(fvec%gsl_vector_complex)
  end subroutine fgsl_vector_complex_free
  subroutine fgsl_vector_complex_c_ptr(res, src) 
    type(c_ptr), intent(in) :: src
    type(fgsl_vector_complex), intent(out) :: res
    res%gsl_vector_complex = src
  end subroutine fgsl_vector_complex_c_ptr
  function fgsl_vector_complex_status(vector_complex)
    type(fgsl_vector_complex), intent(in) :: vector_complex
    logical :: fgsl_vector_complex_status
    fgsl_vector_complex_status = .true.
    if (.not. c_associated(vector_complex%gsl_vector_complex)) fgsl_vector_complex_status = .false.
  end function fgsl_vector_complex_status
  function fgsl_sizeof_vector_complex(w)
    type(fgsl_vector_complex), intent(in) :: w
    integer(fgsl_size_t) :: fgsl_sizeof_vector_complex
    fgsl_sizeof_vector_complex = gsl_aux_sizeof_vector_complex()
  end function fgsl_sizeof_vector_complex

!
! matrices (real)
!
  function fgsl_matrix_init(type)
    real(fgsl_double), intent(in) :: type
    type(fgsl_matrix) :: fgsl_matrix_init
    fgsl_matrix_init%gsl_matrix = fgsl_aux_matrix_double_init()
  end function fgsl_matrix_init
! here again, TARGET needed for actual argument
  function fgsl_matrix_align(array, lda, n, m, fmat)
    integer(fgsl_size_t), intent(in) :: lda, n, m
    real(fgsl_double), dimension(lda, m), target, intent(in) :: array
    type(fgsl_matrix), intent(inout) :: fmat
    integer(fgsl_int) :: fgsl_matrix_align
!
    fgsl_matrix_align = fgsl_aux_matrix_double_align(c_loc(array), lda, &
         n, m, fmat%gsl_matrix)
  end function fgsl_matrix_align
  function fgsl_matrix_pointer_align(ptr, fmat)
    real(fgsl_double), pointer, intent(out) :: ptr(:,:)
    type(fgsl_matrix), intent(in) :: fmat
    integer(fgsl_int) :: fgsl_matrix_pointer_align
!
    real(fgsl_double), pointer :: fp_local(:,:)
    type(c_ptr) :: cp
    integer(fgsl_size_t) :: m, n, lda
    call fgsl_aux_matrix_double_size(fmat%gsl_matrix, lda, n, m)
    cp = gsl_matrix_ptr(fmat%gsl_matrix,0_fgsl_size_t,0_fgsl_size_t)
    call c_f_pointer(cp, fp_local, (/ lda , m /))
    ptr => fp_local(1:n,1:m)
    fgsl_matrix_pointer_align = fgsl_success
  end function fgsl_matrix_pointer_align
  subroutine fgsl_matrix_to_array(result, source)
    real(fgsl_double), intent(inout) :: result(:,:)
    type(fgsl_matrix), intent(in) :: source
!
    integer(fgsl_size_t) :: i, j, kl, m, n, ml, nl, lda
    call fgsl_aux_matrix_double_size(source%gsl_matrix, lda, n, m)
    
    kl = size(result,1)
    nl = min(kl,n)
    kl = size(result,2)
    ml = min(kl,m)
!    write(6, *) 'Number of rows: ', nl, n
!    write(6, *) 'Number of cols: ', ml, m
    do j=1,ml
       do i=1,nl
          result(i,j) = gsl_matrix_get(source%gsl_matrix,j-1,i-1)
       end do
    end do
    do j=1,ml
       do i=nl+1,size(result,1)
          result(i,j) = 0.0_fgsl_double
       end do
    end do
    do j=ml+1,size(result,2)
       do i=1,size(result,1)
          result(i,j) = 0.0_fgsl_double
       end do
    end do
  end subroutine fgsl_matrix_to_array
  subroutine fgsl_matrix_free(fvec)
    type(fgsl_matrix), intent(inout) :: fvec
    call fgsl_aux_matrix_double_free(fvec%gsl_matrix)
  end subroutine fgsl_matrix_free
  subroutine fgsl_matrix_c_ptr(res, src) 
    type(c_ptr), intent(in) :: src
    type(fgsl_matrix), intent(out) :: res
    res%gsl_matrix = src
  end subroutine fgsl_matrix_c_ptr
  function fgsl_matrix_status(matrix)
    type(fgsl_matrix), intent(in) :: matrix
    logical :: fgsl_matrix_status
    fgsl_matrix_status = .true.
    if (.not. c_associated(matrix%gsl_matrix)) fgsl_matrix_status = .false.
  end function fgsl_matrix_status
  function fgsl_sizeof_matrix(w)
    type(fgsl_matrix), intent(in) :: w
    integer(fgsl_size_t) :: fgsl_sizeof_matrix
    fgsl_sizeof_matrix = gsl_aux_sizeof_matrix()
  end function fgsl_sizeof_matrix
!
! matrices (complex)
!
  function fgsl_matrix_complex_init(type)
    complex(fgsl_double_complex), intent(in) :: type
    type(fgsl_matrix_complex) :: fgsl_matrix_complex_init
    fgsl_matrix_complex_init%gsl_matrix_complex = fgsl_aux_matrix_complex_init()
  end function fgsl_matrix_complex_init
! here again, TARGET needed for actual argument
  function fgsl_matrix_complex_align(array, lda, n, m, fmat)
    integer(fgsl_size_t), intent(in) :: lda, n, m
    complex(fgsl_double_complex), dimension(lda, m), target, intent(in) :: array
    type(fgsl_matrix_complex), intent(inout) :: fmat
    integer(fgsl_int) :: fgsl_matrix_complex_align
!
    fgsl_matrix_complex_align = &
         fgsl_aux_matrix_complex_align(c_loc(array), lda, &
         n, m, fmat%gsl_matrix_complex)
  end function fgsl_matrix_complex_align
  function fgsl_matrix_complex_pointer_align(ptr, fmat)
    complex(fgsl_double_complex), pointer, intent(out) :: ptr(:,:)
    type(fgsl_matrix_complex), intent(in) :: fmat
    integer(fgsl_int) :: fgsl_matrix_complex_pointer_align
!
    complex(fgsl_double_complex), pointer :: fp_local(:,:)
    type(c_ptr) :: cp
    integer(fgsl_size_t) :: m, n, lda
    call fgsl_aux_matrix_complex_size(fmat%gsl_matrix_complex, lda, n, m)
    cp = gsl_matrix_complex_ptr(fmat%gsl_matrix_complex,0_fgsl_size_t,0_fgsl_size_t)
    call c_f_pointer(cp, fp_local, (/ lda , m /))
    ptr => fp_local(1:n,1:m)
    fgsl_matrix_complex_pointer_align = fgsl_success
  end function fgsl_matrix_complex_pointer_align
  subroutine fgsl_matrix_complex_to_array(result, source)
    complex(fgsl_double_complex), intent(inout) :: result(:,:)
    type(fgsl_matrix_complex), intent(in) :: source
!
    integer(fgsl_size_t) :: i, j, kl, m, n, ml, nl, lda
    call fgsl_aux_matrix_complex_size(source%gsl_matrix_complex, lda, n, m)
    
    kl = size(result,1)
    nl = min(kl,n)
    kl = size(result,2)
    ml = min(kl,m)
!    write(6, *) 'Number of rows: ', nl, n
!    write(6, *) 'Number of cols: ', ml, m
    do j=1,ml
       do i=1,nl
          result(i,j) = gsl_matrix_complex_get(source%gsl_matrix_complex,j-1,i-1)
       end do
    end do
    do j=1,ml
       do i=nl+1,size(result,1)
          result(i,j) = 0.0_fgsl_double
       end do
    end do
    do j=ml+1,size(result,2)
       do i=1,size(result,1)
          result(i,j) = 0.0_fgsl_double
       end do
    end do
  end subroutine fgsl_matrix_complex_to_array
  subroutine fgsl_matrix_complex_free(fvec)
    type(fgsl_matrix_complex), intent(inout) :: fvec
    call fgsl_aux_matrix_complex_free(fvec%gsl_matrix_complex)
  end subroutine fgsl_matrix_complex_free
  subroutine fgsl_matrix_complex_c_ptr(res, src) 
    type(c_ptr), intent(in) :: src
    type(fgsl_matrix_complex), intent(out) :: res
    res%gsl_matrix_complex = src
  end subroutine fgsl_matrix_complex_c_ptr
  function fgsl_matrix_complex_status(matrix_complex)
    type(fgsl_matrix_complex), intent(in) :: matrix_complex
    logical :: fgsl_matrix_complex_status
    fgsl_matrix_complex_status = .true.
    if (.not. c_associated(matrix_complex%gsl_matrix_complex)) fgsl_matrix_complex_status = .false.
  end function fgsl_matrix_complex_status
  function fgsl_sizeof_matrix_complex(w)
    type(fgsl_matrix_complex), intent(in) :: w
    integer(fgsl_size_t) :: fgsl_sizeof_matrix_complex
    fgsl_sizeof_matrix_complex = gsl_aux_sizeof_matrix_complex()
  end function fgsl_sizeof_matrix_complex
!-*-f90-*-
!
! API: Interpolation
!
  function fgsl_interp_alloc(interp_type, size)
    type(fgsl_interp_type), intent(in) :: interp_type
    integer(fgsl_size_t), intent(in) :: size
    type(fgsl_interp) :: fgsl_interp_alloc
!
    type(c_ptr) :: it
    integer(c_int) :: i
    integer(c_size_t) :: sz
!
    i = interp_type%which
    sz = size
    it = fgsl_aux_interp_alloc(i)
!       write(6, *) 'DEBUG: alloc size is ',sz
       fgsl_interp_alloc%gsl_interp = gsl_interp_alloc(it, sz)
  end function fgsl_interp_alloc
  subroutine fgsl_interp_free(interp)
    type(fgsl_interp), intent(inout) :: interp
!
    call gsl_interp_free(interp%gsl_interp)
  end subroutine fgsl_interp_free
  function fgsl_interp_accel_alloc()
    type(fgsl_interp_accel) :: fgsl_interp_accel_alloc
!
    fgsl_interp_accel_alloc%gsl_interp_accel = gsl_interp_accel_alloc()
  end function fgsl_interp_accel_alloc
  subroutine fgsl_interp_accel_free(acc)
    type(fgsl_interp_accel), intent(inout) :: acc
!
    call gsl_interp_accel_free(acc%gsl_interp_accel)
  end subroutine fgsl_interp_accel_free
  function fgsl_interp_status(interp)
    type(fgsl_interp), intent(in) :: interp
    logical :: fgsl_interp_status
    fgsl_interp_status = .true.
    if (.not. c_associated(interp%gsl_interp)) fgsl_interp_status = .false.
  end function fgsl_interp_status
  function fgsl_interp_accel_status(acc)
    type(fgsl_interp_accel), intent(in) :: acc
    logical :: fgsl_interp_accel_status
    fgsl_interp_accel_status = .true.
    if (.not. c_associated(acc%gsl_interp_accel)) fgsl_interp_accel_status = .false.
  end function fgsl_interp_accel_status
  function fgsl_interp_init(interp, xa, ya, size)
    type(fgsl_interp), intent(inout) :: interp
    integer(fgsl_size_t), intent(in) :: size
    real(fgsl_double), intent(in) :: xa(size), ya(size)
    integer(fgsl_int) :: fgsl_interp_init
    fgsl_interp_init = gsl_interp_init(interp%gsl_interp, xa, ya, size)
  end function fgsl_interp_init
  function fgsl_interp_eval(interp, xa, ya, x, acc)
    type(fgsl_interp), intent(in) :: interp
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) :: xa(:), ya(:), x
    real(fgsl_double) :: fgsl_interp_eval
    fgsl_interp_eval = gsl_interp_eval(interp%gsl_interp, xa, ya, x, &
         acc%gsl_interp_accel)
  end function fgsl_interp_eval
  function fgsl_interp_eval_e(interp, xa, ya, x, acc, y)
    type(fgsl_interp), intent(in) :: interp
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) :: xa(:), ya(:), x
    real(fgsl_double), intent(out) :: y
    integer(fgsl_int) :: fgsl_interp_eval_e
    fgsl_interp_eval_e = gsl_interp_eval_e(interp%gsl_interp, xa, ya, x, &
         acc%gsl_interp_accel, y)
  end function fgsl_interp_eval_e
  function fgsl_interp_eval_integ(interp, xa, ya, a, b, acc)
    type(fgsl_interp), intent(in) :: interp
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) :: xa(:), ya(:), a, b
    real(fgsl_double) :: fgsl_interp_eval_integ
    fgsl_interp_eval_integ = gsl_interp_eval_integ(interp%gsl_interp, xa, ya, a, b, &
         acc%gsl_interp_accel)
  end function fgsl_interp_eval_integ
  function fgsl_interp_eval_integ_e(interp, xa, ya, a, b, acc, result)
    type(fgsl_interp), intent(in) :: interp
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) :: xa(:), ya(:), a, b
    real(fgsl_double), intent(out) :: result
    integer(fgsl_int) :: fgsl_interp_eval_integ_e
    fgsl_interp_eval_integ_e = gsl_interp_eval_integ_e(interp%gsl_interp, xa, ya, a, b, &
         acc%gsl_interp_accel, result)
  end function fgsl_interp_eval_integ_e
  function fgsl_interp_eval_deriv(interp, xa, ya, x, acc)
    type(fgsl_interp), intent(in) :: interp
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) :: xa(:), ya(:), x
    real(fgsl_double) :: fgsl_interp_eval_deriv
    fgsl_interp_eval_deriv = gsl_interp_eval_deriv(interp%gsl_interp, xa, ya, x, &
         acc%gsl_interp_accel)
  end function fgsl_interp_eval_deriv
  function fgsl_interp_eval_deriv_e(interp, xa, ya, x, acc, d)
    type(fgsl_interp), intent(in) :: interp
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) :: xa(:), ya(:), x
    real(fgsl_double), intent(out) :: d
    integer(fgsl_int) :: fgsl_interp_eval_deriv_e
    fgsl_interp_eval_deriv_e = gsl_interp_eval_deriv_e(interp%gsl_interp, xa, ya, x, &
         acc%gsl_interp_accel, d)
  end function fgsl_interp_eval_deriv_e
  function fgsl_interp_eval_deriv2(interp, xa, ya, x, acc)
    type(fgsl_interp), intent(in) :: interp
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) :: xa(:), ya(:), x
    real(fgsl_double) :: fgsl_interp_eval_deriv2
    fgsl_interp_eval_deriv2 = gsl_interp_eval_deriv2(interp%gsl_interp, xa, ya, x, &
         acc%gsl_interp_accel)
  end function fgsl_interp_eval_deriv2
  function fgsl_interp_eval_deriv2_e(interp, xa, ya, x, acc, d2)
    type(fgsl_interp), intent(in) :: interp
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) :: xa(:), ya(:), x
    real(fgsl_double), intent(out) :: d2
    integer(fgsl_int) :: fgsl_interp_eval_deriv2_e
    fgsl_interp_eval_deriv2_e = gsl_interp_eval_deriv2_e(interp%gsl_interp, xa, ya, x, &
         acc%gsl_interp_accel, d2)
  end function fgsl_interp_eval_deriv2_e
  function fgsl_interp_name(interp)
    type(fgsl_interp), intent(in) :: interp
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_interp_name
!
    type(c_ptr) :: name
!
    name = gsl_interp_name(interp%gsl_interp)
    fgsl_interp_name = fgsl_name(name)
  end function fgsl_interp_name
  function fgsl_interp_min_size(interp)
    type(fgsl_interp), intent(in) :: interp
    integer(fgsl_long) :: fgsl_interp_min_size
!
    fgsl_interp_min_size = fgsl_aux_interp_min_size(interp%gsl_interp)
  end function fgsl_interp_min_size
  function fgsl_interp_bsearch(xa, x, index_lo, index_hi)
    real(fgsl_double), intent(in) :: xa(:)
    real(fgsl_double), intent(in) :: x
    integer(fgsl_size_t), intent(in) :: index_lo, index_hi
    integer(fgsl_size_t) ::  fgsl_interp_bsearch
    fgsl_interp_bsearch = gsl_interp_bsearch(xa, x, index_lo, index_hi) + 1
  end function fgsl_interp_bsearch
  function fgsl_interp_accel_find(acc, xa, size, x)
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), dimension(*), intent(in) :: xa
    integer(fgsl_size_t), intent(in) :: size
    real(fgsl_double), intent(in) :: x
    integer(fgsl_size_t) ::  fgsl_interp_accel_find
    fgsl_interp_accel_find = gsl_interp_accel_find(acc%gsl_interp_accel, xa, size, x) + 1
  end function fgsl_interp_accel_find
  function fgsl_spline_alloc(interp_type, size)
    type(fgsl_interp_type), intent(in) :: interp_type
    integer(fgsl_size_t), intent(in) :: size
    type(fgsl_spline) :: fgsl_spline_alloc
!
    type(c_ptr) :: it
    integer(c_int) :: i
    integer(c_size_t) :: sz
!
    i = interp_type%which
    sz = size
    it = fgsl_aux_interp_alloc(i)
    fgsl_spline_alloc%gsl_spline = gsl_spline_alloc(it, sz)
  end function fgsl_spline_alloc
  subroutine fgsl_spline_free(spline)
    type(fgsl_spline), intent(inout) :: spline
!
    call gsl_spline_free(spline%gsl_spline)
  end subroutine fgsl_spline_free
  function fgsl_spline_init(spline, xa, ya, size)
    type(fgsl_spline), intent(inout) :: spline
    integer(fgsl_size_t), intent(in) :: size
    real(fgsl_double), intent(in) :: xa(size), ya(size)
    integer(fgsl_int) :: fgsl_spline_init
    fgsl_spline_init = gsl_spline_init(spline%gsl_spline, xa, ya, size)
  end function fgsl_spline_init
  function fgsl_spline_name(spline)
    type(fgsl_spline), intent(in) :: spline
    character(len=fgsl_strmax) :: fgsl_spline_name
!
    type(c_ptr) :: name
!
    name = gsl_spline_name(spline%gsl_spline)
    fgsl_spline_name = fgsl_name(name)
  end function fgsl_spline_name
  function fgsl_spline_min_size(spline)
    type(fgsl_spline), intent(in) :: spline
    integer(fgsl_long) :: fgsl_spline_min_size
!
    fgsl_spline_min_size = fgsl_aux_spline_min_size(spline%gsl_spline)
  end function fgsl_spline_min_size
  function fgsl_spline_eval(spline, x, acc)
    type(fgsl_spline), intent(in) :: spline
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) ::  x
    real(fgsl_double) :: fgsl_spline_eval
    fgsl_spline_eval = gsl_spline_eval(spline%gsl_spline, x, &
         acc%gsl_interp_accel)
  end function fgsl_spline_eval
  function fgsl_spline_eval_e(spline, x, acc, y)
    type(fgsl_spline), intent(in) :: spline
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) ::  x
    real(fgsl_double), intent(out) :: y
    integer(fgsl_int) :: fgsl_spline_eval_e
    fgsl_spline_eval_e = gsl_spline_eval_e(spline%gsl_spline, x, &
         acc%gsl_interp_accel, y)
  end function fgsl_spline_eval_e
  function fgsl_spline_eval_deriv(spline, x, acc)
    type(fgsl_spline), intent(in) :: spline
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) ::  x
    real(fgsl_double) :: fgsl_spline_eval_deriv
    fgsl_spline_eval_deriv = gsl_spline_eval_deriv(spline%gsl_spline, x, &
         acc%gsl_interp_accel)
  end function fgsl_spline_eval_deriv
  function fgsl_spline_eval_deriv_e(spline, x, acc, y)
    type(fgsl_spline), intent(in) :: spline
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) ::  x
    real(fgsl_double), intent(out) :: y
    integer(fgsl_int) :: fgsl_spline_eval_deriv_e
    fgsl_spline_eval_deriv_e = gsl_spline_eval_deriv_e(spline%gsl_spline, x, &
         acc%gsl_interp_accel, y)
  end function fgsl_spline_eval_deriv_e
  function fgsl_spline_eval_deriv2(spline, x, acc)
    type(fgsl_spline), intent(in) :: spline
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) ::  x
    real(fgsl_double) :: fgsl_spline_eval_deriv2
    fgsl_spline_eval_deriv2 = gsl_spline_eval_deriv2(spline%gsl_spline, x, &
         acc%gsl_interp_accel)
  end function fgsl_spline_eval_deriv2
  function fgsl_spline_eval_deriv2_e(spline, x, acc, y)
    type(fgsl_spline), intent(in) :: spline
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) ::  x
    real(fgsl_double), intent(out) :: y
    integer(fgsl_int) :: fgsl_spline_eval_deriv2_e
    fgsl_spline_eval_deriv2_e = gsl_spline_eval_deriv2_e(spline%gsl_spline, x, &
         acc%gsl_interp_accel, y)
  end function fgsl_spline_eval_deriv2_e
  function fgsl_spline_eval_integ(spline, a, b, acc)
    type(fgsl_spline), intent(in) :: spline
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_spline_eval_integ
    fgsl_spline_eval_integ = gsl_spline_eval_integ(spline%gsl_spline, a, b, &
         acc%gsl_interp_accel)
  end function fgsl_spline_eval_integ
  function fgsl_spline_eval_integ_e(spline, a, b, acc, y)
    type(fgsl_spline), intent(in) :: spline
    type(fgsl_interp_accel), intent(inout) :: acc
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double), intent(out) :: y
    integer(fgsl_int) :: fgsl_spline_eval_integ_e
    fgsl_spline_eval_integ_e = gsl_spline_eval_integ_e(spline%gsl_spline, a, b, &
         acc%gsl_interp_accel, y)
  end function fgsl_spline_eval_integ_e
  function fgsl_spline_status(spline)
    type(fgsl_spline), intent(in) :: spline
    logical :: fgsl_spline_status
    fgsl_spline_status = .true.
    if (.not. c_associated(spline%gsl_spline)) fgsl_spline_status = .false.
  end function fgsl_spline_status
  function fgsl_sizeof_interp(w)
    type(fgsl_interp), intent(in) :: w
    integer(fgsl_size_t) :: fgsl_sizeof_interp
    fgsl_sizeof_interp = gsl_aux_sizeof_interp()
  end function fgsl_sizeof_interp
!-*-f90-*-
!
! API: Permutations and Combinations
!
  function fgsl_permutation_alloc(n)
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_permutation) :: fgsl_permutation_alloc
    fgsl_permutation_alloc%gsl_permutation = gsl_permutation_alloc(n)
  end function fgsl_permutation_alloc
  function fgsl_permutation_calloc(n)
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_permutation) :: fgsl_permutation_calloc
    fgsl_permutation_calloc%gsl_permutation = gsl_permutation_calloc(n)
  end function fgsl_permutation_calloc
  subroutine fgsl_permutation_init(p) 
    type(fgsl_permutation), intent(inout) :: p
    call gsl_permutation_init(p%gsl_permutation)
  end subroutine fgsl_permutation_init
  subroutine fgsl_permutation_free(p) 
    type(fgsl_permutation), intent(inout) :: p
    call gsl_permutation_free(p%gsl_permutation)
  end subroutine fgsl_permutation_free
  function fgsl_permutation_memcpy(dest, src)
    type(fgsl_permutation), intent(inout) :: dest
    type(fgsl_permutation),  intent(in) :: src
    integer(fgsl_int) :: fgsl_permutation_memcpy
    fgsl_permutation_memcpy = gsl_permutation_memcpy(dest%gsl_permutation, src%gsl_permutation)
  end function fgsl_permutation_memcpy
  function fgsl_permutation_get(p, i) 
    type(fgsl_permutation), intent(inout) :: p
    integer(fgsl_size_t), intent(in) :: i
    integer(fgsl_size_t) :: fgsl_permutation_get
    fgsl_permutation_get = gsl_permutation_get(p%gsl_permutation, i)
  end function fgsl_permutation_get
  function fgsl_permutation_swap(p, i, j) 
    type(fgsl_permutation), intent(inout) :: p
    integer(fgsl_size_t), intent(in) :: i, j
    integer(fgsl_int) :: fgsl_permutation_swap
    fgsl_permutation_swap = gsl_permutation_swap(p%gsl_permutation, i, j) 
  end function fgsl_permutation_swap
  function fgsl_permutation_size(p)
    type(fgsl_permutation), intent(in) :: p
    integer(fgsl_size_t) :: fgsl_permutation_size
    fgsl_permutation_size = gsl_permutation_size(p%gsl_permutation)
  end function fgsl_permutation_size
  function fgsl_permutation_data(p)
    type(fgsl_permutation), intent(in) :: p
    integer(fgsl_size_t), pointer :: fgsl_permutation_data(:)
!
    integer(fgsl_size_t) :: size
    type(c_ptr) :: pdata
    size = gsl_permutation_size(p%gsl_permutation)
    pdata = gsl_permutation_data(p%gsl_permutation)
    if (c_associated(pdata)) then
       call c_f_pointer(pdata,fgsl_permutation_data,(/size/))
    else
       nullify(fgsl_permutation_data)
    end if
  end function fgsl_permutation_data
  function fgsl_permutation_valid(p) 
    type(fgsl_permutation), intent(in) :: p
    integer(fgsl_int) :: fgsl_permutation_valid
    fgsl_permutation_valid = gsl_permutation_valid(p%gsl_permutation)
  end function fgsl_permutation_valid
  subroutine fgsl_permutation_reverse(p) 
    type(fgsl_permutation), intent(inout) :: p
    call gsl_permutation_reverse(p%gsl_permutation)
  end subroutine fgsl_permutation_reverse
  function fgsl_permutation_inverse(inv, p)
    type(fgsl_permutation), intent(inout) :: inv
    type(fgsl_permutation),  intent(in) :: p
    integer(fgsl_int) :: fgsl_permutation_inverse
    fgsl_permutation_inverse = gsl_permutation_inverse(inv%gsl_permutation, p%gsl_permutation)
  end function fgsl_permutation_inverse
  function fgsl_permutation_next(p)
    type(fgsl_permutation), intent(in) :: p
    integer(fgsl_int) :: fgsl_permutation_next
    fgsl_permutation_next = gsl_permutation_next(p%gsl_permutation)
  end function fgsl_permutation_next
  function fgsl_permutation_prev(p)
    type(fgsl_permutation), intent(in) :: p
    integer(fgsl_int) :: fgsl_permutation_prev
    fgsl_permutation_prev = gsl_permutation_prev(p%gsl_permutation)
  end function fgsl_permutation_prev
  function fgsl_permute(p, data, stride, n) 
    integer(fgsl_size_t), intent(in) :: p(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double), dimension(:), intent(inout) :: data
    integer(fgsl_int) :: fgsl_permute
    fgsl_permute = gsl_permute(p, data, stride, n)
  end function fgsl_permute
  function fgsl_permute_long(p, data, stride, n) 
    integer(fgsl_size_t), intent(in) :: p(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    integer(fgsl_long), dimension(:), intent(inout) :: data
    integer(fgsl_int) :: fgsl_permute_long
    fgsl_permute_long = gsl_permute_long(p, data, stride, n)
  end function fgsl_permute_long
  function fgsl_permute_inverse(p, data, stride, n) 
    integer(fgsl_size_t), intent(in) :: p(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double), dimension(:), intent(inout) :: data
    integer(fgsl_int) :: fgsl_permute_inverse
    fgsl_permute_inverse = gsl_permute_inverse(p, data, stride, n)
  end function fgsl_permute_inverse
  function fgsl_permute_long_inverse(p, data, stride, n) 
    integer(fgsl_size_t), intent(in) :: p(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    integer(fgsl_long), dimension(:), intent(inout) :: data
    integer(fgsl_int) :: fgsl_permute_long_inverse
    fgsl_permute_long_inverse = gsl_permute_long_inverse(p, data, stride, n)
  end function fgsl_permute_long_inverse
  function fgsl_permute_vector(p,v) 
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector), intent(inout) :: v
    integer(fgsl_int) :: fgsl_permute_vector
    fgsl_permute_vector = &
         gsl_permute_vector(p%gsl_permutation,v%gsl_vector)
  end function fgsl_permute_vector
  function fgsl_permute_vector_inverse(p,v) 
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector), intent(inout) :: v
    integer(fgsl_int) :: fgsl_permute_vector_inverse
    fgsl_permute_vector_inverse = &
         gsl_permute_vector_inverse(p%gsl_permutation,v%gsl_vector)
  end function fgsl_permute_vector_inverse
  function fgsl_permutation_mul(p, pa, pb)
    type(fgsl_permutation), intent(inout) :: p
    type(fgsl_permutation),  intent(in) :: pa, pb
    integer(fgsl_int) :: fgsl_permutation_mul
    fgsl_permutation_mul = gsl_permutation_mul(p%gsl_permutation, &
         pa%gsl_permutation, pb%gsl_permutation)
  end function fgsl_permutation_mul
  function fgsl_permutation_fwrite(stream, p)
    type(fgsl_file), intent(in) :: stream
    type(fgsl_permutation), intent(in) :: p
    integer(fgsl_int) :: fgsl_permutation_fwrite
    fgsl_permutation_fwrite = gsl_permutation_fwrite(stream%gsl_file, p%gsl_permutation)
  end function fgsl_permutation_fwrite
  function fgsl_permutation_fread(stream, p)
    type(fgsl_file), intent(in) :: stream
    type(fgsl_permutation), intent(inout) :: p
    integer(fgsl_int) :: fgsl_permutation_fread
    fgsl_permutation_fread = gsl_permutation_fread(stream%gsl_file, p%gsl_permutation)
  end function fgsl_permutation_fread
  function fgsl_permutation_fprintf(stream, p, format) 
    type(fgsl_file), intent(in) :: stream
    type(fgsl_permutation), intent(in) :: p
    character(kind=fgsl_char, len=*), intent(in) :: format
    integer(fgsl_int) :: fgsl_permutation_fprintf
!
    character(kind=fgsl_char,len=fgsl_strmax) :: lrf
    if (len(trim(format)) < fgsl_strmax) then
       lrf = trim(format) // c_null_char 
       fgsl_permutation_fprintf = &
            gsl_permutation_fprintf(stream%gsl_file, p%gsl_permutation, lrf)
    else
       fgsl_permutation_fprintf = fgsl_failure
    end if
  end function fgsl_permutation_fprintf
  function fgsl_permutation_fscanf(stream, p) 
    type(fgsl_file), intent(in) :: stream
    type(fgsl_permutation), intent(inout) :: p
    integer(fgsl_int) :: fgsl_permutation_fscanf
    fgsl_permutation_fscanf = gsl_permutation_fscanf(stream%gsl_file, p%gsl_permutation)
  end function fgsl_permutation_fscanf
  function fgsl_permutation_linear_to_canonical(q, p)
    type(fgsl_permutation), intent(inout) :: q
    type(fgsl_permutation),  intent(in) :: p
    integer(fgsl_int) :: fgsl_permutation_linear_to_canonical
    fgsl_permutation_linear_to_canonical = &
         gsl_permutation_linear_to_canonical(q%gsl_permutation, p%gsl_permutation)
  end function fgsl_permutation_linear_to_canonical
  function fgsl_permutation_canonical_to_linear(p, q)
    type(fgsl_permutation), intent(inout) :: p
    type(fgsl_permutation),  intent(in) :: q
    integer(fgsl_int) :: fgsl_permutation_canonical_to_linear
    fgsl_permutation_canonical_to_linear = &
         gsl_permutation_canonical_to_linear(p%gsl_permutation, q%gsl_permutation)
  end function fgsl_permutation_canonical_to_linear
  function fgsl_permutation_inversions(p)
    type(fgsl_permutation), intent(in) :: p
    integer(fgsl_size_t) :: fgsl_permutation_inversions
    fgsl_permutation_inversions = gsl_permutation_inversions(p%gsl_permutation)
  end function fgsl_permutation_inversions
  function fgsl_permutation_linear_cycles(p)
    type(fgsl_permutation), intent(in) :: p
    integer(fgsl_size_t) :: fgsl_permutation_linear_cycles
    fgsl_permutation_linear_cycles = gsl_permutation_linear_cycles(p%gsl_permutation)
  end function fgsl_permutation_linear_cycles
  function fgsl_permutation_canonical_cycles(p)
    type(fgsl_permutation), intent(in) :: p
    integer(fgsl_size_t) :: fgsl_permutation_canonical_cycles
    fgsl_permutation_canonical_cycles = gsl_permutation_canonical_cycles(p%gsl_permutation)
  end function fgsl_permutation_canonical_cycles
!  
  function fgsl_combination_alloc(n, k)
    integer(fgsl_size_t), intent(in) :: n, k
    type(fgsl_combination) :: fgsl_combination_alloc
    fgsl_combination_alloc%gsl_combination = gsl_combination_alloc(n, k)
  end function fgsl_combination_alloc
  function fgsl_combination_calloc(n, k)
    integer(fgsl_size_t), intent(in) :: n, k
    type(fgsl_combination) :: fgsl_combination_calloc
    fgsl_combination_calloc%gsl_combination = gsl_combination_calloc(n, k)
  end function fgsl_combination_calloc
  subroutine fgsl_combination_init_first(c) 
    type(fgsl_combination), intent(inout) :: c
    call gsl_combination_init_first(c%gsl_combination)
  end subroutine fgsl_combination_init_first
  subroutine fgsl_combination_init_last(c) 
    type(fgsl_combination), intent(inout) :: c
    call gsl_combination_init_last(c%gsl_combination)
  end subroutine fgsl_combination_init_last
  subroutine fgsl_combination_free(c) 
    type(fgsl_combination), intent(inout) :: c
    call gsl_combination_free(c%gsl_combination)
  end subroutine fgsl_combination_free
  function fgsl_combination_memcpy(dest, src)
    type(fgsl_combination), intent(inout) :: dest
    type(fgsl_combination),  intent(in) :: src
    integer(fgsl_int) :: fgsl_combination_memcpy
    fgsl_combination_memcpy = gsl_combination_memcpy(dest%gsl_combination, src%gsl_combination)
  end function fgsl_combination_memcpy
  function fgsl_combination_get(c, i) 
    type(fgsl_combination), intent(inout) :: c
    integer(fgsl_size_t), intent(in) :: i
    integer(fgsl_size_t) :: fgsl_combination_get
    fgsl_combination_get = gsl_combination_get(c%gsl_combination, i)
  end function fgsl_combination_get
  function fgsl_combination_n(c)
    type(fgsl_combination), intent(in) :: c
    integer(fgsl_size_t) :: fgsl_combination_n
    fgsl_combination_n = gsl_combination_n(c%gsl_combination)
  end function fgsl_combination_n
  function fgsl_combination_k(c)
    type(fgsl_combination), intent(in) :: c
    integer(fgsl_size_t) :: fgsl_combination_k
    fgsl_combination_k = gsl_combination_k(c%gsl_combination)
  end function fgsl_combination_k
  function fgsl_combination_data(c)
    type(fgsl_combination), intent(in) :: c
    integer(fgsl_size_t), pointer :: fgsl_combination_data(:)
!
    integer(fgsl_size_t) :: size
    type(c_ptr) :: cdata
    size = gsl_combination_k(c%gsl_combination)
    cdata = gsl_combination_data(c%gsl_combination)
    if (c_associated(cdata)) then
       call c_f_pointer(cdata,fgsl_combination_data,(/size/))
    else
       nullify(fgsl_combination_data)
    end if
  end function fgsl_combination_data
  function fgsl_combination_valid(c) 
    type(fgsl_combination), intent(in) :: c
    integer(fgsl_int) :: fgsl_combination_valid
    fgsl_combination_valid = gsl_combination_valid(c%gsl_combination)
  end function fgsl_combination_valid
  function fgsl_combination_next(c)
    type(fgsl_combination), intent(in) :: c
    integer(fgsl_int) :: fgsl_combination_next
    fgsl_combination_next = gsl_combination_next(c%gsl_combination)
  end function fgsl_combination_next
  function fgsl_combination_prev(c)
    type(fgsl_combination), intent(in) :: c
    integer(fgsl_int) :: fgsl_combination_prev
    fgsl_combination_prev = gsl_combination_prev(c%gsl_combination)
  end function fgsl_combination_prev
  function fgsl_combination_fwrite(stream, c)
    type(fgsl_file), intent(in) :: stream
    type(fgsl_combination), intent(in) :: c
    integer(fgsl_int) :: fgsl_combination_fwrite
    fgsl_combination_fwrite = gsl_combination_fwrite(stream%gsl_file, c%gsl_combination)
  end function fgsl_combination_fwrite
  function fgsl_combination_fread(stream, c)
    type(fgsl_file), intent(in) :: stream
    type(fgsl_combination), intent(inout) :: c
    integer(fgsl_int) :: fgsl_combination_fread
    fgsl_combination_fread = gsl_combination_fread(stream%gsl_file, c%gsl_combination)
  end function fgsl_combination_fread
  function fgsl_combination_fprintf(stream, c, format) 
    type(fgsl_file), intent(in) :: stream
    type(fgsl_combination), intent(in) :: c
    character(kind=fgsl_char, len=*), intent(in) :: format
    integer(fgsl_int) :: fgsl_combination_fprintf
!
    character(kind=fgsl_char,len=fgsl_strmax) :: lrf
    if (len(trim(format)) < fgsl_strmax) then
       lrf = trim(format) // c_null_char 
       fgsl_combination_fprintf = &
            gsl_combination_fprintf(stream%gsl_file, c%gsl_combination, lrf)
    else
       fgsl_combination_fprintf = fgsl_failure
    end if
  end function fgsl_combination_fprintf
  function fgsl_combination_fscanf(stream, c) 
    type(fgsl_file), intent(in) :: stream
    type(fgsl_combination), intent(inout) :: c
    integer(fgsl_int) :: fgsl_combination_fscanf
    fgsl_combination_fscanf = gsl_combination_fscanf(stream%gsl_file, c%gsl_combination)
  end function fgsl_combination_fscanf

!
  function fgsl_permutation_status(permutation)
    type(fgsl_permutation), intent(in) :: permutation
    logical :: fgsl_permutation_status
    fgsl_permutation_status = .true.
    if (.not. c_associated(permutation%gsl_permutation)) fgsl_permutation_status = .false.
  end function fgsl_permutation_status
  function fgsl_combination_status(combination)
    type(fgsl_combination), intent(in) :: combination
    logical :: fgsl_combination_status
    fgsl_combination_status = .true.
    if (.not. c_associated(combination%gsl_combination)) fgsl_combination_status = .false.
  end function fgsl_combination_status
  function fgsl_sizeof_permutation(p)
    type(fgsl_permutation), intent(in) :: p
    integer(fgsl_size_t) :: fgsl_sizeof_permutation
    fgsl_sizeof_permutation = gsl_aux_sizeof_permutation()
  end function fgsl_sizeof_permutation
  function fgsl_sizeof_combination(c)
    type(fgsl_combination), intent(in) :: c
    integer(fgsl_size_t) :: fgsl_sizeof_combination
    fgsl_sizeof_combination = gsl_aux_sizeof_combination()
  end function fgsl_sizeof_combination


!-*-f90-*-
!
!  API: Sorting
!
  subroutine fgsl_heapsort(array, count, size, compare)
    type(c_ptr) :: array
    integer(fgsl_size_t), intent(in) :: count, size
    interface
       function compare(x, y) bind(c)
         import
         type(c_ptr), value :: x, y
         integer(c_int) :: compare
       end function compare
    end interface
!
    type(c_funptr) :: fp
    fp = c_funloc(compare)
    call gsl_heapsort(array, count, size, fp)
  end subroutine fgsl_heapsort
  function fgsl_heapsort_index(p, array, count, size, compare)
    integer(fgsl_size_t), intent(in) :: count, size
    integer(fgsl_size_t), intent(out) :: p(count)
    type(c_ptr) :: array
    interface
       function compare(x, y) bind(c)
         import
         type(c_ptr), value :: x, y
         integer(c_int) :: compare
       end function compare
    end interface
    integer(fgsl_int) :: fgsl_heapsort_index
!
    type(c_funptr) :: fp
    fp = c_funloc(compare)
    fgsl_heapsort_index = gsl_heapsort_index(p, array, count, size, fp)
  end function fgsl_heapsort_index
  subroutine fgsl_sort_double(data, stride, n)
    real(fgsl_double), intent(inout) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    call gsl_sort(data, stride, n)
  end subroutine fgsl_sort_double
  subroutine fgsl_sort_double_index(p, data, stride, n)
    integer(fgsl_size_t), intent(out) :: p(:)
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    call gsl_sort_index(p, data, stride, n)
  end subroutine fgsl_sort_double_index
  function fgsl_sort_double_smallest(dest, k, src, stride, n)
    integer(fgsl_size_t), intent(in) :: k, stride, n
    real(fgsl_double), intent(out) :: dest(k)
    real(fgsl_double), intent(in) :: src(:)
    integer(fgsl_int) :: fgsl_sort_double_smallest
    fgsl_sort_double_smallest = gsl_sort_smallest(dest, k, src, stride, n)
  end function fgsl_sort_double_smallest
  function fgsl_sort_double_smallest_index(p, k, src, stride, n)
    integer(fgsl_size_t), intent(in) :: k, stride, n
    integer(fgsl_size_t), intent(out) :: p(k)
    real(fgsl_double), intent(in) :: src(:)
    integer(fgsl_int) :: fgsl_sort_double_smallest_index
    fgsl_sort_double_smallest_index = gsl_sort_smallest_index(p, k, &
         src, stride, n)
  end function fgsl_sort_double_smallest_index
  function fgsl_sort_double_largest(dest, k, src, stride, n)
    integer(fgsl_size_t), intent(in) :: k, stride, n
    real(fgsl_double), intent(out) :: dest(k)
    real(fgsl_double), intent(in) :: src(:)
    integer(fgsl_int) :: fgsl_sort_double_largest
    fgsl_sort_double_largest = gsl_sort_largest(dest, k, src, stride, n)
  end function fgsl_sort_double_largest
  function fgsl_sort_double_largest_index(p, k, src, stride, n)
    integer(fgsl_size_t), intent(in) :: k, stride, n
    integer(fgsl_size_t), intent(out) :: p(k)
    real(fgsl_double), intent(in) :: src(:)
    integer(fgsl_int) :: fgsl_sort_double_largest_index
    fgsl_sort_double_largest_index = gsl_sort_largest_index(p, k, &
         src, stride, n)
  end function fgsl_sort_double_largest_index
  subroutine fgsl_sort_long(data, stride, n)
    integer(fgsl_long), intent(inout) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    call gsl_sort_long(data, stride, n)
  end subroutine fgsl_sort_long
  subroutine fgsl_sort_long_index(p, data, stride, n)
    integer(fgsl_size_t), intent(out) :: p(:)
    integer(fgsl_long), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    call gsl_sort_long_index(p, data, stride, n)
  end subroutine fgsl_sort_long_index
  function fgsl_sort_long_smallest(dest, k, src, stride, n)
    integer(fgsl_size_t), intent(in) :: k, stride, n
    integer(fgsl_long), intent(out) :: dest(k)
    integer(fgsl_long), intent(in) :: src(:)
    integer(fgsl_int) :: fgsl_sort_long_smallest
    fgsl_sort_long_smallest = gsl_sort_long_smallest(dest, k, src, stride, n)
  end function fgsl_sort_long_smallest
  function fgsl_sort_long_smallest_index(p, k, src, stride, n)
    integer(fgsl_size_t), intent(in) :: k, stride, n
    integer(fgsl_size_t), intent(out) :: p(k)
    integer(fgsl_long), intent(in) :: src(:)
    integer(fgsl_int) :: fgsl_sort_long_smallest_index
    fgsl_sort_long_smallest_index = gsl_sort_long_smallest_index(p, k, &
         src, stride, n)
  end function fgsl_sort_long_smallest_index
  function fgsl_sort_long_largest(dest, k, src, stride, n)
    integer(fgsl_size_t), intent(in) :: k, stride, n
    integer(fgsl_long), intent(out) :: dest(k)
    integer(fgsl_long), intent(in) :: src(:)
    integer(fgsl_int) :: fgsl_sort_long_largest
    fgsl_sort_long_largest = gsl_sort_long_largest(dest, k, src, stride, n)
  end function fgsl_sort_long_largest
  function fgsl_sort_long_largest_index(p, k, src, stride, n)
    integer(fgsl_size_t), intent(in) :: k, stride, n
    integer(fgsl_size_t), intent(out) :: p(k)
    integer(fgsl_long), intent(in) :: src(:)
    integer(fgsl_int) :: fgsl_sort_long_largest_index
    fgsl_sort_long_largest_index = gsl_sort_long_largest_index(p, k, &
         src, stride, n)
  end function fgsl_sort_long_largest_index
  
  
!-*-f90-*-
!
!  API: Linear algebra support
!
!  Matrices must be entered/read transposed  
!
! LU
!
  function fgsl_linalg_lu_decomp(a, p, signum) 
    type(fgsl_matrix)  :: a
    type(fgsl_permutation) :: p
    integer(fgsl_int) :: signum
    integer(fgsl_int) :: fgsl_linalg_lu_decomp
    fgsl_linalg_lu_decomp = gsl_linalg_lu_decomp(a%gsl_matrix, &
         p%gsl_permutation, signum)
  end function fgsl_linalg_lu_decomp
  function fgsl_linalg_complex_lu_decomp(a, p, signum) 
    type(fgsl_matrix_complex)  :: a
    type(fgsl_permutation) :: p
    integer(fgsl_int) :: signum
    integer(fgsl_int) :: fgsl_linalg_complex_lu_decomp
    fgsl_linalg_complex_lu_decomp = gsl_linalg_complex_lu_decomp( &
         a%gsl_matrix_complex, p%gsl_permutation, signum)
  end function fgsl_linalg_complex_lu_decomp
  function fgsl_linalg_lu_solve(lu, p, b, x) 
    type(fgsl_matrix), intent(in) :: lu
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector), intent(in) :: b
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_lu_solve
    fgsl_linalg_lu_solve = gsl_linalg_lu_solve(lu%gsl_matrix, &
         p%gsl_permutation, b%gsl_vector, x%gsl_vector)
  end function fgsl_linalg_lu_solve
  function fgsl_linalg_complex_lu_solve(lu, p, b, x) 
    type(fgsl_matrix_complex), intent(in) :: lu
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector_complex), intent(in) :: b
    type(fgsl_vector_complex), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_complex_lu_solve
    fgsl_linalg_complex_lu_solve = gsl_linalg_complex_lu_solve( &
         lu%gsl_matrix_complex, p%gsl_permutation, b%gsl_vector_complex, &
         x%gsl_vector_complex)
  end function fgsl_linalg_complex_lu_solve
  function fgsl_linalg_lu_svx(lu, p, x) 
    type(fgsl_matrix), intent(in) :: lu
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_lu_svx
    fgsl_linalg_lu_svx = gsl_linalg_lu_svx(lu%gsl_matrix, &
         p%gsl_permutation, x%gsl_vector)
  end function fgsl_linalg_lu_svx
  function fgsl_linalg_complex_lu_svx(lu, p, x) 
    type(fgsl_matrix_complex), intent(in) :: lu
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector_complex), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_complex_lu_svx
    fgsl_linalg_complex_lu_svx = gsl_linalg_complex_lu_svx(&
         lu%gsl_matrix_complex, p%gsl_permutation, x%gsl_vector_complex)
  end function fgsl_linalg_complex_lu_svx
  function fgsl_linalg_lu_refine (a, lu, p, b, x, residual)
   type(fgsl_matrix), intent(in) :: a, lu
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector), intent(in) :: b
    type(fgsl_vector), intent(inout) :: x
    type(fgsl_vector), intent(inout) :: residual
    integer(fgsl_int) :: fgsl_linalg_lu_refine
    fgsl_linalg_lu_refine = gsl_linalg_lu_refine (a%gsl_matrix, &
         lu%gsl_matrix, p%gsl_permutation, b%gsl_vector, x%gsl_vector, &
         residual%gsl_vector)
  end function fgsl_linalg_lu_refine
  function fgsl_linalg_complex_lu_refine (a, lu, p, b, x, residual)
   type(fgsl_matrix_complex), intent(in) :: a, lu
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector_complex), intent(in) :: b
    type(fgsl_vector_complex), intent(inout) :: x
    type(fgsl_vector_complex), intent(inout) :: residual
    integer(fgsl_int) :: fgsl_linalg_complex_lu_refine
    fgsl_linalg_complex_lu_refine = gsl_linalg_complex_lu_refine( &
         a%gsl_matrix_complex, lu%gsl_matrix_complex, p%gsl_permutation, &
         b%gsl_vector_complex, x%gsl_vector_complex, &
         residual%gsl_vector_complex)
  end function fgsl_linalg_complex_lu_refine
  function fgsl_linalg_lu_invert(lu, p, inverse) 
    type(fgsl_matrix), intent(in) :: lu
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_matrix), intent(inout) :: inverse
    integer(fgsl_int) :: fgsl_linalg_lu_invert
    fgsl_linalg_lu_invert = gsl_linalg_lu_invert(lu%gsl_matrix, &
         p%gsl_permutation, inverse%gsl_matrix)
  end function fgsl_linalg_lu_invert
  function fgsl_linalg_complex_lu_invert(lu, p, inverse) 
    type(fgsl_matrix_complex), intent(in) :: lu
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_matrix_complex), intent(inout) :: inverse
    integer(fgsl_int) :: fgsl_linalg_complex_lu_invert
    fgsl_linalg_complex_lu_invert = gsl_linalg_complex_lu_invert(&
         lu%gsl_matrix_complex, p%gsl_permutation, inverse%gsl_matrix_complex)
  end function fgsl_linalg_complex_lu_invert
  function fgsl_linalg_lu_det(lu, signum) 
    type(fgsl_matrix), intent(in) :: lu
    integer(fgsl_int), intent(in) :: signum
    real(fgsl_double) :: fgsl_linalg_lu_det
    fgsl_linalg_lu_det = gsl_linalg_lu_det(lu%gsl_matrix, signum)
  end function fgsl_linalg_lu_det
  function fgsl_linalg_complex_lu_det(lu, signum) 
    type(fgsl_matrix_complex), intent(in) :: lu
    integer(fgsl_int), intent(in) :: signum
    complex(fgsl_double_complex) :: fgsl_linalg_complex_lu_det
    fgsl_linalg_complex_lu_det = gsl_linalg_complex_lu_det(&
         lu%gsl_matrix_complex, signum)
  end function fgsl_linalg_complex_lu_det
  function fgsl_linalg_lu_lndet(lu) 
    type(fgsl_matrix), intent(in) :: lu
    real(fgsl_double) :: fgsl_linalg_lu_lndet
    fgsl_linalg_lu_lndet = gsl_linalg_lu_lndet(lu%gsl_matrix)
  end function fgsl_linalg_lu_lndet
  function fgsl_linalg_complex_lu_lndet(lu) 
    type(fgsl_matrix_complex), intent(in) :: lu
    real(fgsl_double) :: fgsl_linalg_complex_lu_lndet
    fgsl_linalg_complex_lu_lndet = gsl_linalg_complex_lu_lndet(&
         lu%gsl_matrix_complex)
  end function fgsl_linalg_complex_lu_lndet
  function fgsl_linalg_lu_sgndet(lu, signum) 
    type(fgsl_matrix), intent(in) :: lu
    integer(fgsl_int), intent(in) :: signum
    integer(fgsl_int) :: fgsl_linalg_lu_sgndet
    fgsl_linalg_lu_sgndet = gsl_linalg_lu_sgndet(lu%gsl_matrix, signum)
  end function fgsl_linalg_lu_sgndet
  function fgsl_linalg_complex_lu_sgndet(lu, signum) 
    type(fgsl_matrix_complex), intent(in) :: lu
    integer(fgsl_int), intent(in) :: signum
    complex(fgsl_double_complex) :: fgsl_linalg_complex_lu_sgndet
    fgsl_linalg_complex_lu_sgndet = gsl_linalg_complex_lu_sgndet(&
         lu%gsl_matrix_complex, signum)
  end function fgsl_linalg_complex_lu_sgndet
!
! QR
!
  function fgsl_linalg_qr_decomp (a, tau) 
    type(fgsl_matrix), intent(inout) :: a
    type(fgsl_vector), intent(inout) :: tau
    integer(fgsl_int) :: fgsl_linalg_qr_decomp
    fgsl_linalg_qr_decomp = gsl_linalg_qr_decomp (a%gsl_matrix, tau%gsl_vector) 
  end function fgsl_linalg_qr_decomp
  function fgsl_linalg_qr_solve (qr, tau, b, x) 
    type(fgsl_matrix), intent(in) :: qr 
    type(fgsl_vector), intent(in) :: tau, b
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_qr_solve
    fgsl_linalg_qr_solve = gsl_linalg_qr_solve(qr%gsl_matrix, tau%gsl_vector, &
         b%gsl_vector, x%gsl_vector) 
  end function fgsl_linalg_qr_solve
  function fgsl_linalg_qr_svx (qr, tau, x) 
    type(fgsl_matrix), intent(in) :: qr
    type(fgsl_vector), intent(in) :: tau
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_qr_svx
    fgsl_linalg_qr_svx = gsl_linalg_qr_svx (qr%gsl_matrix, tau%gsl_vector, &
         x%gsl_vector) 
  end function fgsl_linalg_qr_svx
  function fgsl_linalg_qr_lssolve (qr, tau, b, x, residual) 
    type(fgsl_matrix), intent(in) :: qr 
    type(fgsl_vector), intent(in) :: tau, b
    type(fgsl_vector), intent(inout) :: x, residual
    integer(fgsl_int) :: fgsl_linalg_qr_lssolve
    fgsl_linalg_qr_lssolve = gsl_linalg_qr_lssolve (qr%gsl_matrix, tau%gsl_vector, &
         b%gsl_vector, x%gsl_vector, residual%gsl_vector)
  end function fgsl_linalg_qr_lssolve
  function fgsl_linalg_qr_qtvec (qr, tau, v) 
    type(fgsl_matrix), intent(in) :: qr
    type(fgsl_vector), intent(in) :: tau
    type(fgsl_vector), intent(inout) :: v
    integer(fgsl_int) :: fgsl_linalg_qr_qtvec
    fgsl_linalg_qr_qtvec = gsl_linalg_qr_qtvec (qr%gsl_matrix, tau%gsl_vector, &
         v%gsl_vector) 
  end function fgsl_linalg_qr_qtvec
  function fgsl_linalg_qr_qvec (qr, tau, v) 
    type(fgsl_matrix), intent(in) :: qr
    type(fgsl_vector), intent(in) :: tau
    type(fgsl_vector), intent(inout) :: v
    integer(fgsl_int) :: fgsl_linalg_qr_qvec
    fgsl_linalg_qr_qvec = gsl_linalg_qr_qvec (qr%gsl_matrix, tau%gsl_vector, &
         v%gsl_vector) 
  end function fgsl_linalg_qr_qvec
  function fgsl_linalg_qr_qtmat (qr, tau, a) 
    type(fgsl_matrix), intent(in) :: qr
    type(fgsl_vector), intent(in) :: tau
    type(fgsl_matrix), intent(inout) :: a
    integer(fgsl_int) :: fgsl_linalg_qr_qtmat
    fgsl_linalg_qr_qtmat = gsl_linalg_qr_qtmat (qr%gsl_matrix, tau%gsl_vector, &
         a%gsl_matrix) 
  end function fgsl_linalg_qr_qtmat
  function fgsl_linalg_qr_rsolve (qr, b, x) 
    type(fgsl_matrix), intent(in) :: qr 
    type(fgsl_vector), intent(in) :: b
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_qr_rsolve
    fgsl_linalg_qr_rsolve = gsl_linalg_qr_rsolve(qr%gsl_matrix, &
         b%gsl_vector, x%gsl_vector) 
  end function fgsl_linalg_qr_rsolve
  function fgsl_linalg_qr_rsvx (qr, x) 
    type(fgsl_matrix), intent(in) :: qr 
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_qr_rsvx
    fgsl_linalg_qr_rsvx = gsl_linalg_qr_rsvx(qr%gsl_matrix, x%gsl_vector) 
  end function fgsl_linalg_qr_rsvx
  function fgsl_linalg_qr_unpack (qr, tau, q, r) 
    type(fgsl_matrix), intent(in) :: qr
    type(fgsl_vector), intent(in) :: tau
    type(fgsl_matrix), intent(inout) :: q, r
    integer(fgsl_int) :: fgsl_linalg_qr_unpack
    fgsl_linalg_qr_unpack = gsl_linalg_qr_unpack (qr%gsl_matrix, &
         tau%gsl_vector, q%gsl_matrix, r%gsl_matrix) 
  end function fgsl_linalg_qr_unpack
  function fgsl_linalg_qr_qrsolve (q, r, b, x) 
    type(fgsl_matrix), intent(in) :: q, r 
    type(fgsl_vector), intent(in) :: b
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_qr_qrsolve
    fgsl_linalg_qr_qrsolve = gsl_linalg_qr_qrsolve(q%gsl_matrix, &
         r%gsl_matrix, b%gsl_vector, x%gsl_vector) 
  end function fgsl_linalg_qr_qrsolve
  function fgsl_linalg_qr_update (q, r, w, v) 
    type(fgsl_matrix), intent(inout) :: q, r 
    type(fgsl_vector), intent(inout) :: w
    type(fgsl_vector), intent(in) :: v
    integer(fgsl_int) :: fgsl_linalg_qr_update
    fgsl_linalg_qr_update = gsl_linalg_qr_update(q%gsl_matrix, &
         r%gsl_matrix, w%gsl_vector, v%gsl_vector) 
  end function fgsl_linalg_qr_update
  function fgsl_linalg_r_solve (r, b, x) 
    type(fgsl_matrix), intent(in) :: r 
    type(fgsl_vector), intent(in) :: b
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_r_solve
    fgsl_linalg_r_solve = gsl_linalg_r_solve(r%gsl_matrix, &
         b%gsl_vector, x%gsl_vector) 
  end function fgsl_linalg_r_solve
  function fgsl_linalg_r_svx (r, x) 
    type(fgsl_matrix), intent(in) :: r 
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_r_svx
    fgsl_linalg_r_svx = gsl_linalg_r_svx(r%gsl_matrix, x%gsl_vector) 
  end function fgsl_linalg_r_svx
  function fgsl_linalg_qrpt_decomp (a, tau, p, signum, norm) 
    type(fgsl_matrix), intent(inout) :: a
    type(fgsl_vector), intent(inout) :: tau, norm
    type(fgsl_permutation), intent(inout) :: p
    integer(fgsl_int), intent(out) :: signum
    integer(fgsl_int) :: fgsl_linalg_qrpt_decomp
    fgsl_linalg_qrpt_decomp = gsl_linalg_qrpt_decomp (a%gsl_matrix, &
         tau%gsl_vector, p%gsl_permutation, signum, norm%gsl_vector) 
  end function fgsl_linalg_qrpt_decomp
  function fgsl_linalg_qrpt_decomp2 (a, q, r, tau, p, signum, norm) 
    type(fgsl_matrix), intent(in) :: a
    type(fgsl_matrix), intent(inout) :: q, r
    type(fgsl_vector), intent(inout) :: tau, norm
    type(fgsl_permutation), intent(inout) :: p
    integer(fgsl_int), intent(out) :: signum
    integer(fgsl_int) :: fgsl_linalg_qrpt_decomp2
    fgsl_linalg_qrpt_decomp2 = gsl_linalg_qrpt_decomp2 (a%gsl_matrix, &
         q%gsl_matrix, r%gsl_matrix, tau%gsl_vector, p%gsl_permutation, &
         signum, norm%gsl_vector)
  end function fgsl_linalg_qrpt_decomp2
  function fgsl_linalg_qrpt_solve (qr, tau, p, b, x) 
    type(fgsl_matrix), intent(in) :: qr 
    type(fgsl_vector), intent(in) :: tau, b
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_qrpt_solve
    fgsl_linalg_qrpt_solve = gsl_linalg_qrpt_solve(qr%gsl_matrix, &
         tau%gsl_vector, p%gsl_permutation, b%gsl_vector, x%gsl_vector) 
  end function fgsl_linalg_qrpt_solve
  function fgsl_linalg_qrpt_svx (qr, tau, p, x) 
    type(fgsl_matrix), intent(in) :: qr 
    type(fgsl_vector), intent(in) :: tau
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_qrpt_svx
    fgsl_linalg_qrpt_svx = gsl_linalg_qrpt_svx(qr%gsl_matrix, &
         tau%gsl_vector, p%gsl_permutation, x%gsl_vector) 
  end function fgsl_linalg_qrpt_svx
  function fgsl_linalg_qrpt_qrsolve (q, r, p, b, x) 
    type(fgsl_matrix), intent(in) :: q, r 
    type(fgsl_vector), intent(in) :: b
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_qrpt_qrsolve
    fgsl_linalg_qrpt_qrsolve = gsl_linalg_qrpt_qrsolve(q%gsl_matrix, &
         r%gsl_matrix, p%gsl_permutation, b%gsl_vector, x%gsl_vector) 
  end function fgsl_linalg_qrpt_qrsolve
  function fgsl_linalg_qrpt_update (q, r, p, w, v) 
    type(fgsl_matrix), intent(inout) :: q, r 
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector), intent(in) :: v
    type(fgsl_vector), intent(inout) :: w
    integer(fgsl_int) :: fgsl_linalg_qrpt_update
    fgsl_linalg_qrpt_update = gsl_linalg_qrpt_update(q%gsl_matrix, &
         r%gsl_matrix, p%gsl_permutation, w%gsl_vector, v%gsl_vector) 
  end function fgsl_linalg_qrpt_update
  function fgsl_linalg_qrpt_rsolve (qr, p, b, x) 
    type(fgsl_matrix), intent(in) :: qr 
    type(fgsl_vector), intent(in) :: b
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_qrpt_rsolve
    fgsl_linalg_qrpt_rsolve = gsl_linalg_qrpt_rsolve(qr%gsl_matrix, &
         p%gsl_permutation, b%gsl_vector, x%gsl_vector) 
  end function fgsl_linalg_qrpt_rsolve
  function fgsl_linalg_qrpt_rsvx (qr, p, x) 
    type(fgsl_matrix), intent(in) :: qr 
    type(fgsl_permutation), intent(in) :: p
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_qrpt_rsvx
    fgsl_linalg_qrpt_rsvx = gsl_linalg_qrpt_rsvx(qr%gsl_matrix, &
         p%gsl_permutation, x%gsl_vector) 
  end function fgsl_linalg_qrpt_rsvx
!
! SVD
!
  function fgsl_linalg_sv_decomp(a, v, s, work) 
    type(fgsl_matrix), intent(inout)  :: a, v
    type(fgsl_vector), intent(inout) :: s, work
    integer(fgsl_int) :: fgsl_linalg_sv_decomp
    fgsl_linalg_sv_decomp = gsl_linalg_sv_decomp(a%gsl_matrix, &
         v%gsl_matrix, s%gsl_vector, work%gsl_vector)
  end function fgsl_linalg_sv_decomp
  function fgsl_linalg_sv_decomp_mod(a, x, v, s, work) 
    type(fgsl_matrix), intent(inout)  :: a, x, v
    type(fgsl_vector), intent(inout) :: s, work
    integer(fgsl_int) :: fgsl_linalg_sv_decomp_mod
    fgsl_linalg_sv_decomp_mod = gsl_linalg_sv_decomp_mod(a%gsl_matrix, &
         x%gsl_matrix, v%gsl_matrix, s%gsl_vector, work%gsl_vector)
  end function fgsl_linalg_sv_decomp_mod
  function fgsl_linalg_sv_decomp_jacobi(a, v, s) 
    type(fgsl_matrix), intent(inout)  :: a, v
    type(fgsl_vector), intent(inout) :: s
    integer(fgsl_int) :: fgsl_linalg_sv_decomp_jacobi
    fgsl_linalg_sv_decomp_jacobi = gsl_linalg_sv_decomp_jacobi(a%gsl_matrix, &
         v%gsl_matrix, s%gsl_vector)
  end function fgsl_linalg_sv_decomp_jacobi
  function fgsl_linalg_sv_solve(u, v, s, b, x) 
    type(fgsl_matrix), intent(in)  :: u, v
    type(fgsl_vector), intent(in) :: s, b
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_sv_solve
    fgsl_linalg_sv_solve = gsl_linalg_sv_solve(u%gsl_matrix, &
         v%gsl_matrix, s%gsl_vector, b%gsl_vector, x%gsl_vector)
  end function fgsl_linalg_sv_solve
!
! Cholesky
!
  function fgsl_linalg_cholesky_decomp(a) 
    type(fgsl_matrix), intent(inout)  :: a
    integer(fgsl_int) :: fgsl_linalg_cholesky_decomp
    fgsl_linalg_cholesky_decomp = gsl_linalg_cholesky_decomp(a%gsl_matrix)
  end function fgsl_linalg_cholesky_decomp
  function fgsl_linalg_complex_cholesky_decomp(a) 
    type(fgsl_matrix_complex), intent(inout)  :: a
    integer(fgsl_int) :: fgsl_linalg_complex_cholesky_decomp
    fgsl_linalg_complex_cholesky_decomp = &
         gsl_linalg_complex_cholesky_decomp(a%gsl_matrix_complex)
  end function fgsl_linalg_complex_cholesky_decomp
  function fgsl_linalg_cholesky_solve(chol, b, x) 
    type(fgsl_matrix), intent(in) :: chol
    type(fgsl_vector), intent(in) :: b
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_cholesky_solve
    fgsl_linalg_cholesky_solve = gsl_linalg_cholesky_solve(chol%gsl_matrix, &
         b%gsl_vector, x%gsl_vector)
  end function fgsl_linalg_cholesky_solve
  function fgsl_linalg_complex_cholesky_solve(chol, b, x) 
    type(fgsl_matrix_complex), intent(in) :: chol
    type(fgsl_vector_complex), intent(in) :: b
    type(fgsl_vector_complex), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_complex_cholesky_solve
    fgsl_linalg_complex_cholesky_solve = &
         gsl_linalg_complex_cholesky_solve(chol%gsl_matrix_complex, &
         b%gsl_vector_complex, x%gsl_vector_complex)
  end function fgsl_linalg_complex_cholesky_solve
  function fgsl_linalg_cholesky_svx(chol, x) 
    type(fgsl_matrix), intent(in) :: chol
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_cholesky_svx
    fgsl_linalg_cholesky_svx = gsl_linalg_cholesky_svx(chol%gsl_matrix, &
         x%gsl_vector)
  end function fgsl_linalg_cholesky_svx
  function fgsl_linalg_complex_cholesky_svx(chol, x) 
    type(fgsl_matrix_complex), intent(in) :: chol
    type(fgsl_vector_complex), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_complex_cholesky_svx
    fgsl_linalg_complex_cholesky_svx = &
         gsl_linalg_complex_cholesky_svx(chol%gsl_matrix_complex, &
         x%gsl_vector_complex)
  end function fgsl_linalg_complex_cholesky_svx
  function fgsl_linalg_cholesky_invert (chol)
    integer(fgsl_int) :: fgsl_linalg_cholesky_invert
    type(fgsl_matrix), intent(inout) :: chol
    fgsl_linalg_cholesky_invert = gsl_linalg_cholesky_invert(chol%gsl_matrix)
  end function fgsl_linalg_cholesky_invert
!
! Tridiag
!
  function fgsl_linalg_symmtd_decomp(a, tau) 
    type(fgsl_matrix), intent(inout)  :: a
    type(fgsl_vector), intent(inout) :: tau
    integer(fgsl_int) :: fgsl_linalg_symmtd_decomp
    fgsl_linalg_symmtd_decomp = gsl_linalg_symmtd_decomp(a%gsl_matrix, &
         tau%gsl_vector)
  end function fgsl_linalg_symmtd_decomp
  function fgsl_linalg_symmtd_unpack(a, tau, q, diag, subdiag) 
    type(fgsl_matrix), intent(in)  :: a
    type(fgsl_vector), intent(in) :: tau
    type(fgsl_matrix), intent(inout)  :: q
    type(fgsl_vector), intent(inout) :: diag, subdiag
    integer(fgsl_int) :: fgsl_linalg_symmtd_unpack
    fgsl_linalg_symmtd_unpack = gsl_linalg_symmtd_unpack(a%gsl_matrix, &
         tau%gsl_vector, q%gsl_matrix, diag%gsl_vector, subdiag%gsl_vector)
  end function fgsl_linalg_symmtd_unpack
  function fgsl_linalg_symmtd_unpack_t(a, diag, subdiag) 
    type(fgsl_matrix), intent(in)  :: a
    type(fgsl_vector), intent(inout) :: diag, subdiag
    integer(fgsl_int) :: fgsl_linalg_symmtd_unpack_t
    fgsl_linalg_symmtd_unpack_t = gsl_linalg_symmtd_unpack_t(a%gsl_matrix, &
         diag%gsl_vector, subdiag%gsl_vector)
  end function fgsl_linalg_symmtd_unpack_t
  function fgsl_linalg_hermtd_decomp(a, tau) 
    type(fgsl_matrix_complex), intent(inout)  :: a
    type(fgsl_vector_complex), intent(inout) :: tau
    integer(fgsl_int) :: fgsl_linalg_hermtd_decomp
    fgsl_linalg_hermtd_decomp = &
         gsl_linalg_hermtd_decomp(a%gsl_matrix_complex, tau%gsl_vector_complex)
  end function fgsl_linalg_hermtd_decomp
  function fgsl_linalg_hermtd_unpack(a, tau, q, diag, subdiag) 
    type(fgsl_matrix_complex), intent(in)  :: a
    type(fgsl_vector_complex), intent(in) :: tau
    type(fgsl_matrix_complex), intent(inout)  :: q
    type(fgsl_vector), intent(inout) :: diag, subdiag
    integer(fgsl_int) :: fgsl_linalg_hermtd_unpack
    fgsl_linalg_hermtd_unpack = gsl_linalg_hermtd_unpack( &
         a%gsl_matrix_complex, tau%gsl_vector_complex, &
         q%gsl_matrix_complex, diag%gsl_vector, subdiag%gsl_vector)
  end function fgsl_linalg_hermtd_unpack
  function fgsl_linalg_hermtd_unpack_t(a, diag, subdiag) 
    type(fgsl_matrix_complex), intent(in)  :: a
    type(fgsl_vector), intent(inout) :: diag, subdiag
    integer(fgsl_int) :: fgsl_linalg_hermtd_unpack_t
    fgsl_linalg_hermtd_unpack_t = gsl_linalg_hermtd_unpack_t( &
         a%gsl_matrix_complex, diag%gsl_vector, subdiag%gsl_vector)
  end function fgsl_linalg_hermtd_unpack_t
!
! Hessenberg
!
  function fgsl_linalg_hessenberg_decomp(a, tau) 
    type(fgsl_matrix), intent(inout)  :: a
    type(fgsl_vector), intent(inout) :: tau
    integer(fgsl_int) :: fgsl_linalg_hessenberg_decomp
    fgsl_linalg_hessenberg_decomp = gsl_linalg_hessenberg_decomp(a%gsl_matrix, &
         tau%gsl_vector)
  end function fgsl_linalg_hessenberg_decomp
  function fgsl_linalg_hessenberg_unpack(h, tau, u) 
    type(fgsl_matrix), intent(in)  :: h
    type(fgsl_vector), intent(in) :: tau
    type(fgsl_matrix), intent(inout)  :: u
    integer(fgsl_int) :: fgsl_linalg_hessenberg_unpack
    fgsl_linalg_hessenberg_unpack = gsl_linalg_hessenberg_unpack(h%gsl_matrix, &
         tau%gsl_vector, u%gsl_matrix)
  end function fgsl_linalg_hessenberg_unpack
  function fgsl_linalg_hessenberg_unpack_accum(h, tau, v) 
    type(fgsl_matrix), intent(in)  :: h
    type(fgsl_vector), intent(in) :: tau
    type(fgsl_matrix), intent(inout)  :: v
    integer(fgsl_int) :: fgsl_linalg_hessenberg_unpack_accum
    fgsl_linalg_hessenberg_unpack_accum = &
         gsl_linalg_hessenberg_unpack_accum(h%gsl_matrix, &
         tau%gsl_vector, v%gsl_matrix)
  end function fgsl_linalg_hessenberg_unpack_accum
  function fgsl_linalg_hessenberg_set_zero(h) 
    type(fgsl_matrix), intent(inout)  :: h
    integer(fgsl_int) :: fgsl_linalg_hessenberg_set_zero
    fgsl_linalg_hessenberg_set_zero = &
         gsl_linalg_hessenberg_set_zero(h%gsl_matrix)
  end function fgsl_linalg_hessenberg_set_zero
  function fgsl_linalg_hesstri_decomp(a, b, u, v, work) 
    type(fgsl_matrix), intent(inout)  :: a, b, u, v
    type(fgsl_vector), intent(inout) :: work
    integer(fgsl_int) :: fgsl_linalg_hesstri_decomp
    fgsl_linalg_hesstri_decomp = gsl_linalg_hesstri_decomp(a%gsl_matrix, &
         b%gsl_matrix, u%gsl_matrix, v%gsl_matrix, work%gsl_vector)
  end function fgsl_linalg_hesstri_decomp
!
! Bidiag
!
  function fgsl_linalg_bidiag_decomp(a, tau_u, tau_v) 
    type(fgsl_matrix), intent(inout)  :: a
    type(fgsl_vector), intent(inout) :: tau_u, tau_v
    integer(fgsl_int) :: fgsl_linalg_bidiag_decomp
    fgsl_linalg_bidiag_decomp = gsl_linalg_bidiag_decomp(a%gsl_matrix, &
         tau_u%gsl_vector, tau_v%gsl_vector)
  end function fgsl_linalg_bidiag_decomp
  function fgsl_linalg_bidiag_unpack(a, tau_u, u, tau_v, v, diag, &
       superdiag) 
    type(fgsl_matrix), intent(in)  :: a
    type(fgsl_vector), intent(in) :: tau_u, tau_v
    type(fgsl_matrix), intent(inout) :: u, v
    type(fgsl_vector), intent(inout) :: diag, superdiag
    integer(fgsl_int) :: fgsl_linalg_bidiag_unpack
    fgsl_linalg_bidiag_unpack = gsl_linalg_bidiag_unpack(a%gsl_matrix, &
         tau_u%gsl_vector, u%gsl_matrix, tau_v%gsl_vector, v%gsl_matrix, &
         diag%gsl_vector, superdiag%gsl_vector)
  end function fgsl_linalg_bidiag_unpack
  function fgsl_linalg_bidiag_unpack2(a, tau_u, tau_v, v) 
    type(fgsl_matrix), intent(inout)  :: a
    type(fgsl_vector), intent(in) :: tau_u, tau_v
    type(fgsl_matrix), intent(inout) :: v
    integer(fgsl_int) :: fgsl_linalg_bidiag_unpack2
    fgsl_linalg_bidiag_unpack2 = gsl_linalg_bidiag_unpack2(a%gsl_matrix, &
         tau_u%gsl_vector, tau_v%gsl_vector, v%gsl_matrix)
  end function fgsl_linalg_bidiag_unpack2
  function fgsl_linalg_bidiag_unpack_b(a, diag, superdiag) 
    type(fgsl_matrix), intent(in)  :: a
    type(fgsl_vector), intent(inout) :: diag, superdiag
    integer(fgsl_int) :: fgsl_linalg_bidiag_unpack_b
    fgsl_linalg_bidiag_unpack_b = gsl_linalg_bidiag_unpack_b(a%gsl_matrix, &
         diag%gsl_vector, superdiag%gsl_vector)
  end function fgsl_linalg_bidiag_unpack_b
!
! Householder
!
  function fgsl_linalg_householder_transform (v) 
    type(fgsl_vector), intent(inout) :: v
    real(fgsl_double) :: fgsl_linalg_householder_transform
    fgsl_linalg_householder_transform = &
         gsl_linalg_householder_transform(v%gsl_vector)
  end function fgsl_linalg_householder_transform
  function fgsl_linalg_complex_householder_transform (v) 
    type(fgsl_vector), intent(inout) :: v
    complex(fgsl_double_complex) :: fgsl_linalg_complex_householder_transform
    fgsl_linalg_complex_householder_transform = &
         gsl_linalg_complex_householder_transform(v%gsl_vector)
  end function fgsl_linalg_complex_householder_transform
  function fgsl_linalg_householder_hm (tau, v, a)
    real(fgsl_double), intent(in) :: tau
    type(fgsl_vector), intent(in) :: v
    type(fgsl_matrix), intent(inout) :: a
    integer(fgsl_int) :: fgsl_linalg_householder_hm
    fgsl_linalg_householder_hm = &
         gsl_linalg_householder_hm(tau, v%gsl_vector, a%gsl_matrix)
  end function fgsl_linalg_householder_hm
  function fgsl_linalg_complex_householder_hm (tau, v, a)
    complex(fgsl_double_complex), intent(in) :: tau
    type(fgsl_vector_complex), intent(in) :: v
    type(fgsl_matrix_complex), intent(inout) :: a
    integer(fgsl_int) :: fgsl_linalg_complex_householder_hm
!
    type(gsl_complex) :: targ
    targ = tau
    fgsl_linalg_complex_householder_hm = &
         gsl_linalg_complex_householder_hm(targ, v%gsl_vector_complex, &
         a%gsl_matrix_complex)
  end function fgsl_linalg_complex_householder_hm
  function fgsl_linalg_householder_mh (tau, v, a)
    real(fgsl_double), intent(in) :: tau
    type(fgsl_vector), intent(in) :: v
    type(fgsl_matrix), intent(inout) :: a
    integer(fgsl_int) :: fgsl_linalg_householder_mh
    fgsl_linalg_householder_mh = &
         gsl_linalg_householder_mh(tau, v%gsl_vector, a%gsl_matrix)
  end function fgsl_linalg_householder_mh
  function fgsl_linalg_complex_householder_mh (tau, v, a)
    complex(fgsl_double_complex), intent(in) :: tau
    type(fgsl_vector_complex), intent(in) :: v
    type(fgsl_matrix_complex), intent(inout) :: a
    integer(fgsl_int) :: fgsl_linalg_complex_householder_mh
!
    type(gsl_complex) :: targ
    targ = tau
    fgsl_linalg_complex_householder_mh = &
         gsl_linalg_complex_householder_mh(targ, v%gsl_vector_complex, &
         a%gsl_matrix_complex)
  end function fgsl_linalg_complex_householder_mh
  function fgsl_linalg_householder_hv (tau, v, w)
    real(fgsl_double), intent(in) :: tau
    type(fgsl_vector), intent(in) :: v
    type(fgsl_vector), intent(inout) :: w
    integer(fgsl_int) :: fgsl_linalg_householder_hv
    fgsl_linalg_householder_hv = &
         gsl_linalg_householder_hv(tau, v%gsl_vector, w%gsl_vector)
  end function fgsl_linalg_householder_hv
  function fgsl_linalg_complex_householder_hv (tau, v, w)
    complex(fgsl_double_complex), intent(in) :: tau
    type(fgsl_vector_complex), intent(in) :: v
    type(fgsl_vector_complex), intent(inout) :: w
    integer(fgsl_int) :: fgsl_linalg_complex_householder_hv
!
    type(gsl_complex) :: targ
    targ = tau
    fgsl_linalg_complex_householder_hv = &
         gsl_linalg_complex_householder_hv(targ, v%gsl_vector_complex, &
         w%gsl_vector_complex)
  end function fgsl_linalg_complex_householder_hv
  function fgsl_linalg_hh_solve(a, b, x) 
    type(fgsl_matrix), intent(inout) :: a
    type(fgsl_vector), intent(in) :: b
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_hh_solve
    fgsl_linalg_hh_solve = gsl_linalg_hh_solve(a%gsl_matrix, &
         b%gsl_vector, x%gsl_vector)
  end function fgsl_linalg_hh_solve
  function fgsl_linalg_hh_svx(a, x) 
    type(fgsl_matrix), intent(inout) :: a
    type(fgsl_vector), intent(inout) :: x
    integer(fgsl_int) :: fgsl_linalg_hh_svx
    fgsl_linalg_hh_svx = gsl_linalg_hh_svx(a%gsl_matrix, &
         x%gsl_vector)
  end function fgsl_linalg_hh_svx
!
! Tridiagonal
!
  function fgsl_linalg_solve_tridiag(diag, e, f, b, x) 
    type(fgsl_vector), intent(in) :: diag, e, f, b
    type(fgsl_vector), intent(inout) :: x
    integer(c_int) :: fgsl_linalg_solve_tridiag
    fgsl_linalg_solve_tridiag = gsl_linalg_solve_tridiag(diag%gsl_vector, &
         e%gsl_vector, f%gsl_vector, b%gsl_vector, x%gsl_vector)
  end function fgsl_linalg_solve_tridiag
  function fgsl_linalg_solve_symm_tridiag(diag, e, b, x) 
    type(fgsl_vector), intent(in) :: diag, e, b
    type(fgsl_vector), intent(inout) :: x
    integer(c_int) :: fgsl_linalg_solve_symm_tridiag
    fgsl_linalg_solve_symm_tridiag = &
         gsl_linalg_solve_symm_tridiag(diag%gsl_vector, &
         e%gsl_vector, b%gsl_vector, x%gsl_vector)
  end function fgsl_linalg_solve_symm_tridiag
  function fgsl_linalg_solve_cyc_tridiag(diag, e, f, b, x) 
    type(fgsl_vector), intent(in) :: diag, e, f, b
    type(fgsl_vector), intent(inout) :: x
    integer(c_int) :: fgsl_linalg_solve_cyc_tridiag
    fgsl_linalg_solve_cyc_tridiag = gsl_linalg_solve_cyc_tridiag(diag%gsl_vector, &
         e%gsl_vector, f%gsl_vector, b%gsl_vector, x%gsl_vector)
  end function fgsl_linalg_solve_cyc_tridiag
  function fgsl_linalg_solve_symm_cyc_tridiag(diag, e, b, x) 
    type(fgsl_vector), intent(in) :: diag, e, b
    type(fgsl_vector), intent(inout) :: x
    integer(c_int) :: fgsl_linalg_solve_symm_cyc_tridiag
    fgsl_linalg_solve_symm_cyc_tridiag = &
         gsl_linalg_solve_symm_cyc_tridiag(diag%gsl_vector, &
         e%gsl_vector, b%gsl_vector, x%gsl_vector)
  end function fgsl_linalg_solve_symm_cyc_tridiag
  function fgsl_linalg_balance_matrix(a, d) 
    type(fgsl_matrix), intent(inout) :: a
    type(fgsl_vector), intent(inout) :: d
    integer(fgsl_int) :: fgsl_linalg_balance_matrix
    fgsl_linalg_balance_matrix = gsl_linalg_balance_matrix(a%gsl_matrix, &
         d%gsl_vector)
  end function fgsl_linalg_balance_matrix
!-*-f90-*-
!
!  Interfaces: Eigensystem support
!
function fgsl_eigen_symm_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_symm_workspace) :: fgsl_eigen_symm_alloc
  fgsl_eigen_symm_alloc%gsl_eigen_symm_workspace = gsl_eigen_symm_alloc(n)
end function fgsl_eigen_symm_alloc
subroutine fgsl_eigen_symm_free(w) 
  type(fgsl_eigen_symm_workspace) :: w
  call gsl_eigen_symm_free(w%gsl_eigen_symm_workspace)
end subroutine fgsl_eigen_symm_free
function fgsl_eigen_symm(a, eval, w) 
  type(fgsl_matrix), intent(inout) :: a
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_eigen_symm_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_symm
  fgsl_eigen_symm = gsl_eigen_symm(a%gsl_matrix, eval%gsl_vector, &
       w%gsl_eigen_symm_workspace)
end function fgsl_eigen_symm
function fgsl_eigen_symmv_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_symmv_workspace) :: fgsl_eigen_symmv_alloc
  fgsl_eigen_symmv_alloc%gsl_eigen_symmv_workspace = gsl_eigen_symmv_alloc(n)
end function fgsl_eigen_symmv_alloc
subroutine fgsl_eigen_symmv_free(w) 
  type(fgsl_eigen_symmv_workspace) :: w
  call gsl_eigen_symmv_free(w%gsl_eigen_symmv_workspace)
end subroutine fgsl_eigen_symmv_free
function fgsl_eigen_symmv(a, eval, evec, w) 
  type(fgsl_matrix), intent(inout) :: a, evec
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_eigen_symmv_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_symmv
  fgsl_eigen_symmv = gsl_eigen_symmv(a%gsl_matrix, eval%gsl_vector, &
       evec%gsl_matrix, w%gsl_eigen_symmv_workspace)
end function fgsl_eigen_symmv
function fgsl_eigen_herm_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_herm_workspace) :: fgsl_eigen_herm_alloc
  fgsl_eigen_herm_alloc%gsl_eigen_herm_workspace = gsl_eigen_herm_alloc(n)
end function fgsl_eigen_herm_alloc
subroutine fgsl_eigen_herm_free(w) 
  type(fgsl_eigen_herm_workspace) :: w
  call gsl_eigen_herm_free(w%gsl_eigen_herm_workspace)
end subroutine fgsl_eigen_herm_free
function fgsl_eigen_herm(a, eval, w) 
  type(fgsl_matrix_complex), intent(inout) :: a
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_eigen_herm_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_herm
  fgsl_eigen_herm = gsl_eigen_herm(a%gsl_matrix_complex, eval%gsl_vector, &
       w%gsl_eigen_herm_workspace)
end function fgsl_eigen_herm
function fgsl_eigen_hermv_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_hermv_workspace) :: fgsl_eigen_hermv_alloc
  fgsl_eigen_hermv_alloc%gsl_eigen_hermv_workspace = gsl_eigen_hermv_alloc(n)
end function fgsl_eigen_hermv_alloc
subroutine fgsl_eigen_hermv_free(w) 
  type(fgsl_eigen_hermv_workspace) :: w
  call gsl_eigen_hermv_free(w%gsl_eigen_hermv_workspace)
end subroutine fgsl_eigen_hermv_free
function fgsl_eigen_hermv(a, eval, evec, w) 
  type(fgsl_matrix_complex), intent(inout) :: a, evec
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_eigen_hermv_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_hermv
  fgsl_eigen_hermv = gsl_eigen_hermv(a%gsl_matrix_complex, eval%gsl_vector, &
       evec%gsl_matrix_complex, w%gsl_eigen_hermv_workspace)
end function fgsl_eigen_hermv
function fgsl_eigen_nonsymm_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_nonsymm_workspace) :: fgsl_eigen_nonsymm_alloc
  fgsl_eigen_nonsymm_alloc%gsl_eigen_nonsymm_workspace = &
       gsl_eigen_nonsymm_alloc(n)
end function fgsl_eigen_nonsymm_alloc
subroutine fgsl_eigen_nonsymm_free(w) 
  type(fgsl_eigen_nonsymm_workspace) :: w
  call gsl_eigen_nonsymm_free(w%gsl_eigen_nonsymm_workspace)
end subroutine fgsl_eigen_nonsymm_free
subroutine fgsl_eigen_nonsymm_params (compute_t, balance, w) 
  integer(fgsl_int), intent(in) :: compute_t, balance
  type(fgsl_eigen_nonsymm_workspace), intent(inout) :: w
  call gsl_eigen_nonsymm_params (compute_t, balance, &
       w%gsl_eigen_nonsymm_workspace)
end subroutine fgsl_eigen_nonsymm_params
function fgsl_eigen_nonsymm(a, eval, w) 
  type(fgsl_matrix), intent(inout) :: a
  type(fgsl_vector_complex), intent(inout) :: eval
  type(fgsl_eigen_nonsymm_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_nonsymm
  fgsl_eigen_nonsymm = gsl_eigen_nonsymm(a%gsl_matrix, &
       eval%gsl_vector_complex, w%gsl_eigen_nonsymm_workspace)
end function fgsl_eigen_nonsymm
function fgsl_eigen_nonsymm_z(a, eval, z, w) 
  type(fgsl_matrix), intent(inout) :: a, z
  type(fgsl_vector_complex), intent(inout) :: eval
  type(fgsl_eigen_nonsymm_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_nonsymm_z
  fgsl_eigen_nonsymm_z = gsl_eigen_nonsymm_z(a%gsl_matrix, &
       eval%gsl_vector_complex, z%gsl_matrix, w%gsl_eigen_nonsymm_workspace)
end function fgsl_eigen_nonsymm_z
function fgsl_eigen_nonsymmv_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_nonsymmv_workspace) :: fgsl_eigen_nonsymmv_alloc
  fgsl_eigen_nonsymmv_alloc%gsl_eigen_nonsymmv_workspace = &
       gsl_eigen_nonsymmv_alloc(n)
end function fgsl_eigen_nonsymmv_alloc
subroutine fgsl_eigen_nonsymmv_free(w) 
  type(fgsl_eigen_nonsymmv_workspace) :: w
  call gsl_eigen_nonsymmv_free(w%gsl_eigen_nonsymmv_workspace)
end subroutine fgsl_eigen_nonsymmv_free
function fgsl_eigen_nonsymmv(a, eval, evec, w) 
  type(fgsl_matrix), intent(inout) :: a
  type(fgsl_vector_complex), intent(inout) :: eval
  type(fgsl_matrix_complex), intent(inout) :: evec
  type(fgsl_eigen_nonsymmv_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_nonsymmv
  fgsl_eigen_nonsymmv = gsl_eigen_nonsymmv(a%gsl_matrix, &
       eval%gsl_vector_complex, evec%gsl_matrix_complex, &
       w%gsl_eigen_nonsymmv_workspace)
end function fgsl_eigen_nonsymmv
function fgsl_eigen_nonsymmv_z(a, eval, evec, z, w) 
  type(fgsl_matrix), intent(inout) :: a, z
  type(fgsl_vector_complex), intent(inout) :: eval
  type(fgsl_matrix_complex), intent(inout) :: evec
  type(fgsl_eigen_nonsymmv_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_nonsymmv_z
  fgsl_eigen_nonsymmv_z = gsl_eigen_nonsymmv_z(a%gsl_matrix, &
       eval%gsl_vector_complex, evec%gsl_matrix_complex, &
       z%gsl_matrix, w%gsl_eigen_nonsymmv_workspace)
end function fgsl_eigen_nonsymmv_z
function fgsl_eigen_gensymm_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_gensymm_workspace) :: fgsl_eigen_gensymm_alloc
  fgsl_eigen_gensymm_alloc%gsl_eigen_gensymm_workspace = &
       gsl_eigen_gensymm_alloc(n)
end function fgsl_eigen_gensymm_alloc
subroutine fgsl_eigen_gensymm_free(w) 
  type(fgsl_eigen_gensymm_workspace) :: w
  call gsl_eigen_gensymm_free(w%gsl_eigen_gensymm_workspace)
end subroutine fgsl_eigen_gensymm_free
function fgsl_eigen_gensymm(a, b, eval, w) 
  type(fgsl_matrix), intent(inout) :: a, b
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_eigen_gensymm_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_gensymm
  fgsl_eigen_gensymm = gsl_eigen_gensymm(a%gsl_matrix, b%gsl_matrix, &
       eval%gsl_vector, w%gsl_eigen_gensymm_workspace)
end function fgsl_eigen_gensymm
function fgsl_eigen_gensymmv_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_gensymmv_workspace) :: fgsl_eigen_gensymmv_alloc
  fgsl_eigen_gensymmv_alloc%gsl_eigen_gensymmv_workspace = &
       gsl_eigen_gensymmv_alloc(n)
end function fgsl_eigen_gensymmv_alloc
subroutine fgsl_eigen_gensymmv_free(w) 
  type(fgsl_eigen_gensymmv_workspace) :: w
  call gsl_eigen_gensymmv_free(w%gsl_eigen_gensymmv_workspace)
end subroutine fgsl_eigen_gensymmv_free
function fgsl_eigen_gensymmv(a, b, eval, evec, w) 
  type(fgsl_matrix), intent(inout) :: a, b, evec
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_eigen_gensymmv_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_gensymmv
  fgsl_eigen_gensymmv = gsl_eigen_gensymmv(a%gsl_matrix, b%gsl_matrix, &
       eval%gsl_vector, evec%gsl_matrix, w%gsl_eigen_gensymmv_workspace)
end function fgsl_eigen_gensymmv
function fgsl_eigen_genherm_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_genherm_workspace) :: fgsl_eigen_genherm_alloc
  fgsl_eigen_genherm_alloc%gsl_eigen_genherm_workspace = &
       gsl_eigen_genherm_alloc(n)
end function fgsl_eigen_genherm_alloc
subroutine fgsl_eigen_genherm_free(w) 
  type(fgsl_eigen_genherm_workspace) :: w
  call gsl_eigen_genherm_free(w%gsl_eigen_genherm_workspace)
end subroutine fgsl_eigen_genherm_free
function fgsl_eigen_genherm(a, b, eval, w) 
  type(fgsl_matrix_complex), intent(inout) :: a, b
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_eigen_genherm_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_genherm
  fgsl_eigen_genherm = gsl_eigen_genherm(a%gsl_matrix_complex, b%gsl_matrix_complex, &
       eval%gsl_vector, w%gsl_eigen_genherm_workspace)
end function fgsl_eigen_genherm
function fgsl_eigen_genhermv_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_genhermv_workspace) :: fgsl_eigen_genhermv_alloc
  fgsl_eigen_genhermv_alloc%gsl_eigen_genhermv_workspace = &
       gsl_eigen_genhermv_alloc(n)
end function fgsl_eigen_genhermv_alloc
subroutine fgsl_eigen_genhermv_free(w) 
  type(fgsl_eigen_genhermv_workspace) :: w
  call gsl_eigen_genhermv_free(w%gsl_eigen_genhermv_workspace)
end subroutine fgsl_eigen_genhermv_free
function fgsl_eigen_genhermv(a, b, eval, evec, w) 
  type(fgsl_matrix_complex), intent(inout) :: a, b, evec
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_eigen_genhermv_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_genhermv
  fgsl_eigen_genhermv = gsl_eigen_genhermv(a%gsl_matrix_complex, &
       b%gsl_matrix_complex, eval%gsl_vector, evec%gsl_matrix_complex, &
       w%gsl_eigen_genhermv_workspace)
end function fgsl_eigen_genhermv
function fgsl_eigen_gen_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_gen_workspace) :: fgsl_eigen_gen_alloc
  fgsl_eigen_gen_alloc%gsl_eigen_gen_workspace = &
       gsl_eigen_gen_alloc(n)
end function fgsl_eigen_gen_alloc
subroutine fgsl_eigen_gen_free(w) 
  type(fgsl_eigen_gen_workspace) :: w
  call gsl_eigen_gen_free(w%gsl_eigen_gen_workspace)
end subroutine fgsl_eigen_gen_free
subroutine fgsl_eigen_gen_params (compute_s, compute_t, balance, w) 
  integer(fgsl_int), intent(in) :: compute_s, compute_t, balance
  type(fgsl_eigen_gen_workspace), intent(inout) :: w
  call gsl_eigen_gen_params (compute_s, compute_t, balance, &
       w%gsl_eigen_gen_workspace)
end subroutine fgsl_eigen_gen_params
function fgsl_eigen_gen(a, b, alpha, beta, w) 
  type(fgsl_matrix), intent(inout) :: a, b
  type(fgsl_vector_complex), intent(inout) :: alpha
  type(fgsl_vector), intent(inout) :: beta
  type(fgsl_eigen_gen_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_gen
  fgsl_eigen_gen = gsl_eigen_gen(a%gsl_matrix, &
       b%gsl_matrix, alpha%gsl_vector_complex, &
       beta%gsl_vector, w%gsl_eigen_gen_workspace)
end function fgsl_eigen_gen
function fgsl_eigen_gen_qz(a, b, alpha, beta, q, z, w) 
  type(fgsl_matrix), intent(inout) :: a, b, q, z
  type(fgsl_vector_complex), intent(inout) :: alpha
  type(fgsl_vector), intent(inout) :: beta
  type(fgsl_eigen_gen_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_gen_qz
  fgsl_eigen_gen_qz = gsl_eigen_gen_qz(a%gsl_matrix, &
       b%gsl_matrix, alpha%gsl_vector_complex, &
       beta%gsl_vector, q%gsl_matrix, z%gsl_matrix, &
       w%gsl_eigen_gen_workspace)
end function fgsl_eigen_gen_qz
function fgsl_eigen_genv_alloc(n) 
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_eigen_genv_workspace) :: fgsl_eigen_genv_alloc
  fgsl_eigen_genv_alloc%gsl_eigen_genv_workspace = &
       gsl_eigen_genv_alloc(n)
end function fgsl_eigen_genv_alloc
subroutine fgsl_eigen_genv_free(w) 
  type(fgsl_eigen_genv_workspace) :: w
  call gsl_eigen_genv_free(w%gsl_eigen_genv_workspace)
end subroutine fgsl_eigen_genv_free
function fgsl_eigen_genv(a, b, alpha, beta, evec, w) 
  type(fgsl_matrix), intent(inout) :: a, b
  type(fgsl_vector_complex), intent(inout) :: alpha
  type(fgsl_vector), intent(inout) :: beta
  type(fgsl_matrix_complex), intent(inout) :: evec
  type(fgsl_eigen_genv_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_genv
  fgsl_eigen_genv = gsl_eigen_genv(a%gsl_matrix, &
       b%gsl_matrix, alpha%gsl_vector_complex, &
       beta%gsl_vector, evec%gsl_matrix_complex, &
       w%gsl_eigen_genv_workspace)
end function fgsl_eigen_genv
function fgsl_eigen_genv_qz(a, b, alpha, beta, evec, q, z, w) 
  type(fgsl_matrix), intent(inout) :: a, b, q, z
  type(fgsl_vector_complex), intent(inout) :: alpha
  type(fgsl_vector), intent(inout) :: beta
  type(fgsl_matrix_complex), intent(inout) :: evec
  type(fgsl_eigen_genv_workspace) :: w
  integer(fgsl_int) :: fgsl_eigen_genv_qz
  fgsl_eigen_genv_qz = gsl_eigen_genv_qz(a%gsl_matrix, &
       b%gsl_matrix, alpha%gsl_vector_complex, &
       beta%gsl_vector,  evec%gsl_matrix_complex, &
       q%gsl_matrix, z%gsl_matrix, &
       w%gsl_eigen_genv_workspace)
end function fgsl_eigen_genv_qz
function fgsl_eigen_symmv_sort (eval, evec, sort_type)
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_matrix), intent(inout) :: evec 
  integer(fgsl_int), intent(in) :: sort_type
  integer(fgsl_int) :: fgsl_eigen_symmv_sort
  fgsl_eigen_symmv_sort = gsl_eigen_symmv_sort(eval%gsl_vector, &
       evec%gsl_matrix, sort_type)
end function fgsl_eigen_symmv_sort
function fgsl_eigen_hermv_sort (eval, evec, sort_type)
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_matrix_complex), intent(inout) :: evec 
  integer(fgsl_int), intent(in) :: sort_type
  integer(fgsl_int) :: fgsl_eigen_hermv_sort
  fgsl_eigen_hermv_sort = gsl_eigen_hermv_sort(eval%gsl_vector, &
       evec%gsl_matrix_complex, sort_type)
end function fgsl_eigen_hermv_sort
function fgsl_eigen_nonsymmv_sort (eval, evec, sort_type)
  type(fgsl_vector_complex), intent(inout) :: eval
  type(fgsl_matrix_complex), intent(inout) :: evec 
  integer(fgsl_int), intent(in) :: sort_type
  integer(fgsl_int) :: fgsl_eigen_nonsymmv_sort
  fgsl_eigen_nonsymmv_sort = gsl_eigen_nonsymmv_sort( &
       eval%gsl_vector_complex, evec%gsl_matrix_complex, sort_type)
end function fgsl_eigen_nonsymmv_sort
function fgsl_eigen_gensymmv_sort (eval, evec, sort_type)
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_matrix), intent(inout) :: evec 
  integer(fgsl_int), intent(in) :: sort_type
  integer(fgsl_int) :: fgsl_eigen_gensymmv_sort
  fgsl_eigen_gensymmv_sort = gsl_eigen_gensymmv_sort(eval%gsl_vector, &
       evec%gsl_matrix, sort_type)
end function fgsl_eigen_gensymmv_sort
function fgsl_eigen_genhermv_sort (eval, evec, sort_type)
  type(fgsl_vector), intent(inout) :: eval
  type(fgsl_matrix_complex), intent(inout) :: evec 
  integer(fgsl_int), intent(in) :: sort_type
  integer(fgsl_int) :: fgsl_eigen_genhermv_sort
  fgsl_eigen_genhermv_sort = gsl_eigen_genhermv_sort(eval%gsl_vector, &
       evec%gsl_matrix_complex, sort_type)
end function fgsl_eigen_genhermv_sort
function fgsl_eigen_genv_sort (alpha, beta, evec, sort_type)
  type(fgsl_vector_complex), intent(inout) :: alpha
  type(fgsl_vector), intent(inout) :: beta
  type(fgsl_matrix_complex), intent(inout) :: evec 
  integer(fgsl_int), intent(in) :: sort_type
  integer(fgsl_int) :: fgsl_eigen_genv_sort
  fgsl_eigen_genv_sort = gsl_eigen_genv_sort(alpha%gsl_vector_complex, &
       beta%gsl_vector, evec%gsl_matrix_complex, sort_type)
end function fgsl_eigen_genv_sort
!-*-f90-*-
!
!  API: Fast FT support
!
function fgsl_fft_complex_radix2_forward(data, stride, n)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int) :: fgsl_fft_complex_radix2_forward
  fgsl_fft_complex_radix2_forward = gsl_fft_complex_radix2_forward( &
       c_loc(data), stride, n)
end function fgsl_fft_complex_radix2_forward
function fgsl_fft_complex_radix2_transform(data, stride, n, sign)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int), intent(in) :: sign
  integer(fgsl_int) :: fgsl_fft_complex_radix2_transform
  fgsl_fft_complex_radix2_transform = gsl_fft_complex_radix2_transform( &
       c_loc(data), stride, n,sign)
end function fgsl_fft_complex_radix2_transform
function fgsl_fft_complex_radix2_backward(data, stride, n)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int) :: fgsl_fft_complex_radix2_backward
  fgsl_fft_complex_radix2_backward = gsl_fft_complex_radix2_backward( &
       c_loc(data), stride, n)
end function fgsl_fft_complex_radix2_backward
function fgsl_fft_complex_radix2_inverse(data, stride, n)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int) :: fgsl_fft_complex_radix2_inverse
  fgsl_fft_complex_radix2_inverse = gsl_fft_complex_radix2_inverse( &
       c_loc(data), stride, n)
end function fgsl_fft_complex_radix2_inverse
function fgsl_fft_complex_radix2_dif_forward(data, stride, n)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int) :: fgsl_fft_complex_radix2_dif_forward
  fgsl_fft_complex_radix2_dif_forward = gsl_fft_complex_radix2_dif_forward( &
       c_loc(data), stride, n)
end function fgsl_fft_complex_radix2_dif_forward
function fgsl_fft_complex_radix2_dif_transform(data, stride, n, sign)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int), intent(in) :: sign
  integer(fgsl_int) :: fgsl_fft_complex_radix2_dif_transform
  fgsl_fft_complex_radix2_dif_transform = gsl_fft_complex_radix2_dif_transform( &
       c_loc(data), stride, n,sign)
end function fgsl_fft_complex_radix2_dif_transform
function fgsl_fft_complex_radix2_dif_backward(data, stride, n)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int) :: fgsl_fft_complex_radix2_dif_backward
  fgsl_fft_complex_radix2_dif_backward = gsl_fft_complex_radix2_dif_backward( &
       c_loc(data), stride, n)
end function fgsl_fft_complex_radix2_dif_backward
function fgsl_fft_complex_radix2_dif_inverse(data, stride, n)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int) :: fgsl_fft_complex_radix2_dif_inverse
  fgsl_fft_complex_radix2_dif_inverse = gsl_fft_complex_radix2_dif_inverse( &
       c_loc(data), stride, n)
end function fgsl_fft_complex_radix2_dif_inverse
function fgsl_fft_complex_wavetable_alloc(n)
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_fft_complex_wavetable) :: fgsl_fft_complex_wavetable_alloc
  fgsl_fft_complex_wavetable_alloc%gsl_fft_complex_wavetable = &
       gsl_fft_complex_wavetable_alloc(n)
end function fgsl_fft_complex_wavetable_alloc
subroutine fgsl_fft_complex_wavetable_free(w) 
  type(fgsl_fft_complex_wavetable) :: w
  call gsl_fft_complex_wavetable_free(w%gsl_fft_complex_wavetable)
end subroutine fgsl_fft_complex_wavetable_free
function fgsl_fft_complex_workspace_alloc(n)
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_fft_complex_workspace) :: fgsl_fft_complex_workspace_alloc
  fgsl_fft_complex_workspace_alloc%gsl_fft_complex_workspace = &
       gsl_fft_complex_workspace_alloc(n)
end function fgsl_fft_complex_workspace_alloc
subroutine fgsl_fft_complex_workspace_free(w) 
  type(fgsl_fft_complex_workspace) :: w
  call gsl_fft_complex_workspace_free(w%gsl_fft_complex_workspace)
end subroutine fgsl_fft_complex_workspace_free
function fgsl_fft_complex_forward(data, stride, n, wavetable, work)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  type(fgsl_fft_complex_wavetable), intent(in) :: wavetable
  type(fgsl_fft_complex_workspace) :: work
  integer(fgsl_int) :: fgsl_fft_complex_forward
  fgsl_fft_complex_forward = gsl_fft_complex_forward( &
       c_loc(data), stride, n, wavetable%gsl_fft_complex_wavetable, &
       work%gsl_fft_complex_workspace)
end function fgsl_fft_complex_forward
function fgsl_fft_complex_transform(data, stride, n, wavetable, work, sign)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  type(fgsl_fft_complex_wavetable), intent(in) :: wavetable
  type(fgsl_fft_complex_workspace) :: work
  integer(fgsl_int), intent(in) :: sign
  integer(fgsl_int) :: fgsl_fft_complex_transform
  fgsl_fft_complex_transform = gsl_fft_complex_transform( &
       c_loc(data), stride, n, wavetable%gsl_fft_complex_wavetable, &
       work%gsl_fft_complex_workspace, sign)
end function fgsl_fft_complex_transform
function fgsl_fft_complex_backward(data, stride, n, wavetable, work)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  type(fgsl_fft_complex_wavetable), intent(in) :: wavetable
  type(fgsl_fft_complex_workspace) :: work
  integer(fgsl_int) :: fgsl_fft_complex_backward
  fgsl_fft_complex_backward = gsl_fft_complex_backward( &
       c_loc(data), stride, n, wavetable%gsl_fft_complex_wavetable, &
       work%gsl_fft_complex_workspace)
end function fgsl_fft_complex_backward
function fgsl_fft_complex_inverse(data, stride, n, wavetable, work)
  complex(fgsl_double_complex), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  type(fgsl_fft_complex_wavetable), intent(in) :: wavetable
  type(fgsl_fft_complex_workspace) :: work
  integer(fgsl_int) :: fgsl_fft_complex_inverse
  fgsl_fft_complex_inverse = gsl_fft_complex_inverse( &
       c_loc(data), stride, n, wavetable%gsl_fft_complex_wavetable, &
       work%gsl_fft_complex_workspace)
end function fgsl_fft_complex_inverse
function fgsl_fft_real_radix2_transform(data, stride, n)
  real(fgsl_double), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int) :: fgsl_fft_real_radix2_transform
  fgsl_fft_real_radix2_transform = gsl_fft_real_radix2_transform( &
       c_loc(data), stride, n)
end function fgsl_fft_real_radix2_transform
function fgsl_fft_halfcomplex_radix2_inverse(data, stride, n)
  real(fgsl_double), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int) :: fgsl_fft_halfcomplex_radix2_inverse
  fgsl_fft_halfcomplex_radix2_inverse = gsl_fft_halfcomplex_radix2_inverse( &
       c_loc(data), stride, n)
end function fgsl_fft_halfcomplex_radix2_inverse
function fgsl_fft_halfcomplex_radix2_backward(data, stride, n)
  real(fgsl_double), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int) :: fgsl_fft_halfcomplex_radix2_backward
  fgsl_fft_halfcomplex_radix2_backward = gsl_fft_halfcomplex_radix2_backward( &
       c_loc(data), stride, n)
end function fgsl_fft_halfcomplex_radix2_backward
function fgsl_fft_real_wavetable_alloc(n)
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_fft_real_wavetable) :: fgsl_fft_real_wavetable_alloc
  fgsl_fft_real_wavetable_alloc%gsl_fft_real_wavetable = &
       gsl_fft_real_wavetable_alloc(n)
end function fgsl_fft_real_wavetable_alloc
subroutine fgsl_fft_real_wavetable_free(w) 
  type(fgsl_fft_real_wavetable) :: w
  call gsl_fft_real_wavetable_free(w%gsl_fft_real_wavetable)
end subroutine fgsl_fft_real_wavetable_free
function fgsl_fft_halfcomplex_wavetable_alloc(n)
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_fft_halfcomplex_wavetable) :: fgsl_fft_halfcomplex_wavetable_alloc
  fgsl_fft_halfcomplex_wavetable_alloc%gsl_fft_halfcomplex_wavetable = &
       gsl_fft_halfcomplex_wavetable_alloc(n)
end function fgsl_fft_halfcomplex_wavetable_alloc
subroutine fgsl_fft_halfcomplex_wavetable_free(w) 
  type(fgsl_fft_halfcomplex_wavetable) :: w
  call gsl_fft_halfcomplex_wavetable_free(w%gsl_fft_halfcomplex_wavetable)
end subroutine fgsl_fft_halfcomplex_wavetable_free
function fgsl_fft_real_workspace_alloc(n)
  integer(fgsl_size_t), intent(in) :: n
  type(fgsl_fft_real_workspace) :: fgsl_fft_real_workspace_alloc
  fgsl_fft_real_workspace_alloc%gsl_fft_real_workspace = &
       gsl_fft_real_workspace_alloc(n)
end function fgsl_fft_real_workspace_alloc
subroutine fgsl_fft_real_workspace_free(w) 
  type(fgsl_fft_real_workspace) :: w
  call gsl_fft_real_workspace_free(w%gsl_fft_real_workspace)
end subroutine fgsl_fft_real_workspace_free
function fgsl_fft_real_transform(data, stride, n, wavetable, work)
  real(fgsl_double), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  type(fgsl_fft_real_wavetable), intent(in) :: wavetable
  type(fgsl_fft_real_workspace) :: work
  integer(fgsl_int) :: fgsl_fft_real_transform
  fgsl_fft_real_transform = gsl_fft_real_transform( &
       c_loc(data), stride, n, wavetable%gsl_fft_real_wavetable, &
       work%gsl_fft_real_workspace)
end function fgsl_fft_real_transform
function fgsl_fft_halfcomplex_transform(data, stride, n, wavetable, work)
  real(fgsl_double), intent(inout), target :: data(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  type(fgsl_fft_halfcomplex_wavetable), intent(in) :: wavetable
  type(fgsl_fft_real_workspace) :: work
  integer(fgsl_int) :: fgsl_fft_halfcomplex_transform
  fgsl_fft_halfcomplex_transform = gsl_fft_halfcomplex_transform( &
       c_loc(data), stride, n, wavetable%gsl_fft_halfcomplex_wavetable, &
       work%gsl_fft_real_workspace)
end function fgsl_fft_halfcomplex_transform
function fgsl_fft_real_unpack (real_coefficient, complex_coefficient, &
     stride, n)
  real(fgsl_double), intent(in), target :: real_coefficient(*)
  complex(fgsl_double_complex), intent(inout), target :: complex_coefficient(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int) :: fgsl_fft_real_unpack
  fgsl_fft_real_unpack = gsl_fft_real_unpack(c_loc(real_coefficient), &
       c_loc(complex_coefficient), stride, n)
end function fgsl_fft_real_unpack
function fgsl_fft_halfcomplex_unpack (halfcomplex_coefficient, &
     complex_coefficient, stride, n)
  real(fgsl_double), intent(in), target :: halfcomplex_coefficient(*)
  complex(fgsl_double_complex), intent(inout), target :: complex_coefficient(*)
  integer(fgsl_size_t), intent(in) :: stride, n
  integer(fgsl_int) :: fgsl_fft_halfcomplex_unpack
  fgsl_fft_halfcomplex_unpack = gsl_fft_halfcomplex_unpack( &
       c_loc(halfcomplex_coefficient), &
       c_loc(complex_coefficient), stride, n)
end function fgsl_fft_halfcomplex_unpack
!-*-f90-*-
!
! API: Numerical Integration
!
  function fgsl_integration_qng(f, a, b, epsabs, epsrel, result, abserr, neval) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: a, b, epsabs, epsrel
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_size_t), intent(inout) :: neval
    integer(fgsl_int) :: fgsl_integration_qng
    fgsl_integration_qng = gsl_integration_qng(f%gsl_function, a, b, epsabs, epsrel, &
         result, abserr, neval)
  end function fgsl_integration_qng
  function fgsl_integration_workspace_alloc(n) 
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_integration_workspace) :: fgsl_integration_workspace_alloc
    fgsl_integration_workspace_alloc%status = .false.
    fgsl_integration_workspace_alloc%gsl_integration_workspace = &
         gsl_integration_workspace_alloc(n)
    fgsl_integration_workspace_alloc%status = .true.
  end function fgsl_integration_workspace_alloc
  subroutine fgsl_integration_workspace_free(w) 
    type(fgsl_integration_workspace), intent(inout) :: w
    if (c_associated(w%gsl_integration_workspace)) &
         call gsl_integration_workspace_free (w%gsl_integration_workspace)
    w%status = .true.
  end subroutine fgsl_integration_workspace_free
  function fgsl_integration_qag(f, a, b, epsabs, epsrel, limit, key, &
         workspace, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: a, b, epsabs, epsrel
    integer(fgsl_size_t), intent(in) :: limit
    integer(fgsl_int), intent(in) :: key
    type(fgsl_integration_workspace), intent(inout) :: workspace
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_integration_qag
    fgsl_integration_qag = gsl_integration_qag(f%gsl_function, &
         a, b, epsabs, epsrel, limit, key, workspace%gsl_integration_workspace, &
         result, abserr) 
  end function fgsl_integration_qag
  function fgsl_integration_qags(f, a, b, epsabs, epsrel, limit, &
         workspace, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: a, b, epsabs, epsrel
    integer(fgsl_size_t), intent(in) :: limit
    type(fgsl_integration_workspace), intent(inout) :: workspace
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_integration_qags
    fgsl_integration_qags = gsl_integration_qags(f%gsl_function, &
         a, b, epsabs, epsrel, limit, workspace%gsl_integration_workspace, &
         result, abserr) 
  end function fgsl_integration_qags
  function fgsl_integration_qagp(f, pts, npts, epsabs, epsrel, limit, &
         workspace, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: pts(:)
    integer(fgsl_size_t), intent(in) :: npts
    real(fgsl_double), intent(in) :: epsabs, epsrel
    integer(fgsl_size_t), intent(in) :: limit
    type(fgsl_integration_workspace), intent(inout) :: workspace
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_integration_qagp
    fgsl_integration_qagp = gsl_integration_qagp(f%gsl_function, &
         pts, npts, epsabs, epsrel, limit, workspace%gsl_integration_workspace, &
         result, abserr) 
  end function fgsl_integration_qagp
  function fgsl_integration_qagi(f, epsabs, epsrel, limit, &
         workspace, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: epsabs, epsrel
    integer(fgsl_size_t), intent(in) :: limit
    type(fgsl_integration_workspace), intent(inout) :: workspace
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_integration_qagi
    fgsl_integration_qagi = gsl_integration_qagi(f%gsl_function, &
         epsabs, epsrel, limit, workspace%gsl_integration_workspace, &
         result, abserr) 
  end function fgsl_integration_qagi
  function fgsl_integration_qagiu(f, a, epsabs, epsrel, limit, &
         workspace, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: a, epsabs, epsrel
    integer(fgsl_size_t), intent(in) :: limit
    type(fgsl_integration_workspace), intent(inout) :: workspace
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_integration_qagiu
    fgsl_integration_qagiu = gsl_integration_qagiu(f%gsl_function, &
         a, epsabs, epsrel, limit, workspace%gsl_integration_workspace, &
         result, abserr) 
  end function fgsl_integration_qagiu
  function fgsl_integration_qagil(f, b, epsabs, epsrel, limit, &
         workspace, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: b, epsabs, epsrel
    integer(fgsl_size_t), intent(in) :: limit
    type(fgsl_integration_workspace), intent(inout) :: workspace
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_integration_qagil
    fgsl_integration_qagil = gsl_integration_qagil(f%gsl_function, &
         b, epsabs, epsrel, limit, workspace%gsl_integration_workspace, &
         result, abserr) 
  end function fgsl_integration_qagil
  function fgsl_integration_qawc(f, a, b, c, epsabs, epsrel, limit, &
         workspace, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: a, b, c, epsabs, epsrel
    integer(fgsl_size_t), intent(in) :: limit
    type(fgsl_integration_workspace), intent(inout) :: workspace
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_integration_qawc
    fgsl_integration_qawc = gsl_integration_qawc(f%gsl_function, &
         a, b, c, epsabs, epsrel, limit, workspace%gsl_integration_workspace, &
         result, abserr) 
  end function fgsl_integration_qawc
  function fgsl_integration_qaws_table_alloc (alpha, beta, mu, nu) 
    real(fgsl_double), intent(in) :: alpha, beta
    integer(fgsl_int), intent(in) :: mu, nu
    type(fgsl_integration_qaws_table) :: fgsl_integration_qaws_table_alloc
    fgsl_integration_qaws_table_alloc%status = .false.
    fgsl_integration_qaws_table_alloc%gsl_integration_qaws_table = &
         gsl_integration_qaws_table_alloc (alpha, beta, mu, nu) 
    if (c_associated(fgsl_integration_qaws_table_alloc%gsl_integration_qaws_table )) &
         fgsl_integration_qaws_table_alloc%status = .true.
  end function fgsl_integration_qaws_table_alloc
  function fgsl_integration_qaws_table_set (t, alpha, beta, mu, nu) 
    type(fgsl_integration_qaws_table) :: t
    real(fgsl_double), intent(in) :: alpha, beta
    integer(fgsl_int), intent(in) :: mu, nu
    integer(c_int) :: fgsl_integration_qaws_table_set
    fgsl_integration_qaws_table_set = &
         gsl_integration_qaws_table_set (t%gsl_integration_qaws_table, alpha, beta, mu, nu) 
  end function fgsl_integration_qaws_table_set
  subroutine fgsl_integration_qaws_table_free (w) 
    type(fgsl_integration_qaws_table), intent(inout) :: w
    if (c_associated(w%gsl_integration_qaws_table)) &
         call gsl_integration_qaws_table_free (w%gsl_integration_qaws_table)
    w%status = .false.
  end subroutine fgsl_integration_qaws_table_free
  function fgsl_integration_qaws(f, a, b, t, epsabs, epsrel, limit, workspace, &
       result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: a, b, epsabs, epsrel
    type(fgsl_integration_qaws_table), intent(in) :: t
    integer(fgsl_size_t), intent(in) :: limit
    type(fgsl_integration_workspace), intent(inout) :: workspace
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_integration_qaws
    fgsl_integration_qaws = gsl_integration_qaws(f%gsl_function, a, b, &
         t%gsl_integration_qaws_table, epsabs, epsrel, limit, &
         workspace%gsl_integration_workspace, result, abserr)
  end function fgsl_integration_qaws
  function fgsl_integration_qawo_table_alloc(omega, l, sine, n)
    real(fgsl_double), intent(in) :: omega, l
    integer(fgsl_int), intent(in) :: sine
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_integration_qawo_table) :: fgsl_integration_qawo_table_alloc
    fgsl_integration_qawo_table_alloc%status = .false. 
    fgsl_integration_qawo_table_alloc%gsl_integration_qawo_table = &
         gsl_integration_qawo_table_alloc(omega, l, sine, n)
    if (c_associated(fgsl_integration_qawo_table_alloc%gsl_integration_qawo_table)) &
         fgsl_integration_qawo_table_alloc%status = .true. 
  end function fgsl_integration_qawo_table_alloc
  function fgsl_integration_qawo_table_set(t, omega, l, sine) 
    type(fgsl_integration_qawo_table), intent(inout) :: t
    real(fgsl_double), intent(in) :: omega, l
    integer(fgsl_int), intent(in) :: sine
    integer(fgsl_int) :: fgsl_integration_qawo_table_set
    fgsl_integration_qawo_table_set = &
         gsl_integration_qawo_table_set(t%gsl_integration_qawo_table, omega, l, sine)
  end function fgsl_integration_qawo_table_set
  function fgsl_integration_qawo_table_set_length(t, l) 
    type(fgsl_integration_qawo_table), intent(inout) :: t
    real(fgsl_double), intent(in) :: l
    integer(fgsl_int) :: fgsl_integration_qawo_table_set_length
    fgsl_integration_qawo_table_set_length = &
         gsl_integration_qawo_table_set_length(t%gsl_integration_qawo_table, l)
  end function fgsl_integration_qawo_table_set_length
  subroutine fgsl_integration_qawo_table_free (w) 
    type(fgsl_integration_qawo_table), intent(inout) :: w
    if (c_associated(w%gsl_integration_qawo_table)) &
         call gsl_integration_qawo_table_free (w%gsl_integration_qawo_table)
    w%status = .false.
  end subroutine fgsl_integration_qawo_table_free
  function fgsl_integration_qawo (f, a, epsabs, epsrel, limit, workspace, &
         wf, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: a,  epsabs, epsrel
    integer(fgsl_size_t), intent(in) :: limit
    type(fgsl_integration_workspace), intent(inout) :: workspace
    type(fgsl_integration_qawo_table), intent(in) :: wf
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_integration_qawo
    fgsl_integration_qawo = gsl_integration_qawo(f%gsl_function, a, epsabs, epsrel, &
         limit, workspace%gsl_integration_workspace, wf%gsl_integration_qawo_table, &
         result, abserr) 
  end function fgsl_integration_qawo
  function fgsl_integration_qawf(f, a, epsabs, limit, workspace, cyc_workspace, &
       wf, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: a,  epsabs
    integer(fgsl_size_t), intent(in) :: limit
    type(fgsl_integration_workspace), intent(inout) :: workspace, cyc_workspace
    type(fgsl_integration_qawo_table), intent(in) :: wf
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_integration_qawf
    fgsl_integration_qawf = gsl_integration_qawf(f%gsl_function, a, epsabs, limit, &
         workspace%gsl_integration_workspace, cyc_workspace%gsl_integration_workspace, &
         wf%gsl_integration_qawo_table, result, abserr) 
  end function fgsl_integration_qawf
  function fgsl_integration_workspace_status(integration_workspace)
    type(fgsl_integration_workspace), intent(in) :: integration_workspace
    logical :: fgsl_integration_workspace_status
    fgsl_integration_workspace_status = .true.
    if (.not. c_associated(integration_workspace%gsl_integration_workspace)) &
         fgsl_integration_workspace_status = .false.
  end function fgsl_integration_workspace_status
  function fgsl_integration_qaws_table_status(integration_qaws_table)
    type(fgsl_integration_qaws_table), intent(in) :: integration_qaws_table
    logical :: fgsl_integration_qaws_table_status
    fgsl_integration_qaws_table_status = .true.
    if (.not. c_associated(&
         integration_qaws_table%gsl_integration_qaws_table)) &
         fgsl_integration_qaws_table_status = .false.
  end function fgsl_integration_qaws_table_status
  function fgsl_integration_qawo_table_status(integration_qawo_table)
    type(fgsl_integration_qawo_table), intent(in) :: integration_qawo_table
    logical :: fgsl_integration_qawo_table_status
    fgsl_integration_qawo_table_status = .true.
    if (.not. c_associated(&
         integration_qawo_table%gsl_integration_qawo_table)) &
         fgsl_integration_qawo_table_status = .false.
  end function fgsl_integration_qawo_table_status
  function fgsl_sizeof_integration_workspace(w)
    type(fgsl_integration_workspace), intent(in) :: w
    integer(fgsl_size_t) :: fgsl_sizeof_integration_workspace
    fgsl_sizeof_integration_workspace = gsl_aux_sizeof_integration_workspace()
  end function fgsl_sizeof_integration_workspace
  function fgsl_sizeof_integration_qaws_table(w)
    type(fgsl_integration_qaws_table), intent(in) :: w
    integer(fgsl_size_t) :: fgsl_sizeof_integration_qaws_table
    fgsl_sizeof_integration_qaws_table = &
         gsl_aux_sizeof_integration_qaws_table()
  end function fgsl_sizeof_integration_qaws_table
  function fgsl_sizeof_integration_qawo_table(w)
    type(fgsl_integration_qawo_table), intent(in) :: w
    integer(fgsl_size_t) :: fgsl_sizeof_integration_qawo_table
    fgsl_sizeof_integration_qawo_table = &
         gsl_aux_sizeof_integration_qawo_table()
  end function fgsl_sizeof_integration_qawo_table

!-*-f90-*-
!
! API: Random and Quasi-random numbers, distributions
!
  function fgsl_rng_alloc(t) 
    type(fgsl_rng_type), intent(inout) :: t
    type(fgsl_rng) :: fgsl_rng_alloc
    if (.not. c_associated(t%gsl_rng_type)) then
       t%gsl_rng_type = fgsl_aux_rng_assign(t%type)
    end if
    if (c_associated(t%gsl_rng_type)) &
         fgsl_rng_alloc%gsl_rng = gsl_rng_alloc(t%gsl_rng_type)
  end function fgsl_rng_alloc
  subroutine fgsl_rng_set(r, s)
    type(fgsl_rng), intent(inout) :: r
    integer(fgsl_long), intent(in) :: s
    call gsl_rng_set(r%gsl_rng, s)
  end subroutine fgsl_rng_set
  subroutine fgsl_rng_free(r)
    type(fgsl_rng), intent(inout) :: r
    call gsl_rng_free(r%gsl_rng)
  end subroutine fgsl_rng_free
  function fgsl_rng_get(r) 
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_long) :: fgsl_rng_get
    fgsl_rng_get = gsl_rng_get(r%gsl_rng) 
  end function fgsl_rng_get
  function fgsl_rng_uniform(r) 
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double) :: fgsl_rng_uniform
    fgsl_rng_uniform = gsl_rng_uniform(r%gsl_rng) 
  end function fgsl_rng_uniform
  function fgsl_rng_uniform_pos(r) 
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double) :: fgsl_rng_uniform_pos
    fgsl_rng_uniform_pos = gsl_rng_uniform_pos(r%gsl_rng) 
  end function fgsl_rng_uniform_pos
  function fgsl_rng_uniform_int(r, n) 
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_long), intent(in) :: n
    integer(fgsl_long) :: fgsl_rng_uniform_int
    fgsl_rng_uniform_int = gsl_rng_uniform_int(r%gsl_rng, n) 
  end function fgsl_rng_uniform_int
  function fgsl_rng_name(r) 
    type(fgsl_rng), intent(in) :: r
    character(kind=fgsl_char, len=fgsl_strmax) :: fgsl_rng_name
!
    type(c_ptr) :: name
    name = gsl_rng_name(r%gsl_rng)
    fgsl_rng_name = fgsl_name(name) 
  end function fgsl_rng_name
  function fgsl_rng_max(r) 
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_long) :: fgsl_rng_max
    fgsl_rng_max = gsl_rng_max(r%gsl_rng) 
  end function fgsl_rng_max
  function fgsl_rng_min(r) 
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_long) :: fgsl_rng_min
    fgsl_rng_min = gsl_rng_min(r%gsl_rng) 
  end function fgsl_rng_min
! FIXME: state, size and types_setup routines not (yet?) implemented
  function fgsl_rng_env_setup() 
    type(fgsl_rng_type) :: fgsl_rng_env_setup
    fgsl_rng_env_setup%gsl_rng_type = gsl_rng_env_setup() 
  end function fgsl_rng_env_setup
  function fgsl_rng_memcpy(cpy, src) 
    type(fgsl_rng), intent(inout) :: cpy
    type(fgsl_rng), intent(in) :: src
    integer(fgsl_int) :: fgsl_rng_memcpy
    fgsl_rng_memcpy = gsl_rng_memcpy(cpy%gsl_rng, src%gsl_rng) 
  end function fgsl_rng_memcpy
  function fgsl_rng_clone(r) 
    type(fgsl_rng), intent(in) :: r
    type(fgsl_rng) :: fgsl_rng_clone
    fgsl_rng_clone%gsl_rng = gsl_rng_clone(r%gsl_rng) 
  end function fgsl_rng_clone
  function fgsl_rng_fwrite(stream, r)
    type(fgsl_file), intent(in) :: stream
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_int) :: fgsl_rng_fwrite
    fgsl_rng_fwrite = gsl_rng_fwrite(stream%gsl_file, r%gsl_rng)
  end function fgsl_rng_fwrite
  function fgsl_rng_fread(stream, r)
    type(fgsl_file), intent(in) :: stream
    type(fgsl_rng), intent(inout) :: r
    integer(fgsl_int) :: fgsl_rng_fread
    fgsl_rng_fread = gsl_rng_fread(stream%gsl_file, r%gsl_rng)
  end function fgsl_rng_fread
  function fgsl_qrng_alloc(t, d) 
    type(fgsl_qrng_type), intent(in) :: t
    integer(fgsl_int), intent(in) :: d
    type(fgsl_qrng) :: fgsl_qrng_alloc
!
    type(c_ptr) :: type
    type = fgsl_aux_qrng_assign(t%type)
    fgsl_qrng_alloc%gsl_qrng = gsl_qrng_alloc(type, d)
  end function fgsl_qrng_alloc
  subroutine fgsl_qrng_free(r)
    type(fgsl_qrng), intent(inout) :: r
    call gsl_qrng_free(r%gsl_qrng)
  end subroutine fgsl_qrng_free
  subroutine fgsl_qrng_init(r)
    type(fgsl_qrng), intent(inout) :: r
    call gsl_qrng_init(r%gsl_qrng)
  end subroutine fgsl_qrng_init
  function fgsl_qrng_get(q, x)
    type(fgsl_qrng), intent(in) :: q
    real(fgsl_double), intent(out) :: x(:)
    integer(fgsl_int) :: fgsl_qrng_get
    fgsl_qrng_get = gsl_qrng_get(q%gsl_qrng, x)
  end function fgsl_qrng_get
  function fgsl_qrng_name(q) 
    type(fgsl_qrng), intent(in) :: q
    character(kind=fgsl_char, len=fgsl_strmax) :: fgsl_qrng_name
!
    type(c_ptr) :: name
    name = gsl_qrng_name(q%gsl_qrng)
    fgsl_qrng_name = fgsl_name(name) 
  end function fgsl_qrng_name
! FIXME: state and size routines not (yet?) implemented
  function fgsl_qrng_memcpy(cpy, src) 
    type(fgsl_qrng), intent(inout) :: cpy
    type(fgsl_qrng), intent(in) :: src
    integer(fgsl_int) :: fgsl_qrng_memcpy
    fgsl_qrng_memcpy = gsl_qrng_memcpy(cpy%gsl_qrng, src%gsl_qrng) 
  end function fgsl_qrng_memcpy
  function fgsl_qrng_clone(q) 
    type(fgsl_qrng), intent(in) :: q
    type(fgsl_qrng) :: fgsl_qrng_clone
    fgsl_qrng_clone%gsl_qrng = gsl_qrng_clone(q%gsl_qrng) 
  end function fgsl_qrng_clone
!
  function fgsl_ran_gaussian(r, sigma)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_ran_gaussian
    fgsl_ran_gaussian = gsl_ran_gaussian(r%gsl_rng, sigma)
  end function fgsl_ran_gaussian
  function fgsl_ran_gaussian_pdf(x, sigma)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_ran_gaussian_pdf
    fgsl_ran_gaussian_pdf = gsl_ran_gaussian_pdf(x, sigma)
  end function fgsl_ran_gaussian_pdf
  function fgsl_ran_gaussian_ziggurat(r, sigma)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_ran_gaussian_ziggurat
    fgsl_ran_gaussian_ziggurat = gsl_ran_gaussian_ziggurat(r%gsl_rng, sigma)
  end function fgsl_ran_gaussian_ziggurat
  function fgsl_ran_gaussian_ratio_method(r, sigma)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_ran_gaussian_ratio_method
    fgsl_ran_gaussian_ratio_method = gsl_ran_gaussian_ratio_method(r%gsl_rng, sigma)
  end function fgsl_ran_gaussian_ratio_method
  function fgsl_ran_ugaussian(r)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double) :: fgsl_ran_ugaussian
    fgsl_ran_ugaussian = gsl_ran_ugaussian(r%gsl_rng)
  end function fgsl_ran_ugaussian
  function fgsl_ran_ugaussian_pdf(x)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_ran_ugaussian_pdf
    fgsl_ran_ugaussian_pdf = gsl_ran_ugaussian_pdf(x)
  end function fgsl_ran_ugaussian_pdf
  function fgsl_ran_ugaussian_ratio_method(r)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double) :: fgsl_ran_ugaussian_ratio_method
    fgsl_ran_ugaussian_ratio_method = gsl_ran_ugaussian_ratio_method(r%gsl_rng)
  end function fgsl_ran_ugaussian_ratio_method
  function fgsl_cdf_gaussian_p(x, sigma)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_cdf_gaussian_p
    fgsl_cdf_gaussian_p = gsl_cdf_gaussian_p(x, sigma)
  end function fgsl_cdf_gaussian_p
  function fgsl_cdf_gaussian_q(x, sigma)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_cdf_gaussian_q
    fgsl_cdf_gaussian_q = gsl_cdf_gaussian_q(x, sigma)
  end function fgsl_cdf_gaussian_q
  function fgsl_cdf_gaussian_pinv(p, sigma)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_cdf_gaussian_pinv
    fgsl_cdf_gaussian_pinv = gsl_cdf_gaussian_pinv(p, sigma)
  end function fgsl_cdf_gaussian_pinv
  function fgsl_cdf_gaussian_qinv(q, sigma)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_cdf_gaussian_qinv
    fgsl_cdf_gaussian_qinv = gsl_cdf_gaussian_qinv(q, sigma)
  end function fgsl_cdf_gaussian_qinv
  function fgsl_cdf_ugaussian_p(x)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_cdf_ugaussian_p
    fgsl_cdf_ugaussian_p = gsl_cdf_ugaussian_p(x)
  end function fgsl_cdf_ugaussian_p
  function fgsl_cdf_ugaussian_q(x)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_cdf_ugaussian_q
    fgsl_cdf_ugaussian_q = gsl_cdf_ugaussian_q(x)
  end function fgsl_cdf_ugaussian_q
  function fgsl_cdf_ugaussian_pinv(p)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double) :: fgsl_cdf_ugaussian_pinv
    fgsl_cdf_ugaussian_pinv = gsl_cdf_ugaussian_pinv(p)
  end function fgsl_cdf_ugaussian_pinv
  function fgsl_cdf_ugaussian_qinv(q)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double) :: fgsl_cdf_ugaussian_qinv
    fgsl_cdf_ugaussian_qinv = gsl_cdf_ugaussian_qinv(q)
  end function fgsl_cdf_ugaussian_qinv
  function fgsl_ran_gaussian_tail(r, a, sigma)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a, sigma
    real(fgsl_double) :: fgsl_ran_gaussian_tail
    fgsl_ran_gaussian_tail = gsl_ran_gaussian_tail(r%gsl_rng, a, sigma)
  end function fgsl_ran_gaussian_tail
  function fgsl_ran_gaussian_tail_pdf(x, a, sigma)
    real(fgsl_double), intent(in) :: x, a
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_ran_gaussian_tail_pdf
    fgsl_ran_gaussian_tail_pdf = gsl_ran_gaussian_tail_pdf(x, a, sigma)
  end function fgsl_ran_gaussian_tail_pdf
  function fgsl_ran_ugaussian_tail(r, a)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_ran_ugaussian_tail
    fgsl_ran_ugaussian_tail = gsl_ran_ugaussian_tail(r%gsl_rng, a)
  end function fgsl_ran_ugaussian_tail
  function fgsl_ran_ugaussian_tail_pdf(x, a)
    real(fgsl_double), intent(in) :: x, a
    real(fgsl_double) :: fgsl_ran_ugaussian_tail_pdf
    fgsl_ran_ugaussian_tail_pdf = gsl_ran_ugaussian_tail_pdf(x, a)
  end function fgsl_ran_ugaussian_tail_pdf
  subroutine fgsl_ran_bivariate_gaussian(r, sigma_x, sigma_y, rho, x, y) 
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: sigma_x, sigma_y, rho
    real(fgsl_double), intent(out) :: x, y
    call gsl_ran_bivariate_gaussian(r%gsl_rng, sigma_x, sigma_y, rho, x, y) 
  end subroutine fgsl_ran_bivariate_gaussian
  function fgsl_ran_bivariate_gaussian_pdf(x, y, sigma_x, sigma_y, rho) 
    real(fgsl_double), intent(in) :: x, y, sigma_x, sigma_y, rho
    real(fgsl_double) :: fgsl_ran_bivariate_gaussian_pdf
    fgsl_ran_bivariate_gaussian_pdf = &
         gsl_ran_bivariate_gaussian_pdf(x, y, sigma_x, sigma_y, rho)
  end function fgsl_ran_bivariate_gaussian_pdf
  function fgsl_ran_exponential(r, mu)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: mu
    real(fgsl_double) :: fgsl_ran_exponential
    fgsl_ran_exponential = gsl_ran_exponential(r%gsl_rng, mu)
  end function fgsl_ran_exponential
  function fgsl_ran_exponential_pdf(x, mu)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: mu
    real(fgsl_double) :: fgsl_ran_exponential_pdf
    fgsl_ran_exponential_pdf = gsl_ran_exponential_pdf(x, mu)
  end function fgsl_ran_exponential_pdf
  function fgsl_cdf_exponential_p(x, mu)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: mu
    real(fgsl_double) :: fgsl_cdf_exponential_p
    fgsl_cdf_exponential_p = gsl_cdf_exponential_p(x, mu)
  end function fgsl_cdf_exponential_p
  function fgsl_cdf_exponential_q(x, mu)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: mu
    real(fgsl_double) :: fgsl_cdf_exponential_q
    fgsl_cdf_exponential_q = gsl_cdf_exponential_q(x, mu)
  end function fgsl_cdf_exponential_q
  function fgsl_cdf_exponential_pinv(p, mu)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: mu
    real(fgsl_double) :: fgsl_cdf_exponential_pinv
    fgsl_cdf_exponential_pinv = gsl_cdf_exponential_pinv(p, mu)
  end function fgsl_cdf_exponential_pinv
  function fgsl_cdf_exponential_qinv(q, mu)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: mu
    real(fgsl_double) :: fgsl_cdf_exponential_qinv
    fgsl_cdf_exponential_qinv = gsl_cdf_exponential_qinv(q, mu)
  end function fgsl_cdf_exponential_qinv
  function fgsl_ran_laplace(r, a)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_ran_laplace
    fgsl_ran_laplace = gsl_ran_laplace(r%gsl_rng, a)
  end function fgsl_ran_laplace
  function fgsl_ran_laplace_pdf(x, a)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_ran_laplace_pdf
    fgsl_ran_laplace_pdf = gsl_ran_laplace_pdf(x, a)
  end function fgsl_ran_laplace_pdf
  function fgsl_cdf_laplace_p(x, a)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_laplace_p
    fgsl_cdf_laplace_p = gsl_cdf_laplace_p(x, a)
  end function fgsl_cdf_laplace_p
  function fgsl_cdf_laplace_q(x, a)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_laplace_q
    fgsl_cdf_laplace_q = gsl_cdf_laplace_q(x, a)
  end function fgsl_cdf_laplace_q
  function fgsl_cdf_laplace_pinv(p, a)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_laplace_pinv
    fgsl_cdf_laplace_pinv = gsl_cdf_laplace_pinv(p, a)
  end function fgsl_cdf_laplace_pinv
  function fgsl_cdf_laplace_qinv(q, a)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_laplace_qinv
    fgsl_cdf_laplace_qinv = gsl_cdf_laplace_qinv(q, a)
  end function fgsl_cdf_laplace_qinv
  function fgsl_ran_exppow(r, a, b)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_exppow
    fgsl_ran_exppow = gsl_ran_exppow(r%gsl_rng, a, b)
  end function fgsl_ran_exppow
  function fgsl_ran_exppow_pdf(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_exppow_pdf
    fgsl_ran_exppow_pdf = gsl_ran_exppow_pdf(x, a, b)
  end function fgsl_ran_exppow_pdf
  function fgsl_cdf_exppow_p(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_exppow_p
    fgsl_cdf_exppow_p = gsl_cdf_exppow_p(x, a, b)
  end function fgsl_cdf_exppow_p
  function fgsl_cdf_exppow_q(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_exppow_q
    fgsl_cdf_exppow_q = gsl_cdf_exppow_q(x, a, b)
  end function fgsl_cdf_exppow_q
  function fgsl_ran_cauchy(r, a)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_ran_cauchy
    fgsl_ran_cauchy = gsl_ran_cauchy(r%gsl_rng, a)
  end function fgsl_ran_cauchy
  function fgsl_ran_cauchy_pdf(x, a)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_ran_cauchy_pdf
    fgsl_ran_cauchy_pdf = gsl_ran_cauchy_pdf(x, a)
  end function fgsl_ran_cauchy_pdf
  function fgsl_cdf_cauchy_p(x, a)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_cauchy_p
    fgsl_cdf_cauchy_p = gsl_cdf_cauchy_p(x, a)
  end function fgsl_cdf_cauchy_p
  function fgsl_cdf_cauchy_q(x, a)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_cauchy_q
    fgsl_cdf_cauchy_q = gsl_cdf_cauchy_q(x, a)
  end function fgsl_cdf_cauchy_q
  function fgsl_cdf_cauchy_pinv(p, a)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_cauchy_pinv
    fgsl_cdf_cauchy_pinv = gsl_cdf_cauchy_pinv(p, a)
  end function fgsl_cdf_cauchy_pinv
  function fgsl_cdf_cauchy_qinv(q, a)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_cauchy_qinv
    fgsl_cdf_cauchy_qinv = gsl_cdf_cauchy_qinv(q, a)
  end function fgsl_cdf_cauchy_qinv
  function fgsl_ran_rayleigh(r, sigma)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_ran_rayleigh
    fgsl_ran_rayleigh = gsl_ran_rayleigh(r%gsl_rng, sigma)
  end function fgsl_ran_rayleigh
  function fgsl_ran_rayleigh_pdf(x, sigma)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_ran_rayleigh_pdf
    fgsl_ran_rayleigh_pdf = gsl_ran_rayleigh_pdf(x, sigma)
  end function fgsl_ran_rayleigh_pdf
  function fgsl_cdf_rayleigh_p(x, sigma)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_cdf_rayleigh_p
    fgsl_cdf_rayleigh_p = gsl_cdf_rayleigh_p(x, sigma)
  end function fgsl_cdf_rayleigh_p
  function fgsl_cdf_rayleigh_q(x, sigma)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_cdf_rayleigh_q
    fgsl_cdf_rayleigh_q = gsl_cdf_rayleigh_q(x, sigma)
  end function fgsl_cdf_rayleigh_q
  function fgsl_cdf_rayleigh_pinv(p, sigma)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_cdf_rayleigh_pinv
    fgsl_cdf_rayleigh_pinv = gsl_cdf_rayleigh_pinv(p, sigma)
  end function fgsl_cdf_rayleigh_pinv
  function fgsl_cdf_rayleigh_qinv(q, sigma)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_cdf_rayleigh_qinv
    fgsl_cdf_rayleigh_qinv = gsl_cdf_rayleigh_qinv(q, sigma)
  end function fgsl_cdf_rayleigh_qinv
  function fgsl_ran_rayleigh_tail(r, a, sigma)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a, sigma
    real(fgsl_double) :: fgsl_ran_rayleigh_tail
    fgsl_ran_rayleigh_tail = gsl_ran_rayleigh_tail(r%gsl_rng, a, sigma)
  end function fgsl_ran_rayleigh_tail
  function fgsl_ran_rayleigh_tail_pdf(x, a, sigma)
    real(fgsl_double), intent(in) :: x, a
    real(fgsl_double), intent(in) :: sigma
    real(fgsl_double) :: fgsl_ran_rayleigh_tail_pdf
    fgsl_ran_rayleigh_tail_pdf = gsl_ran_rayleigh_tail_pdf(x, a, sigma)
  end function fgsl_ran_rayleigh_tail_pdf
  function fgsl_ran_landau(r)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double) :: fgsl_ran_landau
    fgsl_ran_landau = gsl_ran_landau(r%gsl_rng)
  end function fgsl_ran_landau
  function fgsl_ran_landau_pdf(x)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_ran_landau_pdf
    fgsl_ran_landau_pdf = gsl_ran_landau_pdf(x)
  end function fgsl_ran_landau_pdf
  function fgsl_ran_levy(r, c, alpha)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: c, alpha
    real(fgsl_double) :: fgsl_ran_levy
    fgsl_ran_levy = gsl_ran_levy(r%gsl_rng, c, alpha)
  end function fgsl_ran_levy
  function fgsl_ran_levy_skew(r, c, alpha, beta)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: c, alpha, beta
    real(fgsl_double) :: fgsl_ran_levy_skew
    fgsl_ran_levy_skew = gsl_ran_levy_skew(r%gsl_rng, c, alpha, beta)
  end function fgsl_ran_levy_skew
  function fgsl_ran_gamma(r, a, b)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_gamma
    fgsl_ran_gamma = gsl_ran_gamma(r%gsl_rng, a, b)
  end function fgsl_ran_gamma
  function fgsl_ran_gamma_mt(r, a, b)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_gamma_mt
    fgsl_ran_gamma_mt = gsl_ran_gamma_mt(r%gsl_rng, a, b)
  end function fgsl_ran_gamma_mt
  function fgsl_ran_gamma_pdf(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_gamma_pdf
    fgsl_ran_gamma_pdf = gsl_ran_gamma_pdf(x, a, b)
  end function fgsl_ran_gamma_pdf
  function fgsl_cdf_gamma_p(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gamma_p
    fgsl_cdf_gamma_p = gsl_cdf_gamma_p(x, a, b)
  end function fgsl_cdf_gamma_p
  function fgsl_cdf_gamma_q(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gamma_q
    fgsl_cdf_gamma_q = gsl_cdf_gamma_q(x, a, b)
  end function fgsl_cdf_gamma_q
  function fgsl_cdf_gamma_pinv(p, a, b)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gamma_pinv
    fgsl_cdf_gamma_pinv = gsl_cdf_gamma_pinv(p, a, b)
  end function fgsl_cdf_gamma_pinv
  function fgsl_cdf_gamma_qinv(q, a, b)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gamma_qinv
    fgsl_cdf_gamma_qinv = gsl_cdf_gamma_qinv(q, a, b)
  end function fgsl_cdf_gamma_qinv
  function fgsl_ran_flat(r, a, b)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_flat
    fgsl_ran_flat = gsl_ran_flat(r%gsl_rng, a, b)
  end function fgsl_ran_flat
  function fgsl_ran_flat_pdf(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_flat_pdf
    fgsl_ran_flat_pdf = gsl_ran_flat_pdf(x, a, b)
  end function fgsl_ran_flat_pdf
  function fgsl_cdf_flat_p(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_flat_p
    fgsl_cdf_flat_p = gsl_cdf_flat_p(x, a, b)
  end function fgsl_cdf_flat_p
  function fgsl_cdf_flat_q(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_flat_q
    fgsl_cdf_flat_q = gsl_cdf_flat_q(x, a, b)
  end function fgsl_cdf_flat_q
  function fgsl_cdf_flat_pinv(p, a, b)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_flat_pinv
    fgsl_cdf_flat_pinv = gsl_cdf_flat_pinv(p, a, b)
  end function fgsl_cdf_flat_pinv
  function fgsl_cdf_flat_qinv(q, a, b)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_flat_qinv
    fgsl_cdf_flat_qinv = gsl_cdf_flat_qinv(q, a, b)
  end function fgsl_cdf_flat_qinv
  function fgsl_ran_lognormal(r, zeta, sigma)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: zeta, sigma
    real(fgsl_double) :: fgsl_ran_lognormal
    fgsl_ran_lognormal = gsl_ran_lognormal(r%gsl_rng, zeta, sigma)
  end function fgsl_ran_lognormal
  function fgsl_ran_lognormal_pdf(x, zeta, sigma)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: zeta, sigma
    real(fgsl_double) :: fgsl_ran_lognormal_pdf
    fgsl_ran_lognormal_pdf = gsl_ran_lognormal_pdf(x, zeta, sigma)
  end function fgsl_ran_lognormal_pdf
  function fgsl_cdf_lognormal_p(x, zeta, sigma)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: zeta, sigma
    real(fgsl_double) :: fgsl_cdf_lognormal_p
    fgsl_cdf_lognormal_p = gsl_cdf_lognormal_p(x, zeta, sigma)
  end function fgsl_cdf_lognormal_p
  function fgsl_cdf_lognormal_q(x, zeta, sigma)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: zeta, sigma
    real(fgsl_double) :: fgsl_cdf_lognormal_q
    fgsl_cdf_lognormal_q = gsl_cdf_lognormal_q(x, zeta, sigma)
  end function fgsl_cdf_lognormal_q
  function fgsl_cdf_lognormal_pinv(p, zeta, sigma)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: zeta, sigma
    real(fgsl_double) :: fgsl_cdf_lognormal_pinv
    fgsl_cdf_lognormal_pinv = gsl_cdf_lognormal_pinv(p, zeta, sigma)
  end function fgsl_cdf_lognormal_pinv
  function fgsl_cdf_lognormal_qinv(q, zeta, sigma)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: zeta, sigma
    real(fgsl_double) :: fgsl_cdf_lognormal_qinv
    fgsl_cdf_lognormal_qinv = gsl_cdf_lognormal_qinv(q, zeta, sigma)
  end function fgsl_cdf_lognormal_qinv
  function fgsl_ran_chisq(r, nu)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_ran_chisq
    fgsl_ran_chisq = gsl_ran_chisq(r%gsl_rng, nu)
  end function fgsl_ran_chisq
  function fgsl_ran_chisq_pdf(x, nu)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_ran_chisq_pdf
    fgsl_ran_chisq_pdf = gsl_ran_chisq_pdf(x, nu)
  end function fgsl_ran_chisq_pdf
  function fgsl_cdf_chisq_p(x, nu)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_cdf_chisq_p
    fgsl_cdf_chisq_p = gsl_cdf_chisq_p(x, nu)
  end function fgsl_cdf_chisq_p
  function fgsl_cdf_chisq_q(x, nu)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_cdf_chisq_q
    fgsl_cdf_chisq_q = gsl_cdf_chisq_q(x, nu)
  end function fgsl_cdf_chisq_q
  function fgsl_cdf_chisq_pinv(p, nu)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_cdf_chisq_pinv
    fgsl_cdf_chisq_pinv = gsl_cdf_chisq_pinv(p, nu)
  end function fgsl_cdf_chisq_pinv
  function fgsl_cdf_chisq_qinv(q, nu)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_cdf_chisq_qinv
    fgsl_cdf_chisq_qinv = gsl_cdf_chisq_qinv(q, nu)
  end function fgsl_cdf_chisq_qinv
  function fgsl_ran_fdist(r, nu1, nu2)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: nu1, nu2
    real(fgsl_double) :: fgsl_ran_fdist
    fgsl_ran_fdist = gsl_ran_fdist(r%gsl_rng, nu1, nu2)
  end function fgsl_ran_fdist
  function fgsl_ran_fdist_pdf(x, nu1, nu2)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: nu1, nu2
    real(fgsl_double) :: fgsl_ran_fdist_pdf
    fgsl_ran_fdist_pdf = gsl_ran_fdist_pdf(x, nu1, nu2)
  end function fgsl_ran_fdist_pdf
  function fgsl_cdf_fdist_p(x, nu1, nu2)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: nu1, nu2
    real(fgsl_double) :: fgsl_cdf_fdist_p
    fgsl_cdf_fdist_p = gsl_cdf_fdist_p(x, nu1, nu2)
  end function fgsl_cdf_fdist_p
  function fgsl_cdf_fdist_q(x, nu1, nu2)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: nu1, nu2
    real(fgsl_double) :: fgsl_cdf_fdist_q
    fgsl_cdf_fdist_q = gsl_cdf_fdist_q(x, nu1, nu2)
  end function fgsl_cdf_fdist_q
  function fgsl_cdf_fdist_pinv(p, nu1, nu2)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: nu1, nu2
    real(fgsl_double) :: fgsl_cdf_fdist_pinv
    fgsl_cdf_fdist_pinv = gsl_cdf_fdist_pinv(p, nu1, nu2)
  end function fgsl_cdf_fdist_pinv
  function fgsl_cdf_fdist_qinv(q, nu1, nu2)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: nu1, nu2
    real(fgsl_double) :: fgsl_cdf_fdist_qinv
    fgsl_cdf_fdist_qinv = gsl_cdf_fdist_qinv(q, nu1, nu2)
  end function fgsl_cdf_fdist_qinv
  function fgsl_ran_tdist(r, nu)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_ran_tdist
    fgsl_ran_tdist = gsl_ran_tdist(r%gsl_rng, nu)
  end function fgsl_ran_tdist
  function fgsl_ran_tdist_pdf(x, nu)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_ran_tdist_pdf
    fgsl_ran_tdist_pdf = gsl_ran_tdist_pdf(x, nu)
  end function fgsl_ran_tdist_pdf
  function fgsl_cdf_tdist_p(x, nu)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_cdf_tdist_p
    fgsl_cdf_tdist_p = gsl_cdf_tdist_p(x, nu)
  end function fgsl_cdf_tdist_p
  function fgsl_cdf_tdist_q(x, nu)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_cdf_tdist_q
    fgsl_cdf_tdist_q = gsl_cdf_tdist_q(x, nu)
  end function fgsl_cdf_tdist_q
  function fgsl_cdf_tdist_pinv(p, nu)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_cdf_tdist_pinv
    fgsl_cdf_tdist_pinv = gsl_cdf_tdist_pinv(p, nu)
  end function fgsl_cdf_tdist_pinv
  function fgsl_cdf_tdist_qinv(q, nu)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: nu
    real(fgsl_double) :: fgsl_cdf_tdist_qinv
    fgsl_cdf_tdist_qinv = gsl_cdf_tdist_qinv(q, nu)
  end function fgsl_cdf_tdist_qinv
  function fgsl_ran_beta(r, a, b)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_beta
    fgsl_ran_beta = gsl_ran_beta(r%gsl_rng, a, b)
  end function fgsl_ran_beta
  function fgsl_ran_beta_pdf(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_beta_pdf
    fgsl_ran_beta_pdf = gsl_ran_beta_pdf(x, a, b)
  end function fgsl_ran_beta_pdf
  function fgsl_cdf_beta_p(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_beta_p
    fgsl_cdf_beta_p = gsl_cdf_beta_p(x, a, b)
  end function fgsl_cdf_beta_p
  function fgsl_cdf_beta_q(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_beta_q
    fgsl_cdf_beta_q = gsl_cdf_beta_q(x, a, b)
  end function fgsl_cdf_beta_q
  function fgsl_cdf_beta_pinv(p, a, b)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_beta_pinv
    fgsl_cdf_beta_pinv = gsl_cdf_beta_pinv(p, a, b)
  end function fgsl_cdf_beta_pinv
  function fgsl_cdf_beta_qinv(q, a, b)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_beta_qinv
    fgsl_cdf_beta_qinv = gsl_cdf_beta_qinv(q, a, b)
  end function fgsl_cdf_beta_qinv
  function fgsl_ran_logistic(r, a)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_ran_logistic
    fgsl_ran_logistic = gsl_ran_logistic(r%gsl_rng, a)
  end function fgsl_ran_logistic
  function fgsl_ran_logistic_pdf(x, a)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_ran_logistic_pdf
    fgsl_ran_logistic_pdf = gsl_ran_logistic_pdf(x, a)
  end function fgsl_ran_logistic_pdf
  function fgsl_cdf_logistic_p(x, a)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_logistic_p
    fgsl_cdf_logistic_p = gsl_cdf_logistic_p(x, a)
  end function fgsl_cdf_logistic_p
  function fgsl_cdf_logistic_q(x, a)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_logistic_q
    fgsl_cdf_logistic_q = gsl_cdf_logistic_q(x, a)
  end function fgsl_cdf_logistic_q
  function fgsl_cdf_logistic_pinv(p, a)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_logistic_pinv
    fgsl_cdf_logistic_pinv = gsl_cdf_logistic_pinv(p, a)
  end function fgsl_cdf_logistic_pinv
  function fgsl_cdf_logistic_qinv(q, a)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: a
    real(fgsl_double) :: fgsl_cdf_logistic_qinv
    fgsl_cdf_logistic_qinv = gsl_cdf_logistic_qinv(q, a)
  end function fgsl_cdf_logistic_qinv
  function fgsl_ran_pareto(r, a, b)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_pareto
    fgsl_ran_pareto = gsl_ran_pareto(r%gsl_rng, a, b)
  end function fgsl_ran_pareto
  function fgsl_ran_pareto_pdf(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_pareto_pdf
    fgsl_ran_pareto_pdf = gsl_ran_pareto_pdf(x, a, b)
  end function fgsl_ran_pareto_pdf
  function fgsl_cdf_pareto_p(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_pareto_p
    fgsl_cdf_pareto_p = gsl_cdf_pareto_p(x, a, b)
  end function fgsl_cdf_pareto_p
  function fgsl_cdf_pareto_q(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_pareto_q
    fgsl_cdf_pareto_q = gsl_cdf_pareto_q(x, a, b)
  end function fgsl_cdf_pareto_q
  function fgsl_cdf_pareto_pinv(p, a, b)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_pareto_pinv
    fgsl_cdf_pareto_pinv = gsl_cdf_pareto_pinv(p, a, b)
  end function fgsl_cdf_pareto_pinv
  function fgsl_cdf_pareto_qinv(q, a, b)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_pareto_qinv
    fgsl_cdf_pareto_qinv = gsl_cdf_pareto_qinv(q, a, b)
  end function fgsl_cdf_pareto_qinv
  subroutine fgsl_ran_dir_2d(r, x, y) 
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(out) :: x, y
    call gsl_ran_dir_2d(r%gsl_rng, x, y) 
  end subroutine fgsl_ran_dir_2d
  subroutine fgsl_ran_dir_2d_trig_method(r, x, y) 
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(out) :: x, y
    call gsl_ran_dir_2d_trig_method(r%gsl_rng, x, y) 
  end subroutine fgsl_ran_dir_2d_trig_method
  subroutine fgsl_ran_dir_3d(r, x, y, z) 
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(out) :: x, y, z
    call gsl_ran_dir_3d(r%gsl_rng, x, y, z)
  end subroutine fgsl_ran_dir_3d
  subroutine fgsl_ran_dir_nd(r, n, x) 
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_size_t), intent(in) :: n
    real(fgsl_double), intent(out) :: x
    call gsl_ran_dir_nd(r%gsl_rng, n, x)
  end subroutine fgsl_ran_dir_nd
  function fgsl_ran_weibull(r, a, b)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_weibull
    fgsl_ran_weibull = gsl_ran_weibull(r%gsl_rng, a, b)
  end function fgsl_ran_weibull
  function fgsl_ran_weibull_pdf(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_weibull_pdf
    fgsl_ran_weibull_pdf = gsl_ran_weibull_pdf(x, a, b)
  end function fgsl_ran_weibull_pdf
  function fgsl_cdf_weibull_p(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_weibull_p
    fgsl_cdf_weibull_p = gsl_cdf_weibull_p(x, a, b)
  end function fgsl_cdf_weibull_p
  function fgsl_cdf_weibull_q(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_weibull_q
    fgsl_cdf_weibull_q = gsl_cdf_weibull_q(x, a, b)
  end function fgsl_cdf_weibull_q
  function fgsl_cdf_weibull_pinv(p, a, b)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_weibull_pinv
    fgsl_cdf_weibull_pinv = gsl_cdf_weibull_pinv(p, a, b)
  end function fgsl_cdf_weibull_pinv
  function fgsl_cdf_weibull_qinv(q, a, b)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_weibull_qinv
    fgsl_cdf_weibull_qinv = gsl_cdf_weibull_qinv(q, a, b)
  end function fgsl_cdf_weibull_qinv
  function fgsl_ran_gumbel1(r, a, b)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_gumbel1
    fgsl_ran_gumbel1 = gsl_ran_gumbel1(r%gsl_rng, a, b)
  end function fgsl_ran_gumbel1
  function fgsl_ran_gumbel1_pdf(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_gumbel1_pdf
    fgsl_ran_gumbel1_pdf = gsl_ran_gumbel1_pdf(x, a, b)
  end function fgsl_ran_gumbel1_pdf
  function fgsl_cdf_gumbel1_p(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gumbel1_p
    fgsl_cdf_gumbel1_p = gsl_cdf_gumbel1_p(x, a, b)
  end function fgsl_cdf_gumbel1_p
  function fgsl_cdf_gumbel1_q(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gumbel1_q
    fgsl_cdf_gumbel1_q = gsl_cdf_gumbel1_q(x, a, b)
  end function fgsl_cdf_gumbel1_q
  function fgsl_cdf_gumbel1_pinv(p, a, b)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gumbel1_pinv
    fgsl_cdf_gumbel1_pinv = gsl_cdf_gumbel1_pinv(p, a, b)
  end function fgsl_cdf_gumbel1_pinv
  function fgsl_cdf_gumbel1_qinv(q, a, b)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gumbel1_qinv
    fgsl_cdf_gumbel1_qinv = gsl_cdf_gumbel1_qinv(q, a, b)
  end function fgsl_cdf_gumbel1_qinv
  function fgsl_ran_gumbel2(r, a, b)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_gumbel2
    fgsl_ran_gumbel2 = gsl_ran_gumbel2(r%gsl_rng, a, b)
  end function fgsl_ran_gumbel2
  function fgsl_ran_gumbel2_pdf(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_ran_gumbel2_pdf
    fgsl_ran_gumbel2_pdf = gsl_ran_gumbel2_pdf(x, a, b)
  end function fgsl_ran_gumbel2_pdf
  function fgsl_cdf_gumbel2_p(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gumbel2_p
    fgsl_cdf_gumbel2_p = gsl_cdf_gumbel2_p(x, a, b)
  end function fgsl_cdf_gumbel2_p
  function fgsl_cdf_gumbel2_q(x, a, b)
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gumbel2_q
    fgsl_cdf_gumbel2_q = gsl_cdf_gumbel2_q(x, a, b)
  end function fgsl_cdf_gumbel2_q
  function fgsl_cdf_gumbel2_pinv(p, a, b)
    real(fgsl_double), intent(in) :: p
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gumbel2_pinv
    fgsl_cdf_gumbel2_pinv = gsl_cdf_gumbel2_pinv(p, a, b)
  end function fgsl_cdf_gumbel2_pinv
  function fgsl_cdf_gumbel2_qinv(q, a, b)
    real(fgsl_double), intent(in) :: q
    real(fgsl_double), intent(in) :: a, b
    real(fgsl_double) :: fgsl_cdf_gumbel2_qinv
    fgsl_cdf_gumbel2_qinv = gsl_cdf_gumbel2_qinv(q, a, b)
  end function fgsl_cdf_gumbel2_qinv
  subroutine fgsl_ran_dirichlet(r, k, alpha, theta) 
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_size_t), intent(in) :: k
    real(fgsl_double), intent(in) :: alpha(:)
    real(fgsl_double), intent(out) :: theta(:)
    call gsl_ran_dirichlet(r%gsl_rng, k, alpha, theta) 
  end subroutine fgsl_ran_dirichlet
  function fgsl_ran_dirichlet_pdf(k, alpha, theta) 
    integer(fgsl_size_t), intent(in) :: k
    real(fgsl_double), intent(in) :: alpha(:)
    real(fgsl_double), intent(in) :: theta(:)
    real(fgsl_double) :: fgsl_ran_dirichlet_pdf
    fgsl_ran_dirichlet_pdf = gsl_ran_dirichlet_pdf(k, alpha, theta) 
  end function fgsl_ran_dirichlet_pdf
  function fgsl_ran_dirichlet_lnpdf(k, alpha, theta) 
    integer(fgsl_size_t), intent(in) :: k
    real(fgsl_double), intent(in) :: alpha(:)
    real(fgsl_double), intent(in) :: theta(:)
    real(fgsl_double) :: fgsl_ran_dirichlet_lnpdf
    fgsl_ran_dirichlet_lnpdf = gsl_ran_dirichlet_lnpdf(k, alpha, theta) 
  end function fgsl_ran_dirichlet_lnpdf
  function fgsl_ran_discrete_preproc(k, p)
    integer(fgsl_size_t), intent(in) :: k
    real(fgsl_double), intent(in) :: p(:)
    type(fgsl_ran_discrete_t) :: fgsl_ran_discrete_preproc
    fgsl_ran_discrete_preproc%gsl_ran_discrete_t = &
         gsl_ran_discrete_preproc(k, p)
  end function fgsl_ran_discrete_preproc
  function fgsl_ran_discrete(r, g)
    type(fgsl_rng), intent(in) :: r
    type(fgsl_ran_discrete_t), intent(in) :: g
    integer(fgsl_size_t) :: fgsl_ran_discrete
    fgsl_ran_discrete = gsl_ran_discrete(r%gsl_rng, g%gsl_ran_discrete_t)
  end function fgsl_ran_discrete
  function fgsl_ran_discrete_pdf(k, g) 
    integer(fgsl_size_t), intent(in) :: k
    type(fgsl_ran_discrete_t), intent(in) :: g
    real(fgsl_double) :: fgsl_ran_discrete_pdf
    fgsl_ran_discrete_pdf = gsl_ran_discrete_pdf(k, g%gsl_ran_discrete_t) 
  end function fgsl_ran_discrete_pdf
  subroutine fgsl_ran_discrete_free(g)
    type(fgsl_ran_discrete_t), intent(inout) :: g
    call gsl_ran_discrete_free(g%gsl_ran_discrete_t)
  end subroutine fgsl_ran_discrete_free
  function fgsl_ran_poisson(r, mu)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: mu
    integer(fgsl_int) :: fgsl_ran_poisson
    fgsl_ran_poisson = gsl_ran_poisson(r%gsl_rng, mu)
  end function fgsl_ran_poisson
  function fgsl_ran_poisson_pdf(k, mu)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: mu
    real(fgsl_double) :: fgsl_ran_poisson_pdf
    fgsl_ran_poisson_pdf = gsl_ran_poisson_pdf(k, mu)
  end function fgsl_ran_poisson_pdf
  function fgsl_cdf_poisson_p(k, mu)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: mu
    real(fgsl_double) :: fgsl_cdf_poisson_p
    fgsl_cdf_poisson_p = gsl_cdf_poisson_p(k, mu)
  end function fgsl_cdf_poisson_p
  function fgsl_cdf_poisson_q(k, mu)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: mu
    real(fgsl_double) :: fgsl_cdf_poisson_q
    fgsl_cdf_poisson_q = gsl_cdf_poisson_q(k, mu)
  end function fgsl_cdf_poisson_q
  function fgsl_ran_bernoulli(r, p)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: p
    integer(fgsl_int) :: fgsl_ran_bernoulli
    fgsl_ran_bernoulli = gsl_ran_bernoulli(r%gsl_rng, p)
  end function fgsl_ran_bernoulli
  function fgsl_ran_bernoulli_pdf(k, p)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: p
    real(fgsl_double) :: fgsl_ran_bernoulli_pdf
    fgsl_ran_bernoulli_pdf = gsl_ran_bernoulli_pdf(k, p)
  end function fgsl_ran_bernoulli_pdf
  function fgsl_ran_binomial(r, p, n)
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double), intent(in) :: p
    real(fgsl_double) :: fgsl_ran_binomial
    fgsl_ran_binomial = gsl_ran_binomial(r%gsl_rng, p, n)
  end function fgsl_ran_binomial
  function fgsl_ran_binomial_pdf(k, p, n)
    integer(fgsl_int), intent(in) :: k, n
    real(fgsl_double), intent(in) :: p
    real(fgsl_double) :: fgsl_ran_binomial_pdf
    fgsl_ran_binomial_pdf = gsl_ran_binomial_pdf(k, p, n)
  end function fgsl_ran_binomial_pdf
  function fgsl_cdf_binomial_p(k, p, n)
    integer(fgsl_int), intent(in) :: k, n
    real(fgsl_double), intent(in) :: p
    real(fgsl_double) :: fgsl_cdf_binomial_p
    fgsl_cdf_binomial_p = gsl_cdf_binomial_p(k, p, n)
  end function fgsl_cdf_binomial_p
  function fgsl_cdf_binomial_q(k, p, n)
    integer(fgsl_int), intent(in) :: k, n
    real(fgsl_double), intent(in) :: p
    real(fgsl_double) :: fgsl_cdf_binomial_q
    fgsl_cdf_binomial_q = gsl_cdf_binomial_q(k, p, n)
  end function fgsl_cdf_binomial_q
  subroutine fgsl_ran_multinomial(r, k, nn, p, n)
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_size_t), intent(in) :: k
    integer(fgsl_int), intent(in) :: nn
    real(fgsl_double), intent(in) :: p(:)
    integer(fgsl_int), intent(out) :: n(:)
    call gsl_ran_multinomial(r%gsl_rng, k, nn, p, n)
  end subroutine fgsl_ran_multinomial
  function fgsl_ran_multinomial_pdf(k, p, n)
    integer(fgsl_size_t), intent(in) :: k
    real(fgsl_double), intent(in) :: p(:)
    integer(fgsl_int), intent(in) :: n(:)
    real(fgsl_double) :: fgsl_ran_multinomial_pdf
    fgsl_ran_multinomial_pdf = gsl_ran_multinomial_pdf(k, p, n)
  end function fgsl_ran_multinomial_pdf
  function fgsl_ran_multinomial_lnpdf(k, p, n)
    integer(fgsl_size_t), intent(in) :: k
    real(fgsl_double), intent(in) :: p(:)
    integer(fgsl_int), intent(in) :: n(:)
    real(fgsl_double) :: fgsl_ran_multinomial_lnpdf
    fgsl_ran_multinomial_lnpdf = gsl_ran_multinomial_lnpdf(k, p, n)
  end function fgsl_ran_multinomial_lnpdf
  function fgsl_ran_negative_binomial(r, p, n)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: p, n
    integer(fgsl_int) :: fgsl_ran_negative_binomial
    fgsl_ran_negative_binomial = gsl_ran_negative_binomial(r%gsl_rng, p, n)
  end function fgsl_ran_negative_binomial
  function fgsl_ran_negative_binomial_pdf(k, p, n)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: p, n
    real(fgsl_double) :: fgsl_ran_negative_binomial_pdf
    fgsl_ran_negative_binomial_pdf = gsl_ran_negative_binomial_pdf(k, p, n)
  end function fgsl_ran_negative_binomial_pdf
  function fgsl_cdf_negative_binomial_p(k, p, n)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: p, n
    real(fgsl_double) :: fgsl_cdf_negative_binomial_p
    fgsl_cdf_negative_binomial_p = gsl_cdf_negative_binomial_p(k, p, n)
  end function fgsl_cdf_negative_binomial_p
  function fgsl_cdf_negative_binomial_q(k, p, n)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: p, n
    real(fgsl_double) :: fgsl_cdf_negative_binomial_q
    fgsl_cdf_negative_binomial_q = gsl_cdf_negative_binomial_q(k, p, n)
  end function fgsl_cdf_negative_binomial_q
  function fgsl_ran_pascal(r, p, n)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: p, n
    integer(fgsl_int) :: fgsl_ran_pascal
    fgsl_ran_pascal = gsl_ran_pascal(r%gsl_rng, p, n)
  end function fgsl_ran_pascal
  function fgsl_ran_pascal_pdf(k, p, n)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: p, n
    real(fgsl_double) :: fgsl_ran_pascal_pdf
    fgsl_ran_pascal_pdf = gsl_ran_pascal_pdf(k, p, n)
  end function fgsl_ran_pascal_pdf
  function fgsl_cdf_pascal_p(k, p, n)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: p, n
    real(fgsl_double) :: fgsl_cdf_pascal_p
    fgsl_cdf_pascal_p = gsl_cdf_pascal_p(k, p, n)
  end function fgsl_cdf_pascal_p
  function fgsl_cdf_pascal_q(k, p, n)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: p, n
    real(fgsl_double) :: fgsl_cdf_pascal_q
    fgsl_cdf_pascal_q = gsl_cdf_pascal_q(k, p, n)
  end function fgsl_cdf_pascal_q
  function fgsl_ran_geometric(r, p)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: p
    integer(fgsl_int) :: fgsl_ran_geometric
    fgsl_ran_geometric = gsl_ran_geometric(r%gsl_rng, p)
  end function fgsl_ran_geometric
  function fgsl_ran_geometric_pdf(k, p)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: p
    real(fgsl_double) :: fgsl_ran_geometric_pdf
    fgsl_ran_geometric_pdf = gsl_ran_geometric_pdf(k, p)
  end function fgsl_ran_geometric_pdf
  function fgsl_cdf_geometric_p(k, p)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: p
    real(fgsl_double) :: fgsl_cdf_geometric_p
    fgsl_cdf_geometric_p = gsl_cdf_geometric_p(k, p)
  end function fgsl_cdf_geometric_p
  function fgsl_cdf_geometric_q(k, p)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: p
    real(fgsl_double) :: fgsl_cdf_geometric_q
    fgsl_cdf_geometric_q = gsl_cdf_geometric_q(k, p)
  end function fgsl_cdf_geometric_q
  function fgsl_ran_hypergeometric(r, n1, n2, t) 
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_int), intent(in) :: n1, n2, t
    integer(fgsl_int) :: fgsl_ran_hypergeometric
    fgsl_ran_hypergeometric = gsl_ran_hypergeometric(r%gsl_rng, n1, n2, t)
  end function fgsl_ran_hypergeometric
  function fgsl_ran_hypergeometric_pdf(k, n1, n2, t) 
    integer(fgsl_int), intent(in) :: k
    integer(fgsl_int), intent(in) :: n1, n2, t
    real(fgsl_double)  :: fgsl_ran_hypergeometric_pdf
    fgsl_ran_hypergeometric_pdf = gsl_ran_hypergeometric_pdf(k, n1, n2, t)
  end function fgsl_ran_hypergeometric_pdf
  function fgsl_cdf_hypergeometric_p(k, n1, n2, t) 
    integer(fgsl_int), intent(in) :: k
    integer(fgsl_int), intent(in) :: n1, n2, t
    real(fgsl_double)  :: fgsl_cdf_hypergeometric_p
    fgsl_cdf_hypergeometric_p = gsl_cdf_hypergeometric_p(k, n1, n2, t)
  end function fgsl_cdf_hypergeometric_p
  function fgsl_cdf_hypergeometric_q(k, n1, n2, t) 
    integer(fgsl_int), intent(in) :: k
    integer(fgsl_int), intent(in) :: n1, n2, t
    real(fgsl_double)  :: fgsl_cdf_hypergeometric_q
    fgsl_cdf_hypergeometric_q = gsl_cdf_hypergeometric_q(k, n1, n2, t)
  end function fgsl_cdf_hypergeometric_q
  function fgsl_ran_logarithmic(r, p)
    type(fgsl_rng), intent(in) :: r
    real(fgsl_double), intent(in) :: p
    integer(fgsl_int) :: fgsl_ran_logarithmic
    fgsl_ran_logarithmic = gsl_ran_logarithmic(r%gsl_rng, p)
  end function fgsl_ran_logarithmic
  function fgsl_ran_logarithmic_pdf(k, p)
    integer(fgsl_int), intent(in) :: k
    real(fgsl_double), intent(in) :: p
    real(fgsl_double) :: fgsl_ran_logarithmic_pdf
    fgsl_ran_logarithmic_pdf = gsl_ran_logarithmic_pdf(k, p)
  end function fgsl_ran_logarithmic_pdf
  subroutine fgsl_ran_shuffle(r, base, n, size)
    type(fgsl_rng), intent(in) :: r
    type(c_ptr), intent(in) :: base
    integer(fgsl_size_t), intent(in) :: n, size
    call gsl_ran_shuffle(r%gsl_rng, base, n, size)
  end subroutine fgsl_ran_shuffle
  subroutine fgsl_ran_shuffle_double(r, base, n)
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_size_t), intent(in) :: n
    real(fgsl_double), target, intent(in) :: base(n)
!
    integer(fgsl_size_t) :: size
    type(c_ptr) :: ptr_base
    ptr_base = c_loc(base)
    size = fgsl_sizeof(1.0_fgsl_double)
    call gsl_ran_shuffle(r%gsl_rng, ptr_base, n, size)
  end subroutine fgsl_ran_shuffle_double
  subroutine fgsl_ran_shuffle_size_t(r, base, n)
    type(fgsl_rng), intent(in) :: r
    integer(fgsl_size_t), intent(in) :: n
    integer(fgsl_size_t), target, intent(in) :: base(n)
!
    integer(fgsl_size_t) :: size
    type(c_ptr) :: ptr_base
    ptr_base = c_loc(base)
    size = fgsl_sizeof(1_fgsl_size_t)
    call gsl_ran_shuffle(r%gsl_rng, ptr_base, n, size)
  end subroutine fgsl_ran_shuffle_size_t
  function fgsl_ran_choose(r, dest, k, src, n, size) 
    type(fgsl_rng), intent(in) :: r
    type(c_ptr), intent(in) :: dest, src
    integer(fgsl_size_t), intent(in) :: k, n, size
    integer(fgsl_int) :: fgsl_ran_choose
    fgsl_ran_choose = gsl_ran_choose(r%gsl_rng, dest, k, src, n, size)
  end function fgsl_ran_choose
  subroutine fgsl_ran_sample(r, dest, k, src, n, size) 
    type(fgsl_rng), intent(in) :: r
    type(c_ptr), intent(in) :: dest, src
    integer(fgsl_size_t), intent(in) :: k, n, size
    call gsl_ran_sample(r%gsl_rng, dest, k, src, n, size)
  end subroutine fgsl_ran_sample
! Add-ons
  subroutine fgsl_rng_c_ptr(res, src) 
    type(c_ptr), intent(in) :: src
    type(fgsl_rng), intent(out) :: res
    res%gsl_rng = src
  end subroutine fgsl_rng_c_ptr
  function fgsl_rng_status(rng)
    type(fgsl_rng), intent(in) :: rng
    logical :: fgsl_rng_status
    fgsl_rng_status = .true.
    if (.not. c_associated(rng%gsl_rng)) fgsl_rng_status = .false.
  end function fgsl_rng_status
  function fgsl_qrng_status(qrng)
    type(fgsl_qrng), intent(in) :: qrng
    logical :: fgsl_qrng_status
    fgsl_qrng_status = .true.
    if (.not. c_associated(qrng%gsl_qrng)) fgsl_qrng_status = .false.
  end function fgsl_qrng_status
  function fgsl_ran_discrete_t_status(ran_discrete_t)
    type(fgsl_ran_discrete_t), intent(in) :: ran_discrete_t
    logical :: fgsl_ran_discrete_t_status
    fgsl_ran_discrete_t_status = .true.
    if (.not. c_associated(ran_discrete_t%gsl_ran_discrete_t)) &
         fgsl_ran_discrete_t_status = .false.
  end function fgsl_ran_discrete_t_status


 



 



 



 

!-*-f90-*-
!
! API: Statistics
!
  function fgsl_stats_mean (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double) :: fgsl_stats_mean
    fgsl_stats_mean = gsl_stats_mean (data, stride, n) 
  end function fgsl_stats_mean
  function fgsl_stats_variance (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double) :: fgsl_stats_variance
    fgsl_stats_variance = gsl_stats_variance (data, stride, n) 
  end function fgsl_stats_variance
  function fgsl_stats_variance_m (data, stride, n, mean) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double), intent(in) :: mean
    real(fgsl_double) :: fgsl_stats_variance_m
    fgsl_stats_variance_m = gsl_stats_variance_m (data, stride, n, mean) 
  end function fgsl_stats_variance_m
  function fgsl_stats_sd (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double) :: fgsl_stats_sd
    fgsl_stats_sd = gsl_stats_sd (data, stride, n) 
  end function fgsl_stats_sd
  function fgsl_stats_sd_m (data, stride, n, mean) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double), intent(in) :: mean
    real(fgsl_double) :: fgsl_stats_sd_m
    fgsl_stats_sd_m = gsl_stats_sd_m (data, stride, n, mean) 
  end function fgsl_stats_sd_m
  function fgsl_stats_variance_with_fixed_mean (data, stride, n, mean) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double), intent(in) :: mean
    real(fgsl_double) :: fgsl_stats_variance_with_fixed_mean
    fgsl_stats_variance_with_fixed_mean = &
         gsl_stats_variance_with_fixed_mean (data, stride, n, mean) 
  end function fgsl_stats_variance_with_fixed_mean
  function fgsl_stats_sd_with_fixed_mean (data, stride, n, mean) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double), intent(in) :: mean
    real(fgsl_double) :: fgsl_stats_sd_with_fixed_mean
    fgsl_stats_sd_with_fixed_mean = &
         gsl_stats_sd_with_fixed_mean (data, stride, n, mean) 
  end function fgsl_stats_sd_with_fixed_mean
  function fgsl_stats_absdev (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double) :: fgsl_stats_absdev
    fgsl_stats_absdev = gsl_stats_absdev (data, stride, n) 
  end function fgsl_stats_absdev
  function fgsl_stats_absdev_m (data, stride, n, mean) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double), intent(in) :: mean
    real(fgsl_double) :: fgsl_stats_absdev_m
    fgsl_stats_absdev_m = gsl_stats_absdev_m (data, stride, n, mean) 
  end function fgsl_stats_absdev_m
  function fgsl_stats_skew (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double) :: fgsl_stats_skew
    fgsl_stats_skew = gsl_stats_skew (data, stride, n) 
  end function fgsl_stats_skew
  function fgsl_stats_skew_m_sd (data, stride, n, mean, sd) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double), intent(in) :: mean, sd
    real(fgsl_double) :: fgsl_stats_skew_m_sd
    fgsl_stats_skew_m_sd = gsl_stats_skew_m_sd (data, stride, n, mean, sd) 
  end function fgsl_stats_skew_m_sd
  function fgsl_stats_kurtosis (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double) :: fgsl_stats_kurtosis
    fgsl_stats_kurtosis = gsl_stats_kurtosis (data, stride, n) 
  end function fgsl_stats_kurtosis
  function fgsl_stats_kurtosis_m_sd (data, stride, n, mean, sd) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double), intent(in) :: mean, sd
    real(fgsl_double) :: fgsl_stats_kurtosis_m_sd
    fgsl_stats_kurtosis_m_sd = gsl_stats_kurtosis_m_sd (data, stride, n, mean, sd) 
  end function fgsl_stats_kurtosis_m_sd
  function fgsl_stats_lag1_autocorrelation (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double) :: fgsl_stats_lag1_autocorrelation
    fgsl_stats_lag1_autocorrelation = gsl_stats_lag1_autocorrelation (data, stride, n) 
  end function fgsl_stats_lag1_autocorrelation
  function fgsl_stats_lag1_autocorrelation_m (data, stride, n, mean) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double), intent(in) :: mean
    real(fgsl_double) :: fgsl_stats_lag1_autocorrelation_m
    fgsl_stats_lag1_autocorrelation_m = &
         gsl_stats_lag1_autocorrelation_m (data, stride, n, mean) 
  end function fgsl_stats_lag1_autocorrelation_m
  function fgsl_stats_covariance(data1, stride1, data2, stride2, n) 
    real(fgsl_double), intent(in) :: data1(:), data2(:)
    integer(fgsl_size_t), intent(in) :: stride1, stride2, n
    real(fgsl_double) :: fgsl_stats_covariance
    fgsl_stats_covariance = gsl_stats_covariance(data1, stride1, data2, stride2, n)
  end function fgsl_stats_covariance
  function fgsl_stats_covariance_m(data1, stride1, data2, stride2, n, mean1, mean2) 
    real(fgsl_double), intent(in) :: data1(:), data2(:)
    integer(fgsl_size_t), intent(in) :: stride1, stride2, n
    real(fgsl_double), intent(in) :: mean1, mean2
    real(fgsl_double) :: fgsl_stats_covariance_m
    fgsl_stats_covariance_m = &
         gsl_stats_covariance_m(data1, stride1, data2, stride2, n, mean1, mean2)
  end function fgsl_stats_covariance_m
  function fgsl_stats_correlation(data1, stride1, data2, stride2, n) 
    real(fgsl_double), intent(in) :: data1(:), data2(:)
    integer(fgsl_size_t), intent(in) :: stride1, stride2, n
    real(fgsl_double) :: fgsl_stats_correlation
    fgsl_stats_correlation = gsl_stats_correlation(data1, stride1, data2, stride2, n)
  end function fgsl_stats_correlation
  function fgsl_stats_wmean(w, wstride, data, stride, n) 
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double) :: fgsl_stats_wmean
    fgsl_stats_wmean = gsl_stats_wmean(w, wstride, data, stride, n)
  end function fgsl_stats_wmean
  function fgsl_stats_wvariance(w, wstride, data, stride, n) 
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double) :: fgsl_stats_wvariance
    fgsl_stats_wvariance = gsl_stats_wvariance(w, wstride, data, stride, n)
  end function fgsl_stats_wvariance
  function fgsl_stats_wvariance_m(w, wstride, data, stride, n, mean)
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double), intent(in) :: mean
    real(fgsl_double) :: fgsl_stats_wvariance_m
    fgsl_stats_wvariance_m = gsl_stats_wvariance_m(w, wstride, data, stride, n, mean)
  end function fgsl_stats_wvariance_m
  function fgsl_stats_wsd(w, wstride, data, stride, n) 
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double) :: fgsl_stats_wsd
    fgsl_stats_wsd = gsl_stats_wsd(w, wstride, data, stride, n)
  end function fgsl_stats_wsd
  function fgsl_stats_wsd_m(w, wstride, data, stride, n, mean)
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double), intent(in) :: mean
    real(fgsl_double) :: fgsl_stats_wsd_m
    fgsl_stats_wsd_m = gsl_stats_wsd_m(w, wstride, data, stride, n, mean)
  end function fgsl_stats_wsd_m
  function fgsl_stats_wvariance_with_fixed_mean(w, wstride, data, stride, n, mean)
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double), intent(in) :: mean
    real(fgsl_double) :: fgsl_stats_wvariance_with_fixed_mean
    fgsl_stats_wvariance_with_fixed_mean = &
         gsl_stats_wvariance_with_fixed_mean(w, wstride, data, stride, n, mean)
  end function fgsl_stats_wvariance_with_fixed_mean
  function fgsl_stats_wsd_with_fixed_mean(w, wstride, data, stride, n, mean)
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double), intent(in) :: mean
    real(fgsl_double) :: fgsl_stats_wsd_with_fixed_mean
    fgsl_stats_wsd_with_fixed_mean = &
         gsl_stats_wsd_with_fixed_mean(w, wstride, data, stride, n, mean)
  end function fgsl_stats_wsd_with_fixed_mean
  function fgsl_stats_wabsdev(w, wstride, data, stride, n) 
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double) :: fgsl_stats_wabsdev
    fgsl_stats_wabsdev = gsl_stats_wabsdev(w, wstride, data, stride, n)
  end function fgsl_stats_wabsdev
  function fgsl_stats_wabsdev_m(w, wstride, data, stride, n, mean)
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double), intent(in) :: mean
    real(fgsl_double) :: fgsl_stats_wabsdev_m
    fgsl_stats_wabsdev_m = gsl_stats_wabsdev_m(w, wstride, data, stride, n, mean)
  end function fgsl_stats_wabsdev_m
  function fgsl_stats_wskew(w, wstride, data, stride, n) 
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double) :: fgsl_stats_wskew
    fgsl_stats_wskew = gsl_stats_wskew(w, wstride, data, stride, n)
  end function fgsl_stats_wskew
  function fgsl_stats_wskew_m_sd(w, wstride, data, stride, n, mean, sd)
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double), intent(in) :: mean, sd
    real(fgsl_double) :: fgsl_stats_wskew_m_sd
    fgsl_stats_wskew_m_sd = gsl_stats_wskew_m_sd(w, wstride, data, stride, n, mean, sd)
  end function fgsl_stats_wskew_m_sd
  function fgsl_stats_wkurtosis(w, wstride, data, stride, n) 
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double) :: fgsl_stats_wkurtosis
    fgsl_stats_wkurtosis = gsl_stats_wkurtosis(w, wstride, data, stride, n)
  end function fgsl_stats_wkurtosis
  function fgsl_stats_wkurtosis_m_sd(w, wstride, data, stride, n, mean, sd)
    real(fgsl_double), intent(in) :: w(:), data(:)
    integer(fgsl_size_t), intent(in) :: wstride, stride, n
    real(fgsl_double), intent(in) :: mean, sd
    real(fgsl_double) :: fgsl_stats_wkurtosis_m_sd
    fgsl_stats_wkurtosis_m_sd = &
         gsl_stats_wkurtosis_m_sd(w, wstride, data, stride, n, mean, sd)
  end function fgsl_stats_wkurtosis_m_sd
  function fgsl_stats_max (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double) :: fgsl_stats_max
    fgsl_stats_max = gsl_stats_max (data, stride, n) 
  end function fgsl_stats_max
  function fgsl_stats_min (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double) :: fgsl_stats_min
    fgsl_stats_min = gsl_stats_min (data, stride, n) 
  end function fgsl_stats_min
  subroutine fgsl_stats_minmax (min, max, data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double), intent(out) :: min, max
    call gsl_stats_minmax(min, max, data, stride, n) 
  end subroutine fgsl_stats_minmax
  function fgsl_stats_max_index (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    integer(fgsl_size_t) :: fgsl_stats_max_index
    fgsl_stats_max_index = gsl_stats_max_index (data, stride, n) 
  end function fgsl_stats_max_index
  function fgsl_stats_min_index (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    integer(fgsl_size_t) :: fgsl_stats_min_index
    fgsl_stats_min_index = gsl_stats_min_index (data, stride, n) 
  end function fgsl_stats_min_index
  subroutine fgsl_stats_minmax_index (min_index, max_index, data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    integer(fgsl_size_t), intent(out) :: min_index, max_index
    call gsl_stats_minmax_index(min_index, max_index, data, stride, n) 
  end subroutine fgsl_stats_minmax_index
  function fgsl_stats_median_from_sorted_data (data, stride, n) 
    real(fgsl_double), intent(in) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double) :: fgsl_stats_median_from_sorted_data
    fgsl_stats_median_from_sorted_data = &
         gsl_stats_median_from_sorted_data (data, stride, n) 
  end function fgsl_stats_median_from_sorted_data
  function fgsl_stats_quantile_from_sorted_data (data, stride, n, f) 
    real(fgsl_double), intent(in) :: data(:)
    real(fgsl_double), intent(in) :: f
    integer(fgsl_size_t), intent(in) :: stride, n
    real(fgsl_double) :: fgsl_stats_quantile_from_sorted_data
    fgsl_stats_quantile_from_sorted_data = &
         gsl_stats_quantile_from_sorted_data (data, stride, n, f) 
  end function fgsl_stats_quantile_from_sorted_data
!-*-f90-*-
!
! API: Histograms
!
  function fgsl_histogram_alloc(n)
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_histogram) :: fgsl_histogram_alloc
    fgsl_histogram_alloc%gsl_histogram = gsl_histogram_alloc(n)
  end function fgsl_histogram_alloc
  function fgsl_histogram_set_ranges(h, range, size) 
    type(fgsl_histogram), intent(inout) :: h
    integer(fgsl_size_t), intent(in) :: size
    real(fgsl_double), intent(in) :: range(:)
    integer(fgsl_int) :: fgsl_histogram_set_ranges
    fgsl_histogram_set_ranges = gsl_histogram_set_ranges(h%gsl_histogram, range, size) 
  end function fgsl_histogram_set_ranges
  function fgsl_histogram_set_ranges_uniform(h, xmin, xmax)
    type(fgsl_histogram), intent(inout) :: h
    real(fgsl_double), intent(in) :: xmin, xmax
    integer(fgsl_int) :: fgsl_histogram_set_ranges_uniform
    fgsl_histogram_set_ranges_uniform = &
         gsl_histogram_set_ranges_uniform(h%gsl_histogram, xmin, xmax)
  end function fgsl_histogram_set_ranges_uniform
  subroutine fgsl_histogram_free(h) 
    type(fgsl_histogram), intent(inout) :: h
    call gsl_histogram_free(h%gsl_histogram) 
  end subroutine fgsl_histogram_free
  function fgsl_histogram_memcpy(dest, src) 
    type(fgsl_histogram), intent(inout) :: dest
    type(fgsl_histogram), intent(in) :: src
    integer(fgsl_int) :: fgsl_histogram_memcpy
    fgsl_histogram_memcpy = gsl_histogram_memcpy(dest%gsl_histogram, src%gsl_histogram) 
  end function fgsl_histogram_memcpy
  function fgsl_histogram_clone(src) 
    type(fgsl_histogram), intent(in) :: src
    type(fgsl_histogram) :: fgsl_histogram_clone
    fgsl_histogram_clone%gsl_histogram = gsl_histogram_clone(src%gsl_histogram) 
  end function fgsl_histogram_clone
  function fgsl_histogram_increment(h, x) 
    type(fgsl_histogram), intent(inout) :: h
    real(fgsl_double), intent(in) :: x
    integer(fgsl_int) :: fgsl_histogram_increment
    fgsl_histogram_increment = gsl_histogram_increment(h%gsl_histogram, x) 
  end function fgsl_histogram_increment
  function fgsl_histogram_accumulate(h, x, weight) 
    type(fgsl_histogram), intent(inout) :: h
    real(fgsl_double), intent(in) :: x, weight
    integer(fgsl_int) :: fgsl_histogram_accumulate
    fgsl_histogram_accumulate = gsl_histogram_accumulate(h%gsl_histogram, x, weight) 
  end function fgsl_histogram_accumulate
  function fgsl_histogram_get(h, i)
    type(fgsl_histogram), intent(in) :: h
    integer(fgsl_size_t), intent(in) :: i
    real(fgsl_double) :: fgsl_histogram_get
    fgsl_histogram_get = gsl_histogram_get(h%gsl_histogram, i)
  end function fgsl_histogram_get
  function fgsl_histogram_get_range(h, i, lower, upper) 
    type(fgsl_histogram), intent(in) :: h
    integer(fgsl_size_t), intent(in) :: i
    real(fgsl_double), intent(out) :: lower, upper
    integer(fgsl_int) :: fgsl_histogram_get_range
    fgsl_histogram_get_range = &
         gsl_histogram_get_range(h%gsl_histogram, i, lower, upper) 
  end function fgsl_histogram_get_range
  function fgsl_histogram_max(h)
    type(fgsl_histogram), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram_max
    fgsl_histogram_max = gsl_histogram_max(h%gsl_histogram)
  end function fgsl_histogram_max
  function fgsl_histogram_min(h)
    type(fgsl_histogram), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram_min
    fgsl_histogram_min = gsl_histogram_min(h%gsl_histogram)
  end function fgsl_histogram_min
  function fgsl_histogram_bins(h)
    type(fgsl_histogram), intent(in) :: h
    integer(fgsl_size_t) :: fgsl_histogram_bins
    fgsl_histogram_bins = gsl_histogram_bins(h%gsl_histogram)
  end function fgsl_histogram_bins
  subroutine fgsl_histogram_reset(h) 
    type(fgsl_histogram), intent(inout) :: h
    call gsl_histogram_reset(h%gsl_histogram) 
  end subroutine fgsl_histogram_reset
  function fgsl_histogram_find(h, x, i) 
    type(fgsl_histogram), intent(in) :: h
    real(fgsl_double), intent(in) :: x    
    integer(fgsl_size_t), intent(out) :: i
    integer(fgsl_int) :: fgsl_histogram_find
    fgsl_histogram_find = gsl_histogram_find(h%gsl_histogram, x, i) 
  end function fgsl_histogram_find
  function fgsl_histogram_max_val(h)
    type(fgsl_histogram), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram_max_val
    fgsl_histogram_max_val = gsl_histogram_max_val(h%gsl_histogram)
  end function fgsl_histogram_max_val
  function fgsl_histogram_max_bin(h)
    type(fgsl_histogram), intent(in) :: h
    integer(fgsl_size_t) :: fgsl_histogram_max_bin
    fgsl_histogram_max_bin = gsl_histogram_max_bin(h%gsl_histogram)
  end function fgsl_histogram_max_bin
  function fgsl_histogram_min_val(h)
    type(fgsl_histogram), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram_min_val
    fgsl_histogram_min_val = gsl_histogram_min_val(h%gsl_histogram)
  end function fgsl_histogram_min_val
  function fgsl_histogram_min_bin(h)
    type(fgsl_histogram), intent(in) :: h
    integer(fgsl_size_t) :: fgsl_histogram_min_bin
    fgsl_histogram_min_bin = gsl_histogram_min_bin(h%gsl_histogram)
  end function fgsl_histogram_min_bin
  function fgsl_histogram_mean(h)
    type(fgsl_histogram), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram_mean
    fgsl_histogram_mean = gsl_histogram_mean(h%gsl_histogram)
  end function fgsl_histogram_mean
  function fgsl_histogram_sigma(h)
    type(fgsl_histogram), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram_sigma
    fgsl_histogram_sigma = gsl_histogram_sigma(h%gsl_histogram)
  end function fgsl_histogram_sigma
  function fgsl_histogram_sum(h)
    type(fgsl_histogram), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram_sum
    fgsl_histogram_sum = gsl_histogram_sum(h%gsl_histogram)
  end function fgsl_histogram_sum
  function fgsl_histogram_equal_bins_p(h1, h2)
    type(fgsl_histogram), intent(in) :: h1, h2
    real(fgsl_double) :: fgsl_histogram_equal_bins_p
    fgsl_histogram_equal_bins_p = & 
         gsl_histogram_equal_bins_p(h1%gsl_histogram, h2%gsl_histogram)
  end function fgsl_histogram_equal_bins_p
  function fgsl_histogram_add(h1, h2)
    type(fgsl_histogram), intent(inout) :: h1
    type(fgsl_histogram), intent(in) :: h2
    real(fgsl_double) :: fgsl_histogram_add
    fgsl_histogram_add = & 
         gsl_histogram_add(h1%gsl_histogram, h2%gsl_histogram)
  end function fgsl_histogram_add
  function fgsl_histogram_sub(h1, h2)
    type(fgsl_histogram), intent(inout) :: h1
    type(fgsl_histogram), intent(in) :: h2
    real(fgsl_double) :: fgsl_histogram_sub
    fgsl_histogram_sub = & 
         gsl_histogram_sub(h1%gsl_histogram, h2%gsl_histogram)
  end function fgsl_histogram_sub
  function fgsl_histogram_mul(h1, h2)
    type(fgsl_histogram), intent(inout) :: h1
    type(fgsl_histogram), intent(in) :: h2
    real(fgsl_double) :: fgsl_histogram_mul
    fgsl_histogram_mul = & 
         gsl_histogram_mul(h1%gsl_histogram, h2%gsl_histogram)
  end function fgsl_histogram_mul
  function fgsl_histogram_div(h1, h2)
    type(fgsl_histogram), intent(inout) :: h1
    type(fgsl_histogram), intent(in) :: h2
    real(fgsl_double) :: fgsl_histogram_div
    fgsl_histogram_div = & 
         gsl_histogram_div(h1%gsl_histogram, h2%gsl_histogram)
  end function fgsl_histogram_div
  function fgsl_histogram_scale(h, scale) 
    type(fgsl_histogram), intent(inout) :: h
    real(fgsl_double), intent(in) :: scale
    integer(fgsl_int) :: fgsl_histogram_scale
    fgsl_histogram_scale = gsl_histogram_scale(h%gsl_histogram, scale) 
  end function fgsl_histogram_scale
  function fgsl_histogram_shift(h, offset) 
    type(fgsl_histogram), intent(inout) :: h
    real(fgsl_double), intent(in) :: offset
    integer(fgsl_int) :: fgsl_histogram_shift
    fgsl_histogram_shift = gsl_histogram_shift(h%gsl_histogram, offset) 
  end function fgsl_histogram_shift
  function fgsl_histogram_fwrite(stream, h)
    type(fgsl_file), intent(in) :: stream
    type(fgsl_histogram), intent(in) :: h
    integer(fgsl_int) :: fgsl_histogram_fwrite
    fgsl_histogram_fwrite = gsl_histogram_fwrite(stream%gsl_file, h%gsl_histogram)
  end function fgsl_histogram_fwrite
  function fgsl_histogram_fread(stream, h)
    type(fgsl_file), intent(in) :: stream
    type(fgsl_histogram), intent(inout) :: h
    integer(fgsl_int) :: fgsl_histogram_fread
    fgsl_histogram_fread = gsl_histogram_fread(stream%gsl_file, h%gsl_histogram)
  end function fgsl_histogram_fread
  function fgsl_histogram_fprintf(stream, h, range_format, bin_format) 
    type(fgsl_file), intent(in) :: stream
    type(fgsl_histogram), intent(in) :: h
    character(kind=fgsl_char, len=*), intent(in) :: range_format, bin_format
    integer(fgsl_int) :: fgsl_histogram_fprintf
!
    character(kind=fgsl_char,len=fgsl_strmax) :: lrf, lbf
    if (len(trim(range_format)) < fgsl_strmax .and. &
         len(trim(bin_format)) < fgsl_strmax) then
       lrf = trim(range_format) // c_null_char 
       lbf = trim(bin_format) // c_null_char 
       fgsl_histogram_fprintf = &
            gsl_histogram_fprintf(stream%gsl_file, h%gsl_histogram, lrf, lbf)
    else
       fgsl_histogram_fprintf = fgsl_failure
    end if
  end function fgsl_histogram_fprintf
  function fgsl_histogram_fscanf(stream, h) 
    type(fgsl_file), intent(in) :: stream
    type(fgsl_histogram), intent(inout) :: h
    integer(fgsl_int) :: fgsl_histogram_fscanf
    fgsl_histogram_fscanf = gsl_histogram_fscanf(stream%gsl_file, h%gsl_histogram)
  end function fgsl_histogram_fscanf
  function fgsl_histogram_pdf_alloc(n)
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_histogram_pdf) :: fgsl_histogram_pdf_alloc
    fgsl_histogram_pdf_alloc%gsl_histogram_pdf = gsl_histogram_pdf_alloc(n)
  end function fgsl_histogram_pdf_alloc
  function fgsl_histogram_pdf_init(p, h) 
    type(fgsl_histogram_pdf), intent(inout) :: p
    type(fgsl_histogram), intent(in) :: h
    integer(fgsl_int) :: fgsl_histogram_pdf_init
    fgsl_histogram_pdf_init = gsl_histogram_pdf_init(p%gsl_histogram_pdf, h%gsl_histogram) 
  end function fgsl_histogram_pdf_init
  subroutine fgsl_histogram_pdf_free(p) 
    type(fgsl_histogram_pdf), intent(inout) :: p
    call gsl_histogram_pdf_free(p%gsl_histogram_pdf) 
  end subroutine fgsl_histogram_pdf_free
  function fgsl_histogram_pdf_sample(p, r) 
    type(fgsl_histogram_pdf), intent(in) :: p
    real(fgsl_double), intent(in) :: r
    real(fgsl_double) :: fgsl_histogram_pdf_sample
    fgsl_histogram_pdf_sample = gsl_histogram_pdf_sample(p%gsl_histogram_pdf, r)
  end function fgsl_histogram_pdf_sample
  function fgsl_histogram2d_alloc(nx, ny)
    integer(fgsl_size_t), intent(in) :: nx, ny
    type(fgsl_histogram2d) :: fgsl_histogram2d_alloc
    fgsl_histogram2d_alloc%gsl_histogram2d = gsl_histogram2d_alloc(nx, ny)
  end function fgsl_histogram2d_alloc
  function fgsl_histogram2d_set_ranges(h, xrange, xsize, yrange, ysize) 
    type(fgsl_histogram2d), intent(inout) :: h
    integer(fgsl_size_t), intent(in) :: xsize, ysize
    real(fgsl_double), intent(in) :: xrange(:), yrange(:)
    integer(fgsl_int) :: fgsl_histogram2d_set_ranges
    fgsl_histogram2d_set_ranges = &
         gsl_histogram2d_set_ranges(h%gsl_histogram2d, xrange, xsize, yrange, ysize) 
  end function fgsl_histogram2d_set_ranges
  function fgsl_histogram2d_set_ranges_uniform(h, xmin, xmax, ymin, ymax)
    type(fgsl_histogram2d), intent(inout) :: h
    real(fgsl_double), intent(in) :: xmin, xmax, ymin, ymax
    integer(fgsl_int) :: fgsl_histogram2d_set_ranges_uniform
    fgsl_histogram2d_set_ranges_uniform = &
         gsl_histogram2d_set_ranges_uniform(h%gsl_histogram2d, xmin, xmax, ymin, ymax)
  end function fgsl_histogram2d_set_ranges_uniform
  subroutine fgsl_histogram2d_free(h) 
    type(fgsl_histogram2d), intent(inout) :: h
    call gsl_histogram2d_free(h%gsl_histogram2d) 
  end subroutine fgsl_histogram2d_free
  function fgsl_histogram2d_memcpy(dest, src) 
    type(fgsl_histogram2d), intent(inout) :: dest
    type(fgsl_histogram2d), intent(in) :: src
    integer(fgsl_int) :: fgsl_histogram2d_memcpy
    fgsl_histogram2d_memcpy = &
         gsl_histogram2d_memcpy(dest%gsl_histogram2d, src%gsl_histogram2d) 
  end function fgsl_histogram2d_memcpy
  function fgsl_histogram2d_clone(src) 
    type(fgsl_histogram2d), intent(in) :: src
    type(fgsl_histogram2d) :: fgsl_histogram2d_clone
    fgsl_histogram2d_clone%gsl_histogram2d = gsl_histogram2d_clone(src%gsl_histogram2d) 
  end function fgsl_histogram2d_clone
  function fgsl_histogram2d_increment(h, x, y) 
    type(fgsl_histogram2d), intent(inout) :: h
    real(fgsl_double), intent(in) :: x, y
    integer(fgsl_int) :: fgsl_histogram2d_increment
    fgsl_histogram2d_increment = gsl_histogram2d_increment(h%gsl_histogram2d, x, y) 
  end function fgsl_histogram2d_increment
  function fgsl_histogram2d_accumulate(h, x, y, weight) 
    type(fgsl_histogram2d), intent(inout) :: h
    real(fgsl_double), intent(in) :: x, y, weight
    integer(fgsl_int) :: fgsl_histogram2d_accumulate
    fgsl_histogram2d_accumulate = &
         gsl_histogram2d_accumulate(h%gsl_histogram2d, x, y, weight) 
  end function fgsl_histogram2d_accumulate
  function fgsl_histogram2d_get(h, i, j)
    type(fgsl_histogram2d), intent(in) :: h
    integer(fgsl_size_t), intent(in) :: i, j
    real(fgsl_double) :: fgsl_histogram2d_get
    fgsl_histogram2d_get = gsl_histogram2d_get(h%gsl_histogram2d, i, j)
  end function fgsl_histogram2d_get
  function fgsl_histogram2d_get_xrange(h, i, xlower, xupper) 
    type(fgsl_histogram2d), intent(in) :: h
    integer(fgsl_size_t), intent(in) :: i
    real(fgsl_double), intent(out) :: xlower, xupper
    integer(fgsl_int) :: fgsl_histogram2d_get_xrange
    fgsl_histogram2d_get_xrange = &
         gsl_histogram2d_get_xrange(h%gsl_histogram2d, i, xlower, xupper) 
  end function fgsl_histogram2d_get_xrange
  function fgsl_histogram2d_get_yrange(h, i, ylower, yupper) 
    type(fgsl_histogram2d), intent(in) :: h
    integer(fgsl_size_t), intent(in) :: i
    real(fgsl_double), intent(out) :: ylower, yupper
    integer(fgsl_int) :: fgsl_histogram2d_get_yrange
    fgsl_histogram2d_get_yrange = &
         gsl_histogram2d_get_yrange(h%gsl_histogram2d, i, ylower, yupper) 
  end function fgsl_histogram2d_get_yrange
  function fgsl_histogram2d_xmax(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_xmax
    fgsl_histogram2d_xmax = gsl_histogram2d_xmax(h%gsl_histogram2d)
  end function fgsl_histogram2d_xmax
  function fgsl_histogram2d_xmin(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_xmin
    fgsl_histogram2d_xmin = gsl_histogram2d_xmin(h%gsl_histogram2d)
  end function fgsl_histogram2d_xmin
  function fgsl_histogram2d_nx(h)
    type(fgsl_histogram2d), intent(in) :: h
    integer(fgsl_size_t) :: fgsl_histogram2d_nx
    fgsl_histogram2d_nx = gsl_histogram2d_nx(h%gsl_histogram2d)
  end function fgsl_histogram2d_nx
  function fgsl_histogram2d_ymax(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_ymax
    fgsl_histogram2d_ymax = gsl_histogram2d_ymax(h%gsl_histogram2d)
  end function fgsl_histogram2d_ymax
  function fgsl_histogram2d_ymin(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_ymin
    fgsl_histogram2d_ymin = gsl_histogram2d_ymin(h%gsl_histogram2d)
  end function fgsl_histogram2d_ymin
  function fgsl_histogram2d_ny(h)
    type(fgsl_histogram2d), intent(in) :: h
    integer(fgsl_size_t) :: fgsl_histogram2d_ny
    fgsl_histogram2d_ny = gsl_histogram2d_ny(h%gsl_histogram2d)
  end function fgsl_histogram2d_ny
  subroutine fgsl_histogram2d_reset(h) 
    type(fgsl_histogram2d), intent(inout) :: h
    call gsl_histogram2d_reset(h%gsl_histogram2d) 
  end subroutine fgsl_histogram2d_reset
  function fgsl_histogram2d_find(h, x, y, i, j) 
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double), intent(in) :: x, y    
    integer(fgsl_size_t), intent(out) :: i, j
    integer(fgsl_int) :: fgsl_histogram2d_find
    fgsl_histogram2d_find = gsl_histogram2d_find(h%gsl_histogram2d, x, y, i, j) 
  end function fgsl_histogram2d_find
  function fgsl_histogram2d_max_val(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_max_val
    fgsl_histogram2d_max_val = gsl_histogram2d_max_val(h%gsl_histogram2d)
  end function fgsl_histogram2d_max_val
  subroutine fgsl_histogram2d_max_bin(h, i, j)
    type(fgsl_histogram2d), intent(in) :: h
    integer(fgsl_size_t), intent(out) :: i, j
    call gsl_histogram2d_max_bin(h%gsl_histogram2d, i, j)
  end subroutine fgsl_histogram2d_max_bin
  function fgsl_histogram2d_min_val(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_min_val
    fgsl_histogram2d_min_val = gsl_histogram2d_min_val(h%gsl_histogram2d)
  end function fgsl_histogram2d_min_val
  subroutine fgsl_histogram2d_min_bin(h, i, j)
    type(fgsl_histogram2d), intent(in) :: h
    integer(fgsl_size_t), intent(out) :: i, j
    call gsl_histogram2d_min_bin(h%gsl_histogram2d, i, j)
  end subroutine fgsl_histogram2d_min_bin
  function fgsl_histogram2d_xmean(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_xmean
    fgsl_histogram2d_xmean = gsl_histogram2d_xmean(h%gsl_histogram2d)
  end function fgsl_histogram2d_xmean
  function fgsl_histogram2d_ymean(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_ymean
    fgsl_histogram2d_ymean = gsl_histogram2d_ymean(h%gsl_histogram2d)
  end function fgsl_histogram2d_ymean
  function fgsl_histogram2d_xsigma(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_xsigma
    fgsl_histogram2d_xsigma = gsl_histogram2d_xsigma(h%gsl_histogram2d)
  end function fgsl_histogram2d_xsigma
  function fgsl_histogram2d_ysigma(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_ysigma
    fgsl_histogram2d_ysigma = gsl_histogram2d_ysigma(h%gsl_histogram2d)
  end function fgsl_histogram2d_ysigma
  function fgsl_histogram2d_cov(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_cov
    fgsl_histogram2d_cov = gsl_histogram2d_cov(h%gsl_histogram2d)
  end function fgsl_histogram2d_cov
  function fgsl_histogram2d_sum(h)
    type(fgsl_histogram2d), intent(in) :: h
    real(fgsl_double) :: fgsl_histogram2d_sum
    fgsl_histogram2d_sum = gsl_histogram2d_sum(h%gsl_histogram2d)
  end function fgsl_histogram2d_sum
  function fgsl_histogram2d_equal_bins_p(h1, h2)
    type(fgsl_histogram2d), intent(in) :: h1, h2
    real(fgsl_double) :: fgsl_histogram2d_equal_bins_p
    fgsl_histogram2d_equal_bins_p = & 
         gsl_histogram2d_equal_bins_p(h1%gsl_histogram2d, h2%gsl_histogram2d)
  end function fgsl_histogram2d_equal_bins_p
  function fgsl_histogram2d_add(h1, h2)
    type(fgsl_histogram2d), intent(inout) :: h1
    type(fgsl_histogram2d), intent(in) :: h2
    real(fgsl_double) :: fgsl_histogram2d_add
    fgsl_histogram2d_add = & 
         gsl_histogram2d_add(h1%gsl_histogram2d, h2%gsl_histogram2d)
  end function fgsl_histogram2d_add
  function fgsl_histogram2d_sub(h1, h2)
    type(fgsl_histogram2d), intent(inout) :: h1
    type(fgsl_histogram2d), intent(in) :: h2
    real(fgsl_double) :: fgsl_histogram2d_sub
    fgsl_histogram2d_sub = & 
         gsl_histogram2d_sub(h1%gsl_histogram2d, h2%gsl_histogram2d)
  end function fgsl_histogram2d_sub
  function fgsl_histogram2d_mul(h1, h2)
    type(fgsl_histogram2d), intent(inout) :: h1
    type(fgsl_histogram2d), intent(in) :: h2
    real(fgsl_double) :: fgsl_histogram2d_mul
    fgsl_histogram2d_mul = & 
         gsl_histogram2d_mul(h1%gsl_histogram2d, h2%gsl_histogram2d)
  end function fgsl_histogram2d_mul
  function fgsl_histogram2d_div(h1, h2)
    type(fgsl_histogram2d), intent(inout) :: h1
    type(fgsl_histogram2d), intent(in) :: h2
    real(fgsl_double) :: fgsl_histogram2d_div
    fgsl_histogram2d_div = & 
         gsl_histogram2d_div(h1%gsl_histogram2d, h2%gsl_histogram2d)
  end function fgsl_histogram2d_div
  function fgsl_histogram2d_scale(h, scale) 
    type(fgsl_histogram2d), intent(inout) :: h
    real(fgsl_double), intent(in) :: scale
    integer(fgsl_int) :: fgsl_histogram2d_scale
    fgsl_histogram2d_scale = gsl_histogram2d_scale(h%gsl_histogram2d, scale) 
  end function fgsl_histogram2d_scale
  function fgsl_histogram2d_shift(h, offset) 
    type(fgsl_histogram2d), intent(inout) :: h
    real(fgsl_double), intent(in) :: offset
    integer(fgsl_int) :: fgsl_histogram2d_shift
    fgsl_histogram2d_shift = gsl_histogram2d_shift(h%gsl_histogram2d, offset) 
  end function fgsl_histogram2d_shift
  function fgsl_histogram2d_fwrite(stream, h)
    type(fgsl_file), intent(in) :: stream
    type(fgsl_histogram2d), intent(in) :: h
    integer(fgsl_int) :: fgsl_histogram2d_fwrite
    fgsl_histogram2d_fwrite = gsl_histogram2d_fwrite(stream%gsl_file, h%gsl_histogram2d)
  end function fgsl_histogram2d_fwrite
  function fgsl_histogram2d_fread(stream, h)
    type(fgsl_file), intent(in) :: stream
    type(fgsl_histogram2d), intent(inout) :: h
    integer(fgsl_int) :: fgsl_histogram2d_fread
    fgsl_histogram2d_fread = gsl_histogram2d_fread(stream%gsl_file, h%gsl_histogram2d)
  end function fgsl_histogram2d_fread
  function fgsl_histogram2d_fprintf(stream, h, range_format, bin_format) 
    type(fgsl_file), intent(in) :: stream
    type(fgsl_histogram2d), intent(in) :: h
    character(kind=fgsl_char, len=*), intent(in) :: range_format, bin_format
    integer(fgsl_int) :: fgsl_histogram2d_fprintf
!
    character(kind=fgsl_char,len=fgsl_strmax) :: lrf, lbf
    if (len(trim(range_format)) < fgsl_strmax .and. &
         len(trim(bin_format)) < fgsl_strmax) then
       lrf = trim(range_format) // c_null_char 
       lbf = trim(bin_format) // c_null_char 
       fgsl_histogram2d_fprintf = &
            gsl_histogram2d_fprintf(stream%gsl_file, h%gsl_histogram2d, lrf, lbf)
    else
       fgsl_histogram2d_fprintf = fgsl_failure
    end if
  end function fgsl_histogram2d_fprintf
  function fgsl_histogram2d_fscanf(stream, h) 
    type(fgsl_file), intent(in) :: stream
    type(fgsl_histogram2d), intent(inout) :: h
    integer(fgsl_int) :: fgsl_histogram2d_fscanf
    fgsl_histogram2d_fscanf = gsl_histogram2d_fscanf(stream%gsl_file, h%gsl_histogram2d)
  end function fgsl_histogram2d_fscanf
  function fgsl_histogram2d_pdf_alloc(nx, ny)
    integer(fgsl_size_t), intent(in) :: nx, ny
    type(fgsl_histogram2d_pdf) :: fgsl_histogram2d_pdf_alloc
    fgsl_histogram2d_pdf_alloc%gsl_histogram2d_pdf = gsl_histogram2d_pdf_alloc(nx, ny)
  end function fgsl_histogram2d_pdf_alloc
  function fgsl_histogram2d_pdf_init(p, h) 
    type(fgsl_histogram2d_pdf), intent(inout) :: p
    type(fgsl_histogram2d), intent(in) :: h
    integer(fgsl_int) :: fgsl_histogram2d_pdf_init
    fgsl_histogram2d_pdf_init = gsl_histogram2d_pdf_init(p%gsl_histogram2d_pdf, h%gsl_histogram2d) 
  end function fgsl_histogram2d_pdf_init
  subroutine fgsl_histogram2d_pdf_free(p) 
    type(fgsl_histogram2d_pdf), intent(inout) :: p
    call gsl_histogram2d_pdf_free(p%gsl_histogram2d_pdf) 
  end subroutine fgsl_histogram2d_pdf_free
  function fgsl_histogram2d_pdf_sample(p, r1, r2, x, y) 
    type(fgsl_histogram2d_pdf), intent(in) :: p
    real(fgsl_double), intent(in) :: r1, r2
    real(fgsl_double), intent(out) :: x, y
    integer(fgsl_int) :: fgsl_histogram2d_pdf_sample
    fgsl_histogram2d_pdf_sample = gsl_histogram2d_pdf_sample(p%gsl_histogram2d_pdf, r1, r2, x, y)
  end function fgsl_histogram2d_pdf_sample
  function fgsl_histogram_status(histogram)
    type(fgsl_histogram), intent(in) :: histogram
    logical :: fgsl_histogram_status
    fgsl_histogram_status = .true.
    if (.not. c_associated(histogram%gsl_histogram)) fgsl_histogram_status = .false.
  end function fgsl_histogram_status

!-*-f90-*-
!
! API: Ntuples
! 
  function fgsl_ntuple_create(fname, data, size) 
    character(kind=fgsl_char, len=*), intent(in) :: fname
    type(c_ptr), intent(in) :: data
    integer(fgsl_size_t), intent(in) :: size
    type(fgsl_ntuple) :: fgsl_ntuple_create
!
    character(kind=fgsl_char,len=fgsl_pathmax) :: lname
    if (len(trim(fname)) < fgsl_pathmax) then
       lname = trim(fname) // c_null_char
       fgsl_ntuple_create%gsl_ntuple = &
            gsl_ntuple_create(lname, data, size)
    else
       fgsl_ntuple_create%gsl_ntuple = c_null_ptr
    end if
  end function fgsl_ntuple_create
  function fgsl_ntuple_open(fname, data, size) 
    character(kind=fgsl_char, len=*), intent(in) :: fname
    type(c_ptr), intent(in) :: data
    integer(fgsl_size_t), intent(in) :: size
    type(fgsl_ntuple) :: fgsl_ntuple_open
!
    character(kind=fgsl_char,len=fgsl_pathmax) :: lname
    if (len(trim(fname)) < fgsl_pathmax) then
       lname = trim(fname) // c_null_char
       fgsl_ntuple_open%gsl_ntuple = &
            gsl_ntuple_open(lname, data, size)
    else
       fgsl_ntuple_open%gsl_ntuple = c_null_ptr
    end if
  end function fgsl_ntuple_open
  function fgsl_ntuple_write(ntuple) 
    type(fgsl_ntuple), intent(in) :: ntuple
    integer(fgsl_int) :: fgsl_ntuple_write
    fgsl_ntuple_write = gsl_ntuple_write(ntuple%gsl_ntuple)
  end function fgsl_ntuple_write
  function fgsl_ntuple_bookdata(ntuple) 
    type(fgsl_ntuple), intent(in) :: ntuple
    integer(fgsl_int) :: fgsl_ntuple_bookdata
    fgsl_ntuple_bookdata = gsl_ntuple_write(ntuple%gsl_ntuple)
  end function fgsl_ntuple_bookdata
  function fgsl_ntuple_read(ntuple) 
    type(fgsl_ntuple), intent(inout) :: ntuple
    integer(fgsl_int) :: fgsl_ntuple_read
    fgsl_ntuple_read = gsl_ntuple_read(ntuple%gsl_ntuple)
  end function fgsl_ntuple_read
  function fgsl_ntuple_close(ntuple) 
    type(fgsl_ntuple), intent(inout) :: ntuple
    integer(fgsl_int) :: fgsl_ntuple_close
    fgsl_ntuple_close = gsl_ntuple_close(ntuple%gsl_ntuple)
  end function fgsl_ntuple_close
  function fgsl_ntuple_select_fn_init(func, params)
    interface
       function func(data, params) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: data
         type(c_ptr), value :: params
         integer(c_int) :: func
       end function func
    end interface
    type(c_ptr), intent(in) :: params
    type(fgsl_ntuple_select_fn) :: fgsl_ntuple_select_fn_init
!
    type(c_funptr) :: fp
    fp = c_funloc(func)
    fgsl_ntuple_select_fn_init%gsl_ntuple_select_fn = &
         fgsl_ntuple_select_fn_cinit(fp, params)
  end function fgsl_ntuple_select_fn_init
  function fgsl_ntuple_value_fn_init(func, params)
    interface
       function func(data, params) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: data
         type(c_ptr), value :: params
         real(c_double) :: func
       end function func
    end interface
    type(c_ptr), intent(in) :: params
    type(fgsl_ntuple_value_fn) :: fgsl_ntuple_value_fn_init
!
    type(c_funptr) :: fp
    fp = c_funloc(func)
    fgsl_ntuple_value_fn_init%gsl_ntuple_value_fn = &
         fgsl_ntuple_value_fn_cinit(fp, params)
  end function fgsl_ntuple_value_fn_init
  subroutine fgsl_ntuple_select_fn_free(sfunc)
    type(fgsl_ntuple_select_fn), intent(inout) :: sfunc
    call fgsl_ntuple_select_fn_cfree(sfunc%gsl_ntuple_select_fn)
  end subroutine fgsl_ntuple_select_fn_free
  subroutine fgsl_ntuple_value_fn_free(sfunc)
    type(fgsl_ntuple_value_fn), intent(inout) :: sfunc
    call fgsl_ntuple_value_fn_cfree(sfunc%gsl_ntuple_value_fn)
  end subroutine fgsl_ntuple_value_fn_free
  function fgsl_ntuple_project(h, ntuple, value_func, select_func)
    type(fgsl_histogram), intent(inout) :: h
    type(fgsl_ntuple), intent(in) :: ntuple
    type(fgsl_ntuple_value_fn), intent(in) :: value_func
    type(fgsl_ntuple_select_fn), intent(in) :: select_func
    integer(fgsl_int) :: fgsl_ntuple_project
    fgsl_ntuple_project = gsl_ntuple_project(h%gsl_histogram, &
         ntuple%gsl_ntuple, value_func%gsl_ntuple_value_fn, &
         select_func%gsl_ntuple_select_fn)
  end function fgsl_ntuple_project
! Add-ons
  function fgsl_ntuple_data(ntuple)
    type(fgsl_ntuple), intent(in) :: ntuple
    type(c_ptr) :: fgsl_ntuple_data
    fgsl_ntuple_data = fgsl_aux_ntuple_data(ntuple%gsl_ntuple)
  end function fgsl_ntuple_data
  function fgsl_ntuple_size(ntuple)
    type(fgsl_ntuple), intent(in) :: ntuple
    integer(fgsl_size_t) :: fgsl_ntuple_size
    fgsl_ntuple_size = fgsl_aux_ntuple_size(ntuple%gsl_ntuple)
  end function fgsl_ntuple_size
  function fgsl_ntuple_status(ntuple)
    type(fgsl_ntuple), intent(in) :: ntuple
    logical :: fgsl_ntuple_status
    fgsl_ntuple_status = .true.
    if (.not. c_associated(ntuple%gsl_ntuple)) fgsl_ntuple_status = .false.
  end function fgsl_ntuple_status
  function fgsl_ntuple_value_fn_status(ntuple_value_fn)
    type(fgsl_ntuple_value_fn), intent(in) :: ntuple_value_fn
    logical :: fgsl_ntuple_value_fn_status
    fgsl_ntuple_value_fn_status = .true.
    if (.not. c_associated(ntuple_value_fn%gsl_ntuple_value_fn)) &
         fgsl_ntuple_value_fn_status = .false.
  end function fgsl_ntuple_value_fn_status
  function fgsl_ntuple_select_fn_status(ntuple_select_fn)
    type(fgsl_ntuple_select_fn), intent(in) :: ntuple_select_fn
    logical :: fgsl_ntuple_select_fn_status
    fgsl_ntuple_select_fn_status = .true.
    if (.not. c_associated(ntuple_select_fn%gsl_ntuple_select_fn)) &
         fgsl_ntuple_select_fn_status = .false.
  end function fgsl_ntuple_select_fn_status

!-*-f90-*-
!
! API: Monte Carlo integration
!
  function fgsl_monte_function_init(func, dim, params)
    interface
       function func(x, dim, params) bind(c)
         use, intrinsic :: iso_c_binding
         real(c_double), dimension(*) :: x
         integer(c_size_t), value :: dim
         type(c_ptr), value :: params
         real(c_double) :: func
       end function func
    end interface
    integer(fgsl_size_t), intent(in) :: dim
    type(c_ptr), intent(in) :: params
    type(fgsl_monte_function) :: fgsl_monte_function_init
!
    type(c_funptr) :: fp
    fp = c_funloc(func)
    fgsl_monte_function_init%gsl_monte_function = &
         fgsl_monte_function_cinit(fp, dim, params)
  end function fgsl_monte_function_init
  subroutine fgsl_monte_function_free(func)
    type(fgsl_monte_function), intent(inout) :: func
    call fgsl_monte_function_cfree(func%gsl_monte_function)
  end subroutine fgsl_monte_function_free
  function fgsl_monte_plain_alloc(dim) 
    integer(fgsl_size_t), intent(in) :: dim
    type(fgsl_monte_plain_state) :: fgsl_monte_plain_alloc
    fgsl_monte_plain_alloc%gsl_monte_plain_state = &
         gsl_monte_plain_alloc(dim)
  end function fgsl_monte_plain_alloc
  function fgsl_monte_plain_init(s) 
    type(fgsl_monte_plain_state), intent(in) :: s
    integer(fgsl_int) :: fgsl_monte_plain_init
    fgsl_monte_plain_init = gsl_monte_plain_init(s%gsl_monte_plain_state)
  end function fgsl_monte_plain_init
  function fgsl_monte_plain_integrate(f, xl, xu, dim, calls, r, s, result, abserr) 
    type(fgsl_monte_function), intent(in) :: f
    integer(fgsl_size_t), intent(in) :: dim, calls
    real(fgsl_double), intent(in) :: xl(dim), xu(dim)
    type(fgsl_rng), intent(in) :: r
    type(fgsl_monte_plain_state), intent(in) :: s
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_monte_plain_integrate
    fgsl_monte_plain_integrate = &
         gsl_monte_plain_integrate(f%gsl_monte_function, xl, xu, dim, calls, &
         r%gsl_rng, s%gsl_monte_plain_state, result, abserr)
  end function fgsl_monte_plain_integrate
  subroutine fgsl_monte_plain_free(s)
    type(fgsl_monte_plain_state), intent(inout) :: s
    call gsl_monte_plain_free(s%gsl_monte_plain_state)
  end subroutine fgsl_monte_plain_free
  function fgsl_monte_miser_alloc(dim)
    integer(fgsl_size_t), value :: dim
    type(fgsl_monte_miser_state) :: fgsl_monte_miser_alloc 
    fgsl_monte_miser_alloc%gsl_monte_miser_state = gsl_monte_miser_alloc(dim)
  end function fgsl_monte_miser_alloc
  function fgsl_monte_miser_init(s) 
    type(fgsl_monte_miser_state), intent(in) :: s
    integer(fgsl_int) :: fgsl_monte_miser_init
    fgsl_monte_miser_init = gsl_monte_miser_init(s%gsl_monte_miser_state)
  end function fgsl_monte_miser_init
  function fgsl_monte_miser_integrate(f, xl, xu, dim, calls, r, s, result, abserr) 
    type(fgsl_monte_function), intent(in) :: f
    integer(fgsl_size_t), intent(in) :: dim, calls
    real(fgsl_double), intent(in) :: xl(dim), xu(dim)
    type(fgsl_rng), intent(in) :: r
    type(fgsl_monte_miser_state), intent(in) :: s
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_monte_miser_integrate
    fgsl_monte_miser_integrate = &
         gsl_monte_miser_integrate(f%gsl_monte_function, xl, xu, dim, calls, &
         r%gsl_rng, s%gsl_monte_miser_state, result, abserr)
  end function fgsl_monte_miser_integrate
  subroutine fgsl_monte_miser_free(s)
    type(fgsl_monte_miser_state), intent(inout) :: s
    call gsl_monte_miser_free(s%gsl_monte_miser_state)
  end subroutine fgsl_monte_miser_free
!
  function fgsl_monte_vegas_alloc(dim)
    integer(fgsl_size_t), value :: dim
    type(fgsl_monte_vegas_state) :: fgsl_monte_vegas_alloc 
    fgsl_monte_vegas_alloc%gsl_monte_vegas_state = gsl_monte_vegas_alloc(dim)
  end function fgsl_monte_vegas_alloc
  function fgsl_monte_vegas_init(s) 
    type(fgsl_monte_vegas_state), intent(in) :: s
    integer(fgsl_int) :: fgsl_monte_vegas_init
    fgsl_monte_vegas_init = gsl_monte_vegas_init(s%gsl_monte_vegas_state)
  end function fgsl_monte_vegas_init
  function fgsl_monte_vegas_integrate(f, xl, xu, dim, calls, r, s, result, abserr) 
    type(fgsl_monte_function), intent(in) :: f
    integer(fgsl_size_t), intent(in) :: dim, calls
    real(fgsl_double), intent(in) :: xl(dim), xu(dim)
    type(fgsl_rng), intent(in) :: r
    type(fgsl_monte_vegas_state), intent(in) :: s
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_monte_vegas_integrate
    fgsl_monte_vegas_integrate = &
         gsl_monte_vegas_integrate(f%gsl_monte_function, xl, xu, dim, calls, &
         r%gsl_rng, s%gsl_monte_vegas_state, result, abserr)
  end function fgsl_monte_vegas_integrate
  subroutine fgsl_monte_vegas_free(s)
    type(fgsl_monte_vegas_state), intent(inout) :: s
    call gsl_monte_vegas_free(s%gsl_monte_vegas_state)
  end subroutine fgsl_monte_vegas_free
  function fgsl_monte_vegas_chisq(s) 
    real(fgsl_double) :: fgsl_monte_vegas_chisq
    type(fgsl_monte_vegas_state), intent(in) :: s
    fgsl_monte_vegas_chisq = gsl_monte_vegas_chisq(s%gsl_monte_vegas_state)
  end function fgsl_monte_vegas_chisq
  subroutine fgsl_monte_vegas_runval(s, result, sigma)
    type(fgsl_monte_vegas_state), intent(in) :: s
    real(fgsl_double), intent(out) :: result, sigma
    call gsl_monte_vegas_runval(s%gsl_monte_vegas_state, result, sigma)
  end subroutine fgsl_monte_vegas_runval
  function fgsl_monte_function_status(monte_function)
    type(fgsl_monte_function), intent(in) :: monte_function
    logical :: fgsl_monte_function_status
    fgsl_monte_function_status = .true.
    if (.not. c_associated(monte_function%gsl_monte_function)) &
         fgsl_monte_function_status = .false.
  end function fgsl_monte_function_status
  function fgsl_monte_plain_status(monte_plain)
    type(fgsl_monte_plain_state), intent(in) :: monte_plain
    logical :: fgsl_monte_plain_status
    fgsl_monte_plain_status = .true.
    if (.not. c_associated(monte_plain%gsl_monte_plain_state)) &
         fgsl_monte_plain_status = .false.
  end function fgsl_monte_plain_status
  function fgsl_monte_miser_status(monte_miser)
    type(fgsl_monte_miser_state), intent(in) :: monte_miser
    logical :: fgsl_monte_miser_status
    fgsl_monte_miser_status = .true.
    if (.not. c_associated(monte_miser%gsl_monte_miser_state)) &
         fgsl_monte_miser_status = .false.
  end function fgsl_monte_miser_status
  function fgsl_monte_vegas_status(monte_vegas)
    type(fgsl_monte_vegas_state), intent(in) :: monte_vegas
    logical :: fgsl_monte_vegas_status
    fgsl_monte_vegas_status = .true.
    if (.not. c_associated(monte_vegas%gsl_monte_vegas_state)) &
         fgsl_monte_vegas_status = .false.
  end function fgsl_monte_vegas_status
! add-on routines
! NOTE: in GSL 1.13, accessors were also added to GSL. They're 
! slightly different named and - by necessity - have a differing interface
! For this reason, the functions below are retained.
  subroutine fgsl_monte_miser_setparams(s, estimate_frac, min_calls, &
       min_calls_per_bisection, alpha, dither) 
    type(fgsl_monte_miser_state), intent(inout) :: s
    real(fgsl_double), intent(in) :: estimate_frac, alpha, dither
    integer(fgsl_size_t), intent(in) ::  min_calls, min_calls_per_bisection
    call fgsl_monte_miser_csetparams(s%gsl_monte_miser_state, estimate_frac, &
         min_calls, min_calls_per_bisection, alpha, dither)
  end subroutine fgsl_monte_miser_setparams
  subroutine fgsl_monte_miser_getparams(s, estimate_frac, min_calls, &
       min_calls_per_bisection, alpha, dither) 
    type(fgsl_monte_miser_state), intent(in) :: s
    real(fgsl_double), intent(out) :: estimate_frac, alpha, dither
    integer(fgsl_size_t), intent(out) ::  min_calls, min_calls_per_bisection
    call fgsl_monte_miser_cgetparams(s%gsl_monte_miser_state, estimate_frac, &
         min_calls, min_calls_per_bisection, alpha, dither)
  end subroutine fgsl_monte_miser_getparams
  subroutine fgsl_monte_vegas_setparams(s, result, sigma, chisq, alpha, &
       iterations, stage, mode, verbose, ostream)
    type(fgsl_monte_vegas_state), intent(inout) :: s
    real(fgsl_double), intent(in) :: result, sigma, chisq, alpha
    integer(fgsl_size_t), intent(in) ::  iterations
    integer(fgsl_int), intent(in) :: stage, mode, verbose
    type(fgsl_file), intent(in) :: ostream
    call fgsl_monte_vegas_csetparams(s%gsl_monte_vegas_state, result, &
         sigma, chisq, alpha, iterations, stage, mode, verbose, &
         ostream%gsl_file)
  end subroutine fgsl_monte_vegas_setparams
  subroutine fgsl_monte_vegas_getparams(s, result, sigma, chisq, alpha, &
       iterations, stage, mode, verbose, ostream)
    type(fgsl_monte_vegas_state), intent(in) :: s
    real(fgsl_double), intent(out) :: result, sigma, chisq, alpha
    integer(fgsl_size_t), intent(out) ::  iterations
    integer(fgsl_int), intent(out) :: stage, mode, verbose
    type(fgsl_file), intent(out) :: ostream
    call fgsl_monte_vegas_cgetparams(s%gsl_monte_vegas_state, result, &
         sigma, chisq, alpha, iterations, stage, mode, verbose, &
         ostream%gsl_file)
  end subroutine fgsl_monte_vegas_getparams
  
!-*-f90-*-
!
! API: Simulated annealing
!
  subroutine fgsl_siman_params_init(params, n_tries, iters_fixed_t, &
       step_size, k, t_initial, mu_t, t_min)
    type(fgsl_siman_params_t), intent(inout) :: params
    integer(fgsl_int) :: n_tries, iters_fixed_t
    real(fgsl_double) :: step_size, k, t_initial, mu_t, t_min
    if (.not. associated(params%gsl_siman_params_t)) then
       allocate(params%gsl_siman_params_t)
    end if
    params%gsl_siman_params_t%n_tries = n_tries
    params%gsl_siman_params_t%iters_fixed_t = iters_fixed_t
    params%gsl_siman_params_t%step_size = step_size
    params%gsl_siman_params_t%k = k
    params%gsl_siman_params_t%t_initial = t_initial
    params%gsl_siman_params_t%mu_t = mu_t
    params%gsl_siman_params_t%t_min = t_min
  end subroutine fgsl_siman_params_init
  subroutine fgsl_siman_params_free(params)
    type(fgsl_siman_params_t), intent(inout) :: params
    if (associated(params%gsl_siman_params_t)) then
       deallocate(params%gsl_siman_params_t)
    end if
    params%gsl_siman_params_t => null()
  end subroutine fgsl_siman_params_free
  subroutine fgsl_siman_solve(rng, x0_p, ef, take_step, distance, &
       print_position, copy_func, copy_constructor, destructor, &
       element_size, params) 
    optional :: print_position, copy_func, copy_constructor, destructor, &
         element_size
    type(fgsl_rng), intent(in) :: rng
    type(c_ptr), intent(inout) :: x0_p
    interface
       function ef(xp) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: xp
         real(c_double) :: ef
       end function ef
       subroutine take_step(rng, xp, step_size) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: rng, xp
         real(c_double), value :: step_size
       end subroutine take_step
       function distance(xp, yp) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: xp, yp
         real(c_double) :: distance
       end function distance
       subroutine print_position(xp) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: xp
       end subroutine print_position
       subroutine copy_func(src, dest) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: src, dest
       end subroutine copy_func
       function copy_constructor(xp) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: xp
         type(c_ptr) :: copy_constructor
       end function copy_constructor
       subroutine destructor(xp) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: xp
       end subroutine destructor
    end interface
    integer(fgsl_size_t), intent(in) :: element_size
    type(fgsl_siman_params_t), intent(in) :: params
!
    type(c_funptr) :: p_ef, p_take_step, p_distance, p_print_position, &
         p_copy_func, p_copy_constructor, p_destructor
    p_ef = c_funloc(ef)
    p_take_step = c_funloc(take_step)
    p_distance = c_funloc(distance)
    if (present(print_position)) then
       p_print_position = c_funloc(print_position)
    else
       p_print_position = c_null_funptr
    end if
    if (present(copy_func) .and. present(copy_constructor) &
         .and. present(destructor)) then
       p_copy_func = c_funloc(copy_func)
       p_copy_constructor = c_funloc(copy_constructor)
       p_destructor = c_funloc(destructor)
       if (element_size /= 0) return
    else
       p_copy_func = c_null_funptr
       p_copy_constructor = c_null_funptr
       p_destructor = c_null_funptr
    end if
    if (associated(params%gsl_siman_params_t)) then
       call gsl_siman_solve(rng%gsl_rng, x0_p, p_ef, p_take_step, p_distance, &
            p_print_position, p_copy_func, p_copy_constructor, p_destructor, &
            element_size, params%gsl_siman_params_t)
    end if
  end subroutine fgsl_siman_solve
  function fgsl_siman_params_t_status(siman_params_t)
    type(fgsl_siman_params_t), intent(in) :: siman_params_t
    logical :: fgsl_siman_params_t_status
    fgsl_siman_params_t_status = .true.
    if (.not. associated(siman_params_t%gsl_siman_params_t)) &
         fgsl_siman_params_t_status = .false.
  end function fgsl_siman_params_t_status

!-*-f90-*-
!
! API: Ordinary differential equations
!
  function fgsl_odeiv_system_init(func, dimension, params, jacobian)
    optional :: jacobian
    interface
       function func(t, y, dydt, params) bind(c)
         use, intrinsic :: iso_c_binding
         real(c_double), value :: t
         real(c_double), dimension(*), intent(in) :: y
         real(c_double), dimension(*) :: dydt
         type(c_ptr), value :: params
         integer(c_int) :: func
       end function func
       function jacobian(t, y, dfdy, dfdt, params) bind(c)
         use, intrinsic :: iso_c_binding
         real(c_double), value :: t
         real(c_double), dimension(*), intent(in) :: y
         real(c_double), dimension(*) :: dfdy
         real(c_double), dimension(*) :: dfdt
         type(c_ptr), value :: params
         integer(c_int) :: jacobian
       end function jacobian
    end interface
    integer(fgsl_size_t) :: dimension
    type(c_ptr), intent(in), optional :: params
    type(fgsl_odeiv_system) :: fgsl_odeiv_system_init
!
    type(c_funptr) :: func_loc 
    type(c_funptr) :: jacobian_loc
    type(c_ptr) :: params_loc
! debug
!    integer(c_int) :: status
!    real(c_double) :: y(2), dydt(2)
!    write(6, *) 'Starting init with dimension: ',dimension
!    y = (/1.0_c_double, 0.0_c_double /)
!    status = func(0.0_c_double,y,dydt,params)
!    write(6, *) 'Function output: ',dydt(1:2)
! end debug
    func_loc = c_funloc(func)
    params_loc = c_null_ptr
    jacobian_loc = c_null_funptr
    if (present(jacobian)) jacobian_loc = c_funloc(jacobian)
    if (present(params)) params_loc = params
    fgsl_odeiv_system_init%gsl_odeiv_system = &
         fgsl_odeiv_system_cinit(func_loc, dimension, &
         params_loc, jacobian_loc)
  end function fgsl_odeiv_system_init
  subroutine fgsl_odeiv_system_free(system)
    type(fgsl_odeiv_system), intent(inout) :: system
    call fgsl_odeiv_system_cfree(system%gsl_odeiv_system)
  end subroutine fgsl_odeiv_system_free
  function fgsl_odeiv_step_alloc(t, dim) 
    type(fgsl_odeiv_step_type), intent(in) :: t
    integer(fgsl_size_t), intent(in) :: dim
    type(fgsl_odeiv_step) :: fgsl_odeiv_step_alloc
!
    type(c_ptr) :: step_type
    step_type = fgsl_aux_odeiv_step_alloc(t%which)
    if (c_associated(step_type)) then
       fgsl_odeiv_step_alloc%gsl_odeiv_step = gsl_odeiv_step_alloc(step_type, dim)
    else
       fgsl_odeiv_step_alloc%gsl_odeiv_step = c_null_ptr
    end if
  end function fgsl_odeiv_step_alloc
  function fgsl_odeiv_step_reset(s)
    type(fgsl_odeiv_step), intent(inout) :: s
    integer(fgsl_int) :: fgsl_odeiv_step_reset
    fgsl_odeiv_step_reset = gsl_odeiv_step_reset(s%gsl_odeiv_step)
  end function fgsl_odeiv_step_reset
  subroutine fgsl_odeiv_step_free(s)
    type(fgsl_odeiv_step), intent(inout) :: s
    call gsl_odeiv_step_free(s%gsl_odeiv_step)
  end subroutine fgsl_odeiv_step_free
  function fgsl_odeiv_step_name (s)
    type(fgsl_odeiv_step), intent(in) :: s
    character(kind=fgsl_char, len=fgsl_strmax) :: fgsl_odeiv_step_name
!
    type(c_ptr) :: name
    name = gsl_odeiv_step_name(s%gsl_odeiv_step)
    fgsl_odeiv_step_name = fgsl_name(name)
  end function fgsl_odeiv_step_name
  function fgsl_odeiv_step_order(s)
    type(fgsl_odeiv_step), intent(in) :: s
    integer(fgsl_int) :: fgsl_odeiv_step_order
    fgsl_odeiv_step_order = gsl_odeiv_step_order(s%gsl_odeiv_step)
  end function fgsl_odeiv_step_order
  function fgsl_odeiv_step_apply(s, t, h, y, yerr, dydt_in, dydt_out, dydt)
    type(fgsl_odeiv_step), intent(in) :: s
    real(fgsl_double), intent(in) :: t, h
    real(fgsl_double), intent(inout) :: y(:), yerr(:), dydt_in(:), dydt_out(:)
    type(fgsl_odeiv_system), intent(in) :: dydt
    integer(fgsl_int) :: fgsl_odeiv_step_apply
    fgsl_odeiv_step_apply = gsl_odeiv_step_apply(s%gsl_odeiv_step, t, h, y, yerr, &
         dydt_in, dydt_out, dydt%gsl_odeiv_system)
  end function fgsl_odeiv_step_apply
  function fgsl_odeiv_control_standard_new(eps_abs, eps_rel, a_y, a_dydt)
    real(fgsl_double), intent(in) :: eps_abs, eps_rel, a_y, a_dydt
    type(fgsl_odeiv_control) :: fgsl_odeiv_control_standard_new
    fgsl_odeiv_control_standard_new%gsl_odeiv_control = &
         gsl_odeiv_control_standard_new(eps_abs, eps_rel, a_y, a_dydt)
  end function fgsl_odeiv_control_standard_new
  function fgsl_odeiv_control_y_new(eps_abs, eps_rel)
    real(fgsl_double), intent(in) :: eps_abs, eps_rel
    type(fgsl_odeiv_control) :: fgsl_odeiv_control_y_new
    fgsl_odeiv_control_y_new%gsl_odeiv_control = &
         gsl_odeiv_control_y_new(eps_abs, eps_rel)
  end function fgsl_odeiv_control_y_new
  function fgsl_odeiv_control_yp_new(eps_abs, eps_rel)
    real(fgsl_double), intent(in) :: eps_abs, eps_rel
    type(fgsl_odeiv_control) :: fgsl_odeiv_control_yp_new
    fgsl_odeiv_control_yp_new%gsl_odeiv_control = &
         gsl_odeiv_control_yp_new(eps_abs, eps_rel)
  end function fgsl_odeiv_control_yp_new
  function fgsl_odeiv_control_scaled_new(eps_abs, eps_rel, a_y, a_dydt, scale_abs, dim)
    real(fgsl_double), intent(in) :: eps_abs, eps_rel, a_y, a_dydt
    real(fgsl_double), intent(in) :: scale_abs(:)
    integer(fgsl_size_t), intent(in) :: dim
    type(fgsl_odeiv_control) :: fgsl_odeiv_control_scaled_new
    fgsl_odeiv_control_scaled_new%gsl_odeiv_control = &
         gsl_odeiv_control_scaled_new(eps_abs, eps_rel, a_y, a_dydt, scale_abs, dim)
  end function fgsl_odeiv_control_scaled_new
! FIXME (?) fgsl_odeiv_control_alloc is presently not implemented
  function fgsl_odeiv_control_init(c, eps_abs, eps_rel, a_y, a_dydt)
    type(fgsl_odeiv_control), intent(in) :: c
    real(fgsl_double), intent(in) :: eps_abs, eps_rel, a_y, a_dydt
    integer(fgsl_int) :: fgsl_odeiv_control_init
    fgsl_odeiv_control_init = &
         gsl_odeiv_control_init(c%gsl_odeiv_control, eps_abs, eps_rel, a_y, a_dydt)
  end function fgsl_odeiv_control_init
  subroutine fgsl_odeiv_control_free(c)
    type(fgsl_odeiv_control), intent(inout) :: c
    call gsl_odeiv_control_free(c%gsl_odeiv_control)
  end subroutine fgsl_odeiv_control_free
  function fgsl_odeiv_control_hadjust(c, s, y0, yerr, dydt, h) 
    type(fgsl_odeiv_control), intent(in) :: c
    type(fgsl_odeiv_step), intent(in) :: s
    real(fgsl_double), intent(in) :: y0(:), yerr(:), dydt(:)
    real(fgsl_double), intent(inout) :: h(:)
    integer(fgsl_int) :: fgsl_odeiv_control_hadjust
    fgsl_odeiv_control_hadjust = gsl_odeiv_control_hadjust(c%gsl_odeiv_control, s%gsl_odeiv_step, &
         y0, yerr, dydt, h) 
  end function fgsl_odeiv_control_hadjust
  function fgsl_odeiv_control_name (s)
    type(fgsl_odeiv_step), intent(in) :: s
    character(kind=fgsl_char, len=fgsl_strmax) :: fgsl_odeiv_control_name
!
    type(c_ptr) :: name
    name = gsl_odeiv_control_name(s%gsl_odeiv_step)
    fgsl_odeiv_control_name = fgsl_name(name)
  end function fgsl_odeiv_control_name
  function fgsl_odeiv_evolve_alloc(dim)
    integer(fgsl_size_t), intent(in) :: dim
    type(fgsl_odeiv_evolve) :: fgsl_odeiv_evolve_alloc
    fgsl_odeiv_evolve_alloc%gsl_odeiv_evolve = &
         gsl_odeiv_evolve_alloc(dim)
  end function fgsl_odeiv_evolve_alloc
  function fgsl_odeiv_evolve_apply(e, con, step, dydt, t, t1, h, y) 
    type(fgsl_odeiv_evolve), intent(inout) :: e
    type(fgsl_odeiv_control), intent(inout) :: con
    type(fgsl_odeiv_step), intent(inout) :: step
    type(fgsl_odeiv_system), intent(in) :: dydt
    real(fgsl_double), intent(inout) :: t, h, y(:)
    real(fgsl_double), intent(in) :: t1
    integer(fgsl_int) :: fgsl_odeiv_evolve_apply
!    write(6, *) 'Start evolving'
    fgsl_odeiv_evolve_apply = gsl_odeiv_evolve_apply(e%gsl_odeiv_evolve, &
         con%gsl_odeiv_control, step%gsl_odeiv_step, dydt%gsl_odeiv_system, &
         t, t1, h, y)
  end function fgsl_odeiv_evolve_apply
  function fgsl_odeiv_evolve_reset(s)
    type(fgsl_odeiv_evolve), intent(inout) :: s
    integer(c_int) :: fgsl_odeiv_evolve_reset
    fgsl_odeiv_evolve_reset = gsl_odeiv_evolve_reset(s%gsl_odeiv_evolve)
  end function fgsl_odeiv_evolve_reset
  subroutine fgsl_odeiv_evolve_free(s)
    type(fgsl_odeiv_evolve), intent(inout) :: s
    call gsl_odeiv_evolve_free(s%gsl_odeiv_evolve)
  end subroutine fgsl_odeiv_evolve_free
  function fgsl_odeiv_evolve_status(s)
    type(fgsl_odeiv_evolve), intent(in) :: s
    logical :: fgsl_odeiv_evolve_status
    fgsl_odeiv_evolve_status = .true.
    if (.not. c_associated(s%gsl_odeiv_evolve)) &
         fgsl_odeiv_evolve_status = .false.
  end function fgsl_odeiv_evolve_status
  function fgsl_odeiv_control_status(s)
    type(fgsl_odeiv_control), intent(in) :: s
    logical :: fgsl_odeiv_control_status
    fgsl_odeiv_control_status = .true.
    if (.not. c_associated(s%gsl_odeiv_control)) &
         fgsl_odeiv_control_status = .false.
  end function fgsl_odeiv_control_status
  function fgsl_odeiv_step_status(s)
    type(fgsl_odeiv_step), intent(in) :: s
    logical :: fgsl_odeiv_step_status
    fgsl_odeiv_step_status = .true.
    if (.not. c_associated(s%gsl_odeiv_step)) &
         fgsl_odeiv_step_status = .false.
  end function fgsl_odeiv_step_status
  function fgsl_odeiv_system_status(s)
    type(fgsl_odeiv_system), intent(in) :: s
    logical :: fgsl_odeiv_system_status
    fgsl_odeiv_system_status = .true.
    if (.not. c_associated(s%gsl_odeiv_system)) &
         fgsl_odeiv_system_status = .false.
  end function fgsl_odeiv_system_status
!-*-f90-*-
!
! API: Numerical Differentiation
!
  function fgsl_deriv_central(f, x, h, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: x, h
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_deriv_central
    fgsl_deriv_central = gsl_deriv_central(f%gsl_function, x, h, result, abserr) 
  end function fgsl_deriv_central
  function fgsl_deriv_forward(f, x, h, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: x, h
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_deriv_forward
    fgsl_deriv_forward = gsl_deriv_forward(f%gsl_function, x, h, result, abserr) 
  end function fgsl_deriv_forward
  function fgsl_deriv_backward(f, x, h, result, abserr) 
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: x, h
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_deriv_backward
    fgsl_deriv_backward = gsl_deriv_backward(f%gsl_function, x, h, result, abserr) 
  end function fgsl_deriv_backward
!-*-f90-*-
!
! API: Chebyshev Approximations
!
  function fgsl_cheb_alloc(n)
    integer(fgsl_int), intent(in) :: n
    type(fgsl_cheb_series) :: fgsl_cheb_alloc
    fgsl_cheb_alloc%gsl_cheb_series = gsl_cheb_alloc(n)
  end function fgsl_cheb_alloc
  subroutine fgsl_cheb_free(cs)
    type(fgsl_cheb_series), intent(in) :: cs
    call gsl_cheb_free(cs%gsl_cheb_series)
  end subroutine fgsl_cheb_free
  function fgsl_cheb_init(cs, f, a, b) 
    type(fgsl_cheb_series), intent(inout) :: cs
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: a, b
    integer(fgsl_int) :: fgsl_cheb_init
    fgsl_cheb_init = gsl_cheb_init(cs%gsl_cheb_series, f%gsl_function, a, b) 
  end function fgsl_cheb_init
  function fgsl_cheb_order(cs)
    integer(fgsl_size_t) :: fgsl_cheb_order
    type(fgsl_cheb_series), intent(in) :: cs
    fgsl_cheb_order = gsl_cheb_order(cs%gsl_cheb_series)
  end function fgsl_cheb_order
  function fgsl_cheb_size(cs)
    integer(fgsl_size_t) :: fgsl_cheb_size
    type(fgsl_cheb_series), intent(in) :: cs
    fgsl_cheb_size = gsl_cheb_size(cs%gsl_cheb_series)
  end function fgsl_cheb_size
  function fgsl_cheb_coeffs(cs)
    real(fgsl_double), pointer :: fgsl_cheb_coeffs(:)
    type(fgsl_cheb_series), intent(in) :: cs
    integer(fgsl_size_t) :: isz
    type(c_ptr) :: coeff_ptr
    isz = gsl_cheb_size(cs%gsl_cheb_series)
    coeff_ptr = gsl_cheb_coeffs(cs%gsl_cheb_series)
    if (c_associated(coeff_ptr)) then
       call c_f_pointer(coeff_ptr, fgsl_cheb_coeffs, (/isz/))
    else
       nullify(fgsl_cheb_coeffs)
    end if
  end function fgsl_cheb_coeffs
  function fgsl_cheb_eval(cs, x) 
    type(fgsl_cheb_series), intent(in) :: cs
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_cheb_eval
    fgsl_cheb_eval = gsl_cheb_eval(cs%gsl_cheb_series, x)
  end function fgsl_cheb_eval
  function fgsl_cheb_eval_err(cs, x, result, abserr) 
    type(fgsl_cheb_series), intent(in) :: cs
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_cheb_eval_err
    fgsl_cheb_eval_err = gsl_cheb_eval_err(cs%gsl_cheb_series, x, result, abserr) 
  end function fgsl_cheb_eval_err
  function fgsl_cheb_eval_n(cs, order, x) 
    type(fgsl_cheb_series), intent(in) :: cs
    integer(fgsl_size_t), intent(in) :: order
    real(fgsl_double), intent(in) :: x
    real(fgsl_double) :: fgsl_cheb_eval_n
    fgsl_cheb_eval_n = gsl_cheb_eval_n(cs%gsl_cheb_series, order, x)
  end function fgsl_cheb_eval_n
  function fgsl_cheb_eval_n_err(cs, order, x, result, abserr) 
    type(fgsl_cheb_series), intent(in) :: cs
    integer(fgsl_size_t), intent(in) :: order
    real(fgsl_double), intent(in) :: x
    real(fgsl_double), intent(out) :: result, abserr
    integer(fgsl_int) :: fgsl_cheb_eval_n_err
    fgsl_cheb_eval_n_err = gsl_cheb_eval_n_err(cs%gsl_cheb_series, order, x, result, abserr) 
  end function fgsl_cheb_eval_n_err
  function fgsl_cheb_calc_deriv(deriv, cs) 
    type(fgsl_cheb_series), intent(inout) :: deriv
    type(fgsl_cheb_series), intent(in) :: cs
    integer(fgsl_int) :: fgsl_cheb_calc_deriv
    fgsl_cheb_calc_deriv = gsl_cheb_calc_deriv(deriv%gsl_cheb_series, cs%gsl_cheb_series)
  end function fgsl_cheb_calc_deriv
  function fgsl_cheb_calc_integ(integ, cs) 
    type(fgsl_cheb_series), intent(inout) :: integ
    type(fgsl_cheb_series), intent(in) :: cs
    integer(fgsl_int) :: fgsl_cheb_calc_integ
    fgsl_cheb_calc_integ = gsl_cheb_calc_integ(integ%gsl_cheb_series, cs%gsl_cheb_series)
  end function fgsl_cheb_calc_integ
  function fgsl_cheb_series_status(cheb_series)
    type(fgsl_cheb_series), intent(in) :: cheb_series
    logical :: fgsl_cheb_series_status
    fgsl_cheb_series_status = .true.
    if (.not. c_associated(cheb_series%gsl_cheb_series)) &
         fgsl_cheb_series_status = .false.
  end function fgsl_cheb_series_status
!-*-f90-*-
!
! API: series acceleration
! 
  function fgsl_sum_levin_u_alloc(n) 
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_sum_levin_u_workspace) :: fgsl_sum_levin_u_alloc
    fgsl_sum_levin_u_alloc%gsl_sum_levin_u_workspace = &
         gsl_sum_levin_u_alloc(n) 
  end function fgsl_sum_levin_u_alloc
  function fgsl_sum_levin_u_free(w) 
    type(fgsl_sum_levin_u_workspace), intent(inout) :: w
    integer(fgsl_int) :: fgsl_sum_levin_u_free
    fgsl_sum_levin_u_free = &
         gsl_sum_levin_u_free(w%gsl_sum_levin_u_workspace)
  end function fgsl_sum_levin_u_free
  function fgsl_sum_levin_u_accel(array, array_size, w, sum_accel, abserr) 
    integer(fgsl_size_t), intent(in) :: array_size
    real(fgsl_double), intent(in) :: array(array_size)
    type(fgsl_sum_levin_u_workspace), intent(in) :: w
    real(fgsl_double), intent(out) :: sum_accel, abserr
    integer(fgsl_int) :: fgsl_sum_levin_u_accel
    fgsl_sum_levin_u_accel = gsl_sum_levin_u_accel(array, array_size, &
         w%gsl_sum_levin_u_workspace, sum_accel, abserr)
  end function fgsl_sum_levin_u_accel
  function fgsl_sum_levin_utrunc_alloc(n) 
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_sum_levin_utrunc_workspace) :: fgsl_sum_levin_utrunc_alloc
    fgsl_sum_levin_utrunc_alloc%gsl_sum_levin_utrunc_workspace = &
         gsl_sum_levin_utrunc_alloc(n) 
  end function fgsl_sum_levin_utrunc_alloc
  function fgsl_sum_levin_utrunc_free(w) 
    type(fgsl_sum_levin_utrunc_workspace), intent(inout) :: w
    integer(fgsl_int) :: fgsl_sum_levin_utrunc_free
    fgsl_sum_levin_utrunc_free = &
         gsl_sum_levin_utrunc_free(w%gsl_sum_levin_utrunc_workspace)
  end function fgsl_sum_levin_utrunc_free
  function fgsl_sum_levin_utrunc_accel(array, array_size, w, sum_accel, abserr) 
    integer(fgsl_size_t), intent(in) :: array_size
    real(fgsl_double), intent(in) :: array(array_size)
    type(fgsl_sum_levin_utrunc_workspace), intent(in) :: w
    real(fgsl_double), intent(out) :: sum_accel, abserr
    integer(fgsl_int) :: fgsl_sum_levin_utrunc_accel
    fgsl_sum_levin_utrunc_accel = gsl_sum_levin_utrunc_accel(array, array_size, &
         w%gsl_sum_levin_utrunc_workspace, sum_accel, abserr)
  end function fgsl_sum_levin_utrunc_accel
 
!-*-f90-*-
!
! API: wavelet transforms
! 
  function fgsl_wavelet_alloc(t, k)
    type(fgsl_wavelet_type), intent(in) :: t
    integer(fgsl_size_t), intent(in) :: k
    type(fgsl_wavelet) :: fgsl_wavelet_alloc
!
    type(c_ptr) :: cp
    integer(c_int) :: i
    i = t%which
    cp = fgsl_aux_wavelet_alloc(i)
    fgsl_wavelet_alloc%gsl_wavelet = gsl_wavelet_alloc(cp, k)
  end function fgsl_wavelet_alloc
  function fgsl_wavelet_name(wavelet)
    type(fgsl_wavelet), intent(in) :: wavelet
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_wavelet_name
!
    type(c_ptr) :: name
!
    name = gsl_wavelet_name(wavelet%gsl_wavelet)
    fgsl_wavelet_name = fgsl_name(name)
  end function fgsl_wavelet_name
  subroutine fgsl_wavelet_free(w) 
    type(fgsl_wavelet), intent(inout) :: w
    call gsl_wavelet_free(w%gsl_wavelet)
  end subroutine fgsl_wavelet_free
  function fgsl_wavelet_workspace_alloc(n) 
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_wavelet_workspace) :: fgsl_wavelet_workspace_alloc
    fgsl_wavelet_workspace_alloc%gsl_wavelet_workspace = &
         gsl_wavelet_workspace_alloc(n) 
  end function fgsl_wavelet_workspace_alloc
  subroutine fgsl_wavelet_workspace_free(w) 
    type(fgsl_wavelet_workspace), intent(inout) :: w
    call gsl_wavelet_workspace_free(w%gsl_wavelet_workspace)
  end subroutine fgsl_wavelet_workspace_free
  function fgsl_wavelet_transform(w, data, stride, n, dir, work) 
    type(fgsl_wavelet), intent(in) :: w
    real(fgsl_double), intent(inout) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    integer(fgsl_int), intent(in) :: dir
    type(fgsl_wavelet_workspace), intent(inout) :: work
    integer(fgsl_int) :: fgsl_wavelet_transform
    fgsl_wavelet_transform = gsl_wavelet_transform(w%gsl_wavelet, data, &
         stride, n, dir, work%gsl_wavelet_workspace) 
  end function fgsl_wavelet_transform
  function fgsl_wavelet_transform_forward(w, data, stride, n, work) 
    type(fgsl_wavelet), intent(in) :: w
    real(fgsl_double), intent(inout) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    type(fgsl_wavelet_workspace), intent(inout) :: work
    integer(fgsl_int) :: fgsl_wavelet_transform_forward
    fgsl_wavelet_transform_forward = gsl_wavelet_transform_forward(w%gsl_wavelet, data, &
         stride, n, work%gsl_wavelet_workspace) 
  end function fgsl_wavelet_transform_forward
  function fgsl_wavelet_transform_inverse(w, data, stride, n, work) 
    type(fgsl_wavelet), intent(in) :: w
    real(fgsl_double), intent(inout) :: data(:)
    integer(fgsl_size_t), intent(in) :: stride, n
    type(fgsl_wavelet_workspace), intent(inout) :: work
    integer(fgsl_int) :: fgsl_wavelet_transform_inverse
    fgsl_wavelet_transform_inverse = gsl_wavelet_transform_inverse(w%gsl_wavelet, data, &
         stride, n, work%gsl_wavelet_workspace) 
  end function fgsl_wavelet_transform_inverse
  function fgsl_wavelet2d_transform(w, data, tda, size1, size2, dir, work) 
    type(fgsl_wavelet), intent(in) :: w
    real(fgsl_double), intent(inout) :: data(:)
    integer(fgsl_size_t), intent(in) :: tda, size1, size2
    integer(fgsl_int), intent(in) :: dir
    type(fgsl_wavelet_workspace), intent(inout) :: work
    integer(fgsl_int) :: fgsl_wavelet2d_transform
    fgsl_wavelet2d_transform = gsl_wavelet2d_transform(w%gsl_wavelet, data, &
         tda, size1, size2, dir, work%gsl_wavelet_workspace) 
  end function fgsl_wavelet2d_transform
  function fgsl_wavelet2d_transform_forward(w, data, tda, size1, size2, work) 
    type(fgsl_wavelet), intent(in) :: w
    real(fgsl_double), intent(inout) :: data(:)
    integer(fgsl_size_t), intent(in) :: tda, size1, size2
    type(fgsl_wavelet_workspace), intent(inout) :: work
    integer(fgsl_int) :: fgsl_wavelet2d_transform_forward
    fgsl_wavelet2d_transform_forward = &
         gsl_wavelet2d_transform_forward(w%gsl_wavelet, data, &
         tda, size1, size2, work%gsl_wavelet_workspace) 
  end function fgsl_wavelet2d_transform_forward
  function fgsl_wavelet2d_transform_inverse(w, data, tda, size1, size2, work) 
    type(fgsl_wavelet), intent(in) :: w
    real(fgsl_double), intent(inout) :: data(:)
    integer(fgsl_size_t), intent(in) :: tda, size1, size2
    type(fgsl_wavelet_workspace), intent(inout) :: work
    integer(fgsl_int) :: fgsl_wavelet2d_transform_inverse
    fgsl_wavelet2d_transform_inverse = &
         gsl_wavelet2d_transform_inverse(w%gsl_wavelet, data, &
         tda, size1, size2, work%gsl_wavelet_workspace) 
  end function fgsl_wavelet2d_transform_inverse
  function fgsl_wavelet2d_transform_matrix(w, m, dir, work)
    type(fgsl_wavelet), intent(in) :: w
    type(fgsl_matrix), intent(inout) :: m
    type(fgsl_wavelet_workspace) :: work
    integer(fgsl_int), intent(in) :: dir
    integer(fgsl_int) :: fgsl_wavelet2d_transform_matrix
    fgsl_wavelet2d_transform_matrix = &
         gsl_wavelet2d_transform_matrix(w%gsl_wavelet, m%gsl_matrix, dir, &
         work%gsl_wavelet_workspace)
  end function fgsl_wavelet2d_transform_matrix
  function fgsl_wavelet2d_transform_matrix_forward(w, m, work)
    type(fgsl_wavelet), intent(in) :: w
    type(fgsl_matrix), intent(inout) :: m
    type(fgsl_wavelet_workspace) :: work
    integer(fgsl_int) :: fgsl_wavelet2d_transform_matrix_forward
    fgsl_wavelet2d_transform_matrix_forward = &
         gsl_wavelet2d_transform_matrix_forward(w%gsl_wavelet, m%gsl_matrix, &
         work%gsl_wavelet_workspace)
  end function fgsl_wavelet2d_transform_matrix_forward
  function fgsl_wavelet2d_transform_matrix_inverse(w, m, work)
    type(fgsl_wavelet), intent(in) :: w
    type(fgsl_matrix), intent(inout) :: m
    type(fgsl_wavelet_workspace) :: work
    integer(fgsl_int) :: fgsl_wavelet2d_transform_matrix_inverse
    fgsl_wavelet2d_transform_matrix_inverse = &
         gsl_wavelet2d_transform_matrix_inverse(w%gsl_wavelet, m%gsl_matrix, &
         work%gsl_wavelet_workspace)
  end function fgsl_wavelet2d_transform_matrix_inverse
  function fgsl_wavelet2d_nstransform(w, data, tda, size1, size2, dir, work) 
    type(fgsl_wavelet), intent(in) :: w
    real(fgsl_double), intent(inout) :: data(:)
    integer(fgsl_size_t), intent(in) :: tda, size1, size2
    integer(fgsl_int), intent(in) :: dir
    type(fgsl_wavelet_workspace), intent(inout) :: work
    integer(fgsl_int) :: fgsl_wavelet2d_nstransform
    fgsl_wavelet2d_nstransform = gsl_wavelet2d_nstransform(w%gsl_wavelet, data, &
         tda, size1, size2, dir, work%gsl_wavelet_workspace) 
  end function fgsl_wavelet2d_nstransform
  function fgsl_wavelet2d_nstransform_forward(w, data, tda, size1, size2, work) 
    type(fgsl_wavelet), intent(in) :: w
    real(fgsl_double), intent(inout) :: data(:)
    integer(fgsl_size_t), intent(in) :: tda, size1, size2
    type(fgsl_wavelet_workspace), intent(inout) :: work
    integer(fgsl_int) :: fgsl_wavelet2d_nstransform_forward
    fgsl_wavelet2d_nstransform_forward = &
         gsl_wavelet2d_nstransform_forward(w%gsl_wavelet, data, &
         tda, size1, size2, work%gsl_wavelet_workspace) 
  end function fgsl_wavelet2d_nstransform_forward
  function fgsl_wavelet2d_nstransform_inverse(w, data, tda, size1, size2, work) 
    type(fgsl_wavelet), intent(in) :: w
    real(fgsl_double), intent(inout) :: data(:)
    integer(fgsl_size_t), intent(in) :: tda, size1, size2
    type(fgsl_wavelet_workspace), intent(inout) :: work
    integer(fgsl_int) :: fgsl_wavelet2d_nstransform_inverse
    fgsl_wavelet2d_nstransform_inverse = &
         gsl_wavelet2d_nstransform_inverse(w%gsl_wavelet, data, &
         tda, size1, size2, work%gsl_wavelet_workspace) 
  end function fgsl_wavelet2d_nstransform_inverse
  function fgsl_wavelet2d_nstransform_matrix(w, m, dir, work)
    type(fgsl_wavelet), intent(in) :: w
    type(fgsl_matrix), intent(inout) :: m
    type(fgsl_wavelet_workspace) :: work
    integer(fgsl_int), intent(in) :: dir
    integer(fgsl_int) :: fgsl_wavelet2d_nstransform_matrix
    fgsl_wavelet2d_nstransform_matrix = &
         gsl_wavelet2d_nstransform_matrix(w%gsl_wavelet, m%gsl_matrix, dir, &
         work%gsl_wavelet_workspace)
  end function fgsl_wavelet2d_nstransform_matrix
  function fgsl_wavelet2d_nstransform_matrix_forward(w, m, work)
    type(fgsl_wavelet), intent(in) :: w
    type(fgsl_matrix), intent(inout) :: m
    type(fgsl_wavelet_workspace) :: work
    integer(fgsl_int) :: fgsl_wavelet2d_nstransform_matrix_forward
    fgsl_wavelet2d_nstransform_matrix_forward = &
         gsl_wavelet2d_nstransform_matrix_forward(w%gsl_wavelet, m%gsl_matrix, &
         work%gsl_wavelet_workspace)
  end function fgsl_wavelet2d_nstransform_matrix_forward
  function fgsl_wavelet2d_nstransform_matrix_inverse(w, m, work)
    type(fgsl_wavelet), intent(in) :: w
    type(fgsl_matrix), intent(inout) :: m
    type(fgsl_wavelet_workspace) :: work
    integer(fgsl_int) :: fgsl_wavelet2d_nstransform_matrix_inverse
    fgsl_wavelet2d_nstransform_matrix_inverse = &
         gsl_wavelet2d_nstransform_matrix_inverse(w%gsl_wavelet, m%gsl_matrix, &
         work%gsl_wavelet_workspace)
  end function fgsl_wavelet2d_nstransform_matrix_inverse
  function fgsl_wavelet_status(wavelet)
    type(fgsl_wavelet), intent(in) :: wavelet
    logical :: fgsl_wavelet_status
    fgsl_wavelet_status = .true.
    if (.not. c_associated(wavelet%gsl_wavelet)) fgsl_wavelet_status = .false.
  end function fgsl_wavelet_status
  function fgsl_wavelet_workspace_status(wavelet_workspace)
    type(fgsl_wavelet_workspace), intent(in) :: wavelet_workspace
    logical :: fgsl_wavelet_workspace_status
    fgsl_wavelet_workspace_status = .true.
    if (.not. c_associated(wavelet_workspace%gsl_wavelet_workspace)) &
         fgsl_wavelet_workspace_status = .false.
  end function fgsl_wavelet_workspace_status
  function fgsl_sizeof_wavelet(w)
    type(fgsl_wavelet), intent(in) :: w
    integer(fgsl_size_t) :: fgsl_sizeof_wavelet
    fgsl_sizeof_wavelet = gsl_aux_sizeof_wavelet()
  end function fgsl_sizeof_wavelet
  function fgsl_sizeof_wavelet_workspace(w)
    type(fgsl_wavelet_workspace), intent(in) :: w
    integer(fgsl_size_t) :: fgsl_sizeof_wavelet_workspace
    fgsl_sizeof_wavelet_workspace = gsl_aux_sizeof_wavelet_workspace()
  end function fgsl_sizeof_wavelet_workspace
 
 
 
 
!-*-f90-*-
!
! API: Hankel transforms
! 
  function fgsl_dht_alloc(size) 
    integer(fgsl_size_t), intent(in) :: size
    type(fgsl_dht) :: fgsl_dht_alloc
    fgsl_dht_alloc%gsl_dht = gsl_dht_alloc(size)
  end function fgsl_dht_alloc
  function fgsl_dht_init(t, nu, xmax) 
    type(fgsl_dht), intent(inout) :: t
    real(fgsl_double), intent(in) :: nu, xmax
    integer(fgsl_int) :: fgsl_dht_init
    fgsl_dht_init = gsl_dht_init(t%gsl_dht, nu, xmax) 
  end function fgsl_dht_init
  function fgsl_dht_new(size, nu, xmax) 
    integer(fgsl_size_t), intent(in) :: size
    real(fgsl_double), intent(in) :: nu, xmax
    type(fgsl_dht) :: fgsl_dht_new
    fgsl_dht_new%gsl_dht = gsl_dht_new(size, nu, xmax)
  end function fgsl_dht_new
  subroutine fgsl_dht_free(t)
    type(fgsl_dht), intent(inout) :: t
    call gsl_dht_free(t%gsl_dht)
  end subroutine fgsl_dht_free
  function fgsl_dht_apply(t, f_in, f_out)
    type(fgsl_dht), intent(in) :: t
    real(fgsl_double), intent(in) :: f_in(:)
    real(fgsl_double), intent(out) :: f_out(:)
    integer(fgsl_int) :: fgsl_dht_apply
    fgsl_dht_apply = gsl_dht_apply(t%gsl_dht, f_in, f_out)
  end function fgsl_dht_apply
  function fgsl_dht_x_sample(t, n)
    type(fgsl_dht), intent(in) :: t
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double) :: fgsl_dht_x_sample
    fgsl_dht_x_sample = gsl_dht_x_sample(t%gsl_dht, n)
  end function fgsl_dht_x_sample
  function fgsl_dht_k_sample(t, n)
    type(fgsl_dht), intent(in) :: t
    integer(fgsl_int), intent(in) :: n
    real(fgsl_double) :: fgsl_dht_k_sample
    fgsl_dht_k_sample = gsl_dht_k_sample(t%gsl_dht, n)
  end function fgsl_dht_k_sample
  function fgsl_dht_status(dht)
    type(fgsl_dht), intent(in) :: dht
    logical :: fgsl_dht_status
    fgsl_dht_status = .true.
    if (.not. c_associated(dht%gsl_dht)) fgsl_dht_status = .false.
  end function fgsl_dht_status
!-*-f90-*-
!
! API: Root finding
!
  function fgsl_root_fsolver_alloc(t)
    type(fgsl_root_fsolver_type), intent(in) :: t
    type(fgsl_root_fsolver) :: fgsl_root_fsolver_alloc
    type(c_ptr) :: it
    it = fgsl_aux_fsolver_alloc(t%which)
    if (c_associated(it)) then
       fgsl_root_fsolver_alloc%gsl_root_fsolver = gsl_root_fsolver_alloc(it)
    else
       fgsl_root_fsolver_alloc%gsl_root_fsolver = c_null_ptr
    end if
  end function fgsl_root_fsolver_alloc
  function fgsl_root_fdfsolver_alloc(t)
    type(fgsl_root_fdfsolver_type), intent(in) :: t
    type(fgsl_root_fdfsolver) :: fgsl_root_fdfsolver_alloc
    type(c_ptr) :: it
    it = fgsl_aux_fdfsolver_alloc(t%which)
    if (c_associated(it)) then
       fgsl_root_fdfsolver_alloc%gsl_root_fdfsolver = gsl_root_fdfsolver_alloc(it)
    else
       fgsl_root_fdfsolver_alloc%gsl_root_fdfsolver = c_null_ptr
    end if
  end function fgsl_root_fdfsolver_alloc
  function fgsl_root_fsolver_set(s, f, x_lower, x_upper) 
    type(fgsl_root_fsolver), intent(in) :: s
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: x_lower, x_upper
    integer(fgsl_int) :: fgsl_root_fsolver_set
    fgsl_root_fsolver_set = gsl_root_fsolver_set(s%gsl_root_fsolver, &
         f%gsl_function, x_lower, x_upper)
  end function fgsl_root_fsolver_set
  function fgsl_root_fdfsolver_set(s, fdf, x) 
    type(fgsl_root_fdfsolver), intent(in) :: s
    type(fgsl_function_fdf), intent(in) :: fdf
    real(fgsl_double), intent(in) :: x
    integer(fgsl_int) :: fgsl_root_fdfsolver_set
    fgsl_root_fdfsolver_set = gsl_root_fdfsolver_set(s%gsl_root_fdfsolver, &
         fdf%gsl_function_fdf, x)
  end function fgsl_root_fdfsolver_set
  subroutine fgsl_root_fsolver_free(s) 
    type(fgsl_root_fsolver), intent(inout) :: s
    call gsl_root_fsolver_free(s%gsl_root_fsolver)
  end subroutine fgsl_root_fsolver_free
  subroutine fgsl_root_fdfsolver_free(s) 
    type(fgsl_root_fdfsolver), intent(inout) :: s
    call gsl_root_fdfsolver_free(s%gsl_root_fdfsolver)
  end subroutine fgsl_root_fdfsolver_free
  function fgsl_root_fsolver_name(s) 
    type(fgsl_root_fsolver), intent(in) :: s
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_root_fsolver_name
!
    type(c_ptr) :: name
!
    name = gsl_root_fsolver_name(s%gsl_root_fsolver)
    fgsl_root_fsolver_name = fgsl_name(name)
  end function fgsl_root_fsolver_name
  function fgsl_root_fdfsolver_name(s) 
    type(fgsl_root_fdfsolver), intent(in) :: s
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_root_fdfsolver_name
!
    type(c_ptr) :: name
!
    name = gsl_root_fdfsolver_name(s%gsl_root_fdfsolver)
    fgsl_root_fdfsolver_name = fgsl_name(name)
  end function fgsl_root_fdfsolver_name
  function fgsl_root_fsolver_iterate(s) 
    type(fgsl_root_fsolver), intent(inout) :: s
    integer(fgsl_int) :: fgsl_root_fsolver_iterate
    fgsl_root_fsolver_iterate = gsl_root_fsolver_iterate(s%gsl_root_fsolver)
  end function fgsl_root_fsolver_iterate
  function fgsl_root_fdfsolver_iterate(s) 
    type(fgsl_root_fdfsolver), intent(inout) :: s
    integer(fgsl_int) :: fgsl_root_fdfsolver_iterate
    fgsl_root_fdfsolver_iterate = gsl_root_fdfsolver_iterate(s%gsl_root_fdfsolver)
  end function fgsl_root_fdfsolver_iterate
  function fgsl_root_fsolver_root(s) 
    type(fgsl_root_fsolver), intent(inout) :: s
    real(fgsl_double) :: fgsl_root_fsolver_root
    fgsl_root_fsolver_root = gsl_root_fsolver_root(s%gsl_root_fsolver)
  end function fgsl_root_fsolver_root
  function fgsl_root_fdfsolver_root(s) 
    type(fgsl_root_fdfsolver), intent(inout) :: s
    real(fgsl_double) :: fgsl_root_fdfsolver_root
    fgsl_root_fdfsolver_root = gsl_root_fdfsolver_root(s%gsl_root_fdfsolver)
  end function fgsl_root_fdfsolver_root
  function fgsl_root_fsolver_x_lower(s) 
    type(fgsl_root_fsolver), intent(inout) :: s
    real(fgsl_double) :: fgsl_root_fsolver_x_lower
    fgsl_root_fsolver_x_lower = gsl_root_fsolver_x_lower(s%gsl_root_fsolver)
  end function fgsl_root_fsolver_x_lower
  function fgsl_root_fsolver_x_upper(s) 
    type(fgsl_root_fsolver), intent(inout) :: s
    real(fgsl_double) :: fgsl_root_fsolver_x_upper
    fgsl_root_fsolver_x_upper = gsl_root_fsolver_x_upper(s%gsl_root_fsolver)
  end function fgsl_root_fsolver_x_upper
  function fgsl_root_test_interval(x_lower, x_upper, epsabs, epsrel) 
    real(fgsl_double), intent(in) :: x_lower, x_upper, epsabs, epsrel
    integer(fgsl_int) :: fgsl_root_test_interval
    fgsl_root_test_interval = gsl_root_test_interval(x_lower, x_upper, epsabs, epsrel)
  end function fgsl_root_test_interval
  function fgsl_root_test_delta(x1, x0, epsabs, epsrel) 
    real(fgsl_double), intent(in) :: x1, x0, epsabs, epsrel
    integer(fgsl_int) :: fgsl_root_test_delta
    fgsl_root_test_delta = gsl_root_test_delta(x1, x0, epsabs, epsrel)
  end function fgsl_root_test_delta
  function fgsl_root_test_residual(f, epsabs) 
    real(fgsl_double), intent(in) :: f, epsabs
    integer(fgsl_int) :: fgsl_root_test_residual
    fgsl_root_test_residual = gsl_root_test_residual(f, epsabs) 
  end function fgsl_root_test_residual
!
  function fgsl_root_fsolver_status(s)
    type(fgsl_root_fsolver), intent(in) :: s
    logical :: fgsl_root_fsolver_status
    fgsl_root_fsolver_status = .false.
    if (c_associated(s%gsl_root_fsolver)) &
         fgsl_root_fsolver_status = .true.
  end function fgsl_root_fsolver_status
  function fgsl_root_fdfsolver_status(s)
    type(fgsl_root_fdfsolver), intent(in) :: s
    logical :: fgsl_root_fdfsolver_status
    fgsl_root_fdfsolver_status = .false.
    if (c_associated(s%gsl_root_fdfsolver)) &
         fgsl_root_fdfsolver_status = .true.
  end function fgsl_root_fdfsolver_status
  
!-*-f90-*-
!
! API: Minimization
!
  function fgsl_min_fminimizer_alloc(t) 
    type(fgsl_min_fminimizer_type), intent(in) :: t
    type(fgsl_min_fminimizer) :: fgsl_min_fminimizer_alloc
    type(c_ptr) :: it
    it = fgsl_aux_fminimizer_alloc(t%which)
    if (c_associated(it)) then
       fgsl_min_fminimizer_alloc%gsl_min_fminimizer = gsl_min_fminimizer_alloc(it)
    else
       fgsl_min_fminimizer_alloc%gsl_min_fminimizer = c_null_ptr
    end if
  end function fgsl_min_fminimizer_alloc
  subroutine fgsl_min_fminimizer_free(s)
    type(fgsl_min_fminimizer), intent(inout) :: s
    call gsl_min_fminimizer_free(s%gsl_min_fminimizer)
  end subroutine fgsl_min_fminimizer_free
  function fgsl_min_fminimizer_set(s, f, x_minimum, x_lower, x_upper)
    type(fgsl_min_fminimizer), intent(inout) :: s
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), intent(in) :: x_minimum, x_lower, x_upper
    integer(fgsl_int) :: fgsl_min_fminimizer_set
    fgsl_min_fminimizer_set = gsl_min_fminimizer_set(s%gsl_min_fminimizer, f%gsl_function, &
                              x_minimum, x_lower, x_upper)
  end function fgsl_min_fminimizer_set
  function fgsl_min_fminimizer_set_with_values(s, f, x_minimum, f_minimum, &
           x_lower, f_lower, x_upper, f_upper) 
    type(fgsl_min_fminimizer), intent(inout) :: s
    type(fgsl_function), intent(in) :: f
    real(fgsl_double), value :: x_minimum, f_minimum, &
           x_lower, f_lower, x_upper, f_upper
    integer(fgsl_int) :: fgsl_min_fminimizer_set_with_values
    fgsl_min_fminimizer_set_with_values = gsl_min_fminimizer_set_with_values( &
           s%gsl_min_fminimizer, f%gsl_function, x_minimum, f_minimum, &
           x_lower, f_lower, x_upper, f_upper) 
  end function fgsl_min_fminimizer_set_with_values
  function fgsl_min_fminimizer_iterate(s) 
    type(fgsl_min_fminimizer), value :: s
    integer(fgsl_int) :: fgsl_min_fminimizer_iterate
    fgsl_min_fminimizer_iterate = gsl_min_fminimizer_iterate(s%gsl_min_fminimizer)
  end function fgsl_min_fminimizer_iterate
  function fgsl_min_fminimizer_name(s)
    type(fgsl_min_fminimizer), intent(in) :: s
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_min_fminimizer_name
!
    type(c_ptr) :: name
!
    name = gsl_min_fminimizer_name(s%gsl_min_fminimizer)
    fgsl_min_fminimizer_name = fgsl_name(name)
  end function fgsl_min_fminimizer_name
  function fgsl_min_fminimizer_x_minimum(s) 
    type(fgsl_min_fminimizer), intent(in) :: s
    real(fgsl_double) ::  fgsl_min_fminimizer_x_minimum
    fgsl_min_fminimizer_x_minimum = gsl_min_fminimizer_x_minimum(s%gsl_min_fminimizer)
  end function fgsl_min_fminimizer_x_minimum
  function fgsl_min_fminimizer_x_lower(s) 
    type(fgsl_min_fminimizer), intent(in) :: s
    real(fgsl_double) ::  fgsl_min_fminimizer_x_lower
    fgsl_min_fminimizer_x_lower = gsl_min_fminimizer_x_lower(s%gsl_min_fminimizer)
  end function fgsl_min_fminimizer_x_lower
  function fgsl_min_fminimizer_x_upper(s)
    type(fgsl_min_fminimizer), intent(in) :: s
    real(fgsl_double) ::  fgsl_min_fminimizer_x_upper
    fgsl_min_fminimizer_x_upper = gsl_min_fminimizer_x_upper(s%gsl_min_fminimizer)
  end function fgsl_min_fminimizer_x_upper
  function fgsl_min_fminimizer_f_minimum(s)
    type(fgsl_min_fminimizer), intent(in) :: s
    real(fgsl_double) ::  fgsl_min_fminimizer_f_minimum
    fgsl_min_fminimizer_f_minimum = gsl_min_fminimizer_f_minimum(s%gsl_min_fminimizer)
  end function fgsl_min_fminimizer_f_minimum
  function fgsl_min_fminimizer_f_lower(s) 
    type(fgsl_min_fminimizer), intent(in) :: s
    real(fgsl_double) ::  fgsl_min_fminimizer_f_lower
    fgsl_min_fminimizer_f_lower = gsl_min_fminimizer_f_lower(s%gsl_min_fminimizer)
  end function fgsl_min_fminimizer_f_lower
  function fgsl_min_fminimizer_f_upper(s)
    type(fgsl_min_fminimizer), intent(in) :: s
    real(fgsl_double) ::  fgsl_min_fminimizer_f_upper
    fgsl_min_fminimizer_f_upper = gsl_min_fminimizer_f_upper(s%gsl_min_fminimizer)
  end function fgsl_min_fminimizer_f_upper
  function fgsl_min_test_interval(x_lower, x_upper, epsabs, epsrel) 
    real(fgsl_double), intent(in) :: x_lower, x_upper, epsabs, epsrel
    integer(fgsl_int) :: fgsl_min_test_interval
    fgsl_min_test_interval = gsl_min_test_interval(x_lower, x_upper, epsabs, epsrel)
  end function fgsl_min_test_interval
!
  function fgsl_min_fminimizer_status(s)
    type(fgsl_min_fminimizer), intent(in) :: s
    logical :: fgsl_min_fminimizer_status
    fgsl_min_fminimizer_status = .false.
    if (c_associated(s%gsl_min_fminimizer)) &
         fgsl_min_fminimizer_status = .true.
  end function fgsl_min_fminimizer_status
 
!-*-f90-*-
!
!  API: multi-dimensional root finding
!
  function fgsl_multiroot_function_init(func, ndim, params)
    interface
       function func(x, params, f) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: x, params, f
         integer(c_int) :: func
       end function func
    end interface
    integer(fgsl_size_t), intent(in) :: ndim
    type(c_ptr), intent(in) :: params
    type(fgsl_multiroot_function) :: fgsl_multiroot_function_init
!
    type(c_funptr) :: fp
    fp = c_funloc(func)
    fgsl_multiroot_function_init%gsl_multiroot_function = &
         fgsl_multiroot_function_cinit(fp, ndim, params)
  end function fgsl_multiroot_function_init
  function fgsl_multiroot_function_fdf_init(func, dfunc, fdfunc, ndim, params)
    interface
       function func(x, params, f) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: x, params, f
         integer(c_int) :: func
       end function func
       function dfunc(x, params, df) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: x, params, df
         integer(c_int) :: dfunc
       end function dfunc
       function fdfunc(x, params, f, df) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: x, params, f, df
         integer(c_int) :: fdfunc
       end function fdfunc
    end interface
    integer(fgsl_size_t), intent(in) :: ndim
    type(c_ptr), intent(in) :: params
    type(fgsl_multiroot_function_fdf) :: fgsl_multiroot_function_fdf_init
!
    type(c_funptr) :: fp, dfp, fdfp
    fp = c_funloc(func)
    dfp = c_funloc(dfunc)
    fdfp = c_funloc(fdfunc)
    fgsl_multiroot_function_fdf_init%gsl_multiroot_function_fdf = &
         fgsl_multiroot_function_fdf_cinit(fp, dfp, fdfp, ndim, params)
  end function fgsl_multiroot_function_fdf_init
  subroutine fgsl_multiroot_function_free(fun) 
    type(fgsl_multiroot_function), intent(inout) :: fun
    call fgsl_multiroot_function_cfree(fun%gsl_multiroot_function)
  end subroutine fgsl_multiroot_function_free
  subroutine fgsl_multiroot_function_fdf_free(fun) 
    type(fgsl_multiroot_function_fdf), intent(inout) :: fun
    call fgsl_multiroot_function_fdf_cfree(fun%gsl_multiroot_function_fdf)
  end subroutine fgsl_multiroot_function_fdf_free
  function fgsl_multiroot_fsolver_alloc(t, n) 
    type(fgsl_multiroot_fsolver_type), intent(in) :: t
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_multiroot_fsolver) :: fgsl_multiroot_fsolver_alloc
! 
    type(c_ptr) :: ftype
    ftype = fgsl_aux_multiroot_fsolver_alloc(t%which)
    fgsl_multiroot_fsolver_alloc%gsl_multiroot_fsolver = &
         gsl_multiroot_fsolver_alloc(ftype,n)
  end function fgsl_multiroot_fsolver_alloc
  function fgsl_multiroot_fdfsolver_alloc(t, n) 
    type(fgsl_multiroot_fdfsolver_type), intent(in) :: t
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_multiroot_fdfsolver) :: fgsl_multiroot_fdfsolver_alloc
! 
    type(c_ptr) :: ftype
    ftype = fgsl_aux_multiroot_fdfsolver_alloc(t%which)
    fgsl_multiroot_fdfsolver_alloc%gsl_multiroot_fdfsolver = &
         gsl_multiroot_fdfsolver_alloc(ftype,n)
  end function fgsl_multiroot_fdfsolver_alloc
  subroutine fgsl_multiroot_fsolver_free(s) 
    type(fgsl_multiroot_fsolver), intent(inout) :: s
    call gsl_multiroot_fsolver_free(s%gsl_multiroot_fsolver)
  end subroutine fgsl_multiroot_fsolver_free
  subroutine fgsl_multiroot_fdfsolver_free(s) 
    type(fgsl_multiroot_fdfsolver), intent(inout) :: s
    call gsl_multiroot_fdfsolver_free(s%gsl_multiroot_fdfsolver)
  end subroutine fgsl_multiroot_fdfsolver_free
  function fgsl_multiroot_fsolver_set(s, f, x)
    type(fgsl_multiroot_fsolver), intent(inout) :: s
    type(fgsl_multiroot_function), intent(in)  :: f
    type(fgsl_vector), intent(in)  :: x
    integer(fgsl_int) :: fgsl_multiroot_fsolver_set
    fgsl_multiroot_fsolver_set = gsl_multiroot_fsolver_set(s%gsl_multiroot_fsolver, &
         f%gsl_multiroot_function, x%gsl_vector)
  end function fgsl_multiroot_fsolver_set
  function fgsl_multiroot_fdfsolver_set(s, fdf, x)
    type(fgsl_multiroot_fdfsolver), intent(inout) :: s
    type(fgsl_multiroot_function_fdf), intent(in)  :: fdf
    type(fgsl_vector), intent(in)  :: x
    integer(fgsl_int) :: fgsl_multiroot_fdfsolver_set
    fgsl_multiroot_fdfsolver_set = gsl_multiroot_fdfsolver_set(s%gsl_multiroot_fdfsolver, &
         fdf%gsl_multiroot_function_fdf, x%gsl_vector)
  end function fgsl_multiroot_fdfsolver_set
  function fgsl_multiroot_fsolver_name(s)
    type(fgsl_multiroot_fsolver), intent(in) :: s
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_multiroot_fsolver_name
!
    type(c_ptr) :: name
!
    name = gsl_multiroot_fsolver_name(s%gsl_multiroot_fsolver)
    fgsl_multiroot_fsolver_name = fgsl_name(name)
  end function fgsl_multiroot_fsolver_name
  function fgsl_multiroot_fdfsolver_name(s)
    type(fgsl_multiroot_fdfsolver), intent(in) :: s
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_multiroot_fdfsolver_name
!
    type(c_ptr) :: name
!
    name = gsl_multiroot_fdfsolver_name(s%gsl_multiroot_fdfsolver)
    fgsl_multiroot_fdfsolver_name = fgsl_name(name)
  end function fgsl_multiroot_fdfsolver_name
  function fgsl_multiroot_fsolver_iterate(s)
    type(fgsl_multiroot_fsolver), intent(in) :: s
    integer(fgsl_int) :: fgsl_multiroot_fsolver_iterate
    fgsl_multiroot_fsolver_iterate = &
         gsl_multiroot_fsolver_iterate(s%gsl_multiroot_fsolver)
  end function fgsl_multiroot_fsolver_iterate
  function fgsl_multiroot_fdfsolver_iterate(s)
    type(fgsl_multiroot_fdfsolver), intent(in) :: s
    integer(fgsl_int) :: fgsl_multiroot_fdfsolver_iterate
    fgsl_multiroot_fdfsolver_iterate = &
         gsl_multiroot_fdfsolver_iterate(s%gsl_multiroot_fdfsolver)
  end function fgsl_multiroot_fdfsolver_iterate
  function fgsl_multiroot_fsolver_root(s)
    type(fgsl_multiroot_fsolver), intent(in) :: s
    type(fgsl_vector) :: fgsl_multiroot_fsolver_root
    fgsl_multiroot_fsolver_root%gsl_vector = &
         gsl_multiroot_fsolver_root(s%gsl_multiroot_fsolver)
  end function fgsl_multiroot_fsolver_root
  function fgsl_multiroot_fdfsolver_root(s)
    type(fgsl_multiroot_fdfsolver), intent(in) :: s
    type(fgsl_vector) :: fgsl_multiroot_fdfsolver_root
    fgsl_multiroot_fdfsolver_root%gsl_vector = &
         gsl_multiroot_fdfsolver_root(s%gsl_multiroot_fdfsolver)
  end function fgsl_multiroot_fdfsolver_root
  function fgsl_multiroot_fsolver_f(s)
    type(fgsl_multiroot_fsolver), intent(in) :: s
    type(fgsl_vector) :: fgsl_multiroot_fsolver_f
    fgsl_multiroot_fsolver_f%gsl_vector = &
         gsl_multiroot_fsolver_f(s%gsl_multiroot_fsolver)
  end function fgsl_multiroot_fsolver_f
  function fgsl_multiroot_fdfsolver_f(s)
    type(fgsl_multiroot_fdfsolver), intent(in) :: s
    type(fgsl_vector) :: fgsl_multiroot_fdfsolver_f
    fgsl_multiroot_fdfsolver_f%gsl_vector = &
         gsl_multiroot_fdfsolver_f(s%gsl_multiroot_fdfsolver)
  end function fgsl_multiroot_fdfsolver_f
  function fgsl_multiroot_fsolver_dx(s)
    type(fgsl_multiroot_fsolver), intent(in) :: s
    type(fgsl_vector) :: fgsl_multiroot_fsolver_dx
    fgsl_multiroot_fsolver_dx%gsl_vector = &
         gsl_multiroot_fsolver_dx(s%gsl_multiroot_fsolver)
  end function fgsl_multiroot_fsolver_dx
  function fgsl_multiroot_fdfsolver_dx(s)
    type(fgsl_multiroot_fdfsolver), intent(in) :: s
    type(fgsl_vector) :: fgsl_multiroot_fdfsolver_dx
    fgsl_multiroot_fdfsolver_dx%gsl_vector = &
         gsl_multiroot_fdfsolver_dx(s%gsl_multiroot_fdfsolver)
  end function fgsl_multiroot_fdfsolver_dx
  function fgsl_multiroot_test_delta(dx, x, epsabs, epsrel) 
    type(fgsl_vector), intent(in) :: dx, x
    real(fgsl_double), intent(in) :: epsabs, epsrel
    integer(fgsl_int) :: fgsl_multiroot_test_delta
    fgsl_multiroot_test_delta = gsl_multiroot_test_delta(dx%gsl_vector, &
         x%gsl_vector, epsabs, epsrel)
  end function fgsl_multiroot_test_delta
  function fgsl_multiroot_test_residual(f, epsabs) 
    type(fgsl_vector), intent(in) :: f
    real(fgsl_double), intent(in) :: epsabs
    integer(fgsl_int) :: fgsl_multiroot_test_residual
    fgsl_multiroot_test_residual = &
         gsl_multiroot_test_residual(f%gsl_vector, epsabs)
  end function fgsl_multiroot_test_residual
  function fgsl_multiroot_fsolver_status(s)
    type(fgsl_multiroot_fsolver), intent(in) :: s
    logical :: fgsl_multiroot_fsolver_status
    fgsl_multiroot_fsolver_status = .false.
    if (c_associated(s%gsl_multiroot_fsolver)) &
         fgsl_multiroot_fsolver_status = .true.
  end function fgsl_multiroot_fsolver_status
  function fgsl_multiroot_fdfsolver_status(s)
    type(fgsl_multiroot_fdfsolver), intent(in) :: s
    logical :: fgsl_multiroot_fdfsolver_status
    fgsl_multiroot_fdfsolver_status = .false.
    if (c_associated(s%gsl_multiroot_fdfsolver)) &
         fgsl_multiroot_fdfsolver_status = .true.
  end function fgsl_multiroot_fdfsolver_status
  
!-*-f90-*-
!
!  API: multi-dimensional minimization
!
  function fgsl_multimin_function_init(func, ndim, params)
    interface
       function func(x, params) bind(c)
         import :: c_ptr, c_double
         type(c_ptr), value :: x, params
         real(c_double) :: func
       end function func
    end interface
    integer(fgsl_size_t), intent(in) :: ndim
    type(c_ptr), intent(in) :: params
    type(fgsl_multimin_function) :: fgsl_multimin_function_init
!
    type(c_funptr) :: fp
    fp = c_funloc(func)
    fgsl_multimin_function_init%gsl_multimin_function = &
         fgsl_multimin_function_cinit(fp, ndim, params)
  end function fgsl_multimin_function_init
  function fgsl_multimin_function_fdf_init(func, dfunc, fdfunc, ndim, params)
    interface
       function func(x, params) bind(c)
         import :: c_ptr, c_double
         type(c_ptr), value :: x, params
         real(c_double) :: func
       end function func
       subroutine dfunc(x, params, df) bind(c)
         import :: c_ptr
         type(c_ptr), value :: x, params, df
       end subroutine  dfunc
       subroutine fdfunc(x, params, f, df) bind(c)
         import :: c_ptr, c_double
         real(c_double) :: f
         type(c_ptr), value :: x, params, df
       end subroutine fdfunc
    end interface
    integer(fgsl_size_t), intent(in) :: ndim
    type(c_ptr), intent(in) :: params
    type(fgsl_multimin_function_fdf) :: fgsl_multimin_function_fdf_init
!
    type(c_funptr) :: fp, dfp, fdfp
    fp = c_funloc(func)
    dfp = c_funloc(dfunc)
    fdfp = c_funloc(fdfunc)
    fgsl_multimin_function_fdf_init%gsl_multimin_function_fdf = &
         fgsl_multimin_function_fdf_cinit(fp, dfp, fdfp, ndim, params)
  end function fgsl_multimin_function_fdf_init
  subroutine fgsl_multimin_function_free(fun) 
    type(fgsl_multimin_function), intent(inout) :: fun
    call fgsl_multimin_function_cfree(fun%gsl_multimin_function)
  end subroutine fgsl_multimin_function_free
  subroutine fgsl_multimin_function_fdf_free(fun) 
    type(fgsl_multimin_function_fdf), intent(inout) :: fun
    call fgsl_multimin_function_fdf_cfree(fun%gsl_multimin_function_fdf)
  end subroutine fgsl_multimin_function_fdf_free
  function fgsl_multimin_fminimizer_alloc(t, n) 
    type(fgsl_multimin_fminimizer_type), intent(in) :: t
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_multimin_fminimizer) :: fgsl_multimin_fminimizer_alloc
! 
    type(c_ptr) :: ftype
    ftype = fgsl_aux_multimin_fminimizer_alloc(t%which)
    fgsl_multimin_fminimizer_alloc%gsl_multimin_fminimizer = &
         gsl_multimin_fminimizer_alloc(ftype,n)
  end function fgsl_multimin_fminimizer_alloc
  function fgsl_multimin_fdfminimizer_alloc(t, n) 
    type(fgsl_multimin_fdfminimizer_type), intent(in) :: t
    integer(fgsl_size_t), intent(in) :: n
    type(fgsl_multimin_fdfminimizer) :: fgsl_multimin_fdfminimizer_alloc
! 
    type(c_ptr) :: ftype
    ftype = fgsl_aux_multimin_fdfminimizer_alloc(t%which)
    fgsl_multimin_fdfminimizer_alloc%gsl_multimin_fdfminimizer = &
         gsl_multimin_fdfminimizer_alloc(ftype,n)
  end function fgsl_multimin_fdfminimizer_alloc
  subroutine fgsl_multimin_fminimizer_free(s) 
    type(fgsl_multimin_fminimizer), intent(inout) :: s
    call gsl_multimin_fminimizer_free(s%gsl_multimin_fminimizer)
  end subroutine fgsl_multimin_fminimizer_free
  subroutine fgsl_multimin_fdfminimizer_free(s) 
    type(fgsl_multimin_fdfminimizer), intent(inout) :: s
    call gsl_multimin_fdfminimizer_free(s%gsl_multimin_fdfminimizer)
  end subroutine fgsl_multimin_fdfminimizer_free
  function fgsl_multimin_fminimizer_set(s, f, x, step)
    type(fgsl_multimin_fminimizer), intent(inout) :: s
    type(fgsl_multimin_function), intent(in)  :: f
    type(fgsl_vector), intent(in)  :: x, step
    integer(fgsl_int) :: fgsl_multimin_fminimizer_set
    fgsl_multimin_fminimizer_set = gsl_multimin_fminimizer_set(s%gsl_multimin_fminimizer, &
         f%gsl_multimin_function, x%gsl_vector, step%gsl_vector)
  end function fgsl_multimin_fminimizer_set
  function fgsl_multimin_fdfminimizer_set(s, fdf, x, step, tol)
    type(fgsl_multimin_fdfminimizer), intent(inout) :: s
    type(fgsl_multimin_function_fdf), intent(in)  :: fdf
    type(fgsl_vector), intent(in)  :: x
    real(fgsl_double), intent(in) :: step, tol
    integer(fgsl_int) :: fgsl_multimin_fdfminimizer_set
    fgsl_multimin_fdfminimizer_set = gsl_multimin_fdfminimizer_set(s%gsl_multimin_fdfminimizer, &
         fdf%gsl_multimin_function_fdf, x%gsl_vector, step, tol)
  end function fgsl_multimin_fdfminimizer_set
  function fgsl_multimin_fminimizer_name(s)
    type(fgsl_multimin_fminimizer), intent(in) :: s
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_multimin_fminimizer_name
!
    type(c_ptr) :: name
!
    name = gsl_multimin_fminimizer_name(s%gsl_multimin_fminimizer)
    fgsl_multimin_fminimizer_name = fgsl_name(name)
  end function fgsl_multimin_fminimizer_name
  function fgsl_multimin_fdfminimizer_name(s)
    type(fgsl_multimin_fdfminimizer), intent(in) :: s
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_multimin_fdfminimizer_name
!
    type(c_ptr) :: name
!
    name = gsl_multimin_fdfminimizer_name(s%gsl_multimin_fdfminimizer)
    fgsl_multimin_fdfminimizer_name = fgsl_name(name)
  end function fgsl_multimin_fdfminimizer_name
  function fgsl_multimin_fminimizer_iterate(s)
    type(fgsl_multimin_fminimizer), intent(in) :: s
    integer(fgsl_int) :: fgsl_multimin_fminimizer_iterate
    fgsl_multimin_fminimizer_iterate = &
         gsl_multimin_fminimizer_iterate(s%gsl_multimin_fminimizer)
  end function fgsl_multimin_fminimizer_iterate
  function fgsl_multimin_fdfminimizer_iterate(s)
    type(fgsl_multimin_fdfminimizer), intent(in) :: s
    integer(fgsl_int) :: fgsl_multimin_fdfminimizer_iterate
    fgsl_multimin_fdfminimizer_iterate = &
         gsl_multimin_fdfminimizer_iterate(s%gsl_multimin_fdfminimizer)
  end function fgsl_multimin_fdfminimizer_iterate
  function fgsl_multimin_fminimizer_x(s)
    type(fgsl_multimin_fminimizer), intent(in) :: s
    type(fgsl_vector) :: fgsl_multimin_fminimizer_x
    fgsl_multimin_fminimizer_x%gsl_vector = &
         gsl_multimin_fminimizer_x(s%gsl_multimin_fminimizer)
  end function fgsl_multimin_fminimizer_x
  function fgsl_multimin_fdfminimizer_x(s)
    type(fgsl_multimin_fdfminimizer), intent(in) :: s
    type(fgsl_vector) :: fgsl_multimin_fdfminimizer_x
    fgsl_multimin_fdfminimizer_x%gsl_vector = &
         gsl_multimin_fdfminimizer_x(s%gsl_multimin_fdfminimizer)
  end function fgsl_multimin_fdfminimizer_x
  function fgsl_multimin_fminimizer_minimum(s)
    type(fgsl_multimin_fminimizer), intent(in) :: s
    real(fgsl_double) :: fgsl_multimin_fminimizer_minimum
    fgsl_multimin_fminimizer_minimum = &
         gsl_multimin_fminimizer_minimum(s%gsl_multimin_fminimizer)
  end function fgsl_multimin_fminimizer_minimum
  function fgsl_multimin_fdfminimizer_minimum(s)
    type(fgsl_multimin_fdfminimizer), intent(in) :: s
    real(fgsl_double) :: fgsl_multimin_fdfminimizer_minimum
    fgsl_multimin_fdfminimizer_minimum = &
         gsl_multimin_fdfminimizer_minimum(s%gsl_multimin_fdfminimizer)
  end function fgsl_multimin_fdfminimizer_minimum
  function fgsl_multimin_fdfminimizer_gradient(s)
    type(fgsl_multimin_fdfminimizer), intent(in) :: s
    type(fgsl_vector) :: fgsl_multimin_fdfminimizer_gradient
    fgsl_multimin_fdfminimizer_gradient%gsl_vector = &
         gsl_multimin_fdfminimizer_gradient(s%gsl_multimin_fdfminimizer)
  end function fgsl_multimin_fdfminimizer_gradient
  function fgsl_multimin_fminimizer_size(s)
    type(fgsl_multimin_fminimizer), intent(in) :: s
    real(fgsl_double) :: fgsl_multimin_fminimizer_size
    fgsl_multimin_fminimizer_size = &
         gsl_multimin_fminimizer_size(s%gsl_multimin_fminimizer)
  end function fgsl_multimin_fminimizer_size
  function fgsl_multimin_fdfminimizer_restart(s)
    type(fgsl_multimin_fdfminimizer), intent(in) :: s
    integer(fgsl_int) :: fgsl_multimin_fdfminimizer_restart
    fgsl_multimin_fdfminimizer_restart = &
         gsl_multimin_fdfminimizer_restart(s%gsl_multimin_fdfminimizer)
  end function fgsl_multimin_fdfminimizer_restart
  function fgsl_multimin_test_gradient(g, epsabs) 
    type(fgsl_vector), intent(in) :: g
    real(fgsl_double), intent(in) :: epsabs
    integer(fgsl_int) :: fgsl_multimin_test_gradient
    fgsl_multimin_test_gradient = gsl_multimin_test_gradient(g%gsl_vector, epsabs)
  end function fgsl_multimin_test_gradient
  function fgsl_multimin_test_size(size, epsabs) 
    real(fgsl_double), intent(in) :: size, epsabs
    integer(fgsl_int) :: fgsl_multimin_test_size
    fgsl_multimin_test_size = &
         gsl_multimin_test_size(size, epsabs)
  end function fgsl_multimin_test_size
  function fgsl_multimin_fminimizer_status(s)
    type(fgsl_multimin_fminimizer), intent(in) :: s
    logical :: fgsl_multimin_fminimizer_status
    fgsl_multimin_fminimizer_status = .false.
    if (c_associated(s%gsl_multimin_fminimizer)) &
         fgsl_multimin_fminimizer_status = .true.
  end function fgsl_multimin_fminimizer_status
  function fgsl_multimin_fdfminimizer_status(s)
    type(fgsl_multimin_fdfminimizer), intent(in) :: s
    logical :: fgsl_multimin_fdfminimizer_status
    fgsl_multimin_fdfminimizer_status = .false.
    if (c_associated(s%gsl_multimin_fdfminimizer)) &
         fgsl_multimin_fdfminimizer_status = .true.
  end function fgsl_multimin_fdfminimizer_status
  
!-*-f90-*-
!
!  API: Fitting
!
  function fgsl_fit_linear(x, xstride, y, ystride, n, c0, c1, &
       cov00, cov01, cov11, sumsq)
    real(fgsl_double), intent(in) :: x(:), y(:)
    integer(fgsl_size_t), intent(in) :: xstride, ystride, n
    real(fgsl_double), intent(out) :: c0, c1, cov00, cov01, cov11, sumsq
    integer(fgsl_int) :: fgsl_fit_linear
    fgsl_fit_linear = gsl_fit_linear(x, xstride, y, ystride, n, c0, c1, &
       cov00, cov01, cov11, sumsq)
  end function fgsl_fit_linear
  function fgsl_fit_wlinear(x, xstride, w, wstride, y, ystride, n, c0, c1, &
       cov00, cov01, cov11, chisq)
    real(fgsl_double), intent(in) :: x(:), y(:), w(:)
    integer(fgsl_size_t), intent(in) :: xstride, ystride, wstride, n
    real(fgsl_double), intent(out) :: c0, c1, cov00, cov01, cov11, chisq
    integer(fgsl_int) :: fgsl_fit_wlinear
    fgsl_fit_wlinear = gsl_fit_wlinear(x, xstride, w, wstride, y, ystride, n, c0, c1, &
       cov00, cov01, cov11, chisq)
  end function fgsl_fit_wlinear
  function fgsl_fit_linear_est(x, c0, c1, cov00, cov01, cov11, y, y_err)
    real(fgsl_double), intent(in) :: x, c0, c1, cov00, cov01, cov11
    real(fgsl_double), intent(out) ::  y, y_err
    integer(fgsl_int) :: fgsl_fit_linear_est
    fgsl_fit_linear_est = gsl_fit_linear_est(x, c0, c1, cov00, cov01, cov11, y, y_err)
  end function fgsl_fit_linear_est
  function fgsl_fit_mul(x, xstride, y, ystride, n, c1, &
        cov11, sumsq)
    real(fgsl_double), intent(in) :: x(:), y(:)
    integer(fgsl_size_t), intent(in) :: xstride, ystride, n
    real(fgsl_double), intent(out) :: c1,  cov11, sumsq
    integer(fgsl_int) :: fgsl_fit_mul
    fgsl_fit_mul = gsl_fit_mul(x, xstride, y, ystride, n, c1, &
        cov11, sumsq)
  end function fgsl_fit_mul
  function fgsl_fit_wmul(x, xstride, w, wstride, y, ystride, n, c1, &
        cov11, chisq)
    real(fgsl_double), intent(in) :: x(:), y(:), w(:)
    integer(fgsl_size_t), intent(in) :: xstride, ystride, wstride, n
    real(fgsl_double), intent(out) :: c1, cov11, chisq
    integer(fgsl_int) :: fgsl_fit_wmul
    fgsl_fit_wmul = gsl_fit_wmul(x, xstride, w, wstride, y, ystride, n, c1, &
        cov11, chisq)
  end function fgsl_fit_wmul
  function fgsl_fit_mul_est(x, c1, cov11, y, y_err)
    real(fgsl_double), intent(in) :: x, c1, cov11
    real(fgsl_double), intent(out) ::  y, y_err
    integer(fgsl_int) :: fgsl_fit_mul_est
    fgsl_fit_mul_est = gsl_fit_mul_est(x, c1, cov11, y, y_err)
  end function fgsl_fit_mul_est
  function fgsl_multifit_linear_alloc(n, p) 
    integer(fgsl_size_t), intent(in) :: n, p
    type(fgsl_multifit_linear_workspace) :: fgsl_multifit_linear_alloc
    fgsl_multifit_linear_alloc%gsl_multifit_linear_workspace = &
         gsl_multifit_linear_alloc(n, p) 
  end function fgsl_multifit_linear_alloc
  subroutine fgsl_multifit_linear_free(w) 
    type(fgsl_multifit_linear_workspace), intent(inout) :: w
    call gsl_multifit_linear_free(w%gsl_multifit_linear_workspace) 
  end subroutine fgsl_multifit_linear_free
  function fgsl_multifit_linear(x, y, c, cov, chisq, work)
    type(fgsl_matrix), intent(in) :: x
    type(fgsl_vector), intent(in) :: y
    type(fgsl_vector), intent(inout) :: c
    type(fgsl_matrix), intent(inout) :: cov
    type(fgsl_multifit_linear_workspace), intent(inout) :: work
    real(fgsl_double), intent(inout) :: chisq
    integer(fgsl_int) :: fgsl_multifit_linear
    fgsl_multifit_linear = gsl_multifit_linear(x%gsl_matrix, y%gsl_vector, &
         c%gsl_vector, cov%gsl_matrix, chisq, work%gsl_multifit_linear_workspace)
  end function fgsl_multifit_linear
  function fgsl_multifit_linear_svd(x, y, tol, rank, c, cov, chisq, work)
    type(fgsl_matrix), intent(in) :: x
    type(fgsl_vector), intent(in) :: y
    real(fgsl_double), intent(in) :: tol
    integer(fgsl_size_t), intent(inout) :: rank
    type(fgsl_vector), intent(inout) :: c
    type(fgsl_matrix), intent(inout) :: cov
    type(fgsl_multifit_linear_workspace), intent(inout) :: work
    real(fgsl_double), intent(inout) :: chisq
    integer(fgsl_int) :: fgsl_multifit_linear_svd
    fgsl_multifit_linear_svd = gsl_multifit_linear_svd(x%gsl_matrix, y%gsl_vector, tol, rank, &
         c%gsl_vector, cov%gsl_matrix, chisq, work%gsl_multifit_linear_workspace)
  end function fgsl_multifit_linear_svd
  function fgsl_multifit_wlinear(x, w, y, c, cov, chisq, work)
    type(fgsl_matrix), intent(in) :: x
    type(fgsl_vector), intent(in) :: w, y
    type(fgsl_vector), intent(inout) :: c
    type(fgsl_matrix), intent(inout) :: cov
    type(fgsl_multifit_linear_workspace), intent(inout) :: work
    real(fgsl_double), intent(inout) :: chisq
    integer(fgsl_int) :: fgsl_multifit_wlinear
    fgsl_multifit_wlinear = gsl_multifit_wlinear(x%gsl_matrix, w%gsl_vector, y%gsl_vector, &
         c%gsl_vector, cov%gsl_matrix, chisq, work%gsl_multifit_linear_workspace)
  end function fgsl_multifit_wlinear
  function fgsl_multifit_wlinear_svd(x, w, y, tol, rank, c, cov, chisq, work)
    type(fgsl_matrix), intent(in) :: x
    type(fgsl_vector), intent(in) :: w, y
    real(fgsl_double), intent(in) :: tol
    integer(fgsl_size_t), intent(inout) :: rank
    type(fgsl_vector), intent(inout) :: c
    type(fgsl_matrix), intent(inout) :: cov
    type(fgsl_multifit_linear_workspace), intent(inout) :: work
    real(fgsl_double), intent(inout) :: chisq
    integer(fgsl_int) :: fgsl_multifit_wlinear_svd
    fgsl_multifit_wlinear_svd = gsl_multifit_wlinear_svd(x%gsl_matrix, w%gsl_vector, &
         y%gsl_vector, tol, rank, c%gsl_vector, cov%gsl_matrix, chisq, &
         work%gsl_multifit_linear_workspace)
  end function fgsl_multifit_wlinear_svd
  function fgsl_multifit_linear_est(x, c, cov, y, y_err)
    type(fgsl_vector), intent(in) :: x, c
    type(fgsl_matrix), intent(in) :: cov
    real(fgsl_double), intent(inout) :: y, y_err
    integer(fgsl_int) :: fgsl_multifit_linear_est
    fgsl_multifit_linear_est = gsl_multifit_linear_est(x%gsl_vector, c%gsl_vector, &
         cov%gsl_matrix, y, y_err)
  end function fgsl_multifit_linear_est
  function fgsl_multifit_linear_residuals(x, y, c, r)
    type(fgsl_matrix), intent(in) :: x
    type(fgsl_vector), intent(in) :: y, c
    type(fgsl_vector), intent(inout) :: r
    integer(fgsl_int) :: fgsl_multifit_linear_residuals
    fgsl_multifit_linear_residuals = gsl_multifit_linear_residuals( &
         x%gsl_matrix, y%gsl_vector, c%gsl_vector, r%gsl_vector)
  end function fgsl_multifit_linear_residuals
  function fgsl_multifit_status(multifit)
    type(fgsl_multifit_linear_workspace), intent(in) :: multifit
    logical :: fgsl_multifit_status
    fgsl_multifit_status = .true.
    if (.not. c_associated(multifit%gsl_multifit_linear_workspace)) fgsl_multifit_status = .false.
  end function fgsl_multifit_status
!-*-f90-*-
!
!  API: non-linear least-squares fitting
!
  function fgsl_multifit_function_init(func, ndim, p, params)
    interface
       function func(x, params, f) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: x, params, f
         integer(c_int) :: func
       end function func
    end interface
    integer(fgsl_size_t), intent(in) :: ndim, p
    type(c_ptr), intent(in) :: params
    type(fgsl_multifit_function) :: fgsl_multifit_function_init
!
    type(c_funptr) :: fp
    fp = c_funloc(func)
    fgsl_multifit_function_init%gsl_multifit_function = &
         fgsl_multifit_function_cinit(fp, ndim, p, params)
  end function fgsl_multifit_function_init
  function fgsl_multifit_function_fdf_init(func, dfunc, fdfunc, ndim, p, params)
    interface
       function func(x, params, f) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: x, params, f
         integer(c_int) :: func
       end function func
       function dfunc(x, params, df) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: x, params, df
         integer(c_int) :: dfunc
       end function dfunc
       function fdfunc(x, params, f, df) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: x, params, f, df
         integer(c_int) :: fdfunc
       end function fdfunc
    end interface
    integer(fgsl_size_t), intent(in) :: ndim, p
    type(c_ptr), intent(in) :: params
    type(fgsl_multifit_function_fdf) :: fgsl_multifit_function_fdf_init
!
    type(c_funptr) :: fp, dfp, fdfp
    fp = c_funloc(func)
    dfp = c_funloc(dfunc)
    fdfp = c_funloc(fdfunc)
    fgsl_multifit_function_fdf_init%gsl_multifit_function_fdf = &
         fgsl_multifit_function_fdf_cinit(fp, dfp, fdfp, ndim, p, params)
  end function fgsl_multifit_function_fdf_init
  subroutine fgsl_multifit_function_free(fun) 
    type(fgsl_multifit_function), intent(inout) :: fun
    call fgsl_multifit_function_cfree(fun%gsl_multifit_function)
  end subroutine fgsl_multifit_function_free
  subroutine fgsl_multifit_function_fdf_free(fun) 
    type(fgsl_multifit_function_fdf), intent(inout) :: fun
    call fgsl_multifit_function_fdf_cfree(fun%gsl_multifit_function_fdf)
  end subroutine fgsl_multifit_function_fdf_free
  function fgsl_multifit_fsolver_alloc(t, n, p) 
    type(fgsl_multifit_fsolver_type), intent(in) :: t
    integer(fgsl_size_t), intent(in) :: n, p
    type(fgsl_multifit_fsolver) :: fgsl_multifit_fsolver_alloc
! 
    type(c_ptr) :: ftype
    ftype = fgsl_aux_multifit_fsolver_alloc(t%which)
    fgsl_multifit_fsolver_alloc%gsl_multifit_fsolver = &
         gsl_multifit_fsolver_alloc(ftype, n, p)
  end function fgsl_multifit_fsolver_alloc
  function fgsl_multifit_fdfsolver_alloc(t, n, p) 
    type(fgsl_multifit_fdfsolver_type), intent(in) :: t
    integer(fgsl_size_t), intent(in) :: n, p
    type(fgsl_multifit_fdfsolver) :: fgsl_multifit_fdfsolver_alloc
! 
    type(c_ptr) :: ftype
    ftype = fgsl_aux_multifit_fdfsolver_alloc(t%which)
    fgsl_multifit_fdfsolver_alloc%gsl_multifit_fdfsolver = &
         gsl_multifit_fdfsolver_alloc(ftype, n, p)
  end function fgsl_multifit_fdfsolver_alloc
  subroutine fgsl_multifit_fsolver_free(s) 
    type(fgsl_multifit_fsolver), intent(inout) :: s
    call gsl_multifit_fsolver_free(s%gsl_multifit_fsolver)
  end subroutine fgsl_multifit_fsolver_free
  subroutine fgsl_multifit_fdfsolver_free(s) 
    type(fgsl_multifit_fdfsolver), intent(inout) :: s
    call gsl_multifit_fdfsolver_free(s%gsl_multifit_fdfsolver)
  end subroutine fgsl_multifit_fdfsolver_free
  function fgsl_multifit_fsolver_set(s, f, x)
    type(fgsl_multifit_fsolver), intent(inout) :: s
    type(fgsl_multifit_function), intent(in)  :: f
    type(fgsl_vector), intent(in)  :: x
    integer(fgsl_int) :: fgsl_multifit_fsolver_set
    fgsl_multifit_fsolver_set = gsl_multifit_fsolver_set(s%gsl_multifit_fsolver, &
         f%gsl_multifit_function, x%gsl_vector)
  end function fgsl_multifit_fsolver_set
  function fgsl_multifit_fdfsolver_set(s, fdf, x)
    type(fgsl_multifit_fdfsolver), intent(inout) :: s
    type(fgsl_multifit_function_fdf), intent(in)  :: fdf
    type(fgsl_vector), intent(in)  :: x
    integer(fgsl_int) :: fgsl_multifit_fdfsolver_set
    fgsl_multifit_fdfsolver_set = gsl_multifit_fdfsolver_set(s%gsl_multifit_fdfsolver, &
         fdf%gsl_multifit_function_fdf, x%gsl_vector)
  end function fgsl_multifit_fdfsolver_set
  function fgsl_multifit_fsolver_name(s)
    type(fgsl_multifit_fsolver), intent(in) :: s
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_multifit_fsolver_name
!
    type(c_ptr) :: name
!
    name = gsl_multifit_fsolver_name(s%gsl_multifit_fsolver)
    fgsl_multifit_fsolver_name = fgsl_name(name)
  end function fgsl_multifit_fsolver_name
  function fgsl_multifit_fdfsolver_name(s)
    type(fgsl_multifit_fdfsolver), intent(in) :: s
    character(kind=fgsl_char,len=fgsl_strmax) :: fgsl_multifit_fdfsolver_name
!
    type(c_ptr) :: name
!
    name = gsl_multifit_fdfsolver_name(s%gsl_multifit_fdfsolver)
    fgsl_multifit_fdfsolver_name = fgsl_name(name)
  end function fgsl_multifit_fdfsolver_name
  function fgsl_multifit_fsolver_iterate(s)
    type(fgsl_multifit_fsolver), intent(in) :: s
    integer(fgsl_int) :: fgsl_multifit_fsolver_iterate
    fgsl_multifit_fsolver_iterate = &
         gsl_multifit_fsolver_iterate(s%gsl_multifit_fsolver)
  end function fgsl_multifit_fsolver_iterate
  function fgsl_multifit_fdfsolver_iterate(s)
    type(fgsl_multifit_fdfsolver), intent(in) :: s
    integer(fgsl_int) :: fgsl_multifit_fdfsolver_iterate
    fgsl_multifit_fdfsolver_iterate = &
         gsl_multifit_fdfsolver_iterate(s%gsl_multifit_fdfsolver)
  end function fgsl_multifit_fdfsolver_iterate
  function fgsl_multifit_fsolver_position(s)
    type(fgsl_multifit_fsolver), intent(in) :: s
    type(fgsl_vector) :: fgsl_multifit_fsolver_position
    fgsl_multifit_fsolver_position%gsl_vector = &
         gsl_multifit_fsolver_position(s%gsl_multifit_fsolver)
  end function fgsl_multifit_fsolver_position
  function fgsl_multifit_fdfsolver_position(s)
    type(fgsl_multifit_fdfsolver), intent(in) :: s
    type(fgsl_vector) :: fgsl_multifit_fdfsolver_position
    fgsl_multifit_fdfsolver_position%gsl_vector = &
         gsl_multifit_fdfsolver_position(s%gsl_multifit_fdfsolver)
  end function fgsl_multifit_fdfsolver_position
  function fgsl_multifit_fdfsolver_dx(s)
    type(fgsl_multifit_fdfsolver), intent(in) :: s
    type(fgsl_vector) :: fgsl_multifit_fdfsolver_dx
    fgsl_multifit_fdfsolver_dx%gsl_vector = &
         gsl_multifit_fdfsolver_dx(s%gsl_multifit_fdfsolver)
  end function fgsl_multifit_fdfsolver_dx
  function fgsl_multifit_fdfsolver_f(s)
    type(fgsl_multifit_fdfsolver), intent(in) :: s
    type(fgsl_vector) :: fgsl_multifit_fdfsolver_f
    fgsl_multifit_fdfsolver_f%gsl_vector = &
         gsl_multifit_fdfsolver_f(s%gsl_multifit_fdfsolver)
  end function fgsl_multifit_fdfsolver_f
  function fgsl_multifit_fdfsolver_jac(s)
    type(fgsl_multifit_fdfsolver), intent(in) :: s
    type(fgsl_matrix) :: fgsl_multifit_fdfsolver_jac
    fgsl_multifit_fdfsolver_jac%gsl_matrix = &
         gsl_multifit_fdfsolver_jac(s%gsl_multifit_fdfsolver)
  end function fgsl_multifit_fdfsolver_jac
  function fgsl_multifit_test_delta(dx, x, epsabs, epsrel) 
    type(fgsl_vector), intent(in) :: dx, x
    real(fgsl_double), intent(in) :: epsabs, epsrel
    integer(fgsl_int) :: fgsl_multifit_test_delta
    fgsl_multifit_test_delta = gsl_multifit_test_delta(dx%gsl_vector, x%gsl_vector, &
         epsabs, epsrel)
  end function fgsl_multifit_test_delta
  function fgsl_multifit_test_gradient(g, epsabs) 
    type(fgsl_vector), intent(in) :: g
    real(fgsl_double), intent(in) :: epsabs
    integer(fgsl_int) :: fgsl_multifit_test_gradient
    fgsl_multifit_test_gradient = gsl_multifit_test_gradient(g%gsl_vector, epsabs)
  end function fgsl_multifit_test_gradient
  function fgsl_multifit_gradient(j, f, g) 
    type(fgsl_matrix), intent(in) :: j
    type(fgsl_vector), intent(in) :: f
    type(fgsl_vector), intent(inout) :: g
    integer(fgsl_int) :: fgsl_multifit_gradient
    fgsl_multifit_gradient = gsl_multifit_gradient(j%gsl_matrix, f%gsl_vector, &
         g%gsl_vector)
  end function fgsl_multifit_gradient
  function fgsl_multifit_covar(j, epsrel, covar) 
    type(fgsl_matrix), intent(in) :: j
    real(fgsl_double), intent(in) :: epsrel
    type(fgsl_matrix), intent(inout) :: covar
    integer(fgsl_int) :: fgsl_multifit_covar
    fgsl_multifit_covar = gsl_multifit_covar(j%gsl_matrix, epsrel, &
         covar%gsl_matrix)
  end function fgsl_multifit_covar
  function fgsl_multifit_fsolver_status(s)
    type(fgsl_multifit_fsolver), intent(in) :: s
    logical :: fgsl_multifit_fsolver_status
    fgsl_multifit_fsolver_status = .false.
    if (c_associated(s%gsl_multifit_fsolver)) &
         fgsl_multifit_fsolver_status = .true.
  end function fgsl_multifit_fsolver_status
  function fgsl_multifit_fdfsolver_status(s)
    type(fgsl_multifit_fdfsolver), intent(in) :: s
    logical :: fgsl_multifit_fdfsolver_status
    fgsl_multifit_fdfsolver_status = .false.
    if (c_associated(s%gsl_multifit_fdfsolver)) &
         fgsl_multifit_fdfsolver_status = .true.
  end function fgsl_multifit_fdfsolver_status
  
!-*-f90-*-
!
! API: Basis splines
!
function fgsl_bspline_alloc(k, nbreak)
  integer(fgsl_size_t), intent(in) :: k, nbreak
  type(fgsl_bspline_workspace) :: fgsl_bspline_alloc
!
  fgsl_bspline_alloc%gsl_bspline_workspace = gsl_bspline_alloc(k, nbreak)
end function fgsl_bspline_alloc
subroutine fgsl_bspline_free (w)
  type(fgsl_bspline_workspace), intent(inout) :: w
  call gsl_bspline_free (w%gsl_bspline_workspace)
end subroutine fgsl_bspline_free
function fgsl_bspline_deriv_alloc(k)
  integer(fgsl_size_t), intent(in) :: k
  type(fgsl_bspline_deriv_workspace) :: fgsl_bspline_deriv_alloc
!
  fgsl_bspline_deriv_alloc%gsl_bspline_deriv_workspace =  &
       gsl_bspline_deriv_alloc(k)
end function fgsl_bspline_deriv_alloc
subroutine fgsl_bspline_deriv_free (w)
  type(fgsl_bspline_deriv_workspace), intent(inout) :: w
  call gsl_bspline_deriv_free (w%gsl_bspline_deriv_workspace)
end subroutine fgsl_bspline_deriv_free
function fgsl_bspline_knots(breakpts, w)
  integer(fgsl_int) :: fgsl_bspline_knots
  type(fgsl_vector), intent(in) :: breakpts
  type(fgsl_bspline_workspace), intent(inout) :: w
  fgsl_bspline_knots = gsl_bspline_knots(breakpts%gsl_vector, &
       w%gsl_bspline_workspace)
end function fgsl_bspline_knots
function fgsl_bspline_knots_uniform(a, b, w)
  integer(fgsl_int) :: fgsl_bspline_knots_uniform
  real(fgsl_double), intent(in) :: a, b
  type(fgsl_bspline_workspace), intent(inout) :: w
  fgsl_bspline_knots_uniform = gsl_bspline_knots_uniform(a, b, &
       w%gsl_bspline_workspace)
end function fgsl_bspline_knots_uniform
function fgsl_bspline_eval(x, b, w)
  integer(fgsl_int) :: fgsl_bspline_eval
  real(fgsl_double), intent(in) :: x
  type(fgsl_vector), intent(inout) :: b
  type(fgsl_bspline_workspace), intent(inout) :: w
  fgsl_bspline_eval = gsl_bspline_eval(x, b%gsl_vector, &
       w%gsl_bspline_workspace)
end function fgsl_bspline_eval
function fgsl_bspline_eval_nonzero(x, bk, istart, iend, w)
  integer(fgsl_int) :: fgsl_bspline_eval_nonzero
  real(fgsl_double), intent(in) :: x
  type(fgsl_vector), intent(inout) :: bk
  integer(fgsl_size_t), intent(inout) :: istart, iend
  type(fgsl_bspline_workspace), intent(inout) :: w
  fgsl_bspline_eval_nonzero = gsl_bspline_eval_nonzero(x, bk%gsl_vector, &
       istart, iend, w%gsl_bspline_workspace)
end function fgsl_bspline_eval_nonzero
function fgsl_bspline_deriv_eval(x, nderiv, db, w, dw)
  integer(fgsl_int) :: fgsl_bspline_deriv_eval
  real(fgsl_double), intent(in) :: x
  integer(fgsl_size_t), intent(in) :: nderiv
  type(fgsl_matrix), intent(inout) :: db
  type(fgsl_bspline_workspace), intent(inout) :: w
  type(fgsl_bspline_deriv_workspace), intent(inout) :: dw
  fgsl_bspline_deriv_eval = gsl_bspline_deriv_eval(x, nderiv, db%gsl_matrix, &
       w%gsl_bspline_workspace, dw%gsl_bspline_deriv_workspace)
end function fgsl_bspline_deriv_eval
function fgsl_bspline_deriv_eval_nonzero(x, nderiv, db, istart, iend, w, dw)
  integer(fgsl_int) :: fgsl_bspline_deriv_eval_nonzero
  real(fgsl_double), intent(in) :: x
  integer(fgsl_size_t), intent(in) :: nderiv
  type(fgsl_matrix), intent(inout) :: db
  integer(fgsl_size_t), intent(inout) :: istart, iend
  type(fgsl_bspline_workspace), intent(inout) :: w
  type(fgsl_bspline_deriv_workspace), intent(inout) :: dw
  fgsl_bspline_deriv_eval_nonzero = gsl_bspline_deriv_eval_nonzero(x, nderiv, &
       db%gsl_matrix, istart, iend, w%gsl_bspline_workspace, &
       dw%gsl_bspline_deriv_workspace)
end function fgsl_bspline_deriv_eval_nonzero
function fgsl_bspline_ncoeffs(w)
  integer(fgsl_size_t) :: fgsl_bspline_ncoeffs
  type(fgsl_bspline_workspace), intent(inout) :: w
  fgsl_bspline_ncoeffs = gsl_bspline_ncoeffs(w%gsl_bspline_workspace)
end function fgsl_bspline_ncoeffs
function fgsl_bspline_greville_abscissa(i, w) 
  real(fgsl_double) :: fgsl_bspline_greville_abscissa
  integer(fgsl_size_t) :: i
  type(fgsl_bspline_workspace), intent(in) :: w
  fgsl_bspline_greville_abscissa = gsl_bspline_greville_abscissa(i, w%gsl_bspline_workspace)
end function fgsl_bspline_greville_abscissa
    
!-*-f90-*-
!
!  API: IEEE
!
  subroutine fgsl_ieee_fprintf_float(stream, x) 
    type(fgsl_file), intent(in) :: stream
    real(fgsl_float) :: x
    call gsl_ieee_fprintf_float(stream%gsl_file, x)
  end subroutine fgsl_ieee_fprintf_float
  subroutine fgsl_ieee_fprintf_double(stream, x) 
    type(fgsl_file), intent(in) :: stream
    real(fgsl_double) :: x
    call gsl_ieee_fprintf_double(stream%gsl_file, x)
  end subroutine fgsl_ieee_fprintf_double
  subroutine fgsl_ieee_printf_float(x) 
    real(fgsl_float) :: x
    call gsl_ieee_printf_float(x)
  end subroutine fgsl_ieee_printf_float
  subroutine fgsl_ieee_printf_double(x) 
    real(fgsl_double) :: x
    call gsl_ieee_printf_double(x)
  end subroutine fgsl_ieee_printf_double
  subroutine fgsl_ieee_env_setup()
    call gsl_ieee_env_setup()
  end subroutine fgsl_ieee_env_setup
! FIXME (maybe): some routines from the include file have not been implemented. 
end module fgsl
