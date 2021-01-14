/*
 *  Created on: 02 Mar 2016
 */
#pragma once
//#include "noise_models.h"
#include "build_lorentzian.h"
//#include "interpol.h"
//#include "function_rot.h"

using Eigen::VectorXi;
using Eigen::VectorXd;

VectorXd model_MS_Global_a1nl_etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);
VectorXd model_MS_Global_a1n_etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);
VectorXd model_MS_Global_a1l_etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);
VectorXd model_MS_Global_a1etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic_v2(const VectorXd& params, const VectorXi& params_length, const VectorXd& x); // Added on 10 Feb 2020
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic_v3(const VectorXd& params, const VectorXi& params_length, const VectorXd& x); // Added on 10 Feb 2020
VectorXd model_MS_Global_a1acta3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);
VectorXd model_MS_Global_a1etaa3_Harvey1985(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);
VectorXd model_MS_Global_a1acta3_Harvey1985(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);

VectorXd model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);
VectorXd model_MS_Global_a1etaa3_AppWidth_HarveyLike_v2(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);

VectorXd model_MS_local_basic(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);
VectorXd model_MS_local_Hnlm(const VectorXd& params, const VectorXi& params_length, const VectorXd& x); // Added on 14 Feb 2020 : Handles local fit without fitting inclination... instead used Hnlm

VectorXd model_Test_Gaussian(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);
VectorXd model_Harvey_Gaussian(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);

VectorXd model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v2(const VectorXd& params, const VectorXi& params_length, const VectorXd& x); 
VectorXd model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3(const VectorXd& params, const VectorXi& params_length, const VectorXd& x); 
VectorXd model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3(const VectorXd& params, const VectorXi& params_length, const VectorXd& x);

bool debug(const VectorXd& model, const long double Hl, const long double fl, const long double a1, const long double eta, const long double a3,
               const long double asym, const long double Wl, const long double el, const long double step, const double inclination, const VectorXd& ratios,
               const long double trunc_c, bool exit_c=true);
bool debug_solver(const VectorXd& x, const VectorXd& fl1_all, const VectorXd& fl0_all, const int el, const long double delta0l, 
                  const long double DPl, const long double alpha_g, const long double q_star, const long double sigma_p_l1);
