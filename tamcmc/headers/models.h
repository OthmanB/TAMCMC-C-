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

//double lin_interpol(VectorXd x, VectorXd y, double x_int);
VectorXd model_MS_Global_a1nl_etaa3_HarveyLike(VectorXd params, VectorXi params_length, VectorXd x);
VectorXd model_MS_Global_a1n_etaa3_HarveyLike(VectorXd params, VectorXi params_length, VectorXd x);
VectorXd model_MS_Global_a1l_etaa3_HarveyLike(VectorXd params, VectorXi params_length, VectorXd x);
VectorXd model_MS_Global_a1etaa3_HarveyLike(VectorXd params, VectorXi params_length, VectorXd x);
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic(VectorXd params, VectorXi params_length, VectorXd x);
VectorXd model_MS_Global_a1acta3_HarveyLike(VectorXd params, VectorXi params_length, VectorXd x);
VectorXd model_MS_Global_a1etaa3_Harvey1985(VectorXd params, VectorXi params_length, VectorXd x);
VectorXd model_MS_Global_a1acta3_Harvey1985(VectorXd params, VectorXi params_length, VectorXd x);

VectorXd model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1(VectorXd params, VectorXi params_length, VectorXd x);
VectorXd model_MS_Global_a1etaa3_AppWidth_HarveyLike_v2(VectorXd params, VectorXi params_length, VectorXd x);

VectorXd model_MS_local_basic(VectorXd params, VectorXi params_length, VectorXd x);

VectorXd model_Test_Gaussian(VectorXd params, VectorXi params_length, VectorXd x);
VectorXd model_Harvey_Gaussian(VectorXd params, VectorXi params_length, VectorXd x);
