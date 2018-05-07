#include <math.h>
#include <Eigen/Dense>
#include "function_rot.h"
# include <iostream>
# include <iomanip>
#include <string>

using Eigen::VectorXd;


VectorXd build_l_mode_a1l_etaa3(const VectorXd x_l, double H_l, double fc_l, double f_s1, double f_s2, double eta, double a3, double asym, double gamma_l, const int l, VectorXd V);
VectorXd build_l_mode_a1etaa3(const VectorXd x_l, double H_l, double fc_l, double f_s, double eta, double a3, double asym, double gamma_l, const int l, VectorXd V);
VectorXd build_l_mode_a1acta3(const VectorXd x_l, double H_l, double fc_l, double f_s, double eta, double a3, double b, double alpha, double asym, double gamma_l, const int l, VectorXd V);

VectorXd optimum_lorentzian_calc_a1l_etaa3(const VectorXd x, VectorXd y, double H_l, double fc_l, double f_s1, double f_s2, double eta, double a3, double asym, double gamma_l, const int l, VectorXd V, double step, const double c);
VectorXd optimum_lorentzian_calc_a1etaa3(const VectorXd x, VectorXd y, double H_l, double fc_l, double f_s, double eta, double a3, double asym, double gamma_l, const int l, VectorXd V, double step, const double c);
VectorXd optimum_lorentzian_calc_a1acta3(const VectorXd x, VectorXd y, double H_l, double fc_l, double f_s, double eta, double a3, 
		double b, double alpha, double asym, double gamma_l, const int l, VectorXd V, double step, const double c);
