/*
 * priors_dictionary.h
 *
 *  Contains the headers for the function that define the priors.
 *  Currently set priors: 
 *       - Uniform
 *       - Gaussian
 *       - Jeffrey
 *       - Uniform-Gaussian
 *       - Gaussian-Uniform
 *       - Gaussian-Uniform-Gaussian
 *
 *  Created on: 09 Apr 2016
 *      Author: obenomar
 */
#pragma once
#include <iostream>
#include <iomanip>
#include <math.h>
#include <Eigen/Dense>
 
 long double logP_gaussian(long double mean, long double sigma, long double x);
 long double logP_gaussian_truncated(long double mean, long double sigma, long double x, long double bmin, long double bmax);
 long double logP_uniform(long double b_min, long double b_max, long double x);
 long double logP_uniform_abs(long double b_min, long double b_max, long double x);
 long double logP_uniform_cos(long double b_min, long double b_max, long double x);
 long double logP_multivariate_gaussian( Eigen::VectorXd mean, Eigen::MatrixXd Matrix, Eigen::VectorXd x);
 long double logP_jeffrey(long double hmin, long double hmax, long double h);
long double logP_jeffrey_abs(long double hmin, long double hmax, long double h);
 long double logP_uniform_gaussian( long double b_min, long double b_max, 
			long double sigma, long double x);
 long double logP_gaussian_uniform( long double b_min, long double b_max, 
			long double sigma, long double x);
 long double logP_gaussian_uniform_gaussian( long double b_min, long double b_max, 
			long double sigma1, long double sigma2, long double x);
