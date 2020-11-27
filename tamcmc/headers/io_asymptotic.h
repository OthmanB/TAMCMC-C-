#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "data.h" // contains the structure Data
//#include "string_handler.h"
#include "io_models.h"
#include "function_rot.h"
#include "io_ms_global.h"
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

Input_Data build_init_asymptotic(const MCMC_files inputs_MS_global, const bool verbose, const double resol);