#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>

#include "data.h" // contains the structure Data
#include "string_handler.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

MCMC_files read_MCMC_file_evolved_Global(const std::string cfg_model_file, const bool verbose);
Input_Data build_init_evolved_Global(const MCMC_files inputs_MS_global, const bool verbose);
int fatalerror_msg_io_evolved_Global(const std::string varname, const std::string param_type, const std::string syntax_vals, const std::string example_vals);
