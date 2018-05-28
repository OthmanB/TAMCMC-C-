#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>

#include "data.h" // contains the structure Data
#include "string_handler.h"
#include "io_models.h" // Contains the function that create the final vector

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

MCMC_files read_MCMC_file_MS_Global(const std::string cfg_model_file, const bool verbose);
Input_Data build_init_MS_Global(const MCMC_files inputs_MS_global, const bool verbose);
short int set_noise_params(Input_Data *Noise_in, const MatrixXd noise_s2, const VectorXd noise_params);
short int fatalerror_msg_io_MS_Global(const std::string varname, const std::string param_type, const std::string syntax_vals, const std::string example_vals);
