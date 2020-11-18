/*
 * Config_presets.h
 *
 * Header file that contains all kind of class/structures
 * used to tune the configuration using preconfigured 
 * step-by-step analysis
 * 
 *  Created on: 29 May 2016
 *      Author: obenomar
 */

#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "config.h"
#include "string_handler.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;


class Config_presets{

	private:
		std::string cfg_file_presets; // The configuration file that defines the presets (master program, equivalent to iterative_MCMC.pro in my IDL program)
		std::vector<int> Nt_learn0;
		long dN_mixing0;
	public:
		bool force_manual_config;
		std::string manual_config_file;
		std::string cfg_models_dir;
		std::string cfg_out_dir;
		std::vector<std::string> processing;
		std::vector<std::string> core_out; // specify the core of the filename for the restoring files/data files
		std::vector<std::string> core_in; // specify the core of the filename for the restoring files/data files
		VectorXi Nsamples;
		VectorXd c0;
		VectorXi restore;
		std::vector< std::vector<std::string> > table_ids;

		int current_id_ind; // index pointing to the current object id that we have to process
		int current_process_ind; // index pointing to the current processing step that we have to execute
		int current_slice_ind; // The index pointing to the current subdataset range (e.g. frequency range for io_local models)
		int Nslices; // Total number of subdataset ranges (ie, the max value that current_slice_ind can take)
		int first_slice_ind;
		int last_slice_ind;
		
		int first_id_ind; // Index of the first object id that we have to process
		int first_process_ind; // Index of the first processing step that we have to execute
		int last_id_ind; // Index of the last object id that we have to process
		int last_process_ind; // Index of the last processing step that we have to execute

		Config_presets(std::string cfg_file_in, Config *cfg0);
		void read_cfg_file_presets();
		void apply_presets(Config *cfg, const bool verbose);

		bool isdir(const std::string pathname);
		std::string shell_exec(const std::string cmd);
		bool generate_dir_tree(const std::string rootdir, const std::vector<std::string> subdirs);
		void generate_default_dirtree();
		int msg_handler(const std::string file, const std::string error_type, const std::string fct_name, const std::string arguments, const bool fatal);
		std::string format_line(const std::string str);

};
