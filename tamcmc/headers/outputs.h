/*
 * Outputs.h
 *
 * Contains all kind of functions
 * used to monitor and diagnose the Tempered MCMC with
 * the adaptive scheme of Atchade 2005 (MALA algorithm)
 * 
 *  Created on: 24 Mar 2016
 *      Author: obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "matrices.h"
#include "config.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

//using namespace std;


///////// The Superclass //////////
class Outputs{
	private:
		int Nbuffer; // How many samples are we keeping before saving
		int Nchains;
		int Nvars;
		int Ncons;
		int Ndata_obs; // Number of points in the data
		int Ndata_obs_limit; // For models with greater number of points... WE DO NOT ALLOW to save the models: CRASH THE CODE
        	int Nsamples; // Total number of requested samples
        	bool get_statcriteria;
		bool get_proposal_params;
		bool get_params; // if 1 save all chains parameters
		bool get_parallel_tempering_params;
		bool get_models;
		std::string file_format; // text, binary or HDF5 whenever possible. For the moment only plain text files (ascii) are possible

		bool erase_old_files; // If set to 1, any previous output file will be erased before proceeding. Turn this to 0 if you want to append an existing file
	
		std::string output_dir;
		// -- This is the output filenames for text or binary format -- //
		// --               The extension is not required            -- //
		// --           It is automatically set to .txt or .bin      -- //
		std::string params_txtbin_fileout;
		std::string proposal_txtbin_fileout;
		std::string parallel_tempering_txtbin_fileout;
		std::string model_txtbin_fileout;
		std::string stat_criteria_txtbin_fileout;
		std::string acceptance_txtbin_fileout;
		std::string file_ext; // The extension of the files: txt or bin
		// ------------------------------------------------------------ //
		// --       This is the output filename for HDF5 format      -- //
		// --      In HDF5 all outputs are save in the same file     -- //
		//std::string hdf5_fileout;

		bool do_restore_variables, do_restore_proposal, do_restore_last_index;
		std::string restore_dir;
		std::string restore_file_in; // file that contains data required in order to load an incomplete run
		std::string restore_file_out; // file that contains data required in order to save an incomplete run

	public:
		struct Buffer_proposals{
			/*
 			 * This is a buffer which contains all values that are saved
 		 	 * for the proposal law
			*/
			int counts; // How many samples are written so far? Should be iteratively updated
			int Ncopy; // How many times we have copied the buffer in the files?
			double target_acceptance;
			MatrixXd sigmas; // (i, Nvars)
			MatrixXd **mus; // (i, Nchains, Nvars)
			MatrixXd ***covarmats; // (i, Nchains, Nvars, Nvars)
			std::vector<std::string> vars_names;
			MatrixXd Pmoves;
			std::vector< std::vector<bool> > moveds; // A matrix of booleans
		};
		struct Buffer_parallel_tempering{
			/*
			 * This is a buffer which contains all values that are saved
			 * for the parallel tempering parameters
			*/
			int counts; // How many samples are written so far? Should be iteratively updated
			int Ncopy; // How many times we have copied the buffer in the files?
			VectorXd Tcoefs; // List of temperatures (of course invariant towards counts)
			std::vector<bool> attempt_mixing; // Tells us if we actually tried to mix the chains (should be the term i%dN_mixing == 0 in MALA::execute())
			VectorXi chain0s; // by definition chain1 = chain0 + 1... no need to save it then
			//MatrixXd Pswitchs;
			VectorXd Pswitchs;
			//std::vector<std::vector<bool> > switcheds; // A matrix of booleans
			std::vector<bool> switcheds; // A vector of booleans
		};
		struct Buffer_params{
			/*
			 * This is a buffer which contains all values that are saved
			 * for the models
			*/		
			int counts; // How many samples are written so far? Should be iteratively updated
			int Ncopy; // How many times we have copied the buffer in the files?
			MatrixXd **vars; // (i, Nchains, Nvars)
			VectorXd cons;
			std::vector<std::string> vars_names;
			std::vector<std::string> cons_names;
			std::vector<std::string> params_names;
			VectorXi plength;
			VectorXi relax;
		};
		struct Buffer_models{
			/*
			 * This is a buffer which contains all values that are saved
			 * for the models
			*/
		
			int counts; // How many samples are written so far? Should be iteratively updated
			int Ncopy; // How many times we have copied the buffer in the files?
			MatrixXd **models; // (i, Nchains, Ndata_obs)
		};
		struct Buffer_stat_criteria{
			/*
			 * This is a buffer which contains all values that are saved
			 * for the likelihoods, priors and posteriors
			*/
			int counts; // How many samples are written so far? Should be iteratively updated
			int Ncopy; // How many times we have copied the buffer in the files?
			MatrixXd logLikelihoods;
			MatrixXd logPriors;
			MatrixXd logPosteriors;	
		};
		struct Buffer_acceptance{
			//VectorXd xaxis;
			int Ncopy;
			double xaxis;
			VectorXd acceptance_rate;
		};
		struct Buffer_restore{
			int counts;
			int Ncopy;
			int Nsamples_sofar; // number of samples sofar written on files
			//int Naverage;
			VectorXd sigmas; // (Nvars)
			MatrixXd mus; // (Nchains, Nvars)
			MatrixXd **covarmats; // (Nchains, Nvars, Nvars)
			MatrixXd vars; // (Nchains, Nvars)
			// Averaged values of Naverage ==> Might give more robust estimates of the covariance matrix
			VectorXd sigmas_mean;
			MatrixXd mus_mean;
			MatrixXd **covarmats_mean;
			MatrixXd vars_mean;
		};

		Outputs(Config *cfg, VectorXd Tchains); // The constructor
		~Outputs(); // The destructor
		void destroy_3dMatrix(MatrixXd** m3d, const int depth);
		void destroy_4dMatrix(MatrixXd*** m4d, int depth, int Nch);

		// ---- Procedure that write the buffer in plain ASCII ----
		void write_txt_prop_params(long, bool); // writes the class Buffer_proposal into a file
		void write_txt_params(long, bool); // writes the class Buffer_params into a file
		void write_txt_parallel_temp_params(long, bool); // writes the class Buffer_parallel_tempering into a file
		void write_txt_models(long, bool);
		void write_txt_stat_criteria(long, bool);
		void write_txt_acceptance();
		void write_buffer_restore();
		// ---- Procedure that write the buffer in binary ----
		void write_bin_prop_params(long Nrest);
		void write_bin_params(long Nrest);
		void write_bin_parallel_temp_params(long Nrest);
		void write_bin_models(long Nrest);
		void write_bin_stat_criteria(long Nrest);

		// ---- Procedures for deciding when the buffer should be written ----
		void update_buffer_stat_criteria(VectorXd logLikelihood, VectorXd logPrior, VectorXd logPosterior); 	
		void update_buffer_models(MatrixXd models_in);
		//void update_buffer_acceptance(std::vector<bool>);  // There is not need for such a function because the acceptance is written when the 'proposal' buffer is full
		void update_buffer_params(MatrixXd params_in);
		void update_buffer_ptempering(bool tempted_mixing, int chain_A, long double Probaswitch, bool bool_switched);
		void update_buffer_proposals(VectorXd sigma_chains, MatrixXd mu_chains, VectorXd Pmove_chains, 
					     std::vector<bool>, MatrixXd **covarmat_chains);
		void update_buffer_restore(VectorXd sigma_chains, MatrixXd mu_chains, MatrixXd **covarmat_chains, MatrixXd params_in);

		void init_buffer_stat_criteria();
		void init_buffer_models();
		void init_buffer_acceptance();
		void init_buffer_params(VectorXd cons_in, std::vector<std::string> vars_names,
					std::vector<std::string> cons_names, std::vector<std::string> params_names, VectorXi relax, VectorXi plength);
		void init_buffer_parallel_tempering(VectorXd Tchains);
		void init_buffer_proposals(std::vector<std::string> vars_names);
		void init_buffer_restore(long init_i);

		bool file_exists (const std::string& name); // Test whether a file exist or not
		long count_accepted_vals(std::vector< std::vector<bool> > moves_chain, long m, long Nbuffer); // count the number of accepted positions
		MatrixXd reject_rate(std::vector< std::vector<bool> > moveds, long Nbuffer, long Ncopy, int Nchains, long Nsamples_init);  // Determine the rejection rate for the chains. line 0 has a x-axis. line 1 the rejection rate
		int get_Nbuffer();

		Buffer_proposals buf_proposal;
		Buffer_parallel_tempering buf_parallel_temp;
		Buffer_models buf_models;
		Buffer_stat_criteria buf_stat_crit;
		Buffer_params buf_params;
		Buffer_acceptance buf_acceptance;
		Buffer_restore buf_restore;
};



