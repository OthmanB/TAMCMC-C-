/*
 * io_ms_global.cpp
 *
 * Contains all kind of functions
 * used to arrange/handle the strings for the MS_Global models
 * 
 *  Created on: 20 Jun 2016
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "data.h" // contains the structure Data
//#include "string_handler.h"
#include "io_ms_global.h"
#include "io_models.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

MCMC_files read_MCMC_file_MS_Global(const std::string cfg_model_file, const bool verbose){

	int i, out, nl, el, cpt;
	std::vector<int> pos;
	std::string line0, char0, char1, char2;
	std::vector<std::string> word, tmp;
	std::ifstream cfg_session;
        MatrixXd tmpXd;

	MCMC_files iMS_global;


    cpt=0;
    i=0;
    out=0;
    nl=200; // maximum number of lines
    el=4; // maximum degree of the modes
    //cfg_session.open(modeling.cfg_model_file.c_str());
    cfg_session.open(cfg_model_file.c_str());
    if (cfg_session.is_open()) {

	char0="#"; 
	std::getline(cfg_session, line0);
	if(verbose == 1) {std::cout << "  - Global parameters:" << std::endl;}
	iMS_global.freq_range.resize(2);
	iMS_global.els.resize(400);
	iMS_global.freqs_ref.resize(400);
	while((out < 3) && (!cfg_session.eof())){

		line0=strtrim(line0);
		char0=strtrim(line0.substr(0, 1));
		char1=line0.substr(1, 1);
		if (char0 == "#" && char1 == "K"){
			word=strsplit(line0, "= \t");
			iMS_global.ID=strtrim(strtrim(word[1]));
			if(verbose == 1) {std::cout << "           ID=" << iMS_global.ID << std::endl;}
		}
		if (char0 == "!" && char1 != "!"){
			word=strsplit(line0, " ");
			iMS_global.Dnu=str_to_dbl(word[1]);
			if(verbose == 1) {std::cout << "           Dnu =" << iMS_global.Dnu << std::endl;}
		}
		if (char0 == "!" && char1 == "!"){
			word=strsplit(line0, " ");
			iMS_global.C_l=str_to_dbl(word[1]);
			if(verbose == 1) {std::cout << "           C_l =" << iMS_global.C_l << std::endl;}
		}
		if (char0 == "*"){
			word=strsplit(line0, " ");
			iMS_global.freq_range[0]=str_to_dbl(word[1]);
			iMS_global.freq_range[1]=str_to_dbl(word[2]);
			if(verbose == 1) {std::cout << "           freq_range = [" << iMS_global.freq_range[0] << " , " << iMS_global.freq_range[1] << "]" << std::endl;}
		}
		if ((char0 != "#") && (char0 != "!") && (char0 != "*")){
			word=strsplit(line0, " ");
			if(strtrim(word[0]) == "p" || strtrim(word[0]) == "g" || strtrim(word[0]) == "co"){
				iMS_global.param_type.push_back(strtrim(word[0]));
				iMS_global.els[cpt]=str_to_int(word[1]);
                iMS_global.freqs_ref[cpt]=str_to_dbl(word[2]);
                if (word.size() >=4){
                    iMS_global.relax_freq.push_back(str_to_bool(word[3]));
                } else{
                    std::cout << "Warning: The .model file does not specify if the frequency f=" << word[2] << " is fixed/free! ==> Using default condition (free parameter) " << std::endl;
                    iMS_global.relax_freq.push_back(1);  // By default, we fit frequencies
                }
                if (word.size() >=5){
                    iMS_global.relax_gamma.push_back(str_to_bool(word[4]));
                } else{
                    std::cout << "Warning: The .model file does not specify if the Width of the mode at frequency f=" << word[2] << " is fixed/free! ==> Using default condition (free parameter) " << std::endl;
                    iMS_global.relax_gamma.push_back(1);
                }
                //std::cout << "H" << std::endl;
                if (word.size() >=6){
                    iMS_global.relax_H.push_back(str_to_bool(word[5]));
                } else{
                    std::cout << "Warning: The .model file does not specify if the Height of the mode at frequency f=" << word[2] << " is fixed/free! ==> Using default condition (free parameter) " << std::endl;
                    iMS_global.relax_H.push_back(1);
                }
                cpt=cpt+1;
			} else{
				std::cout << "Problem with the keyword that identifies the type of mode" << std::endl;
				std::cout << "The type of mode should be either 'p', 'g' or 'co'" << std::endl;
				std::cout << "THe program will exit now" << std::endl;
				exit(EXIT_FAILURE);
			}
			
		}
		if (char0 == "#"){
			out=out+1;
		}	
		i=i+1;
		char0="";
		char1="";
		std::getline(cfg_session, line0);
	}
	iMS_global.freqs_ref.conservativeResize(iMS_global.relax_freq.size());
	iMS_global.els.conservativeResize(iMS_global.relax_freq.size());
	if(verbose == 1) {
		std::cout << " - List of relaxed parameters:" << std::endl;
		std::cout << "          l    nu       relax(nu)   relax(W)  relax(H)" << std::endl;
		for(int k=0; k<iMS_global.freqs_ref.size();k++){
			std::cout << "              " << iMS_global.els[k] << "       " << iMS_global.freqs_ref[k] << "       " << iMS_global.relax_freq[k] << "  "
				  << iMS_global.relax_gamma[k] << "      " << iMS_global.relax_H[k] << std::endl;
		}
	}
    
	// -------------------------------------
	
	i=0;
	cpt=0;
	iMS_global.hyper_priors.resize(10);
	if(verbose == 1) {std::cout << " - Hyper priors:" << std::endl;}
	while ((out < 4) && !cfg_session.eof()){ // the priors, until we reach the next # symbol
			std::getline(cfg_session, line0);
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){
				iMS_global.hyper_priors[i]=str_to_dbl(line0);
				cpt=cpt+1;
			} else{
				 out=out+1;
			}
			i=i+1;
	  }
	  iMS_global.hyper_priors.conservativeResize(cpt);
	  if(verbose == 1) {
		std::cout << iMS_global.hyper_priors.transpose() << std::endl;
	  }

	i=0;
	cpt=0;
	iMS_global.eigen_params.resize(200,6);
	std::getline(cfg_session, line0);
	if(verbose == 1) {std::cout << " - Initial guesses and frequency priors:" << std::endl;}
	while ((out < 5) && !cfg_session.eof() ){ //the priors, until we reach the 5th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){
				word=strsplit(line0, " \t");
				iMS_global.eigen_params.row(i)=str_to_Xdarr(line0, " \t");
				cpt=cpt+1;
		 	} else{
				out=out+1;
			}
			i=i+1;
			std::getline(cfg_session, line0);
	}
	iMS_global.eigen_params.conservativeResize(cpt, 6);
	if(verbose == 1) {
		std::cout << iMS_global.eigen_params << std::endl;
	}
	// -------------------------------------

	i=0;
	cpt=0;
	iMS_global.noise_params.resize(10);
	while ( (out < 6) && !cfg_session.eof() ){  // the priors, until we reach the 6th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){		
				word=strsplit(line0," \t");
          			for(int j=0; j<word.size(); j++){
          			 	iMS_global.noise_params[cpt+j]=str_to_dbl(word[j]);
          			}
				cpt=cpt+word.size();
		 	} else{
				 out=out+1;
			}
			i=i+1;
			std::getline(cfg_session, line0);
	  }
    // Ensure a compatibility with the 3 Harvey profile + White noise standard
      if(cpt < 10){
        tmpXd=iMS_global.noise_params.segment(0,cpt);
        iMS_global.noise_params.setConstant(-1);
        iMS_global.noise_params.segment(10-cpt, cpt)=tmpXd;
      }
	if(verbose == 1) {
		std::cout << " - Noise inputs:" << std::endl;
		std::cout << iMS_global.noise_params.transpose() << std::endl;
	}
	// -------------------------------------	

	i=0;
	cpt=0;
	iMS_global.noise_s2.resize(10, 3);
	while ((out < 7) && !cfg_session.eof() ){ //the priors, until we reach the 7th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){
				word=strsplit(line0," \t");
				iMS_global.noise_s2.row(i)=str_to_Xdarr(line0, " \t");
				cpt=cpt+1;
		 	} else{
				out=out+1;
			}
			i=i+1;
			std::getline(cfg_session, line0);
	}
    if(i < 11){
        tmpXd=iMS_global.noise_s2.block(0,0,cpt,iMS_global.noise_s2.cols());
        iMS_global.noise_s2.setConstant(-1);
        iMS_global.noise_s2.block(10-cpt, 0, tmpXd.rows(), tmpXd.cols())=tmpXd;
    }
    if(verbose == 1) {
		std::cout << " - Noise information extracted during step s2:" << std::endl;
		std::cout << iMS_global.noise_s2 << std::endl;
    }
    //	exit(EXIT_SUCCESS);
	// -------------------------------------	
	
	// -------- The common parameters follow ------
	VectorXd a;
	i=0;
	cpt=0;
	iMS_global.modes_common.resize(20,5);
	iMS_global.modes_common.setConstant(-9999); // up to 10 variables and 4 prior parameters
	while ( (out < 9) && !cfg_session.eof() ){ // the initial values for the common parameters + priors, until we reach the 9th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			word=strsplit(line0," \t");
			if (char0 != "#"){
				word=strsplit(line0," \t");
				iMS_global.common_names.push_back(strtrim(word[0]));
				iMS_global.common_names_priors.push_back(strtrim(word[1]));
 				word.erase(word.begin()); // erase the slot containing the keyword name
				word.erase(word.begin()); // erase the slot containing the type of prior/switch
				a=arrstr_to_Xdarrdbl(word);
				
				for(int k=0; k<a.size();k++){
					iMS_global.modes_common(cpt, k)=a[k];
				}
				cpt=cpt+1;
		 	} else{
				out=out+1;
			}
			i=i+1;
			std::getline(cfg_session, line0);
	}

	iMS_global.modes_common.conservativeResize(iMS_global.common_names.size(), 5);

	if(verbose == 1) {
		std::cout << " - Common parameters for modes:" << std::endl;
		for(i=0; i<iMS_global.common_names.size();i++){
			std::cout << "  " << iMS_global.common_names[i] << "  ";
			std::cout << "  " << iMS_global.common_names_priors[i] << "  ";
			std::cout << "  " << iMS_global.modes_common.row(i) << std::endl;
		}
	}
    
	// -------------------------------------
   cfg_session.close();
   } else {
   		std::cout << "Unable to open the file: " << cfg_model_file << std::endl;
   		std::cout << "Check that the file exist and that the path is correct" << std::endl;
   		std::cout << "Cannot proceed" << std::endl;
   		std::cout << "The program will exit now" << std::endl;
   		exit(EXIT_FAILURE);
   }
   return iMS_global;
}

Input_Data build_init_MS_Global(const MCMC_files inputs_MS_global, const bool verbose, const double resol){

	const long double pi = 3.141592653589793238L;
	const double G=6.667e-8;
	const double Teff_sun= 5777; 
	const double Dnu_sun=135.1;
	const double numax_sun=3150.;
	const double R_sun=6.96342e5; //in km
	const double M_sun=1.98855e30; //in kg
	const double rho_sun=M_sun*1e3/(4*pi*std::pow(R_sun*1e5,3)/3); //in g.cm-3
	const int Nmax_prior_params=4; // The maximum number of parameters for the priors. Should be 4 in all my code

	double rho=pow(inputs_MS_global.Dnu/Dnu_sun,2.) * rho_sun;
	double Dnl=0.75, trunc_c=-1;

	// All Default booleans
	bool do_a11_eq_a12=1, do_avg_a1n=1, do_amp=0;
	bool bool_a1sini=0, bool_a1cosi=0;

	int lmax, en, Ntot, p0;
	double tol=1e-2, tmp;
	VectorXi pos_el, pos_relax0, els_eigen, Nf_el(4), plength;
	VectorXd extra_priors, tmpXd;
	std::vector<int> pos_relax;
	std::vector<double> f_inputs, h_inputs, w_inputs, f_priors_min, f_priors_max, f_el;;
	std::vector<bool> f_relax, h_relax, w_relax; 
	std::vector<int> rf_el, rw_el, rh_el;
	
	std::string tmpstr_h, tmpstr;
	
	Input_Data Snlm_in, Vis_in, Inc_in, Noise_in, freq_in, height_in, width_in; // This is by block, each category of parameters		
	Input_Data all_in; // The final structure of parameters, using the standards of my code
	IO_models io_calls; // function dictionary that is used to initialise, create and add parameters to the Input_Data structure

	// Flatening and ordering of all the inputs/relax variables
	lmax=inputs_MS_global.els.maxCoeff();

	// -- Initialisation of structures --
    // --- Look for common instruction That must be run before the setup ---------
	all_in.model_fullname=" "; // Default is an empty string
	//all_in.prior_fullname="prior_MS_Global"; // Default set of prior
    for(int i=0; i<inputs_MS_global.common_names.size(); i++){
        if(inputs_MS_global.common_names[i] == "model_fullname" ){ // This defines if we assume S11=S22 or not (the executed model will be different)
        	all_in.model_fullname=inputs_MS_global.common_names_priors[i];
            if(all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike_Classic"){
            	//Previously corresponding to average_a1nl     bool    1    1 
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            }
            if(all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike"){
            	//Previously corresponding to average_a1nl     bool    1    1 
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            }
            if(all_in.model_fullname == "model_MS_Global_a1etaa3_Harvey1985"){
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;            	
            }
            if(all_in.model_fullname == "model_MS_Global_a1acta3_HarveyLike"){
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            }
            if(all_in.model_fullname == "model_MS_Global_a1acta3_Harvey1985"){
         		do_a11_eq_a12=1;
           		do_avg_a1n=1;          	
           	}
           	if(all_in.model_fullname == "model_MS_Global_a1l_etaa3_HarveyLike"){
           		//Previously corresponding to average_a1nl     bool    0    1 
           	    do_a11_eq_a12=0;
           		do_avg_a1n=1;
           	}
        	if(all_in.model_fullname == "model_MS_Global_a1n_etaa3_HarveyLike"){
           		//Previously corresponding to average_a1nl     bool    1    0            	    
           		do_a11_eq_a12=1;
        		do_avg_a1n=0;
            }
        	if(all_in.model_fullname == "model_MS_Global_a1nl_etaa3_HarveyLike"){
        		//Previously corresponding to average_a1nl     bool    0    0
        		do_a11_eq_a12=0;            		
        		do_avg_a1n=0;
            }
        }
        if(inputs_MS_global.common_names[i] == "fit_squareAmplitude_instead_Height" ){ 
        		do_amp=1;
        	if(inputs_MS_global.common_names_priors[i] != "bool"){
				fatalerror_msg_io_MS_Global("fit_squareAmplitude_instead_Height", "bool", "[0/1]  ", "1 " );
			} else{
				do_amp=inputs_MS_global.modes_common(i,0);
				std::cout << "Using do_amp = " << do_amp << std::endl;
			}

        }
 
    }
    if(all_in.model_fullname == " "){
    	std::cout << "Model name empty. Cannot proceed. Check that the .model file contains the model_fullname variable." << std::endl;
    	exit(EXIT_FAILURE);
    }
 	
 	io_calls.initialise_param(&Vis_in, lmax, Nmax_prior_params, -1, -1);
	io_calls.initialise_param(&Inc_in, 1, Nmax_prior_params, -1, -1);
	
	// -----------------------------------------------------------------
	// ------------ Handling Frequencies/Widths/Heights ----------------
	// -----------------------------------------------------------------
	Nf_el.setZero();
	for(int el=0; el<lmax+1; el++){
		els_eigen.resize(inputs_MS_global.eigen_params.col(0).size());
		for(int i=0; i<els_eigen.size();i++){
			els_eigen[i]=inputs_MS_global.eigen_params(i,0);
		}
		// --- sublist of relax taken from the first frequency list at the given el----
		pos_relax0=where_int(inputs_MS_global.els, el); // Get all the positions for the modes of degree el in the relax tab
		for(int i=0; i<pos_relax0.size(); i++){ // fill temporary variables for the given el 
			f_el.push_back(inputs_MS_global.freqs_ref[pos_relax0[i]]);
			rf_el.push_back(inputs_MS_global.relax_freq[pos_relax0[i]]);
			rw_el.push_back(inputs_MS_global.relax_gamma[pos_relax0[i]]);
			rh_el.push_back(inputs_MS_global.relax_H[pos_relax0[i]]);
		}
		
		pos_el=where_int(els_eigen, el); // find where we have modes of degree el in the eigen vectors
		Nf_el[el]=pos_el.size(); // Get the Number of frequencies of a given el... this will be used in plength
		for(int en=0; en<pos_el.size(); en++){
			f_inputs.push_back(inputs_MS_global.eigen_params(pos_el[en],1));
			f_priors_min.push_back(inputs_MS_global.eigen_params(pos_el[en],2));
			f_priors_max.push_back(inputs_MS_global.eigen_params(pos_el[en],3));
			if(el ==0){
				w_inputs.push_back(inputs_MS_global.eigen_params(pos_el[en],4)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE WIDTH ARE INTERPOLATED USING l=0
				h_inputs.push_back(inputs_MS_global.eigen_params(pos_el[en],5)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE HEIGHT ARE SCALED USING VISBILITIES				
			}

			pos_relax=where_dbl(f_el, inputs_MS_global.eigen_params(pos_el[en],1), tol); // Get the (unique) position in the relax tab for the seeked frequency

			if(pos_relax[0] != -1 && pos_relax.size() == 1){ // We impose uniqueness of the solution, within the tolerance
				f_relax.push_back(rf_el[pos_relax[0]]);
				if(el ==0){
					w_relax.push_back(rw_el[pos_relax[0]]);
					h_relax.push_back(rh_el[pos_relax[0]]);
				}
			} else{
				std::cout << "Error when preparing the vector of parameters using the MCMC file" << std::endl;
				std::cout << "The uniqueness of the frequency " << inputs_MS_global.eigen_params(pos_el[en],1) << " is not respected" << std::endl;
				std::cout << "     - pos_relax.size() =" << pos_relax.size() << std::endl;
				std::cout << "     - pos_relax[0] =" << pos_relax[0] << std::endl;
				if(pos_relax.size() != 1){
					std::cout << "Check that you do not have duplicated values in your MCMC file " << std::endl;

				}
				if(pos_relax[0] == -1){
					std::cout << "Check that you do not the same values in the relax list and in the eigen list of your MCMC file " << std::endl;
				}
				std::cout << "The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		// Empty the temporary variables for the relax and the freqs_ref at a given el
		f_el.resize(0);
		rf_el.resize(0);
		rw_el.resize(0);
		rh_el.resize(0);
	}
	

	// ------------------------------------------------------------------------------------------
	// ------------------------------- Handling the Common parameters ---------------------------
	// ------------------------------------------------------------------------------------------
	if( do_amp){
		std::cout << "   ===> Requested to fit squared amplitudes instead of Height... Converting height inputs into A_squared = pi*Height*Width..." << std::endl;
		tmpstr_h="Amplitude_l";
		for(int i=0; i<h_inputs.size(); i++){
			h_inputs[i]=pi*w_inputs[i]*h_inputs[i]; 
		}
	} else{
		tmpstr_h="Height_l";
	}
    // Set default value of priors for Height Width and frequency
	io_calls.initialise_param(&height_in, h_relax.size(), Nmax_prior_params, -1, -1);
	io_calls.initialise_param(&width_in, w_relax.size(), Nmax_prior_params, -1, -1); 
	io_calls.initialise_param(&freq_in, f_relax.size(), Nmax_prior_params, -1, -1); 

	tmpXd.resize(4);
	tmpXd << 1, 100000., -9999., -9999.; // default hmin and hmax for the Jeffreys prior
	for(int i=0; i<h_inputs.size(); i++){
		if(h_relax[i]){
			io_calls.fill_param(&height_in, tmpstr_h, "Jeffreys", h_inputs[i], tmpXd, i, 0);	
		} else{
			io_calls.fill_param(&height_in, tmpstr_h, "Fix", h_inputs[i], tmpXd, i, 0);			
		}
	}
	tmpXd << 0.1, 50., -9999., -9999.; // default hmin and hmax for the Jeffreys prior
	for(int i=0; i<w_inputs.size(); i++){
		if(w_relax[i]){ 
			io_calls.fill_param(&width_in, "Width_l", "Jeffreys", w_inputs[i], tmpXd, i, 0);	
		} else{
			io_calls.fill_param(&width_in, "Width_l", "Fix", w_inputs[i], tmpXd, i, 0);			
		}
	}

	// --- Default setup for frequencies ---
	for(int i=0; i<f_inputs.size(); i++){
		if(f_relax[i]){
			tmpXd << f_priors_min[i], f_priors_max[i], 0.01*inputs_MS_global.Dnu, 0.01*inputs_MS_global.Dnu; // default parameters for a GUG prior on frequencies
			io_calls.fill_param(&freq_in, "Frequency_l", "GUG", f_inputs[i], tmpXd, i, 0);	
		} else{
			io_calls.fill_param(&freq_in, "Frequency_l", "Fix", f_inputs[i], tmpXd, i, 0);			
		}
	}

	// ----- Switch between the models that handle averaging over n,l or both -----
   if(do_a11_eq_a12 == 1 && do_avg_a1n == 1){
        io_calls.initialise_param(&Snlm_in, 6, Nmax_prior_params, -1, -1);
    }  
    if(do_a11_eq_a12 == 0 && do_avg_a1n == 1){
        io_calls.initialise_param(&Snlm_in, 7, Nmax_prior_params, -1, -1);
    }
    if(do_a11_eq_a12 == 1 && do_avg_a1n == 0){
		if (Nf_el[1] == Nf_el[2]){
			io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1], Nmax_prior_params, -1, -1);
		} else {
			std::cout << "When considering a11=a22"<< std::endl;
			std::cout <<" You must have as many l=1 than l=2" << std::endl;
			std::cout <<" Here we have: Nf(l=1) = " << Nf_el[1] << " and Nf(l=2) = " << Nf_el[2] << std::endl;
			std::cout <<" Check the .model file" << std::endl;
			std::cout <<" The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
    }
    if(do_a11_eq_a12 == 0 && do_avg_a1n == 0){
    	io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1]+Nf_el[2], Nmax_prior_params, -1, -1);
    }
    
	// -------------- Set Extra_priors ----------------	
	extra_priors.resize(2);
	extra_priors[0]=1; // By default, we apply a smoothness condition
	extra_priors[1]=2.; // By default, the smoothness coeficient is 2 microHz
	// ------------------------------------------------
	
	for(int i=0; i<inputs_MS_global.common_names.size(); i++){
		// --- Common parameters than can be run during setup ---
		if(inputs_MS_global.common_names[i] == "freq_smoothness" || inputs_MS_global.common_names[i] == "Freq_smoothness"){
			if(inputs_MS_global.common_names_priors[i] != "bool"){
				fatalerror_msg_io_MS_Global("freq_smoothness", "bool", "[0/1]     [Smoothness coeficient in microHz]", "1     2.0" );
			} else{
				extra_priors[0]=inputs_MS_global.modes_common(i,0);
				extra_priors[1]=inputs_MS_global.modes_common(i,1);
			}
		}
		if(inputs_MS_global.common_names[i] == "trunc_c"){
			if(inputs_MS_global.common_names_priors[i] != "Fix"){
				fatalerror_msg_io_MS_Global("trunc_c", "Fix", "[Truncation parameter]", "20" );
			} else{
				trunc_c=inputs_MS_global.modes_common(i,0); // This value is added to the input vector at the end of this function.
			}
		}	

		// --- Frequencies ---
		if(inputs_MS_global.common_names[i] == "Frequency" || inputs_MS_global.common_names[i] == "frequency"){ 
			if(inputs_MS_global.common_names_priors[i] == "GUG" || inputs_MS_global.common_names_priors[i] == "Uniform"){
				for(int p0=0; p0<f_inputs.size(); p0++){
					if(f_relax[p0]){
						if(inputs_MS_global.common_names_priors[i] == "GUG"){
							tmpXd << f_priors_min[p0], f_priors_max[p0], inputs_MS_global.modes_common(i,3), inputs_MS_global.modes_common(i,4);
						} else{
							tmpXd << f_priors_min[p0], f_priors_max[p0], -9999, -9999; 						
						}
						io_calls.fill_param(&freq_in, "Frequency_l", inputs_MS_global.common_names_priors[i], f_inputs[p0], tmpXd, p0, 0);	
					} else{
						io_calls.fill_param(&freq_in, "Frequency_l", "Fix", f_inputs[p0], tmpXd, p0, 0);			
					}
				}
			} else{
				fatalerror_msg_io_MS_Global(inputs_MS_global.common_names[i], "GUG or Uniform", "", "" );
			}
		}

		// --- Height or Amplitude ---
		if(inputs_MS_global.common_names[i] == "height" || inputs_MS_global.common_names[i] == "Height" || 
		   inputs_MS_global.common_names[i] == "amplitude" || inputs_MS_global.common_names[i] == "Amplitude"){

			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
					fatalerror_msg_io_MS_Global(inputs_MS_global.common_names[i], "Fix_Auto", "", "" );
				}
			for(p0=0; p0<h_inputs.size(); p0++){
				if(h_relax[p0]){
					io_calls.fill_param(&height_in, tmpstr_h,  inputs_MS_global.common_names_priors[i], h_inputs[p0],  inputs_MS_global.modes_common.row(i), p0, 1);	
				} else{
					io_calls.fill_param(&height_in, tmpstr_h,  "Fix", h_inputs[p0],  inputs_MS_global.modes_common.row(i), p0, 1);		
				}
			}
		}
		// -- Mode Width ---
		if(inputs_MS_global.common_names[i] == "width" || inputs_MS_global.common_names[i] == "Width"){
				if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
					//fatalerror_msg_io_MS_Global(inputs_MS_global.common_names[i], "Fix_Auto", "", "" );
                    tmpstr="Jeffreys";
                    tmpXd << resol, inputs_MS_global.Dnu/3., -9999., -9999.;
                    std::cout << "Fix_Auto requested for Widths... For all free Widths, the prior will be with this syntax:" << std::endl;
                    std::cout << "          " << std::left << std::setw(15) << tmpstr << " [Spectrum Resolution]   [Deltanu / 3]   -9999    -9999" << std::endl;
                    std::cout << "          " << "Resolution: " << resol << std::endl;
                } else{
                    tmpstr=inputs_MS_global.common_names_priors[i];
                    tmpXd=inputs_MS_global.modes_common.row(i);
                }
			for(p0=0; p0<w_inputs.size(); p0++){
				if(w_relax[p0]){
					io_calls.fill_param(&width_in, "Width_l", tmpstr, w_inputs[p0],  tmpXd, p0, 0);
				} else{
					io_calls.fill_param(&width_in, "Width_l",  "Fix", w_inputs[p0],  inputs_MS_global.modes_common.row(i), p0, 1);		
				}
			}
		}			
		// --- Splittings and asymetry ---
		if(inputs_MS_global.common_names[i] == "splitting_a1" || inputs_MS_global.common_names[i] == "Splitting_a1"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("splitting_a1", "Fix_Auto", "", "" );
			}
			p0=0;
			io_calls.fill_param(&Snlm_in, "Splitting_a1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
        	if(do_a11_eq_a12 == 0 && do_avg_a1n == 1){ // Case a1(l), we put the same prior for a1(2) than a1(1)
		    	std::cout << "do_a11_eq_a12 == 0 && do_avg_a1n == 1" << std::endl;
		    	p0=6;
        	    io_calls.fill_param(&Snlm_in, Snlm_in.inputs_names[0], Snlm_in.priors_names[0], Snlm_in.inputs[0], inputs_MS_global.modes_common.row(i), p0, 1);
        	}
        	if(do_a11_eq_a12 == 1 && do_avg_a1n == 0){ // Case a1(n), we put the same prior for a1(n)
				std::cout << "do_a11_eq_a12 == 1 && do_avg_a1n == 0" << std::endl;
				for(int kk=0; kk<Nf_el[1]; kk++){
					p0=6+kk;
					io_calls.fill_param(&Snlm_in, Snlm_in.inputs_names[0], Snlm_in.priors_names[0], Snlm_in.inputs[0], inputs_MS_global.modes_common.row(i), p0, 1); 
				}
				p0=0;
				io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, inputs_MS_global.modes_common.row(i), p0, 1);
        	}
        	if(do_a11_eq_a12 == 0 && do_avg_a1n == 0){ // Case a1(n,l), we put the same prior for a1(n,l)
				std::cout << "do_a11_eq_a12 == 0 && do_avg_a1n == 0" << std::endl;
				for(int kk=0; kk<Nf_el[1]+Nf_el[2]; kk++){
					p0=6+kk;
					io_calls.fill_param(&Snlm_in, Snlm_in.inputs_names[0], Snlm_in.priors_names[0], Snlm_in.inputs[0], inputs_MS_global.modes_common.row(i), p0, 1); 
				}
				p0=0;
				io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, inputs_MS_global.modes_common.row(i), p0, 1); // Erase values of the default Splitting_a1 block
        	}
		}
		if(inputs_MS_global.common_names[i] == "asphericity_eta"|| inputs_MS_global.common_names[i] == "Asphericity_eta"){  
			Snlm_in.inputs_names[1]="Asphericity_eta";
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){ // Case where the centrifugal force is added but FIXED.
				Snlm_in.priors_names[1]="Fix"; //In all cases the centrifugal force is fixed
				Snlm_in.relax[1]=0;
				if(inputs_MS_global.modes_common(i,0) == 1){
  					if (Snlm_in.inputs[0] != -9999){ 
						//eta=(4./3.)*!pi*Dnl*(a1_init*1d-6)^2/(rho*G)
						Snlm_in.inputs[1]=(4./3.)*pi*Dnl*pow(Snlm_in.inputs[0]*1e-6,2.)/(rho*G);
						std::cout << " -------------" << std::endl;
						std::cout << "Asphericity_eta given with Fix_Auto ==> Setting the asphericity to the centrifugal force" << std::endl;
						std::cout << "      eta=" << Snlm_in.inputs[1] << std::endl;
						std::cout << " -------------" << std::endl;
					} else{
						std::cout << "Warning: the keyword 'asphericity_eta' must appear after the keyword splitting_a1" << std::endl;
						std::cout << "         This because the splitting_a1 is used to define the initial value of asphericity" << std::endl;
						std::cout << "         Edit the .MCMC file accordingly" << std::endl;
						std::cout << "The program will exit now" << std::endl;
						exit(EXIT_FAILURE);
					}
				} else{
					Snlm_in.inputs[1]=0;
				}
			} else{ // Case where the centrifugal force is user-defined: Could be free or fixed
				p0=1;
				io_calls.fill_param(&Snlm_in, "Asphericity_eta", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			}
		}

		if(inputs_MS_global.common_names[i] == "splitting_a3" || inputs_MS_global.common_names[i] == "Splitting_a3"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("splitting_a3", "Fix_Auto", "", "" );
			}
			p0=2;
			io_calls.fill_param(&Snlm_in, "Splitting_a3", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
		}

		if(inputs_MS_global.common_names[i] == "asymetry" || inputs_MS_global.common_names[i] == "Asymetry"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("asymetry", "Fix_Auto", "", "" );
			}
			p0=5;
			io_calls.fill_param(&Snlm_in, "Lorentzian_asymetry", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
		}
		// --- Dealing with visibilities ---
		if(inputs_MS_global.common_names[i] == "visibility_l1" || inputs_MS_global.common_names[i] == "Visibility_l1"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("visibility_l1", "Fix_Auto", "", "" );
			}
			if(lmax >= 1){
				p0=0;
				io_calls.fill_param(&Vis_in, "Visibility_l1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
			} else{
				std::cout << "Warning: lmax=" << lmax << " but keyword 'visibility_l1' detected" << std::endl;
				std::cout << "         This visibilitiy input will be ignored" << std::endl;
				std::cout << "         Proceeding..." << std::endl;
			}
		}
		if(inputs_MS_global.common_names[i] == "visibility_l2" || inputs_MS_global.common_names[i] == "Visibility_l2"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("visibility_l2", "Fix_Auto", "", "" );
			}
			if(lmax >= 2){
				p0=1;
				io_calls.fill_param(&Vis_in, "Visibility_l2", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
			} else{
				std::cout << "Warning: lmax=" << lmax << " but keyword 'visibility_l2' detected" << std::endl;
				std::cout << "         This visibilitiy input will be ignored" << std::endl;
				std::cout << "         Proceeding..." << std::endl;
			}
		}
		if(inputs_MS_global.common_names[i] == "visibility_l3" || inputs_MS_global.common_names[i] == "Visibility_l3"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("visibility_l3", "Fix_Auto", "", "" );
			}
			if(lmax >= 3){
				p0=2;
				io_calls.fill_param(&Vis_in, "Visibility_l3", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
			} else{
				std::cout << "Warning: lmax=" << lmax << " but keyword 'visibility_l3' detected" << std::endl;
				std::cout << "         This visibilitiy input will be ignored" << std::endl;
				std::cout << "         Proceeding..." << std::endl;
			}
		}
		// --- Finally the inclination ---
		if(inputs_MS_global.common_names[i] == "inclination" || inputs_MS_global.common_names[i] == "Inclination"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("inclination", "Fix_Auto", "", "" );
			}
			if( inputs_MS_global.modes_common(i,0) >= 90){ // Avoid some nan when computing the prior
				    tmp=89.99999; 
			} else{
				tmp=inputs_MS_global.modes_common(i,0);
			}
			p0=0;
			io_calls.fill_param(&Inc_in, "Inclination", inputs_MS_global.common_names_priors[i], tmp, inputs_MS_global.modes_common.row(i), p0, 1);
		}
		
		if(inputs_MS_global.common_names[i] == "sqrt(splitting_a1).cosi"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("sqrt(splitting_a1).cosi", "Fix_Auto", "", "" );
			}
			p0=3;
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).cosi", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);

            if(do_a11_eq_a12 == 0 || do_avg_a1n == 0){ // In that case, we put the same prior for a1(2) than a1(1)
				std::cout << "Warning: do_a11_eq_a12=0  and/or do_avg_a1n=0 (models *_a1n* or *_a1l*) was requested but is not available when fitting sqrt(a1).cosi and sqrt(a1).sini" << std::endl;
				std::cout << "         You must modify the code accordingly if you want to implement a11 and a12 in that scenario" << std::endl;
				std::cout << "         The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
            }
            bool_a1cosi=1;
		}
		
		if(inputs_MS_global.common_names[i] == "sqrt(splitting_a1).sini"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("sqrt(splitting_a1).sini", "Fix_Auto", "", "" );
			}
			p0=4;
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).sini", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
			
            if(do_a11_eq_a12 == 0 || do_avg_a1n == 0){ // In that case, we put the same prior for a1(2) than a1(1)
				std::cout << "Warning: do_a11_eq_a12=0 and/or do_avg_a1n=0 (models *_a1n* or *_a1l*) was requested but is not available when fitting sqrt(a1).cosi and sqrt(a1).sini" << std::endl;
				std::cout << "         You must modify the code accordingly if you want to implement a11 and a12 in that scenario" << std::endl;
				std::cout << "         The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
            }
            bool_a1sini=1;
		}
	}

if(bool_a1cosi != bool_a1sini){ // Case when one of the projected splitting quantities is missing ==> Problem
	std::cout << "Warning: Both 'sqrt(splitting_a1).sini' and 'sqrt(splitting_a1).cosi' keywords must appear" << std::endl;
	std::cout << "         It is forbidden to use only one of them" << std::endl;
	std::cout << "         Edit the .MODEL file accordingly" << std::endl;
	std::cout << "The program will exit now" << std::endl;
	exit(EXIT_FAILURE);
}
if((bool_a1cosi == 0) && (bool_a1sini == 0)){ // Case where Inclination and Splitting_a1 are supposed to be used 
	if( (all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike") || (all_in.model_fullname == "model_MS_Global_a1etaa3_Harvey1985")){
		std::cout << "Warning: Splitting_a1 and Inclination keywords detected while requested model is model_MS_Global_a1etaa3_HarveyLike... ADAPTING THE VARIABLE FOR ALLOWING THE FIT TO WORK" << std::endl;
		
		std::cout << "         Replacing variables splitting_a1 and inclination by sqrt(splitting_a1).cos(i) and sqrt(splitting_a1).sin(i)..." << std::endl;

		if( Inc_in.priors_names[0]=="Fix" && Snlm_in.priors_names[0]=="Fix"){
			std::cout << "         - Inclination and splitting_a1 were 'Fixed' ==> Fixing sqrt(splitting_a1).cos(i) and sqrt(splitting_a1).sin(i)" << std::endl; 
			
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).cosi", "Fix", sqrt(Snlm_in.inputs[0])*cos(Inc_in.inputs[0]*pi/180.), Snlm_in.priors.col(0), 3, 0);
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).sini", "Fix", sqrt(Snlm_in.inputs[0])*sin(Inc_in.inputs[0]*pi/180.), Snlm_in.priors.col(0), 4, 0);
		} else{
			std::cout << "	       - Inclination and/or Splitting_a1 requested as a free parameter ==> Freeing both sqrt(splitting_a1).cos(i) and sqrt(splitting_a1).sin(i)" << std::endl;
			std::cout << "         - Setting the prior to the one requested for the Splitting_a1 assuming non-informative prior on inclination: " << Snlm_in.priors_names[0] << std::endl;
			std::cout << "         - Value of the prior are set for cosi=sini=1, ie the max is set by sqrt(splitting_a1)" << std::endl;
			std::cout << "           Beware that this may not be what you want!" << std::endl;

			Snlm_in.priors(1,0)=std::sqrt(Snlm_in.priors(1,0)); // Convert the max boundary of splitting_a1 into the max boundary of sqrt(splitting_a1)
			
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).cosi", Snlm_in.priors_names[0], sqrt(Snlm_in.inputs[0])*cos(Inc_in.inputs[0]*pi/180.), Snlm_in.priors.col(0), 3, 0);
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).sini", Snlm_in.priors_names[0], sqrt(Snlm_in.inputs[0])*sin(Inc_in.inputs[0]*pi/180.), Snlm_in.priors.col(0), 4, 0);
		}

		if(Snlm_in.inputs[3] < 1e-2){ // Avoid issues with the edge of uniform priors
			Snlm_in.inputs[3]=1e-2;
		}
		if(Snlm_in.inputs[4] < 1e-2){ // Avoid issues with the edge of uniform priors
			Snlm_in.inputs[4]=1e-2;
		}

        std::cout << "            Slots for the variables Inclination and Splitting are forced to be FIXED" << std::endl;
        std::cout << "            Be aware that may lead to unwished results if you are not careful about the used model" << std::endl;
        io_calls.fill_param(&Inc_in, "Empty", "Fix", 0, Inc_in.priors.col(0), 0, 1); // Note that inputs_MS_global.modes_common.row(0) is not used... just dummy
        io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, Snlm_in.priors.col(0), 0, 1); // Note that inputs_MS_global.modes_common.row(0) is not used... just dummy
	}
}
if((bool_a1cosi == 1) && (bool_a1sini ==1)){
	if (all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike_Classic"){
		std::cout << "Warning: We cannot use " << all_in.model_fullname << " with variables sqrt(splitting_a1).cosi and sqrt(splitting_a1).sini..." << std::endl;
		std::cout << "         Use Inclination and Splitting_a1 instead for the model this model"  << std::endl;
		std::cout << "	       Alternatively, you can use 'model_MS_Global_a1etaa3_HarveyLike with sqrt(splitting_a1).cosi and sqrt(splitting_a1).sini keywords'" <<std::endl;
		std::cout << "         The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
    std::cout << "    NOTICE: sqrt(splitting_a1).sini and sqrt(splitting_a1).cosi superseed inclination and splitting" << std::endl;
    std::cout << "            Slots for the variables Inclination and Splitting are therefore forced to be FIXED" << std::endl;
    std::cout << "            Be aware that may lead to unwished results if you are not careful about the used model" << std::endl;
   	
   	io_calls.fill_param(&Inc_in, "Empty", "Fix", 0, Inc_in.priors.col(0), 0, 1); // Note that inputs_MS_global.modes_common.row(0) is not used... just dummy   	
   	io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, Snlm_in.priors.col(0), 0, 1); // "Splitting_a1" default values are erased
}
	// ----------------------------------------------------
	// ---------------- Handling noise --------------------
	// ----------------------------------------------------
	io_calls.initialise_param(&Noise_in, 10, Nmax_prior_params, -1, -1);
	set_noise_params(&Noise_in, inputs_MS_global.noise_s2, inputs_MS_global.noise_params); 	// Set defaults

	// -------------------------------------------------------------------------------------
	// ------------- Sticking everything together in a Input_Data structure --------------
	// -------------------------------------------------------------------------------------
	
	plength.resize(11);
	plength[0]=h_inputs.size(); plength[1]=lmax                  ; plength[2]=Nf_el[0];
	plength[3]=Nf_el[1]       ; plength[4]=Nf_el[2]		   ; plength[5]=Nf_el[3];
	plength[6]=Snlm_in.inputs.size(); plength[7]=w_inputs.size() ; plength[8]=Noise_in.inputs.size(); 
	plength[9]=Inc_in.inputs.size();
	plength[10]=2; // This is trunc_c and do_amp;
	
	io_calls.initialise_param(&all_in, plength.sum(), Nmax_prior_params, plength, extra_priors);
	
	// --- Put the Height or Amplitudes---
	p0=0;
	io_calls.add_param(&all_in, &height_in, p0); // height_in may contain either height or amplitudes depending on specified keywords
	// --- Put the Visibilities ---
	p0=all_in.plength[0];
	io_calls.add_param(&all_in, &Vis_in, p0);
	
	// --- Put the Frequencies ---
	p0=all_in.plength[0] + all_in.plength[1];	
	io_calls.add_param(&all_in, &freq_in, p0);

	// --- Put the Snlm (splittings and asymetry) ---
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5];
	io_calls.add_param(&all_in, &Snlm_in, p0);
	
	// --- Put the Width ---
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] +  all_in.plength[3]  + all_in.plength[4] + all_in.plength[5] + all_in.plength[6];
	io_calls.add_param(&all_in, &width_in, p0);
		
	// --- Put the Noise ---
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] +  all_in.plength[3]  + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7];
	io_calls.add_param(&all_in, &Noise_in, p0);
	
	// --- Put the Inclination ---
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8];
	io_calls.add_param(&all_in, &Inc_in, p0);
	
	// --- Add trunc_c that controls the truncation of the Lorentzian ---
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8] + all_in.plength[9];
	io_calls.fill_param(&all_in, "Truncation parameter", "Fix", trunc_c, inputs_MS_global.modes_common.row(0), p0, 1);
	if (all_in.inputs[p0] <= 0){
		std::cout << "Warning: trunc_c <= 0. This is forbidden. Setting default to 10000. (No truncation)" << std::endl;
		all_in.inputs[p0]=10000.; // In case of a non-sense value for c, we use Full-Lorentzian as default
	}
	// -- Add the Amplitude switch --
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8] + all_in.plength[9] + 1;
	io_calls.fill_param(&all_in, "Switch for fit of Amplitudes or Heights", "Fix", do_amp, inputs_MS_global.modes_common.row(0), p0,1);
		
			
	if(verbose == 1){
		std::cout << " ----------------- Configuration summary -------------------" << std::endl;
		std::cout << "Model Name = " << all_in.model_fullname << std::endl;
		if(all_in.extra_priors[0] == 0){ 
			std::cout << "   freq_smoothness is set to 0 ==> NO smoothness condition on frequencies" << std::endl;
		} else{
			std::cout << "   freq_smoothness is set to 1 ==> APPLIES a smoothness condition on frequencies" << std::endl;
			std::cout << "   smoothness coeficient as specified by the user (of defined by default): " << all_in.extra_priors[1] << " microHz" << std::endl;
		}
		std::cout << " -----------------------------------------------------------" << std::endl;
		std::cout << " ---------- Configuration of the input vectors -------------" << std::endl;
		std::cout << "    Table of inputs " << std::endl;
		io_calls.show_param(all_in, 1);
		std::cout << " -----------------------------------------------------------" << std::endl;

		std::cout << "    The fit includes:" << std::endl;
		std::cout << "          - " << all_in.plength[0] << "  Heights parameters "<< std::endl; 
		std::cout << "          - " << all_in.plength[1] << "  Visibility parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[2] << "  l=0 Frequencies parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[3] << "  l=1 Frequencies parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[4] << "  l=2 Frequencies parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[5] << "  l=3 Frequencies parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[6] << "  Splitting parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[7] << "  Width parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[8] << "  Noise parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[9] << "  Inclination parameters" << std::endl; 
		std::cout << std::endl << "          - " << "Total Number of Parameters: " << all_in.plength.sum() << std::endl; 
		std::cout << " -----------------------------------------------------------" << std::endl;
	}
	
	//std::cout << "Exiting test " << std::endl;
	//exit(EXIT_SUCCESS);

return all_in;
}


short int set_noise_params(Input_Data *Noise_in, const MatrixXd noise_s2, const VectorXd noise_params){
/*
 *
 * A function that prepares the Noise_in Data structure using 
 * inputs of the .model file regrouped inside the noise_s2 matrixXd
 *
*/

	(*Noise_in).inputs_names[0]="Harvey-Noise_H"; (*Noise_in).inputs_names[3]="Harvey-Noise_H"; (*Noise_in).inputs_names[6]="Harvey-Noise_H";
	(*Noise_in).inputs_names[1]="Harvey-Noise_tc"; (*Noise_in).inputs_names[4]="Harvey-Noise_tc"; (*Noise_in).inputs_names[7]="Harvey-Noise_tc";
	(*Noise_in).inputs_names[2]="Harvey-Noise_p"; (*Noise_in).inputs_names[5]="Harvey-Noise_p"; (*Noise_in).inputs_names[8]="Harvey-Noise_p";
	(*Noise_in).inputs_names[9]="White_Noise_N0";

	(*Noise_in).priors_names[0]="Fix"; (*Noise_in).priors_names[3]="Fix"; (*Noise_in).priors_names[6]="Gaussian";
	(*Noise_in).priors_names[1]="Fix"; (*Noise_in).priors_names[4]="Fix"; (*Noise_in).priors_names[7]="Gaussian";
	(*Noise_in).priors_names[2]="Fix"; (*Noise_in).priors_names[5]="Fix"; (*Noise_in).priors_names[8]="Gaussian";
	(*Noise_in).priors_names[9]="Gaussian";

	(*Noise_in).relax[0]=0; (*Noise_in).relax[3]=0; (*Noise_in).relax[6]=1;
	(*Noise_in).relax[1]=0; (*Noise_in).relax[4]=0; (*Noise_in).relax[7]=1;
	(*Noise_in).relax[2]=0; (*Noise_in).relax[5]=0; (*Noise_in).relax[8]=1;
	(*Noise_in).relax[9]=1;
	(*Noise_in).inputs=noise_params;

	// Handle cases with negative H or tc ==> Harvey is Fix to 0 (no Harvey) <==> case of simulations with white noise
	if(((*Noise_in).inputs[0] <= 0) || ((*Noise_in).inputs[1] <= 0) || ((*Noise_in).inputs[2] <= 0)){
		(*Noise_in).priors_names[0]="Fix"; (*Noise_in).priors_names[1]="Fix"; (*Noise_in).priors_names[2]="Fix";
		(*Noise_in).relax[0]=0;            (*Noise_in).relax[1]=0;            (*Noise_in).relax[2]=0;
		(*Noise_in).inputs[0]=0;		    (*Noise_in).inputs[1]=0;		    (*Noise_in).inputs[2]=1;
	}
	if(((*Noise_in).inputs[3] <= 0) || ((*Noise_in).inputs[4] <= 0) || ((*Noise_in).inputs[5] <= 0)){
		(*Noise_in).priors_names[3]="Fix"; (*Noise_in).priors_names[4]="Fix"; (*Noise_in).priors_names[5]="Fix";
		(*Noise_in).relax[3]=0;            (*Noise_in).relax[4]=0;            (*Noise_in).relax[5]=0;
		(*Noise_in).inputs[3]=0;		    (*Noise_in).inputs[4]=0;		    (*Noise_in).inputs[5]=1;
	}
	if(((*Noise_in).inputs[6] <= 0) || ((*Noise_in).inputs[7] <= 0) || ((*Noise_in).inputs[8] <= 0)){
		(*Noise_in).priors_names[6]="Fix"; (*Noise_in).priors_names[7]="Fix"; (*Noise_in).priors_names[8]="Fix";
		(*Noise_in).relax[6]=0;            (*Noise_in).relax[7]=0;            (*Noise_in).relax[8]=0;
		(*Noise_in).inputs[6]=0;		    (*Noise_in).inputs[7]=0;		    (*Noise_in).inputs[8]=1;
	}
	// --- Center of the Gaussian -----
	(*Noise_in).priors(0,6)=noise_s2(6,0); 
	(*Noise_in).priors(0,7)=noise_s2(7,0);
	(*Noise_in).priors(0,8)=noise_s2(8,0);
	(*Noise_in).priors(0,9)=noise_s2(9,0);
	// --- sigma of the Gaussian
	(*Noise_in).priors(1,6)=(noise_s2(6,1) + noise_s2(6,2))*3./2;
	(*Noise_in).priors(1,7)=(noise_s2(7,1) + noise_s2(7,2))*3./2;
	if(noise_s2(8,1) !=0){ // If p has given errors then set uncertainty to 3*error
		(*Noise_in).priors(1,8)=(noise_s2(8,1) + noise_s2(8,2))*3./2;
	} else{ // If p is given with null-errors then set uncertainty to 0.1*p
		(*Noise_in).priors(1,8)=(*Noise_in).priors(0,8)*0.1;
	}
	//(*Noise_in).priors(1,9)=(noise_s2(9,1) + noise_s2(9,2))*10./2;
	(*Noise_in).priors(1,9)=(noise_s2(9,1) + noise_s2(9,2));
	if((((*Noise_in).priors(1,6)/(*Noise_in).priors(0,6)) <= 0.05) && (*Noise_in).priors_names[6] != "Fix"){ // If the given relative uncertainty on H3 is smaller than 5%
		std::cout << "Warning: The relative uncertainty on the Height of the high-frequency Harvey profile" << std::endl;
		std::cout << "         is smaller than 5% in the .MCMC file. This is too small" << std::endl;
		std::cout << "         ==> Relative uncertainty forced to be of 5%" << std::endl;
		(*Noise_in).priors(1,6)=(*Noise_in).priors(0,6)*0.05;
		std::cout << "	       Resuming..." << std::endl;
	}
	if((((*Noise_in).priors(1,7)/(*Noise_in).priors(0,7)) <= 0.005) && (*Noise_in).priors_names[7] != "Fix"){ // If the given relative uncertainty on tc3 is smaller than 5%
		std::cout << "Warning: The relative uncertainty on the timescale of the high-frequency Harvey profile" << std::endl;
		std::cout << "         is smaller than 0.5% in the .MCMC file. This is too small" << std::endl;
		std::cout << "         ==> Relative uncertainty forced to be of 0.5%" << std::endl;
		(*Noise_in).priors(1,7)=(*Noise_in).priors(0,7)*0.005;
		std::cout << "	       Resuming..." << std::endl;
	}
	if((((*Noise_in).priors(1,8)/(*Noise_in).priors(0,8)) <= 0.05) && (*Noise_in).priors_names[8] != "Fix"){ // If the given relative uncertainty on p is smaller than 5%
		std::cout << "Warning: The relative uncertainty on the power of the high-frequency Harvey profile" << std::endl;
		std::cout << "         is smaller than 5% in the .MCMC file. This is too small" << std::endl;
		std::cout << "         ==> Relative uncertainty forced to be of 5%" << std::endl;
		(*Noise_in).priors(1,8)=(*Noise_in).priors(0,8)*0.05;
		std::cout << "	       Resuming..." << std::endl;
	}
	if((((*Noise_in).priors(1,9)/(*Noise_in).priors(0,9)) <= 0.0005) && (*Noise_in).priors_names[9] != "Fix"){ // If the given relative uncertainty on N0 is smaller than 5%
		std::cout << "Warning: The relative uncertainty on the Height of the high-frequency Harvey profile" << std::endl;
		std::cout << "         is smaller than 0.05% in the .MCMC file. This is too small" << std::endl;
		std::cout << "         ==> Relative uncertainty forced to be of 0.05%" << std::endl;
		(*Noise_in).priors(1,9)=(*Noise_in).priors(0,9)*0.0005;
		std::cout << "	       Resuming..." << std::endl;
	}
	return 0;
}

short int fatalerror_msg_io_MS_Global(const std::string varname, const std::string param_type, const std::string syntax_vals, const std::string example_vals){
/*
* Function that handle error messages and warnings for io_MS_Global
*/
	if(param_type == "Fix_Auto"){
		std::cout << "         Warning: Fix_Auto is not implemented for that parameter" << std::endl;
		std::cout << "         		You must choose the prior yourself" << std::endl;
	} else{
		std::cout << "         Warning: " << varname << " should always be defined as '" << param_type << "'" << std::endl;
	}
	if(syntax_vals != ""){
		std::cout << "                  The syntax should be as follow: " << varname << "    " << param_type << "    " << syntax_vals << std::endl; 
	}
	if(example_vals !=""){
		std::cout << "                  Example: " << varname << "    " << param_type << "    " << example_vals << std::endl; 
	}
	std::cout << "         The program will exit now" << std::endl;
	exit(EXIT_FAILURE);
	return -1;
}

