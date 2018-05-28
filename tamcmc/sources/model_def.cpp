/*
 * model_def.cpp
 *
 * Contains the methods associated to the class 'Model'
 * The header file of that class is in encapsulators.h
 * Note that all models should be put in models.cpp
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include "model_def.h"
#include "models.h"   // This contains the models that can be used
#include "likelihoods.h"  // This contains the likelihoods that can be used
#include "priors_calc.h"  // This contains the priors that can be used
//#include "config.h"
//#include "data.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;


Model_def::Model_def(Config *config, VectorXd Tcoefs, bool verbose){
	
	bool error;
	double warning_thld;
	
	Nmodels=config->MALA.Nchains;
	model_fct_name=config->modeling.model_fct_name;
	likelihood_fct_name=config->modeling.likelihood_fct_name;
	prior_fct_name=config->modeling.prior_fct_name;
	model_fct_name_switch=config->modeling.model_fct_name_switch;
	likelihood_fct_name_switch=config->modeling.likelihood_fct_name_switch;
	prior_fct_name_switch=config->modeling.prior_fct_name_switch;
	
	params_names=config->modeling.inputs.inputs_names;
	priors_params_names=config->modeling.inputs.priors_names;
	priors_params_names_switch=config->modeling.inputs.priors_names_switch;
	relax=config->modeling.inputs.relax;
	plength=config->modeling.inputs.plength;
	extra_priors=config->modeling.inputs.extra_priors;

	likelihood_params=config->modeling.likelihood_params;

	Nparams=plength.sum();
	params.resize(Nmodels, Nparams);

	for (int m=0; m<Nmodels; m++){
		params.row(m)=config->modeling.inputs.inputs;
	}
	priors_params=config->modeling.inputs.priors;
	
	Pmove.setZero(Nmodels); // initialize the value of the move probabilty for all models
	moved.assign(Nmodels, 0); // initialize to False the move action indicator for all models

	swaped=0;
	Pswap=0;

	// Additional quantities from the random generator (for debug only)
	comparator_MH.setZero(Nmodels); // initialize the value of the move probabilty for all models
	comparator_PT=0; // initialize the value of the move probabilty for all models

	// --- initialise vectors of variables/constant ------
	Nvars=0; // initialize the number of variables to 0
	Ncons=0; // initialize the number of constant to 0
	vars.resize(Nmodels, Nparams); // set to maximum possible size. This must be 2D
	cons.resize(Nparams);
	index_to_relax.resize(Nparams);
	vars_names.resize(Nparams); // set to maximum possible size this is 1D
	cons_names.resize(Nparams);


	if( relax.sum() <= 1){
		std::cout << "Less than one parameter free! No minimization possible!" << std::endl;
		std::cout << "You need to set relax to 1 for at least two parameters" << std::endl;
	}
	for (int i=0; i<Nparams;i++){
		if(relax[i] == 1){
			for(int mod=0; mod<Nmodels; mod++){
				vars(mod,Nvars)=params(mod,i);
			}
			vars_names[Nvars]=params_names[i];
			index_to_relax[Nvars]=i; // used to quickly refresh params when calculating all the models
			Nvars=Nvars+1;
		} else{
			cons[Ncons]=params(0,i);
			cons_names[Ncons]=params_names[i];
			Ncons=Ncons+1;
		}
	}

	vars.conservativeResize(Nmodels, Nvars); // Resize by keeping the elements inside
	vars_names.resize(Nvars); // This is a vector... we can do a normal resize
	cons.conservativeResize(Ncons); // This is also a vector, not a VectorXd or a MatrixXd
	cons_names.resize(Ncons); // again a vector...
	index_to_relax.resize(Nvars); // vector with position of params that vary

	// ------ Restoring last values on the restore file if requested -------
	error=0;
	if(config->outputs.do_restore_variables == 1 || config->outputs.do_restore_proposal == 1){
		std::cout << "   [4] Restore of old variable requested..." << std::endl;
		std::cout << "        ...Performing consistency checks for restored variables..." << std::endl;
		if(config->restored_vals.Nvars != Nvars){
			std::cout << "    - Inconsistency in the number of variables: The config.cfg describes a model with Nvars=" << Nvars << " while the restored files suggest Nvars= " << config->restored_vals.Nvars << std::endl; 
		error=1;
		}
		if(config->restored_vals.Nchains != Nmodels){
			std::cout << "    - Inconsistency in the number of chains: The config.cfg describes a model with Nchains=" << Nmodels << " while the restored files suggest Nchains= " << config->restored_vals.Nchains << std::endl; 
		error=1;
		}
		if(error == 1){
			std::cout << "The program will stop now" << std::endl;
			exit(EXIT_FAILURE);
		} else {
			std::cout << "Tests passed..." << std::endl;
			if (config->outputs.do_restore_variables == 1){
				std::cout << "       Setting up the MCMC using the restored variables..." << std::endl;
				vars=config->restored_vals.vars;
				for(int m=0; m<Nmodels; m++){				
					update_params_with_vars(m);
				}
			}
		}
	}

	// Initialize the model, logLikelihood, logPrior, logPosterior variables
	model.resize(Nmodels, (*config).data.data.Nx);
	logLikelihood.resize(Nmodels);
	logPrior.resize(Nmodels);
	logPosterior.resize(Nmodels);

	Data data_in=config->data.data;
	bool empty_container=0; // Used to know whether we fill MatrixXd/VectorXd for all Nmodels
	if(empty_container == 0){
		for(int m=0; m<Nmodels; m++){
			generate_model(&data_in, m, Tcoefs); // No need of the returned value
		} 
	}

	if(verbose == 1){	
		warning_thld=5000.;
		std::cout << "Checking that there is not problem with the likelihood or the priors..." << std::endl;
		std::cout << "      - logLikelihood[" << 0 << "]=" << logLikelihood[0]  << std::endl;
		std::cout << "      - logPrior[" << 0 << "]=" << logPrior[0]  << std::endl;
		std::cout << "      - logPosterior[" << 0 << "]=" << logPosterior[0]  << std::endl;
		if((logLikelihood[0] == INFINITY) || (logLikelihood[0] == -INFINITY)){
			std::cout << " --------------------------------- WARNING ----------------------------------" << std::endl;
			std::cout << "                  The initial logLikelihood[0] is INFINITY "<< std::endl;
			std::cout << "    Your initial vector of parameter must have a FINITE initial likelihood   " << std::endl;
			std::cout << "                         The program will exit now " << std::endl;
			std::cout << " ---------------------------------------------------------------------------" << std::endl;
		}
		if(std::abs(logPrior[0]) > warning_thld){
			std::cout << " --------------------------------- WARNING ----------------------------------" << std::endl;
			std::cout << "                    The initial logPrior[0] >" << warning_thld << std::endl;
			std::cout << "  The penalisation arising from the priors is exceeding the warning threshold " << std::endl;
			std::cout << "   Check that your initial vector of parameter is compatible with the priors" << std::endl;
			std::cout << "                         The program will exit now " << std::endl;
			std::cout << " ---------------------------------------------------------------------------" << std::endl;
		}
	
		std::cout << "   - Variables " << std::endl;
		for(int i=0; i<Nvars; i++){
			std::cout << vars_names[i]  << " | "  << vars(0,i) << std::endl;
		}
		std::cout << "   - Constant " << std::endl;
		for(int i=0; i<Ncons; i++){
			std::cout << cons_names[i] << " | " << cons[i] << std::endl;
		}
		std::cout << " ------------- " << std::endl;
	}

}

Model_def::Model_def(){ // The Empty constructor (overload)

	// Empty because no initialisation of internal variable is made when calling this constructor
}

Model_def::~Model_def(){ // The destructor

	// Empty because no pointer are declared into the class ==> destruction handled by the compiler
}

/* 
* 
* Allows explicit call of call_model function without the need of setting internal class variables
* To use when model_def class has to be used only for a few model computation.
* Do not use this on large loop... it will be slow
*/
VectorXd Model_def::call_model_explicit(Data *data_struc, const VectorXi plength0, const VectorXd params0, const int model_case){


	params.resize(1, params0.size());
	params.row(0)=params0;
	model_fct_name_switch=model_case;
	plength=plength0;
	return call_model(data_struc, 0);

}

VectorXd Model_def::call_model(Data *data_struc, int m){

	VectorXd fail(data_struc->Nx);
	fail.setZero();

	switch(model_fct_name_switch){
		case 0: // model_Test_Gaussian
		  return model_Test_Gaussian(params.row(m), plength, (*data_struc).x);
		  break;
		case 1: // model_Harvey_Gaussian
		  return model_Harvey_Gaussian(params.row(m), plength, (*data_struc).x);
		  break;
		case 2: // model_MS_Global with a1, eta (asphericity), a3, asymetry, Generalized Harvey function. 
		  return model_MS_Global_a1etaa3_HarveyLike(params.row(m), plength, (*data_struc).x);
		  break;
		case 3: // model_MS_Global with a1, eta (asphericity), a3, asymetry, Generalized Harvey function. Inclination and splitting are not fitted directly. sqrt(a1).cos(i) and sqrt(a1).sin(i) instead
		  return model_MS_Global_a1etaa3_HarveyLike_Classic(params.row(m), plength, (*data_struc).x);
		  break;
		case 4: // model_MS_Global with a1, eta (asphericity), a3, asymetry, Original Harvey function 
		  return model_MS_Global_a1etaa3_Harvey1985(params.row(m), plength, (*data_struc).x);
		  break;
		case 5: // model_MS_Global with a1, magb and magalfa (asphericity), a3, asymetry, Original Harvey function 
		  std::cout << "Obselete model that is not anymore supported: model_MS_Global_a1acta3_HarveyLike" << std::endl;
		  std::cout << "The program will exit now" << std::endl;
		  exit(EXIT_FAILURE);
		  return fail;
		  //return model_MS_Global_a1acta3_Harvey1985(params.row(m), plength, (*data_struc).x);
		  break;
        	case 6:// model_MS_Global with a1(l=1), a1(2), a1(3)=(a1(1)+a1(2))/2, eta (asphericity), a3, asymetry, Generalized Harvey function
            	  return model_MS_Global_a1l_etaa3_HarveyLike(params.row(m), plength, (*data_struc).x);
		  break;
        	case 7:// model_MS_Global with a1(n, l=1)=a1(n, l=2), a1(3)=(a1(1)+a1(2))/2, eta (asphericity), a3, asymetry, Generalized Harvey function
            	  return model_MS_Global_a1n_etaa3_HarveyLike(params.row(m), plength, (*data_struc).x);
		  break;
        	case 8:// model_MS_Global with a1(n, l), a1(3)=(a1(n,1)+a1(n,2))/2, eta (asphericity), a3, asymetry, Generalized Harvey function
            	  return model_MS_Global_a1nl_etaa3_HarveyLike(params.row(m), plength, (*data_struc).x);
          break;
		default:
		  std::cout << " Problem in model_def.cpp! " << std::endl;
		  std::cout << " model_fct_names_switch = " << model_fct_name_switch << std::endl;
		  std::cout << " This value is not associated to any known case statement " << std::endl;
		  std::cout << " Keywords with valid statements so far:" << std::endl;
		  std::cout << "    - 'Test_Gaussian' (For Debug only)" << std::endl;
		  std::cout << "    - 'model_Harvey_Gaussian'" << std::endl;
		  std::cout << "    - 'model_MS_Global_a1etaa3_HarveyLike'" << std::endl;
		  std::cout << "    - 'model_MS_Global_a1etaa3_HarveyLike_Classic'" << std::endl;
		  std::cout << "    - 'model_MS_Global_a1etaa3_Harvey1985'" << std::endl;

          std::cout << "    - 'model_MS_Global_a1l_etaa3_HarveyLike'" << std::endl;
		  std::cout << "    - 'model_MS_Global_a1n_etaa3_HarveyLike'" << std::endl;
          std::cout << "    - 'model_MS_Global_a1nl_etaa3_HarveyLike'" << std::endl;
		  std::cout << " The program will exit now" << std::endl;
		  exit(EXIT_FAILURE);
	}


	return fail; // This might never be used but avoid warnings from the compiler
}

long double Model_def::call_likelihood(Data *data_struc, int m, VectorXd Tcoefs){
/* 
 * call the likelihood using its name. Before returning the log(Likelihood), it saves it into 'logLikelihood'
*/
	double p;
	long double logL;

	switch(likelihood_fct_name_switch){
		case 0: // chi(2,2p)
		  p=likelihood_params; // to construct the chi(2,2p) statistics, we need p
		  logL=likelihood_chi22p((*data_struc).y, model.row(m), p); // 'Model_def::model' must be calculated prior to the execution of this function
		  return logL/Tcoefs[m];
		  break;
		case 1: // chi_square
		  logL=likelihood_chi_square((*data_struc).y, model.row(m), (*data_struc).sigma_y); // Here no need of likelihood_params
		  return logL/Tcoefs[m];
		  break;
		default:
		  std::cout << " Problem in model_def.cpp! " << std::endl;
		  std::cout << " likelihood_fct_names_switch = " << likelihood_fct_name_switch << std::endl;
		  std::cout << " This value is not associated to any known case statement " << std::endl;
		  std::cout << " Keywords with valid statements so far:" << std::endl;
		  std::cout << "    - 'chi(2,2p)'" << std::endl;
		  std::cout << "    - 'chi_square'" << std::endl;
		  std::cout << " The program will exit now" << std::endl;
		  exit(EXIT_FAILURE);
	}

	return 0; // This might never be used but avoid warnings from the compiler
}

long double Model_def::call_prior(Data *data_struc, int m){
/* 
 * call the prior using its name. Before returning the log(prior), we save it into 'logPrior'
*/
	//bool passed=0;
	double p;

	switch(prior_fct_name_switch){
		case 0: // model_Test_Gaussian
		  return priors_Test_Gaussian(params.row(m), plength, priors_params, priors_params_names_switch);
		  break;
		case 1: // model_Harvey_Gaussian
		  return priors_Harvey_Gaussian(params.row(m), plength, priors_params, priors_params_names_switch);
		  break;
		case 2: // model_MS_Global
		  return priors_MS_Global(params.row(m), plength, priors_params, priors_params_names_switch, extra_priors);
		  break;
		default:
		  std::cout << " Problem in model_def.cpp! " << std::endl;
		  std::cout << " prior_fct_name_switch = " << prior_fct_name_switch << std::endl;
		  std::cout << " This value is not associated to any known case statement " << std::endl;
		  std::cout << " Keywords with valid statements so far:" << std::endl;
		  std::cout << "    - 'Test_Gaussian' (For Debug only)" << std::endl;
		  std::cout << "    - 'model_Harvey_Gaussian'" << std::endl;
		  std::cout << "    - 'MS_Global'" << std::endl;
		  std::cout << " The program will exit now" << std::endl;
		  exit(EXIT_FAILURE);
	}

	return 0; // This might never be used but avoid warnings from the compiler
}

long double Model_def::generate_model(Data *data_struc, long m, VectorXd Tcoefs){
/*
 * call successively call_model, call_likelihood and call_prior and then calculates the logPosterior. This is also returned.
*/
	model.row(m)=call_model(data_struc, m);
	logLikelihood[m]=call_likelihood(data_struc, m, Tcoefs); // logLikelihood saved at element m of the vector
	logPrior[m]=call_prior(data_struc, m); // logprior saved at element m of the vector
	logPosterior[m]= logLikelihood[m] + logPrior[m]; // logPosterior saved at element m of the vector
	
return logPosterior[m];
}


void Model_def::update_params_with_vars(long m){

	for( int i=0; i<index_to_relax.size(); i++){
		//std::cout << "i=" << i << std::endl;
		params(m,index_to_relax[i])=vars(m,i);
	}
	

}


