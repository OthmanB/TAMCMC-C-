/*
 * priors_calc.cpp
 *
 *  Contains the function that calculate the penalization
 *  according to presets. So far, the only priors preset is for
 *  a main sequence star: priors_MS_Global()
 * 
 *  Created on: 02 Mar 2016
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "math.h"
#include "priors_calc.h"
#include "stats_dictionary.h"
#include "derivatives_handler.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

long double priors_MS_Global(const VectorXd params, const VectorXi params_length, const MatrixXd priors_params, const VectorXi priors_names_switch, const VectorXd extra_priors){

	long double f=0;

	const int smooth_switch=extra_priors[0];
	const double scoef=extra_priors[1];
	const double a3ova1_limit=extra_priors[2];
	const int impose_normHnlm=extra_priors[3];
	//const int numax=extra_priors[3]; 'Prior on numax, only applied if non-zero Not used.
	
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 3 for a global MS model (a1,a2,a3)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3; // Total number of modes

	double Dnu, d02, a1, alfa, b, fmax, Q11, max_b, el, em;
	Deriv_out frstder, scdder;

	// Apply the priors as defined in the configuration defined by the user and read by 'io_MS_global.cpp'
	f=f + apply_generic_priors(params, priors_params, priors_names_switch);

	// ----- Add a positivity condition on visibilities ------
	for(int i=Nmax; i<=Nmax+lmax; i++){
		if(params[i] < 0){
			f=-INFINITY; 
		}
	}
	// ----- Add a positivity condition on inclination -------
	// The prior could return values -90<i<90. We want it to give only 0<i<90
	//f=f+logP_uniform(0., 90., params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise]);
	switch(impose_normHnlm){ 
		case 1: // Case specific to model_MS_Global_a1etaa3_HarveyLike_Classic_v2
			//l=1: m=0, m=+/-1
			f=f+logP_uniform(0, 1.+1e-10, params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise] + 2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+1]); // The sum must be positive
			//l=2: m=0, m=+/-1, m=+/-2
 			// The sum must be positive
			f=f+logP_uniform(0, 1+1e-10, params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+2]+ 
								   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+3]+ 
								   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+4]
							); // The sum must be positive
			//l=3: m=0, m=+/-1, m=+/-2, m=+/-3
 			// The sum must be positive
			f=f+logP_uniform(0, 1+1e-10, params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+5]+ 
								   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+6]+ 
								   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+7]+
								   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+8]
								   ); // The sum must be positive

		break;
		case 2:
			std::cout << "priors_MS_Global: impose_normHnlm=2 YET TO BE IMPLEMENTED!" << std::endl;
			exit(EXIT_SUCCESS);
		break;
	}

	// Determine the large separation
	frstder=Frstder_adaptive_reggrid(params.segment(Nmax+lmax, Nfl0)); // First derivative of fl0 gives Dnu
	Dnu=frstder.deriv.sum();
	
	// Apply a prior on the d02
	//std::cout << "--- d02 --" << std::endl;
	if(Nfl0 == Nfl2){
		for(int i=0; i<Nfl0; i++){
			d02=params[Nmax+lmax+i] - params[Nmax+lmax+Nfl0+Nfl1+i];
			f=f+logP_gaussian_uniform( 0, Dnu/3., 0.015*Dnu, d02); // This is mainly for F stars
		}
	}
	// Set the smootheness condition handled by derivatives_handler.cpp
	switch(smooth_switch){
			case 1: // Case with smoothness
				//std::cout << " ------- Frequency derivatives ------" << std::endl;	
				if(Nfl0 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax, Nfl0)); // The l=0 frequencies
				}
				//std::cout << "--- Fl0 --" << std::endl;
				for(int i=0; i<Nfl0; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
				}
			
				if(Nfl1 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0, Nfl1)); // The l=1 frequencies
				}
				//std::cout << "--- Fl1 --" << std::endl;
				for(int i=0; i<Nfl1; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
				}
	
				if(Nfl2 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0+Nfl1, Nfl2)); // The l=2 frequencies
				}
				//std::cout << "--- Fl2 --" << std::endl;
				for(int i=0; i<Nfl2; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
	
				}
	
				if(Nfl3 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0+Nfl1+Nfl2, Nfl3)); // The l=3 frequencies
				}
				//std::cout << "--- Fl3 --" << std::endl;
				for(int i=0; i<Nfl3; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
				}
			  	break;
	}
	// Prior on a3/a1 ratio. a3 << a1 is enforced here by putting a3ova1_limit
	if(std::abs(params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+2]/params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3]) >= a3ova1_limit){
		f=-INFINITY;
	}

/*	std::cout << "a3 =" << params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+2] << std::endl;
	std::cout << "a1 ="	<< params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3] << std::endl;
	std::cout << "ratio =" << params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+2]/params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3] << std::endl;
	std::cout << "a3ova1_limit = " << a3ova1_limit << std::endl; 
	std::cout << "after a3ova1_limit " << f << std::endl;
	std::cout << "extra_priors = " << extra_priors << std::endl;
*/
	// Implement securities to avoid unphysical quantities that might lead to NaNs
	if(params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+4] < 0){ // Impose that the power coeficient of magnetic effect is positive
		f=-INFINITY;
	}
	if(priors_names_switch[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+3] != 0){
		if( (params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+3] < 0) || // Harvey profile height
		    (params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+4] < 0) || // Harvey profile tc
		    (params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+5] < 0) ){ // Harvey profile p
				f=-INFINITY;
		}
	}
	if(priors_names_switch[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+6] != 0){
		if( (params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+6] < 0) || // Harvey profile height
		    (params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+7] < 0) || // Harvey profile tc
		    (params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+8] < 0) ){// Harvey profile p 
				f=-INFINITY;
		}
	}
	if((priors_names_switch[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+9] != 0) && (params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+9] < 0)){
		f=-INFINITY;
	}
	//exit(EXIT_SUCCESS);
	
return f;
} 


long double priors_local(const VectorXd params, const VectorXi params_length, const MatrixXd priors_params, const VectorXi priors_names_switch, const VectorXd extra_priors){

	long double f=0;

	//const int smooth_switch=extra_priors[0];
	const double a3ova1_limit=extra_priors[2];
	const int impose_normHnlm=extra_priors[3];
	//const int numax=extra_priors[3]; 'Prior on numax, only applied if non-zero Not used.
	
	const int Nmax=params_length[0]; // Number of Heights
	const int Nvis=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 3 for a global MS model (a1,a2,a3)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
	const int lmax=3; // lmax is hardcoded here

	int pos0;
	VectorXi Nf_el(4);  
	Nf_el[0]=Nfl0; Nf_el[1]=Nfl1; Nf_el[2]=Nfl2; Nf_el[3]=Nfl3; 

	//double Dnu, d02, scoef, a1, alfa, b, fmax, Q11, max_b, el, em;
	
 	// Apply the priors as defined in the configuration defined by the user and read by 'io_MS_global.cpp'
	f=f + apply_generic_priors(params, priors_params, priors_names_switch);
	// ----- Add a positivity condition on inclination -------
	// The prior could return values -90<i<90. We want it to give only 0<i<90
	//f=f+logP_uniform(0., 90., params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise]);


// ALL THIS IS NOT NECESSARY BECAUSE WE DO NOT DEAL WITH VISIBILITIES... NO NORMALISATION TO 1 REQUIRED
// THEREFORE NO NEED TO EXCLUDE VALUES BEYOND 1.
/*	switch(impose_normHnlm){ 
		case 1: // Case specific to model_MS_Global_a1etaa3_HarveyLike_Classic_v2
			//l=1: m=0, m=+/-1
			std::cout << "priors_MS_local: impose_normHnlm=1 Not required!" << std::endl;
			std::cout << "                 Likely bug in io_local.cpp ... Exiting now" << std::endl; 
			exit(EXIT_SUCCESS);	
			break;
		case 2:
			std::cout << "In case 2" << std::endl;
			std::cout << "lmax=" << lmax << std::endl;
			for(int el=1; el<=lmax; el++){
					for(int en=0; en<Nf_el[el]; en++){
						for(int em=0; em<=el; em++){
							std::cout << "el=" << el << "   ==> " <<  Nf_el.segment(0,el+1)  << std::endl;
							pos0=Nf_el.segment(0,el+1).sum() + (el+1)*en;
							std::cout << pos0 << std::endl;
							switch(el){
								case 1: //l=1: m=0, m=+/-1
									std::cout << "In l=1: m=0, m=+/-1" << std::endl;
									f=f+logP_uniform(0, 1, params[pos0] + 2*params[pos0+1]); // The sum must be positive
									break;
								case 2: //l=2: m=0, m=+/-1, m=+/-2
									std::cout << "In l=2: m=0, m=+/-1, m=+/-2" << std::endl;
									f=f+logP_uniform(0, 1, params[pos0] + 2*params[pos0+1] + 2*params[pos0+2]); // The sum must be positive
									break;
								case 3:
									std::cout << "In l=3: m=0, m=+/-1, m=+/-2, m=+/-3" << std::endl;
									f=f+logP_uniform(0, 1, params[pos0]+ 2*params[pos0+1]+ 2*params[pos0+2] + 2*params[pos0+3]); // The sum must be positive
									break;
							}
						}
					}
				}

			break;
	}
*/
	// Prior on a3/a1 ratio. a3 << a1 is enforced here by putting a3ova1_limit
	if (params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3] !=0){ // If there is a a1 to be fitted
		if(std::abs(params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+2]/params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3]) >= a3ova1_limit){
			f=-INFINITY;
		}
	} else{
		if ((params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+3] !=0) &&  (params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+4] !=0)){ // If the sqrt(a1).cosi and sqrt(a1).sin(i) are not 0
			//std::cout << "Case sqrt(a1).cos and sin" << std::endl;
			if(std::abs(params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+2]/(pow(params[Nmax + Nvis + Nfl0+Nfl1+Nfl2+Nfl3+3],2)+ pow(params[Nmax + Nvis + Nfl0+Nfl1+Nfl2+Nfl3+4],2))) >= a3ova1_limit){
				f=-INFINITY;
			}
		}
	}
	//std::cout << f << std::endl;
	
	// Prior in inclination positivity (if relevant, ie if it is directly fitted)
	//std::cout << " priors_names_switch[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] = " << priors_names_switch[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] << std::endl;
	//std::cout << " params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] =" << params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] << std::endl;
	if((priors_names_switch[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] != 0) && (params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] < 0)){
		f=-INFINITY;
	}
	//std::cout << f << std::endl;

	// BELOW: COMMENTED ON 13 AUG 2020 AS THESE LINES MAKE NO SENSE TO ME AND INDUCE A CRASH IN SOME OCCASION. I SUSPECT THAT THIS WAS POSITIVITY CONDITION ON INC
	// WAS REPLACE BY THE LINE ABOVE BUT KEPT HERE IN CASE I FINALY UNDERSTAND THOSE LINES
	//if((priors_names_switch[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+9] != 0) && (params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+9] < 0)){
	//	f=-INFINITY;
	//}

	//exit(EXIT_SUCCESS);
	//std::cout << "Priors done" << std::endl;
return f;
} 







//long double priors_Test_Gaussian(const VectorXd params, const VectorXi param_length, const MatrixXd priors_params, const std::vector<std::string> priors_names){
long double priors_Test_Gaussian(const VectorXd params, const VectorXi params_length, const MatrixXd priors_params, const VectorXi priors_names_switch){

	long double f=0;

	f=f + apply_generic_priors(params, priors_params, priors_names_switch);

	//// More rigid strategy... How much time is spent on the cases ??
	//for(i=0; i<params.size(); i++){
	//	f=f + logP_uniform(priors_params(0, i), priors_params(1, i), params[i]);
	//}
	//exit(EXIT_SUCCESS);
return f;
} 

//long double priors_Test_Gaussian(const VectorXd params, const VectorXi param_length, const MatrixXd priors_params, const std::vector<std::string> priors_names){
long double priors_Harvey_Gaussian(const VectorXd params, const VectorXi params_length, const MatrixXd priors_params, const VectorXi priors_names_switch){

	long double f=0;

	f=f + apply_generic_priors(params, priors_params, priors_names_switch);

return f;
} 



long double apply_generic_priors(const VectorXd params, const MatrixXd priors_params, const VectorXi priors_names_switch){
/*
 * This function apply the generic priors, that are defined in the configuration file and that have been translated into
 * The constant priors_params (the parameters that define the priors) and into priors_names_switch (translation of 
 * priors_names into integers in order to be handled by the switch / case statement.
 *
*/
	int i;
	long double pena=0;
	for(i=0; i<params.size(); i++){
		switch(priors_names_switch[i]){
			case 0: // No Prior Applied --> The parameter is Fixed or for some reason the prior is set to NONE
			  //std::cout << "CASE 0" << std::endl;
			  //std::cout << "[" << i << "] FIX OR NONE, penalizationa=" << pena << std::endl;
			  break;
			case 1: // Uniform Prior
			  pena=pena + logP_uniform(priors_params(0, i), priors_params(1, i), params[i]);
			  //std::cout << "CASE 1" << std::endl;
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i) << std::endl;
			  //std::cout << "    pena=" << pena  << std::endl;
			  break;
			case 2: // Gaussian Prior
			  pena=pena + logP_gaussian(priors_params(0, i), priors_params(1, i),params[i]);
			  //std::cout << "CASE 2" << std::endl;
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i) << std::endl;
			  //std::cout << "    pena=" << pena  << std::endl;
			  break;
			case 3: // Multivariate Gaussian Prior
			  std::cout << "The Multivariate Gaussian Prior IS NOT HANDLED" << std::endl;
		      std::cout << "The program will exit now" << std::endl;
			  exit(EXIT_FAILURE);
			  //pena=pena + logP_multivariate_gaussian( priors_params(0, i), Matrix, params[i]);
			  break;
			case 4: // Jeffreys Prior
			  //std::cout << "CASE 4" << std::endl;
			  pena=pena + logP_jeffrey(priors_params(0, i), priors_params(1, i), params[i]);
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i) << std::endl;
			  //std::cout << "    pena=" << pena  << std::endl;	  
			  break;
			case 5: // Uniform-Gaussian Prior
			  //std::cout << "CASE 5" << std::endl;
			  pena=pena + logP_uniform_gaussian( priors_params(0, i), priors_params(1, i), priors_params(2, i), params[i]); // bmin, bmax, sigma, x
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i);
			  //std::cout << "  " << "priors_params(2, " << i << ")=" << priors_params(2, i) << std::endl;
			  //std::cout << "    pena=" << pena << std::endl;			  
			  break;
			case 6: // Gaussian-Uniform Prior
			  //std::cout << "CASE 6" << std::endl;
			  pena=pena + logP_gaussian_uniform( priors_params(0, i), priors_params(1, i), priors_params(2, i), params[i]); // bmin, bmax, sigma, x
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i);
			  //std::cout << "  " << "priors_params(2, " << i << ")=" << priors_params(2, i) << std::endl;
			  //std::cout << "    pena=" << pena << std::endl;
			  break;
			case 7: // Gaussian-Uniform-Gaussian Prior
			  //std::cout << "CASE 7" << std::endl;
			  pena=pena + logP_gaussian_uniform_gaussian( priors_params(0, i), priors_params(1, i), priors_params(2, i),  priors_params(3, i), params[i]); // bmin, bmax, sigma1, sigma2, x
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i);
			  //std::cout << "  " << "priors_params(2, " << i << ")=" << priors_params(2, i) << std::endl;
			  //std::cout << "  " << "priors_params(3, " << i << ")=" << priors_params(3, i) << std::endl;
			  //std::cout << "    pena=" << pena << std::endl;
			  break;
			case 8: // Uniform Prior on abs(params[i])
			  //std::cout << "CASE 8" << std::endl;
			  pena=pena + logP_uniform_abs(priors_params(0, i), priors_params(1, i), params[i]);
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i);
			  //std::cout << "    pena=" << pena << std::endl;
			  break;
			case 9: // Uniform Prior on cos(params[i]*pi/180)
			  //std::cout << "CASE 9" << std::endl;
			  pena=pena + logP_uniform_cos(priors_params(0, i), priors_params(1, i), params[i]);
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i);
			  //std::cout << "    pena=" << pena << std::endl;
			  break;
			case 10: // Jeffreys Prior on abs(param[i])
			  //std::cout << "CASE 4" << std::endl;
			  pena=pena + logP_jeffrey_abs(priors_params(0, i), priors_params(1, i), params[i]);
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i) << std::endl;
			  //std::cout << "    pena=" << pena  << std::endl;	  
			  break;
			case 11: // Case Auto: No Prior Applied here --> The user must set his own prior
			  //std::cout << "CASE 10" << std::endl;
			  //std::cout << "[" << i << "] Auto, pena=" << pena << std::endl;
			  break;
			default:
			  std::cout << " Problem in priors_calc.cpp! " << std::endl;
			  std::cout << " priors_names_switch[" << i << "]=" << priors_names_switch[i] << std::endl;
		          std::cout << " This value is not associated to any known case statement " << std::endl;
			  std::cout << " The program will exit now" << std::endl;
			  exit(EXIT_FAILURE);
		}
	}

return pena;
}


