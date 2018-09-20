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

	//const long double pi = 3.141592653589793238462643383279502884L;
	//const long double G=6.667e-8;
	//const long double Teff_sun= 5777; 
	//const long double Dnu_sun=135.1;
	//const long double numax_sun=3150.;
	//const long double R_sun=6.96342e5; //in km
	//const long double M_sun=1.98855e30; //in kg
	//const long double rho_sun=M_sun*1e3/(4L*pi*std::pow(R_sun*1e5,3L)/3L); //in g.cm-3

	const int smooth_switch=extra_priors[0];
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
	
	double Dnu, d02, scoef, a1, alfa, b, fmax, Q11, max_b, el, em;
	Deriv_out frstder, scdder;

	//std::cout << "Initial f=" << f << std::endl;	
	
	// Apply the priors as defined in the configuration defined by the user and read by 'io_MS_global.cpp'
	f=f + apply_generic_priors(params, priors_params, priors_names_switch);
	//std::cout << "After generic priors ==> f=" << f << std::endl;	

	// ----- Add a positivity condition on visibilities ------
	for(int i=Nmax; i<=Nmax+lmax; i++){
		if(params[i] < 0){
			f=-INFINITY; 
		}
	}
	// ----- Add a positivity condition on inclination -------
	// The prior on cos(i) is returning values -90<i<90. We want it to give only 0<i<90
	f=f+logP_uniform(0., 90., params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise]);
	
	// Determine the large separation
	frstder=Frstder_adaptive_reggrid(params.segment(Nmax+lmax, Nfl0)); // First derivative of fl0 gives Dnu
	Dnu=frstder.deriv.sum();
	
	// Apply a prior on the d02
	//std::cout << "--- d02 --" << std::endl;
	if(Nfl0 == Nfl2){
		for(int i=0; i<Nfl0; i++){
			d02=params[Nmax+lmax+i] - params[Nmax+lmax+Nfl0+Nfl1+i];
			//if (abs(d02) <= Dnu/3.){
			f=f+logP_gaussian_uniform( 0, Dnu/3., 0.015*Dnu, d02); // This is mainly for F stars
			//}
			//std::cout << "  " << d02;
		}
	}
	//std::cout << std::endl;
	//std::cout << "After prior on d02 ==> f=" << f << std::endl;	
		
	
	// Set the smootheness condition handled by derivatives_handler.cpp
	//scoef=2.;
	//std::cout << extra_priors << std::endl;
	//exit(EXIT_SUCCESS);
	switch(smooth_switch){
			case 1: // Case with smoothness
				scoef=extra_priors[1];
				//std::cout << " ------- Frequency derivatives ------" << std::endl;	
				if(Nfl0 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax, Nfl0)); // The l=0 frequencies
				}
				//std::cout << "--- Fl0 --" << std::endl;
				for(int i=0; i<Nfl0; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
					//std::cout << "  " << scdder.deriv[i];
				}
				//std::cout << std::endl;
				//std::cout << "After scd der on fl0 ==> f=" << f << std::endl;	
			
				if(Nfl1 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0, Nfl1)); // The l=1 frequencies
				}
				//std::cout << "--- Fl1 --" << std::endl;
				for(int i=0; i<Nfl1; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
					//std::cout << "  " << scdder.deriv[i];	
				}
				//std::cout << std::endl;
				//std::cout << "After scd der on fl1 ==> f=" << f << std::endl;	
	
				if(Nfl2 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0+Nfl1, Nfl2)); // The l=2 frequencies
				}
				//std::cout << "--- Fl2 --" << std::endl;
				for(int i=0; i<Nfl2; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
					//std::cout << "  " << scdder.deriv[i];
	
				}
				//std::cout << std::endl;
				//std::cout << "After scd der on fl2 ==> f=" << f << std::endl;	
	
				if(Nfl3 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0+Nfl1+Nfl2, Nfl3)); // The l=3 frequencies
				}
				//std::cout << "--- Fl3 --" << std::endl;
				for(int i=0; i<Nfl3; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
					//std::cout << "  " << scdder.deriv[i];
				}
				//std::cout << std::endl;
				//std::cout << "After scd der on fl3 ==> f=" << f << std::endl;	
			  	break;
	}

		
	// If requested, we impose a user-defined prior on the intensity coeficient for the magnetic effect on rotational splitting
	if(priors_names_switch[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+3] == 10){
		a1=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3];
		alfa=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+4];
		b=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+3];
		fmax=priors_params.block(1, Nmax+lmax, 1, Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3).maxCoeff();
		el=1.; em=1.;
		Q11=(el*(el+1) - 3.*pow(em,2))/( (2*el - 1)*(2*el + 3) );
		max_b=1.1*std::abs( a1/(Q11 * pow(fmax*1e-3,alfa)) ); //Beware, it is an hyperparameter
		f=f+logP_jeffrey(0.05, max_b, std::abs(b));
		
		std::cout << "This part of the code has to be checked" << std::endl;
		std::cout << "a1="<<a1 << std::endl;
		std::cout << "alfa="<<alfa << std::endl;
		std::cout << "b="<<b << std::endl;
		std::cout << "Q11="<<Q11 << std::endl;
		std::cout << "fmax="<<fmax << std::endl;
		std::cout << "max_b="<<max_b << std::endl;
		std::cout << "Review these values and validate that everything is alright" << std::endl;
		std::cout << "The program will stop now"<< std::endl;
		exit(EXIT_SUCCESS);	
	}
	//std::cout << "After user defined prior ==> f=" << f << std::endl;
		
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
	//std::cout << "After securities ==> f=" << f << std::endl;	  
	
	/*std::cout << " ---- Noise slot ----" << std::endl;
	std::cout << "H_harvey=" << params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+3] << std::endl;
	std::cout << "White noise=" << params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+9] << std::endl;
	std::cout << "Inc=" << params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+10] << std::endl;
	
	exit(EXIT_SUCCESS);
	*/
return f;
} 

long double priors_evolved_Global(const VectorXd params, const VectorXi params_length, const MatrixXd priors_params, const VectorXi priors_names_switch, const VectorXd extra_priors){

	long double f=0;

	//const long double pi = 3.141592653589793238462643383279502884L;
	//const long double G=6.667e-8;
	//const long double Teff_sun= 5777; 
	//const long double Dnu_sun=135.1;
	//const long double numax_sun=3150.;
	//const long double R_sun=6.96342e5; //in km
	//const long double M_sun=1.98855e30; //in kg
	//const long double rho_sun=M_sun*1e3/(4L*pi*std::pow(R_sun*1e5,3L)/3L); //in g.cm-3

	const int smooth_switch=extra_priors[0];
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
	const int Ntrunc=params_length[10]; // Truncation parameter
	const int NHl1=params_length[11]; // Heights l=1
	const int NWl1=params_length[12]; // Widths l=1
	const int Nsplit1=params_length[13]; // Splitting l=1
	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	
	double Dnu, d02, scoef, a1, alfa, b, fmax, Q11, max_b, el, em;
	Deriv_out frstder, scdder;

	//std::cout << "Initial f=" << f << std::endl;	
	
	// Apply the priors as defined in the configuration defined by the user and read by 'io_MS_global.cpp'
	f=f + apply_generic_priors(params, priors_params, priors_names_switch);
	//std::cout << "After generic priors ==> f=" << f << std::endl;	

	// ----- Add a positivity condition on visibilities ------
	for(int i=Nmax; i<=Nmax+lmax; i++){
		if(params[i] < 0){
			f=-INFINITY; 
		}
	}
	// ----- Add a positivity condition on inclination -------
	// The prior on cos(i) is returning values -90<i<90. We want it to give only 0<i<90
	f=f+logP_uniform(0., 90., params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise]);
	
	// Determine the large separation
	frstder=Frstder_adaptive_reggrid(params.segment(Nmax+lmax, Nfl0)); // First derivative of fl0 gives Dnu
	Dnu=frstder.deriv.sum();
	
	// Apply a prior on the d02
	//std::cout << "--- d02 --" << std::endl;
	if(Nfl0 == Nfl2){
		for(int i=0; i<Nfl0; i++){
			d02=params[Nmax+lmax+i] - params[Nmax+lmax+Nfl0+Nfl1+i];
			//if (abs(d02) <= Dnu/3.){
			f=f+logP_gaussian_uniform( 0, Dnu/3., 0.015*Dnu, d02); // This is mainly for F stars
			//}
			//std::cout << "  " << d02;
		}
	}
	//std::cout << std::endl;
	//std::cout << "After prior on d02 ==> f=" << f << std::endl;	
		
	
	// Set the smootheness condition handled by derivatives_handler.cpp
	//scoef=2.;
	//std::cout << extra_priors << std::endl;
	//exit(EXIT_SUCCESS);
	switch(smooth_switch){
			case 1: // Case with smoothness
				scoef=extra_priors[1];
				//std::cout << " ------- Frequency derivatives ------" << std::endl;	
				if(Nfl0 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax, Nfl0)); // The l=0 frequencies
				}
				//std::cout << "--- Fl0 --" << std::endl;
				for(int i=0; i<Nfl0; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
					//std::cout << "  " << scdder.deriv[i];
				}
				//std::cout << std::endl;
				//std::cout << "After scd der on fl0 ==> f=" << f << std::endl;	
			
				if(Nfl1 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0, Nfl1)); // The l=1 frequencies
				}
				//std::cout << "--- Fl1 --" << std::endl;
				//for(int i=0; i<Nfl1; i++){
				//	f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
				//	//std::cout << "  " << scdder.deriv[i];	
				//}
				//std::cout << std::endl;
				//std::cout << "After scd der on fl1 ==> f=" << f << std::endl;	
	
				if(Nfl2 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0+Nfl1, Nfl2)); // The l=2 frequencies
				}
				//std::cout << "--- Fl2 --" << std::endl;
				for(int i=0; i<Nfl2; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
					//std::cout << "  " << scdder.deriv[i];
	
				}
				//std::cout << std::endl;
				//std::cout << "After scd der on fl2 ==> f=" << f << std::endl;	
	
				if(Nfl3 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0+Nfl1+Nfl2, Nfl3)); // The l=3 frequencies
				}
				//std::cout << "--- Fl3 --" << std::endl;
				for(int i=0; i<Nfl3; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
					//std::cout << "  " << scdder.deriv[i];
				}
				//std::cout << std::endl;
				//std::cout << "After scd der on fl3 ==> f=" << f << std::endl;	
			  	break;
	}

		
		
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
	//std::cout << "After securities ==> f=" << f << std::endl;	  
	
	/*std::cout << " ---- Noise slot ----" << std::endl;
	std::cout << "H_harvey=" << params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+3] << std::endl;
	std::cout << "White noise=" << params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+9] << std::endl;
	std::cout << "Inc=" << params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+10] << std::endl;
	
	exit(EXIT_SUCCESS);
	*/
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
			  //std::cout << "[" << i << "] FIX OR NONE, pena=" << pena << std::endl;
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


