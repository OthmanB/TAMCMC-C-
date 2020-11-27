# include <iostream>
# include <iomanip>
#include <fstream>
# include <Eigen/Dense>
#include "models.h"
#include "noise_models.h"
#include "../../external/ARMM/solver_mm.h"
#include "../../external/ARMM/bump_DP.h"
//#include <cmath>

using Eigen::VectorXd;
using Eigen::VectorXi;

double lin_interpol(VectorXd x, VectorXd y, double x_int);
VectorXd amplitude_ratio(int l, double beta);

VectorXd model_MS_Global_a1l_etaa3_HarveyLike(VectorXd params, VectorXi params_length, VectorXd x){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax.
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */
    
    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;

    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1(l=1),eta,a3, magb, magalfa, asym, a1(l=2))
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    
    double inclination;
    double trunc_c;
    
    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());
    
    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a11, a12,eta,a3, asym;
    
    int Nharvey;
    long cpt;
    
    //std::ofstream file_out;
    
    /*
     -------------------------------------------------------
     ------- Gathering information about the modes ---------
     -------------------------------------------------------
     */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    
	inclination=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise];
	
    // Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
    ratios_l0.setOnes();
    if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
        ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
    if(lmax >=2){
        Vl2=std::abs(params[Nmax+1]);
        ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;
        
    }
    if(lmax >=3){
        Vl3=std::abs(params[Nmax+2]);
        ratios_l3=amplitude_ratio(3, inclination);
    }
    
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
    Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);
    
    a11=std::abs(params[Nmax + lmax + Nf]);
    a12=std::abs(params[Nmax + lmax + Nf+6]);
    eta=params[Nmax + lmax + Nf + 1];
    a3=params[Nmax + lmax + Nf + 2];
    asym=params[Nmax+lmax + Nf + 5];

    
    model_final.setZero();
    
    //std::cout << "model_MS_Global_a1l_etaa3_HarveyLike" << std::endl;
    //std::cout << "a11=" << a11 << std::endl;
    //std::cout << "a12=" << a12 << std::endl;
    /* -------------------------------------------------------
     --------- Computing the models for the modes  ---------
     -------------------------------------------------------
     */
    cpt=0;
    for(long n=0; n<Nmax; n++){
        
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
			//std::cout << "[0] do conversion" << std::endl;
		} else{
			Hl0=std::abs(params[n]);
		}		
        model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl0, fl0, a11, a12, eta, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        
        if(lmax >=1){
 			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			if(do_amp){
				Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
				//std::cout << "[1] do conversion" << std::endl;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
           model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl1, fl1, a11, a12, eta, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        }
        if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
				//std::cout << "[2] do conversion" << std::endl;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl2, fl2, a11, a12, eta, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        }
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl3, fl3, a11, a12, eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        }
    }
    //std::cin.ignore();
    
    /* -------------------------------------------------------
     ------- Gathering information about the noise ---------
     -------------------------------------------------------
     */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
    
    /* -------------------------------------------------------
     ---------- Computing the mode of the noise ------------
     -------------------------------------------------------
     */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    
    return model_final;
}

VectorXd model_MS_Global_a1n_etaa3_HarveyLike(VectorXd params, VectorXi params_length, VectorXd x){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax.
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */
    
    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	    
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be Nmax+6 for a ALL global MS model (<a1>0=0,eta,a3, magb, magalfa, asym, a1(n))
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

    double inclination, trunc_c;
    
    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());
    
    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise), a11, a12; //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta,a3, asym, a1_l1, a1_l2;
    
    int Nharvey;
    long cpt;
    
    //std::ofstream file_out;
    
    /*
     -------------------------------------------------------
     ------- Gathering information about the modes ---------
     -------------------------------------------------------
     */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    inclination=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise];

   
    // Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
    ratios_l0.setOnes();
    if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
        ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
    if(lmax >=2){
        Vl2=std::abs(params[Nmax+1]);
        ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;
        
    }
    if(lmax >=3){
        Vl3=std::abs(params[Nmax+2]);
        ratios_l3=amplitude_ratio(3, inclination);
    }
    
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
    Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);
    
    a11=params.segment(Nmax + lmax + Nf + 6, Nmax);
    a12=a11;
    eta=params[Nmax + lmax + Nf + 1];
    a3=params[Nmax + lmax + Nf + 2];
    asym=params[Nmax+lmax + Nf + 5];
    
    model_final.setZero();
    
    //std::cout << "model_MS_Global_a1n_etaa3_HarveyLike" << std::endl;
    //std::cout << "a11=" << a11.transpose() << std::endl;
    //std::cout << "a12=" << a12.transpose() << std::endl;
    /* -------------------------------------------------------
     --------- Computing the models for the modes  ---------
     -------------------------------------------------------
     */
    cpt=0;
    for(long n=0; n<Nmax; n++){
        
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
			//std::cout << "[0] do conversion" << std::endl;
		} else{
			Hl0=std::abs(params[n]);
		}		
		a1_l1=std::abs(a11[n]);
		a1_l2=std::abs(a12[n]);
        model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl0, fl0, a1_l1, a1_l2, eta, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);

        if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			if(do_amp){
				Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
				//std::cout << "[1] do conversion" << std::endl;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl1, fl1, a1_l1, a1_l2, eta, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        }
        if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
				//std::cout << "[2] do conversion" << std::endl;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl2, fl2, a1_l1, a1_l2, eta, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        }
        //std::cout << "lmax=" << lmax << std::endl;
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl3, fl3, a1_l1, a1_l2, eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        }
    }
    //std::cin.ignore();
    
    /* -------------------------------------------------------
     ------- Gathering information about the noise ---------
     -------------------------------------------------------
     */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
    
    /* -------------------------------------------------------
     ---------- Computing the mode of the noise ------------
     -------------------------------------------------------
     */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    
    return model_final;
}


VectorXd model_MS_Global_a1nl_etaa3_HarveyLike(VectorXd params, VectorXi params_length, VectorXd x){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax.
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */
    
    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	    
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 2*Nmax+6 for a ALL global MS model (<a1>0=0,eta,a3, magb, magalfa, asym, a1(n,l))
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    
    double inclination, trunc_c;
    
    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());
    
    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise), a11, a12;; //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta,a3, asym, a1_l1, a1_l2;
    
    int Nharvey;
    long cpt;
    
    //std::ofstream file_out;
    
    /*
     -------------------------------------------------------
     ------- Gathering information about the modes ---------
     -------------------------------------------------------
     */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    inclination=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise];
  
  
    // Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
    ratios_l0.setOnes();
    if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
        ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
    if(lmax >=2){
        Vl2=std::abs(params[Nmax+1]);
        ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;
        
    }
    if(lmax >=3){
        Vl3=std::abs(params[Nmax+2]);
        ratios_l3=amplitude_ratio(3, inclination);
    }
    
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
    Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);
    
    a11=params.segment(Nmax + lmax + Nf + 6, Nmax);
    a12=params.segment(Nmax + lmax + Nf + 6 + Nmax, Nmax);
    eta=params[Nmax + lmax + Nf + 1];
    a3=params[Nmax + lmax + Nf + 2];
    
    asym=params[Nmax+lmax + Nf + 5];
 
    model_final.setZero();
    
    //std::cout << "model_MS_Global_a1nl_etaa3_HarveyLike" << std::endl;
    //std::cout << "a11=" << a11.transpose() << std::endl;
    //std::cout << "a12=" << a12.transpose() << std::endl;
    /* -------------------------------------------------------
     --------- Computing the models for the modes  ---------
     -------------------------------------------------------
     */
    cpt=0;
    for(long n=0; n<Nmax; n++){
        
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
			//std::cout << "[0] do conversion" << std::endl;
		} else{
			Hl0=std::abs(params[n]);
		}		
		a1_l1=std::abs(a11[n]);
		a1_l2=std::abs(a12[n]);
		
        model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl0, fl0, a1_l1, a1_l2, eta, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        
        if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			if(do_amp){
				Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
				//std::cout << "[1] do conversion" << std::endl;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl1, fl1, a1_l1, a1_l2, eta, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        }
        if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
				//std::cout << "[2] do conversion" << std::endl;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl2, fl2, a1_l1, a1_l2, eta, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        }
        if(lmax >=3){
 			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
           model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl3, fl3, a1_l1, a1_l2, eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        }
    }
    //std::cin.ignore();
    
    /* -------------------------------------------------------
     ------- Gathering information about the noise ---------
     -------------------------------------------------------
     */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
    
    /* -------------------------------------------------------
     ---------- Computing the mode of the noise ------------
     -------------------------------------------------------
     */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    
    return model_final;
}

VectorXd model_MS_Global_a1etaa3_HarveyLike(VectorXd params, VectorXi params_length, VectorXd x){
	/* Model of the power spectrum of a Main sequence solar-like star
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
	 *          Size MUST be 0 otherwise (this check is not made in this function)
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

	double inclination;
	

	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta,a3, asym;

	int Nharvey;
	long cpt;

	//std::ofstream file_out;

	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	
  	a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);
    
    inclination=std::atan(params[Nmax + lmax + Nf+4]/params[Nmax + lmax + Nf+3]); 
    inclination=inclination*180./pi;

    //std::cout << "a1 = " << a1 << std::endl;
    //std::cout << "inclination = " << inclination << std::endl;

	// Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
	ratios_l0.setOnes();
	if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
		ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
	if(lmax >=2){
		Vl2=std::abs(params[Nmax+1]);
		ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;

	}
	if(lmax >=3){
		Vl3=std::abs(params[Nmax+2]);
		ratios_l3=amplitude_ratio(3, inclination);
	}

	fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
	Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);

	//a1=std::abs(params[Nmax + lmax + Nf]);
	eta=params[Nmax + lmax + Nf + 1];
	a3=params[Nmax + lmax + Nf + 2];
	asym=params[Nmax+lmax + Nf + 5];
	
	model_final.setZero();
	
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	cpt=0;
	for(long n=0; n<Nmax; n++){
		
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);
			
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
		} else{
			Hl0=std::abs(params[n]);
		}		

		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);

		if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			if(do_amp){
				//Hl1=std::abs(params[n]/(pi*Wl1));
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
		}
		if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
                //Hl2=std::abs(params[n]/(pi*Wl2));
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
		}
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				//Hl3=std::abs(params[n]/(pi*Wl3));
                Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
		}		
	}
    //std::cin.ignore();

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
	//std::cout << "Before computing noise model" << std::endl;
	noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
	Nharvey=(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
		
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
	
	//std::cout << "End test" << std::endl;
    //exit(EXIT_SUCCESS);

	return model_final;
}



VectorXd model_MS_Global_a1etaa3_Harvey1985(VectorXd params, VectorXi params_length, VectorXd x){
	/* Model of the power spectrum of a Main sequence solar-like star
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
	 *          Size MUST be 0 otherwise (this check is not made in this function)
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

	double inclination, trunc_c;

	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta,a3, asym;

	int Nharvey;
	long cpt;

	//std::ofstream file_out;

	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];

	inclination=std::atan(params[Nmax + lmax + Nf+4]/params[Nmax + lmax + Nf+3]); 
	inclination=inclination*180./pi;
	a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);

	std::cout << "a1 = " << a1 << std::endl;
	std::cout << "inclination = " << inclination << std::endl;

	// Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
	ratios_l0.setOnes();
	if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
		ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
	if(lmax >=2){
		Vl2=std::abs(params[Nmax+1]);
		ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;

	}
	if(lmax >=3){
		Vl3=std::abs(params[Nmax+2]);
		ratios_l3=amplitude_ratio(3, inclination);
	}

	fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
	Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);

	//a1=std::abs(params[Nmax + lmax + Nf]);
	eta=params[Nmax + lmax + Nf + 1];
	a3=params[Nmax + lmax + Nf + 2];
	asym=params[Nmax+lmax + Nf + 5];
	std::cout << "a1=" << a1 << std::endl;
	std::cout << "eta=" << eta << std::endl;
	std::cout << "a3=" << a3 << std::endl;
	std::cout << "asym=" << asym << std::endl;
	exit(EXIT_SUCCESS);
	
	model_final.setZero();
	
	//std::cout << "Here 2 " << std::endl;
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	cpt=0;
	for(long n=0; n<Nmax; n++){
						
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);
			
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
			//std::cout << "[0] do conversion" << std::endl;
		} else{
			Hl0=std::abs(params[n]);
		}		
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);

		if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			if(do_amp){
				Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
				//std::cout << "[1] do conversion" << std::endl;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
		}
		if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
				//std::cout << "[2] do conversion" << std::endl;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
		}
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
		}		
	}
    //std::cin.ignore();

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
	noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
	//noise_params.array().abs();
	Nharvey=(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
		
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey1985(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background

	return model_final;
}

// Last modifcation: on 05/02/2020: Handling the possibility that the user asked for Heights Hlm to be fitted instead of inclination
// This might be a temporary fix for testing. If successfull, a full derivation into a specific model might be better
// Added on 05/02/2020: Handling the possibility that the user asked for Heights Hlm to be fitted instead of inclination
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic(VectorXd params, VectorXi params_length, VectorXd x){
	/* Model of the power spectrum of a Main sequence solar-like star
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
	 *          Size MUST be 0 otherwise (this check is not made in this function)
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

	double inclination, trunc_c;

	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta,a3, asym;

	int Nharvey;
	long cpt;

	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];

    inclination=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise]; // This version of the code take inclination as direct argument

	// Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
	ratios_l0.setOnes();
	if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
		ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
	if(lmax >=2){
		Vl2=std::abs(params[Nmax+1]);
		ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;
	}
	if(lmax >=3){
	   Vl3=std::abs(params[Nmax+2]);
	   ratios_l3=amplitude_ratio(3, inclination);
    }
 
	fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
	Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);

	a1=std::abs(params[Nmax + lmax + Nf]);
	eta=params[Nmax + lmax + Nf + 1];
	a3=params[Nmax + lmax + Nf + 2];
	asym=params[Nmax+lmax + Nf + 5];
	
	model_final.setZero();
	
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	cpt=0;
	for(long n=0; n<Nmax; n++){
				
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);	
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
			//std::cout << "[0] do conversion" << std::endl;
		} else{
			Hl0=std::abs(params[n]);
		}		
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
		
		if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			if(do_amp){
				Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
				//std::cout << "[1] do conversion" << std::endl;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
			//model_final=model_final + model_l1;
 		}
		if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
				//std::cout << "[2] do conversion" << std::endl;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
		}
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
		}		
	}

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
	//std::cout << "Before computing noise model" << std::endl;
	noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
	Nharvey=(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
		
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
	

	return model_final;
}


// Added on 10 Feb 2020
// Alternative version  of the legacy function that is not dealing with the inclination
// but instead fit the mode heights, given a visibility
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic_v2(VectorXd params, VectorXi params_length, VectorXd x){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */

    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Equal to the  number of heights for m components: Sum(l*(l+1))_l={1,2,3}
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

    double trunc_c;

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta,a3, asym;

    int Nharvey;
    long cpt;

    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
    if(lmax >=2){
        Vl2=std::abs(params[Nmax+1]);
        //std::cout << "Vl2 " << Vl2 << std::endl;

    }
    if(lmax >=3){
        Vl3=std::abs(params[Nmax+2]);
    }

    // l=1 Heights are defined by two parameters 
    ratios_l0[0]=1;

    ratios_l1[0]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+1]);  // m=-1
    ratios_l1[1]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise]); // m=0
    ratios_l1[2]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+1]); // m=+1 

    // l=2 Heights are defined by three parameters  
    ratios_l2[0]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+4]);  // m=-2
    ratios_l2[1]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+3]); // m=-1
    ratios_l2[2]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+2]); // m=0 
    ratios_l2[3]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+3]); // m=+1
    ratios_l2[4]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+4]); // m=+2 
 
    // l=3 Heights are defined by four parameters 
    ratios_l3[0]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+8]);  // m=-3
    ratios_l3[1]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+7]); // m=-2
    ratios_l3[2]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+6]); // m=-1 
    ratios_l3[3]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+5]);  // m=0
    ratios_l3[4]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+6]); // m=+1
    ratios_l3[5]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+7]); // m=+2 
    ratios_l3[6]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+8]); // m=+3
 
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
    Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);

    a1=std::abs(params[Nmax + lmax + Nf]);
    eta=params[Nmax + lmax + Nf + 1];
    a3=params[Nmax + lmax + Nf + 2];
    asym=params[Nmax+lmax + Nf + 5];
    
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */
    cpt=0;
    for(long n=0; n<Nmax; n++){
                
        fl0=fl0_all[n];
        Wl0=std::abs(Wl0_all[n]);   
        if(do_amp){
            Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
            //std::cout << "[0] do conversion" << std::endl;
        } else{
            Hl0=std::abs(params[n]);
        }       
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        
        if(lmax >=1){
            fl1=params[Nmax+lmax+Nfl0+n];
            Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
            if(do_amp){
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
                //std::cout << "[1] do conversion" << std::endl;
            } else{
                Hl1=std::abs(params[n]*Vl1);
            }               
            model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
            //model_final=model_final + model_l1;
        }
        if(lmax >=2){
            fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
            Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
            if(do_amp){
                Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
                //std::cout << "[2] do conversion" << std::endl;
            } else{
                Hl2=std::abs(params[n]*Vl2);
            }   
            model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        }
        if(lmax >=3){
            fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
            Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
            if(do_amp){
                Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
            } else{
                Hl3=std::abs(params[n]*Vl3);            
            }       
            model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        }     
    }

    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise backgound
    

    //exit(EXIT_SUCCESS);

    return model_final;
}


// Added on 10 Feb 2020
// Alternative version  of the legacy function that is not dealing with the inclination
// but instead fit the mode heights directly. Visibilities are not considered here.
// WARNING: This model is not really suitable for a global fit due to the large set of parameters that might
//          Arise from letting all heights free for all m E ([l*(l+1)-1]/2 + 1) * Ntotal ~ 100 free parameters for a global fit with l=3
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic_v3(VectorXd params, VectorXi params_length, VectorXd x){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */

    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Equal to the  number of heights for m components: Sum(l*(l+1))_l={1,2,3}
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

    double trunc_c;

    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
    VectorXd Hl0(1), Hl1(3), Hl2(5), Hl3(7);
    double fl0, fl1, fl2, fl3, Wl0, Wl1, Wl2, Wl3, a1,eta,a3, asym;
    int Nharvey;
    long cpt, pos0;

    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];

    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
    Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);

    a1=std::abs(params[Nmax + lmax + Nf]);
    eta=params[Nmax + lmax + Nf + 1];
    a3=params[Nmax + lmax + Nf + 2];
    asym=params[Nmax+lmax + Nf + 5];
    
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */
    cpt=0;
    for(long n=0; n<Nmax; n++){
                
        fl0=fl0_all[n];
        Wl0=std::abs(Wl0_all[n]);   
        if(do_amp){
            Hl0[0]=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
            //std::cout << "[0] do conversion" << std::endl;
        } else{
            Hl0[0]=std::abs(params[n]);
        }       
        model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl0, fl0, a1, eta, a3, asym, Wl0, 0, step, trunc_c);
 
        if(lmax >=1){
            fl1=params[Nmax+lmax+Nfl0+n];
            Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
            pos0=2*n;
            Hl1[0]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];
            Hl1[1]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0];
            Hl1[2]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];
            if(do_amp){
                //tmp=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0, 2);
                Hl1=Hl1/(pi*Wl1);
            } 
            Hl1=Hl1.cwiseAbs();
            model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl1, fl1, a1, eta, a3,asym, Wl1, 1, step, trunc_c);
        }
        if(lmax >=2){
            fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
            Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
            pos0=3*n;
            Hl2[0]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+2];
            Hl2[1]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];     
            Hl2[2]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0];
            Hl2[3]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];
            Hl2[4]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+2];
            if(do_amp){
                //Hl2=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0, 3)/(pi*Wl2);
                Hl2=Hl2/(pi*Wl2);
           }
            Hl2=Hl2.cwiseAbs();
            model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl2, fl2, a1, eta, a3,asym, Wl2, 2, step, trunc_c);
        }
        if(lmax >=3){
            fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
            Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
            pos0=4*n;
            Hl3[0]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+3];
            Hl3[1]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+2];
            Hl3[2]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];     
            Hl3[3]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0];
            Hl3[4]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];
            Hl3[5]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+2];
            Hl3[6]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+3];
            if(do_amp){
               //Hl3=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0, 4)/(pi*Wl3);
                Hl3=Hl3/(pi*Wl3);
            }
            Hl3=Hl3.cwiseAbs();
            model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl3, fl3, a1, eta, a3, asym, Wl3, 3, step, trunc_c);
        }  
 /*       std::cout << "[" << n << "]" << std::endl;
        std::cout << "fl0[" << n << "] = " << fl0 << std::endl;
        std::cout << "Wl0[" << n << "] = " << Wl0 << std::endl;
        std::cout << "Hl0[" << n << "] = " << Hl0 << std::endl;
         std::cout << "---------" << std::endl;
        std::cout << "fl1[" << n << "] = " << fl1 << std::endl;
        std::cout << "Wl1[" << n << "] = " << Wl1 << std::endl;
        std::cout << "Hl1[" << n << "] = " << Hl1 << std::endl;
        std::cout << "---------" << std::endl;
        std::cout << "fl2[" << n << "] = " << fl2 << std::endl;
        std::cout << "Wl2[" << n << "] = " << Wl2 << std::endl;
        std::cout << "Hl2[" << n << "] = " << Hl2 << std::endl;
        std::cout << "---------" << std::endl;
        if(lmax>=3){
            std::cout << "fl3[" << n << "] = " << fl3 << std::endl;
            std::cout << "Wl3[" << n << "] = " << Hl3 << std::endl;
            std::cout << "Hl3[" << n << "] = " << Hl3 << std::endl;
            std::cout << "---------" << std::endl;
        }
        std::cout << "a1 = " << Wl0 << std::endl;
        std::cout << "eta = " << Wl1 << std::endl;
        std::cout << "a3 = " << Wl2 << std::endl;

       for(long indd=0; indd<model_final.size(); indd++){
            if (std::isfinite(model_final[indd]) == false){
                std::cout << indd << " Not finite" << std::endl;
                std::cout << "      " << model_final[indd] << std::endl;
                exit(EXIT_SUCCESS);
            }
        }  
*/
    }
    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    
    //exit(EXIT_SUCCESS);

    return model_final;
}
// ------------------------------
// ------------------------------
// ------------------------------

VectorXd model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1(VectorXd params, VectorXi params_length, VectorXd x){
	/* Model of the power spectrum of a Main sequence solar-like star
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 * Width a following the Appourchaux et al. 2014, 566, 20 and Appourchaux et al. 2016, 595, C2 (Corrigendum) relation.
	 * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
	 *          Size MUST be 0 otherwise (this check is not made in this function)
	 * The v1 version has the particularity to not explicitly impose numax as a parameter. Instead, it is 
	 * calculated at each iteration by a weighted average of the frequencies, with the weight being the heights
	 * This model has therefore 5 parameters for the widths
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
	
	double inclination;

	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta,a3, asym;
	double numax, Htot, lnGamma0, lnLorentz;
	double e;
	
	int Nharvey;
	long cpt;
	
	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	inclination=std::atan(params[Nmax + lmax + Nf+4]/params[Nmax + lmax + Nf+3]); 
	inclination=inclination*180./pi;
	a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);

    //std::cout << "a1 = " << a1 << std::endl;
    //std::cout << "inclination = " << inclination << std::endl;

	// Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
	ratios_l0.setOnes();
	if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
		ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
	if(lmax >=2){
		Vl2=std::abs(params[Nmax+1]);
		ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;

	}
	if(lmax >=3){
		Vl3=std::abs(params[Nmax+2]);
		ratios_l3=amplitude_ratio(3, inclination);
	}

	fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
	
	// ----- Widths ----
	numax=0.;
	Htot=0.;
	for(long n=0; n<Nmax; n++){
		numax=numax+params[n]*params[Nmax + lmax + n]; // Adding up the l=0: H(l=0, n)*nu(l=0, n)
		Htot=Htot+params[n];
		if(lmax>=1){
			numax=numax+params[n]*Vl1*params[Nmax+lmax+Nfl0+n]; // Adding up the l=1: H(l=1, n)*nu(l=1, n)
			Htot=Htot+params[n]*Vl1;
		}
		if(lmax>=2){
			numax=numax+params[n]*Vl2*params[Nmax+lmax+Nfl0+Nfl1+n]; // Adding up the l=2: H(l=2, n)*nu(l=2, n)
			Htot=Htot+params[n]*Vl2;
		}
		if(lmax>=3){
			numax=numax+params[n]*Vl3*params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n]; // Adding up the l=3: H(l=3, n)*nu(l=3, n)
			Htot=Htot+params[n]*Vl3;
		}
	}
	numax=numax/Htot;
	
	std::cout << "numax :" << numax << std::endl;
	// -----------------
	//a1=std::abs(params[Nmax + lmax + Nf]);
	eta=params[Nmax + lmax + Nf + 1];
	a3=params[Nmax + lmax + Nf + 2];
	asym=params[Nmax+lmax + Nf + 5];
	
	model_final.setZero();
	
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	cpt=0;
	for(long n=0; n<Nmax; n++){
		fl0=fl0_all[n];
		lnGamma0=params[Nmax + lmax + Nf + Nsplit+1] * log(fl0/numax) + log(params[Nmax+lmax+Nf+Nsplit+2]);
		e=2.*log(fl0/params[Nmax+lmax+Nf+Nsplit+0]) / log(params[Nmax+lmax+Nf+Nsplit+3]/numax);
		lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+4])/(1. + pow(e,2));		
		Wl0=exp(lnGamma0 + lnLorentz);

		//std::cout << "Wl0=" << Wl0 << std::endl;
		 
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
		} else{
			Hl0=std::abs(params[n]);
		}		

		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
		
		if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+1] * log(fl1/numax) + log(params[Nmax+lmax+Nf+Nsplit+2]);
			e=2.*log(fl1/params[Nmax+lmax+Nf+Nsplit+0]) / log(params[Nmax+lmax+Nf+Nsplit+3]/numax);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+4])/(1. + pow(e,2));		
			Wl1=exp(lnGamma0 + lnLorentz);

			//std::cout << "Wl1=" << Wl1 << std::endl;

			if(do_amp){
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
		}
		if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+1] * log(fl2/numax) + log(params[Nmax+lmax+Nf+Nsplit+2]);
			e=2.*log(fl2/params[Nmax+lmax+Nf+Nsplit+0]) / log(params[Nmax+lmax+Nf+Nsplit+3]/numax);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+4])/(1. + pow(e,2));		
			Wl2=exp(lnGamma0 + lnLorentz);

			//std::cout << "Wl2=" << Wl2 << std::endl;

			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
		}
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+1] * log(fl3/numax) + log(params[Nmax+lmax+Nf+Nsplit+2]);
			e=2.*log(fl3/params[Nmax+lmax+Nf+Nsplit+0]) / log(params[Nmax+lmax+Nf+Nsplit+3]/numax);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+4])/(1. + pow(e,2));		
			Wl3=exp(lnGamma0 + lnLorentz);

			//std::cout << "Wl3=" << Wl3 << std::endl;

			if(do_amp){
              Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
		}		
	}
    //std::cin.ignore();

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
	//std::cout << "Before computing noise model" << std::endl;
	noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
	Nharvey=(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
		
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
	
	//std::cout << "End test" << std::endl;
       //exit(EXIT_SUCCESS);

	return model_final;
}


VectorXd model_MS_Global_a1etaa3_AppWidth_HarveyLike_v2(VectorXd params, VectorXi params_length, VectorXd x){
	/* Model of the power spectrum of a Main sequence solar-like star
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 * Width a following the Appourchaux et al. 2014, 566, 20 and Appourchaux et al. 2016, 595, C2 (Corrigendum) relation.
	 * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
	 *          Size MUST be 0 otherwise (this check is not made in this function)
	 * The v2 version has the particularity to explicitly use numax as a parameter. 
	 * This model has therefore 6 parameters for the widths contrary to the 5 parameters of the v1 version
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
	
	double inclination;
	

	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta,a3, asym;
	double Htot, lnGamma0, lnLorentz;
	double e;
	
	int Nharvey;
	long cpt;
	

	//std::ofstream file_out;

	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	inclination=std::atan(params[Nmax + lmax + Nf+4]/params[Nmax + lmax + Nf+3]); 
	inclination=inclination*180./pi;
	a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);

    //std::cout << "a1 = " << a1 << std::endl;
    //std::cout << "inclination = " << inclination << std::endl;

	// Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
	ratios_l0.setOnes();
	if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
		ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
	if(lmax >=2){
		Vl2=std::abs(params[Nmax+1]);
		ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;

	}
	if(lmax >=3){
		Vl3=std::abs(params[Nmax+2]);
		ratios_l3=amplitude_ratio(3, inclination);
	}

	fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
	
	// -----------------
	//a1=std::abs(params[Nmax + lmax + Nf]);
	eta=params[Nmax + lmax + Nf + 1];
	a3=params[Nmax + lmax + Nf + 2];
	asym=params[Nmax+lmax + Nf + 5];
	
	model_final.setZero();
	
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	cpt=0;
	for(long n=0; n<Nmax; n++){
		fl0=fl0_all[n];
		lnGamma0=params[Nmax + lmax + Nf + Nsplit+2] * log(fl0/params[Nmax + lmax + Nf + Nsplit + 0]) + log(params[Nmax+lmax+Nf+Nsplit+3]);
		e=2.*log(fl0/params[Nmax+lmax+Nf+Nsplit+1]) / log(params[Nmax+lmax+Nf+Nsplit+4]/params[Nmax + lmax +Nf + Nsplit + 0]);
		lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+5])/(1. + pow(e,2));		
		
		Wl0=exp(lnGamma0 + lnLorentz);

		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
		} else{
			Hl0=std::abs(params[n]);
		}		

		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);

		//std::cout << "Wl0=" << Wl0 << std::endl;
		
		if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+2] * log(fl1/params[Nmax + lmax + Nf + Nsplit + 0]) + log(params[Nmax+lmax+Nf+Nsplit+3]);
			e=2.*log(fl1/params[Nmax+lmax+Nf+Nsplit+1]) / log(params[Nmax+lmax+Nf+Nsplit+4]/params[Nmax + lmax +Nf + Nsplit + 0]);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+5])/(1. + pow(e,2));		
			Wl1=exp(lnGamma0 + lnLorentz);

			//std::cout << "Wl1=" << Wl1 << std::endl;

			if(do_amp){
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
		}
		if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+2] * log(fl2/params[Nmax + lmax + Nf + Nsplit + 0]) + log(params[Nmax+lmax+Nf+Nsplit+3]);
			e=2.*log(fl2/params[Nmax+lmax+Nf+Nsplit+1]) / log(params[Nmax+lmax+Nf+Nsplit+4]/params[Nmax + lmax +Nf + Nsplit + 0]);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+5])/(1. + pow(e,2));		
			Wl2=exp(lnGamma0 + lnLorentz);

			//std::cout << "Wl2=" << Wl2 << std::endl;

			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
		}
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+2] * log(fl3/params[Nmax + lmax + Nf + Nsplit + 0]) + log(params[Nmax+lmax+Nf+Nsplit+3]);
			e=2.*log(fl3/params[Nmax+lmax+Nf+Nsplit+1]) / log(params[Nmax+lmax+Nf+Nsplit+4]/params[Nmax + lmax +Nf + Nsplit + 0]);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+5])/(1. + pow(e,2));		
			Wl3=exp(lnGamma0 + lnLorentz);

			if(do_amp){
                Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
		}		
	}
    //std::cin.ignore();

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
	//std::cout << "Before computing noise model" << std::endl;
	noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
	Nharvey=(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
		
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
	
	//std::cout << "End test" << std::endl;
        //exit(EXIT_SUCCESS);

	return model_final;
}


// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// //////////////////////////////   Models for local fit \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

VectorXd model_MS_local_basic(VectorXd params, VectorXi params_length, VectorXd x){

	/* Model for a local fit of the power spectrum of a solar-like star
	 * Assumptions of constant splitting over the fit window make it more suitable for a fit of a MS star
	 * but RGB might be handled also provided that the user focus on single mode fitting (ie, avoids fitting groups of different l)
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	
	const int Nmax=params_length[0]; // Total Number of Heights
	const int Nvis=params_length[1]; // number of visibilities
	const int Nfl0=params_length[2]; // number of l=0 frequencies, heights and widths
	const int Nfl1=params_length[3]; // number of l=1 frequencies, heights and widths
	const int Nfl2=params_length[4]; // number of l=2 frequencies, heights and widths
	const int Nfl3=params_length[5]; // number of l=3 frequencies, heights and widths
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a MS_local (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // Total number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Might be only the white noise for a local model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const double trunc_c=params[Nmax+Nvis+Nf+Nsplit+Nwidth+Nnoise+Ninc];
	const bool do_amp=params[Nmax+Nvis+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
	
	double inclination;
	
	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta,a3, asym;

	int Nharvey;
	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	inclination=std::atan(params[Nmax + Nvis + Nf+4]/params[Nmax + Nvis + Nf+3]); 
	inclination=inclination*180./pi;
	a1=pow(params[Nmax + Nvis + Nf+3],2)+ pow(params[Nmax + Nvis + Nf+4],2);

    //std::cout << "a1 = " << a1 << std::endl;
    //std::cout << "inclination = " << inclination << std::endl;

	// Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
	ratios_l0.setOnes();
	if(Nfl1 >=1){
//        Vl1=std::abs(params[Nmax]);
		ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
	if(Nfl2 >=1){
//		Vl2=std::abs(params[Nmax+1]);
		ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;

	}
	if(Nfl3 >=1){
//		Vl3=std::abs(params[Nmax+2]);
		ratios_l3=amplitude_ratio(3, inclination);
	}

	//a1=std::abs(params[Nmax + Nvis + Nf]);
	eta=params[Nmax + Nvis + Nf + 1];
	a3=params[Nmax + Nvis + Nf + 2];
	asym=params[Nmax+Nvis + Nf + 5];
	
	model_final.setZero();
	
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	for(long n=0; n<Nfl0; n++){		
		fl0=params[Nmax + Nvis + n];
		Wl0=std::abs(params[Nmax + Nvis + Nf + Nsplit + n ]);		
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
		} else{
			Hl0=std::abs(params[n]);
		}		
		//std::cout << "fl0 = " << fl0 << "      Hl0 = " << Hl0 << "      Wl0 = " << Wl0 << std::endl;
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
	}
	for(long n=0; n<Nfl1; n++){		
		fl1=params[Nmax + Nvis + Nfl0 + n];
		Wl1=std::abs(params[Nmax + Nvis + Nf + Nsplit + Nfl0 + n ]);	
		if(do_amp){
    		Hl1=std::abs(params[Nfl0 + n]/(pi*Wl1));
		} else{
			Hl1=std::abs(params[Nfl0 + n]);
		}				
		//std::cout << "fl1 = " << fl1 << "      Hl1 = " << Hl1 << "      Wl1 = " << Wl1 << std::endl;
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
	}
	
	for(long n=0; n<Nfl2; n++){		
		fl2=params[Nmax + Nvis + Nfl0 + Nfl1 + n];
		Wl2=std::abs(params[Nmax+ Nvis + Nf + Nsplit + Nfl0 + Nfl1 + n ]);	
		if(do_amp){
			Hl2=std::abs(params[Nfl0 + Nfl1 + n]/(pi*Wl2));
		} else{
			Hl2=std::abs(params[Nfl0 + Nfl1 + n]);
		}	
		//std::cout << "fl2 = " << fl2 << "      Hl2 = " << Hl2 << "      Wl2 = " << Wl2 << std::endl;
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
	}
	for(long n=0; n<Nfl3; n++){		
        //std::cout << "BEFORE l=3" << std::endl;
		fl3=params[Nmax + Nvis + Nfl0 + Nfl1 + Nfl2 + n];
        //std::cout << "fl3 = " << fl3 << std::endl;
		Wl3=std::abs(params[Nmax + Nvis + Nf + Nsplit + Nfl0 + Nfl1 + Nfl2 + n ]);
        //std::cout << "Wl3 = " << Wl3 << std::endl;
		if(do_amp){
        	Hl3=std::abs(params[Nfl0 + Nfl1 + Nfl2 + n]/(pi*Wl3));
		} else{
			Hl3=std::abs(params[Nfl0 + Nfl1 + Nfl2 + n]);			
		}		
		//std::cout << "fl3 = " << fl3 << "      Hl3 = " << Hl3 << "      Wl3 = " << Wl3 << std::endl;
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        //std::cout << "AFTER l=3" << std::endl;
	}

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
    //std::cout << "Before noise_params" << std::endl;
	noise_params=params.segment(Nmax+Nvis+Nf+Nsplit+Nwidth, Nnoise);
	Nharvey=0; //(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
	//std::cout << "End test" << std::endl;
    //exit(EXIT_SUCCESS);

	return model_final;
}


VectorXd model_MS_local_Hnlm(VectorXd params, VectorXi params_length, VectorXd x){

    /* Model for a local fit of the power spectrum of a solar-like star
     * Assumptions of constant splitting over the fit window make it more suitable for a fit of a MS star
     * but RGB might be handled also provided that the user focus on single mode fitting (ie, avoids fitting groups of different l)
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * This model differs from model_MS_local_basic() by the fact that it fits the heights Hnlm instead than the inclination.
     * Therefore it compares to the model_MS_Global_a1etaa3_HarveyLike_Classic_v3() but for the local class of models
     * This means that the it preserve the symetry H(n,l,+m) = H(n,l, -m). It has to be used jointly with extra_priors[3]=2
     * Which imposes Sum(Hnlm)_{m=-l, m+l} = 1 in priors_local()
     */

    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    
    const int Nmax=params_length[0]; // Total Number of Heights
    const int Nvis=params_length[1]; // number of visibilities
    const int Nfl0=params_length[2]; // number of l=0 frequencies, heights and widths
    const int Nfl1=params_length[3]; // number of l=1 frequencies, heights and widths
    const int Nfl2=params_length[4]; // number of l=2 frequencies, heights and widths
    const int Nfl3=params_length[5]; // number of l=3 frequencies, heights and widths
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a MS_local (a1,eta,a3, magb, magalfa, asym)
    const int Nwidth=params_length[7]; // Total number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Might be only the white noise for a local model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+Nvis+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+Nvis+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    
    int pos0;
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd noise_params(Nnoise);
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Wl0, Wl1, Wl2, Wl3, a1,eta,a3, asym;
    VectorXd Hl0(1), Hl1(3), Hl2(5), Hl3(7);
    int Nharvey;
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    a1=std::abs(params[Nmax + Nvis + Nf]);
    eta=params[Nmax + Nvis + Nf + 1];
    a3=params[Nmax + Nvis + Nf + 2];
    asym=params[Nmax+Nvis + Nf + 5];
    
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */
    for(long n=0; n<Nfl0; n++){     
        fl0=params[Nmax + Nvis + n];
        Wl0=std::abs(params[Nmax + Nvis + Nf + Nsplit + n ]);       
        if(do_amp){
            Hl0[0]=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
        } else{
            Hl0[0]=std::abs(params[n]);
        }       
        //std::cout << "fl0 = " << fl0 << "      Hl0 = " << Hl0 << "      Wl0 = " << Wl0 << std::endl;
        model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl0, fl0, a1, eta, a3, asym, Wl0, 0, step, trunc_c);
    }
    for(long n=0; n<Nfl1; n++){     
        fl1=params[Nmax + Nvis + Nfl0 + n];
        Wl1=std::abs(params[Nmax + Nvis + Nf + Nsplit + Nfl0 + n ]);    
        pos0=2*n;
        Hl1[0]=params[Nfl0 + pos0+1];
        Hl1[1]=params[Nfl0 + pos0];
        Hl1[2]=params[Nfl0 + pos0+1];
        if(do_amp){
                Hl1=Hl1/(pi*Wl1);
        }
        Hl1=Hl1.cwiseAbs(); 
        model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl1, fl1, a1, eta, a3,asym, Wl1, 1, step, trunc_c);
    }    
    for(long n=0; n<Nfl2; n++){     
        fl2=params[Nmax + Nvis + Nfl0 + Nfl1 + n];
        Wl2=std::abs(params[Nmax+ Nvis + Nf + Nsplit + Nfl0 + Nfl1 + n ]);  
        pos0=3*n;
        Hl2[0]=params[Nfl0 + Nfl1 + pos0+2];
        Hl2[1]=params[Nfl0 + Nfl1 + pos0+1];     
        Hl2[2]=params[Nfl0 + Nfl1 + pos0];
        Hl2[3]=params[Nfl0 + Nfl1 + pos0+1];
        Hl2[4]=params[Nfl0 + Nfl1 + pos0+2];
        if(do_amp){
            Hl2=Hl2/(pi*Wl2);
        }
        Hl2=Hl2.cwiseAbs();
        model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl2, fl2, a1, eta, a3,asym, Wl2, 2, step, trunc_c);
    }
    for(long n=0; n<Nfl3; n++){     
        fl3=params[Nmax + Nvis + Nfl0 + Nfl1 + Nfl2 + n];
        Wl3=std::abs(params[Nmax + Nvis + Nf + Nsplit + Nfl0 + Nfl1 + Nfl2 + n ]);
        pos0=4*n;
        Hl3[0]=params[Nfl0 + Nfl1 + Nfl2 + pos0+3];
        Hl3[1]=params[Nfl0 + Nfl1 + Nfl2 + pos0+2];
        Hl3[2]=params[Nfl0 + Nfl1 + Nfl2 + pos0+1];     
        Hl3[3]=params[Nfl0 + Nfl1 + Nfl2 + pos0];
        Hl3[4]=params[Nfl0 + Nfl1 + Nfl2 + pos0+1];
        Hl3[5]=params[Nfl0 + Nfl1 + Nfl2 + pos0+2];
        Hl3[6]=params[Nfl0 + Nfl1 + Nfl2 + pos0+3];
        if(do_amp){
            Hl3=Hl3/(pi*Wl3);
       }
        Hl3=Hl3.cwiseAbs();
        model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl3, fl3, a1, eta, a3, asym, Wl3, 3, step, trunc_c);
    }

    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    noise_params=params.segment(Nmax+Nvis+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=0; //(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    
    //std::cout << "End test" << std::endl;
    //exit(EXIT_SUCCESS);

    return model_final;
}


////////////////////////////////////
// Models for asymptotic fitting ///
////////////////////////////////////

VectorXd model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v2(VectorXd params, VectorXi params_length, VectorXd x){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Width a following the Appourchaux et al. 2014, 566, 20 and Appourchaux et al. 2016, 595, C2 (Corrigendum) relation.
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     * The v2 version has the particularity to explicitly impose numax as a parameter. This is because we cannot create a v1
     * version unless we include so iterative process to evaluate numax from a weighted average of the frequencies. 
     * This would lead to uncessary complications and increase model computation time
     * This model has therefore 5 parameters for the widths. 
     * NOTE THAT IN THE RGB_asympt model, the v2 version is not implemented
     */

    const double step=x[2]-x[1]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1, ie, visibilities
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of parameters to describe the l=1 mixed modes: delta0l, DPl, alpha_g, q, sigma_p, sigma_g, sigma_m
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of zones to describe the rotation profile (2 for rot_core and rot_renv + 1 eta + 1 a3 + 1 asym )
    const int Nwidth=params_length[7]; // number of parameters for the l=0 widths. Should be 5 here as it uses the Appourchaux profile
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    
    double inclination;

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), Hl0_all(Nmax), noise_params(Nnoise), fl1_all, Wl1_all, Hl1p_all, Hl1_all,a1_l1, a1_l2(Nfl2), a1_l3(Nfl3); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta,a3, asym;
    double numax, Htot, lnGamma0, lnLorentz;
    double e;
    
    int Nharvey;
    long cpt;
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    inclination=params[Nmax + lmax + Nf+Nsplit + Nwidth + Nnoise]; 
    std::cout << "inclination =" << inclination << std::endl;
    //inclination=std::atan(params[Nmax + lmax + Nf+4]/params[Nmax + lmax + Nf+3]); 
    //inclination=inclination*180./pi;
    //a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);

    // Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
    ratios_l0.setOnes();
    if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
        ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
    if(lmax >=2){
        Vl2=std::abs(params[Nmax+1]);
        ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;

    }
    if(lmax >=3){
        Vl3=std::abs(params[Nmax+2]);
        ratios_l3=amplitude_ratio(3, inclination);
    }

    // --- Preparing profiles for l=0 modes ---
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths and mixed modes determination 
 
    for (int n=0; n<Nmax;n++)
    {
        lnGamma0=params[Nmax + lmax + Nf + Nsplit+2] * log(fl0_all[n]/params[Nmax + lmax + Nf + Nsplit + 0]) + log(params[Nmax+lmax+Nf+Nsplit+3]);
        e=2.*log(fl0_all[n]/params[Nmax+lmax+Nf+Nsplit+1]) / log(params[Nmax+lmax+Nf+Nsplit+4]/params[Nmax + lmax +Nf + Nsplit + 0]);
        lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+5])/(1. + pow(e,2));      
        Wl0_all[n]=exp(lnGamma0 + lnLorentz);
    }

   if(do_amp){
        Hl0_all=params.segment(0, Nmax).cwiseProduct(Wl0_all.cwiseInverse() / pi); // A^2/(pi.Gamma)
        Hl0_all=Hl0_all.cwiseAbs();
    } else{
        Hl0_all=params.segment(0, Nmax).cwiseAbs();
    }       

    // --------------
    // --------------
    // --------------
    // --------------
    // ---- Mixed modes handling ----
    // --------------
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen_m(seed); 
  
    const double delta0l=params[Nmax + lmax + Nfl0];
    const double DPl=params[Nmax + lmax + Nfl0 + 1];
    const double alpha_g=params[Nmax + lmax + Nfl0 + 2];
    const double q_star=params[Nmax + lmax + Nfl0 + 3];
    const double sigma_p_l1=params[Nmax + lmax + Nfl0 + 4];
    const double sigma_g_l1=params[Nmax + lmax + Nfl0 + 5];
    const double sigma_m_l1=params[Nmax + lmax + Nfl0 + 6];
    const double rot_env=params[Nmax + lmax + Nf];
    const double rot_core=params[Nmax + lmax + Nf+1];

    // -------- DEBUG OF l=1 mixed modes -----
    std::cout << "delta0l =" << delta0l << std::endl;
    std::cout << "DPl =" << DPl << std::endl;
    std::cout << "alpha_g =" << alpha_g << std::endl;
    std::cout << "q_star =" << q_star << std::endl;
    std::cout << "sigma_p =" << sigma_p_l1 << std::endl;
    std::cout << "sigma_g =" << sigma_g_l1 << std::endl;
    std::cout << "sigma_m =" << sigma_m_l1 << std::endl;
    std::cout << "rot_env =" << rot_env << std::endl;
    std::cout << "rot_core =" << rot_core << std::endl;
    // ---------------------------------------

    std::normal_distribution<double> distrib_m(0.,sigma_m_l1);

    double r, tmp;
    VectorXi posOK;  
    VectorXd ksi_pg, h1_h0_ratio;
    Data_eigensols freqs_l1;
    
    freqs_l1=solve_mm_asymptotic_O2from_l0(fl0_all, 1, delta0l, DPl, alpha_g, q_star, sigma_p_l1, step, true, false); // note that we use the true data resolution (step) for optimising computation
    // Filter solutions that endup at frequencies higher/lower than the nu_l0 because we will need to extrapolate height/widths otherwise...
    posOK=where_in_range(freqs_l1.nu_m, fl0_all.minCoeff(), fl0_all.maxCoeff(), false);
    fl1_all.resize(posOK.size());
    
    for (int i=0; i<posOK.size();i++)
    {
        fl1_all[i]=freqs_l1.nu_m[posOK[i]];
        if (sigma_m_l1 !=0) // If requested, we add a random gaussian qty to the mixed mode solution
        {
            r = distrib_m(gen_m);
            fl1_all[i]=fl1_all[i]+r;
        }
    }
    
    // Generating widths profiles for l=1 modes using the ksi function
    ksi_pg=ksi_fct2(fl1_all, freqs_l1.nu_p, freqs_l1.nu_g, freqs_l1.dnup, freqs_l1.dPg, q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
    h1_h0_ratio=h_l_rgb(ksi_pg); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019)
    
    Hl1p_all.resize(fl1_all.size());
    for (int i=0; i<fl1_all.size();i++)
    {
        tmp=lin_interpol(fl0_all, params.segment(0, Nmax), fl1_all[i]); // interpolate Hl0 to fl1 positions
        Hl1p_all[i]=tmp;
    }
    Hl1_all=h1_h0_ratio.cwiseProduct(Hl1p_all*Vl1);
    Wl1_all=gamma_l_fct2(ksi_pg, fl1_all, fl0_all, Wl0_all, h1_h0_ratio, 1); // generate the mixed modes widths
    // Generating splittings with a two-zone averaged rotation rates
    a1_l1=dnu_rot_2zones(ksi_pg, rot_env, rot_core);
    a1_l2.setConstant(rot_env);
    a1_l3.setConstant(rot_env);
    eta=params[Nmax + lmax + Nf + 2];
    a3=params[Nmax + lmax + Nf + 3];
    asym=params[Nmax+lmax + Nf + 4];
  
    // --------------
    // --------------
    // --------------
    // --------------
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */
    cpt=0;
    for(long n=0; n<Nfl0; n++){
        fl0=fl0_all[n];
        Wl0=Wl0_all[n];
        Hl0=Hl0_all[n];
        //std::cout << "Wl0=" << Wl0 << std::endl;
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, 0, eta, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
    }
    for (long n=0; n<freqs_l1.nu_m.size(); n++){
        fl1=freqs_l1.nu_m[n];
        Wl1=Wl1_all[n];
        Hl1=Hl1_all[n];              
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1_l1[n], eta, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
    }
    for(long n=0; n<Nfl2; n++){ 
        fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
        lnGamma0=params[Nmax + lmax + Nf + Nsplit+2] * log(fl2/params[Nmax + lmax + Nf + Nsplit + 0]) + log(params[Nmax+lmax+Nf+Nsplit+3]);
        e=2.*log(fl2/params[Nmax+lmax+Nf+Nsplit+1]) / log(params[Nmax+lmax+Nf+Nsplit+4]/params[Nmax + lmax +Nf + Nsplit + 0]);
        lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+5])/(1. + pow(e,2));      
        Wl2=exp(lnGamma0 + lnLorentz);
        //std::cout << "Wl2=" << Wl2 << std::endl;
        if(do_amp){
            Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
        } else{
            Hl2=std::abs(params[n]*Vl2);
        }   
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1_l2[n], eta, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
    }
    for(long n=0; n<Nfl3; n++){ 
        fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
        lnGamma0=params[Nmax + lmax + Nf + Nsplit+2] * log(fl3/params[Nmax + lmax + Nf + Nsplit + 0]) + log(params[Nmax+lmax+Nf+Nsplit+3]);
        e=2.*log(fl3/params[Nmax+lmax+Nf+Nsplit+1]) / log(params[Nmax+lmax+Nf+Nsplit+4]/params[Nmax + lmax +Nf + Nsplit + 0]);
        lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+5])/(1. + pow(e,2));      
        Wl3=exp(lnGamma0 + lnLorentz);
        //std::cout << "Wl3=" << Wl3 << std::endl;
        if(do_amp){
            Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
        } else{
            Hl3=std::abs(params[n]*Vl3);            
        }       
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1_l3[n], eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
    }       
    //std::cin.ignore();

    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    std::cout << "fl0    / Wl0     / Hl0" << std::endl;
    for (int i =0; i<fl0_all.size(); i++){
        std::cout << fl0_all[i]  << "    "  << Wl0_all[i]  << "    " << Hl0_all[i] << std::endl;
    }
    std::cout << "fl1    / Wl1     / Hl1    / a1_l1" << std::endl;
    for (int i =0; i<fl1_all.size(); i++){
        std::cout << fl1_all[i]  << "    "  << Wl1_all[i]  << "    " << Hl1_all[i]  << "    " << a1_l1[i]  << std::endl;
    }
    std::cout << "a1_l2 = " << a1_l2.transpose() << std::endl;
    
    std::cout << "a1_l3 = " << a1_l3.transpose() << std::endl;

    std::cout << "End test" << std::endl;
    exit(EXIT_SUCCESS);

    return model_final;
}


/*
void asymptotic_mm_fromtemplate(VectorXd params, VectorXi params_length, VectorXd x){

    int seed=(unsigned)time(NULL);
    srand(seed);

    std::random_device rd;
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> distrib(0 , 1);
    
    // ----- Constants ------
    const double PI = 3.141592653589793238462643;
    const double G=6.667e-8;
    const double Teff_sun=5777;
    const double Dnu_sun=135.1;
    const double numax_sun=3150;
    const double R_sun=6.96342e5;
    const double M_sun=1.98855e30;
    const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);
    // ----------------------
    const std::string cpath=getcwd(NULL, 0);
    
    const double lmax=3;
    const double Nmax_pm=6; // Number of radial order to inject on each side of numax
    double nmax_spread=params[16];
    
    double xmin, xmax, a1;
    double H, tau, p, N0;
    
    Cfg_synthetic_star cfg_star;
    Params_synthetic_star params_out;
    // ------------------------------------
    
    // -----------
    cfg_star.Teff_star=-1;
    
    cfg_star.numax_star=params[0];
    cfg_star.maxHNR_l0=params[1];
    cfg_star.Gamma_max_l0=params[2];
    cfg_star.H0_spread=params[3];
    
    cfg_star.Dnu_star=params[4];
    cfg_star.epsilon_star=params[5];
    cfg_star.delta0l_percent_star=params[6];
    cfg_star.beta_p_star=params[7];
    
    cfg_star.DPl_star=params[8];                
    cfg_star.alpha_g_star=params[9];
    cfg_star.q_star=params[10];

    cfg_star.rot_env_input=params[11];
    cfg_star.rot_ratio_input=-1;
    cfg_star.rot_core_input=params[12];
 
    cfg_star.Vl.resize(3);
    cfg_star.Vl << params[13], params[14], params[15];
 
    cfg_star.noise_params_harvey_like.resize(8);
    cfg_star.noise_params_harvey_like <<  params[19], params[20] , params[21] , params[22] , params[23] , params[24] , params[25]  , params[26];    //[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0]
    
    cfg_star.fmin=cfg_star.numax_star -Nmax_pm*cfg_star.Dnu_star;
    cfg_star.fmax=cfg_star.numax_star +(Nmax_pm+2)*cfg_star.Dnu_star;
    cfg_star.output_file_rot=cpath + "/external/ARMM/star_params.rot";
    cfg_star.filetemplate = template_file;
    cfg_star.nmax_star=cfg_star.numax_star/cfg_star.Dnu_star - cfg_star.epsilon_star;
    cfg_star.alpha_p_star=cfg_star.beta_p_star/cfg_star.nmax_star;
    cfg_star.sigma_m=params[17];
    cfg_star.sigma_p=params[18];

    if (std::abs(nmax_spread) > 0)
    {
        try
        {
            xmin=cfg_star.nmax_star*(1. - std::abs(nmax_spread)/100.);
            xmax=cfg_star.nmax_star*(1. + std::abs(nmax_spread)/100.);
            cfg_star.nmax_star=xmin + (xmax-xmin)*distrib(gen);
            
        }
        catch(...)
        {
            std::cout << "Error debug info:" << std::endl;
            std::cout << "cfg_star.nmax: " << cfg_star.nmax_star << std::endl;
            std::cout << "nmax_spread: " << nmax_spread << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    //inc_star=inc_rad*180./PI;
    
    // b. Generate the mode profiles and frequencies
    params_out=make_synthetic_asymptotic_star(cfg_star);

    if (cfg_star.Dnu_star <= 15){
        std::cout << "    Model with small Dnu ==> many mixed modes. This might be long to find the solutions..." << std::endl; 
    }

    tau=params[20] * pow(cfg_star.numax_star*1e-6,params[21]) + params[22]; // Granulation timescale (in seconds)
    H=params[17] * pow(cfg_star.numax_star*1e-6,params[18]) + params[19]; // Granulation Amplitude
    H=H/tau ; //This is due to the used definition for the Harvey profile (conversion from Hz to microHz)
    tau=tau/1000. ; //conversion in ksec
    p=params[23];// power law:  MUST BE CLOSE TO 2
    N0=params[24];
}
*/


VectorXd model_Harvey_Gaussian(VectorXd params, VectorXi params_length, VectorXd x){
/*
 * A model in which we fit a Gaussian + Harvey-Like profile.
 * Parameters are assumed to be in that order: Maximum Height, variance/width, central frequency, White noise
*/
	VectorXd nu0(x.size()), model_final(x.size()), noise_params;
	int Nharvey;

	// ------ Setting the Gaussian -------
	model_final= -0.5 * (x - nu0.setConstant(params[2])).array().square() /pow(std::abs(params[1]),2);
	model_final= std::abs(params[0])*model_final.array().exp();
	// ----------------------------------

	// ---- Setting the Noise model -----
	Nharvey=1;
	noise_params=params.segment(3, 4); // pick the 3 elements, begining from the index 3
	//model_final=harvey1985(noise_params.array().abs(), x, model_final, Nharvey);
	model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey);
	// ----------------------------------

	//std::cout << noise_params << std::endl;
	//exit(EXIT_SUCCESS);

	return model_final;
}

VectorXd model_Harvey1985_Gaussian(VectorXd params, VectorXi params_length, VectorXd x){
/*
 * A model in which we fit a Gaussian + Harvey-Like profile.
 * Parameters are assumed to be in that order: Maximum Height, variance/width, central frequency, White noise
*/
	VectorXd nu0(x.size()), model_final(x.size()), noise_params;
	int Nharvey;

	// ------ Setting the Gaussian -------
	model_final= -0.5 * (x - nu0.setConstant(params[2])).array().square() /pow(std::abs(params[1]),2);
	model_final= std::abs(params[0])*model_final.array().exp();
	// ----------------------------------

	// ---- Setting the Noise model -----
	Nharvey=1;
	noise_params=params.segment(3, 4); // pick the 3 elements, begining from the index 3
	model_final=harvey1985(noise_params.array().abs(), x, model_final, Nharvey);
	//model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey);
	// ----------------------------------

	//std::cout << noise_params << std::endl;
	//exit(EXIT_SUCCESS);

	return model_final;
}


VectorXd model_Test_Gaussian(VectorXd params, VectorXi params_length, VectorXd x){
/*
 * A model in which we fit a Gaussian + White noise.
 * Parameters are assumed to be in that order: Maximum Height, variance/width, central frequency, White noise
*/
	VectorXd nu0(x.size()), model_final(x.size()), tmp(x.size());

	tmp.setConstant(params[3]);

	model_final= -0.5 * (x - nu0.setConstant(params[2])).array().square() /pow(params[1],2);
	model_final= params[0]*model_final.array().exp();
	model_final= model_final + tmp;
	return model_final;
}


