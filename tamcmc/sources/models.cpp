# include <iostream>
# include <iomanip>
#include <fstream>
# include <Eigen/Dense>
#include "models.h"
#include "noise_models.h"
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

	//std::ofstream file_out;

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
		fl3=params[Nmax + Nvis + Nfl0 + Nfl1 + Nfl2 + n];
		Wl3=std::abs(params[Nmax + Nvis + Nf + Nsplit + Nfl0 + Nfl1 + Nfl2 + n ]);
		if(do_amp){
        	Hl3=std::abs(params[Nfl0 + Nfl1 + Nfl2 + n]/(pi*Wl3));
		} else{
			Hl3=std::abs(params[Nfl0 + Nfl1 + Nfl2 + n]);			
		}		
		//std::cout << "fl3 = " << fl3 << "      Hl3 = " << Hl3 << "      Wl3 = " << Wl3 << std::endl;
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
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


