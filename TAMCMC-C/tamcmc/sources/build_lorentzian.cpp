/*
 * build_lorentzian.cpp
 *
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>
#include "build_lorentzian.h"
# include <iostream>
# include <iomanip>

using Eigen::VectorXd;

VectorXd build_l_mode_a1l_etaa3(const VectorXd x_l, double H_l, double fc_l, double f_s1, double f_s2, double eta, double a3, double asym, double gamma_l, const int l, VectorXd V){
    /*
     * This model IS WITHOUT THE ASSUMPTION S11=S22. It includes:
     *      - Asymetry of Lorentzian asym
     *      - splitting a1(l=1) and a1(l=2). ASSUMES a1(l=3) =  (a1(1) + a1(2))/2.
     *      - an Asphericity parameter eta
     *      - latitudinal effect a3(l=2) only. We consider a3(l=3)=0
     */
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double Qlm, clm, f_s;
    
    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
            if(l == 1){
                clm=0; // a3=0 for l=1 BY DEFINITION
                f_s=f_s1;
            }
            if(l == 2){
                clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
                f_s=f_s2;
            }
            if(l == 3){
                clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
                f_s=(f_s1 + f_s2)/2.; // APPROXIMATION
            }
            profile=(x_l - tmp.setConstant(fc_l*(1. + eta*Qlm) + m*f_s + clm*a3)).array().square();
            profile=4*profile/pow(gamma_l,2);
        } else{
            profile=(x_l - tmp.setConstant(fc_l)).array().square();
            profile=4*profile/pow(gamma_l,2);
        }
        if(asym == 0){ //Model with no asymetry
            result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
        } else{
            tmp.setConstant(1);
            asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
            result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
        }
    }
    
    return result;
}

VectorXd build_l_mode_a1etaa3(const VectorXd x_l, double H_l, double fc_l, double f_s, double eta, double a3, double asym, double gamma_l, const int l, VectorXd V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - an Asphericity parameter eta
 *      - latitudinal effect a3
*/
	const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
	double Qlm, clm;

	result.setZero();
	for(int m=-l; m<=l; m++){
		if(l != 0){
			Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
			if(l == 1){
				clm=m; // a3 for l=1
			}
			if(l == 2){
				clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
			}
			if(l == 3){
				clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
			}
			profile=(x_l - tmp.setConstant(fc_l*(1. + eta*Qlm) + m*f_s + clm*a3)).array().square();
			profile=4*profile/pow(gamma_l,2);
		} else{
			profile=(x_l - tmp.setConstant(fc_l)).array().square();
			profile=4*profile/pow(gamma_l,2);
		}
		if(asym == 0){ //Model with no asymetry
			result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
		} else{
			tmp.setConstant(1);
			asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
			result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
		}
	}

return result;
}

VectorXd build_l_mode_asym_act(const VectorXd x_l, double H_l, double fc_l, double f_s, double eta, double a3, double b, double alpha, double asym, double gamma_l, const int l, VectorXd V){
/*
 * ---- OBSELETE FUNCTION ----
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - centrifugal force effect eta.nu.a1^2 (rotation-induced oblateness)
 *      - latitudinal effect a3
 *      - effect of magnetic field of the form b.nu^alpha
 * ---------------------------
*/
	const long Nxl=x_l.size();
        VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
	double Qlm, clm;

	result.setZero();
	for(int m=-l; m<=l; m++){
		if(l != 0){
			Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
			if(l == 1){
				clm=m; // a3 for l=1
			}
			if(l == 2){
				clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
			}
			if(l == 3){
				clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
			}
		
			profile=(x_l - tmp.setConstant(fc_l + Qlm*(eta*fc_l*pow(f_s,2) + b*pow(fc_l*1e-3,alpha)) + m*f_s + clm*a3) ).array().square();
			profile=4*profile/pow(gamma_l,2);
		} else{
			profile=(x_l - tmp.setConstant(fc_l)).array().square();
			profile=4*profile/pow(gamma_l,2);
		}
		if(asym == 0){ //Model with no asymetry
			result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
		} else{
			tmp.setConstant(1);
			asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
			//std::cout << "asymetry=" << asymetry.transpose() <<std::endl;
			//tmp2=((tmp.setConstant(1) + profile)).cwiseInverse();
			result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
			//std::cout << "in build_lorentzian" <<std::endl;
		}
	}

//std::cout << "result=" << result.transpose() <<std::endl;
return result;
}

VectorXd optimum_lorentzian_calc_a1l_etaa3(const VectorXd x, VectorXd y, double H_l, double fc_l, double f_s1, double f_s2, double eta, double a3, double asym, double gamma_l, const int l, VectorXd V, double step, const double c){
    /*
     function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
     that contains the lorentzian model.
     BEWARE: USES build_l_mode_a1l_etaa3() ==> Asphericity is a linear term in nu
     */
    //const double c=20.;
    double pmin, pmax, f_s;
    long imin, imax;
    VectorXd m0, x_l;
    //VectorXd mall(x.size());
    switch(l){
        case 0:
            f_s=0.;
        case 1:
            f_s=f_s1;
        case 2:
            f_s=f_s2;
        case 3:
            f_s=(f_s1 + f_s2)/2.;
    }
    if(gamma_l >= 1 && f_s >= 1){
        if(l != 0){
            pmin=fc_l - c*(l*f_s + gamma_l);
            pmax=fc_l + c*(l*f_s + gamma_l);
        } else{
            pmin=fc_l - c*gamma_l*2.2;
            pmax=fc_l + c*gamma_l*2.2;
        }
    }
    if(gamma_l <= 1 && f_s >= 1){
        if(l !=0){
            pmin=fc_l -c*(l*f_s + 1);
            pmax=fc_l + c*(l*f_s + 1);
        } else{
            pmin=fc_l - c*2.2;
            pmax=fc_l + c*2.2;
        }
    }
    if(gamma_l >= 1 && f_s <= 1){
        if(l != 0){
            pmin=fc_l - c*(l + gamma_l);
            pmax=fc_l + c*(l + gamma_l);
        } else{
            pmin=fc_l - c*2.2*gamma_l;
            pmax=fc_l + c*2.2*gamma_l;
        }
    }
    if(gamma_l <= 1 && f_s <= 1){
        if(l !=0){
            pmin=fc_l - c*(l+1);
            pmax=fc_l + c*(l+1);
        } else{
            pmin=fc_l -c*2.2;
            pmax=fc_l +c*2.2;
        }
    }
    // ---------- Handling boundaries ----------
    //     case when the proposed values lead to (pmax - step) < min(x)....
    if( (pmax - step) < x.head(1)(0)){
        pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
    }
    //     case when the proposed values lead to (pmin + step) >= max(x)...
    if( (pmin + step) >= x.tail(1)(0)){
        pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
    }
    
    imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
    imax=ceil((pmax-x.head(1)(0))/step);
    
    if(imin < 0){imin=0;} // ensure that we always stay within the vector x
    if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x
    if(imax-imin <= 0){
        std::cout << "Warning imax -imin <= 0 : imin=" << imin << "   imax=" << imax << std::endl;
        std::cout << " - pmin=" << pmin << "   pmax=" << pmax << std::endl;
        std::cout << " - step=" << step << std::endl;
        std::cout << " --------" << std::endl;
        std::cout << " - l=" << l << std::endl;
        std::cout << " - fc_l=" << fc_l << std::endl;
        std::cout << " - H_l=" << H_l << std::endl;
        std::cout << " - gamma_l=" << gamma_l << std::endl;
        std::cout << " - f_s=" << f_s << std::endl;
        std::cout << " - eta=" << eta << std::endl;
        std::cout << " - a3=" << a3 << std::endl;
        std::cout << " - asym=" << asym << std::endl;
        std::cout << " --------" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
    //std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
    
    //std::cout << "step=" << step << std::endl;
    //std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
    //std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();
    
    x_l=x.segment(imin, imax-imin);
    
    m0=build_l_mode_a1l_etaa3(x_l, H_l, fc_l, f_s1, f_s2, eta, a3, asym, gamma_l, l, V);
    //mall.setZero();
    //mall.segment(imin, imax-imin)=m0;
    y.segment(imin, imax-imin)= y.segment(imin, imax-imin) + m0;
    return y;
}


VectorXd optimum_lorentzian_calc_a1etaa3(const VectorXd x, VectorXd y, double H_l, double fc_l, double f_s, double eta, double a3, double asym, double gamma_l, const int l, VectorXd V, double step, const double c){
/*
	function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
	that contains the lorentzian model.
	BEWARE: USES build_l_mode_a1etaa3() ==> Asphericity is a linear term in nu
*/
	//const double c=20.;
	double pmin, pmax;
	long imin, imax;
	VectorXd m0, x_l;
	//VectorXd mall(x.size());

	if(gamma_l >= 1 && f_s >= 1){
		if(l != 0){
			pmin=fc_l - c*(l*f_s + gamma_l);
			pmax=fc_l + c*(l*f_s + gamma_l);
		} else{
			pmin=fc_l - c*gamma_l*2.2;
			pmax=fc_l + c*gamma_l*2.2;
		}
	}
	if(gamma_l <= 1 && f_s >= 1){
		if(l !=0){
			pmin=fc_l -c*(l*f_s + 1);
			pmax=fc_l + c*(l*f_s + 1);
		} else{
			pmin=fc_l - c*2.2;
			pmax=fc_l + c*2.2;
		}
	}
	if(gamma_l >= 1 && f_s <= 1){
		if(l != 0){
			pmin=fc_l - c*(l + gamma_l);
			pmax=fc_l + c*(l + gamma_l);
		} else{
			pmin=fc_l - c*2.2*gamma_l;
			pmax=fc_l + c*2.2*gamma_l;
		}
	}
	if(gamma_l <= 1 && f_s <= 1){
		if(l !=0){
			pmin=fc_l - c*(l+1);
			pmax=fc_l + c*(l+1);
		} else{
			pmin=fc_l -c*2.2;
			pmax=fc_l +c*2.2;
		}
	}
	// ---------- Handling boundaries ----------
	//     case when the proposed values lead to (pmax - step) < min(x)....
	if( (pmax - step) < x.head(1)(0)){ 
		pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
	}
	//     case when the proposed values lead to (pmin + step) >= max(x)...
	if( (pmin + step) >= x.tail(1)(0)){
		pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
	}

	imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
	imax=ceil((pmax-x.head(1)(0))/step);

	if(imin < 0){imin=0;} // ensure that we always stay within the vector x
	if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x
	if(imax-imin <= 0){
		std::cout << "Warning imax -imin <= 0 : imin=" << imin << "   imax=" << imax << std::endl;
		std::cout << " - pmin=" << pmin << "   pmax=" << pmax << std::endl;
		std::cout << " - step=" << step << std::endl;
		std::cout << " --------" << std::endl;
		std::cout << " - l=" << l << std::endl;
		std::cout << " - fc_l=" << fc_l << std::endl;
		std::cout << " - H_l=" << H_l << std::endl;
		std::cout << " - gamma_l=" << gamma_l << std::endl;
		std::cout << " - f_s=" << f_s << std::endl;
		std::cout << " - eta=" << eta << std::endl;
		std::cout << " - a3=" << a3 << std::endl;
		std::cout << " - asym=" << asym << std::endl;
		std::cout << " --------" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	//std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
	//std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
	
	//std::cout << "step=" << step << std::endl;
	//std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
	//std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();

	x_l=x.segment(imin, imax-imin);
 
	m0=build_l_mode_a1etaa3(x_l, H_l, fc_l, f_s, eta, a3, asym, gamma_l, l, V);
	//mall.setZero();
	//mall.segment(imin, imax-imin)=m0;
    	y.segment(imin, imax-imin)= y.segment(imin, imax-imin) + m0;
return y;
}

VectorXd optimum_lorentzian_calc_a1acta3(const VectorXd x, VectorXd y, double H_l, double fc_l, double f_s, double eta, double a3, 
		double b, double alpha, double asym, double gamma_l, const int l, VectorXd V, double step, const double c){
/*
	function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
	that contains the lorentzian model.
	BEWARE: USES build_l_mode_asym_act() ==> Include mode asymetry and effect of activity

*/
	//const double c=20.;
	double pmin, pmax;
	long imin, imax;
	VectorXd m0, x_l;
	//VectorXd mall(x.size());

	if(gamma_l >= 1 && f_s >= 1){
		if(l != 0){
			pmin=fc_l - c*(l*f_s + gamma_l);
			pmax=fc_l + c*(l*f_s + gamma_l);
		} else{
			pmin=fc_l - c*gamma_l*2.2;
			pmax=fc_l + c*gamma_l*2.2;
		}
	}
	if(gamma_l <= 1 && f_s >= 1){
		if(l !=0){
			pmin=fc_l -c*(l*f_s + 1);
			pmax=fc_l + c*(l*f_s + 1);
		} else{
			pmin=fc_l - c*2.2;
			pmax=fc_l + c*2.2;
		}
	}
	if(gamma_l >= 1 && f_s <= 1){
		if(l != 0){
			pmin=fc_l - c*(l + gamma_l);
			pmax=fc_l + c*(l + gamma_l);
		} else{
			pmin=fc_l - c*2.2*gamma_l;
			pmax=fc_l + c*2.2*gamma_l;
		}
	}
	if(gamma_l <= 1 && f_s <= 1){
		if(l !=0){
			pmin=fc_l - c*(l+1);
			pmax=fc_l + c*(l+1);
		} else{
			pmin=fc_l -c*2.2;
			pmax=fc_l +c*2.2;
		}
	}
	// ---------- Handling boundaries ----------
	//     case when the proposed values lead to (pmax - step) < min(x)....
	if( (pmax - step) < x.head(1)(0)){ 
		pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
	}
	//     case when the proposed values lead to (pmin + step) >= max(x)...
	if( (pmin + step) >= x.tail(1)(0)){
		pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
	}

	imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
	imax=ceil((pmax-x.head(1)(0))/step);

	if(imin < 0){imin=0;} // ensure that we always stay within the vector x
	if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x
	if(imax-imin <= 0){
		std::cout << "Warning imax -imin <= 0 : imin=" << imin << "   imax=" << imax << std::endl;
		std::cout << " - pmin=" << pmin << "   pmax=" << pmax << std::endl;
		std::cout << " - step=" << step << std::endl;
		std::cout << " --------" << std::endl;
		std::cout << " - l=" << l << std::endl;
		std::cout << " - fc_l=" << fc_l << std::endl;
		std::cout << " - H_l=" << H_l << std::endl;
		std::cout << " - gamma_l=" << gamma_l << std::endl;
		std::cout << " - f_s=" << f_s << std::endl;
		std::cout << " - eta=" << eta << std::endl;
		std::cout << " - a3=" << a3 << std::endl;
		std::cout << " - b=" << b << std::endl;
		std::cout << " - alpha=" << alpha << std::endl;
		std::cout << " - asym=" << asym << std::endl;
		std::cout << " --------" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	//std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
	//std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
	
	//std::cout << "step=" << step << std::endl;
	//std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
	//std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();

	x_l=x.segment(imin, imax-imin);
 
	//m0=build_l_mode(x_l, H_l, fc_l, f_s, a2, a3, gamma_l, l, V);
	m0=build_l_mode_asym_act(x_l, H_l, fc_l, f_s, eta, a3, b, alpha, asym, gamma_l, l, V);
	//mall.setZero();
	//mall.segment(imin, imax-imin)=m0;
        y.segment(imin, imax-imin)= y.segment(imin, imax-imin) + m0;
return y;
}


