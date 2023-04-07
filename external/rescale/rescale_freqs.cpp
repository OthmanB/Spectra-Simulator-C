#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include "../../linspace.h"
#include "../../linfit.h"
#include "decompose_nu.h"
#include "rescale_freqs.h"
#include "data.h"
using Eigen::VectorXd;

Freq_modes rescale_freqs(const double Dnu_star, const double epsilon_star, const Freq_modes freqs_ref, const VectorXd& d0l_star){
    /*
        A rescaling function decompose the frequency by identifying all of the terms of the asymptotic
        in order to isolate the residual error term. Then it uses this residual error (with a proper rescaling) to generate a new set 
        of frequencies following the asymptotic with the desired Dnu_star, epsilon_star, dl1_star (d01), dl2_star (d02), dl3_star (d13).
        Finally, a consistency check is made by re-decomposing the rescaled frequency so that you can check that the rescaling process
        worked as expected
    */
   // Declarations
   const bool verbose = false;
   const double Cfactor = 0.25;
   int l;
   double Dnu_ref, epsilon_ref, n0_ref, n0_star, e;
   VectorXd tmp, n_ref, n_star, fcoefs, fcoefs_star;
   Freq_modes freqs_star;
   Data_asympt_p asymp_l1,asymp_l2,asymp_l3;
   freqs_star.error_status=false; // By default, no error
   // Initial check
   if(freqs_ref.fl0.size() == 0){
        std::cout << "Warning: You must at least set freqs_ref.fl0 !" << std::endl;
        std::cout << "         Nothing to do" << std::endl;
        freqs_star.fl0.resize(1); freqs_star.fl0.setConstant(-1);
        freqs_star.fl1.resize(1); freqs_star.fl1.setConstant(-1);
        freqs_star.fl2.resize(1); freqs_star.fl2.setConstant(-1);
        freqs_star.fl3.resize(1); freqs_star.fl3.setConstant(-1);
        return freqs_star;
   }
    // Compute Dnu and epsilon
    n_ref=linspace(0, freqs_ref.fl0.size()-1, freqs_ref.fl0.size());
    fcoefs=linfit(n_ref, freqs_ref.fl0);
    Dnu_ref=fcoefs[0];
    epsilon_ref=modf(fcoefs[1]/fcoefs[0], &n0_ref);

    // Rescaling and epsilon shifting l=0
    tmp.resize(freqs_ref.fl0.size());
    tmp.setConstant(- epsilon_ref + epsilon_star);
    freqs_star.fl0=(freqs_ref.fl0/Dnu_ref + tmp) * Dnu_star;
    n_star=linspace(0, freqs_star.fl0.size()-1, freqs_star.fl0.size()); // First array without knowing n0
    fcoefs_star=linfit(n_star, freqs_star.fl0);
    //Dnu_l=fcoefs_star[0];
    e=modf(fcoefs_star[1]/fcoefs_star[0], &n0_star);
    n_star=linspace(n0_star, n0_star + freqs_star.fl0.size()-1, freqs_star.fl0.size()); // Recompute n_star now that we know n0  
    // l=1
    l=1;
    if(freqs_ref.fl1.size() != 0){
        if (d0l_star[0] <= -9999){
            std::cout << "Error: fl1 and dl1_star must be jointly set! Otherwise do not pass them as arguments" <<  std::endl;
            std::cout << "       Cannot rescale" << std::endl;
            freqs_star.error_status=true;
            //exit(EXIT_FAILURE);
        }
        // Rescaling the fl_ref:
        //       1. Extract all of asymtptotic elements. In particular nl and O2_l
        asymp_l1=decompose_nu_nl(1, freqs_ref.fl0, freqs_ref.fl1, Cfactor, verbose);
        //       2. Reconstruct a new relation using the requested Dnu_star and epsilon_star along with nl and O2_l*Dnu_star/Dnu_ref
        tmp.resize(asymp_l1.n.size()); tmp.setConstant((epsilon_star + 1.*l/2)*Dnu_star + d0l_star[0]);
        freqs_star.fl1=asymp_l1.n*Dnu_star + tmp + asymp_l1.O2_term*Dnu_star/Dnu_ref;
    } else{
        freqs_star.fl1.resize(1);
        freqs_star.fl1.setConstant(-1);
    }
    // l=2
    l=2;
    if(freqs_ref.fl2.size() != 0){
        if (d0l_star[1] <= -9999){
            std::cerr << "Error: fl2 and dl2_star must be jointly set! Otherwise do not pass them as arguments" <<  std::endl;
            std::cerr << "       Cannot rescale" << std::endl;
            freqs_star.error_status=true;            
            //exit(EXIT_FAILURE);
        }
        // Rescaling the fl_ref:
        //       1. Extract all of asymtptotic elements. In particular nl and O2_l
        asymp_l2=decompose_nu_nl(2, freqs_ref.fl0, freqs_ref.fl2, Cfactor, verbose);
        //       2. Reconstruct a new relation using the requested Dnu_star and epsilon_star along with nl and O2_l*Dnu_star/Dnu_ref
        tmp.resize(asymp_l2.n.size()); tmp.setConstant((epsilon_star + 1.*l/2)*Dnu_star + d0l_star[1]);
        freqs_star.fl2=asymp_l2.n*Dnu_star + tmp + asymp_l2.O2_term*Dnu_star/Dnu_ref;
    } else{
        freqs_star.fl2.resize(1);
        freqs_star.fl2.setConstant(-1);
    }
    // l=3
    l=3;
    if(freqs_ref.fl3.size() != 0){
        if (d0l_star[2] <= -9999){
            std::cerr << "Error: fl2 and dl2_star must be jointly set! Otherwise do not pass them as arguments" <<  std::endl;
            std::cerr << "       Cannot rescale" << std::endl;
           freqs_star.error_status=true;            
            //exit(EXIT_FAILURE);
        }
        // Rescaling the fl_ref:
        //       1. Extract all of asymtptotic elements. In particular nl and O2_l
        asymp_l3=decompose_nu_nl(3, freqs_ref.fl0, freqs_ref.fl3, Cfactor, verbose);
        //       2. Reconstruct a new relation using the requested Dnu_star and epsilon_star along with nl and O2_l*Dnu_star/Dnu_ref
        tmp.resize(asymp_l3.n.size()); tmp.setConstant((epsilon_star + 1.*l/2)*Dnu_star + d0l_star[2]);
        freqs_star.fl3=asymp_l3.n*Dnu_star + tmp + asymp_l3.O2_term*Dnu_star/Dnu_ref;
    } else{
        freqs_star.fl3.resize(1);
        freqs_star.fl3.setConstant(-1);
    }
    return freqs_star;
}