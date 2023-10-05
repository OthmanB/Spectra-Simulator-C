// -------------------
// ---- Functions adapted from the solver_mm.py function ----
/* This contains all the functions for solving the asymptotic relation
# of the mixed modes, as they have been tested during their development
# All this arise from reading few papers from Benoit Mosser and 
# The PhD thesis from Charlotte Gehand:
# https://arxiv.org/pdf/1203.0689.pdf (Mosser paper on mixed modes)
# https://arxiv.org/pdf/1004.0449.pdf (older Mosser paper on scaling relations for gaussian_width, Amp etc.. - 2010 -)
# https://arxiv.org/pdf/1011.1928.pdf (The universal pattern introduced with the curvature - Fig. 3 - )
# https://arxiv.org/pdf/1411.1082.pdf
# https://tel.archives-ouvertes.fr/tel-02128409/document

# Examples and tests function have been built using asymptotic relations in the python original code.
# But note that they should be applicable to ANY set of value of:
#	 nu_p(nu), nu_g(nu), Dnu_p(nu) and DPl(nu) (meaning handling glitches)
*/
// ------------------
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif
//#include "solver_mm.h"

#include "data_solver.h"
#include "string_handler.h"
#include "derivatives_handler.h"
#include "interpol.h"
#include "linfit.h"

using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::VectorXd;

Eigen::VectorXd removeDuplicates(const Eigen::VectorXd& nu_m_all, double tol) {
    Eigen::VectorXd uniqueVec;
    
    for(int i = 0; i < nu_m_all.size(); i++) {
        bool isDuplicate = false;
        
        for(int j = 0; j < uniqueVec.size(); j++) {
            if(std::abs(nu_m_all[i] - uniqueVec[j]) <= tol) {
                isDuplicate = true;
                break;
            }
        }
        
        if(!isDuplicate) {
            uniqueVec.conservativeResize(uniqueVec.size() + 1);
            uniqueVec[uniqueVec.size() - 1] = nu_m_all[i];
        }
    }
    
    return uniqueVec;
}

// Function that detects sign changes
// If the sign went from + to - tag it with a -1
// If the sign went from - to + tag it with a +1
// If there is no change of sign tag it with a 0
// 0 is dealt as a zone of change of sign as well. eg. if we pass from 0 to 2, then the result is +1
// Inputs:
//    - x: input vector for which we want to know sign changes
//    - return_indices: if true (default), the function returns positions at which the sign changed
//						if false, it returns a vector of size(x)-1 with the tags for the sign changes (or not)
VectorXi sign_change(const VectorXd& x, bool return_indices=true)
{
	//bool bool_tmp;
	VectorXi s(x.size()-1), pos_s(x.size()-1);
	s.setConstant(0); // Vector of tags for sign changes
	pos_s.setConstant(-1); // Vector of indices
	long i, j=0; 
	for (i = 0; i<x.size()-1; i++)
	{
		if ((   x[i+1]>=0 && x[i] >=0   )|| (  x[i+1]<=0 && x[i]   )) // No sign change case or values are constant at 0
		{
			s[i]=0;
		}
		if (  (x[i+1]>=0 && x[i] <0) || (x[i+1]>0 && x[i] <=0)  ) // Sign change case from - to +
		{
			s[i]=1;
			pos_s[j]=i;
			j=j+1;
		}
		if (  (x[i+1]<=0 && x[i] >0) || (x[i+1]<=0 && x[i] >=0)   ) // Sign change case from + to -
		{
			s[i]=-1;
			pos_s[j]=i;
			j=j+1;
		}
	}
	pos_s.conservativeResize(j); // Discard slots that where not used in the vector
	if (return_indices == true)
	{
		return pos_s;
	} else
	{
		return s;
	}
}


VectorXd pnu_fct(const VectorXd& nu, const long double nu_p)
{
	if (nu.size() == 0){
		std::cout << "Vector nu in pnu_fct is of size 0. Cannot pursue." << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		VectorXd tmp(nu.size()), pnu(nu.size());
		tmp.setConstant(nu_p);
		pnu=nu - tmp;
		return pnu;
	}
}

long double pnu_fct(const long double nu, const long double nu_p)
{
	long double pnu=nu-nu_p;
	return pnu;
}

VectorXd gnu_fct(const VectorXd& nu, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q)
{
	const long double pi = 3.141592653589793238L;
	VectorXd X(nu.size()), gnu(nu.size()), tmp(nu.size());
	tmp.setConstant(nu_g);

	X=pi * (nu.cwiseInverse() - tmp.cwiseInverse())*1e6 / DPl;
	tmp=q*X.array().tan(); //.tan();
	gnu=Dnu_p * tmp.array().atan()/pi; //.atan() / pi;
	return gnu;
}

long double gnu_fct(const long double nu, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q)
{
	const long double pi = 3.141592653589793238L;
	long double X, gnu;
	
	X=pi * (1./nu - 1./nu_g)*1e6 / DPl;
	gnu=Dnu_p * atan(q*tan(X))/pi; //.atan() / pi;
	return gnu;
}

/* A small function that generate a serie of p modes using the asymptotic relation
# at the second order as per defined in Mosser et al. 2018, equation 22 (https://www.aanda.org/articles/aa/pdf/2018/10/aa32777-18.pdf)
# delta0l and alpha and nmax must be set, a
# Note that we have the following relationship between D0 and delta0l:
#			delta0l=-l(l+1) D0 / Dnu_p
# Such that delta0l=-l(l+1) gamma / 100, if gamma is in % of Dnu_p
# r: Allows you to add an extra term to Dnu_o
*/
long double asympt_nu_p(const long double Dnu_p, const int np, const long double epsilon, const int l, 
	const long double delta0l, const long double alpha, const long double nmax, long double r=0)
{

/*
	std::cout << " np=" << np << std::endl;
	std::cout << " epsilon=" << epsilon << std::endl;
	std::cout << " l=" << l << std::endl;
	std::cout << " delta0l=" << delta0l << std::endl;
	std::cout << " alpha=" << alpha << std::endl;
	std::cout << " nmax=" << nmax << std::endl;
	std::cout << " Dnu_p=" << Dnu_p << std::endl;
*/	
	long double nu_p=(np + epsilon + l/2. + delta0l + alpha*std::pow(np - nmax, 2) / 2)*Dnu_p;
	if (nu_p < 0.0)
	{
		std::cout << " WARNING: NEGATIVE FREQUENCIES DETECTED: IMPOSING POSITIVITY" << std::endl;
		std::cout << " nu_p: " << nu_p << std::endl;
		std::cout << " Cannot pursue " << std::endl;
		exit(EXIT_FAILURE);
	}
	return nu_p+r;
}


/* A small function that generate a serie of p modes based on a shifting of a series of l=0 modes
# and on the asymptotic relation. This effectively allow to account for 2nd order terms of p modes
# delta0l : small spacing 
# r: Allows you to add an extra term to Dnu_p
# WARNING: COULD BE SOME PROBLEMS ON THE EDGES... IF PROBLEM ARE FOUND, WE MIGHT NEED REPLACE THIS ALGO BY 
#          AN INTERPOLATION AT NP... BUT MORE COSTLY
*/

long double asympt_nu_p_from_l0(const VectorXd& nu_l0, const long double Dnu_p, const int np, const long double epsilon, const int l, 
	const long double delta0l, long double r=0)
{
	long double nu_p_l;

	const int np_min0=int(floor(nu_l0.minCoeff()/Dnu_p - epsilon));
	const int np_max0=int(floor(nu_l0.maxCoeff()/Dnu_p - epsilon));

	if (np >= np_min0 && np < np_max0){ // If we are inside the nu_l0 range, we make the shifting
		nu_p_l= nu_l0[np-np_min0] + l/2.*Dnu_p + delta0l;
	} else{
		nu_p_l=(np + epsilon + l/2.)*Dnu_p + delta0l;
	}
	if (nu_p_l < 0.0)
	{
		std::cout << " WARNING: NEGATIVE FREQUENCIES DETECTED: IMPOSING to nu_p_l + r=0" << std::endl;
		std::cout << " nu_p_l: " << nu_p_l << std::endl;
		 nu_p_l=0;
		 r=0;
	}
	return nu_p_l+r;
}

VectorXd asympt_nu_p_from_l0_Xd(const VectorXd& nu_l0, const long double Dnu_p, const int l, const long double delta0l, long double fmin=-1, long double fmax=-1)
{
	VectorXd nu_l0_long, nu_l1_long, nu_l1, tmp;
	VectorXi posOK;

	if (fmin == -1){
		fmin=0;
	}
	if(fmax == -1){
		fmax=nu_l0.maxCoeff() + 10*Dnu_p;
	}
	// Extend the vector by extrapolating edges... this to avoid bad behavior at the edges (missing frequencies)
	nu_l0_long.resize(nu_l0.size() + 6);
	nu_l1_long.resize(nu_l0.size() + 6);

	nu_l0_long[0]=nu_l0.minCoeff() - 3*Dnu_p;
	nu_l0_long[1]=nu_l0.minCoeff() - 2*Dnu_p;
	nu_l0_long[2]=nu_l0.minCoeff() - Dnu_p;
	for (int k=0; k<nu_l0.size();k++){
		nu_l0_long[k+3]=nu_l0[k];
	}
	nu_l0_long[nu_l0_long.size()-3]=nu_l0.maxCoeff() + Dnu_p;
	nu_l0_long[nu_l0_long.size()-2]=nu_l0.maxCoeff() + 2*Dnu_p;
	nu_l0_long[nu_l0_long.size()-1]=nu_l0.maxCoeff() + 3*Dnu_p;

	tmp.resize(nu_l0_long.size());
	tmp.setConstant(l/2.*Dnu_p + delta0l);
	nu_l1_long=nu_l0_long + tmp;

    posOK=where_in_range(nu_l1_long, fmin, fmax, 0); // Remove frequencies out of the requested range
    nu_l1.resize(posOK.size());
    if (posOK[0] != -1){
    	for (int k=0;k<posOK.size(); k++){
    		nu_l1[k]=nu_l1_long[posOK[k]];
    	}
    }
	else{
		std::cout << "Serious issue when computing frequencies with solver_mm::asympt_nu_p_from_l0_Xd()" << std::endl;
		std::cout << "No frequency found in the specified range (" << fmin << "," << fmax << "). Debug required" << std::endl;
		exit(EXIT_FAILURE);
	}
	return nu_l1;
}

/* 
	Compute the asymptotic relation for the g modes.
	r: an optional parameter that can be added to the Period (e.g. a random quantity)	
*/
long double asympt_nu_g(const long double DPl, const int ng, const long double alpha, long double r=0)
{
	const long double Pl=(ng + alpha)*DPl;
	return 1e6/(Pl+r);
}

/*
This the main function that solves the mixed mode asymptotic relation
which is of the type p(nu) = g(nu)
This solver specifically solve the case:
      nu - nu_p = Dnu*arctan(q tan(1/(nu DPl) - 1/(nu_g*DPl)))
      It tries to find the intersect between p(nu) and g(nu)
      using an interpolation and over a range of frequency such that nu is in [numin, numax]
Parameters:
	- Mandatory: 
	     nu_p (double) : frequency of a given p-mode (in microHz)
	     nu_g (double): frequency of a given g-mode (in microHz)
	     Dnu_p (double): Large separation for p modes (in microHz)
	     DP1 (double): Period spacing for g modes (in seconds)
	     q (double): Coupling term (no unit, should be between 0 and 1)
	- Optional:
		numin (double): Minimum frequency considered for the solution (in microHz)
		numax (double): Maximum frequency considered for the solution (in microHz)
		resol (double): Base resolution for the interpolated base function. The interpolation may miss solutions 
		       if this is set too low. Typically, the resolution parameter should be higher than the
		       spectral resolution of your spectrum so that all resolved modes should be found.
		       This is also used for creating the nu axis for visualisation (in microHz).
		factor (double): Define how fine will be the new tiny grid used for performing the interpolation. This is important
				to avoid extrapolation (which is forbiden and will result in crash of the code). Typically, the default
				value factor=0.05 can compute mixed modes for frequency down to 80microHz. Going below requires a smaller factor

		returns_axis: If True, returns nu, pnu and gnu (see optional reutrns below). Mainly for debug
Returns a structure with:
	nu_m: An array with all solutions that match p(nu) = g(nu)
	nu (optional): The frequency axis used as reference for finding the intersection
	pnu (optional): The curve for p(nu)
	gnu (optional): The curve g(nu)
*/
Data_coresolver solver_mm(const long double nu_p, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q, 
	const long double numin, const long double numax, const long double resol, const bool returns_axis=false, const bool verbose=false, const long double factor=0.05)
{
	//const int Nmmax=500; // Number of maximum mixed modes solutions that can be found 

	int i, s_ok;
	long double range_min, range_max, nu_m_proposed, ratio,  ysol_pnu, ysol_gnu;
	Data_coresolver results;
	VectorXi idx;
	VectorXd nu, pnu, gnu, nu_local, pnu_local, gnu_local, nu_m, ysol_all;

	if (nu_g >= numin && nu_g <= numax){
	// Generate a frequency axis that has a fixed resolution and that span from numin to numax
	nu = Eigen::VectorXd::LinSpaced(long((numax-numin)/resol), numin, numax);

	// Function p(nu) describing the p modes
	pnu=pnu_fct(nu, nu_p);
	// Function g(nu) describing the g modes 
	gnu=gnu_fct(nu, nu_g, Dnu_p, DPl, q);
	
	/* Find when p(nu) = g(nu) by looking for solution of p(nu) - g(nu) = 0
	#     Method 1: Direct Interpolation... Works only for single solutions ==> Not used here
	#int_fct = interpolate.interp1d(pnu - gnu, nu)
	#nu_m=int_fct(0)
	#     Method 2: (a) Find indices close to sign changes for p(nu) - g(nu)
	#               (b) Then perform an iterative interpolation in narrow ranges
	#                   near the approximate solutions. How narrow is the range is defined
	#					by the resolution parameter resol, which in this case can be view
	#					as the minimum precision.
	*/
	idx=sign_change(pnu-gnu);
	s_ok=0;
	nu_m.resize(idx.size());
	for (long ind=0; ind<idx.size();ind++)
	{
		//std::cout << "idx[" << ind << "] =" << idx[ind] << std::endl;
		// Define a small local range around each of the best solutions
		range_min=nu[idx[ind]] - 2*resol;
		range_max=nu[idx[ind]] + 2*resol;
		// Redefine nu, pnu and gnu for that local range
		
		nu_local = Eigen::VectorXd::LinSpaced(long((range_max-range_min)/(resol*factor)), range_min, range_max);
		pnu_local=pnu_fct(nu_local, nu_p);
		gnu_local=gnu_fct(nu_local, nu_g, Dnu_p, DPl, q);	

		// Perform the interpolation on the local range and append the solution to the nu_m list
		nu_m_proposed=lin_interpol(pnu_local - gnu_local, nu_local, 0);
		try
		{	
			ysol_gnu=gnu_fct(nu_m_proposed, nu_g, Dnu_p, DPl, q);
			ysol_pnu=pnu_fct(nu_m_proposed, nu_p);
		}
		catch(...)
		{
			std::cout << "Interpolation issue detected. Debuging information:" << std::endl;
			std::cout << "    nu_p: " <<  nu_p << std::endl;
			std::cout << "    nu_g: " <<  nu_g << std::endl;
			std::cout << "    Dnu_p: "<< Dnu_p << std::endl;
			std::cout << "    DPl: "<< DPl << std::endl;
			std::cout << "    q: " << q << std::endl;
			std::cout << "    numin: "<< numin << std::endl;
			std::cout << "    numax: "<< numax << std::endl;
			std::cout << "    resol:"<< resol << std::endl;
			std::cout << "    factor:"<< factor << std::endl;
			std::cout << " ------------" << std::endl;
			std::cout << "range_min/max: "<< range_min << range_max << std::endl;
			std::cout << "  nu_local: "<< nu_local << std::endl;
			std::cout << "  pnu_local: "<< pnu_local << std::endl;
			std::cout << "  gnu_local: "<< gnu_local << std::endl;
			std::cout << " ------------" << std::endl;
			std::cout << " int_fct  ==>  nu_local      /   pnu_local - gnu_local : " << std::endl;
			for (i=0; i<nu_local.size(); i++)
			{
				std::cout << "    " <<  nu_local[i]<< pnu_local[i]-gnu_local[i] << std::endl;
			}
			exit(EXIT_FAILURE);
		}
		ratio=ysol_gnu/ysol_pnu;
		if (verbose == true)
		{	std::cout << "-------"<< std::endl;
			std::cout << "nu_m:"<<  nu_m_proposed<< std::endl;
			std::cout << "Ratio:"<<  ratio << std::endl;
		}
		// Sometimes, the interpolator mess up due to the limits of validity for the atan function
		// The way to keep real intersection is to verify after interpolation that we really
		// have p(nu_m_proposed) = g(nu_m_proposed). We then only keeps solutions that satisfy
		// a precision criteria of 0.1%.
		if ((ratio >= 0.999) && (ratio <= 1.001))
		{
			nu_m[s_ok]=nu_m_proposed;
			s_ok=s_ok+1;
		}
	}
	nu_m.conservativeResize(s_ok);
	ysol_all=gnu_fct(nu_m, nu_g, Dnu_p, DPl, q);
	}
	if (returns_axis == true){
		results.nu_m=nu_m;
		results.ysol=ysol_all;
		results.nu=nu;
		results.pnu=pnu;
		results.gnu=gnu;
		return results;
	}
	else{
		results.nu_m=nu_m;
		results.ysol=ysol_all;
		return results;
	}
}

Data_eigensols solve_mm_asymptotic_O2p(const long double Dnu_p, const long double epsilon, const int el, const long double delta0l, const long double alpha_p, 
	const long double nmax, const long double DPl, const long double alpha, const long double q, const long double sigma_p, 
	const long double fmin, const long double fmax, const long double resol, bool returns_pg_freqs, bool verbose)
{
	const bool returns_axis=true;
	const int Nmmax=100000; //Ngmax+Npmax;
	const double tol=2*resol; // Tolerance while searching for double solutions of mixed modes
	
	unsigned seed_p = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen_p(seed_p); //, gen_g(seed_g);
	std::normal_distribution<double> distrib_p(0.,sigma_p);
	
	bool success;
	int np, ng, np_min, np_max, ng_min, ng_max;
	double nu_p, nu_g, Dnu_p_local, DPl_local; // Dnu_p_local and DPl_local are important if modes do not follow exactly the asymptotic relation.
	double fact=0.04;  // Default factor
	double r;

	VectorXd nu_p_all, nu_g_all, nu_m_all(Nmmax);

	Data_eigensols nu_sols;
	Deriv_out deriv_p, deriv_g;

	// Use fmin and fmax to define the number of pure p modes and pure g modes to be considered
	np_min=int(floor(fmin/Dnu_p - epsilon - el/2 - delta0l));
	np_max=int(ceil(fmax/Dnu_p - epsilon - el/2 - delta0l));

	np_min=int(floor(np_min - alpha_p*std::pow(np_min - nmax, 2) /2.));
	np_max=int(ceil(np_max + alpha_p*std::pow(np_max - nmax, 2) /2.)); // CHECK THIS DUE TO - -

	ng_min=int(floor(1e6/(fmax*DPl) - alpha));
	ng_max=int(ceil(1e6/(fmin*DPl) - alpha));

	if (np_min <= 0)
	{
		np_min=1;
	}
	if (fmin <= 150) // overrides of the default factor in case fmin is low
	{
		fact=0.01;
	}
	if (fmin <= 50)
	{
		fact=0.005;
	}
	// Handling the p and g modes, randomized or not
	nu_p_all.resize(np_max-np_min);
	nu_g_all.resize(ng_max-ng_min);
	for (int np=np_min; np<np_max; np++)
	{
		nu_p=asympt_nu_p(Dnu_p, np, epsilon, el, delta0l, alpha_p, nmax);
		if (sigma_p == 0)
		{
			nu_p=asympt_nu_p(Dnu_p, np, epsilon, el, delta0l, alpha_p, nmax);
		} else{
			r = distrib_p(gen_p);
			nu_p=asympt_nu_p(Dnu_p, np, epsilon, el, delta0l, alpha_p, nmax, r);
		}		
		nu_p_all[np-np_min]=nu_p;
	}

	for (int ng=ng_min; ng<ng_max;ng++)
	{
		nu_g=asympt_nu_g(DPl, ng, alpha);
		nu_g_all[ng-ng_min]=nu_g;
	}
	deriv_p=Frstder_adaptive_reggrid(nu_p_all);
	deriv_g.deriv.resize(nu_g_all.size());
	deriv_g.deriv.setConstant(DPl);
	std::vector<double> filteredVec;
	#pragma omp parallel for collapse(2)
	for (size_t np = 0; np < nu_p_all.size(); np++) {
		for (size_t ng = 0; ng < nu_g_all.size(); ng++) {
			double nu_p = nu_p_all[np];
			double nu_g = nu_g_all[ng];
			double Dnu_p_local = Dnu_p * (1.0 + alpha_p * (np + np_min - nmax));
			double DPl_local = DPl;

			Data_coresolver sols_iter = solver_mm(nu_p, nu_g, Dnu_p_local, DPl_local, q, nu_p - 1.75 * Dnu_p, nu_p + 1.75 * Dnu_p, resol, returns_axis, verbose, fact);
			if (sols_iter.nu_m.size() > 0) {
				for (int i = 0; i < sols_iter.nu_m.size(); i++) {
					if (sols_iter.nu_m[i] >= fmin && sols_iter.nu_m[i] <= fmax) {
						#pragma omp critical
						{
							filteredVec.push_back(sols_iter.nu_m[i]);
						}
					}
				}
			}
		}
	}
	
	std::sort(filteredVec.begin(), filteredVec.end());
	filteredVec.erase(std::unique(filteredVec.begin(), filteredVec.end(), [tol](double a, double b) {
		return std::abs(a - b) <= tol;
	}), filteredVec.end());

	nu_m_all.resize(filteredVec.size());
	std::copy(filteredVec.begin(), filteredVec.end(), nu_m_all.data());
	nu_m_all.resize(filteredVec.size());
	//#pragma omp parallel for
	//for (int i = 0; i < filteredVec.size(); i++) {
	//	nu_m_all[i] = filteredVec[i];
	//}

	if (returns_pg_freqs == true)
	{
		nu_sols.nu_m=nu_m_all;
		nu_sols.nu_p=nu_p_all;
		nu_sols.nu_g=nu_g_all;
		nu_sols.dnup=deriv_p.deriv;
		nu_sols.dPg=deriv_g.deriv;
		return nu_sols;
	} else
	{
		nu_sols.nu_m=nu_m_all;
		return nu_sols;
	}
}



// This function uses solver_mm to find solutions from a spectrum
// of pure p modes and pure g modes following the asymptotic relations at the second order for p modes and the first order for g modes
//
//	nu_l0_in: Frequencies for the l=0 modes. Used to derive nu_l1_p and therefore Dnu and epsilon
//	el: Degree of the mode
//	delta0l: first order shift related to core structure (and to D0)
//	alpha_p: Second order shift relate to the mode curvature
//	nmax: radial order at numax
//	DPl: average Period spacing of the g modes
//	alpha: phase offset for the g modes
//	q: coupling strength
//	sigma_p: standard deviation controling the randomisation of individual p modes. Set it to 0 for no spread
//	sigma_g: standard deviation controling the randomisation of individial g modes. Set it to 0 for no spread
//  resol: Control the grid resolution. Might be set to the resolution of the spectrum
//  returns_pg_freqs: If true, returns the values for calculated p and g modes
//  verbose: If true, print the solution on screen 
Data_eigensols solve_mm_asymptotic_O2from_l0(const VectorXd& nu_l0_in, const int el, const long double delta0l, 
    const long double DPl, const long double alpha, const long double q, const long double sigma_p, 
	const long double resol, bool returns_pg_freqs, bool verbose, const long double freq_min, const long double freq_max)
{

	const bool returns_axis=true;
	const int Nmmax=100000; //Ngmax+Npmax;
	const double tol=2*resol; // Tolerance while searching for double solutions of mixed modes
	unsigned seed_p = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen_p(seed_p); //, gen_g(seed_g);
	std::normal_distribution<double> distrib_p(0.,sigma_p);
	
	int np, ng, ng_min, ng_max;//, np_min, np_max, attempts;
	double nu_p, nu_g, Dnu_p, epsilon, Dnu_p_local, DPl_local, fmin, fmax; // Dnu_p_local and DPl_local are important if modes does not follow exactly the asymptotic relation.
	double fact=0.04;  // Default factor

	VectorXi test;
	VectorXd fit, nu_p_all, nu_g_all, nu_m_all(Nmmax), results(Nmmax);	

	Data_coresolver sols_iter;
	Data_eigensols nu_sols;
	Deriv_out deriv_p, deriv_g;

	const Eigen::VectorXd tmp = Eigen::VectorXd::LinSpaced(nu_l0_in.size(), 0, nu_l0_in.size()-1);
    
	fit=linfit(tmp, nu_l0_in); // fit[0] is the slope ==> Dnu and fit[1] is the ordinate at origin ==> fit[1]/fit[0] = epsilon
	Dnu_p=fit[0];
	epsilon=fit[1]/fit[0];
	epsilon=epsilon - floor(epsilon);
	fmin=nu_l0_in.minCoeff() - Dnu_p;
	fmax=nu_l0_in.maxCoeff() + Dnu_p;

	nu_m_all.setConstant(-9999);
	if (fmin < 0){
		fmin=0;
	}

	ng_min=int(floor(1e6/(fmax*DPl) - alpha));
	ng_max=int(ceil(1e6/(fmin*DPl) - alpha));

	if (fmin <= 150) // overrides of the default factor in case fmin is low
	{
		fact=0.01;
	}
	if (fmin <= 50)
	{
		fact=0.005;
	}

	// Handling the p and g modes, randomized or not
	nu_g_all.resize(ng_max-ng_min);

	// Step of extrapolating edges to avoid egdes effect when shifting l=0 frequencies to generate l=1 p modes
	nu_p_all=asympt_nu_p_from_l0_Xd(nu_l0_in, Dnu_p, el, delta0l, fmin, fmax);
	
	for (int ng=ng_min; ng<ng_max;ng++)
	{
		nu_g=asympt_nu_g(DPl, ng, alpha);
		nu_g_all[ng-ng_min]=nu_g;
	}
	
	deriv_p=Frstder_adaptive_reggrid(nu_p_all);
	deriv_g.deriv.resize(nu_g_all.size());
	deriv_g.deriv.setConstant(DPl);
	std::vector<double> filteredVec;
	#pragma omp parallel for collapse(2)
	for (size_t np = 0; np < nu_p_all.size(); np++) {
		for (size_t ng = 0; ng < nu_g_all.size(); ng++) {
			double nu_p = nu_p_all[np];
			double nu_g = nu_g_all[ng];
			// This is the local Dnu_p which differs from the average Dnu_p because of the curvature. The solver needs basically d(nu_p)/dnp , which is Dnu if O2 terms are 0.
			Dnu_p_local=deriv_p.deriv[np]; 
			DPl_local=DPl; // The solver needs here d(nu_g)/dng. Here we assume no core glitches so that it is the same as DPl. 	
			
			Data_coresolver sols_iter = solver_mm(nu_p, nu_g, Dnu_p_local, DPl_local, q, nu_p - 1.75 * Dnu_p, nu_p + 1.75 * Dnu_p, resol, returns_axis, verbose, fact);
			if (sols_iter.nu_m.size() > 0) {
				for (int i = 0; i < sols_iter.nu_m.size(); i++) {
					if (sols_iter.nu_m[i] >= freq_min && sols_iter.nu_m[i] <= freq_max) {
						#pragma omp critical
						{
							filteredVec.push_back(sols_iter.nu_m[i]);
						}
					}
				}
			}
		}
	}
	std::sort(filteredVec.begin(), filteredVec.end());
	filteredVec.erase(std::unique(filteredVec.begin(), filteredVec.end(), [tol](double a, double b) {
		return std::abs(a - b) <= tol;
	}), filteredVec.end());

	nu_m_all.resize(filteredVec.size());
	std::copy(filteredVec.begin(), filteredVec.end(), nu_m_all.data());
	nu_m_all.resize(filteredVec.size());

	if (returns_pg_freqs == true)
	{
		nu_sols.nu_m=nu_m_all;
		nu_sols.nu_p=nu_p_all;
		nu_sols.nu_g=nu_g_all;
		nu_sols.dnup=deriv_p.deriv;
		nu_sols.dPg=deriv_g.deriv;
		return nu_sols;
	} else
	{
		nu_sols.nu_m=nu_m_all;
		return nu_sols;
	}
}


// This function uses solver_mm to find solutions from a spectrum
// of pure p modes and pure g modes following the asymptotic relations at the second order for p modes and the first order for g modes
//
// CONTRARY TO *froml0, this function takes nu(l) as in input directly. HOWEVER NOTE THAT 
// IN THIS FORM, IT USES linfit() TO COMPUTE Dnu from nu(l). THIS IS NOT OPTIMAL AS Dnu(l=0) 
// IS USUALLY SLIGHTLY DIFFERENT THAN Dnu(l)
// ANOTHER DIFFERENCE IS THAT THIS FUNCTION SEARCH SOLUTION OVER +/- 1.75*Dnu_p_local 
// WHILE THE OTHER SIMILAR FUNCTIONS (O2 of froml0) SEARCH OVER A CONSTANT WINDOW Dnu_p
// BOTH ARE ACCEPTABLE BUT WILL LEAD TO SLIGHTLY DIFFERENT SET OF SOLUTIONS
//	nu_p_all: Frequencies for the l modes.
//	el: Degree of the mode
//	alpha_p: Second order shift relate to the mode curvature
//	nmax: radial order at numax
//	DPl: average Period spacing of the g modes
//	alpha: phase offset for the g modes
//	q: coupling strength
//	sigma_p: standard deviation controling the randomisation of individual p modes. Set it to 0 for no spread
//	sigma_g: standard deviation controling the randomisation of individial g modes. Set it to 0 for no spread
//  resol: Control the grid resolution. Might be set to the resolution of the spectrum
//  returns_pg_freqs: If true, returns the values for calculated p and g modes
//  verbose: If true, print the solution on screen 
Data_eigensols solve_mm_asymptotic_O2from_nupl(const VectorXd& nu_p_all, const int el, //const long double delta0l, 
    const long double DPl, const long double alpha, const long double q, const long double sigma_p, 
	const long double resol, bool returns_pg_freqs, bool verbose, const long double freq_min, const long double freq_max)
{

	const bool returns_axis=true;
	const int Nmmax=100000; //Ngmax+Npmax;
	const double tol=2*resol; // Tolerance while searching for double solutions of mixed modes
	unsigned seed_p = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen_p(seed_p); //, gen_g(seed_g);
	std::normal_distribution<double> distrib_p(0.,sigma_p);

	int np, ng, ng_min, ng_max;//, np_min, np_max, attempts;
	double nu_p, nu_g, Dnu_p, epsilon, Dnu_p_local, DPl_local, fmin, fmax; // Dnu_p_local and DPl_local are important if modes does not follow exactly the asymptotic relation.
	double fact=0.04;  // Default factor

	VectorXi test;
	VectorXd fit, nu_g_all, nu_m_all(Nmmax);	

	Data_coresolver sols_iter;
	Data_eigensols nu_sols;
	Deriv_out deriv_p, deriv_g;

	const Eigen::VectorXd tmp = Eigen::VectorXd::LinSpaced(nu_p_all.size(), 0, nu_p_all.size()-1);
    
	// There is a problem here because nu_p_all IS NOT for l=0
	// This function need to be updated to include Dnu_p as input.
	fit=linfit(tmp, nu_p_all); // fit[0] is the slope ==> Dnu and fit[1] is the ordinate at origin ==> fit[1]/fit[0] = epsilon
	Dnu_p=fit[0];

	fmin=nu_p_all.minCoeff() - Dnu_p; // Range for setting the number of g modes 
	fmax=nu_p_all.maxCoeff() + Dnu_p;

	nu_m_all.setConstant(-9999);
	if (fmin < 0){
		fmin=0;
	}

	ng_min=int(floor(1e6/(fmax*DPl) - alpha));
	ng_max=int(ceil(1e6/(fmin*DPl) - alpha));

	if (fmin <= 150) // overrides of the default factor in case fmin is low
	{
		fact=0.01;
	}
	if (fmin <= 50)
	{
		fact=0.005;
	}

	// Handling the p and g modes, randomized or not
	nu_g_all.resize(ng_max-ng_min);

	for (int ng=ng_min; ng<ng_max;ng++)
	{
		nu_g=asympt_nu_g(DPl, ng, alpha);
		nu_g_all[ng-ng_min]=nu_g;
	}	
	deriv_p=Frstder_adaptive_reggrid(nu_p_all);
	deriv_g.deriv.resize(nu_g_all.size());
	deriv_g.deriv.setConstant(DPl);
	
	std::vector<double> filteredVec;
	#pragma omp parallel for collapse(2)
	for (size_t np = 0; np < nu_p_all.size(); np++) {
		for (size_t ng = 0; ng < nu_g_all.size(); ng++) {
			double nu_p = nu_p_all[np];
			double nu_g = nu_g_all[ng];
			// This is the local Dnu_p which differs from the average Dnu_p because of the curvature. The solver needs basically d(nu_p)/dnp , which is Dnu if O2 terms are 0.
			Dnu_p_local=deriv_p.deriv[np]; 
			DPl_local=DPl; // The solver needs here d(nu_g)/dng. Here we assume no core glitches so that it is the same as DPl. 	
			Data_coresolver sols_iter = solver_mm(nu_p, nu_g, Dnu_p_local, DPl_local, q, nu_p - 1.75 * Dnu_p_local, nu_p + 1.75 * Dnu_p_local, resol, returns_axis, verbose, fact);
			if (sols_iter.nu_m.size() > 0) {
				for (int i = 0; i < sols_iter.nu_m.size(); i++) {
					if (sols_iter.nu_m[i] >= freq_min && sols_iter.nu_m[i] <= freq_max) {
						#pragma omp critical
						{
							filteredVec.push_back(sols_iter.nu_m[i]);
						}
					}
				}
			}
		}
	}
	std::sort(filteredVec.begin(), filteredVec.end());
	filteredVec.erase(std::unique(filteredVec.begin(), filteredVec.end(), [tol](double a, double b) {
		return std::abs(a - b) <= tol;
	}), filteredVec.end());

	nu_m_all.resize(filteredVec.size());
	std::copy(filteredVec.begin(), filteredVec.end(), nu_m_all.data());
	nu_m_all.resize(filteredVec.size());
	
	
	if (returns_pg_freqs == true)
	{
		nu_sols.nu_m=nu_m_all;
		nu_sols.nu_p=nu_p_all;
		nu_sols.nu_g=nu_g_all;
		nu_sols.dnup=deriv_p.deriv;
		nu_sols.dPg=deriv_g.deriv;
		return nu_sols;
	} else
	{
		nu_sols.nu_m=nu_m_all;
		return nu_sols;
	}
}