// -------------------
// ---- Functions adapted from the bump_DP.py function ----
// This contains all the function that describes the Bumped period spacing
// and how they can be used to derived mode amplitudes from inertia ratio
// as they have been developed and tested. This arises from the following publications:
// https://arxiv.org/pdf/1509.06193.pdf (Inertia and ksi relation)
// https://www.aanda.org/articles/aa/pdf/2015/08/aa26449-15.pdf (Eq. 17 for the rotation - splitting relation)
// https://arxiv.org/pdf/1401.3096.pdf (Fig 13 and 14 for the evolution - rotation relation in SG and RGB) 
// https://arxiv.org/pdf/1505.06087.pdf (Eq. 3 for determining log(g) using Teff and numax, used for getting the evolution stage in conjonction with Fig. 13 from above)
// https://iopscience.iop.org/article/10.1088/2041-8205/781/2/L29/pdf (Inertia and Height relation)
// https://arxiv.org/pdf/1707.05989.pdf (Fig.5, distribution of surface rotation for 361 RGB stars)

// ------------------
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "version_solver.h"
#include "data_solver.h"
#include "string_handler.h"
//#include "interpol.h" // 
#include "../../interpol.h"
//#include "noise_models.h" // get the harvey_1985 function IF STANDALONE OR TAMCMC
#include "../../noise_models.h" // IF inside the spectrum simulator
#include "solver_mm.h"
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

/*
# the ksi function as defined in equation 14 of Mosser+2017 (https://arxiv.org/pdf/1509.06193.pdf)
# Inputs: 
# 	nu: The freqency axis (microHz)
# 	nu_p : Frequenc of a p-mode (microHz)
#	nu_g : Frequency of a g-mode (microHz)
#	Dnu_p: Large separation. DOES NOT NEED TO BE A CONSTANT (fonction of np or nu to account for glitches)
#   DPl: Period Spacing (seconds)
#	q : Coupling term (no unit)
*/
VectorXd ksi_fct1(const VectorXd& nu, const long double nu_p, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q)
{
	

	const long double pi = 3.141592653589793238L;
	VectorXd cos_upterm, cos_downterm, front_term, tmp(nu.size()), tmp2(nu.size()), ksi(nu.size());

	tmp.setConstant(1./nu_g);
	cos_upterm=pi * 1e6 * (nu.cwiseInverse() - tmp)/DPl;

	tmp.setConstant(nu_p);
	cos_downterm=pi * (nu - tmp) /Dnu_p;
	front_term= 1e-6 * nu.array().square() * DPl / (q * Dnu_p); // relation accounting for units in Hz and in seconds

	//ksi=1./(1. + front_term * cos_upterm.array().cos().square()/cos_downterm.array().cos().square());
	tmp2=cos_upterm.array().cos().square()/cos_downterm.array().cos().square();

	tmp.setConstant(1);
	ksi=front_term.cwiseProduct(tmp2);
	ksi=(tmp + ksi).cwiseInverse();

/*
	std::cout << "---- Debug ksi_fct1 ----" << std::endl;
	std::cout << "nu_p = " << nu_p << std::endl;
	std::cout << "nu_g = " << nu_g << std::endl;
	std::cout << "cos_upterm =" << cos_upterm << std::endl;
	std::cout << "cos_downterm =" << cos_downterm << std::endl;
	std::cout << "front_term =" << front_term << std::endl;
	std::cout << "front * cos^2/sin^2 =" << front_term.cwiseProduct(tmp2) << std::endl;
	std::cout << "ksi = " << ksi << std::endl;
	std::cout << "-----" << std::endl;
	std::cout << "-----" << std::endl;
	exit(EXIT_SUCCESS);
*/
	return ksi;
}

/*
# Variant of ksi_fct that deals with arrays for nu_p, nu_g, Dnu_p, DPl
# This requires that nu_p and Dnu_p have the same dimension
# Also nu_g and DPl must have the same dimension
# Additional parameter:
#  - norm-method: When set to "fast", normalise by the max of the ksi_pg calculated at the
#				   Frequencies nu given by the user
#				   When set to 'exact', normalise by the max of a heavily interpolated function
#				   of ksi_pg. Allows a much a higher precision, but will be slower
#				   This could be usefull as in case of low ng, the norm is badly estimated in
#				   "fast" mode. Then we need to use a more continuous function to evaluate the norm
*/
VectorXd ksi_fct2(const VectorXd& nu, const VectorXd& nu_p, const VectorXd& nu_g, const VectorXd& Dnu_p, const VectorXd& DPl, const long double q, const std::string norm_method="fast")
{
	const int Lp=nu_p.size();
	const int Lg=nu_g.size();
	const long double resol=1e6/(4*365.*86400.); // Fix the grid resolution to 4 years (converted into microHz)

	VectorXd ksi_tmp, ksi_pg(nu.size()), nu4norm, ksi4norm;

	int Ndata;
	long double norm_coef, fmin,fmax;

	ksi_pg.setConstant(nu.size());
	ksi_pg.setZero();
//#pragma omp parallel for default(shared) private(ksi_tmp)
	for (int np=0; np<Lp;np++)
	{
		for (int ng=0; ng<Lg; ng++)
		{
			ksi_tmp=ksi_fct1(nu, nu_p[np], nu_g[ng], Dnu_p[np], DPl[ng], q);
			ksi_pg=ksi_pg + ksi_tmp;
		}
	}
	if (norm_method == "fast"){
		norm_coef=ksi_pg.maxCoeff();
	}
	else
	{ // We build a very resolved 'continuous function of the frequency to calculate the norm'
		if (nu_p.minCoeff() >= nu_g.minCoeff()){
			fmin=nu_g.minCoeff();
		} else{
			fmin=nu_p.minCoeff();
		}
		if (nu_p.maxCoeff() >= nu_g.maxCoeff()){
			fmax=nu_p.maxCoeff();
		} else{
			fmax=nu_g.maxCoeff();
		}
		Ndata=int((fmax-fmin)/resol);
		nu4norm=linspace(fmin, fmax, Ndata);
		ksi4norm.resize(nu4norm.size());
		ksi4norm.setZero();
//#pragma omp parallel for default(shared) private(ksi_tmp)
		for (int np=0; np<Lp; np++){
			for (int ng=0; ng<Lg; ng++){
				ksi_tmp=ksi_fct1(nu4norm, nu_p[np], nu_g[ng], Dnu_p[np], DPl[ng], q);
				ksi4norm=ksi4norm + ksi_tmp;
			}	
		}
		norm_coef=ksi4norm.maxCoeff();
	}
	ksi_pg=ksi_pg/norm_coef;
	// Ensuring that round off error don't lead to values higher than 1...
	for(int i=0; i<ksi_pg.size(); i++){
		if (ksi_pg[i]>1){
			ksi_pg[i]=1;
		}
	}
	return ksi_pg;
}

VectorXd gamma_l_fct2(const VectorXd& ksi_pg, const VectorXd& nu_m, const VectorXd& nu_p_l0, const VectorXd& width_l0, const VectorXd& hl_h0_ratio, const int el, const long double factor=1.0)
{
	long double width0_at_l;
	VectorXd width_l(ksi_pg.size());

	if (  (nu_p_l0.size() != width_l0.size()) || (ksi_pg.size() != nu_m.size()) )
	{
		std::cout << "Inconsistency between the size of the Width and l=0 frequency array or between ksi_pg and nu_m arrays" << std::endl;
		std::cout << "Cannot pursue. The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	} else
	{
		// Perform the interpolation
		for (int i=0; i<ksi_pg.size(); i++)
		{
			width0_at_l=lin_interpol(nu_p_l0, width_l0, nu_m[i]);
			width_l[i]=width0_at_l * (1. - factor*ksi_pg[i])/ std::sqrt(hl_h0_ratio[i]);
		}
		// ---- DEBUG LINES ----
		//std::cout << "   DEBUG FOR gamma_l_fct2..." << std::endl;
		//std::cout << "ksi_pg     ,   width_l    ,   hl_h0_ratio" << std::endl;
		//for (int i=0; i<ksi_pg.size(); i++)
		//{
		//	std::cout << ksi_pg[i] << "    " << width_l[i] << "    "  << hl_h0_ratio[i] << std::endl;
		//}
	}
	return width_l;
}

VectorXd h_l_rgb(const VectorXd& ksi_pg, const long double factor=1.0)
{
	const double tol=1e-5;
	VectorXi pos;
	VectorXd tmp(ksi_pg.size()), hl_h0;

	tmp.setConstant(1);
	hl_h0=tmp - factor*ksi_pg;
	hl_h0=hl_h0.array().sqrt();
	pos=where_dbl(hl_h0, 0, tol);
	if(pos[0] != -1)
	{
		for (int i=0;i<pos.size();i++)
		{
			hl_h0[pos[i]] = 1e-10;
		}
	}
	// --- DEBUG LINES ---
	//std::cout << "   DEBUG FOR h_l_rgb..." << std::endl;
	//std::cout << "  ksi_pg       hl_l0 " << std::endl;
	//for (int i=0; i<ksi_pg.size(); i++)
	//{
	//	std::cout << ksi_pg[i] << "    " << hl_h0[i] << std::endl;
	//}
	return hl_h0;
}

// Put here the code for reading template files that contain heights and width profiles
template_file read_templatefile(const std::string file, const bool ignore_errors=true){

	const std::string delimiter="=";
	const int Ncols=3;
	const int Nrows=200;
	
	size_t pos;
	int cpt;	
	//std::vector<std::string> tmp;

	std::string key;

	std::string line0, char0, char1;
	std::vector<std::string> word, tmp;
	std::ifstream tmpfile_session;

	template_file template_data;
	
	tmpfile_session.open(file.c_str());
    if (tmpfile_session.is_open() == false) {
 		std::cout << "Unable to open the file: " << file << std::endl;
   		std::cout << "Check that the file exist and that the path is correct" << std::endl;
   		std::cout << "Cannot proceed" << std::endl;
   		std::cout << "The program will exit now" << std::endl;
   		exit(EXIT_FAILURE);
    } 

	char0="#"; 
	while(char0 == "#" && (!tmpfile_session.eof())) // Ignore comments in the header
	{
		std::getline(tmpfile_session, line0);
		line0=strtrim(line0);
		char0=strtrim(line0.substr(0, 1));
	}

	// Initialise to out of range values
	template_data.ID_ref=" ";
	template_data.numax_ref=-1;
	template_data.Dnu_ref=-1;
	template_data.epsilon_ref=-1;

	while (char0 != "#" && (!tmpfile_session.eof())) // Look for keywords for numax, Dnu and ID
	{
		pos=line0.find(delimiter);
		key=line0.substr(0, pos); // What is before the delimiter is used as a key
		if (key == "ID_ref")
		{
			template_data.ID_ref=str_to_dbl(line0.substr(pos + delimiter.length(), std::string::npos));
		}
		if (key == "numax_ref")
		{	
			template_data.numax_ref=str_to_dbl(line0.substr(pos + delimiter.length(), std::string::npos));
		}
		if (key == "Dnu_ref")
		{
			template_data.Dnu_ref=str_to_dbl(line0.substr(pos + delimiter.length(), std::string::npos));
		}
		if (key == "epsilon_ref")
		{
			template_data.epsilon_ref=str_to_dbl(line0.substr(pos + delimiter.length(), std::string::npos));
		}
		std::getline(tmpfile_session, line0);
		line0=strtrim(line0);
		char0=strtrim(line0.substr(0, 1));
	}
	if (ignore_errors == false)
	{
		if (template_data.numax_ref == -1 || template_data.Dnu_ref == -1 || template_data.epsilon_ref == -1)
		{
			std::cout <<"Error: Could not find at least one of the keywords defining the global pulsation parameters" << std::endl;
			std::cout <<"Check that the following keywords are present in the template file: " << file  << std::endl;
			std::cout <<"       The program will exit now"  << std::endl;
			exit(EXIT_FAILURE);
		}
		if (template_data.ID_ref == " ")
		{
			std::cout <<"Warning: ID_ref is not set"  << std::endl;
			std::cout <<"         This does not prevent the code to run, but may result in a more difficult tracking of the used reference mode profiles in the future" << std::endl;
		}
	}
	// Process the table that should be marked by a "#" first
	template_data.data_ref.resize(Nrows, Ncols);
	cpt=0;
	while (!tmpfile_session.eof())
	{ 
		std::getline(tmpfile_session, line0); // read each line of the table
		template_data.data_ref.row(cpt)=str_to_Xdarr(line0, " \t");
		cpt=cpt+1;
	}
	template_data.data_ref.conservativeResize(cpt,Ncols);
	return template_data;
}

Data_2vectXd width_height_load_rescale(const VectorXd& nu_star, const long double Dnu_star, const long double numax_star, const std::string file)
{
	int Nref, Nstar;
	long double n_at_numax_ref, height_ref_at_numax, gamma_ref_at_numax, epsilon_star, n_at_numax_star, w_tmp, h_tmp;
	VectorXd tmp, tmp_ref, nu_ref, height_ref, gamma_ref, en_list_ref, en_list_star, w_star, h_star;
	Data_2vectXd out;

	template_file template_data;

	template_data= read_templatefile(file);

	nu_ref=template_data.data_ref.col(0); //[:,0] // CHECKS NEEDED REALLY COL OR IS IT ROW?
	height_ref=template_data.data_ref.col(1);//[:,1]
	gamma_ref=template_data.data_ref.col(2);//[:,2]

	height_ref_at_numax=lin_interpol(nu_ref, height_ref, template_data.numax_ref);
	gamma_ref_at_numax=lin_interpol(nu_ref, gamma_ref, template_data.numax_ref);
	n_at_numax_ref=template_data.numax_ref/template_data.Dnu_ref - template_data.epsilon_ref;

	Nref=nu_ref.size();
	tmp.resize(Nref);
	en_list_ref.resize(Nref);
	tmp.setConstant(template_data.epsilon_ref);
	en_list_ref=nu_ref/template_data.Dnu_ref - tmp; // This list will be monotonic
	// ------------------------------------------------------------------------------------
	// Rescaling using the base frequencies given above for the Sun
	epsilon_star=0;
	for (int i=0; i< nu_star.size(); i++){
		epsilon_star=epsilon_star + ( std::fmod(nu_star[i]/Dnu_star, 1));
	}
	epsilon_star=epsilon_star/nu_star.size();

	n_at_numax_star=numax_star/Dnu_star - epsilon_star;
	
	Nstar=nu_star.size();
	tmp.resize(Nstar);
	en_list_star.resize(Nstar);

	tmp.setConstant(epsilon_star);
	en_list_star=nu_star/Dnu_star - tmp;

	tmp_ref.resize(en_list_ref.size());
	tmp_ref.setConstant(n_at_numax_ref);
			
	tmp_ref=en_list_ref - tmp_ref;

	w_star.resize(en_list_star.size());
	h_star.resize(en_list_star.size());

	for (int en=0; en<en_list_star.size(); en++){
		w_tmp=lin_interpol(tmp_ref, gamma_ref/gamma_ref_at_numax, en_list_star[en] - n_at_numax_star);
		h_tmp=lin_interpol(tmp_ref, height_ref/height_ref_at_numax, en_list_star[en] - n_at_numax_star); 
		w_star[en]=w_tmp;
		h_star[en]=h_tmp;
	}
	h_star=h_star/h_star.maxCoeff(); // Normalise so that that HNR(numax)=1 if white noise N0=1

	out.vecXd1=w_star;
	out.vecXd2=h_star;
	return out;
}


// Simple way of computing the core rotation from the surface rotation. Used if we want uniform 
// distribution of rotation in the envelope and a uniform population of core-to-envelope ratios 
// 	 (1) rot_envelope: average rotation in the envelope
//	 (2) core2envelope_star: average rotation in the core 
Data_rot2zone rot_2zones_v2(const long double rot_envelope, const long double core2envelope_star, std::string output_file_rot=" ")
{

	Data_rot2zone rot2data;

	// Determine the core rotation
	const long double rot_core=core2envelope_star * rot_envelope;
	//std::ostringstream strg;
	std::ofstream outfile;
	if (output_file_rot != " ")
	{
		outfile.open(output_file_rot.c_str());
		if (outfile.is_open()){
			//outfile << strg.str().c_str();
			outfile << "#Average envelope rotation (microHz) /  Average core rotation  (microHz)\n";
			outfile << rot_envelope <<  "  " << rot_core;
			outfile.close();
		}
	}

	rot2data.rot_env=rot_envelope;
	rot2data.rot_core=rot_core;

	return rot2data;
}

// Simple way of computing the core rotation from the surface rotation. Used if we want uniform 
// distribution of rotation in the envelope and a uniform population of core rotation
// 	 (1) rot_envelope: average rotation in the envelope
//	 (2) rot_core: average rotation in the core 
Data_rot2zone rot_2zones_v3(const long double rot_envelope, const long double rot_core, std::string output_file_rot=" ")
{
	Data_rot2zone rot2data;
	std::ofstream outfile;

	if (output_file_rot != " ")
	{
		outfile.open(output_file_rot.c_str());
		if (outfile.is_open()){
			outfile << "#Average envelope rotation (microHz) /  Average core rotation  (microHz)\n";
			outfile << rot_envelope <<  "  " << rot_core;
			outfile.close();
		}
	}
	rot2data.rot_env=rot_envelope;
	rot2data.rot_core=rot_core;

	return rot2data;
}

// Function that determine the rotation in the envelope, here approximated to be the surface rotation.
// Inspired by the surface rotation from Ceillier et al. 2017 (https://arxiv.org/pdf/1707.05989.pdf), Fig. 5
// They have a skewed distribution going for ~30 days to ~160 days with a peak around 60 days. 
// For simplicity, I just insert a truncated gaussian distribution with rotation between 30  and 90 and a median of 60.
// The truncation happens at sigma. Values are given in days
// Returns: 
//	rot_s: rotation frequency in microHz
long double rot_envelope(long double med=60., long double sigma=3.)
{
	const long double var=30./sigma; // in days
	
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen(seed); 
	std::normal_distribution<double> distrib(med,var);
	
	long double period_s=distrib(gen);
	long double rot_s;

	if (period_s < med - var*sigma)
	{
		period_s=med-var*sigma;
	}
	if (period_s > med + var*sigma)
	{
		period_s=med+var*sigma;
	}
	rot_s=1e6/(86400.*period_s);

	return rot_s;
}

/*
# Splitting of modes assuming a two-zone (or two-zone averaged) rotation profile
# rot_envelope and rot_core must be in Hz units (or micro, nano, etc...)
# ksi_pg: The ksi function that describes the degree of mixture between p and g modes in function of the more frequency
# rot_envelope: average rotation in the envelope. Must be a scalar
# rot_core: average rotation in the core. Must be a scalar
# Returns:
#	dnu_rot: A vector of same size as ksi_pg
*/
VectorXd dnu_rot_2zones(const VectorXd& ksi_pg, const long double rot_envelope, const long double rot_core)
{
	//VectorXd rc(ksi_pg.size()), re(ksi_pg.size());
	VectorXd re(ksi_pg.size());

	//rc.setConstant(rot_core/2);
	re.setConstant(rot_envelope);

	return ksi_pg*(rot_core/2 - rot_envelope) + re;
}	

//Assumptions: nu_max is derived from the requested Dnu_star parameter using the relation from Stello+2009. 
//	Dnu ~ 0.263*numax^0.77 (no uncertainty implemented here)
long double numax_from_stello2009(const long double Dnu_star, const long double spread)
{
	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	
	// Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	const long double beta0=0.263; // according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	const long double beta1=0.77; // according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	
	long double  nu_max=std::pow(10, std::log10(Dnu_star/beta0)/beta1);
//	std::cout <<"nu_max=", nu_max)	
	if (std::abs(spread)>0) // Add a unifrom spread around numax, spread must be given in fraction (eg 5% is 0.05)
	{	
		std::uniform_real_distribution<double> distrib(nu_max*(1.-std::abs(spread)),nu_max*(1. + std::abs(spread)));
		nu_max=distrib(gen);
	} else
	{
		std::cout << "Spread argument in bump_DP.cpp::numax_from_stello2009() is 0 or negative... Ignoring it" << std::endl;
	}
//	std::cout <<"spread: ", spread)
//	std::cout <<"New nu_max=", nu_max)
	return nu_max;
}


/*
# The main function that generate a set of parameters used to generate Lorentzian profiles FOR SIMULATIONS IN THE C++ SIMULATOR
# Assumptions: 
#     - The frequencies of l=0, l=2 and l=3 modes follow exactly the asymtptotic relation for the p mdoes
#	  - The frequencies of the l=1 modes follow exactly the asymptotitc relation for the mixed modes
#	  - Widths of l=0, 2, 3 modes are rescaled using the synthetic relation from Appourchaux+2014 applied to the solar profile
#	  - Widths of l=1 mixed modes are defined using the ksi function, scaled using l=0 modes. Thus l=0 modes fixes the upper
#		limit of the mode width
#	  - Heights of l=0 modes are rescaled using the measured heights of the solar modes
#	  - Heights of l=1 mixed modes are set based on equipartition of energy accounting for damping and excitation assuming no radiative pressure
#		GrosJean+201X showed that this was valid for not-too-evolved RGB and for SG.
#	  - Bolometric visibilities in height are assumed to be 1, 1.5, 0.5, 0.07 for l=0,1,2,3 respectively
#	  - Splitting is implemented in different ways... see the various models
# 	  Warning for widths and heights: if fmin/fmax is too large, you may have an error as an extrapolation 
#									  would be required, which I forbid by code. An fmin/fmax that englobes
#									  10 radial orders should not pose any issue.
# Input: 
#	Dnu_star: Large separation for the p modes
#   epsilon_star: phase offset term for the asymtptotic relation of the p modes
#   delta0l_star, alpha_p_star, nmax_star: Instead of D0_star, these parameters can be used to create 2nd order asymptotic p modes (see Mosser+2018, Eq.22) 
#   DPl_star: The period spacing for l=1 g modes 
#	alpha_star: The phase offset term for the asymptotic relation of the g modes
#   q_star: Coupling coeficient between p and g modes
#	fmin: Minimum frequency for the modes that should be included in the calculations
#   fmax: Maximum frequency for the modes that should be included in the calculations
# Outputs:
#	nu_lx: Frequencies of the l=x modes. x is between 0 and 3
#	nu_p_l1: Base p modes frequencies used to build the frequencies for the l=1 mixed modes
#	nu_g_l1: Base p modes frequencies used to build the frequencies for the l=1 mixed modes
#	width_lx: Widths of the l=x modes. x is between 0 and 3
#   height_lx: Heights of the l=x modes. x is between 0 and 3 
*/

Params_synthetic_star make_synthetic_asymptotic_star(Cfg_synthetic_star cfg_star)
{
	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	//std::uniform_real_distribution<double> distrib(xmin,xmax);
	std::uniform_real_distribution<double> distrib(0 , 1);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen_m(seed); 
	std::normal_distribution<double> distrib_m(0.,cfg_star.sigma_m);

	//std::cout << "cfg_star.use_nu_nl = " << cfg_star.use_nu_nl << std::endl;
	//exit(EXIT_SUCCESS);
	int en, el, np_min, np_max;
	long double r, tmp, resol, c, xmin ,xmax, delta0l_star;
	VectorXi posOK;
	VectorXd tmpXd, noise_params_harvey1985(4), noise_l0, hmax_l0, Dnu_p, DPl, ksi_pg, h1_h0_ratio;
		//nu_p_l1, nu_g_l1, 
	VectorXd nu_l0, nu_m_l1, nu_l2, nu_l3, 
		height_l0, height_l1, height_l2, height_l3, height_l1p, 
		width_l0, width_l1, width_l2, width_l3,
		a1_l1, a1_l2, a1_l3; // Simulate a single harvey profile
	int ng, ng_min, ng_max;  // used only if cfg_star.use_nu_nl = true
	VectorXd nu_p_l1, nu_g_l1; // used only if cfg_star.use_nu_nl = true

	Data_2vectXd width_height_l0;
	Data_rot2zone rot2data;
	Params_synthetic_star params_out;
	Data_eigensols freqs;

	//std::cout << cfg_star.numax_star << std::endl;
	//exit(EXIT_SUCCESS);
	//Defining what should be Hmax_l0 in order to get the desired HNR
	//                   04           1         2           3            4           5          6       7
	//noise_params_harvey_like=[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0] // 
	noise_params_harvey1985[0] = cfg_star.noise_params_harvey_like[0] * std::pow(cfg_star.numax_star*1e-6,cfg_star.noise_params_harvey_like[1]) + cfg_star.noise_params_harvey_like[2]; // Granulation Amplitude
	noise_params_harvey1985[1] = cfg_star.noise_params_harvey_like[3] * std::pow(cfg_star.numax_star*1e-6, cfg_star.noise_params_harvey_like[4]) + cfg_star.noise_params_harvey_like[5]; // Granulation timescale (in seconds)
	noise_params_harvey1985[0] = noise_params_harvey1985[0]/noise_params_harvey1985[1];
	noise_params_harvey1985[1]= noise_params_harvey1985[1]/1000.;

	noise_params_harvey1985[2]=cfg_star.noise_params_harvey_like[6];
	noise_params_harvey1985[3]=cfg_star.noise_params_harvey_like[7];

	// Fix the resolution to 4 years (converted into microHz)
	resol=1e6/(4*365.*86400.);

	if (cfg_star.use_nu_nl == false){ // Do all the stuff related to computing frequencies only if requested
		// ----- l=0 modes -----
		// This section generate l=0 modes following the asymptotic relation of p modes, and make
		// rescaled width and height profiles for the star using the solar width and height profiles

		// Use fmin and fmax to define the number of pure p modes and pure g modes to be considered
		np_min=int(floor(cfg_star.fmin/cfg_star.Dnu_star - cfg_star.epsilon_star));
		np_max=int(ceil(cfg_star.fmax/cfg_star.Dnu_star - cfg_star.epsilon_star));
		np_min=int(floor(np_min - cfg_star.alpha_p_star*std::pow(np_min - cfg_star.nmax_star, 2) /2.));
		np_max=int(ceil(np_max + cfg_star.alpha_p_star*std::pow(np_max - cfg_star.nmax_star, 2) /2.));  // The minus plus is there because (np_max - nmax_star)^2 is always positive
		
		if (np_min < 1)
		{
			np_min=1;
		}
	//	std::cout << "np_min =" << np_min << std::endl;
	//	std::cout << "np_max =" << np_max << std::endl;
		nu_l0.resize(np_max-np_min);
		for (en=np_min; en<np_max; en++)
		{
			tmp=asympt_nu_p(cfg_star.Dnu_star, en, cfg_star.epsilon_star, 0, 0, cfg_star.alpha_p_star, cfg_star.nmax_star);
			nu_l0[en-np_min]=tmp;
		}
	} else{
		nu_l0=cfg_star.nu_nl.row(0).segment(0, cfg_star.Nf_el[0]);
	}
	width_height_l0=width_height_load_rescale(nu_l0, cfg_star.Dnu_star, cfg_star.numax_star, cfg_star.filetemplate); // Function that ensure that Hmax and Wplateau is at numax
	width_l0=width_height_l0.vecXd1;
	height_l0=width_height_l0.vecXd2;
	noise_l0.resize(nu_l0.size());
	noise_l0.setZero();
	noise_l0=harvey1985(noise_params_harvey1985, nu_l0, noise_l0, 1); // Iterate on Noise_l0 to update it by putting the noise profile with one harvey profile
		
	c=1; // This is the ratio of HNR between the reference star and the target simulated star: maxHNR_l0/maxHNR_ref.
	hmax_l0=cfg_star.maxHNR_l0*noise_l0*c;
	height_l0=height_l0.cwiseProduct(hmax_l0); // height_l0 being normalised to 1 on width_height_load_rescale, getting the desired hmax_l0 requires just to multiply height_l0 by hmax_l0

	if (std::abs(cfg_star.H0_spread) > 0)
	{
		try
		{
			for (int i=0; i<height_l0.size();i++)
			{
				xmin=height_l0[i]*(1. - std::abs(cfg_star.H0_spread)/100.);
				xmax=height_l0[i]*(1. + std::abs(cfg_star.H0_spread)/100.);
				height_l0[i]=xmin + (xmax-xmin)*distrib(gen);
			}
		}
		catch(...)
		{
			std::cout << "Error debug info:" << std::endl;
			std::cout << "nu_l0 = " << nu_l0 << std::endl;
			std::cout << "hmax_l0 = " << hmax_l0 << std::endl;
			std::cout << "Height_l0: " << height_l0 << std::endl;
			std::cout << "cfg_star.H0_spread: " << cfg_star.H0_spread << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	width_l0=width_l0*cfg_star.Gamma_max_l0;
	// ------- l=1 modes ------
		// Use the solver to get mixed modes
	if (cfg_star.use_nu_nl == false){ // Do all the stuff related to computing frequencies only if requested
		el=1;
		delta0l_star=-el*(el + 1) * cfg_star.delta0l_percent_star / 100.;
		freqs=solve_mm_asymptotic_O2p(cfg_star.Dnu_star, cfg_star.epsilon_star, el, delta0l_star, cfg_star.alpha_p_star, cfg_star.nmax_star, cfg_star.DPl_star, 
									cfg_star.alpha_g_star, cfg_star.q_star, cfg_star.sigma_p, cfg_star.fmin, cfg_star.fmax, resol, true, false);

		// Filter solutions that endup at frequencies higher/lower than the nu_l0 because we will need to extrapolate height/widths otherwise...
		posOK=where_in_range(freqs.nu_m, nu_l0.minCoeff(), nu_l0.maxCoeff(), false);
		nu_m_l1.resize(posOK.size());
		for (int i=0; i<posOK.size();i++)
		{
			nu_m_l1[i]=freqs.nu_m[posOK[i]];
			if (cfg_star.sigma_m !=0) // If requested, we add a random gaussian qty to the mixed mode solution
			{
				r = distrib_m(gen_m);
				nu_m_l1[i]=nu_m_l1[i]+r;
			}
		}
		
		// Generating widths profiles for l=1 modes using the ksi function
		Dnu_p=freqs.dnup;
		DPl=freqs.dPg; 
		ksi_pg=ksi_fct2(nu_m_l1, freqs.nu_p, freqs.nu_g, Dnu_p, DPl, cfg_star.q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
	} else{
		nu_m_l1=cfg_star.nu_nl.row(1).segment(0, cfg_star.Nf_el[1]);
		//std::cout << "nu_m_l1=" << nu_m_l1.transpose() << std::endl;
		if (cfg_star.Dnu_star == -1 || cfg_star.DPl_star == -1 || cfg_star.alpha_g_star == -1 || cfg_star.q_star == -1){
			std::cout << "Error: cfg_star.use_nu_l is set to true " << std::endl;
			std::cout << "       You must provide in the cfg_star structure: " << std::endl;
			std::cout << "             - nu_nl in a matrix form" << std::endl;
			std::cout << "             - Dnu_star" << std::endl;
			std::cout << "             - DPl_star" << std::endl;
			std::cout << "             - alpha_g_star" << std::endl;
			std::cout << "             - q_star" << std::endl;
			exit(EXIT_FAILURE);
		}
		// get epsilon_p
		//std::cout << "epsilon " << std::endl;
		cfg_star.epsilon_star=0;
		for (int c=0; c<cfg_star.nu_nl.rows(); c++){
			cfg_star.epsilon_star=cfg_star.epsilon_star + cfg_star.nu_nl(0,c)/cfg_star.Dnu_star  - floor(cfg_star.nu_nl(0,c)/cfg_star.Dnu_star);
		}
		cfg_star.epsilon_star=cfg_star.epsilon_star/cfg_star.nu_nl.rows();
		std::cout << "     epsilon = " << cfg_star.epsilon_star << std::endl;
		el=1;
		delta0l_star=-el*(el + 1) * cfg_star.delta0l_percent_star / 100.;
		np_min=int(floor(cfg_star.fmin/cfg_star.Dnu_star - cfg_star.epsilon_star - el/2 - delta0l_star));
		np_max=int(ceil(cfg_star.fmax/cfg_star.Dnu_star - cfg_star.epsilon_star - el/2 - delta0l_star));
		ng_min=int(floor(1e6/(cfg_star.fmax*cfg_star.DPl_star) - cfg_star.alpha_g_star));
		ng_max=int(ceil(1e6/(cfg_star.fmin*cfg_star.DPl_star) - cfg_star.alpha_g_star));
		nu_p_l1.resize(np_max-np_min);
		//std::cout << "delta0l_star = " << delta0l_star << std::endl;
		for (en=np_min; en<np_max;en++){
			//nu_p=asympt_nu_p(Dnu_p, en, cfg_star.epsilon_star, 1, 0, 0, 0);	
			nu_p_l1[en-np_min]=asympt_nu_p(cfg_star.Dnu_star, en, cfg_star.epsilon_star, 1, delta0l_star, 0, 0);
		}
		//std::cout << "nu_p_l1 = " << nu_p_l1.transpose() << std::endl;
		nu_g_l1.resize(ng_max-ng_min);
		for (ng=ng_min; ng<ng_max;ng++){
			nu_g_l1[ng-ng_min]=asympt_nu_g(cfg_star.DPl_star, ng, cfg_star.alpha_g_star);
		}
		//std::cout << "nu_g_l1 = " << nu_g_l1.transpose() << std::endl;
		//std::cout << " ---" << std::endl;
		Dnu_p.resize(nu_p_l1.size());
		DPl.resize(nu_g_l1.size());
		Dnu_p.setConstant(cfg_star.Dnu_star);
		DPl.setConstant(cfg_star.DPl_star);
		ksi_pg=ksi_fct2(nu_m_l1, nu_p_l1, nu_g_l1, Dnu_p, DPl, cfg_star.q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
	}
	h1_h0_ratio=h_l_rgb(ksi_pg, cfg_star.Hfactor); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019). Hfactor added on May 2, 2022
	
	height_l1p.resize(nu_m_l1.size());
	for (int i=0; i<nu_m_l1.size();i++)
	{
		tmp=lin_interpol(nu_l0, height_l0, nu_m_l1[i]);
		height_l1p[i]=tmp;
	}

	height_l1p=height_l1p*cfg_star.Vl[0];
	height_l1=h1_h0_ratio.cwiseProduct(height_l1p);
	width_l1=gamma_l_fct2(ksi_pg, nu_m_l1, nu_l0, width_l0, h1_h0_ratio, el, cfg_star.Wfactor); // Wfactor added on May 2, 2022
	
	// Generating splittings with a two-zone averaged rotation rates
	if (cfg_star.rot_env_input >=0)
	{
		cfg_star.Teff_star=-1;
		if (cfg_star.rot_ratio_input > 0)
		{
			rot2data=rot_2zones_v2(cfg_star.rot_env_input, cfg_star.rot_ratio_input, cfg_star.output_file_rot); //rot_env, rot_c
		}
		else
		{
			rot2data=rot_2zones_v3(cfg_star.rot_env_input, cfg_star.rot_core_input, cfg_star.output_file_rot); //rot_env, rot_c
		}
	} else{
		// Determine the envelope rotation
		cfg_star.rot_env_input=rot_envelope(); // Warning: This function has optional arguments
		if (cfg_star.Teff_star >=0)
		{
			std::cout << " The option of determining the rotation rate through Teff_star is not anymore supported in this C++ only code version." << std::endl;
			std::cout << " For that purpose, please use the older version that integrates calls to python to handle mixed modes" << std::endl;
			std::cout << " The program will be terminated" << std::endl;
			exit(EXIT_SUCCESS);
			//rot_env, rot_c, rot_env_true, rot_c_true=rot_2zones_v1(rot_env_input, cfg_star.numax_star, Teff_star, sigma_logg=0.1, randomize_core_rot=True, output_file_rot=output_file_rot)
		}
	}

	a1_l1=dnu_rot_2zones(ksi_pg, rot2data.rot_env, rot2data.rot_core);
	// ------- l=2 modes -----
	if (cfg_star.use_nu_nl == false){ // Do all the stuff related to computing frequencies only if requested
		el=2;
		delta0l_star=-el*(el + 1) * cfg_star.delta0l_percent_star / 100.;
		nu_l2.resize(np_max-np_min);
		for (int en=np_min; en< np_max;en++)
		{
			tmp=asympt_nu_p(cfg_star.Dnu_star, en, cfg_star.epsilon_star, el, delta0l_star, cfg_star.alpha_p_star, cfg_star.nmax_star);
			nu_l2[en-np_min]=tmp;
		}
	} else{
		nu_l2=cfg_star.nu_nl.row(2).segment(0, cfg_star.Nf_el[2]);
		//std::cout << "nu_l2 = " << nu_l2.transpose() << std::endl;
	}
	// Filter solutions that endup at frequencies higher/lower than the nu_l0 because we will need to extrapolate height/widths otherwise...
	posOK=where_in_range(nu_l2, nu_l0.minCoeff(), nu_l0.maxCoeff(), false);
	tmpXd=nu_l2;
	nu_l2.resize(posOK.size());
	height_l2.resize(posOK.size());
	width_l2.resize(posOK.size());
	for (int i=0; i<posOK.size();i++)
	{
		nu_l2[i]=tmpXd[posOK[i]];
		tmp=lin_interpol(nu_l0, height_l0, nu_l2[i]);
		height_l2[i]=tmp*cfg_star.Vl[1];
		tmp=lin_interpol(nu_l0, width_l0, nu_l2[i]);
		width_l2[i]=tmp;		
	}
	// Assume that the l=2 modes are only sensitive to the envelope rotation
	a1_l2.resize(nu_l2.size());
	a1_l2.setConstant(rot2data.rot_env);
	
	// ------ l=3 modes ----
	if (cfg_star.use_nu_nl == false){ // Do all the stuff related to computing frequencies only if requested
		el=3;
		delta0l_star=-el*(el + 1) * cfg_star.delta0l_percent_star / 100.;
		nu_l3.resize(np_max-np_min);
		for (int en=np_min; en<np_max;en++)
		{
			tmp=asympt_nu_p(cfg_star.Dnu_star, en, cfg_star.epsilon_star, el, delta0l_star, cfg_star.alpha_p_star, cfg_star.nmax_star);
			nu_l3[en-np_min]=tmp;
		}
	} else{
		nu_l3=cfg_star.nu_nl.row(3).segment(0, cfg_star.Nf_el[3]);
		//std::cout << "nu_l3 = " << nu_l3 << std::endl;
	}
	posOK=where_in_range(nu_l3, nu_l0.minCoeff(), nu_l0.maxCoeff(), false);
	tmpXd=nu_l3;
	nu_l3.resize(posOK.size());
	height_l3.resize(posOK.size());
	width_l3.resize(posOK.size());
	for (int i=0; i<posOK.size();i++)
	{
		nu_l3[i]=tmpXd[posOK[i]];
		tmp=lin_interpol(nu_l0, height_l0, nu_l3[i]);
		height_l3[i]=tmp*cfg_star.Vl[2];
		tmp=lin_interpol(nu_l0, width_l0, nu_l3[i]);
		width_l3[i]=tmp;		
	}
	// Assume that the l=3 modes are only sensitive to the envelope rotation
	a1_l3.resize(nu_l3.size());
	a1_l3.setConstant(rot2data.rot_env);//=numpy.repeat(rot_env, len(nu_l3))

	params_out.nu_l0=nu_l0;
	if (cfg_star.use_nu_nl == false){ // Do all the stuff related to computing frequencies only if requested
		params_out.nu_p_l1=freqs.nu_p;
		params_out.nu_g_l1=freqs.nu_g;
	} else{
		params_out.nu_p_l1=nu_p_l1;
		params_out.nu_g_l1=nu_g_l1;
	}
	params_out.nu_l0=nu_l0;
	params_out.nu_m_l1=nu_m_l1;
	params_out.nu_l2=nu_l2;
	params_out.nu_l3=nu_l3;
	params_out.width_l0=width_l0;
	params_out.width_l1=width_l1;
	params_out.width_l2=width_l2;
	params_out.width_l3=width_l3;
	params_out.height_l0=height_l0;
	params_out.height_l1=height_l1;
	params_out.height_l2=height_l2;
	params_out.height_l3=height_l3;
	params_out.a1_l1=a1_l1;
	params_out.a1_l2=a1_l2;
	params_out.a1_l3=a1_l3;
	/*
	std::cout << "nu_l0  = " << nu_l0.transpose() << std::endl;
	std::cout << "nu_l1  = " << nu_m_l1.transpose() << std::endl;
	std::cout << "nu_l2  = " << nu_l2.transpose() << std::endl;
	std::cout << "nu_l3  = " << nu_l3.transpose() << std::endl;	
	std::cout << "width_l0  = " << width_l0 << std::endl;
	std::cout << "width_l1  = " << width_l1 << std::endl;
	std::cout << "width_l2  = " << width_l2 << std::endl;
	std::cout << "width_l3  = " << width_l3 << std::endl;
	std::cout << "height_l0 = " << height_l0 << std::endl;
	std::cout << "height_l1 = " << height_l1 << std::endl;
	std::cout << "height_l2 = " << height_l2 << std::endl;
	std::cout << "height_l3 = " << height_l3 << std::endl;
	std::cout << "a1_l1     = " << a1_l1 << std::endl;
	std::cout << "a1_l2     = " << a1_l2 << std::endl;
	std::cout << "a1_l3     = " << a1_l3 << std::endl;
	*/
	return params_out;
}


Cfg_synthetic_star test_make_synthetic_asymptotic_star_sg(void){
	// Define global Pulsation parameters
	int el=1;

	// Set the current path variables
	std::string cpath=getcwd(NULL, 0);
	
	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params_out;

	cfg_star.output_file_rot=" ";

	cfg_star.Dnu_star=55;
	cfg_star.epsilon_star=0.1;
	cfg_star.delta0l_percent_star=1./100;
	cfg_star.beta_p_star=0.00;

	cfg_star.DPl_star=350.;
	cfg_star.alpha_g_star=0.; // Parameters for g modes that follow exactly the asymptotic relation of g modes for a star with radiative core
	cfg_star.q_star=0.15;

	cfg_star.rot_env_input=30.;
	cfg_star.rot_ratio_input=5.;
	cfg_star.rot_core_input=-1;
	cfg_star.output_file_rot=cpath + "/test.rot";

	cfg_star.maxHNR_l0=5.;

	// Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, 0); // Second argument is the random spread
	cfg_star.fmin=cfg_star.numax_star - 6*cfg_star.Dnu_star;
	cfg_star.fmax=cfg_star.numax_star + 6*cfg_star.Dnu_star;

	cfg_star.nmax_star=cfg_star.numax_star/cfg_star.Dnu_star - cfg_star.epsilon_star;
	cfg_star.alpha_p_star=cfg_star.beta_p_star/cfg_star.nmax_star;

	cfg_star.noise_params_harvey_like.resize(8);
	cfg_star.Vl.resize(3);

	cfg_star.noise_params_harvey_like <<  1., -2. , 0. , 1. ,-1. , 0. , 2.  ,1.;
	cfg_star.Vl << 1.5,0.5, 0.07;
	
	cfg_star.Gamma_max_l0=1;
	cfg_star.Teff_star=-1;
	cfg_star.H0_spread=0;
	//cfg_star.filetemplate=cpath + "/templates/11771760.template";
	cfg_star.filetemplate=cpath + "/templates/Sun.template";

	cfg_star.sigma_p=0;
	cfg_star.sigma_m=0;

	params_out=make_synthetic_asymptotic_star(cfg_star);

	std::cout << " ----- FINAL DIAGNOSTICS ------" << std::endl;
	std::cout << "[en]   type   l    nu_l        w_l         h_l         a1_l" << std::endl;
	for (int en=0; en<params_out.nu_l0.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "p   0    " << params_out.nu_l0[en] <<  "       " << params_out.width_l0[en] <<  "       "  << params_out.height_l0[en] <<    "   0  " << std::endl;
	}
	for (int en=0; en<params_out.nu_p_l1.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "p   1    " << params_out.nu_p_l1[en] <<  "       " << "          -             " <<  "       "  << "          -             " <<  "   -  " << std::endl;
	}
	for (int en=0; en<params_out.nu_g_l1.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "g   1    " << params_out.nu_g_l1[en] <<  "       " << "          -             " <<  "       "  << "          -             " <<  "   -  " << std::endl;
	}
	for (int en=0; en<params_out.nu_m_l1.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "m   1    " << params_out.nu_m_l1[en] <<  "       " << params_out.width_l1[en] <<  "       "  << params_out.height_l1[en] <<  "     " << params_out.a1_l1[en] << std::endl;
	}
	for (int en=0; en<params_out.nu_l2.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "p   2    " << params_out.nu_l2[en] <<  "       " << params_out.width_l2[en] <<  "       "  << params_out.height_l2[en] <<  "    "  << params_out.a1_l2[en] << std::endl;
	}
	for (int en=0; en<params_out.nu_l3.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "p   3    " << params_out.nu_l3[en] <<  "       " << params_out.width_l3[en] <<  "       "  << params_out.height_l3[en] <<  "   "  << params_out.a1_l3[en] << std::endl;
	}
	std::cout << " ------------------------------" << std::endl;
	return cfg_star;
}


Cfg_synthetic_star test_make_synthetic_asymptotic_star_rgb(void){
	// Define global Pulsation parameters
	int el=1;

	// Set the current path variables
	std::string cpath=getcwd(NULL, 0);
	
	Cfg_synthetic_star cfg_star;
	Params_synthetic_star params_out;

	cfg_star.output_file_rot=" ";

	cfg_star.Dnu_star=15;
	cfg_star.epsilon_star=0.1;
	cfg_star.delta0l_percent_star=1./100;
	cfg_star.beta_p_star=0.00;

	cfg_star.DPl_star=80.;
	cfg_star.alpha_g_star=0.; // Parameters for g modes that follow exactly the asymptotic relation of g modes for a star with radiative core
	cfg_star.q_star=0.15;

	cfg_star.rot_env_input=30.;
	cfg_star.rot_ratio_input=5.;
	cfg_star.rot_core_input=-1;
	cfg_star.output_file_rot=cpath + "/test.rot";

	cfg_star.maxHNR_l0=5.;

	// Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	cfg_star.numax_star=numax_from_stello2009(cfg_star.Dnu_star, 0); // Second argument is the random spread
	cfg_star.fmin=cfg_star.numax_star - 6*cfg_star.Dnu_star;
	cfg_star.fmax=cfg_star.numax_star + 6*cfg_star.Dnu_star;

	cfg_star.nmax_star=cfg_star.numax_star/cfg_star.Dnu_star - cfg_star.epsilon_star;
	cfg_star.alpha_p_star=cfg_star.beta_p_star/cfg_star.nmax_star;

	cfg_star.noise_params_harvey_like.resize(8);
	cfg_star.Vl.resize(3);

	cfg_star.noise_params_harvey_like <<  1., -2. , 0. , 1. ,-1. , 0. , 2.  ,1.;
	cfg_star.Vl << 1.5,0.5, 0.07;
	
	cfg_star.Gamma_max_l0=1;
	cfg_star.Teff_star=-1;
	cfg_star.H0_spread=0;
	//cfg_star.filetemplate=cpath + "/templates/11771760.template";
	cfg_star.filetemplate=cpath + "/templates/Sun.template";

	cfg_star.sigma_p=0.0*cfg_star.Dnu_star;
	cfg_star.sigma_m=0.05*cfg_star.Dnu_star;

	params_out=make_synthetic_asymptotic_star(cfg_star);

	std::cout << " ----- FINAL DIAGNOSTICS ------" << std::endl;
	std::cout << "[en]   type   l    nu_l        w_l         h_l         a1_l" << std::endl;
	for (int en=0; en<params_out.nu_l0.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "p   0    " << params_out.nu_l0[en] <<  "       " << params_out.width_l0[en] <<  "       "  << params_out.height_l0[en] <<    "   0  " << std::endl;
	}
	for (int en=0; en<params_out.nu_p_l1.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "p   1    " << params_out.nu_p_l1[en] <<  "       " << "          -             " <<  "       "  << "          -             " <<  "   -  " << std::endl;
	}
	for (int en=0; en<params_out.nu_g_l1.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "g   1    " << params_out.nu_g_l1[en] <<  "       " << "          -             " <<  "       "  << "          -             " <<  "   -  " << std::endl;
	}
	for (int en=0; en<params_out.nu_m_l1.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "m   1    " << params_out.nu_m_l1[en] <<  "       " << params_out.width_l1[en] <<  "       "  << params_out.height_l1[en] <<  "     " << params_out.a1_l1[en] << std::endl;
	}
	for (int en=0; en<params_out.nu_l2.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "p   2    " << params_out.nu_l2[en] <<  "       " << params_out.width_l2[en] <<  "       "  << params_out.height_l2[en] <<  "    "  << params_out.a1_l2[en] << std::endl;
	}
	for (int en=0; en<params_out.nu_l3.size(); en++ )
	{
		std::cout << "[" << en << "]   " << "p   3    " << params_out.nu_l3[en] <<  "       " << params_out.width_l3[en] <<  "       "  << params_out.height_l3[en] <<  "   "  << params_out.a1_l3[en] << std::endl;
	}
	std::cout << " ------------------------------" << std::endl;
	return cfg_star;
}

/*
int main(void){
	Cfg_synthetic_star cfg_star;
	//cfg_star=test_make_synthetic_asymptotic_star_sg();
	cfg_star=test_make_synthetic_asymptotic_star_rgb();
}
*/