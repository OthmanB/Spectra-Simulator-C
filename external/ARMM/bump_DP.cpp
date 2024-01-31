/** 
 * @file bump_DP.cpp
 * @brief Functions computing the full profile of the mixed modes as well as artificial asteroseismic parameters for a RGB stars.
 *
 * This file contains all the functions that describe the Bumped period spacing and how they can be used to derive mode amplitudes from inertia ratio. The functions have been developed and tested based on the following publications:
 * - [Inertia and ksi relation](https://arxiv.org/pdf/1509.06193.pdf)
 * - [Eq. 17 for the rotation - splitting relation](https://www.aanda.org/articles/aa/pdf/2015/08/aa26449-15.pdf)
 * - [Fig 13 and 14 for the evolution - rotation relation in SG and RGB](https://arxiv.org/pdf/1401.3096.pdf)
 * - [Eq. 3 for determining log(g) using Teff and numax, used for getting the evolution stage in conjunction with Fig. 13 from above](https://arxiv.org/pdf/1505.06087.pdf)
 * - [Inertia and Height relation](https://iopscience.iop.org/article/10.1088/2041-8205/781/2/L29/pdf)
 * - [Fig.5, distribution of surface rotation for 361 RGB stars](https://arxiv.org/pdf/1707.05989.pdf)
 */

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "version_solver.h"
#include "data_solver.h"
#include "string_handler.h"
#include "interpol.h"
#include "noise_models.h" // get the harvey_1985 function
#include "solver_mm.h"
#include "bump_DP.h"
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

/**
 * @brief Calculates the ksi function as defined in equation 14 of Mosser+2017 (https://arxiv.org/pdf/1509.06193.pdf).
 *
 * This function calculates the zeta function using the given parameters and equations.
 *
 * @param nu The frequency axis in microHz.
 * @param nu_p The frequency of a p-mode in microHz.
 * @param nu_g The frequency of a g-mode in microHz.
 * @param Dnu_p The large separation, which can be a function of np or nu to account for glitches.
 * @param DPl The period spacing in seconds.
 * @param q The coupling term, without any unit.
 * @return The ksi vector which is the local zeta function.
 */
VectorXd ksi_fct1(const VectorXd& nu, const long double nu_p, const long double nu_g, const long double Dnu_p, const long double DPl, const long double q)
{
	

	const long double pi = M_PI;
	VectorXd cos_upterm, cos_downterm, front_term, tmp(nu.size()), tmp2(nu.size()), ksi(nu.size());

	tmp.setConstant(1./nu_g);
	cos_upterm=pi * 1e6 * (nu.cwiseInverse() - tmp)/DPl;

	tmp.setConstant(nu_p);
	cos_downterm=pi * (nu - tmp) /Dnu_p;
	front_term= 1e-6 * nu.array().square() * DPl / (q * Dnu_p); // relation accounting for units in Hz and in seconds
	tmp2=cos_upterm.array().cos().square()/cos_downterm.array().cos().square();
	tmp.setConstant(1);
	ksi=front_term.cwiseProduct(tmp2);
	ksi=(tmp + ksi).cwiseInverse();
	return ksi;
}

/**
 * @brief Calculates the zeta function for mixed modes in precise mode.
 *
 * This function calculates the ksi_pg vector using the ksi_fct1 function for each combination of nu_p and nu_g.
 *
 * @param nu The frequency vector.
 * @param nu_p The list of p modes in a vector form.
 * @param nu_g The list of g modes in a vector form.
 * @param Dnu_p The local values of the large separation, through a vector.
 * @param DPl The local values of the period spacing in the form of a vector.
 * @param q The coupling strength q between p and g modes.
 * @param norm_method The method for normalizing the ksi_pg vector.
 *                    Possible values: "fast" or "exact". In the exact mode, the ksi_fct2_precise function is used.
 * @return The ksi_pg vector which is the zeta function intervening in mixed modes.
 */
Eigen::VectorXd ksi_fct2(const Eigen::VectorXd& nu, const Eigen::VectorXd& nu_p, const Eigen::VectorXd& nu_g, const Eigen::VectorXd& Dnu_p, const Eigen::VectorXd& DPl, const long double q, const std::string norm_method)
{
    const int Lp = nu_p.size();
    const int Lg = nu_g.size();
    const long double resol = 1e6 / (4 * 365. * 86400.); // Fix the grid resolution to 4 years (converted into microHz)
    Eigen::VectorXd ksi_tmp, ksi_pg(nu.size()), nu_highres, ksi_highres;
    int Ndata;
    long double norm_coef, fmin, fmax;
    ksi_pg.setZero();
    if (norm_method == "fast"){
		for (int np = 0; np < Lp; np++){
			for (int ng = 0; ng < Lg; ng++){
				ksi_tmp = ksi_fct1(nu, nu_p[np], nu_g[ng], Dnu_p[np], DPl[ng], q);
				ksi_pg = ksi_pg + ksi_tmp;
			}
		}
        norm_coef = ksi_pg.maxCoeff();
		ksi_pg = ksi_pg / norm_coef;
        // Ensuring that round-off errors don't lead to values higher than 1...
        for (int i = 0; i < ksi_pg.size(); i++){
            if (ksi_pg[i] > 1){
                ksi_pg[i] = 1;
            }
        }
    } else{
        ksi_pg=ksi_fct2_precise(nu, nu_p, nu_g, Dnu_p, DPl, q);
    }
    return ksi_pg;
}

/**
 * @brief Calculates the zeta function for mixed modes in precise mode.
 *
 * This function calculates the ksi_pg vector using the ksi_fct1 function for each combination of nu_p and nu_g.
 * It also performs parallel computation using OpenMP for better performance.
 *
 * @param nu The frequency vector.
 * @param nu_p The list of p modes in a  vector form.
 * @param nu_g The list of g modes in a vector form.
 * @param Dnu_p The local values of the large separation, through a vector.
 * @param DPl The local values of the period spacing in the form of a vector.
 * @param q The coupling strength q between p and g modes.
 * @return The ksi_pg vector which is the zeta function intervening in mixed modes.
 */
Eigen::VectorXd ksi_fct2_precise(const Eigen::VectorXd& nu, const Eigen::VectorXd& nu_p, const Eigen::VectorXd& nu_g, const Eigen::VectorXd& Dnu_p, const Eigen::VectorXd& DPl, const long double q)
{
	// A slightly more optimized version since 17 Sept 2023. 10% performance increase + proper omp implementation
    const int Lp = nu_p.size();
    const int Lg = nu_g.size();
    const long double resol = 1e6 / (4 * 365. * 86400.); // Fix the grid resolution to 4 years (converted into microHz)
    const long double fmin = (nu_p.minCoeff() >= nu_g.minCoeff()) ? nu_g.minCoeff() : nu_p.minCoeff();
    const long double fmax = (nu_p.maxCoeff() >= nu_g.maxCoeff()) ? nu_p.maxCoeff() : nu_g.maxCoeff();
    const int Ndata = int((fmax - fmin) / resol);
    const Eigen::VectorXd nu_highres = Eigen::VectorXd::LinSpaced(Ndata, fmin, fmax);
    Eigen::VectorXd ksi_pg(nu.size());//, ksi_tmp(nu.size());
    ksi_pg.setZero();
    #pragma omp parallel for shared(ksi_pg)
    for (int np = 0; np < Lp; np++)
    {
        Eigen::VectorXd ksi2_local(nu.size());
        ksi2_local.setZero();
        for (int ng = 0; ng < Lg; ng++){           
            // The function of interest
            ksi2_local += ksi_fct1(nu, nu_p[np], nu_g[ng], Dnu_p[np], DPl[ng], q); //ksi_fct1(nu, nu_p, nu_g, Dnu_p, DPl, q, np, ng);
        }
        #pragma omp critical
        {
            ksi_pg += ksi2_local;
        }  
    }
    Eigen::VectorXd ksi_highres(nu_highres.size());
    ksi_highres.setZero();
    #pragma omp parallel for shared(ksi_highres)
    for (int np = 0; np < Lp; np++)
    {
        Eigen::VectorXd ksi2_highres_local(nu_highres.size());
        ksi2_highres_local.setZero();
        for (int ng = 0; ng < Lg; ng++){           
            // The function of interest
            ksi2_highres_local += ksi_fct1(nu_highres, nu_p[np], nu_g[ng], Dnu_p[np], DPl[ng], q);//ksi_fct1(nu_highres, nu_p, nu_g, Dnu_p, DPl, q, np, ng);
        }
        #pragma omp critical
        {
            ksi_highres += ksi2_highres_local;
        }
    }
    //#pragma omp barrier
    const long double norm_coef = ksi_highres.maxCoeff();
    ksi_pg=ksi_pg/norm_coef;
    
	#pragma omp parallel for
    for (int i = 0; i < nu.size(); i++)
    {
        if (ksi_pg[i] > 1)
        {
            ksi_pg[i] = 1;
        }
    }  
    return ksi_pg;
}


/**
 * @brief Calculates the mode width of mixed modes based on the given ksi_pg, nu_m, nu_p_l0, width_l0, hl_h0_ratio, el, and a factor.
 *
 * This function calculates the width_l vector by performing interpolation using the given ksi_pg, nu_m, nu_p_l0, width_l0, hl_h0_ratio, el, and factor.
 * The width0_at_l is calculated using linear interpolation between nu_p_l0 and width_l0 for each element in nu_m.
 * The width_l vector is then calculated by multiplying width0_at_l with (1. - factor*ksi_pg[i]) and dividing by the square root of hl_h0_ratio[i].
 *
 * @param ksi_pg The ksi_pg vector.
 * @param nu_m The nu_m vector.
 * @param nu_p_l0 The nu_p_l0 vector.
 * @param width_l0 The width_l0 vector.
 * @param hl_h0_ratio The hl_h0_ratio vector.
 * @param el The el value.
 * @param factor The factor value.
 * @return The width_l vector.
 *
 * @note The size of nu_p_l0 and width_l0 must be the same.
 * @note The size of ksi_pg and nu_m must be the same.
 * @note If the sizes are inconsistent, the function will print an error message and exit the program.
 */
VectorXd gamma_l_fct2(const VectorXd& ksi_pg, const VectorXd& nu_m, const VectorXd& nu_p_l0, const VectorXd& width_l0, const VectorXd& hl_h0_ratio, const int el, const long double factor)
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
	}
	return width_l;
}

/**
 * @brief Calculates the H(l)/H(l=0) ratio based on the given ksi_pg and a factor.
 *
 * This function calculates the h_l_rgb vector by subtracting factor*ksi_pg from a vector of ones, and then taking the square root of each element.
 * If any element in the resulting vector is smaller than a given tolerance, it is set to a value close to 0 at the machine precision level to avoid division by 0 or infinity.
 *
 * @param ksi_pg The ksi_pg vector.
 * @param factor The factor value.
 * @return The h_l_rgb vector.
 */
VectorXd h_l_rgb(const VectorXd& ksi_pg, const long double factor)
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
	return hl_h0;
}

/**
 * @brief Reads a template file and extracts relevant information.
 *
 * This function reads a template file and extracts the ID reference, numax reference, Dnu reference, epsilon reference, and data reference from the file. The file should be formatted with key-value pairs separated by a delimiter. Comments in the header are ignored.
 *
 * @param file The path to the template file.
 * @param ignore_errors Flag indicating whether to ignore errors if the required keywords are not found in the file.
 * @return A template_file object containing the extracted information.
 */
template_file read_templatefile(const std::string file, const bool ignore_errors){

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

/**
 * @brief Loads and rescales data from a file based on reference values.
 *
 * This function loads data from a file and rescales it based on reference values. The rescaling is done using interpolation and the extracted information from the template file. The rescaled data is then returned as a Data_2vectXd object.
 *
 * @param nu_star The reference frequencies.
 * @param Dnu_star The reference Dnu value.
 * @param numax_star The reference numax value.
 * @param file The path to the template file.
 * @return A Data_2vectXd object containing the rescaled data.
 */
Data_2vectXd width_height_load_rescale(const VectorXd& nu_star, const long double Dnu_star, const long double numax_star, const std::string file)
{
	int Nref, Nstar;
	long double n_at_numax_ref, height_ref_at_numax, gamma_ref_at_numax, epsilon_star, n_at_numax_star, w_tmp, h_tmp;
	VectorXd tmp, tmp_ref, nu_ref, height_ref, gamma_ref, en_list_ref, en_list_star, w_star, h_star;
	Data_2vectXd out;

	template_file template_data;

	template_data= read_templatefile(file);

	nu_ref=template_data.data_ref.col(0); //[:,0] 
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


/**
 * @brief Computes the core rotation from the surface rotation.
 *
 * This function computes the core rotation based on the given surface rotation and the average core-to-envelope ratio. It is used to achieve a uniform distribution of rotation in the envelope and a uniform population of core-to-envelope ratios.
 *
 * @param rot_envelope The average rotation in the envelope (in microHz).
 * @param core2envelope_star The average rotation in the core.
 * @param output_file_rot The path to the output file to save the computed rotation values.
 * @return A Data_rot2zone object containing the computed rotation values.
 */
Data_rot2zone rot_2zones_v2(const long double rot_envelope, const long double core2envelope_star, std::string output_file_rot)
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

/**
 * @brief Computes the core rotation from the surface rotation.
 *
 * This function computes the core rotation based on the given surface rotation and the average core-to-envelope ratio. It is used to achieve a uniform distribution of rotation in the envelope and a uniform population of core-to-envelope ratios.
 *
 * @param rot_envelope The average rotation in the envelope (in microHz).
 * @param rot_core The average rotation in the core (in microHz).
 * @param output_file_rot The path to the output file to save the computed rotation values.
 * @return A Data_rot2zone object containing the computed rotation values.
 */
Data_rot2zone rot_2zones_v3(const long double rot_envelope, const long double rot_core, std::string output_file_rot)
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

/**
 * @brief Computes the rotation in the envelope based on a truncated Gaussian distribution.
 *
 * This function determines the rotation in the envelope based on a truncated Gaussian distribution. The distribution is inspired by the surface rotation from Ceillier et al. 2017 (https://arxiv.org/pdf/1707.05989.pdf), Fig. 5, which has a skewed distribution ranging from ~30 days to ~160 days with a peak around 60 days. For simplicity, a truncated Gaussian distribution is used with rotation values between 30 and 90 days and a median of 60 days. The truncation happens at sigma. The values are given in days.
 *
 * @param med The median of the distribution (in days).
 * @param sigma The standard deviation of the distribution (in days).
 * @return The rotation frequency in microHz.
 */
long double rot_envelope(long double med, long double sigma)
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

/**
 * @brief Computes the frequency splitting of modes assuming a two-zone rotation profile.
 *
 * This function calculates the frequency splitting of modes based on a two-zone rotation profile. The rotation rates in the envelope and core are given as input parameters. The function also takes the ksi function, which describes the degree of mixture between p and g modes as a function of frequency.
 *
 * @param ksi_pg The ksi function that describes the degree of mixture between p and g modes as a function of frequency.
 * @param rot_envelope The average rotation in the envelope. Must be a scalar.
 * @param rot_core The average rotation in the core. Must be a scalar.
 * @return A VectorXd object containing the frequency splitting values.
 */
VectorXd dnu_rot_2zones(const VectorXd& ksi_pg, const long double rot_envelope, const long double rot_core)
{
	VectorXd re(ksi_pg.size());
	re.setConstant(rot_envelope);

	return ksi_pg*(rot_core/2 - rot_envelope) + re;
}	


/**
 * @brief Computes the value of numax based on the Dnu_star parameter using the relation from Stello+2009.
 *
 * This function calculates the value of numax based on the Dnu_star parameter using the relation from Stello+2009. The relation is given by Dnu ~ 0.263 * numax^0.77. The function also allows for adding a uniform spread around numax if the spread parameter is provided.
 *
 * @param Dnu_star The Dnu_star parameter.
 * @param spread The spread around numax, given as a fraction (e.g., 5% is 0.05). Default value is 0.
 * @return The value of numax.
 */
long double numax_from_stello2009(const long double Dnu_star, const long double spread)
{
	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	
	// Define the frequency range for the calculation by (1) getting numax from Dnu and (2) fixing a range around numax
	const long double beta0=0.263; // according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	const long double beta1=0.77; // according to Stello+2009, we have Dnu_p ~ 0.263*numax^0.77 (https://arxiv.org/pdf/0909.5193.pdf)
	
	long double  nu_max=std::pow(10, std::log10(Dnu_star/beta0)/beta1);
	if (std::abs(spread)>0) // Add a unifrom spread around numax, spread must be given in fraction (eg 5% is 0.05)
	{	
		std::uniform_real_distribution<double> distrib(nu_max*(1.-std::abs(spread)),nu_max*(1. + std::abs(spread)));
		nu_max=distrib(gen);
	} else
	{
		std::cout << "Spread argument in bump_DP.cpp::numax_from_stello2009() is 0 or negative... Ignoring it" << std::endl;
	}
	return nu_max;
}


/**
 * @brief Generates a set of mode parameters for simulating an evolved star based on the asymptotic relations.
 *
 * This function generates a set of mode parameters for simulating an evolved star based on the asymptotic relations. The function takes a structure `cfg_star` as input, which contains all the important parameters for creating the mode parameters. The mode parameters include frequencies, widths, heights, and a-coefficients for each mode of degree l=0, l=1, l=2, and l=3.
 *
 * @param cfg_star A structure that contains all the important parameters for creating the mode parameters, such as (see cfg_star structure itself for the full list), 
 *                 - Dnu_star: Large separation for the p modes.
 *                 - epsilon_star: Phase offset term for the asymptotic relation of the p modes.
 *                 - delta0l_star, alpha_p_star, nmax_star: Parameters for creating 2nd order asymptotic p modes.
 *                 - DPl_star: Period spacing for l=1 g modes.
 *                 - alpha_star: Phase offset term for the asymptotic relation of the g modes.
 *                 - q_star: Coupling coefficient between p and g modes.
 *                 - fmin: Minimum frequency for the modes that should be included in the calculations.
 *                 - fmax: Maximum frequency for the modes that should be included in the calculations.
 * @param legacynoise An optional boolean that is defining how the input noise is dealt with. If legacynoise = True (default), then the fit will use the old definition inspired from Kallinger2014 and will convert that into vectors suitable for the harvey_1985() function. Otherwise, it expects a set of Harvey-like parameters that it will use directly using the harvey_like() function. This latest solution is prefered as it can be jointly used with the pure derivation of these parameters from the Kallinger2014 relations. In that case, the conversion happens outside this function, which becomes "noise-agnostic"
 * @return A structure `params_out` that contains the mode parameters for simulating the evolved star.
 *         - nu_lx: Frequencies of the l=x modes, where x is between 0 and 3.
 *         - nu_p_l1: Base p mode frequencies used to build the frequencies for the l=1 mixed modes.
 *         - nu_g_l1: Base g mode frequencies used to build the frequencies for the l=1 mixed modes.
 *         - width_lx: Widths of the l=x modes, where x is between 0 and 3.
 *         - height_lx: Heights of the l=x modes, where x is between 0 and 3.
 *         - aj_lx: The a-coefficients of order j for each mode of degree l=x.
 */
Params_synthetic_star make_synthetic_asymptotic_star(Cfg_synthetic_star cfg_star, const bool legacynoise)
{
	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> distrib(0 , 1);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen_m(seed); 
	std::normal_distribution<double> distrib_m(0.,cfg_star.sigma_m);


	int en, el, np_min, np_max;
	long double r, tmp, resol, c, xmin ,xmax, delta0l_star;
	VectorXi posOK;
	VectorXd tmpXd, noise_params_harvey1985(4), noise_l0, hmax_l0, Dnu_p, DPl, ksi_pg, h1_h0_ratio;
		//nu_p_l1, nu_g_l1, 
	VectorXd nu_l0, nu_m_l1, nu_l2, nu_l3, 
		height_l0, height_l1, height_l2, height_l3, height_l1p, 
		width_l0, width_l1, width_l2, width_l3,
		a1_l1, a1_l2, a1_l3, // Average rotations
		a2_l1, a2_l2, a2_l3, a4_l2, a4_l3, a6_l3, // Activity or Magnetic effects
		a3_l2, a3_l3, a5_l3; // Latitudinal differential rotation effects

	Data_2vectXd width_height_l0;
	Data_rot2zone rot2data;
	Params_synthetic_star params_out;
	Data_eigensols freqs;

	// Allow to debug content of cfg_star
	//displayCfgSyntheticStar(cfg_star); 
	if(legacynoise == true){
		//noise_params_harvey_like=[A_Pgran ,  B_Pgran , C_Pgran   ,  A_taugran ,  B_taugran  , C_taugran    , p      N0] // 
		noise_params_harvey1985[0] = cfg_star.noise_params_harvey_like[0] * std::pow(cfg_star.numax_star*1e-6,cfg_star.noise_params_harvey_like[1]) + cfg_star.noise_params_harvey_like[2]; // Granulation Amplitude
		noise_params_harvey1985[1] = cfg_star.noise_params_harvey_like[3] * std::pow(cfg_star.numax_star*1e-6, cfg_star.noise_params_harvey_like[4]) + cfg_star.noise_params_harvey_like[5]; // Granulation timescale (in seconds)
		noise_params_harvey1985[0] = noise_params_harvey1985[0]/noise_params_harvey1985[1];
		noise_params_harvey1985[1]= noise_params_harvey1985[1]/1000.;

		noise_params_harvey1985[2]=cfg_star.noise_params_harvey_like[6];
		noise_params_harvey1985[3]=cfg_star.noise_params_harvey_like[7];
	}
	
	// Fix the resolution to 4 years (converted into microHz)
	resol=1e6/(4*365.*86400.);
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
	nu_l0.resize(np_max-np_min);
	for (en=np_min; en<np_max; en++)
	{
		tmp=asympt_nu_p(cfg_star.Dnu_star, en, cfg_star.epsilon_star, 0, 0, cfg_star.alpha_p_star, cfg_star.nmax_star);
		nu_l0[en-np_min]=tmp;
	}
	width_height_l0=width_height_load_rescale(nu_l0, cfg_star.Dnu_star, cfg_star.numax_star, cfg_star.filetemplate); // Function that ensure that Hmax and Wplateau is at numax
	width_l0=width_height_l0.vecXd1;
	height_l0=width_height_l0.vecXd2;
	noise_l0.resize(nu_l0.size());
	noise_l0.setZero();
	if(legacynoise == true){
		noise_l0=harvey1985(noise_params_harvey1985, nu_l0, noise_l0, 1); // Iterate on Noise_l0 to update it by putting the noise profile with one harvey profile
	} else{
		int Nharvey=(cfg_star.noise_params_harvey_like.size() -1)/3; // Assumes that we have Nharvey + white noise
		noise_l0=harvey_like(cfg_star.noise_params_harvey_like, nu_l0, noise_l0, Nharvey); 
	}	
	
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
	el=1;
	delta0l_star=-el*(el + 1) * cfg_star.delta0l_percent_star / 100.;	
	freqs=solve_mm_asymptotic_O2p(cfg_star.Dnu_star, cfg_star.epsilon_star, el, delta0l_star, cfg_star.alpha_p_star, cfg_star.nmax_star, cfg_star.DPl_star, 
								  cfg_star.alpha_g_star, cfg_star.q_star, cfg_star.sigma_p, cfg_star.fmin, cfg_star.fmax, resol, true, false);
	const u_int8_t neverfail = 1;
	// Case where we return p modes only to avoid any crash further in the code... 
	// But this is a MS star, so we Show a Warning about this
	if (freqs.nu_m.size() == 0){ 
		std::cout << " ---------------------------------------------------------" << std::endl;
		std::cout << "    WARNING:  THE COMPUTED STAR IS NOT A RGB OR A SUBGIANT" << std::endl;
		std::cout << "    This is because Dnu and DP1 are both high leading to no mixed modes around numax" << std::endl;
		std::cout << "    The input parameters are:" << std::endl;
		std::cout << "        - cfg_star.Dnu_star     = " << cfg_star.Dnu_star << std::endl;
		std::cout << "        - cfg_star.epsilon_star = " << cfg_star.epsilon_star << std::endl;
		std::cout << "        - delta0l_star          = " << delta0l_star << std::endl;
		std::cout << "        - cfg_star.alpha_p_star = " << cfg_star.alpha_p_star << std::endl;
		std::cout << "        - cfg_star.nmax_star    = " << cfg_star.nmax_star << std::endl;
		std::cout << "        - cfg_star.DPl_star     = " << cfg_star.DPl_star << std::endl;
		std::cout << "        - cfg_star.alpha_g_star = " << cfg_star.alpha_g_star << std::endl;
		std::cout << "        - cfg_star.q_star       = " << cfg_star.q_star << std::endl;
		std::cout << "        - cfg_star.sigma_p      = " <<  cfg_star.sigma_p << std::endl;
		std::cout << "        - cfg_star.fmin         = " <<  cfg_star.fmin << std::endl;
		std::cout << "        - cfg_star.fmax         = " <<  cfg_star.fmax << std::endl;
		std::cout << "        - resol                 = " <<  resol << std::endl;
		std::cout << "     With resulting p and g mode frequencies:" << std::endl;
		std::cout << "        - nu_p : " << freqs.nu_p.transpose() << std::endl;
		std::cout << "        - nu_g : " << freqs.nu_g << std::endl;
		if (neverfail == 0){ // Case where we exit when facing this edge case
		    std::cerr << "     To avoid this, provide smaller Dnu or DP1 so that the density of g modes around" <<std::endl;
			std::cerr << "     the p modes expected to be visible (within [fmin,fmax]) is sufficient to ensure accurate computation" << std::endl;
			std::cerr << "     EXITING..." << std::endl;
			exit(EXIT_FAILURE);
		}
		if (neverfail == 1){ // Case where we go straight to the end of the program and skip the computation. outputs will be empty so that you can skip the whole star generation
			std::cerr << "     PURSUING BUT SKIPPING THE STAR COMPUTATION" << std::endl;
			std::cerr << "     The output structure will be mostly empty and the star should be entirely skipped..." << std::endl;
			params_out.failed=true;
			return params_out;
		}
		if (neverfail == 2){ // Case where the star is computed as if it was a MS star
			freqs.nu_m=freqs.nu_p;
			std::cout << "     PURSUING..." << std::endl;
		}
		std::cout << " ---------------------------------------------------------" << std::endl;
	}
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
	h1_h0_ratio=h_l_rgb(ksi_pg, cfg_star.Hfactor); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019). Hfactor Added on May 2, 2022
	height_l1p.resize(nu_m_l1.size());
	for (int i=0; i<nu_m_l1.size();i++)
	{
		tmp=lin_interpol(nu_l0, height_l0, nu_m_l1[i]);
		height_l1p[i]=tmp;
	}
	
	height_l1p=height_l1p*cfg_star.Vl[0];
	height_l1=h1_h0_ratio.cwiseProduct(height_l1p);
	width_l1=gamma_l_fct2(ksi_pg, nu_m_l1, nu_l0, width_l0, h1_h0_ratio, el, cfg_star.Wfactor); //Wfactor Added on May 2, 2022
	
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
	el=2;
	delta0l_star=-el*(el + 1) * cfg_star.delta0l_percent_star / 100.;
	nu_l2.resize(np_max-np_min);
	for (int en=np_min; en< np_max;en++)
	{
		tmp=asympt_nu_p(cfg_star.Dnu_star, en, cfg_star.epsilon_star, el, delta0l_star, cfg_star.alpha_p_star, cfg_star.nmax_star);
		nu_l2[en-np_min]=tmp;
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
	el=3;
	delta0l_star=-el*(el + 1) * cfg_star.delta0l_percent_star / 100.;
	nu_l3.resize(np_max-np_min);
	for (int en=np_min; en<np_max;en++)
	{
		tmp=asympt_nu_p(cfg_star.Dnu_star, en, cfg_star.epsilon_star, el, delta0l_star, cfg_star.alpha_p_star, cfg_star.nmax_star);
		nu_l3[en-np_min]=tmp;
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
	
	//-----  ADDED ON 13 Sept -----
	// Implementation of the latitudinal differential rotation for outer layers
	// This uses the new substructure env_lat_dif_rot that has all its values initialised to 0 
	// by default. This dummy value is used to identify the different scenarios.
	// Implementation of a3 and a5 is possible in various ways. 
	// However, we recommend using only two situations due to the quality of the available data
	// (unless your a3,a5 values come from eg a rotational model). that is Inside cfg_star.env_lat_dif_rot Either:
	//      - a3_l2, a3_l3 and a5_l3 are set to 0, then this is a case without differential rotation
	//      - a3_l2 = a3_l3 set to some values ~ few percent of a1. And keep a5_l3 =0
	a3_l2.resize(nu_l2.size());
	a3_l3.resize(nu_l3.size());
	a5_l3.resize(nu_l3.size());
	a3_l2.setConstant(cfg_star.env_lat_dif_rot.a3_l2);
	a3_l3.setConstant(cfg_star.env_lat_dif_rot.a3_l3);
	a5_l3.setConstant(cfg_star.env_lat_dif_rot.a5_l3);
	// Implementation of asphericity parameters. As for the differential rotation, the recommendation is 
	// to keep it to 0 (default value), unless you have a model of eg activity and magnetic effects
	a2_l1.resize(nu_m_l1.size());
	a2_l2.resize(nu_l2.size());
	a2_l3.resize(nu_l3.size());
	a4_l2.resize(nu_l2.size());
	a4_l3.resize(nu_l3.size());
	a6_l3.resize(nu_l3.size());
	a2_l1.setConstant(cfg_star.env_aspher.a2_l1);
	a2_l2.setConstant(cfg_star.env_aspher.a2_l2);
	a2_l3.setConstant(cfg_star.env_aspher.a2_l3);
	a4_l2.setConstant(cfg_star.env_aspher.a4_l2);
	a4_l3.setConstant(cfg_star.env_aspher.a4_l3);
	a6_l3.setConstant(cfg_star.env_aspher.a6_l3);
	
	// ----- 

	params_out.nu_l0=nu_l0;
	params_out.nu_p_l1=freqs.nu_p;
	params_out.nu_g_l1=freqs.nu_g;
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
	params_out.a2_l1=a2_l1;
	params_out.a2_l2=a2_l2;
	params_out.a2_l3=a2_l3;
	params_out.a3_l2=a3_l2;
	params_out.a3_l3=a3_l3;
	params_out.a4_l2=a4_l2;
	params_out.a4_l3=a4_l3;
	params_out.a5_l3=a5_l3;
	params_out.a6_l3=a6_l3;

	return params_out;
}


/**
 * @brief Displays the values of the parameters in the Cfg_synthetic_star structure.
 *
 * This function displays the values of the parameters in the Cfg_synthetic_star structure. The parameters include Teff_star, numax_star, Dnu_star, epsilon_star, and other parameters related to the synthetic star simulation.
 *
 * @param cfg The Cfg_synthetic_star structure containing the parameters.
 */
void displayCfgSyntheticStar(const Cfg_synthetic_star& cfg) {
    std::cout << "Teff_star: " << cfg.Teff_star << std::endl;
    std::cout << "numax_star: " << cfg.numax_star << std::endl;
    std::cout << "Dnu_star: " << cfg.Dnu_star << std::endl;
    std::cout << "epsilon_star: " << cfg.epsilon_star << std::endl;
    std::cout << "delta0l_percent_star: " << cfg.delta0l_percent_star << std::endl;
    std::cout << "beta_p_star: " << cfg.beta_p_star << std::endl;
    std::cout << "alpha_p_star: " << cfg.alpha_p_star << std::endl;
    std::cout << "nmax_star: " << cfg.nmax_star << std::endl;
    std::cout << "DPl_star: " << cfg.DPl_star << std::endl;
    std::cout << "alpha_g_star: " << cfg.alpha_g_star << std::endl;
    std::cout << "q_star: " << cfg.q_star << std::endl;
    std::cout << "fmin: " << cfg.fmin << std::endl;
    std::cout << "fmax: " << cfg.fmax << std::endl;
    std::cout << "maxHNR_l0: " << cfg.maxHNR_l0 << std::endl;
    std::cout << "noise_params_harvey_like: " << cfg.noise_params_harvey_like << std::endl;
    std::cout << "Gamma_max_l0: " << cfg.Gamma_max_l0 << std::endl;
    std::cout << "rot_env_input: " << cfg.rot_env_input << std::endl;
    std::cout << "rot_ratio_input: " << cfg.rot_ratio_input << std::endl;
    std::cout << "rot_core_input: " << cfg.rot_core_input << std::endl;
    std::cout << "env_lat_dif_rot.a3_l2: " << cfg.env_lat_dif_rot.a3_l2 << std::endl;
    std::cout << "env_lat_dif_rot.a3_l3: " << cfg.env_lat_dif_rot.a3_l3 << std::endl;
    std::cout << "env_lat_dif_rot.a5_l3: " << cfg.env_lat_dif_rot.a5_l3 << std::endl;
    std::cout << "env_aspher.a2_l1: " << cfg.env_aspher.a2_l1 << std::endl;
    std::cout << "env_aspher.a2_l2: " << cfg.env_aspher.a2_l2 << std::endl;
    std::cout << "env_aspher.a2_l3: " << cfg.env_aspher.a2_l3 << std::endl;
    std::cout << "env_aspher.a4_l2: " << cfg.env_aspher.a4_l2 << std::endl;
    std::cout << "env_aspher.a4_l3: " << cfg.env_aspher.a4_l3 << std::endl;
    std::cout << "env_aspher.a6_l3: " << cfg.env_aspher.a6_l3 << std::endl;
    std::cout << "output_file_rot: " << cfg.output_file_rot << std::endl;
    std::cout << "Vl: " << cfg.Vl.transpose() << std::endl;
    std::cout << "H0_spread: " << cfg.H0_spread << std::endl;
    std::cout << "filetemplate: " << cfg.filetemplate << std::endl;
    std::cout << "sigma_p: " << cfg.sigma_p << std::endl;
    std::cout << "sigma_m: " << cfg.sigma_m << std::endl;
    std::cout << "Hfactor: " << cfg.Hfactor << std::endl;
    std::cout << "Wfactor: " << cfg.Wfactor << std::endl;
    //std::cout << "inclination: " << cfg.inclination << std::endl;
    //std::cout << "nu_nl: " << cfg.nu_nl << std::endl;
    //std::cout << "Nf_el: " << cfg.Nf_el.transpose() << std::endl;
    //std::cout << "use_nu_nl: " << cfg.use_nu_nl << std::endl;
}