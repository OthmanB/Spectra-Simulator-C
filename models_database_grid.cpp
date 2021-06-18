/*
 * models_database.cpp
 *
 * Header file that contains all kind of methods
 * used to generate models for the pulsation/noise
 * 
 *  Created on: 20 Apr 2016
 *      Author: obenomar
 */
# include <iostream>
# include <iomanip>
#include <fstream>
# include <Eigen/Dense>
#include "models_database_grid.h"
#include "interpol.h"
#include "linfit.h"
#include "noise_models.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

/*
void generate_cfg_asymptotic_act_asym_Hgauss(VectorXd input_params, std::string file_out_modes, std::string file_out_noise){

	// ----- Constants ------
	const double PI = 3.141592653589793238462643;
	const double G=6.667e-8;
	const double Teff_sun=5777;
	const double Dnu_sun=135.1;
	const double numax_sun=3150;
	const double R_sun=6.96342e5;
	const double M_sun=1.98855e30;
	const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);

        //std::cout << "      generate_cfg_asymptotic_act_asym_Hgauss" << std::endl;

	VectorXd Visibilities(4);
	Visibilities << 1, 1.5, 0.5, 0.08;
	// ----------------------

	double ks=2; // controls the width of the gaussian for heights

	//std::cout << "input_params=" << input_params.transpose() << std::endl;
	//std::cout << "input_params.size()=" << input_params.size() << std::endl;
	// ------- Deploy the parameters ------
	double numax=input_params[0];
	double Dnu=input_params[1];
	double epsilon=input_params[2];
	double D0=input_params[3];
	double Max_Height=input_params[4];
	double Width=input_params[5];
	double lmax=input_params[6];
	int    Nmax=input_params[7];
	double a1=input_params[8];
	double a3=input_params[9];
	double b=input_params[10];
	double alfa=input_params[11];
	double beta_asym=input_params[12];
	double inc=input_params[13];	
	VectorXd input_noise(9);
	if((input_params.size() - 14) == 9){
		input_noise=input_params.segment(14, 9);
	} else{ // If we did not provide 3 lines of 3 elements, then this means that we have 2 lines of 3 elements + 1 line with one elements (7 parameters)
		input_noise.segment(0, 7)=input_params.segment(14, 7);
		input_noise[7]=-2;
		input_noise[8]=-2;
	}
	// ------------------------------------

    //std::cout << "           - Variables..." << std::endl;

	// --------- Variables ---------
	int k;
	double el, en, n_at_numax, height, eta;
	VectorXd en_list(Nmax);
	MatrixXd nu(int(lmax+1), Nmax), h(int(lmax+1), Nmax), w(int(lmax+1), Nmax), s_a1(int(lmax+1), Nmax), s_eta(int(lmax+1), Nmax), 
		 s_a3(int(lmax+1), Nmax), s_asym(int(lmax+1), Nmax), s_b(int(lmax+1), Nmax), s_alfa(int(lmax+1), Nmax), i(int(lmax+1), Nmax), 
		 mode_params(int(lmax+1)*Nmax, 11), noise_params(3,3);
	// -----------------------------

	if((Nmax % 2) !=0){
		std::cout << "Nmax must be a odd number. Please change that value accordingly" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	k=0;
	n_at_numax=numax/Dnu - epsilon + el/2 + el*(el+1)*D0/Dnu;
	for(en=-Nmax/2 + 1; en<=Nmax/2; en++){
		en_list[k]=floor(n_at_numax) + en;
		k=k+1;
	}	

	// Define the centrifugal force effect	
	eta=(4./3)*PI * pow( a1*1e-6 ,2) / (G * rho_sun) * pow(Dnu_sun/Dnu,2);
	//std::cout << "Fixing centrifugal term eta = " << eta << std::endl;
	
	// Create a list of frequencies, Height, Width, Splitting, Centrifugal terms, latitudinal terms and stellar inclination

    //std::cout << "           - Parameters..." << std::endl;
    
	for(el=0; el<=lmax; el++){
		for(en=0; en<Nmax; en++){
			if(el == 0 || el == 1){
				nu(el,en)=(en_list[en] + epsilon + el/2)*Dnu - el*(el + 1)*D0;
			}
			if(el == 2 || el == 3){
				nu(el,en)=(en_list[en]-1 + epsilon + el/2)*Dnu - el*(el + 1)*D0;
			}
			if(el > 3){
				std::cout << "Generating frequencies with degree higher than l=3 is not implemented" << std::endl;
				std::cout << "Please set lmax<4" << std::endl;
				std::cout << "The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
			}
			
			height=Max_Height*exp(-0.5 * pow( (nu(el, en) - numax)/(ks*Dnu), 2));
			h(el, en)=height*Visibilities[el];

			w(el,en)=Width;

			s_a1(el,en)=a1;
			s_eta(el,en)=eta;
			s_a3(el,en)=a3;
			s_b(el,en)=b;
			s_alfa(el,en)=alfa;
			s_asym(el,en)=beta_asym;
			i(el,en)=inc;
		}
	}

    //std::cout << "           - mode_params preparation..." << std::endl;

	// Summarizing the information into a parameter table
	for(el=0; el<=lmax; el++){
		for(k=el*Nmax; k<(el+1)*Nmax; k++){
			mode_params(k,0)=el;
			mode_params(k,1)=nu(el , k-el*Nmax);
			mode_params(k,2)=h(el , k-el*Nmax);
			mode_params(k,3)=w(el , k-el*Nmax);
			mode_params(k,4)=s_a1(el , k-el*Nmax);
			mode_params(k,5)=s_eta(el , k-el*Nmax);
			mode_params(k,6)=s_a3(el , k-el*Nmax);
			mode_params(k,7)=s_b(el , k-el*Nmax);
			mode_params(k,8)=s_alfa(el , k-el*Nmax);
			mode_params(k,9)=s_asym(el , k-el*Nmax);
			mode_params(k,10)=i(el , k-el*Nmax);
		}
	}
	

	// A FUNCTION THAT WRITES THE PARAMETERS
	write_star_mode_params_act_asym(mode_params, file_out_modes);

    //std::cout << "           - Noise..." << std::endl;

	// Defining the noise profile parameters
	// Note about the noise: -1 means that it is ignored. -2 mean that the value is irrelevant
	for(int e=0; e<3; e++){
		for(int k=0; k<3; k++){
			noise_params(e, k)=input_noise(3*e + k);
		}
	}
	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(noise_params, file_out_noise);

    //std::cout << "           - Exit" << std::endl;


}
*/


/*
*  This procedure :
*		- Use a set of frequencies from models (given in input_model)
*		- Use a reference star to generate/rescale Widths and Heights profiles. Visibilities are calculated there. The input format is a simple Matrix-formated table
*               - Implement a1, inclination
*		- Implement Harvey profile noise (+ White Noise) as defined by Karoff et al. 2010
*               - MUST BE MODIFIED IF YOU WISH: a2, a3, asymmetry
*/
void generate_cfg_from_refstar_HWscaled(VectorXd input_params, Model_data input_model, std::string file_ref_star, std::string file_out_modes, std::string file_out_noise){

	// ----- Constants ------
	const double PI = 3.141592653589793238462643;
	const double G=6.667e-8;
	const double Teff_sun=5777;
	const double Dnu_sun=135.1;
	const double numax_sun=3150;
	const double R_sun=6.96342e5;
	const double M_sun=1.98855e30;
	const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);

	const std::string delimiter=" ";
	const bool verbose_data=0;

	const int lmax=3, Nmax_min=10; // We simulate l=0,1,2,3

        //std::cout << "      generate_cfg_from_refstar_HWscaled" << std::endl;

	VectorXd Visibilities(4);
	Visibilities.setOnes();
	Data_Nd data_ref;
	//Visibilities << 1, 1.5, 0.5, 0.08;

	// ----- Used Constant and variables of the reference star -----
	const int pref_freq=1, pref_h=2, pref_w=3, pref_n=11; // positions in the reference table for frequencies, heights, widths and the local noise
	int Nmax, k, el, en;
	double dnu_ref, numax_ref, norm;
	std::vector<double> pos;
	VectorXi pos2;
	VectorXd  tmp, tmp2, ones, x_ref, height_ref, fit;
	MatrixXd nu_ref, w_ref, h_ref, a_ref, hnr_ref;
	// --------- Variables used for the simulated star ---------
	std::string varname;
	int Nmax_star;
	double c, dnu_star, numax_star;
	VectorXd x, gamma, height, hnr, noise_star;
	MatrixXd nu, h, w, s_a1, s_eta, s_a3, s_asym, s_b, s_alfa, inclination, mode_params;
	// -----------------------------

	// ------- Deploy the parameters of the simulation given by input_params ------
	double Mind=input_params[0];
	double maxHNR=input_params[1];
	double Gamma_Hmax_star=input_params[2];
	double a1=input_params[3];
	double inc=input_params[4];
	MatrixXd input_noise(3,3);
	
	input_noise.row(0)=input_params.segment(5,3);
	input_noise.row(1)=input_params.segment(8,3);
	input_noise(2, 0)=input_params[11];
	input_noise(2, 1)=-2;
	input_noise(2, 2)=-2;

	//std::cout << "input_noise=" << input_noise << std::endl;

	// -------- Extracting numax_star and Dnu_star -------
	varname="dnu";
	pos2=where_strXi(input_model.labels_params, varname);
	dnu_star=input_model.params[pos2[0]];
	varname="numax";
	pos2=where_strXi(input_model.labels_params, varname);
	numax_star=input_model.params[pos2[0]];
	// ------------------------------------

	std::cout << "      dnu_star=" << dnu_star << std::endl;
	std::cout << "      numax_star=" << numax_star << std::endl;

	// ------- Read the parameters of the reference star ----------
	//std::cout << file_ref_star << std::endl;
	data_ref=read_data_ascii_Ncols(file_ref_star, delimiter, verbose_data);

	//std::cout << "data_ref=" << std::endl;
	//std::cout << data_ref.data << std::endl;
	//exit(EXIT_SUCCESS);

	// Sorting usefull data, perform requirement checks and compute mode HNR	
	if (data_ref.data.col(0).maxCoeff() < lmax){
		std::cout << " --------- Reference star inputs: Requirements not fullfilled --------- " << std::endl;
		std::cout << " This spectrum simulator requires reference widths, height and noise level" << std::endl;
		std::cout << " for modes of degree l=0,1,2,3. The current lmax=" << data_ref.data.col(0).maxCoeff() << std::endl;
		std::cout << " is unsuficient. Please give valid reference parameters for the reference star" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_SUCCESS);
	}
	for(el=0; el<lmax+1;el++){
		pos=where(data_ref.data.col(0), "=", el, 0);
		if(el == 0){
			if(pos.size() >= Nmax_min){
				Nmax=pos.size()+2; // We add 2 because of boundary conditions: nu=0 ==> w=h=a=hnr=0. Same at nu=10000.
				std::cout << "Nmax=" << Nmax << std::endl;
				nu_ref.resize(lmax+1, Nmax);
				w_ref.resize(lmax+1, Nmax);
				h_ref.resize(lmax+1, Nmax);
				a_ref.resize(lmax+1, Nmax);
				hnr_ref.resize(lmax+1, Nmax);
			} else{
				std::cout << " --------- Reference star inputs: Requirements not fullfilled --------- " << std::endl;
				std::cout << " The number of radial order should be at least Nmax=" << Nmax_min << std::endl;
				std::cout << " Please add modes for the reference star " << std::endl;
				std::cout << " The program will exit now" << std::endl;
				exit(EXIT_SUCCESS);
			}
		} else{
			if(pos.size() != Nmax-2){
				std::cout << " --------- Reference star inputs: Requirements not fullfilled --------- " << std::endl;
				std::cout << " The number of radial order for each degree should be the same " << std::endl;
				std::cout << " Please check your reference parameters for the reference star " << std::endl;
				std::cout << " The program will exit now" << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		//std::cout << "pos.size()=" << pos.size() << std::endl;
		// ----- Imposing boundary conditions ------
		nu_ref(el,0)=0;
		w_ref(el,0)=0;
		h_ref(el,0)=0;
		a_ref(el,0)=0;
		hnr_ref(el,0)=0;
		nu_ref(el,Nmax-1)=10000.; // Value at which we should not expect to detect pulsations anytime soon...
		w_ref(el,Nmax-1)=50;
		h_ref(el,Nmax-1)=0.0001;
		a_ref(el,Nmax-1)=0;
		hnr_ref(el,Nmax-1)=0;

		// ----------------------------------------
		en=1;
		for(int ii=0; ii<pos.size();ii++){
			nu_ref(el,en)=data_ref.data(pos[ii],pref_freq);
			w_ref(el,en)=data_ref.data(pos[ii],pref_w);
			h_ref(el,en)=data_ref.data(pos[ii],pref_h);
			a_ref(el,en)=sqrt(PI*h_ref(el,en)*w_ref(el,en)); // mode amplitude... used for determining numax(el)
			hnr_ref(el,en)=h_ref(el,en)/data_ref.data(pos[ii],pref_n);
			//std::cout << "(" << el << "," << en << ")" << data_ref.data(pos[ii],pref_freq) << " " << data_ref.data(pos[ii],pref_w) << " " << data_ref.data(pos[ii],pref_h) << " " << sqrt(PI*h_ref(el,en)*w_ref(el,en)) << " " << h_ref(el,en)/data_ref.data(pos[ii],pref_n) << std::endl;
			en=en+1;
		}
	}

	//std::cout << "hnr_ref=" << std::endl;
	//std::cout << hnr_ref << std::endl;
	// Getting numax
	numax_ref=0;
	norm=0;
	for(el=0;el<lmax+1;el++){
		pos=where(a_ref.row(el), "=", a_ref.row(el).maxCoeff(), 0);
		numax_ref=numax_ref + nu_ref(el, pos[0]) * a_ref.row(el).maxCoeff(); // Weighted mean
		norm=norm + a_ref.row(el).maxCoeff();
		//std::cout << "numax(l="<< el << ")=" << nu_ref(el, pos[0]) <<std::endl;
	}
	numax_ref=numax_ref/norm;
	//std::cout << "mean numax=" << numax_ref << std::endl;

	// Getting \Delta\nu
	ones.resize(Nmax-2);
	ones.setOnes();

	std::cout << "      numax_ref=" << numax_ref << std::endl;
	//std::cout << "Chosen numax_ref=" << numax_ref << std::endl;

	//std::cout << nu_ref << std::endl; ///numax_ref- ones << std::endl;
	//exit(EXIT_SUCCESS);

	tmp=nu_ref.row(0).segment(1,Nmax-2);
	pos=where(tmp/numax_ref-ones, ">=", -0.25, 0);// AND nu_l0/numax_star-1 le 0.2);
	tmp2.resize(pos.size());
	for(int i=0; i<pos.size(); i++){tmp2[i]=tmp[pos[i]];}
	//std::cout << "tmp2=" << tmp2 << std::endl;	
	ones.setOnes(tmp2.size());
	pos=where(tmp2/numax_ref-ones, "<=", 0.25, 0);
	tmp.resize(pos.size());
	for(int i=0; i<pos.size(); i++){tmp[i]=tmp2[pos[i]];}
	//std::cout << "tmp=" << tmp << std::endl;
	if(tmp2.size() <= 3){
		std::cout << "Warning: small vector of frequencies detected for Mind=" << Mind << std::endl;
		std::cout << "         using the full dataset of frequencies for calculating Deltanu" << std::endl;
		//exit(EXIT_SUCCESS);
		tmp2.resize(nu_ref.col(0).size());
		tmp2=nu_ref.col(0);
	}

	tmp.setLinSpaced(tmp2.size(), 0, tmp2.size()-1);
	//std::cout << "fitting" << std::endl;
	fit=linfit(tmp, tmp2); // Use a linear fit of the modes around numax to get Delta\nu
	dnu_ref=fit[0];
	std::cout << "      dnu_ref=" << dnu_ref << std::endl;

	// --------- Variables size allocation ---------
	Nmax_star=10000.; // large value to begin with
	for(el=0; el<lmax+1; el++){
		//std::cout << "el=" << el << std::endl;
		pos=where(input_model.freqs.col(0), "=", el, 0);
		if (pos.size() < Nmax_star){ Nmax_star=pos.size();} // Keep the minimal number of radial order from the list of modes (ensure that Nmax_star is same for all l)
		//std::cout << "Nmax_star=" << Nmax_star << std::endl;
	}
	//std::cout << "Nmax_star=" << Nmax_star << std::endl;

	nu.resize(lmax+1, Nmax_star), h.resize(lmax+1, Nmax_star), w.resize(lmax+1, Nmax_star), s_a1.resize(lmax+1, Nmax_star), s_eta.resize(lmax+1, Nmax_star), 
		 s_a3.resize(lmax+1, Nmax_star), s_asym.resize(lmax+1, Nmax_star), s_b.resize(lmax+1, Nmax_star), s_alfa.resize(lmax+1, Nmax_star), inclination.resize(lmax+1, Nmax_star), 
		 mode_params.resize((lmax+1)*Nmax_star, 11);
	// Initialize variables for second order effects to neutral values (e.g. 0)
	s_eta.setZero(); //Neglect the centrifugal effects
	s_a3.setZero(); //Neglect the first order effect of the latitudinal rotation
	s_asym.setZero(); //Neglect the mode asymmetry
	s_b.setZero(); // Do not consider any additional frequency-dependent star distorsion
	s_alfa.setOnes(); // Do not consider any additional frequency-dependent star distorsion
	// -----------------------------

	// ---------- Setting up the variables for the simulated star ----------
	for(el=0; el<lmax+1; el++){
		//std::cout << "el=" << el << std::endl;
		// ------- Organise frequencies in a (l,n) table ---------
		pos=where(input_model.freqs.col(0), "=", el, 0);
		//for(int i=0;i<pos.size();i++){
		for(int i=0;i<Nmax_star;i++){
			//std::cout << "[" << i << "/" << Nmax_star-1 << "] " << "l=" << input_model.freqs(pos[i],0) << "pos=" << pos[i]  << " input_model.freqs(2,pos[i])= " << input_model.freqs(pos[i],2) << std::endl;
			nu(el,i)=input_model.freqs(pos[i],2);
		}
		// -------------------------------------------------------

		// ----- Perform rescalings -------
		x_ref.resize(Nmax);
		for(int en=0; en<Nmax;en++){
			x_ref[en]=(nu_ref(el,en)-numax_ref)/dnu_ref;
			//std::cout << x_ref[en] << std::endl;
		}

		//std::cout << "l=" << el << " Before interpols" << std::endl;
		x.resize(Nmax_star);
		for(int en=0; en<Nmax_star;en++){
			x[en]=(nu(el,en)-numax_star)/dnu_star;
			//std::cout << x[en] << std::endl;
		}
		/*
		std::cout << " x_ref      w_ref " << std::endl;
		for(int en=0; en<Nmax; en++){		
			std::cout << en << "  " << x_ref[en] << "  " <<  w_ref(el,en) << std::endl;
		}
		*/
		gamma.resize(Nmax_star);
		//std::cout << " x      gamma " << std::endl;
		for(int en=0; en<Nmax_star; en++){
			gamma[en]=lin_interpol(x_ref, w_ref.row(el), x[en]);
			//std::cout << en << "  " <<  x[en] << "  " << gamma[en] << std::endl;
		}
		/*
		std::cout << " x_ref      h_ref " << std::endl;
		for(int en=0; en<Nmax; en++){		
			std::cout << en << "  " << x_ref[en] << "  " << h_ref(el,en) << std::endl;
		}
		*/
		hnr.resize(Nmax_star);
		//std::cout << " x      h " << std::endl;
		for(int en=0; en<Nmax_star; en++){
			hnr[en]=lin_interpol(x_ref, hnr_ref.row(el), x[en]);
			//std::cout << en << "  " <<  x[en] << "  " << hnr[en] << std::endl;
		}
		if(el==0){ 
			c=maxHNR/hnr.maxCoeff(); // we will use l=0 to compute the new maxHNR
			//std::cout << "c=" << c << std::endl;
		}
		noise_star=harvey_1985(input_noise, nu.row(el)); // the local noise level of the simulated star

		// We solve: max(HNR_star) = c max(HNR_ref) ==> c=max(HNR_star)/max(HNR_ref)
		// We then use c to rescale the heights hstar(nu) noting that: c* href(nu)/Nref(nu) = hstar(nu)/Nstar(nu)
		height.resize(Nmax_star);
		for (int en=0; en<Nmax_star; en++){
			height[en]=hnr[en] * c * noise_star[en];
			//std::cout << height[en]/noise_star[en]/maxHNR << std::endl;
		}
		
		pos=where(height, "=", height.maxCoeff(),0);
		gamma=gamma * Gamma_Hmax_star/gamma(pos[0]); // Rescaling of the Widths in the y-axis. Gamma_star(pos[0]) is the width at Hmax
		// --------------------------

		// --------------- Other parameters setup -----------------
		h.row(el)=height;
		w.row(el)=gamma;
		s_a1.row(el).setConstant(a1);
		inclination.row(el).setConstant(inc);
	}

	// --------------------------------------------------------------------

	// ---------- Summarizing the information into suitable inputs for the writting function -------------
	for(el=0; el<lmax+1; el++){		
		for(k=el*Nmax_star; k<(el+1)*Nmax_star; k++){
			mode_params(k,0)=el;
			mode_params(k,1)=nu(el , k-el*Nmax_star);
			mode_params(k,2)=h(el , k-el*Nmax_star);
			mode_params(k,3)=w(el , k-el*Nmax_star);
			mode_params(k,4)=s_a1(el , k-el*Nmax_star);
			mode_params(k,5)=s_eta(el , k-el*Nmax_star);
			mode_params(k,6)=s_a3(el , k-el*Nmax_star);
			mode_params(k,7)=s_b(el , k-el*Nmax_star);
			mode_params(k,8)=s_alfa(el , k-el*Nmax_star);
			mode_params(k,9)=s_asym(el , k-el*Nmax_star);
			mode_params(k,10)=inclination(el , k-el*Nmax_star);
		}
	}

	// A FUNCTION THAT WRITES THE PARAMETERS
	write_star_mode_params_act_asym(mode_params, file_out_modes);

   	//std::cout << "           - Noise..." << std::endl;

	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(input_noise, file_out_noise);

    	//std::cout << "           - Exit" << std::endl;
}


/*
*  This procedure :
*		- Use a set of frequencies from models (given in input_model)g++ -O3 -I ../eigen -I ../ -fopenmp -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lgsl -lgslcblas
*		- Use a reference star to generate/rescale Widths and Heights profiles. Visibilities are calculated there. The input format is a simple Matrix-formated table
*               - Implement a1, a3, inclination
*		- Implement A Single Harvey profile noise which scales with numax (+ a White Noise). The definition of the Harvey parameters are as defined by Karoff et al. 2010
*                 Recommended coefficients for the scaling are Pgran = A numax^B + C with A=10^-4 and B=-2, C=0. t_gran = A numax^B + C with A=1 and B=-1 and C=0
*/
void generate_cfg_from_refstar_HWscaled_GRANscaled(VectorXd input_params, Model_data input_model, std::string file_ref_star, std::string file_out_modes, std::string file_out_noise){

	// ----- Constants ------
	const double PI = 3.141592653589793238462643;
	const double G=6.667e-8;
	const double Teff_sun=5777;
	const double Dnu_sun=135.1;
	const double numax_sun=3150;
	const double R_sun=6.96342e5;
	const double M_sun=1.98855e30;
	const double rho_sun=M_sun*1e3/(4*PI*pow(R_sun*1e5,3)/3);

	const std::string delimiter=" ";
	const bool verbose_data=0;

	const int lmax=3, Nmax_min=10; // We simulate l=0,1,2,3

        //std::cout << "      generate_cfg_from_refstar_HWscaled_GRANscaled" << std::endl;

	VectorXd Visibilities(4);
	Visibilities.setOnes();
	Data_Nd data_ref;
	//Visibilities << 1, 1.5, 0.5, 0.08;

	// ----- Used Constant and variables of the reference star -----
	const int pref_freq=1, pref_h=2, pref_w=3, pref_n=11; // positions in the reference table for frequencies, heights, widths and the local noise
	int Nmax, k, el, en;
	double dnu_ref, numax_ref, norm;
	std::vector<double> pos;
	VectorXi pos2;
	VectorXd  tmp, tmp2, ones, x_ref, height_ref, fit;
	MatrixXd nu_ref, w_ref, h_ref, a_ref, hnr_ref;
	// --------- Variables used for the simulated star ---------
	std::string varname;
	int Nmax_star;
	double c, dnu_star, numax_star;
	VectorXd x, gamma, height, hnr, noise_star;
	MatrixXd nu, h, w, s_a1, s_eta, s_a3, s_asym, s_b, s_alfa, inclination, mode_params;
	// -----------------------------

	// ------- Deploy the parameters of the simulation given by input_params ------
	double Mind=input_params[0];
	double maxHNR=input_params[1];
	double Gamma_Hmax_star=input_params[2];
	double a1=input_params[3];
	double inc=input_params[4];
	double H, tau, p; // Coefficient for the noise
	MatrixXd input_noise(2,3);
	
	// -------- Extracting numax_star and Dnu_star -------
	varname="dnu";
	pos2=where_strXi(input_model.labels_params, varname);
	dnu_star=input_model.params[pos2[0]];
	varname="numax";
	pos2=where_strXi(input_model.labels_params, varname);
	numax_star=input_model.params[pos2[0]];
	// ------------------------------------

	std::cout << "      dnu_star=" << dnu_star << std::endl;
	std::cout << "      numax_star=" << numax_star << std::endl;
	
	// --------- Scaling the noise according to the given parameters --------
	//std::cout << input_params.segment(5, 7) << std::endl;

	tau=input_params[8] * pow(numax_star*1e-6,input_params[9]) + input_params[10]; // Granulation timescale (in seconds)

	H=input_params[5] * pow(numax_star*1e-6,input_params[6]) + input_params[7]; // Granulation Amplitude
	//std::cout << "P=" << H << std::endl;
	H=H/tau ; //This is due to the used definition for the Harvey profile (conversion from Hz to microHz)
	tau=tau/1000. ; //conversion in ksec
	p=2;// Fixed power law
	input_noise(0,0)=H;
	input_noise(0,1)=tau;
	input_noise(0,2)=p; 
	input_noise(1, 0)=input_params[11]; // White noise
	input_noise(1, 1)=-2;
	input_noise(1, 2)=-2;

	//std::cout << "H=" << H << std::endl;
	//std::cout << "tau=" << tau << std::endl;
	//std::cout << "p=" << p << std::endl;
	//exit(EXIT_SUCCESS);

	// ------- Read the parameters of the reference star ----------
	data_ref=read_data_ascii_Ncols(file_ref_star, delimiter, verbose_data);

	// Sorting usefull data, perform requirement checks and compute mode HNR	
	if (data_ref.data.col(0).maxCoeff() < lmax){
		std::cout << " --------- Reference star inputs: Requirements not fullfilled --------- " << std::endl;
		std::cout << " This spectrum simulator requires reference widths, height and noise level" << std::endl;
		std::cout << " for modes of degree l=0,1,2,3. The current lmax=" << data_ref.data.col(0).maxCoeff() << std::endl;
		std::cout << " is unsuficient. Please give valid reference parameters for the reference star" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_SUCCESS);
	}
	for(el=0; el<lmax+1;el++){
		pos=where(data_ref.data.col(0), "=", el, 0);
		if(el == 0){
			if(pos.size() >= Nmax_min){
				Nmax=pos.size()+2; // We add 2 because of boundary conditions: nu=0 ==> w=h=a=hnr=0. Same at nu=10000.
				nu_ref.resize(lmax+1, Nmax);
				w_ref.resize(lmax+1, Nmax);
				h_ref.resize(lmax+1, Nmax);
				a_ref.resize(lmax+1, Nmax);
				hnr_ref.resize(lmax+1, Nmax);
			} else{
				std::cout << " --------- Reference star inputs: Requirements not fullfilled --------- " << std::endl;
				std::cout << " The number of radial order should be at least Nmax=" << Nmax_min << std::endl;
				std::cout << " Please add modes for the reference star " << std::endl;
				std::cout << " The program will exit now" << std::endl;
				exit(EXIT_SUCCESS);
			}
		} else{
			if(pos.size() != Nmax-2){
				std::cout << " --------- Reference star inputs: Requirements not fullfilled --------- " << std::endl;
				std::cout << " The number of radial order for each degree should be the same " << std::endl;
				std::cout << " Please check your reference parameters for the reference star " << std::endl;
				std::cout << " The program will exit now" << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		// ----- Imposing boundary conditions ------
		nu_ref(el,0)=0;
		w_ref(el,0)=0;
		h_ref(el,0)=0;
		a_ref(el,0)=0;
		hnr_ref(el,0)=0;
		nu_ref(el,Nmax-1)=10000.; // Value at which we should not expect to detect pulsations anytime soon...
		w_ref(el,Nmax-1)=50;
		h_ref(el,Nmax-1)=0.0001;
		a_ref(el,Nmax-1)=0;
		hnr_ref(el,Nmax-1)=0;

		// ----------------------------------------
		//for(en=1; en<Nmax-1;en++){
		en=1;
		for(int ii=0; ii<pos.size();ii++){
			nu_ref(el,en)=data_ref.data(pos[ii],pref_freq);
			w_ref(el,en)=data_ref.data(pos[ii],pref_w);
			h_ref(el,en)=data_ref.data(pos[ii],pref_h);
			a_ref(el,en)=sqrt(PI*h_ref(el,en)*w_ref(el,en)); // mode amplitude... used for determining numax(el)
			hnr_ref(el,en)=h_ref(el,en)/data_ref.data(pos[ii],pref_n);
			en=en+1;
		}
	}

	//std::cout << "hnr_ref=" << std::endl;
	//std::cout << hnr_ref << std::endl;
	// Getting numax
	numax_ref=0;
	norm=0;
	for(el=0;el<lmax+1;el++){
		pos=where(a_ref.row(el), "=", a_ref.row(el).maxCoeff(), 0);
		numax_ref=numax_ref + nu_ref(el, pos[0]) * a_ref.row(el).maxCoeff(); // Weighted mean
		norm=norm + a_ref.row(el).maxCoeff();
		//std::cout << "numax(l="<< el << ")=" << nu_ref(el, pos[0]) <<std::endl;
	}
	numax_ref=numax_ref/norm;
	//std::cout << "mean numax=" << numax_ref << std::endl;

	// Getting \Delta\nu
	ones.resize(Nmax-2);
	ones.setOnes();

	std::cout << "      numax_ref=" << numax_ref << std::endl;
	//std::cout << "Chosen numax_ref=" << numax_ref << std::endl;

	//std::cout << nu_ref << std::endl; ///numax_ref- ones << std::endl;
	//exit(EXIT_SUCCESS);

	tmp=nu_ref.row(0).segment(1,Nmax-2);
	pos=where(tmp/numax_ref-ones, ">=", -0.25, 0);// AND nu_l0/numax_star-1 le 0.2);
	tmp2.resize(pos.size());
	for(int i=0; i<pos.size(); i++){tmp2[i]=tmp[pos[i]];}
	ones.setOnes(tmp2.size());
	pos=where(tmp2/numax_ref-ones, "<=", 0.25, 0);
	tmp.resize(pos.size());
	for(int i=0; i<pos.size(); i++){tmp[i]=tmp2[pos[i]];}
	//std::cout << "tmp=" << tmp << std::endl;
	if(tmp2.size() <= 3){
		std::cout << "Warning: small vector of frequencies detected for Mind=" << Mind << std::endl;
		std::cout << "         using the full dataset of frequencies for calculating Deltanu" << std::endl;
		//exit(EXIT_SUCCESS);
		tmp2.resize(nu_ref.col(0).size());
		tmp2=nu_ref.col(0);
	}

	tmp.setLinSpaced(tmp2.size(), 0, tmp2.size()-1);
	//std::cout << "fitting" << std::endl;
	fit=linfit(tmp, tmp2); // Use a linear fit of the modes around numax to get Delta\nu
	dnu_ref=fit[0];
	std::cout << "      dnu_ref=" << dnu_ref << std::endl;

	// --------- Variables size allocation ---------
	Nmax_star=10000.; // large value to begin with
	for(el=0; el<lmax+1; el++){
		//std::cout << "el=" << el << std::endl;
		pos=where(input_model.freqs.col(0), "=", el, 0);
		if (pos.size() < Nmax_star){ Nmax_star=pos.size();} // Keep the minimal number of radial order from the list of modes (ensure that Nmax_star is same for all l)
		//std::cout << "Nmax_star=" << Nmax_star << std::endl;
	}
	//std::cout << "Nmax_star=" << Nmax_star << std::endl;

	nu.resize(lmax+1, Nmax_star), h.resize(lmax+1, Nmax_star), w.resize(lmax+1, Nmax_star), s_a1.resize(lmax+1, Nmax_star), s_eta.resize(lmax+1, Nmax_star), 
		 s_a3.resize(lmax+1, Nmax_star), s_asym.resize(lmax+1, Nmax_star), s_b.resize(lmax+1, Nmax_star), s_alfa.resize(lmax+1, Nmax_star), inclination.resize(lmax+1, Nmax_star), 
		 mode_params.resize((lmax+1)*Nmax_star, 11);
	// Initialize variables for second order effects to neutral values (e.g. 0)
	s_eta.setZero(); //Neglect the centrifugal effects
	s_a3.setZero(); //Neglect the first order effect of the latitudinal rotation
	s_asym.setZero(); //Neglect the mode asymmetry
	s_b.setZero(); // Do not consider any additional frequency-dependent star distorsion
	s_alfa.setOnes(); // Do not consider any additional frequency-dependent star distorsion
	// -----------------------------

	// ---------- Setting up the variables for the simulated star ----------
	for(el=0; el<lmax+1; el++){
		//std::cout << "el=" << el << std::endl;
		// ------- Organise frequencies in a (l,n) table ---------
		pos=where(input_model.freqs.col(0), "=", el, 0);
		//for(int i=0;i<pos.size();i++){
		for(int i=0;i<Nmax_star;i++){
			//std::cout << "[" << i << "/" << Nmax_star-1 << "] " << "l=" << input_model.freqs(pos[i],0) << "pos=" << pos[i]  << " input_model.freqs(2,pos[i])= " << input_model.freqs(pos[i],2) << std::endl;
			nu(el,i)=input_model.freqs(pos[i],2);
		}
		// -------------------------------------------------------

		// ----- Perform rescalings -------
		x_ref.resize(Nmax);
		for(int en=0; en<Nmax;en++){
			x_ref[en]=(nu_ref(el,en)-numax_ref)/dnu_ref;
			//std::cout << x_ref[en] << std::endl;
		}

		//std::cout << "l=" << el << " Before interpols" << std::endl;
		x.resize(Nmax_star);
		for(int en=0; en<Nmax_star;en++){
			x[en]=(nu(el,en)-numax_star)/dnu_star;
			//std::cout << x[en] << std::endl;
		}
		/*
		std::cout << " x_ref      w_ref " << std::endl;
		for(int en=0; en<Nmax; en++){		
			std::cout << en << "  " << x_ref[en] << "  " <<  w_ref(el,en) << std::endl;
		}
		*/
		gamma.resize(Nmax_star);
		//std::cout << " x      gamma " << std::endl;
		for(int en=0; en<Nmax_star; en++){
			gamma[en]=lin_interpol(x_ref, w_ref.row(el), x[en]);
			//std::cout << en << "  " <<  x[en] << "  " << gamma[en] << std::endl;
		}
		/*
		std::cout << " x_ref      h_ref " << std::endl;
		for(int en=0; en<Nmax; en++){		
			std::cout << en << "  " << x_ref[en] << "  " << h_ref(el,en) << std::endl;
		}
		*/
		hnr.resize(Nmax_star);
		//std::cout << " x      h " << std::endl;
		for(int en=0; en<Nmax_star; en++){
			hnr[en]=lin_interpol(x_ref, hnr_ref.row(el), x[en]);
			//std::cout << en << "  " <<  x[en] << "  " << hnr[en] << std::endl;
		}
		if(el==0){ 
			c=maxHNR/hnr.maxCoeff(); // we will use l=0 to compute the new maxHNR
			//std::cout << "c=" << c << std::endl;
		}
		noise_star=harvey_1985(input_noise, nu.row(el)); // the local noise level of the simulated star

		// We solve: max(HNR_star) = c max(HNR_ref) ==> c=max(HNR_star)/max(HNR_ref)
		// We then use c to rescale the heights hstar(nu) noting that: c* href(nu)/Nref(nu) = hstar(nu)/Nstar(nu)
		height.resize(Nmax_star);
		for (int en=0; en<Nmax_star; en++){
			height[en]=hnr[en] * c * noise_star[en];
			//std::cout << height[en]/noise_star[en]/maxHNR << std::endl;
		}
		
		pos=where(height, "=", height.maxCoeff(),0);
		gamma=gamma * Gamma_Hmax_star/gamma(pos[0]); // Rescaling of the Widths in the y-axis. Gamma_star(pos[0]) is the width at Hmax
		// --------------------------

		// --------------- Other parameters setup -----------------
		h.row(el)=height;
		w.row(el)=gamma;
		s_a1.row(el).setConstant(a1);
		inclination.row(el).setConstant(inc);
	}

	// --------------------------------------------------------------------

	// ---------- Summarizing the information into suitable inputs for the writting function -------------
	for(el=0; el<lmax+1; el++){		
		for(k=el*Nmax_star; k<(el+1)*Nmax_star; k++){
			mode_params(k,0)=el;
			mode_params(k,1)=nu(el , k-el*Nmax_star);
			mode_params(k,2)=h(el , k-el*Nmax_star);
			mode_params(k,3)=w(el , k-el*Nmax_star);
			mode_params(k,4)=s_a1(el , k-el*Nmax_star);
			mode_params(k,5)=s_eta(el , k-el*Nmax_star);
			mode_params(k,6)=s_a3(el , k-el*Nmax_star);
			mode_params(k,7)=s_b(el , k-el*Nmax_star);
			mode_params(k,8)=s_alfa(el , k-el*Nmax_star);
			mode_params(k,9)=s_asym(el , k-el*Nmax_star);
			mode_params(k,10)=inclination(el , k-el*Nmax_star);
		}
	}

	// A FUNCTION THAT WRITES THE PARAMETERS
	write_star_mode_params_act_asym(mode_params, file_out_modes);

   	//std::cout << "           - Noise..." << std::endl;

	// A FUNCTION THAT WRITES THE Noise
	write_star_noise_params(input_noise, file_out_noise);

    	//std::cout << "           - Exit" << std::endl;
}


 
