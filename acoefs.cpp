/*
Date: 16 Nov 2021
Functions that handle the acoefficients
It can create nu(n,l,m) from acoeffs for l<=3, j={1,2,3,4,5,6}
Or it can decompose into aj coefficients
Adapted from acoefs.py from the acoefs_check project (see github)
*/
#include <math.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

long double Hslm_Ritzoller1991(const int s,const int l, const int m){
	const int L=l*(l+1);
	long double Hsm;

	if (s == 0){
		Hsm=1;
	}
	if(s == 1){
		Hsm=2*m;
	}
	if(s == 2){
		Hsm=6*pow(m,2) - 2*L;
	}
	if (s == 3){
		Hsm=20*pow(m,3) - 4*(3*L-1)*m;
	}
	if (s == 4){
		Hsm=70*pow(m,4) - 10*(6*L-5)*pow(m,2) + 6*L*(L-2);
	}
	if (s == 5){
		Hsm=252*pow(m,5)-140*(2*L-3)*pow(m,3) + (20*L*(3*L-10) + 48)*m;
	}
	if (s == 6){
		Hsm=924*pow(m,6) - 420*pow(m,4) * (3*L-7) + 84*pow(m,2) *(5*pow(L,2) - 25*L + 14) - 20*L*(pow(L,2) - 8*L + 12);
	}
	if (s > 6){
		std::cout << "Warning in Hslm_Ritzoller1991: s>6 not supported" << std::endl;
		return -1;
	}
	return Hsm;
}

long double Pslm(const int s,const int l,const int m){
	// These take the Ritzoller1991 coefficients and normalise them by Pslm(l)=l
	// As per specified in Schou, JCD, Thompson 1994
	// so basically we solved Pslm(l)=l=c*Hslm and find c
	long double H, c, Ps;

	if (s==0){
		Ps=l;
	}
	if (s==1){
		Ps=m;
	}
	if (s==2){
		if (l>0){
			Ps=(3*pow(m,2) -l*(l+1))/(2*l-1);
		} else{
			Ps=0;
		}
	}
	if (s==3){
		if (l>1){
			Ps=(5*pow(m,3) - (3*l*(l+1)-1)*m)/((l-1)*(2*l-1));
		} else{
			Ps=0;
		}
	}
	if (s==4){
		H=(35*pow(m,4) - 5*(6*l*(l+1)-5)*pow(m,2)) + 3*l*(l+1)*(l*(l+1)-2);
		c=2*(l-1)*(2*l-1)*(2*l-3);
		if (c !=0){
			Ps=H/c;
		} else{
			Ps=0;
		}
	}
	if (s==5){
		H=Hslm_Ritzoller1991(s,l,m);
		c=8*(4*pow(l,4) - 20*pow(l,3) + 35*pow(l,2) - 25*l + 6);
		if (c !=0){
			Ps=H/c;
		} else{
			Ps=0;
		}
	}
	if (s==6){
		H=Hslm_Ritzoller1991(s,l,m);
		c=64*pow(l,5) - 480*pow(l,4) + 1360*pow(l,3) - 1800*pow(l,2) + 1096*l - 240;
		if (c !=0){
			Ps=H/c;
		} else{
			Ps=0;
		}
	}
	if (s>6){
		std::cout << "Warning in Pslm(): s>6 not supported" << std::endl;
		std::cout << "The program will return 0" << std::endl;
		return 0;
	}
	return Ps;
}


VectorXd Tnlm(VectorXd& nu_nlm, const int l){
	VectorXd tnlm;
	if (l != 0 && l<=3){
		tnlm.resize(l);
	}
	switch (l){
		case 1:
			//Tn11
			tnlm[0]=(nu_nlm[2]- nu_nlm[0])/2;
			break;
		case 2:
			//Tn21
			tnlm[0]=(nu_nlm[3]-nu_nlm[1])/2;
			//Tn22
			tnlm[1]=(nu_nlm[4]-nu_nlm[0])/4;
			break;
		case 3:
			//Tn31
			tnlm[0]=(nu_nlm[4]-nu_nlm[2])/2;
			//Tn32
			tnlm[1]=(nu_nlm[5]-nu_nlm[1])/4;
			//Tn33
			tnlm[2]=(nu_nlm[6]-nu_nlm[0])/6;
			break;
		default:
			tnlm.resize(1);
			tnlm[0]=0;
			std::cout << " Warning in Tnlm: Please enter a value of 0<l<4" << std::endl;
			std::cout << "     The function will return [0]" << std::endl;
	}
	return tnlm;	
}

VectorXd Snlm(VectorXd& nu_nlm, const int l){
	VectorXd snlm;
	if (l != 0 && l<=3){
		snlm.resize(l);
	}
	switch (l){
		case 1:
			//Sn11
			snlm[0]=(nu_nlm[0] + nu_nlm[2])/2 - nu_nlm[1];
			break;
		case 2:
			//Sn21
			snlm[0]=(nu_nlm[1] + nu_nlm[3])/2 - nu_nlm[2];
			//Sn22
			snlm[1]=(nu_nlm[0] + nu_nlm[4])/2 - nu_nlm[2];
			break;
		case 3:
			//Sn31
			snlm[0]=(nu_nlm[2] + nu_nlm[4])/2 - nu_nlm[3];
			//Sn32
			snlm[1]=(nu_nlm[1] + nu_nlm[5])/2 - nu_nlm[3];
			//Sn33
			snlm[2]=(nu_nlm[0] + nu_nlm[6])/2 - nu_nlm[3];
			break;
		default:
			snlm.resize(1);
			snlm[0]=0;
			std::cout << " Please enter a value of 0<l<4" << std::endl;
			std::cout << "     The function will return [0]" << std::endl;
	}
	return snlm;
}

VectorXd nunlm_from_acoefs(const long double nunl0, const int l, 
	const long double a1, const long double a2, const long double a3, const long double a4, const long double a5, const long double a6){
	// This function compute nu_nlm from a series of a-coeficient and provided the central frequency without splitting nunl0
	VectorXd nu_nlm(2*l+1);
	for (int m=-l; m<=l; m++){
		nu_nlm[m+l]=nunl0 + a1*Pslm(1,l,m) + a2*Pslm(2,l,m) + a3*Pslm(3,l,m) + a4*Pslm(4,l,m)+ a5*Pslm(5,l,m) + a6*Pslm(6,l,m);
	}
	return nu_nlm;
}


VectorXd eval_acoefs(const int l, VectorXd& nu_nls){ // We expect nu_nls=[nu(-l), nu(-l+1), ... , nu[0], nu[1], ... nu(l)]
	// Function that gets the splitted frequencies of a given mode (n,l) and determines the analytical a-coefficients
	// from a1 to a6 and for l<=3
	// More details on the equations that lead to this are on Benomar+2021 and in acoefs.py on github acoefs_check project
	long double Num_a1, Den_a1;
	const long double A=-14;//, B=0, D=0;// A0=-14, A1=-3, B0=1,B1=-2/3, C0=-15,C1=-4, D0=18,D1=14/3;
	long double C, E, F, G;
	long double Pi21, Pi22, Pi23, Pi41, Pi42, Pi43, Pi61, Pi62, Pi63;

	VectorXd aj(6), tnlm, snlm;
	aj.setZero();

	if (l ==0){
		std::cout <<"There is no a1 coefficient for l==0" << std::endl;
		std::cout << "The function will return 0 for all coefficients " << std::endl;
		return aj;
	}
	if (l>3){
		std::cout << "an for l>3 not implemented. Should you need it, better to solve this algorithmically using equation A3-6 from Schou, JCD, Thompson, 1994" << std::endl;
		std::cout << "The program will return 0 for all coefficients" << std::endl;
		return aj;
	}
	snlm=Snlm(nu_nls, l);
	tnlm=Tnlm(nu_nls, l);
	switch (l){
		case 1:
			aj[0]=tnlm[0]; //(nu_nls[2]- nu_nls[0])/2;
			aj[1]=snlm[0]/3; //((nu_nls[0] + nu_nls[2])/2 - nu_nls[1])/3;
			break;
		case 2:
			// a1 Term:
			Num_a1=tnlm[0]+4*tnlm[1]; //T21-2*T22*Pslm(3,2,1)/Pslm(3,2,2);
			Den_a1=5; //Pslm(1,2,1) - Pslm(3,2,1)*Pslm(1,2,2)/Pslm(3,2,2);
			aj[0]=Num_a1/Den_a1;
			// a2 term:
			aj[1]=(2*snlm[1] - snlm[0])/7;
			// a3 term:
			aj[2]=(tnlm[1] - tnlm[0])/5;
			// a4 term:
			aj[3]=(snlm[1] - 4*snlm[0])/70.;
			break;
		case 3:
			// Odds terms
			// We have Tnlm = Sum [a_{2j-1} P_{2j-1}] with j=[1,M/2]
			// a1 term:
			aj[0]=tnlm[0]/14 + 2*tnlm[1]/7 + 9*tnlm[2]/14;
			// a3 term: 
			aj[2]=-tnlm[0]/9 - 2*tnlm[1]/9 + tnlm[2]/3;
			// a5 term: 
			aj[4] = tnlm[2]/42 + 5*tnlm[0]/126 - 4*tnlm[1]/63;
			// Even Terms
			// We have Snlm = Sum[a2j (P2j - P2j(0))] with j=[1, M/2]
			// a2 term:
			aj[1]=(-15*snlm[0] + 25*snlm[2])/126;
			// a4 term:
			aj[3]=13*(snlm[0] - 7*snlm[1] + 3*snlm[2])/1001;
			// a6 term:
			aj[5]=(15*snlm[0] - 6*snlm[1] + snlm[2])/1386;
			break;
		default:
			std::cout << "Error eval_acoefs(): l must be between 1 and 3. Your input is l=" << l << std::endl;
			std::cout << "The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
			break;
		}
	return aj;
}
