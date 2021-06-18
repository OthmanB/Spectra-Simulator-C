/*
 * format.cpp
 *
 *  Functions that format long or str into
 *  a suitable format for read/write
 * 
 *  Created on: 10 Oct 2017
 *      Author: obenomar
 */

//#include <math.h>
//#include <Eigen/Dense>
# include <iostream>
# include <iomanip>
# include "format.h"
# include "ioproc.h" // contains the string handlers

//# include <vector>

using Eigen::VectorXd;
using Eigen::MatrixXd;


std::string identifier2chain(long identifier){

	std::string out;

	if  (identifier < 10) { out="000000" + strtrim(lng_to_str(identifier)); }
	if ((identifier >= 10) && (identifier < 100)) { out="00000" + strtrim(lng_to_str(identifier));}
	if ((identifier >= 100) && (identifier < 1000)) { out="0000" + strtrim(lng_to_str(identifier));}
	if ((identifier >= 1000) && (identifier < 10000)) { out="000" + strtrim(lng_to_str(identifier));}
	if ((identifier >= 10000) && (identifier  < 100000)) { out= "00" + strtrim(lng_to_str(identifier)); }
	if ((identifier >= 100000) && (identifier< 1000000)) { out="0" + strtrim(lng_to_str(identifier));}
	if (identifier  >  10000000) { 
		std::cout << "Warning: This cannot handle greater number than 99999" << std::endl;
		std::cout << "Pursuing will lead to an improper naming of the models" << std::endl;
		std::cout << "Please update identifier2chain in order to handle greater numbers" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	}
	return out;
}

std::string identifier2chain(std::string identifier){

	std::string out;

	if  (str_to_dbl(identifier) < 10) { out="000000" + strtrim(identifier); }
	if ((str_to_dbl(identifier) >= 10) && (str_to_dbl(identifier)    < 100)) { out="00000" + strtrim(identifier);}
	if ((str_to_dbl(identifier) >= 100) && (str_to_dbl(identifier)   < 1000)) { out="0000" + strtrim(identifier);}
	if ((str_to_dbl(identifier) >= 1000) && (str_to_dbl(identifier)  < 10000)) { out="000" + strtrim(identifier);}
	if ((str_to_dbl(identifier)  >= 10000) && (str_to_dbl(identifier)  < 100000)) { out= "00" + strtrim(identifier); }
	if ((str_to_dbl(identifier) >= 100000) && (str_to_dbl(identifier)< 1000000)) { out="0" + strtrim(identifier);}
	if (str_to_dbl(identifier)  >  10000000) { 
		std::cout << "Warning: This cannot handle greater number than 9999999" << std::endl;
		std::cout << "Pursuing will lead to an improper naming of the models" << std::endl;
		std::cout << "Please update identifier2chain in order to handle greater numbers" << std::endl;
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	}
	return out;
}

/* 
 * Function that format properly the name of the frequency file
 * produced by ADIPLS
 * Used to read ADIPLS frequency outputs 
*/
std::string format_freqname(std::string id){

	std::string delimiters, out="";
	std::vector<std::string> split;
	long identifier;
	delimiters=".";
	split=strsplit(id, delimiters);
	identifier=str_to_lng(split[0]);

	if  (identifier < 10) { out="00000" + id;}
	if ((identifier >= 10) && (identifier < 100)) { out="0000" + id;}
	if ((identifier >= 100) && (identifier < 1000)) { out="000" + id;}
	if ((identifier >= 1000) && (identifier < 10000)) { out="00" + id;}
	if ((identifier >= 10000) && (identifier  < 100000)) { out= "0" + id;}
	if ((identifier >= 100000) && (identifier  < 1000000)) { out= id;}
	if (identifier  >  1000000) {
		std::cout << "Warning: This cannot handle greater number than 999999" << std::endl;
		std::cout << "Pursuing will lead to an improper naming of the models" << std::endl;
		std::cout << "Please update format_freqname if you think you need to handle greater numbers" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	//std::cout << split[0] << std::endl;
	//std::cout << identifier << std::endl;
	//std::cout << out << std::endl;
	//exit(EXIT_SUCCESS);

	return out;
}



