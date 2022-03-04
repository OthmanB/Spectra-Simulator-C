/*
 * combi.cpp
 *
 *  Function that handle combination generation
 *  and read/write into files
 * 
 *  Created on: 10 Oct 2017
 *      Author: obenomar
 */

# include <iostream>
# include <iomanip>
#include <fstream>
# include <Eigen/Dense>
# include <vector>
#include <string>
# include "data.h"
# include "ioproc.h"
# include "format.h"

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;


/*
* Function that generate all possible combinations provided that the number
* of variables does not exceed 10. This function could certainly be better written
* by the use of combinatory function... but I could not find the proper algorithm
* so that I went for the simplest solution
*/
MatrixXd define_all_combinations(const MatrixXd Values, const VectorXi Nvalues, const int Nparams){

	if(Nparams > 10){
		std::cout << "This function calculates combinations up to maximum Nparams=10" << std::endl;
		std::cout << "Your value of Nparams is: ' + strtrim(Nparams,2) + ' which is greater than 10" << std::endl;
		std::cout << "Update the program to handle larger number of parameters. The program will exit now" << std::endl;
		exit(EXIT_SUCCESS);
	}
	long Ncombi=Nvalues.prod();
	MatrixXd allcombi(Ncombi, Nparams);
	long z=0;	

	std::cout << "     Total number of combination to be processed:" << Ncombi << std::endl;
	switch(Nparams){
		case 1: 
		   //std::cout << "Values=" << Values << std::endl;
		   allcombi.col(0)=Values;
		   //std::cout << "out" << std::endl;
		   break;
		case 2:
		   for(long i=0; i<Nvalues[0];i++){
		   	for(long j=0; j<Nvalues[1]; j++){
		 		//std::cout << "(" << i << "," << j <<")  <===>" << z << std::endl;
				allcombi.row(z) << Values(i,0), Values(j,1);
				//std::cout << Values(i,0) << "  " << Values(j,1) << std::endl;
				z=z+1.;
			}
		    }
		    break;
		case 3:
		   for(long i=0; i<Nvalues[0];i++){
			for(long j=0; j<Nvalues[1]; j++){
				for(long k=0; k<Nvalues[2]; k++){
					allcombi.row(z) << Values(i,0), Values(j,1), Values(k,2);
					z=z+1.;
				}
			}
		   }
		   break;
		case 4:		
		   for(long i=0; i<Nvalues[0];i++){
			for(long j=0; j<Nvalues[1]; j++){
				for(long k=0; k<Nvalues[2]; k++){
					for(long l=0; l<Nvalues[3]; l++){
						allcombi.row(z) << Values(i,0), Values(j,1), Values(k,2),Values(l,3);
						z=z+1.;
					}
				}
			}
		   }
		   break;
		case 5:
		   for(long i=0; i<Nvalues[0];i++){
			for(long j=0; j<Nvalues[1]; j++){
				for(long k=0; k<Nvalues[2]; k++){
					for(long l=0; l<Nvalues[3]; l++){
						for(long m=0; m<Nvalues[4]; m++){
							allcombi.row(z) << Values(i,0), Values(j,1), Values(k,2),Values(l,3), Values(m,4);
							z=z+1.;
						}
					}
				}
			}
		   }
		   break;
		case 6:
		   for(long i=0; i<Nvalues[0];i++){
			for(long j=0; j<Nvalues[1]; j++){
				for(long k=0; k<Nvalues[2]; k++){
					for(long l=0; l<Nvalues[3]; l++){
						for(long m=0; m<Nvalues[4]; m++){
							for(long n=0; n<Nvalues[5]; n++){
								allcombi.row(z) << Values(i,0), Values(j,1), Values(k,2),Values(l,3), Values(m,4), Values(n,5);
								z=z+1.;
							}
						}
					}
				}
			}
		   }
		   break;
		case 7:
		   for(long i=0; i<Nvalues[0];i++){
			for(long j=0; j<Nvalues[1]; j++){
				for(long k=0; k<Nvalues[2]; k++){
					for(long l=0; l<Nvalues[3]; l++){
						for(long m=0; m<Nvalues[4]; m++){
							for(long n=0; n<Nvalues[5]; n++){
								for(long o=0; o<Nvalues[6]; o++){
									allcombi.row(z) << Values(i,0), Values(j,1), Values(k,2),Values(l,3), Values(m,4), Values(n,5), Values(o,6);
									z=z+1.;
								}
							}
						}
					}
				}
			}
		   }
		   break;
		case 8:
		   for(long i=0; i<Nvalues[0];i++){
			for(long j=0; j<Nvalues[1]; j++){
				for(long k=0; k<Nvalues[2]; k++){
					for(long l=0; l<Nvalues[3]; l++){
						for(long m=0; m<Nvalues[4]; m++){
							for(long n=0; n<Nvalues[5]; n++){
								for(long o=0; o<Nvalues[6]; o++){
									for(long p=0; p<Nvalues[7]; p++){
										allcombi.row(z) << Values(i,0), Values(j,1), Values(k,2),Values(l,3), 
													   Values(m,4), Values(n,5), Values(o,6), Values(p,7);
										z=z+1.;
									}
								}
							}
						}
					}
				}
			}
		   }
		   break;
		case 9:
		   for(long i=0; i<Nvalues[0];i++){
			for(long j=0; j<Nvalues[1]; j++){
				for(long k=0; k<Nvalues[2]; k++){
					for(long l=0; l<Nvalues[3]; l++){
						for(long m=0; m<Nvalues[4]; m++){
							for(long n=0; n<Nvalues[5]; n++){
								for(long o=0; o<Nvalues[6]; o++){
									for(long p=0; p<Nvalues[7]; p++){
										for(long q=0; q<Nvalues[8]; q++){
											allcombi.row(z) << Values(i,0), Values(j,1), Values(k,2),Values(l,3), 
											 			   Values(m,4), Values(n,5), Values(o,6), Values(p,7), 
											 			   Values(q,7);
											z=z+1.;
										}
									}
								}
							}
						}
					}
				}
			}
		   }
		   break;
		case 10:
		   for(long i=0; i<Nvalues[0];i++){
			for(long j=0; j<Nvalues[1]; j++){
				for(long k=0; k<Nvalues[2]; k++){
					for(long l=0; l<Nvalues[3]; l++){
						for(long m=0; m<Nvalues[4]; m++){
							for(long n=0; n<Nvalues[5]; n++){
								for(long o=0; o<Nvalues[6]; o++){
									for(long p=0; p<Nvalues[7]; p++){
										for(long q=0; q<Nvalues[8]; q++){
											for(long r=0; r<Nvalues[9]; r++){
												allcombi.row(z) << Values(i,0), Values(j,1), Values(k,2),Values(l,3), 
											 				   Values(m,4), Values(n,5), Values(o,6), Values(p,7), 
											 				   Values(q,7), Values(r,8);
												z=z+1.;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		   }
		   break;
	}
return allcombi;
}



long read_id_allcombi(std::string file_combi){

	std::string lastline;
	std::vector<std::string> vals_last;

	lastline=read_lastline_ascii(file_combi);
	vals_last=strsplit(strtrim(lastline), " ");

	std::cout << "lastline=" << lastline << std::endl;
	for(int i=0; i<vals_last.size(); i++){
		std::cout << "vals_last[" << i << "]=" << vals_last[i] << std::endl;	
	}
	return str_to_lng(vals_last[0]);
}

std::string write_allcombi(MatrixXd allcombi, VectorXd cte_params, Config_Data cfg, std::string fileout, bool erase_old_file, long iter, long id0, 
		    std::vector<std::string> cte_names, std::vector<std::string> var_names, std::vector<std::string> param_names){
	
	int Nchars, precision;
	std::string id_str;
	VectorXd input_params, var_params;
	std::ofstream outfile;

	Nchars = 17;
	precision = 5;

	if(erase_old_file == 1 && iter == 0) {
		outfile.open(fileout.c_str()); // write a new file
	} else{
		outfile.open(fileout.c_str(), std::ios::app); // append
	}
	if(outfile.is_open()){
		//std::cout << "File opened" << std::endl;
		if(erase_old_file == 1 && iter == 0) { // Write Header only if we do not erase the old file AND this is the first execution of the function
			outfile << "model_name= " << cfg.model_name << std::endl;
			outfile << " --------------------------" << std::endl;
			outfile << "  List of all combinations " << std::endl;
			outfile << " --------------------------" << std::endl;
			outfile << "#" << std::setw(7) << "id  ";
			for(int s=0; s<param_names.size(); s++){
				outfile << std::setw(Nchars) << param_names[s];
			}
			outfile << std::endl;
		} 

		for(int i=0; i<allcombi.rows(); i++){
			id_str=identifier2chain(i + id0); // The identifier corresponds to the index of the current process + the initial id0
			outfile << std::setw(7) << id_str;
			//std::cout << "id_str=" << id_str << std::endl;
			
			var_params=allcombi.row(i).transpose();
			input_params=order_input_params(cte_params, var_params, cte_names, var_names, param_names);
			//std::cout << "input_params=" << input_params << std::endl;
			for(int j=0; j<input_params.size(); j++){			
				outfile << std::setw(Nchars) << std::setprecision(precision) << input_params(j);	
			}
			outfile << std::endl;
		}
	outfile.close();
	}  
	else {
		std::cout << " Unable to open file " << fileout << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	return id_str;
}


