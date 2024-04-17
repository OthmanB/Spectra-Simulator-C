/**
 * @file plots_diags.cpp
 * @brief Plot the model using Gnuplot.
 *
 * Handling plots using gnuplot
 *
 * @date 06 Oct 2016
 * @author obenomar
 */

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "plots_diags.h"
#include "io_star_params.h"
#include "ioproc.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;


void gnuplt_model(VectorXd x, VectorXd y, VectorXd model, double scoef1, double scoef2, std::string file_model){
/* 
 * Simple function that handle the plots of the model.
 * Here the data are read in a temporary file that is created here. It will contain:
 *    - Header with the smoothing coefficient scoef1 and scoef2
 *    - col1 the frequency
 *    - col2 the noisy spectrum
 *    - col3 a smooth of the noisy spectrum at a level scoef1
 *    - col4 a smooth of the noisy spectrum at a level scoef2
 *    - col5 the asymptotic spectrum
 * If ps is set to 1, then the plot is written in file_out. Otherwise, it is shown
 * on the screen.
*/

    bool ps;
    std::string filename_txt_spectrum;
    std::stringstream ss;
    std::string xlabel, ylabel, xunit, yunit;
    
 	ps=1;
	xlabel="Frequency";
	xunit="(uHz)";
	ylabel="Power";
	yunit="(ppm^2/microHz)";

	Gnuplot gn;
	
	// Create a temporary file with the data that we wish to plot
	filename_txt_spectrum="graph_tmp.txt";
    write_spectrum(x, y, model, scoef1, scoef2, filename_txt_spectrum); // This function is inside the file 'write_star_params.cpp'

    	if (ps == 0) {
		//gn << "set term X \n";
    	gn << "set term xterm \n";
    	} else {
		gn << "set term post eps enhanced color font 'Times-Bold, 15'\n";
		gn << "set out '" + file_model + "'\n"; //+ ".eps'\n";
   	 }
   	 // Setup the common plot configuration
 	   gn << "set datafile commentschars '#!'\n";
 	   gn << "set autoscale \n"; // scale axis automatically
 	   gn << "unset label \n";
 	   gn << "set key default \n"; // restore the default position for the legends
 	   gn << "set key top right \n"; // legend will be on bottom right

  	  // The plot for the likelihood of all parallel chains
       	gn << "set xlabel '" + xlabel + " " + xunit + "'" << std::endl;
       	gn << "set ylabel '" + ylabel + " " + yunit + "'" << std::endl;
	gn << "plot '" + filename_txt_spectrum +  "' using 1:2 with lines title 'scoef="+dbl_to_str(scoef1) +"'" + 
		    ",'" + filename_txt_spectrum +  "' using 1:3 with lines title 'scoef="+dbl_to_str(scoef2) +"'" +
		    ",'" + filename_txt_spectrum +  "' using 1:5 with lines title 'Model'";
  	gn << std::endl;
		
    gn << "set term xterm" << std::endl;


}
