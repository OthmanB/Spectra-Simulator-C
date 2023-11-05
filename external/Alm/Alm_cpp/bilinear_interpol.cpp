#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <gsl/gsl_interp2d.h>
#include "data.h"
#include "string_handler.h"

using namespace std;

// Load grid data from a gzip file (compressed using gzip) containing a single txt file with data as follow:
// Line 1: The number of rows in the grid followed by the number of columns in the grid
// Line 2: The x-coordinates for each row (in ascending order)
// Column 1: The y-coordinates for each column (in ascending order)
// Lines [3,Nx] - Columns [2,Ny] : The z-values for each point in the grid, listed in row-major order (i.e., the first n_cols values are the z-values for the first row, the next n_cols values are the z-values for the second row, and so on)
GridData loadGridData(const string& filename) {
    string line;
    string s;
    vector<std::string> v, sp;
    std::vector<double> ztmp;
    GridData data;
    // Define input stream
    ifstream file_in(filename, ios_base::in | ios_base::binary);
    if (!file_in) {
        throw runtime_error("Failed to open file: " + filename);
    }

    // Define filter chain for gzip-compressed file
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::gzip_decompressor()); // for gzip-compressed file
    in.push(file_in);
    while (getline(in, line)) {
        //std::cout << line << std::endl;
        s=strtrim(line); // remove white spaces at the end and begining
        v=strsplit(s, "="); // check the variable
        sp=strsplit(s, " \t");
        if (v[0] == "x"){
            // Read the first line which is the x-axis
            data.n_cols=sp.size();
            for(size_t i=0; i<data.n_cols; i++){
                data.x.push_back(str_to_dbl(sp[i]));
            }
        }
        if (v[0] == "y"){
            // Read the first line which is the x-axis
            data.n_rows=sp.size();
            for(size_t i=0; i<data.n_rows; i++){
                data.y.push_back(str_to_dbl(sp[i]));
            }
        } // x,y are done at that stage. We proceed with z. These are named in the data file just for the user
        //if (v[0] == "z" || (xdone == true && ydone == true)){
        if (v[0] == "z"){
            continue;
        } else{
            ztmp.clear();
            if(v[0] != "x" && v[0] != "y"){
                // The other lines are lines of z-values
                for(size_t i=0;i<data.n_cols;i++){
                    ztmp.push_back(str_to_dbl(sp[i]));
                }
                data.z.push_back(ztmp); // push the line into the 2D vector
                ztmp.clear();
            }
        }
    }  

    // Process the last line, if data are found on it. This allows to account for case where the user did not add a \n at the last line of the data file
	if (getline(in,line)){
        s=strtrim(line);
        sp=strsplit(s, " \t");
        // The other lines are lines of z-values
        //ztmp.clear();
        ztmp.clear();
        for(size_t i=0;i<data.n_cols;i++){
            ztmp.push_back(str_to_dbl(sp[i]));
        }
        data.z.push_back(ztmp); // push the line into the 2D vector
        ztmp.clear();
	}
    return data;
}
 
GridData4gsl flatten_grid(const GridData& data){
    GridData4gsl newGrid;
    newGrid.nx= data.x.size();
    newGrid.ny = data.y.size();
    newGrid.x= new double[newGrid.nx];
    newGrid.y= new double[newGrid.ny];
    newGrid.z = new double[newGrid.nx * newGrid.ny];
    for (int i = 0; i < newGrid.nx; i++) {
        newGrid.x[i] = data.x[i];
    }
    for (int j = 0; j < newGrid.ny; j++) {
        newGrid.y[j] = data.y[j];
    }
    for (int j = 0; j < newGrid.ny; j++) {
        for (int i = 0; i < newGrid.nx; i++) {
            newGrid.z[j * newGrid.nx + i] = data.z[j][i];
        }
    }
    return newGrid; 
}

gsl_interp2d* init_2dgrid(const GridData4gsl& data_flatten){
    // Initialize interp2d object
    //const gsl_interp2d_type* T = gsl_interp2d_bilinear;//gsl_interp2d_bicubic;
    const gsl_interp2d_type* T = gsl_interp2d_bicubic;
    gsl_interp2d* interp = gsl_interp2d_alloc(T, data_flatten.nx, data_flatten.ny);
    gsl_interp2d_init(interp, data_flatten.x, data_flatten.y, data_flatten.z, data_flatten.nx, data_flatten.ny);
    return interp;
}

// Interpolation for most cases (no within large loop calls)
double interpolate_core(gsl_interp2d* interp, const GridData4gsl& data_flatten, double x, double y){
    double result;
    gsl_interp_accel* xacc = gsl_interp_accel_alloc();
    gsl_interp_accel* yacc = gsl_interp_accel_alloc();
    
    gsl_interp2d_eval_e(interp, data_flatten.x, data_flatten.y, data_flatten.z, x, y, xacc, yacc, &result);    
    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);
    return result;
}

// Interpolation for most cases (no within large loop calls)
double interpolate(const GridData& data, double x, double y){
    GridData4gsl data_flatten=flatten_grid(data);
    gsl_interp2d* interp=init_2dgrid(data_flatten);
    double r=interpolate_core(interp, data_flatten, x, y);
    // Free resources
    delete[] data_flatten.x;
    delete[] data_flatten.y;
    delete[] data_flatten.z;
    gsl_interp2d_free(interp);
    // Evaluate interpolated value
    return r;
}
