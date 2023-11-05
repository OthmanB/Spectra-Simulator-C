#include "gtest/gtest.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include "activity.h"
#include "data.h"
#include "make_grids.h"
#include "do_grids.h"
#include "gzip_compress.h"
#include "bilinear_interpol.h"
#include "Alm_interpol.h"

using namespace std;

int rem_grids(const std::string dir_to_rem){
    int rem;
    const std::string fileroot = dir_to_rem+"/A"; // dir_to_rem is eg 'grids'. 
    std::string filename;
    for (size_t l=1; l<=3; l++){
        for (size_t m=0; m<2*l+1;m++){
            filename=fileroot + std::to_string(l) + std::to_string(m) + ".gz";
            rem=remove(filename.c_str());
            if (rem !=0){
                std::cerr << "Error code " + std::to_string(rem) + "Failed to remove the file\n";
                return rem;
            } 
        }
    }
    rem=rmdir(dir_to_rem.c_str());
    if (rem !=0){
        std::cerr << "Error code " + std::to_string(rem) + "Failed to remove the directory\n";
        return rem;
    }
    return 1;
}

void show_matrix_Alm(GridData_Alm grid, int m){
        int l=(grid.Alm_grid.size()-1)/2;
        // Show the matrix for debug
        std::cout << std::setw(22);
        for (int i=0; i<grid.theta.size(); i++){
            std::cout << grid.theta[i] << std::setw(14);
        }
        std::cout << std::endl << std::setw(21) << "i=  ";
        for (int i=0; i<grid.theta.size(); i++){
            std::cout << i << std::setw(14);
        }
        std::cout << std::endl;
        for (int j=0; j<grid.delta.size(); j++){
            for (int i=0; i<grid.theta.size(); i++){
                if (i==0){
                    std::cout << grid.delta[j] << "  j=" << j << "  ";
                }
                std::cout << grid.Alm_grid[m+l][j][i] << std::setw(14);
            }
            std::cout << std::endl;
        }
}

void show_matrix(GridData grid){
        int l=(grid.z.size()-1)/2;
        // Show the matrix for debug
        std::cout << std::setw(22);
        for (int i=0; i<grid.x.size(); i++){
            std::cout << grid.x[i] << std::setw(14);
        }
        std::cout << std::endl << std::setw(21) << "i=  ";
        for (int i=0; i<grid.x.size(); i++){
            std::cout << i << std::setw(14);
        }
        std::cout << std::endl;
        for (int j=0; j<grid.y.size(); j++){
            for (int i=0; i<grid.x.size(); i++){
                if (i==0){
                    std::cout << grid.y[j] << "  j=" << j << "  ";
                }
                std::cout << grid.z[j][i] << std::setw(14);
            }
            std::cout << std::endl;
        }
}

double getmax(std::vector<double> x){
    double max=0;
    for(size_t i=0;i<x.size();i++){
        if (x[i] > max){
            max=x[i];
        }
    }
    return max;
}

double getmin(std::vector<double> x){
    double min=getmax(x);
    for(size_t i=0;i<x.size();i++){
        if (x[i] < min){
            min=x[i];
        }
    }
    return min;
}

int random_int(int min, int max) {
  std::random_device rd;  // obtain a random number from hardware
  std::mt19937 gen(rd()); // seed the generator
  
  // create a distribution
  std::uniform_int_distribution<> distr(min, max);
  
  // generate and return a random integer within the range [min, max]
  return distr(gen);
}

double random_double(double min, double max) {
  std::random_device rd;  // obtain a random number from hardware
  std::mt19937 gen(rd()); // seed the generator

  // create a distribution
  std::uniform_real_distribution<> distr(min, max);

  // generate and return a random double within the range [min, max]\n  
  return distr(gen);
}


// ---- Tests for gzip_compress ----
TEST(CompressionTests, TestCompressGzip) {
    // Test compress_gzip function by creating a test file,
    // compressing it and verifying that the compressed file exists
    std::string source = "testfile.txt";
    std::string destination = "testfile.txt.gz";

    // Create a test file with some content
    std::ofstream outfile(source);
    outfile << "Hello this is a test file" << std::endl;
    outfile.close();

    // Compress the test file
    compress_gzip(source, destination);

    // Verify that the compressed file exists
    std::ifstream infile(destination);
    bool exists = infile.good();
    infile.close();
    ASSERT_TRUE(exists);

    // Clean up
    std::remove(source.c_str());
    std::remove(destination.c_str());
}

TEST(FileIOTests, TestWriteToFile) {
    // Test writeToFile function by writing some data to a file
    // and verifying that the file contains the same data

    std::vector<double> x = {1, 2, 3};
    std::vector<double> y = {4, 5, 6};
    std::vector< std::vector<double> > z = {{7, 8, 9}, {10, 11, 12}};

    std::string filename = "testdata.txt";
    writeToFile(x, y, z, filename);

    std::ifstream infile(filename);
    std::string line;
    std::getline(infile, line);
    ASSERT_EQ(line, "x=1 2 3 ");
    std::getline(infile, line);
    ASSERT_EQ(line, "y=4 5 6 ");
    std::getline(infile, line);
    ASSERT_EQ(line, "z=");
    std::getline(infile, line);
    ASSERT_EQ(line, "7 8 9 ");
    std::getline(infile, line);
    ASSERT_EQ(line, "10 11 12 ");

    infile.close();
    // Clean up
    std::remove(filename.c_str());
}

// -----------

// Tests for make_grids.cpp
TEST(make_Alm_grid, validInputs) {
    const int l = 1;
    const double theta_min = 0.;
    const double theta_max = M_PI/2;
    const double delta_min = 0.0;
    const double delta_max = M_PI/4;
    const double resol = 0.1;
    const std::string ftype = "gate";

    GridData_Alm grid = make_Alm_grid(l, theta_min, theta_max, delta_min, delta_max, resol, ftype, false);

    const double max_theta = getmax(grid.theta);
    const double max_delta = getmax(grid.delta);

    // Check the basic structure
    EXPECT_EQ(grid.Ntheta, ceil((theta_max - theta_min) / resol));
    EXPECT_EQ(grid.Ndelta, ceil((delta_max - delta_min) / resol));
    EXPECT_EQ(max_theta, theta_max);
    EXPECT_EQ(max_delta, delta_max);
    EXPECT_EQ(grid.l, l);
    EXPECT_EQ(grid.Alm_grid[0].size(), grid.Ndelta);
    EXPECT_EQ(grid.Alm_grid[0][0].size(), grid.Ntheta);
    EXPECT_EQ(grid.Alm_grid.size(), 2*grid.l+1);

    // Check the inner values: Are they put in the right location in the 3D array?
    // In particular, we want to be sure that the theta and delta axis are not reversed and properly addressed
    // For this test we do some property-based test
    const int Niter=20; //We perform 20 random reading tests in the array
    int m,i0,j0;
    double Alm_norm, theta0, delta0, expected;
    std::cout << "iter=";
    for(size_t iter=0; iter<Niter; iter++){
        std::cout << iter+1 << "...";
        m=random_int(-l, l);
        i0=random_int(0, grid.Ntheta-1);
        j0=random_int(0, grid.Ndelta-1);
        // We compute the Alm with the Alm function before comparing it with the tabulated value
        expected=Alm(l, m, grid.theta[i0], grid.delta[j0], ftype); 
        if (ftype == "gauss"){
            Alm_norm= gauss_filter_cte(theta0, delta0);
        }
        else{
            Alm_norm=1;
        }
        expected=expected/Alm_norm;
       // Comparison
        EXPECT_DOUBLE_EQ(expected, grid.Alm_grid[m+l][j0][i0]);
    }
    std::cout << "End." << std::endl;

    // Cleanup
    int rem;
    std::string fileroot = "grids/A";
    std::string filename;
    for (size_t l=1; l<=3; l++){
        for (size_t m=0; m<2*l+1;m++){
            filename=fileroot + std::to_string(l) + std::to_string(m) + ".gz";
            rem=remove(filename.c_str());
            if (rem !=0){
                std::cerr << "Error code " + std::to_string(rem) + "Failed to remove the file\n";
            } 
        }
    }
    rem=rmdir("grids");
    if (rem !=0){
        std::cerr << "Error code " + std::to_string(rem) + "Failed to remove the directory\n";
    }

}

TEST(make_Alm_grid, invalidInputs) {
    const int l = 1;
    double theta_min = 0;
    double theta_max = M_PI + M_PI/10;
    double delta_min = 0.0;
    double delta_max = M_PI/4;
    double resol = 0.05;
    std::string ftype = "gate";
    EXPECT_EXIT(make_Alm_grid(l, theta_min, theta_max, delta_min, delta_max, resol, ftype, false), 
        ::testing::ExitedWithCode(EXIT_FAILURE), "Error: theta_min must be > 0 and theta_max must < Pi with also theta_max>theta_min.*");
    
    theta_min = 0.;
    theta_max = M_PI/2;
    delta_min = 0.0;
    delta_max = M_PI + M_PI/10;
    EXPECT_EXIT(make_Alm_grid(l, theta_min, theta_max, delta_min, delta_max, resol, ftype, false), 
        ::testing::ExitedWithCode(EXIT_FAILURE), "Error: delta_min must be > 0 and delta_max must < Pi with also delta_max>delta_min.*");

    theta_min = 0.1;
    theta_max = M_PI/4;
    delta_min = 0.05;
    delta_max = M_PI/3;
    resol=-0.1;
    EXPECT_EXIT(make_Alm_grid(l, theta_min, theta_max, delta_min, delta_max, resol, ftype, false), 
        ::testing::ExitedWithCode(EXIT_FAILURE), "Error: resol must be > 0.*");

    ftype = "invalid";
    EXPECT_EXIT(make_Alm_grid(l, theta_min, theta_max, delta_min, delta_max, resol, ftype, false), 
        ::testing::ExitedWithCode(EXIT_FAILURE), "Error: ftype must be either 'gate', 'triangle' or 'gauss'.*");
    
}

// Tests for do_grids.cpp
TEST(MakeGridsTest, ValidInputs) {
    // test with valid input values and filter type
    char* argv[] = {const_cast<char*>("prog_name"),
            const_cast<char*>("--resol=0.1"),
            const_cast<char*>("--lmax=3"), 
            const_cast<char*>("--ftype=gate"),  
            const_cast<char*>("--verbose=false")};
    int argc = sizeof(argv) / sizeof(argv[0]);

    EXPECT_EQ(do_grids(argc, argv), 0);
    rem_grids("grids/");
}

TEST(MakeGridsTest, InvalidInputs) {
    // test with invalid input values for resol and filter type
    char* argv[] = {const_cast<char*>("prog_name"), 
            const_cast<char*>("--resol=-0.1"), 
            const_cast<char*>("--lmax=3"), 
            const_cast<char*>("--ftype=invalid_type"),  
            const_cast<char*>("--verbose=false")};
    int argc = sizeof(argv) / sizeof(argv[0]);

    EXPECT_EQ(do_grids(argc, argv), 1);
}

TEST(MakeGridsTest, SaveGrid) {
    // test the ability to save the generated grid data
    char* argv[] = {const_cast<char*>("prog_name"), 
            const_cast<char*>("--outdir=grid_dir"), 
            const_cast<char*>("--resol=0.1"), 
            const_cast<char*>("--lmax=3"), 
            const_cast<char*>("--ftype=gate"),  
            const_cast<char*>("--verbose=false")};
    int argc = sizeof(argv) / sizeof(argv[0]);

    EXPECT_EQ(do_grids(argc, argv), 0);
    rem_grids("grid_dir/");
}

TEST(MakeGridsTest, UnexpectedError) {
    // test for unexpected error such as insufficient disk space
    char* argv[] = {const_cast<char*>("prog_name"), 
            const_cast<char*>("--outdir=/"), 
            const_cast<char*>("--resol=0.1"), 
            const_cast<char*>("--lmax=3"), 
            const_cast<char*>("--ftype=gate"),  
            const_cast<char*>("--verbose=false")};
    int argc = sizeof(argv) / sizeof(argv[0]);

    EXPECT_EQ(do_grids(argc, argv), 1);
}


// A function that generates the reference array of data
// This is the same content as the file grid_data_test.txt.gz but in plain text
GridData gen_refdata(){
    GridData data_ref;
    std::vector<double> x = { 0.0, 0.5, 1.0, 1.5, 2.};
    std::vector<double> y = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
    std::vector<std::vector<double>> z = { { 1.0, 2.0, 3.0, 4.0, 5.0 }, 
                                           { 2.0, 3.0, 4.0, 5.0, 6.0 }, 
                                           { 3.0, 4.0, 5.0, 6.0, 7.0},
                                           { 4.0, 5.0, 6.0, 7.0, 8.0},
                                           { 5.0, 6.0, 7.0, 8.0, 9.0},
                                           { 6.0, 7.0, 8.0, 9.0, 10.0},
                                        };
    data_ref.x=x;
    data_ref.y=y;
    data_ref.z=z;
    
    return data_ref;
}

// ------ Interpolation functions test -----

// Test that loadGridData correctly reads a file containing grid data
TEST(loadGridData, correctGridData) {
    // load grid data and check values
    GridData data = loadGridData("../../../data/Alm_grids_CPP/basic_tests/grid_data_test.txt.gz");
    // load the hardcoded reference data for the test
    GridData data_ref = gen_refdata();
    EXPECT_EQ(data.n_rows, data_ref.y.size());
    EXPECT_EQ(data.n_cols, data_ref.x.size());
    for(size_t i=0; i<data.x.size();i++){
        EXPECT_EQ(data.x[i], data_ref.x[i]);
    }
    for(size_t i=0; i<data.y.size();i++){
        EXPECT_EQ(data.y[i], data_ref.y[i]);
    }
    for(size_t j=0; j<data.y.size();j++){
        for(size_t i=0; i<data.x.size();i++){
            EXPECT_EQ(data.z[j][i], data_ref.z[j][i]);
        }
    }
    //show_matrix(data);
}

// Test that interpolate performs correctly the bilinear interpolation for a given set of input values
TEST(interpolate, correctInterpolation) {
    // load the hard coded reference grid data
    GridData data_ref=gen_refdata();

    // interpolate at the nodes of the grids: These should return the same values as the nodes
    for(size_t j=0; j<data_ref.y.size();j++){
        for(size_t i=0; i<data_ref.x.size();i++){
            double c = interpolate(data_ref, data_ref.x[i], data_ref.y[j]);
            EXPECT_NEAR(c, data_ref.z[j][i], 1e-6);
        }
    }
    //show_matrix(data_ref);
}
 
 
// -------- Test for the Alm interpolation ----
TEST(Alm_interpTest, testValidInputs_gate_case) {
    // directory with test grid files: achieves errors smaller than 3e-2 (test on 1000 samples)
   //const std::string grid_dir = "../../../data/Alm_grids_CPP/test_grids/"; 
    //
    // directory with test grid files: achieves errors guaranteed to be lower than 9e-3 (test on 1000 samples)
    const std::string grid_dir = "../../../data/Alm_grids_CPP/2deg_grids/"; 
    //
    // directory with test grid files: achieves 2e-3 errors (test on 500 samples)
    //const std::string grid_dir = "../../../data/Alm_grids_CPP/1deg_grids/"; 
    //
    const std::string ftype="gate";
    const int lmax=3;
    const int Niter=50; // total number of random tests
    // Test valid inputs and check if the function returns the expected output
    int l; // some input value
    int m; // some input value
    long double theta0; // some input value
    long double delta; // some input value
    long double interp_output;
    long double expected_output;
    // Property based testing to verify that the interpolation is precise enough on a grid
    for(int iter=0; iter<Niter;iter++){
        l=random_int(1,lmax);
        m=random_int(-l,l);
        theta0=random_double(0, M_PI/2);
        delta=random_double(0, M_PI/4);
        interp_output = Alm_interp(l, m, theta0, delta, ftype, grid_dir); // approximation function
        expected_output= Alm(l, m, theta0, delta, ftype); // full integral computation
        // Compare. The tolerance depends on the sparsity of the grid and of the interpolation algo
        std::cout << "["<<iter<<"]         l = " << l;
        std::cout << "         m = " << m;
        std::cout << "    theta0 = " << theta0;
        std::cout << "    delta  = " << delta << std::endl;
        
        std::cout << "    interp =" << interp_output << 
            "     expected = " << expected_output << 
            "    diff    =" << 1-interp_output/expected_output << std::endl;
        if (expected_output >= 1e-4){ // If the output is expected to be >> 0 then tolerance is stict 
            //EXPECT_NEAR(expected_output, interp_output, 6e-3); // Pass Threshold for 1 deg grids
            EXPECT_NEAR(expected_output, interp_output, 1.5e-2);    // Pass Threshold for 2 deg grids
            //EXPECT_NEAR(expected_output, interp_output, 3e-2);  // Pass Threshold for 5 deg grids
        } else{ // If the output is expected to be very small, we don't need a high tolerance
            EXPECT_NEAR(expected_output, interp_output, 1e-1);  
        }
    }
}


// The grid approach has a problem in gauss case
/*
TEST(Alm_interpTest, testValidInputs_gauss_case) {
    // directory with test grid files: achieves errors smaller than 3e-2 (test on 1000 samples)
    //const std::string grid_dir = "../../../data/Alm_grids_CPP/test_grids/"; // 5 deg resolution grid
    //
    // directory with test grid files: achieves errors guaranteed to be lower than 1.5e-2 (test on 1000 samples)
    const std::string grid_dir = "../../../data/Alm_grids_CPP/2deg_grids/"; 
    //
    // directory with test grid files: achieves 1.5e-3 errors (test on 1000 samples)
    //const std::string grid_dir = "../../../data/Alm_grids_CPP/1deg_grids/"; 
    //
    const std::string ftype="gauss";
    const int lmax=3;
    const int Niter=50; // total number of random tests
    // Test valid inputs and check if the function returns the expected output
    int l; // some input value
    int m; // some input value
    long double theta0; // some input value
    long double delta; // some input value
    long double interp_output;
    long double expected_output;
    // Property based testing to verify that the interpolation is precise enough on a grid
    for(int iter=0; iter<Niter;iter++){
        l=random_int(1,lmax);
        m=random_int(-l,l);
        theta0=random_double(0, M_PI/2);
        delta=random_double(0, M_PI/4);
        interp_output = Alm_interp(l, m, theta0, delta, ftype, grid_dir); // approximation function
        expected_output= Alm(l, m, theta0, delta, ftype); // full integral computation
        // Compare. The tolerance depends on the sparsity of the grid and of the interpolation algo
        std::cout << "["<<iter<<"]         l = " << l;
        std::cout << "         m = " << m;
        std::cout << "    theta0 = " << theta0;
        std::cout << "    delta  = " << delta << std::endl;
        std::cout << "    interp =" << interp_output << 
            "     expected = " << expected_output << 
            "    diff    =" << 1-interp_output/expected_output << std::endl;
        // If the output is expected to be >> 0 then tolerance is strict 
        // Also in the gauss case, there is some larger imprecision around pi/2 and 0 
        // due to the normalisation error. So we do not threat these edges on the same way
        // The threshold is 2.5 degres from the edges
        if (expected_output >= 1e-4 || theta0 < 0.05 || theta0 > 1.5){ 
            EXPECT_NEAR(expected_output, interp_output, 0.1);  // Pass Threshold for 5 deg grids
        } else{ // If the output is expected to be very small, we don't need a high tolerance
            EXPECT_NEAR(expected_output, interp_output, 1);    // 
        }
    }
}
*/

TEST(Alm_interpTest, testValidInputs_triangle_case) {
    // directory with test grid files: achieves errors smaller than 3e-2 (test on 1000 samples)
    //const std::string grid_dir = "../../../data/Alm_grids_CPP/test_grids/"; 
    //
    // directory with test grid files: achieves errors guaranteed to be lower than 1.5e-2 (test on 1000 samples)
    const std::string grid_dir = "../../../data/Alm_grids_CPP/2deg_grids/"; 
    //
    // directory with test grid files: achieves 6e-3 errors (test on 500 samples)
    //const std::string grid_dir = "../../../data/Alm_grids_CPP/1deg_grids/"; 
    //
    const std::string ftype="triangle";
    const int lmax=3;
    const int Niter=50; // total number of random tests
    // Test valid inputs and check if the function returns the expected output
    int l; // some input value
    int m; // some input value
    long double theta0; // some input value
    long double delta; // some input value
    long double interp_output;
    long double expected_output;
    // Property based testing to verify that the interpolation is precise enough on a grid
    for(int iter=0; iter<Niter;iter++){
        l=random_int(1,lmax);
        m=random_int(-l,l);
        theta0=random_double(0, M_PI/2);
        delta=random_double(0, M_PI/4);
        interp_output = Alm_interp(l, m, theta0, delta, ftype, grid_dir); // approximation function
        expected_output= Alm(l, m, theta0, delta, ftype); // full integral computation
        // Compare. The tolerance depends on the sparsity of the grid and of the interpolation algo
        std::cout << "["<<iter<<"]         l = " << l;
        std::cout << "         m = " << m;
        std::cout << "    theta0 = " << theta0;
        std::cout << "    delta  = " << delta << std::endl;
        
        std::cout << "    interp =" << interp_output << 
            "     expected = " << expected_output << 
            "    diff    =" << 1-interp_output/expected_output << std::endl;
        if (expected_output >= 1e-4){ // If the output is expected to be >> 0 then tolerance is stict 
            //EXPECT_NEAR(expected_output, interp_output, 6e-3); // Pass Threshold for 1 deg grids
            EXPECT_NEAR(expected_output, interp_output, 1.5e-2);    // Pass Threshold for 2 deg grids
            //EXPECT_NEAR(expected_output, interp_output, 3e-2);  // Pass Threshold for 5 deg grids
        } else{
             EXPECT_NEAR(expected_output, interp_output, 1e-1);    // Pass Threshold for 2 deg grids           
        }
    }
}

//

// -------- Test for the Alm fast interpolation ----
TEST(LoadAllDataTest, invalidDir){
    const std::string grid_dir = "../../../data/Alm_grids_CPP/test_grids/";
    const std::string ftype = "no_exist";
    GridData_Alm_fast result = loadAllData(grid_dir, ftype);
    EXPECT_TRUE(result.error);
}


TEST(LoadAllDataTest, validInput){
    const std::string grid_dir = "../../../data/Alm_grids_CPP/test_grids/";
    const std::string ftype = "gate";
    GridData_Alm_fast result = loadAllData(grid_dir, ftype);
    EXPECT_FALSE(result.error);
    EXPECT_GT(result.A10.z.size(), 0);
    EXPECT_GT(result.A11.z.size(), 0);
    EXPECT_GT(result.A20.z.size(), 0);
    EXPECT_GT(result.A21.z.size(), 0);
    EXPECT_GT(result.A22.z.size(), 0);
    EXPECT_GT(result.A30.z.size(), 0);
    EXPECT_GT(result.A31.z.size(), 0);
    EXPECT_GT(result.A32.z.size(), 0);
    EXPECT_GT(result.A33.z.size(), 0);
    /*
    if(result.error == false){
        std::cout << "   -----  A10 ----- " << std::endl;
        show_matrix(result.A10);
        std::cout << "   -----  A11 ----- " << std::endl;
        show_matrix(result.A11);
        std::cout << "   -----  A20 ----- " << std::endl;
        show_matrix(result.A20);
        std::cout << "   -----  A21 ----- " << std::endl;
        show_matrix(result.A21);
        std::cout << "   -----  A22 ----- " << std::endl;
        show_matrix(result.A22);
        std::cout << "   -----  A30 ----- " << std::endl;
        show_matrix(result.A30);
        std::cout << "   -----  A31 ----- " << std::endl;
        show_matrix(result.A31);
        std::cout << "   -----  A32 ----- " << std::endl;
        show_matrix(result.A32);
        std::cout << "   -----  A33 ----- " << std::endl;
        show_matrix(result.A33);
    }
    */
}


TEST(Alm_interp4iterTest, testValidInputs_gate_case) {
    // directory with test grid files: achieves errors smaller than 3e-2 (test on 1000 samples)
   //const std::string grid_dir = "../../../data/Alm_grids_CPP/test_grids/"; 
    //
    // directory with test grid files: achieves errors guaranteed to be lower than 9e-3 (test on 1000 samples)
    const std::string grid_dir = "../../../data/Alm_grids_CPP/2deg_grids/"; 
    //
    // directory with test grid files: achieves 2e-3 errors (test on 500 samples)
    //const std::string grid_dir = "../../../data/Alm_grids_CPP/1deg_grids/"; 
    //
    const std::string ftype="gate";
    const int lmax=3;
    const int Niter=50; // total number of random tests
    // Test valid inputs and check if the function returns the expected output
    int l; // some input value
    int m; // some input value
    long double theta0; // some input value
    long double delta; // some input value
    long double interp_output;
    long double expected_output;
    //
    // We prepare the grids...
    GridData_Alm_fast grids=loadAllData(grid_dir, ftype);
    //
    // Property based testing to verify that the interpolation is precise enough on a grid
    for(int iter=0; iter<Niter;iter++){
        l=random_int(1,lmax);
        m=random_int(-l,l);
        theta0=random_double(0, M_PI/2);
        delta=random_double(0, M_PI/4);
        interp_output = Alm_interp_iter(l, m, theta0, delta, ftype, grids); // approximation function
        expected_output= Alm(l, m, theta0, delta, ftype); // full integral computation
        // Compare. The tolerance depends on the sparsity of the grid and of the interpolation algo
        std::cout << "["<<iter<<"]         l = " << l;
        std::cout << "         m = " << m;
        std::cout << "    theta0 = " << theta0;
        std::cout << "    delta  = " << delta << std::endl;
        
        std::cout << "    interp =" << interp_output << 
            "     expected = " << expected_output << 
            "    diff    =" << 1-interp_output/expected_output << std::endl;
        if (expected_output >= 1e-4){ // If the output is expected to be >> 0 then tolerance is stict 
            //EXPECT_NEAR(expected_output, interp_output, 6e-3); // Pass Threshold for 1 deg grids
            EXPECT_NEAR(expected_output, interp_output, 1.5e-2);    // Pass Threshold for 2 deg grids
            //EXPECT_NEAR(expected_output, interp_output, 3e-2);  // Pass Threshold for 5 deg grids
        } else{ // If the output is expected to be very small, we don't need a high tolerance
            EXPECT_NEAR(expected_output, interp_output, 1e-1);  
        }
    }
}

TEST(Alm_interp4iterTest, testValidInputs_triangle_case) {
    // directory with test grid files: achieves errors smaller than 3e-2 (test on 1000 samples)
   //const std::string grid_dir = "../../../data/Alm_grids_CPP/test_grids/"; 
    //
    // directory with test grid files: achieves errors guaranteed to be lower than 9e-3 (test on 1000 samples)
    const std::string grid_dir = "../../../data/Alm_grids_CPP/2deg_grids/"; 
    //
    // directory with test grid files: achieves 2e-3 errors (test on 500 samples)
    //const std::string grid_dir = "../../../data/Alm_grids_CPP/1deg_grids/"; 
    //
    const std::string ftype="triangle";
    const int lmax=3;
    const int Niter=50; // total number of random tests
    // Test valid inputs and check if the function returns the expected output
    int l; // some input value
    int m; // some input value
    long double theta0; // some input value
    long double delta; // some input value
    long double interp_output;
    long double expected_output;
    //
    // We prepare the grids...
    GridData_Alm_fast grids=loadAllData(grid_dir, ftype);
    //
    // Property based testing to verify that the interpolation is precise enough on a grid
    for(int iter=0; iter<Niter;iter++){
        l=random_int(1,lmax);
        m=random_int(-l,l);
        theta0=random_double(0, M_PI/2);
        delta=random_double(0, M_PI/4);
        interp_output = Alm_interp_iter(l, m, theta0, delta, ftype, grids); // approximation function
        expected_output= Alm(l, m, theta0, delta, ftype); // full integral computation
        // Compare. The tolerance depends on the sparsity of the grid and of the interpolation algo
        std::cout << "["<<iter<<"]         l = " << l;
        std::cout << "         m = " << m;
        std::cout << "    theta0 = " << theta0;
        std::cout << "    delta  = " << delta << std::endl;
        
        std::cout << "    interp =" << interp_output << 
            "     expected = " << expected_output << 
            "    diff    =" << 1-interp_output/expected_output << std::endl;
        if (expected_output >= 1e-4){ // If the output is expected to be >> 0 then tolerance is stict 
            //EXPECT_NEAR(expected_output, interp_output, 6e-3); // Pass Threshold for 1 deg grids
            EXPECT_NEAR(expected_output, interp_output, 1.5e-2);    // Pass Threshold for 2 deg grids
            //EXPECT_NEAR(expected_output, interp_output, 3e-2);  // Pass Threshold for 5 deg grids
        } else{ // If the output is expected to be very small, we don't need a high tolerance
            EXPECT_NEAR(expected_output, interp_output, 1e-1);  
        }
    }
}

TEST(Alm_interp4iter_initialisedTest, testValidInputs_gate_case) {
    // directory with test grid files: achieves errors smaller than 3e-2 (test on 1000 samples)
   //const std::string grid_dir = "../../../data/Alm_grids_CPP/test_grids/"; 
    //
    // directory with test grid files: achieves errors guaranteed to be lower than 9e-3 (test on 1000 samples)
    //const std::string grid_dir = "../../../data/Alm_grids_CPP/2deg_grids/"; 
    //
    // directory with test grid files: achieves 2e-3 errors (test on 500 samples)
    const std::string grid_dir = "../../../data/Alm_grids_CPP/1deg_grids/"; 
    //
    const std::string ftype="gate";
    const int lmax=3;
    const int Niter=50; // total number of random tests
    // Test valid inputs and check if the function returns the expected output
    int l; // some input value
    int m; // some input value
    long double theta0; // some input value
    long double delta; // some input value
    long double interp_output;
    long double expected_output;
    //
    // We prepare the grids...
    GridData_Alm_fast grids=loadAllData(grid_dir, ftype);
    // Pre-initialisation of the grid into gsl : Flattening + gsl init
    gsl_funcs funcs_data;
    funcs_data.flat_grid_A10=flatten_grid(grids.A10);
    funcs_data.flat_grid_A11=flatten_grid(grids.A11);
    funcs_data.flat_grid_A20=flatten_grid(grids.A20);
    funcs_data.flat_grid_A21=flatten_grid(grids.A21);
    funcs_data.flat_grid_A22=flatten_grid(grids.A22);
    funcs_data.flat_grid_A30=flatten_grid(grids.A30);
    funcs_data.flat_grid_A31=flatten_grid(grids.A31);
    funcs_data.flat_grid_A32=flatten_grid(grids.A32);
    funcs_data.flat_grid_A33=flatten_grid(grids.A33);
    funcs_data.interp_A10=init_2dgrid(funcs_data.flat_grid_A10);
    funcs_data.interp_A11=init_2dgrid(funcs_data.flat_grid_A11);
    funcs_data.interp_A20=init_2dgrid(funcs_data.flat_grid_A20);
    funcs_data.interp_A21=init_2dgrid(funcs_data.flat_grid_A21);
    funcs_data.interp_A22=init_2dgrid(funcs_data.flat_grid_A22);
    funcs_data.interp_A30=init_2dgrid(funcs_data.flat_grid_A30);
    funcs_data.interp_A31=init_2dgrid(funcs_data.flat_grid_A31);
    funcs_data.interp_A32=init_2dgrid(funcs_data.flat_grid_A32);
    funcs_data.interp_A33=init_2dgrid(funcs_data.flat_grid_A33);
    
    // Property based testing to verify that the interpolation is precise enough on a grid
    for(int iter=0; iter<Niter;iter++){
        l=random_int(1,lmax);
        m=random_int(-l,l);
        theta0=random_double(0, M_PI/2);
        delta=random_double(0, M_PI/4);
        interp_output = Alm_interp_iter_preinitialised(l, m, theta0, delta, ftype, funcs_data); // approximation function
        expected_output= Alm(l, m, theta0, delta, ftype); // full integral computation
        // Compare. The tolerance depends on the sparsity of the grid and of the interpolation algo
        std::cout << "["<<iter<<"]         l = " << l;
        std::cout << "         m = " << m;
        std::cout << "    theta0 = " << theta0;
        std::cout << "    delta  = " << delta << std::endl;
        
        std::cout << "    interp =" << interp_output << 
            "     expected = " << expected_output << 
            "    diff    =" << 1-interp_output/expected_output << std::endl;
        if (expected_output >= 1e-4){ // If the output is expected to be >> 0 then tolerance is stict 
            //EXPECT_NEAR(expected_output, interp_output, 6e-3); // Pass Threshold for 1 deg grids
            EXPECT_NEAR(expected_output, interp_output, 1.5e-2);    // Pass Threshold for 2 deg grids
            //EXPECT_NEAR(expected_output, interp_output, 3e-2);  // Pass Threshold for 5 deg grids
        } else{ // If the output is expected to be very small, we don't need a high tolerance
            EXPECT_NEAR(expected_output, interp_output, 1e-1);  
        }
    }
}

TEST(Alm_interp4iter_initialisedTest, testValidInputs_triangle_case) {
    // directory with test grid files: achieves errors smaller than 3e-2 (test on 1000 samples)
   //const std::string grid_dir = "../../../data/Alm_grids_CPP/test_grids/"; 
    //
    // directory with test grid files: achieves errors guaranteed to be lower than 9e-3 (test on 1000 samples)
    //const std::string grid_dir = "../../../data/Alm_grids_CPP/2deg_grids/"; 
    //
    // directory with test grid files: achieves 2e-3 errors (test on 500 samples)
    const std::string grid_dir = "../../../data/Alm_grids_CPP/1deg_grids/"; 
    //
    const std::string ftype="triangle";
    const int lmax=3;
    const int Niter=50; // total number of random tests
    // Test valid inputs and check if the function returns the expected output
    int l; // some input value
    int m; // some input value
    long double theta0; // some input value
    long double delta; // some input value
    long double interp_output;
    long double expected_output;
    //
    // We prepare the grids...
    GridData_Alm_fast grids=loadAllData(grid_dir, ftype);
    // Pre-initialisation of the grid into gsl : Flattening + gsl init
    gsl_funcs funcs_data;
    funcs_data.flat_grid_A10=flatten_grid(grids.A10);
    funcs_data.flat_grid_A11=flatten_grid(grids.A11);
    funcs_data.flat_grid_A20=flatten_grid(grids.A20);
    funcs_data.flat_grid_A21=flatten_grid(grids.A21);
    funcs_data.flat_grid_A22=flatten_grid(grids.A22);
    funcs_data.flat_grid_A30=flatten_grid(grids.A30);
    funcs_data.flat_grid_A31=flatten_grid(grids.A31);
    funcs_data.flat_grid_A32=flatten_grid(grids.A32);
    funcs_data.flat_grid_A33=flatten_grid(grids.A33);
    funcs_data.interp_A10=init_2dgrid(funcs_data.flat_grid_A10);
    funcs_data.interp_A11=init_2dgrid(funcs_data.flat_grid_A11);
    funcs_data.interp_A20=init_2dgrid(funcs_data.flat_grid_A20);
    funcs_data.interp_A21=init_2dgrid(funcs_data.flat_grid_A21);
    funcs_data.interp_A22=init_2dgrid(funcs_data.flat_grid_A22);
    funcs_data.interp_A30=init_2dgrid(funcs_data.flat_grid_A30);
    funcs_data.interp_A31=init_2dgrid(funcs_data.flat_grid_A31);
    funcs_data.interp_A32=init_2dgrid(funcs_data.flat_grid_A32);
    funcs_data.interp_A33=init_2dgrid(funcs_data.flat_grid_A33);
    // Property based testing to verify that the interpolation is precise enough on a grid
    for(int iter=0; iter<Niter;iter++){
        l=random_int(1,lmax);
        m=random_int(-l,l);
        theta0=random_double(0, M_PI/2);
        delta=random_double(0, M_PI/4);
        interp_output = Alm_interp_iter_preinitialised(l, m, theta0, delta, ftype, funcs_data); // approximation function
        expected_output= Alm(l, m, theta0, delta, ftype); // full integral computation
        // Compare. The tolerance depends on the sparsity of the grid and of the interpolation algo
        std::cout << "["<<iter<<"]         l = " << l;
        std::cout << "         m = " << m;
        std::cout << "    theta0 = " << theta0;
        std::cout << "    delta  = " << delta << std::endl;
        
        std::cout << "    interp =" << interp_output << 
            "     expected = " << expected_output << 
            "    diff    =" << 1-interp_output/expected_output << std::endl;
        if (expected_output >= 1e-4){ // If the output is expected to be >> 0 then tolerance is stict 
            //EXPECT_NEAR(expected_output, interp_output, 6e-3); // Pass Threshold for 1 deg grids
            EXPECT_NEAR(expected_output, interp_output, 1.5e-2);    // Pass Threshold for 2 deg grids
            //EXPECT_NEAR(expected_output, interp_output, 3e-2);  // Pass Threshold for 5 deg grids
        } else{ // If the output is expected to be very small, we don't need a high tolerance
            EXPECT_NEAR(expected_output, interp_output, 1e-1);  
        }
    }
}
