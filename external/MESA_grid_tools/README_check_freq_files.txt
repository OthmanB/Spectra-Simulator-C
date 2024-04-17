Purpose:
     (1) Read a MESA summary file as defined by Kuldeep
     (2) Goes through a directory and check if freq files exist for each Model ID of the summary file
     (3) Report missing freq files in an output file

How to compile:
g++ -O3 -I ../eigen -I ../ -fopenmp -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lgsl -lgslcblas ../format.cpp ../write_star_params.cpp ../ioproc.cpp check_freq_files.cpp -o check.out
