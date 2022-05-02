g++ -O3 -std=c++11 -o test.prg main.cpp bump_DP.cpp solver_mm.cpp interpol.cpp string_handler.cpp derivatives_handler.cpp noise_models.cpp linfit.cpp linspace.cpp -I$EIGEN3_INCLUDE_DIR -L/usr/local/lib -I/usr/local/include -lboost_system -lboost_filesystem -Xpreprocessor -fopenmp -lomp
#-Wall
