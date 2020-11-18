g++ -Ofast -std=c++11 -I$EIGEN3_INCLUDE_DIR -L/usr/local/lib -I/usr/local/include -L/usr/local/Cellar/boost/1.74.0/lib -I/usr/local/Cellar/boost/1.74.0/include -fopenmp -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lgsl -lgslcblas artificial_spectrum.cpp io_star_params.cpp build_lorentzian.cpp function_rot.cpp plots_diags.cpp  iterative_artificial_spectrum.cpp models_database.cpp  random_JB.cpp string_handler.cpp noise_models.cpp bump_DP.cpp solver_mm.cpp interpol.cpp derivatives_handler.cpp -o ./sim.out


#g++ -Ofast -std=c++11 -I$EIGEN3_INCLUDE_DIR -L/usr/local/lib -I/usr/local/include artificial_spectrum.cpp io_star_params.cpp build_lorentzian.cpp function_rot.cpp plots_diags.cpp  iterative_artificial_spectrum.cpp models_database.cpp  random_JB.cpp string_handler.cpp noise_models.cpp bump_DP.cpp solver_mm.cpp interpol.cpp derivatives_handler.cpp -o ./sim.out


