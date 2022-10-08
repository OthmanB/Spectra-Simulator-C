cmake_minimum_required(VERSION 3.4.1)

project(specsim VERSION 1.0.4)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_USE_MULTITHREADED TRUE)
option(WITH_OPENMP "OpenMP Library" ON)
option(BUILD_ON_DALMA "Building for NYU DALMA Supercomputer" OFF)

if(NOT BUILD_ON_DALMA)
	if(CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_CXX_FLAGS "-O2 -fopenmp")
	endif()
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		set(CMAKE_CXX_FLAGS "-O2 -fopenmp=libomp")
	endif()
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
		set(CMAKE_CXX_FLAGS "-O2")
	endif()
else()
	message(STATUS "Building for NYU DALMA Supercomputer...")
	message(WARNING "   - Remember to load the modules required for DALMA before running cmake")
	message(STATUS  "   - Adding optimisation flags for DALMA...")
	# SSE / AVX switch
	option(USE_AVX "Build with AVX support" ON)
	if(USE_AVX)
	  set(CMAKE_CXX_FLAGS "-fopenmp -O3 -mavx2 -mfma -march=haswell")
	  #set(CMAKE_CXX_FLAGS "-O3 -mavx2 -mfma -march=haswell")
	else()
	  set(CMAKE_CXX_FLAGS "-O3 -msse4.2 -march=westmere")
	endif()
endif()

#Bring the headers into the project
set(HEADERS  #_COMMON
    artificial_spectrum.h
    data.h
    gnuplot-iostream.h
    models_database.h
    derivatives_handler.h
    interpol.h
    linfit.h
    linspace.h
    plots_diags.h
    version.h
    build_lorentzian.h
    function_rot.h
    io_star_params.h
    io_star_common.h
    noise_models.h
    ioproc.h
    external/integrate/activity.h
    combi.h
    stellar_models.h
    models_database_grid.h
    format.h
    external/ARMM/bump_DP.h
    external/ARMM/solver_mm.h
    external/ARMM/data_solver.h
    acoefs.h
    )

#Bring the core program sources into the project
set(SOURCES_MAIN  #_COMMON
	iterative_artificial_spectrum.cpp
	artificial_spectrum.cpp
    models_database.cpp
    derivatives_handler.cpp
    interpol.cpp
    linfit.cpp
    linspace.cpp
    plots_diags.cpp
    build_lorentzian.cpp
    function_rot.cpp
    io_star_params.cpp
    io_star_common.cpp
    noise_models.cpp
    ioproc.cpp
    random_JB.cpp
    combi.cpp
    stellar_models.cpp
    models_database_grid.cpp
    format.cpp
    acoefs.cpp
    external/integrate/activity.cpp
    external/ARMM/bump_DP.cpp
    external/ARMM/solver_mm.cpp
    )

set(SOURCES_TEST_LORENTZIAN
    test_fct/lorentzian_tests.cpp
    function_rot.cpp
    ioproc.cpp
    build_lorentzian.cpp
    linspace.cpp
    acoefs.cpp
    external/integrate/activity.cpp
    )


# Check and Load dependencies
find_package(Gnuplot REQUIRED)
if(NOT GNUPLOT_FOUND)
	message(FATAL_ERROR "Gnuplot is required for rendering. Cannot pursue. Please install gnuplot")
endif()

if (WITH_OPENMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
	    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	else()
		message(FATAL_ERROR "OpenMP not found. Use -DWITH_OPENMP=OFF to build without OpenMP")
	endif()
endif()

#Look for eigen and explicitly specify to use it. EIGEN3_INCLUDE_DIR Must be specified in the bashrc
include_directories( "$ENV{EIGEN3_INCLUDE_DIR}" )
include_directories()

#Dependencies that are not portable (BOOST)
set(Boost_USE_MULTITHREADED=OFF)
find_package(Boost COMPONENTS system filesystem iostreams REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

add_executable(${CMAKE_PROJECT_NAME} ${HEADERS} ${SOURCES_MAIN})
target_link_libraries(${CMAKE_PROJECT_NAME}  ${Boost_SYSTEM_LIBRARY_RELEASE} ${Boost_FILESYSTEM_LIBRARY_RELEASE} ${Boost_IOSTREAMS_LIBRARY_RELEASE})

add_executable(test_lorentzians ${HEADERS} ${SOURCES_TEST_LORENTZIAN})
target_link_libraries(${CMAKE_PROJECT_NAME}  ${Boost_SYSTEM_LIBRARY_RELEASE} ${Boost_FILESYSTEM_LIBRARY_RELEASE} ${Boost_IOSTREAMS_LIBRARY_RELEASE})



