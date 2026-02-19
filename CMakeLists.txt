cmake_minimum_required(VERSION 3.5)

project(specsim VERSION 1.1)

# Check for macOS
if(APPLE)
    set(MACOS 1)
else()
    set(MACOS 0)
endif()

# Check for Linux
if(CMAKE_COMPILER_IS_GNUCXX)
    set(LINUX 1)
else()
    set(LINUX 0)
endif()

if(NOT CMAKE_COMPILER_IS_GNUCXX AND NOT MACOS)
  message(FATAL_ERROR "This program works only on Mac or Linux")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_USE_MULTITHREADED TRUE)
option(WITH_OPENMP "OpenMP Library" ON)
option(BUILD_ON_DALMA "Building for NYU DALMA Supercomputer" OFF)

# Configure the header file
# configure_file(rgen_config.h.in rgen_config.h)

if(NOT BUILD_ON_DALMA)
	if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")
		#set(CMAKE_CXX_FLAGS "-O2 -fopenmp")
	endif()
	if (APPLE)
    set(CMAKE_CXX_FLAGS "-O2")
		#set(CMAKE_CXX_FLAGS "-O2 -fopenmp=libomp")
	endif()
else()
	message(STATUS "Building for NYU DALMA Supercomputer...")
	message(WARNING "   - Remember to load the modules required for DALMA before running cmake")
	message(STATUS  "   - Adding optimisation flags for DALMA...")
	# SSE / AVX switch
	option(USE_AVX "Build with AVX support" ON)
	if(USE_AVX)
	  set(CMAKE_CXX_FLAGS "-fopenmp -O3 -mavx2 -mfma -march=haswell")
	else()
	  set(CMAKE_CXX_FLAGS "-O3 -msse4.2 -march=westmere")
	endif()
endif()


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
	    logging.cpp
	    rng.cpp
	    combi.cpp
    stellar_models.cpp
    models_database_grid.cpp
    format.cpp
    acoefs.cpp
    external/Alm/Alm_cpp/activity.cpp
    external/ARMM/bump_DP.cpp
    external/ARMM/solver_mm.cpp
    external/rescale/rescale_freqs.cpp
    external/rescale/decompose_nu.cpp
    external/rescale/data.h
    )

set(SOURCES_TEST_LORENTZIAN
    test_fct/lorentzian_tests.cpp
    function_rot.cpp
    ioproc.cpp
    build_lorentzian.cpp
    linspace.cpp
    acoefs.cpp
    external/Alm/Alm_cpp/activity.cpp
    )


# Check and Load dependencies
# NOTE: gnuplot is only required at *runtime* when plots are requested (doplots=1).
# Some environments (CI, minimal containers) can build and run the simulator with doplots=0
# without having the gnuplot executable installed. Therefore, we do not hard-require it here.
find_package(Gnuplot)
if(NOT GNUPLOT_FOUND)
	message(WARNING "Gnuplot executable not found. Plotting will fail if doplots=1. Install gnuplot or run with doplots=0")
endif()

if (WITH_OPENMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
	    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
      if(WITH_OPENMP)
        string(APPEND CMAKE_CXX_FLAGS_RELEASE " -fopenmp")
      endif()
  endif()
endif()


#
include_directories()


# Dependencies that are not portable (BOOST)
#
# Note: CMake >= 4 removed FindBoost. Some Boost distributions (e.g. Homebrew)
# ship only `BoostConfig.cmake` without per-component `boost_<component>Config.cmake`
# files. In that case, `find_package(Boost COMPONENTS ...)` fails.
#
# To keep builds working across environments, we locate Boost headers and the
# required component libraries manually.

find_path(BOOST_INCLUDE_DIR boost/filesystem.hpp)
if(NOT BOOST_INCLUDE_DIR)
	message(FATAL_ERROR "Boost headers not found (expected boost/filesystem.hpp). Install Boost or set CMAKE_PREFIX_PATH.")
endif()

find_library(BOOST_FILESYSTEM_LIBRARY NAMES boost_filesystem)
find_library(BOOST_PROGRAM_OPTIONS_LIBRARY NAMES boost_program_options)
find_library(BOOST_IOSTREAMS_LIBRARY NAMES boost_iostreams)
find_library(BOOST_SYSTEM_LIBRARY NAMES boost_system)

if(NOT BOOST_FILESYSTEM_LIBRARY)
	message(FATAL_ERROR "Boost.Filesystem library not found (libboost_filesystem).")
endif()
if(NOT BOOST_PROGRAM_OPTIONS_LIBRARY)
	message(FATAL_ERROR "Boost.ProgramOptions library not found (libboost_program_options).")
endif()
if(NOT BOOST_IOSTREAMS_LIBRARY)
	message(FATAL_ERROR "Boost.Iostreams library not found (libboost_iostreams).")
endif()

set(BOOST_LIBRARIES
	${BOOST_FILESYSTEM_LIBRARY}
	${BOOST_PROGRAM_OPTIONS_LIBRARY}
	${BOOST_IOSTREAMS_LIBRARY}
)
if(BOOST_SYSTEM_LIBRARY)
	list(APPEND BOOST_LIBRARIES ${BOOST_SYSTEM_LIBRARY})
endif()

include_directories(${BOOST_INCLUDE_DIR})

add_executable(${CMAKE_PROJECT_NAME}  ${SOURCES_MAIN})
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
if(MACOS)
  find_package(Eigen3 REQUIRED NO_MODULE)
  target_compile_definitions(${CMAKE_PROJECT_NAME}  PRIVATE MACOS)
  target_link_libraries(${CMAKE_PROJECT_NAME}  Eigen3::Eigen  ${BOOST_LIBRARIES} "-framework CoreFoundation")
elseif(LINUX)
  #Look for eigen and explicitly specify to use it. EIGEN3_INCLUDE_DIR Must be specified in the bashrc
  include_directories( "$ENV{EIGEN3_INCLUDE_DIR}" )
  target_compile_definitions(${CMAKE_PROJECT_NAME}  PRIVATE LINUX)
  target_link_libraries(${CMAKE_PROJECT_NAME}  ${BOOST_LIBRARIES}  rt)
endif()
