cmake_minimum_required(VERSION 3.5)

project(ARMMSolver VERSION 1.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_USE_MULTITHREADED TRUE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include/ -O2 -fopenmp")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

#Bring the core program sources into the project
set(SOURCES_MAIN  #_COMMON
    do_solve.cpp
    bump_DP.cpp
    solver_mm.cpp
    interpol.cpp
    string_handler.cpp
    derivatives_handler.cpp
    noise_models.cpp
    linfit.cpp
    readparams_job.cpp
    )

#Bring the core program sources into the project
set(SOURCES_SCAN  
    scanrange_solver.cpp
    bump_DP.cpp
    solver_mm.cpp
    interpol.cpp
    string_handler.cpp
    derivatives_handler.cpp
    noise_models.cpp
    linfit.cpp
    readparams_job.cpp
    )

#Bring the core program sources into the project
set(SOURCES_MAKE_STAR  
configure_make_star.cpp
make_star.cpp
bump_DP.cpp
solver_mm.cpp
interpol.cpp
string_handler.cpp
derivatives_handler.cpp
noise_models.cpp
linfit.cpp
readparams_job.cpp
writeparams_job.cpp
)

if(APPLE)
    find_package(OpenMP REQUIRED COMPONENTS C)
    find_package(Eigen3 REQUIRED)
else()
    find_package(OpenMP REQUIRED)
endif()


#Look for eigen and explicitly specify to use it. EIGEN3_INCLUDE_DIR Must be specified in the bashrc
include_directories( "$ENV{EIGEN3_INCLUDE_DIR}" )

#Dependencies that are not portable (BOOST)
find_package(Boost COMPONENTS system filesystem iostreams program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES_MAIN} ${DIAG_SRC})
target_link_libraries(${CMAKE_PROJECT_NAME} ${Boost_LIBRARIES})

add_executable(scan_range_solver ${SOURCES_SCAN} ${DIAG_SRC})
target_link_libraries(scan_range_solver ${Boost_LIBRARIES})

add_executable(make_star ${SOURCES_MAKE_STAR} ${DIAG_SRC})
target_link_libraries(make_star ${Boost_LIBRARIES})
