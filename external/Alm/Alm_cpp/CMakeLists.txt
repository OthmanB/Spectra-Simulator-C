cmake_minimum_required(VERSION 3.5)

project(Alm VERSION 2.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_USE_MULTITHREADED TRUE)
set(CMAKE_CXX_FLAGS "-I/usr/include/  -O2")
set(CMAKE_C_FLAGS "-I/usr/include/  -O2")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../../bin)

#Bring the core program sources into the project
set(SOURCES_MAIN  #_COMMON
    main.cpp
    activity.cpp
    Alm_interpol.cpp
    bilinear_interpol.cpp
    string_handler.cpp
    GaussLegendre1D.hpp
    GaussLegendre2D.hpp
    )

set(SOURCES_INTEGRAND  
    get_filter.cpp
    activity.cpp
    string_handler.cpp
    GaussLegendre1D.hpp
    GaussLegendre2D.hpp
    )

set(SOURCES_GRIDMAKER 
    main_do_grids.cpp
    do_grids.cpp
    make_grids.cpp
    linspace.cpp  # linspace is still required here because it is using linspace_vec (not Eigen-based vectors)
    gzip_compress.cpp
    activity.cpp
    string_handler.cpp
    GaussLegendre1D.hpp
    GaussLegendre2D.hpp
    )

set(SOURCES_TMP
    gzip_compress.cpp
    Alm_interpol.cpp
    bilinear_interpol.cpp
    string_handler.cpp
    )

#Look for eigen and explicitly specify to use it. EIGEN3_INCLUDE_DIR Must be specified in the bashrc
include_directories( "$ENV{EIGEN3_INCLUDE_DIR}" )

#Dependencies that are not portable (BOOST)
find_package(Boost REQUIRED COMPONENTS iostreams program_options filesystem)
find_package(GSL REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES_MAIN})
target_link_libraries(${CMAKE_PROJECT_NAME} ${Boost_LIBRARIES}  ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})

add_executable(Filter_Alm ${SOURCES_INTEGRAND})
target_link_libraries(Filter_Alm Boost::boost)

add_executable(GridMaker ${SOURCES_GRIDMAKER})
target_link_libraries(GridMaker ${Boost_LIBRARIES})

## Add your source files here
#add_executable(Alm_interpol ${SOURCES_TMP})
## Link the required libraries
#target_link_libraries(Alm_interpol ${Boost_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})

