cmake_minimum_required(VERSION 3.4.1)

project(Alm VERSION 1.3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Boost_USE_MULTITHREADED TRUE)
set(CMAKE_CXX_FLAGS "-I/usr/include/  -O2")
set(CMAKE_C_FLAGS "-I/usr/include/  -O2")

#Bring the headers into the project
#set(HEADERS  #_COMMON
#	
#   )

#Bring the core program sources into the project
set(SOURCES_MAIN  #_COMMON
    main.cpp
    activity.cpp
    linspace.cpp
    GaussLegendre1D.hpp
    GaussLegendre2D.hpp
    )

#Look for eigen and explicitly specify to use it. EIGEN3_INCLUDE_DIR Must be specified in the bashrc
include_directories( "$ENV{EIGEN3_INCLUDE_DIR}" )

#Dependencies that are not portable (BOOST)
#find_package(Boost COMPONENTS system filesystem iostreams math REQUIRED)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

add_executable(${CMAKE_PROJECT_NAME} ${HEADERS} ${SOURCES_MAIN} ${DIAG_SRC})
#target_link_libraries(${CMAKE_PROJECT_NAME}  ${Boost_SYSTEM_LIBRARY_RELEASE} ${Boost_MATH_LIBRARY_RELEASE} ${Boost_FILESYSTEM_LIBRARY_RELEASE} ${Boost_IOSTREAMS_LIBRARY_RELEASE} Boost::boost)
target_link_libraries(${CMAKE_PROJECT_NAME} Boost::boost)
