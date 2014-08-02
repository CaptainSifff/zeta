cmake_minimum_required(VERSION 2.6)
project(zeta CXX)
SET(CMAKE_CXX_FLAGS_RELEASE "-O -march=corei7 -mtune=corei7 -Wconversion -pipe")
#SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ftree-vectorizer-verbose=2 -march=corei7 -g")
SET(CMAKE_CXX_FLAGS_DEBUG "-march=core2 -gstabs+ -Wall -Wabi -Wconversion")
SET(CMAKE_CXX_FLAGS_PROFILED "-O -march=core2 -pg -Wconversion")

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_CONFIG_H -I. -I..")

IF(CMAKE_COMPILER_IS_GNUCXX)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

add_executable(zeta main.cpp)

install(TARGETS zeta RUNTIME DESTINATION bin)
