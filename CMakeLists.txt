cmake_minimum_required(VERSION 3.1)

project(minnow)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# We require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(KSW_FLAGS "-march=native -DHAVE_KALLOC")
set(RBF_CPP_FLAGS "-mbmi2 -msse4 -pthread -std=c++11 -W -Wall -Wextra -Wpointer-arith -Wunused -Wwrite-strings -openmp -Wno-unknown-pragmas -Wno-unused-function -O3")
#set(WARN_ALL_THINGS "-fdiagnostics-color=always -Wall -Wcast-align -Wcast-qual -Wconversion -Wctor-dtor-privacy -Wdisabled-optimization -Wdouble-promotion -Wextra -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wno-sign-conversion -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wpedantic -Wredundant-decls -Wshadow -Wstrict-aliasing=1 -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unknown-pragmas -Wuseless-cast") 

#set(WARN_ALL_THINGS "-fdiagnostics-color=always -Wall -Wcast-align -Wcast-qual -Wconversion -Wctor-dtor-privacy -Wdisabled-optimization -Wdouble-promotion -Wduplicated-branches -Wduplicated-cond -Wextra -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wno-sign-conversion -Wnoexcept -Wnull-dereference -Wold-style-cast -Woverloaded-virtual -Wpedantic -Wredundant-decls -Wrestrict -Wshadow -Wstrict-aliasing=1 -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unknown-pragmas -Wuseless-cast") 


#set(OPT_FLAGS "-O3 -g -ffast-math -fPIC -funroll-loops -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing")
#set(OPT_FLAGS "-fprofile-generate -O3 -DNDEBUG -funroll-loops -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing")
#set(OPT_FLAGS "-fprofile-use -fprofile-correction -O3 -DNDEBUG -funroll-loops -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing")
set(DEBUG_FLAGS "-pg -g -gstabs -no-pie")

set(RBF_CPP_FLAGS "${KSW_FLAGS} ${RBF_CPP_FLAGS} ${OPT_FLAGS}")

set(CMAKE_CXX_FLAGS ${RBF_CPP_FLAGS})
message("${CMAKE_CXX_FLAGS}")

#find_package(JeMalloc)
#    if(JEMALLOC_FOUND)
#      include_directories(SYSTEM ${JEMALLOC_INCLUDE_DIRS})
#    endif()
#
#
#include(ExternalProject)

set(GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

message("Build system will NOT compile libgff")
##
# Let us check the sha sum of our pacakges if we have the right tools
##
#set(SHASUM ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_shasum.sh)

#message("Build system will compile libgff")
#message("==================================================================")
#ExternalProject_Add(libgff
#    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
#    DOWNLOAD_COMMAND curl -k -L https://github.com/COMBINE-lab/libgff/archive/v1.1.tar.gz -o libgff.tgz &&
#                     ${SHASUM} 34597ecd3718e6d6e7a1dd4350220664be2c1cd992070d1218564a20d7aa5697  libgff.tgz &&
#                     tar -xzvf libgff.tgz 
#    ##
#    #URL https://github.com/COMBINE-lab/libgff/archive/v1.1.tar.gz
#    #DOWNLOAD_NAME libff.tgz
#    #URL_HASH SHA1=37b3147d78391d1fabbe6a0df313fbf516abbc6f
#    #TLS_VERIFY FALSE
#    ##
#    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff-1.1
#    #UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
#    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
#    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff-1.1/build
#    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_SOURCE_DIR}/external/install
#)
#ExternalProject_Add_Step(libgff makedir
#  COMMAND mkdir -p <SOURCE_DIR>/build 
#  COMMENT "Make build directory"
#  DEPENDEES download 
#  DEPENDERS configure)


###
#
# Done building external dependencies.
#
###

set (CPACK_SOURCE_IGNORE_FILES
"/external/"
".git/")


include_directories(
  include
#  external/install/include
)

link_directories(
  lib
#  external/install/lib
)
add_subdirectory(src)
