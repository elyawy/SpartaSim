cmake_minimum_required(VERSION 3.4...3.18)
project(Sparta)

add_subdirectory(pybind11)
ADD_LIBRARY(LibsModule 
    SpartaSRC/Block.cpp
    SpartaSRC/DistTable.cpp
    SpartaSRC/errorMsg.cpp
    SpartaSRC/euclidean_distance.cpp
    SpartaSRC/Event.cpp
    SpartaSRC/FastZip.cpp
    SpartaSRC/geometricDist.cpp
    SpartaSRC/lavaletteDist.cpp
    SpartaSRC/LinkedList.cpp
    SpartaSRC/logFile.cpp
    SpartaSRC/MSA.cpp
    SpartaSRC/MSA_decomposer.cpp
    SpartaSRC/Node.cpp
    SpartaSRC/Parameters.cpp
    SpartaSRC/param_priors.cpp
    SpartaSRC/poissonDist.cpp
    SpartaSRC/read_seqs.cpp
    SpartaSRC/readTree.cpp
    SpartaSRC/simulator.cpp
    SpartaSRC/SpartaABC_options.cpp
    SpartaSRC/summ_stats_wrapper.cpp
    SpartaSRC/tree.cpp
    SpartaSRC/treeIt.cpp
    SpartaSRC/treeUtil.cpp
    SpartaSRC/zipfDist.cpp
    SpartaSRC/ConversionUtils.cpp
    SpartaSRC/someUtil.cpp
    SpartaSRC/alphabet.cpp
    SpartaSRC/talRandom.cpp
)

pybind11_add_module(Sparta src/main.cpp)
target_link_libraries(Sparta PRIVATE LibsModule)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(Sparta
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall --std=c++11 -O3 -fPIC" )
set( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -O3 -fPIC" )
