cmake_minimum_required(VERSION 3.13)
project(SIMPQS)
enable_language(C)
set(CMAKE_C_STANDARD 11)
#include_directories(/pari/bin/include/pari)
#include_directories(pari)
#add_executable(PariTrys ecm_basic.gp.c ../master/masterMain.c ../worker/workerMain.c ../utils.c ../utils.h)
include_directories(.)
set(SOURCES
        utils/utils.h
        utils/utils.c
        utils/gmp_patch.c
        utils/gmp_patch.h SIMPQS.h)
add_executable(SIMQS_Worker ${SOURCES} worker/workerMain.c worker/sievingSIMPQS.c worker/sievingSIMPQS.h factorizerQuick.c factorizerQuick.h)
add_executable(SIMQS_Master ${SOURCES} master/masterMain.c CONFIGURATION.h)

set(CMAKE_C_FLAGS "-Wall -lgmp")
#target_link_libraries(PariTrys /pari)
#include(pari/bin/include/pari/pari.h)

