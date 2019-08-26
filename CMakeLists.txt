cmake_minimum_required(VERSION 3.13)
project(SIMPQS)
enable_language(C)
set(CMAKE_C_STANDARD 11)
#include_directories(/pari/bin/include/pari)
#include_directories(pari)
#add_executable(PariTrys ecm_basic.gp.c ../master/masterMain.c ../worker/workerMain.c ../utils.c ../utils.h)
include_directories(.)
set(SOURCES_COMMON
        utils/utils.h
        utils/utils.c
        utils/gmp_patch.c
        utils/gmp_patch.h
        SIMPQS.h
        CONFIGURATION.h)
set(SOURCES_WORKER
        factorization/factorizerQuick.c
        factorization/factorizerQuick.h
        worker/sievingSIMPQS.c
        worker/sievingSIMPQS.h
        )
set(SOURCES_MASTER
        master/masterMain.c
        )

add_executable(SIMQS_Worker ${SOURCES_COMMON} ${SOURCES_WORKER})
add_executable(SIMQS_Master ${SOURCES_COMMON} ${SOURCES_MASTER})

set(CMAKE_C_FLAGS "-Wall -lgmp -pthread -lm -ggdb")
#target_link_libraries(PariTrys /pari)
#include(pari/bin/include/pari/pari.h)

