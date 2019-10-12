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
        CONFIGURATION.h utils/precomputations.c utils/precomputations.h)
set(SOURCES_FACTORIZE factorization/factorizerQuick.h factorization/factorizerQuick.c)
set(SOURCES_MATRIX matrix/matrix.c matrix/matrix.h)
set(SOURCES_WORKER
        worker/workerMain.c
        worker/sievingSIMPQS.c
        worker/sievingSIMPQS.h
        ${SOURCES_FACTORIZE}
        ${SOURCES_MATRIX}
        worker/reports.c worker/report.h)
set(SOURCES_MASTER
        master/masterMain.c
        )
set(GCC_FLAGS
        "-Wall -Wextra -lgmp -lmpfr -pthread -lm -ggdb -O2"
        )
add_executable(SIMQS_Worker ${SOURCES_COMMON} ${SOURCES_WORKER})
add_executable(SIMQS_Master ${SOURCES_COMMON} ${SOURCES_MASTER})
add_executable(TRIAL_DIVISION_MAIN_TEST ${SOURCES_FACTORIZE} ${SOURCES_COMMON} )
add_executable(TRIAL_MULTIPLE_POLYNOMIAL ${SOURCES_COMMON} worker/report.h)
add_executable(TRIAL try.c)
set(CMAKE_C_FLAGS ${GCC_FLAGS})
#target_link_libraries(PariTrys /pari)
#include(pari/bin/include/pari/pari.h)

