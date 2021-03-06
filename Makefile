CC=gcc
#libs
CFLAGS=-lgmp -lmpfr -pthread -lm -Wall     #TODO OLD
#add header in path of compilation
CFLAGS+= -I .

###GDB MAX COMPATIBILITY
CFLAGS+= -ggdb -Wall
##GCC OPTIMIZATION
#CFLAGS+= -O2
CFLAGS+=  -D _LARGEFILE64_SOURCE -D _GNU_SOURCE     #dep for fcntl and seek64

#TODO APP BASIC CONFIGS #######################################
###FIXED TIMEOUT TODO COMMENT FOR ADAPTIVE TIMEOUT
#CFLAGS+= -D TIMEOUTFIX
###USE DEFAULT (MACRO DEFINED) FOR TX CONFIGURATION TODO DECOMMENT
#CFLAGS+= -D TX_CONFIG_FIXED
###########################################################

################		TESTS EXTRA CONFIG	#############
TESTCONFIG_MACRO =
TESTCONFIG_MACRO+= -D TEST_QUIET_PRINT		#DISABLE SOME PRINT DURING TEST
### PCKPAYLOAD WILL BE READED FROM A FILE NAMED pcksize,
### RELATED MACRO WILL BE SETTED IF THE FILE IS FOUNDED
#_PCKPAYLOAD=$(shell cat pcksize)
#ifneq ($(_PCKPAYLOAD),)						#TRUE IF NOT EMPTY=> FOUNDED DATA IN FILE
#TESTCONFIG_MACRO+= -D PCKPAYLOADSIZE=$(_PCKPAYLOAD)	#set pcksize from file
#endif
###TODO NEXT MACRO TRIGGER GCC ATOMIC VERSION OF SENDER/RECEIVER
#TESTCONFIG_MACRO+= -D GCCATOMICS		#decomment for gcc atomics trigger on test
#############################################################

##basic source files vars
WORKER=$(shell find worker -iname "*.c" )
CFILES_COMMON=$(shell find utils -iname "*.c" )
CFILES_COMMON+=$(shell find utils -iname "*.h" )
CFILES_COMMON+=SIMPQS.h CONFIGURATION.h
CFILES_FACTORIZE= factorization/factorizerQuick.h factorization/factorizerQuick.c matrix/matrix.c matrix/matrix.h
HEADERS=$(shell find -iname "*.h" )
worker.o:  $(WORKER) $(CFILES_COMMON) $(CFILES_FACTORIZE) $(HEADERS)
	$(CC) -o $@ $(WORKER) $(CFLAGS)  $(CFILES_COMMON) $(CFILES_FACTORIZE)
print:
	@echo ${TESTCONFIG_MACRO} $(_EXIT_SUCCESS)
clean:
	rm worker.o


