CC = g++
ARCH = $(shell getconf LONG_BIT)
OPTS = -O3
#/home/hchao/usr/local/lib,include have NLOPT, GALIB, BOOST libraries and header files.

#the following is for armadillo
ARMA_INCLUDE=-I$(HOME)/usr/include
ARMA_LIBRARY=-L$(HOME)/usr/lib64 -larmadillo
#the following is for simtrack
SIMTRACK_INC=-I./simtrack
SIMTRACK_LIB=-lm -lgsl -lgslcblas -L./simtrack -lsimtrack$(ARCH)
#the following is for GAlib
GA_INCLUDE=-I/home/hchao/usr/local/include
GA_LIBRARY=-L/home/hchao/usr/local/lib -lm -lga
NLOPT_INCLUDE=-I/home/hchao/usr/local/include
NLOPT_LIBRARY=-L/home/hchao/usr/local/lib -lm -lnlopt
BOOST_INCLUDE=-I/home/hchao/usr/local/include
BOOST_LIBRARY=-L/home/hchao/usr/local/lib 

#all: nlopt2 loadQ_track Qerrgen
all: run62_omp

#basic math, random number, distribution generation
myroutine.o: myroutine.cpp myroutine.h
	$(CC) $(OPTS) -c $<
#complmentary error function, space charge kick
myroutine1.o: myroutine1.cpp myroutine1.h myroutine.cpp myroutine.h Faddeeva.cc Faddeeva.hh
	$(CC) $(OPTS) -c $<
#ring operations
myroutine2.o: myroutine2.cpp myroutine2.h myroutine.cpp myroutine.h
	$(CC) $(OPTS) -c $< $(ARMA_INCLUDE)
Faddeeva.o: Faddeeva.cc Faddeeva.hh
	$(CC) $(OPTS) -c $<

run62_omp.o: run62_omp.cpp 
	$(CC) $(OPTS) -fopenmp -c $< $(SIMTRACK_INC)
run62_omp: run62_omp.o FODO4.o myroutine.o myroutine1.o myroutine2.o Faddeeva.o simtrack/libsimtrack$(ARCH).a
	$(CC) -fopenmp $^ $(SIMTRACK_INC) $(SIMTRACK_LIB) $(ARMA_INCLUDE) $(ARMA_LIBRARY) -o $@

#Qerrgen: Qerrgen.cpp FODO4.o myroutine.o myroutine2.o simtrack/libsimtrack$(ARCH).a
#	$(CC) $^ $(SIMTRACK_INC) $(SIMTRACK_LIB) $(ARMA_INCLUDE) $(ARMA_LIBRARY) -o $@
#nlopt_obj1.o: nlopt_obj1.cpp nlopt_obj1.h
#	$(CC) $(OPTS) -c $< $(ARMA_INCLUDE) $(SIMTRACK_INC) $(BOOST_INCLUDE) #include the last to fix split_regex problem

#nlopt2.o: nlopt2.cpp 
#	$(CC) $(OPTS) -c $< $(SIMTRACK_INC) $(ARMA_INCLUDE) $(NLOPT_INCLUDE) 
#nlopt2: nlopt2.o nlopt_obj1.o FODO4.o myroutine.o myroutine2.o simtrack/libsimtrack$(ARCH).a
#	$(CC) $^ $(SIMTRACK_INC) $(SIMTRACK_LIB) $(NLOPT_INCLUDE) $(NLOPT_LIBRARY) $(ARMA_INCLUDE) $(ARMA_LIBRARY) -lboost_regex -o $@

FODO4.o: FODO4.cpp FODO4.h 
	$(CC) $(OPTS) -c $< $(SIMTRACK_INC)


clean:
	rm -f *.o 
lib:
	cd simtrack; make
libclean:
	cd simtrack; make clean

.PHONY: clean lib libclean
