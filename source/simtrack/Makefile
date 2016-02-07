ARCH = $(shell getconf LONG_BIT)
GSL_INC= $(HOME)/usr/local/gsl/include

OBJS=smath.o tpsa.o elepass.o tracking.o physics.o interface.o
all: $(OBJS) libsimtrack$(ARCH).a

tpsa.o: tpsa.cpp tpsa.h
	g++ -O3 -c $<
smath.o: smath.cpp smath.h
	g++ -O3 -c $< -I$(GSL_INC)
elepass.o: elepass.cpp elepass.h global.h smath.h tpsa.h
	g++ -O3 -c $<
tracking.o: tracking.cpp tracking.h global.h elepass.h smath.h
	g++ -O3 -c $<
physics.o: physics.cpp physics.h global.h elepass.h tracking.h
	g++ -O3 -c $<
interface.o: interface.cpp interface.h global.h elepass.h tracking.h
	g++ -O3 -c $<


libsimtrack$(ARCH).a: $(OBJS)
	ar -r $@ $(OBJS)

clean:
	rm -f *.o libsimtrack$(ARCH).a

