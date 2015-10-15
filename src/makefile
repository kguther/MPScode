CC=icpc
LFLAGS=-llapacke -llapack -lblas -lgfortran
DCFLAGS=-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

all: #main.o parameters.o network.o Qoperator.o Qstate.o
	$(CC) -o MPScode main.cpp parameters.cpp $(LFLAGS) $(DCFLAGS) -Wall -g

main.o: main.cpp network.h Qstate.h Qoperator.h arraycreation.h arrayprocessing.h
	$(CC) -c main.cpp $(LFLAGS) $(DCFLAGS)

parameters.o: parameters.cpp parameters.h
	$(CC) -c parameters.cpp

network.o: network.cpp network.h arraycreation.h arrayprocessing.h Qstate.h Qoperator.h
	$(CC) -c network.cpp $(LFLAGS) $(DCFLAGS)
