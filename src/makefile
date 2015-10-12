CC=icpc
LFLAGS=-llapacke -llapack -lblas -lgfortran
DCFLAGS=-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

all: main.o parameters.o
	$(CC) -o MPScode main.cpp parameters.cpp $(LFLAGS) $(DCFLAGS) -Wall -g

main.o: main.cpp
	$(CC) -c main.cpp $(LFLAGS) $(DCFLAGS)

parameters.o: parameters.cpp
	$(CC) -c parameters.cpp
