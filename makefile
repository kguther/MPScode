CC=icpc
LFLAGS=-llapacke -llapack -lblas -lgfortran
DCFLAGS=-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

all: main.o parameters.o
	$(CC) -o MPScode main.o parameters.o $(LFLAGS) $(DCFLAGS) -Wall

main.o: main.cpp
	$(CC) -c main.cpp

parameters.o: parameters.cpp
	$(CC) -c parameters.cpp
