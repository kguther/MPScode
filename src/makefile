CC=icpc
LFLAGS=-llapacke -llapack -lblas -lgfortran
DCFLAGS=-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP
APCK_INCLUDE=-I /usr/include/arpack++
SUPERLU_INCLUDE=-I /usr/include/superlu
SUPERLU_LIB=/usr/lib/libsuperlu.a
APCK_LIB=/usr/lib/libarpack.a /usr/lib/libarpack++.a
CFLAGS=-Wall -O3
OBJECTS=main.o parameters.o network.o arrayprocessing.o Qsystem.o optHMatrix.o mps.o
TEMPLATES=arraycreation.h mpo.h pContraction.h
NAME=MPScode

all: $(OBJECTS)
	$(CC) -o $(NAME) $(OBJECTS) $(APCK_LIB) $(SUPERLU_LIB) $(LFLAGS) $(DCFLAGS) $(CFLAGS) 

$(OBJECTS): %.o: %.cpp $(TEMPLATES)
	$(CC) -c $(DCFLAGS) $(CFLAGS) $(APCK_INCLUDE) $< -o $@

clean:
	\rm *.o
