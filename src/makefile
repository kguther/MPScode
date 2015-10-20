CC=icpc
LFLAGS=-llapacke -llapack -lblas -lgfortran
DCFLAGS=-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP
APCK_INCLUDE=-I /usr/include/arpack++
SUPERLU_INCLUDE=-I /usr/include/superlu
SUPERLU_LIB=/usr/lib/libsuperlu.a
APCK_LIB=/usr/lib/libarpack.a /usr/lib/libarpack++.a
CFLAGS=-Wall -g
OBJECTS=main.o parameters.o network.o arrayprocessing.o
NAME=MPScode

all: $(OBJECTS)
	$(CC) -o $(NAME) $(OBJECTS) $(SUPERLU_LIB) $(APCK_LIB) $(LFLAGS) $(DCFLAGS) $(CFLAGS) 

$(OBJECTS): %.o: %.cpp
	$(CC) -c $(DCFLAGS) $(CFLAGS) $(SUPERLU_INCLUDE) $(APCK_INCLUDE) $< -o $@

clean:
	\rm *.o
