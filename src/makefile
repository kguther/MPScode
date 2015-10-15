CC=icpc
LFLAGS=-llapacke -llapack -lblas -lgfortran
DCFLAGS=-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP
CFLAGS=-Wall -g
OBJECTS=main.o parameters.o network.o arrayprocessing.o
NAME=MPScode

all: $(OBJECTS)
	$(CC) -o $(NAME) $(OBJECTS) $(LFLAGS) $(DCFLAGS) $(CFLAGS)

$(OBJECTS): %.o: %.cpp
	$(CC) -c $(DCFLAGS) $(CFLAGS) $< -o $@

clean:
	\rm *.o
