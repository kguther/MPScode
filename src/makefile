CC=icpc
APCK_LIB=/usr/lib/libarpack.a /usr/lib/libarpack++.a
MKL_PATH=/opt/intel/compilers_and_libaries_2016.0.109/linux/mkl/lib/intel64/
MKL_INCLUDE_PATH=/opt/intel/compilers_and_libaries_2016.0.109/linux/mkl/include
LFLAGS=-L$(MKL_PATH) -lmkl_gf_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -lgfortran
DCFLAGS=-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP
CFLAGS=-Wall -qopenmp -O3 $(MKL_INCLUDE) 
APCK_INCLUDE=-I/usr/include/arpack++
MKL_INCLUDE=-I$(MKL_INCLUDE_PATH)
OBJECTS=network.o arrayprocessing.o Qsystem.o optHMatrix.o mps.o globalMeasurement.o baseMeasurement.o iterativeMeasurement.o overlap.o stateArray.o network_enrichment.o projector.o quantumNumber.o dimensionTable.o basisQNOrderMatrix.o blockHMatrix.o localHSpaces.o localMeasurementSeries.o simulation.o problemOperators.o
TEMPLATES=arraycreation.h mpo.h pContraction.h tmpContainer.h localMpo.h
NAME=vMPScode

all: main.o  $(OBJECTS)
	$(CC) -o $(NAME) $(OBJECTS) main.o $(APCK_LIB) $(CFLAGS) $(LFLAGS)


$(OBJECTS): %.o: %.cpp $(TEMPLATES)
	$(CC) -c $(CFLAGS) $(APCK_INCLUDE) $< -o $@

main.o: main.cpp problemOperators.h $(TEMPLATES)
	$(CC) -c $(CFLAGS) $(APCK_INCLUDE) main.cpp -o main.o

clean:
	\rm *.o
