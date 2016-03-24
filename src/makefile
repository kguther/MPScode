CC=icpc
APCK_LIB=/usr/lib/libarpack.a /usr/lib/libarpack++.a
MKL_PATH=/opt/intel/compilers_and_libraries_2016.0.109/linux/mkl/lib/intel64/
MKL_INCLUDE_PATH=/opt/intel/compilers_and_libraries_2016.0.109/linux/mkl/include
LPKFLAGS=-llapacke -llapack -lblas
MKLFLAGS=-lmkl_gf_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm
USE_MKL=1
ifdef USE_MKL
LBFLAGS=$(MKLFLAGS)
DMKL=-DUSE_MKL
NAME=scMPScode
endif
ifndef USE_MKL
LBFLAGS=$(LPKFLAGS)
DMKL=-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP
NAME=holleMPScode
endif
LFLAGS=$(MPI_LINK_FLAGS) -L$(MKL_PATH) $(LBFLAGS) -lgfortran
DCFLAGS=$(DMKL)
MPI_COMPILE_FLAGS=$(shell mpic++ --showme:compile)
MPI_LINK_FLAGS=$(shell mpic++ --showme:link)
CFLAGS=-Wall -fopenmp -O3 $(MKL_INCLUDE)
APCK_INCLUDE=-I/usr/include/arpack++
MKL_INCLUDE=-I$(MKL_INCLUDE_PATH)
OBJECTS=network.o arrayprocessing.o Qsystem.o optHMatrix.o mps.o globalMeasurement.o baseMeasurement.o iterativeMeasurement.o overlap.o stateArray.o network_enrichment.o projector.o quantumNumber.o dimensionTable.o basisQNOrderMatrix.o blockHMatrix.o localHSpaces.o localMeasurementSeries.o simulation.o problemOperators.o interface.o main.o parameters.o exactGroundState.o baseTensor.o
TEMPLATES=arraycreation.h mpo.h pContraction.h tmpContainer.h localMpo.h

all: $(OBJECTS)
	$(CC) -o $(NAME) $(OBJECTS) $(APCK_LIB) $(CFLAGS) $(LFLAGS) $(DCFLAGS)

$(OBJECTS): %.o: %.cpp $(TEMPLATES)
	$(CC) -c $(CFLAGS) $(APCK_INCLUDE) $(MPI_COMPILE_FLAGS) $(DCFLAGS) $< -o $@

clean:
	\rm *.o
