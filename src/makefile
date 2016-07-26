CC=icpc
APCK_LIB=/usr/lib/libarpack.a /usr/lib/libarpack++.a
MKL_PATH=/opt/intel/compilers_and_libraries_2016.0.109/linux/mkl/lib/intel64/
MKL_INCLUDE_PATH=/opt/intel/compilers_and_libraries_2016.0.109/linux/mkl/include
LPKFLAGS=-llapacke -llapack -lblas
MKLFLAGS=-lmkl_gf_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm
LBFLAGS=$(MKLFLAGS)
DMKL=-DUSE_MKL
NAME=scMPScode
DLPCK=-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP
HOLLE_NAME=holleMPScode
LFLAGS=$(MPI_LINK_FLAGS) -L$(MKL_PATH) -lgfortran
MPI_COMPILE_FLAGS=$(shell mpic++ --showme:compile)
MPI_LINK_FLAGS=$(shell mpic++ --showme:link)
CFLAGS=-Wall -std=c++11 -fopenmp -O3 $(MKL_INCLUDE)
CFLAGS_DEBUG=-Wall -std=c++11 -fopenmp -g $(MKL_INCLUDE)
APCK_INCLUDE=-I/usr/include/arpack++
MKL_INCLUDE=-I$(MKL_INCLUDE_PATH)
DEBUG_NAME=dMPScode
OBJECTS=network.o arrayprocessing.o Qsystem.o optHMatrix.o mps.o globalMeasurement.o baseMeasurement.o iterativeMeasurement.o overlap.o stateArray.o network_enrichment.o projector.o quantumNumber.o dimensionTable.o basisQNOrderMatrix.o blockHMatrix.o localHSpaces.o localMeasurementSeries.o simulation.o problemOperators.o interface.o main.o parameters.o exactGroundState.o verifyQN.o pseudoQuantumNumber.o siteQNOrderMatrix.o truncation.o heisenbergChain.o
OBJECTS_DEBUG=$(OBJECTS:.o=.dbg)
SOURCE=$(OBJECTS:.o=.cpp)
OBJECTS_HOLLE=$(OBJECTS:.o=.hll)
DEPENDENCIES=.depend

-include $(DEPENDENCIES)

fuchur: $(OBJECTS)
	$(CC) -o $(NAME) $(OBJECTS) $(APCK_LIB) $(CFLAGS) $(LFLAGS) $(MKLFLAGS) $(DMKL)

debug: $(OBJECTS_DEBUG)
	$(CC) -o $(DEBUG_NAME) $(OBJECTS_DEBUG) $(APCK_LIB) $(CFLAGS_DEBUG) $(LFLAGS) $(LPKFLAGS) $(DLPCK)

holle: $(OBJECTS_HOLLE)
	$(CC) -o $(HOLLE_NAME) $(OBJECTS_HOLLE) $(APCK_LIB) $(CFLAGS) $(LFLAGS) $(LPKFLAGS) $(DLPCK)

dep: $(SOURCE)
	$(CC) $(APCK_INCLUDE) $(MPI_COMPILE_FLAGS) -MM $(SOURCE) > $(DEPENDENCIES)

$(OBJECTS): %.o: %.cpp 
	$(CC) -c $(CFLAGS) $(APCK_INCLUDE) $(MPI_COMPILE_FLAGS) $(DMKL) $< -o $@

$(OBJECTS_DEBUG): %.dbg: %.cpp 
	$(CC) -c $(CFLAGS_DEBUG) $(APCK_INCLUDE) $(MPI_COMPILE_FLAGS) $(DLPCK) $< -o $@

$(OBJECTS_HOLLE): %.hll: %.cpp
	$(CC) -c $(CFLAGS) $(APCK_INCLUDE) $(MPI_COMPILE_FLAGS) $(DLPCK) $< -o $@

clean:
	\rm *.o
	\rm *.hll
	\rm *.dbg
