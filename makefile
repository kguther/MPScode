CC=icpc
APCK_LIB=/usr/lib/libarpack.a /usr/lib/x86_64-linux-gnu/libarpack++.a
MKL_PATH=/opt/intel/compilers_and_libraries_2016.0.109/linux/mkl/lib/intel64/
MKL_INCLUDE_PATH=/opt/intel/compilers_and_libraries_2016.0.109/linux/mkl/include
INCL_SRC_PATH=-I src/headers -I src/headers/systems
LPKFLAGS=-llapacke -llapack -lblas -liomp5 -lpthread
MKLFLAGS=-lmkl_gf_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm
LBFLAGS=$(MKLFLAGS)
DMKL=-DUSE_MKL -DREAL_MPS_ENTRIES
NAME=scMPScode
DLPCK=-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP
DLPCK=-DREAL_MPS_ENTRIES
HOLLE_NAME=holleMPScode
LFLAGS=$(MPI_LINK_FLAGS) -L$(MKL_PATH) -lgfortran
MPI_COMPILE_FLAGS=$(shell mpic++ --showme:compile)
MPI_LINK_FLAGS=$(shell mpic++ --showme:link)
CFLAGS=-Wall -std=c++11 -qopenmp -O3
CFLAGS_DEBUG=-Wall -std=c++11 -qopenmp -g
APCK_INCLUDE=-I/usr/include/arpack++
MKL_INCLUDE=-I$(MKL_INCLUDE_PATH)
DEBUG_NAME=dMPScode
OBJECTS_LIB=src/parameters/parameters.o src/network/network.o src/auxiliary/arrayprocessing.o src/network/optimizer/optHMatrix.o src/network/mps.o src/network/measurement/globalMeasurement.o src/network/measurement/baseMeasurement.o src/network/measurement/iterativeMeasurement.o src/network/overlap.o src/network/stateArray.o src/QN/network_enrichment.o src/network/optimizer/projector.o src/QN/quantumNumber.o src/parameters/dimensionTable.o src/QN/basisQNOrderMatrix.o src/network/optimizer/blockHMatrix.o src/parameters/localHSpaces.o src/network/measurement/localMeasurementSeries.o src/QN/exactGroundState.o src/QN/verifyQN.o src/QN/pseudoQuantumNumber.o src/QN/siteQNOrderMatrix.o src/QN/truncation.o src/auxiliary/linalgWrapper.o
OBJECTS=$(OBJECTS_LIB) src/frontend/Qsystem.o src/frontend/simulation.o src/frontend/problemOperators.o src/frontend/interface.o src/frontend/main.o src/frontend/heisenbergChain.o
OBJECTS_DEBUG=$(OBJECTS:.o=.dbg)
OBJECTS_LIB_HOLLE=$(OBJECTS_LIB:.o=.hll)
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
	$(CC) $(APCK_INCLUDE) $(MPI_COMPILE_FLAGS) $(INCL_SRC_PATH) -MM -std=c++11  $(SOURCE) > $(DEPENDENCIES)

library: $(OBJECTS_LIB_HOLLE)
	ar rcs src/lib/libvmps.a $(OBJECTS_LIB_HOLLE)

$(OBJECTS): %.o: %.cpp 
	$(CC) -c $(CFLAGS) $(APCK_INCLUDE) $(MPI_COMPILE_FLAGS) $(DMKL) $(INCL_SRC_PATH) $(MKL_INCLUDE) $< -o $@

$(OBJECTS_DEBUG): %.dbg: %.cpp 
	$(CC) -c $(CFLAGS_DEBUG) $(APCK_INCLUDE) $(MPI_COMPILE_FLAGS) $(DLPCK) $(INCL_SRC_PATH) $< -o $@

$(OBJECTS_HOLLE): %.hll: %.cpp
	$(CC) -c $(CFLAGS) $(APCK_INCLUDE) $(DLPCK) $(MPI_COMPILE_FLAGS) $(INCL_SRC_PATH) $< -o $@

cleanHolle:
	rm $(OBJECTS_HOLLE)
cleanFuchur:
	rm $(OBJECTS)
cleanDebug:
	rm $(OBJECTS_DEBUG)
