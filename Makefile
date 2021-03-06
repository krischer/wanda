NETCDF_LIB = /apps/monch/netcdf/4.3.1/intel/14.0.1/nonparallel/lib
EXODUS_LIB = /apps/monch/exodus/6.02/intel/14.0.1/nonparallel/lib
HDF5_LIB   = /apps/monch/hdf5/1.8.12/intel/14.0.1/lib 

EXODUS_INC = /apps/monch/exodus/6.02/intel/14.0.1/nonparallel/include
NETCDF_INC = /apps/monch/netcdf/4.3.1/intel/14.0.1/nonparallel/include

CXX      = mpic++ 
CC       = icc
CXXFLAGS = -O3 -std=c++11 -isystem$(NETCDF_INC) -isystem$(EXODUS_INC) -DwithMPI -shared
CFLAGS   = -O3

LDFLAGS = -L$(EXODUS_LIB) -lexodus -L$(NETCDF_LIB) -lnetcdf_c++4 -lnetcdf -L$(HDF5_LIB) -lhdf5_hl -lhdf5 -lz

OBJS_GENERIC     = \
	./src/main.o \
	./src/wendy.o \
	./src/kdtree.o
								 
./src/%.o: ./src/%.c
	@$(CC) 						$(CFLAGS) -c -o $@ $<
	@echo "CC  $<"

./src/%.o: ./src/%.cpp
	@$(CXX) $(CXXFLAGS) -c -o $@ $<
	@echo "CXX $<"

###################
all: main

main: $(OBJS_GENERIC)
	@$(CXX) $(LDFLAGS) -o testWendy $(OBJS_GENERIC) $(LDFLAGS)
	@echo "Linking ... $<"
	
###################	
clean:
	$(RM) ./src/*.o ./testWendy
