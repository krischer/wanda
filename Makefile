NETCDF_LIB = /apps/monch/netcdf/4.3.1/intel/14.0.1/nonparallel/lib
HDF5_LIB   = /apps/monch/hdf5/1.8.12/intel/14.0.1/lib 

NETCDF_INC = /apps/monch/netcdf/4.3.1/intel/14.0.1/nonparallel/include

CXX      = mpic++ 
CC       = icc
CXXFLAGS = -O3 -std=c++11 -fopenmp -isystem $(NETCDF_INC) -DwithMPI 
CFLAGS   = -O3

LDFLAGS = -fopenmp -L$(NETCDF_LIB) -lnetcdf_c++4 -lnetcdf -L$(HDF5_LIB) -lhdf5_hl -lhdf5 -lz

OBJS_GENERIC     = \
	./main.o \
	./wendy.o \
	./kdtree.o
								 
###################
all: main

main: $(OBJS_GENERIC)
	$(CXX) $(LDFLAGS) -o testWendy $(OBJS_GENERIC) $(LDFLAGS)
	
###################	
clean:
	$(RM) *.o
