CXX = icpc
CXXFLAGS = -c -I/home/anandps/local/adolc/include -I/home/anandps/local/petsc_3.4.5/include -I/home/anandps/local/mpich/include  -std=c++11 #-debug -g -traceback
#CXXFLAGS = -c -I$(HOME)/local/adolc/include -g3 -std=c++11

#CXXFLAGS = -c -I/home/anandps/local/adolc/include -I/home/anandps/local/petsc_3.4.5/include -I/home/anandps/local/mpich/include  #-debug

OBJECTS = petsc.o residual.o utils.o io.o metrics.o main.o

BIN = bin/
SRC = src/


single:$(OBJECTS)
	$(CXX) -o mujo $(OBJECTS) -L$(HOME)/local/adolc/lib64 -ladolc -L$(HOME)/local/petsc_3.4.5/lib -lpetsc -L/home/anandps/local/mpich/lib -lmpich
	mv mujo $(BIN).

%.o:$(SRC)%.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm *.o bin/mujo

#all: pre-build main-build
#pre-build:
#	swig -c++ -python api.i
#main-build: $(OBJ)
#	$(CXX) $(CFLAGS) $(OBJ)
#	icpc -shared *.o -o _api.so -openmp -L/home/anandps/local/adolc/lib64 -ladolc
#clean:
#	rm -f core *.o *~ *.so *.cxx
