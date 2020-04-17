#
# Makefile for timing examples
#

CXX = mpicxx
LD  = mpicxx
SILO_LIB = /Users/lily/libs/silo/lib
SILO_INC = /Users/lily/libs/silo/include
LDFLAGS = 
CXXINCLUDE = 
CXXFLAGS = -O3 -I$(SILO_INC) -std=c++11
LDFLAGS = -L$(SILO_LIB) -lsilo

OBJECTS = main.o parameters_setup.o Wave_Solve1D.o Darrays.o
EXEC = ex.x

.PHONY: clean

exe: $(EXEC)


$(EXEC): $(OBJECTS)
	$(LD) $(LDFLAGS) $(OBJECTS) -o $(EXEC)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -c $<

clean:
	rm -f $(OBJECTS) $(EXEC)

start:
	rm -f $(OBJECTS) $(EXEC)
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -c $<
	mpirun -np 2 ./$(EXEC)

