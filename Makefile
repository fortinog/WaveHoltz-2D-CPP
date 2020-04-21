# Set direction of source and header files, and build direction
SRC_DIR    = src
INC_DIR    = include
BUILD_DIR  = .build
LDFLAGS    = ""
LDPATH     = SRC_DIR

# list all source and object files
SRC_FILES  = $(wildcard $(SRC_DIR)/*.cpp)
SRC_FILES += $(wildcard $(SRC_DIR)/*/*.cpp)
OBJ_FILES  = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_FILES))
INC        = -I$(INC_DIR)

# CXX        = g++
CXX        = mpicxx
CXXFLAGS   = -O3 

# Rule for linking main
main: $(OBJ_FILES)
	$(CXX) -o $@ $(OBJ_FILES) $(LDFLAGS) -L$(LDPATH)

# Rule for building all src files
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC) 
	@$(CXX) -MM $< -MP -MT $@ -MF $(@:.o=.d) $(INC) 

.PHONY: all clean sandwich para

# use 'make clean' to remove object files and executable
clean:
	rm -fr $(BUILD_DIR) 
	rm -f  main 

debug: CXXFLAGS += -DDEBUG -g
debug: main

# IN CASE OF EMERGENCY
sandwich:
	make clean
	make
	./main
para:
	make clean
	make
	mpirun -n 2 ./main

# include the dependency files
-include $(OBJ_FILES:.o=.d)
