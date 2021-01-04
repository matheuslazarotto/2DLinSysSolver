# Folder with '.cpp' and '.hpp' files (definitions)
VPATH=./src

# Compilation flags
CXXFLAGS=-Wall -I$(VPATH) #-Werror -O3

# Linked libraries
LINKFLAGS = -lm -lgsl -lgslcblas

# Files to be included
DEPENDENCIES=main.cpp vector_aux.cpp ode_sys_solv.cpp auxf.cpp

# Executable
TARGETS=2dlinsyssolver

build: $(TARGETS)

$(TARGETS): $(DEPENDENCIES)
	$(CXX) $(CXXFLAGS) -o $(TARGETS) $^ $(LINKFLAGS)

.PHONY: clean
clean:
	-rm -f $(TARGETS)
