CC = g++
FLAGS = -O3 -std=c++17 -fopenmp
DEBUG_FLAGS = -g -std=c++17
INCLUDE = -I ../include
MAIN = main.cpp
EXE = TICG-engine

# Python configuration - uses current environment's python
PYTHON ?= python
PYBIND_INCLUDES = $(shell $(PYTHON) -m pybind11 --includes)
PYTHON_EXT_SUFFIX = $(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")

# macOS requires -undefined dynamic_lookup for Python extensions
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    PYBIND_LDFLAGS = -undefined dynamic_lookup
else
    PYBIND_LDFLAGS =
endif

all: $(EXE)

random_mars.o: random_mars.cpp random_mars.h
	$(CC) $(FLAGS) -c $<

Grid.o: Grid.cpp Grid.h
	$(CC) $(FLAGS) $(INCLUDE) -c $<

Sim.o: Sim.cpp Sim.h
	$(CC) $(FLAGS) $(INCLUDE) -c $<

$(EXE): $(MAIN) Grid.o Sim.o random_mars.o
	$(CC) $(FLAGS) $(INCLUDE) -o $@ $^

debug2: $(MAIN) Grid.o Sim.o random_mars.o
	$(CC) $(DEBUG_FLAGS) $(INCLUDE) -o TICG-engine $^

debug:
	g++ $(DEBUG_FLAGS) -I ../include -c Grid.cpp
	g++ $(DEBUG_FLAGS) -I ../include -c Sim.cpp
	g++ $(DEBUG_FLAGS) -c random_mars.cpp
	g++ $(DEBUG_FLAGS) -I ../include -o TICG-debug $(MAIN) Grid.o Sim.o random_mars.o

pybind:
	$(CC) -O3 -shared -std=c++17 $(INCLUDE) -fPIC $(PYBIND_INCLUDES) $(PYBIND_LDFLAGS) pybind_Sim.cpp -o scribe_engine$(PYTHON_EXT_SUFFIX)

clean:
	rm -f *.o

get_contacts: Grid.o Sim.o random_mars.o
	$(CC) $(FLAGS) $(INCLUDE) get_contacts.cpp -o $@ $^
