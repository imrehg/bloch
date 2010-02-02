
CC=g++
CXXFLAGS=-c -g -Wall -O3 -fopenmp -DPARALLEL_MODE_OMP
LDFLAGS=
SOURCES=comb_FAST_SWEEP.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=comb_FAST_SWEEP

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CXXFLAGS) $< -o $@

clean:
	rm -rf *o ${EXECUTABLE}
