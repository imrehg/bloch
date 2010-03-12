CC=g++
CXXFLAGS=-O3 -llapack -lblas -fopenmp 
LATLAS=/usr/lib/atlas/liblapack.a /usr/lib/atlas/libblas.a -latlas
LDFLAGS=-o
MAIN=comb_SWEEP_D1_main.cpp
D1=comb_SWEEP_D1.cpp
O=$(MAIN:.cpp=.o)
O2=$(D1:.cpp=.o)
EXECUTABLE=comb

all: $(O) $(O2) final
	
$(O): atom.cpp comb.h $(MAIN)
	$(CC) -c $(MAIN) $(CXXFLAGS)

$(O2): atom.cpp comb.h $(D1)
	$(CC) -c $(D1) $(CXXFLAGS)
final:
	$(CC) $(LDFLAGS) $(EXECUTABLE) $(O) $(O2) $(CXXFLAGS)
