CC=g++
CXXFLAGS=-O3 -fopenmp
LDFLAGS=-fopenmp
MAIN=comb_SWEEP_D1_main.cpp
D1=comb_SWEEP_D1.cpp
O=$(MAIN:.cpp=.o)
O2=$(D1:.cpp=.o)
EXECUTABLE=comb

all: $(O) $(O2) final
	
$(O): atom.cpp comb.h $(MAIN)
	$(CC) -c $(CXXFLAGS) $(MAIN)

$(O2): atom.cpp comb.h $(D1)
	$(CC) -c $(CXXFLAGS) $(D1)
final:
	$(CC) $(LDFLAGS) -o $(EXECUTABLE) $(O) $(O2) 
