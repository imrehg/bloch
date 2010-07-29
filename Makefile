CC=g++
CXXFLAGS=-O3 -llapack -lblas -fopenmp 
LATLAS=/usr/lib/atlas/liblapack.a /usr/lib/atlas/libblas.a -latlas
LDFLAGS=-o
MAIN=comb_SWEEP_main.cpp
D1=comb_SWEEP.cpp
MAIN_CW=CW_SWEEP_main.cpp
CW=CW_SWEEP.cpp
O=$(MAIN:.cpp=.o)
O2=$(D1:.cpp=.o)
O3=$(MAIN_CW:.cpp=.o)
O4=$(CW:.cpp=.o)
EXECUTABLE=comb
EXECUTABLE_CW=CW

comb: $(O) $(O2) final
	
$(O): atom.cpp comb.h $(MAIN)
	$(CC) -c $(MAIN) $(CXXFLAGS)

$(O2): atom.cpp comb.h $(D1)
	$(CC) -c $(D1) $(CXXFLAGS)
final:
	$(CC) $(LDFLAGS) $(EXECUTABLE) $(O) $(O2) $(CXXFLAGS)

CW:  $(O3) $(O4) final_CW

$(O3): atom.cpp CW.h $(MAIN_CW)
	$(CC) -c $(MAIN_CW) $(CXXFLAGS)

$(O4): atom.cpp CW.h $(CW)
	$(CC) -c $(CW) $(CXXFLAGS)

final_CW:
	$(CC) $(LDFLAGS) $(EXECUTABLE_CW) $(O3) $(O4) $(CXXFLAGS)


clean:
	rm *.o
