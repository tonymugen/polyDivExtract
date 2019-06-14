#################################################
#
# Project makefile
#
#################################################

CXX = g++
AXTOBJ = parseAXT.o
TSTOUT = test
DIVSITES = divSites
CXXFLAGS = -O3 -march=native -std=c++11

all : $(TSTOUT) $(DIVSITES)
.PHONY : all

$(TSTOUT) : test.cpp $(AXTOBJ)
	$(CXX) test.cpp $(AXTOBJ) -o $(TSTOUT) $(CXXFLAGS)

$(DIVSITES) : divSites.cpp utilities.hpp $(AXTOBJ)
	$(CXX) divSites.cpp $(AXTOBJ) -o $(DIVSITES) $(CXXFLAGS)

$(AXTOBJ) : parseAXT.cpp parseAXT.hpp
	$(CXX) -c parseAXT.cpp $(CXXFLAGS)

.PHONY : clean
clean:
	-rm *.o $(TSTOUT) $(DIVSITES)

