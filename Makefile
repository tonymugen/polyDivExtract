#################################################
#
# Project makefile
#
#################################################

CXX = g++
AXTOBJ = parseAXT.o
VCFOBJ = parseVCF.o
TSTOUT = test
DIVSITES = divSites
CXXFLAGS = -O3 -march=native -std=c++11

all : $(DIVSITES) $(TSTOUT)
.PHONY : all

$(TSTOUT) : test.cpp $(AXTOBJ) $(VCFOBJ)
	$(CXX) test.cpp $(AXTOBJ) $(VCFOBJ) -o $(TSTOUT) $(CXXFLAGS)

$(DIVSITES) : divSites.cpp utilities.hpp $(AXTOBJ)
	$(CXX) divSites.cpp $(AXTOBJ) -o $(DIVSITES) $(CXXFLAGS)

$(AXTOBJ) : parseAXT.cpp parseAXT.hpp
	$(CXX) -c parseAXT.cpp $(CXXFLAGS)

$(VCFOBJ) : parseAXT.cpp parseAXT.hpp parseVCF.cpp parseVCF.hpp
	$(CXX) -c parseVCF.cpp $(CXXFLAGS)

.PHONY : clean
clean:
	-rm *.o $(TSTOUT) $(DIVSITES)

