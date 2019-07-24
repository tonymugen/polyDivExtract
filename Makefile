#################################################
#
# Project makefile
#
#################################################

CXX = g++
AXTOBJ = parseAXT.o
VCFOBJ = parseVCF.o
DIVSITES = divSites
POLYSITES = polySites
CXXFLAGS = -O3 -march=native -std=c++11

all : $(DIVSITES) $(POLYSITES)
.PHONY : all

$(POLYSITES) : polySites.cpp utilities.hpp $(AXTOBJ) $(VCFOBJ)
	$(CXX) polySites.cpp $(AXTOBJ) $(VCFOBJ) -o $(POLYSITES) $(CXXFLAGS)

$(DIVSITES) : divSites.cpp utilities.hpp $(AXTOBJ)
	$(CXX) divSites.cpp $(AXTOBJ) -o $(DIVSITES) $(CXXFLAGS)

$(AXTOBJ) : parseAXT.cpp parseAXT.hpp
	$(CXX) -c parseAXT.cpp $(CXXFLAGS)

$(VCFOBJ) : parseAXT.cpp parseAXT.hpp parseVCF.cpp parseVCF.hpp
	$(CXX) -c parseVCF.cpp $(CXXFLAGS)

.PHONY : clean
clean:
	-rm *.o $(POLYSITES) $(DIVSITES)

