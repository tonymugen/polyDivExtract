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
SORT = fastaSort
CXXFLAGS = -O3 -march=native -std=c++11

all : $(DIVSITES) $(POLYSITES) $(SORT)
.PHONY : all

$(SORT) : fastaSort.cpp utilities.hpp
	$(CXX) fastaSort.cpp -o $(SORT) $(CXXFLAGS)

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
	-rm *.o $(POLYSITES) $(DIVSITES) $(SORT)

