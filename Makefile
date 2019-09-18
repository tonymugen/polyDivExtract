#################################################
#
# Project makefile
#
#################################################

CXX = g++
AXTOBJ = parseAXT.o
VCFOBJ = parseVCF.o
FFOBJ = ffExtract.o
DIVSITES = divSites
POLYSITES = polySites
SORT = fastaSort
GFFS = getFFsites
CXXFLAGS = -O3 -march=native -std=c++11

all : $(DIVSITES) $(POLYSITES) $(SORT) $(GFFS)
.PHONY : all

$(GFFS) : getFFsites.cpp utilities.hpp $(FFOBJ)
	$(CXX) getFFsites.cpp $(FFOBJ) -o $(GFFS) $(CXXFLAGS)

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

$(FFOBJ) : ffExtract.cpp ffExtract.hpp
	$(CXX) -c ffExtract.cpp $(CXXFLAGS)

.PHONY : clean
clean:
	-rm *.o $(POLYSITES) $(DIVSITES) $(SORT) $(GFFS)

