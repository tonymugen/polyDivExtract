#################################################
#
# Project makefile
#
#################################################

CXX = g++
AXTOBJ = parseAXT.o
TSTOUT = test
CXXFLAGS = -O3 -march=native -std=c++11

$(TSTOUT) : test.cpp $(AXTOBJ)
	$(CXX) test.cpp $(AXTOBJ) -o $(TSTOUT) $(CXXFLAGS)

$(AXTOBJ) : parseAXT.cpp parseAXT.hpp
	$(CXX) -c parseAXT.cpp $(CXXFLAGS)

.PHONY : clean
clean:
	-rm *.o $(TSTOUT)

