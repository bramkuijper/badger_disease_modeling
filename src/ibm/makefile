EXE=xbadger
CPPMAIN=badger_spatial.cpp
CPPIND=individual.cpp
CPPINDH=individual.hpp

all : $(EXE)

$(EXE) : $(CPPMAIN)  $(CPPIND) $(CPPINDH)
	g++ -Wall -O3 -o $(EXE) $(CPPMAIN) $(CPPIND)


.PHONY : 
	clean


clean :
	rm -rf $(EXE)

