EXE=xbadger
CPPMAIN=badger_spatial.cpp
CPPIND=individual.cpp
CPPINDH=individual.hpp

all : $(EXE)

$(EXE) : $(CPPMAIN)  $(CPPIND) $(CPPINDH)
	g++ -ggdb3 -std=c++11 -Wall -O3 -o $(EXE) $(CPPMAIN) $(CPPIND)


.PHONY : 
	clean
	cleanall


clean :
	rm -rf $(EXE)

cleanall :
	rm -rf $(EXE)
	rm -rf sim*

