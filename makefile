# **************************************************
# Variables to control Makefile operation
CXX= -g++ -fopenmp -static-libstdc++
CXXFLAGS= -Wall -g -O3
#*****************************************************
# Targets needed to bring the executable up to date
all: program 

program:folder main.o TissueBacteria.o Bacteria.o Nodes.o SignalCommon.o ConfigParser.o commonData.o GeoVector.o ranum2.o Diffusion2D.o Grid.o TissueGrid.o driver.o Fungi.o growthFunctions.o hyphaeSegment.o constants.o 
	$(CXX) $(CXXFLAGS) -o program main.o TissueBacteria.o Bacteria.o Nodes.o SignalCommon.o ConfigParser.o commonData.o GeoVector.o ranum2.o Diffusion2D.o Grid.o TissueGrid.o driver.o Fungi.o growthFunctions.o hyphaeSegment.o constants.o 
folder: 
		mkdir -p ./DataOutput ./Animation
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

Bacteria.o: Bacteria.cpp
	$(CXX) $(CXXFLAGS) -c Bacteria.cpp

TissueBacteria.o: TissueBacteria.cpp
	$(CXX) $(CXXFLAGS) -c TissueBacteria.cpp

Nodes.o: Nodes.cpp
	$(CXX) $(CXXFLAGS) -c Nodes.cpp

ConfigParser.o: ConfigParser.cpp
	$(CXX) $(CXXFLAGS) -c ConfigParser.cpp
	
commonData.o: commonData.cpp
	$(CXX) $(CXXFLAGS) -c commonData.cpp

SignalCommon.o: SignalCommon.cpp
	$(CXX) $(CXXFLAGS) -c SignalCommon.cpp
	
GeoVector.o: GeoVector.cpp
	$(CXX) $(CXXFLAGS) -c GeoVector.cpp

ranum2.o: ranum2.cpp
	$(CXX) $(CXXFLAGS) -c ranum2.cpp

Diffusion2D.o: Diffusion2D.cpp
	$(CXX) $(CXXFLAGS) -c Diffusion2D.cpp
	
Grid.o: Grid.cpp
	$(CXX) $(CXXFLAGS) -c Grid.cpp

TissueGrid.o: TissueGrid.cpp
	$(CXX) $(CXXFLAGS) -c TissueGrid.cpp

driver.o: driver.cpp
	$(CXX) $(CXXFLAGS) -c driver.cpp

Fungi.o: Fungi.cpp
	$(CXX) $(CXXFLAGS) -c Fungi.cpp

growthFunctions.o: growthFunctions.cpp
	$(CXX) $(CXXFLAGS) -c growthFunctions.cpp

hyphaeSegment.o: hyphaeSegment.cpp
	$(CXX) $(CXXFLAGS) -c hyphaeSegment.cpp

constants.o: constants.cpp
	$(CXX) $(CXXFLAGS) -c constants.cpp

clean: wipe
		rm -rf strain_vec.txt stress_vec.txt *o program

wipe:
