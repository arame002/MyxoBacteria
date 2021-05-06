# -*- Makefile -*-

# target: dependencies
#	action

driver: driver.o growthFunctions.o constants.o
	g++ -g driver.o growthFunctions.o constants.o -o driver -l m 

driver.o: driver.cpp growthFunctions.h hyphaeSegment.h constants.h
	g++ -g -c driver.cpp

growthFunctions.o: growthFunctions.cpp hyphaeSegment.h
	g++ -g -c growthFunctions.cpp

constants.o: constants.cpp
	g++ -g -c constants.cpp

clean: 
	rm -f *.o driver