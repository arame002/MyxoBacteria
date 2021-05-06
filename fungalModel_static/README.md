# Fungal Growth Model - Static & Simple

## Description of .cpp Files: 
### driver.cpp
- *Purpose:* 
  - Main script used to run fungal growth simulation. 
  - Will generate a very simple mycelia structure composed of 8 hyphae segments.
  - The information for each segment is stored in a vector of class instances.
- *Output:* When the file is ran, it will print out the endpoints and angles of each segment.

Important:
- The i^th segment is considered a tip if
`hy[i].can_extend == true`
- The location of the tip is given by
`hy[i].x2`
&
`hy[i].y2`

### constants.cpp
- *Purpose:* Contains the fixed parameters used throughout the code.
- *Note:* In this case, there are only 3 parameters listed, but in the dynamic version of the code there are many.


### growthFunctions.cpp
- *Purpose:* Contains functions that are used to help the fungal structure grow.
- *Note:* In this case, there is only one function, but in the dynamic version of the code there are many.
  - The function here (`add_new_hyphae`) generates the geometry for the newly formed segments and classifies the relationship with neighboring segments.

## Description of .h Files:
### constants.h
Header file for constants.cpp

### growthFunctions.h
Header file for growthFunctions.cpp

### hyphaeSegment.h
Defines the class structure for each hyphae segment:
- starting endpoints of segment (`x1` & `y1`)
- ending endpoints of segment (`x2` & `y2`)
- angle of segment (`angle`)
- whether or not a segment can branch and extend (`can_branch` & `can_extend`)
- neighbor information: who it extends/branched from (`from_who`), who it extends to (`extend_to`), and who it branches to (`branch_to`)



