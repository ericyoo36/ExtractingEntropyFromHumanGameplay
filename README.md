# CPSC530-CounterStrikeEntropy
Project for CPSC530 to generate entropy from Counter-Strike: Global Offensive.  The mouse data we used exists in newcoordinates1.txt, newcoordinates2.txt, newcoordinates3.txt, and newcoordinatesConcat.txt, the last of which represents the three coordinates files concatenated.

# Instructions to run
Compile with regular java compilation, make sure to include the .jar files from our repo in the compilation process.

To run our program, execute EntropyProcessor with commandline arguments as such: "EntropyProcessor [input filename] [confidence value 0.0 - 1.0]"

To run the ahk script, download AutoHotKey and compile using AHK, then run the .exe file it generates and press f12 to start and stop recording coordinates.  It will save a file called "coordinates.txt" in whatever folder the executable is in once you stop recording.
