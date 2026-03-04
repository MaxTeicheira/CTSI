#!/bin/bash

#compile all the .c file in use to .o object file.

g++ -c src/Detector.cpp
g++ -c src/CTSI.cpp
g++ -c src/main.cpp
g++ -c src/Event.cpp
g++ -c src/MyRandom.cpp
g++ -c src/Preamplification.cpp
g++ -c src/Simulation.cpp

g++ -o ctsi.exe main.o CTSI.o Detector.o Event.o MyRandom.o Preamplification.o Simulation.o

chmod a+x ctsi.exe

