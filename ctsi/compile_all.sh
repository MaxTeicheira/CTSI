#!/bin/bash

#compile all the .c file in use to .o object file.

g++ -c Detector.cpp
g++ -c CTSI.cpp
g++ -c main.cpp
g++ -c Event.cpp
g++ -c MyRandom.cpp
g++ -c Preamplification.cpp
g++ -c Simulation.cpp

g++ -o ctsi.exe main.o CTSI.o Detector.o Event.o MyRandom.o Preamplification.o Simulation.o

chmod a+x ctsi.exe

