#!/bin/bash
export LIBS=$(root-config --glibs)
export CC=$(root-config --cflags)
c++ $CC $LIBS sim.C particle.cxx particleType.cxx resonanceType.cxx -o sim.exe
