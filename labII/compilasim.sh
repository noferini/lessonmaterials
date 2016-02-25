#!/bin/bash
export LIBS=$(root-config --glibs)
export CC=$(root-config --cflags)
export LD=-Wl,--no-as-needed
c++ $LD $CC $LIBS sim.C particle.cxx particleType.cxx resonanceType.cxx -o sim.exe
