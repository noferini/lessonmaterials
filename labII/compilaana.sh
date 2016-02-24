#!/bin/bash
export LIBS=$(root-config --glibs)
export CC=$(root-config --cflags)
c++ $CC $LIBS analyze.C particle.cxx particleType.cxx resonanceType.cxx -o ana.exe
