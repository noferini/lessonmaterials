#!/bin/bash
export LIBS=$(root-config --glibs)
export CC=$(root-config --cflags)
export LD=-Wl,--no-as-needed
c++ -I../ $LD $CC $LIBS analyze3D.C ../particle.cxx ../particleType.cxx ../resonanceType.cxx -o ana.exe
