#!/bin/bash
c++ -I$ROOTSYS/include -L$ROOTSYS/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic analyze.C ../particle.cxx ../particleType.cxx ../resonanceType.cxx -o ana.exe
