#include<stdio.h>
#include "particleType.h"

//ClassImp(particleType)

particleType::particleType(const char *name,double mass,int charge):
  fName(name),
  fMass(mass),
  fCharge(charge)
{
}

void particleType::Print() const{
  printf("%-20s mass=%8.6f        %4i\n",GetParticleName(),GetMass(),GetCharge());
}
