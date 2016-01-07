#include<stdio.h>
#include "resonanceType.h"

//ClassImp(resonanceType)

resonanceType::resonanceType(const char *name,double mass,int charge,double width):
  particleType(name,mass,charge),
  fWidth(width)
{
}

void resonanceType::Print() const{
  printf("%-20s mass=%8.6f        %4i     (width=%5.2f)\n",GetParticleName(),GetMass(),GetCharge(),GetWidth());
}
