#ifndef RESONANCETYPE_H
#define RESONANCETYPE_H

/////////////////////////////////////////////////
//                                             //
//      Resonance Type Class                   //
//                                             //
/////////////////////////////////////////////////


#include "particleType.h"

class resonanceType : public particleType
{
 public:
  resonanceType(const char *name,double mass,int charge,double width=0.0);
  ~resonanceType(){};

  int  IsResonance() const {return int(1);};

  double GetWidth() const {return fWidth;};

  void Print() const;


 private:
  const double fWidth;

  //  ClassDef(resonanceType,1) 

};
#endif
