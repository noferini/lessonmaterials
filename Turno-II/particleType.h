#ifndef PARTICLETYPE_H
#define PARTICLETYPE_H

/////////////////////////////////////////////////
//                                             //
//      Particle Type Class                    //
//                                             //
/////////////////////////////////////////////////


class particleType
{
 public:
  particleType():fName(" "), fMass(0.), fCharge(0){};
  particleType(const char *name,double mass,int charge);
  virtual ~particleType(){};

  const char *GetParticleName() const {return fName;};
  double GetMass() const {return fMass;};
  int GetCharge() const {return fCharge;};

  virtual int  IsResonance() const {return int(0);};

  virtual void Print() const;
  
 private:
  const char *fName;
  const double fMass;
  const int fCharge;

  // ClassDef(particleType,1) 
};
#endif
