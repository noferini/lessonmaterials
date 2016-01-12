#ifndef PARTICLE_H
#define PARTICLE_H

/////////////////////////////////////////////////
//                                             //
//      Particle Class                         //
//                                             //
/////////////////////////////////////////////////


class particleType;  //questa e' una "forward declaration"

class particle //: public TObject
{
 public:
  //membri pubblici
  particle();
 particle(const particle& p);
  virtual ~particle();
  particle(int iparticle,double px=0.0,double py=0.0,double pz=0.0);
  particle(const char *name,double px=0.0,double py=0.0,double pz=0.0);

  static int AddParticleType(const char *name,double mass,int charge,double width=-1.0);
  static void PrintParticleType();

  int GetParticleType() const {return fIparticle;};
  void Print() const;

  void ChangeParticleType(int iparticle);
  void ChangeParticleType(const char *name);

  particle& operator=(const particle &value);
  particle& operator+=(const particle &value);
  particle operator+(const particle &value) const;

  double GetPx()const {return fPx;}
  double GetPy()const {return fPy;}
  double GetPz()const {return fPz;}
  double GetMass() const;
  double GetEnergy() const;

  int Decay2body(particle &dau1,particle &dau2) const;
  static int Decay2body(particle &dau1,particle &dau2,float mass,float px=0,float py=0,float pz=0);
  int Decay3body(particle &dau1,particle &dau2,particle &dau3) const;

  double InvMass(particle & other)const;
  double InvMass(particle & other,particle & other2)const;
  
  void SetP(double px,double py,double pz){fPx=px,fPy=py,fPz=pz;};
  
  static const int fMaxNumParticleType=10; //reso public per funzionare in root
 private:

  void Boost(double bx, double by, double bz);
  static int FindParticle(const char *name);

  static int fNparticleType;
  static particleType *fParticleType[fMaxNumParticleType];

  double fPx,fPy,fPz;  
  int fIparticle;
  
//  ClassDef(particle,1)

};

#endif
