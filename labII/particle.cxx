#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include "particle.h"
#include "particleType.h"
#include "resonanceType.h"

//ClassImp(particle)
int particle::fNparticleType = 0;
particleType *particle::fParticleType[particle::fMaxNumParticleType];

particle::particle(){
  fIparticle =0;
  fPx = 0;
  fPy = 0;
  fPz = 0;
}

particle::particle(const particle& p)//:TObject(p)
{
  fIparticle = p.fIparticle;
  fPx = p.fPx;
  fPy = p.fPy;
  fPz = p.fPz;
}

particle::particle(int iparticle,double px,double py, double pz):
  fPx(px),
  fPy(py),
  fPz(pz)
{
  if(iparticle < fNparticleType && iparticle >=0){
    fIparticle = iparticle;
  }
  else{
    printf("Particle %d doesn't exist in the stack\n",iparticle);
    fIparticle = -1;    
  }
}

particle::particle(const char *name,double px,double py, double pz):
  fPx(px),
  fPy(py),
  fPz(pz)
{
  int ip=FindParticle(name);
  if(ip != -1){
    fIparticle = ip;
  }
  else{
    printf("Particle \"%s\" doesn't exist in the stack\n",name);
    fIparticle = -1;    
  }
}


particle::~particle(){
}

int particle::FindParticle(const char *name){

  for(int i=0;i<fNparticleType;i++){
    const char *currentType = fParticleType[i]->GetParticleName();
    
    int k=0;
    while((currentType[k] == name[k]) && currentType[k] != '\0' && name[k] != '\0'){
      k++;
    }

    if(currentType[k] == name[k]) return i;    
  }

  return -1; // if no match
}

int particle::AddParticleType(const char *name,double mass,int charge,double width){
  if(fNparticleType < fMaxNumParticleType){
    int ip = FindParticle(name);
    if(ip != -1){
      printf("A particle with this name (\"%s\") already exists in the stack (nothing done)\n",name);
      return 2;
    }
    if(width > 0) fParticleType[fNparticleType] = new resonanceType(name,mass,charge,width);
    else fParticleType[fNparticleType] = new particleType(name,mass,charge);

    fNparticleType++;
  }
  else{
    printf("Stack is full because you have already inserted %i particles (nothing done)\n",fMaxNumParticleType);
    return 1;
  }

  return 0;
}

void particle::PrintParticleType(){
  if(fNparticleType){
    printf("%-20s                    %5s\n","Particle Name","charge");
    for(int i=0;i<fNparticleType;i++){
      if(!fParticleType[i]->IsResonance()) 
	fParticleType[i]->Print();
      else
	((resonanceType *) fParticleType[i])->Print();
    }
    printf("\n");
  }
  else
    printf("No particle types defined\n");
}

void particle:: Print() const{
  if(fIparticle!=-1)
    printf("type=%i) %-20s p=(%7.3f,%7.3f,%7.3f)\n",fIparticle,fParticleType[fIparticle]->GetParticleName(),fPx,fPy,fPz);
  else{
    printf("particle not valid!!!        ");
    printf("p=(%7.3f,%7.3f,%7.3f)\n",fPx,fPy,fPz);
  }
}

particle& particle::operator=(const particle &value){
  if (&value==this)return *this; //protezione contro autoassegnamenti
  fPx = value.fPx;
  fPy = value.fPy;
  fPz = value.fPz;
  fIparticle = value.fIparticle;
  return *this;
}

particle& particle::operator+=(const particle & value){
  fPx += value.fPx;
  fPy += value.fPy;
  fPz += value.fPz;

  return *this;
}


particle  particle::operator+(const particle &value) const{
  particle result(*this);  
  result += value;  
  return result;
}

void particle::ChangeParticleType(int iparticle){
  if(iparticle >=0 && iparticle < fNparticleType) fIparticle = iparticle;
  else printf("Particle %d doesn't exist\n",iparticle);
}

void particle::ChangeParticleType(const char *name){
  int ip=FindParticle(name);
  if(ip != -1){
    fIparticle = ip;
  }
  else{
    printf("Particle \"%s\" doesn't exist in the stack\n",name);
  }

}

// extra methods

double particle::GetMass() const {
  if(fIparticle > -1){
    return fParticleType[fIparticle]->GetMass();
  }
  else{
    return 0.0;
  }
}

double particle::GetEnergy() const {
  double mass = GetMass();
  return sqrt(fPx*fPx + fPy*fPy + fPz*fPz + mass*mass); 
}

void particle::Boost(double bx, double by, double bz)
{

  double energy = GetEnergy();

  //Boost this Lorentz vector
  double b2 = bx*bx + by*by + bz*bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx*fPx + by*fPy + bz*fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

  fPx += gamma2*bp*bx + gamma*bx*energy;
  fPy += gamma2*bp*by + gamma*by*energy;
  fPz += gamma2*bp*bz + gamma*bz*energy;
}

int particle::Decay2body(particle &dau1,particle &dau2) const {
  if(GetMass() == 0.0){
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }
  
  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if(fIparticle > -1 && fParticleType[fIparticle]->IsResonance()){ // add width effect
    // gaussian random numbers
    float x1, x2, w, y1;//, y2;
    
    double invnum = 1./RAND_MAX;
    do {
      x1 = 2.0 * rand()*invnum - 1.0;
      x2 = 2.0 * rand()*invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    // y2 = x2 * w;

    massMot += ((resonanceType *) fParticleType[fIparticle])->GetWidth() * y1;

  }

  if(massMot < massDau1 + massDau2){
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }
  
  double pout = sqrt((massMot*massMot - (massDau1+massDau2)*(massDau1+massDau2))*(massMot*massMot - (massDau1-massDau2)*(massDau1-massDau2)))/massMot*0.5;

  double norm = 6.283/RAND_MAX;

  double phi = rand()*norm;
  double theta = rand()*norm*0.5 - 1.57075;
  dau1.SetP(pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta));
  dau2.SetP(-pout*sin(theta)*cos(phi),-pout*sin(theta)*sin(phi),-pout*cos(theta));

  double energy = sqrt(fPx*fPx + fPy*fPy + fPz*fPz + massMot*massMot);

  double bx = fPx/energy;
  double by = fPy/energy;
  double bz = fPz/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);

  return 0;
}


double particle::InvMass(particle & other) const{
  double energy = GetEnergy() + other.GetEnergy();
  double p2 = (fPx+other.GetPx())*(fPx+other.GetPx()) + (fPy+other.GetPy())*(fPy+other.GetPy()) + (fPz+other.GetPz())*(fPz+other.GetPz());
  
  return sqrt(energy*energy - p2);

}

int particle::Decay2body(particle &dau1,particle &dau2,float mass,float px,float py,float pz) {
  if(mass == 0.0){
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }
  
  double massMot = mass;
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if(massMot < massDau1 + massDau2){
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }
  
  double pout = sqrt((massMot*massMot - (massDau1+massDau2)*(massDau1+massDau2))*(massMot*massMot - (massDau1-massDau2)*(massDau1-massDau2)))/massMot*0.5;

  double norm = 6.283/RAND_MAX;

  double phi = rand()*norm;
  double theta = rand()*norm*0.5 - 1.57075;
  dau1.SetP(pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta));
  dau2.SetP(-pout*sin(theta)*cos(phi),-pout*sin(theta)*sin(phi),-pout*cos(theta));

  double energy = sqrt(px*px + py*py + pz*pz + massMot*massMot);

  double bx = px/energy;
  double by = py/energy;
  double bz = pz/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);

  return 0;
}

int particle::Decay3body(particle &dau1,particle &dau2,particle &dau3) const {
  if(GetMass() == 0.0){
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }
  
  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();
  double massDau3 = dau3.GetMass();

  double massLimit2 = (massMot-massDau2)*(massMot-massDau2);

  double invnum = 1./RAND_MAX;

  if(fIparticle > -1 && fParticleType[fIparticle]->IsResonance()){ // add width effect
    // gaussian random numbers
    float x1, x2, w, y1;//, y2;
    do {
      x1 = 2.0 * rand()*invnum - 1.0;
      x2 = 2.0 * rand()*invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    // y2 = x2 * w;

    massMot += ((resonanceType *) fParticleType[fIparticle])->GetWidth() * y1;

  }

  if(massMot < massDau1 + massDau2 + massDau3){
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }

  double xran = rand()*invnum;
  double mass13 = -1;
  int counter = 0;
  while(mass13*mass13/massLimit2 > xran){ // to assure a dalitz homogenous plot
    xran = rand()*invnum;
    double mass12 = (massDau1 + massDau2)*(massDau1 + massDau2);
    mass12 += ((massMot-massDau3)*(massMot-massDau3) - mass12) * xran;
    mass12 = sqrt(mass12);
    
    // perform decay mass3 and mass12 and then mass12 decay
    double pout = sqrt((massMot*massMot - (massDau3+mass12)*(massDau3+mass12))*(massMot*massMot - (massDau3-mass12)*(massDau3-mass12)))/massMot*0.5;
    
    double norm = 6.283/RAND_MAX;
    
    double phi = rand()*norm;
    double theta = rand()*norm*0.5 - 1.57075;
    dau3.SetP(pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta));
    particle::Decay2body(dau1,dau2,mass12,pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta));
    
    
    mass13 = dau3.InvMass(dau1);
    xran = rand()*invnum;
    counter ++;
    if(counter > 20) printf("Some problems in performing decay -> counter = %i (m12 =%f, m13=%f)\n",counter,mass12,mass13);
  }

  double energy = sqrt(fPx*fPx + fPy*fPy + fPz*fPz + massMot*massMot);

  double bx = fPx/energy;
  double by = fPy/energy;
  double bz = fPz/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);
  dau3.Boost(bx,by,bz);

  return 0;
}
