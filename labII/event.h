#ifndef EVENT_H
#define EVENT_H

#include "TClonesArray.h"
/////////////////////////////////////////////////
//                                             //
//      Particle Class                         //
//                                             //
/////////////////////////////////////////////////

#include"particle.h"

class event
{
 public:
  //membri pubblici
 event() {fPart = new TClonesArray("particle");};

 virtual ~event() {};
 const particle *GetParticle(Int_t i) const {if(i >= 0 && i < GetEntries()) return (particle *) fPart->At(i); else return NULL;};
 Int_t GetEntries() const {return fPart->GetEntries();};
 void SetParticle(Int_t i,particle part) {TClonesArray &ar = *fPart;new(ar[i]) particle(part);};

 void Reset(){for(Int_t i=GetEntries()-1;i >= 0;i--) if(GetParticle(i)) fPart->RemoveAt(i);};

 private:

  TClonesArray *fPart;

  ClassDef(event,1)

};

#endif
