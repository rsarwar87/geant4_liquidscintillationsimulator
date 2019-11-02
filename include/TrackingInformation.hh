#ifndef TrackInformation_h
#define TrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class TrackInformation : public G4VUserTrackInformation 
{
public:
  TrackInformation();
  TrackInformation(const G4Track* aTrack);
  TrackInformation(const TrackInformation* aTrackInfo);
  virtual ~TrackInformation();
   
  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);

  TrackInformation& operator =(const TrackInformation& right);
  
  void SetSourceTrackInformation(const G4Track* aTrack);
  virtual void Print() const;

public:
  inline G4int GetTrackingStatus() const {return fTrackingStatus;}
  inline void  SetTrackingStatus(G4int i) {fTrackingStatus = i;}
  inline G4int GetSourceTrackID() const {return fSourceTrackID;}

private:
  // Information of the primary track at the primary vertex
  G4int                 fOriginalTrackID;  // Track ID of primary particle
  G4ParticleDefinition* fParticleDefinition;
  G4ThreeVector         fOriginalPosition;
  G4ThreeVector         fOriginalMomentum;
  G4double              fOriginalEnergy;
  G4double              fOriginalTime;

  G4int                 fTrackingStatus;
  // trackingStatus = 1 : primary or secondary track which has not yet reached 
  //                      to calorimeter
  //                = 0 : track which or ancester of which has reached to 
  //                      calorimeter
  //                = 2 : track or its ancester had once reached to calorimeter 
  //                      and then escaped from it
  // Information of the track which reached to the calorimeter boundary at the 
  // boundary surface. 
  // This information is valid only for trackingStatus = 0 or 2
  G4int                 fSourceTrackID;
  G4ParticleDefinition* fSourceDefinition;
  G4ThreeVector         fSourcePosition;
  G4ThreeVector         fSourceMomentum;
  G4double              fSourceEnergy;
  G4double              fSourceTime;
public:
  G4String                 fParentType;
  G4int                 fParentID;
  G4int                 fID;
};

extern G4ThreadLocal G4Allocator<TrackInformation> * aTrackInformationAllocator;
//extern  G4Allocator<TrackInformation> aTrackInformationAllocator;

inline void* TrackInformation::operator new(size_t)
{
  if(!aTrackInformationAllocator)
    aTrackInformationAllocator = new G4Allocator<TrackInformation>;
  return (void*)aTrackInformationAllocator->MallocSingle();
}

inline void TrackInformation::operator delete(void *aTrackInfo)
{ aTrackInformationAllocator->FreeSingle((TrackInformation*)aTrackInfo);}

#endif