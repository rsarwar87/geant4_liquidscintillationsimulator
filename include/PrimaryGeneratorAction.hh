//******************************************************************************
// PrimaryGeneratorAction.hh
//
// This class is a class derived from G4VUserPrimaryGeneratorAction for
// constructing the process used to generate incident particles.
//
// 1.00 JMV, LLNL, JAN-2007:  First version.
//******************************************************************************ss
//
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <iomanip>
#include <random>
#include "G4DataVector.hh"
#include "G4ParticleGun.hh"
#include "G4SPSPosDistribution.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "MultipleSource.hh"
#include "Randomize.hh"
#include "SponFIss_FIF.hh"
#include "SponFiss.hh"
#include "SponFissIsotope.hh"
#include "vector"

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction();

  G4ParticleGun* GetParticleGun() { return fParticleGun; };

  static G4String name;
  static int energy;
  static int mode;
  static bool gamma;
  static bool neutron;
  static bool mono;
  static bool AmLi;
  static bool beam;
  static bool Co;
  static bool sfif;

 public:
  void GeneratePrimaries(G4Event* anEvent);

 private:
  G4ParticleGun* fParticleGun;
  static G4Mutex aMutex;
  SponFiss* iso;
  // SponFiss_FF* fif;
  G4SPSPosDistribution* posDist;

  MultipleSource* fissionSource;
  G4double time;
  G4double totalIntensity;
  G4int nisotopes;

  std::random_device rd;
  static std::discrete_distribution<> intensity_dis;
  static std::discrete_distribution<> leak_dis;
  static double decay_time;
};

#endif
