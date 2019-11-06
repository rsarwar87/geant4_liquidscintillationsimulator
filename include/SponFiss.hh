/*
 * SponFiss.hh
 *
 *  Created on: 31 May 2017
 *      Author: sarwarr
 */

#ifndef INCLUDE_SPONFISS_HH_
#define INCLUDE_SPONFISS_HH_
#define FISSION_NEW
#define USEFREYA

#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4LLNLFission.hh"
#include "G4Neutron.hh"
#include "G4PrimaryVertex.hh"
#include "G4SPSPosDistribution.hh"
#include "Randomize.hh"
#ifdef FISSION_NEW
#include "fissionEvent.h"
#else
#include "G4fissionEvent.hh"
#endif

#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4LLNLFission.hh"
#include "G4Neutron.hh"
#include "G4PrimaryVertex.hh"
#include "MyRun.hh"
#include "Randomize.hh"
#ifdef FISSION_NEW
#include "fissionEvent.h"
#else
#include "G4fissionEvent.hh"
#endif
class G4Event;

class SponFiss {
 public:
  SponFiss();
  SponFiss(G4int iso,
           G4SPSPosDistribution* pos /*, PrimaryGeneratorAction *run*/);
  ~SponFiss();

 public:
  void GeneratePrimaryVertex(G4Event* anEvent, G4double time = 0, int mode = 3);
  // Set the verbosity level.
  void SetVerbosity(G4int verb) { verbosityLevel = verb; };

 private:
  G4int isotope;
  G4ThreeVector particle_polarization;
  G4ParticleDefinition* neutron_definition;
  G4ParticleDefinition* photon_definition;
  G4SPSPosDistribution* posDist;
  // Verbosity
  G4int verbosityLevel;
  // PrimaryGeneratorAction fRun;
};

#endif /* INCLUDE_SPONFISS_HH_ */
