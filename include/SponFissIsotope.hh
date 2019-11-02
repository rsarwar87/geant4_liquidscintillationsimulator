/*
 * SponFissIsotope.hh
 *
 *  Created on: 30 May 2017
 *      Author: sarwarr
 */

#ifndef WORKSPACE_INCLUDE_SPONFISSISOTOPE_HH_
#define WORKSPACE_INCLUDE_SPONFISSISOTOPE_HH_

#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4LLNLFission.hh"
#include "G4Neutron.hh"
#include "G4PrimaryVertex.hh"
#include "Randomize.hh"
#include "SingleSource.hh"
#ifdef FISSION_NEW
#include "fissionEvent.h"
#else
#include "G4fissionEvent.hh"
#endif

class G4Event;

class SponFissIsotope : public SingleSource {
 public:
  SponFissIsotope();
  SponFissIsotope(G4int iso);
  ~SponFissIsotope();

 public:
  void GeneratePrimaryVertex(G4Event* anEvent);
  // Set the verbosity level.
  void SetVerbosity(G4int verb) { verbosityLevel = verb; };

  static bool angular_correlation; 
 private:
  G4int isotope;
  G4ThreeVector particle_polarization;
  G4ParticleDefinition* neutron_definition;
  G4ParticleDefinition* photon_definition;
  G4SPSPosDistribution* posDist;

  // Verbosity
  G4int verbosityLevel;
};

#endif /* WORKSPACE_INCLUDE_SPONFISSISOTOPE_HH_ */
