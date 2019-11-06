/*
 * EMPhysics.hh
 *
 *  Created on: 25 May 2017
 *      Author: sarwarr
 */

#ifndef WORKSPACE_INCLUDE_EMPHYSICS_HH_
#define WORKSPACE_INCLUDE_EMPHYSICS_HH_

#include "G4ios.hh"
#include "globals.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eplusAnnihilation.hh"

class EMPhysics : public G4VPhysicsConstructor {
 public:
  EMPhysics(const G4String& name = "EM");
  virtual ~EMPhysics();

 public:
  // This method will be invoked in the Construct() method.
  // each particle type will be instantiated
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();
};

#endif /* WORKSPACE_INCLUDE_EMPHYSICS_HH_ */
