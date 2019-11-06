/*
 * G4FissLib_new.hh
 *
 *  Created on: 30 May 2017
 *      Author: sarwarr
 */

#ifndef WORKSPACE_INCLUDE_G4FISSLIB_NEW_HH_
#define WORKSPACE_INCLUDE_G4FISSLIB_NEW_HH_

#include "G4HadronicInteraction.hh"
#include "G4NeutronHPChannel.hh"
#include "G4NeutronHPThermalBoost.hh"
#include "globals.hh"

class G4FissLib_new : public G4HadronicInteraction {
 public:
  G4FissLib_new();
  ~G4FissLib_new();

  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
                                 G4Nucleus& aTargetNucleus);
  const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const;

 private:
  G4double* xSec;
  G4NeutronHPChannel* theFission;
  G4String dirName;
  G4int numEle;
};

#endif /* WORKSPACE_INCLUDE_G4FISSLIB_NEW_HH_ */
