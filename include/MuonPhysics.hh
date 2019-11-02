/*
 * MuonPhysics.hh
 *
 *  Created on: 25 May 2017
 *      Author: sarwarr
 */

#ifndef WORKSPACE_INCLUDE_MUONPHYSICS_HH_
#define WORKSPACE_INCLUDE_MUONPHYSICS_HH_


#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuIonisation.hh"
#include "G4hIonisation.hh"

#include "G4MuonMinusCapture.hh"

class MuonPhysics: public G4VPhysicsConstructor
{
public:

	MuonPhysics(const G4String& name = "muon");
	virtual ~MuonPhysics();

	// This method will be invoked in the Construct() method.
	// each particle type will be instantiated
	virtual void ConstructParticle();

	// This method will be invoked in the Construct() method.
	// each physics process will be instantiated and
	// registered to the process manager of each particle type
	virtual void ConstructProcess();

};

#endif
