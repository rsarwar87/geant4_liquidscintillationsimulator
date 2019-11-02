/*
 * GeneralPhysics.hh
 *
 *  Created on: 25 May 2017
 *      Author: sarwarr
 */

#ifndef WORKSPACE_INCLUDE_GENERALPHYSICS_HH_
#define WORKSPACE_INCLUDE_GENERALPHYSICS_HH_



#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

class GeneralPhysics: public G4VPhysicsConstructor
{
public:

	GeneralPhysics(const G4String& name = "general");
	virtual ~GeneralPhysics();

	// This method will be invoked in the Construct() method.
	// each particle type will be instantiated
	virtual void ConstructParticle();

	// This method will be invoked in the Construct() method.
	// each physics process will be instantiated and
	// registered to the process manager of each particle type
	virtual void ConstructProcess();

};



#endif /* WORKSPACE_INCLUDE_GENERALPHYSICS_HH_ */
