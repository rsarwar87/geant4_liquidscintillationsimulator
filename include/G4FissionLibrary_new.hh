/*
 * G4FissionLibrary_new.hh
 *
 *  Created on: 30 May 2017
 *      Author: sarwarr
 */

#ifndef WORKSPACE_INCLUDE_G4FISSIONLIBRARY_NEW_HH_
#define WORKSPACE_INCLUDE_G4FISSIONLIBRARY_NEW_HH_


#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4DynamicParticleVector.hh"
#include "G4HadFinalState.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4ParticleHPNames.hh"
#include "G4ParticleHPParticleYield.hh"
#include "G4ParticleHPFissionERelease.hh"
#include "G4ParticleHPEnergyDistribution.hh"
#include "G4ParticleHPPhotonDist.hh"
#include "G4ParticleHPAngular.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4LLNLFission.hh"
#include "fissionEvent.h"

class G4FissionLibrary_new: public G4ParticleHPFinalState
{
public:

	G4FissionLibrary_new();
	~G4FissionLibrary_new();

	//void Init (G4double A, G4double Z, G4String & dirName, G4String &);
	void Init(G4double A, G4double Z, G4int M, G4String & dirName, G4String &,
			G4ParticleDefinition*);
	G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack);
	G4ParticleHPFinalState * New();

private:
	fissionEvent* fe;
	G4int theIsotope; // used to call fissionEvent constructor
	G4int theZ;
	G4int theA;
	G4IonTable* theTableOfIons;
	G4double targetMass;
	void SampleMult(const G4HadProjectile & theTrack, G4int* nPrompt,
			G4int* gPrompt, G4double eKinetic);
	inline G4ParticleHPFissionERelease * GetEnergyRelease()
	{
		return &theEnergyRelease;
	}
	G4ParticleHPParticleYield theFinalStateNeutrons;
	G4ParticleHPEnergyDistribution thePromptNeutronEnDis;
	G4ParticleHPEnergyDistribution theDelayedNeutronEnDis;
	G4ParticleHPAngular theNeutronAngularDis;

	G4ParticleHPPhotonDist theFinalStatePhotons;
	G4ParticleHPFissionERelease theEnergyRelease;
};

#endif /* WORKSPACE_INCLUDE_G4FISSIONLIBRARY_NEW_HH_ */
