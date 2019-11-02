/*
 * SponFiss.hh
 *
 *  Created on: 27 Feb 2018
 *      Author: sarwarr
 */

#ifndef INCLUDE_SPONFISS_FF_HH_
#define INCLUDE_SPONFISS_FF_HH_

#include "G4DynamicParticle.hh"
#include "globals.hh"
#include "G4PrimaryVertex.hh"
#include "G4Event.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4SPSPosDistribution.hh"
#include <fstream>
#include <iostream>
#include <string>


class SponFiss_FF
{
public:
	SponFiss_FF(){};
	SponFiss_FF(G4SPSPosDistribution *pos);
	~SponFiss_FF();

public:
	void GeneratePrimaryVertex(G4Event* anEvent, G4double time = 0, int mode = 3);
	// Set the verbosity level.
	void SetVerbosity(G4int verb)
	{
		verbosityLevel = verb;
	}
	;
	static std::string filename;
private:
	G4int isotope;
	static std::ifstream *bfile;
	bool fopen = false;
	void update_angular(std::vector<G4ThreeVector> &fe);
	G4ThreeVector particle_polarization;
	G4ParticleDefinition* neutron_definition;
	G4ParticleDefinition* photon_definition;
	G4SPSPosDistribution* posDist;
	G4int verbosityLevel;
};



#endif /* INCLUDE_SPONFISS_HH_ */
