/*
 * SingleSource.hh
 *
 *  Created on: 30 May 2017
 *      Author: sarwarr
 */

#ifndef WORKSPACE_INCLUDE_SINGLESOURCE_HH_
#define WORKSPACE_INCLUDE_SINGLESOURCE_HH_



#include "G4VPrimaryGenerator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
//
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSRandomGenerator.hh"

class SingleSource: public G4VPrimaryGenerator
{
public:
	SingleSource();
	~SingleSource();
	virtual void GeneratePrimaryVertex(G4Event *evt);
	//
	G4SPSPosDistribution* GetPosDist()
	{
		return posGenerator;
	}
	;
	G4SPSAngDistribution* GetAngDist()
	{
		return angGenerator;
	}
	;
	G4SPSEneDistribution* GetEneDist()
	{
		return eneGenerator;
	}
	;
	G4SPSRandomGenerator* GetBiasRndm()
	{
		return biasRndm;
	}
	;

	// Set the verbosity level.
	void SetVerbosity(G4int);

	// Set the particle species
	void SetParticleDefinition(G4ParticleDefinition * aParticleDefinition);
	inline G4ParticleDefinition * GetParticleDefinition()
	{
		return particle_definition;
	}
	;

	inline void SetParticleCharge(G4double aCharge)
	{
		particle_charge = aCharge;
	}
	;

	// Set polarization
	inline void SetParticlePolarization(G4ThreeVector aVal)
	{
		particle_polarization = aVal;
	}
	;
	inline G4ThreeVector GetParticlePolarization()
	{
		return particle_polarization;
	}
	;

	// Set Time.
	inline void SetParticleTime(G4double aTime)
	{
		particle_time = aTime;
	}
	;
	inline G4double GetParticleTime()
	{
		return particle_time;
	}
	;

	inline void SetNumberOfParticles(G4int i)
	{
		NumberOfParticlesToBeGenerated = i;
	}
	;
	//
	inline G4int GetNumberOfParticles()
	{
		return NumberOfParticlesToBeGenerated;
	}
	;
	inline G4ThreeVector GetParticlePosition()
	{
		return particle_position;
	}
	;
	inline G4ThreeVector GetParticleMomentumDirection()
	{
		return particle_momentum_direction;
	}
	;
	inline G4double GetParticleEnergy()
	{
		return particle_energy;
	}
	;

private:

	G4SPSPosDistribution* posGenerator;
	G4SPSAngDistribution* angGenerator;
	G4SPSEneDistribution* eneGenerator;
	G4SPSRandomGenerator* biasRndm;
	//
	// Other particle properties
	G4int NumberOfParticlesToBeGenerated;
	G4ParticleDefinition * particle_definition;
	G4ParticleMomentum particle_momentum_direction;
	G4double particle_energy;
	G4double particle_charge;
	G4ThreeVector particle_position;
	G4double particle_time;
	G4ThreeVector particle_polarization;
	G4double particle_weight;

	// Verbosity
	G4int verbosityLevel;

};


#endif /* WORKSPACE_INCLUDE_SINGLESOURCE_HH_ */
