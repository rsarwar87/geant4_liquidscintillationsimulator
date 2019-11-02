/*
 * SingleSource.cc
 *
 *  Created on: 30 May 2017
 *      Author: sarwarr
 */

#include "SingleSource.hh"
#include <cmath>
#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "Randomize.hh"

SingleSource::SingleSource() {
  // Initialise all variables
  // Position distribution Variables

  NumberOfParticlesToBeGenerated = 1;
  particle_definition = G4Geantino::GeantinoDefinition();
  G4ThreeVector zero;
  particle_momentum_direction = G4ParticleMomentum(1, 0, 0);
  particle_energy = 1.0 * MeV;
  particle_position = zero;
  particle_time = 0.0;
  particle_polarization = zero;
  particle_charge = 0.0;
  particle_weight = 1.0;

  biasRndm = new G4SPSRandomGenerator();
  posGenerator = new G4SPSPosDistribution();
  posGenerator->SetBiasRndm(biasRndm);
  angGenerator = new G4SPSAngDistribution();
  angGenerator->SetPosDistribution(posGenerator);
  angGenerator->SetBiasRndm(biasRndm);
  eneGenerator = new G4SPSEneDistribution();
  eneGenerator->SetBiasRndm(biasRndm);

  // verbosity
  verbosityLevel = 0;
}

SingleSource::~SingleSource() {}

void SingleSource::SetVerbosity(int vL) {
  verbosityLevel = 0;
  posGenerator->SetVerbosity(verbosityLevel);
  angGenerator->SetVerbosity(verbosityLevel);
  eneGenerator->SetVerbosity(verbosityLevel);
  G4cout << "Verbosity Set to: " << verbosityLevel << G4endl;
}

void SingleSource::SetParticleDefinition(
    G4ParticleDefinition* aParticleDefinition) {
  particle_definition = aParticleDefinition;
  particle_charge = particle_definition->GetPDGCharge();
}

void SingleSource::GeneratePrimaryVertex(G4Event* evt) {
  if (particle_definition == NULL) return;
  if (verbosityLevel > 1)
    G4cout << " NumberOfParticlesToBeGenerated: "
           << NumberOfParticlesToBeGenerated << G4endl;

  // Position stuff
  particle_position = posGenerator->GenerateOne();

  // create a new vertex
  G4PrimaryVertex* vertex =
      new G4PrimaryVertex(particle_position, particle_time);

  for (G4int i = 0; i < NumberOfParticlesToBeGenerated; i++) {
    // Angular stuff
    particle_momentum_direction = angGenerator->GenerateOne();
    // Energy stuff
    particle_energy = eneGenerator->GenerateOne(particle_definition);

    if (verbosityLevel >= 2)
      G4cout << "Creating primaries and assigning to vertex" << G4endl;
    // create new primaries and set them to the vertex
    G4double mass = particle_definition->GetPDGMass();
    G4double energy = particle_energy + mass;
    G4double pmom = std::sqrt(energy * energy - mass * mass);
    G4double px = pmom * particle_momentum_direction.x();
    G4double py = pmom * particle_momentum_direction.y();
    G4double pz = pmom * particle_momentum_direction.z();

    if (verbosityLevel > 1) {
      G4cout << "Particle name: " << particle_definition->GetParticleName()
             << G4endl;
      G4cout << "       Energy: " << particle_energy << G4endl;
      G4cout << "     Position: " << particle_position << G4endl;
      G4cout << "    Direction: " << particle_momentum_direction << G4endl;
    }
    G4PrimaryParticle* particle =
        new G4PrimaryParticle(particle_definition, px, py, pz);
    particle->SetMass(mass);
    particle->SetCharge(particle_charge);
    particle->SetPolarization(particle_polarization.x(),
                              particle_polarization.y(),
                              particle_polarization.z());
    // Set bweight equal to the multiple of all non-zero weights
    particle_weight = biasRndm->GetBiasWeight();
    // pass it to primary particle
    particle->SetWeight(particle_weight);

    vertex->SetPrimary(particle);
  }
  // now pass the weight to the primary vertex. CANNOT be used here!
  //  vertex->SetWeight(particle_weight);
  evt->AddPrimaryVertex(vertex);
  if (verbosityLevel > 1) G4cout << " Primary Vetex generated !" << G4endl;
}

