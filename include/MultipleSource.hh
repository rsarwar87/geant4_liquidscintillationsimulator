/*
 * MultipleSource.hh
 *
 *  Created on: 30 May 2017
 *      Author: sarwarr
 */

#ifndef WORKSPACE_INCLUDE_MULTIPLESOURCE_HH_
#define WORKSPACE_INCLUDE_MULTIPLESOURCE_HH_

#include <vector>
#include "globals.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "SingleSource.hh"

class MultipleSource : public G4VPrimaryGenerator {
  //
 public:
  MultipleSource();
  MultipleSource(SingleSource* src, G4double);
  ~MultipleSource();

  void GeneratePrimaryVertex(G4Event*);

  G4int GetNumberofSource() { return G4int(sourceVector.size()); };
  void ListSource();
  void SetCurrentSourceto(G4int);
  void SetCurrentSourceIntensity(G4double);
  SingleSource* GetCurrentSource() { return currentSource; };
  G4int GetCurrentSourceIndex() { return currentSourceIdx; };
  G4double GetCurrentSourceIntensity() {
    return sourceIntensity[currentSourceIdx];
  };
  void ClearAll();
  void AddaSource(G4double);
  void AddaSource(SingleSource*, G4double);
  void DeleteaSource(G4int);

  // Set the verbosity level.
  void SetVerbosity(G4int i) { currentSource->SetVerbosity(i); };

  // Set if multiple vertex per event.
  void SetMultipleVertex(G4bool av) { multiple_vertex = av; };

  // set if flat_sampling is applied in multiple source case

  void SetFlatSampling(G4bool av) {
    flat_sampling = av;
    normalised = false;
  };

  // Set the particle species
  void SetParticleDefinition(G4ParticleDefinition* aParticleDefinition) {
    currentSource->SetParticleDefinition(aParticleDefinition);
  };

  G4ParticleDefinition* GetParticleDefinition() {
    return currentSource->GetParticleDefinition();
  };

  void SetParticleCharge(G4double aCharge) {
    currentSource->SetParticleCharge(aCharge);
  };

  // Set polarization
  void SetParticlePolarization(G4ThreeVector aVal) {
    currentSource->SetParticlePolarization(aVal);
  };
  G4ThreeVector GetParticlePolarization() {
    return currentSource->GetParticlePolarization();
  };

  // Set Time.
  void SetParticleTime(G4double aTime) {
    currentSource->SetParticleTime(aTime);
  };
  G4double GetParticleTime() { return currentSource->GetParticleTime(); };

  void SetNumberOfParticles(G4int i) {
    currentSource->SetNumberOfParticles(i);
  };
  //
  G4int GetNumberOfParticles() {
    return currentSource->GetNumberOfParticles();
  };
  G4ThreeVector GetParticlePosition() {
    return currentSource->GetParticlePosition();
  };
  G4ThreeVector GetParticleMomentumDirection() {
    return currentSource->GetParticleMomentumDirection();
  };
  G4double GetParticleEnergy() { return currentSource->GetParticleEnergy(); };

 private:
  void IntensityNormalization();

 private:
  G4bool multiple_vertex;
  G4bool flat_sampling;
  G4bool normalised;
  G4int currentSourceIdx;
  SingleSource* currentSource;
  std::vector<SingleSource*> sourceVector;
  std::vector<G4double> sourceIntensity;
  std::vector<G4double> sourceProbability;
};

#endif /* WORKSPACE_INCLUDE_MULTIPLESOURCE_HH_ */
