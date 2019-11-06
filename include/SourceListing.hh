/*
 * SourceListing.hh
 *
 *  Created on: 15 Mar 2016
 *      Author: sarwarr
 */

#ifndef INCLUDE_SOURCELISTING_HH_
#define INCLUDE_SOURCELISTING_HH_

#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include "globals.hh"
using namespace std;
class SourceListing {
 public:
  SourceListing();
  virtual ~SourceListing();
  static SourceListing *Instance();

  G4int GetTotalNSource() { return TotalNSource; };
  G4double GetTotalActivity() { return TotalActivity; };
  G4String *GetSourceNames() { return SourceNames; };
  G4double *GetSourceIntensities() { return SourceIntensities; };
  G4String *GetParticleNames() { return ParticleNames; };
  vector<G4double> *GetMultiplicityDistribution() {
    return MultiplicityDistribution;
  };
  vector<G4double> *GetWattSpectrum() { return WattSpectrum; };
  G4bool *IsWatt() { return Watt; };
  G4int *GetNEnergyBins() { return NEnergyBins; };
  vector<G4double> *GetEnergyBins() { return EnergyBins; };
  vector<G4double> *GetIntensityBins() { return IntensityBins; };
  G4int *GetIntervelType() { return IntervelType; };
  G4bool UseDefaultSource() { return !end; };

  bool parse_line(vector<string> *tok);

 private:
  static SourceListing *fgInstance;

  G4int TotalNSource;
  G4double TotalActivity;
  G4String SourceNames[10];
  G4double SourceIntensities[10];
  G4String ParticleNames[10];
  G4double MultiplicityDistSize[10];
  vector<G4double> MultiplicityDistribution[10];
  vector<G4double> WattSpectrum[10];
  G4bool Watt[10];
  G4int NEnergyBins[10];
  vector<G4double> EnergyBins[10];
  vector<G4double> IntensityBins[10];
  G4int IntervelType[10];
  G4bool end;
  ifstream fdSourceListing;
};

#endif /* INCLUDE_SOURCELISTING_HH_ */
