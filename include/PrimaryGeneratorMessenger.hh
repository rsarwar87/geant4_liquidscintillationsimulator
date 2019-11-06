/*
 * PrimaryGeneratorMessenger.hh
 *
 *  Created on: 15 Mar 2016
 *      Author: sarwarr
 */

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorMessenger : public G4UImessenger {
 public:
  PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
  virtual ~PrimaryGeneratorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

 private:
  PrimaryGeneratorAction* fOpNoviceAction;
  G4UIdirectory* fGunDir;
  G4UIcmdWithADoubleAndUnit* fPolarCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
