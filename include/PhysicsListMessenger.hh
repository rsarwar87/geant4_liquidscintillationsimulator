
#ifndef NeutronHPMessenger_h
#define NeutronHPMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class NeutronHPphysics;
class G4UIdirectory;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class NeutronHPMessenger : public G4UImessenger {
 public:
  NeutronHPMessenger(NeutronHPphysics*);
  ~NeutronHPMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

 private:
  NeutronHPphysics* fNeutronPhysics;

  G4UIdirectory* fPhysDir;
  G4UIcmdWithABool* fThermalCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
