

#ifndef Run_h
#define Run_h 1

#include <map>
#include "G4Run.hh"
#include "G4VProcess.hh"
#include "SteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run {
 public:
  Run(DetectorConstruction *);
  ~Run();

 public:
  void CountProcesses(const G4VProcess *process);
  void ParticleCount(G4String, G4double);
  void SumTrackLength(G4int, G4int, G4double, G4double, G4double, G4double);

  void SetPrimary(G4ParticleDefinition *particle, G4double energy);

  /*void UpdateSource(int n_nu, int g_nu);
  void UpdateSource(int id, int nu, G4double energy);
   void UpdateSource(int id, int nu, int ang, int energy);*/

  void EndOfRun();

  virtual void Merge(const G4Run *);
  virtual void RecordEvent(const G4Event *evt);

  void Reset();

  // get methods

  const G4int *GetParticleIncidentEnergy(int pidx = 0) const {
    return fIncidentEnergy[pidx];
  }
  const G4int *GetParticleEntryEnergy(int pidx = 0) const {
    return fEntryEnergy[pidx];
  }
  const G4int *GetParticleDeposit(int pidx = 0) const {
    return fParticleDeposit[pidx];
  }
  const G4int *GetElectronDeposited(int pidx = 0) const {
    return fElectronDeposited[pidx];
  }
  const G4int *GetProtonDeposited(int pidx = 0) const {
    return fProtonDeposited[pidx];
  }
  const G4int *GetElectronProduced(int pidx = 0) const {
    return fElectronProduced[pidx];
  }
  const G4int *GetProtonProduced(int pidx = 0) const {
    return fProtonProduced[pidx];
  }
  const G4int *GetOphotonProduced(int pidx = 0) const {
    return fOphotonProduced[pidx];
  }
  const G4int *GetPMTDistribution(int pidx = 0) const {
    return fOphotonDeposited[pidx];
  }
  const G4int *GetCerenkovDistribution(int pidx = 0) const {
    return fCerenkovProduced[pidx];
  }

  const G4int *GetLightResponse(int pidx = 0) const {
    return fEnergyTime[pidx][1];
  }
  const G4int *GetLightHistogram(int pidx = 0) const {
    return fCountTime[pidx][1];
  }
  const G4int *GetPMTResponse(int pidx = 0) const { return fEnergyTime[pidx][0]; }
  const G4int *GetPMTHistogram(int pidx = 0) const {
    return fCountTime[pidx][0];
  }

  const G4int GetTotalRegisterPartile(int idx) const {
    return fEventRegistered[idx];
  }

  const G4int *GetGammaSourceSpec() const { return spec_theroy[0]; }
  const G4int *GetNeutronSourceSpec() const { return spec_theroy[1]; }

  const G4int *GetGammaSourceMultiplicity() const { return multi_theroy[0]; }
  const G4int *GetNeutronSourceMultiplicity() const { return multi_theroy[1]; }
  const G4int *GetJointSourceMultiplicity() const { return multi_theroy[2]; }

  const G4int *GetGammaDetectedMultiplicity() const {
    return multi_detected[0];
  }
  const G4int *GetNeutronDetectedMultiplicity() const {
    return multi_detected[1];
  }
  const G4int *GetJointDetectedMultiplicity() const {
    return multi_detected[2];
  }
  const G4int *GetGammaDetectedCXMultiplicity() const {
    return multi_detectedcx[0];
  }
  const G4int *GetNeutronDetectedCXMultiplicity() const {
    return multi_detectedcx[1];
  }
  const G4int *GetJointDetectedCXMultiplicity() const {
    return multi_detectedcx[2];
  }

  G4int *GetRossi(int id) const { return rossi_[id]; }
  G4int **GetTime() const { return _time; }
  G4int *GetRossiCX(int id) const { return rossi_cx[id]; }
  G4int *DetectedSpectrum() const { return n_spec; }

  G4int **GetAngularContour() const { return angular_contour; }
  G4int **GetAngularPlot() const { return angular_plot; }

  G4int **GetAngularCCXontour() const { return angular_contourcx; }
  G4int **GetAngularCXPlot() const { return angular_plotcx; }
  static G4Mutex aMutex;

  G4int gEventNumber = 0;

 private:
  struct ParticleData {
    ParticleData() : fCount(0), fEmean(0.), fEmin(0.), fEmax(0.) {}
    ParticleData(G4int count, G4double ekin, G4double emin, G4double emax)
        : fCount(count), fEmean(ekin), fEmin(emin), fEmax(emax) {}
    G4int fCount;
    G4double fEmean;
    G4double fEmin;
    G4double fEmax;
  };

 private:
  struct entry {
    int didx = -1;
    int pid = -1;
  };
  G4int **fIncidentEnergy;
  G4int **fEntryEnergy;
  G4int **fParticleDeposit;
  G4int **fElectronProduced;
  G4int **fElectronDeposited;
  G4int **fProtonProduced;
  G4int **fProtonDeposited;
  G4int **fOphotonProduced;
  G4int **fCerenkovProduced;
  G4int **fOphotonDeposited;

  G4int *spec_theroy[2];
  G4int *n_spec;

  G4int *multi_theroy[3];
  G4int *multi_detected[3];
  G4int *multi_detectedcx[3];

  G4int *rossi_[3];
  G4int **_time;
  G4int *rossi_cx[3];

  G4int **angular_plot;
  G4int **angular_contour;
  G4int **angular_plotcx;
  G4int **angular_contourcx;

  G4int fEventNumber;
  G4int *fEventRegistered;

  G4int ***fEnergyTime;
  G4int ***fCountTime;

  std::ofstream fout;
  G4String filename;

 private:
  DetectorConstruction *fDetector;
  G4ParticleDefinition *fParticle;
  G4double fEkin;

  std::map<G4String, G4int> fProcCounter;
  std::map<G4String, ParticleData> fParticleDataMap;

  G4int fNbStep1, fNbStep2;
  G4double fTrackLen1, fTrackLen2;
  G4double fTime1, fTime2;

  void ProcessCoincidence(std::map<double, entry> storage,
                          std::vector<int> &angular,
                          std::vector<int> &angularcx, int &multi,
                          int &multi_cx, int *rossi, int *rossicxs);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
