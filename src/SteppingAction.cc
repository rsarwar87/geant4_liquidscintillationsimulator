

#include "SteppingAction.hh"
#include "G4EventManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4SDManager.hh"
#include "G4StepPoint.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackStatus.hh"
#include "G4UnitsTable.hh"
#include "G4VPhysicalVolume.hh"
#include "HistoManager.hh"
#include "MyRun.hh"
#include "TrackingAction.hh"
#include "TrackingInformation.hh"

#include "G4RunManager.hh"
#include "RunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal SteppingAction* SteppingAction::fgInstance = 0;

__thread RecodedEvent* SteppingAction::crEvent = NULL;
__thread RecodedParticle* SteppingAction::crParticle = NULL;

__thread G4int SteppingAction::_cnnt = 0;
__thread G4int SteppingAction::_cnnt2 = -1;
__thread G4double SteppingAction::_eng_l = -1;
__thread G4double SteppingAction::_tm_l = 1000000;

int SteppingAction::gwidth = 20;
bool SteppingAction::spec = false;
SteppingAction::SteppingAction(G4String fn, TrackingAction* TrAct)
    : G4UserSteppingAction(), fTrackingAction(TrAct) {
  // Reset();
  filename = fn;
  fgInstance = this;
  fout.open("SteppingAction", std::ios::out | std::ios::trunc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() {
  fout.close();
  fgInstance = 0;
}

void SteppingAction::Reset() {
  // release particle buffers
  for (int i = 0; i < 2; i++)
    if (rParticle[i].size() > 0) {
      for (int j = 0; j < rParticle[i].size(); j++) {
        if (rParticle[i].at(j) != NULL) delete rParticle[i].at(j);
        rParticle[i].at(j) = NULL;
      }
      rParticle[i].clear();
    }
  // release event buffers
  if (rEvent.size() > 0) {
    for (int j = 0; j < rEvent.size(); j++) {
      rEvent.at(j)->delete_class();
      if (rEvent.at(j) != NULL) delete rEvent.at(j);
      rEvent.at(j) = NULL;
    }
    rEvent.clear();
  }
  crEvent = NULL;
  crParticle = NULL;
  _cnnt = 0;
  _cnnt2 = -1;
  _eng_l = 1;
  _tm_l = 1000000;
  return;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction* SteppingAction::Instance() { return fgInstance; }

void SteppingAction::StackParticle(const G4Step* step,
                                   const G4StepPoint* point) {
  G4Track* track = step->GetTrack();
  if ((_eng_l != step->GetPreStepPoint()->GetKineticEnergy()) &&
      (_cnnt2 != track->GetTrackID())) {
    TrackInformation* trackInfo =
        (TrackInformation*)(track->GetUserInformation());
    _cnnt++;

    _cnnt2 = track->GetTrackID();
    _tm_l = step->GetTrack()->GetGlobalTime();
    ParticleName =
        track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
    _eng_l = point->GetKineticEnergy();
    crParticle =
        new RecodedParticle(_cnnt, ParticleName, _eng_l, track->GetGlobalTime(),
                            point->GetMomentumDirection());
    // if (ParticleName != "neutron")
    // crParticle->Pout(_cnnt*1000+track->GetParentID() );
    trackInfo->fID = _cnnt;
    if (track->GetParentID() == 0) trackInfo->fParentType = ParticleName;
    track->SetUserInformation(trackInfo);
    // if (ParticleName != "neutron")
    {
      /*G4cout << ParticleName << " " << trackInfo->fParentID  << " " <<
                      trackInfo->fParentType  << " " << trackInfo->fID << " "
                      << track->GetKineticEnergy()/keV<< G4endl;*/
      // track->GetDynamicParticle()->GetPrimaryParticle()->Print();
    }

    rParticle[ParticleName == "neutron" ? 1 : 0].push_back(crParticle);
  }
}
void SteppingAction::UserSteppingAction(const G4Step* step) {
  // count processes
  //

  const G4StepPoint* endPoint = step->GetPostStepPoint();
  const G4VProcess* process = endPoint->GetProcessDefinedStep();
  Run* run =
      static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->CountProcesses(process);

  //
  // collect information on the first particle
  if (step->GetTrack()->GetTrackID() == 1) {
    G4double ekin = endPoint->GetKineticEnergy();
    G4double trackl = step->GetTrack()->GetTrackLength();
    G4double time = step->GetTrack()->GetLocalTime();
    fTrackingAction->UpdateTrackInfo(ekin, trackl, time);
    G4AnalysisManager::Instance()->FillH1(7, ekin);
  }
  G4Track* track = step->GetTrack();

  // collect information on the particles generated and populate the buffers
  TrackInformation* trackInfo =
      (TrackInformation*)(track->GetUserInformation());
  if ((track->GetParentID() == 0)) {
    // if (step->GetPreStepPoint()->GetKineticEnergy()/keV < 300)
    // track->SetTrackStatus(fStopAndKill);
    StackParticle(step, step->GetPreStepPoint());

  } else if (trackInfo->GetTrackingStatus() > 0) {
    if ((track->GetDynamicParticle()
             ->GetParticleDefinition()
             ->GetParticleName() == "gamma") &&
        _eng_l != track->GetParentID()) {
      StackParticle(step, step->GetPreStepPoint());
      _eng_l = track->GetParentID();
    }
  }
  ParticleName =
      track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  if (step->GetTrack()->GetVolume()->GetName() != "Scintillator") return;

  //	if (trackInfo->GetTrackingStatus() > 0)
  {
    if (trackInfo->fID > _cnnt + 1) {
      // G4cout << " Particle not found " << track->GetKineticEnergy()/keV << "
      // " << trackInfo->fParentID << " " << _cnnt << G4endl;
      return;
    } else if (trackInfo->fID != crParticle->particleid) {
      RecodedParticle* tmp = NULL;

      for (int j = 0; j < 2; j++) {
        int sz = rParticle[j].size();
        int found = -1;
        for (int i = 0; i < sz; i++) {
          if (rParticle[j].at(i)->particleid == trackInfo->fID) {
            found = i;
            break;
          }
        }
        if (found != -1) {
          tmp = rParticle[j].at(found);
          break;
        }
      }
      if (tmp != NULL) {
        crParticle = tmp;
      } else {
        // G4cout << " Particle not found " << track->GetKineticEnergy()/keV <<
        // " " << trackInfo->fParentID  << G4endl;
        return;
      }
    }
  }

  // initialize the boarder process function.
  static G4ThreadLocal G4OpBoundaryProcess* boundary = NULL;
  if (!boundary) {
    G4ProcessManager* pm =
        step->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for (i = 0; i < nprocesses; i++) {
      if ((*pv)[i]->GetProcessName() == "OpBoundary") {
        boundary = (G4OpBoundaryProcess*)(*pv)[i];
        break;
      }
    }
  }

  // collect information on the detectors in which the simulation is taking
  // place
  int idx_det = step->GetTrack()->GetVolume()->GetCopyNo();
  if ((crEvent == NULL) || (crEvent->detectorID != idx_det)) {
    int sz = rEvent.size();
    int found = -1;
    for (int i = 0; i < sz; i++) {
      if (rEvent.at(i)->detectorID == idx_det) {
        found = i;
        break;
      }
    }
    if (found != -1)
      crEvent = rEvent.at(found);
    else {
      crEvent = new RecodedEvent(idx_det, crParticle,
                                 step->GetTrack()->GetGlobalTime());
      rEvent.push_back(crEvent);
    }
  }

  // check if particle is new in the detector
  crEvent->CheckPart(crParticle);

  // G4cout << idx_det<<G4endl;
  // for optical photons, detect the boundary absorption os the particles and
  // count them
  if (track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() ==
      "opticalphoton") {
    TrackInformation* trackInfo2 =
        (TrackInformation*)(track->GetUserInformation());

    if (trackInfo2->fParentID > _cnnt + 1) {
      return;
    } else if (trackInfo2->fParentID != crParticle->particleid) {
      RecodedParticle* tmp = NULL;

      for (int j = 0; j < 2; j++) {
        int sz = rParticle[j].size();
        int found = -1;
        for (int i = 0; i < sz; i++) {
          if (rParticle[j].at(i)->particleid == trackInfo2->fParentID) {
            found = i;
            break;
          }
        }
        if (found != -1) {
          tmp = rParticle[j].at(found);
          break;
        }
      }
      if (tmp != NULL) {
        crParticle = tmp;
        crEvent->CheckPart(crParticle);
      } else
        return;
    }

    G4OpBoundaryProcessStatus boundaryStatus = boundary->GetStatus();
    if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
      double idx = 0;
      if (boundaryStatus == Detection) {
        crEvent->recordProduction(
            0, step->GetPreStepPoint()->GetTotalEnergy() * 425.4885,
            step->GetPostStepPoint()->GetGlobalTime());
      }
    }
    return;
  }

  // identify compton and elastic scattering of photons and neutrons,
  // respectively.
  if (ParticleName == "neutron" || ParticleName == "gamma") {
    /* G4cout << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
     << "\t"; G4cout <<
     step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()  <<
     "\t"; G4cout << ParticleName<< "\t" << step->GetDeltaEnergy()  << "\t";
     G4cout <<   step->GetPostStepPoint()->GetTotalEnergy() << "\t";
     G4cout <<   step->GetPreStepPoint()->GetKineticEnergy()  << "\t";
     G4cout <<   step->GetPostStepPoint()->GetKineticEnergy() << "\n";*/
    if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() ==
        "compt") {
      // G4cout << "compt" << G4endl;
      // PrimaryDeltaE =  step->GetPostStepPoint()->GetTotalEnergy();
      G4double deng = step->GetPreStepPoint()->GetTotalEnergy() -
                      step->GetPostStepPoint()->GetTotalEnergy();
      crEvent->recordReaction(0, deng,
                              step->GetPostStepPoint()->GetGlobalTime(),
                              step->GetPreStepPoint()->GetTotalEnergy());
      /*G4double ekin = endPoint->GetKineticEnergy();
      G4double trackl = step->GetTrack()->GetTrackLength();
      G4double time = step->GetTrack()->GetLocalTime();
      fTrackingAction->UpdateTrackInfo(ekin, trackl, time);
      G4AnalysisManager::Instance()->FillH1(7, ekin);*/
    }
    if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() ==
        "hadElastic") {
      // G4cout << "hadElastic" << G4endl;
      G4double deng = step->GetPreStepPoint()->GetKineticEnergy() -
                      step->GetPostStepPoint()->GetKineticEnergy();
      crEvent->recordReaction(1, deng,
                              step->GetPostStepPoint()->GetGlobalTime(),
                              step->GetPreStepPoint()->GetKineticEnergy());

      /*G4double ekin = endPoint->GetKineticEnergy();
      G4double trackl = step->GetTrack()->GetTrackLength();
      G4double time = step->GetTrack()->GetLocalTime();
      fTrackingAction->UpdateTrackInfo(ekin, trackl, time);
      G4AnalysisManager::Instance()->FillH1(7, ekin);*/
    }
  }
  if (ParticleName == "e-") {
    /* G4cout << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
     << "\t";

     G4cout <<
     step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()  <<
     "\t";

     G4cout << ParticleName<< "\t" << step->GetDeltaEnergy()  << "\t";

     G4cout <<  step->GetPostStepPoint()->GetTotalEnergy() << "\t";
     G4cout <<   step->GetTotalEnergyDeposit() << "\t";
     G4cout <<   step->GetPreStepPoint()->GetKineticEnergy()  << "\t";
     G4cout <<   step->GetPostStepPoint()->GetKineticEnergy() << "\n";*/
    if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() ==
        "eIoni") {
      G4double deng = step->GetPreStepPoint()->GetKineticEnergy() -
                      step->GetPostStepPoint()->GetKineticEnergy();
      crEvent->recordReaction(2, deng);
      /*if (fElectronCounter == 0 || update)
       {
       fElectronEnergy += step->GetPreStepPoint()->GetKineticEnergy();
       update = false;
       }
       fElectronCounter++;*/
    }
  }
  if (ParticleName == "proton") {
    /*G4cout << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
     << "\t";

     G4cout <<
     step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()  <<
     "\t";

     G4cout << ParticleName<< "\t" << step->GetDeltaEnergy()  << "\t";

     G4cout <<   step->GetPostStepPoint()->GetTotalEnergy() << "\t";
     G4cout <<   step->GetTotalEnergyDeposit() << "\t";
     G4cout <<   step->GetPreStepPoint()->GetKineticEnergy()  << "\t";
     G4cout <<  step->GetPostStepPoint()->GetKineticEnergy() << "\n";*/
    if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() ==
        "hIoni") {
      // if (step->GetPreStepPoint()->GetKineticEnergy() ==
      // step->GetPostStepPoint()->GetKineticEnergy() ) return;
      G4double deng = step->GetPreStepPoint()->GetKineticEnergy() -
                      step->GetPostStepPoint()->GetKineticEnergy();
      crEvent->recordReaction(3, deng);
    }
  }

  const std::vector<const G4Track*>* secondaries =
      step->GetSecondaryInCurrentStep();
  // G4cout << secondaries->size() << "\n";
  if (secondaries->size() > 0) {
    // fout << "\n " <<
    // G4OpticalPhoton::OpticalPhotonDefinition()->GetParticleName()<< "\t";
    for (unsigned int i = 0; i < secondaries->size(); ++i) {
      if (secondaries->at(i)->GetParentID() > 0) {
        if (secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() ==
            G4OpticalPhoton::OpticalPhotonDefinition()) {
          if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() ==
              "Scintillation") {
            crEvent->recordProduction(1,
                                      step->GetPostStepPoint()->GetGlobalTime(),
                                      secondaries->at(i)->GetKineticEnergy());
          }
          if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() ==
              "Cerenkov") {
            crEvent->recordProduction(5,
                                      secondaries->at(i)->GetKineticEnergy());
          }
        } else if (secondaries->at(i)
                           ->GetDynamicParticle()
                           ->GetParticleDefinition()
                           ->GetParticleName() == "neutron" &&
                   ParticleName != "neutron" &&
                   track->GetDynamicParticle()
                           ->GetParticleDefinition()
                           ->GetParticleName() != "proton") {
          // G4cout << secondaries->at(i)->GetVolume()->GetName()<< " " <<
          // ParticleName << " " <<
          // secondaries->at(i)->GetCreatorProcess()->GetProcessName()  <<
          // G4endl; G4cout <<
          // secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()->GetParticleName()<<
          // "\n";

          crEvent->recordProduction(3, secondaries->at(i)->GetKineticEnergy());
        } else if (secondaries->at(i)
                           ->GetDynamicParticle()
                           ->GetParticleDefinition()
                           ->GetParticleName() == "proton" &&
                   ParticleName != "proton" &&
                   track->GetDynamicParticle()
                           ->GetParticleDefinition()
                           ->GetParticleName() != "proton") {
          // G4cout << secondaries->at(i)->GetVolume()->GetName()<< " " <<
          // ParticleName << " " <<
          // secondaries->at(i)->GetCreatorProcess()->GetProcessName()  <<
          // G4endl; G4cout <<
          // secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()->GetParticleName()<<
          // "\n";
          crEvent->recordProduction(2, secondaries->at(i)->GetKineticEnergy());
        } else if (secondaries->at(i)
                           ->GetDynamicParticle()
                           ->GetParticleDefinition()
                           ->GetParticleName() == "e-" &&
                   ParticleName != "e-" &&
                   track->GetDynamicParticle()
                           ->GetParticleDefinition()
                           ->GetParticleName() != "e-") {
          crEvent->recordProduction(4, secondaries->at(i)->GetKineticEnergy());
        } else if (secondaries->at(i)
                       ->GetDynamicParticle()
                       ->GetParticleDefinition()
                       ->GetParticleName() == "gamma") {
          /*      G4cout <<
             G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << "
             Gamma detected " << secondaries->at(i)->GetKineticEnergy()/MeV << "
             " << secondaries->at(i)->GetParentID() << G4endl; G4cout << "Parent
             " <<  ParticleName << "\t" <<
             step->GetPostStepPoint()->GetTotalEnergy()/keV
                        << "MeV \t"<<
             step->GetPreStepPoint()->GetKineticEnergy()/keV  << " keV \t" <<
                          step->GetPostStepPoint()->GetKineticEnergy() << " keV
             \t"
                          <<
             step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
             << " keV \t"<< G4endl;
                //secondaries->at(i)->SetTrackStatus(fStopAndKill);
                G4cout <<
             step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
             << " " <<
                        step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessType()
             << " " <<
             step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessSubType()
             << G4endl;
                */
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

