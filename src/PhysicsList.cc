
#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "HadronElasticPhysicsHP.hh"
#include "NeutronHPphysics.hh"

#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"

#include "G4IonINCLXXPhysics.hh"
#include "G4IonPhysics.hh"

#include "G4StoppingPhysics.hh"
#include "GammaNuclearPhysics.hh"

#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4DecayPhysics.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4LossTableManager.hh"
#include "G4MesonConstructor.hh"
#include "G4OpticalPhysics.hh"
#include "G4ProcessManager.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4ShortLivedConstructor.hh"
#include "PhysListEmStandard.hh"

#include "G4Cerenkov.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpMieHG.hh"
#include "G4OpRayleigh.hh"
#include "G4Scintillation.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EMPhysics.hh"
#include "GeneralPhysics.hh"
#include "MuonPhysics.hh"

G4ThreadLocal G4int PhysicsList::fMaxNumPhotonStep = 20;
G4ThreadLocal G4Cerenkov* PhysicsList::fCerenkovProcess = 0;
G4ThreadLocal G4Scintillation* PhysicsList::fScintillationProcess = 0;
G4ThreadLocal G4OpAbsorption* PhysicsList::fAbsorptionProcess = 0;
G4ThreadLocal G4OpRayleigh* PhysicsList::fRayleighScatteringProcess = 0;
G4ThreadLocal G4OpMieHG* PhysicsList::fMieHGScatteringProcess = 0;
G4ThreadLocal G4OpBoundaryProcess* PhysicsList::fBoundaryProcess = 0;

PhysicsList::PhysicsList() : G4VModularPhysicsList() {
  G4int verb = 0;
  SetVerboseLevel(verb);

  // add new units
  //
  new G4UnitDefinition("millielectronVolt", "meV", "Energy", 1.e-3 * eV);
  new G4UnitDefinition("mm2/g", "mm2/g", "Surface/Mass", mm2 / g);
  new G4UnitDefinition("um2/mg", "um2/mg", "Surface/Mass", um * um / mg);

  // Neutron Physics
  RegisterPhysics(new NeutronHPphysics("neutronHP"));

  // RegisterPhysics(new HadronElasticPhysicsHP(verb));

  RegisterPhysics(new G4HadronPhysicsQGSP_BIC_HP(verb));

  // Ion Physics
  RegisterPhysics(new G4IonPhysics(verb));
  ////RegisterPhysics( new G4IonINCLXXPhysics(verb));

  // stopping Particles
  RegisterPhysics(new G4StoppingPhysics(verb));

  // Gamma-Nuclear Physics
  // EM physics
  RegisterPhysics(new ElectromagneticPhysics());

  // Decay
  // RegisterPhysics(new G4DecayPhysics());

  // Radioactive decay
  // RegisterPhysics(new G4RadioactiveDecayPhysics());

  defaultCutValue = 1.0 * mm;

  // General Physics
  // RegisterPhysics(new GeneralPhysics("general"));

  // EM Physics
  RegisterPhysics(new EMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(new MuonPhysics("muon"));

  // Optical Physics
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  RegisterPhysics(opticalPhysics);

  opticalPhysics->SetWLSTimeProfile("delta");
  opticalPhysics->SetScintillationByParticleType(true);
  opticalPhysics->SetScintillationYieldFactor(1);
  opticalPhysics->SetScintillationExcitationRatio(0.0);

  opticalPhysics->SetTrackSecondariesFirst(kCerenkov, true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation, true);  //*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle() {
  G4BosonConstructor pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts() {
  // SetCutValue(0 * mm, "proton");
  /*SetCutValue(10 * km, "e-");
   SetCutValue(10 * km, "e+");
   SetCutValue(10 * km, "gamma");*/
}
/*
void PhysicsList::ConstructProcess()
{
        AddTransportation();
        ConstructDecay();
        ConstructEM();
        ConstructOp();
 }*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
#include "G4Decay.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructDecay()
{
        // Add Decay Process
        G4Decay* theDecayProcess = new G4Decay();
        //auto particleIterator = GetParticleIterator();
        aparticleIterator->reset();
        while ((*particleIterator)())
        {
                G4ParticleDefinition* particle = particleIterator->value();
                G4ProcessManager* pmanager = particle->GetProcessManager();
                if (theDecayProcess->IsApplicable(*particle))
                {
                        pmanager->AddProcess(theDecayProcess);
                        // set ordering for PostStepDoIt and AtRestDoIt
                        pmanager->SetProcessOrdering(theDecayProcess,
idxPostStep); pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
                }
        }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MuMultipleScattering.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuBremsstrahlung.hh"
#include "G4MuIonisation.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructEM()
{
//	auto particleIterator = GetParticleIterator();
        particleIterator->reset();
        while ((*particleIterator)())
        {
                G4ParticleDefinition* particle = particleIterator->value();
                G4ProcessManager* pmanager = particle->GetProcessManager();
                G4String particleName = particle->GetParticleName();

                if (particleName == "gamma")
                {
                        // gamma
                        // Construct processes for gamma
                        pmanager->AddDiscreteProcess(new G4GammaConversion());
                        pmanager->AddDiscreteProcess(new G4ComptonScattering());
                        pmanager->AddDiscreteProcess(new
G4PhotoElectricEffect());

                }
                else if (particleName == "e-")
                {
                        //electron
                        // Construct processes for electron
                        pmanager->AddProcess(new G4eMultipleScattering(), -1, 1,
1); pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);
                        pmanager->AddProcess(new G4eBremsstrahlung(), -1, 3, 3);

                }
                else if (particleName == "e+")
                {
                        //positron
                        // Construct processes for positron
                        pmanager->AddProcess(new G4eMultipleScattering(), -1, 1,
1); pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);
                        pmanager->AddProcess(new G4eBremsstrahlung(), -1, 3, 3);
                        pmanager->AddProcess(new G4eplusAnnihilation(), 0, -1,
4);

                }
                else if (particleName == "mu+" || particleName == "mu-")
                {
                        //muon
                        // Construct processes for muon
                        pmanager->AddProcess(new G4MuMultipleScattering(), -1,
1, 1); pmanager->AddProcess(new G4MuIonisation(), -1, 2, 2);
                        pmanager->AddProcess(new G4MuBremsstrahlung(), -1, 3,
3); pmanager->AddProcess(new G4MuPairProduction(), -1, 4, 4);

                }
                else
                {
                        if ((particle->GetPDGCharge() != 0.0)
                                        && (particle->GetParticleName() !=
"chargedgeantino")
                                        && !particle->IsShortLived())
                        {
                                // all others charged particles except geantino
                                pmanager->AddProcess(new
G4hMultipleScattering(), -1, 1, 1); pmanager->AddProcess(new G4hIonisation(),
-1, 2, 2);
                        }
                }
        }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4Threading.hh"

void PhysicsList::ConstructOp()
{
        int fVerboseLevel = 0;
        fCerenkovProcess = new G4Cerenkov("Cerenkov");
        fCerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotonStep);
        fCerenkovProcess->SetMaxBetaChangePerStep(10.0);
        fCerenkovProcess->SetTrackSecondariesFirst(true);
        fScintillationProcess = new G4Scintillation("Scintillation");
        fScintillationProcess->SetScintillationYieldFactor(1.);
        fScintillationProcess->SetTrackSecondariesFirst(true);
        fAbsorptionProcess = new G4OpAbsorption();
        fRayleighScatteringProcess = new G4OpRayleigh();
        fMieHGScatteringProcess = new G4OpMieHG();
        fBoundaryProcess = new G4OpBoundaryProcess();

        fCerenkovProcess->SetVerboseLevel(fVerboseLevel);
        fScintillationProcess->SetVerboseLevel(fVerboseLevel);
        fAbsorptionProcess->SetVerboseLevel(fVerboseLevel);
        fRayleighScatteringProcess->SetVerboseLevel(fVerboseLevel);
        fMieHGScatteringProcess->SetVerboseLevel(fVerboseLevel);
        fBoundaryProcess->SetVerboseLevel(fVerboseLevel);

        // Use Birks Correction in the Scintillation process
        if (G4Threading::IsMasterThread())
        {
                G4EmSaturation* emSaturation =
                                G4LossTableManager::Instance()->EmSaturation();
                fScintillationProcess->AddSaturation(emSaturation);
        }

//	auto particleIterator = GetParticleIterator();
        particleIterator->reset();
        while ((*particleIterator)())
        {
                G4ParticleDefinition* particle = particleIterator->value();
                G4ProcessManager* pmanager = particle->GetProcessManager();
                G4String particleName = particle->GetParticleName();
                if (fCerenkovProcess->IsApplicable(*particle))
                {
                        pmanager->AddProcess(fCerenkovProcess);
                        pmanager->SetProcessOrdering(fCerenkovProcess,
idxPostStep);
                }
                if (fScintillationProcess->IsApplicable(*particle))
                {
                        pmanager->AddProcess(fScintillationProcess);
                        pmanager->SetProcessOrderingToLast(fScintillationProcess,
                                        idxAtRest);
                        pmanager->SetProcessOrderingToLast(fScintillationProcess,
                                        idxPostStep);
                }
                if (particleName == "opticalphoton")
                {
                        G4cout << " AddDiscreteProcess to OpticalPhoton " <<
G4endl; pmanager->AddDiscreteProcess(fAbsorptionProcess);
                        pmanager->AddDiscreteProcess(fRayleighScatteringProcess);
                        pmanager->AddDiscreteProcess(fMieHGScatteringProcess);
                        pmanager->AddDiscreteProcess(fBoundaryProcess);
                }
        }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
