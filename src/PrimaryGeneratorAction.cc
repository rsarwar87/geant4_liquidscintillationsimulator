
#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int PrimaryGeneratorAction::energy = 500;
int PrimaryGeneratorAction::mode = 3;
double PrimaryGeneratorAction::decay_time = 0;
G4String PrimaryGeneratorAction::name = "neutron";
bool PrimaryGeneratorAction::gamma = false;
bool PrimaryGeneratorAction::neutron = true;
bool PrimaryGeneratorAction::mono = false;
bool PrimaryGeneratorAction::Co = false;
bool PrimaryGeneratorAction::beam = false;
bool PrimaryGeneratorAction::AmLi = false;
bool PrimaryGeneratorAction::sfif = false;
G4Mutex PrimaryGeneratorAction::aMutex = G4MUTEX_INITIALIZER;

std::discrete_distribution<> PrimaryGeneratorAction::intensity_dis;
std::discrete_distribution<> PrimaryGeneratorAction::leak_dis;

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {
  std::vector<double> intensity = {.965, .035};
  std::vector<double> leak = {.994, .006};

  fParticleGun = new G4ParticleGun(1);
  fParticleGun->SetParticleDefinition(
      G4ParticleTable::GetParticleTable()->FindParticle(name));
  fParticleGun->SetParticleEnergy(energy * keV);
  fParticleGun->SetParticlePosition(DetectorConstruction::source_position);
  fParticleGun->SetParticleTime(0.0 * ns);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  G4cout << DetectorConstruction::source_position << G4endl;
  std::discrete_distribution<> tmp_leak(leak.begin(), leak.end());
  std::discrete_distribution<> tmp_intensity(intensity.begin(),
                                             intensity.end());

  intensity_dis = tmp_intensity;
  leak_dis = tmp_leak;

  // Specify isotopic composition and fission rates in fissions/sec
  time = 0;  // set to 0 initially
  G4ThreeVector* center = &DetectorConstruction::source_position;
  G4double radius = .2 * cm;

  // Specify isotopic composition and fission rates in fissions/sec
  G4DataVector* isotopeList = new G4DataVector(2);
  isotopeList->operator[](0) = 92238;

  G4DataVector* intensityList = new G4DataVector(2);
  intensityList->operator[](0) = 2.368;

  time = 0;  // set to 0 initially

  SponFissIsotope* currentIsotope;
  G4SPSPosDistribution* posDist2;

  totalIntensity = 0.;
  nisotopes = isotopeList->size();
  for (G4int i = 0; i < nisotopes; i++) {
    currentIsotope =
        new SponFissIsotope(static_cast<G4int>(isotopeList->operator[](i)));
    posDist2 = currentIsotope->GetPosDist();
    posDist2->SetPosDisType("Volume");
    posDist2->SetPosDisShape("Sphere");
    posDist2->SetCentreCoords(*center);
    posDist2->SetRadius(radius);
    totalIntensity += intensityList->operator[](i);
    if (i == 0)
      fissionSource =
          new MultipleSource(currentIsotope, intensityList->operator[](i));
    else
      fissionSource->AddaSource(currentIsotope, intensityList->operator[](i));
    fissionSource->SetVerbosity(2);
  }

  posDist = new G4SPSPosDistribution();
  posDist->SetPosDisType("Volume");
  posDist->SetPosDisShape("Sphere");
  posDist->SetCentreCoords(*center);
  posDist->SetRadius(radius);
  iso = new SponFiss(98252, posDist);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fParticleGun;
  delete posDist;
  delete iso;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  if (mono) {
    // fParticleGun->SetParticleEnergy((G4UniformRand()+.3 ) * MeV);
    fParticleGun->SetParticleEnergy(energy * keV);
    if (!beam) {
      G4ThreeVector direction;
      direction.setRThetaPhi(1.0, std::acos(G4UniformRand() * 2 - 1),
                             (G4UniformRand() * 2 - 1) * 180 * deg);
      fParticleGun->SetParticleMomentumDirection(direction);
    }
    fParticleGun->GeneratePrimaryVertex(anEvent);
  } else if (Co) {
    fParticleGun->SetParticleEnergy(1121 * keV);
    G4ThreeVector direction;
    direction.setRThetaPhi(1.0, std::acos(G4UniformRand() * 2 - 1),
                           (G4UniformRand() * 2 - 1) * 180 * deg);
    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    fParticleGun->SetParticleEnergy(1333 * keV);
    direction.setRThetaPhi(1.0, std::acos(G4UniformRand() * 2 - 1),
                           (G4UniformRand() * 2 - 1) * 180 * deg);
    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->GeneratePrimaryVertex(anEvent);
  } else if (AmLi) {
    fParticleGun->SetParticleEnergy((G4UniformRand() + .3) * MeV);
    G4ThreeVector direction;
    direction.setRThetaPhi(1.0, std::acos(G4UniformRand() * 2 - 1),
             (G4UniformRand() * 2 - 1) * 180 * deg);
    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->GeneratePrimaryVertex(anEvent);
  } else if (sfif) {
    G4MUTEXLOCK(&aMutex);
    static SponFiss_FF* fif = new SponFiss_FF(posDist);
    fif->GeneratePrimaryVertex(anEvent);
    G4MUTEXUNLOCK(&aMutex);
  } else {
    decay_time += 1 / 331000;
    G4MUTEXLOCK(&aMutex);
    iso->GeneratePrimaryVertex(anEvent, decay_time * ns, mode);
    G4MUTEXUNLOCK(&aMutex);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
