//******************************************************************************
// SponFiss.cc
//
//******************************************************************************
//
#include "SponFIss_FIF.hh"
#include <cmath>
#include <string>
#include <vector>
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4INCLThreeVector.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "SingleSource.hh"

std::ifstream* SponFiss_FF::bfile = NULL;
std::string SponFiss_FF::filename = "dump.bin";
//----------------------------------------------------------------------------//
SponFiss_FF::SponFiss_FF(G4SPSPosDistribution* pos) {
  neutron_definition = G4Neutron::Neutron();
  photon_definition = G4Gamma::Gamma();

  // verbosity
  verbosityLevel = 0;
  // fRun = run;
  posDist = pos;

  if (bfile == NULL) {
    // bfile = &DetectorConstruction::bfile;
    bfile = new std::ifstream(filename, std::ios::in | std::ios::binary);

    if (!bfile->is_open()) {
      G4cout << "FIF Failed: " << filename << G4endl;
      exit(-1);
    }
    G4cout << "FIF initialized" << G4endl;
  }
}

//----------------------------------------------------------------------------//
SponFiss_FF::~SponFiss_FF() {}

void SponFiss_FF::update_angular(std::vector<G4ThreeVector>& fe) {
  return;
  for (int n1 = 0; n1 < fe.size(); n1++) {
    double u1 = fe[n1].getX(), v1 = fe[n1].getY(), w1 = fe[n1].getZ();
    for (int n2 = n1 + 1; n2 < fe.size(); n2++) {
      double u2 = fe[n2].getX(), v2 = fe[n2].getY(), w2 = fe[n2].getZ();
      double scalar_product = u1 * u2 + v1 * v2 + w1 * w2;
      // cout << scalar_product << endl;
      int bin_index = (int)(98 * (scalar_product + 1) / 2);
      DetectorConstruction::hist[0][bin_index]++;
      if (n2 - n1 < 10) DetectorConstruction::hist[n2 - n1][bin_index]++;
    }
  }
}
//----------------------------------------------------------------------------//
void SponFiss_FF::GeneratePrimaryVertex(G4Event* anEvent, G4double time,
                                        int mode) {
retry:
  if (bfile->eof()) {
    bfile->clear();
    bfile->seekg(0, std::ios::beg);
    return;
  }
  // Generate a spontaneous fission using the fission library and emit
  // the neutrons and gamma-rays

  short nPrompt = 0;
  short* _nPrompt = new short[1];

  bfile->read((char*)_nPrompt, sizeof(char) * 2);
  if (bfile->eof()) goto retry;
  // fRun->fRun->UpdateSource(nPrompt, gPrompt);

  G4double rota = (G4UniformRand() * 180 * deg);

  nPrompt = *_nPrompt;
  if (verbosityLevel > 1) {
    G4cout << " nPrompt: " << nPrompt << G4endl;
  }
  delete _nPrompt;
  if (nPrompt < 1) return;

  std::vector<G4ThreeVector> dir;
  // Position
  G4ThreeVector sampled_particle_position =
      DetectorConstruction::source_position;
  // create a new vertex
  G4PrimaryVertex* vertex = new G4PrimaryVertex(sampled_particle_position, 0.);

  G4double mom, momx, momy, momz, eng;

  if (verbosityLevel >= 2)
    G4cout << "Creating primaries and assigning to vertex" << G4endl;
  G4DynamicParticle* it;
  // Build neutrons

  if (PrimaryGeneratorAction::neutron && nPrompt > 0)
    for (G4int i = nPrompt - 1; i > -1; i--) {
      char* tmp = new char[4 * 4];
      bfile->read((char*)tmp, sizeof(float) * 4);
      if (bfile->eof()) goto retry;
      float* fl = (float*)tmp;

      eng = fl[0] / 1000;

      G4INCL::ThreeVector tmp_rot = G4INCL::ThreeVector(fl[1], fl[2], fl[3]);
      tmp_rot.rotate(rota, G4INCL::ThreeVector(0, 0, 1));
      momx = tmp_rot.getX();
      momy = tmp_rot.getY();
      momz = tmp_rot.getZ();

      dir.push_back(G4ThreeVector(momx, momy, momz));
      // continue;
      if (eng > 19.9) eng = 19;

      it = new G4DynamicParticle();
      it->SetDefinition(neutron_definition);
      it->SetKineticEnergy(eng);
      mom = it->GetTotalMomentum();
      delete it;

      delete[] tmp;

      G4PrimaryParticle* particle = new G4PrimaryParticle(
          neutron_definition, mom * momx, mom * momy, mom * momz, eng * MeV);
      particle->SetMass(neutron_definition->GetPDGMass());
      particle->SetCharge(neutron_definition->GetPDGCharge());
      particle->SetPolarization(particle_polarization.x(),
                                particle_polarization.y(),
                                particle_polarization.z());

      if (verbosityLevel > 1) {
        G4cout << "Particle name: " << particle->GetG4code()->GetParticleName()
               << G4endl;
        G4cout << "     Momentum: " << particle->GetMomentum() << G4endl;
        G4cout << "     Position: " << vertex->GetPosition() << G4endl;
        G4cout << "     Energy:  " << particle->GetKineticEnergy() << G4endl;
      }
      vertex->SetPrimary(particle);
    }
  anEvent->AddPrimaryVertex(vertex);
  // if (dir.size() > 2) update_angular(dir);
  if (dir.size() > 2)
    for (int n1 = 0; n1 < dir.size(); n1++) {
      double u1 = dir.at(n1).getX(), v1 = dir.at(n1).getY(),
             w1 = dir.at(n1).getZ();
      for (int n2 = n1 + 1; n2 < dir.size(); n2++) {
        double u2 = dir.at(n2).getX(), v2 = dir.at(n2).getY(),
               w2 = dir.at(n2).getZ();

        double scalar_product = u1 * u2 + v1 * v2 + w1 * w2;

        int bin_index = (int)(98 * (scalar_product + 1) / 2);
        DetectorConstruction::hist[0][bin_index]++;
        if (n2 - n1 < 10) DetectorConstruction::hist[n2 - n1][bin_index]++;
      }
    }

  if (verbosityLevel > 1) G4cout << " Primary Vetex generated !" << G4endl;
}
