#include "Trajectory.hh"
#include "G4TrajectoryPoint.hh"
#include "G4Trajectory.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ThreeVector.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4Polymarker.hh"
#include "globals.hh"

G4ThreadLocal G4Allocator<Trajectory>* TrajectoryAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Trajectory::Trajectory()
  :G4Trajectory(),fWls(false),fDrawit(false),fForceNoDraw(false),fForceDraw(false)
{
  fParticleDefinition=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Trajectory::Trajectory(const G4Track* aTrack)
  :G4Trajectory(aTrack),fWls(false),fDrawit(false)
{
  fParticleDefinition=aTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Trajectory::Trajectory(Trajectory &right)
  :G4Trajectory(right),fWls(right.fWls),fDrawit(right.fDrawit)
{
  fParticleDefinition=right.fParticleDefinition;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Trajectory::~Trajectory() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Trajectory::DrawTrajectory() const
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (!pVVisManager) return;

  G4Polyline trajectoryLine;

  for (G4int i = 0; i < GetPointEntries() ; i++) {
    G4VTrajectoryPoint* aTrajectoryPoint = GetPoint(i);
    const std::vector<G4ThreeVector>* auxiliaries
      = aTrajectoryPoint->GetAuxiliaryPoints();
    if (auxiliaries) {
      for (size_t iAux = 0; iAux < auxiliaries->size(); ++iAux) {
        const G4ThreeVector pos((*auxiliaries)[iAux]);
          trajectoryLine.push_back(pos);
      }
    }
    const G4ThreeVector pos(aTrajectoryPoint->GetPosition());
      trajectoryLine.push_back(pos);
  }

    G4Colour colour;

    if(fParticleDefinition == "opticalphoton")
		colour = G4Colour(1., 0., 0., 100.); //red
    else if(fParticleDefinition == "neutron")
        colour = G4Colour(0.,0.,1.); // blue
    else if(fParticleDefinition == "proton")
        colour = G4Colour(1.,0.,1); //grey
    else if(fParticleDefinition == "gamma")
        colour = G4Colour(0.,1.,0.); //lime
    else if(fParticleDefinition == "e-")
        colour = G4Colour(1.,1.,0.); // yellow
    else
	{
		//G4cout << fParticleDefinition << G4endl;
		colour = G4Colour(1., 1., 1., 1.);
	}
    G4VisAttributes trajectoryLineAttribs(colour);
    trajectoryLine.SetVisAttributes(&trajectoryLineAttribs);
    pVVisManager->Draw(trajectoryLine);

}
