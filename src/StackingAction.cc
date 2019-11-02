//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file StackingAction.cc
/// \brief Implementation of the StackingAction class
//
// $Id: StackingAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"
#include "MyRun.hh"

#include "G4RunManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "TrackingInformation.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
:
		G4UserStackingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(
		const G4Track* aTrack)
{
	//keep primary particle
	if (aTrack->GetParentID() == 0)
	{
		TrackInformation* trackInfo;
		trackInfo = new TrackInformation(aTrack);
		trackInfo->SetTrackingStatus(1);
		trackInfo->fParentID = 0;
		G4Track* theTrack = (G4Track*)aTrack;
		theTrack->SetUserInformation(trackInfo);
		return fUrgent;
	}

	//count secondary particles
	G4String name = aTrack->GetDefinition()->GetParticleName();
	G4double energy = aTrack->GetKineticEnergy();
	Run* run =
			static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
	run->ParticleCount(name, energy);

	//kill all secondaries
	if (aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
	{ // particle is optical photon
		if (aTrack->GetParentID() > 0)
		{ // particle is secondary
			if (aTrack->GetCreatorProcess()->GetProcessName()
					== "Scintillation")
				fScintillationCounter++;
			if (aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov")
				fCerenkovCounter++;
			return fWaiting;
		}
	}
	if (aTrack->GetDefinition()->GetParticleName() == "gamma")
	{ // particle is optical photon
		if (aTrack->GetParentID() > 0)
		{ // particle is secondary
			if (aTrack->GetKineticEnergy()/keV < 10)
				return fKill;
			else if (aTrack->GetCreatorProcess()->GetProcessName()
					== "nCapture")
			{
			//if (aTrack->GetVolume()->GetName() == "Scintillator")
				G4Track* theTrack = (G4Track*)aTrack;
				//theTrack->SetParentID(0);
				TrackInformation* trackInfo
				= (TrackInformation*)(aTrack->GetUserInformation());
				trackInfo->SetTrackingStatus(1);
				theTrack->SetUserInformation(trackInfo);
				return fUrgent;
			}
			else 
			{
				//G4cout << " Stacking " << aTrack->GetCreatorProcess()->GetProcessName() << " " <<  aTrack->GetKineticEnergy() << G4endl;
				G4Track* theTrack = (G4Track*)aTrack;
				//theTrack->SetParentID(0);
				TrackInformation* trackInfo
				= (TrackInformation*)(aTrack->GetUserInformation());
				trackInfo->SetTrackingStatus(2);
				theTrack->SetUserInformation(trackInfo);
				return fUrgent;
			}
		}
	}
	if (energy*MeV > 20*MeV) return fKill; 
	return fUrgent;
}

void StackingAction::NewStage()
{
	//G4cout << "Number of Scintillation photons produced in this event : "
	//       << fScintillationCounter << G4endl;
	//G4cout << "Number of Cerenkov photons produced in this event : "
	//       << fCerenkovCounter << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingAction::PrepareNewEvent()
{
	fScintillationCounter = 0;
	fCerenkovCounter = 0;
	fCount = 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
