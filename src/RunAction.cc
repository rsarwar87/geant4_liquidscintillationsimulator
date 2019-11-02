

#include "RunAction.hh"
#include "MyRun.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <iomanip>
#include <vector>
#include <string.h>

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double RunAction::cutOu = 50;
RunAction::RunAction(G4String fn, DetectorConstruction* det,
		PrimaryGeneratorAction* prim) :
		G4UserRunAction(), fDetector(det), fPrimary(prim), fRun(0), fHistoManager(
				0)
{
	// Book predefined histograms
	fHistoManager = new HistoManager();
	fTimer = new G4Timer;
	filename = fn;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
	delete fTimer;
	delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{
	fRun = new Run(fDetector);
	return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::multiplicity_spec(const G4Run* aRun)
{
	std::ofstream fout;
	fout.open(filename + "_spec_" + std::to_string(PrimaryGeneratorAction::energy), 
						std::ios::out | std::ios::trunc);
	int detected = fRun->GetTotalRegisterPartile(0);
	fout << "\n\n\n\n number of event = " << aRun->GetNumberOfEvent() << "\n";
	fout << "\n number of detected event = " << detected;
	
	
	for (int pid = 0; pid < 2; pid++)	
	{
		if (!PrimaryGeneratorAction::neutron && (pid == 1)) continue;
		if (!PrimaryGeneratorAction::gamma && (pid == 0)) continue;
		int sp = 16;
		fout << ((pid==0) ? "\n\n\n\nGamma" : "\n\nNeutron" ) <<" Distribution \n";
		fout << std::setw(sp) << "Energy [keV]" << "\t"
					<< std::setw(sp) << "Incident [keV]" << "\t" 
					<< std::setw(sp) << "EntryEnergy [keV]" << "\t" 
					<< std::setw(sp) << "Particle Dep [kev]" << "\t"
					<< std::setw(sp) << "Proton Dep [kev]" << "\t" 
					<< std::setw(sp) << "Electron Dep [kev]" << "\t" 
					<< std::setw(sp) << "Optical Dep [kev]" << "\t"
					<< std::setw(sp) << "Time [ns]" << "\t" 
					<< std::setw(sp) << "LightRes" << "\t" 
					<< std::setw(sp) << "LightHistogram" << "\t"
					<< std::setw(sp) << "PMTLResponse" << "\t"
					<< std::setw(sp) << "PMTHistogaram\n" ;
			
			const G4int *fIncidentEnergy = fRun->GetParticleIncidentEnergy(pid);
			const G4int *fEntryEnergy = fRun->GetParticleEntryEnergy(pid);
			const G4int *fParticleDeposit = fRun->GetParticleDeposit(pid);
			const G4int *fElectronDeposited = fRun->GetElectronDeposited(pid);
			const G4int *fProtonDeposited = fRun->GetProtonDeposited(pid);
			const G4int *fOphotonDeposited = fRun->GetPMTDistribution(pid);
			
			const G4int *fLightResponse = fRun->GetLightResponse(pid);
			const G4int *fLightHistogram = fRun->GetLightHistogram(pid);
			const G4int *fPMTResponse = fRun->GetPMTResponse(pid);
			const G4int *fPMTHistogram = fRun->GetPMTHistogram(pid);
			
			for (uint i = 0; i < 4500; i++)
				fout << std::setw(sp) << i*2+1 << "\t" 
					<< std::setw(sp) << fIncidentEnergy[i] << "\t" 
					<< std::setw(sp) << fEntryEnergy[i]  << "\t" 
					<< std::setw(sp) << fParticleDeposit[i] << "\t"
					<< std::setw(sp) << fProtonDeposited[i] << "\t" 
					<< std::setw(sp) << fElectronDeposited[i] << "\t" 
					<< std::setw(sp) << fOphotonDeposited[i] << "\t"
					/*<< std::setw(sp) << i << "\t" 
					<< std::setw(sp) << fLightResponse[i] << "\t" 
					<< std::setw(sp) << fLightHistogram[i] << "\t"
					<< std::setw(sp) << fPMTResponse[i] << "\t"
					<< std::setw(sp) << fPMTHistogram[i]*/ 
					<< "\n" ;
			
	}
	fout.close();
	G4cout << "printing theretical specs \n\n" << G4endl;
}
void RunAction::multiplicity_th()
{
	std::ofstream fd;
	fd.open(filename + "_coin_" + std::to_string(PrimaryGeneratorAction::energy), 
						std::ios::out | std::ios::trunc);
	{
		const G4int * g_th = fRun->GetGammaSourceMultiplicity();
		const G4int * n_th = fRun->GetNeutronSourceMultiplicity();
		const G4int * j_th = fRun->GetJointSourceMultiplicity();
		fd << "Theoretical Multiplicity: Gamma\n";
		for (int i = 0; i < 32; i++)
			fd << "\t" << g_th[i];
		fd << "\nTheoretical Multiplicity: Neutron\n";
		for (int i = 0; i < 32; i++)
			fd << "\t" << n_th[i];
		fd << "\nTheoretical Multiplicity: Joint\n";
		for (int i = 0; i < 64; i++)
			fd << "\t" << j_th[i];
	}
	{
		const G4int *g_th = fRun->GetGammaDetectedMultiplicity();
		const G4int *n_th = fRun->GetNeutronDetectedMultiplicity();
		const G4int *j_th = fRun->GetJointDetectedMultiplicity();

		const G4int *g_ac = fRun->GetGammaDetectedCXMultiplicity();
		const G4int *n_ac = fRun->GetNeutronDetectedCXMultiplicity();
		const G4int *j_ac = fRun->GetJointDetectedCXMultiplicity();
		
		fd << "\n\nDetected Multiplicity: Gamma\n";
		for (int i = 0; i < 32; i++)
			fd << "\t" << g_th[i];
		fd << "\nDetected Multiplicity: Neutron\n";
		for (int i = 0; i < 32; i++)
			fd << "\t" << n_th[i];
		fd << "\nDetected Multiplicity: Joint\n";
		for (int i = 0; i < 64; i++)
			fd << "\t" << j_th[i];
		
		fd << "\n\nXT Corrected Multiplicity: Gamma\n";
		for (int i = 0; i < 32; i++)
			fd << "\t" << g_ac[i];
		fd << "\nXT Corrected Multiplicity: Neutron\n";
		for (int i = 0; i < 32; i++)
			fd << "\t" << n_ac[i];
		fd << "\nXT Corrected Multiplicity: Joint\n";
		for (int i = 0; i < 64; i++)
			fd << "\t" << j_ac[i];
	}
	
	
	G4int* n_spec = fRun->DetectedSpectrum();
	
	fd << "\n\n\nSpec\n";
	for (int k = 0; k < 500; k++)
		fd << "\t" << n_spec[k];
	fd << "\n";
	
	fd << "\n\n\nRossi\n";
	for (int i = 0; i < 3; i++)
	{
		G4int* Rossi = fRun->GetRossi(i);
		for (int k = 0; k < 500; k++)
			fd << "\t" << Rossi[k];
		fd << "\n";
	}
	fd << "\n\n\nRossi CX\n";
	for (int i = 0; i < 3; i++)
	{
		G4int* Rossi_cx = fRun->GetRossiCX(i);
		for (int k = 0; k < 500; k++)
			fd << "\t" << Rossi_cx[k];
		fd << "\n";
	}
	
	G4int** tm = fRun->GetTime();
	fd << "\n\n\nTR\n";
	for (int i = 0; i < DETECTOR_COUNT; i++)
	{
		for (int k = 0; k < 500; k++)
			fd << "\t" << tm[i][k];
		fd << "\n";
	}
	
	G4int **hist = fRun->GetAngularPlot();
	G4int **hist3 = fRun->GetAngularContour();

		fd << "\n totals ==>";
		for (int i = 0; i < 16; i++)
			fd << hist[0][i] << "\t";

		fd << "\n singlets ==>";
		for (int i = 0; i < 16; i++)
			fd << hist[1][i] << "\t";


		fd << "\n couplets ==>";
		for (int i = 0; i < 16; i++)
			fd << hist[2][i] << "\t";


		fd << "\n triplets ==>";
		for (int i = 0; i < 16; i++)
			fd << hist[3][i] << "\t";


		fd << "\n quarts ==>";
		for (int i = 0; i < 16; i++)
			fd << hist[4][i] << "\t";

		fd << "\n pentlets ==>";
		for (int i = 0; i < 16; i++)
			fd << hist[5][i] << "\t";


		fd << "\n\n octopus ==>\n";
		for (int i = 0; i < 16; i++) 
		{
			for (int j = 0; j < 16; j++)
				fd << hist3[i][j] << " ";
			fd << "\n";
		}
		

		
		
	hist = fRun->GetAngularCXPlot();
	hist3 = fRun->GetAngularCCXontour();

		fd << "\n totals CX ==>";
		for (int i = 0; i < 16; i++)
			fd << hist[0][i] << "\t";

		fd << "\n singlets CX==>";
		for (int i = 0; i < 16; i++)
			fd << hist[1][i] << "\t";


		fd << "\n couplets CX==>";
		for (int i = 0; i < 16; i++)
			fd << hist[2][i] << "\t";


		fd << "\n triplets CX==>";
		for (int i = 0; i < 16; i++)
			fd << hist[3][i] << "\t";


		fd << "\n quarts CX==>";
		for (int i = 0; i < 16; i++)
			fd << hist[4][i] << "\t";

		fd << "\n pentlets CX==>";
		for (int i = 0; i < 16; i++)
			fd << hist[5][i] << "\t";


		fd << "\n\n octopus CX==>\n";
		for (int i = 0; i < 16; i++) 
		{
			for (int j = 0; j < 16; j++)
				fd << hist3[i][j] << " ";
			fd << "\n";
		}
	fd << "\n\n angu_th ==>\n";
	for (int i = 0; i < 100; i++)
    {
       for (int j = 0; j < 10; j++)
          fd << DetectorConstruction::hist[j][i] << " ";
       fd << "\n";
    }
	fd.close();
	G4cout << "printing theretical multiplicity \n\n" << G4endl;
}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
	G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
	fTimer->Start();
	// save Rndm status
	G4RunManager::GetRunManager()->SetRandomNumberStore(false);
	if (isMaster)
		G4Random::showEngineStatus();

	// keep run condition
	if (fPrimary)
	{
		G4ParticleDefinition* particle =
		fPrimary->GetParticleGun()->GetParticleDefinition();
		G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
		fRun->SetPrimary(particle, energy);
	}

	//histograms
	//
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	if (analysisManager->IsActive())
	{
		analysisManager->OpenFile();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
	if (isMaster)
		fRun->EndOfRun();

	//save histograms
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	if (analysisManager->IsActive())
	{
		analysisManager->Write();
		analysisManager->CloseFile();
	}

	// show Rndm status
	if (isMaster)
		G4Random::showEngineStatus();

	if (IsMaster())
	{
		fTimer->Stop();
		G4cout << "Global result with " << fRun->gEventNumber << " "
				<< *fTimer << G4endl;
	}
	else
	{
		G4cout << "Local thread result with " << fRun->gEventNumber << " "
		<< *fTimer << G4endl;
		return;
	}

	fTimer->Stop();

	G4cout << "printing files \n\n\n" << G4endl;
	multiplicity_th();
	multiplicity_spec(aRun);

	G4cout << "Graceful Exit \n\n\n";
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
