
#include <map>
#include "MyRun.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "SteppingAction.hh"
#include "G4Threading.hh"
#include <RunAction.hh>
#include "G4Timer.hh"

#include "G4Run.hh"
#include "EventAction.hh"
#include <fstream>
#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include <iomanip>
#include "PrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include <iomanip>
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det) :
		G4Run(), fDetector(det), fParticle(0), fEkin(0.), fNbStep1(0), fNbStep2(
				0), fTrackLen1(0.), fTrackLen2(0.), fTime1(0.), fTime2(0.)
{
	Reset();
	filename = "EventAction_" + PrimaryGeneratorAction::name + "_"
			+ std::to_string(PrimaryGeneratorAction::energy);
	fout.open(filename /*+ PrimaryGeneratorAction::name*/,
			std::ios::out | std::ios::trunc);
	fout << "\n\n" << std::setw(8) << "EventId" << "  " << std::setw(14)
			<< "   GlobalTime" << "  " << std::setw(5) << "DetID" << "  "
			<< std::setw(5) << "ParID" << "  " << std::setw(10)
			<< " PrimaryEnrg" << "  " << std::setw(10) << "PrimaryDepo" << "  "
			<< std::setw(6) << "PriCnt" << "  " << "  " << std::setw(10)
			<< "ProtnEnrg" << "  " << std::setw(10) << "ProtonDepo" << "  "
			<< std::setw(6) << "ProCnt" << "  " << std::setw(10) << "ElectEnrg"
			<< "  " << std::setw(10) << "ElectDepo" << "  " << std::setw(6)
			<< "EleCnt" << "  " << std::setw(10) << "OpticEmi" << "  "
			<< std::setw(6) << "OptCnt" << "  " << std::setw(10)
			<< "PMTDetected" << "  " << std::setw(6) << "PMTCnt" << "  "
			<< std::setw(10) << "CerenkovEmm" << "  " << std::setw(10)
			<< "CerenkovCnt \n";

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Mutex Run::aMutex = G4MUTEX_INITIALIZER;
Run::~Run()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{
	fParticle = particle;
	fEkin = energy;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(const G4VProcess* process)
{
	G4String procName = process->GetProcessName();
	std::map<G4String, G4int>::iterator it = fProcCounter.find(procName);
	if (it == fProcCounter.end())
	{
		fProcCounter[procName] = 1;
	}
	else
	{
		fProcCounter[procName]++;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleCount(G4String name, G4double Ekin)
{
	std::map<G4String, ParticleData>::iterator it = fParticleDataMap.find(name);
	if (it == fParticleDataMap.end())
	{
		fParticleDataMap[name] = ParticleData(1, Ekin, Ekin, Ekin);
	}
	else
	{
		ParticleData& data = it->second;
		data.fCount++;
		data.fEmean += Ekin;
		//update min max
		G4double emin = data.fEmin;
		if (Ekin < emin)
			data.fEmin = Ekin;
		G4double emax = data.fEmax;
		if (Ekin > emax)
			data.fEmax = Ekin;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumTrackLength(G4int nstep1, G4int nstep2, G4double trackl1,
		G4double trackl2, G4double time1, G4double time2)
{
	fNbStep1 += nstep1;
	fNbStep2 += nstep2;
	fTrackLen1 += trackl1;
	fTrackLen2 += trackl2;
	fTime1 += time1;
	fTime2 += time2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
	const Run* localRun = static_cast<const Run*>(run);

	//primary particle info
	//
	fParticle = localRun->fParticle;
	fEkin = localRun->fEkin;

	// accumulate sums
	//
	fNbStep1 += localRun->fNbStep1;
	fNbStep2 += localRun->fNbStep2;
	fTrackLen1 += localRun->fTrackLen1;
	fTrackLen2 += localRun->fTrackLen2;
	fTime1 += localRun->fTime1;
	fTime2 += localRun->fTime2;

	//map: processes count
	std::map<G4String, G4int>::const_iterator itp;
	for (itp = localRun->fProcCounter.begin();
			itp != localRun->fProcCounter.end(); ++itp)
	{

		G4String procName = itp->first;
		G4int localCount = itp->second;
		if (fProcCounter.find(procName) == fProcCounter.end())
		{
			fProcCounter[procName] = localCount;
		}
		else
		{
			fProcCounter[procName] += localCount;
		}
	}

	//map: created particles count
	std::map<G4String, ParticleData>::const_iterator itn;
	for (itn = localRun->fParticleDataMap.begin();
			itn != localRun->fParticleDataMap.end(); ++itn)
	{

		G4String name = itn->first;
		const ParticleData& localData = itn->second;
		if (fParticleDataMap.find(name) == fParticleDataMap.end())
		{
			fParticleDataMap[name] = ParticleData(localData.fCount,
					localData.fEmean, localData.fEmin, localData.fEmax);
		}
		else
		{
			ParticleData& data = fParticleDataMap[name];
			data.fCount += localData.fCount;
			data.fEmean += localData.fEmean;
			G4double emin = localData.fEmin;
			if (emin < data.fEmin)
				data.fEmin = emin;
			G4double emax = localData.fEmax;
			if (emax > data.fEmax)
				data.fEmax = emax;
		}
	}

	gEventNumber += localRun->gEventNumber;
	G4cout << "local event count = " << localRun->gEventNumber << " " 
			<< "global event count = " << gEventNumber << G4endl;

	for (uint k = 0; k <  2; k++)
		for (uint i = 0; i < 5000; i++)
		{
			fLightResponse[k][i] += localRun->fLightResponse[k][i];
			fLightHistogram[k][i] += localRun->fLightHistogram[k][i];
			fPMTResponse[k][i] += localRun->fPMTResponse[k][i];
			fPMTHistogram[k][i] += localRun->fPMTHistogram[k][i];
			fIncidentEnergy[k][i] += localRun->fIncidentEnergy[k][i];
			fParticleDeposit[k][i] += localRun->fParticleDeposit[k][i];
			fElectronDeposited[k][i] += localRun->fElectronDeposited[k][i];
			fElectronProduced[k][i] += localRun->fElectronProduced[k][i];
			fProtonProduced[k][i] += localRun->fProtonProduced[k][i];
			fProtonDeposited[k][i] += localRun->fProtonDeposited[k][i];
			fOphotonProduced[k][i] += localRun->fOphotonProduced[k][i];
			fCerenkovProduced[k][i] += localRun->fCerenkovProduced[k][i];
			fOphotonDeposited[k][i] += localRun->fOphotonDeposited[k][i];
			fEntryEnergy[k][i] += localRun->fEntryEnergy[k][i];
			
			
			
		}
	for (int i = 0; i < 3; i++)
		for (int k = 0; k < 500; k++)
		{
			rossi_[i][k] += localRun->rossi_[i][k];
			rossi_cx[i][k] += localRun->rossi_cx[i][k];
		}
	for (int k = 0; k < 500; k++)
		n_spec[k] += localRun->n_spec[k];
	
	for (int i = 0; i < DETECTOR_COUNT; i++)
		for (int k = 0; k < 500; k++)
			_time[i][k] += localRun->_time[i][k];
		
	for (uint i = 0; i < 5000; i++)
	{
		spec_theroy[0][i] += localRun->spec_theroy[0][i];
		spec_theroy[1][i] += localRun->spec_theroy[1][i];
	}
	
	for (int ii = 0; ii < 20; ii++)
	{
		fEventRegistered[ii] += localRun->fEventRegistered[ii];
	}

	for(int i  = 0; i  < 16; i ++)
	for(int ii = 0; ii < 16; ii++)
	{
		angular_plot[i][ii] += localRun->angular_plot[i][ii];
		angular_contour[i][ii] += localRun->angular_contour[i][ii];
		angular_plotcx[i][ii] += localRun->angular_plotcx[i][ii];
		angular_contourcx[i][ii] += localRun->angular_contourcx[i][ii];
	}
	
	for (int ii = 0; ii < 32; ii++)
	{
		multi_theroy[0][ii] += localRun->multi_theroy[0][ii];
		multi_theroy[1][ii] += localRun->multi_theroy[1][ii];
		multi_theroy[2][ii] += localRun->multi_theroy[2][ii];
		multi_detected[0][ii] += localRun->multi_detected[0][ii];
		multi_detected[1][ii] += localRun->multi_detected[1][ii];
		multi_detected[2][ii] += localRun->multi_detected[2][ii];
		multi_detectedcx[0][ii] += localRun->multi_detectedcx[0][ii];
		multi_detectedcx[1][ii] += localRun->multi_detectedcx[1][ii];
		multi_detectedcx[2][ii] += localRun->multi_detectedcx[2][ii];
	}

	
	G4Run::Merge(run);
}



void Run::Reset()
{

	fPMTResponse = new G4int*[ 2];
	fLightResponse = new G4int*[ 2];
	fPMTHistogram = new G4int*[ 2];
	fLightHistogram = new G4int*[ 2];
	fOphotonDeposited = new G4int*[ 2];
	fCerenkovProduced = new G4int*[ 2];

	fEntryEnergy= new G4int*[2];
	fOphotonProduced = new G4int*[ 2];
	fProtonDeposited = new G4int*[ 2];
	fParticleDeposit = new G4int*[ 2];
	fProtonProduced = new G4int*[ 2];
	fElectronProduced = new G4int*[2];
	fElectronDeposited = new G4int*[2];
	fIncidentEnergy = new G4int*[2];

	
	_time = new G4int*[DETECTOR_COUNT];
	for (int i = 0; i < DETECTOR_COUNT; i++)
		_time[i] = (G4int *) calloc(500, sizeof(G4int));
	for (int i = 0; i < 2; i++)
	{
		fPMTResponse[i] = (G4int *) calloc(5000, sizeof(G4int));
		fLightResponse[i] = (G4int *) calloc(5000, sizeof(G4int));
		fPMTHistogram[i] = (G4int *) calloc(5000, sizeof(G4int));
		fLightHistogram[i] = (G4int *) calloc(5000, sizeof(G4int));
		
		fOphotonDeposited[i] = (G4int *) calloc(5000, sizeof(G4int));
		fParticleDeposit[i] = (G4int *) calloc(5000, sizeof(G4int));
		fProtonDeposited[i] = (G4int *) calloc(5000, sizeof(G4int));
		fCerenkovProduced[i] = (G4int *) calloc(5000, sizeof(G4int));
		fOphotonProduced[i] = (G4int *) calloc(5000, sizeof(G4int));
		fProtonProduced[i] = (G4int *) calloc(5000, sizeof(G4int));
		fElectronProduced[i] = (G4int *) calloc(5000, sizeof(G4int));
		fElectronDeposited[i] = (G4int *) calloc(5000, sizeof(G4int));
		fIncidentEnergy[i] = (G4int *) calloc(5000, sizeof(G4int));
		fEntryEnergy[i] = (G4int *) calloc(5000, sizeof(G4int));
		
	}

	rossi_[0] = (G4int *) calloc(500, sizeof(G4int));
	rossi_[1] = (G4int *) calloc(500, sizeof(G4int));
	rossi_[2] = (G4int *) calloc(500, sizeof(G4int));
	rossi_cx[0] = (G4int *) calloc(500, sizeof(G4int));
	rossi_cx[1] = (G4int *) calloc(500, sizeof(G4int));
	rossi_cx[2] = (G4int *) calloc(500, sizeof(G4int));
	n_spec = (G4int *) calloc(500, sizeof(G4int));
	
	fEventRegistered = (G4int *) calloc(64, sizeof(G4int));
	multi_detectedcx[0] = (G4int *) calloc(64, sizeof(G4int));
	multi_detectedcx[1] = (G4int *) calloc(64, sizeof(G4int));
	multi_detectedcx[2] = (G4int *) calloc(64, sizeof(G4int));
	multi_detected[0] = (G4int *) calloc(64, sizeof(G4int));
	multi_detected[1] = (G4int *) calloc(64, sizeof(G4int));
	multi_detected[2] = (G4int *) calloc(64, sizeof(G4int));
	multi_theroy[0] = (G4int *) calloc(64, sizeof(G4int));
	multi_theroy[1] = (G4int *) calloc(64, sizeof(G4int));
	multi_theroy[2] = (G4int *) calloc(64, sizeof(G4int));

	angular_plot = new int*[16];
	angular_contour =  new int*[16];
	angular_plotcx = new int*[16];
	angular_contourcx =  new int*[16];

	for (int i = 0; i < 16; i++)
	{
		angular_contour[i] = (G4int *) calloc(16, sizeof(G4int));
		angular_plot[i] = (G4int *) calloc(16, sizeof(G4int));
		angular_contourcx[i] = (G4int *) calloc(16, sizeof(G4int));
		angular_plotcx[i] = (G4int *) calloc(16, sizeof(G4int));
	}


	spec_theroy[0] = (G4int *) calloc(5000, sizeof(G4int));
	spec_theroy[1] = (G4int *) calloc(5000, sizeof(G4int));
	fEventNumber = 0;
	return;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
	G4int prec = 5, wid = prec + 2;
	G4int dfprec = G4cout.precision(prec);
	//fout.precision(dfprec);

	//run condition
	G4String Particle = fParticle->GetParticleName();
	G4cout << "\n The run is " << numberOfEvent << " " << Particle << " of "
			<< G4BestUnit(fEkin, "Energy") << G4endl;
	//fout << "\n The run is " << numberOfEvent << " " << Particle << " of "
	//		<< G4BestUnit(fEkin, "Energy") << G4endl;

	if (numberOfEvent == 0)
	{
		G4cout.precision(dfprec);
		//fout.precision(dfprec);
		return;
	}

	//frequency of processes
	//
	G4cout << "\n Process calls frequency :" << G4endl;
	//fout << "\n Process calls frequency :" << G4endl;
	G4int survive = 0;
	std::map<G4String, G4int>::iterator it;
	for (it = fProcCounter.begin(); it != fProcCounter.end(); it++)
	{
		G4String procName = it->first;
		G4int count = it->second;
		G4cout << "\t" << procName << " \t = " << count << G4endl;
		//fout << "\t" << procName << " \t = " << count << G4endl;
		if (procName == "Transportation")
			survive = count;
	}
	G4cout << G4endl;
	//fout << G4endl;

	if (survive > 0)
	{
		G4cout << "\n Nb of incident particles surviving: " << survive << G4endl;
		//fout << "\n Nb of incident particles surviving: " << survive << G4endl;
	}

	// total track length of incident neutron 
	//
	G4cout << "\n Parcours of incident neutron:";
	//fout << "\n Parcours of incident neutron:";

	G4double meanCollision1 = (G4double) fNbStep1 / numberOfEvent;
	G4double meanCollision2 = (G4double) fNbStep2 / numberOfEvent ;
	G4double meanCollisTota = meanCollision1 + meanCollision2;

	G4cout << "\n   nb of collisions    E>1*eV= " << meanCollision1
			<< "      E<1*eV= " << meanCollision2 << "       total= "
			<< meanCollisTota;
	/*fout << "\n   nb of collisions    E>1*eV= " << meanCollision1
			<< "      E<1*eV= " << meanCollision2 << "       total= "
			<< meanCollisTota;*/

	G4double meanTrackLen1 = fTrackLen1 / numberOfEvent;
	G4double meanTrackLen2 = fTrackLen2 / numberOfEvent;
	G4double meanTrackLtot = meanTrackLen1 + meanTrackLen2;

	G4cout << "\n   track length        E>1*eV= "
			<< G4BestUnit(meanTrackLen1, "Length") << "  E<1*eV= "
			<< G4BestUnit(meanTrackLen2, "Length") << "   total= "
			<< G4BestUnit(meanTrackLtot, "Length");
	/*fout << "\n   track length        E>1*eV= "
			<< G4BestUnit(meanTrackLen1, "Length") << "  E<1*eV= "
			<< G4BestUnit(meanTrackLen2, "Length") << "   total= "
			<< G4BestUnit(meanTrackLtot, "Length");*/

	G4double meanTime1 = fTime1 / numberOfEvent;
	G4double meanTime2 = fTime2 / numberOfEvent;
	G4double meanTimeTo = meanTime1 + meanTime2;

	G4cout << "\n   time of flight      E>1*eV= "
			<< G4BestUnit(meanTime1, "Time") << "  E<1*eV= "
			<< G4BestUnit(meanTime2, "Time") << "   total= "
			<< G4BestUnit(meanTimeTo, "Time") << G4endl;
	/*fout << "\n   time of flight      E>1*eV= " << G4BestUnit(meanTime1, "Time")
			<< "  E<1*eV= " << G4BestUnit(meanTime2, "Time") << "   total= "
			<< G4BestUnit(meanTimeTo, "Time") << G4endl;*/

	//particles count
	//
	G4cout << "\n List of generated particles:" << G4endl;
	//fout << "\n List of generated particles:" << G4endl;

	std::map<G4String, ParticleData>::iterator itn;
	for (itn = fParticleDataMap.begin(); itn != fParticleDataMap.end(); itn++)
	{
		G4String name = itn->first;
		ParticleData data = itn->second;
		G4int count = data.fCount;
		G4double eMean = data.fEmean / count;
		G4double eMin = data.fEmin;
		G4double eMax = data.fEmax;

		G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
				<< "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
				<< "\t( " << G4BestUnit(eMin, "Energy") << " --> "
				<< G4BestUnit(eMax, "Energy") << ")" << G4endl;
		/*fout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
				<< "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
				<< "\t( " << G4BestUnit(eMin, "Energy") << " --> "
				<< G4BestUnit(eMax, "Energy") << ")" << G4endl;*/
	}

	//normalize histograms
	////G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	////G4double factor = 1./numberOfEvent;
	////analysisManager->ScaleH1(3,factor);

	//remove all contents in fProcCounter, fCount
	fProcCounter.clear();
	fParticleDataMap.clear();

	//restore default format
	G4cout.precision(dfprec);
	//fout.precision(dfprec);
}
void Run::RecordEvent(const G4Event* evt)
{

	gEventNumber++;
	
	if (gEventNumber % 1000 == 0)
		G4cout << "NPS: " << gEventNumber << " Neutron: " << multi_detected[1][1] 
					<< " " << multi_detected[1][2] << " " << multi_detected[1][3] << G4endl;
	
	SteppingAction* SA = SteppingAction::Instance();
	std::vector<RecodedParticle*>* nContainer = SA->GetRecodedNeutrons();
	std::vector<RecodedParticle*>* pContainer = SA->GetRecodedPhotons();
	/*G4cout << "====================" << G4endl;
	G4cout << "Event = " <<  gEventNumber << G4endl;
	G4cout << "number of neutron = " <<  nContainer->size() << G4endl;
	G4cout << "number of  photon = " <<  pContainer->size() << G4endl;*/
	int cnt_g = 0, cnt_n = 0;
	for (int i = 0; i < nContainer->size(); i++) 
	{
		//nContainer->at(i)->Pout(i);
		if ((nContainer->at(i)->IncidentEnergy/keV) < 9555)
			spec_theroy[1][(int)(nContainer->at(i)->IncidentEnergy/keV/2)]++;
		cnt_n++;
		
	}
	for (int i = 0; i < pContainer->size(); i++) 
	{
		//pContainer->at(i)->Pout(i);
		if ((pContainer->at(i)->IncidentEnergy/keV) < 9555)
			spec_theroy[0][(int)(pContainer->at(i)->IncidentEnergy/keV/2)]++;
		cnt_g++;
	}
	if (cnt_g < 33) multi_theroy[0][cnt_g]++;
	if (cnt_n < 33) multi_theroy[1][cnt_n]++;
	if (cnt_g+cnt_n < 33) multi_theroy[2][cnt_g+cnt_n]++;
			
	std::vector<RecodedEvent*>* eContainer = SA->GetRecodedEvents();
	
	std::vector<RecodedEvent*> dPhotons;
	std::vector<RecodedEvent*> dNeutrons;
	std::vector<RecodedEvent*> dJoint;
	int idpn = 0;
	if (eContainer->size()>0){
		
		for (int i = 0; i < eContainer->size(); i++) 
		{
			RecodedEvent* tmp = eContainer->at(i);
			int detID = 0;
			//tmp->Print();
			idpn = tmp->ParticleType();
			if ((tmp->depoEnergy[idpn][PMT] == 0) && (tmp->reacEnergy[idpn][PMT] > 0)) continue;
			int fPartType = (tmp->name == "neutron" ? 1 : 0);
			if ((tmp->depoEnergy[idpn][PMT] == 0) && (tmp->depoEnergy[idpn][PMT] > tmp->reacEnergy[idpn][fPartType])) continue;
			if (!tmp->IsValid()) continue;
			
			
			if (fPartType == 1) dNeutrons.push_back(tmp);
			else dPhotons.push_back(tmp);
			dJoint.push_back(tmp);
			
			if ((tmp->particledef.at(0).IncidentEnergy/keV) < 9555)
				fIncidentEnergy[fPartType][(int) (tmp->ptr_particledef->IncidentEnergy/keV/2)]++;
			if ((tmp->particledef.at(0).EntryEnergy/keV) < 9555)
				fEntryEnergy[fPartType][(int) (tmp->ptr_particledef->EntryEnergy/keV/2)]++;
			
			
				if (tmp->reacEnergy[idpn][fPartType]/keV < 9555)
					fParticleDeposit[fPartType ][(int) (tmp->reacEnergy[idpn][fPartType]/keV/2)]++;

			
			if (tmp->reacEnergy[idpn][2]/keV < 9555)
				fElectronDeposited[fPartType][(int) (tmp->reacEnergy[idpn][2]/keV/2)]++;
			
			if (tmp->depoEnergy[idpn][electron]/keV < 9555)
				fElectronProduced[fPartType][(int) (tmp->depoEnergy[idpn][electron]/keV/2)]++;
			
			if (tmp->depoEnergy[idpn][proton]/keV > 0 && tmp->depoEnergy[idpn][proton]/keV < 9555)
				fProtonProduced[fPartType][(int) (tmp->depoEnergy[idpn][proton]/keV/2)]++;
			
			if (tmp->reacEnergy[idpn][hIoni]/keV > 0 && tmp->reacEnergy[idpn][hIoni]/keV < 9555)
				fProtonDeposited[fPartType][(int) (tmp->reacEnergy[idpn][hIoni]/keV/2)]++;
			
			if (tmp->depoEnergy[idpn][optical]/keV > 0 && tmp->depoEnergy[idpn][optical]/keV < 9555)
				fOphotonProduced[fPartType][(int) (tmp->depoEnergy[idpn][optical]/keV/2)]++;
			if (tmp->depoEnergy[idpn][Cerenkov]/keV > 0 && tmp->depoEnergy[idpn][Cerenkov]/keV < 9555)
				fCerenkovProduced[fPartType][(int) (tmp->depoEnergy[idpn][Cerenkov]/keV/2)]++;
			
			if (tmp->depoEnergy[idpn][PMT]/keV > 0 && tmp->depoEnergy[idpn][PMT]/keV  < 9555)
				fOphotonDeposited[fPartType][(int) (tmp->depoEnergy[idpn][PMT]/keV /2)]++;
			
			
			fEventRegistered[detID]++;
			fEventNumber++;
		}
	}
	
	
	std::map<double, entry> stNeutron;
	std::map<double, entry> stJoint;
	for (int i = 0; i < dNeutrons.size(); i++)
	{
		double t = dNeutrons.at(i)->GetTriggerTime(RunAction::cutOu)/ns;
		//G4cout << i <<" y " <<  t <<  " " << dNeutrons.at(i)->depoEnergy[PMT]/keV <<G4endl;
		if (t < 0) continue;
		entry tmp_;
		tmp_.didx = dNeutrons.at(i)->detectorID;
		tmp_.pid = dNeutrons.at(i)->particledef.at(0).particleid;
		stNeutron[t] = tmp_;
		stJoint[t] = tmp_;
		
		if (t < 500)
			if (i > 0 || PrimaryGeneratorAction::mono )
			n_spec[(int)t]++;
		
	}
	
	//if (dNeutrons.size() > 0) G4cout << dNeutrons.size() << " " << stNeutron.size() << G4endl;
	std::map<double, entry> stPhoton;
	for (int i = 0; i < dPhotons.size(); i++)
	{
		double t = dPhotons.at(i)->GetTriggerTime(RunAction::cutOu)/ns;
		if (t < 0) continue;
		entry tmp_;
		tmp_.didx = dPhotons.at(i)->detectorID;
		tmp_.pid = dPhotons.at(i)->particledef.at(0).particleid;
		stPhoton[t] = tmp_;
		stJoint[t] = tmp_;
	}
	
	int multi_n = 0, multi_g = 0, multi_j = 0;
	int multi_nxc = 0, multi_gxc = 0, multi_jxc = 0;
	std::vector<int> angular, ang, an;
	std::vector<int> angularcx;
	
	if ((stNeutron.size() > 0))
		if (stNeutron.size() == 1) {multi_n++;multi_nxc++;}
		else
			ProcessCoincidence(stNeutron, angular, angularcx, multi_n, multi_nxc, rossi_[1], rossi_cx[1]);
	
	if ((stPhoton.size() > 0))
		if (stPhoton.size() == 1) {multi_g++;multi_gxc++;}
		else
			ProcessCoincidence(stPhoton, ang, an, multi_g, multi_gxc, rossi_[0], rossi_cx[0]);
	an.clear(); ang.clear();
	if ((stJoint.size() > 0) ) 
		if (stJoint.size() == 1) {multi_j++;multi_jxc++;}
		else
			ProcessCoincidence(stJoint, ang, an, multi_j, multi_jxc, rossi_[2], rossi_cx[2]);
		
	multi_detected[1][multi_n]++;
	multi_detected[0][multi_g]++;
	multi_detected[2][multi_j]++;
	multi_detectedcx[1][multi_nxc]++;
	multi_detectedcx[0][multi_gxc]++;
	multi_detectedcx[2][multi_jxc]++;
		
		

	if (angular.size() > 1)
	{
		int base = angular.at(0)+1;
		int shift = 8 - base;
		int v = 0;
		for (int i = 1; i < angular.size(); i++)
		{
			int val = angular.at(i) + shift + 1;
			if (val < 1 ) val += 15;
			else if (val > 15) val -= 15;
				
			angular_plot[0][val]++;
			angular_plot[i][val]++;
			if (i == 1) v = val;
			else if (i == 2) angular_contour[v][val]++;
		}
	}
	if (angularcx.size() > 1)
	{
		int base = angularcx.at(0)+1;
		int shift = 8 - base;
		int v = 0;
		for (int i = 1; i < angularcx.size(); i++)
		{
			int val = angularcx.at(i) + shift + 1;
			if (val < 1 ) val += 15;
			else if (val > 15) val -= 15;
				
			angular_plotcx[0][val]++;
			angular_plotcx[i][val]++;
			if (i == 1) v = val;
			else if (i == 2) angular_contourcx[v][val]++;
		}
	}
	
	
	SteppingAction::Instance()->Reset();
	//G4cout << "\n";

	G4Run::RecordEvent(evt);
}
void Run::ProcessCoincidence(std::map<double, entry> storage, std::vector<int>& angular, 
				std::vector<int>& angularcx, int &multi, int &multi_cx, int* rossi, int* rossicx)
{
	std::map<double, entry>::iterator it = storage.begin();
	double tim = -1;
	bool cx_map[64] = {false};
	for (; it != storage.end(); it++)
	{				
		entry en = it->second;
		//G4cout << "ProcessCoincidence tim " << tim << G4endl; 
		int deltaT = (int)(it->first - tim); 
		if (tim == -1) 
		{
			tim = it->first; 
			multi++; 
			multi_cx++; 
			angular.push_back(en.didx);
			angularcx.push_back(en.didx);
			cx_map[en.pid] = true; 
		}
		else if (deltaT  < SteppingAction::gwidth) 
		{
			//G4cout << "ProcessCoincidence multi t " << en.didx << G4endl;
			multi++;
			angular.push_back(en.didx);
			rossi[deltaT]++;
			if (!cx_map[en.pid]) {
				multi_cx++;
				//G4cout << "ProcessCoincidence cx multi_cx " << multi_cx << G4endl;
				angularcx.push_back(en.didx);
				rossicx[deltaT]++;
				//G4cout << "ProcessCoincidence cx angular" << angularcx.size() << G4endl;
				//G4cout << "ProcessCoincidence cs rossi" << (int) (it->first - tim) << G4endl;
			}
			else
				_time[en.didx][deltaT]++;
			cx_map[en.pid] = true; 
			//G4cout << "ProcessCoincidence angular" << angular.size() << G4endl;
			//G4cout << "ProcessCoincidence rossi" << (int) (it->first - tim) << G4endl;
		} 
		else if(deltaT < 500)
		{
			//G4cout << "ProcessCoincidence rossi t " << x << G4endl;
			rossi[deltaT]++;
			if (!cx_map[en.pid]) {
				rossicx[deltaT]++;
			}
			else
				_time[en.didx][deltaT]++;
			cx_map[en.pid] = true; 
		}
	}
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
