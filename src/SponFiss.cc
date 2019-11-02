//******************************************************************************
// SponFiss.cc
//
// 1.00 JMV, LLNL, OCT-2010: version compatible with Geant 4.9.3.
//******************************************************************************
//
#include "SponFiss.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorConstruction.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include <cmath>
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "SingleSource.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGeneratorAction.hh"

//----------------------------------------------------------------------------//
SponFiss::SponFiss()
{
}

//----------------------------------------------------------------------------//
SponFiss::SponFiss(G4int iso,
		G4SPSPosDistribution *pos/*, PrimaryGeneratorAction *run*/)
{
	neutron_definition = G4Neutron::Neutron();
	photon_definition = G4Gamma::Gamma();

	// verbosity
	verbosityLevel = 0;
	//fRun = run;
	isotope = iso;
	posDist = pos;
}

//----------------------------------------------------------------------------//
SponFiss::~SponFiss()
{
}

//----------------------------------------------------------------------------//
void SponFiss::GeneratePrimaryVertex(G4Event* anEvent, G4double time, int mode)
{
	// Generate a spontaneous fission using the fission library and emit
	// the neutrons and gamma-rays

	fissionEvent* fe = new fissionEvent(isotope, 0, -1., 0., 0);
	fe->setCf252Option(2, 0);
	fe->setCorrelationOption(mode);
	if (3 == fe->getCorrelationOption())
	{
		int err_len = 1000;
		char* error_message = new char[err_len];
		fe->getFREYAerrors(&err_len, error_message);
		if (err_len>1)
		{
			G4ExceptionDescription ed;
			ed << "Call to new fissionEvent("
			<< "isotope=" << isotope << ", "
			<< "time=" << time << ", "
			<< "nubar=-1." << ", "
			<< "eng=0." << ", "
			<< "0) failed with error message from FREYA: "
			<< G4endl
			<< error_message;
			delete [] error_message;
			G4Exception("G4FissionLibrary_new::SampleMult", "freya001", FatalException,
					ed);
		}
		delete [] error_message;
	}
	G4int nPrompt, gPrompt;
	nPrompt = fe->getNeutronNu();
	gPrompt = fe->getPhotonNu();
	//fRun->fRun->UpdateSource(nPrompt, gPrompt);
	
	//nPrompt = nPrompt != 0 ? 1 : 0;


	if (verbosityLevel > 1)
	{
		G4cout << " nPrompt: " << nPrompt  << "   gPrompt: " << gPrompt
		<< G4endl;
	}

	// Position
	//posDist = GetPosDist();
	G4ThreeVector sampled_particle_position = DetectorConstruction::source_position;

	// create a new vertex
	G4PrimaryVertex* vertex = new G4PrimaryVertex(sampled_particle_position,
			0.);

	G4double mom, momx, momy, momz, eng;

	if (verbosityLevel >= 2)
		G4cout << "Creating primaries and assigning to vertex" << G4endl;

	for(int n1=0; n1<nPrompt; n1++) {
         double u1 = fe->getNeutronDircosu(n1), v1 = fe->getNeutronDircosv(n1), w1 = fe->getNeutronDircosw(n1);
         for(int n2=n1+1; n2<nPrompt; n2++) 
		 {
            double u2 = fe->getNeutronDircosu(n2), v2 = fe->getNeutronDircosv(n2), w2 = fe->getNeutronDircosw(n2);
            double scalar_product = u1*u2+v1*v2+w1*w2;
        
            int bin_index = (int) (98*(scalar_product+1)/2);
            DetectorConstruction::hist[0][bin_index]++;
	        if (n2 - n1 < 10)
              DetectorConstruction::hist[n2-n1][bin_index]++;
        }
    }
	
	G4DynamicParticle* it;
	// Build neutrons
	if (PrimaryGeneratorAction::neutron && nPrompt > 0)
	for (G4int i = nPrompt-1; i > -1; i--)
	{
		it = new G4DynamicParticle();
		it->SetDefinition(neutron_definition);
			eng = fe->getNeutronEnergy(i);
		if (eng > 19.9)
			eng = 19;
		it->SetKineticEnergy(eng);
		mom = it->GetTotalMomentum();

		momx = mom * fe->getNeutronDircosu(i);
		momy = mom * fe->getNeutronDircosv(i);
		momz = mom * fe->getNeutronDircosw(i);

		G4PrimaryParticle* particle = new G4PrimaryParticle(neutron_definition,
				momx, momy, momz, eng * MeV);
		//particle->SetMomentumDirection(G4ThreeVector(1., 0., 0.));
		particle->SetMass(neutron_definition->GetPDGMass());
		particle->SetCharge(neutron_definition->GetPDGCharge());
		particle->SetPolarization(particle_polarization.x(),
				particle_polarization.y(), particle_polarization.z());


		if (verbosityLevel > 1)
		{
			//G4cout << "Particle name: "
			//		<< particle->GetG4code()->GetParticleName() << G4endl;
			//G4cout << "     Momentum: " << particle->GetMomentum() << G4endl;
			//G4cout << "     Position: " << vertex->GetPosition() << G4endl;
			G4cout << "     Energy:  " << particle->GetKineticEnergy() << G4endl;
		}

			if (fe->getNeutronAge(i) != -1)
				particle->SetProperTime(fe->getNeutronAge(i) * ns);
			else
				particle->SetProperTime(0 * ns);

		vertex->SetPrimary(particle);
			/*	fRun - fRun->UpdateSource(1, i, eng);
			double u1 = fe->getNeutronDircosu(i), v1 = fe->getNeutronDircosv(i),
					w1 = fe->getNeutronDircosw(i);
			for (int n2 = i + 1; n2 < nPrompt; n2++)
			{
				double u2 = fe->getNeutronDircosu(n2), v2 =
						fe->getNeutronDircosv(n2), w2 = fe->getNeutronDircosw(
						n2);
				double scalar_product = u1 * u2 + v1 * v2 + w1 * w2;

				int andex = (int) (100 * (scalar_product + 1) / 2);
				fRun - fRun->UpdateSource(1, i, andex, 0);

			}
			 */
			delete it;
	}


	// Build gammas
	if (PrimaryGeneratorAction::gamma)
	for (G4int i = 0; i < gPrompt; i++)
	{
		it = new G4DynamicParticle();
		it->SetDefinition(photon_definition);
		eng = fe->getPhotonEnergy(i);
			if (eng > 19.9)
				eng = 19;
		it->SetKineticEnergy(eng);
		mom = it->GetTotalMomentum();

		momx = mom * fe->getPhotonDircosu(i);
		momy = mom * fe->getPhotonDircosv(i);
		momz = mom * fe->getPhotonDircosw(i);

		G4PrimaryParticle* particle = new G4PrimaryParticle(photon_definition,
				momx, momy, momz, eng * MeV);
		particle->SetMass(photon_definition->GetPDGMass());
		particle->SetCharge(photon_definition->GetPDGCharge());
		particle->SetPolarization(particle_polarization.x(),
				particle_polarization.y(), particle_polarization.z());

		if (verbosityLevel > 1)
		{
			G4cout << "Particle name: "
					<< particle->GetG4code()->GetParticleName() << G4endl;
			G4cout << "     Momentum: " << particle->GetMomentum() << G4endl;
			G4cout << "     Position: " << vertex->GetPosition() << G4endl;
		}

			if (fe->getPhotonAge(i) != -1)
				particle->SetProperTime(fe->getPhotonAge(i) * ns);
			else
				particle->SetProperTime(0 * ns);
			vertex->SetPrimary(particle);

			/*fRun - fRun->UpdateSource(1, i, eng);
			double u1 = fe->getPhotonDircosu(i), v1 = fe->getPhotonDircosv(i),
					w1 = fe->getPhotonDircosw(i);
			for (int n2 = i + 1; n2 < gPrompt; n2++)
			{
				double u2 = fe->getNeutronDircosu(n2), v2 =
						fe->getNeutronDircosv(n2), w2 = fe->getNeutronDircosw(
						n2);
				double scalar_product = u1 * u2 + v1 * v2 + w1 * w2;

				int andex = (int) (100 * (scalar_product + 1) / 2);
				fRun - fRun->UpdateSource(1, i, andex, 0);

			 }*/

		delete it;
	}
	delete fe;
	vertex->SetT0(time);
//  G4cout << "         Time: "<<vertex->GetT0()/second << G4endl;

	anEvent->AddPrimaryVertex(vertex);
	if (verbosityLevel > 1)
		G4cout << " Primary Vetex generated !" << G4endl;
	}
