
#define FISSION_NEW
#define USEFREYA
#define G4MULTITHREADED
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#include <PhysicsList.hh>
#include <DetectorConstruction.hh>

#include <ActionInitialization.hh>

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
#include "SourceListing.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " OpNovice [-m macro ] [-u UIsession] [-t nThreads] [-r seed] [-c cutoff] [-mono beam]"
           << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 20 ) {
    PrintUsage();
    return 1;
  }
  SourceListing *SL;
  //SL = new SourceListing();
  G4String macro;
  G4String session;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif

  G4long myseed = 345354;
  for ( G4int i=1; i<argc; i=i+2 ) 
  {
		if      ( G4String(argv[i]) == "-m" ) macro   = argv[i+1];
		else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
		else if ( G4String(argv[i]) == "-r" ) myseed  = atoi(argv[i+1]);
		else if (G4String(argv[i]) == "-c")
			RunAction::cutOu = atoi(argv[i + 1]);
		else if (G4String(argv[i]) == "-l")
		{
			DetectorConstruction::nb_cryst = 15;
			DetectorConstruction::ring_R1 = 26.25 * cm;
		}
		else if (G4String(argv[i]) == "-w")
		{
			DetectorConstruction::ring_W1 = atoi(argv[i+1])/10 * cm;
		}
		else if (G4String(argv[i]) == "-j")
		{
			PrimaryGeneratorAction::gamma = true;
			PrimaryGeneratorAction::neutron = true;
		}
		else if (G4String(argv[i]) == "-mono") 
		{
			PrimaryGeneratorAction::mono = true;
			if ( atoi(argv[i + 1]) == 1 ) 
				PrimaryGeneratorAction::beam = true;

		}
        else if (G4String(argv[i]) == "-AmLi")
        {
            PrimaryGeneratorAction::AmLi = true;
            PrimaryGeneratorAction::neutron = true;
			PrimaryGeneratorAction::energy = 4121;
        }
        else if (G4String(argv[i]) == "-cmod")
        {
            PrimaryGeneratorAction::sfif = true;
			std::string s1(argv[i+1]);
			SponFiss_FF::filename = s1;
        }
        else if (G4String(argv[i]) == "-Co")
        {
            PrimaryGeneratorAction::Co = true;
            PrimaryGeneratorAction::neutron = false;
            PrimaryGeneratorAction::gamma = true;
			PrimaryGeneratorAction::name = "gamma";
			PrimaryGeneratorAction::mono = false;
			PrimaryGeneratorAction::beam = false;
			PrimaryGeneratorAction::energy = 1121;
        }
		else if (G4String(argv[i]) == "-g")
		{
			PrimaryGeneratorAction::gamma = true;
			PrimaryGeneratorAction::neutron = false;
			PrimaryGeneratorAction::energy = atoi(argv[i + 1]);
            PrimaryGeneratorAction::name = "gamma";

		}
		else if (G4String(argv[i]) == "-n")
		{
			PrimaryGeneratorAction::gamma = false;
			PrimaryGeneratorAction::neutron = true;
			PrimaryGeneratorAction::energy = atoi(argv[i + 1]);
		}
		else if (G4String(argv[i]) == "-mode")
		{
			PrimaryGeneratorAction::mode = atoi(argv[i + 1]);
		}
		else if (G4String(argv[i]) == "-gw")
		{
			SteppingAction::gwidth = atoi(argv[i + 1]);
		}
		else if (G4String(argv[i]) == "-lancs")
		{
			DetectorConstruction::lancs = true;
		}

#ifdef G4MULTITHREADED
		else if ( G4String(argv[i]) == "-t" ) {
                    nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
        }
#endif
		else {
			std::cout << argv[i] << std::endl;
			PrintUsage();
			return 1;
		}
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  // Seed the random number generator manually
  G4Random::setTheSeed(myseed);
	DetectorConstruction* dDet = new DetectorConstruction();
  // Set mandatory initialization classes
  //
  // Detector construction
	runManager->SetUserInitialization(dDet);
  // Physics list
  runManager-> SetUserInitialization(new PhysicsList());
  // User action initialization
	runManager->SetUserInitialization(new ActionInitialization(dDet));

  // Initialize G4 kernel
  //
  runManager->Initialize();

#ifdef G4VIS_USE
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro.size() ) {
     // Batch mode
     G4String command = "/control/execute ";
     UImanager->ApplyCommand(command+macro);
  }
  else // Define UI session for interactive mode
  {
#ifdef G4UI_USE
     G4UIExecutive * ui = new G4UIExecutive(argc,argv,session);
#ifdef G4VIS_USE
     UImanager->ApplyCommand("/control/execute vis.mac");
#else
     UImanager->ApplyCommand("/control/execute OpNovice.in");
#endif
     if (ui->IsGUI())
        UImanager->ApplyCommand("/control/execute gui.mac");
     ui->SessionStart();
     delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
