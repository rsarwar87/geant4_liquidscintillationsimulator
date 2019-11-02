
#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4StepPoint.hh"
#include "globals.hh"
#include <map>
#include <fstream>
#include <vector>
#include "RecordedEvent.hh"
class TrackingAction;
#define DETECTOR_COUNT 16

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction : public G4UserSteppingAction
{
  public:
	SteppingAction(G4String, TrackingAction*);
	~SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
	static SteppingAction* Instance();
	// method from the base class
	void Reset();
	
	
	std::vector<RecodedEvent*>* GetRecodedEvents() { return &rEvent; };
	std::vector<RecodedParticle*>* GetRecodedNeutrons() { return &rParticle[1]; };
	std::vector<RecodedParticle*>* GetRecodedPhotons() { return &rParticle[0]; };
	static bool spec;
	static int gwidth;
private:
	
	G4String ParticleName;
	void StackParticle(const G4Step* step, const G4StepPoint * point);
	static G4ThreadLocal RecodedEvent* crEvent;
	static G4ThreadLocal RecodedParticle* crParticle;
	
	
	static G4ThreadLocal G4int _cnnt;
    static G4ThreadLocal G4int _cnnt2;
    static G4ThreadLocal G4double _eng_l, _tm_l;
	
	std::ofstream fout;
	TrackingAction* fTrackingAction;
	G4String filename;
	static G4ThreadLocal SteppingAction* fgInstance;

	std::vector<RecodedEvent*> rEvent;
	std::vector<RecodedParticle*> rParticle[2];

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
