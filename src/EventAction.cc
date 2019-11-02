
#include "EventAction.hh"
#include <DetectorConstruction.hh>
#include <RunAction.hh>
#include <SteppingAction.hh>
#include <iomanip>
#include <vector>
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

// EventAction* EventAction::fgInstance = 0;
// int EventAction::cnt = 0;
// G4Mutex EventAction::aMutex = G4MUTEX_INITIALIZER;

EventAction::EventAction(G4String fn) : G4UserEventAction() {
  // fgInstance = this;
}

EventAction::~EventAction() {
  // fgInstance = 0;
}
// static access method
void EventAction::BeginOfEventAction(const G4Event* event) {
  // G4MUTEXLOCK(&aMutex);
  // G4cout << "EventAction::BeginOfEventAction" << G4endl;
  // G4MUTEXUNLOCK(&aMutex);
  return;
}

void EventAction::EndOfEventAction(const G4Event* event) {
  // G4cout << "EventAction::EndOfEventAction" << G4endl;

  return;
}
