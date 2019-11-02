/*
 * EventAction.h
 *
 *  Created on: 14 Mar 2016
 *      Author: sarwarr
 */

#ifndef INCLUDE_EVENTACTION_HH_
#define INCLUDE_EVENTACTION_HH_

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <fstream>
#include <map>
#include <vector>
#include "G4Threading.hh"
#include "G4AutoLock.hh"

class EventAction : public G4UserEventAction
{
public:
	EventAction(G4String fn= "EventAction" );
	virtual ~EventAction();
// static access method
	//static EventAction* Instance();
	virtual void BeginOfEventAction(const G4Event* event);
	virtual void EndOfEventAction(const G4Event* event);



private:
	//static EventAction* fgInstance;
	int cnt = 0;
	//static G4Mutex aMutex;
};



#endif

