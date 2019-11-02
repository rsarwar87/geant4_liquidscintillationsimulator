
#ifndef RecodedEvent_h
#define RecodedEvent_h 1

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
enum ptype {
  PMT = 0,
  optical = 1,
  proton = 2,
  neutron = 3,
  electron = 4,
  Cerenkov = 5,
  holder = 6
};
enum rtype { compt = 0, hadElastic = 1, eIoni = 2, hIoni = 3, alpha = 4 };
class RecodedParticle {
 public:
  RecodedParticle(int pid, std::string nm, G4double eng, G4double tm,
                  G4ThreeVector dir) {
    particleid = pid;
    name = nm;
    direction = dir;
    IncidentEnergy = eng;
    Time = tm;
  }
  std::string Pout(int j) {
    G4cout << "Particle No " << particleid << G4endl;
    G4cout << "    Particle Type " << name << G4endl;
    G4cout << "    IncidentEnergy " << IncidentEnergy / MeV << G4endl;
    G4cout << "    EntryEnergy " << EntryEnergy / MeV << G4endl;
    G4cout << "    Time " << Time << G4endl;
    G4cout << "    Direction " << direction << G4endl;
  }

  int particleid = -1;
  std::string name = "";
  G4ThreeVector direction;
  G4double IncidentEnergy = 0;
  G4double EntryEnergy = 0;
  G4double Time = 0;
  bool del = true;
};

class RecodedEvent {
 public:
  RecodedEvent(int det, RecodedParticle *nue, G4double time) {
    detectorID = det;
    particledef.push_back(*nue);
    masterTime = time;
    triggermap = new std::map<G4double, G4double>[100];
    for (int i = 0; i < 7; i++) {
      intervaltimeEnergy[i] = (G4double *)calloc(10000, sizeof(G4double));
      intervaltimeCounter[i] = (G4int *)calloc(10000, sizeof(G4int));
    }

    for (int i = 0; i < 100; i++)
      for (int j = 0; j < 7; j++) {
        firstDepo[i][j] = -1;
        firstInteraction[i][j] = -1;
      }
    for (int i = 0; i < 100; i++) triggermap[i].clear();
  };
  void delete_class() {
    for (int i = 0; i < 100; i++) triggermap[i].clear();
    delete[] triggermap;
    for (int i = 0; i < 7; i++) {
      if (intervaltimeEnergy[i] != NULL) free(intervaltimeEnergy[i]);
      if (intervaltimeCounter[i] != NULL) free(intervaltimeCounter[i]);
    }
  }
  std::map<G4double, G4double> *triggermap;
  int detectorID = -1;
  std::vector<RecodedParticle> particledef;
  RecodedParticle *ptr_particledef = NULL;
  G4double masterTime = -1;
  G4double triggerTime = -1;
  std::string name = "";
  G4double firstDepo[100][7] = {{-1, -1, -1, -1, -1, -1, -1}};
  G4double firstInteraction[100][7] = {{-1, -1, -1, -1, -1, -1, -1}};
  int nNeutron = 0;
  int nPhoton = 0;
  int nUnique = 0;
  G4double reacCounter[100][5] = {{0}};
  G4double reacEnergy[100][5] = {{0}};
  G4double depoCounter[100][7] = {{0}};
  G4double depoEnergy[100][7] = {{0}};
  G4double *intervaltimeEnergy[7] = {NULL};
  G4int *intervaltimeCounter[7] = {NULL};
  G4int idp = 0, idn = 0, idx = 0;
  G4int pCount = 0;
  bool IsValid() {
    ParticleType();
    // return true;
    if (name == "neutron" &&
        reacEnergy[idx][1] / keV < 1 * depoEnergy[idx][0] / keV)
      return false;  // 4 1.5

    if (name == "neutron" &&
        ptr_particledef->EntryEnergy / keV < reacEnergy[idx][1] / keV)
      return false;

    // if (name == "neutron" && reacEnergy[idx][1]/keV < 3.5*(cutoff)) return
    // false; //4 1.5

    return true;
  }

  int ParticleType() {
    /*ptr_particledef = &particledef.at(0);
    name = ptr_particledef->name;
    return 0;*/
    G4double engp = 0, engn = 0;
    if (ptr_particledef != NULL) return idx;
    for (int i = 0; i < particledef.size(); i++) {
      if (reacEnergy[i][0] > engp) {
        engp = reacEnergy[i][0];
        idp = i;
      }
      if (reacEnergy[i][1] > engn) {
        engn = reacEnergy[i][1];
        idn = i;
      }
    }
    if (engp > engn * 2) {
      ptr_particledef = &particledef.at(idp);
      idx = idp;
    } else {
      ptr_particledef = &particledef.at(idn);
      idx = idn;
    }

    name = ptr_particledef->name;
    // if (name != "neutron") G4cout << ptr_particledef->name << G4endl;
    /*if (!IsValid())
    {
            Print();
            name = "gamma";
            ptr_particledef->name = name;
    }*/
    return idx;
  }

  G4double GetTriggerTime(G4double cutoff) {
    if (!IsValid()) return -1;
    triggerTime = -1;
    if (depoEnergy[idx][0] / keV < cutoff) return triggerTime;
    // if (name == "gamma" && reacEnergy[idx][0]/keV < (cutoff)) return
    // triggerTime;

    if (triggermap[idx].size() < 1) return triggerTime;
    std::map<G4double, G4double>::iterator it = triggermap[idx].begin();
    for (; it != triggermap[idx].end(); it++)
      if (it->second / keV > cutoff) {
        triggerTime = it->first;
        break;
      }
    // G4cout << "GetTriggerTime cutoff depoEnergy " << cutoff << " " <<
    // depoEnergy[idx][0]/keV << " " << triggerTime  << G4endl;
    return triggerTime;
  };

  bool hasTriggered(G4double ftrigger, G4double gwidth, G4double cutoff = 0) {
    triggerTime = GetTriggerTime(cutoff);
    if (triggerTime < 0) return false;
    return ((ftrigger + gwidth) < triggerTime);
  }

  bool CheckPart(RecodedParticle *crParticle) {
    // return;
    if (particledef.size() > 95) return false;
    for (int i = 0; i < particledef.size(); i++) {
      pCount = i;
      if (particledef.at(i).particleid == crParticle->particleid) return true;
    }
    particledef.push_back(*crParticle);
    pCount = particledef.size() - 1;
    return false;
  }

  void recordProduction(int _idx, G4double eng, G4double time = 0) {
    if (_idx < 2) {
      int id = 0;
      if (firstDepo[pCount][_idx] == -1)
        firstDepo[pCount][_idx] = time;
      else
        id = (int)(time * 10 - firstDepo[pCount][_idx] * 10);
      if (id > -1 && id < 10000) {
        intervaltimeCounter[_idx][id]++;
        intervaltimeEnergy[_idx][id] += eng;
      }
    }

    depoEnergy[pCount][_idx] += eng;
    depoCounter[pCount][_idx]++;

    if (_idx == 0) triggermap[pCount][time] = depoEnergy[pCount][_idx];
  }
  void recordReaction(int _idx, G4double eng, G4double time = 0,
                      G4double eng2 = 0) {
    reacEnergy[pCount][_idx] += eng;
    reacCounter[pCount][_idx]++;

    if (_idx < 2) {
      if (firstInteraction[pCount][_idx] == -1) {
        firstInteraction[pCount][_idx] = time;
        particledef.at(pCount).EntryEnergy = eng2;
      }
    }
  }

  void Print() {
    ParticleType();
    if (depoEnergy[idx][PMT] == 0) return;

    G4cout << " \\\\\\\\\\\\\\ " << G4endl;
    G4cout << "Det No " << detectorID << G4endl;
    G4cout << "    masterTime " << G4BestUnit(masterTime, "Time") << G4endl;
    G4cout << "Number of Particle" << particledef.size() << G4endl;
    for (int i = 0; i < particledef.size(); i++) {
      G4cout << " ----------------- " << G4endl;
      G4cout << "    Particle No " << particledef.at(i).particleid
             << (idx == i ? " Accepted " : " ") << particledef.at(i).name
             << G4endl;
      G4cout << "    Particle I Energy "
             << G4BestUnit(particledef.at(i).IncidentEnergy, "Energy")
             << G4endl;
      G4cout << "    Particle E Energy "
             << G4BestUnit(particledef.at(i).EntryEnergy, "Energy") << G4endl;
      G4cout << "    firstDepo[i][0] " << G4BestUnit(firstDepo[i][0], "Time")
             << G4endl;
      G4cout << "    firstDepo[i][1] " << G4BestUnit(firstDepo[i][1], "Time")
             << G4endl;
      G4cout << "    depoEnergy[i][0] "
             << G4BestUnit(depoEnergy[i][0], "Energy") << G4endl;
      G4cout << "    depoEnergy[i][1] "
             << G4BestUnit(depoEnergy[i][1], "Energy") << G4endl;
      G4cout << "    firstInteraction[i][0] "
             << G4BestUnit(firstInteraction[i][0], "Time") << G4endl;
      G4cout << "    firstInteraction[i][1] "
             << G4BestUnit(firstInteraction[i][1], "Time") << G4endl;
      G4cout << "    reacEnergy[i][0] "
             << G4BestUnit(reacEnergy[i][0], "Energy") << G4endl;
      G4cout << "    reacEnergy[i][1] "
             << G4BestUnit(reacEnergy[i][1], "Energy") << G4endl;
      G4cout << " ----------------- " << G4endl;
    }
    G4cout << " \\\\\\\\\\\\\\ " << G4endl;
  }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
