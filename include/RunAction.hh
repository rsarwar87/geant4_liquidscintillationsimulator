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
/// \file RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id: RunAction.hh 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Timer.hh"
#include "globals.hh"
#include "SteppingAction.hh"
#include <map>
#include <vector>
#include <armadillo>
using namespace arma;
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#define MAX_MULTIPLICITY 64
#define MAX_ARRAY 512
class DetectorConstruction;
class Run;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction: public G4UserRunAction
{
public:
	RunAction(G4String, DetectorConstruction*, PrimaryGeneratorAction*);
	~RunAction();

public:
	virtual G4Run* GenerateRun();
	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);

	static double cutOu;
private:
	DetectorConstruction* fDetector;
	PrimaryGeneratorAction* fPrimary;
	Run* fRun;
	HistoManager* fHistoManager;
	G4Timer* fTimer;
	G4String filename;
	void multiplicity_spec(const G4Run* aRun);
	void multiplicity_th();

	void multiplicity(std::vector<double>* trig_time,
			std::vector<std::vector<double>>* det_idx,
			std::vector<std::vector<G4int>> *dix, int whch);

	void loadBar(int x, int n, int r, int w, int nps);
	void unique1(mat &fields, cube &fields2, int *counter, int val, int ipt);
	void unique2(cube &fields, int cnt, int ipt, int tag, int val,
			double energy, bool rosi);
	uint issue_npar_id(int pid, int nps);
	void initialize_clear(bool init = false);
	mat normalize(mat fields);
	mat calc_moments(mat fields);
	int cnt_particle(cube fields, int cnt, int ipt);
	bool update_base(mat &summary, mat &xt_db);
	int cnt_XT(cube fields, int cnt, int limit, cube fields2, int cnt2,
			int ipt);
	int find_idx(cube &fields, int cnt, int val, int ipt);
	int find_idx2(cube &fields, int id, int val, int ipt);
	void fynamY_calc();
	void fynamY(double tme, double energy, int pid);
	void Print_Multi(string fn, bool VDU, bool rossi);

	cube unique_particle_in_cell, unique_cell_in_particle;
	mat unique_cell_history, unique_particle_history;
	mat summary_, summary_accidental, xt_db_, xt_db_accidental;
	mat rosi_distribution_;

	double new_time;
	double old_time;
	vector<uint>* npar_list;
	int cnt_npar[2], cnt_ncell[2];
	double _cuttoff = 0;
	bool started;
	double _det_window;
	double time_start = 0;
	double tg_time;
	int stage_tm = 0;

	int section_var[2][MAX_ARRAY];
	int64_t fyn_var_cnt[MAX_ARRAY];
	int sec_var_size;
	double fyn_sum;
	double fyn_time[MAX_ARRAY];
	double fyn_mean[MAX_ARRAY];
	int fyn_reminder[MAX_ARRAY];
	double fyn_mean_auto[MAX_ARRAY];
	double fyn_mean_math[MAX_ARRAY];
	double fyn_var_auto[MAX_ARRAY];
	double fyn_var_math[MAX_ARRAY];
	double fyn_Y_auto[MAX_ARRAY];
	double fyn_Y_math[MAX_ARRAY];
	unsigned long int fyn_var_dis_size[MAX_ARRAY];
	int8_t *fyn_var_dist[MAX_ARRAY];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
