
#include <DetectorConstruction.hh>
#include "G4Box.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"

#include "G4NistElementBuilder.hh"
#include "G4NistManager.hh"
#include "G4NistMaterialBuilder.hh"

#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Scintillation.hh"
#include "G4SolidStore.hh"

#include "G4Isotope.hh"
#include "G4MaterialTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Sphere.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DetectorConstruction::hist[10][100] = {{0}};
G4int DetectorConstruction::nb_cryst = 8;
G4double DetectorConstruction::ring_R1 = 20.25 * cm;
G4double DetectorConstruction::ring_W1 = 0 * cm;
bool DetectorConstruction::lancs = false;

G4ThreeVector DetectorConstruction::source_position;
DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(), fCheckOverlaps(true) {
  fExpHall_x = 6 * m;
  fExpHall_y = 4 * m;
  fExpHall_z = 4.3 * m;
  fTank_x = fTank_z = 91.76 * cm;
  fTank_y = 133.6 * cm;
  fBubble_x = fBubble_z = 58.56 * cm;
  fBubble_y = 121.8 * cm;
  buldge_thickness = 0.2 * cm;
  inner_linning_thickness = 0.05 * cm;
  outter_linning_thickness = 0.05 * cm;
  outter_steel = 0.8 * cm;
  polly_thick[0] = 0.5 * cm;
  polly_thick[1] = 7.5 * cm;
  polly_thick[2] = 5 * cm;
  polly_thick[3] = 2.5 * cm;
  det_place[0] = 64.04 * cm;
  det_place[1] = 24.436 * cm;
  det_place[2] = 12.218 * cm;
  det_place[3] = 0 * cm;
  det_place[4] = -12.218 * cm;
  det_place[5] = -24.436 * cm;

  lsd_place[0] = 60.06 * cm;
  lsd_place[1] = 22.53 * cm;
  lsd_place[2] = 7.5 * cm;
  lsd_place[3] = 0 * cm;
  lsd_place[4] = -7.5 * cm;
  lsd_place[5] = -22.53 * cm;

  source_height = source_radius = 0;
  visAttHide = new G4VisAttributes(false);

  /*bfile.open("dump.bin", std::ios::in | std::ios::binary);

  if (!bfile.is_open())
  {
          G4cout << "FIF Failed"<<G4endl;
          exit(-1);
  }
  G4cout << "FIF initialized"<<G4endl;*/

  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {}

void DetectorConstruction::DefineMaterials() {
  G4NistManager* nist = G4NistManager::Instance();
  G4double a;  // atomic mass
  G4double z;  // atomic number

  G4int polyPMMA = 1;
  G4int nC_PMMA = 3 + 2 * polyPMMA;
  G4int nH_PMMA = 6 + 2 * polyPMMA;
  G4int polyeth = 1;
  G4int nC_eth = 2 * polyeth;
  G4int nH_eth = 4 * polyeth;
  //
  // ------------ Generate & Add Material Properties Table ------------
  //
  G4double photonEnergy[] = {
      2.034 * eV, 2.068 * eV, 2.103 * eV, 2.139 * eV, 2.177 * eV, 2.216 * eV,
      2.256 * eV, 2.298 * eV, 2.341 * eV, 2.386 * eV, 2.433 * eV, 2.481 * eV,
      2.532 * eV, 2.585 * eV, 2.640 * eV, 2.697 * eV, 2.757 * eV, 2.820 * eV,
      2.885 * eV, 2.954 * eV, 3.026 * eV, 3.102 * eV, 3.181 * eV, 3.265 * eV,
      3.353 * eV, 3.446 * eV, 3.545 * eV, 3.649 * eV, 3.760 * eV, 3.877 * eV,
      4.002 * eV, 4.136 * eV};

  const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

  G4double refractiveIndex1[] = {
      1.3435, 1.344,  1.3445, 1.345,  1.3455, 1.346,  1.3465, 1.347,
      1.3475, 1.348,  1.3485, 1.3492, 1.35,   1.3505, 1.351,  1.3518,
      1.3522, 1.3530, 1.3535, 1.354,  1.3545, 1.355,  1.3555, 1.356,
      1.3568, 1.3572, 1.358,  1.3585, 1.359,  1.3595, 1.36,   1.3608};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

  G4double absorption[] = {
      3.448 * m,  4.082 * m,  6.329 * m,  9.174 * m,  12.346 * m, 13.889 * m,
      15.152 * m, 17.241 * m, 18.868 * m, 20.000 * m, 26.316 * m, 35.714 * m,
      45.455 * m, 47.619 * m, 52.632 * m, 52.632 * m, 55.556 * m, 52.632 * m,
      52.632 * m, 47.619 * m, 45.455 * m, 41.667 * m, 37.037 * m, 33.333 * m,
      30.000 * m, 28.500 * m, 27.000 * m, 24.500 * m, 22.000 * m, 19.500 * m,
      17.500 * m, 14.500 * m};

  G4double wls_Energy[] = {2.00 * eV, 2.87 * eV, 2.90 * eV, 3.47 * eV};
  const G4int wlsnum = sizeof(wls_Energy) / sizeof(G4double);

  G4double AbsFiber[] = {9.00 * m, 9.00 * m, 0.1 * mm, 0.1 * mm};
  assert(sizeof(AbsFiber) == sizeof(wls_Energy));

  // fiber
  G4double EmissionFib[] = {1.0, 1.0, 0.0, 0.0};
  assert(sizeof(EmissionFib) == sizeof(wls_Energy));
  G4double RefractiveIndexFiber[] = {1.60, 1.60, 1.60, 1.60};
  assert(sizeof(RefractiveIndexFiber) == sizeof(wls_Energy));
  /////////////////////////////////////////////////////////
  //***Elements
  G4String symbol;
  G4double density;
  G4int Z, A, n_iso;

  fN = new G4Element(symbol = "N", symbol = "N", n_iso = 2);
  fiN = new G4Isotope(symbol = "N", Z = 7, A = 14);
  fN->AddIsotope(fiN, 99.6 * perCent);
  G4Isotope* fiN15 = new G4Isotope(symbol = "N", Z = 7, A = 15);
  fN->AddIsotope(fiN15, 99.6 * perCent);

  fO = new G4Element(symbol = "O", symbol = "O", n_iso = 3);
  fiO = new G4Isotope(symbol = "O", Z = 8, A = 16);
  G4Isotope* fiO17 = new G4Isotope(symbol = "O", Z = 8, A = 17);
  G4Isotope* fiO18 = new G4Isotope(symbol = "O", Z = 8, A = 18);
  fO->AddIsotope(fiO, 99.76 * perCent);
  fO->AddIsotope(fiO17, 0.04 * perCent);
  fO->AddIsotope(fiO18, 0.2 * perCent);

  fH = new G4Element(symbol = "H", symbol = "H", n_iso = 2);
  fiH = new G4Isotope(symbol = "H", Z = 1, A = 1);
  G4Isotope* fiH2 = new G4Isotope(symbol = "H", Z = 1, A = 2);
  fH->AddIsotope(fiH, 99.98 * perCent);
  fH->AddIsotope(fiH2, 0.02 * perCent);

  fC = new G4Element(symbol = "C", symbol = "C", n_iso = 2);
  fiC = new G4Isotope(symbol = "C", Z = 6, A = 12);
  G4Isotope* fiC2 = new G4Isotope(symbol = "C", Z = 6, A = 13);
  fC->AddIsotope(fiC, 98.9 * perCent);
  fC->AddIsotope(fiC2, 1.1 * perCent);

  fCa = new G4Element(symbol = "Ca", symbol = "Ca", n_iso = 6);
  fiCa = new G4Isotope(symbol = "Ca", Z = 20, A = 40);
  G4Isotope* fiCa2 = new G4Isotope(symbol = "Ca", Z = 20, A = 42);
  G4Isotope* fiCa3 = new G4Isotope(symbol = "Ca", Z = 20, A = 43);
  G4Isotope* fiCa4 = new G4Isotope(symbol = "Ca", Z = 20, A = 44);
  G4Isotope* fiCa5 = new G4Isotope(symbol = "Ca", Z = 20, A = 46);
  G4Isotope* fiCa6 = new G4Isotope(symbol = "Ca", Z = 20, A = 48);
  fCa->AddIsotope(fiCa, 96.941 * perCent);
  fCa->AddIsotope(fiCa2, .647 * perCent);
  fCa->AddIsotope(fiCa3, .135 * perCent);
  fCa->AddIsotope(fiCa4, 2.086 * perCent);
  fCa->AddIsotope(fiCa5, .004 * perCent);
  fCa->AddIsotope(fiCa6, .187 * perCent);

  fS = new G4Element(symbol = "S", symbol = "S", n_iso = 4);
  fiS = new G4Isotope(symbol = "S", Z = 16, A = 32);
  G4Isotope* fiS2 = new G4Isotope(symbol = "S", Z = 16, A = 33);
  G4Isotope* fiS3 = new G4Isotope(symbol = "S", Z = 16, A = 34);
  G4Isotope* fiS4 = new G4Isotope(symbol = "S", Z = 16, A = 36);
  fS->AddIsotope(fiS, 94.99 * perCent);
  fS->AddIsotope(fiS2, .75 * perCent);
  fS->AddIsotope(fiS3, 4.25 * perCent);
  fS->AddIsotope(fiS4, .01 * perCent);

  fLi = new G4Element(symbol = "L", symbol = "Li", n_iso = 2);
  fiLi = new G4Isotope(symbol = "Li", Z = 3, A = 6);
  G4Isotope* fiLi2 = new G4Isotope(symbol = "Li", Z = 3, A = 7);
  fLi->AddIsotope(fiLi, 7.59 * perCent);
  fLi->AddIsotope(fiLi2, 92.41 * perCent);

  fAm = new G4Element("Am", "Am", z = 92., a = 241 * g / mole);

  fiU233 = new G4Isotope(symbol = "U", Z = 92, A = 234);
  fiU235 = new G4Isotope(symbol = "U", Z = 92, A = 235);
  fiU236 = new G4Isotope(symbol = "U", Z = 92, A = 236);
  fiU238 = new G4Isotope(symbol = "U", Z = 92, A = 238);
  fU = new G4Element(symbol = "U", symbol = "U", n_iso = 4);
  fU->AddIsotope(fiU233, .149 * perCent);
  fU->AddIsotope(fiU235, U235p * perCent);
  fU->AddIsotope(fiU236, .197 * perCent);
  fU->AddIsotope(fiU238, (100 - .149 - .197 - U235p) * perCent);
  fUOx = new G4Material("UOx", density = 12.69 * g / cm3, 2);
  fUOx->AddElement(fU, 84.5 * perCent);
  fUOx->AddElement(fO, 15.5 * perCent);

  Z = 92;                                        // 94;
  fiPu238 = new G4Isotope(symbol = "Pu", Z, A);  // = 238);
  fiPu239 = new G4Isotope(symbol = "Pu", Z, A);  // = 239);
  fiPu240 = new G4Isotope(symbol = "Pu", Z, A);  // = 240);
  fiPu241 = new G4Isotope(symbol = "Pu", Z, A);  // = 241);
  fiPu242 = new G4Isotope(symbol = "Pu", Z, A);  // = 242);
  fPu = new G4Element(symbol = "U", symbol = "U", n_iso = 5);
  fPu->AddIsotope(fiPu238, .149 * perCent);
  fPu->AddIsotope(fiPu239, U235p * perCent);
  fPu->AddIsotope(fiPu240, .197 * perCent);
  fPu->AddIsotope(fiPu241, (100 - .149 - .197 - U235p) * perCent);
  fPu->AddIsotope(fiPu242, .2 * perCent);

  Z = 92;                                        // 94;
  fiCm242 = new G4Isotope(symbol = "Cm", Z, A);  // = 242);
  fiCm244 = new G4Isotope(symbol = "Cm", Z, A);  // = 244);
  fiCm248 = new G4Isotope(symbol = "Cm", Z, A);  // = 248);
  fiCm246 = new G4Isotope(symbol = "Cm", Z, A);  // = 246);
  fCm = new G4Element(symbol = "Cm", symbol = "Cm", n_iso = 4);
  fCm->AddIsotope(fiCm242, 22 * perCent);
  fCm->AddIsotope(fiCm244, 70 * perCent);
  fCm->AddIsotope(fiCm246, 5 * perCent);
  fCm->AddIsotope(fiCm248, 3 * perCent);

  fPu_src = new G4Material("POx", density = 12.69 * g / cm3, 3);
  fPu_src->AddElement(fPu, 80.0 * perCent);
  fPu_src->AddElement(fCm, 4.5 * perCent);
  fPu_src->AddElement(fO, 15.5 * perCent);

  fiCf249 = new G4Isotope(symbol = "Cf", Z, A);  // = 249);
  fiCf250 = new G4Isotope(symbol = "Cf", Z, A);  // = 250);
  fiCf251 = new G4Isotope(symbol = "Cf", Z, A);  // = 251);
  fiCf252 = new G4Isotope(symbol = "Cf", Z, A);  // = 252);
  fiCf254 = new G4Isotope(symbol = "Cf", Z, A);  // = 254);
  fiCf256 = new G4Isotope(symbol = "Cf", Z, A);  // = 254);
  fCf = new G4Element(symbol = "Cf", symbol = "Cf", n_iso = 6);
  fCf->AddIsotope(fiCf249, 3.411 * perCent);
  fCf->AddIsotope(fiCf250, 8.702 * perCent);
  fCf->AddIsotope(fiCf251, 2.6 * perCent);
  fCf->AddIsotope(fiCf252, 85.273 * perCent);
  fCf->AddIsotope(fiCf254, .004 * perCent);
  fCf->AddIsotope(fiCf256, .01 * perCent);

  fCf_src = new G4Material("Cf-252", density = 12.69 * g / cm3, 2);
  fCf_src->AddElement(fCf, 96.625 * perCent);
  fCf_src->AddElement(fCm, (100 - 96.625) * perCent);

  fAmLi = new G4Material("AmLi", density = 12 * g / cm3, 2);
  fAmLi->AddElement(fAm, 1);
  fAmLi->AddElement(fLi, 3);
  //***Materials
  // Aluminum
  fAl = new G4Material("Al", z = 13., a = 26.98 * g / mole,
                       density = 2.7 * g / cm3);

  // Vacuum
  fVacuum = new G4Material("Vacuum", z = 1., a = 1.01 * g / mole,
                           density = universe_mean_density, kStateGas,
                           0.1 * kelvin, 1.e-19 * pascal);
  // Steel
  fSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  // concrete
  fConcrete = nist->FindOrBuildMaterial("G4_CONCRETE");
  // plaster
  fplaster = new G4Material("Plaster", density = 0.85 * g / cm3, 3);
  fplaster->AddElement(fO, 4);
  fplaster->AddElement(fCa, 1);
  fplaster->AddElement(fS, 1);
  // source
  // wood
  fWood = new G4Material("Wood", density = 0.90 * g / cm3, 4);
  fWood->AddElement(fC, 50 * perCent);
  fWood->AddElement(fO, 42 * perCent);
  fWood->AddElement(fH, 7 * perCent);
  fWood->AddElement(fN, 1 * perCent);
  // lead
  fLead = nist->FindOrBuildMaterial("G4_LEAD_OXIDE");

  // water
  /*fwater = new G4Material("Water", density = 1.0 * g / cm3, 2);
   fwater->AddElement(fH, 2);
   fwater->AddElement(fO, 1);*/

  fwater = new G4Material("Water_ts", 1.000 * g / cm3, 2, kStateLiquid,
                          300 * kelvin, 150 * bar);
  // fwater->AddElement(fH, 2);
  fwater->AddElement(fO, 1);
  G4Element* H = new G4Element("TS_H_of_Water", "H", 1., 1.0079 * g / mole);
  fwater->AddElement(H, 2);
  // fwater->AddElement(O, 1);
  fwater->GetIonisation()->SetMeanExcitationEnergy(78.0 * eV);

  G4double energy_water[] = {
      1.56962 * eV, 1.58974 * eV, 1.61039 * eV, 1.63157 * eV, 1.65333 * eV,
      1.67567 * eV, 1.69863 * eV, 1.72222 * eV, 1.74647 * eV, 1.77142 * eV,
      1.7971 * eV,  1.82352 * eV, 1.85074 * eV, 1.87878 * eV, 1.90769 * eV,
      1.93749 * eV, 1.96825 * eV, 1.99999 * eV, 2.03278 * eV, 2.06666 * eV,
      2.10169 * eV, 2.13793 * eV, 2.17543 * eV, 2.21428 * eV, 2.25454 * eV,
      2.29629 * eV, 2.33962 * eV, 2.38461 * eV, 2.43137 * eV, 2.47999 * eV,
      2.53061 * eV, 2.58333 * eV, 2.63829 * eV, 2.69565 * eV, 2.75555 * eV,
      2.81817 * eV, 2.88371 * eV, 2.95237 * eV, 3.02438 * eV, 3.09999 * eV,
      3.17948 * eV, 3.26315 * eV, 3.35134 * eV, 3.44444 * eV, 3.54285 * eV,
      3.64705 * eV, 3.75757 * eV, 3.87499 * eV, 3.99999 * eV, 4.13332 * eV,
      4.27585 * eV, 4.42856 * eV, 4.59258 * eV, 4.76922 * eV, 4.95999 * eV,
      5.16665 * eV, 5.39129 * eV, 5.63635 * eV, 5.90475 * eV, 6.19998 * eV};
  const G4int numentries_water = sizeof(energy_water) / sizeof(G4double);
  G4double mie_water[] = {
      167024.4 * m, 158726.7 * m, 150742 * m,   143062.5 * m, 135680.2 * m,
      128587.4 * m, 121776.3 * m, 115239.5 * m, 108969.5 * m, 102958.8 * m,
      97200.35 * m, 91686.86 * m, 86411.33 * m, 81366.79 * m, 76546.42 * m,
      71943.46 * m, 67551.29 * m, 63363.36 * m, 59373.25 * m, 55574.61 * m,
      51961.24 * m, 48527.00 * m, 45265.87 * m, 42171.94 * m, 39239.39 * m,
      36462.50 * m, 33835.68 * m, 31353.41 * m, 29010.30 * m, 26801.03 * m,
      24720.42 * m, 22763.36 * m, 20924.88 * m, 19200.07 * m, 17584.16 * m,
      16072.45 * m, 14660.38 * m, 13343.46 * m, 12117.33 * m, 10977.70 * m,
      9920.416 * m, 8941.407 * m, 8036.711 * m, 7202.470 * m, 6434.927 * m,
      5730.429 * m, 5085.425 * m, 4496.467 * m, 3960.210 * m, 3473.413 * m,
      3032.937 * m, 2635.746 * m, 2278.907 * m, 1959.588 * m, 1675.064 * m,
      1422.710 * m, 1200.004 * m, 1004.528 * m, 833.9666 * m, 686.1063 * m};

  assert(sizeof(mie_water) == sizeof(energy_water));

  // gforward, gbackward, forward backward ratio
  G4double mie_water_const[3] = {0.99, 0.99, 0.8};
  // gforward, gbackward, forward backward ratio
  G4MaterialPropertiesTable* matH2O = new G4MaterialPropertiesTable();

  matH2O->AddProperty("RINDEX", photonEnergy, refractiveIndex1, nEntries)
      ->SetSpline(true);
  matH2O->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries)
      ->SetSpline(true);
  matH2O->AddProperty("MIEHG", energy_water, mie_water, numentries_water)
      ->SetSpline(true);
  matH2O->AddConstProperty("MIEHG_FORWARD", mie_water_const[0]);
  matH2O->AddConstProperty("MIEHG_BACKWARD", mie_water_const[1]);
  matH2O->AddConstProperty("MIEHG_FORWARD_RATIO", mie_water_const[2]);
  // fwater->SetMaterialPropertiesTable(matH2O);
  // fwater->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  // HE3

  G4Isotope* He3 = new G4Isotope("He3", z = 2, 3, a = 235.01 * g / mole);
  G4Element* eHe3 = new G4Element("He3Det", "He3", 1);
  eHe3->AddIsotope(He3, 100. * perCent);
  fHe3 = new G4Material("Plaster", density = 0.0495 * kg / m3, 1);
  fHe3->AddElement(eHe3, 1);
  fCd = nist->FindOrBuildMaterial("G4_Cd");
  fPb = nist->FindOrBuildMaterial("G4_Pb");

  // Pyrex
  fPyrex = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
  fPyrex->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
  fPyrex->SetMaterialPropertiesTable(matH2O);
  // Glass
  fGlass = new G4Material("Glass", density = 1.032 * g / cm3, 2);
  fGlass->AddElement(fC, 91.533 * perCent);
  fGlass->AddElement(fH, 8.467 * perCent);
  fGlass->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
  fGlass->SetMaterialPropertiesTable(matH2O);

  // Pstyrene
  fPstyrene = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
  fPstyrene->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
  fPstyrene->SetMaterialPropertiesTable(matH2O);
  // Double cladding(flourinated polyethylene)
  fPethylene2 = new G4Material("Pethylene2", density = 1400 * kg / m3, 2);
  fPethylene2->AddElement(fH, nH_eth);
  fPethylene2->AddElement(fC, nC_eth);
  G4double RefractiveIndexPethylene2[] = {1.42, 1.42, 1.42, 1.42};
  assert(sizeof(RefractiveIndexPethylene2) == sizeof(wls_Energy));
  G4MaterialPropertiesTable* Pethylene2Properties =
      new G4MaterialPropertiesTable();
  Pethylene2Properties->AddProperty("RINDEX", wls_Energy,
                                    RefractiveIndexPethylene2, wlsnum);
  Pethylene2Properties->AddProperty("ABSLENGTH", wls_Energy, AbsFiber, wlsnum);
  fPethylene2->SetMaterialPropertiesTable(Pethylene2Properties);
  // Cladding(polyethylene)
  fPethylene1 = new G4Material("Pethylene1", density = 1200 * kg / m3, 2);
  fPethylene1->AddElement(fH, nH_eth);
  fPethylene1->AddElement(fC, nC_eth);
  G4double RefractiveIndexPethylene1[] = {1.49, 1.49, 1.49, 1.49};
  assert(sizeof(RefractiveIndexPethylene1) == sizeof(wls_Energy));
  G4MaterialPropertiesTable* Pethylene1Properties =
      new G4MaterialPropertiesTable();
  Pethylene1Properties->AddProperty("RINDEX", wls_Energy,
                                    RefractiveIndexPethylene1, wlsnum);
  Pethylene1Properties->AddProperty("ABSLENGTH", wls_Energy, AbsFiber, wlsnum);
  fPethylene1->SetMaterialPropertiesTable(Pethylene1Properties);
  // Air
  fAir = nist->FindOrBuildMaterial("G4_AIR");
  G4double refractiveIndex2[] = {
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};

  G4MaterialPropertiesTable* matAir = new G4MaterialPropertiesTable();
  matAir->AddProperty("RINDEX", photonEnergy, refractiveIndex2, nEntries);
  fAir->SetMaterialPropertiesTable(matAir);

  // photocathod + scintillator
  // Fiber(PMMA)
  fPMMA = new G4Material("PMMA", density = 1190 * kg / m3, 3);
  fPMMA->AddElement(fH, 52);
  fPMMA->AddElement(fC, 43);
  fPMMA->AddElement(fO, 18);
  /*ffiber_mat->AddProperty("RINDEX",wls_Energy,RefractiveIndexFiber,wlsnum);
   ffiber_mat->AddProperty("WLSABSLENGTH",wls_Energy,AbsFiber,wlsnum);
   ffiber_mat->AddProperty("WLSCOMPONENT",wls_Energy,EmissionFib,wlsnum);
   ffiber_mat->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
   fPMMA->SetMaterialPropertiesTable(ffiber_mat);*/
  fscintillator = new G4Material("Scintillator", density = 0.959 * g / cm3, 2,
                                 kStateLiquid);
  fscintillator->AddElement(fH, 5);
  fscintillator->AddElement(fC, 4);
  // fscintillator->AddElement(fO, 6);
  fPyrex->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  G4double scintillator_Energy[] = {3.2204 * eV, 3.0996 * eV, 2.9876 * eV,
                                    2.9173 * eV, 2.8700 * eV, 2.8114 * eV,
                                    2.3616 * eV};
  const G4int scintillatornum = sizeof(scintillator_Energy) / sizeof(G4double);

  G4double SCY_Energy[201] = {
      0.01000 * MeV, 0.1000 * MeV, 0.2000 * MeV, 0.3000 * MeV, 0.4000 * MeV,
      0.5000 * MeV,  0.6000 * MeV, 0.7000 * MeV, 0.8000 * MeV, 0.9000 * MeV,
      1 * MeV,       1.100 * MeV,  1.200 * MeV,  1.300 * MeV,  1.400 * MeV,
      1.500 * MeV,   1.600 * MeV,  1.700 * MeV,  1.800 * MeV,  1.900 * MeV,
      2 * MeV,       2.100 * MeV,  2.200 * MeV,  2.300 * MeV,  2.400 * MeV,
      2.500 * MeV,   2.600 * MeV,  2.700 * MeV,  2.800 * MeV,  2.900 * MeV,
      3 * MeV,       3.100 * MeV,  3.200 * MeV,  3.300 * MeV,  3.400 * MeV,
      3.500 * MeV,   3.600 * MeV,  3.700 * MeV,  3.800 * MeV,  3.900 * MeV,
      4 * MeV,       4.100 * MeV,  4.200 * MeV,  4.300 * MeV,  4.400 * MeV,
      4.500 * MeV,   4.600 * MeV,  4.700 * MeV,  4.800 * MeV,  4.900 * MeV,
      5 * MeV,       5.100 * MeV,  5.200 * MeV,  5.300 * MeV,  5.400 * MeV,
      5.500 * MeV,   5.600 * MeV,  5.700 * MeV,  5.800 * MeV,  5.900 * MeV,
      6 * MeV,       6.100 * MeV,  6.200 * MeV,  6.300 * MeV,  6.400 * MeV,
      6.500 * MeV,   6.600 * MeV,  6.700 * MeV,  6.800 * MeV,  6.900 * MeV,
      7 * MeV,       7.100 * MeV,  7.200 * MeV,  7.300 * MeV,  7.400 * MeV,
      7.500 * MeV,   7.600 * MeV,  7.700 * MeV,  7.800 * MeV,  7.900 * MeV,
      8 * MeV,       8.100 * MeV,  8.200 * MeV,  8.300 * MeV,  8.400 * MeV,
      8.500 * MeV,   8.600 * MeV,  8.700 * MeV,  8.800 * MeV,  8.900 * MeV,
      9 * MeV,       9.100 * MeV,  9.200 * MeV,  9.300 * MeV,  9.400 * MeV,
      9.500 * MeV,   9.600 * MeV,  9.700 * MeV,  9.800 * MeV,  9.900 * MeV,
      10 * MeV,      10.10 * MeV,  10.20 * MeV,  10.30 * MeV,  10.40 * MeV,
      10.50 * MeV,   10.60 * MeV,  10.70 * MeV,  10.80 * MeV,  10.90 * MeV,
      11 * MeV,      11.10 * MeV,  11.20 * MeV,  11.30 * MeV,  11.40 * MeV,
      11.50 * MeV,   11.60 * MeV,  11.70 * MeV,  11.80 * MeV,  11.90 * MeV,
      12 * MeV,      12.10 * MeV,  12.20 * MeV,  12.30 * MeV,  12.40 * MeV,
      12.50 * MeV,   12.60 * MeV,  12.70 * MeV,  12.80 * MeV,  12.90 * MeV,
      13 * MeV,      13.10 * MeV,  13.20 * MeV,  13.30 * MeV,  13.40 * MeV,
      13.50 * MeV,   13.60 * MeV,  13.70 * MeV,  13.80 * MeV,  13.90 * MeV,
      14 * MeV,      14.10 * MeV,  14.20 * MeV,  14.30 * MeV,  14.40 * MeV,
      14.50 * MeV,   14.60 * MeV,  14.70 * MeV,  14.80 * MeV,  14.90 * MeV,
      15 * MeV,      15.10 * MeV,  15.20 * MeV,  15.30 * MeV,  15.40 * MeV,
      15.50 * MeV,   15.60 * MeV,  15.70 * MeV,  15.80 * MeV,  15.90 * MeV,
      16 * MeV,      16.10 * MeV,  16.20 * MeV,  16.30 * MeV,  16.40 * MeV,
      16.50 * MeV,   16.60 * MeV,  16.70 * MeV,  16.80 * MeV,  16.90 * MeV,
      17 * MeV,      17.10 * MeV,  17.20 * MeV,  17.30 * MeV,  17.40 * MeV,
      17.50 * MeV,   17.60 * MeV,  17.70 * MeV,  17.80 * MeV,  17.90 * MeV,
      18 * MeV,      18.10 * MeV,  18.20 * MeV,  18.30 * MeV,  18.40 * MeV,
      18.50 * MeV,   18.60 * MeV,  18.70 * MeV,  18.80 * MeV,  18.90 * MeV,
      19 * MeV,      19.10 * MeV,  19.20 * MeV,  19.30 * MeV,  19.40 * MeV,
      19.50 * MeV,   19.60 * MeV,  19.70 * MeV,  19.80 * MeV,  19.90 * MeV,
      20 * MeV};

  G4double SCY_Electron[201] = {
      123,      1230,     2460,     3690.000, 4920,     6150,     7380.000,
      8610,     9840,     11070,    12300,    13530.00, 14760.00, 15990,
      17220,    18450,    19680,    20910.00, 22140,    23370,    24600,
      25830,    27060.00, 28290.00, 29520.00, 30750,    31980,    33210,
      34440,    35670.00, 36900,    38130,    39360,    40590,    41820.00,
      43050,    44280,    45510,    46740,    47970.00, 49200,    50430.00,
      51660,    52890,    54120.00, 55350,    56580.00, 57810,    59040.00,
      60270.00, 61500,    62730.00, 63960,    65190.00, 66420,    67650,
      68880,    70110,    71340.00, 72570,    73800,    75030,    76260,
      77490.00, 78720,    79950,    81180,    82410,    83640.00, 84870,
      86100,    87330,    88560,    89790.00, 91020,    92250,    93480,
      94710,    95940.00, 97170,    98400,    99630,    100860.0, 102090.0,
      103320,   104550,   105780,   107010.0, 108240.0, 109470,   110700,
      111930,   113160.0, 114390.0, 115620,   116850,   118080.0, 119310.0,
      120540.0, 121770,   123000,   124230,   125460.0, 126690.0, 127920.0,
      129150,   130380,   131610,   132840,   134070,   135300,   136530,
      137760,   138990,   140220,   141450,   142680,   143910,   145140,
      146370,   147600,   148830,   150060,   151290,   152520.0, 153750,
      154980,   156210,   157440,   158670.0, 159900,   161130,   162360,
      163590,   164820.0, 166050,   167280,   168510,   169740,   170970.0,
      172200,   173430,   174660,   175890,   177120.0, 178350,   179580,
      180810,   182040,   183270.0, 184500,   185730,   186960,   188190,
      189420.0, 190650,   191880,   193110,   194340,   195570.0, 196800,
      198030.0, 199260,   200490,   201720.0, 202950,   204180.0, 205410,
      206640,   207870.0, 209100,   210330.0, 211560,   212790,   214020.0,
      215250,   216480.0, 217710,   218940,   220170.0, 221400,   222630.0,
      223860,   225090,   226320.0, 227550,   228780.0, 230010,   231240,
      232470.0, 233700,   234930.0, 236160,   237390,   238620.0, 239850,
      241080.0, 242310,   243540,   244770.0, 246000};
  G4double SCY_Proton[201] = {
      4.55700250495120, 58.2718603538135, 144.245466277281, 257.110173580758,
      396.079060189878, 560.388231958463, 749.296148795582, 962.082970512382,
      1198.04992181157, 1456.51867585940, 1736.83075589639, 2038.34695435885,
      2360.44676899881, 2702.52785550500, 3064.00549614188, 3444.31208393803,
      3842.89662196899, 4259.22423729255, 4692.77570910807, 5143.04701072333,
      5609.54886492490, 6091.80631235992, 6589.35829254840, 7101.75723715646,
      7628.56867517177, 8169.37084963287, 8723.75434557418, 9291.32172885864,
      9871.68719557923, 10464.4762317202, 11069.3252827778, 11685.8814330487,
      12313.8020943040, 12952.7547035729, 13602.4164297713, 14262.4738889143,
      14932.6228676632, 15612.5680549620, 16302.0227815277, 17000.7087669635,
      17708.3558742720, 18424.7018715531, 19149.4922006745, 19882.4797527122,
      20623.4246499620, 21372.0940343304, 22128.2618619176, 22891.7087036116,
      23662.2215515183, 24439.5936310549, 25223.6242185431, 26014.1184641396,
      26810.8872199485, 27613.7468731633, 28422.5191840923, 29237.0311289233,
      30057.1147470900, 30882.6069931048, 31713.3495927271, 32549.1889033415,
      33389.9757784206, 34235.5654359550, 35085.8173307332, 35940.5950303589,
      36799.7660948961, 37663.2019600374, 38530.7778236897, 39402.3725358796,
      40277.8684918803, 41157.1515284664, 42040.1108232038, 42926.6387966873,
      43816.6310176390, 44709.9861107837, 45606.6056674201, 46506.3941586093,
      47409.2588509032, 48315.1097245396, 49223.8593940304, 50135.4230310745,
      51049.7182897258, 51966.6652337514, 52886.1862661157, 53808.2060605272,
      54732.6514949899, 55659.4515872976, 56588.5374324169, 57519.8421417026,
      58453.3007838910, 59388.8503278209, 60326.4295868300, 61265.9791647792,
      62207.4414036557, 63150.7603327097, 64095.8816190799, 65042.7525198630,
      65991.3218355865, 66941.5398650431, 67893.3583614461, 68846.7304898684,
      69801.6107859273, 70757.9551156769, 71715.7206366763, 72674.8657601950,
      73635.3501145255, 74597.1345093695, 75560.1809012655, 76524.4523600290,
      77489.9130361739, 78456.5281292883, 79424.2638573346, 80393.0874268487,
      81362.9670040114, 82333.8716865654, 83305.7714765555, 84278.6372538648,
      85252.4407505272, 86227.1545257903, 87202.7519419080, 88179.2071406427,
      89156.4950204536, 90134.5912143536, 91113.4720684138, 92093.1146208971,
      93073.4965820024, 94054.5963142013, 95036.3928131508, 96018.8656891640,
      97001.9951492235, 97985.7619795204, 98970.1475285051, 99955.1336904338,
      100940.702889396, 101926.838063811, 102913.522651376, 103900.740574453,
      104888.476225893, 105876.714455260, 106865.440555471, 107854.640249821,
      108844.299679387, 109834.405390805, 110824.944324400, 111815.903802666,
      112807.271519083, 113799.035527259, 114791.184230395, 115783.706371051,
      116776.591021216, 117769.827572670, 118763.405727624, 119757.315489636,
      120751.547154793, 121746.091303149, 122740.938790419, 123736.080739908,
      124731.508534683, 125727.213809971, 126723.188445780, 127719.424559734,
      128715.914500124, 129712.650839152, 130709.626366385, 131706.834082388,
      132704.267192554, 133701.919101110, 134699.783405293, 135697.853889711,
      136696.124520850, 137694.589441760, 138693.242966883, 139692.079577039,
      140691.093914558, 141690.280778552, 142689.635120330, 143689.152038941,
      144688.826776852, 145688.654715751, 146688.631372473, 147688.752395046,
      148689.013558850, 149689.410762892, 150689.940026186, 151690.597484246,
      152691.379385669, 153692.282088835, 154693.302058686, 155694.435863615,
      156695.680172431, 157697.031751429, 158698.487461529, 159700.044255513,
      160701.699175334, 161703.449349508, 162705.291990578, 163707.224392659,
      164709.243929050, 165711.348049912, 166713.534280027, 167715.800216608,
      168718.143527181};

  G4double scintillator_SCINT[] = {0.1, 0.65, 0.75, 1.0, 0.8, 0.7, 0.1};
  assert(sizeof(scintillator_SCINT) == sizeof(scintillator_Energy));
  G4double scintillator_RIND[] = {1.59, 1.58, 1.58, 1.57, 1.56, 1.55, 1.54};
  assert(sizeof(scintillator_RIND) == sizeof(scintillator_Energy));
  G4double scintillator_ABSL[] = {35. * cm, 35. * cm, 35. * cm, 35. * cm,
                                  35. * cm, 35. * cm, 35. * cm};
  assert(sizeof(scintillator_ABSL) == sizeof(scintillator_Energy));
  fScintillator_mt = new G4MaterialPropertiesTable();
  fScintillator_mt->AddProperty("FASTCOMPONENT", scintillator_Energy,
                                scintillator_SCINT, scintillatornum);
  fScintillator_mt->AddProperty("SLOWCOMPONENT", scintillator_Energy,
                                scintillator_SCINT, scintillatornum);
  fScintillator_mt->AddProperty("RINDEX", scintillator_Energy,
                                scintillator_RIND, scintillatornum);
  fScintillator_mt->AddProperty("ABSLENGTH", scintillator_Energy,
                                scintillator_ABSL, scintillatornum);
  fScintillator_mt->AddProperty("ELECTRONSCINTILLATIONYIELD", SCY_Energy,
                                SCY_Electron, 200);
  fScintillator_mt->AddProperty("PROTONSCINTILLATIONYIELD", SCY_Energy,
                                SCY_Proton, 200);
  // fScintillator_mt->AddConstProperty("SCINTILLATIONYIELD", 12300. / MeV);
  fScintillator_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
  fScintillator_mt->AddConstProperty("FASTTIMECONSTANT", 3.5 * ns);
  fScintillator_mt->AddConstProperty("SLOWTIMECONSTANT", 32. * ns);
  fScintillator_mt->AddConstProperty("YIELDRATIO", 1);
  fscintillator->SetMaterialPropertiesTable(fScintillator_mt);

  cryst_mat = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
  G4double rIndexPstyrene[] = {1.5, 1.5, 1.5, 1.5};
  assert(sizeof(rIndexPstyrene) == sizeof(wls_Energy));
  G4double absorption1[] = {2. * cm, 2. * cm, 2. * cm, 2. * cm};
  assert(sizeof(absorption1) == sizeof(wls_Energy));
  G4double scintilFast[] = {0.00, 0.00, 1.00, 1.00};
  assert(sizeof(scintilFast) == sizeof(wls_Energy));
  fcryst_mat = new G4MaterialPropertiesTable();
  fcryst_mat->AddProperty("RINDEX", wls_Energy, rIndexPstyrene, wlsnum);
  fcryst_mat->AddProperty("ABSLENGTH", wls_Energy, absorption1, wlsnum);
  fcryst_mat->AddProperty("FASTCOMPONENT", wls_Energy, scintilFast, wlsnum);
  cryst_mat->SetMaterialPropertiesTable(fcryst_mat);
  // fPMMA->SetMaterialPropertiesTable(fcryst_mat);
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {
  G4Box* expHall_box =
      new G4Box("World", fExpHall_x / 2, fExpHall_y / 2, fExpHall_z / 2);

  G4LogicalVolume* expHall_log =
      new G4LogicalVolume(expHall_box, fAir, "World", 0, 0, 0);

  G4VPhysicalVolume* expHall_phys = new G4PVPlacement(
      0, G4ThreeVector(), expHall_log, "World", 0, false, 0, fCheckOverlaps);

  // expHall_log->SetVisAttributes(visAttHide);
  // Room
  /* DefineRoom(expHall_log);

   //Tank
   DefineTank(expHall_log);*/

  // Chamber
  // DefineChamber(expHall_log);

  if (!lancs)
    DefineORNL(expHall_log, 15);
  else
    DefineLANCS(expHall_log);

  // source_position = G4ThreeVector(0 * cm, 0 * cm, 0*cm);
  G4cout << "EJ309 created" << G4endl;
  fPBox = expHall_phys;
  return expHall_phys;
}

void DetectorConstruction::DefineLANCS(G4LogicalVolume* expHall_log) {
  DefineRoom(expHall_log);

  // Tank
  DefineTank(expHall_log);
  new G4PVPlacement(
      0, G4ThreeVector(-140. * cm, 7. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 1),                               // its logical volume
      "Detector",                                          // its name
      expHall_log,                                         // its mother  volume
      false,            // no boolean operation
      3,                // copy number
      fCheckOverlaps);  // checking overlaps*/

  new G4PVPlacement(
      0, G4ThreeVector(-140. * cm, -7. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 2),  // its logical volume
      "Detector",             // its name
      expHall_log,            // its mother  volume
      false,                  // no boolean operation
      4,                      // copy number
      fCheckOverlaps);        // checking overlaps*/
  new G4PVPlacement(
      0, G4ThreeVector(-138. * cm, 20. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 3),  // its logical volume
      "Detector",             // its name
      expHall_log,            // its mother  volume
      false,                  // no boolean operation
      3,                      // copy number
      fCheckOverlaps);        // checking overlaps*/

  new G4PVPlacement(
      0,
      G4ThreeVector(-138. * cm, -20. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 4),                              // its logical volume
      "Detector",                                         // its name
      expHall_log,                                        // its mother  volume
      false,            // no boolean operation
      4,                // copy number
      fCheckOverlaps);  // checking overlaps*/
  new G4PVPlacement(
      0, G4ThreeVector(-136. * cm, 33. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 5),  // its logical volume
      "Detector",             // its name
      expHall_log,            // its mother  volume
      false,                  // no boolean operation
      3,                      // copy number
      fCheckOverlaps);        // checking overlaps*/

  new G4PVPlacement(
      0,
      G4ThreeVector(-136. * cm, -33. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 6),                              // its logical volume
      "Detector",                                         // its name
      expHall_log,                                        // its mother  volume
      false,            // no boolean operation
      4,                // copy number
      fCheckOverlaps);  // checking overlaps*/

  new G4PVPlacement(
      0, G4ThreeVector(-130. * cm, 45. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 7),  // its logical volume
      "Detector",             // its name
      expHall_log,            // its mother  volume
      false,                  // no boolean operation
      3,                      // copy number
      fCheckOverlaps);        // checking overlaps*/

  new G4PVPlacement(
      0,
      G4ThreeVector(-130. * cm, -45. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 8),                              // its logical volume
      "Detector",                                         // its name
      expHall_log,                                        // its mother  volume
      false,            // no boolean operation
      4,                // copy number
      fCheckOverlaps);  // checking overlaps*/

  new G4PVPlacement(
      0, G4ThreeVector(-125. * cm, 55. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 9),  // its logical volume
      "Detector",             // its name
      expHall_log,            // its mother  volume
      false,                  // no boolean operation
      3,                      // copy number
      fCheckOverlaps);        // checking overlaps*/

  new G4PVPlacement(
      0,
      G4ThreeVector(-125. * cm, -55. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 10),                             // its logical volume
      "Detector",                                         // its name
      expHall_log,                                        // its mother  volume
      false,            // no boolean operation
      4,                // copy number
      fCheckOverlaps);  // checking overlaps*/

  new G4PVPlacement(
      0, G4ThreeVector(-120. * cm, 67. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 11),  // its logical volume
      "Detector",              // its name
      expHall_log,             // its mother  volume
      false,                   // no boolean operation
      3,                       // copy number
      fCheckOverlaps);         // checking overlaps*/

  new G4PVPlacement(
      0,
      G4ThreeVector(-120. * cm, -67. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 12),                             // its logical volume
      "Detector",                                         // its name
      expHall_log,                                        // its mother  volume
      false,            // no boolean operation
      4,                // copy number
      fCheckOverlaps);  // checking overlaps*/
  new G4PVPlacement(
      0,
      G4ThreeVector(-117. * cm, -80. * cm, -110.5 * cm),  // rotation,position
      DefineEJ309(false, 13),                             // its logical volume
      "Detector",                                         // its name
      expHall_log,                                        // its mother  volume
      false,            // no boolean operation
      4,                // copy number
      fCheckOverlaps);  // checking overlaps*/

  G4RotationMatrix rotm = G4RotationMatrix();
  rotm.rotateY(-90 * deg);
  new G4PVPlacement(
      G4Transform3D(rotm, G4ThreeVector(-120. * cm, 30. * cm,
                                        -80.5 * cm)),  // rotation,position
      DefineEJ309(false, 14),                          // its logical volume
      "Detector",                                      // its name
      expHall_log,                                     // its mother  volume
      false,                                           // no boolean operation
      4,                                               // copy number
      fCheckOverlaps);                                 // checking overlaps*/
  new G4PVPlacement(
      G4Transform3D(rotm, G4ThreeVector(-120. * cm, -30. * cm,
                                        -80.5 * cm)),  // rotation,position
      DefineEJ309(false, 0),                           // its logical volume
      "Detector",                                      // its name
      expHall_log,                                     // its mother  volume
      false,                                           // no boolean operation
      4,                                               // copy number
      fCheckOverlaps);                                 // checking overlaps*/
                                                       // Chamber
  // DefineChamber(expHall_log);
  // DefineORNL(expHall_log, 15);
}

void DetectorConstruction::DefineORNL(G4LogicalVolume* logicWorld, int nd) {
  // 		Floot and cieling
  G4double Ceilingz = 40 * cm;
  G4Box* solidCeiling =
      new G4Box("CeilingG", 0.5 * fExpHall_x, 0.5 * fExpHall_y,
                0.5 * Ceilingz);  // its size
  G4LogicalVolume* logicCeiling = new G4LogicalVolume(solidCeiling, fConcrete,
                                                      "CeilingLV");  // its name
  // logicCeiling->SetVisAttributes(visAttHide);
  new G4PVPlacement(
      0,                                               // no rotation
      G4ThreeVector(0, 0, 175 * cm + 0.5 * Ceilingz),  // at (0,0,0)
      logicCeiling,                                    // its logical volume
      "Floor",                                         // its name
      logicWorld,                                      // its mother  volume
      false,                                           // no boolean operation
      0,                                               // copy number
      fCheckOverlaps);                                 // checking overlaps

  new G4PVPlacement(
      0,                                                // no rotation
      G4ThreeVector(0, 0, -175 * cm - 0.5 * Ceilingz),  // at (0,0,0)
      logicCeiling,                                     // its logical volume
      "Ceiling",                                        // its name
      logicWorld,                                       // its mother  volume
      false,                                            // no boolean operation
      0,                                                // copy number
      fCheckOverlaps);                                  // checking overlaps

  // Table
  G4double l_table = 150 * cm, h_table = 2.5 * cm, b_table = 120 * cm;
  G4double l_tleg = 10 * cm, h_tleg = 100 * cm, b_tleg = 10 * cm;
  G4Box* solid_tablebase = new G4Box("T_base", 0.5 * l_table, 0.5 * b_table,
                                     0.5 * h_table);  // its size
  G4LogicalVolume* logic_tablebase =
      new G4LogicalVolume(solid_tablebase, fAl, solid_tablebase->GetName());
  G4Box* solid_TableLeg = new G4Box("T_leg", 0.5 * l_tleg, 0.5 * b_tleg,
                                    0.5 * h_tleg);  // its size
  G4LogicalVolume* logicTableLeg =
      new G4LogicalVolume(solid_TableLeg, fAl, solid_TableLeg->GetName());
  G4double table_center = -175 * cm + h_tleg + h_table;
  new G4PVPlacement(
      0,                                                  // no rotation
      G4ThreeVector(0, 0, table_center - 0.5 * h_table),  // at (0,0,0)
      logic_tablebase,                                    // its logical volume
      logic_tablebase->GetName(),                         // its name
      logicWorld,                                         // its mother  volume
      false,            // no boolean operation
      0,                // copy number
      fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
      0,  // no rotation
      G4ThreeVector(l_table / 2 - l_tleg / 2, -b_table / 2 + b_tleg / 2,
                    table_center - h_table - 0.5 * h_tleg),  // at (0,0,0)
      logicTableLeg,             // its logical volume
      logicTableLeg->GetName(),  // its name
      logicWorld,                // its mother  volume
      false,                     // no boolean operation
      0,                         // copy number
      fCheckOverlaps);           // checking overlaps

  new G4PVPlacement(
      0,  // no rotation
      G4ThreeVector(l_table / 2 - l_tleg / 2, b_table / 2 - b_tleg / 2,
                    table_center - h_table - 0.5 * h_tleg),  // at (0,0,0)
      logicTableLeg,             // its logical volume
      logicTableLeg->GetName(),  // its name
      logicWorld,                // its mother  volume
      false,                     // no boolean operation
      0,                         // copy number
      fCheckOverlaps);           // checking overlaps
  new G4PVPlacement(
      0,  // no rotation
      G4ThreeVector(-l_table / 2 + l_tleg / 2, -b_table / 2 + b_tleg / 2,
                    table_center - h_table - 0.5 * h_tleg),  // at (0,0,0)
      logicTableLeg,             // its logical volume
      logicTableLeg->GetName(),  // its name
      logicWorld,                // its mother  volume
      false,                     // no boolean operation
      0,                         // copy number
      fCheckOverlaps);           // checking overlaps

  new G4PVPlacement(
      0,  // no rotation
      G4ThreeVector(-l_table / 2 + l_tleg / 2, b_table / 2 - b_tleg / 2,
                    table_center - h_table - 0.5 * h_tleg),  // at (0,0,0)
      logicTableLeg,             // its logical volume
      logicTableLeg->GetName(),  // its name
      logicWorld,                // its mother  volume
      false,                     // no boolean operation
      0,                         // copy number
      fCheckOverlaps);           // checking overlaps

  new G4PVPlacement(
      0,  // no rotation
      G4ThreeVector(l_tleg / 2, b_table / 2 - b_tleg / 2,
                    table_center - h_table - 0.5 * h_tleg),  // at (0,0,0)
      logicTableLeg,             // its logical volume
      logicTableLeg->GetName(),  // its name
      logicWorld,                // its mother  volume
      false,                     // no boolean operation
      0,                         // copy number
      fCheckOverlaps);           // checking overlaps
  new G4PVPlacement(
      0,  // no rotation
      G4ThreeVector(-l_tleg / 2, b_table / 2 - b_tleg / 2,
                    table_center - h_table - 0.5 * h_tleg),  // at (0,0,0)
      logicTableLeg,             // its logical volume
      logicTableLeg->GetName(),  // its name
      logicWorld,                // its mother  volume
      false,                     // no boolean operation
      0,                         // copy number
      fCheckOverlaps);           // checking overlaps

  new G4PVPlacement(
      0,  // no rotation
      G4ThreeVector(l_tleg / 2, -b_table / 2 + b_tleg / 2,
                    table_center - h_table - 0.5 * h_tleg),  // at (0,0,0)
      logicTableLeg,             // its logical volume
      logicTableLeg->GetName(),  // its name
      logicWorld,                // its mother  volume
      false,                     // no boolean operation
      0,                         // copy number
      fCheckOverlaps);           // checking overlaps
  new G4PVPlacement(
      0,  // no rotation
      G4ThreeVector(-l_tleg / 2, -b_table / 2 + b_tleg / 2,
                    table_center - h_table - 0.5 * h_tleg),  // at (0,0,0)
      logicTableLeg,             // its logical volume
      logicTableLeg->GetName(),  // its name
      logicWorld,                // its mother  volume
      false,                     // no boolean operation
      0,                         // copy number
      fCheckOverlaps);           // checking overlaps

  // detector ring
  //
  // ring
  //
  // define one ring as an "envelope" made of air. This will be filled
  // by the crystals (dautgher volumes). The logical volume of the ring
  // (which contains all daughters) will be then placed many times
  //

  G4double cryst_dX = 10 * cm, cryst_dY = 10 * cm, cryst_dZ = 40 * cm;

  G4int nb_rings = 1;

  // The inner radius of the ring is accomodated such that
  // the nb_cryst crystals do not overlap. The approach presented
  // here is pratically equivalent to ring_R1*twopi = nb_cryst*cryst_dY,
  // i.e. the crystals fill entirely the inner circumference.
  // The outer radius is calculated accordingly, by adding the thickness
  // of the crystals and an extra term 1/cos. The computation might
  // screw up if there are too-few crystals.
  //
  G4double dPhi = twopi / nb_cryst, half_dPhi = 0.5 * dPhi;
  G4double cosdPhi = std::cos(half_dPhi);
  G4double tandPhi = std::tan(half_dPhi);
  //
  G4double ring_R2 = (ring_R1 + cryst_dZ) / cosdPhi;
  G4Tubs* solidRing = new G4Tubs("Ring",          // name
                                 ring_R1,         // inner radius
                                 ring_R2,         // outer radius
                                 0.5 * cryst_dX,  // height
                                 0.,              // start angle
                                 twopi);          // spanning angle

  // det support
  G4double l_sup = 3.8 * cm, h_sup = 30 * cm, b_sup = 3.8 * cm;
  G4Box* solid_support = new G4Box("Support", 0.5 * l_sup, 0.5 * b_sup,
                                   0.5 * h_sup);  // its size
  G4Box* solid_support2 =
      new G4Box("Support", 0.5 * (l_sup - 1), 0.5 * (b_sup - 1),
                0.5 * h_sup);  // its size
  G4LogicalVolume* logic_support =
      new G4LogicalVolume(solid_support, fAl, solid_support->GetName());
  G4LogicalVolume* logic_support2 =
      new G4LogicalVolume(solid_support2, fAir, solid_support->GetName());
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0 * cm),            // at (0,0,0)
                    logic_support2, logic_support->GetName(),  // its name
                    logic_support,  // its mother  volume
                    false, 0, fCheckOverlaps);
  G4Tubs* solidRingSupport = new G4Tubs("SRing",      // name
                                        ring_R1,      // inner radius
                                        ring_R2,      // outer radius
                                        0.5 * l_sup,  // height
                                        0.,           // start angle
                                        twopi);       // spanning angle

  G4LogicalVolume* logicRing = new G4LogicalVolume(solidRing,  // its solid
                                                   fAir,       // its material
                                                   "Ring");    // its name
  G4LogicalVolume* logicRingSup =
      new G4LogicalVolume(solidRingSupport,  // its solid
                          fAir,              // its material
                          "RingSupport");    // its name
  for (G4int icrys = 0; icrys < nb_cryst; icrys++) {
    G4double phi = icrys * dPhi;
    // create a rotation matrix... (identity, by defauly)
    G4RotationMatrix rotm = G4RotationMatrix();
    //... and apply rotations. Notice that all rotations are
    // referred to the mother volume.
    rotm.rotateY(90 * deg);
    rotm.rotateZ(phi);
    // Calculate position with respect to the reference frame
    // of the mother volume
    G4ThreeVector uz = G4ThreeVector(std::cos(phi), std::sin(phi), 0.);
    G4ThreeVector position = (ring_R1 + 0.5 * cryst_dZ) * uz;

    G4Transform3D transform = G4Transform3D(rotm, position);

    // Place the crystal with the appropriate transformation
    new G4PVPlacement(transform,                  // rotation,position
                      DefineEJ309(false, icrys),  // its logical volume
                      "Detector",                 // its name
                      logicRing,                  // its mother  volume
                      false,                      // no boolean operation
                      icrys,                      // copy number
                      fCheckOverlaps);            // checking overlaps

    new G4PVPlacement(transform,          // rotation,position
                      logic_support,      // its logical volume
                      "DetectorSupport",  // its name
                      logicRingSup,       // its mother  volume
                      false,              // no boolean operation
                      icrys,              // copy number
                      fCheckOverlaps);    // checking overlaps
  }

  new G4PVPlacement(
      0,                                                        // no rotation
      G4ThreeVector(0 * cm, 0 * cm, table_center + l_sup / 2),  // at (0,0,0)
      logicRingSup,             // its logical volume
      logicRingSup->GetName(),  // its name
      logicWorld,               // its mother  volume
      false,                    // no boolean operation
      0,                        // copy number
      fCheckOverlaps);          // checking overlaps

  new G4PVPlacement(
      0,  // no rotation
      G4ThreeVector(0 * cm, 0 * cm,
                    table_center + l_sup + cryst_dX / 2),  // at (0,0,0)
      logicRing,                                           // its logical volume
      logicRing->GetName(),                                // its name
      logicWorld,                                          // its mother  volume
      false,            // no boolean operation
      0,                // copy number
      fCheckOverlaps);  // checking overlaps

  // lead
  G4double pb_rad1 = 20 * cm, pb_dz = 20.0 * cm;
  G4double pb_rad2 = 19.5 * cm;
  G4Tubs* pb_tube1 =
      new G4Tubs("pb_shield", pb_rad2, pb_rad1, pb_dz / 2, 0., twopi);
  G4LogicalVolume* volPb =
      new G4LogicalVolume(pb_tube1, fPb, pb_tube1->GetName());
  new G4PVPlacement(
      0,                                                        // no rotation
      G4ThreeVector(0 * cm, 0 * cm, table_center + pb_dz / 2),  // at (0,0,0)
      volPb,             // its logical volume
      volPb->GetName(),  // its name
      logicWorld,        // its mother  volume
      false,             // no boolean operation
      0,                 // copy number
      fCheckOverlaps);   // checking overlaps

  // source
  G4double src_rad1 = .2 * cm, can1_dz = 1 * cm;
  G4double src_rad2 = .1 * cm, can2_dz = 0.7 * cm;
  G4Tubs* logic_can =
      new G4Tubs("logic_can", 0, src_rad1, can1_dz / 2, 0., twopi);
  G4Tubs* logic_src = new G4Tubs("SRC", 0, src_rad2, (can2_dz) / 2, 0., twopi);

  G4LogicalVolume* volSrccan =
      new G4LogicalVolume(logic_can, fAl, logic_can->GetName());
  G4LogicalVolume* volSrc =
      new G4LogicalVolume(logic_src, fCf_src, logic_src->GetName());

  new G4PVPlacement(0,                                 // no rotation
                    G4ThreeVector(0 * cm, 0 * cm, 0),  // at (0,0,0)
                    volSrc,                            // its logical volume
                    volSrccan->GetName(),              // its name
                    volSrccan,                         // its mother  volume
                    false,                             // no boolean operation
                    0,                                 // copy number
                    fCheckOverlaps);                   // checking overlaps

  // support stack
  if (ring_W1 == 0 * cm) {
    G4double tube1_rad1 = 2.65 * cm, tube1_dz = 3.0 * cm;
    G4double tube2_rad1 = 1.55 * cm, tube2_dz = 1.7 * cm;
    G4double tube1_rad2 = 2.6 * cm, tube1_dz2 = 2.95 * cm;
    G4double tube2_rad2 = 1.5 * cm, tube2_dz2 = 1.65 * cm;
    G4double l_sup = 20 * cm, h_sup = 3.8 * cm, b_sup = 3.8 * cm;
    G4Tubs* logic_tube1 =
        new G4Tubs("logic_tube1", 0, tube1_rad1, tube1_dz / 2, 0., twopi);
    G4Tubs* logic_tube2 =
        new G4Tubs("logic_tube2", 0, tube2_rad1, (tube2_dz) / 2, 0., twopi);
    G4Tubs* logic_tube1i =
        new G4Tubs("logic_tube2i", 0, tube1_rad2, tube1_dz2 / 2, 0., twopi);
    G4Tubs* logic_tube2i =
        new G4Tubs("logic_tube2i", 0, tube2_rad2, (tube2_dz2) / 2, 0., twopi);
    G4LogicalVolume* volTube1 =
        new G4LogicalVolume(logic_tube1, fAl, logic_tube1->GetName());
    G4LogicalVolume* volTube2 =
        new G4LogicalVolume(logic_tube2, fAl, logic_tube2->GetName());
    G4LogicalVolume* volTube1i =
        new G4LogicalVolume(logic_tube1i, fAir, logic_tube1i->GetName());
    G4LogicalVolume* volTube2i =
        new G4LogicalVolume(logic_tube2i, fAir, logic_tube2i->GetName());
    new G4PVPlacement(0,                                          // no rotation
                      G4ThreeVector(0 * cm, 0 * cm, -.025 * cm),  // at (0,0,0)
                      volTube1i,             // its logical volume
                      volTube1i->GetName(),  // its name
                      volTube1,              // its mother  volume
                      false,                 // no boolean operation
                      0,                     // copy number
                      fCheckOverlaps);       // checking overlaps

    new G4PVPlacement(0,                                          // no rotation
                      G4ThreeVector(0 * cm, 0 * cm, -.025 * cm),  // at (0,0,0)
                      volTube2i,             // its logical volume
                      volTube2i->GetName(),  // its name
                      volTube2,              // its mother  volume
                      false,                 // no boolean operation
                      0,                     // copy number
                      fCheckOverlaps);       // checking overlaps

    G4Box* solid_support = new G4Box("Support", 0.5 * l_sup, 0.5 * b_sup,
                                     0.5 * h_sup);  // its size
    G4Box* solid_support2 =
        new G4Box("Support", 0.5 * (l_sup), 0.5 * (b_sup - 1),
                  0.5 * (h_sup - 1));  // its size
    G4LogicalVolume* logic_support =
        new G4LogicalVolume(solid_support, fAl, solid_support->GetName());
    G4LogicalVolume* logic_support2 =
        new G4LogicalVolume(solid_support2, fAir, solid_support->GetName());

    new G4PVPlacement(0, G4ThreeVector(0, 0, 0 * cm),            // at (0,0,0)
                      logic_support2, logic_support->GetName(),  // its name
                      logic_support,  // its mother  volume
                      false, 0, fCheckOverlaps);

    new G4PVPlacement(
        0,                                                        // no rotation
        G4ThreeVector(0 * cm, 0 * cm, table_center + h_sup / 2),  // at (0,0,0)
        logic_support,             // its logical volume
        logic_support->GetName(),  // its name
        logicWorld,                // its mother  volume
        false,                     // no boolean operation
        0,                         // copy number
        fCheckOverlaps);           // checking overlaps
    new G4PVPlacement(
        0,  // no rotation
        G4ThreeVector(0 * cm, 0 * cm,
                      table_center + h_sup + tube1_dz / 2),  // at (0,0,0)
        volTube1,             // its logical volume
        volTube1->GetName(),  // its name
        logicWorld,           // its mother  volume
        false,                // no boolean operation
        0,                    // copy number
        fCheckOverlaps);      // checking overlaps
    new G4PVPlacement(0,      // no rotation
                      G4ThreeVector(0 * cm, 0 * cm,
                                    table_center + h_sup + tube1_dz +
                                        tube2_dz / 2),  // at (0,0,0)
                      volTube2,                         // its logical volume
                      volTube2->GetName(),              // its name
                      logicWorld,                       // its mother  volume
                      false,                            // no boolean operation
                      0,                                // copy number
                      fCheckOverlaps);                  // checking overlaps

    new G4PVPlacement(0,  // no rotation
                      G4ThreeVector(0 * cm, 0 * cm,
                                    table_center + h_sup + tube1_dz + tube2_dz +
                                        can1_dz / 2),  // at (0,0,0)
                      volSrccan,                       // its logical volume
                      volSrccan->GetName(),            // its name
                      logicWorld,                      // its mother  volume
                      false,                           // no boolean operation
                      0,                               // copy number
                      fCheckOverlaps);                 // checking overlaps
    source_position =
        G4ThreeVector(0 * cm, 0 * cm,
                      table_center + h_sup + tube1_dz + tube2_dz + can1_dz / 2);
  } else {
    G4double water_dz = (-63 * cm - table_center);
    G4Tubs* logic_can =
        new G4Tubs("water_can", 0, ring_W1, water_dz, 0., twopi);
    G4LogicalVolume* volWater =
        new G4LogicalVolume(logic_can, fwater, logic_can->GetName());

    new G4PVPlacement(0,                                      // no rotation
                      G4ThreeVector(0 * cm, 0 * cm, 0 * cm),  // at (0,0,0)
                      volSrccan,             // its logical volume
                      volSrccan->GetName(),  // its name
                      volWater,              // its mother  volume
                      false,                 // no boolean operation
                      0,                     // copy number
                      fCheckOverlaps);       // checking overlaps
    new G4PVPlacement(0,                     // no rotation
                      G4ThreeVector(0 * cm, 0 * cm, -63 * cm),  // at (0,0,0)
                      volWater,             // its logical volume
                      volWater->GetName(),  // its name
                      logicWorld,           // its mother  volume
                      false,                // no boolean operation
                      0,                    // copy number
                      fCheckOverlaps);      // checking overlaps
    source_position = G4ThreeVector(0 * cm, 0, -63 * cm);
  }

  ;
}
void DetectorConstruction::SurfaceProperties(G4LogicalVolume* fHousing_log2,
                                             G4LogicalVolume* fPhotocath_log) {
  G4double ephoton[] = {2. * eV, 3.47 * eV};
  const G4int num = sizeof(ephoton) / sizeof(G4double);

  //**Scintillator housing properties
  G4double reflectivity[] = {1.35, 1.40};
  assert(sizeof(reflectivity) == sizeof(ephoton));
  G4double efficiency[] = {0.0, 0.0};
  assert(sizeof(efficiency) == sizeof(ephoton));
  G4MaterialPropertiesTable* scintHsngPT = new G4MaterialPropertiesTable();
  scintHsngPT->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
  scintHsngPT->AddProperty("EFFICIENCY", ephoton, efficiency, num);
  G4OpticalSurface* OpScintHousingSurface = new G4OpticalSurface(
      "HousingSurface", unified, polished, dielectric_metal);
  OpScintHousingSurface->SetMaterialPropertiesTable(scintHsngPT);

  //**Photocathode surface properties
  G4double photocath_EFF[] = {.25, .25};  // Enables 'detection' of photons
  assert(sizeof(photocath_EFF) == sizeof(ephoton));
  G4double photocath_ReR[] = {1.92, 1.92};
  assert(sizeof(photocath_ReR) == sizeof(ephoton));
  G4double photocath_ImR[] = {1.69, 1.69};
  assert(sizeof(photocath_ImR) == sizeof(ephoton));
  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY", ephoton, photocath_EFF, num);
  photocath_mt->AddProperty("REALRINDEX", ephoton, photocath_ReR, num);
  photocath_mt->AddProperty("IMAGINARYRINDEX", ephoton, photocath_ImR, num);
  G4OpticalSurface* photocath_opsurf = new G4OpticalSurface(
      "photocath_opsurf", glisur, polished, dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

  //**Create logical skin surfaces
  new G4LogicalSkinSurface("ScintSurface", fHousing_log2,
                           OpScintHousingSurface);
  new G4LogicalSkinSurface("photocath_surf", fPhotocath_log, photocath_opsurf);
}

G4LogicalVolume* DetectorConstruction::DefineEJ309(bool upright, int copyNo) {
  G4double scin_frame_dx = 10 * cm, scin_frame_dy = 10 * cm,
           scin_frame_dz = 40 * cm;
  G4double scin_endcap_dxy = 10 * cm, scin_endcap_dz = 12 * cm;
  G4double scin_room_chamber_dxy = 9 * cm, scin_room_dz = 11 * cm;
  G4double pyrex_glass_rad = 4.5 * cm, pyrex_glass_dz = 1 * cm;
  G4double glass_frame_rad1 = 5 * cm, glass_frame_rad2 = 4.5 * cm,
           glass_frame_dz = 0.75 * cm;
  G4double glass_rad1 = 4.5 * cm, glass_dz = 0.75 * cm, glass_dz1 = 0.5 * cm;
  G4double tube1_rad1 = 4.5 * cm, tube1_dz = 5.0 * cm;
  G4double tube2_rad1 = 3.0 * cm, tube2_dz = 10.0 * cm;
  G4double tube3_rad1 = 3.0 * cm, tube3_dz = 4.5 * cm;
  // Solids
  G4Box* solidDetectorFrame =
      new G4Box("DetectorFrame", scin_frame_dx / 2, scin_frame_dy / 2,
                scin_frame_dz / 2);  // detectorBlock
  G4Box* solidScintillatorChamber =
      new G4Box("ScintillatorChamber", scin_endcap_dxy / 2, scin_endcap_dxy / 2,
                scin_endcap_dz / 2);  // AlEndCap

  G4Tubs* SolidGlassFrame =
      new G4Tubs("GlassFrame", glass_frame_rad2, glass_frame_rad1,
                 glass_frame_dz / 2, 0., twopi);  // detectorBlock
  // G4IntersectionSolid IntersFrames("IntersFrames", SolidGlassFrame,
  // solidDetectorFrame);

  G4Tubs* SolidglassBlock1 =
      new G4Tubs("GlassBlock1", 0, glass_rad1, glass_dz1 / 2, 0., twopi);
  G4Tubs* SolidglassBlock2 =
      new G4Tubs("GlassBlock2", 0, glass_rad1, glass_dz / 2, 0., twopi);
  G4Cons* solidPMTCone = new G4Cons(
      "PMTCone", 0, tube1_rad1,  //   //G4double  pRmin1, G4double  pRmax1,
      0, tube2_rad1, (1.75), 0, 2 * twopi / 3);
  G4Tubs* solidPMTTopTube = new G4Tubs("PMTTopTube", 0, tube3_rad1,
                                       (tube2_dz + tube3_dz) / 2, 0., twopi);
  G4Tubs* solidPMTBottomTube =
      new G4Tubs("PMTBottomTube", 0, tube1_rad1, (tube1_dz) / 2, 0., twopi);

  // logics
  G4LogicalVolume* logicDetectorFrame = new G4LogicalVolume(
      solidDetectorFrame, fAir, solidDetectorFrame->GetName());  // its name
  G4LogicalVolume* logicScintillatorChamber = new G4LogicalVolume(
      solidScintillatorChamber, fAl, solidScintillatorChamber->GetName());
  // G4LogicalVolume* logicDetECap1 = new G4LogicalVolume(solidDetECap, fAl,
  // "DetEndCap");

  G4Box* solidScintillationChamber = NULL;
  G4LogicalVolume* logicScintillationChamber = NULL;
  G4Box* solidEmptyChamber = NULL;
  G4LogicalVolume* logicEmptyChamber = NULL;

  G4LogicalVolume* logicGlassFrame =
      new G4LogicalVolume(SolidGlassFrame, fAl,
                          SolidGlassFrame->GetName());  // its name
  // G4LogicalVolume* logicGlassFrame = new G4LogicalVolume(&IntersFrames, fAl,
  // "GlassFrame");            //its name
  G4LogicalVolume* logicDetGlass1 = new G4LogicalVolume(
      SolidglassBlock1, fPMMA, SolidglassBlock1->GetName());  // its name
  G4LogicalVolume* logicDetGlass2 =
      new G4LogicalVolume(SolidglassBlock2, fPMMA, SolidglassBlock2->GetName());
  G4LogicalVolume* logicPMTCone =
      new G4LogicalVolume(solidPMTCone, fAl, solidPMTCone->GetName());
  G4LogicalVolume* logicPMTSTube =
      new G4LogicalVolume(solidPMTTopTube, fAl,
                          solidPMTTopTube->GetName());  // its name
  G4LogicalVolume* logicPMTBTube = new G4LogicalVolume(
      solidPMTBottomTube, fAl, solidPMTBottomTube->GetName());

  if (upright) {
    solidScintillationChamber = new G4Box(
        "Scintillator", scin_room_chamber_dxy / 2, (scin_room_chamber_dxy) / 2,
        scin_room_dz * 0.6 / 2);  // detectorchamber
    logicScintillationChamber =
        new G4LogicalVolume(solidScintillationChamber, fscintillator,
                            solidScintillationChamber->GetName());
    solidEmptyChamber = new G4Box(
        "EmptyChamber", scin_room_chamber_dxy / 2, (scin_room_chamber_dxy) / 2,
        scin_room_dz * 0.4 / 2);  // detectorchamberEmpty
    logicEmptyChamber = new G4LogicalVolume(solidEmptyChamber, fAir,
                                            solidEmptyChamber->GetName());
    SurfaceProperties(logicScintillatorChamber, logicDetGlass1);

    new G4PVPlacement(0, G4ThreeVector(0, 0, -2.2 * cm),  // at (0,0,0)3,3
                      logicScintillationChamber,
                      solidScintillationChamber->GetName(),  // its name
                      logicScintillatorChamber,  // its mother  volume
                      false, copyNo, fCheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 3.3 * cm),  // at (0,0,0)
                      logicEmptyChamber,
                      solidEmptyChamber->GetName(),  // its name
                      logicScintillatorChamber,      // its mother  volume
                      false, copyNo, fCheckOverlaps);
  } else {
    solidEmptyChamber = new G4Box("EmptyChamber", scin_room_chamber_dxy / 2,
                                  (scin_room_chamber_dxy * 0.4) / 2,
                                  scin_room_dz / 2);  // detectorchamberEmpty
    logicEmptyChamber = new G4LogicalVolume(solidEmptyChamber, fAir,
                                            solidEmptyChamber->GetName());
    solidScintillationChamber =
        new G4Box("Scintillator", scin_room_chamber_dxy / 2,
                  (scin_room_chamber_dxy * 0.6) / 2,
                  scin_room_dz / 2);  // detectorchamber
    logicScintillationChamber =
        new G4LogicalVolume(solidScintillationChamber, fscintillator,
                            solidScintillationChamber->GetName());

    SurfaceProperties(logicScintillatorChamber, logicDetGlass1);
    new G4PVPlacement(0, G4ThreeVector(0, -1.8 * cm, 0),  // at (0,0,0)3,3
                      logicScintillationChamber,
                      logicScintillationChamber->GetName(),  // its name
                      logicScintillatorChamber,  // its mother  volume
                      false, copyNo, fCheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 2.7 * cm, 0),  // at (0,0,0)
                      logicEmptyChamber,
                      logicEmptyChamber->GetName(),  // its name
                      logicScintillatorChamber,      // its mother  volume
                      false, copyNo, fCheckOverlaps);
  }

  new G4PVPlacement(0, G4ThreeVector(0, 0, 5.75 * cm),          // at (0,0,0)
                    logicDetGlass1, logicDetGlass1->GetName(),  // its name
                    logicScintillatorChamber,  // its mother  volume
                    false, copyNo, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0, 0, -14 * cm),  // at (0,0,0)
                    logicScintillatorChamber,
                    logicScintillatorChamber->GetName(),  // its name
                    logicDetectorFrame,                   // its mother  volume
                    false, copyNo, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0, 0, -7.625 * cm),          // at (0,0,0)
                    logicGlassFrame, logicGlassFrame->GetName(),  // its name
                    logicDetectorFrame,  // its mother  volume
                    false, copyNo, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0, 0, -7.625 * cm),        // at (0,0,0)
                    logicDetGlass2, logicDetGlass2->GetName(),  // its name
                    logicDetectorFrame,  // its mother  volume
                    false, copyNo, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0, 0, -4.75 * cm),       // at (0,0,0)
                    logicPMTBTube, logicPMTBTube->GetName(),  // its name
                    logicDetectorFrame,  // its mother  volume
                    false, copyNo, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5 * cm),      // at (0,0,0)
                    logicPMTCone, logicPMTCone->GetName(),  // its name
                    logicDetectorFrame,  // its mother  volume
                    false, copyNo, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 8.5 * cm),         // at (0,0,0)
                    logicPMTSTube, logicPMTSTube->GetName(),  // its name
                    logicDetectorFrame,  // its mother  volume
                    false, copyNo, fCheckOverlaps);
  /*new G4PVPlacement(0,
                  G4ThreeVector(-h_sup / 2 + (scin_frame_dx + h_sup) / 2, 0, -5
     * cm), //at (0,0,0) logic_support, logic_support->GetName(), //its name
                  logicDetectorFrame,              //its mother  volume
                  false, copyNo, fCheckOverlaps);*/

  return logicDetectorFrame;
}

void DetectorConstruction::DefineChamber(G4LogicalVolume* logicWorld) {
  // Buldge
  G4Box* geom_empty = new G4Box("Buldge", 0.5 * fBubble_x, 0.5 * fBubble_y,
                                0.5 * fBubble_z);  // its size
  G4LogicalVolume* logicEmpty =
      new G4LogicalVolume(geom_empty,                    // its solid
                          fAir, geom_empty->GetName());  // its name

  visAttOthers = new G4VisAttributes(true, G4Colour(0., 0., 0., 1.));
  visAttOthers->SetForceSolid(true);
  logicEmpty->SetVisAttributes(visAttOthers);
  // Buldgethickness
  G4double X = fBubble_x + buldge_thickness * 2;
  G4double Z = fBubble_z + buldge_thickness * 2;
  G4Box* geom_empty_frame = new G4Box("BuldgeFrame", 0.5 * (X), 0.5 * fBubble_y,
                                      0.5 * Z);  // its size
  G4LogicalVolume* logicEmptyFrame =
      new G4LogicalVolume(geom_empty_frame, fAl, geom_empty_frame->GetName());
  visAttOthers = new G4VisAttributes(true, G4Colour(1., 0., 0., .3));
  visAttOthers->SetForceSolid(true);
  // logicEmptyFrame->SetVisAttributes(visAttOthers);
  new G4PVPlacement(0,                         // no rotation
                    G4ThreeVector(0, 0, 0.0),  // at (0,0,0)
                    logicEmpty,                // its logical volume
                    logicEmpty->GetName(),     // its name
                    logicEmptyFrame,           // its mother  volume
                    false,                     // no boolean operation
                    0,                         // copy number
                    fCheckOverlaps);           // checking overlaps

  // Buldgethickness
  X += inner_linning_thickness * 2;
  Z += inner_linning_thickness * 2;
  G4Box* geam_inner_linning = new G4Box("InnerLinning", 0.5 * (X),
                                        0.5 * fBubble_y, 0.5 * Z);  // its size
  G4LogicalVolume* logicInnerLinning = new G4LogicalVolume(
      geam_inner_linning, fCd, geam_inner_linning->GetName());
  new G4PVPlacement(0,                           // no rotation
                    G4ThreeVector(0, 0, 0.0),    // at (0,0,0)
                    logicEmptyFrame,             // its logical volume
                    logicEmptyFrame->GetName(),  // its name
                    logicInnerLinning,           // its mother  volume
                    false,                       // no boolean operation
                    0,                           // copy number
                    fCheckOverlaps);             // checking overlaps
  visAttOthers = new G4VisAttributes(true, G4Colour(0., 1., 0., .3));
  visAttOthers->SetForceSolid(true);
  // logicInnerLinning->SetVisAttributes(visAttOthers);
  // Buldgethickness
  X += polly_thick[0] * 2;
  Z += polly_thick[0] * 2;
  G4Box* geam_polly1 = new G4Box("Polly1", 0.5 * (X), 0.5 * fBubble_y,
                                 0.5 * Z);  // its size
  G4LogicalVolume* logicPolly1 =
      new G4LogicalVolume(geam_polly1, fAir, geam_polly1->GetName());
  new G4PVPlacement(0,                             // no rotation
                    G4ThreeVector(0, 0, 0.0),      // at (0,0,0)
                    logicInnerLinning,             // its logical volume
                    logicInnerLinning->GetName(),  // its name
                    logicPolly1,                   // its mother  volume
                    false,                         // no boolean operation
                    0,                             // copy number
                    fCheckOverlaps);               // checking overlaps
  visAttOthers = new G4VisAttributes(true, G4Colour(0., 0., 1., .3));
  visAttOthers->SetForceSolid(true);
  // logicPolly1->SetVisAttributes(visAttOthers);
  // detectors

  X += (polly_thick[1] * 2 + polly_thick[2] * 2);
  Z += (polly_thick[1] * 2 + polly_thick[2] * 2);
  G4Box* geam_polly2 = new G4Box("Polly2", 0.5 * (X), 0.5 * fBubble_y,
                                 0.5 * Z);  // its size
  G4LogicalVolume* logicPolly2 =
      new G4LogicalVolume(geam_polly2, fAir, geam_polly2->GetName());
  new G4PVPlacement(0,                         // no rotation
                    G4ThreeVector(0, 0, 0.0),  // at (0,0,0)
                    logicPolly1,               // its logical volume
                    logicPolly1->GetName(),    // its name
                    logicPolly2,               // its mother  volume
                    false,                     // no boolean operation
                    0,                         // copy number
                    fCheckOverlaps);           // checking overlaps
  int nrow = 2;
  double gap = 15;
  double shift = gap * nrow / 2 - 5;
  if (nrow == 1) {
    gap = 0;
    shift = 0;
  }
  for (int i = 0; i < nrow; i++) {
    new G4PVPlacement(0,
                      G4ThreeVector(36.28 * cm, (-shift + gap * i) * cm,
                                    lsd_place[1]),     // at (0,0,0)
                      DefineEJ309(false, 1 + i * 16),  // its logical volume
                      "LSD", logicPolly2, true, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(36.28 * cm, (0 - shift + gap * i) * cm,
                                    lsd_place[2]),     // at (0,0,0)
                      DefineEJ309(false, 2 + i * 16),  // its logical volume
                      "LSD", logicPolly2, true, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(36.28 * cm, (0 - shift + gap * i) * cm,
                                    lsd_place[5]),     // at (0,0,0)
                      DefineEJ309(false, 3 + i * 16),  // its logical volume
                      "LSD", logicPolly2, true, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(36.28 * cm, (0 - shift + gap * i) * cm,
                                    lsd_place[4]),     // at (0,0,0)
                      DefineEJ309(false, 4 + i * 16),  // its logical volume
                      "LSD", logicPolly2, true, 0, fCheckOverlaps);

    new G4PVPlacement(0,
                      G4ThreeVector(-36.28 * cm, (0 - shift + gap * i) * cm,
                                    lsd_place[1]),     // at (0,0,0)
                      DefineEJ309(false, 5 + i * 16),  // its logical volume
                      "LSD", logicPolly2, true, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(-36.28 * cm, (0 - shift + gap * i) * cm,
                                    lsd_place[2]),     // at (0,0,0)
                      DefineEJ309(false, 6 + i * 16),  // its logical volume
                      "LSD", logicPolly2, true, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(-36.28 * cm, (0 - shift + gap * i) * cm,
                                    lsd_place[4]),     // at (0,0,0)
                      DefineEJ309(false, 7 + i * 16),  // its logical volume
                      "LSD", logicPolly2, true, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(-36.28 * cm, (0 - shift + gap * i) * cm,
                                    lsd_place[5]),     // at (0,0,0)
                      DefineEJ309(false, 8 + i * 16),  // its logical volume
                      "LSD", logicPolly2, true, 0, fCheckOverlaps);

    new G4PVPlacement(0,
                      G4ThreeVector(lsd_place[1], (0 - shift + gap * i) * cm,
                                    36.28 * cm),       // at (0,0,0)
                      DefineEJ309(false, 9 + i * 16),  // its logical volume
                      "LSD", logicPolly2, false, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(lsd_place[2], (0 - shift + gap * i) * cm,
                                    36.28 * cm),        // at (0,0,0)
                      DefineEJ309(false, 10 + i * 16),  // its logical volume
                      "LSD", logicPolly2, false, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(lsd_place[5], (0 - shift + gap * i) * cm,
                                    36.28 * cm),        // at (0,0,0)
                      DefineEJ309(false, 11 + i * 16),  // its logical volume
                      "LSD", logicPolly2, false, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(lsd_place[4], (0 - shift + gap * i) * cm,
                                    36.28 * cm),        // at (0,0,0)
                      DefineEJ309(false, 12 + i * 16),  // its logical volume
                      "LSD", logicPolly2, false, 0, fCheckOverlaps);

    new G4PVPlacement(0,
                      G4ThreeVector(lsd_place[1], (0 - shift + gap * i) * cm,
                                    -36.28 * cm),       // at (0,0,0)
                      DefineEJ309(false, 14 + i * 16),  // its logical volume
                      "LSD", logicPolly2, false, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(lsd_place[2], (0 - shift + gap * i) * cm,
                                    -36.28 * cm),       // at (0,0,0)
                      DefineEJ309(false, 15 + i * 16),  // its logical volume
                      "LSD", logicPolly2, false, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(lsd_place[5], (0 - shift + gap * i) * cm,
                                    -36.28 * cm),       // at (0,0,0)
                      DefineEJ309(false, 16 + i * 16),  // its logical volume
                      "LSD", logicPolly2, false, 0, fCheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(lsd_place[4], (0 - shift + gap * i) * cm,
                                    -36.28 * cm),       // at (0,0,0)
                      DefineEJ309(false, 13 + i * 16),  // its logical volume
                      "LSD", logicPolly2, false, 0, fCheckOverlaps);
  }

  X += outter_linning_thickness * 2;
  Z += outter_linning_thickness * 2;
  G4Box* geam_Outer = new G4Box("OuterLinning", 0.5 * (X), 0.5 * fBubble_y,
                                0.5 * Z);  // its size
  G4LogicalVolume* logicOuterLinning =
      new G4LogicalVolume(geam_Outer, fCd, geam_Outer->GetName());
  new G4PVPlacement(0,                         // no rotation
                    G4ThreeVector(0, 0, 0.0),  // at (0,0,0)
                    logicPolly2,               // its logical volume
                    logicPolly2->GetName(),    // its name
                    logicOuterLinning,         // its mother  volume
                    false,                     // no boolean operation
                    0,                         // copy number
                    fCheckOverlaps);           // checking overlaps

  X += polly_thick[3] * 2;
  Z += polly_thick[3] * 2;
  G4Box* geam_polly4 = new G4Box("Polly4", 0.5 * (X), 0.5 * fBubble_y,
                                 0.5 * Z);  // its size
  G4LogicalVolume* logicPolly4 =
      new G4LogicalVolume(geam_polly4, fAir, geam_polly4->GetName());
  new G4PVPlacement(0,                             // no rotation
                    G4ThreeVector(0, 0, 0.0),      // at (0,0,0)
                    logicOuterLinning,             // its logical volume
                    logicOuterLinning->GetName(),  // its name
                    logicPolly4,                   // its mother  volume
                    false,                         // no boolean operation
                    0,                             // copy number
                    fCheckOverlaps);               // checking overlaps
  // logicPolly4->SetVisAttributes(visAttOthers);

  X += outter_steel * 2;
  Z += outter_steel * 2;
  G4Box* geam_steel =
      new G4Box("Steel", 0.5 * (X), 0.5 * fTank_y, 0.5 * Z);  // its size
  G4LogicalVolume* logicSteel =
      new G4LogicalVolume(geam_steel, fSteel, geam_steel->GetName());
  new G4PVPlacement(0,                         // no rotation
                    G4ThreeVector(0, 0, 0.0),  // at (0,0,0)
                    logicPolly4,               // its logical volume
                    logicPolly4->GetName(),    // its name
                    logicSteel,                // its mother  volume
                    false,                     // no boolean operation
                    0,                         // copy number
                    fCheckOverlaps);           // checking overlaps
  visAttOthers = new G4VisAttributes(true, G4Colour(0., 1., 0., .3));
  visAttOthers->SetForceSolid(true);
  // logicInnerLinning->SetVisAttributes(visAttOthers);
  // logicOuterLinning->SetVisAttributes(visAttOthers);
  visAttOthers = new G4VisAttributes(true, G4Colour(1., 0, 1., .3));
  visAttOthers->SetForceSolid(true);
  logicSteel->SetVisAttributes(visAttOthers);

  G4double Y = 5.1 * cm;
  X = Z = 90.16 * cm;
  G4Box* geam_plate =
      new G4Box("Plate", 0.5 * (X), 0.5 * Y, 0.5 * Z);  // its size
  G4Box* geam_plateCd =
      new G4Box("CdPlate", 0.5 * (X), 0.5 * 0.05, 0.5 * Z);  // its size
  G4LogicalVolume* logicPlateCd =
      new G4LogicalVolume(geam_plateCd, fCd, geam_plate->GetName());
  G4LogicalVolume* logicPlateN =
      new G4LogicalVolume(geam_plate, fPethylene1, geam_plate->GetName());
  new G4PVPlacement(0,                                        // no rotation
                    G4ThreeVector(0, -0.5 * 0.05 * cm, 0.0),  // at (0,0,0)
                    logicPlateCd,                       // its logical volume
                    logicPlateCd->GetName(),            // its name
                    logicPlateN,                        // its mother  volume
                    false,                              // no boolean operation
                    0,                                  // copy number
                    fCheckOverlaps);                    // checking overlaps
  new G4PVPlacement(0,                                  // no rotation
                    G4ThreeVector(0, 2.525 * cm, 0.0),  // at (0,0,0)
                    logicPlateCd,                       // its logical volume
                    logicPlateCd->GetName(),            // its name
                    logicPlateN,                        // its mother  volume
                    false,                              // no boolean operation
                    0,                                  // copy number
                    fCheckOverlaps);                    // checking overlaps
  G4LogicalVolume* logicPlateF =
      new G4LogicalVolume(geam_plate, fPethylene1, geam_plate->GetName());
  new G4PVPlacement(0,                                 // no rotation
                    G4ThreeVector(0, -2.525 * cm, 0),  // at (0,0,0)
                    logicPlateCd,                      // its logical volume
                    logicPlateCd->GetName(),           // its name
                    logicPlateF,                       // its mother  volume
                    false,                             // no boolean operation
                    0,                                 // copy number
                    fCheckOverlaps);                   // checking overlaps
  new G4PVPlacement(0,                                 // no rotation
                    G4ThreeVector(0, 0.5 * 0.05 * cm, 0.0),  // at (0,0,0)
                    logicPlateCd,             // its logical volume
                    logicPlateCd->GetName(),  // its name
                    logicPlateF,              // its mother  volume
                    false,                    // no boolean operation
                    0,                        // copy number
                    fCheckOverlaps);          // checking overlaps77

  new G4PVPlacement(0,                                   // no rotation
                    G4ThreeVector(0, -63.45 * cm, 0.0),  // at (0,0,0)
                    logicPlateN,                         // its logical volume
                    logicPlateN->GetName(),              // its name
                    logicSteel,                          // its mother  volume
                    false,                               // no boolean operation
                    0,                                   // copy number
                    fCheckOverlaps);                     // checking overlaps
  new G4PVPlacement(0,                                   // no rotation
                    G4ThreeVector(0, 63.45 * cm, 0.0),   // at (0,0,0)
                    logicPlateF,                         // its logical volume
                    logicPlateF->GetName(),              // its name
                    logicSteel,                          // its mother  volume
                    false,                               // no boolean operation
                    0,                                   // copy number
                    fCheckOverlaps);                     // checking overlaps

  new G4PVPlacement(0,                       // no rotation
                    G4ThreeVector(0, 0, 0),  // at (0,0,0)
                    logicSteel,              // its logical volume
                    logicSteel->GetName(),   // its name
                    logicWorld,              // its mother  volume
                    false,                   // no boolean operation
                    0,                       // copy number
                    fCheckOverlaps);         // checking overlaps

  return;
}

G4LogicalVolume* DetectorConstruction::DefineHE3(bool upright, int copyNo) {
  G4Box* solidDetectorFrame = new G4Box(
      "DetectorFrame", 5 * cm / 2, fBubble_y / 2, 5 * cm / 2);  // detectorBlock
  G4LogicalVolume* logicDetectorFrame = new G4LogicalVolume(
      solidDetectorFrame, fAir, solidDetectorFrame->GetName());

  G4RotationMatrix rotm = G4RotationMatrix();
  rotm.rotateX(-90 * deg);

  G4Tubs* SolidAlFrame =
      new G4Tubs("AlFrame", 0, 5 * cm / 2, fBubble_y / 2, 0., twopi);
  G4LogicalVolume* logicAlFrame =
      new G4LogicalVolume(SolidAlFrame,                   // its solid
                          fAl, SolidAlFrame->GetName());  // its name
  G4Tubs* SolidHe3 =
      new G4Tubs("He3Avtive", 0, 4.8 * cm / 2, fBubble_y / 2, 0., twopi);
  G4LogicalVolume* logicHe3 =
      new G4LogicalVolume(SolidHe3,                    // its solid
                          fHe3, SolidHe3->GetName());  // its name

  new G4PVPlacement(0,                         // no rotation
                    G4ThreeVector(0, 0, 0.0),  // at (0,0,0)
                    logicHe3,                  // its logical volume
                    logicHe3->GetName(),       // its name
                    logicAlFrame,              // its mother  volume
                    false,                     // no boolean operation
                    copyNo,                    // copy number
                    fCheckOverlaps);           // checking overlaps

  new G4PVPlacement(G4Transform3D(rotm, G4ThreeVector(0, 0, 0)),  // at (0,0,0)
                    logicAlFrame,             // its logical volume
                    logicAlFrame->GetName(),  // its name
                    logicDetectorFrame,       // its mother  volume
                    false,                    // no boolean operation
                    copyNo,                   // copy number
                    fCheckOverlaps);          // checking overlaps

  return logicDetectorFrame;
}

void DetectorConstruction::DefineTank(G4LogicalVolume* logicWorld) {
  // Source
  G4double source_rad1 = 0.4 * cm, source_rad2 = 0.25 * cm, source_dz1 = 1 * cm,
           source_dz2 = 0.75 * cm;
  G4Tubs* SolidSourceCF =
      new G4Tubs("CfSource", 0, source_rad2, source_dz2 / 2, 0., twopi);
  G4Tubs* SolidSourceSt =
      new G4Tubs("CfCapsule", 0, source_rad1, source_dz1 / 2, 0., twopi);
  G4LogicalVolume* logicSourceCF =
      new G4LogicalVolume(SolidSourceCF,                       // its solid
                          fCf_src, SolidSourceCF->GetName());  // its name
  G4LogicalVolume* logicSourceSt =
      new G4LogicalVolume(SolidSourceSt, fSteel,
                          SolidSourceSt->GetName());  // its name

  visAttOthers = new G4VisAttributes(true, G4Colour(1., 0., 0., 1.));
  visAttOthers->SetForceSolid(true);
  logicSourceCF->SetVisAttributes(visAttOthers);
  new G4PVPlacement(0,                         // no rotation
                    G4ThreeVector(0, 0, 0.0),  // at (0,0,0)
                    logicSourceCF,             // its logical volume
                    logicSourceCF->GetName(),  // its name
                    logicSourceSt,             // its mother  volume
                    false,                     // no boolean operation
                    0,                         // copy number
                    fCheckOverlaps);           // checking overlaps
  visAttOthers = new G4VisAttributes(true, G4Colour(0., 1., 1., .3));
  visAttOthers->SetForceSolid(true);
  // Water
  G4double waterX = 100 * cm;
  G4double waterY = 100 * cm;
  G4double waterZ = 100 * cm;
  G4Box* solidWater = new G4Box("Water", 0.5 * waterX, 0.5 * waterY,
                                0.5 * waterZ);  // its size
  G4LogicalVolume* logicWater =
      new G4LogicalVolume(solidWater, fwater,
                          solidWater->GetName());  // its name
  logicWater->SetVisAttributes(visAttOthers);
  new G4PVPlacement(0,                                        // no rotation
                    G4ThreeVector(-49.55 * cm, 0, -10 * cm),  // at (0,0,0)
                    logicSourceSt,             // its logical volume
                    logicSourceSt->GetName(),  // its name
                    logicWater,                // its mother  volume
                    false,                     // no boolean operation
                    0,                         // copy number
                    fCheckOverlaps);           // checking overlaps

  visAttOthers = new G4VisAttributes(true, G4Colour(1., 1., 1., .2));
  visAttOthers->SetForceSolid(true);
  // FiberGlass1
  G4Box* solidFiberGlass1 =
      new G4Box("FiberGlass1", 0.5 * (102) * cm, 0.5 * (102) * cm,
                0.5 * (101) * cm);  // its size
  G4LogicalVolume* logicFiberGlass1 = new G4LogicalVolume(
      solidFiberGlass1, fPyrex, solidFiberGlass1->GetName());  // its name
  logicFiberGlass1->SetVisAttributes(visAttOthers);
  new G4PVPlacement(0,                              // no rotation
                    G4ThreeVector(0, 0, 0.5 * cm),  // at (0,0,0)
                    logicWater,                     // its logical volume
                    logicWater->GetName(),          // its name
                    logicFiberGlass1,               // its mother  volume
                    false,                          // no boolean operation
                    0,                              // copy number
                    fCheckOverlaps);                // checking overlaps

  // FiberAir1
  G4Box* solidAir1 = new G4Box("Air1", 0.5 * (130) * cm, 0.5 * (130) * cm,
                               0.5 * (130) * cm);  // its size
  G4LogicalVolume* logicAir1 =
      new G4LogicalVolume(solidAir1, fAir,
                          solidAir1->GetName());  // its name
  G4PhysicalVolumeStore::GetInstance()->Register(new G4PVPlacement(
      0,                                                // no rotation
      G4ThreeVector(0, 0, 0.5 * (102 - 130 - 1) * cm),  // at (0,0,0)
      logicFiberGlass1,                                 // its logical volume
      logicFiberGlass1->GetName(),                      // its name
      logicAir1,                                        // its mother  volume
      false,                                            // no boolean operation
      0,                                                // copy number
      fCheckOverlaps));                                 // checking overlaps
  // FiberGlass2

  G4Box* solidFiberGlass2 =
      new G4Box("FiberGlass2", 0.5 * 132 * cm, 0.5 * 132 * cm,
                0.5 * 131 * cm);  // its size
  G4LogicalVolume* logicFiberGlass2 = new G4LogicalVolume(
      solidFiberGlass2, fPyrex, solidFiberGlass2->GetName());  // its name
  logicFiberGlass2->SetVisAttributes(visAttOthers);
  new G4PVPlacement(0,                              // no rotation
                    G4ThreeVector(0, 0, 0.5 * cm),  // at (0,0,0)
                    logicAir1,                      // its logical volume
                    logicAir1->GetName(),           // its name
                    logicFiberGlass2,               // its mother  volume
                    false,                          // no boolean operation
                    0,                              // copy number
                    fCheckOverlaps);                // checking overlaps
  // FiberAir2
  G4Box* solidAir2 = new G4Box("Air2", 0.5 * 154 * cm, 0.5 * 154 * cm,
                               0.5 * 154 * cm);  // its size
  G4LogicalVolume* logicAir2 =
      new G4LogicalVolume(solidAir2, fAir,
                          solidAir2->GetName());  // its name
  new G4PVPlacement(
      0,                                                // no rotation
      G4ThreeVector(0, 0, 0.5 * (132 - 154 - 1) * cm),  // at (0,0,0)
      logicFiberGlass2,                                 // its logical volume
      logicFiberGlass2->GetName(),                      // its name
      logicAir2,                                        // its mother  volume
      false,                                            // no boolean operation
      0,                                                // copy number
      fCheckOverlaps);                                  // checking overlaps
  // Steel
  G4Box* solidSteel = new G4Box("Steel", 0.5 * 160 * cm, 0.5 * 160 * cm,
                                0.5 * 160 * cm);  // its size
  G4LogicalVolume* logicSteel =
      new G4LogicalVolume(solidSteel, fSteel,
                          solidSteel->GetName());  // its name
  visAttOthers = new G4VisAttributes(true, G4Colour(1., 1., 0., .2));
  visAttOthers->SetForceSolid(true);
  logicSteel->SetVisAttributes(visAttOthers);
  new G4PVPlacement(0,                       // no rotation
                    G4ThreeVector(0, 0, 0),  // at (0,0,0)
                    logicAir2,               // its logical volume
                    logicAir2->GetName(),    // its name
                    logicSteel,              // its mother  volume
                    false,                   // no boolean operation
                    0,                       // copy number
                    fCheckOverlaps);         // checking overlaps

  new G4PVPlacement(0,                                       // no rotation
                    G4ThreeVector(-12.5 * cm, 0, -95 * cm),  // at (0,0,0)
                    logicSteel,             // its logical volume
                    solidSteel->GetName(),  // its name
                    logicWorld,             // its mother  volume
                    false,                  // no boolean operation
                    0,                      // copy number
                    fCheckOverlaps);        // checking overlaps

  buldge_thickness = source_rad2;
  buldge_thickness = source_dz2 / 2;
  source_position = G4ThreeVector(-49.55 * cm, 0, -10 * cm) +
                    G4ThreeVector(-12.5 * cm, 0, -95 * cm) +
                    G4ThreeVector(0, 0, 0.5 * (132 - 154 - 1) * cm) +
                    G4ThreeVector(0, 0, 0.5 * cm) +
                    G4ThreeVector(0, 0, 0.5 * (102 - 130 - 1) * cm) +
                    G4ThreeVector(0, 0, 0.5 * cm);
}

void DetectorConstruction::DefineRoom(G4LogicalVolume* logicWorld) {
  // Side door
  G4double sideroomX = 6.0 * m;
  G4double sideroomY = 0.5 * m;
  G4double roomY = 3.5 * m;

  // Plaster
  G4double plasterY = 30 * cm;
  G4Box* solidPlaster = new G4Box("Plaster", 0.5 * sideroomX, 0.5 * plasterY,
                                  0.5 * roomY);  // its size
  G4LogicalVolume* logicPlaster = new G4LogicalVolume(solidPlaster, fConcrete,
                                                      "Plaster");  // its name
  logicPlaster->SetVisAttributes(visAttHide);
  new G4PVPlacement(0,                                 // no rotation
                    G4ThreeVector(0, (180) * cm, 0),   // at (0,0,0)
                    logicPlaster,                      // its logical volume
                    "Plaster1",                        // its name
                    logicWorld,                        // its mother  volume
                    false,                             // no boolean operation
                    0,                                 // copy number
                    fCheckOverlaps);                   // checking overlaps
  new G4PVPlacement(0,                                 // no rotation
                    G4ThreeVector(0, -(180) * cm, 0),  // at (0,0,0)
                    logicPlaster,                      // its logical volume
                    "Plaster2",                        // its name
                    logicWorld,                        // its mother  volume
                    false,                             // no boolean operation
                    0,                                 // copy number
                    fCheckOverlaps);                   // checking overlaps

  G4double frontroomX = 0.3 * m;
  G4double frontroomY = 3.3 * m;
  G4Box* solidOutside = new G4Box("OutRoom", 0.5 * frontroomX, 0.5 * frontroomY,
                                  0.5 * roomY);  // its size
  G4LogicalVolume* logicOutside =
      new G4LogicalVolume(solidOutside, fConcrete,
                          "OutsideRoom");  // its name
  // logicOutside->SetVisAttributes(visAttHide);
  new G4PVPlacement(0,                                   // no rotation
                    G4ThreeVector(227 * cm, 0 * cm, 0),  // at (0,0,0)
                    logicOutside,                        // its logical volume
                    "OutsideWall",                       // its name
                    logicWorld,                          // its mother  volume
                    false,                               // no boolean operation
                    0,                                   // copy number
                    fCheckOverlaps);                     // checking overlaps
  new G4PVPlacement(
      0,                                              // no rotation
      G4ThreeVector(-(127.5 + 72.5 / 2) * cm, 0, 0),  // at (0,0,0)
      logicOutside,                                   // its logical volume
      "Outside",                                      // its name
      logicWorld,                                     // its mother  volume
      false,                                          // no boolean operation
      0,                                              // copy number
      fCheckOverlaps);                                // checking overlaps

  // Ceiling & Floor
  G4double Ceilingz = 40 * cm;
  G4Box* solidCeiling =
      new G4Box("CeilingG", 0.5 * fExpHall_x, 0.5 * fExpHall_y,
                0.5 * Ceilingz);  // its size
  G4LogicalVolume* logicCeiling = new G4LogicalVolume(solidCeiling, fConcrete,
                                                      "CeilingLV");  // its name
  logicCeiling->SetVisAttributes(visAttHide);
  new G4PVPlacement(
      0,                                               // no rotation
      G4ThreeVector(0, 0, 175 * cm + 0.5 * Ceilingz),  // at (0,0,0)
      logicCeiling,                                    // its logical volume
      "Floor",                                         // its name
      logicWorld,                                      // its mother  volume
      false,                                           // no boolean operation
      0,                                               // copy number
      fCheckOverlaps);                                 // checking overlaps
  new G4PVPlacement(
      0,                                                // no rotation
      G4ThreeVector(0, 0, -175 * cm - 0.5 * Ceilingz),  // at (0,0,0)
      logicCeiling,                                     // its logical volume
      "Ceiling",                                        // its name
      logicWorld,                                       // its mother  volume
      false,                                            // no boolean operation
      0,                                                // copy number
      fCheckOverlaps);                                  // checking overlaps
}

