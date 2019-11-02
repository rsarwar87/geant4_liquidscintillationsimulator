
#ifndef OpNoviceDetectorConstruction_h
#define OpNoviceDetectorConstruction_h 1

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
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4Box.hh"
#include "G4Cache.hh"
#include "G4Cons.hh"
#include "G4IntersectionSolid.hh"
#include "G4Isotope.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VisAttributes.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

 public:
  virtual G4VPhysicalVolume* Construct();

  G4VPhysicalVolume* GetWorld() { return fPBox; };

  G4double GetSize() { return fExpHall_x * fExpHall_y * fExpHall_z; };
  G4Material* GetMaterial() { return fwater; };
  static G4int nb_cryst;
  static G4int hist[10][100];
  static G4double ring_R1, ring_W1;
  static bool lancs;
  static G4ThreeVector source_position;

 private:
  void DefineMaterials();
  G4LogicalVolume* DefineEJ309(bool upright, int copyNo);
  void DefineORNL(G4LogicalVolume* logicWorld, int nd = 15);
  G4LogicalVolume* DefineHE3(bool upright, int copyNo);
  void DefineTank(G4LogicalVolume* logicWorld);
  void DefineRoom(G4LogicalVolume* logicWorld);
  void DefineChamber(G4LogicalVolume* logicWorld);
  void SurfaceProperties(G4LogicalVolume* fHousing_log2,
                         G4LogicalVolume* fPhotocath_log);
  void DefineLANCS(G4LogicalVolume* expHall_log);

  G4VPhysicalVolume* fPBox;
  G4double fExpHall_x;
  G4double fExpHall_y;
  G4double fExpHall_z;

  G4double fBubble_x;
  G4double fBubble_y;
  G4double fBubble_z;

  G4bool fCheckOverlaps;

  G4double fTank_x;
  G4double fTank_y;
  G4double fTank_z;

  // Materials & Elements
  G4Isotope* fiN;
  G4Isotope* fiO;
  G4Isotope* fiH;
  G4Isotope* fiAm;
  G4Isotope* fiLi;
  G4Isotope* fiC;
  G4Isotope* fiCa;
  G4Isotope* fiS;
  G4Isotope* fiPu238;
  G4Isotope* fiPu239;
  G4Isotope* fiPu240;
  G4Isotope* fiPu241;
  G4Isotope* fiPu242;
  G4Isotope* fiCm242;
  G4Isotope* fiCm244;
  G4Isotope* fiCm248;
  G4Isotope* fiCm246;
  G4Isotope* fiCf250;
  G4Isotope* fiCf251;
  G4Isotope* fiCf249;
  G4Isotope* fiCf252;
  G4Isotope* fiCf254;
  G4Isotope* fiCf256;
  G4Isotope* fiU233;
  G4Isotope* fiU235;
  G4Isotope* fiU236;
  G4Isotope* fiU238;

  double U235p = 20.1;

  G4Element* fN;
  G4Element* fO;
  G4Element* fCf;
  G4Element* fAm;
  G4Element* fLi;
  G4Material* fCd;
  G4Material* fHe3;
  G4Element* fC;
  G4Element* fCa;
  G4Element* fU;
  G4Element* fS;
  G4Element* fPu;
  G4Element* fCm;
  G4Element* fH;
  G4Material* fAl;
  G4Material* fUOx;
  G4Material* fPb;
  G4Material* fAir;
  G4Material* fVacuum;
  G4Material* fAmLi;
  G4Material* fGlass;
  G4Material* fPstyrene;
  G4Material* fPMMA;
  G4Material* fPethylene1;
  G4Material* fPethylene2;
  G4Material* fConcrete;
  G4Material* fSteel;
  G4Material* fwater;
  G4Material* fPyrex;
  G4Material* fplaster;
  G4Material* fCf_src;
  G4Material* fPu_src;
  G4Material* fU_src;
  G4Material* fWood;
  G4Material* fmu_metal;
  G4Material* fscintillator;
  G4Material* fLead;
  G4Material* cryst_mat;

  G4MaterialPropertiesTable* fScintillator_mt;
  G4MaterialPropertiesTable* fcryst_mat;
  G4MaterialPropertiesTable* ffiber_mat;
  G4double source_radius;
  G4double source_height;
  G4double buldge_thickness, inner_linning_thickness, outter_linning_thickness;
  G4double outter_steel;
  G4double polly_thick[4];
  G4double det_place[6];
  G4double lsd_place[6];

  G4VisAttributes* visAttHide;
  G4VisAttributes* visAttOthers;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*OpNoviceDetectorConstruction_h*/
