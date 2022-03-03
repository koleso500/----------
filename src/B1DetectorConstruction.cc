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
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "MakeLeaf.hh"






B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }



B1DetectorConstruction::~B1DetectorConstruction()
{ }



G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 200*cm, env_sizeZ = 200*cm;


  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //


  //
  // targetW
  //
  G4Material* W = nist->FindOrBuildMaterial("G4_W");
  G4ThreeVector pos1 = G4ThreeVector(0, 0, (-1000.+5./2.)*mm);

  G4Box* targetW =
    new G4Box("targetW",
    5./2*mm, 5./2*mm, 2./2.*mm);

  G4LogicalVolume* targetWL =
    new G4LogicalVolume(targetW,         //its solid
                        W,               //its material
                        "targetW");      //its name

  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    targetWL,             //its logical volume
                    "targetW",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //
  // targetCu
  //
  G4Material* Cu = nist->FindOrBuildMaterial("G4_Cu");
  G4ThreeVector pos2 = G4ThreeVector(0, 0, -1000*mm);

  G4Box* targetCu =
    new G4Box("targetCu",
    5./2*mm, 5./2*mm, 5./2*mm);

  G4LogicalVolume* targetCuL =
    new G4LogicalVolume(targetCu,         //its solid
                        Cu,          //its material
                        "targetCu");           //its name

  new G4PVPlacement(0,                       //no rotation
                    pos2,                    //at position
                    targetCuL,             //its logical volume
                    "targetCu",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


//upper part primary collimator
G4ThreeVector posPC = G4ThreeVector(0, 0, (-1000+6/2)*mm);
G4Tubs* PCUW = new G4Tubs("PCUW", 10.*mm, 40.*mm, 10.*mm, 0.*deg, 360.*deg);
G4LogicalVolume* PCUWL =
  new G4LogicalVolume(PCUW,         //its solid
                      W,          //its material
                      "PCUW");  // Trapezoid shape


new G4PVPlacement(0,                       //no rotation
                  posPC,                    //at position
                  PCUWL,             //its logical volume
                  "PCUW",                //its name
                  logicWorld,                //its mother  volume
                  false,                   //no boolean operation
                  0,                       //copy number
                  checkOverlaps);


//lower part primary collimator
posPC = G4ThreeVector(0, 0, (-1000+6/2+10+60/2)*mm);
G4Tubs* PCLW = new G4Tubs("PCLW", 0, 40.*mm, 60./2*mm, 0.*deg, 360.*deg);
G4LogicalVolume* PCLWL = new G4LogicalVolume(PCLW,         //its solid
                                        W,          //its material
                                        "PCLW");  // Trapezoid shape


new G4PVPlacement(0,                       //no rotation
                  posPC,                    //at position
                  PCLWL,             //its logical volume
                  "PCLW",                //its name
                  logicWorld,                //its mother  volume
                  false,                   //no boolean operation
                  0,                       //copy number
                  checkOverlaps);

G4Material* Vacuum = nist->FindOrBuildMaterial("G4_Galactic");
G4Cons* InPCLW = new G4Cons ("InPCLW", 0., 10*mm, 0., (76.*std::tan(14.*deg))*mm, 30.*mm,  0.*deg, 360.*deg);
G4LogicalVolume* InPCLWL = new G4LogicalVolume(InPCLW,         //its solid
                                        Vacuum,          //its material
                                        "InPCLW");


  new G4PVPlacement(0,                       //no rotation
                    posPC,                    //at position
                    InPCLWL,             //its logical volume
                    "InPCLW",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);


//vacuum window
//а черт его знает, пока без него

G4ThreeVector centre = G4ThreeVector (0,0, -1025*mm);

//ionization chamber

  G4Material* KAPTON = nist->FindOrBuildMaterial("G4_KAPTON");
    G4Tubs* ICTubeW = new G4Tubs("ionizationChamberTube", 0., 3.75*2.54*10.*mm, 0.005*25.4*mm, 0.*deg, 360.*deg);
    G4Tubs* ICTubeP = new G4Tubs("ionizationChamberTube", 0., 3.75*2.54*10.*mm, 0.002*25.4*mm, 0.*deg, 360.*deg);


    // W1
    centre= G4ThreeVector(0.,0.,(-1090+15+148.35)*mm);
    G4LogicalVolume *PCUTubeW1LV = new G4LogicalVolume(ICTubeW, KAPTON, "ionizationChamberTubeW1LV");
    new G4PVPlacement(0, centre, PCUTubeW1LV,"ionizationChamberTubeW1PV", logicWorld, false, 0, checkOverlaps);
    // P1
    centre= G4ThreeVector(0.,0.,(-1090+15+150.73)*mm);
    G4LogicalVolume *PCUTubeP1LV = new G4LogicalVolume(ICTubeP, KAPTON, "ionizationChamberTubeP1LV");
    new G4PVPlacement(0, centre, PCUTubeP1LV, "ionizationChamberTubeP1PV", logicWorld, false, 0, checkOverlaps);
    // W2
    centre= G4ThreeVector(0.,0.,(-1090+15+155.5)*mm);
    G4LogicalVolume *PCUTubeW2LV = new G4LogicalVolume(ICTubeW, KAPTON, "ionizationChamberTubeW2LV");
    new G4PVPlacement(0, centre, PCUTubeW2LV, "ionizationChamberTubeW2PV", logicWorld, false, 0, checkOverlaps);

    // P2
    centre= G4ThreeVector(0.,0.,(-1090+15+153.12)*mm);
    G4LogicalVolume *PCUTubeP2LV = new G4LogicalVolume(ICTubeP, KAPTON, "ionizationChamberTubeP2LV");
    new G4PVPlacement(0, centre, PCUTubeP2LV, "ionizationChamberTubeP2PV", logicWorld, false, 0, checkOverlaps);

    // W3
    centre= G4ThreeVector(0.,0.,(-1090+15+162.65)*mm);
    G4LogicalVolume *PCUTubeW3LV = new G4LogicalVolume(ICTubeW, KAPTON, "ionizationChamberTubeW3LV");
    new G4PVPlacement(0, centre, PCUTubeW3LV, "ionizationChamberTubeW3PV", logicWorld, false, 0, checkOverlaps);

    // P3
    centre= G4ThreeVector(0.,0.,(-1090+15+157.88)*mm);
    G4LogicalVolume *PCUTubeP3LV = new G4LogicalVolume(ICTubeP, KAPTON, "ionizationChamberTubeP3LV");
    new G4PVPlacement(0, centre, PCUTubeP3LV, "ionizationChamberTubeP3PV", logicWorld, false, 0, checkOverlaps);

    // P4
    centre= G4ThreeVector(0.,0.,(-1090+15+160.27)*mm);
    G4LogicalVolume *PCUTubeP4LV = new G4LogicalVolume(ICTubeP, KAPTON, "ionizationChamberTubeP4LV");
    new G4PVPlacement(0, centre, PCUTubeP4LV, "ionizationChamberTubeP4PV", logicWorld, false, 0, checkOverlaps);






//mlc
//centre= G4ThreeVector(0.,0.*mm,(-880)*mm);
//G4Box* box1 = new G4Box("jaw1box", 90.*mm, 90.*mm, 78./2.*mm);
//G4LogicalVolume *boxl1 = new G4LogicalVolume (box1, W, "jaw1boxl");
//new G4PVPlacement (0, centre, boxl1, "jaw1boxp", logicWorld, false,0, checkOverlaps);


//centre= G4ThreeVector(0.,0.*mm,(-790)*mm);
//G4Box* box2 = new G4Box("jaw1box", 90.*mm, 90.*mm, 78./2.*mm);
//G4LogicalVolume *boxl2 = new G4LogicalVolume (box2, W, "jaw1boxl");
//new G4PVPlacement (0, centre, boxl2, "jaw1boxp", logicWorld, false,0, checkOverlaps);
    




//phantom
G4double A, Z;
A = 1.01*g/mole;
G4Element* elH = new G4Element ("Hydrogen","H",Z = 1.,A);


A = 16.00*g/mole;
G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);

G4double d= 1.00*g/cm3; //density
G4int natoms, ncomponents;
G4Material* PMMA = new G4Material("Water",d,ncomponents=2);
PMMA->AddElement(elH, natoms=2);
PMMA->AddElement(elO, natoms=1);

G4Box* phantom = new G4Box ("phantom", 20.*cm, 20.*cm, 20.*cm);
G4LogicalVolume *phantoml = new G4LogicalVolume (phantom, PMMA, "phantoml");
new G4PVPlacement (0, G4ThreeVector(0,0,0), phantoml, "phantomp", logicWorld, false, 0, checkOverlaps);

//MakeLeaf(G4ThreeVector(0 * mm, 0 * mm, -870 * mm), logicWorld);    
//MakeLeaf(G4ThreeVector(0 * mm, 0 * mm, -780 * mm), logicWorld); 
//что-то
    for (int i=1;i<114;i++){
        //G4Threevector = G4ThreeVector(0, i, 0);
        MakeLeaf(G4ThreeVector(0, 0,-685 -i*2 * CLHEP::mm), logicWorld);
    }
//  fScoringVolume = targetCuL;

  //
  //always return the physical World
  //
//G4PVPlacement* MakeLeaf(G4ThreeVector&& pos, G4LogicalVolume* mother_vol)
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

