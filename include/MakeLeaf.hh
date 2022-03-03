#ifndef __MAKE_LEAF_HH__
#define __MAKE_LEAF_HH__

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

G4PVPlacement* MakeLeaf(G4ThreeVector&& pos, G4LogicalVolume* mother_vol)
{
    G4NistManager* nist = G4NistManager::Instance();
    auto wolfram = nist->FindOrBuildMaterial("G4_W");
    auto air = nist->FindOrBuildMaterial("G4_AIR");

    G4Box* leaf_solid = new G4Box("Leaf", 9 * CLHEP::cm, 9 * CLHEP::cm, 0.1 * CLHEP::cm);
    G4LogicalVolume* leaf_logic = new G4LogicalVolume(leaf_solid, wolfram, "Leaf");
    G4PVPlacement* leaf_physical = new G4PVPlacement(
      0,
      pos,
      leaf_logic,
      "Leaf",
      mother_vol,
      0,
      0,
      false);

    G4Box* hole_solid = new G4Box("Hole", 2 * CLHEP::cm, 2 * CLHEP::cm, 0.1 * CLHEP::cm);
    G4LogicalVolume* hole_logic = new G4LogicalVolume(hole_solid, air, "Hole");
    G4PVPlacement* hole_physical = new G4PVPlacement(
      0,
      pos,
      hole_logic,
      "Leaf",
      mother_vol,
      0,
      0,
      false);

    return leaf_physical;
}

#endif
