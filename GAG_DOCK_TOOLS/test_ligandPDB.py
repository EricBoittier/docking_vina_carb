from LigandPDB import *
from ProteinPDB import *
import imolecule
import os
path = "/Users/ericboittier/Desktop/docking/"
pdbs = os.listdir(path+"protein-ligand-complexes")
ligands = os.listdir(path+"ligands")

def drawLigand(filename):
    print(imolecule.draw(filename, size=(500, 500)))

def isAxialorEq(b):
    if b:
        print("Is axial")
    else:
        print("Is equitorial")



for ligand in ligands:
    if ligand.__contains__("LIGAND"):
        print("")
        print(ligand)
        ligand = LigandPDB("/Users/ericboittier/Desktop/docking/ligands/"+ligand)
        print(ligand.getPsi())


        # for atom, label in ligand.atom_labels.items():
        #     if label == "C1":
        #         try:
        #             isAxialorEq(ligand.isAxial(atom))
        #         except:
        #             pass



