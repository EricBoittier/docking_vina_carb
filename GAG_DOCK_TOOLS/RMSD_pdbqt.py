import numpy as np
import os
import re


"""
Script to calculate the RMSD of two identical ligands in .pdbqt format after docking with Vina-Carb.
"""
class RMSD_pdbqt(object):
    """

    """
    def __init__(self, filename):
        super(RMSD_pdbqt, self).__init__()
        self.file = open(filename, "r")
        self.file_lines = self.file.readlines()
        self.atoms = {}
        self.setAtoms()

    def setAtoms(self):
        for line in self.file_lines:
            if line.startswith("ATOM"):
                split = line.split()
                self.atoms[split[1]] = [split[2], float(split[6]), float(split[7]), float(split[8])]

    def getAtoms(self):
        return self.atoms


basic = ["C1", "C2", "C3", "C4", "C5", "O5"]
sulfur = ["S", "S2", "S3", "S4"]


def calculateRMSD(ligand1, ligand2, name):
    """
    Calculates RMSD of two ligands in pdbqt format...

    :param ligand1: RMSD_pdbqt(object)
    :param ligand2: RMSD_pdbqt(object)
    :return: RMSD of ligands 1 and 2 (float)
    """
    ligand1atoms = ligand1.getAtoms()
    ligand2atoms = ligand2.getAtoms()
    sumOfDifSqr = 0

    for id in range(len(ligand1atoms)):
        if len(ligand2atoms) > 0 and ligand1atoms[str(id+1)][0] in basic:
            try:
                for value in difSqr(ligand2atoms[str(id+1)][1:], ligand1atoms[str(id+1)][1:]):
                    sumOfDifSqr += value
            except KeyError:
                print("Key Error: {}".format(name))

    return np.sqrt(sumOfDifSqr)

def difSqr(predicted, actual):
    predicted = np.array(predicted)
    actual = np.array(actual)
    x = predicted - actual
    return x * x

