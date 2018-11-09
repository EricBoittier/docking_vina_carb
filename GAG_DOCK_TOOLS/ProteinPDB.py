import os
from LigandPDB import *
from FindCarbohydrateName import *


class Count:
    def __init__(self, start):
        self.count = start

    def forward(self):
        self.count += 1
        return self.count

    def getCount(self):
        return self.count


class ProteinPDB(PDB):
    '''
    A class to handle a .pdb file and extract information from it
    '''

    def __init__(self, filename):
        super().__init__(filename)
        self.resolution = "-"
        self.protein_name = "-"
        self.description = "-"
        self.pH = "-"
        self.heteroatoms = []
        self.DOI = ""
        self.ligands = list(nx.connected_component_subgraphs(self.graph))
        self.all_ligands_list = []
        temp = []
        for ligand in self.ligands:
            if len(ligand.nodes) > 1:
                ligand = list(ligand.nodes)
                catch = True
                for atomID in ligand:
                    try:
                        if not self.atoms[atomID].isHeteroAtom():
                            catch = False
                    except KeyError:
                        pass
                if catch:
                    temp.append(ligand)
        self.ligands = temp
        self.findInfo()
        self.known_residues = []
        self.residue_atoms = {}
        self.findKnownResidues()
        self.makeAllLigandsList()

    def makeAllLigandsList(self):
        for ligand in self.ligands:
            for id in ligand:
                self.all_ligands_list.append(id)

    def findKnownResidues(self):
        for line in self.heteroatoms:
            if not self.known_residues.__contains__(line):
                self.known_residues.append(line)
                self.residue_atoms[line.split()[1]] = line

    def findInfo(self):
        for line in self.lines:
            if line.startswith("COMPND   2 MOLECULE"):
                self.protein_name = line.split(" ", 5)[5].replace(
                    ";", "").strip("\n").capitalize()
            if line.startswith("REMARK 200  PH  "):
                self.pH = line.split(":")[-1]
                if self.pH.__contains__("NULL"):
                    self.pH = "-"
            if line.startswith("JRNL        DOI"):
                self.DOI = line.split()[2].strip("\n")
            if line.startswith("HEADER"):
                self.description = line.split(
                    " ", 1)[1][0:-40].strip(" ").capitalize()
            if line.startswith("REMARK   2 RESOLUTION."):
                self.resolution = line.split()[3]
            if line.startswith("HETATM"):
                if line.split()[3] != "MSE":
                    self.heteroatoms.append(line)

    def getInfo(self):
        return [str(self.filename.strip("PDBs/")), self.resolution, self.protein_name, self.description, self.DOI, self.pH]

    def getLigand(self, int):
        return self.ligands[int]

    def getLigands(self):
        return list(self.ligands)

    def getProteinName(self):
        return self.protein_name

    def getKnownResidues(self):
        return self.known_residues

    def getHeteroAtoms(self):
        return self.heteroatoms

    def saveProteinWithoutLigands(self, save_path):
        print(self.filename)
        with open("{}/{}_WithoutLigands.pdb".format(save_path, self.filename), "w") as file:
            for line in self.lines:
                try:
                    if line.split()[1] not in self.all_ligands_list:
                        file.write(line)
                except IndexError:
                    # This is ok, it just means the line is only one word long, e.g. END, so it's valid
                    file.write(line)

    def saveLigands(self, load_path):
        count = 1
        for ligand in self.ligands:
            print("saving {}/{}_LIGAND_{}.pdb".format(load_path, self.filename, count))
            with open("{}/{}_LIGAND_{}.pdb".format(load_path, self.filename, count), "w") as file:
                connections = []
                for id in ligand:
                    for line in self.lines:
                        try:
                            if line.split()[1] == id:
                                if line.split()[0] == "HETATM":
                                    file.write(line)
                                else:
                                    connections.append(line)
                        except IndexError:
                            # This is ok, it just means the line is only one word long, e.g. END, do it's valid
                            pass

                for connection in connections:
                    file.write(connection)

            count += 1






