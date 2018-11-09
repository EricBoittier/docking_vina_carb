import os

class editPDBQT(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.filename = filepath.split("/")[-1]
        self.path = ""
        for path in filepath.split("/")[:-1]:
            self.path += path + "/"

        self.file = open(filepath, "r")
        self.lines = self.file.readlines()


    def make_flexible_glycosidic_linkages(self):
        newfile = open(self.path+"flexible_phi_psi_"+self.filename, "w")
        for line in self.lines:
            if line.__contains__("between atoms") and self.is_glycosidic_linkage(line):
                s = line[0:12]+"A"+line[14:]
                newfile.write(s)
            if line.__contains__("between atoms") and not self.is_glycosidic_linkage(line):
                s = line[0:12]+"I"+line[14:]
                newfile.write(s)
            else:
                newfile.write(line)


    def is_glycosidic_linkage(self, line):
        split = line.split()
        atom1 = split[-1].split("_")[0]
        atom2 = split[-3].split("_")[0]

        if len(atom1) == 2 and len(atom2) == 2:
            if atom1[0] == "C" and atom2[0] == "O" or atom1[0] == "C" and atom2[0] == "O":
                if atom1 == "O1" or atom2 == "O1":
                    return True
                elif atom1[1] != atom2[1]:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

editPDBQT("/Volumes/Eric/Desktop/docking/ligands/pdbqts/1XMN_LIGAND_2_withHs.pdb.pdbqt").make_flexible_glycosidic_linkages()

path_to_run = "/Users/ericboittier/Desktop/docking_november/run/"

for directory in os.listdir(path_to_run):
    if directory.__contains__("LIGAND"):
        for file in os.listdir(path_to_run+directory):
            if file.__contains__(".mol2") and not file.__contains__("flexible"):
                print(file)
                editPDBQT(path_to_run+directory+"/"+file).make_flexible_glycosidic_linkages()
                os.remove(path_to_run+directory+"/"+file)





