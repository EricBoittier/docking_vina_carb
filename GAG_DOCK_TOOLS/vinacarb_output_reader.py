import os

class PDBQTOutput:
    def __init__(self, result, bonds, atoms):
        self.result = result
        self.bonds = bonds
        self.atoms = atoms

    def get_atoms(self):
        return self.atoms

    def get_bonds(self):
        return self.bonds

    def get_result(self):
        return self.result


class File:
    """
    A class to handle files
    """
    def __init__(self, filename):
        self.pdbqt_files = {}
        self.filename = filename
        self.file_lines = []
        self.read_file()
        self.make_pqbqt_files()

    def filename_print(self):
        print(self.filename)

    def set_filename(self, string):
        self.filename = string

    def get_filename(self):
        return self.filename

    def get_file(self):
        return self.file_lines

    def read_file(self):
        with open(self.filename, 'r') as f:
            self.file_lines = f.readlines()

    def make_pqbqt_files(self):
        root_index = []
        torsion_index = []
        endroot_index = []
        result_index = []
        for x in range(len(self.file_lines)):
            if self.file_lines[x].__contains__("ENDROOT"):
                endroot_index.append(x)
            elif self.file_lines[x].__contains__("ROOT"):
                root_index.append(x)
            elif self.file_lines[x].__contains__("active torsions"):
                torsion_index.append(x)
            elif self.file_lines[x].__contains__("VINA RESULT"):
                result_index.append(x)

        count = 0
        while count < len(result_index):
            self.pdbqt_files[count] = PDBQTOutput(self.file_lines[result_index[count]],
                                                  self.file_lines[torsion_index[count]:root_index[count]],
                                                  self.file_lines[root_index[count]:endroot_index[count]])
            count+=1

        count = 0
        while count < len(result_index):
            with open(self.filename[0:-6]+str(count)+".pdbqt", "w") as f:
                # f.write("#")
                # for z in self.pdbqt_files[count].get_result():
                #     f.write(z)
                for x in self.pdbqt_files[count].get_bonds():
                    f.write(x)
                for y in self.pdbqt_files[count].get_atoms():
                    f.write(y)
                f.write("ENDROOT \n")
                f.write("TORSDOF "+str(len(self.pdbqt_files[count].get_bonds())-2)+"\n")
            count+=1


    def get_pdbqt_files(self):
        return self.pdbqt_files
#
# if __name__ == '__main__':
#     import sys
#     filename = sys.argv[1]
#     file = File(filename)
#     print(file.get_pdbqt_files())

p = "/Users/ericboittier/Desktop/GAGTOOLS/docking_data_processing/docking_full_run_1"
for y in os.listdir(p):
    if y.__contains__("LIGAND"):
        for x in os.listdir(p+"/"+y):
            if x.__contains__("_out."):
                f = File(p+"/"+y+"/"+x)
                f.get_pdbqt_files()