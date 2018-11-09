import os


class Config_editor(object):
    """docstring for Config_editor"""

    def __init__(self, path, filename, coeff, cut_off):
        super(Config_editor, self).__init__()
        self.path, self.filename = path, filename
        self.filepath = path + "/" + filename
        self.file = open(self.filepath)
        self.filelines = self.file.readlines()
        self.chi_co = coeff
        self.chi_cut = cut_off
        self.makeNewConfig()

    def makeNewConfig(self):
        newfile = open(self.path + "{}_{}_".format(self.chi_co,
                                                   self.chi_cut) + self.filename, "w")
        newfile.write("chi_coeff = {}\n".format(self.chi_co))
        newfile.write("chi_cutoff = {}\n".format(self.chi_cut))
        for line in self.filelines:
            if not line.__contains__("chi_"):
                newfile.write(line)

path = "/Users/ericboittier/Desktop/docking_rigid/"
for dir in os.listdir(path):
    if dir.__contains__("LIGAND"):
        for file in os.listdir(path + dir):
            if file.__contains__("conf") and len(file.split("_")) == 3:
                Config_editor(
                    path + dir + "/", file, 0, 12)
