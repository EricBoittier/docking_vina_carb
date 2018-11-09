class Atom(object):
    """
    An Object that holds the data of a line for a pdb
    file.

    Example of a line:

    HETATM 3891  O3  GC4 A 516  24.953  25.336  -2.509  1.00 20.48  O

    """

    def __init__(self, line):
        super(Atom, self).__init__()
        self.line = line
        self.kind = line[0:5].strip()
        self.id = line[6:10].strip()
        self.atomtype = line[12:15].strip()
        self.ligandtype = line[17:20].strip()
        self.chain = line[21].strip()
        self.ligandID = line[22:25].strip()
        self.x = float(line[30:37].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:53].strip())
        self.connections = []
        self.atomname = line[76:77].strip()

    def getID(self):
        return self.id

    def getLigandType(self):
        return self.ligandtype

    def getLigandID(self):
        return self.ligandID

    def getAtomType(self):
        return self.atomtype

    def getAtomename(self):
        return self.atomname

    def getXYZ(self):
        return [self.x, self.y, self.z]

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def getZ(self):
        return self.z

    def isAtom(self, atomtype):
        return self.atomtype.__contains__(atomtype)

    def getKind(self):
        return self.kind

    def isHeteroAtom(self):
        return self.kind == "HETATM"
