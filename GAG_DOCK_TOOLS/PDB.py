from Atom import *
import networkx as nx
import os


class PDB(object):
    """docstring for LigandPDB"""

    def __init__(self, filename):
        super(PDB, self).__init__()
        self.filepath = filename
        self.filename = filename
        self.lines = []

        with open(filename, "r") as f:
            self.lines = f.readlines()

        self.filename = self.filename.split("/")[-1]
        self.filename = self.filename[0:-4]

        self.atoms = {}

        self.ring_names = []

        self.graph = nx.Graph()

        self.setAtomsAndConnections()

        self.connections = nx.to_dict_of_lists(self.graph)
        self.connections = {k: v for k,v in self.connections.items() if v}


        self.rings = nx.cycle_basis(self.graph)

        self.nameRings()

    def ringHash(self, ring):
        hash = 0
        for atom in ring:
            hash += int(atom) * 3
        return hash


    def getRings(self):
        return self.rings

    def getRingNames(self):
        return self.ring_names


    def getNeighbours(self, atomID):
        self.graph.__getitem__(atomID)

    def setAtomsAndConnections(self):
        for line in self.lines:
            if line.split()[0].__contains__("AT") and len(line.split()) > 10:
                atom = Atom(line)
                self.graph.add_node(atom.getID())
                self.atoms[line.split()[1]] = Atom(line)
            elif line.startswith("CONECT"):
                split = line.split()
                for connection in split[2:]:
                    self.graph.add_edge(split[1], connection)

    def getAllConnections(self):
        return self.connections

    def getConnections(self, atomID):
        return self.connections[atomID]

    def getFilename(self):
        return self.filename

    def getAtoms(self):
        return self.atoms

    def getLines(self):
        return self.lines

    def getGraph(self):
        return self.graph

    def nameRings(self):
        for ring in self.rings:
            try:
                self.ring_names.append(
                    (self.ringHash(ring), self.atoms[ring[0]].getLigandType()))
            except KeyError:
                pass

    def getFilepath(self):
        return self.filepath

    def getFilename(self):
        return self.filename

    def getRingName(self, ring):
        for ring_from_list, name in self.ring_names:
            if ring == ring_from_list:
                return name

    def getRingNamefromAtom(self, atom):
        return self.getRingName(self.getRing(atom))

    def getRing(self, atom):
        for ring in self.rings:
            if atom in ring:
                return self.ringHash(ring)

    def rename(self, new_name):
        with open(new_name, "w") as file:
            print("writing")
            for line in self.lines:
                file.write(line)
            #os.remove(self.filepath)