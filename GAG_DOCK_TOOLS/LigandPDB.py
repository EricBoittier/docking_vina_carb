from PDB import *
import os
from abbreviationToSugarName import *
import networkx as nx
import matplotlib.pyplot as plt
from XYZdihedral import dihedral
from mpl_toolkits.mplot3d import Axes3D
import numpy as np



class LigandPDB(PDB):
    """docstring for LigandPDB"""

    def __init__(self, filename):
        super().__init__(filename)
        self.filename = filename
        self.atom_labels = {}
        self.name = ""
        self.linkages = {}
        self.atom_positions = {}
        self.phi_angles = []
        self.psi_angles = []
        self.addLabels()
        self.addPositions()
        self.findPsi()
        self.findPhi()
        self.addLinkages()
        self.sugarGraph = nx.DiGraph()

        self.ring_names = []
        self.renameRings()

        self.findName()
        self.atomsToCompare = {}
        self.findAtomsToCompare()



        self.xMin = self.yMin = self.zMin = 10000000000
        self.zMax = self.xMax = self.yMax = -10000000000
        self.centreX = self.centreY = self.centreZ = self.sizeX = self.sizeY = self.sizeZ = 0
        self.findADbox()

    def getFilename(self):
        return self.filename

    def isAxialorEq(self, b):
        if b:
            return "ax"
        else:
            return "eq"

    def getAtomLabels(self):
        return self.atom_labels

    def renameRings(self):
        for ring in self.rings:
            anomer = ""
            try:
                for atom in ring:
                    if self.getAtomLabels()[atom] == "C1":
                        try:
                            anomer = self.isAxialorEq(self.isAxial(atom))
                        except:
                            pass
                self.ring_names.append(
                    (self.ringHash(ring), self.atoms[ring[0]].getLigandType(), anomer))
            except KeyError:
                pass

    def isAxial(self, atom):
        test_label = atom
        test = self.atoms[atom].getXYZ()
        neighbour_not_in_ring = ""
        neighbour_not_in_ring_label = ""

        in_ring = ""
        in_ring_label = ""

        in_rings_neighbour = ""
        in_rings_neighbour_label = ""

        try:

            for neighbour in self.connections[atom]:
                if not self.isAtomInRing(neighbour):
                    neighbour_not_in_ring = self.atoms[neighbour].getXYZ()
                    neighbour_not_in_ring_label = neighbour
                elif self.isAtomInRing(neighbour):
                    in_ring_label = neighbour
                    in_ring = self.atoms[neighbour].getXYZ()

            for neighbour in self.connections[in_ring_label]:
                if self.isAtomInRing(neighbour) and neighbour != atom:
                    in_rings_neighbour = self.atoms[neighbour].getXYZ()
                    in_rings_neighbour_label = neighbour
                    break

        except KeyError:
            pass

        return abs(dihedral(neighbour_not_in_ring, test, in_ring, in_rings_neighbour)) < 90



    def findAtomsToCompare(self):
        """

        """
        for ring in self.rings:
            ringhash = self.ringHash(ring)
            self.atomsToCompare[ringhash] = []
            for atom in ring:
                self.atomsToCompare[ringhash].append(self.atoms[atom])
                """
                Breadth first search

                """
                for start in self.connections[atom]:
                    if start not in ring[0]:
                        explored = []
                        queue = [start]
                        levels = {start: 0}
                        visited = [ring[0], start]
                        while queue:
                            node = queue.pop(0)
                            explored.append(node)
                            neighbours = self.connections[node]

                            for neighbour in neighbours:
                                if neighbour not in visited:
                                    queue.append(neighbour)
                                    visited.append(neighbour)
                                    levels[neighbour] = levels[node] + 1

                        for once_visited in explored:
                            if self.atom_labels[once_visited][0] == "S" or self.atom_labels[once_visited][0] == "N":
                                if self.atoms[once_visited] not in self.atomsToCompare[ringhash]:
                                    self.atomsToCompare[ringhash].append(self.atoms[once_visited])


    def getAtomsToCompare(self):
        return self.atomsToCompare


    def findPsi(self):
        """

        C1 - O - Cx - Cx-1

        Psi is defined as the dihedral angle between C1, the oxygen of the
        glycosidic bond, CX (where x is the ring carbon number connected to the
        glycosidic oxygen), and CX-1 is the carbon numbered 1 less than the CX.

        :return:
        """
        c1 = []
        glycO = []
        cx = []
        cxtake1 = []
        linkage = ""
        C1 = ""
        CX = ""

        psi_angles = []

        for posC1 in self.connections:
            if self.isC1(posC1):
                c1 = self.XYZ(posC1)
                C1 = posC1
                for pos_glycO in self.connections[posC1]:
                    if self.isGylcosidic(pos_glycO):
                        glycO = self.XYZ(pos_glycO)
                        """
                        Now the idea is, if you find CX, find Carbon X-1
                        """
                        for poscx in self.connections[pos_glycO]:
                            if poscx != posC1:
                                cx = self.XYZ(poscx)
                                CX = poscx

                                """
                                Cases for (C4, C3), (C3, C2), (C2, C1)
                                """
                                #C4, C3
                                if self.isC4(poscx):
                                    for connection in self.connections[poscx]:
                                        if self.isC3(connection):
                                            posCx = connection
                                            cxtake1 = self.XYZ(connection)
                                            cx = self.XYZ(poscx)
                                            linkage = "1-4"

                                #C3, C2
                                if self.isC3(poscx):
                                    for connection in self.connections[poscx]:
                                        if self.isC2(connection):
                                            posCx = connection
                                            cxtake1 = self.XYZ(connection)
                                            linkage = "1-3"

                                #C2, C1
                                if self.isC2(poscx):
                                    for connection in self.connections[poscx]:
                                        if self.isC1(connection):
                                            posCx = connection
                                            self.cxtake1 = self.XYZ(connection)
                                            linkage = "1-2"

                                if cxtake1 != []:
                                    name = "{}-({})-{}".format\
                                        (abbreviationToSugarName(self.getRingNamefromAtom(C1)), linkage,
                                         abbreviationToSugarName(self.getRingNamefromAtom(CX)))

                                    anomer1 = self.isAxialorEq(self.isAxial(C1))
                                    anomer2 = self.isAxialorEq(self.isAxial(poscx))
                                    linkage = "{}{}-{}{}".format(anomer1, linkage.split("-")[0], anomer2, linkage.split("-")[1])
                                    psi_angles.append(
                                        [linkage, pos_glycO, dihedral(c1, glycO, cx, cxtake1), [C1, CX], name])

        self.psi_angles = psi_angles

    def findPhi(self):
        """
        O5 - C1 - O - Cx

        Psi is defined as the dihedral angle between the ring oxygen, carbon 1,
        the oxygen of the glycosidic bond,CX (where x is the ring carbon number
        connected to the glycosidic oxygen).

        :return:
        """
        o5 = []
        c1 = []
        glycO = []
        cx = []
        linkage = ""

        phi_angles = []

        O5 = ""
        CX = ""

        for posO5 in self.connections:
            if self.isO5(posO5):
                o5 = self.XYZ(posO5)
                O5 = posO5
                for pos_c1 in self.connections[posO5]:
                    if self.isC1(pos_c1):
                        c1 = self.XYZ(pos_c1)
                        for pos_glycO in self.connections[pos_c1]:
                            if self.isGylcosidic(pos_glycO):
                                glycO = self.XYZ(pos_glycO)
                                for pos_cx in self.connections[pos_glycO]:
                                    if pos_cx != pos_c1:
                                        if self.isC1(pos_cx):
                                            linkage = "1-1"
                                            cx = self.XYZ(pos_cx)
                                            CX = pos_cx
                                        elif self.isC2(pos_cx):
                                            linkage = "1-2"
                                            cx = self.XYZ(pos_cx)
                                            CX = pos_cx
                                        elif self.isC3(pos_cx):
                                            linkage = "1-3"
                                            cx = self.XYZ(pos_cx)
                                            CX = pos_cx
                                        elif self.isC4(pos_cx):
                                            linkage = "1-4"
                                            cx = self.XYZ(pos_cx)
                                            CX = pos_cx

                                        name = "{}-({})-{}".format\
                                            (abbreviationToSugarName(self.getRingNamefromAtom(O5)), linkage,
                                             abbreviationToSugarName(self.getRingNamefromAtom(CX)))

                                        anomer1 = self.isAxialorEq(self.isAxial(pos_c1))
                                        anomer2 = self.isAxialorEq(self.isAxial(pos_cx))
                                        linkage = "{}{}-{}{}".format(anomer1, linkage.split("-")[0], anomer2, linkage.split("-")[1])

                                        phi_angles.append(
                                            [linkage, pos_glycO, dihedral(o5, c1, glycO, cx), [O5, CX], catchNameExceptions(name)])
        self.phi_angles = phi_angles

    def XYZ(self, atomID):
        return self.atoms[atomID].getXYZ()

    def getPsi(self):
        return self.psi_angles

    def getPhi(self):
        return self.phi_angles

####################################################################

        # Atom Label Logic

    def isConnected(self, atomID1, atomID2):
        return self.graph.has_edge(atomID1, atomID2)

    def isC1(self, atomID):
        is_carbon5 = self.isC5(atomID)
        inRing = self.isAtomInRing(atomID)

        if not is_carbon5 and inRing:
            for neighbours in self.connections[atomID]:
                if self.isO5(neighbours):
                    return True

        try:
            if self.atom_labels[atomID] == "C1":
                return True
        except KeyError:
            return False

    def isC2(self, atomID):
        if self.atoms[atomID].isAtom("C") and self.isAtomInRing(atomID):
            for connections in self.connections[atomID]:
                if self.isC1(connections):
                    return True
        return False

    def isC3(self, atomID):
        if self.atoms[atomID].getAtomType()[0] == "C" and self.isAtomInRing(atomID):
            for connections in self.connections[atomID]:
                if self.isC4(connections) and not self.isC5(atomID):
                    return True
        return False

    def isC4(self, atomID):
        for neighbours in self.getConnections(atomID):
            if self.isC5(neighbours) and self.atoms[atomID].isAtom("C") and self.isAtomInRing(atomID):
                return True
        try:
            if self.atom_labels[atomID] == "C4":
                return True
        except KeyError:
            return False

    def isO5(self, atomID):
        if len(self.getConnections(atomID))\
                == 2 and self.atoms[atomID].isAtom("O") \
                and self.isAtomInRing(atomID):
            return True
        else:
            return False

    def isGylcosidic(self, atomID):
        if len(self.connections[atomID]) == 2 and self.atoms[atomID].isAtom("O") and not self.isAtomInRing(atomID):
            if self.isAtomInRing(self.connections[atomID][0]) and self.isAtomInRing(self.connections[atomID][1]):
                return True
        return False

    def isC5(self, atomID):
        if self.atoms[atomID].isAtom("C") and self.isAtomInRing(atomID):
            for connections in self.connections[atomID]:
                if self.atoms[connections].isAtom("C") and not self.isAtomInRing(connections):
                    return True
        try:
            if self.atom_labels[atomID] == "C5":
                return True
        except KeyError:
            return False

    def isAtomInRing(self, atomID):
        for ring in self.rings:
            for atom in ring:
                if atom == atomID:
                    return True
        else:
            return False

    def getLabel(self, atomID):
        if self.isC1(atomID):
            return "C1"
        elif self.isC2(atomID):
            return "C2"
        elif self.isC3(atomID):
            return "C3"
        elif self.isC4(atomID):
            return "C4"
        elif self.isC5(atomID):
            return "C5"
        elif self.isO5(atomID):
            return "O5"
        elif self.isGylcosidic(atomID):
            return "Link"
        else:
            return self.atoms[atomID].getAtomType()

    def getLabels(self):
        return self.atom_labels

    def addLabels(self):
        for atomIDs in self.atoms:
            self.atom_labels[atomIDs] = self.getLabel(atomIDs)

    def addPosition(self, atomID, xyz):
        try:
            self.atom_positions[atomID] = [xyz, "{}".format(self.getLabel(atomID))]
        except KeyError:
            pass

    def addPositions(self):
        for atomIDs in self.atoms:
            self.addPosition(atomIDs, self.atoms[atomIDs].getXYZ())

    def getPositions(self):
        return self.atom_positions

    def addLinkages(self):
        for entry in self.psi_angles:
            self.linkages[self.getRing(entry[3][0])] = entry[0]

    def getName(self):
        return self.name

    def findName(self):
        for linkage in self.psi_angles:
            self.sugarGraph.add_node(self.getRing(linkage[3][0]))
            self.sugarGraph.add_node(self.getRing(linkage[3][1]))
            self.sugarGraph.add_edge(self.getRing(linkage[3][0]), self.getRing(linkage[3][1]))

        name = ""
        start = ""
        end = ""
        for node in self.sugarGraph.nodes:
            if len(list(self.sugarGraph.successors(node))) == 0:
                end = node
            elif len(list(self.sugarGraph.predecessors(node))) == 0:
                start = node

        shortestpath = nx.shortest_path(self.sugarGraph, start, end)

        for c in shortestpath:
            n = self.getAnomericName(c)
            name += abbreviationToSugarName(n[0])
            try:
                name += "-({}{})-".format(n[1], self.linkages[c])
            except KeyError:
                pass

        self.name = catchNameExceptions(name)

    def getAnomericName(self, ring):
        for ring_from_list, name, a in self.ring_names:
            if ring == ring_from_list:
                return [name, a]

    def getSNFGname(self):
        for linkage in self.psi_angles:
            self.sugarGraph.add_node(self.getRing(linkage[3][0]))
            self.sugarGraph.add_node(self.getRing(linkage[3][1]))
            self.sugarGraph.add_edge(self.getRing(linkage[3][0]), self.getRing(linkage[3][1]))

        name = ""
        start = ""
        end = ""
        for node in self.sugarGraph.nodes:
            if len(list(self.sugarGraph.successors(node))) == 0:
                end = node
            elif len(list(self.sugarGraph.predecessors(node))) == 0:
                start = node

        shortestpath = nx.shortest_path(self.sugarGraph, start, end)

        for c in shortestpath:
            name += getSNFGname(self.getRingName(c))

        return name

    def findADbox(self):
        for atom in self.atom_positions.values():
            if atom[0][0] < self.xMin:
                self.xMin = atom[0][0]

            if atom[0][0] > self.xMax:
                self.xMax = atom[0][0]

            if atom[0][1] < self.yMin:
                self.yMin = atom[0][1]

            if atom[0][1] > self.yMax:
                self.yMax = atom[0][1]

            if atom[0][2] < self.zMin:
                self.zMin = atom[0][2]

            if atom[0][2] > self.zMax:
                self.zMax = atom[0][2]

        self.centreX = int((self.xMin+self.xMax)/2)
        self.centreY = int((self.yMin+self.yMax)/2)
        self.centreZ = int((self.zMin+self.zMax)/2)
        self.sizeX = int(abs(self.xMax-self.xMin))
        self.sizeY = int(abs(self.yMax-self.yMin))
        self.sizeZ = int(abs(self.zMax-self.zMin))

    def getADBox(self):
        '''

        the bounds of the ligand plus 5 angstroms (2.5 angstroms either side of the box) to allow for
        flexibility of hydrogen bonding residues.

        :return: array size 6 of ints corresponding to centre x, y, z and corresponding sizes.
        '''
        return [self.centreX, self.centreY, self.centreZ, self.sizeX+5, self.sizeY+5, self.sizeZ+5]

    def makeConfig(self, path, ADbox):
        text = "center_x = {}\ncenter_y = {}\ncenter_z = {}\nsize_x = {}\nsize_y = {}\nsize_z" \
               " = {}\nenergy_range = 12\nexhaustiveness = 80\nnum_modes = 100".format(ADbox[0],ADbox[1],ADbox[2],
                                                                                       ADbox[3],ADbox[4],ADbox[5])
        with open(path+"/"+self.filename+".conf", "w") as f:
            f.write(text)


############################################################################################################

    def draw3DGraph(self):
        # nx.draw(self.graph, with_labels=True)
        # plt.show()

        # 3D network plot
        G = self.graph
        pos = self.atom_positions

        def colour_helper_function(label):
            if label != "":
                if label[0] == "C":
                    return "tab:gray"
                elif label[0] == "O" or label[0] == "L":
                    return "r"
                else:
                    return "b"
            else:
                return "b"

        def opacity_helper_function(label):
            if label != "":
                return 0.8
            else:
                return 0.1

        with plt.style.context(('ggplot')):

            fig = plt.figure(figsize=(10,7))
            ax = Axes3D(fig)

            # Loop on the pos dictionary to extract the x,y,z coordinates of each node
            for key, value in pos.items():
                xi = value[0][0]
                yi = value[0][1]
                zi = value[0][2]

                # Scatter plot
                ax.scatter(xi, yi, zi, c=colour_helper_function(value[1]), edgecolors='k', alpha=opacity_helper_function(value[1]))
                ax.text(xi, yi, zi, value[1])

            # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
            # Those two points are the extrema of the line to be plotted
            for i,j in enumerate(G.edges()):

                x = np.array((pos[j[0]][0][0], pos[j[1]][0][0]))
                y = np.array((pos[j[0]][0][1], pos[j[1]][0][1]))
                z = np.array((pos[j[0]][0][2], pos[j[1]][0][2]))


            # Plot the connecting lines
                ax.plot(x, y, z, c='black', alpha=0.5)

            # plt.title("{} Phi: {} Psi: {}".format(CarbohydrateName(self.PDBLigand).getName(),
            #                                                 self.GT.findPhi()[0][2], self.GT.findPsi()[0][2]))

            plt.axis('off')
            plt.show()

