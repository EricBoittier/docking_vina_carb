from LigandPDB import *
import numpy as np



def calculateRMSD(ligand1, ligand2):
    lis = []
    ligand1_atoms = ligand1.getAtomsToCompare()
    ligand2_atoms = ligand2.getAtomsToCompare()
    sumOfDifSqr = 0
    count = 0

    for items in ligand2_atoms:
        print(items)
    print("")
    for items in ligand1_atoms:
        print(items)

    # print("Comparisons:")
    for ring1hash, ring1 in ligand1_atoms.items():
        print(ring1hash)
        print(ligand2_atoms)
        for ring_1_atom in ring1:
            for ring_2_atom in ligand2_atoms[ring1hash]:
                if ring_1_atom.getID() == ring_2_atom.getID() and ring_1_atom.getAtomType() == ring_2_atom.getAtomType():
                    print("{} {}".format(ring_1_atom.getID(), ring_2_atom.getID()))
                    print("{} {}".format(ring_1_atom.getAtomType(), ring_2_atom.getAtomType()))

                    if not (ring_1_atom.getXYZ(), ring_2_atom.getXYZ()) or not (ring_2_atom.getXYZ(),
                                                                            ring_1_atom.getXYZ()) in lis:
                        lis.append((ring_1_atom.getXYZ(), ring_2_atom.getXYZ()))


    for atom1, atom2 in lis:
        for value in difSqr(atom1, atom2):
            print(value)
            sumOfDifSqr += value
            count += 1
        print("")

    print("count {}".format(count))
    print((count/3)/5/2)

    return np.sqrt(sumOfDifSqr / count)


def difSqr(predicted, actual):
    print("{} {}".format(predicted, actual))
    predicted = np.array(predicted)
    actual = np.array(actual)
    x = predicted - actual
    return x * x


def getRMSD(pdb_string_1, pdb_string_2):
    try:
        referenceLigand = RMSD(pdb_string_1)
        comparisonLigand = RMSD(pdb_string_2)
        return "RMSD = {}\n".format(calculateRMSD(referenceLigand, comparisonLigand))
    except IndexError:
        pass


if __name__ == '__main__':
    import sys
    try:
        referenceLigand = LigandPDB(sys.argv[1])
        comparisonLigand = LigandPDB(sys.argv[1])
        print("RMSD = {}\n".format(calculateRMSD(
            referenceLigand, comparisonLigand)))
    except IndexError:
        pass

