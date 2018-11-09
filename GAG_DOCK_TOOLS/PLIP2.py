import sys
import os
sys.path.insert(0, '/Volumes/Eric/plip-stable')
from plip.modules.preparation import PDBComplex


non_GAG_residues = ["NI", "CH2", "SO4", "TDG", 'CAC', 'ASG', 'CA', 'EPE', 'IPA', 'NA', 'PO4', 'GOL', 'A3P', 'CIT',
                    'MPD', 'ACY', '0G6', 'ZN', 'FMT', 'PEG', 'MG', 'E64', 'NO3', 'ACT', 'JHM', 'TLA', 'DTT', 'PCA',
                    'CO3', 'MN', 'PA5', 'FAD', 'THJ', 'SIA', 'PT', 'K', 'AMP',  'TL',  'MRD', 'DMJ', 'CS', 'ACP', 'CO',
                    'EDO', 'LDA', 'ACE', 'RET', 'PLM', 'C8E', 'BME', 'HG',  'FOR',  'BEK', 'ADE', 'NHE',  'BEZ', 'HEM',
                    'OS', 'NH2',  'CD', 'HEZ', 'ASO',  'VO4',  'FSM', 'PG4',  'BEN',  '6PG', 'NPO', '2PO', 'ACH', 'PEP',
                    'NO2', 'UNX',  '293', 'IFL', 'XX6',  'XX7', 'AGG', 'MPT', 'PGE', 'FE', '1PG', 'GBL','MES', 'A46',
                    'AZI', 'AVE', 'AVF', 'AVD',  'BCD', 'OH', 'CYS']

charged_pos = ["ARG", "HIS", "LYS"]
charged_neg = ["ASP", "GLU"]
polar_uncharged = ["SER", "THR", "ASN", "GLN"]
special = ["CYS", "SEC", "GLY", "PRO"]
hydrophobic = ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"]

amino_acids = charged_neg + charged_neg + polar_uncharged + special + hydrophobic

print(amino_acids)

class InteractionCounter(object):
    def __init__(self):
        self.amino_acids = ["ARG", "HIS", "LYS", "ASP", "GLU", "SER", "THR", "ASN", "GLN", "CYS", "SEC", "GLY", "PRO",
                            "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"]

        self.amino_acids.sort()

        self.HbondIDs = []
        self.SaltbridgeIDs = []

        self.salt_bridges = {}
        self.hbonds = {}
        for aa in self.amino_acids:
            self.salt_bridges[aa] = 0
            self.hbonds[aa] = 0

    def addHbond(self, residue, ID):
        if residue in amino_acids and ID not in self.HbondIDs:
            self.HbondIDs.append(ID)
            self.hbonds[residue] += 1

    def addSaltBridge(self, residue, ID):
        if residue in amino_acids and ID not in self.SaltbridgeIDs:
            self.SaltbridgeIDs.append(ID)
            self.salt_bridges[residue] += 1


    def make_line_for_CSV(self):
        s = ""
        for aa in amino_acids:
            s += str(self.salt_bridges[aa]) + ","
        for aa in amino_acids:
            s += str(self.hbonds[aa]) + ","
        s = s[:-1] + "\n"
        return s




"""
This code is for a program that loops through a folder of pdb files are writes the interactions of each residue to
a file.

"""

class PLIP_DATA(object):
    """
    A class that takes in a dictionary of interaction data and returns a line for the .csv file
    """
    def __init__(self, interaction_dictionary):
        self.interaction_dictionary = interaction_dictionary
        #print(self.interaction_dictionary)
        self.ligand_residues = self.interaction_dictionary.lig_members
        self.ligand_id_chain_code = []

        for res_name, chain, id in self.ligand_residues:

            self.ligand_id_chain_code.append(str(id)+chain)

        self.saltbridges_list = self.interaction_dictionary.saltbridge_lneg
        self.hbond_ligand_donor_list = self.interaction_dictionary.all_hbonds_ldon
        self.hbond_protein_donor_list = self.interaction_dictionary.all_hbonds_pdon

        """
        Dictionaries where the key is the ligand ID and containing a list of interaction objects from PLIP
        """
        self.saltbridges_dict = {}
        self.hbond_ligand_donor_dict = {}
        self.hbond_protein_donor_dict = {}

        self.fill_hbond_protein_donor_dict()
        self.fill_hbond_ligand_donor_dict()
        self.fill_saltbridge_dict()

        print(self.saltbridges_dict.keys())


        self.csv_strings = []

        for name, chain, resid in self.ligand_residues:
            interaction_count = InteractionCounter()
            print(resid)

            if resid in self.saltbridges_dict.keys():
                for saltbridge in self.saltbridges_dict[resid]:
                    print(saltbridge)
                    interaction_count.addSaltBridge(saltbridge.restype, saltbridge.resnr)

            if resid in self.hbond_protein_donor_dict.keys():
                for hbonds in self.hbond_protein_donor_dict[resid]:
                    print(hbonds)
                    interaction_count.addHbond(hbonds.restype, hbonds.resnr)

            if resid in self.hbond_ligand_donor_dict.keys():
                for hbonds in self.hbond_ligand_donor_dict[resid]:
                    print(hbonds)
                    interaction_count.addHbond(hbonds.restype, hbonds.resnr)

            print(interaction_count.make_line_for_CSV())


            #self.get_saltbridge_residue_counts(resid[2])
            #self.csv_strings.append(self.make_csv_string(resid[2]))



    def get_csv_string(self):
        return self.csv_strings


    def make_csv_string(self, residue_id):
        """

        :param residue_id: residue id to construct string
        :return: csv string in the form: int type, protein res, ligand atom type, distance

        """

        string = ""


        #for dict in [self.saltbridges_dict, self.hbond_ligand_donor_dict, self.hbond_protein_donor_dict]:
        smallest_distance = 1000
        second_smallest_distance = 1000

        smallest = 0
        second_smallest = 0

        if residue_id in self.saltbridges_dict.keys():
            for saltbridge in self.saltbridges_dict[residue_id]:


                if saltbridge.distance < smallest_distance:

                    second_smallest_distance = smallest_distance
                    second_smallest = smallest

                    smallest_distance = saltbridge.distance
                    smallest = saltbridge

                elif saltbridge.distance < second_smallest_distance:

                    second_smallest_distance = saltbridge.distance
                    second_smallest = saltbridge

            try:
                string += "{},{},{},{},{},{},".format(smallest.negative.fgroup, smallest.restype, smallest.distance,
                                                      second_smallest.negative.fgroup, second_smallest.restype,
                                                      second_smallest.distance)
            except:
                try:
                    string += "{},{},{},none,none,none,".format(smallest.negative.fgroup, smallest.restype, smallest.distance)
                except:
                    string += "none,none,none,none,none,none,"
        else:
            string+="none,none,none,none,none,none,"

        smallest_distance = 1000
        second_smallest_distance = 1000

        smallest = 0
        second_smallest = 0

        if residue_id in self.hbond_ligand_donor_dict.keys():
            for saltbridge in self.hbond_ligand_donor_dict[residue_id]:


                if saltbridge.distance_ah < smallest_distance:

                    second_smallest_distance = smallest_distance
                    second_smallest = smallest

                    smallest_distance = saltbridge.distance_ah
                    smallest = saltbridge

                elif saltbridge.distance_ah < second_smallest_distance:

                    second_smallest_distance = saltbridge.distance_ah
                    second_smallest = saltbridge

            try:
                string += "{},{},{},{},{},{},".format(smallest.atype, smallest.restype, smallest.distance_ah,
                                                      second_smallest.atype, second_smallest.restype,
                                                      second_smallest.distance_ah)
            except:
                try:
                    string += "{},{},{},none,none,none,".format(smallest.atype, smallest.restype, smallest.distance_ah)
                except:
                    string += "none,none,none,none,none,none,"
        else:
            string+="none,none,none,none,none,none,"

        smallest_distance = 1000
        second_smallest_distance = 1000

        smallest = 0
        second_smallest = 0

        if residue_id in self.hbond_protein_donor_dict.keys():
            for saltbridge in self.hbond_protein_donor_dict[residue_id]:


                if saltbridge.distance_ah < smallest_distance:

                    second_smallest_distance = smallest_distance
                    second_smallest = smallest

                    smallest_distance = saltbridge.distance_ah
                    smallest = saltbridge

                elif saltbridge.distance_ah < second_smallest_distance:

                    second_smallest_distance = saltbridge.distance_ah
                    second_smallest = saltbridge

            try:
                string += "{},{},{},{},{},{},".format(smallest.atype, smallest.restype, smallest.distance_ah,
                                                      second_smallest.atype, second_smallest.restype,
                                                      second_smallest.distance_ah)
            except:
                try:
                    string += "{},{},{},none,none,none,".format(smallest.atype, smallest.restype, smallest.distance_ah)
                except:
                    string += "none,none,none,none,none,none,"
        else:
            string+="none,none,none,none,none,none,"


        return string[:-1]+"\n" # to avoid including the last comma


    def fill_saltbridge_dict(self):
        for saltbridge in self.saltbridges_list:
            if saltbridge.resnr_l not in self.saltbridges_dict.keys():
                self.saltbridges_dict[saltbridge.resnr_l] = []
            self.saltbridges_dict[saltbridge.resnr_l].append(saltbridge)

    def fill_hbond_ligand_donor_dict(self):
        for hbond in self.hbond_ligand_donor_list:
            if hbond.type == 'strong' and hbond.resnr_l not in self.hbond_ligand_donor_dict.keys():
                self.hbond_ligand_donor_dict[hbond.resnr_l] = []
            if hbond.type == 'strong':
                self.hbond_ligand_donor_dict[hbond.resnr_l].append(hbond)

    def fill_hbond_protein_donor_dict(self):
        for hbond in self.hbond_protein_donor_list:
            if hbond.type == 'strong' and hbond.resnr_l not in self.hbond_protein_donor_dict.keys():
                self.hbond_protein_donor_dict[hbond.resnr_l] = []
            if hbond.type == 'strong':
                self.hbond_protein_donor_dict[hbond.resnr_l].append(hbond)



path = "/Volumes/Eric/Desktop/docking/protein-ligand-complexes/"

file = open("GAG_interactions.txt", "w")

interactions_csv = open("GAG_interactions.csv", "w")
interactions_csv.write("Salt Bridge 1 atom type, Res, "
                       "distance, Salt Bridge 2 atom type, Res, distance, HB ligand 1 atom type, Res, distance, "
                       "HB ligand 2 atom type, Res, distance, HB protein 1 atom type, Res, distance, HB protein 2 "
                       "atom type, Res, distance,")
count = 1

# loops through the files in the ../../pdb path
for ligand in os.listdir(path):
    if ligand.__contains__(".pdb"):
        my_mol = PDBComplex()
        file.write("LIGAND\n")
        print("{} {} out of {}".format(ligand, count, len(os.listdir(path))))
        count += 1
        file.write(ligand)
        file.write("\n")
        # try:
        my_mol.load_pdb(path+ligand) # Load the PDB file into PLIP class
        my_mol.analyze()
        residues = list(my_mol.interaction_sets.items())
        for res, interaction in residues:
            if not res.split(":")[-1] in non_GAG_residues:
                plip = PLIP_DATA(my_mol.interaction_sets[res])
                for string in plip.get_csv_string():
                    interactions_csv.write(string)
                    break


        #     print(res)
        #     file.write("RESIDUE \n")
        #     file.write(res)
        #     file.write("\n")
        #     file.write("INTERACTIONS \n")
        #     my_interactions = my_mol.interaction_sets[res]
        #     print(my_interactions.__dict__)
        #     for key, value in my_interactions.__dict__.items():
        #         print("{} {}".format(key, value))
        #
        #     for interaction in my_interactions.all_itypes:
        #
        #         print(interaction.restype)
        #
        #         print(str(interaction).split("(")[0])
        #
        #         try:
        #             print(interaction.distance)
        #         except:
        #             pass
        #         try:
        #             print(interaction.distance_ad)
        #         except:
        #             pass
        #         file.write("{} {} \n".format(str(interaction).split("(")[0], interaction.restype))
        # # except:
        # #     pass

