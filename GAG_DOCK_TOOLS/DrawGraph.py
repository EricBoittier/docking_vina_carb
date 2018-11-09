import os
import matplotlib.pyplot as plt
import numpy as np
from FindCarbohydrateName import *
from FindGlycosidicTorsions import *

class DrawGraph(object):
    """docstring for DrawGraph"""
    def __init__(self, PDBLigand):
        super(DrawGraph, self).__init__()
        self.PDBLigand = PDBLigand
        self.graph = self.PDBLigand.getGraph()
        self.atom_positions = self.PDBLigand.getPositions()
        
