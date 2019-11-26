import numpy as np
from pyrosetta import *
from Bio import pairwise2
from utils import *


class Alignment(object):
    """Class for storing alignment information."""
    def __init__(self, tpose, mpose, tselect=[], mselect=[]):
        """tselect and mselect should be residue selectors (vectors of booleans)"""
        self.target = tpose
        self.mobile = mpose
        self.target_residues = res_selector_to_size_list(tselect)
        self.mobile_residues = res_selector_to_size_list(mselect)
        self.initiate_sequences()
        self.atoms = ['N','C','CA','O']

    def set_atoms(self, atoms):
        self.atoms = atoms

    def initiate_sequences(self):
        """Take list of residue numbers and return a sequence"""
        self.target_sequence = ''
        for resnum in self.target_residues:
            self.target_sequence += self.target.sequence(resnum, resnum)
        self.mobile_sequence = ''
        for resnum in self.mobile_residues:
            self.mobile_sequence += self.mobile.sequence(resnum, resnum)

    def set_target_residues(self, tselect):
        """Setter for target residues. Takes res selector (vector of booleans)."""
        self.target_residues = res_selector_to_size_list(tselect)
        self.initiate_sequences()

    def set_mobile_residues(self, mselect):
        """Setter for mobile residues. Takes res selector (vector of booleans)."""
        self.mobile_residues = res_selector_to_size_list(mselect)
        self.initiate_sequences()

    def create_shell(self, shell, focus_residue_list):
        focus_residues = ''
        for resi in focus_residue_list:
            focus_residues += (str(resi) + ',')

        focus_residue_selector = rosetta.core.select.\
                residue_selector.ResidueIndexSelector(focus_residues)
        selector = rosetta.core.select.residue_selector.\
                NeighborhoodResidueSelector(focus_residue_selector,\
                shell, include_focus_in_subset=True)

        self.target_residues =\
                res_selector_to_size_list(selector.apply(self.target))
        self.mobile_residues =\
                res_selector_to_size_list(selector.apply(self.mobile))

        print(self.target_residues)
        print(self.mobile_residues)


class Transformation(object):
    """Class for storing rotation & transformation information"""

    def __init__(self,rotation,translation):
        self.rotation = rotation
        self.translation = translation


def get_superimpose_transformation(P1, P2):
    '''Get the superimpose transformation that transfoms a list of
    points P1 to another list of points P2.
    From XingJie Pan'''
    if len(P1) != len(P2):
        raise Exception("Sets to be superimposed must have same number of points.")

    com1 = np.mean(P1, axis=0)
    com2 = np.mean(P2, axis=0)

    R = np.dot(np.transpose(np.array(P1) - com1), np.array(P2) - com2)
    V, S, W = np.linalg.svd(R)

    if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
        V[:, -1] = -V[:, -1]

    M = np.transpose(np.array(np.dot(V, W)))

    return M, com2 - np.dot(M, com1)


def sequence_align(seq1, seq2):
    return pairwise2.align.localxs(seq1, seq2, -0.1, -0.1)


init()
tpose = pose_from_pdb('8cho.pdb')
mpose = pose_from_pdb('1qjg.pdb')
align = Alignment(tpose, mpose)
align.create_shell(10,[38])
