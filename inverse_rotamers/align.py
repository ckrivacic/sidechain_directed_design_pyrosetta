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
        self.set_target_sequence()
        self.set_mobile_sequence()
        self.atoms = ['N','C','CA','O']

    def set_atoms(self, atoms):
        self.atoms = atoms

    def set_target_sequence(self):
        """Take list of residue numbers and return a sequence"""
        self.target_sequence = ''
        for resnum in self.target_residues:
            self.target_sequence += self.target.sequence(resnum, resnum)

    def set_mobile_sequence(self):
        self.mobile_sequence = ''
        for resnum in self.mobile_residues:
            self.mobile_sequence += self.mobile.sequence(resnum, resnum)

    def set_target_residues(self, tselect):
        """Setter for target residues. Takes res selector (vector of booleans)."""
        self.target_residues = res_selector_to_size_list(tselect)
        self.set_target_sequence()

    def set_mobile_residues(self, mselect):
        """Setter for mobile residues. Takes res selector (vector of booleans)."""
        self.mobile_residues = res_selector_to_size_list(mselect)
        self.set_mobile_sequence()

    def create_shell(self, shell, focus_residue_list, mobile_focus_list=None):
        focus_residues = ''
        for resi in focus_residue_list:
            focus_residues += (str(resi) + ',')
        if mobile_focus_list:
            mobile_focus_residues = ''
            for resi in mobile_focus_list:
                mobile_focus_residues += (str(resi) + ',')
        else:
            mobile_focus_residues = focus_residues

        focus_residue_selector = rosetta.core.select.\
                residue_selector.ResidueIndexSelector(focus_residues)
        selector = rosetta.core.select.residue_selector.\
                NeighborhoodResidueSelector(focus_residue_selector,\
                shell, include_focus_in_subset=True)

        self.set_target_residues(selector.apply(self.target))
        self.set_mobile_residues(selector.apply(self.mobile))

    def match_align(self, match_only=True):
        """
        Figure out which residue numbers to align based on sequence alignment.
        """
        alignments = pairwise2.align.localxs(self.target_sequence,\
                self.mobile_sequence, -0.1, -0.1)
        #for alignment in alignments:
        tseq = alignments[0][0]
        mseq = alignments[0][1]
        tmatch = []
        mmatch = []
        t_iter = self.target_residues.__iter__()
        m_iter = self.mobile_residues.__iter__()
        target_align_residues = []
        mobile_align_residues = []
        for i in range(len(tseq)):
            t = tseq[i]
            m = mseq[i]
            if t == '-':
                next(m_iter)
            elif m == '-':
                next(t_iter)
            elif match_only and t != m:
                next(t_iter)
                next(m_iter)
            else:
                target_align_residues.append(next(t_iter))
                mobile_align_residues.append(next(m_iter))
        superimpose_poses_by_residues(self.target, target_align_residues,\
                self.mobile, mobile_align_residues)


class Transformation(object):
    """Class for storing rotation & transformation information"""

    def __init__(self,rotation,translation):
        self.rotation = rotation
        self.translation = translation


def get_backbone_points(pose, residues):
    '''Get backbone points for residues in a pose.'''
    points = []

    for res in residues:
        for atom in ['N', 'CA', 'C']:
            points.append(xyzV_to_np_array(pose.residue(res).xyz(atom)))

    return points


def superimpose_poses_by_residues(pose_source, residues_source, pose_target, residues_target):
    '''Superimpose residues in a source pose into residues in a target pose.
    Only backbone atoms are used for the superimposition.
    '''
    assert(len(residues_source) == len(residues_target))

    # Get the points to be superimposed

    points_source = get_backbone_points(pose_source, residues_source)
    points_target = get_backbone_points(pose_target, residues_target)

    # Get the rigid body transformation

    M, t = get_superimpose_transformation(points_source, points_target)

    # Transform the source pose

    pose_source.apply_transform_Rx_plus_v(np_array_to_xyzM(M), 
            np_array_to_xyzV(t))


def calc_backbone_RMSD(pose1, residues1, pose2, residues2):
    '''Calculate backbone RMSD between two poses for specific positions.'''
    assert(len(residues1) == len(residues2))

    def RMSD(points1, poinsts2):
        '''Calcualte RMSD between two lists of numpy points.'''
        diff = [points1[i] - poinsts2[i] for i in range(len(points1))]
        return np.sqrt(sum(np.dot(d, d) for d in diff) / len(diff))

    points1 = get_backbone_points(pose1, residues1)
    points2 = get_backbone_points(pose2, residues2)

    return RMSD(points1, points2)


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
align.match_align()
