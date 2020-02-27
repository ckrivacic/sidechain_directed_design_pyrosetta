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
        self.atoms = ['N','C','CA']
        self.bb_rmsd = None

    def set_atoms(self, atoms):
        self.atoms = atoms

    def set_target_sequence(self):
        """Take list of residue numbers and return a sequence"""
        self.target_sequence = ''
        target_residues = []
        if len(self.target_residues) < 1:
            for i in range(0, self.target.size()):
                target_residues.append(i+1)
            self.target_residues = target_residues
        for resnum in self.target_residues:
            self.target_sequence += self.target.sequence(resnum, resnum)

    def set_mobile_sequence(self):
        self.mobile_sequence = ''
        mobile_residues = []
        if len(self.mobile_residues) < 1:
            for i in range(0, self.mobile.size()):
                mobile_residues.append(i+1)
            self.mobile_residues = mobile_residues
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

    def create_shell(self, shell, focus_residue_list,
            mobile_focus_list=None, protein_only=True, tchain=None,
            mchain=None):
        print('Creating shell for alignment of size ' + str(shell))
        focus_residues = ''
        for resi in focus_residue_list:
            focus_residues += (str(resi) + ',')
        focus_residue_selector = rosetta.core.select.\
                residue_selector.ResidueIndexSelector(focus_residues)
        if mobile_focus_list:
            mobile_focus_residues = ''
            for resi in mobile_focus_list:
                mobile_focus_residues += (str(resi) + ',')
            mobile_focus_selector = \
                    rosetta.core.select.residue_selector.\
                    ResidueIndexSelector(mobile_focus_residues)
        else:
            mobile_focus_selector = focus_residue_selector

        target_selector = rosetta.core.select.residue_selector.\
                NeighborhoodResidueSelector(focus_residue_selector,\
                shell, include_focus_in_subset=True)
        mobile_selector = rosetta.core.select.residue_selector.\
                NeighborhoodResidueSelector(mobile_focus_selector,\
                shell, include_focus_in_subset=True)

        if protein_only:
            print('Selecting only protein residues')
            property_selector = rosetta.core.select.residue_selector.\
                    ResiduePropertySelector(rosetta.core.chemical.ResidueProperty.PROTEIN)
            target_selector = rosetta.core.select.residue_selector.\
                    AndResidueSelector(target_selector,
                            property_selector)
            mobile_selector = rosetta.core.select.residue_selector.\
                    AndResidueSelector(mobile_selector,
                            property_selector)

        if tchain:
            print('tchain is {}'.format(tchain))
            tchain_selector = rosetta.core.select.residue_selector.\
                    ChainSelector(tchain)
            target_selector = rosetta.core.select.residue_selector.\
                    AndResidueSelector(target_selector, tchain_selector)

        if mchain:
            print('mchain is {}'.format(mchain))
            mchain_selector = rosetta.core.select.residue_selector.\
                    ChainSelector(mchain)
            mobile_selector = rosetta.core.select.residue_selector.\
                    AndResidueSelector(mobile_selector, mchain_selector)


        self.set_target_residues(target_selector.apply(self.target))
        self.set_mobile_residues(mobile_selector.apply(self.mobile))

    def align_sequences(self, alignments, ind, match_only=True):

        oneletter =\
                ['A','R','N','D','B','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']

        tseq = alignments[ind][0]
        mseq = alignments[ind][1]
        tmatch = []
        mmatch = []
        t_iter = self.target_residues.__iter__()
        m_iter = self.mobile_residues.__iter__()
        target_align_residues = []
        mobile_align_residues = []
        self.mismatches = []

        for i in range(len(tseq)):
            t = tseq[i]
            m = mseq[i]
            if t != '-' and t in oneletter:
                tres = next(t_iter)
            if m != '-' and m in oneletter:
                mres = next(m_iter)
                if t != '-':
                    if t == m:
                        target_align_residues.append(tres)
                        mobile_align_residues.append(mres)
                    elif t != m:
                        self.mismatches.append(Mismatch(tres, mres))
                        if not match_only:
                            target_align_residues.append(tres)
                            mobile_align_residues.append(mres)
        superimpose_poses_by_residues(self.mobile, mobile_align_residues,\
                self.target, target_align_residues)

        bb_rmsd = calc_backbone_RMSD(self.mobile, mobile_align_residues,
                self.target, target_align_residues)

        return bb_rmsd

    def match_align(self, match_only=True):
        """
        Figure out which residue numbers to align based on sequence alignment.
        """

        print('Running match_align')

        # Note about alignment parameters: Used to have a small -0.1
        # penalty for gap creation & extension. Got rid of this in favor
        # of getting better alignments, but the downside is that for
        # cases where you want to know which residues are mutants of
        # each other, you don't get that info anymore.
        alignments = pairwise2.align.globalxs(self.target_sequence,\
                self.mobile_sequence, 0, 0)

        #for alignment in alignments:
        if not alignments:
            return False
        print(alignments)
        best_rmsd = 9999
        i = 0
        for i in range(0, len(alignments)):
            bb_rmsd = self.align_sequences(alignments, i,
                    match_only=match_only)
            if bb_rmsd < best_rmsd:
                best_rmsd = bb_rmsd
                best_i = i
        self.bb_rmsd = self.align_sequences(alignments, best_i,
                match_only=match_only)

        return True


class Mismatch(object):
    """Class for storing mismatch info."""
    def __init__(self, target_res, mobile_res):
        self.target_ = target_res
        self.mobile_ = mobile_res

    @property
    def target(self):
        return self.target_

    @target.setter
    def target(self, target_res):
        self.target_ = target_res

    @property
    def mobile(self):
        return self.mobile_

    @mobile.setter
    def mobile(self, mobile_res):
        self.mobile_ = mobile_res


def get_shell_selector(shell, focus_selector, include_focus=True):
    shell_selector = \
            rosetta.core.select.residue_selector.NeighborhoodResidueSelector(
                    focus_selector, float(shell),
                    include_focus_in_subset=include_focus
                    )
    return shell_selector



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


def test_compare(pdb1, pdb1_focus_resnum, pdb1_focus_chain, pdb2,
        pdb2_focus_resnum, pdb2_focus_chain, shell):
    from benchmark_utils import pose_from_rcsb
    pose1 = pose_from_rcsb(pdb1, prefix='temp/')
    pose1reslist = [pose1.pdb_info().pdb2pose(pdb1_focus_chain,
            pdb1_focus_resnum)]
    pose2 = pose_from_rcsb(pdb2, prefix='temp/')
    pose2reslist = [pose2.pdb_info().pdb2pose(pdb2_focus_chain,
            pdb2_focus_resnum)]
    aligner = Alignment(pose1, pose2)
    aligner.create_shell(shell, pose1reslist,
            mobile_focus_list=pose2reslist)
    aligner.match_align()
    return aligner

def test_compare2(shell):
    from benchmark_utils import pose_from_rcsb
    pose2 = pose_from_file('temp/1nhs.cif.gz')
    pose1 = pose_from_file('temp/1npx.cif.gz')
    pose2reslist = [pose2.pdb_info().pdb2pose('A',
            41)]
    pose1reslist = [pose1.pdb_info().pdb2pose('A',
            41)]
    aligner = Alignment(pose1, pose2)
    aligner.create_shell(shell, pose1reslist,
            mobile_focus_list=pose2reslist)
    aligner.match_align()
    return aligner

def sequence_align(seq1, seq2):
    return pairwise2.align.localxs(seq1, seq2, -0.1, -0.1)
