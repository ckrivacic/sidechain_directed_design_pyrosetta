'''
Utility file for creating a constraints file from a PDB and a residue type.
After downloading, gets coordinates from the motif atoms on the of interest.

Usage:
    python3 create_constraints.py <pdbid> <resnum> <chain> <comma_separated_atoms>


Downloads from https://files.rcsb.org/download/PDBID.pdb
'''

import sys, os, wget
from utils import *
from pyrosetta import *
from numeric import *
from align import Alignment
#from klab import docopt


def constraints_from_pose(reference_pose, res_dict):
    """
    Make constraints from a pose given the resnum and list of atoms.
    res_dict should have the following structure:
    {resnum:[atom_list]}
    ex:
    {38:['N','C','CA']}
    Note: This is intended to return the constraints, which should then be
    applied to a DIFFERENT pose.
    """
    func = rosetta.core.scoring.constraints.BoundFunc(0, 0.05, 0.4, "coordcst")
    coordinate_constraints = []
    for resi in res_dict:
        for atom in res_dict[resi]:
            xyzV = reference_pose.residue(resi).xyz(atom) 
            fixed_pt = reference_pose.atom_tree().root().atom_id()
            atomno = reference_pose.residue(resi).atom_index(atom)
            atom_id = rosetta.core.id.AtomID(atomno, resi)
            coordinate_constraint = \
                    rosetta.core.scoring.constraints.CoordinateConstraint(\
                    atom_id, fixed_pt, xyzV, func
                    )
            coordinate_constraints.append(coordinate_constraint)
    
    return coordinate_constraints


def constrain_mutant_to_wt(mutant_pose, wt_pose, focus_residues):
    aligner = Alignment(mutant_pose, wt_pose)
    aligner.create_shell(10, focus_residues)
    aligner.match_align()
    alignment_dict = {}
    for focus_residue in focus_residues:
        alignment_dict[focus_residue] = ['N','C','CA']
    csts = constraints_from_pose(aligner.mobile, alignment_dict)
    for cst in csts:
        print(cst)
        aligner.target.add_constraint(cst)


init()
mutant_pose = pose_from_pdb('1qjg.pdb')
wt_pose = pose_from_pdb('8cho.pdb')
constrain_mutant_to_wt(mutant_pose, wt_pose, [38])
