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
from inverse_rotamers import *
import time
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


def pose_from_rcsb(pdbid):
    if not os.path.exists(pdbid + '.pdb'):
        url = 'https://files.rcsb.org/download/' + pdbid + '.pdb'
        wget.download(url, pdbid + '.pdb')

    return pose_from_file(pdbid + '.pdb')


def prepare_pdbid_for_modeling(wt_pdbid, mut_pdbid, focus_resnum, focus_restype):
    wt_pose = pose_from_rcsb(wt_pdbid)
    mut_pose = pose_from_rcsb(mut_pdbid)
    designable, repackable = choose_designable_residues(mut_pose,
            [focus_resnum])
    task_factory = setup_task_factory(mut_pose, designable, repackable,
            motif_dict={focus_resnum:focus_restype}, layered_design=False,
            prepare_focus=True)
    constrain_mutant_to_wt(mut_pose, wt_pose, [focus_resnum])
    return mut_pose, designable, repackable, task_factory

if __name__=='main':
    init()
    pose, designable, repackable, task_factory = prepare_pdbid_for_modeling('4s0w','1cv1',111,'V')
    start_time = time.time()
    loopmodeler = get_loop_modeler(pose, designable, repackable, 111, task_factory=task_factory,
            fast=True, mover='ngk', resbuffer=4)
    end_time = time.time() - start_time
    print('total time elapsed: ',end_time)
    pose.dump_file('out.pdb')
