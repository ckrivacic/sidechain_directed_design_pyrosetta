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


def constrain_mutant_to_wt(mutant_pose, wt_pose, mut_focus_residues,
        wt_focus_residues, constrain=True):
    aligner = Alignment(mutant_pose, wt_pose)
    aligner.create_shell(10, mut_focus_residues, mobile_focus_list=wt_focus_residues)
    aligner.match_align()
    alignment_dict = {}
    for focus_residue in mut_focus_residues:
        alignment_dict[focus_residue] = ['N','C','CA']
    if constrain:
        csts = constraints_from_pose(aligner.mobile, alignment_dict)
        for cst in csts:
            aligner.target.add_constraint(cst)

    return aligner


def pose_from_rcsb(pdbid, prefix=None):
    if prefix:
        path = prefix + '/' + pdbid
    else:
        path = pdbid
    if not os.path.exists(path):
        url = 'https://files.rcsb.org/download/' + pdbid + '.pdb'
        wget.download(url, path + '.pdb')
    pyrosetta.toolbox.cleanATOM(path + '.pdb')
    pose = rosetta.core.import_pose.get_pdb_and_cleanup(path + '.clean.pdb')

    return pose


def pose_from_pdbredo(pdbid, prefix):
    path = os.path.join(prefix, pdbid[1:3], pdbid,\
            pdbid + '_final_tot.pdb')
    if not os.path.exists(path[:-4] + '.clean.pdb'):
        toolbox.cleanATOM(path)
    pose = rosetta.core.import_pose.get_pdb_and_cleanup(path[:-4] \
            + '.clean.pdb')


def prepare_pdbid_for_modeling(wt_pdbid, mut_pdbid, motif_dict,
        mut_focus, wt_focus, prefix=None, constrain=True):
    '''
    Given 2 pdbs and a motif dict, prepare them for modeling (get designable
    residues and task factory).
    Motif dict should have the following structure:
    {int(resnum):one_letter_resname}
    Ex:
    {4:'E',5:'V'}
    '''
    wt_pose = pose_from_pdbredo(wt_pdbid,
            prefix='/netapp/home/krivacic/pdb_redo')
    mut_pose = pose_from_pdbredo(mut_pdbid,
            prefix='/netapp/home/krivacic/pdb_redo')
    focus_list = [key for key in motif_dict]
    designable, repackable = choose_designable_residues(mut_pose,
            focus_list)
    task_factory = setup_task_factory(mut_pose, designable, repackable,
            motif_dict=motif_dict, layered_design=False,
            prepare_focus=True)
    aligner = constrain_mutant_to_wt(mut_pose, wt_pose, [mut_focus],
            [wt_focus], constrain=constrain)
    return designable, repackable, task_factory, aligner

if __name__=='main':
    init()
    pose, designable, repackable, task_factory = prepare_pdbid_for_modeling('4s0w','1cv1',111,'V')
    start_time = time.time()
    loopmodeler = get_loop_modeler(pose, designable, repackable, 111, task_factory=task_factory,
            fast=True, mover='ngk', resbuffer=4)
    end_time = time.time() - start_time
    print('total time elapsed: ',end_time)
    pose.dump_file('out.pdb')
