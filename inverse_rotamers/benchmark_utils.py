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


def constraints_from_pose(reference_pose, target_pose, res_dict):
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
            xyzV = reference_pose.residue(resi[0]).xyz(atom)
            fixed_pt = reference_pose.atom_tree().root().atom_id()
            atomno = target_pose.residue(resi[1]).atom_index(atom)
            atom_id = rosetta.core.id.AtomID(atomno, resi[1])
            coordinate_constraint = \
                    rosetta.core.scoring.constraints.CoordinateConstraint(\
                    atom_id, fixed_pt, xyzV, func
                    )
            coordinate_constraints.append(coordinate_constraint)
    
    return coordinate_constraints


def constrain_mutant_to_wt(mutant_pose, wt_pose, mut_focus_residues,
        wt_focus_residues, constrain=True, shell=6):
    aligner = Alignment(mutant_pose, wt_pose)

    # Get chain info (this is a roundabout way of getting this and
    # should probably be refactored at some point.)
    # Initially this just used the Rosetta chain, but (older?) versions
    # of PyRosetta do not seem to support that.
    tchain = mutant_pose.pdb_info().chain(mut_focus_residues[0])
    mchain = wt_pose.pdb_info().chain(wt_focus_residues[0])

    # Create a selection shell
    aligner.create_shell(shell, mut_focus_residues,
            mobile_focus_list=wt_focus_residues,
            tchain=tchain,
            mchain=mchain)
    aligner.match_align()
    if constrain:
        alignment_dict = {}
        assert(len(mut_focus_residues) == len(wt_focus_residues))
        for i in range(0, len(mut_focus_residues)):
            alignment_dict[(wt_focus_residues[i],mut_focus_residues[i])] = ['N','C','CA']
        csts = constraints_from_pose(aligner.mobile, aligner.target, alignment_dict)
        for cst in csts:
            aligner.target.add_constraint(cst)

    return aligner


def pose_from_rcsb(pdbid, prefix=None):
    if prefix:
        path = prefix + '/' + pdbid
    else:
        path = pdbid
    if not os.path.exists(path + '.pdb'):
        url = 'https://files.rcsb.org/download/' + pdbid + '.pdb'
        wget.download(url, path + '.pdb')
    pyrosetta.toolbox.cleanATOM(path + '.pdb')
    pose = rosetta.core.import_pose.get_pdb_and_cleanup(path + '.clean.pdb')

    return pose

def pose_from_netapp_pdb(pdbid,\
        prefix='/netapp/database/pdb/remediated/mmCIF/', suffix='.cif.gz'):
    path = os.path.join(prefix, pdbid[1:3], pdbid + suffix)
    print('opening file ' + path)
    pose = pose_from_file(path)
    return pose


def custom_open(pdbid,\
        prefix='/wynton/home/kortemme/krivacic/intelligent_design/sidechain_directed_design_pyrosetta/backrub_pointmutant_benchmark/benchmark_pdbs',
        suffix='.pdb'):
    path = os.path.join(prefix, pdbid + suffix)
    print('opening file ' + path)
    pose = pose_from_file(path)
    return pose


def pose_from_pdbredo(pdbid, prefix):
    path = os.path.join(prefix, pdbid[1:3], pdbid,\
            pdbid + '_final_tot.pdb')
    if not os.path.exists(path[:-4] + '.clean.pdb'):
        toolbox.cleanATOM(path)
    print(path[:-4] + '.clean.pdb')
    pose = rosetta.core.import_pose.get_pdb_and_cleanup(path[:-4] \
            + '.clean.pdb')
    return pose


def prepare_pdbids_for_modeling(wt_pdbid, mut_pdbid, focus_mismatch_list,
        prefix='/wynton/home/kortemme/krivacic/intelligent_design/sidechain_directed_design_pyrosetta/backrub_pointmutant_benchmark/benchmark_pdbs', constrain=True, shell=6.0):
    '''
    Update info below.

    #Given 2 pdbs and a motif dict, prepare them for modeling (get designable
    #residues and task factory).
    #Motif dict should have the following structure:
    #{int(resnum):one_letter_resname}
    #Ex:
    #{4:'E',5:'V'}
    '''
    wt_pose = custom_open(wt_pdbid, prefix=prefix)
    '''
    wt_pose = pose_from_pdbredo(wt_pdbid,
            prefix=prefix)
    '''
    mut_pose = custom_open(mut_pdbid, prefix=prefix)
    '''
    mut_pose = pose_from_pdbredo(mut_pdbid,
            prefix=prefix)
    '''

    mut_pair = MutantPair(mut_pose, wt_pose, focus_mismatch_list, 
            shell=shell, cst=constrain)

    focus_list = [x.target for x in focus_mismatch_list]
    designable, repackable = choose_designable_residues(mut_pair.mut.pose,
            focus_list, dshell=shell, pshell=shell+(shell/2.0))

    task_factory = setup_task_factory(mut_pair.mut.pose, designable, repackable,
            motif_dict=mut_pair.motif_dict, layered_design=False,
            prepare_focus=True)

    return designable, repackable, task_factory, mut_pair


class TestPose(object):
    """
    Class for storing information about a pose to be tested compared to mutant. 
    Information includes focus residues, surrounding mutations, etc.
    """
    def __init__(self, pose):
        self.pose_ = pose
        self.focus_ = None
        self.mut_positions_ = None

    @property
    def pose(self):
        return self.pose_

    @pose.setter
    def pose(self, pose):
        self.pose_ = pose

    @property
    def focus(self):
        return self.focus_

    @focus.setter
    def focus(self, focus):
        if not self.pose_.residue(focus):
            raise ValueError("Focus residue is not valid for this pose.")
        self.focus_ = focus

    @property
    def mobile(self):
        return self.mobile_

    @mobile.setter
    def mobile(self, boolean):
        self.mobile_ = boolean

    @property
    def mismatches(self):
        """
        Mutations should be Mismatch object from align.py
        """
        return self.mut_positions_

    @mismatches.setter
    def mismatches(self, muts):
        mismatches = []
        print('Parsing mismatches')
        for mut in muts:
            if self.mobile:
                mismatches.append(mut.mobile)
            else:
                mismatches.append(mut.target)
        self.mut_positions_ = mismatches


class MutantPair(object):
    """Stores pair of mutants to be tested"""
    def __init__(self, mut_pose, wt_pose, focus_mismatches, shell=6, cst=False):
        self.mut_ = TestPose(mut_pose)
        self.mut_.mobile = False
        self.wt_ = TestPose(wt_pose)
        self.wt_.mobile = True

        target_focus_list = [x.target for x in focus_mismatches]
        mobile_focus_list = [x.mobile for x in focus_mismatches]
        self.aligner_ = constrain_mutant_to_wt(self.mut_.pose,
                self.wt_.pose,
                target_focus_list, mobile_focus_list, constrain=cst,
                shell=shell)

        self.mut_.mismatches = self.aligner_.mismatches
        self.wt_.mismatches = self.aligner_.mismatches

        motif_dict = {}
        for mismatch in self.aligner_.mismatches:
            wt_res_type = self.wt_.pose.residue(mismatch.mobile).name1()
            motif_dict[mismatch.target] = str(wt_res_type)
        self.motif_dict_ = motif_dict

    @property
    def mut(self):
        return self.mut_

    @mut.setter
    def mut(self, pose):
        self.mut_ = pose

    @property
    def wt(self):
        return self.wt_

    @wt.setter
    def wt(self, pose):
        self.wt_ = pose

    @property
    def aligner(self):
        return self.aligner_

    @aligner.setter
    def aligner(self, aligner):
        self.aligner_ = aligner

    @property
    def motif_dict(self):
        return self.motif_dict_

    @motif_dict.setter
    def motif_dict(self, motif_dict):
        self.motif_dict_ = motif_dict


def finish_io(temp, final, suffix=''):
    from shutil import copyfile
    pdbfile = 'pdbs_' + suffix + '.tar.gz'
    cmd = 'tar -czvf ' + os.path.join(temp, pdbfile) + ' --directory=' + temp + ' .'
    print('RUNNING COMMAND ' + cmd)
    os.system(cmd)
    if not os.path.exists(final):
        os.makedirs(final, exist_ok=True)
    copyfile(os.path.join(temp, pdbfile),
            os.path.join(final, pdbfile))
# For testing use the following command
# stuff = prepare_pdbids_for_modeling('8cho','1qj4',mm, prefix='/home/ckrivacic/pdb-redo')
