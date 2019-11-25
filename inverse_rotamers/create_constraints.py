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
#from klab import docopt


def create_constraints(pdbid, res_dict):
    '''
    Make constraints from a pdb given a resnum, chain and list of atoms.
    res_dict should have the following structure:
    {resnum+chain:[atom_list]}
    ex:
    {'38A':['N','CA','CB'],'56A':['N','CA','CB']}
    '''
    if not os.path.exists(pdbid + '.pdb'):
        url = 'https://files.rcsb.org/download/' + pdbid + '.pdb'
        wget.download(url, pdbid + '.pdb')

    def parse_line(line):
        split = line.split()
        if len(split) >= 9:
            atom_name = split[2]
            resn = split[3]
            chain = split[4]
            resnum = split[5]
            xyz = [split[6],split[7],split[8]]

            return atom_name, chain, resnum, xyz
        else:
            return None, None, None, None

    constraints_list = []
    f = open(pdbid + '.pdb')
    for line in f:
        if line.split()[0] == 'ATOM':
            atom_name, chain, resnum, xyz = parse_line(line)
            resi = resnum + chain
            if resi in res_dict:
                if atom_name in res_dict[resi]:
                    print(xyz)
                    constraints_list.append('CoordinateConstraint ' + atom_name
                            + ' ' + str(resnum) + ' ' + atom_name + ' 1 ' + 
                            str(xyz[0]) + ' ' + str(xyz[1]) + ' ' + str(xyz[2]))
                    print(constraints_list)


def align_poses_and_constrain(mpose, tpose, res_dict, shell=None):
    '''
    Align 2 poses and create constraints for a target residue and atom list.
    res_dict should have the following structure:
    {resnum:[atom_list]}
    ex:
    {'38':['N','CA','CB'],'56A':['N','CA','CB']}
    '''
    #aligner = rosetta.protocols.stepwise.modeler.align.StepWisePoseAligner(mpose)
    focus_residues = ''
    for key in res_dict:
        focus_residues += key

    aligner =\
            rosetta.protocols.stepwise.modeler.align.StepWisePoseAligner(tpose)
    if shell:
        focus_residue_selector =\
                rosetta.core.select.residue_selector.ResidueIndexSelector(focus_residues)
        selector =\
                rosetta.core.select.residue_selector.NeighborhoodResidueSelector(focus_residue_selector,
                10.0, include_focus_in_subset=False)

        mresidues = intlist_to_vector1_size(res_selector_to_size_list(selector.apply(mpose)))
        tresidues = intlist_to_vector1_size(res_selector_to_size_list(selector.apply(tpose)))
    else:
        mresidues = intlist_to_vector1_size([n for n in range(1, mpose.size()+1)])
        tresidues = intlist_to_vector1_size([n for n in range(1, tpose.size()+1)])

    partition_res = rosetta.protocols.stepwise.modeler.\
            figure_out_root_partition_res(tpose, tresidues)
    if len(partition_res) == 0:
        print('no partition res found')
        partition_res.append(tpose.fold_tree().root())
        
    #aligner.set_root_partition_res(intlist_to_vector1_size(tresidues))
    print("partition residues:", partition_res)

    aligner.set_root_partition_res(partition_res)

    full_sequence = rosetta.core.pose.full_model_info.const_full_model_info(tpose).full_sequence()
    print(full_sequence)
    print(rosetta.core.pose.rna.remove_bracketed(full_sequence))
    print('------------')

    aligner.apply(mpose)


    

#create_constraints('8cho',{'38A':['N','CA','CB']})
init()
chemical_manager = rosetta.utility.SingletonBase_core_chemical_ChemicalManager_t
#chemical_manager = rosetta.core.chemical.ChemicalManager()
rsd_set = chemical_manager.get_instance().residue_type_set('fa_standard')
tpose = rosetta.core.import_pose.get_pdb_with_full_model_info('8cho.pdb',rsd_set)
mpose = rosetta.core.import_pose.get_pdb_with_full_model_info('1qjg.pdb',rsd_set)
mpose = mpose.split_by_chain(1)
align_poses_and_constrain(mpose, tpose, {'38':['N','CA','CB']}, shell=10)
#tpose = rosetta.core.import_pose.get_pdb_with_full_model_info('1jvm.pdb',rsd_set)
#mpose = rosetta.core.import_pose.get_pdb_with_full_model_info('3hpl.pdb',rsd_set)
#align_poses_and_constrain(mpose, tpose, {'71':['N','CA','CB']}, shell=None)
