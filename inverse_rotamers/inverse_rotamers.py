from pyrosetta import *
#import pyrosetta.rosetta
import numpy as np
from utils import *
from parsers import parse_restraints
from numeric import *
from test_utils import plot_3d


"""
This is a test file to get the basics of inverse rotamers & constraint
generation figured out. The tasks are as follows:

    1. Import a pose. (trivial)
    2. Import constraints on a residue's functional group. (partially
    done, need to generate these constraints from relative constraints.)
    3. Generate inverse rotamers (TryRotamers, then align them probably)
    (done)
    4. Figure out what closest rotamer to bb is (done)
    5. Define bb constraints based on closest rotamer
    6. Sample backbone conformations
    7. Repack/minimize/TBD

Will start with 1 and 2, then potentially update this file to perform
the rest of the tasks. TBD. For now this is just to figure out inverse
rotamers.
"""


init("-ignore_unrecognized_res -extrachi_cutoff 0 -ex1 -ex2 -ex3 -ex4")


def nearest_backbone_rmsd(rotamer, nearest_residue,
        alignment_atoms):
    """
    Measure backbone RMSD between a rotamer and the nearest residue on
    the design protein.
    """

    distances = np.array([])
    for atom in alignment_atoms:
        rot_xyz = xyzV_to_np_array(rotamer.xyz(atom))
        near_xyz = xyzV_to_np_array(nearest_residue.xyz(atom))
        distances = np.append(distances,euclidean_distance(rot_xyz,near_xyz))

    return np.sqrt((distances**2).mean())


class ConstrainToInvRot(object):

    def __init__(self):
        self.alignment_atoms = ['N', 'CA', 'CB']
        self.pdb_path = 'test_inputs/8cho.pdb'
        # constraints = 'test_inputs/test_constraints.cst'
        self.constraints = 'test_inputs/8cho_cst_D.cst'

        self.pose = pose_from_pdb(self.pdb_path)

        dummy_pose = rosetta.core.pose.Pose()
        self.global_residue_type_set = dummy_pose.residue_type_set_for_pose()
        # self.residue_type = global_residue_type_set.get_representative_type_base_name("ASP")

        self.restraints = parse_restraints(self.constraints)

    def create_inverse_rotamers(self, resn):
        """
        Creates a set of inverse rotamers that are overlaid onto a set of
        coordinate restraints. 

        In the future, the positioning information coming from the motif
        database will be relative. To address this, I'll need 2
        transformation steps; one that aligns the motif coordinates onto my
        protein, then another to create the inverse rotamers.
        """
        restype = self.global_residue_type_set.get_representative_type_base_name(resn)
        print("creating rotamer set")
        rotamer_set = bb_independent_rotamers_extra_chi(restype)

        xyzarray = []
        restrained_atoms = []
        for restraint in self.restraints:
            xyzarray.append(restraint.coord)
            restrained_atoms.append(restraint.atom_name)

        xyzarray = np.array(xyzarray)

        for residue in rotamer_set:
            # Just to look @ before and after of full residue
            # full_res_initial = []
            # for atom in residue.atoms():
            #    full_res_initial.append(np.array(atom.xyz()))
            # full_res_initial = np.array(full_res_initial)

            residue_restrained_xyz = []
            for atom in restrained_atoms:
                residue_restrained_xyz.append(np.array(residue.xyz(atom)))
            residue_restrained_xyz = np.array(residue_restrained_xyz)

            rotation, translation = \
                    get_transformation(residue_restrained_xyz,xyzarray)

            # The following transformation should be saved in rotamer_set,
            # so that's ultimately what we'll return.
            residue.apply_transform_Rx_plus_v(np_array_to_xyzM(rotation),np_array_to_xyzV(translation))


            # check to see new residue overlay
            # residue_new_xyz = []
            # for atom in residue.atoms():
            #    residue_new_xyz.append(np.array(atom.xyz()))
            # residue_new_xyz = np.array(residue_new_xyz)
            # plot_3d(full_res_initial,xyzarray,residue_new_xyz)

        self.rotamer_set = rotamer_set

    def choosing_func(self,invrot_rmsds):
        """
        I can expand this function later to make it more sophisticated,
        but for now it will just choose the inverse rotamer with the
        lowest RMSD.
        """
        best_rmsd = None
        for rmsd in invrot_rmsds:
            if not best_rmsd or rmsd < best_rmsd:
                best_rmsd = rmsd
                best_invrot = invrot_rmsds[rmsd]
        
        return best_invrot, best_rmsd

    def choose_rotamer(self):

        invrot_rmsds = {}
        for rotamer in self.rotamer_set:
            bb_array = []
            bb_indices = rotamer.all_bb_atoms()
            for index in bb_indices:
                bb_array.append(xyzV_to_np_array(rotamer.xyz(index)))
            bb_array = np.array(bb_array)
            # I realize I'm doing xyzV to np and then back, but not sure if
            # I can find the centroid of an xyzV easily, and this works for
            # now. Consider changing if I can.
            bb_centroid = np_array_to_xyzV([sum(x)/len(x) for x in
                zip(*bb_array)])

            nearest_residue =\
                rosetta.protocols.protein_interface_design.movers.\
                find_nearest_residue_to_coord(self.pose, bb_centroid, 1)

            rmsd = nearest_backbone_rmsd(rotamer,
                    self.pose.residue(nearest_residue), self.alignment_atoms)
            invrot_rmsds[rmsd] = (rotamer, nearest_residue)

        chosen_rotamer, chosen_rmsd = self.choosing_func(invrot_rmsds)
        print(chosen_rmsd)


cst_test = ConstrainToInvRot()
rotamer_set = cst_test.create_inverse_rotamers('ASP')
cst_test.choose_rotamer()
