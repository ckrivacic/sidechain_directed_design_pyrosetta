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


class ConstrainToInvRot(object):

    def __init__(self):
        self.alignment_atoms = ['N', 'CA', 'CB']
        self.pdb_path = 'test_inputs/8cho.pdb'
        # constraints = 'test_inputs/test_constraints.cst'
        self.constraints = 'test_inputs/8cho_cst_E.cst'

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

            self.nearest_residue =\
                rosetta.protocols.protein_interface_design.movers.\
                find_nearest_residue_to_coord(self.pose, bb_centroid, 1)

            rmsd = backbone_rmsd(rotamer,
                    self.pose.residue(self.nearest_residue), self.alignment_atoms)
            invrot_rmsds[rmsd] = (rotamer, self.nearest_residue)

        self.inverse_rotamer, best_rmsd = self.choosing_func(invrot_rmsds)
        print(best_rmsd)
        print(self.inverse_rotamer)
    
    def make_constraints_from_inverse_rotamer(self):
        """
        Create constraints on the pose backbone to the nearest inverse
        rotamer.
        """
        if not hasattr(self, 'inverse_rotamer'):
            print("No rotamer chosen. Maybe something went "\
            "wrong with ConstrainToInvRot.choose_rotamer.")
            return
        if not hasattr(self,'nearest_residue'):
            print("Nearest residue not chosen. Maybe something "\
            "went wrong with ConstrainToInvRot.choose_rotamer.")
            return

        inverse_rotamers = rosetta.std.list_std_shared_ptr_const_core_conformation_Residue_std_allocator_std_shared_ptr_const_core_conformation_Residue_t()
        inverse_rotamers.push_back(self.inverse_rotamer[0])
        """
        cst_generator =\
        rosetta.protocols.forge.constraints.InverseRotamersRCG(self.nearest_residue,
                self.nearest_residue, inverse_rotamers)
        print(self.pose.constraint_set())
        cst_generator.generate_remodel_constraints(self.pose)
        print(self.pose.constraint_set())

        # csts =
        # rosetta.protocols.toolbox.match_enzdes_util.constrain_pose_res_to_invrots()
        """

        sfxn = rosetta.core.scoring.constraints.BoundFunc(0, 0.05, 0.4,
                "invrot")
        seqpos = rosetta.utility.vector1_unsigned_long()
        seqpos.append(self.inverse_rotamer[1])
        #seqpos = [38]
        ambiguous_constraint = \
            rosetta.protocols.toolbox.match_enzdes_util.constrain_pose_res_to_invrots(inverse_rotamers, seqpos, self.pose, sfxn)


        #print(ambiguous_constraint)
        csts = rosetta.utility.vector1_std_shared_ptr_const_core_scoring_constraints_Constraint_t()
        csts.append(ambiguous_constraint)
        self.sfxn = create_score_function("ref2015_cst")
        #score_manager = rosetta.core.scoring.ScoreTypeManager()
        #score_term = score_manager.score_type_from_name("coordinate_constraint")
        #sfxn.set_weight(score_term, 1.0)
        self.sfxn(self.pose)
        #self.pose.add_constraints(csts)
        self.pose.add_constraint(ambiguous_constraint)
        #print(self.pose.constraint_set())

def minimize_pose(pose, residues_bb_movable, residues_sc_movable):
    mm = rosetta.core.kinematics.MoveMap()
    mm.set_chi(False)
    mm.set_bb(False)

    for i in residues_bb_movable:
        mm.set_bb(i, True)
        mm.set_chi(i, True)

    for i in residues_sc_movable:
        mm.set_chi(i, True)

    min_opts = rosetta.core.optimization.MinimizerOptions(
            "lbfgs_armijo_nonmonotone", 0.01, True )


    minmover = rosetta.protocols.minimization_packing.MinMover()
    minmover.movemap(mm)
    minmover.min_options(min_opts)
    minmover.score_function(self.sfxn)

    minmover.apply(pose)
    
    pose.dump_file('out.pdb')


def setup_movemap(residues_bb_movable, residues_sc_movable):
    mm = rosetta.core.kinematics.MoveMap()
    mm.set_chi(False)
    mm.set_bb(False)

    for i in residues_bb_movable:
        mm.set_bb(i, True)
        mm.set_chi(i, True)
    
    for i in residues_sc_movable:
        mm.set_chi(i, True)
    
    return mm


def setup_restrained_sfxn(restraint_types, weights):
    sfxn = create_score_function("ref2015_cst")
    score_manager = rosetta.core.scoring.ScoreTypeManager()
    for i in range(0, len(restraint_types)):
        score_term = score_manager.score_type_from_name(restraint_types[i])
        sfxn.set_weight(score_term, weights[i])

    return sfxn


def fast_relax(pose, residues_bb_movable, residues_sc_movable):
    '''Fast relax the pose'''
    mm = setup_movemap(residues_bb_movable, residues_sc_movable)
    sfxn = setup_restrained_sfxn(['coordinate_constraint'],[2.0])

    fast_relax_rounds = 5
    fast_relax = rosetta.protocols.relax.FastRelax(sfxn, fast_relax_rounds)
    fast_relax.set_movemap(mm) 
    
    fast_relax.apply(pose)
    pose.dump_file('out.pdb')


def fast_design(pose, residues_bb_movable, residues_sc_movable,
        resfile):
    '''Run fast design on the pose'''
    mm = setup_movemap(residues_bb_movable, residues_sc_movable)
    sfxn = setup_restrained_sfxn(['coordinate_constraint'],[2.0])


cst_test = ConstrainToInvRot()
rotamer_set = cst_test.create_inverse_rotamers('GLU')
cst_test.choose_rotamer()
cst_test.make_constraints_from_inverse_rotamer()

bb_movable = [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
        46, 47]

fast_relax(cst_test.pose,bb_movable, [])
