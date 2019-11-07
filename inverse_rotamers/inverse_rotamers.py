from pyrosetta import *
#import pyrosetta.rosetta
import numpy as np
from utils import *
from parsers import parse_restraints
from numeric import *
from test_utils import *
from itertools import compress
from mover_utils import *


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


init("-ignore_unrecognized_res -extrachi_cutoff 0 -ex1 -ex2 -out:overwrite " +\
        "-lh:db_path=/home/ckrivacic/rosetta/database/loophash_db/ " +\
        "-lh:loopsizes 6 -out:level 400")


class ConstrainToInvRot(object):

    def __init__(self):
        self.alignment_atoms = ['N', 'CA', 'CB']
        self.pdb_path = 'test_inputs/8cho_clean_relaxed.pdb'
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

        #i = 0
        for residue in rotamer_set:
            # dump_invrot(residue, 'nontransformed_invrots/invrot_' + str(i) +
            #         '.pdb')
            # Just to look @ before and after of full residue
            # full_res_initial = []
            # for atom in residue.atoms():
            #    full_res_initial.append(np.array(atom.xyz()))
            # full_res_initial = np.array(full_res_initial)

            residue_restrained_xyz = []
            residue_all_xyz = []
            for atom in restrained_atoms:
                residue_restrained_xyz.append(np.array(residue.xyz(atom)))
            residue_restrained_xyz = np.array(residue_restrained_xyz)
            for index in range(1, residue.natoms() + 1):
                residue_all_xyz.append(np.array(residue.xyz(index)))
            residue_all_xyz = np.array(residue_all_xyz)
            residue_average = [sum(x)/len(x) for x in zip(*residue_all_xyz)]

            rotation, translation = \
                    get_transformation(residue_restrained_xyz, xyzarray,
                            average=residue_average)

            # The following transformation should be saved in rotamer_set,
            # so that's ultimately what we'll return.
            residue.apply_transform_Rx_plus_v(np_array_to_xyzM(rotation),np_array_to_xyzV(translation))
            # dump_invrot(residue, 'transformed_invrots/trans_invrot_' + str(i) + '.pdb')
            # i += 1


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
        # i = 0
        for rmsd in invrot_rmsds:
            # dump_invrot(invrot_rmsds[rmsd][0], 'inverse_rotamers/invrot_' + str(i) + '.pdb')
            if not best_rmsd or rmsd < best_rmsd:
                best_rmsd = rmsd
                best_invrot = invrot_rmsds[rmsd]
            # i += 1

        dump_invrot(best_invrot[0],'inverse_rotamers/best.pdb')
        return best_invrot, best_rmsd

    def choose_rotamer(self):

        invrot_rmsds = {}
        i = 0
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


def fast_relax(pose, residues_bb_movable, residues_sc_movable, selectors=True):
    '''Fast relax the pose'''
    if selectors==False:
        mm = setup_movemap(residues_bb_movable, residues_sc_movable)
    else:
        mm = setup_movemap_from_resselectors(residues_bb_movable,
                residues_sc_movable)
    sfxn = setup_restrained_sfxn(['coordinate_constraint'],[2.0])

    fast_relax_rounds = 5
    fast_relax = rosetta.protocols.relax.FastRelax(sfxn, fast_relax_rounds)
    fast_relax.set_movemap(mm) 
    
    fast_relax.apply(pose)
    pose.dump_file('out.pdb')



def fast_design(pose, designable_selector, repackable_selector,
        movemap=None, task_factory=None):
    '''Run fast design on the pose'''
    mm = setup_movemap_from_resselectors(designable_selector,
            repackable_selector)
    sfxn = setup_restrained_sfxn(['coordinate_constraint'],[1])


    fastdesign = rosetta.protocols.denovo_design.movers.FastDesign()
    fastdesign.set_movemap(mm)
    fastdesign.set_scorefxn(sfxn)
    fastdesign.set_up_default_task_factory()
    fastdesign.set_task_factory(task_factory)
    fastdesign.apply(pose)

    pose.dump_file('out.pdb')

def model_loops(pose, designable_selector, repackable_selector,
        focus_residue, movemap=None, task_factory=None, 
        mover='ngk', fast=False, resbuffer=3):
    '''Run loop modeler on the pose (default to NGK)
    Available movers: 
        - NGK
        - Loophash KIC'''
    mm = setup_movemap_from_resselectors(designable_selector,
            repackable_selector)
    sfxn = setup_restrained_sfxn(['coordinate_constraint'],[1.0])

    loopmodeler = rosetta.protocols.loop_modeler.LoopModeler()
    if mover=='ngk':
        loopmodeler.setup_kic_config()
        loops = generate_loops_from_res_selector(pose, designable_selector,
                focus_residue, resbuffer=resbuffer)
    elif mover=='lhk':
        '''A note on LHK: You need to mutate focus residues to their motif
        residue before running.'''
        assert(resbuffer >= 4)
        loopmodeler.setup_loophash_kic_config(True, str(focus_residue))
        loops = generate_loops_simple(pose, focus_residue, resbuffer)


    loopmodeler.set_loops(loops)
    #loopmodeler.set_cen_scorefxn(sfxn)
    loopmodeler.set_fa_scorefxn(sfxn)

    if fast:
        loopmodeler.centroid_stage().mark_as_test_run()
        loopmodeler.fullatom_stage().mark_as_test_run()

    if task_factory:
        loopmodeler.set_task_factory(task_factory)

    loopmodeler.apply(pose)
    pose.dump_file('out.pdb')



cst_test = ConstrainToInvRot()
rotamer_set = cst_test.create_inverse_rotamers('GLU')
cst_test.choose_rotamer()
cst_test.make_constraints_from_inverse_rotamer()

designable, repackable = choose_designable_residues(cst_test.pose, [38])
task_factory = setup_task_factory(cst_test.pose, designable, repackable,
        motif_dict={38:'E'},layered_design=False, prepare_focus=True)

bb_movable = [i for i in range(1,cst_test.pose.size() + 1)]
sc_movable = []
#fast_relax(cst_test.pose, bb_movable, sc_movable, selectors=False)
#print(cst_test.pose.constraint_set())
#print(designable)
#loops = generate_loops_from_res_selector(cst_test.pose, designable, 38)
#fast_design(cst_test.pose, designable, repackable, task_factory=task_factory)

model_loops(cst_test.pose, designable, repackable, 38,
        #task_factory=task_factory, 
        fast=False, mover='lhk', resbuffer=4)

