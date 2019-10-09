import numpy as np
from pyrosetta import *

class Transformation(object):
    """Class for storing rotation & transformation information"""

    def __init__(self,rotation,translation):
        self.rotation = rotation
        self.translation = translation


# I will need a function to calculate the rotation and translation
# vectors to align my inverse rotamers to target coordinates. Let's
# start with that, then we'll write the functions to figure out which
# atoms go where.
def get_transformation(template_coordinate_set,
        target_coordinate_set):

    """Function that returns rotation and translation vectors to align
    two vectors of equal dimensions"""

    # coordinate sets assumed to be of the form
    # [[x1,y1,z1],[x2,y2,z2],...[xn,yn,zn]]
    template_average = [sum(x)/len(x) for x in
            zip(*template_coordinate_set)]
    target_average = [sum(x)/len(x) for x in zip(*target_coordinate_set)]

    template_zeroed = template_coordinate_set - template_average
    target_zeroed = target_coordinate_set - target_average
    matrix = np.dot(template_zeroed.T,target_zeroed)

    U, s, Vh = np.linalg.svd(matrix)
    Id = np.array([[1,0,0],
                  [0,1,0],
                  [0,0,np.sign(np.linalg.det(matrix))]])

    rotation = np.dot(Vh.T, np.dot(Id, U.T))
    translation = target_average - np.dot(template_average,rotation.T)

    return rotation, translation


def apply_transformation(Transformation, template_coordinate_set):
    return np.dot(template_coordinate_set, Transformation.rotation.T) +\
            Transformation.translation


def simple_align(mobile, target):
    # Handles getting and applying transformations so I don't have to
    # keep track of which one is template and which one is target
    # through multiple functions
    t = Transformation(*get_transformation(mobile, target))
    return apply_transformation(t, mobile)


def bb_independent_rotamers_extra_chi(rot_restype):
    """
    The bb_independent_rotamers function in Rosetta ignores -ex1 and
    -ex2 flags in the command line, and naively generates rotamers only
    for phi 180 and psi 180.
    This function tries to generate a complete rotamer set for a residue
    type.
    """
    firstres = rosetta.core.conformation.Residue(rot_restype, True)
    dummy_pose = rosetta.core.pose.Pose()
    dummy_pose.append_residue_by_jump(firstres, 0)
    
    if rot_restype.is_polymer():
        rosetta.core.pose.add_lower_terminus_type_to_pose_residue(dummy_pose,
                1)
        rosetta.core.pose.add_upper_terminus_type_to_pose_residue(dummy_pose,
                1)
    
    
    # Set up sfxn and dummy task
    dummy_sfxn = rosetta.core.scoring.ScoreFunction()
    dummy_sfxn(dummy_pose)
    dummy_task = rosetta.core.pack.task.TaskFactory.create_packer_task(dummy_pose)

    # options for dummy task
    dummy_task.initialize_from_command_line()
    dummy_task.initialize_extra_rotamer_flags_from_command_line()
    dummy_task.show_residue_task(1)
    dummy_task.nonconst_residue_task(1).restrict_to_repacking()
    dummy_task.nonconst_residue_task(1).or_include_current(False)
    dummy_task.nonconst_residue_task(1).or_fix_his_tautomer(True)
    dummy_png = rosetta.core.pack.create_packer_graph(dummy_pose,\
            dummy_sfxn, dummy_task)

    rotset = rosetta.core.pack.rotamer_set.RotamerSetFactory.create_rotamer_set(dummy_pose)
    rotset.set_resid(1)
    rotset.build_rotamers(dummy_pose, dummy_sfxn, dummy_task, dummy_png)

    to_return = []
    for i in range(1, rotset.num_rotamers() + 1):
        print(i)


def test_rotamer_gen(resn):
    dummy_pose = rosetta.core.pose.Pose()
    global_residue_type_set = dummy_pose.residue_type_set_for_pose()
    restype = global_residue_type_set.get_representative_type_base_name(resn)
    bb_independent_rotamers_extra_chi(restype)

init('-ex1 -ex2 -ex3 -ex4')
test_rotamer_gen("ASP")

