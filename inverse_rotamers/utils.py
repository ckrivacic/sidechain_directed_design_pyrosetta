import numpy as np
from numeric import *
from pyrosetta import *

class Transformation(object):
    """Class for storing rotation & transformation information"""

    def __init__(self,rotation,translation):
        self.rotation = rotation
        self.translation = translation


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


# I will need a function to calculate the rotation and translation
# vectors to align my inverse rotamers to target coordinates. Let's
# start with that, then we'll write the functions to figure out which
# atoms go where.
def get_transformation(template_coordinate_set,
        target_coordinate_set, average=None):

    """
    DEPRECIATED; XingJie's function above works w/o the 'average' variable.

    Function that returns rotation and translation vectors to align
    two vectors of equal dimensions. Can give this function an average to use
    for calculating rotation matrix. 
    The purpose of the "average" variable is to rotate the object around a
    center that is outside the coordinate sets. For instance if aligning a part
    of a mobile residue onto part of a template residue, 'average' should be the center
    of the mobile residue."""

    # coordinate sets assumed to be of the form
    # [[x1,y1,z1],[x2,y2,z2],...[xn,yn,zn]]
    template_average = [sum(x)/len(x) for x in
            zip(*template_coordinate_set)]
    # May need a different average to calculate rotation if aligning a subset
    # of something.
    if average:
        rotation_average = average
    else:
        rotation_average = template_average

    target_average = [sum(x)/len(x) for x in zip(*target_coordinate_set)]

    template_zeroed = template_coordinate_set - rotation_average
    target_zeroed = target_coordinate_set - target_average
    matrix = np.dot(template_zeroed.T,target_zeroed)

    U, s, Vh = np.linalg.svd(matrix)
    if (np.linalg.det(U) * np.linalg.det(Vh)) < 0.0:
        U[:,-1] = -U[:,-1]
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


def list_to_str(l):
    return ','.join(list(str(i) for i in l))


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
    # dummy_task.show_residue_task(1)
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
        rot = firstres.clone()
        for j in range(1, firstres.nchi() + 1):
            rot.set_chi(j, rotset.rotamer(i).chi(j))
            dummy_pose.set_chi(j, 1, rotset.rotamer(i).chi(j))
            #dummy_pose.dump_file('test_outputs/testout_' + str(i) + '.pdb')
        to_return.append(rot)
    return to_return#, dummy_pose


def bb_independent_rotamers_extra_chi_phipsi(rot_restype):
    """
    The bb_independent_rotamers function in Rosetta ignores -ex1 and
    -ex2 flags in the command line, and naively generates rotamers only
    for phi 180 and psi 180.
    This function tries to generate a complete rotamer set for a residue
    type.
    """
    firstres = rosetta.core.conformation.Residue(rot_restype, True)
    dummy_pose = rosetta.core.pose.Pose()
    global_residue_type_set = dummy_pose.residue_type_set_for_pose()
    gly = global_residue_type_set.get_representative_type_base_name('GLY')
    glyres = rosetta.core.conformation.Residue(gly, True)
    dummy_pose.append_residue_by_jump(glyres, 0)
    dummy_pose.append_residue_by_bond(firstres,
            build_ideal_geometry=True)
    dummy_pose.append_residue_by_bond(glyres, build_ideal_geometry=True)

    to_return = []
    for psi in np.linspace(-180., 180., 40):
        for phi in np.linspace(-180., 180., 40):

            #if rot_restype.is_polymer():
            #    rosetta.core.pose.add_lower_terminus_type_to_pose_residue(dummy_pose,
            #            1)
            #    rosetta.core.pose.add_upper_terminus_type_to_pose_residue(dummy_pose,
            #            1)

            dummy_pose.set_psi(2, psi)
            dummy_pose.set_phi(2, phi)
            # Set up sfxn and dummy task
            dummy_sfxn = rosetta.core.scoring.ScoreFunction()
            dummy_sfxn(dummy_pose)
            dummy_task = rosetta.core.pack.task.TaskFactory.create_packer_task(dummy_pose)

            # options for dummy task
            dummy_task.initialize_from_command_line()
            dummy_task.initialize_extra_rotamer_flags_from_command_line()
            # dummy_task.show_residue_task(1)
            dummy_task.nonconst_residue_task(2).restrict_to_repacking()
            dummy_task.nonconst_residue_task(2).or_include_current(False)
            dummy_task.nonconst_residue_task(2).or_fix_his_tautomer(True)
            dummy_png = rosetta.core.pack.create_packer_graph(dummy_pose,\
                    dummy_sfxn, dummy_task)

            rotset = rosetta.core.pack.rotamer_set.RotamerSetFactory.create_rotamer_set(dummy_pose)
            rotset.set_resid(2)
            rotset.build_rotamers(dummy_pose, dummy_sfxn, dummy_task, dummy_png)

            for i in range(1, rotset.num_rotamers() + 1):
                rot = firstres.clone()
                for j in range(1, firstres.nchi() + 1):
                    rot.set_chi(j, rotset.rotamer(i).chi(j))
                    dummy_pose.set_chi(j, 2, rotset.rotamer(i).chi(j))
                    dummy_pose.dump_file('test_outputs/testout_phi' +
                            str(phi) + '_psi' + str(psi) + '_' + str(i) + '.pdb')
                to_return.append(rot)
    return to_return, dummy_pose


def apply_rotamer(pose, rotamer):
    pose.set_torsion(rotamer.id, rotamer.value)


def test_rotamer_gen(resn):
    dummy_pose = rosetta.core.pose.Pose()
    global_residue_type_set = dummy_pose.residue_type_set_for_pose()
    restype = global_residue_type_set.get_representative_type_base_name(resn)
    bb_independent_rotamers_extra_chi(restype)


def make_constraints_from_inverse_rotamer(inverse_rotamer, resid, pose,
        atoms=['N','C','CA']):
    func = rosetta.core.scoring.constraints.BoundFunc(0, 0.05, 0.4,
            "invrot")
    coordinate_constraints = []
    for atom in atoms:
        xyzV = inverse_rotamer.xyz(atom)
        fixed_pt = pose.atom_tree().root().atom_id()
        atomno = pose.residue(resid).atom_index(atom)
        atom_id = rosetta.core.id.AtomID(atomno, resid)
        coordinate_constraint = \
        rosetta.core.scoring.constraints.CoordinateConstraint(
                atom_id, fixed_pt, xyzV, func
                )
        coordinate_constraints.append(coordinate_constraint)

    return coordinate_constraints


def res_selector_to_size_list(resselector):
    size_list = []
    for i, boolean in enumerate(resselector):
        if boolean == True:
            size_list.append(int(i + 1))

    return intlist_to_vector1_size(size_list)


def distance_pdb(pose1, res1_pdb, res1_chain, pose2, res2_pdb, res2_chain):
    xyz1 = pose1.residue(pose1.pdb_info().pdb2pose(res1_chain, res1_pdb)).xyz('CA')
    xyz2 = pose2.residue(pose2.pdb_info().pdb2pose(res2_chain, res2_pdb)).xyz('CA')
    return euclidean_distance(xyz1, xyz2)


def distance_rosetta(pose1, res1, pose2, res2):
    xyz1 = pose1.residue(res1).xyz('CA')
    xyz2 = pose2.residue(res2).xyz('CA')
    return euclidean_distance(xyz1, xyz2)


def oneletter_to_threeletter(aa):
    one_to_three = {
            'a':'ala', 'r':'arg', 'n':'asn', 'd':'asp', 'c':'cys', 'e':'glu',
            'q':'gln', 'g':'gly', 'h':'his', 'i':'ile', 'l':'leu', 'k':'lys',
            'm':'met', 'f':'phe', 'p':'pro', 's':'ser', 't':'thr', 'w':'trp',
            'y':'tyr', 'v':'val'
            }
    return one_to_three[aa.lower()]


def chemical_aa_from_oneletter(aa):
    one_to_three = {
            'a':rosetta.core.chemical.AA.aa_ala,
            'r':rosetta.core.chemical.AA.aa_arg,
            'n':rosetta.core.chemical.AA.aa_asn,
            'd':rosetta.core.chemical.AA.aa_asp,
            'c':rosetta.core.chemical.AA.aa_cys,
            'e':rosetta.core.chemical.AA.aa_glu,
            'q':rosetta.core.chemical.AA.aa_gln,
            'g':rosetta.core.chemical.AA.aa_gly,
            'h':rosetta.core.chemical.AA.aa_his,
            'i':rosetta.core.chemical.AA.aa_ile,
            'l':rosetta.core.chemical.AA.aa_leu,
            'k':rosetta.core.chemical.AA.aa_lys,
            'm':rosetta.core.chemical.AA.aa_met,
            'f':rosetta.core.chemical.AA.aa_phe,
            'p':rosetta.core.chemical.AA.aa_pro,
            's':rosetta.core.chemical.AA.aa_ser,
            't':rosetta.core.chemical.AA.aa_thr,
            'w':rosetta.core.chemical.AA.aa_trp,
            'y':rosetta.core.chemical.AA.aa_tyr,
            'v':rosetta.core.chemical.AA.aa_val
            }
    return one_to_three[aa.lower()]
#init('-extrachi_cutoff 0 -ex1 -ex2')
#test_rotamer_gen("ASP")
