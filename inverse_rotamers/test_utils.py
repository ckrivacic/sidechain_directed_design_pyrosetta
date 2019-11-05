from numeric import xyz_from_3d_array
from pyrosetta import *

def plot_3d(mobile,target,mobile_moved):
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt

    fig = plt.Figure()
    ax = plt.axes(projection='3d')
    ax.scatter(*xyz_from_3d_array(mobile),c='blue')
    ax.scatter(*xyz_from_3d_array(target),c='green')
    ax.scatter(*xyz_from_3d_array(mobile_moved),c='red')
    plt.show()


def dump_invrot(inverse_rotamer, filename):
    dummy_pose = rosetta.core.pose.Pose()
    dummy_pose.append_residue_by_jump(inverse_rotamer.clone(), 1)

    dummy_pose.dump_file(filename)
