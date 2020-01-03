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


def plot(list_of_tuples, xlabel='x', ylabel='y', title='title', groups=None,
        zorders=None, vline=None, unitline=False):
    '''Plotting function. Input should be a list of tuples (x, y). For ex:
    [(x1, y1), (x2, y2)]'''

    from matplotlib import pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    if groups is not None and zorders is not None:
        for data, group, order in zip(list_of_tuples, groups, zorders):
            x, y = data
            ax.scatter(x, y, label=group, zorder=order)

    elif groups is not None:
        for data, group in zip(list_of_tuples, groups):
            x, y = data
            ax.scatter(x, y, label=group)

    elif groups is None and zorders is None:
        for data in list_of_tuples:
            x, y = data
            ax.scatter(x, y)

    if unitline:
        ax.plot(x, x)
    if vline:
        ax.axvline(vline)

    plt.legend(loc=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    #plt.show()

    return plt, fig, ax
