import numpy as np
import math
import pyrosetta
from pyrosetta import rosetta


def xyzV_to_np_array(xyz):
    return np.array([xyz.x, xyz.y, xyz.z])

def np_array_to_xyzV(a):
    return rosetta.numeric.xyzVector_double_t(a[0], a[1], a[2])

def intlist_to_vector1_size(a):
    vector = rosetta.utility.vector1_unsigned_long()
    for item in a:
        vector.append(item)
    return vector

def xyzM_to_np_array(M):
    return np.array([[M.xx, M.xy, M.xz],
                     [M.yx, M.yy, M.yz],
                     [M.zx, M.zy, M.zz]])

def np_array_to_xyzM(a):
    return rosetta.numeric.xyzMatrix_double_t.rows(
            a[0][0], a[0][1], a[0][2],
            a[1][0], a[1][1], a[1][2],
            a[2][0], a[2][1], a[2][2])

def mult_np_transformation(T1, T2):
    '''Multiply two numpy rigid body transformations'''
    M1, v1 = T1
    M2, v2 = T2
    
    return np.dot(M1, M2), np.dot(M1, v2) + v1

def inverse_np_transformation(T):
    '''Inverse an numpy rigid body transformation.'''
    M, v = T
    
    invM = np.linalg.inv(M)   
    return invM, - np.dot(invM, v)


def RMSD(points1, poinsts2):
    '''Calcualte RMSD between two lists of numpy points.'''
    diff = [points1[i] - poinsts2[i] for i in range(len(points1))]
    return np.sqrt(sum(np.dot(d, d) for d in diff) / len(diff))


def xyz_from_3d_array(array):
    """
    Takes a 3-dimensional numpy array and returns lists of x, y, and z
    coordinates.
    """
    x = array[:,0]
    y = array[:,1]
    z = array[:,2]

    return x,y,z


def xyz_to_array(xyz):
    """
    Convert a list of strings representing a 3D coordinate to floats and return
    the coordinate as a ``numpy`` array.
    """
    return np.array([float(x) for x in xyz])


def euclidean_distance(xyz1, xyz2):
    """
    Simple function for calculating euclidean distance between two points.
    """
    dist = [(a - b)**2 for a,b in zip(xyz1, xyz2)]
    return math.sqrt(sum(dist))


def backbone_rmsd(rotamer, residue,
        alignment_atoms):
    """
    Measure backbone RMSD between a rotamer and the nearest residue on
    the design protein.
    """

    distances = np.array([])
    for atom in alignment_atoms:
        rot_xyz = xyzV_to_np_array(rotamer.xyz(atom))
        near_xyz = xyzV_to_np_array(residue.xyz(atom))
        distances = np.append(distances,euclidean_distance(rot_xyz,near_xyz))

    return np.sqrt((distances**2).mean())
