from numeric import xyz_to_array
from scipy.spatial.distance import euclidean

def parse_restraints(path):
    restraints = []
    parsers = {
            'CoordinateConstraint': CoordinateConstraint,
            }

    with open(path) as file:
        for line in file:
            if not line.strip(): continue
            if line.startswith('#'): continue

            tokens = line.split()
            key, args = tokens[0], tokens[1:]

            if key not in parsers:
                raise IOError("Cannot parse '{0}' restraints.".format(key))

            restraint = parsers[key](args)
            restraints.append(restraint)

    return restraints


class CoordinateConstraint(object):

    def __init__(self, args):
        self.metric = 'dist'
        self.atom_name = args[0]
        self.atom_names = [args[0]]
        self.residue_id = int(args[1])
        self.residue_ids = [self.residue_id]
        self.atom = self.atom_name, self.residue_id
        self.coord = xyz_to_array(args[4:7])

    def distance_from_ideal(self, atom_xyzs):
        return euclidean(self.coord, atom_xyzs[self.atom])
