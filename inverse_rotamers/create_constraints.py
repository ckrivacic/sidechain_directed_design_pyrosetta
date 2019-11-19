'''
Utility file for creating a constraints file from a PDB and a residue type.
After downloading, gets coordinates from the motif atoms on the of interest.

Usage:
    python3 create_constraints.py <pdbid> <resnum> <chain> <comma_separated_atoms>


Downloads from https://files.rcsb.org/download/PDBID.pdb
'''

import sys, os, wget
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


create_constraints('8cho',{'38A':['N','CA','CB']})
