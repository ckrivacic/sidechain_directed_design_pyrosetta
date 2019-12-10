from create_constraints import prepare_pdbid_for_modeling
from inverse_rotamers import *
import time, sys, re
import pandas as pd


def import_t4_dataframe(path):
    df = pd.read_csv(path)
    df['mutant_dict'] = [[] for x in range(len(df))]
    for i, row in df.iterrows():
        mutlist = []
        for mutation in row['mutant'].split(','):
            mutlist.append(mutation.strip())
        #df.set_value(i,'mutant',mutlist)
        df.at[i, 'mutant'] = mutlist
        resdict = {}
        for mutation in mutlist:
            wt, resnum, mut = re.split('([0-9]+)',mutation)
            wt = wt.upper()
            mut = mut.upper()
            resnum = int(resnum)
            resdict[resnum] = wt
        df.at[i, 'mutant_dict'] = [resdict]
    return df


def import_benchmark_dataframe(path):
    df = pd.read_csv(path)
    return df

if __name__=='__main__':
    df = import_benchmark_dataframe(sys.argv[1])
    print(df)
    '''
    Going to need to get all focus residues on their own line, so we can calc
    bb rmsd separately.
    '''
    task_num = 1
    init('-ignore_unrecognized_res')
    wt_pdb = '3fa0'
    #df = import_t4_dataframe('t4_inputs/t4_lysozyme.csv')
    #pattern = re.compile('\w\d{1,3}\w')
    row = df.loc[task_num]
    wt_pdbid = row['wt']
    mut_pdbid = row['mutant']
    out_dict = row.to_dict()
    focus_res = int(row['mut_res'])
    motif_dict = {focus_res:row['wt_restype']}
    designable, repackable, task_factory, aligner = \
            prepare_pdbid_for_modeling(wt_pdbid, mut_pdbid, motif_dict,
                    focus_res, int(row['wt_res']), constrain=row['constrain'])

    out_dict['pre_rmsd'] = aligner.bb_rmsd
    out_dict['pre_dist'] = distance_rosetta(aligner.target, focus_res,
            aligner.mobile, row['wt_res'])
    print(out_dict)
    '''
    for i, row in df.iterrows():
        if len(row['mutant']) > 1:
            print(row['pdbid'])
            continue
        mut_pdb = row['pdbid']
        print('MUT PDB-+-----------------------------', mut_pdb)
        for focus_residue in row['mutant_dict'][0]:
            print('FOCUS RESIDUE-----------------', focus_residue)
            mut_pose, designable, repackable, task_factory, bb_rmsd =\
                    prepare_pdbid_for_modeling(wt_pdb, mut_pdb,\
                            row['mutant_dict'][0], focus_residue,
                            prefix='t4_inputs')
            loopmodeler = get_loop_modeler(mut_pose, designable, repackable,
                    focus_residue, task_factory=task_factory, fast=True,
                    mover='ngk', resbuffer=4)
            start_time = time.time()
            loopmodeler.apply(mut_pose)
            elapsed = time.time() - start_time
            print('ELAPSED TIME = ',elapsed)
            mut_pose.dump_file('t4_outputs/ngk_' + mut_pdb + '_resi_' +\
                    str(focus_residue) + '.pdb')
    '''
