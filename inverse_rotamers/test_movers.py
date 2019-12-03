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


if __name__=='__main__':
    init('-ignore_unrecognized_res')
    wt_pdb = '3fa0'
    df = import_t4_dataframe('t4_inputs/t4_lysozyme.csv')
    pattern = re.compile('\w\d{1,3}\w')
    print(df)
    for i, row in df.iterrows():
        if len(row['mutant']) > 1:
            print(row['pdbid'])
            continue
        mut_pdb = row['pdbid']
        print('MUT PDB-+-----------------------------', mut_pdb)
        for focus_residue in row['mutant_dict'][0]:
            print('FOCUS RESIDUE-----------------', focus_residue)
            mut_pose, designable, repackable, task_factory =\
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
    pose, designable, repackable, task_factory = prepare_pdbid_for_modeling('4s0w','1cv1',111,'V')
    start_time = time.time()
    loopmodeler = get_loop_modeler(pose, designable, repackable, 111, task_factory=task_factory,
            fast=True, mover='ngk', resbuffer=4)
    loopmodeler.apply(pose)
    end_time = time.time() - start_time
    print('total time elapsed: ',end_time)
    pose.dump_file('out.pdb')
    '''
