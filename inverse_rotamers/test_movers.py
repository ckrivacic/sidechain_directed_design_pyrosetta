from create_constraints import prepare_pdbid_for_modeling
from inverse_rotamers import *
import time, sys, re
import pandas as pd


def import_t4_dataframe(path):
    df = pd.read_csv(path)
    for i, row in df.iterrows():
        mutlist = []
        for mutation in row['mutant'].split(','):
            mutlist.append(mutation.strip())
        df.set_value(i,'mutant',mutlist)
    return df


if __name__=='__main__':
    print('hi')
    wt_pdb = '3fa0'
    df = import_t4_dataframe('t4_inputs/t4_lysozyme.csv')
    print(df)
    pattern = re.compile('\w\d{1,3}\w')
    for i, row in df.iterrows():
        mut_pdb = row['pdbid']
        for mutation in row['mutant']:
            wt, resnum, mut = re.split('([0-9]+)',mutation)


    '''
    init()
    pose, designable, repackable, task_factory = prepare_pdbid_for_modeling('4s0w','1cv1',111,'V')
    start_time = time.time()
    loopmodeler = get_loop_modeler(pose, designable, repackable, 111, task_factory=task_factory,
            fast=True, mover='ngk', resbuffer=4)
    loopmodeler.apply(pose)
    end_time = time.time() - start_time
    print('total time elapsed: ',end_time)
    pose.dump_file('out.pdb')
    '''
