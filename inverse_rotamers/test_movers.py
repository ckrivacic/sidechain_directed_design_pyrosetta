#$ -S /netapp/home/krivacic/.python36/bin/python3
#$ -l mem_free=4G
#$ -cwd
import sys
sys.path.insert(1,
        '/netapp/home/krivacic/intelligent_design/sidechain_directed_design_pyrosetta/inverse_rotamers/')
from create_constraints import prepare_pdbid_for_modeling
from utils import distance_rosetta
from inverse_rotamers import *
import time, re, os, pickle
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
    mover = sys.argv[2]
    df = import_benchmark_dataframe(sys.argv[1])
    print(df)
    '''
    Going to need to get all focus residues on their own line, so we can calc
    bb rmsd separately.
    '''
    #task_num = 1 # make sure to subtract 1 from SGE TASK ID for the real thing
    task_num = int(os.environ['SGE_TASK_ID']) - 1
    num_models = 20
    row_num = task_num//num_models
    init('-ignore_unrecognized_res')
    #wt_pdb = '3fa0'
    #df = import_t4_dataframe('t4_inputs/t4_lysozyme.csv')
    #pattern = re.compile('\w\d{1,3}\w')
    row = df.loc[row_num]

    outdir = mover + str(row['constrain'])
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    default_sfxn = create_score_function('ref2015')
    wt_pdbid = row['wt']
    mut_pdbid = row['mutant']
    out_dict = row.to_dict()
    focus_res = int(row['mut_res'])
    motif_dict = {focus_res:row['wt_restype']}

    designable, repackable, task_factory, aligner = \
            prepare_pdbid_for_modeling(wt_pdbid, mut_pdbid, motif_dict,
                    focus_res, int(row['wt_res']), constrain=row['constrain'])

    if not os.path.exists(outdir + '/aligned/'):
        os.mkdir(outdir + '/aligned/')
    aligner.mobile.dump_scored_pdb(outdir + '/aligned/' + wt_pdbid +
            '_' + str(task_num) + '.pdb', default_sfxn)
    aligner.target.dump_scored_pdb(outdir + '/aligned/' + mut_pdbid +
            '_' + str(task_num) + '.pdb', default_sfxn)

    out_dict['pre_rmsd'] = aligner.bb_rmsd
    out_dict['pre_dist'] = distance_rosetta(aligner.target, focus_res,
            aligner.mobile, row['wt_res'])
    print(out_dict)

    if mover == 'ngk':
        loopmodeler = get_loop_modeler(aligner.target, designable, repackable,
                focus_res, task_factory=task_factory, fast=False,
                mover='ngk', resbuffer=4)
        start_time = time.time()
        loopmodeler.apply(aligner.target)
        elapsed = time.time() - start_time
        out_dict['elapsed_time'] = elapsed
        aligner.target.remove_constraints()

        aligner.match_align()
        out_dict['post_rmsd'] = aligner.bb_rmsd
        out_dict['post_dist'] = distance_rosetta(aligner.target,
            focus_res, aligner.mobile, row['wt_res'])

        if not os.path.exists(outdir + '/ngk_fast/'):
            os.mkdir(outdir + '/ngk/')

        aligner.target.dump_scored_pdb(outdir + '/ngk/' + mut_pdbid +
                '_' + str(task_num) + '_ngk.pdb', default_sfxn)
        
        relaxer = fast_relax(aligner.target, designable, repackable,
                selectors=True)
        relaxer.apply(aligner.target)

        aligner.match_align()
        out_dict['post_rmsd_relaxed'] = aligner.bb_rmsd
        out_dict['post_dist_relaxed'] = distance_rosetta(aligner.target,
            focus_res, aligner.mobile, row['wt_res'])

        aligner.target.dump_scored_pdb(outdir + '/ngk/' + mut_pdbid +
                '_' + str(task_num) + '_relaxed.pdb', default_sfxn)
        print(default_sfxn(aligner.target))
        print(out_dict)
        with open(outdir + '/results_task_' + str(task_num), 'wb') as f:
            pickle.dump(out_dict, f)


    elif mover == 'fastdesign':
        fastdesign = fast_design(aligner.target, designable, repackable,
                task_factory=task_factory)    
        start_time = time.time()
        fastdesign.apply(aligner.target)
        elapsed = time.time() - start_time
        out_dict['elapsed_time'] = elapsed
        aligner.target.remove_constraints()

        aligner.match_align()
        out_dict['post_rmsd'] = aligner.bb_rmsd
        out_dict['post_dist'] = distance_rosetta(aligner.target, 
                focus_res, aligner.mobile, row['wt_res'])

        if not os.path.exists(outdir + '/fastdesign/'):
            os.mkdir(outdir + '/fastdesign/')

        aligner.target.dump_scored_pdb(outdir + '/fastdesign/' +
                mut_pdbid + '_' + str(task_num) + '_fastdesign.pdb',
                default_sfxn)

        relaxer = fast_relax(aligner.target, designable, repackable,
                selectors=True)
        relaxer.apply(aligner.target)

        aligner.match_align()
        out_dict['post_rmsd_relaxed'] = aligner.bb_rmsd
        out_dict['post_dist_relaxed'] = distance_rosetta(aligner.target,
            focus_res, aligner.mobile, row['wt_res'])

        aligner.target.dump_scored_pdb(outdir + '/fastdesign/' +
                mut_pdbid + '_' + str(task_num) + '_fastdesign_relaxed.pdb',
                default_sfxn)
        out_dict['final_score'] = default_sfxn(aligner.target)
        print(default_sfxn(aligner.target))
        print(out_dict)
        with open(outdir + '/results_task_' + str(task_num), 'wb') as f:
            pickle.dump(out_dict, f)
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
