#$ -S /netapp/home/krivacic/.python36/bin/python3
#$ -l mem_free=4G
#$ -cwd
import sys
sys.path.insert(1,
        '/netapp/home/krivacic/intelligent_design/sidechain_directed_design_pyrosetta/inverse_rotamers/')
from benchmark_utils import *
from align import Mismatch
from utils import distance_rosetta
from inverse_rotamers import *
import time, re, os, pickle
import pandas as pd

'''
Please put path to pdbredo as $PDBREDO in your bashrc
Usage:
    test_movers.py <input_df> <output_dir> <mover> <submission_number> [br]
'''

def import_benchmark_dataframe(path):
    df = pd.read_csv(path)
    return df

def import_backrub_dataframe(path):
    df = pd.read_csv(path, sep='\t')
    df.rename(columns={'pdb1':'wt','chain1':'wt_chain','resnum1':'wt_resnum','restype1':'wt_restype','pdb2':'mutant','chain2':'mut_chain','resnum2':'mut_resnum','restype2':'mut_restype'},
            inplace=True)
    return df


if __name__=='__main__':
    pdbredo_directory = '/netapp/home/krivacic/pdb_redo'
    shell=6
    task_num = int(os.environ['SGE_TASK_ID']) - 0
    #task_num = 0 # make sure to subtract 1 from SGE TASK ID for the real thing
    num_models = 250
    row_num = task_num//num_models

    mover = sys.argv[3]

    offset = (int(sys.argv[4]) - 1) * num_models * 100
    task_num += offset

    # backrub variable tells us if we're reading from the backrub
    # benchmark dataframe, NOT whether we're using the backrub mover
    backrub = True if (len(sys.argv) > 5 and sys.argv[5] == 'br') else False
    if backrub:
        df = import_backrub_dataframe(sys.argv[1])
    else:
        df = import_benchmark_dataframe(sys.argv[1])
    print(df)
    '''
    Going to need to get all focus residues on their own line, so we can calc
    bb rmsd separately.
    '''
    init('-ignore_unrecognized_res')
    row = df.loc[row_num]


    default_sfxn = create_score_function('ref2015')
    wt_pdbid = row['wt'].lower()
    mut_pdbid = row['mutant'].lower()

    constrain = 'constrained' if (task_num%2 == 0) else 'unconstrained'
    outdir = os.path.join(sys.argv[2], wt_pdbid + '_' + mut_pdbid, constrain)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    wt_pose = pose_from_pdbredo(wt_pdbid, pdbredo_directory)
    mut_pose = pose_from_pdbredo(mut_pdbid, pdbredo_directory)

    if backrub:
        wtfocus = wt_pose.pdb_info().pdb2pose(row['wt_chain'],
                row['wt_resnum'])
        mutfocus = mut_pose.pdb_info().pdb2pose(row['mut_chain'],
                row['mut_resnum'])
        focus = Mismatch(mutfocus, wtfocus)
    else:
        focus = Mismatch(int(row['mut_res']), int(row['wt_res']))
    #mut_pair = MutantPair(mut_pose, wt_pose, [focus], shell=shell)

    out_dict = row.to_dict()
    ##focus_res = int(row['mut_res'])
    ##motif_dict = {focus_res:row['wt_restype']}
    if constrain == 'constrained':
        cst = True
    else:
        cst = False

    designable, repackable, task_factory, mut_pair = \
            prepare_pdbids_for_modeling(wt_pdbid, mut_pdbid, [focus],
                    constrain=cst)


    #if not os.path.exists(outdir + '/aligned/'):
    #    os.mkdir(outdir + '/aligned/')
    ##aligner.mobile.dump_scored_pdb(outdir + '/aligned/' + wt_pdbid +
    #        '_' + str(task_num) + '.pdb', default_sfxn)
    #aligner.target.dump_scored_pdb(outdir + '/aligned/' + mut_pdbid +
    #        '_' + str(task_num) + '.pdb', default_sfxn)

    out_dict['pre_rmsd'] = mut_pair.aligner.bb_rmsd
    out_dict['pre_dist'] = distance_rosetta(mut_pair.aligner.target,
            focus.target,
            mut_pair.aligner.mobile, focus.mobile)
    print(out_dict)


    if mover == 'ngk':
        loopmodeler = get_loop_modeler(mut_pair.aligner.target, designable, repackable,
                focus.target, task_factory=task_factory, fast=False,
                mover='ngk', resbuffer=4)
        start_time = time.time()
        loopmodeler.apply(mut_pair.aligner.target)
        elapsed = time.time() - start_time
        out_dict['elapsed_time'] = elapsed
        mut_pair.aligner.target.remove_constraints()

        mut_pair.aligner.match_align()
        out_dict['post_rmsd'] = mut_pair.aligner.bb_rmsd
        out_dict['post_dist'] = distance_rosetta(mut_pair.aligner.target,
            focus.target, mut_pair.aligner.mobile, focus.mobile)
        out_dict['post_score'] = default_sfxn(mut_pair.aligner.target)

        #aligner.target.dump_scored_pdb(outdir + '/ngk/' + mut_pdbid +
        #        '_' + str(task_num) + '_ngk.pdb', default_sfxn)
        
        relaxer = fast_relax(mut_pair.aligner.target, designable, repackable,
                selectors=True)
        relaxer.apply(mut_pair.aligner.target)

        mut_pair.aligner.match_align()
        out_dict['post_rmsd_relaxed'] = mut_pair.aligner.bb_rmsd
        out_dict['post_dist_relaxed'] = distance_rosetta(mut_pair.aligner.target,
            focus.target, mut_pair.aligner.mobile, focus.mobile)

        mut_pair.aligner.target.dump_scored_pdb(outdir + '/' + wt_pdbid
                + '_' + mut_pdbid +
                '_' + str(task_num) + '_relaxed.pdb', default_sfxn)
        out_dict['final_score'] = default_sfxn(mut_pair.aligner.target)
        print(default_sfxn(mut_pair.aligner.target))
        print(out_dict)
        with open(outdir + '/results_task_' + str(task_num), 'wb') as f:
            pickle.dump(out_dict, f)


    if mover == 'ngkf':
        loopmodeler = get_loop_modeler(mut_pair.aligner.target, designable, repackable,
                focus.target, task_factory=task_factory, fast=True,
                mover='ngk', resbuffer=4)
        start_time = time.time()
        loopmodeler.apply(mut_pair.aligner.target)
        elapsed = time.time() - start_time
        out_dict['elapsed_time'] = elapsed
        mut_pair.aligner.target.remove_constraints()

        mut_pair.aligner.match_align()
        out_dict['post_rmsd'] = mut_pair.aligner.bb_rmsd
        out_dict['post_dist'] = distance_rosetta(mut_pair.aligner.target,
            focus.target, mut_pair.aligner.mobile, focus.mobile)
        out_dict['post_score'] = default_sfxn(mut_pair.aligner.target)

        #aligner.target.dump_scored_pdb(outdir + '/ngk/' + mut_pdbid +
        #        '_' + str(task_num) + '_ngk.pdb', default_sfxn)
        
        relaxer = fast_relax(mut_pair.aligner.target, designable, repackable,
                selectors=True)
        relaxer.apply(mut_pair.aligner.target)

        mut_pair.aligner.match_align()
        out_dict['post_rmsd_relaxed'] = mut_pair.aligner.bb_rmsd
        out_dict['post_dist_relaxed'] = distance_rosetta(mut_pair.aligner.target,
            focus.target, mut_pair.aligner.mobile, focus.mobile)

        mut_pair.aligner.target.dump_scored_pdb(outdir +'/' + wt_pdbid +
                '_' + mut_pdbid +
                '_' + str(task_num) + '_relaxed.pdb', default_sfxn)
        out_dict['final_score'] = default_sfxn(mut_pair.aligner.target)
        print(default_sfxn(mut_pair.aligner.target))
        print(out_dict)
        with open(outdir + '/results_task_' + str(task_num), 'wb') as f:
            pickle.dump(out_dict, f)


    elif mover == 'fastdesign':
        fastdesign = fast_design(mut_pair.aligner.target, designable, repackable,
                task_factory=task_factory)    
        start_time = time.time()
        fastdesign.apply(mut_pair.aligner.target)
        elapsed = time.time() - start_time
        out_dict['elapsed_time'] = elapsed
        mut_pair.aligner.target.remove_constraints()

        mut_pair.aligner.match_align()
        out_dict['post_rmsd'] = mut_pair.aligner.bb_rmsd
        out_dict['post_dist'] = distance_rosetta(mut_pair.aligner.target, 
            focus.target, mut_pair.aligner.mobile, focus.mobile)
        out_dict['post_score'] = default_sfxn(mut_pair.aligner.target)

        ##aligner.target.dump_scored_pdb(outdir + '/fastdesign/' +
        #        mut_pdbid + '_' + str(task_num) + '_fastdesign.pdb',
        #        default_sfxn)

        relaxer = fast_relax(mut_pair.aligner.target, designable, repackable,
                selectors=True)
        relaxer.apply(mut_pair.aligner.target)

        mut_pair.aligner.match_align()
        out_dict['post_rmsd_relaxed'] = mut_pair.aligner.bb_rmsd
        out_dict['post_dist_relaxed'] = distance_rosetta(mut_pair.aligner.target,
            focus.target, mut_pair.aligner.mobile, focus.mobile)

        mut_pair.aligner.target.dump_scored_pdb(outdir + '/' +
                wt_pdbid + '_' + mut_pdbid + '_' + str(task_num) + '_fastdesign_relaxed.pdb',
                default_sfxn)
        out_dict['final_score'] = default_sfxn(mut_pair.aligner.target)
        print(default_sfxn(mut_pair.aligner.target))
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
