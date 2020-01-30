#! /wynton/home/kortemme/krivacic/software/anaconda3/bin/python3
##$ -S /wynton/home/kortemme/krivacic/software/python38/bin/python3
#$ -l mem_free=4G
#$ -cwd
import sys
sys.path.insert(1, "/wynton/home/kortemme/krivacic/intelligent_design/sidechain_directed_design_pyrosetta/inverse_rotamers/")
from benchmark_utils import *
from align import Mismatch
from utils import distance_rosetta
from inverse_rotamers import *
import time, re, os, pickle
import pandas as pd
import gzip
from copy import deepcopy

'''
Please put path to pdbredo as $PDBREDO in your bashrc
Usage:
    test_movers.py <input_df> <output_dir> <mover> <submission_number> [br]
'''
def smart_open(path, options):
    try:
        if path.endswith('.gz'): 
            return gzip.open(path, options)
        else: 
            return open(path, options)
    except:
        print('could not open file')

def annotate_pdb(pdb_path, out_dict):
    f = open(pdb_path, 'a')
    for key in out_dict:
        if type(out_dict[key]) != type('asdf'):
            f.write('\nmetric_' + key + ' ')
        else:
            f.write('\n' + key + ' ')
        f.write(str(out_dict[key]))
    f.close()

def import_benchmark_dataframe(path):
    df = pd.read_csv(path)
    return df

def import_backrub_dataframe(path):
    df = pd.read_csv(path, sep='\t')
    df.rename(columns={'pdb1':'wt','chain1':'wt_chain','resnum1':'wt_resnum','restype1':'wt_restype','pdb2':'mutant','chain2':'mut_chain','resnum2':'mut_resnum','restype2':'mut_restype'},
            inplace=True)
    return df


if __name__=='__main__':
    pdbredo_directory = '/wynton/home/kortemme/krivacic/pdb_redo'
    shell=6
    task_num = int(os.environ['SGE_TASK_ID']) - 1
    #task_num = 0 # make sure to subtract 1 from SGE TASK ID for the real thing
    num_models = 150

    mover = sys.argv[3]

    #offset = (int(sys.argv[4]) - 1) * num_models * 100
    #task_num += offset
    row_num = task_num//2

    # backrub variable tells us if we're reading from the backrub
    # benchmark dataframe, NOT whether we're using the backrub mover
    backrub = True if (len(sys.argv) > 4 and sys.argv[4] == 'br') else False
    if backrub:
        df = import_backrub_dataframe(sys.argv[1])
    else:
        df = import_benchmark_dataframe(sys.argv[1])
    print(df)
    '''
    Going to need to get all focus residues on their own line, so we can calc
    bb rmsd separately.
    '''
    init('-ignore_unrecognized_res -pdb_gz')
    row = df.loc[row_num]


    default_sfxn = create_score_function('ref2015')
    wt_pdbid = row['wt'].lower()
    mut_pdbid = row['mutant'].lower()

    constrain = 'constrained' if (task_num%2 == 0) else 'unconstrained'
    outdir_final = os.path.join(sys.argv[2], mut_pdbid + '_' + wt_pdbid, constrain)
    if not os.path.exists(outdir_final):
        os.makedirs(outdir_final, exist_ok=True)

    if 'TMPDIR' in os.environ:
        os_tmp = os.environ['TMPDIR']
    else:
        os_tmp = os.path.join('/scratch',os.environ['USER'])
    outdir_temp = os.path.join(os_tmp, str(task_num))
    if not os.path.exists(outdir_temp):
        os.makedirs(outdir_temp, exist_ok=True)

    # Try to get poses from pdbredo, otherwise download from rcsb.
    try:
        wt_pose = pose_from_pdbredo(wt_pdbid, pdbredo_directory)
    except:
        wt_pose = pose_from_rcsb(wt_pdbid, os_tmp)

    try:
        mut_pose = pose_from_pdbredo(mut_pdbid, pdbredo_directory)
    except:
        mut_pose = pose_from_rcsb(mut_pdbid, os_tmp)

    if backrub:
        wtfocus = wt_pose.pdb_info().pdb2pose(row['wt_chain'],
                row['wt_resnum'])
        mutfocus = mut_pose.pdb_info().pdb2pose(row['mut_chain'],
                row['mut_resnum'])
        focus = Mismatch(mutfocus, wtfocus)
    else:
        focus = Mismatch(int(row['mut_res']), int(row['wt_res']))
    #mut_pair = MutantPair(mut_pose, wt_pose, [focus], shell=shell)

    output_data = []

    for jobnum in range(0, num_models):

        try:
            out_dict = deepcopy(row.to_dict())
            ##focus_res = int(row['mut_res'])
            ##motif_dict = {focus_res:row['wt_restype']}
            if constrain == 'constrained':
                cst = True
            else:
                cst = False

            designable, repackable, task_factory, mut_pair = \
                    prepare_pdbids_for_modeling(wt_pdbid, mut_pdbid, [focus],
                            constrain=cst)

            out_dict['pre_rmsd'] = mut_pair.aligner.bb_rmsd
            out_dict['pre_dist'] = distance_rosetta(mut_pair.aligner.target,
                    focus.target,
                    mut_pair.aligner.mobile, focus.mobile)
            out_dict['mover'] = mover
            out_dict['decoy_number'] = jobnum


            if mover == 'ngk':
                modeler = get_loop_modeler(mut_pair.aligner.target, designable, repackable,
                        focus.target, task_factory=task_factory, fast=False,
                        mover='ngk', resbuffer=4)
            elif mover == 'ngkf':
                modeler = get_loop_modeler(mut_pair.aligner.target, designable, repackable,
                        focus.target, task_factory=task_factory, fast=True,
                        mover='ngk', resbuffer=4)
            elif mover == 'lhk':
                for key in mut_pair.motif_dict:
                    aatype = chemical_aa_from_oneletter(mut_pair.motif_dict[key])
                    mut_res = rosetta.protocols.simple_moves.MutateResidue(key,
                            aatype)
                    mut_res.apply(mut_pair.aligner.target)
                modeler = get_loop_modeler(mut_pair.aligner.target,
                        designable, repackable, focus.target,
                        task_factory=task_factory, fast=False,
                        mover='lhk', resbuffer=4)
            elif mover == 'fastdesign':
                modeler = fast_design(mut_pair.aligner.target, designable, repackable,
                        task_factory=task_factory)    
            elif mover == 'br':
                modeler = get_backrub_protocol(mut_pair.motif_dict,
                        shell=shell, kt=1.6, task_factory=task_factory,
                        ntrials=10000, stride=10000)


            start_time = time.time()
            modeler.apply(mut_pair.aligner.target)
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
            pdb_path = os.path.join(outdir_temp, wt_pdbid + '_' + mut_pdbid + '_' +
                    str(jobnum) + '_' + mover + '.pdb.gz')
            out_dict['path'] = pdb_path
            mut_pair.aligner.target.dump_scored_pdb(pdb_path, default_sfxn)
            
            relaxer = fast_relax(mut_pair.aligner.target, designable, repackable,
                    selectors=True)
            relaxer.apply(mut_pair.aligner.target)

            mut_pair.aligner.match_align()
            out_dict['post_rmsd_relaxed'] = mut_pair.aligner.bb_rmsd
            out_dict['post_dist_relaxed'] = distance_rosetta(mut_pair.aligner.target,
                focus.target, mut_pair.aligner.mobile, focus.mobile)

            pdb_path_rel = os.path.join(outdir_temp, wt_pdbid + '_' + mut_pdbid + '_' +
                    str(jobnum) + '_' + mover + '_relaxed.pdb.gz')
            out_dict['path_relaxed'] = pdb_path_rel
            mut_pair.aligner.target.dump_scored_pdb(pdb_path_rel, default_sfxn)
            out_dict['final_score'] = default_sfxn(mut_pair.aligner.target)
            print(default_sfxn(mut_pair.aligner.target))
            print(out_dict)
            #with open(outdir_temp + '/results_' + str(jobnum) + '.pkl', 'wb') as f:
            #    pickle.dump(out_dict, f)
            output_data.append(out_dict)

            annotate_pdb(pdb_path, out_dict)
            annotate_pdb(pdb_path_rel, out_dict)
        except:
            with open(os.path.join(outdir_temp, 'errors.txt'),'a') as f:
                f.write('Job number ' + str(jobnum) + ' failed \n')
    df_out = pd.DataFrame.from_records(output_data)
    print(df_out)
    with open(outdir_final + '/results.pkl', 'wb') as f:
        pickle.dump(df_out, f)
    finish_io(outdir_temp, outdir_final)
