#$ -S /netapp/home/krivacic/.python36/bin/python3
#$ -l mem_free=4G
#$ -cwd
import sys
sys.path.insert(1,
        '/netapp/home/krivacic/intelligent_design/sidechain_directed_design_pyrosetta/inverse_rotamers/')
from benchmark_utils import *
from utils import distance_rosetta
from inverse_rotamers import *
import time, re, os, pickle
import pandas as pd


def import_benchmark_dataframe(path):
    df = pd.read_csv(path)
    return df


if __name__=='__main__':
    pdbredo_directory = '/netapp/home/krivacic/pdb_redo'
    shell=6
    task_num = int(os.environ['SGE_TASK_ID']) - 1
    #task_num = 1 # make sure to subtract 1 from SGE TASK ID for the real thing
    num_models = 250
    row_num = task_num//num_models

    mover = sys.argv[2]
    df = import_benchmark_dataframe(sys.argv[1])
    print(df)
    '''
    Going to need to get all focus residues on their own line, so we can calc
    bb rmsd separately.
    '''
    init('-ignore_unrecognized_res')
    row = df.loc[row_num]

    outdir = 'funnel_test/' + mover + '_' + str(row['constrain'])
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    default_sfxn = create_score_function('ref2015')
    wt_pdbid = row['wt']
    wt_pose = pose_from_pdbredo(wt_pdbid, prefix=pdbredo_directory)
    mut_pdbid = row['mutant']
    mut_pose = pose_from_pdbredo(mut_pdbid, prefix=pdbredo_directory)
    
    focus = Mismatch(int(row['mut_res']), int(row['wt_res']))
    mut_pair = MutantPair(mut_pose, wt_pose, [focus], shell=shell)

    out_dict = row.to_dict()
    ##focus_res = int(row['mut_res'])
    ##motif_dict = {focus_res:row['wt_restype']}

    designable, repackable = choose_designable_residues(mut_pose, [focus])
    task_factory = setup_task_factory(mut_pair.mut, mut_pair.wt,
            designable, repackable, motif_dict=mut_pair.motif_dict,
            layered_design=False, prepare_focus=True)

    if not os.path.exists(outdir + '/aligned/'):
        os.mkdir(outdir + '/aligned/')
    #aligner.mobile.dump_scored_pdb(outdir + '/aligned/' + wt_pdbid +
    #        '_' + str(task_num) + '.pdb', default_sfxn)
    #aligner.target.dump_scored_pdb(outdir + '/aligned/' + mut_pdbid +
    #        '_' + str(task_num) + '.pdb', default_sfxn)

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
        out_dict['post_score'] = default_sfxn(aligner.target)

        if not os.path.exists(outdir + '/ngk/'):
            os.mkdir(outdir + '/ngk/')

        #aligner.target.dump_scored_pdb(outdir + '/ngk/' + mut_pdbid +
        #        '_' + str(task_num) + '_ngk.pdb', default_sfxn)
        
        relaxer = fast_relax(aligner.target, designable, repackable,
                selectors=True)
        relaxer.apply(aligner.target)

        aligner.match_align()
        out_dict['post_rmsd_relaxed'] = aligner.bb_rmsd
        out_dict['post_dist_relaxed'] = distance_rosetta(aligner.target,
            focus_res, aligner.mobile, row['wt_res'])

        aligner.target.dump_scored_pdb(outdir + '/ngk/' + mut_pdbid +
                '_' + str(task_num) + '_relaxed.pdb', default_sfxn)
        out_dict['final_score'] = default_sfxn(aligner.target)
        print(default_sfxn(aligner.target))
        print(out_dict)
        with open(outdir + '/results_task_' + str(task_num), 'wb') as f:
            pickle.dump(out_dict, f)


    if mover == 'ngkf':
        loopmodeler = get_loop_modeler(aligner.target, designable, repackable,
                focus_res, task_factory=task_factory, fast=True,
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
        out_dict['post_score'] = default_sfxn(aligner.target)

        if not os.path.exists(outdir + '/ngk/'):
            os.mkdir(outdir + '/ngk/')

        #aligner.target.dump_scored_pdb(outdir + '/ngk/' + mut_pdbid +
        #        '_' + str(task_num) + '_ngk.pdb', default_sfxn)
        
        relaxer = fast_relax(aligner.target, designable, repackable,
                selectors=True)
        relaxer.apply(aligner.target)

        aligner.match_align()
        out_dict['post_rmsd_relaxed'] = aligner.bb_rmsd
        out_dict['post_dist_relaxed'] = distance_rosetta(aligner.target,
            focus_res, aligner.mobile, row['wt_res'])

        aligner.target.dump_scored_pdb(outdir + '/ngk/' + mut_pdbid +
                '_' + str(task_num) + '_relaxed.pdb', default_sfxn)
        out_dict['final_score'] = default_sfxn(aligner.target)
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
        out_dict['post_score'] = default_sfxn(aligner.target)

        if not os.path.exists(outdir + '/fastdesign/'):
            os.mkdir(outdir + '/fastdesign/')

        #aligner.target.dump_scored_pdb(outdir + '/fastdesign/' +
        #        mut_pdbid + '_' + str(task_num) + '_fastdesign.pdb',
        #        default_sfxn)

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
