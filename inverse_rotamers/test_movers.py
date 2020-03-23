#! /wynton/home/kortemme/krivacic/software/anaconda36/bin/python3
#$ -l mem_free=4G
#$ -cwd
"""
Usage: test_movers <input_df> <output_dir> <mover> [options]

Options:
    --br  Parse the backrub dataframe
    --cen_temp_cycles=NUM  How many centroid temp cycles (for loop
    modeler movers)?  [default:]
    --fa_temp_cycles=NUM  How many fa temp cycles (for loop modeler
    movers)?  [default:]
    --fast  For loop modelers, turn on test run
    --task=INT  Choose a specific task to execute (debuggin purposes
    only)
    --nomutate  Don't mutate before applying mover (testing only)
"""
import sys
import docopt
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
    args = docopt.docopt(__doc__)
    print(args)
    denom = 8
    pdbredo_directory = '/wynton/home/kortemme/krivacic/pdb_redo'
    shell=10
    if args['--task']:
        task_num = int(args['--task'])
    else:
        task_num = int(os.environ['SGE_TASK_ID']) - 1
    #task_num = 518 # make sure to subtract 1 from SGE TASK ID for the real thing
    num_models = 50

    mover = args['<mover>']

    #offset = (int(sys.argv[4]) - 1) * num_models * 100
    #task_num += offset
    row_num = task_num//denom

    # backrub variable tells us if we're reading from the backrub
    # benchmark dataframe, NOT whether we're using the backrub mover
    if args['--br']:
        df = import_backrub_dataframe(args['<input_df>'])
    else:
        df = import_benchmark_dataframe(args['<input_df>'])
    print(df)

    row = df.loc[row_num]
    wt_pdbid = row['wt'].lower()
    mut_pdbid = row['mutant'].lower()

    constrain = 'constrained' if (task_num%2 == 0) else 'unconstrained'
    outdir_final = os.path.join(args['<output_dir>'], mut_pdbid + '_' + wt_pdbid, constrain)
    pklout = outdir_final + '/results_' + str(task_num%denom) + '.pkl'
    if os.path.exists(pklout):
        print('Job already completed. Exiting.')
        sys.exit()

    if not os.path.exists(outdir_final):
        os.makedirs(outdir_final, exist_ok=True)

    '''
    Going to need to get all focus residues on their own line, so we can calc
    bb rmsd separately.
    '''
    init("-ignore_unrecognized_res -ex1 -ex2 " +\
            "-lh:db_path=/wynton/home/kortemme/krivacic/rosetta/database/loophash_db/ " +\
            "-lh:loopsizes 6 -pdb_gz "+ \
            "-total_threads 1")


    default_sfxn = create_score_function('ref2015')

    if 'TMPDIR' in os.environ:
        os_tmp = os.environ['TMPDIR']
    elif args['--task']:
        os_tmp = os.path.join(os.environ['HOME'], 'temp')
    else:
        os_tmp = os.path.join('/scratch',os.environ['USER'])
    outdir_temp = os.path.join(os_tmp, str(task_num))
    if not os.path.exists(outdir_temp):
        os.makedirs(outdir_temp, exist_ok=True)

    # Try to get poses from pdbredo, otherwise download from rcsb.
    wt_pose = custom_open(wt_pdbid,
            prefix='/wynton/home/kortemme/krivacic/intelligent_design/sidechain_directed_design_pyrosetta/backrub_pointmutant_benchmark/benchmark_pdbs',
            suffix='.pdb')
    wt_score = default_sfxn(wt_pose)

    mut_pose = custom_open(mut_pdbid,
            prefix='/wynton/home/kortemme/krivacic/intelligent_design/sidechain_directed_design_pyrosetta/backrub_pointmutant_benchmark/benchmark_pdbs',
            suffix='.pdb')
    mut_score = default_sfxn(mut_pose)

    if args['--br']:
        wtfocus = wt_pose.pdb_info().pdb2pose(row['wt_chain'],
                row['wt_resnum'])
        print('WTFOCUS')
        print(wtfocus)
        mutfocus = mut_pose.pdb_info().pdb2pose(row['mut_chain'],
                row['mut_resnum'])
        print('MUTFOCUS')
        print(mutfocus)
        focus = Mismatch(mutfocus, wtfocus)
        print('MISMATCH TARGET')
        print(focus.target)
        print('MISMATCH MOBILE')
        print(focus.mobile)
    else:
        focus = Mismatch(int(row['mut_res']), int(row['wt_res']))
    #mut_pair = MutantPair(mut_pose, wt_pose, [focus], shell=shell)

    print('-----------------------------------')
    print('COMPARING THE FOLLOWING 2 PDBIDs:')
    print('{}, resi {} chain {}, rosetta res {}'.format(wt_pdbid,
        row['wt_resnum'], row['wt_chain'], wtfocus))
    print('{}, resi {} chain {}, rosetta res {}'.format(mut_pdbid,
        row['mut_resnum'], row['mut_chain'], mutfocus))
    print('TASK NUMBER IS ' + str(task_num))
    print('ROW NUMBER IS ' + str(row_num))
    print('-----------------------------------')

    # Minimize poses before going any further
    mut_chain_selector = rosetta.core.select.residue_selector.ChainSelector(
            row['mut_chain']
    )
    mut_minimizer = get_minimizer(res_selector_to_size_list(mut_chain_selector.apply(mut_pose)),[])
    print('Minimizing mut')
    mut_minimizer.apply(mut_pose)
    mut_score_minimized = default_sfxn(mut_pose)
    wt_chain_selector = rosetta.core.select.residue_selector.ChainSelector(
            row['wt_chain']
    )
    wt_minimizer = get_minimizer(res_selector_to_size_list(wt_chain_selector.apply(wt_pose)), [])
    print('Minimizing wt')
    wt_minimizer.apply(wt_pose)
    wt_score_minimized = default_sfxn(mut_pose)

    output_data = []

    for jobnum in range(0, num_models):
        print('Running job {}'.format(jobnum))

        out_dict = deepcopy(row.to_dict())
        #wt_pose_copy = deepcopy(wt_pose)
        wt_pose_copy = custom_open(wt_pdbid,
                prefix='/wynton/home/kortemme/krivacic/intelligent_design/sidechain_directed_design_pyrosetta/backrub_pointmutant_benchmark/benchmark_pdbs',
                suffix='.pdb')
        wt_minimizer.apply(wt_pose_copy)
        if not os.path.exists(os.path.join(outdir_final,
            wt_pdbid + '_min.pdb.gz')):
            wt_pose_copy.dump_scored_pdb(os.path.join(outdir_final,
                wt_pdbid + '_min.pdb.gz'), default_sfxn)
        #mut_pose_copy = deepcopy(mut_pose)
        mut_pose_copy = custom_open(mut_pdbid,
                prefix='/wynton/home/kortemme/krivacic/intelligent_design/sidechain_directed_design_pyrosetta/backrub_pointmutant_benchmark/benchmark_pdbs',
                suffix='.pdb')
        mut_minimizer.apply(mut_pose_copy)
        if not os.path.exists(os.path.join(outdir_final,
            mut_pdbid + '_min.pdb.gz')):
            wt_pose_copy.dump_scored_pdb(os.path.join(outdir_final,
                mut_pdbid + '_min.pdb.gz'), default_sfxn)
        ##focus_res = int(row['mut_res'])
        ##motif_dict = {focus_res:row['wt_restype']}
        if constrain == 'constrained':
            cst = True
        else:
            cst = False

        designable, repackable, task_factory, mut_pair = \
                prepare_pdbids_for_modeling(wt_pose_copy, mut_pose_copy, [focus],
                        constrain=cst, shell=shell)

        out_dict['pre_score_wt'] = wt_score
        out_dict['pre_score_mut'] = mut_score
        out_dict['pre_score_wt_min'] = wt_score_minimized
        out_dict['pre_score_mut_min'] = mut_score_minimized

        print('MOTIFS HERE')
        print(mut_pair.motif_dict)

        out_dict['pre_rmsd'] = mut_pair.aligner.bb_rmsd
        out_dict['pre_dist'] = distance_rosetta(mut_pair.aligner.target,
                focus.target,
                mut_pair.aligner.mobile, focus.mobile)
        out_dict['mover'] = mover
        out_dict['decoy_number'] = jobnum
        print(out_dict)


        if not args['--nomutate']:
            for key in mut_pair.motif_dict:
                aatype = oneletter_to_threeletter(mut_pair.motif_dict[key]).upper()
                mut_res = rosetta.protocols.simple_moves.MutateResidue(key,
                        aatype)
                print('MADE MUT_RES OBJ')
                mut_res.apply(mut_pair.aligner.target) # apply() or make_mutation()?
                print('APPLIED MUT_RES')

            print('Minimizing mut after mutation')
            mut_minimizer.apply(mut_pair.aligner.target)

        if args['--fast']:
            fast=True
        else:
            fast=False

        # Get modeler object for corresponding mover type
        if mover == 'ngk':
            modeler = get_loop_modeler(mut_pair.aligner.target, designable, repackable,
                    focus.target, task_factory=task_factory, fast=fast,
                    mover='ngk', resbuffer=4)
        elif mover == 'ngkf':
            modeler = get_loop_modeler(mut_pair.aligner.target, designable, repackable,
                    focus.target, task_factory=task_factory, fast=True,
                    mover='ngk', resbuffer=4)
        elif mover == 'lhk':
            modeler = get_loop_modeler(mut_pair.aligner.target,
                    designable, repackable, focus.target,
                    task_factory=task_factory, fast=fast,
                    mover='lhk', resbuffer=4)
        elif mover == 'fastdesign':
            modeler = fast_design(mut_pair.aligner.target, designable, repackable,
                    task_factory=task_factory)    
        elif mover == 'br':
            modeler = get_backrub_protocol(mut_pair.motif_dict,
                    shell=shell, kt=1.6, task_factory=task_factory,
                    ntrials=100, stride=100)
        elif mover == 'jacobi':
            modeler = get_jacobi_refiner(mut_pair.aligner.target,
                    focus.target, resbuffer=4)
        elif mover=='ngk_jacobi' or mover=='lhk_jacobi':
            modeler = Modeler()
            loopmodeler = get_loop_modeler(mut_pair.aligner.target, designable, repackable,
                    focus.target, task_factory=task_factory, fast=fast,
                    mover=mover.split('_')[0], resbuffer=4)
            jacobi = get_jacobi_refiner(mut_pair.aligner.target,
                    focus.target, resbuffer=4)
            repackable_post_jacobi = [x for x in range(focus.target - 4,
                focus.target + 5)]
            pack_mover = get_pack_rotamers(mut_pair.aligner.target,
                    repackable_post_jacobi, default_sfxn)
            modeler.add_modeler(loopmodeler)
            modeler.add_modeler(jacobi)
            modeler.add_modeler(pack_mover)


        # Set mover options for loop modelers
        if mover in ['lhk','ngk','ngkf']:
            if args['--cen_temp_cycles']:
                print('0---dfaklfjalsk; HIHI HIHI HI')
                modeler.centroid_stage().set_temp_cycles(int(args['--cen_etmp_cycles']),True)
            if args['--fa_temp_cycles']:
                modeler.fullatom_stage().set_temp_cycles(int(args['--fa_temp_cycles']),True)
                print(modeler.fullatom_stage().get_temp_cycles())

        start_time = time.time()
        #mut_pair.aligner.target.dump_pdb('TEST_OUT.pdb')
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
        print(out_dict)

        annotate_pdb(pdb_path, out_dict)
        annotate_pdb(pdb_path_rel, out_dict)

    df_out = pd.DataFrame.from_records(output_data)
    print(df_out)
    print('Job finished, transferring files...')
    with open(pklout, 'wb') as f:
        pickle.dump(df_out, f)
    finish_io(outdir_temp, outdir_final, str(task_num%denom))
