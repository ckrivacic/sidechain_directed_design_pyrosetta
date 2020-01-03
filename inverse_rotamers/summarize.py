import pandas as pd
import numpy as np
import sys, os, glob
import pickle
from test_utils import plot

'''
Input: directory of directories with constrained and unconstrained
subdirectories. Should look like this:

    input_directory/pdb_pair/[constrained or unconstrained]/*.pdb

For each subdirectory, summarize:
    - % designs improved from pre_dist or pre_rmsd

'''

def summarize(input_dir, summary='mean'):
    avgs = []
    for pair in os.listdir(input_dir):
        curr_dir = os.path.join(input_dir, pair)

        if os.path.isdir(curr_dir):

            for cst in os.listdir(curr_dir):
                data_list = []
                folder = os.path.join(curr_dir, cst)

                for filename in glob.glob(folder + '/results*'):
                    with open(filename, 'rb') as f:
                        data = pickle.load(f)
                    data_list.append(data)

                df = pd.DataFrame(data_list)
                avgs_dict = {}

                if df.shape[0] > 0:
                    for col in df.columns:
                        try:
                            if summary == 'mean':
                                print('hey what')
                                avgs_dict[col + '_sum'] = df[col].mean()
                            elif summary == 'low_score':
                                avgs_dict[col + '_sum'] = \
                                        df[col].loc[[df['final_score'].idxmin()]]
                        except:
                            avgs_dict[col + '_sum'] = df[col][0]

                    if cst == 'constrained':
                        avgs_dict['constrained'] = True
                    else:
                        avgs_dict['constrained'] = False

                    avgs_dict['path'] = folder
                    
                    avgs.append(avgs_dict)

    return pd.DataFrame(avgs)


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    input_dir = sys.argv[1]
    df = summarize(input_dir, summary='low_score')
    data1 = (df[df['constrained']==True]['pre_rmsd_sum'],
        df[df['constrained']==True]['post_rmsd_relaxed_sum'])
    data2 = (df[df['constrained']==False]['pre_rmsd_sum'],
            df[df['constrained']==False]['post_rmsd_relaxed_sum'])
    data = [data1, data2]
    groups = ['constrained', 'unconstrained']
    plt = plot(data, groups=groups, xlabel='rmsd pre-mover',
            ylabel='rmsd post-mover', title='rmsd comparison', unitline=True)
    plt.show()
