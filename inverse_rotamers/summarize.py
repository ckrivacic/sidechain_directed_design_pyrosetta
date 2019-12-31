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

def summarize(input_dir):
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
                            avgs_dict[col + '_avg'] = df[col].mean()
                        except:
                            avgs_dict[col + '_avg'] = df[col][0]
                    if type(avgs_dict[col + '_avg']) != type('abcd'):
                        if np.isna(avgs_dict[col + '_avg']):
                            print(df[col])
                        if np.isnull(avgs_dict[col + '_avg']):
                            print(df[col])
                        if avgs_dict[col + '_avg'] == 'NaN':
                            print(df[col])

                    if cst == 'constrained':
                        avgs_dict['constrained'] = True
                    else:
                        avgs_dict['constrained'] = False

                    avgs_dict['path'] = folder
                    
                    avgs.append(avgs_dict)

    return pd.DataFrame(avgs)


if __name__=='__main__':
    from matplotlib import pyplot as plt
    input_dir = sys.argv[1]
    df = summarize(input_dir)
    data1 = (df[df['constrained']==True]['pre_rmsd_avg'],
        df[df['constrained']==True]['post_rmsd_relaxed_avg'])
    data2 = (df[df['constrained']==False]['pre_rmsd_avg'],
            df[df['constrained']==False]['post_rmsd_relaxed_avg'])
    data = [data1, data2]
    groups = ['constrained', 'unconstrained']
    #plt = plot(data, groups=groups)

    list_of_tuples = data
    zorders = None
    title = 'title'
    xlabel = 'x'
    ylabel = 'y'

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for datum, group in zip(data, groups):
        x, y = datum
        print(x)
        ax.scatter(x, y, label=group)

    plt.legend(loc=1)
    plt.xlabel('else')
    plt.ylabel('okay')
    plt.title('titular')

    plt.show()
