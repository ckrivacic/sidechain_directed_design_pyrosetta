import pandas as pd
import numpy as np
import sys, os, glob
import pickle
from test_utils import plot
from numeric import euclidean_distance

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

                #df = pd.DataFrame(data_list)
                #if os.path.exists(folder + '/results.pkl'):
                for filename in glob.glob(folder + '/results*'):
                    with open(filename, 'rb') as f:
                        df = pickle.load(f)

                    avgs_dict = {}
                    if df.shape[0] > 0:
                        for col in df.columns:
                            try:
                                if summary == 'mean':
                                    avgs_dict[col + '_sum'] = df[col].mean()
                                elif summary == 'low_score':
                                    avgs_dict[col + '_sum'] = \
                                            df[col].loc[[df['final_score'].idxmin()]]
                                elif summary == 'median':
                                    avgs_dict[col + '_sum'] = df[col].median()
                                elif summary == 'percent_improved':
                                    if col.split('_')[0] == 'pre':
                                        avgs_dict[col + '_sum'] = df[col][0]
                                    else:
                                        subname = col.split('_')[1]
                                        if subname == 'dist' or subname == 'rmsd':
                                            fraction = df[df[col] < df['pre_' +
                                                subname]][col].size / df['pre_' +
                                                        subname].size
                                            if fraction == 1.0:
                                                print(df[curr_dir][0])
                                            avgs_dict[col + '_sum'] = 100 * fraction
                                    #avgs_dict[col + '_sum'] = df[df[]]
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
    #mid = 'dist'
    if len(sys.argv) < 3:
        mid = input('Analyze distance (dist) or RMSD (rmsd)? ')
    if len(sys.argv) < 4:
        summary = input('Enter data analysis method (mean, median, low_score, percent_improved): ')
    if len(sys.argv) > 3:
        mid = sys.argv[2]
        summary = sys.argv[3]
    relaxed = False
    x = 'pre_' + mid + '_sum'
    if not relaxed:
        y = 'post_' + mid + '_sum'
    if relaxed:
        y = 'post_' + mid + '_relaxed_sum'
    #summary = 'mean'

    input_dir = sys.argv[1]
    df = summarize(input_dir, summary=summary)
    df_cst = df[df['constrained']==True]
    df_uncst = df[df['constrained']==False]
    data1 = (df_cst[x],
        df_cst[y])
    data2 = (df_uncst[x],
            df_uncst[y])
    data = [data1, data2]
    groups = ['constrained', 'unconstrained']
    if summary != 'percent_improved':
        unitline = True
    else:
        unitline = False
    title = mid + ' comparison summarized by ' + summary
    plt, fig, ax = plot(data, groups=groups, xlabel=x,
            ylabel=y, title=title, unitline=unitline)


    def on_pick(event):
        if not hasattr(event, 'ind'):
            return True
        ind = event.ind
        #for a, b in enumerate(ind):
        data = event.artist.get_offsets()
        if event.artist.get_label() == 'constrained':
            pick_df = df_cst
        else:
            pick_df = df_uncst

        xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
        if len(ind) > 1:
            low_dist = None
            low_index = None
            for i in ind:
                dist = euclidean_distance(data[i], (xmouse,ymouse))
                if (not low_dist) or (dist < low_dist):
                    low_dist = dist
                    low_index = i
        else:
            low_index = ind[0]

        pathlist = pick_df.iloc[low_index]['path'].split('/')
        #pathlist = pick_df['path'][low_index].split('/')
        if pathlist[-1] == 'constrained' or pathlist[-1] == 'unconstrained':
            pathlist = pathlist[:-1]
        path = os.path.abspath('/'.join(pathlist))
        #cmd = 'show_my_designs ' + path + '/*' # Need pdb files

        relaxed_str = 'y' if relaxed else 'n'
        cmd = 'python3.7 ' +\
        os.environ['HOME'] + '/intelligent_design/sidechain_directed_design/inverse_rotamers/plot_funnels.py '\
        + path + ' ' + relaxed_str

        os.system(cmd)

    #fig.canvas.mpl_connect('button_press_event', on_pick)
    fig.canvas.mpl_connect('pick_event', on_pick)
    plt.show()
