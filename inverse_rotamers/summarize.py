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

def summarize(input_dir, summary='mean', force=False, relaxed=False,
        by='dist', threshold=None):
    if threshold and summary=='percent_improved':
        outpath = os.path.join(input_dir,
                'summary_of_{}_by_{}_thresh_{}.pkl'\
                        .format(by, summary, threshold))
    else:
        outpath = os.path.join(input_dir,
            'summary_of_{}_by_{}.pkl'.format(by, summary))
    if not force:
        if os.path.exists(outpath):
            with open(outpath, 'rb') as f:
                return pickle.load(f)

    avgs = []
    for pair in os.listdir(input_dir):
        curr_dir = os.path.join(input_dir, pair)

        if os.path.isdir(curr_dir):

            for cst in os.listdir(curr_dir):
                data_list = []
                folder = os.path.join(curr_dir, cst)

                #df = pd.DataFrame(data_list)
                #if os.path.exists(folder + '/results.pkl'):
                # Open all resulst pkl files
                for filename in glob.glob(folder + '/results*'):
                    with open(filename, 'rb') as f:
                        try:
                            df = pickle.load(f)
                            data_list.append(df)
                        except:
                            print('Error opening file ', filename)
                avgs_dict = {}
                if len(data_list) > 0:
                    df = pd.concat(data_list, ignore_index=True)
                    for col in df.columns:
                        if not df.dtypes[col] == np.object:
                            avgs_dict[col + '_std'] = df[col].std
                            if summary == 'mean':
                                avgs_dict[col + '_sum'] = df[col].mean()
                            elif summary == 'low_score':
                                ycol = 'final_score' if relaxed else 'post_score'
                                avgs_dict[col + '_sum'] = \
                                        df.loc[df[ycol].idxmin(),
                                                col]
                            elif summary == 'median':
                                avgs_dict[col + '_sum'] = df[col].median()
                            elif summary == 'structure':
                                ycol = 'post_{}_relaxed'.format(by) if relaxed else\
                                    'post_{}'.format(by)
                                df['delta_{}'.format(by)] = df[ycol] - df['pre_{}'.format(by)]
                                avgs_dict[col + '_sum'] = \
                                        df.loc[df['delta_{}'.format(by)].idxmin(),
                                                col]
                            elif summary == 'percent_improved':
                                # Columns pre-mover need no change
                                if col.split('_')[0] == 'pre':
                                    avgs_dict[col + '_sum'] = df[col][0]
                                else:
                                    splt = col.split('_')
                                    if len(splt) > 1:
                                        subname = col.split('_')[1]
                                        if not threshold:
                                            threshold = 0.0
                                        if subname == 'dist' or subname == 'rmsd':
                                            fraction = df[df[col] < (1 -
                                                threshold) * df['pre_' +
                                                subname]][col].size / df['pre_' +
                                                        subname].size
                                            avgs_dict[col + '_sum'] = 100 * fraction
                                        else:
                                            avgs_dict[col + '_sum'] =\
                                                df.loc[0,col]
                                    else:
                                        avgs_dict[col + '_sum'] =\
                                                df.loc[0,col]
                        else:
                            avgs_dict[col + '_sum'] = df[col][0]

                            if cst == 'constrained':
                                avgs_dict['constrained'] = True
                            else:
                                avgs_dict['constrained'] = False

                            avgs_dict['path'] = folder
                            
                    avgs.append(avgs_dict)
    avgs = pd.DataFrame(avgs)
    with open(outpath, 'wb') as f:
        pickle.dump(avgs, f)
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
    df = summarize(input_dir, summary=summary, relaxed=relaxed, by=mid)
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
            ylabel=y, title=title, unitline=unitline, markersize=4)


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
        #cmd = 'python3.7 ' +\
        cmd = 'python3 ' +\
        os.environ['HOME'] + '/intelligent_design/sidechain_directed_design_pyrosetta/inverse_rotamers/plot_funnels.py '\
        + path + ' ' + relaxed_str

        os.system(cmd)

    #fig.canvas.mpl_connect('button_press_event', on_pick)
    fig.canvas.mpl_connect('pick_event', on_pick)
    plt.show()
