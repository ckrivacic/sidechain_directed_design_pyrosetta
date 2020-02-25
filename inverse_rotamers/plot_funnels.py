import pickle, sys, os, glob
import pandas as pd
from matplotlib import pyplot as plt
from numeric import euclidean_distance

if __name__=='__main__':
    #by = 'final_score'
    #by = 'post_dist_relaxed'
    #by = 'post_dist'
    #by = 'post_rmsd'
    #by = 'post_rmsd_relaxed'

    directory = sys.argv[1]
    true = directory + '/constrained'
    false = directory + '/unconstrained'
    #outdir = sys.argv[2]
    relaxed = sys.argv[2]
    if relaxed == 'y':
        cols = ['post_dist_relaxed','post_rmsd_relaxed']
    elif relaxed == 'n':
        cols = ['post_dist', 'post_rmsd','elapsed_time']

    #difference = 139

    for by in cols:
        data_cst = []
        for filename in glob.glob(true + '/results*'):
            f = open(filename, 'rb')
            data_cst.append(pickle.load(f))
            f.close()

        data_uncst = []
        for filename in glob.glob(false + '/results*'):
            f = open(filename, 'rb')
            data_uncst.append(pickle.load(f))
            f.close()

        cst_list = []
        for filename in glob.glob(true + '/results*.pkl'):
            with open(filename, 'rb') as f:
                df = pickle.load(f)
                cst_list.append(df)
        df_cst = pd.concat(cst_list)

        uncst_list = []
        for filename in glob.glob(false + '/results*.pkl'):
            with open(filename, 'rb') as f:
                df = pickle.load(f)
                uncst_list.append(df)
        df_uncst = pd.concat(uncst_list)

        #print(df_cst)
        #print(df_uncst)

        for wt in df_cst['wt'].unique():
            for mut in df_cst['mutant'].unique():
                ycol = 'final_score' if relaxed=='y' else 'post_score'
                to_plot = []
                groups = []
                zorders = []

                # Add constrained data
                x = df_cst[(df_cst['wt'] == wt) &
                        (df_cst['mutant'] == mut)][by]
                y = df_cst[(df_cst['wt'] == wt) &
                        (df_cst['mutant'] == mut)][ycol]
                if len(x) > 0 and len(y) > 0:

                    to_plot.append((x, y))
                    groups.append(mut + ' to ' + wt + ' constrained')
                    zorder = 2
                    zorders.append(zorder)

                # Add unconstrained data
                x = df_uncst[(df_uncst['wt'] == wt) &
                        (df_uncst['mutant'] == mut)][by]
                y = df_uncst[(df_uncst['wt'] == wt) &
                        (df_uncst['mutant'] == mut)][ycol]
                if len(x) > 0 and len(y) > 0:

                    to_plot.append((x, y))
                    groups.append(mut + ' to ' + wt + ' unconstrained')
                    zorder = 1
                    zorders.append(zorder)

                # Create plot
                if len(to_plot) > 0:

                    def on_pick(event):
                        if not hasattr(event, 'ind'):
                            return True
                        ind = event.ind
                        data = event.artist.get_offsets()
                        if event.artist.get_label().endswith('unconstrained'):
                            pick_df = df_uncst
                            cst = 'unconstrained'
                        else:
                            pick_df = df_cst
                            cst = 'constrained'

                        xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
                        if len(ind) > 1:
                            low_dist = None
                            low_index = None
                            for i in ind:
                                dist = euclidean_distance(data[i],
                                        (xmouse,ymouse))
                                if (not low_dist) or (dist < low_dist):
                                    low_dist = dist
                                    low_index = i
                        else:
                            low_index = ind[0]

                        row = pick_df.iloc[low_index]
                        pdb_file = row['path'].split('/')[-1]
                        print(pdb_file)
                        subfolder = '_'.join(pdb_file.split('_')[0:2])
                        pdb_path = os.path.join(directory, cst, pdb_file)
                        print(pdb_path)


                    fig = plt.figure()
                    ax = fig.add_subplot(1, 1, 1, picker=True)

                    for data, group, order in zip(to_plot, groups, zorders):
                        x, y = data
                        ax.scatter(x, y, label=group, zorder=order, picker=True)

                    if by=='post_dist_relaxed' or by=='post_dist':
                        #pre_dist = df_cst[(df_cst['mutant'] == mut) &
                        #    (df_cst['wt'] == wt)]['pre_dist']
                        #ax.axvline(pre_dist[0])
                        pre_dist = df['pre_dist'][0]
                        ax.axvline(pre_dist)
                    elif by=='post_rmsd_relaxed' or by=='post_rmsd':
                        #ax.axvline(df_cst[(df_cst['mutant'] == mut) &
                        #    (df_cst['wt'] == wt)]['pre_rmsd'][0])
                        pre_rmsd = df['pre_rmsd'][0]
                        ax.axvline(pre_rmsd)
                    
                    plt.legend(loc=1)
                    plt.xlabel(by)
                    fig.canvas.mpl_connect('pick_event', on_pick)
                    plt.ylabel('Total score')
                    #plt.savefig(outdir + '/fastdesign_' + by + '_' + groups[0] +
                            #'_' + groups[1] + '.png')
                    plt.show()
