import pickle, sys, os, glob
import pandas as pd
from matplotlib import pyplot as plt


if __name__=='__main__':
    #by = 'final_score'
    #by = 'post_dist_relaxed'
    #by = 'post_dist'
    #by = 'post_rmsd'
    #by = 'post_rmsd_relaxed'

    directory = sys.argv[1]
    true = directory + '/constrained'
    false = directory + '/unconstrained'
    outdir = sys.argv[2]
    relaxed = sys.argv[3]
    if relaxed == 'y':
        cols = ['post_dist_relaxed','post_rmsd_relaxed']
    elif relaxed == 'n':
        cols = ['post_dist', 'post_rmsd']

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

        df_cst = pd.DataFrame(data_cst)
        df_uncst = pd.DataFrame(data_uncst)
        print(df_cst)
        print(df_uncst)

        for wt in df_cst['wt'].unique():
            for mut in df_cst['mutant'].unique():
                to_plot = []
                groups = []
                zorders = []

                # Add constrained data
                x = df_cst[(df_cst['wt'] == wt) &
                        (df_cst['mutant'] == mut)][by]
                y = df_cst[(df_cst['wt'] == wt) &
                        (df_cst['mutant'] == mut)]['final_score']
                if len(x) > 0 and len(y) > 0:

                    to_plot.append((x, y))
                    groups.append(mut + ' to ' + wt + ' constrained')
                    zorder = 2
                    zorders.append(zorder)

                # Add unconstrained data
                x = df_uncst[(df_uncst['wt'] == wt) &
                        (df_uncst['mutant'] == mut)][by]
                y = df_uncst[(df_uncst['wt'] == wt) &
                        (df_uncst['mutant'] == mut)]['final_score']
                if len(x) > 0 and len(y) > 0:

                    to_plot.append((x, y))
                    groups.append(mut + ' to ' + wt + ' unconstrained')
                    zorder = 1
                    zorders.append(zorder)

                # Create plot
                if len(to_plot) > 0:
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 1, 1)

                    for data, group, order in zip(to_plot, groups, zorders):
                        x, y = data
                        ax.scatter(x, y, label=group, zorder=order)

                    if by=='post_dist_relaxed' or by=='post_dist':
                        pre_dist = df_cst[(df_cst['mutant'] == mut) &
                            (df_cst['wt'] == wt)]['pre_dist']
                        ax.axvline(pre_dist[0])
                    elif by=='post_rmsd_relaxed' or by=='post_rmsd':
                        ax.axvline(df_cst[(df_cst['mutant'] == mut) &
                            (df_cst['wt'] == wt)]['pre_rmsd'][0])
                    
                    plt.legend(loc=1)
                    plt.xlabel(by)
                    plt.ylabel('Total score')
                    #plt.savefig(outdir + '/fastdesign_' + by + '_' + groups[0] +
                            #'_' + groups[1] + '.png')
                    plt.show()
