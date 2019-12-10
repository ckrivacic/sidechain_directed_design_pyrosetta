import pickle, sys, os, glob
import pandas as pd
from matplotlib import pyplot as plt


if __name__=='__main__':
    #by = 'final_score'
    #by = 'post_dist_relaxed'
    #by = 'post_dist'
    #by = 'post_rmsd'
    #by = 'post_rmsd_relaxed'

    true = sys.argv[1]
    false = sys.argv[2]
    outdir = sys.argv[3]
    relaxed = sys.argv[4]
    if relaxed == 'y':
        cols = ['post_dist_relaxed','post_rmsd_relaxed']
    elif relaxed == 'n':
        cols = ['post_dist', 'post_rmsd']

    #difference = 139

    for by in cols:
        all_data = []
        for filename in glob.glob(true + '/results*'):
            f = open(filename, 'rb')
            all_data.append(pickle.load(f))
            f.close()

        for filename in glob.glob(false + '/results*'):
            f = open(filename, 'rb')
            all_data.append(pickle.load(f))
            f.close()

        df = pd.DataFrame(all_data)
        print(df)
        for wt in df['wt'].unique():
            for mut in df['mutant'].unique():
                to_plot = []
                groups = []
                zorders = []
                for bul in [True, False]:
                    x = df[(df['constrain'] == bul) & (df['wt'] == wt) &
                            (df['mutant'] == mut)][by]
                    y = df[(df['constrain'] == bul) & (df['wt'] == wt) &
                            (df['mutant'] == mut)]['final_score']
                    if len(x) > 0 and len(y) > 0:

                        to_plot.append((x, y))
                        constrained = ' constrained' if bul else ' unconstrained'
                        groups.append(mut + ' to ' + wt + constrained)
                        zorder = 2 if bul else 1
                        zorders.append(zorder)

                # Create plot
                if len(to_plot) > 0:
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 1, 1)

                    for data, group, order in zip(to_plot, groups, zorders):
                        x, y = data
                        ax.scatter(x, y, label=group, zorder=order)
                    
                    plt.legend(loc=1)
                    plt.xlabel(by)
                    plt.ylabel('Total score')
                    plt.savefig(outdir + '/fastdesign_' + by + '_' + groups[0] +
                            '_' + groups[1] + '.png')
                    plt.show()
