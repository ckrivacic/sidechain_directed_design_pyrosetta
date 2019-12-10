import pickle, sys, os, glob
import pandas as pd
from matplotlib import pyplot as plt


if __name__=='__main__':
    #by = 'final_score'
    #by = 'post_dist_relaxed'
    #by = 'post_dist'
    #by = 'post_rmsd'
    by = 'post_rmsd_relaxed'

    true = sys.argv[1]
    false = sys.argv[2]

    difference = 139

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
    to_plot = []
    groups = []
    for wt in df['wt'].unique():
        for mut in df['mutant'].unique():
            x = df[(df['constrain'] == False) & (df['wt'] == wt) &
                    (df['mutant'] == mut)][by]
            y = df[(df['constrain'] == True) & (df['wt'] == wt) &
                    (df['mutant'] == mut)][by]
            if len(x) > 0 and len(y) > 0:

                to_plot.append((x, y))
                groups.append(mut + ' to ' + wt)

    # Create plot
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    xmin_all = 9999.9
    xmax_all = -99999999.0
    ymax_all = -99999999.0
    ymin_all = 9999.9
    for data, group in zip(to_plot, groups):
        x, y = data
        xmin = min(x)
        ymin = min(y)
        xmax = max(x)
        ymax = max(y)
        ax.scatter(x, y, label=group)
        if xmin < xmin_all:
            xmin_all = xmin
        if ymin < ymin_all:
            ymin_all = ymin
        if xmax > xmax_all:
            xmax_all = xmax
        if ymax > ymax_all:
            ymax_all = ymax

    maximum = max(xmax, ymax)
    minimum = min(xmin, ymin)
    
    ax.plot([minimum, maximum], [minimum, maximum], 'k-')
    plt.legend(loc=2)
    plt.xlabel('No constraints')
    plt.ylabel('N, C, CA constraints')
    plt.show()
