import pickle, sys, os, glob
import pandas as pd
from matplotlib import pyplot as plt


if __name__=='__main__':
    by = 'post_dist_relaxed'

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

    for data, group in zip(to_plot, groups):
        x, y = data
        ax.scatter(x, y, label=group)
    
    ax.plot([0.25, 2.25] , [0.25, 2.25], 'k-')
    plt.legend(loc=2)
    plt.show()
