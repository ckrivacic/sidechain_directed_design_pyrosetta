"""
Usage:
    bar_plots <folders>... [options]

Options:
    --by=BY  Analyze by distance or rmsd  [default: dist]
    --sum=SUM  How to summarize data of individual folders (can be one of mean,
    median, low_score, or percent_improved)  [default: percent_improved]
    --relaxed, -r  Choose to use the post-relaxed structures/data
    --cst=CST  Constrained, unconstrained, or both?  [default: both]

"""


from summarize import summarize
import docopt
import pandas as pd
import numpy as np

if __name__ == "__main__":
    from matplotlib import pyplot as plt

    args = docopt.docopt(__doc__)
    print(args)
    mid = args['--by']
    summary = args['--sum']
    relaxed = args['--relaxed']
    x = 'pre_' + mid + '_sum'
    if not relaxed:
        y = 'post_' + mid + '_sum'
    if relaxed:
        y = 'post_' + mid + '_relaxed_sum'

    bins = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.6, 2.0, 10.0]
    labels = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.6, 2.0]
    ind = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots()

    i = 0
    for input_dir in args['<folders>']:
        folder_label = input_dir.split('/')[-3]
        df = summarize(input_dir, summary=summary)

        yvals_cst = []
        stderr_cst = []
        yvals_uncst = []
        stderr_uncst = []

        df['binned'] = pd.cut(df[x], bins=bins, labels=labels)
        df.to_csv('test_out.csv')
        print(df)
        
        for label in labels:
            #print(label)
            cst_mean = df[(df['constrained']==True) &
                (df['binned']==label)][y].mean()
            yvals_cst.append(cst_mean)
            stderr_cst.append(df[(df['constrained']==True) &
                (df['binned']==label)][y].std())

            yvals_uncst.append(df[(df['constrained']==False) &
                (df['binned']==label)][y].mean())
            stderr_uncst.append(df[(df['constrained']==False) &
                (df['binned']==label)][y].std())

        print(yvals_cst)
        if args['--cst'] == 'cst' or args['--cst'] == 'both':
            ax.bar(ind + (2 * i) * width / (len(args['<folders>'])),
                    yvals_cst, width/(len(args['<folders>'])), yerr=stderr_cst,
                    label=folder_label + '_cst')
        if args['--cst'] == 'uncst' or args['--cst'] == 'both':
            ax.bar(ind + (2 * i + 1) * width / (len(args['<folders>'])),                
                    yvals_uncst, width/(len(args['<folders>'])), yerr=stderr_uncst,
                    label=folder_label + '_uncst')

        i += 1


    ax.set_ylabel('Mean percent improved decoys')
    ax.set_xlabel('Binned distance of starting structure to WT')
    ax.set_title('Sampling methods broken down by distance of point mutant to WT')
    ax.set_xticks(ind)
    ax.set_xticklabels([x for x in labels])
    ax.legend()

    fig.tight_layout()
    plt.show()
