"""
Usage:
    bar_plots <folders>... [options]

Options:
    --by=BY  Analyze by distance or rmsd  [default: dist]
    --sum=SUM  How to summarize data of individual folders (can be one of mean,
    median, low_score, or percent_improved)  [default: percent_improved]
    --relaxed, -r  Choose to use the post-relaxed structures/data
    --cst=CST  Constrained, unconstrained, or both?  [default: both]
    --bin=BIN  Plot examples from a particular bin.  [default: ]
    --force, -f  Force to re-analyze if analysis pkl file is found.

"""


from summarize import summarize
import docopt
import pandas as pd
import numpy as np


def row_index_from_pdbs(mut, wt, df):
    row = df.index[(df['mutant_sum'] == mut) & (df['wt_sum'] == wt)]
    return row.tolist()


if __name__ == "__main__":
    from matplotlib import pyplot as plt

    args = docopt.docopt(__doc__)
    print(args)
    mid = args['--by']
    summary = args['--sum']
    relaxed = args['--relaxed']
    force = args['--force']
    x = 'pre_' + mid + '_sum'
    if not relaxed:
        y = 'post_' + mid + '_sum'
    if relaxed:
        y = 'post_' + mid + '_relaxed_sum'

    bins = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 10.0]
    labels = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0]
    ind = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots()

    # i keeps track of which folder we're on for formatting purposes
    i = 0
    ticklabels = []
    for input_dir in args['<folders>']:
        folder_label = input_dir.split('/')[-3]
        df = summarize(input_dir, summary=summary, force=force,
                relaxed=relaxed)

        yvals_cst = []
        stderr_cst = []
        yvals_uncst = []
        stderr_uncst = []

        df['binned'] = pd.cut(df[x], bins=bins, labels=labels)
    
        if not args['--bin']:
            ind = np.arange(len(labels))
            for label in labels:
                cst_mean = df[(df['constrained']==True) &
                    (df['binned']==label)][y].mean()
                yvals_cst.append(cst_mean)
                stderr_cst.append(df[(df['constrained']==True) &
                    (df['binned']==label)][y].std())

                yvals_uncst.append(df[(df['constrained']==False) &
                    (df['binned']==label)][y].mean())
                stderr_uncst.append(df[(df['constrained']==False) &
                    (df['binned']==label)][y].std())

            if args['--cst'] == 'cst' or args['--cst'] == 'both':
                ax.bar(ind + (2 * i) * width / (len(args['<folders>'])),
                        yvals_cst, width/(len(args['<folders>'])), yerr=stderr_cst,
                        label=folder_label + '_cst')
            if args['--cst'] == 'uncst' or args['--cst'] == 'both':
                ax.bar(ind + (2 * i + 1) * width / (len(args['<folders>'])),                
                        yvals_uncst, width/(len(args['<folders>'])), yerr=stderr_uncst,
                        label=folder_label + '_uncst')

            ax.set_ylabel('Mean percent improved decoys')
            ax.set_xlabel('Binned distance of starting structure to WT')
            ax.set_title('Sampling methods broken down by distance of point mutant to WT')
            ax.set_xticks(ind)
            ax.set_xticklabels([x for x in labels])

        else:
            import random
            num_examples = 50
            label = float(args['--bin'])
            df_cst = df[df['constrained']==True]
            df_uncst = df[df['constrained']==False]
            df_cst = df_cst[df_cst['binned']==label]
            print('{} examples chosen out of {}'.format(num_examples,
                len(df_cst.index)))
            df_uncst = df_uncst[df_uncst['binned']==label]

            yvals_cst = []
            yvals_uncst = []
            if i == 0:
                indices = random.sample(list(df_cst.index.values),min(num_examples,
                    len(df_cst.index)))

                for index in indices:
                    yval_cst = df_cst.loc[index, y] - df_cst.loc[index, x]

                    wt = df_cst.loc[index, 'wt_sum']
                    mut = df_cst.loc[index, 'mutant_sum']
                    uncst_row = df_uncst.index[(df_uncst['wt_sum']==wt) &
                            (df_uncst['mutant_sum']==mut)].tolist()
                    if len(uncst_row) > 0:
                        uncst_row = uncst_row[0]
                        yval_uncst = df_uncst.loc[uncst_row, y] - df_uncst.loc[uncst_row, x]
                        yvals_uncst.append(yval_uncst)
                        yvals_cst.append(yval_cst)
                        ticklabels.append(mut + '_' + wt)
            else:
                for tick in ticklabels:
                    mut = tick.split('_')[0]
                    wt = tick.split('_')[1]
                    cst_row = row_index_from_pdbs(mut, wt, df_cst)
                    if len(cst_row) > 0:
                        cst_row = cst_row[0]
                        yval_cst = df_cst.loc[cst_row, y] - df_cst.loc[cst_row, x]
                        yvals_cst.append(yval_cst)
                    else:
                        yvals_cst.append(0)

                    uncst_row = row_index_from_pdbs(mut, wt, df_uncst)
                    if len(uncst_row) > 0:
                        uncst_row = uncst_row[0]
                        yval_uncst = df_uncst.loc[uncst_row, y] - df_uncst.loc[uncst_row, x]
                        yvals_uncst.append(yval_uncst)
                    else:
                        yvals_uncst.append(0)


            assert(len(yvals_cst) == len(yvals_uncst))
            ind = np.arange(len(yvals_cst))

            if args['--cst'] == 'cst' or args['--cst'] == 'both':
                ax.bar(ind + (2 * i) * width/(len(args['<folders>'])),
                        yvals_cst, width / (len(args['<folders>'])),
                        label=folder_label + '_cst')

            if args['--cst'] == 'uncst' or args['--cst'] == 'both':
                ax.bar(ind + (2 * i + 1) * width/(len(args['<folders>'])),
                        yvals_uncst, width / (len(args['<folders>'])),
                        label=folder_label + '_uncst')

            ax.set_ylabel('{} delta-{}'.format(args['--sum'], args['--by']))
            ax.set_xlabel('PDB pair (mut_wt)')
            ax.set_xticks(ind)
            ax.set_xticklabels(ticklabels)
            ax.set_title('{} examples from bin {}'.format(num_examples,
                label))
            plt.xticks(rotation=70)



        i += 1


    ax.legend()
    fig.tight_layout()
    plt.show()
