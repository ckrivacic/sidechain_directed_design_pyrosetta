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
    --binby=BINBY  What to bin by (distance, rmsd, relsurf)  
    [default:]

"""


from summarize import summarize
import docopt
import pandas as pd
import numpy as np
from itertools import cycle


def row_index_from_pdbs(mut, wt, df):
    row = df.index[(df['mutant_sum'] == mut) & (df['wt_sum'] == wt)]
    return row.tolist()


if __name__ == "__main__":
    from matplotlib import pyplot as plt

    name_dict = {'fastdesign':'FastRelax ','ngkf':'Next-gen KIC (fast) ',
            'br':'Backrub (100 trials) '}
    colors = cycle(['navy', 'cornflowerblue', 'darkgreen',
            'lightgreen', 'firebrick', 'lightcoral'])

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

    if args['--binby'] == 'relsurf_sum':
        bins = [0, 0.02, 0.04, 0.06, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
                1.0]
        labels = [0, 0.02, 0.04, 0.06, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
                0.5]
        labels_str = ['0-0.02', '0.02-0.04', '0.04-0.06', '0.06-0.1',
                '0.1-0.15', '0.15-0.2', '0.2-0.25', '0.25-0.3',
                '0.3-0.4', '0.4-0.5', '>0.5']
    else:
        bins = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0,10.0]
        labels = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]
        labels_str = ['0-0.1', '0.1-0.2', '0.2-0.4', '0.4-0.6',
                '0.6-0.8', '0.8-1.0', '>1.0']
    '''
    else:
        bins = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 10.0]
        labels = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5]
        labels_str = ['0-0.1', '0.1-0.2', '0.2-0.4', '0.4-0.6',
                '0.6-0.8', '0.8-1.0', '1.0-1.5', '>1.5']
    '''
    
    ind = np.arange(len(labels))
    width = 0.45
    fig, ax = plt.subplots()
    #fig = plt.figure(figsize=(10,4.8))

    # i keeps track of which folder we're on for formatting purposes
    i = 0
    ticklabels = []
    pdbs = []
    #pdbs = ['1a6g_1a6m','1b9k_1kyf','1bu7_1jme','1ffr_1edq','1hvf_1hve','1i4o_1kmc','1qop_1k8y','1qul_1quj']
    #pdbs = ['1ccp_3ccp','1ec0_1hvs','5eaa_1qit','7pti_1aal']
    for input_dir in args['<folders>']:
        folder_label = name_dict[input_dir.split('/')[-3]]
        df = summarize(input_dir, summary=summary, force=force,
                relaxed=relaxed, by=mid)

        yvals_cst = []
        stderr_cst = []
        yvals_uncst = []
        stderr_uncst = []

        print(df.columns)
        if args['--binby']:
            df['binned'] = pd.cut(df[args['--binby']], bins=bins,
                    labels=labels)
        else:
            df['binned'] = pd.cut(df[x], bins=bins, labels=labels)
            #df = df[df['binned'] != 1.0]
    
        by_dict = {'dist':'distance', 'rmsd':'RMSD'}
        sum_dict = {'low_score':'Lowest scoring', 'mean':'Mean',
                'median':'Median', 
                'percent_improved':'Percentage of improved',
                'structure':'Closest structure '}
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
                        yvals_cst, width/(len(args['<folders>'])), #yerr=stderr_cst,
                        label=folder_label + '(constrained)',
                        color = next(colors))
            if args['--cst'] == 'uncst' or args['--cst'] == 'both':
                ax.bar(ind + (2 * i + 1) * width / (len(args['<folders>'])),                
                        yvals_uncst, width/(len(args['<folders>'])), #yerr=stderr_uncst,
                        label=folder_label + '(unconstrained)',
                        color=next(colors))

            ax.set_ylabel('Mean percentage of improved decoys')
            ax.set_xlabel('Binned {} (Å) of starting structure α-carbon to WT'\
                     .format(by_dict[args['--by']]))
            ax.set_title('Sampling methods broken down by {} of point mutant to WT'\
                    .format(by_dict[args['--by']]))
            ax.set_xticks(ind + width/(len(args['<folders>'])/2))
            ax.set_xticklabels(labels_str)

        else:
            import random
            num_examples = 10
            label = float(args['--bin'])
            df_cst = df[df['constrained']==True]
            df_uncst = df[df['constrained']==False]
            df_cst = df_cst[df_cst['binned']==label]
            print('{} examples chosen out of {}'.format(num_examples,
                len(df_cst.index)))
            df_uncst = df_uncst[df_uncst['binned']==label]

            yvals_cst = []
            yvals_uncst = []
            if len(pdbs) < 1:
                indices = random.sample(list(df_cst.index.values),min(num_examples,
                    len(df_cst.index)))

                for index in indices:
                    if not args['--sum'] == 'percent_improved':
                        yval_cst = df_cst.loc[index, y] - df_cst.loc[index, x]
                    else:
                        yval_cst = df_cst.loc[index, y]

                    wt = df_cst.loc[index, 'wt_sum']
                    wtchain = df_cst.loc[index, 'wt_chain_sum']
                    wtresnum = df_cst.loc[index, 'wt_resnum_sum']
                    mut = df_cst.loc[index, 'mutant_sum']
                    mutchain = df_cst.loc[index, 'mut_chain_sum']
                    mutresnum = df_cst.loc[index, 'mut_resnum_sum']

                    uncst_row = df_uncst.index[(df_uncst['wt_sum']==wt) &
                            (df_uncst['mutant_sum']==mut)].tolist()
                    if len(uncst_row) > 0:
                        uncst_row = uncst_row[0]
                        if not args['--sum'] == 'percent_improved':
                            yval_uncst = df_uncst.loc[uncst_row, y] - df_uncst.loc[uncst_row, x]
                        else:
                            yval_uncst = df_uncst.loc[uncst_row, y]
                        yvals_uncst.append(yval_uncst)
                        yvals_cst.append(yval_cst)
                        ticklabels.append('{} ({}{}),\n{} ({}{})'\
                                .format(mut, mutchain, mutresnum, wt,
                                    wtchain, wtresnum))
                        pdbs.append('_'.join([mut, wt]))
                print(ticklabels)
            else:
                for tick in pdbs:
                    mut = tick.split('_')[0]
                    wt = tick.split('_')[1]
                    cst_row = row_index_from_pdbs(mut.upper(),
                            wt.upper(), df_cst)
                    print(cst_row)
                    if len(cst_row) > 0:

                        cst_row = cst_row[0]

                        # Block only relevant when custom pdbs are used
                        '''
                        wt = df_cst.loc[cst_row, 'wt_sum']
                        wtchain = df_cst.loc[cst_row, 'wt_chain_sum']
                        wtresnum = df_cst.loc[cst_row, 'wt_resnum_sum']
                        mut = df_cst.loc[cst_row, 'mutant_sum']
                        mutchain = df_cst.loc[cst_row, 'mut_chain_sum']
                        mutresnum = df_cst.loc[cst_row, 'mut_resnum_sum']
                        ticklabels.append('{} ({}{}),\n{} ({}{})'\
                                .format(mut, mutchain, mutresnum, wt,
                                    wtchain, wtresnum))
                        '''
                        
                        if not args['--sum'] == 'percent_improved':
                            yval_cst = df_cst.loc[cst_row, y] - df_cst.loc[cst_row, x]
                        else:
                            yval_cst = df_cst.loc[cst_row, y]
                        yvals_cst.append(yval_cst)
                    else:
                        yvals_cst.append(0)

                    uncst_row = row_index_from_pdbs(mut.upper(),
                            wt.upper(), df_uncst)
                    if len(uncst_row) > 0:
                        uncst_row = uncst_row[0]
                        if not args['--sum'] == 'percent_improved':
                            yval_uncst = df_uncst.loc[uncst_row, y] - df_uncst.loc[uncst_row, x]
                        else:
                            yval_uncst = df_uncst.loc[uncst_row, y]
                        yvals_uncst.append(yval_uncst)
                    else:
                        yvals_uncst.append(0)


            assert(len(yvals_cst) == len(yvals_uncst))
            ind = np.arange(len(yvals_cst))

            if args['--cst'] == 'cst' or args['--cst'] == 'both':
                ax.bar(ind + (2 * i) * width/(len(args['<folders>'])),
                        yvals_cst, width / (len(args['<folders>'])),
                        label=folder_label + '(constrained)',
                        color=next(colors))

            if args['--cst'] == 'uncst' or args['--cst'] == 'both':
                ax.bar(ind + (2 * i + 1) * width/(len(args['<folders>'])),
                        yvals_uncst, width / (len(args['<folders>'])),
                        label=folder_label + '(unconstrained)',
                        color=next(colors))

            ax.set_ylabel('{} Δ-{}'.format(sum_dict[args['--sum']],
                by_dict[args['--by']]))
            ax.set_xlabel('Mutant-wildtype pair')
            ax.set_xticks(ind + width/(len(args['<folders>'])/2))
            ax.set_xticklabels(ticklabels)
            ax.set_title('{} examples from bin {}'.format(num_examples,
                label))
            plt.xticks(rotation=70)



        i += 1


    ax.legend()
    fig.tight_layout()
    plt.show()
