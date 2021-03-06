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
    [default: dist]
    --bins=INT  How many bins  [default: 10]
    --restype=1letterAA  Specify a restype or comma-separated list of
    resypes to analyze.  [default: ]
    --delta, -d  Report y-values as deltas
    --pthresh=FLOAT  Between 0 and 1. How much better must a dist or
    rmsd be before it's considered in the percent_improved metric?
    [default: 0.0]
    --absthresh=FLOAT  Absolute threshold for considering percent
    improvement  [default: 0.0]
    --violin, -v  Make a violin plot instead of a bar plot

"""


from summarize import summarize
import matplotlib.patches as mpatches
import docopt
import pandas as pd
import numpy as np
from itertools import cycle


def row_index_from_pdbs(mut, wt, df):
    row = df.index[(df['mutant_sum'] == mut) & (df['wt_sum'] == wt)]
    return row.tolist()


if __name__ == "__main__":
    from matplotlib import pyplot as plt

    name_dict = {'fastdesign':'FastRelax ','ngk':'Next-gen KIC (fast) ',
            'br':'Backrub (100 trials) ', 'lhk':'Loophash KIC (fast)',
            'jacobi_refine':'Jacobi refinement',
            'ngk_jacobi':'NGK->Jacobi refine->pack'}
    colors = cycle(['dimgray','darkgray', 'purple', 'thistle', 'navy', 'cornflowerblue', 'darkgreen',
            'lightgreen', 'firebrick', 'lightcoral', 'darkgoldenrod',
            'gold', 'saddlebrown', 'sandybrown'])

    args = docopt.docopt(__doc__)
    print(args)
    mid = args['--by']
    summary = args['--sum']
    relaxed = args['--relaxed']
    force = args['--force']
    delta = args['--delta']
    pthresh = float(args['--pthresh'])
    absthresh = float(args['--absthresh'])
    x = 'pre_' + mid + '_sum'
    if not relaxed:
        y = 'post_' + mid + '_sum'
    if relaxed:
        y = 'post_' + mid + '_relaxed_sum'

    '''
    else:
        bins = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 10.0]
        labels = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5]
        labels_str = ['0-0.1', '0.1-0.2', '0.2-0.4', '0.4-0.6',
                '0.6-0.8', '0.8-1.0', '1.0-1.5', '>1.5']
    '''
    df = summarize(args['<folders>'][0], summary=summary, force=force,
            relaxed=relaxed, by=mid)

    if args['--binby'] == 'rmsd':
        args['--binby'] = 'pre_rmsd_sum'
    if args['--binby'] == 'relsurf_sum':
        bins = [0, 0.05,10.0]
        labels = [0,0.05]
        labels_str = ['Buried','Exposed']

    elif args['--binby'] == 'dist':
        bins = [0, 0.1, 0.2, 0.4, 0.6,10.0]
        labels = [0, 0.1, 0.2, 0.4, 0.6]
        '''
        labels_str = ['0-0.1', '0.1-0.2', '0.2-0.4', '0.4-0.6',
                '0.6-0.8', '0.8-1.0', '1.0-1.2', '1.2-1.5', '>1.5']
        '''
        labels_str = ['0-0.1','0.1-0.2','0.2-0.4','0.4-0.6','>0.6']

    elif args['--binby'] == 'pre_rmsd_sum':
        bins = [0, 0.1, 0.2, 0.4, 0.6, 10]
        labels = [0, 0.1, 0.2, 0.4, 0.6]
        labels_str = ['0-0.1', '0.1-0.2', '0.2-0.4', '0.4-0.6', '>0.6']

    elif args['--binby'] == 'restype':
        labels = ['R','H','K','D','E','N','Q','S','T','C','G','P','A','I','L','M','F','W','Y','V']
        labels_str = labels

    elif args['--binby'] == 'resgroup':
        label_dict = {'R':'Charged','H':'Aro','K':'Charged','D':'Charged','E':'Charged',
        'N':'Polar','Q':'Polar','S':'Polar','T':'Polar',
        'C':'P, G, C', 'G':'P, G, C', 'P': 'P, G, C', 
        'A':'Hydrophobic', 'I':'Hydrophobic',
        'L':'Hydrophobic','M':'Hydrophobic','F':'Aro','W':'Aro','Y':'Aro','V':'Hydrophobic'}
        labels = ['Charged', 'Polar', 'P, G, C', 'Hydrophobic', 'Aro']
        labels_str = labels
    elif args['--binby'] == 'ss':
        label_dict = {'G': 'Helix', 'H': 'Helix', 'I': 'Helix','T': 'Helix', 
                'E':'β-sheet', 'B':'β-sheet', 
                'S': 'Loop', 'C': 'Loop', '_':'Loop'}
        labels = ['Helix', 'β-sheet', 'Loop']
        labels_str = labels

    elif args['--binby'] == 'dssp':
        labels = ['G', 'H', 'I', 'T', 'E', 'B', 'S', 'C']
        labels_str = labels
    elif args['--binby'] == 'none':
        labels=[0]
        labels_str = [0]

    else:
        stop = df[args['--binby']].max()
        start = df[args['--binby']].min()
        step = (stop - start) / float(args['--bins'])
        bins = np.arange(start, stop, step)
        print(bins)
        labels = bins[0:-1]
        labels_str = []
        for b in range(0, len(bins)-1):
            labels_str.append('{0:.3g}-{1:.3g}'.format(bins[b],bins[b+1]))
            if args['--bin']:
                if float(args['--bin']) < bins[b+1] and\
                        float(args['--bin']) > bins[b]:
                    args['--bin'] = bins[b]


    ind = np.arange(len(labels))
    width = 0.45
    fig, ax = plt.subplots()
    fig.set_size_inches(10,6)
    #fig = plt.figure(figsize=(10,4.8))

    # i keeps track of which folder we're on for formatting purposes
    i = 0
    ticklabels = []
    pdbs = []
    #pdbs = ['1a6g_1a6m','1b9k_1kyf','1bu7_1jme','1ffr_1edq','1hvf_1hve','1i4o_1kmc','1qop_1k8y','1qul_1quj']
    #pdbs = ['1ccp_3ccp','1ec0_1hvs','5eaa_1qit','7pti_1aal']
    pdbs=['1bri_1a2p', '1djc_3blm', '1kjh_1mt7']

    violin_labels = []
    def add_label(violin, label):
        color = violin["bodies"][0].get_facecolor().flatten()
        violin_labels.append((mpatches.Patch(color=color), label))

    for input_dir in args['<folders>']:
        folder_label = name_dict[input_dir.split('/')[-3]]
        df = summarize(input_dir, summary=summary, force=force,
                relaxed=relaxed, by=mid, threshold=pthresh,
                absthresh=absthresh)

        yvals_cst = []
        stderr_cst = []
        bincount_cst = []

        yvals_uncst = []
        stderr_uncst = []
        bincount_uncst = []

        print(df.columns)

        if args['--binby']=='none':
            df['binned']=0
        elif args['--binby'] and args['--binby']!='dist' and\
                args['--binby']!='restype' and args['--binby'] !=\
                'resgroup' and args['--binby'] != 'dssp' and\
                args['--binby'] != 'ss':
            df['binned'] = pd.cut(df[args['--binby']], bins=bins,
                    labels=labels)
        elif args['--binby']=='restype':
            df['binned'] = df['wt_restype_sum']
        elif args['--binby']=='resgroup':
            df['binned'] = df['wt_restype_sum'].map(label_dict)
        elif args['--binby']=='dssp':
            df['binned'] = df['ss_sum']
        elif args['--binby']=='ss':
            print(df['ss_sum'].unique())
            df['binned'] = df['ss_sum'].map(label_dict)
        else:
            df['binned'] = pd.cut(df[x], bins=bins, labels=labels)
            #df = df[df['binned'] != 1.0]
    
        by_dict = {'dist':'distance', 'rmsd':'RMSD'}
        sum_dict = {'low_score':'Lowest scoring', 'mean':'Mean',
                'median':'Median', 
                'percent_improved':'Percentage of improved',
                'structure':'Closest structure '}

        polar_only=False
        if polar_only:
            df = df[df['wt_restype_sum'].isin(['N','E','D','Q','K','R'])]
        big_only=False
        if big_only:
            df = df[df['wt_restype_sum'].isin(['F','Y','W','H'])]

        if args['--restype']:
            restypes = args['--restype'].split(',')
            print('Restricting to restype(s) ', args['--restype'])
            df = df[df['wt_restype_sum'].isin(restypes)]
            #df_uncst = df_uncst[df_uncst['wt_restype_sum']==args['--restype']]


        if args['--bin']:
            '''
            PLOT EXAMPLE PDBS
            '''
            import random
            num_examples = 10
            if args['--binby']=='restype':
                label = args['--bin']
            else:
                label = float(args['--bin'])
            df_cst = df[df['constrained']==True]
            df_uncst = df[df['constrained']==False]
            df_cst = df_cst[df_cst['binned']==label]
            print('Choosing {} examples out of {}'.format(num_examples,
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
                        ticklabels.append('{} ({}{}) to \n{} ({}{})'\
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
                        wt = df_cst.loc[cst_row, 'wt_sum']
                        wtchain = df_cst.loc[cst_row, 'wt_chain_sum']
                        wtresnum = df_cst.loc[cst_row, 'wt_resnum_sum']
                        mut = df_cst.loc[cst_row, 'mutant_sum']
                        mutchain = df_cst.loc[cst_row, 'mut_chain_sum']
                        mutresnum = df_cst.loc[cst_row, 'mut_resnum_sum']
                        ticklabels.append('{} ({}{}),\n{} ({}{})'\
                                .format(mut, mutchain, mutresnum, wt,
                                    wtchain, wtresnum))
                        
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
            width_split = 2 if args['--cst']=='both' else 1
            ax.set_xticks(ind + (width/width_split))
            ax.set_xticklabels(ticklabels)
            ax.set_title('{} examples from bin {}'.format(num_examples,
                label))
            plt.xticks(rotation=70)


        elif args['--violin']:
            '''
            VIOLIN PLOT OF SUMMARY DATA
            '''
            for label in labels:
                df_cst_binned = df[(df['constrained']==True) &
                        (df['binned']==label)]
                df_uncst_binned = df[(df['constrained']==False) &
                        (df['binned']==label)]

                nans = [float('nan'), float('nan')]

                cst_vals = df_cst_binned[y].values
                if len(cst_vals)==0:
                    yvals_cst.append(nans)
                else:
                    yvals_cst.append(cst_vals)
                bincount_cst.append(len(cst_vals))

                uncst_vals = df_uncst_binned[y].values
                if len(uncst_vals)==0:
                    yvals_uncst.append(nans)
                else:
                    yvals_uncst.append(uncst_vals)
                bincount_uncst.append(len(uncst_vals))

            widthmodifier = 0.15 if args['--cst']=='both' else 0.3
            widthmodifier=0.3
            widthlist = widthmodifier*len(args['<folders>']) * width/(len(args['<folders>']))
            if args['--cst']=='cst' or args['--cst']=='both':
                vplt = ax.violinplot(yvals_cst, ind + (2*i)*width /
                        (len(args['<folders>'])),
                        widths=widthlist, 
                        showmedians=True,
                        showmeans=True)
                # Coloring
                edgecolor = next(colors)
                for pc in vplt['bodies']:
                    pc.set_edgecolor(edgecolor)
                for partname in ('cbars','cmins','cmaxes','cmeans','cmedians'):
                    vp = vplt[partname]
                    vp.set_edgecolor(edgecolor)
                if args['--cst']=='both':
                    facecolor = edgecolor
                else:
                    facecolor = next(colors)
                for pc in vplt['bodies']:
                    pc.set_facecolor(facecolor)

                vplt['cmedians'].set_linestyles('dashed')

                add_label(vplt, folder_label + ' (constrained)')

                # Label number of bins
                for j, v in enumerate(yvals_cst):
                    ax.text(ind[j] + (2*i) * width /
                            len(args['<folders>']) - (0.5 *
                                width)/len(args['<folders>']),
                            0,
                            bincount_cst[j])

            if args['--cst']=='uncst' or args['--cst']=='both':
                vplt = ax.violinplot(yvals_uncst, ind + (2*i + 1)*width /
                        (len(args['<folders>'])),
                        widths=widthlist, 
                        showmedians=True,
                        showmeans=True)
                # Coloring
                edgecolor = next(colors)
                for pc in vplt['bodies']:
                    pc.set_edgecolor(edgecolor)
                for partname in ('cbars','cmins','cmaxes','cmeans','cmedians'):
                    vp = vplt[partname]
                    vp.set_edgecolor(edgecolor)
                if args['--cst']=='both':
                    facecolor = edgecolor
                else:
                    facecolor = next(colors)
                for pc in vplt['bodies']:
                    pc.set_facecolor(facecolor)

                vplt['cmedians'].set_linestyles('dashed')

                add_label(vplt, folder_label + ' (unconstrained)')

                # Label number of bins
                for j, v in enumerate(yvals_cst):
                    ax.text(ind[j] + (2*i + 1) * width /
                            len(args['<folders>']) - (0.5 *
                                width)/len(args['<folders>']),
                            0,
                            bincount_cst[j])
            ax.set_ylabel('Mean percentage of improved decoys')
            ax.set_xlabel('Binned {} (Å) of starting structure α-carbon to WT'\
                     .format(by_dict[args['--by']]))
            ax.set_title('Sampling methods broken down by {} of point mutant to WT'\
                    .format(by_dict[args['--by']]))
            width_split = 2 if args['--cst']=='both' else 1
            ax.set_xticks(ind + (width/width_split))
            ax.set_xticklabels(labels_str)

        else:
            '''
            BAR PLOT SUMMARY OF ALL DATA
            '''
            for label in labels:
                df_cst_binned = df[(df['constrained']==True) &
                        (df['binned']==label)]
                df_uncst_binned = df[(df['constrained']==False) &
                        (df['binned']==label)]

                if delta:
                    deltas_cst = df_cst_binned[y] - df_cst_binned[x]
                    deltas_uncst = df_uncst_binned[y] -\
                        df_uncst_binned[x]

                    #cst_mean = deltas_cst.mean()
                    cst_mean = deltas_cst.median()
                    stderr_cst_val = deltas_cst.std()
                    #uncst_mean = deltas_uncst.mean()
                    uncst_mean = deltas_uncst.median()
                    stderr_uncst_val = deltas_uncst.std()


                else:
                    cst_mean = df_cst_binned[y].mean()
                    stderr_cst_val = df_cst_binned[y].std()
                    uncst_mean = df_uncst_binned[y].mean()
                    stderr_uncst_val = df_uncst_binned[y].std()

                yvals_cst.append(cst_mean)
                stderr_cst.append(stderr_cst_val)
                bincount_cst.append(len(df_cst_binned.index))

                yvals_uncst.append(uncst_mean)
                stderr_uncst.append(stderr_uncst_val)
                bincount_uncst.append(len(df_uncst_binned.index))

            if args['--cst'] == 'cst' or args['--cst'] == 'both':
                ax.bar(ind + (2 * i) * width / (len(args['<folders>'])),
                        yvals_cst, width/(len(args['<folders>'])), #yerr=stderr_cst,
                        label=folder_label + '(constrained)',
                        color = next(colors))
                for j, v in enumerate(yvals_cst):
                    if np.isnan(v):
                        v = 0
                    if args['--sum']=='percent_improved':
                        v += 0.05
                    ax.text(ind[j] + (2*i) * width /
                            len(args['<folders>']) - (0.5 *
                                width)/len(args['<folders>']),
                            v ,
                            bincount_cst[j])
            if args['--cst'] == 'uncst' or args['--cst'] == 'both':
                ax.bar(ind + (2 * i + 1) * width / (len(args['<folders>'])),                
                        yvals_uncst, width/(len(args['<folders>'])), #yerr=stderr_uncst,
                        label=folder_label + '(unconstrained)',
                        color=next(colors))
                for j, v in enumerate(yvals_uncst):
                    if np.isnan(v):
                        v = 0
                    if args['--sum']=='percent_improved':
                        v += 0.05
                    ax.text(ind[j] + (2*i + 1) * width /
                            len(args['<folders>']) - (0.5 *
                                width)/len(args['<folders>']),
                            v,
                            bincount_uncst[j])

            ax.set_ylabel('Mean percentage of improved decoys')
            ax.set_xlabel('Binned {} (Å) of starting structure α-carbon to WT'\
                     .format(by_dict[args['--by']]))
            ax.set_title('Sampling methods broken down by {} of point mutant to WT'\
                    .format(by_dict[args['--by']]))
            width_split = 2 if args['--cst']=='both' else 1
            ax.set_xticks(ind + (width/width_split))
            ax.set_xticklabels(labels_str)


        i += 1

    if args['--violin']:
        #ax.legend(*zip(*violin_labels),loc=2)
        pass
    else:
        ax.legend()

    fig.tight_layout()
    plt.show()
