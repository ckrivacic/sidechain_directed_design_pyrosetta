"""
Usage:
    heatmap.py <folders>... [options]

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
import sys

args = docopt.docopt(__doc__)
if args['--binby'] == 'resgroup':
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

def row_index_from_pdbs(mut, wt, df):
    row = df.index[(df['mutant_sum'] == mut) & (df['wt_sum'] == wt)]
    return row.tolist()


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    import seaborn as sns

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
    #fig = plt.figure(figsize=(10,4.8))

    # i keeps track of which folder we're on for formatting purposes
    i = 0
    ticklabels = []
    pdbs = []
    #pdbs = ['1a6g_1a6m','1b9k_1kyf','1bu7_1jme','1ffr_1edq','1hvf_1hve','1i4o_1kmc','1qop_1k8y','1qul_1quj']
    #pdbs = ['1ccp_3ccp','1ec0_1hvs','5eaa_1qit','7pti_1aal']

    violin_labels = []
    def add_label(violin, label):
        color = violin["bodies"][0].get_facecolor().flatten()
        violin_labels.append((mpatches.Patch(color=color), label))

    df_cst_list = []
    df_uncst_list = []
    for input_dir in args['<folders>']:
        print('Parsing dataframes from {}'.format(input_dir))
        folder_label = name_dict[input_dir.split('/')[-3]]
        df = summarize(input_dir, summary=summary, force=force,
                relaxed=relaxed, by=mid, threshold=pthresh,
                absthresh=absthresh)

        df_cst = df[(df['constrained']==True)]
        df_cst_list.append(df_cst)
        df_uncst = df[(df['constrained']==False)]
        df_uncst_list.append(df_uncst)
        i += 1

    df_cst_final = pd.concat(df_cst_list, ignore_index=True)
    df_cst_final['pdbs'] = df_cst_final['mutant_sum'] + '_' + df_cst_final['wt_sum']
    df_uncst_final = pd.concat(df_uncst_list, ignore_index=True)
    df_uncst_final['pdbs'] = df_uncst_final['mutant_sum'] + '_' + df_uncst_final['wt_sum']

    if args['--binby']=='resgroup':
        df_cst_final['binned'] = df_cst_final['wt_restype_sum'].map(label_dict)
    elif args['--binby']=='ss':
        print(df['ss_sum'].unique())
        df_cst_final['binned'] = df_cst_final['ss_sum'].map(label_dict)

    def reshape_dataframe(df, metric):
        pivot = pd.pivot_table(df,
                values=metric,index=['pdbs','binned'],columns=['mover_sum'],
                dropna=True)
        return pivot

    cst_pivot = reshape_dataframe(df_cst_final,'post_dist_sum')
    #uncst_pivot = reshape_dataframe(df_uncst_final,'post_dist_sum')
    '''
    good_rows = []    
    for idx, row in cst_pivot.iterrows():
        print(row)
        if not row.isnull().values.any():
            good_rows.append(row)
    '''
    cst_pivot = cst_pivot.dropna(how='any')
    
    #cst_pivot.stack().reset_index()
    #print(cst_pivot.iloc[:, df.columns.get_level_values(0)=='binned'])
    print(cst_pivot.xs('binned',level='Col',axis=1))
    species = cst_pivot.pop('binned')
    lut = dict(zip(species.unique(), 'rgb'))
    print(lut)
    row_colors = species.map(lut)

    cst_pivot.dropna(how='any')
    ax = sns.clustermap(cst_pivot, row_colors=row_colors)
    plt.show()
