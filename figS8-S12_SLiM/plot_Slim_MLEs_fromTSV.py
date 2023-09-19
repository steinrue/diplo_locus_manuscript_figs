#!/usr/bin/env python -W ignore::DeprecationWarning
"""collect LRs on each replicate & plot MLE distribution"""
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

import sys, time, os, re, pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


s_regex = re.compile(r's(\.\d+)')
H_regex = re.compile(r'_h(\.?\d+)_')


def annotate_outliers(DF, ax, which_s, cond_list, outlier_bounds=None, note_outliers: bool=False):
    if outlier_bounds is not None:
        # check that it's a tuple of two numbers
        assert (len(outlier_bounds) == 2) and isinstance(outlier_bounds[0], float), \
            f'Wrong format for outlier bounds: {outlier_bounds}. Must be a touple of two ascending floats.'
        lower, upper = outlier_bounds
        top_outliers = DF[(DF[f"shat"] > upper)]
        bottom_outliers = DF[(DF[f"shat"] < lower)]
    else:
        # collect everything with abs. larger than 0.1
        top_outliers = DF[DF.shat > 0.1]
        bottom_outliers = DF[DF.shat < -0.1]
        DF_subset = DF[(DF.shat >= -0.1) & (DF.shat <= 0.1)]
        lower, upper = DF_subset.shat.min(), DF_subset.shat.max()
    # set y range
    ax.set_ylim(lower, upper)

    ymin, ymax = ax.get_ylim()
    # annotate their counts on the plot
    if note_outliers:
        ## find the x coords of all the vertical lines, which are the central axes of each box
        # print( [set(l.get_xdata()) for l in ax.get_lines()])
        all_x_coords = [l.get_xdata()[0] for l in ax.get_lines() if len(set(l.get_xdata())) == 1]
        all_x_coords = sorted(list(set(all_x_coords)))
        ### sanity check
        assert len(all_x_coords) == len(cond_list)
        ### match them with H & s; might be the stupidest way but as long as it works...
        matched_xCoords = []
        for i, cond in enumerate(cond_list):
            s = float(s_regex.findall(cond)[0])
            H = float(H_regex.findall(cond)[0])
            matched_xCoords.append([s, H, all_x_coords[i]])
        matched_xCoords = pd.DataFrame(matched_xCoords, columns=[which_s, 'h', 'xCoord'])
        ### assign H & s to these numbers
        ## top
        if top_outliers.shape[0] > 0:
            top_outliers_counts = pd.pivot_table(top_outliers, values='rep',
                                                 index=[which_s, 'h'],
                                                 aggfunc=len).reset_index()
            top_outliers_counts = top_outliers_counts.rename(columns={'rep': 'counts'})
            # merge to get their xcoords
            top_outliers_counts = top_outliers_counts.merge(matched_xCoords, how='left',
                                                            on=[which_s, 'h'])
            # now annotate
            # print(top_outliers_counts)
            if top_outliers_counts.shape[0] > 5:
                noise = np.random.normal(0, 2.5e-4, top_outliers_counts.shape[0])
                noise = np.abs(noise)
            else:
                noise = [0] * top_outliers_counts.shape[0]
            for i, row in top_outliers_counts.iterrows():
                # row: true_s, H, xCoord, count
                tag = row.counts.astype(int)
                # print(tag)
                this_x = row.xCoord
                ax.annotate(text=str(tag), xy=(this_x, ymax),
                            xytext=(this_x - noise[i]*(-1)**i, ymax + 0.14*(ymax - ymin) + noise[i]*(-1)**i ),
                            arrowprops=dict(arrowstyle="->", color="#66666622"), fontsize='x-small',
                            bbox=dict(boxstyle="round", pad=0.15, edgecolor='#aaaaaa66', fc="#ffffff88"), ha='center', va='center',
                            color='black', annotation_clip=True)
        # bottom
        if bottom_outliers.shape[0] > 0:
            bottom_outliers_counts = pd.pivot_table(bottom_outliers, values='rep',
                                                    index=[which_s, 'h'],
                                                    aggfunc=len).reset_index()
            bottom_outliers_counts = bottom_outliers_counts.rename(columns={'rep': 'counts'})
            # print(bottom_outliers_counts, bottom_outliers)
            # merge to get their xcoords
            bottom_outliers_counts = bottom_outliers_counts.merge(matched_xCoords, how='left',
                                                                  on=[which_s, 'h'])
            labels = r''
            # now annotate
            for i, row in bottom_outliers_counts.iterrows():
                # row: s, h, xCoord, count
                tag = row.counts.astype(int)
                # print(tag)
                this_x = row.xCoord
                ax.annotate(text=str(tag), xy=(this_x, ymin),
                            xytext=(this_x, ymin - 0.14*(ymax - ymin)),  # + noise[i]*(-1)**i
                            arrowprops=dict(arrowstyle="->", color="#66666622"), fontsize='x-small',
                            bbox=dict(boxstyle="round", pad=0.125, edgecolor='#aaaaaa66', fc='#ffffff88'), ha='center', va='center',
                            color='black', annotation_clip=True)
    return ax, top_outliers, bottom_outliers


def main():
    ABSL_PATH = '/gpfs/data/steinruecken-lab/XiaohengStuff/Diffusion_spect/Simulations/singleLocus/'  #
    s2_scoreFile = sys.argv[1]
    if not os.path.exists(s2_scoreFile):
        raise FileNotFoundError(s2_scoreFile)

    stdVar = sys.argv[2]
    if stdVar == "sfs":
        sv_insert = '_stdVarSFS'
        vars2_slist = ['.0', '.001', '.002', '.005']  #, '.02', '.05', '.01'
        vars1_slist = ['.0002', '.0004', '.0006', '.0008', '.001']  #, '.01', '.02', '.05'
    elif stdVar == "new":
        sv_insert = '_newSV'
        vars2_slist = ['.0', '.001', '.002', '.003', '.004', '.005']  #, '.01', '.02', '.05'
        vars1_slist = ['.0002', '.0004', '.0006', '.0008', '.001']  #, '.01', '.02', '.05'
    else:
        return False

    fig_prefix = sys.argv[3]

    # unless otherwise specified
    on_off_grid = 'off'

    h_list = ['0', '.5', '1', '5']
    vars2_Cond_list = [f's{s}_h{h}{sv_insert}'
                       for s in vars2_slist
                       for h in ['.5', '0', '1']
                       if f's{s}_h{h}' not in ['s.0_h0', 's.0_h1']]
    # vars1_Cond_list = [f's{s}_h5{sv_insert}' for s in vars1_slist]  # + [f's.0_h.5{sv_insert}']

    H_colors = {".5": plt.cm.tab10(0),  # blue
                "0.5": plt.cm.tab10(0),  # in case it doesn't recognize
                "0": plt.cm.tab10(1),  # orange
                "0.0": plt.cm.tab10(1),  # orange
                "1": plt.cm.tab10(2),  # green
                "1.0": plt.cm.tab10(2),  # green
                "5": plt.cm.tab10(3),  # red
                "5.0": plt.cm.tab10(3),  # red
                "Inf": plt.cm.tab10(4),
                "fix_s1": plt.cm.tab10(5)}  # purple

    Dom = {".5": "additive", "0.5": "additive",
           '0': "complete recessive", '0.0': "complete recessive",
           "1": "complete dominance", "1.0": "complete dominance",
           "5": "over-dominance", "5.0": "over-dominance",
           "Inf": r"$s_{AA}=0$"}

    def _express_H(H):
        if H == "Inf":
            return Dom[H]
        else:
            return f'h = {float(H):g}, {Dom[H]}'


    # read scores
    ## header: 'Condition', 's2', 'h', 'rep', 'position', 'ongrid_shat', 'shat', 'MLR'
    s2_scores = pd.read_csv(s2_scoreFile, sep="\t")
    s2_scores = s2_scores[s2_scores.Condition.isin(vars2_Cond_list)]
    # print(s2_scores.h.unique(), type(s2_scores.h.unique()[0]))
    # s1_scores = pd.read_csv(s1_scoreFile, sep="\t")
    # s1_scores = s1_scores[s1_scores.Condition.isin(vars1_Cond_list)]
    # s1_scores['s1'] = s1_scores.s2.astype(float) * s1_scores.h.astype(float)
    # s1_scores.shat = s1_scores.shat.astype(float) * s1_scores.h.astype(float)

    # now plot MLEs
    box_fig, ax1 = plt.subplots(figsize=(6, 3.5))
    ax1 = sns.boxplot(x='s2', y='shat', hue=s2_scores.h.astype(str), data=s2_scores, palette=H_colors,  # "Pastel1",
                      fliersize=0.5, orient="v", linewidth=0.5, ax=ax1)

    # aux lines
    ax1.grid(visible=True, axis='y', lw=0.3)
    x_left, x_right = ax1.get_xlim()
    s2_list_2plot = sorted(s2_scores.s2.astype(float).unique())
    ax1.hlines(y=s2_list_2plot, xmin=x_left, xmax=x_right,
               colors=['darkred'] * len(s2_list_2plot), linestyle='dashed', lw=0.5)  #
    # annotate outliers
    ax1, s2_top_outliers, s2_bottom_outliers = annotate_outliers(s2_scores, ax1, 's2', vars2_Cond_list,
                                                                 outlier_bounds=(-0.015, 0.015), note_outliers=True)

    # double check axis labels
    ax1.set_ylabel(None)
    ax1.set_title('Inferred $~\widehat{s}_{AA}$', loc='left', fontsize='medium') #, fontweight='bold'
    ax1.set_xlabel(r'True $s_{AA}$', fontsize='medium')
    ax1.tick_params(axis='both', labelsize="x-small")
    # get legend handles
    handles_L, labels_L = ax1.get_legend_handles_labels()
    labels_L = [_express_H(h) for h in labels_L]
    ax1.legend(handles_L, labels_L, framealpha=0, fontsize='x-small', frameon=False,
                   title=r'Dominance', title_fontsize='small', loc='best')
    ymin, ymax = ax1.get_ylim()
    print(ymin, ymax)

    plt.tight_layout()
    # figname = f'{ABSL_PATH}DL-1D_200kb_t9n40_noS001_center10kb_{on_off_grid}Grid-maxLLR_loc_boxplots.png'

    figname1 = f'{ABSL_PATH}{fig_prefix}_200kb_t9x500gen_n40{sv_insert}_vars2_{on_off_grid}Grid_MLEs_boxplots.png'
    # figname = f'{ABSL_PATH}DL-1D_200kb_t9n40_{on_off_grid}Grid-maxLLR_loc_violinplots.png' #
    box_fig.savefig(figname1)


if __name__ == '__main__':
    main()
