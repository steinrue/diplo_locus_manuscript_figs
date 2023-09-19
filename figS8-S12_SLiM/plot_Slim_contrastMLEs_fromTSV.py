#!/usr/bin/env python -W ignore::DeprecationWarning
"""collect LRs on each replicate & plot distances to the true locus"""
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
        top_outliers = DF[(DF.shat > upper)]
        bottom_outliers = DF[(DF.shat < lower)]
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
        # assert len(all_x_coords) == len(cond_list)
        assert len(all_x_coords) == len(DF.s2.unique()) * len(DF.h.unique()), \
            f'len(all_x_coords) = {len(all_x_coords)}, len(DF.s2.unique()) = {len(DF.s2.unique())},\n' \
            f' len(DF.h.unique()) = {len(DF.h.unique())}, len(cond_list) = {len(cond_list)}'
        ### match them with H & s; might be the stupidest way but as long as it works...
        matched_xCoords = []
        # print(f'len(all_x_coords) = {len(all_x_coords)}, len(cond_list) = {len(cond_list)}\n{cond_list}')
        # for i, cond in enumerate(cond_list):
        #     s = float(s_regex.findall(cond)[0])
        #     H = float(H_regex.findall(cond)[0])
        counter = 0
        for s in DF.s2.unique():
            for H in DF.h.unique():
                matched_xCoords.append([s, H, all_x_coords[counter]])
                counter += 1
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
                            xytext=(this_x - noise[i]*(-1)**i, ymax + 0.15*(ymax - ymin) + noise[i]*(-1)**i ),
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
                            xytext=(this_x, ymin - 0.2*(ymax - ymin)),  # + noise[i]*(-1)**i
                            arrowprops=dict(arrowstyle="->", color="#66666622"), fontsize='x-small',
                            bbox=dict(boxstyle="round", pad=0.125, edgecolor='#aaaaaa66', fc='#ffffff88'), ha='center', va='center',
                            color='black', annotation_clip=True)
    return ax, top_outliers, bottom_outliers


def main():
    ABSL_PATH = '/gpfs/data/steinruecken-lab/XiaohengStuff/Diffusion_spect/Simulations/singleLocus/'  #
    s2_scoreFile = sys.argv[1]
    if not os.path.exists(s2_scoreFile):
        raise FileNotFoundError(s2_scoreFile)
    s1_scoreFile = sys.argv[2]
    if not os.path.exists(s1_scoreFile):
        raise FileNotFoundError(s1_scoreFile)

    stdVar = sys.argv[3]
    if stdVar == "sfs":
        sv_insert = '_stdVarSFS'
        vars2_slist = ['.0', '.001', '.002', '.005']
        vars1_slist = ['.0002', '.0004', '.0006', '.0008', '.001']
    elif stdVar == "new":
        sv_insert = '_newSV'
        vars2_slist = ['.0', '.001', '.002', '.003', '.004', '.005']
        vars1_slist = ['.0', '.0002', '.0004', '.0006', '.0008', '.001']
    else:
        return False

    fig_prefix = sys.argv[4]

    # unless otherwise specified
    on_off_grid = 'off'

    # h_list = ['0', '.5', '1', '5']
    vars2_Cond_list = [f's{s}_h{h}{sv_insert}'
                       for s in vars2_slist
                       for h in ['.5', '0', '1']
                       if f's{s}_h{h}' not in ['s.0_h0', 's.0_h1']]
    vars1_Cond_list = [f's{s}_h5{sv_insert}' for s in vars1_slist]  # + [f's.0_h.5{sv_insert}']

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
    s1_scores = pd.read_csv(s1_scoreFile, sep="\t")
    s1_scores = s1_scores[s1_scores.Condition.isin(vars1_Cond_list)]
    s1_scores['s1'] = s1_scores.s2.astype(float) * s1_scores.h.astype(float)
    s1_scores.shat = s1_scores.shat.astype(float) * s1_scores.h.astype(float)

    # now plot MLEs
    box_fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(7, 3.5), sharey=True,
                                       gridspec_kw={'width_ratios': (0.5, 0.25)})

    if 'selected' in s2_scoreFile:
        ax1 = sns.boxplot(x='s2', y='shat', hue=s2_scores.h.astype(str), data=s2_scores, palette=H_colors,  # "Pastel1",
                          fliersize=0.5, orient="v", linewidth=0.5, ax=ax1)
        ax1.set_title(r'A) Inferred $~\widehat{s}_{AA}$', loc='left', fontsize='medium') #, fontweight='bold'
        # annotate outliers
        ax1, s2_top_outliers, s2_bottom_outliers = annotate_outliers(s2_scores, ax1, 's2', vars2_Cond_list,
                                                                     outlier_bounds=(-0.01, 0.015), note_outliers=True)
    else:
        assert 'max' in s2_scoreFile, s2_scoreFile

        ax1 = sns.boxplot(x='s2', y=s2_scores.shat.abs(), hue=s2_scores.h.astype(str), data=s2_scores,
                          palette=H_colors, fliersize=0.5, orient="v", linewidth=0.5, ax=ax1 )
        ax1.set_title(r'A) Inferred $~|\widehat{s}_{AA}|$ of top LR sites', loc='left', fontsize='medium') #, fontweight='bold'
        # annotate outliers
        s2_scores.shat = s2_scores.shat.abs()
        ax1, s2_top_outliers, s2_bottom_outliers = annotate_outliers(s2_scores, ax1, 's2', vars2_Cond_list,
                                                                     outlier_bounds=(-0.005, 0.015), note_outliers=True)
    # aux lines
    ax1.grid(visible=True, axis='y', lw=0.3)
    x_left, x_right = ax1.get_xlim()
    s2_list_2plot = sorted(s2_scores.s2.astype(float).unique())
    print(s2_list_2plot)
    ax1.hlines(y=s2_list_2plot, xmin=x_left, xmax=x_right,
               colors=['darkred'] * len(s2_list_2plot), linestyle='dashed', lw=0.5)  #

    # double check axis labels
    ax1.set_ylabel(None)
    ax1.set_xlabel(r'True $s_{AA}$', fontsize='medium')
    ax1.tick_params(axis='both', labelsize="x-small")
    # get legend handles
    handles_L, labels_L = ax1.get_legend_handles_labels()
    labels_L = [_express_H(h) for h in labels_L]
    ax1.legend(handles_L, labels_L, framealpha=0, fontsize='x-small', frameon=False,
               loc='best')  # title=r'Dominance', title_fontsize='small',
    ymin, ymax = ax1.get_ylim()
    print(ymin, ymax)

    # move on to s1
    # aux lines
    ax2.grid(visible=True, axis='y', lw=0.3)
    if 'selected' in s1_scoreFile:
        ax2 = sns.boxplot(x='s1', y='shat', hue=s1_scores.h.astype(str), data=s1_scores, palette=H_colors,  # "Pastel1",
                          fliersize=0.5, orient="v", linewidth=0.5, ax=ax2)
        ax2.set_title(r'B) Inferred $~\widehat{s}_{Aa}$', loc='left', fontsize='medium') #, fontweight='bold'
        # annotate outliers
        ax2, s1_top_outliers, s2_bottom_outliers = annotate_outliers(s1_scores, ax2, 's1', vars1_Cond_list,
                                                                     outlier_bounds=(-0.015, 0.015), note_outliers=True)
    else:
        assert 'max' in s1_scoreFile, s1_scoreFile
        ax2 = sns.boxplot(x='s2', y=s1_scores.shat.abs(), hue=s1_scores.h.astype(str), data=s1_scores,
                          palette=H_colors, fliersize=0.5, orient="v", linewidth=0.5, ax=ax2)
        ax2.set_title(r'B) Inferred $~|\widehat{s}_{Aa}|$ of top LR sites', loc='left', fontsize='medium') #, fontweight='bold'
        # annotate outliers
        s1_scores.shat = s1_scores.shat.abs()
        ax2, s1_top_outliers, s1_bottom_outliers = annotate_outliers(s1_scores, ax2, 's1', vars1_Cond_list,
                                                                     outlier_bounds=(-0.005, 0.015), note_outliers=True)
    x_left, x_right = ax2.get_xlim()
    s1_list_2plot = sorted(s1_scores.s1.unique().astype(float).tolist())
    print(s1_list_2plot)
    ax2.hlines(y=s1_list_2plot, xmin=x_left, xmax=x_right,
               colors=['darkred'] * len(s1_list_2plot), linestyle='dashed', lw=0.5)  #
    # outlier_bounds=None,
    # double check axis labels
    ax2.set_ylabel(None)
    ax2.set_xlabel(r'True $s_{Aa}$', fontsize='medium')
    handles_R, labels_R = ax2.get_legend_handles_labels()
    labels_R = [_express_H(h) for h in labels_R]
    # print(labels_R)
    # some more tick label configuring
    # xticks = ax2.get_xticklabels()
    # print(xticks)
    ax2.set_xticklabels(labels=[f'{s1:.3f}' for s1 in s1_list_2plot])
    ax2.tick_params(axis='both', labelsize="x-small")
    ax2.set_ylim(ymin, ymax)

    # add shared legend & title <- nope
    # handles, labels = list(handles_L) + list(handles_R), list(labels_L) + list(labels_R)
    ax2.legend(handles_R, labels_R, framealpha=0, fontsize='x-small', frameon=False,
               loc='best')  # title='Dominance', title_fontsize='small',

    # box_fig.suptitle() <-- no need for title
    plt.tight_layout()

    figname1 = f'{ABSL_PATH}{fig_prefix}_200kb_t9x500gen_n40{sv_insert}_{on_off_grid}Grid_MLEs_contrastBoxes.png'
    # figname = f'{ABSL_PATH}DL-1D_200kb_t9n40_{on_off_grid}Grid-maxLLR_loc_violinplots.png' #
    box_fig.savefig(figname1, dpi=350)


if __name__ == '__main__':
    main()
