#!/usr/bin/env python -W ignore::DeprecationWarning
"""plot distances from maxLR sites to the true locus"""
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import sys, os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def read_scores(sel_scoreFile, max_scoreFile, which_s, H_pool, varS_Cond_list):
    # sanity check
    if not os.path.exists(sel_scoreFile):
        raise FileNotFoundError(sel_scoreFile)
    if not os.path.exists(max_scoreFile):
        raise FileNotFoundError(max_scoreFile)
    # read scores
    ## header: 'Condition', 's2', 'h', 'rep', 'position', 'ongrid_shat', 'shat', 'MLR'
    sel_scores = pd.read_csv(sel_scoreFile, sep="\t")
    sel_scores = sel_scores[['Condition', 's2', 'h', 'rep', 'position', 'shat', 'MLR']]
    sel_scores.position = sel_scores.position.astype(int)

    max_scores = pd.read_csv(max_scoreFile, sep="\t")[['Condition', 's2', 'h', 'rep', 'position', 'shat', 'MLR']]
    max_scores.position = max_scores.position.astype(int)

    print(sel_scores.shape, max_scores.shape)

    sel_scores = sel_scores[sel_scores.Condition.isin(varS_Cond_list) &
                            sel_scores.h.astype(str).isin(H_pool)]

    max_scores = max_scores[max_scores.Condition.isin(varS_Cond_list) &
                            max_scores.h.astype(str).isin(H_pool)]

    # integrate the two df
    print(sel_scores.shape, max_scores.shape)
    merged_scores = sel_scores.merge(max_scores, how='inner', on=['Condition', 's2', 'h', 'rep'],
                                     # left_on=['Condition', 's2', 'h', 'rep', 'position', 'shat', 'MLR'],
                                     # right_on=['Condition', 's2', 'h', 'rep', 'position', 'shat', 'MLR'],
                                     suffixes=['_true', '_maxLR'])
    if which_s == 's1':
        merged_scores['s1'] = merged_scores.s2.astype(float) * merged_scores.h.astype(float)
    else:
        assert which_s == 's2', f'Unrecognized `which_s`: {which_s}'

    # sanity check
    assert which_s in merged_scores.columns, f'\"{which_s}\" not among the columns: {merged_scores.columns}'
    merged_scores['true_s'] = merged_scores[which_s].apply('{:.3f}'.format)

    print(merged_scores.shape)

    return merged_scores


def _express_H(H):
    Dom = {".5": "additive", "0.5": "additive",
           '0': "complete recessive", '0.0': "complete recessive",
           "1": "complete dominance", "1.0": "complete dominance",
           "5": "over-dominance", "5.0": "over-dominance",
           "Inf": r"$s_{AA}=0$"}
    if H == "Inf":
        return Dom[H]
    elif float(H) > 1:
        return f'h = {float(H):g}, {Dom[H]}\n' + r'(grouped by $s_{Aa}$ values)'
    else:
        return f'h = {float(H):g}, {Dom[H]}'


def plot_dist_boxplots(merged_scores, fig_prefix, H_colors, var_labels, fig_title,
                       add_dots=False, point_colors=None):
    # now we plot
    fig, bp = plt.subplots(figsize=(7, 3.5), sharey=True)

    # set up frame
    sns.set_theme(style='whitegrid')  # , rc=Theme_params
    bp.grid(visible=True, axis='x', lw=0.3)
    bp.axhline(0, xmin=0, xmax=1, lw=0.8)
    # plot
    bp = sns.boxplot(data=merged_scores, x='true_s', y='dist', hue='H',
                     linewidth=0.5, palette=H_colors, dodge=True,  # facealpha=0.7,
                     fliersize=0, orient="v", ax=bp)

    if add_dots:
        assert point_colors is not None, 'Need color map for stripplot'
        bp = sns.swarmplot(data=merged_scores, x='true_s', y='dist', hue='H',
                           palette=point_colors, size=1.2, dodge=True, orient='v', ax=bp)
    # axes and legend
    bp.set_ylabel('Distance to Target (kbp)', fontsize="small")
    ## transform (bp --> kbp) ticklabels
    yticks = bp.get_yticks()
    bp.set_yticks(yticks, labels=[f'{p/1e3:g}' for p in yticks])
    bp.set_xlabel(f'{var_labels[0]} ( or t{var_labels[1][1:]} for h=5)', fontsize="small")
    bp.tick_params(axis='both', labelsize="x-small")
    handles, labels = bp.get_legend_handles_labels()
    # print(handles, labels)
    # maybe use shorter labels
    # labels = [_express_H(h) for h in labels]
    # because both violin & dot plots generate legend
    n_legendlabs = int(len(handles) // 2)
    lg = bp.legend(handles[:n_legendlabs], labels[:n_legendlabs], loc="best", title=r'Dominance $h$',  # "lower right"
                   # loc="best", bbox_to_anchor=(0.5, 0), borderaxespad=6, facecolor='None',
                   edgecolor='None', markerscale=2, framealpha=0.6, ncol=2,
                   fontsize="xx-small", title_fontsize="x-small", frameon=True)
    for t in lg.get_texts():
        t.set_va('center')  # baseline_'bottom'
    bp.set_title(fig_title, fontsize='small')
    fig.tight_layout()
    # figname = f'{fig_prefix}_Dist_200kb_t9x500gen_n40_varS1+varS2_Unif_offGrid-maxLLR_boxplots.png'
    figname = f'{fig_prefix}_Dist_200kb_t9x500gen_n40_varS1+varS2_Unif_offGrid-maxLLR_boxplots.pdf'
    fig.savefig(figname, dpi=500, transparent=True)
    print(figname)


def main():
    fig_prefix = sys.argv[1]
    vars2_slist = ['.0', '.001', '.002', '.003', '.004', '.005']
    vars1_slist = ['.0', '.0002', '.0004', '.0006', '.0008', '.001']
    # sel_scoreFile_s2 = 'DL-varS2_selected-sites_NotLost_minMAF.05_200kb_t9x500gen_n40_3h6s_Concatenated_offGrid_Dist.tsv.gz'
    # max_scoreFile_s2 = 'DL-varS2_maxLR-sites_NotLost_minMAF.05_200kb_t9x500gen_n40_3h6s_Concatenated_offGrid_Dist.tsv.gz'
    # sel_scoreFile_s1 = 'DL-varS1_selected-sites_NotLost_minMAF.05_200kb_t9x500gen_n40_2h6s_Concatenated_offGrid_Dist.tsv.gz'
    # max_scoreFile_s1 = 'DL-varS1_maxLR-sites_NotLost_minMAF.05_200kb_t9x500gen_n40_2h6s_Concatenated_offGrid_Dist.tsv.gz'
    sel_scoreFile_s2, max_scoreFile_s2, sel_scoreFile_s1, max_scoreFile_s1 = sys.argv[2:]

    # some plotting parameters
    # H_2plot = [['0.0', '0.5', '1.0'], ['5.0']]
    H_2plot = [['0.0', '0.5', '1.0'], ['5']]
    # true_s_2plot = [vars2_slist, vars1_slist]
    H_colors = {".5": plt.cm.tab20(0*2+1),  # blue
                "0.5": plt.cm.tab20(0*2+1),  # in case it doesn't recognize
                "0": plt.cm.tab20(1*2+1),  # orange
                "0.0": plt.cm.tab20(1*2+1),  # orange
                "1": plt.cm.tab20(2*2+1),  # green
                "1.0": plt.cm.tab20(2*2+1),  # green
                "5": plt.cm.tab20(3*2+1),  # red
                "5.0": plt.cm.tab20(3*2+1),  # red
                "Inf": plt.cm.tab20(4*2+1),  # purple
                "fix_s1": plt.cm.tab20(5*2+1)}  # brown

    dot_colors = {".5": plt.cm.tab20(0*2),  # blue
                  "0.5": plt.cm.tab20(0*2),  # in case it doesn't recognize
                  "0": plt.cm.tab20(1*2),  # orange
                  "0.0": plt.cm.tab20(1*2),  # orange
                  "1": plt.cm.tab20(2*2),  # green
                  "1.0": plt.cm.tab20(2*2),  # green
                  "5": plt.cm.tab20(3*2),  # red
                  "5.0": plt.cm.tab20(3*2),  # red
                  "Inf": plt.cm.tab20(4*2),  # purple
                  "fix_s1": plt.cm.tab20(5*2)}  # brown

    var_labels = [r"True $s_{AA}$", r"True $s_{Aa}$"]

    # vars2_Cond_list = [f's{s}_h{h}' for s in vars2_slist
    #                    for h in ['.5', '0', '1'] if f's{s}_h{h}' not in ['s.0_h0', 's.0_h1']]
    vars2_Cond_list = [f's{s}_h{h}' for s in vars2_slist for h in ['.5', '0', '1']]
    # vars1_Cond_list = [f's{s}_h5' for s in vars1_slist] + [f's.0_h.5']
    vars1_Cond_list = [f's{s}_h5' for s in vars1_slist]

    # read var_s2
    merged_scores_s2 = read_scores(sel_scoreFile_s2, max_scoreFile_s2, 's2', H_2plot[0], vars2_Cond_list)

    # read var_s1
    merged_scores_s1 = read_scores(sel_scoreFile_s1, max_scoreFile_s1, 's1', H_2plot[1], vars1_Cond_list)

    # now merge s1 & s2
    merged_scores = merged_scores_s2.merge(merged_scores_s1, suffixes=['_varS2', '_varS1'], how='outer',
                                           on=['true_s', 'h', 'rep', 'position_true', 'position_maxLR']
                                           )
    # get distance
    merged_scores['dist'] = merged_scores.position_maxLR.astype(int) - merged_scores.position_true.astype(int)
    # just checkin'
    print("Range of dist:", merged_scores.dist.min(), merged_scores.dist.max())
    # for plotting convenience
    merged_scores['H'] = merged_scores.h.apply('{:g}'.format).astype(str)

    # now plot
    fig_title = 'Distance of highest LLR site from selection target'
    plot_dist_boxplots(merged_scores, fig_prefix, H_colors, var_labels, fig_title,
                       add_dots=True, point_colors=dot_colors)


if __name__ == '__main__':
    main()
