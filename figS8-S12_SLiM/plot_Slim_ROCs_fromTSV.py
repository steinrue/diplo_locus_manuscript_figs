#!/usr/bin/env python -W ignore::DeprecationWarning
"""collect LRs on each replicate & plot distances to the true locus"""
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

import sys, time, os, re, pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def get_ROCs(neut_data, sel_data):
    # make sure things are sorted, default is ascending
    neut_data = np.sort(neut_data)
    sel_data = np.sort(sel_data)
    TP, FP, TN, FN = [], [], [], []
    for thr in range(len(neut_data)):
        TP.append(np.sum(sel_data >= thr))
        FP.append(np.sum(neut_data >= thr))
        TN.append(np.sum(neut_data < thr))
        FN.append(np.sum(sel_data < thr))
    return TP, FP, TN, FN


def plot_ROCs_from_DF(ax, DF, which_s, true_s_2plot):
    ax.set_xlabel('False Positive Rate', fontsize="small")
    ax.set_ylabel('True Positive Rate', fontsize="small")
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.04, 1.04)
    ax.tick_params(labelsize=8)
    ax.axhline(0, 0, 1, ls=':', c='darkgray', lw=1)
    ax.axvline(0, 0, 1, ls=':', c='darkgray', lw=1)
    ax.axline((0, 0), (1, 1), ls=":", lw=0.5, c="0.5")

    line_colors = plt.cm.Dark2(range(len(true_s_2plot)))

    for idx, true_s in enumerate(true_s_2plot):
        # s_pair_key = (f'{float(H) * float(true_s):g}', true_s)
        neut_pool = np.array(DF[DF[which_s].astype(float) == 0].MLR.astype(float))
        sel_pool = np.array(DF[DF[which_s].astype(float) == float(true_s)].MLR.astype(float))
        # print(which_s, true_s, len(neut_pool), len(sel_pool))
        TP, FP, TN, FN = get_ROCs(neut_pool, sel_pool)
        # print(true_s, [TP[:3], FP[:3], TN[:3], FN[:3]])
        if np.any((np.array(TP) + np.array(FN)) == 0):
            TPR = np.zeros(np.array(TP).shape)
            denom = np.array(TP) + np.array(FN)
            indice = np.where(denom > 0)
            TPR[indice] = np.array(FP)[indice] / (np.array(FP)[indice] + np.array(TN)[indice])
            print(true_s, indice, TPR[:5])
        else:
            TPR = np.array(TP) / (np.array(TP) + np.array(FN))

        if np.any((np.array(FP) + np.array(TN)) == 0):
            FPR = np.zeros(np.array(FP).shape)
            denom = np.array(FP) + np.array(TN)
            indice = np.where(denom > 0)
            FPR[indice] = np.array(FP)[indice] / (np.array(FP)[indice] + np.array(TN)[indice])
            print(true_s, indice, FPR[:5])
        else:
            FPR = np.array(FP) / (np.array(FP) + np.array(TN))
        ax.plot(FPR, TPR, color=line_colors[idx], label=f'{float(true_s):.3f}')

    return ax


s_regex = re.compile(r's(\.\d+)')
H_regex = re.compile(r'_h(\.?\d+)_')


def plot_varS2_rocs_fromDF_wrapper(figname_prefix, DF, H_2plot, true_s_2plot, panel_labs, which_s: str, group: str):
    roc_figname = f'{figname_prefix}_ROCs.png'
    # three cols of ROCs
    roc_fig, roc_ax = plt.subplots(nrows=1, ncols=len(H_2plot), sharex=True, sharey=True,
                                   figsize=(2.25 * len(H_2plot), 2.75))

    for i, h in enumerate(H_2plot):
        DF_temp = DF[(DF.h.astype(float) == float(h))]  # | (DF.s2 == 0)
        print(i, h, DF_temp.shape)
        roc_ax[i] = plot_ROCs_from_DF(roc_ax[i], DF_temp,
                                      which_s=which_s, true_s_2plot=true_s_2plot)
        roc_ax[i].set_title(f'{panel_labs[i]})\t$h={h}$', loc='left', fontsize="medium", fontweight='bold')

    handles, labels = roc_ax[0].get_legend_handles_labels()
    # print(handles, labels)
    if which_s == 's2':
        var_label = r"True $s_{AA}$"
    else:
        assert which_s == 's1', which_s
        var_label = r"True $s_{Aa}$"

    plt.legend(handles, labels, loc="lower right",  # loc="best", bbox_to_anchor=(0.5, 0), borderaxespad=6,
               # 'center left',outsidencol=len(true_s_2plot), # mode='expand',
               edgecolor='None', facecolor='None', markerscale=1,
               fontsize="x-small", title_fontsize="small", frameon=False, title=var_label)
    roc_fig.suptitle(group, fontsize='medium', y=0.925)  #
    plt.tight_layout()  # h_pad=1.5
    roc_fig.savefig(roc_figname, dpi=350, transparent=True)
    print(roc_figname)


def plot_bothS_rocs_fromDF_wrapper(figname_prefix, s2_DF, s1_DF, H_2plot, true_s_2plot, panel_labs, suptitle):
    roc_figname = f'{figname_prefix}_ROCs.png'
    # sanity check
    assert len(H_2plot) == 2, H_2plot
    assert len(true_s_2plot) == 2, true_s_2plot
    num_panels = len(H_2plot[0]) + len(H_2plot[1])
    # prepare two sets of panel labels
    var_label = [r"True $s_{AA}$", r"True $s_{Aa}$"]
    DFs = [s2_DF, s1_DF]
    s_keys = ['s2', 's1']
    # three cols of ROCs
    roc_fig, roc_ax = plt.subplots(nrows=1, ncols=num_panels, sharex=True, sharey=True,
                                   figsize=(2.25 * num_panels, 2.75))
    fig_i_counter = 0
    for batch_i, thisS_H_2plot in enumerate(H_2plot):
        which_s = s_keys[batch_i]
        for h_i, h in enumerate(thisS_H_2plot):
            DF_temp = DFs[batch_i][(DFs[batch_i].h.astype(float) == float(h))]  # | (DF.s2 == 0)
            print(batch_i, h, DF_temp.shape)
            # print(DF_temp.h.unique(), DF_temp.s2.value_counts())
            roc_ax[fig_i_counter] = plot_ROCs_from_DF(roc_ax[fig_i_counter], DF_temp,
                                                      which_s=which_s, true_s_2plot=true_s_2plot[batch_i])
            roc_ax[fig_i_counter].set_title(f'{panel_labs[fig_i_counter]})\t$h={h}$',
                                            loc='left', fontsize="medium", fontweight='bold')
            if h_i + 1 < len(thisS_H_2plot):
                # remove legend
                roc_ax[fig_i_counter].legend().remove()
            else:
                handles, labels = roc_ax[fig_i_counter].get_legend_handles_labels()
                roc_ax[fig_i_counter].legend(handles, labels, loc="lower right",
                                             # loc="best", bbox_to_anchor=(0.5, 0), borderaxespad=6,
                                             edgecolor='None', facecolor='None', markerscale=1,
                                             fontsize="xx-small", title_fontsize="x-small", frameon=False,
                                             title=var_label[batch_i])
            # increment
            fig_i_counter += 1

    roc_fig.suptitle(suptitle, fontsize='medium', y=0.925)  #
    plt.tight_layout()  # h_pad=1.5
    roc_fig.savefig(roc_figname, dpi=350, transparent=True)
    print(roc_figname)


def main():
    # ABSL_PATH = '/gpfs/data/steinruecken-lab/XiaohengStuff/Diffusion_spect/Simulations/singleLocus/'  #
    s2_scoreFile = sys.argv[1]
    if not os.path.exists(s2_scoreFile):
        raise FileNotFoundError(s2_scoreFile)
    s1_scoreFile = sys.argv[2]
    if not os.path.exists(s1_scoreFile):
        raise FileNotFoundError(s1_scoreFile)
    stdVar = sys.argv[3]
    if stdVar == "sfs":
        sv_insert = '_stdVarSFS'
        vars2_slist = ['.0', '.001', '.002', '.005']  # , '.01', '.02', '.05'
        vars1_slist = ['.0002', '.0004', '.0006', '.0008', '.001']  # , '.01', '.02', '.05'
        simGroup = 'Selection on standing variation with f > 0.1, restart from after msprime'
    elif stdVar == "new":
        sv_insert = '_newSV'
        vars2_slist = ['.0', '.001', '.002', '.003', '.004', '.005']  # , '.01', '.02', '.05'
        vars1_slist = ['.0', '.0002', '.0004', '.0006', '.0008', '.001']  # , '.01', '.02', '.05'
        simGroup = 'Selection on standing variation'
    else:
        return False

    fig_prefix = sys.argv[4]

    vars2_Cond_list = [f's{s}_h{h}{sv_insert}'
                       for s in vars2_slist
                       for h in ['.5', '0', '1']
                       if f's{s}_h{h}' not in ['s.0_h0', 's.0_h1']]
    vars1_Cond_list = [f's{s}_h5{sv_insert}' for s in vars1_slist] + [f's.0_h.5{sv_insert}']

    if 'select' in s2_scoreFile:
        sites = 'LLR of selection target'
    else:
        assert 'max' in s2_scoreFile, s2_scoreFile
        sites = 'LLR at sites with maximal LLR'

    # read scores
    ## header: 'Condition', 's2', 'h', 'rep', 'position', 'ongrid_shat', 'shat', 'MLR'
    s2_scores = pd.read_csv(s2_scoreFile, sep="\t")
    # print(s2_scores.Condition.unique())
    s2_scores = s2_scores[s2_scores.Condition.isin(vars2_Cond_list)]
    # print(s2_scores.Condition.unique())
    s1_scores = pd.read_csv(s1_scoreFile, sep="\t")
    # print(s2_scores.Condition.unique())
    s1_scores = s1_scores[s1_scores.Condition.isin(vars1_Cond_list)]
    # print(s1_scores.Condition.value_counts())
    s1_scores['s1'] = s1_scores.s2.astype(float) * s1_scores.h.astype(float)
    # print(s1_scores.s1.value_counts())
    s1_scores.shat = s1_scores.shat.astype(float) * s1_scores.h.astype(float)

    # now plot ROCs
    plot_varS2_rocs_fromDF_wrapper(f'{fig_prefix}_200kb_t9x500gen_n40{sv_insert}_off-Grid',
                                   s2_scores, H_2plot=['0', '0.5', '1'],
                                   true_s_2plot=s2_scores.s2.unique(), panel_labs=["A", "B", "C"],
                                   which_s='s2', group=simGroup + '. ' + sites)
    # clear the slate
    plt.clf()
    plot_bothS_rocs_fromDF_wrapper(f'{fig_prefix}_200kb_t9x500gen_n40{sv_insert}_varS2+varS1_off-Grid',
                                   s2_scores, s1_scores, H_2plot=[['0', '0.5', '1'], ['5']],
                                   true_s_2plot=[s2_scores.s2.unique(), s1_scores.s1.unique()],
                                   panel_labs=["A", "B", "C", "D"],
                                   suptitle=simGroup + '. ' + sites)


if __name__ == '__main__':
    main()
