import numpy as np
import pandas as pd
import pickle, sys
import seaborn as sns
import pathlib
import matplotlib.pyplot as plt
import bz2




MY_PALETTE = [plt.cm.tab10(0), plt.cm.tab10(2)]
ONSET_RESULTS_PICKLE_FILE = pathlib.Path ("../supp_figS8_onset/final_lls.pickle.bz2")




def loadPickledResults (pickleFile):

    # unpickle
    pickleData = pickle.load (bz2.open (pickleFile, "rb"))
    numSimReps = pickleData['numSimReps']
    constMLEs = pickleData['constMLEs']
    onsetMLEs = pickleData['onsetMLEs']

    # then build a dataframes for seaborn
    mles = {
        'onset' : onsetMLEs,
        'constant' : constMLEs,
    }
    frames = {}
    for (thisSimScenario, thisMLEs) in mles.items():
        thisColumns = ['onsetScenario', 'trueS2', 'rep', 's2Hat']
        if (thisSimScenario == 'onset'):
            thisColumns.append ('tHat')
        thisFrame = None
        # fill it
        for (k,v) in thisMLEs.items():
            (thisScenario, thisS2String) = k
            thisDict = {
                'onsetScenario' : np.repeat (thisScenario, numSimReps),
                'trueS2' : np.repeat (float(thisS2String), numSimReps),
                'rep' : np.arange (numSimReps),
                'baseLL' : v[:,0],
                'maxLL' : v[:,1],
                's2Hat' : v[:,2],
            }
            if (thisSimScenario == 'onset'):
                thisDict['tHat'] = v[:,3]

            # and extend the big frame
            thisDF = pd.DataFrame (thisDict)
            if (thisFrame is None):
                thisFrame = thisDF
            else:
                thisFrame = pd.concat([thisFrame, thisDF], axis=0)
        
        # store the new frame
        frames[thisSimScenario] = thisFrame

    # then put the frame with the onset results together
    preFrame = frames['onset']
    onsetFrame = preFrame[preFrame['onsetScenario'] == 'onset']

    return onsetFrame


def boxplotOnset (ax, onsetFrame):
    ax = sns.boxplot(data=onsetFrame, x="trueS2", y="tHat", fliersize=0.5, orient="v", ax=ax)

    ax.grid(visible=True, axis='y', lw=0.3)
    # add aux lines
    trueT = 2000
    ax.axline ((0,trueT), slope=0, color='C3', linestyle='dashed', lw=1)
    ax.tick_params(axis='both', labelsize="small")
    ax.set_xlabel(r'True $s_{AA}$', fontsize="small")
    ax.set_ylabel(r'Inferred $\hat{t}_{o}$', fontsize="small")

    return ax


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


def plot_ROCs(ax, LRs_varS2, H, true_s2s):
    # ax.set_title("ROC", fontsize="small")
    # ax.set_title('C)', loc='left', fontsize="medium") #, fontstyle='bold'
    ax.set_xlabel('False Positive Rate', fontsize="small")
    ax.set_ylabel('True Positive Rate', fontsize="small")
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.04, 1.04)
    ax.tick_params(labelsize=8)
    ax.axhline(0, 0, 1, ls=':', c='darkgray', lw=1)
    ax.axvline(0, 0, 1, ls=':', c='darkgray', lw=1)
    ax.axline((0, 0), (1, 1), ls=":", lw=0.5, c="0.5")

    line_colors = plt.cm.Dark2(range(len(true_s2s)))

    for idx, true_s2 in enumerate(true_s2s):
        s_pair_key = (f'{float(H) * float(true_s2):g}', true_s2)
        TP, FP, TN, FN = get_ROCs(LRs_varS2[H][("0", "0")], LRs_varS2[H][s_pair_key])  # [:,0][:,0]
        TPR = np.array(TP) / (np.array(TP) + np.array(FN))
        FPR = np.array(FP) / (np.array(FP) + np.array(TN))
        ax.plot(FPR, TPR, color=line_colors[idx], label=str(true_s2))

    return ax


def _simH_to_llH(simH):
    if simH == ".5":
        return "0.5"
    elif simH == "Inf":
        return simH
    else:
        return simH


def manual_MLEdict_to_pd(MLEdict):
    DF = pd.DataFrame(columns=['H', 'true_s1', 'true_s2', 'rep', 'ongrid_shat', 'offgrid_shat'])
    # print(MLEdict.keys())
    for H, subdict in MLEdict.items():
        print(H, subdict.keys(), type(H))
        for s_pair, MLEs in subdict.items():
            s1, s2 = s_pair
            try:
                ongrid_shat, offgrid_shat = zip(*MLEs)
            except Exception as e:
                print(e)
                print(MLEs.shape)
                sys.exit()
            df = pd.DataFrame.from_dict({'ongrid_shat': ongrid_shat, 'offgrid_shat': offgrid_shat})
            df['rep'] = df.index
            df['H'] = H
            df['true_s1'] = s1
            df['true_s2'] = s2
            DF = pd.concat([DF, df], axis=0)
    # print(DF.head())
    return DF


def plot_MLEs(ax, metaS2hats, H_2plot, true_s_2plot):
    # print(sHat_varS.keys())  # , llInit: str
    meltedDF_shat = manual_MLEdict_to_pd(metaS2hats)
    meltedDF_shat = meltedDF_shat.rename(columns={'offgrid_shat': 's2hat'})
    # print(meltedDF_shat.H.unique())
    # remove hits on the boundary (when MLE takes -8 or -9)
    print('Before removing boundary & anomaly hits, shape = ', meltedDF_shat.shape)
    meltedDF_shat = meltedDF_shat[meltedDF_shat.s2hat > -8]
    print('After removing boundary & anomaly hits, shape = ', meltedDF_shat.shape)
    # get the subset of good sampling schemes
    # DF_subset = meltedDF_shat[meltedDF_shat.H.isin([H_2plot]) &
    DF_subset = meltedDF_shat[meltedDF_shat.H.isin(H_2plot) &
                              meltedDF_shat.true_s2.astype(str).isin(true_s_2plot)].copy()
    # turn true_s2 to numeric
    DF_subset['true_s2'] = pd.to_numeric(DF_subset.loc[:, 'true_s2'])
    # print(DF_subset.iloc[:10,:]), errors='coerce'

    ax = sns.boxplot(data=DF_subset, x="true_s2", y="s2hat", hue="H",
                     fliersize=0.5, orient="v", ax=ax, palette=MY_PALETTE)
    ax.grid(visible=True, axis='y', lw=0.3)
    ax.hlines(y=list(map(float, true_s_2plot)), xmin=-0.5, xmax=len(true_s_2plot) - 0.4,
              colors=['black'] * len(true_s_2plot), linestyle='dashed', lw=0.5)
    # ax.set_ylabel(None)
    ax.set_ylabel(r'Inferred $\widehat{s}_{AA}$', fontsize="small")
    # ax.set_title(r'Inferred $\widehat{s}_{AA}$', fontsize="small")
    ax.set_ylim(-0.005, 0.015)
    ax.set_xlim(-0.5, len(true_s_2plot) - 0.4)
    ax.set_xlabel('True $s_{AA}$', fontsize="small")
    # axD.legend(loc=2, title='Dominance $h$', frameon=False) #'best'
    ax.get_legend().remove()
    return ax


def main():

    plt.rcParams['font.size'] = 10

    figname = sys.argv[1]
    Arrangement = [['A', 'B'],
                   ['C', 'D']]
    fig, ax = plt.subplot_mosaic(Arrangement, figsize=(5.5, 4),
                                 gridspec_kw={'width_ratios': [2, 4], 'height_ratios': [1, 1]})
    sns.set_theme(style="whitegrid")
    # H_2plot = ["0.5"]  # "0", , "1"
    H_2plot = ["0.5", "1.0"]  # "0", , "1"
    true_s_2plot = ("0", "0.001", "0.002", "0.003", "0.004", "0.005")  # , "0.01"

    LRs_varS2, LLs_varS2, s2Hat_varS2 = {}, {}, {}
    for H in H_2plot:
        if (H == "0.5"):
            short_H = ".5"
        elif (H == "1.0"):
            short_H = "1"
        else:
            short_H = H
        with open(
                (
                    "../supp_figS1-S7_DLsims/simulations/"
                    f"varS2_fromSeed87235_simInitSV_h{short_H}x5s_MAF.05"
                    f"_t9n40_500reps_LLfixH{H}_LLinitUnif_51xgeomGrid75e-2_offGridMLE.pkl"
                ), 
                'rb') as pk:
            (LRs_varS2[H], s2Hat_varS2[H], LLs_varS2[H]) = pickle.load(pk)
        pk.close()

    # we also need the onsetFrame
    onsetFrame = loadPickledResults (ONSET_RESULTS_PICKLE_FILE)

    # panel A: ROC h=0.5
    ax['A'] = plot_ROCs(ax['A'], LRs_varS2, H="0.5", true_s2s=true_s_2plot)
    ax['A'].set_title('a', loc='left', fontsize="medium", fontweight='bold')
    ax['A'].legend(loc="best",  # 'center left', ncol=2, bbox_to_anchor=(1, 0.5),
                   edgecolor='None', facecolor='None',
                   fontsize="xx-small", title_fontsize="x-small", frameon=False, title="True $s_{AA}$")

    # panel B: boxplots h=0.5
    # ax['B'] = plot_MLEs(ax['B'], s2Hat_varS2, "0.5", true_s_2plot)
    ax['B'] = plot_MLEs(ax['B'], s2Hat_varS2, ["0.5", "1.0"], true_s_2plot)
    ax['B'].tick_params(labelsize=8)
    ax['B'].set_title('b', loc='left', fontsize="medium", fontweight='bold')  # , fontstyle='bold'
    ax['B'].legend(["h=0.5", "h=1"], loc="best", edgecolor='None', facecolor='None',
                   fontsize="xx-small", title_fontsize="x-small", frameon=False)
    leg = ax['B'].get_legend()
    leg.legend_handles[0].set_facecolor (MY_PALETTE[0])
    leg.legend_handles[1].set_facecolor (MY_PALETTE[1])

    # panel C: ROC h=1.0
    ax['C'] = plot_ROCs(ax['C'], LRs_varS2, H="1.0", true_s2s=true_s_2plot)
    ax['C'].set_title('c', loc='left', fontsize="medium", fontweight='bold')
    ax['C'].legend(loc="best",  # 'center left', ncol=2, bbox_to_anchor=(1, 0.5),
                   edgecolor='None', facecolor='None',
                   fontsize="xx-small", title_fontsize="x-small", frameon=False, title="True $s_{AA}$")

    # # panel D: boxplots onset
    ax['D'] = boxplotOnset (ax['D'], onsetFrame)
    ax['D'].tick_params(labelsize=8)
    ax['D'].set_title('d', loc='left', fontsize="medium", fontweight='bold')  # , fontstyle='bold'

    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.975, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.6)

    # fig.tight_layout()
    # fig.savefig(figname, dpi=500, transparent=True)
    fig.savefig(figname)


if __name__ == '__main__':
    main()
