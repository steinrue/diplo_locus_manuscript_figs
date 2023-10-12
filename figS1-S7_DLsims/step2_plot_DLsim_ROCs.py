"""Collect DL simulation pkls, compute and plot LR-based ROCs
Usage:
python %prog <generic_s2_pkl> <generic_s1_pkl> <sim_init> <fig_prefix>
"""
import pickle, sys, re, subprocess
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# local import
from diplo_locus.utility import _reformat_LL_DF_to_matrix


def _simH_to_llH(simH):
    if simH == ".5":
        return "0.5"
    elif simH == "Inf":
        return simH
    else:
        return simH


def _llH_to_simH(llH):
    if llH == "0.5":
        return ".5"
    else:
        return llH


def get_ROCs(neut_data, sel_data):
    # make sure things are sorted, default is ascending
    neut_data = np.sort(neut_data[~np.isnan(neut_data)])
    sel_data = np.sort(sel_data[~np.isnan(sel_data)])
    TP, FP, TN, FN = [], [], [], []
    for thr in range(len(neut_data)):
        if thr == np.nan:
            continue
        TP.append(np.sum(sel_data >= thr))
        FP.append(np.sum(neut_data >= thr))
        TN.append(np.sum(neut_data < thr))
        FN.append(np.sum(sel_data < thr))
    return TP, FP, TN, FN


def plot_ROCs(ax, LRs_varS2, H, true_s2s, which_s: str):
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

    for idx, true_s in enumerate(true_s2s):
        if which_s == 's2':
            s_pair_key = (f'{float(H) * float(true_s):g}', true_s)
        else:
            assert which_s == 's1', which_s
            try:
                s_pair_key = (true_s, f'{float(true_s) / float(H):g}')
            except Exception as e:
                print(e)
                # raise ValueError
                s_pair_key = (true_s, '0')

        neut_pool = LRs_varS2[H][_llH_to_simH(H)][("0", "0")]
        sel_pool = LRs_varS2[H][_llH_to_simH(H)][s_pair_key]
        # sanity check:
        assert (np.isnan(neut_pool).sum() < len(neut_pool)) and (np.isnan(sel_pool).sum() < len(sel_pool)), \
            f'LRs_varS2[{H}][{_llH_to_simH(H)}][{s_pair_key}] = {LRs_varS2[H][_llH_to_simH(H)][s_pair_key]}'

        TP, FP, TN, FN = get_ROCs(neut_pool, sel_pool)  # [:,0][:,0]

        # if np.any()
        # print(len(TP), len(FP), len(TN), len(FN))
        if np.any((np.array(TP) + np.array(FN)) == 0) or np.any(np.array(FP) + np.array(TN) == 0):
            passed_idx = np.where((np.array(TP) + np.array(FN) > 0) & (np.array(FP) + np.array(TN) > 0) )[0].astype(int)
            print(len(passed_idx))
            if len(passed_idx) > 0:
                nnTP, nnFP, nnTN, nnFN = TP[passed_idx], FP[passed_idx], TN[passed_idx], FN[passed_idx]
                # print(len(nnTP), len(nnFP), len(nnTN), len(nnFN))
                TPR = np.array(nnTP) / (np.array(nnTP) + np.array(nnFN))
                FPR = np.array(nnTP) / (np.array(nnTP) + np.array(nnFN))
            else:
                print(f"skip H={H}, {which_s}={true_s}, pair={s_pair_key}")
                # print(TP, FP, TN, FN)
                print(H, LRs_varS2[H].keys())
                print(_llH_to_simH(H), LRs_varS2[H][_llH_to_simH(H)].keys())
                print(LRs_varS2[H][_llH_to_simH(H)][("0", "0")])
                print(f'LRs_varS2[{H}][{_llH_to_simH(H)}][{s_pair_key}] = {LRs_varS2[H][_llH_to_simH(H)][s_pair_key]}' )
                # continue
                sys.exit()
        else:
            TPR = np.array(TP) / (np.array(TP) + np.array(FN))
            FPR = np.array(FP) / (np.array(FP) + np.array(TN))
        # TPR = np.divide(np.array(TP), np.array(TP) + np.array(FN),
        #                 where=(np.array(TP) + np.array(FN) != 0), out=np.zeros_like(np.array(TP)), dtype=int)
        # FPR = np.divide(np.array(FP), np.array(FP) + np.array(TN),
        #                 where=(np.array(FP) + np.array(TN) != 0), out=np.zeros_like(np.array(FP)), dtype=int)
        ax.plot(FPR, TPR, color=line_colors[idx], label=str(true_s))

    return ax


def _express_init(init):
    # simulation inits
    if 'Freq.' in init:
        return 'freq. 0.01'
    elif 'SV' in init:
        return 'standing var.'
    # LL comp inits
    elif init == 'Unif':
        return 'uniform initial'
    elif init == 'Freq':
        return 'initial freq. 0.01'
    else:
        raise ValueError(f'Can\'t recognize init={init}.')


def plot_rocs_wrapper(figname_prefix, LRs_varS2, H_2plot, true_s_2plot, panel_labs, sim_init, LL_init, which_s: str):
    roc_figname = f'{figname_prefix}_ROCs.png'
    # three cols of ROCs
    roc_fig, roc_ax = plt.subplots(nrows=1, ncols=len(H_2plot), sharex=True, sharey=True,
                                   figsize=(2.25 * len(H_2plot), 2.75))

    for i, h in enumerate(H_2plot):
        roc_ax[i] = plot_ROCs(roc_ax[i], LRs_varS2, H=h, true_s2s=true_s_2plot, which_s=which_s)
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
    if len(H_2plot) < 3:
        roc_fig.suptitle(f'Simulation initialized with {_express_init(sim_init)}.\n'
                         f'Likelihood computed with {_express_init(LL_init)}.',
                         fontsize='medium', ha='center', va='top', y=0.925)  #
    else:
        roc_fig.suptitle(f'Simulation initialized with {_express_init(sim_init)}. '
                         f'Likelihood computed with {_express_init(LL_init)}.',
                         fontsize='medium', ha='center', va='top', y=0.925)  #
    plt.tight_layout()  # h_pad=1.5
    roc_fig.savefig(roc_figname, dpi=350, transparent=True)
    print(roc_figname)


def main():
    s2_pkl_name = sys.argv[1]
    s1_pkl_name = sys.argv[2]
    simInit = sys.argv[3]
    # llInit = sys.argv[4]
    figname_prefix = sys.argv[4]
    # prepare to extract
    trueH_s2_regex = re.compile(r'h(\.?[0-9]+)x\ds')
    trueH_s1_regex = re.compile(r'simH(\.?[0-9]+|\w*)x\ds')
    llH_regex = re.compile(r'LLfixH(\d?.?\d+)_')
    simInit_regex = re.compile(r'simInit([a-z|A-Z]+\.?[0-9]*)_')
    llInit_regex = re.compile(r'LLinit([a-z|A-Z]+)_')
    fixS2_regex = re.compile(r'LLfixS2-(\d+)_')

    # load data
    metaLRs, metaLLs, metaS_hats = {'s2': {'Unif': {}, 'Freq': {}}, 's1': {'Unif': {}, 'Freq': {}}}, \
                                   {'s2': {'Unif': {}, 'Freq': {}}, 's1': {'Unif': {}, 'Freq': {}}}, \
                                   {'s2': {'Unif': {}, 'Freq': {}}, 's1': {'Unif': {}, 'Freq': {}}}
    keynames = ('s2', 's1')
    llInits = []
    ## load s2 first
    for i, pkl_bash_name in enumerate([s2_pkl_name, s1_pkl_name]):
        key_s = keynames[i]
        pkl_output = subprocess.check_output(f"ls {pkl_bash_name}", shell=True).decode()
        # print(pkl_output)
        pkl_list = pkl_output.split('\n')

        if len(pkl_list) == 0:
            print('pkl files doesn\'t exist.')
            sys.exit()
        # now we continue
        simInits_temp = []
        for pkl_name in pkl_list:
            if pkl_name == '':
                continue
            try:
                simInit = simInit_regex.findall(pkl_name)[0]
                llInit = llInit_regex.findall(pkl_name)[0]
            except IndexError:
                print(pkl_name, pkl_list)
                sys.exit()
            simInits_temp.append(simInit)
            llInits.append(llInit)
            assert len(set(simInits_temp)) == 1, simInits_temp
            # assert len(set(llInits)) == 1, llInits

            if key_s == 's1':
                simH = trueH_s1_regex.findall(pkl_name)[0]
                if 'LLfixH' in pkl_name:
                    llH = llH_regex.findall(pkl_name)[0]
                    # print(f'Loading {key_s} simulations from simInit = {simInit},'
                    #     f'simH = {simH}, ll_init = {llInit}, fixed ll H = {llH}.')
                elif 'LLfixS2' in pkl_name:
                    fixS2 = fixS2_regex.findall(pkl_name)[0]
                    # print(f'Loading {key_s} simulations from simInit = {simInit},'
                    #       f'simH = {simH}, ll_init = {llInit}, fixed ll s2 = {fixS2}.')
                    if int(fixS2) == 0:
                        llH = "Inf"
                    else:
                        # assert fixS2.isdecimal(), f'Cannot parse {fixS2} as numeric.'
                        print('s2 fixed; H value vary by s1')
                        llH = "fix_s1"
                else:
                    raise ValueError(f'Cannot recognize filename: {pkl_name}')
            else:
                simH = trueH_s2_regex.findall(pkl_name)[0]
                llH = llH_regex.findall(pkl_name)[0]
                if float(llH) >= 1 or float(llH) == 0:
                    llH = str(int(float(llH)))
            if llH != _simH_to_llH(simH):
                continue
            print(f'Loading {key_s} simulations from simInit = {simInit},'
                    f'simH = {simH}, ll_init = {llInit}, llH = {llH}.')

            # initialize
            if llH not in metaLRs[key_s][llInit]:
                metaLRs[key_s][llInit][llH], metaS_hats[key_s][llInit][llH], metaLLs[key_s][llInit][llH] = {}, {}, {}
            # each will be a dictionary with llH-decided s_pairs as key
            with open(pkl_name, 'rb') as pk:
                temp_LRs, temp_Shats, temp_LLs = pickle.load(pk)
                # metaLRs[key_s][llInit][llH], metaLLs[key_s][llInit][llH], metaS_hats[key_s][llInit][llH] = pickle.load(pk)
            pk.close()
            metaLRs[key_s][llInit][llH][simH], metaS_hats[key_s][llInit][llH][simH], metaLLs[key_s][llInit][llH][
                simH] = temp_LRs, temp_Shats, temp_LLs
            # print(llH, simH)
            # print(metaLRs[key_s][llInit].keys(), llH, metaLRs[key_s][llInit][llH].keys())  #, metaLRs[key_s][llInit][llH][simH].keys()
            assert ('0', '0') in metaLRs[key_s][llInit][llH][simH], metaLRs[key_s][llInit][llH][simH].keys()
            # sanity check
            for s_pair, scores in temp_LRs.items():
                assert (np.isnan(scores).sum() < len(scores)), \
                    f'temp_LRs[{llH}][{simH}][{s_pair}] = {scores}.\n{pkl_name}'

    # done loading file. Now plot
    H_2plot = (["0.5", "0", "1"], ["5", "Inf"])  # "0.5", "0", "1",".5",
    # H_2plot = ([0.5, 0., 1.], [5., "Inf"])  # "0.5", "0", "1",
    Dom = {".5": "additive", "0.5": "additive", '0': "complete recessive", "1": "complete dominance", "5": "over-dominance",
           "Inf": r"$s_{AA}=0$"}
    true_s_2plot = (["0", "0.001", "0.002", "0.003", "0.004", "0.005"],  #
                    ["0", "0.001", "0.002", "0.003", "0.004", "0.005"])  # , "0.01", "0.02"

    # just give it a shot?
    # plot ROCs (fig S1--3)
    llInits = list(set(llInits))
    for i, llInit in enumerate(llInits):
        plot_rocs_wrapper(f'{figname_prefix}_ll{llInit}_0-1H', metaLRs['s2'][llInit], H_2plot[0], true_s_2plot[0],
                          panel_labs=['A', 'B', 'C'], sim_init=simInit, LL_init=llInit, which_s='s2')  # , 'D'
        # clean the slate
        plt.clf()

        plot_rocs_wrapper(f'{figname_prefix}_ll{llInit}_5-InfH', metaLRs['s1'][llInit], H_2plot[1], true_s_2plot[1],
                          panel_labs=['A', 'B'], sim_init=simInit, LL_init=llInit, which_s='s1')  # , 'D', 'C'
        # clean the slate
        plt.clf()


if __name__ == '__main__':
    main()
