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


# prepare to extract
trueH_s2_regex = re.compile(r'h(\.?[0-9]+)x\ds')
llH_regex = re.compile(r'LLfixH(\d?.?\d+)_')
simInit_regex = re.compile(r'simInit([a-z|A-Z]+\.?[0-9]*)_')
llInit_regex = re.compile(r'LLinit([a-z|A-Z]+)_')


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


def load_pkl_LRs(freq_pkl_name, sv_pkl_name):
    # load data. key order: simInit, llinit, H
    metaLRs = {"sim_freq": {'Unif': {}, 'Freq': {}}, "sim_sv": {'Unif': {}, 'Freq': {}}}
    # metaLLs = {"sim_freq": {'Unif': {}, 'Freq': {}}, "sim_sv": {'Unif': {}, 'Freq': {}}}

    keynames = ("sim_freq", "sim_sv")
    llInits = []
    ## load data
    for i, pkl_bash_name in enumerate([freq_pkl_name, sv_pkl_name]):
        key_simInit = keynames[i]
        pkl_output = subprocess.check_output(f"ls {pkl_bash_name}", shell=True).decode()
        # print(pkl_output)
        pkl_list = pkl_output.split('\n')

        if len(pkl_list) == 0:
            print('pkl files doesn\'t exist.')
            sys.exit()
        # now we continue
        for pkl_name in pkl_list:
            if pkl_name == '':
                continue
            try:
                llInit = llInit_regex.findall(pkl_name)[0]
            except IndexError:
                print(pkl_name, pkl_list)
                sys.exit()
            llInits.append(llInit)

            simH = trueH_s2_regex.findall(pkl_name)[0]
            llH = llH_regex.findall(pkl_name)[0]
            if float(llH) >= 1 or float(llH) == 0:
                llH = str(int(float(llH)))
            # only consider sim_H == LL_H, and only read ('0', '0')
            if llH != _simH_to_llH(simH):
                continue
            print(f'Loading {key_simInit} simulations from simH = {simH}, ll_init = {llInit}, llH = {llH}.')

            # initialize
            if llH not in metaLRs[key_simInit][llInit]:
                metaLRs[key_simInit][llInit][llH] = {}
                # metaLLs[key_simInit][llInit][llH] = {}
            # each will be a dictionary with llH-decided s_pairs as key
            with open(pkl_name, 'rb') as pk:
                temp_LRs, temp_Shats, temp_LLs = pickle.load(pk)
            pk.close()
            assert ('0', '0') in temp_LRs, temp_LRs.keys()

            metaLRs[key_simInit][llInit][llH][simH] = temp_LRs[('0', '0')]
            # metaLLs[key_simInit][llInit][llH][simH] = temp_LLs[('0', '0')]
    return metaLRs


from scipy.stats import chi2


def main():
    # let's only look at H == 0.5 here
    freq_pkl_name = sys.argv[1]
    sv_pkl_name = sys.argv[2]
    figname = sys.argv[3]

    # load file
    metaLRs = load_pkl_LRs(freq_pkl_name, sv_pkl_name)

    # done loading file. Now plot
    init_colors = plt.cm.Dark2(range(3))
    init_shapes = ['o', '+', 'x']
    init_sizes = [1, 2, 2]
    init_labels = ['sim: Freq. 0.01,      LR: Freq. 0.01',
                   'sim: Freq. 0.01,      LR: Uniform',
                   'sim: Standing var., LR: Uniform']

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))
    # set up
    ax.set_ylabel(r'$-log_{10}p_{\chi^2(1)}$', fontsize="medium")
    ax.set_xlabel(r'$-log_{10}$ quantile', fontsize="small")
    ax.tick_params(axis='both', labelsize="x-small")
    # ref x=y line
    ax.axline((0, 0), (1, 1), ls="--", lw=0.5, color="grey")
    # loop through inits
    H = "0.5"
    for i, LRs in enumerate([metaLRs["sim_freq"]['Freq'], metaLRs["sim_freq"]["Unif"], metaLRs["sim_sv"]['Unif']]):
        # now we plotLRs[H]
        mlrs = LRs[H][_llH_to_simH(H)]
        n_reps = len(mlrs)
        pval_fromLRs = chi2.sf(mlrs, df=1)
        pval_fromLRs.sort()
        # print(pval_fromLRs[:5], pval_fromLRs[n_reps//2 - 2: n_reps//2 + 2], pval_fromLRs[-3:])
        y_coords = -np.log10(pval_fromLRs)
        x_coords = -np.log10(np.arange(1, n_reps + 1)/n_reps)
        ax.plot(x_coords, y_coords, init_shapes[i], ms=init_sizes[i], c=init_colors[i],
                ls="", label=init_labels[i])
    # make it square
    ax.set_xlim( ax.get_ylim() )

    ax.legend(title='Initial Conditions', markerscale=1.5,
              loc="best", edgecolor='None', facecolor='None',
              fontsize="xx-small", title_fontsize="x-small")

    plt.tight_layout()
    fig.savefig(figname, dpi=500, transparent=True)


if __name__ == '__main__':
    main()
