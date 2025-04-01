import numpy as np
import pandas as pd
import pickle, sys
import seaborn as sns
import matplotlib.pyplot as plt
import bz2
import pickle
from scipy.stats import chi2


def _reformat_LL_DF_to_matrix(onGrid_LLs):
    # assume the DF has both row names and column names
    loci_names = list(onGrid_LLs.ID)
    s_pairs = [colname for colname in onGrid_LLs.columns if colname != "ID"]
    # make sure they're numbers
    s_pairs = [list(map(float, pair.strip("()").split(','))) for pair in s_pairs]
    s1_list, s2_list = map(lambda x: np.array(sorted(list(set(x)))),
                           zip(*s_pairs))
    LLmatrix = np.array(onGrid_LLs.loc[:, onGrid_LLs.columns != "ID"], dtype=float)  # first column is ID
    # make sure s_pairs have tupples
    s_pairs = [tuple(pair) for pair in s_pairs]
    return np.array(LLmatrix), np.array(loci_names), s1_list, s2_list, s_pairs, len(s_pairs)


def plot_ASIP(ax, llFile, leftmost=False):
    # llFile = 'ex1_ASIP-init_51x51linGrid1e-1_LLmatrices.table'
    LLfile = pd.read_csv(llFile, sep="\t", comment="#")
    LLmatrix, loci_name, s1_list, s2_list, s_pairs, num_pairs = _reformat_LL_DF_to_matrix(LLfile)
    asip = LLmatrix[1, :].reshape((len(s2_list), len(s1_list)))
    # ### on-grid max
    max_idx = np.unravel_index(np.argmax(asip), asip.shape)
    snp_peak = asip[max_idx] / np.log(10)
    s1_0, s2_0 = s1_list[max_idx[1]], s2_list[max_idx[0]]

    # ct = ax.contour(s1_list, s2_list, asip / np.log(10), levels=40, cmap='copper')  # magma cividis
    ct = ax.contour(s1_list, s2_list, asip / np.log(10), levels=20, cmap='copper')  # magma cividis
    ax.clabel(ct, ct.levels[1::2], inline=True, fontsize='x-small')
    ax.set_xlabel('$s_{Aa}$')
    ax.get_xaxis().set_ticks([-0.075, 0, 0.075], ["-0.75", "0.00", "0.75"])
    if leftmost:
        ax.set_ylabel('$s_{AA}$')
        ax.get_yaxis().set_ticks([-0.075, 0, 0.075], ["-0.75", "0.00", "0.75"])
    else:
        # ax.get_yaxis().set_ticks([])
        ax.get_yaxis().set_visible(False)
    ax.plot(s1_0, s2_0, 'x')

    ax.annotate('(%.3f, %.3f)\n$log_{10}\mathcal{L}_{max}=$\n%.4f' % (s1_0, s2_0, snp_peak),
                (s1_0 + 0.005, s2_0 - 0.002), fontsize="xx-small", ha='left', va='top',
                xytext=(s1_0, s2_0), textcoords='offset points')
    # add reference lines
    ax.axhline(0, 0, 1, ls='--', c='darkgray', lw=1)
    ax.axvline(0, 0, 1, ls='--', c='darkgray', lw=1)
    ax.axline((0, 0), (0.01, .01), ls='--', c='darkgray', lw=1)
    ax.axline((0, 0), (0.005, .01), ls='--', c='darkgray', lw=1)
    return ax


def main():

    # plt.rcParams['font.size'] = 10
    plt.rcParams['font.size'] = 8
    thisLabelSize = 7

    figname = sys.argv[1]
    # Arrangement = [['A', 'B'],
    #                ['C', 'C'],
    #                ['D', 'D']]
    Arrangement = [['A', 'B', 'C']]
    # fig, ax = plt.subplot_mosaic(Arrangement, figsize=(6.5, 6.5),
    #                              gridspec_kw={'width_ratios': [1, 1], 'height_ratios': [2, 1.5, 1.5]})
    fig, ax = plt.subplot_mosaic(Arrangement, figsize=(6.5, 2.8),
                                 gridspec_kw={'width_ratios': [1, 1, 1]})
    sns.set_theme(style="whitegrid")

    # panel A: ASIP -- earliest
    thisTime = 17000
    llFile = f'ex1_ASIP_t{thisTime}_-init_51x51linGrid1e-1_LLmatrices.table'
    ax['A'] = plot_ASIP(ax['A'], llFile, leftmost=True)
    ax['A'].tick_params(labelsize=thisLabelSize)
    ax['A'].set_title(fr'A) $t_0=${thisTime}', loc='left', fontsize="small")
    ax['A'].set_aspect('equal')


    # panel A: ASIP -- middle
    thisTime = 15000
    # thisTime = 17000
    llFile = f'ex1_ASIP_t{thisTime}_-init_51x51linGrid1e-1_LLmatrices.table'
    ax['B'] = plot_ASIP(ax['B'], llFile, leftmost=False)
    ax['B'].tick_params(labelsize=thisLabelSize)
    ax['B'].set_title(fr'B) $t_0=${thisTime}', loc='left', fontsize="small")
    ax['B'].set_aspect('equal')

    # panel A: ASIP -- latest
    thisTime = 13105
    # thisTime = 17000
    llFile = f'ex1_ASIP_t{thisTime}_-init_51x51linGrid1e-1_LLmatrices.table'
    ax['C'] = plot_ASIP(ax['C'], llFile, leftmost=False)
    ax['C'].tick_params(labelsize=thisLabelSize)
    ax['C'].set_title(fr'C) $t_0=${thisTime}', loc='left', fontsize="small")
    ax['C'].set_aspect('equal')

    fig.suptitle ('$\it{ASIP}$ in ancient horses', fontsize="medium")

    plt.subplots_adjust(left=0.08, bottom=0.1, right=0.98, top=0.88, wspace=0.03, hspace=0.12)
    # plt.tight_layout()
    # fig.savefig(figname, dpi=500, transparent=True)
    fig.savefig(figname, dpi=600)


if __name__ == '__main__':
    main()
