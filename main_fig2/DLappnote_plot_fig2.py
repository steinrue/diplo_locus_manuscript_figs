import numpy as np
import pandas as pd
import pickle, sys
import seaborn as sns
import matplotlib.pyplot as plt
import bz2
import pickle


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


def plot_ASIP(ax):
    # llFile = 'ex1_ASIP-init_51x51linGrid1e-1_LLmatrices.table'
    # llFile = '../supp_ASIP/ex1_ASIP_t17000_-init_51x51linGrid1e-1_LLmatrices.table'
    llFile = '../supp_figS16_ASIP/ex1_ASIP_t13105_-init_51x51linGrid1e-1_LLmatrices.table'
    LLfile = pd.read_csv(llFile, sep="\t", comment="#")
    LLmatrix, loci_name, s1_list, s2_list, s_pairs, num_pairs = _reformat_LL_DF_to_matrix(LLfile)
    asip = LLmatrix[1, :].reshape((len(s2_list), len(s1_list)))
    # ### on-grid max
    max_idx = np.unravel_index(np.argmax(asip), asip.shape)
    snp_peak = asip[max_idx] / np.log(10)
    s1_0, s2_0 = s1_list[max_idx[1]], s2_list[max_idx[0]]
    ### find neut
    # neut_i, neut_j = list(s2_list).index(0), list(s1_list).index(0)
    # snp_neut = asip[neut_i, neut_j] / np.log(10)
    # read maxLR
    # mlrFile = '../../examples/example_output/Ex1_ASIP-init_50x50linGrid5e-2_off-grid_maxLLs.txt'
    # maxLR = pd.read_csv(mlrFile, sep="\t", comment="%")
    # maxLR.rename(columns={'# locus': 'locus'})

    # ct = ax.contour(s1_list, s2_list, asip / np.log(10), levels=40, cmap='copper')  # magma cividis
    ct = ax.contour(s1_list, s2_list, asip / np.log(10), levels=35, cmap='copper')  # magma cividis
    ax.clabel(ct, ct.levels[1::2], inline=True, fontsize='x-small')
    ax.set_xlabel('$s_{Aa}$')
    ax.set_ylabel('$s_{AA}$')
    ax.set_title('$\it{ASIP}$ in ancient horses', fontsize="small")  # $log_{10}$likelihood for
    # ax.set_title('A)', loc='left', fontsize="medium") #, fontstyle='bold'
    ax.plot(s1_0, s2_0, 'x')

    # ax.annotate('$log_{10}\mathcal{L}_{max}=$\n%.4f\n(%.3f, %.3f)' %
    ax.annotate('(%.3f, %.3f)\n$log_{10}\mathcal{L}_{max}=$\n%.4f' % (s1_0, s2_0, snp_peak),
                (s1_0 + 0.005, s2_0 - 0.002), fontsize="x-small", ha='left', va='top',
                xytext=(s1_0, s2_0), textcoords='offset points')
    # add reference lines
    ax.axhline(0, 0, 1, ls='--', c='darkgray', lw=1)
    ax.axvline(0, 0, 1, ls='--', c='darkgray', lw=1)
    ax.axline((0, 0), (0.01, .01), ls='--', c='darkgray', lw=1)
    ax.axline((0, 0), (0.005, .01), ls='--', c='darkgray', lw=1)
    return ax


def plot_LCT(ax):
    ## load data
    ukfile = 'UK_v54.1_1240K_from4500_rs4988235_MAF.05_51x51linGrid2e-1_LLmatrices.table.gz'
    ukLLs = pd.read_csv(ukfile, sep="\t", comment="#")  # , header=[0], index_col=0
    ukLLs["ID"] = ukLLs.index
    ukLLmatrix, snp_names, s1_list, s2_list, s_pairs, num_pairs = _reformat_LL_DF_to_matrix(ukLLs)
    ## get things in shape
    surface = np.array(ukLLmatrix[0, :].reshape((len(s2_list)), len(s1_list)))
    # ### on-grid max
    max_idx = np.unravel_index(np.argmax(surface), surface.shape)
    snp_peak = surface[max_idx] / np.log(10)
    s1_0, s2_0 = s1_list[max_idx[1]], s2_list[max_idx[0]]
    # print(max_idx, s1_0, s2_0, snp_peak)
    ### plot contour
    ct2 = ax.contour(s1_list, s2_list, surface / np.log(10), levels=40, cmap='copper')  # magma cividis
    ax.clabel(ct2, ct2.levels, inline=True, fontsize='x-small')  # [::1]
    ax.set_xlabel('$s_{Aa}$')
    ax.set_xlim(-0.2, 0.2)
    ax.set_ylim(-0.2, 0.2)
    ax.set_ylabel('$s_{AA}$')
    ax.set_title(r'$\it{LCT}$ rs4988235 in GB', fontsize="small")
    # ax.set_title('B)', loc='left', fontsize="medium") #, fontstyle='oblique'

    ax.plot(s1_0, s2_0, 'x')
    ax.annotate('(%.3f, %.3f)\n$log_{10}\mathcal{L}_{max}=$\n%.4f' % (s1_0, s2_0, snp_peak),
                (s1_0, s2_0), xytext=(s1_0, s2_0), ha='left', va='top', fontsize='x-small',
                textcoords='offset points', annotation_clip=True)
    # add reference lines
    ax.axhline(0, 0, 1, ls='--', c='darkgray', lw=1)
    ax.axvline(0, 0, 1, ls='--', c='darkgray', lw=1)
    ax.axline((0, 0), (0.01, .01), ls='--', c='darkgray', lw=1)
    ax.axline((0, 0), (0.005, .01), ls='--', c='darkgray', lw=1)
    return ax


def _split_locus_name(locus_name):
    split_names = locus_name.split("_")
    if len(split_names) > 3:
        Chr, physPos = split_names[0], split_names[1]
        rsID = '_'.join(split_names[2:])
    elif len(split_names) == 3:
        Chr, physPos, rsID = split_names
    else:
        return False
    return int(Chr), int(physPos), rsID


from scipy.stats import chi2


def read_single_file(filename):
    # chrDF = pd.read_csv(filename, sep="\t", comment="#", header=None)  #
    chrDF = pd.read_csv(filename, sep="\t", comment="#")  #
    # print(chrDF.head())
    if len(chrDF.columns) == 7:
        chrDF.columns = ['ID', 'ongrid_s2hat', 'ongrid_maxLogLikelihood',
                         's2hat', 'maxLogLikelihood', 'MLR', 'chi2_p']
    else:  # just read from file
        with open(filename, 'r') as scores:
            for l in scores:
                if l.startswith("ID"):
                    header = l.strip().split("\t")
                    assert len(header) == len(chrDF.columns), f'len(header) = {len(header)};' \
                                                              f'len(chrDF.columns) = {len(chrDF.columns)}.'
                    chrDF.columns = header
        print(chrDF.head())
        # sys.exit()
    # split the locus name <-- currently included in the file already!
    if 'Chr' not in chrDF.columns or 'physPos' not in chrDF.columns or 'rsID' not in chrDF.columns:
        chrDF['Chr'], chrDF['physPos'], chrDF['rsID'] = zip(*chrDF['ID'].apply(_split_locus_name))

    # if only 3 columns, only ongrid max are available
    if chrDF.shape[1] <= 5:
        # then there's no interpolated value
        assert 'maxLogLikelihood' not in chrDF.columns, f'Chr{c} maxLL file header: {chrDF.columns}'
        # rename
        chrDF = chrDF.rename(columns={'ongrid_s2hat': 's2_hat'})
        # reorder
        chrDF = chrDF[["Chr", "physPos", "rsID", "MLR", "s2_hat"]]
        # get p-value
        chrDF['chi2_pval'] = chrDF['MLR'].apply(chi2.sf, df=1)
    elif 'chi2_p' in chrDF.columns:
        assert 'MLR' in chrDF.columns, f'chrDF.columns = {chrDF.columns}'
        assert 's2hat' in chrDF.columns, f'chrDF.columns = {chrDF.columns}'
        # filter -8s & -9s
        chrDF = chrDF[(chrDF.s2hat != -8) & (chrDF.s2hat != -9)]
        # rename
        chrDF = chrDF.rename(columns={'s2hat': 's2_hat', 'chi2_p': 'chi2_pval'})
        # reorder
        chrDF = chrDF[["Chr", "physPos", "rsID", "MLR", "s2_hat", "chi2_pval"]]
    else:
        # there must exist off-grid MLEs
        assert 's2hat' in chrDF.columns, f'chrDF.columns = {chrDF.columns}'
        # filter -8s & -9s
        chrDF = chrDF[(chrDF.s2hat != -8) & (chrDF.s2hat != -9)]
        # rename
        chrDF = chrDF.rename(columns={'s2hat': 's2_hat'})
        # reorder
        chrDF = chrDF[["Chr", "physPos", "rsID", "MLR", "s2_hat"]]
        # get p-value
        chrDF['chi2_pval'] = chrDF['MLR'].apply(chi2.sf, df=1)
    # done
    return chrDF


# TODO annotate top SNP
def plot_chr_manhattan(ax, filename, analysisType, chrom, cutoff):
    # read file
    if (analysisType == "lct"):
        dt = read_single_file(filename)
        # now plot
        ## make sure we only have one chrom
        assert dt['Chr'].value_counts().shape[0] == 1, \
            f'Only accept data from the same chrom. Now chrom in {dt.Chr.value_counts()}'
        theseTicks = [0, 50, 100, 150, 200, 250]
    elif (analysisType == "dsim"):
        # we have it pickled
        pickleData = pickle.load (bz2.open (filename, "rb"))
        commonPositions = pickleData["commonPositions"]
        fisherPValues = pickleData["fisherPValues"]
        # need some different names, so the script works
        thisChr = f"chr{chrom}"
        dt = pd.DataFrame ({
            'physPos' : commonPositions[thisChr],
            'chi2_pval' : fisherPValues[thisChr],
        })
        theseTicks = [0, 5, 10, 15, 20]
    else:
        assert (False), "Only valif analysis-types are ['lct', 'dmel']."
    
    dt['log10_p'] = -np.log10(np.array(dt.chi2_pval.copy(), dtype=float))
    # ax = sns.scatterplot(data=dt, x="physPos", y="log10_p", ax=ax)
    ax.plot(dt.physPos / 1e6, dt.log10_p, ".", markersize=2, rasterized=True,label='_nolegend_')
    # plot cutoff
    ax.axhline(y=-np.log10(cutoff / dt.shape[0]), ls='dashed', lw=0.5, c='black')
    if (analysisType == "lct"):
        # annotate top SNP
        top = dt.iloc[np.argmax(dt.log10_p), :]
        print(top)  # top['rsID']
        ax.annotate('rs4988235', xy=(top['physPos']/1e6, top['log10_p']),
                    xytext=(3 + top['physPos']/1e6, top['log10_p']),
                    fontsize=8, ha='left', va='top', family='monospace',
                    annotation_clip=True, bbox=dict(boxstyle='round', alpha=0.4,
                                                    pad=0.2, lw=0.5, fc='w', edgecolor='black'))
    ax.set_xticks(ticks=theseTicks)
    ax.set_xlabel(f'Positions on chromosome {chrom} (Mbp)', fontsize="small")
    ax.set_ylabel(r'$-\log_{10}p_{\chi^2(1)}$', fontsize="small")
    ax.legend(["BF threshold"], fontsize="small")
    return ax


def main():

    plt.rcParams['font.size'] = 10

    figname = sys.argv[1]
    Arrangement = [['A', 'B'],
                   ['C', 'C'],
                   ['D', 'D']]
    fig, ax = plt.subplot_mosaic(Arrangement, figsize=(6.5, 6.5),
                                 gridspec_kw={'width_ratios': [1, 1], 'height_ratios': [2, 1.5, 1.5]})
    sns.set_theme(style="whitegrid")

    # panel A: ASIP
    ax['A'] = plot_ASIP(ax['A'])
    ax['A'].tick_params(labelsize=8)
    ax['A'].set_title('a', loc='left', fontsize="medium", fontweight='bold')

    # panel B: LCT
    ax['B'] = plot_LCT(ax['B'])
    ax['B'].tick_params(labelsize=8)
    ax['B'].set_title('b', loc='left', fontsize="medium", fontweight='bold')

    # panel C: chr2
    filename = '../supp_figS17-S18_LCT/v54_chr2_UK_1240K_from4500_MAF.05_fixH.5_51x5e-1geomGrid_off-grid_maxLLs.txt'
    ax['C'] = plot_chr_manhattan(ax['C'], filename, "lct", chrom=2, cutoff=0.05)
    ## de-spine
    sns.despine(ax=ax['C'])
    ax['C'].tick_params(labelsize=8)
    ax['C'].set_title('c', loc='left', fontsize="medium", fontweight='bold')  # , fontstyle='bold'

    # panel D: dsim
    filename = '../supp_figS19-S21_dsim/results/analyzed_F0-F60_fisher.pickle.bz2'
    ax['D'] = plot_chr_manhattan(ax['D'], filename, "dsim", chrom="2L", cutoff=0.05)
    ## de-spine
    sns.despine(ax=ax['D'])
    ax['D'].tick_params(labelsize=8)
    ax['D'].set_title('d', loc='left', fontsize="medium", fontweight='bold')  # , fontstyle='bold'

    plt.tight_layout()
    # fig.savefig(figname, dpi=500, transparent=True)
    fig.savefig(figname, dpi=600)


if __name__ == '__main__':
    main()
