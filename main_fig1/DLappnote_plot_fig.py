import numpy as np
import pandas as pd
import pickle, sys
import seaborn as sns
import matplotlib.pyplot as plt


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
    # llFile = '../../examples/example_output/Ex1_ASIP-init_50x50linGrid5e-2_LLmatrices.table'
    llFile = 'ASIP-init_200x200xlinGrid5e-2_LLmatrices.table'
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

    ct = ax.contour(s1_list, s2_list, asip / np.log(10), levels=40, cmap='copper')  # magma cividis
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
    # ukfile = 'England_v54_selectLoci_1240K_4500-0BP_noSG_noRelatives_noContam_gen27_Ne5000_101x101linGrid5e-1_Unif_LLsurfaces.table.gz'
    # ukfile = 'rs4988235_v54_England_since4500_2D_50x5e-1linGrid_LLmatrices.table'
    ukfile = 'UK_noSG_1240K_noDG_noFam_strictPASS_from4500_to0_rs4988235_minK2_missing0.05_MAF.05_200x200linGrid.2_LLmatrices.table.gz'
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
    ax.set_title(r'$\it{LCT}/\it{MCM6}$ rs4988235 in UK', fontsize="small")
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
    DF_subset = meltedDF_shat[meltedDF_shat.H.isin(H_2plot) &
                              meltedDF_shat.true_s2.astype(str).isin(true_s_2plot)].copy()
    # turn true_s2 to numeric
    DF_subset['true_s2'] = pd.to_numeric(DF_subset.loc[:, 'true_s2'])
    # print(DF_subset.iloc[:10,:]), errors='coerce'

    ax = sns.boxplot(data=DF_subset, x="true_s2", y="s2hat", hue="H",
                     fliersize=0.5, orient="v", ax=ax)
    ax.grid(visible=True, axis='y', lw=0.3)
    ax.hlines(y=list(map(float, true_s_2plot)), xmin=-0.5, xmax=len(true_s_2plot) - 0.4,
              colors=['black'] * len(true_s_2plot), linestyle='dashed', lw=0.5)
    # ax.set_ylabel('Inferred $\widehat{s}_{AA}$')
    ax.set_ylabel(None)
    ax.set_title(r'Inferred $\widehat{s}_{AA}$', fontsize="small")
    # ax.set_ylim(-0.0025, 0.01)
    ax.set_xlim(-0.5, len(true_s_2plot) - 0.4)
    ax.set_xlabel('True $s_{AA}$', fontsize="small")
    # axD.legend(loc=2, title='Dominance $h$', frameon=False) #'best'
    ax.get_legend().remove()
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
    chrDF = pd.read_csv(filename, sep="\t", comment="#", header=None)  #
    # print(chrDF.head())
    if len(chrDF.columns) == 7:
        chrDF.columns = ['# locus', 'ongrid_s2hat', 'ongrid_maxLogLikelihood',
                         's2hat', 'maxLogLikelihood', 'MLR', 'chi2_p']
    else:  # just read from file
        with open(filename, 'r') as scores:
            for l in scores:
                if l.startswith("# locus"):
                    header = l.strip().split("\t")
                    assert len(header) == len(chrDF.columns), f'len(header) = {len(header)};' \
                                                              f'len(chrDF.columns) = {len(chrDF.columns)}.'
                    chrDF.columns = header
        print(chrDF.head())
        # sys.exit()
    # split the locus name <-- currently included in the file already!
    if 'Chr' not in chrDF.columns or 'physPos' not in chrDF.columns or 'rsID' not in chrDF.columns:
        chrDF['Chr'], chrDF['physPos'], chrDF['rsID'] = zip(*chrDF['# locus'].apply(_split_locus_name))

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
def plot_chr_manhattan(ax, filename, chrom, cutoff):
    # read file
    dt = read_single_file(filename)
    # now plot
    ## make sure we only have one chrom
    assert dt['Chr'].value_counts().shape[0] == 1, \
        f'Only accept data from the same chrom. Now chrom in {dt.Chr.value_counts()}'
    dt['log10_p'] = -np.log10(np.array(dt.chi2_pval.copy(), dtype=float))
    # ax = sns.scatterplot(data=dt, x="physPos", y="log10_p", ax=ax)
    ax.plot(dt.physPos / 1e6, dt.log10_p, ".", markersize=2)
    # plot cutoff
    ax.axhline(y=-np.log10(cutoff / dt.shape[0]), ls='dashed', lw=0.5, c='black')
    # annotate top SNP
    top = dt.iloc[np.argmax(dt.log10_p), :]
    print(top)  # top['rsID']
    ax.annotate('rs4988235', xy=(top['physPos']/1e6, top['log10_p']),
                xytext=(3 + top['physPos']/1e6, top['log10_p']),
                fontsize=8, ha='left', va='top', family='monospace',
                annotation_clip=True, bbox=dict(boxstyle='round', alpha=0.4,
                                                pad=0.2, lw=0.5, fc='w', edgecolor='black'))
    ax.set_xticks(ticks=[0, 50, 100, 150, 200, 250])
    ax.set_xlabel(f'Positions on chromosome {chrom} (Mbp)', fontsize="small")
    ax.set_ylabel(r'$-\log_{10}p_{\chi^2(1)}$', fontsize="small")
    return ax


def main():
    figname = sys.argv[1]
    Arrangement = [['A', 'B', 'B', 'B'],
                   ['C', 'C', 'C', 'C'],
                   ['D', 'D', 'E', 'E']]
    fig, ax = plt.subplot_mosaic(Arrangement, figsize=(7.5, 7.5),
                                 gridspec_kw={'width_ratios': [9, 1, 5, 5], 'height_ratios': [1.2, 1, 2]})
    sns.set_theme(style="whitegrid")
    H_2plot = ["0.5"]  # "0", , "1"
    true_s_2plot = ("0", "0.001", "0.002", "0.003", "0.004", "0.005")  # , "0.01"

    LRs_varS2, LLs_varS2, s2Hat_varS2 = {}, {}, {}
    for H in H_2plot:
        if H == "0.5":
            short_H = ".5"
        else:
            short_H = H
        with open(
                f'LLs_fromSeed87235_simInitSV_h{short_H}x7s_MAF.05_t9n40_500reps_LLfixH{H}_LLinitUnif_51xgeomGrid75e-2_offGridMLE.pkl',
                'rb') as pk:
            (LRs_varS2[H], s2Hat_varS2[H], LLs_varS2[H]) = pickle.load(pk)
        pk.close()

    # panel A: ROC
    ax['A'] = plot_ROCs(ax['A'], LRs_varS2, H="0.5", true_s2s=true_s_2plot)
    ax['A'].set_title('A)', loc='left', fontsize="medium", fontweight='bold')
    ax['A'].legend(loc="best",  # 'center left', ncol=2, bbox_to_anchor=(1, 0.5),
                   edgecolor='None', facecolor='None',
                   fontsize="xx-small", title_fontsize="x-small", frameon=False, title="True $s_{AA}$")

    # panel B: boxplots
    ax['B'] = plot_MLEs(ax['B'], s2Hat_varS2, H_2plot, true_s_2plot)
    ax['B'].tick_params(labelsize=8)
    ax['B'].set_title('B)', loc='left', fontsize="medium", fontweight='bold')  # , fontstyle='bold'

    # panel C: chr2
    filename = 'v54_UK-chr2_1240K_from4500_to0_minK2_minObs0.01_MAF.05_fixH.5_51x5e-1geomGrid_off-grid_maxLLs.txt'
    ax['C'] = plot_chr_manhattan(ax['C'], filename, chrom=2, cutoff=0.05)
    ## de-spine
    sns.despine(ax=ax['C'])
    ax['C'].tick_params(labelsize=8)
    ax['C'].set_title('C)', loc='left', fontsize="medium", fontweight='bold')  # , fontstyle='bold'

    # panel D: LCT
    ax['D'] = plot_LCT(ax['D'])
    ax['D'].tick_params(labelsize=8)
    ax['D'].set_title('D)', loc='left', fontsize="medium", fontweight='bold')

    # panel E: ASIP
    ax['E'] = plot_ASIP(ax['E'])
    ax['E'].tick_params(labelsize=8)
    ax['E'].set_title('E)', loc='left', fontsize="medium", fontweight='bold')

    # with open(f"t10n40_every1000gen_strongerFilter_4h8s2x500reps_MLRs-s2Hats_50x50geomGrid_initUnif.pkl", "rb") as pk:
    # metaLRs, metaLLs, metaS2hats = {}, {}, {}

    plt.tight_layout()
    fig.savefig(figname, dpi=500, transparent=True)


if __name__ == '__main__':
    main()
