"""Only plot boxplots; Two panels: one with s2_grid + H in (0,1), another with h=5, s1_grid, and LL_fix_s2=0
Add annotations for initial conditions as fig title for clarity"""
import numpy as np
import pandas as pd
import pickle, sys, re, subprocess
import seaborn as sns
import matplotlib.pyplot as plt


# from diplo_locus.utility import _reformat_LL_DF_to_matrix

def _simH_to_llH(simH):
    if simH == ".5":
        return "0.5"
    elif simH == "Inf":
        return simH
    else:
        return simH


def _llH_to_simH(llH):
    # return '%.1f' % (float(llH))
    if llH == "0.5":
        return ".5"
    else:
        return llH


def manual_MLEdict_to_pd(MLEdict, same_H=False):
    DF = pd.DataFrame(columns=['llH', 'simH', 'true_s1', 'true_s2', 'rep', 'ongrid_shat', 'offgrid_shat'])
    # print(MLEdict.keys())
    for LL_H, subdict in MLEdict.items():
        # print(LL_H, subdict.keys(), type(LL_H))
        for simH, subsubdict in subdict.items():
            if same_H:
                if _simH_to_llH(simH) != LL_H:
                    continue

            for s_pair, MLEs in subsubdict.items():
                s1, s2 = s_pair
                ongrid_shat, offgrid_shat = zip(*MLEs)
                df = pd.DataFrame.from_dict({'ongrid_shat': ongrid_shat, 'offgrid_shat': offgrid_shat})
                df['rep'] = df.index
                # if LL_H.isdecimal():
                #     df['llH'] = float(LL_H)
                # else:
                df['llH'] = LL_H
                if simH.isdecimal():
                    df['simH'] = float(simH)
                else:
                    df['simH'] = simH
                df['true_s1'] = s1
                df['true_s2'] = s2
                DF = pd.concat([DF, df], axis=0)
    # print(DF.head())
    return DF


def plot_MLE_boxplot_varS(ax, sHat_varS, H_2plot, true_s_2plot, which_s: str,
                          absolute_MLEs: bool, colors=None, collect_outliers=False,
                          outlier_bounds=None, note_outliers=False):
    # print(sHat_varS.keys())  # , llInit: str
    meltedDF_shat = manual_MLEdict_to_pd(sHat_varS, same_H=True)
    # for plotting purpose, only keep llH == simH
    # meltedDF_shat = meltedDF_shat[meltedDF_shat.llH == meltedDF_shat.simH]
    meltedDF_shat['H'] = meltedDF_shat['llH']
    meltedDF_shat = meltedDF_shat.drop(columns=['llH', 'simH'])

    if which_s == "s2":
        meltedDF_shat = meltedDF_shat.rename(columns={'offgrid_shat': 's2hat'})
        meltedDF_shat['s1hat'] = 0.
        # print(meltedDF_shat.H.unique())
        for H in meltedDF_shat.H.unique():
            temp = meltedDF_shat[meltedDF_shat.H == H].copy()
            if str(H).isdecimal():
                temp.loc[:, 's1hat'] = temp.s2hat * float(H)
            elif float(H) == 0.5:
                print(f'meltedDF_shat[meltedDF_shat.H.astype(float) == float(H)].shape = {temp.shape}')
                temp.loc[:, 's1hat'] = temp.s2hat * float(H)
            else:
                raise ValueError('cannot process this H value for now:', H)
            meltedDF_shat[meltedDF_shat.H == H] = temp
    else:
        assert which_s == 's1', which_s
        meltedDF_shat = meltedDF_shat.rename(columns={'offgrid_shat': 's1hat'})
        meltedDF_shat['s2hat'] = 0.  # meltedDF_shat['s1hat'], , 0.
        # print(meltedDF_shat.H.unique())
        for H in meltedDF_shat.H.unique():
            temp = meltedDF_shat[meltedDF_shat.H == H].copy()
            if float(H) == 0.5:
                print(f'meltedDF_shat[meltedDF_shat.H.astype(float) == float(H)].shape = {temp.shape}')
            if str(H).isdecimal():
                assert float(H) > 0, H
                temp.loc[:, 's2hat'] = temp.s1hat / float(H)
            elif str(H).lower() == "inf":
                temp.loc[:, 's2hat'] = 0
            else:
                raise ValueError('cannot process this H value for now:', H)
            meltedDF_shat[meltedDF_shat.H == H] = temp
    # remove hits on the boundary (when MLE takes -8 or -9)
    print('Before removing boundary & anomaly hits, shape = ', meltedDF_shat.shape)
    meltedDF_shat = meltedDF_shat[meltedDF_shat[f"{which_s}hat"] > -8].copy()
    print('After removing boundary & anomaly hits, shape = ', meltedDF_shat.shape)
    # get the subset of good sampling schemes
    DF_subset = meltedDF_shat[meltedDF_shat.H.isin(H_2plot) &
                              meltedDF_shat[f"true_{which_s}"].astype(str).isin(true_s_2plot)].copy()

    # take abs if needed
    if absolute_MLEs:
        print('Showing absolute values of MLEs')
        DF_subset[f"{which_s}hat"] = DF_subset[f"{which_s}hat"].abs()
    else:
        print('Showing signed values of MLEs')

    if colors is None:
        ax = sns.boxplot(data=DF_subset, x=f"true_{which_s}", y=f"{which_s}hat", hue="H",
                         palette="tab10", order=true_s_2plot,
                         fliersize=0.5, orient="v", linewidth=0.5, ax=ax)
    else:
        ax = sns.boxplot(data=DF_subset, x=f"true_{which_s}", y=f"{which_s}hat", hue="H",
                         palette=colors, order=true_s_2plot,
                         fliersize=0.5, orient="v", linewidth=0.5, ax=ax)

    if collect_outliers:
        if outlier_bounds is not None:
            # check that it's a tuple of two numbers
            assert (len(outlier_bounds) == 2) and isinstance(outlier_bounds[0], float), \
                f'Wrong format for outlier bounds: {outlier_bounds}. Must be a touple of two ascending floats.'
            lower, upper = outlier_bounds
            top_outliers = DF_subset[(DF_subset[f"{which_s}hat"] > upper)]
            bottom_outliers = DF_subset[(DF_subset[f"{which_s}hat"] < lower)]
            # set y range
            ax.set_ylim(lower, upper)
        else:
            # collect everything with abs. larger than 0.1
            top_outliers = DF_subset[DF_subset[f"{which_s}hat"] > 0.1]
            bottom_outliers = DF_subset[DF_subset[f"{which_s}hat"] < -0.1]
        ymin, ymax = ax.get_ylim()
        # annotate their counts on the plot
        if note_outliers:
            ## find the x coords of all the vertical lines, which are the central axes of each box
            # print( [set(l.get_xdata()) for l in ax.get_lines()])
            all_x_coords = [l.get_xdata()[0] for l in ax.get_lines() if len(set(l.get_xdata())) == 1]
            all_x_coords = sorted(list(set(all_x_coords)))
            ### sanity check (we should have len(Hs) x len(s) boxes)
            assert len(all_x_coords) == (len(H_2plot) * len(true_s_2plot)), f'\nlen(all_x_coords) = {len(all_x_coords)},' \
                                                                            f'len(H_2plot) = {len(H_2plot)},' \
                                                                            f'len(true_s_2plot) = {len(true_s_2plot)},' \
                                                                            f'{all_x_coords}'
            ### match them with H & s; might be the stupidest way but as long as it works...
            matched_xCoords = []
            counter = 0
            for i, s in enumerate(true_s_2plot):
                for j, H in enumerate(H_2plot):
                    assert (i * len(H_2plot) + j) == counter, f'i * len(H_2plot) + j = {i} * {len(H_2plot)} + {j}\n' \
                                                              f'counter = {counter}'
                    matched_xCoords.append([s, H, all_x_coords[counter]])
                    counter += 1
            matched_xCoords = pd.DataFrame(matched_xCoords, columns=[f'true_{which_s}', 'H', 'xCoord'])
            ### assign H & s to these numbers
            ## top
            if top_outliers.shape[0] > 0:
                top_outliers_counts = pd.pivot_table(top_outliers, values='rep',
                                                     index=[f'true_{which_s}', 'H'],
                                                     aggfunc=len).reset_index()
                top_outliers_counts = top_outliers_counts.rename(columns={'rep': 'counts'})
                # merge to get their xcoords
                top_outliers_counts = top_outliers_counts.merge(matched_xCoords, how='left',
                                                                on=[f'true_{which_s}', 'H'])
                # now annotate
                # print(top_outliers_counts)
                if top_outliers_counts.shape[0] > 5:
                    noise = np.random.normal(0, 2.5e-4, top_outliers_counts.shape[0])
                    noise = np.abs(noise)
                else:
                    noise = [0] * top_outliers_counts.shape[0]
                for i, row in top_outliers_counts.iterrows():
                    # row: true_s, H, xCoord, count
                    tag = row.counts
                    # print(tag)
                    this_x = row.xCoord
                    ax.annotate(text=str(tag), xy=(this_x, ymax),
                                xytext=(this_x - noise[i]*(-1)**i, ymax + 0.145*(ymax - ymin) + noise[i]*(-1)**i ),
                                arrowprops=dict(arrowstyle="->", color="#66666622"), fontsize='x-small',
                                bbox=dict(boxstyle="round", pad=0.15, edgecolor='#aaaaaa66', fc="#ffffff88"), ha='center', va='center',
                                color='black', annotation_clip=True)
            # bottom
            if bottom_outliers.shape[0] > 0:
                bottom_outliers_counts = pd.pivot_table(bottom_outliers, values='rep',
                                                        index=[f'true_{which_s}', 'H'],
                                                        aggfunc=len).reset_index()
                bottom_outliers_counts = bottom_outliers_counts.rename(columns={'rep': 'counts'})
                # merge to get their xcoords
                bottom_outliers_counts = bottom_outliers_counts.merge(matched_xCoords, how='left',
                                                                      on=[f'true_{which_s}', 'H'])
                labels = r''
                # now annotate
                for i, row in bottom_outliers_counts.iterrows():
                    # row: true_s, H, xCoord, count
                    tag = row.counts
                    # print(tag)
                    this_x = row.xCoord
                    ax.annotate(text=str(tag), xy=(this_x, ymin),
                                xytext=(this_x, ymin - 0.14*(ymax - ymin)),  # + noise[i]*(-1)**i
                                arrowprops=dict(arrowstyle="->", color="#66666622"), fontsize='x-small',
                                bbox=dict(boxstyle="round", pad=0.125, edgecolor='#aaaaaa66', fc='#ffffff88'), ha='center', va='center',
                                color='black', annotation_clip=True)
    else:
        empty_frame = pd.DataFrame([], columns=meltedDF_shat.columns)
        top_outliers, bottom_outliers = empty_frame, empty_frame
    outlier_df = pd.concat([top_outliers, bottom_outliers], axis=0, ignore_index=True)

    ax.grid(visible=True, axis='y', lw=0.3)
    # add aux lines
    left_x, right_x = ax.get_xlim()
    ax.hlines(y=list(map(float, true_s_2plot)), xmin=left_x, xmax=right_x,  # label=true_s_2plot,
              colors=['black'] * len(true_s_2plot), linestyle='dashed', lw=0.5)
    ax.tick_params(axis='both', labelsize="x-small")
    return ax, outlier_df


# prepare to extract
trueH_s2_regex = re.compile(r'h(\.?[0-9]+)x\ds')
trueH_s1_regex = re.compile(r'simH(\.?[0-9]+|\w*)x\ds')
llH_regex = re.compile(r'LLfixH(\d?.?\d+)_')
simInit_regex = re.compile(r'simInit([a-z|A-Z]+\.?[0-9]*)_')
llInit_regex = re.compile(r'LLinit([a-z|A-Z]+)_')
fixS2_regex = re.compile(r'LLfixS2-(\d+)_')


def load_pkl_MLEs(s2_pkl_name, s1_pkl_name):
    # load data
    metaS_hats = {'s2': {'Unif': {}, 'Freq': {}}, 's1': {'Unif': {}, 'Freq': {}}}
    keynames = ('s2', 's1')
    llInits = []
    simInits = []
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
            simInits.append(simInits_temp[0])

            if key_s == 's1':
                simH = trueH_s1_regex.findall(pkl_name)[0]
                if 'LLfixH' in pkl_name:
                    llH = llH_regex.findall(pkl_name)[0]
                    # print(f'Loading {key_s} simulations from simInit = {simInit},'
                    #       f'simH = {simH}, ll_init = {llInit}, fixed ll H = {llH}.')
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
                # print(f'Loading {key_s} simulations from simInit = {simInit},'
                #       f'simH = {simH}, ll_init = {llInit}, llH = {llH}.')
                if float(llH) >= 1 or float(llH) == 0:
                    llH = str(int(float(llH)))

            # initialize
            if llH not in metaS_hats[key_s][llInit]:
                metaS_hats[key_s][llInit][llH] = {}
            # each will be a dictionary with llH-decided s_pairs as key
            with open(pkl_name, 'rb') as pk:
                temp_LRs, temp_Shats, temp_LLs = pickle.load(pk)
                # metaLRs[key_s][llInit][llH], metaLLs[key_s][llInit][llH], metaS_hats[key_s][llInit][llH] = pickle.load(pk)
            pk.close()
            metaS_hats[key_s][llInit][llH][simH] = temp_Shats
            # print(llH, simH)
            # print(metaLRs[key_s][llInit].keys(), llH, metaLRs[key_s][llInit][llH].keys())  #, metaLRs[key_s][llInit][llH][simH].keys()
            assert ('0', '0') in metaS_hats[key_s][llInit][llH][simH], metaS_hats[key_s][llInit][llH][simH].keys()
        # out of for-loop, move on to next
    return metaS_hats, list(set(simInits)), list(set(llInits))


def main():
    s2_pkl_name = sys.argv[1]
    s1_pkl_name = sys.argv[2]
    simInit = sys.argv[3]
    # llInit = sys.argv[4]
    figname_prefix = sys.argv[4]
    # done loading file. Now plot
    H_2plot = (["0", "0.5", "1"], ["5", "Inf"])  # "0.5", "0", "1",".5",
    # H_2plot = ([0.5, 0., 1.], [5., "Inf"])  # "0.5", "0", "1",
    H_colors = {".5": plt.cm.tab10(0),  # blue
                "0.5": plt.cm.tab10(0),  # in case it doesn't recognize
                "0": plt.cm.tab10(1),  # orange
                "1": plt.cm.tab10(2),  # green
                "5": plt.cm.tab10(3),  # red
                "Inf": plt.cm.tab10(4),
                "fix_s1": plt.cm.tab10(5)}  # purple
    Dom = {".5": "additive", "0.5": "additive", '0': "complete recessive", "1": "complete dominant", "5": "over-dominant",
           "Inf": r"$s_{AA}=0$, over-dominant"}
    true_s_2plot = (["0", "0.001", "0.002", "0.003", "0.004", "0.005"],  #
                    ["0", "0.001", "0.002", "0.003", "0.004", "0.005"])  # , "0.01", "0.02"

    metaS_hats, simInits, llInits = load_pkl_MLEs(s2_pkl_name, s1_pkl_name)
    assert len(simInits) == 1, simInits

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

    def _express_H(H):
        if H == "Inf":
            return Dom[H]
        else:
            return f'h = {float(H):g}, {Dom[H]}'

    # llInits = list(set(llInits))
    for i, llInit in enumerate(llInits):
        fig, (left_ax, right_ax) = plt.subplots(ncols=2, gridspec_kw={'width_ratios': (0.5, 0.3)},
                                                figsize=(7, 3.5))  # , sharey=TruemetaS_hats
        if 'max' in s2_pkl_name:
            # take absolute values
            s2_absTag = True
        else:
            s2_absTag = False
        left_ax, s2_outliers = plot_MLE_boxplot_varS(left_ax, metaS_hats['s2'][llInit], H_2plot[0], true_s_2plot[0],
                                                     which_s='s2', absolute_MLEs=s2_absTag, colors=H_colors,
                                                     collect_outliers=True, note_outliers=True,
                                                     outlier_bounds=(-0.01, 0.015))  #
        # double check axis labels
        left_ax.set_ylabel(None)
        left_ax.set_title(r'Inferred $~\widehat{s}_{AA}$', loc='center', fontsize='small')
        left_ax.set_title('A)', loc='left', fontsize='medium', fontweight='bold') #
        left_ax.set_xlabel(r'True $s_{AA}$', fontsize='medium')
        left_ax.tick_params(axis='both', labelsize="x-small")
        # get legend handles
        handles_L, labels_L = left_ax.get_legend_handles_labels()

        labels_L = [_express_H(h) for h in labels_L]
        # print(labels_L)
        left_ax.legend(handles_L, labels_L, framealpha=0, fontsize='x-small', frameon=False,
                       title=None, loc='best')  # , title_fontsize='small'
        ymin, ymax = left_ax.get_ylim()
        print(ymin, ymax)
        if 'max' in s1_pkl_name:
            # take absolute values
            s1_absTag = True
        else:
            s1_absTag = False
        right_ax, s1_outliers = plot_MLE_boxplot_varS(right_ax, metaS_hats['s1'][llInit], H_2plot[1], true_s_2plot[1],
                                                      which_s='s1', absolute_MLEs=s1_absTag, colors=H_colors,
                                                      collect_outliers=True, note_outliers=True,
                                                      outlier_bounds=(-0.01, 0.015))  #
        # double check axis labels
        right_ax.set_ylabel(None)
        right_ax.set_title(r'Inferred $~\widehat{s}_{Aa}$', loc='center', fontsize='small')
        right_ax.set_title('B)', loc='left', fontsize='medium', fontweight='bold') #
        right_ax.set_xlabel(r'True $s_{Aa}$', fontsize='medium')
        right_ax.tick_params(axis='both', labelsize="x-small")
        # get legend handles
        handles_R, labels_R = right_ax.get_legend_handles_labels()
        labels_R = [_express_H(h) for h in labels_R]
        # print(labels_R)
        right_ax.set_ylim(ymin, ymax)

        # add shared legend & title <- nope
        # handles, labels = list(handles_L) + list(handles_R), list(labels_L) + list(labels_R)
        right_ax.legend(handles_R, labels_R, framealpha=0, fontsize='x-small', frameon=False,
                        title=None, loc='best')  # , title_fontsize='small'
        fig.suptitle(f'Simulation initialized w. {_express_init(simInit)}.  Likelihood computed w. {_express_init(llInit)}.',
                     fontsize='medium', ha='center', va='top', y=0.9) #
        plt.tight_layout()
        fig.savefig(f'{figname_prefix}_LLinit{llInit}_MLE_twoBoxes.png', dpi=350, transparent=True)

        # save outliers
        outliers = s2_outliers.merge(s1_outliers, how="outer", on=["H", "true_s1", "true_s2", "rep", "s1hat", "s2hat", "ongrid_shat"])
        # only write outliers when they exist...
        if outliers.shape[0] > 0:
            outliers['simH'] = outliers.true_s1.astype(float) / outliers.true_s2.astype(float)
            print(outliers.head())
            outliers.to_csv(f'{figname_prefix}_LLinit{llInit}_MLE_outliers.tsv', sep="\t", na_rep='NA', index_label='rep')


if __name__ == '__main__':
    main()
