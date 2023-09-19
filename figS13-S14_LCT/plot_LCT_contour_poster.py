import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def load_surface(tablename, ID_col: str = "ID"):
    """return the matrix and its rolnames & colnames"""
    onGrid_LLs = pd.read_csv(tablename, delimiter="\t", comment="#")  # , index_col=0, header=[0]
    # assume the DF has both row names and column names
    if ID_col in list(onGrid_LLs.columns):
        loci_names = onGrid_LLs[ID_col].to_list()
        # s_pairs = list(onGrid_LLs.columns[1:])
        s_pairs = [colname for colname in onGrid_LLs.columns if colname != ID_col]
    else:
        loci_names = list(onGrid_LLs.index)
        s_pairs = list(onGrid_LLs.columns)
    # make sure they're numbers
    # print(s_pairs)
    try:
        s_pairs = [list(map(float, pair.strip("()").split(','))) for pair in s_pairs]
    except Exception as e:
        print(e)
        print(s_pairs[:5])
        print(list(onGrid_LLs.columns).index(ID_col), len(s_pairs))
        print(s_pairs[-5:])
        sys.exit()
    s1_list, s2_list = map(lambda x: np.array(sorted(list(set(x)))),
                           zip(*s_pairs))
    LLcontainer = np.array(onGrid_LLs.loc[:, onGrid_LLs.columns != ID_col])
    # make sure s_pairs have tupples
    s_pairs = [tuple(pair) for pair in s_pairs]
    return np.array(LLcontainer), np.array(loci_names), s1_list, s2_list, s_pairs, len(s_pairs)


# blegh
def _read_all_float(num: str):
    num = num.strip()
    if num == '':
        raise ValueError(f'num = \"{num}\"')
    elif num[0] == '-':
        return -1. * float(num[1:])
    else:
        assert num.isnumeric(), f'\"{num}\".isnumeric() == False'
        return float(num)


def main():
    # read commandline args
    inputfile, outprefix = sys.argv[1:]
    # SNPs = pd.read_csv(inputfile, delimiter="\t", header=0, index_col=0, comment="#")
    LLcontainer, loci_names, s1_list, s2_list, s_pairs, num_pairs = load_surface(inputfile, ID_col="ID")

    for i, snp in enumerate(loci_names):
        if 'rs4988235' not in snp:
            continue

        surface = np.array(LLcontainer[i, :].reshape((len(s2_list)), len(s1_list)))

        # crop it to -0.2 to 0.2
        which_x_idx = np.where((s2_list >= -0.15) & (s2_list <= 0.2))[0]
        which_y_idx = np.where((s1_list >= -0.15) & (s1_list <= 0.2))[0]
        s2_list_zoomin = s2_list[which_x_idx]
        s1_list_zoomin = s1_list[which_y_idx]
        surface_zoomin = (surface[which_x_idx, :])[:, which_y_idx]
        print(surface_zoomin.shape, len(s2_list_zoomin))

        # find neut
        neut_i, neut_j = list(s2_list_zoomin).index(0), list(s1_list_zoomin).index(0)
        snp_neut = surface_zoomin[neut_i, neut_j] / np.log(10)
        max_idx = np.unravel_index(np.argmax(surface_zoomin), surface_zoomin.shape)
        snp_peak = surface_zoomin[max_idx] / np.log(10)
        s1_0, s2_0 = s1_list_zoomin[max_idx[1]], s2_list_zoomin[max_idx[0]]
        print(snp, max_idx, s1_0, s2_0, snp_peak)
        # now let's plot
        figname = f'{outprefix}_{snp}_contour.png'
        fig, ax = plt.subplots(figsize=(4.5, 4))  #
        ct = ax.contour(s1_list_zoomin, s2_list_zoomin, surface_zoomin / np.log(10), levels=70, cmap='coolwarm')  # copperinfernoReds_r
        ax.clabel(ct, ct.levels[1::2], inline=True, fontsize='x-small')
        # title
        mlr = (snp_peak - snp_neut) * 2  # <-- this is already log10
        ax.set_title(r"$LCT$/$MCM6$, rs4988235" + ', ' +  # '\n' +
                     r"$log_{10}\mathcal{L.R.}_{max}\approx$" + f'{mlr:.3f}', fontsize='small')
        # ax.set_axis_on()
        ax.set_xlim(-0.15, 0.2)
        ax.set_xlabel('$s_{Aa}$', fontsize='small')
        # print(ax.get_xticklabels())
        print(ax.get_xticks())
        ax.set_ylim(-0.15, 0.2)
        ax.set_ylabel('$s_{AA}$', fontsize='small')
        # print(ax.get_yticklabels())
        print(ax.get_yticks())
        ax.tick_params(axis='both', labelsize="x-small", width=0.5, direction='inout')
        # add secondary axes
        ax_x = ax.twiny()
        ax_x.set_xlabel(r'$\sigma_{Aa}$', fontsize='small', color='#360000')
        ax_x.set_xlim(-0.15, 0.2)
        # ax_x.set_xticks(ax.get_xticks())
        ax_x.set_xticklabels([f'{s * 1e4:g}' for s in ax.get_xticks()])
        ax_x.tick_params(axis='x', labelsize="x-small", color='#360000', labelcolor='#360000')

        ax_y = ax.twinx()
        ax_y.set_ylabel(r'$\sigma_{AA}$', fontsize='small', color='#360000', rotation=270)
        ax_y.set_ylim(-0.15, 0.2)
        # ax_y.set_yticks(ax.get_yticks())
        # ax_y.set_yticklabels([f'{_read_all_float(lab.get_text())*1e4:g}' for lab in ax.get_yticklabels()])
        ax_y.set_yticklabels([f'{s * 1e4:g}' for s in ax.get_yticks()])
        ax_y.tick_params(axis='y', labelsize="x-small", color='#360000', labelcolor='#360000')

        # annotate on-grid max
        ax.plot(s1_0, s2_0, 'x', c='#360000')
        # ax.set_title(f'%s on chr%s, %s\nLR=%.4f', size='small')
        # ax.annotate('$log_{10}\mathcal{L}_{max}=$%.4f\n(%.3f, %.3f)' % (s1_0, s2_0, snp_peak),
        ax.annotate(r'$(\hat{s}_{Aa}, \hat{s}_{AA})\approx$' + '(%.3f, %.3f)\n' % (s1_0, s2_0) +
                    r'$log_{10}\mathcal{L}_{max}=$' + '%.4f' % snp_peak,
                    (s1_0, s2_0), xytext=(s1_0, s2_0-0.01), ha='center', va='top', fontsize='x-small',
                    bbox=dict(facecolor='w', alpha=0.75, pad=0.5, linewidth=0.),
                    annotation_clip=True)  # textcoords='offset points',
        # add axis line
        ax.axvline(x=0, ymin=0, ymax=1, ls='--', color='#888888', lw=0.5)
        ax.axhline(y=0, xmin=0, xmax=1, ls='--', color='#888888', lw=0.5)
        ax.axline((0, 0), (0.01, .01), ls='-.', c='darkgray', lw=0.7)
        ax.axline((0, 0), (0.005, .01), ls='-.', c='darkgray', lw=0.7)
        # annotate h
        ax.annotate('Complete dominance\n($h=1$)', xy=(.15, .15), xytext=(0.17, 0.12), color='gray',
                    fontsize='xx-small', ha='center', va='top', #textcoords='offset points',
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.4', color='gray', lw=0.75),
                    bbox=dict(boxstyle="round", alpha=0.4, pad=0.25, lw=0., fc='w'), annotation_clip=True)

        ax.annotate('Semi-dominance\n(additive, $h=0.5$)', xy=(.07, .14), xytext=(-0.02, 0.1), color='gray',
                    fontsize='xx-small', ha='right', va='top', #textcoords='offset points',
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.25', color='gray', lw=0.75),
                    bbox=dict(boxstyle="round", alpha=0.4, pad=0.25, lw=0., fc='w'), annotation_clip=True)

        ax.annotate('Complete recessive\n($h=0$)', xy=(0, -0.1), xytext=(0.12, -0.1), color='gray',
                    fontsize='xx-small', ha='center', va='bottom', #textcoords='offset points',
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2', color='gray', lw=0.75),
                    bbox=dict(boxstyle="round", alpha=0.4, pad=0.25, lw=0., fc='w'), annotation_clip=True)
        # adjust border
        ax.spines['top'].set_color('#360000')
        ax_x.spines['top'].set_color('#360000')
        ax.spines['right'].set_color('#360000')
        ax_y.spines['right'].set_color('#360000')

        for this_ax in [ax, ax_x, ax_y]:
            for side in this_ax.spines.keys():
                this_ax.spines[side].set_linewidth(0.3)

        # ax.set_xlim(-0.1, 0.102)
        fig.tight_layout()
        plt.savefig(figname, dpi=500, transparent=True)


if __name__ == '__main__':
    main()
