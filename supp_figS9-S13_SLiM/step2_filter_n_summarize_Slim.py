#!/usr/bin/env python -W ignore::DeprecationWarning
"""collect LRs of 1) selected site and 2) site with highest MLR on each replicate
Usage:
python %prog <minMAF> <name_prefix> [s1]
"""
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

import sys, time, os, re
import numpy as np
import pandas as pd

rep_regex = re.compile(r'rep([0-9]+)-')
pos_regex = re.compile(r'-([0-9]{4,})$')
int_regex = re.compile(r'\s(\d+)\s')
frac_regex = re.compile(r'(0?\.\d+)')
trueH_regex = re.compile(r'_h(\.?[0-9]+)_')


def _split_rep_ID(IDstring):
    """return tuple (rep, pos)"""
    (rep, pos) = IDstring.split('_')
    return rep, pos


def _extract_rep_selPos_fromCounts(filename):
    """read ms-formatted allele count, extract the pos of selected locus from formatter"""
    with open(filename, 'r') as counts:
        for line_idx, line in enumerate(counts):
            if line_idx == 1:
                # t0, sel_pos, tk = int_regex.findall(line)
                ints = [n for n in int_regex.findall(line) if n != '']
                assert (len(ints) == 3) or (
                            len(ints) == 4), f'Extract {len(ints)} digit strings ({ints}) from line {line_idx + 1}:' \
                                             f'{line}\n{filename}'
                sel_pos = int(ints[1])
            elif line_idx > 1:
                break
            else:
                continue
    counts.close()
    return sel_pos


# need testing
def _extract_metaInfo_fromCounts(filename):
    """read ms-formatted allele count, extract info from formatter
    Basically an expanded `_extract_rep_selPos_fromCounts()`
    """
    with open(filename, 'r') as counts:
        for line_idx, line in enumerate(counts):
            if line_idx == 1:
                # t0, sel_pos, tk = int_regex.findall(line)
                ints = int_regex.findall(line)
                assert len(ints) == 3, f'Extract {len(ints)} digit strings from line {line_idx + 1}: {line}'
                sel_pos = int(ints[1])
            elif line_idx == 2:
                Traj_idx = int_regex.findall(line)
                Traj_idx = list(map(int, Traj_idx))
            elif line_idx == 3:
                Traj_ct = int_regex.findall(line)
                Traj_ct = list(map(int, Traj_ct))
            elif line_idx > 4:
                break
            else:
                continue
    counts.close()
    # sanity check
    assert len(Traj_idx) == len(Traj_ct), f'# generations: {len(Traj_idx)}, # recorded counts: {len(Traj_ct)}'
    Traj = zip(Traj_idx, Traj_ct)
    return sel_pos, Traj


from diplo_locus.utility import parse_allele_counts


def read_slim_counts_and_filter(filename: str, masks, minMAF: float = 0):
    """read ms-formatted allele count, apply filters, and return list of positions
    `masks` is a list (or other iterables) of strings
    """
    sel_pos = _extract_rep_selPos_fromCounts(filename)
    # outdated interface
    # pos_names, samples, sampleSizes = parse_allele_counts(filename, minMAF=minMAF)
    pos_names, samples, sampleSizes = parse_allele_counts(filename)
    # the above 3 already passed poolMAF filter (if any)
    # initialize mask
    passed_rows = np.array(np.zeros(len(pos_names)), dtype=bool)
    # now filter more
    sample_freqs = samples / sampleSizes
    # filter by MAF <-- just to doublecheck
    # THIS IS ACTUALLY NECESSARY, SINCE parse_allele_counts() doesn't do it anymore
    pool_freqs = samples.sum(axis=1) / sampleSizes.sum(axis=1)
    passed_freqs = (pool_freqs > minMAF)
    passed_rows |= passed_freqs
    # "notAllFixed" filter is built-in in slim's parsing script
    for mask in masks:
        if mask == "2morePoly":
            # if mask.endswith("morePoly"):
            #     P = int(int_regex.findall(mask)[0])
            is_poly = (sample_freqs > 0) & (sample_freqs < 1)
            num_poly = is_poly.sum(axis=1)
            passed_rows &= (num_poly >= 2)
        elif mask in ["1morePoly", "notMono", "hasPoly"]:  # more synonyms to add (?)
            is_poly = (sample_freqs > 0) & (sample_freqs < 1)
            num_poly = is_poly.sum(axis=1)
            passed_rows &= (num_poly > 0)
        elif mask in ['', ' ']:
            continue
        else:
            print("Unrecognized mask:", mask)
            # no filter, so do nothing

    pos_list = list(map(int, pos_names[passed_rows]))

    return pos_list, sel_pos


def extract_MLRs_for_poslist(mlrfilename, pos_list, sel_pos, on_off_grid, cand_range=None, seq_len=None):
    """read maxLLs.txt file, extract MLR and MLE of selected sites and maxLLR site
    `cand_range` is a tuple of two ints for positions
    return two 4-element lists, `SelPos_pdRow` & `maxPos_pdRow`: [position, ongrid_shat, shat, MLR]
    """
    scores = pd.read_csv(mlrfilename, sep="\t", comment="#")
    # getting MLR
    if on_off_grid == 'on':  # 'ongrid_maxLogLikelihood': 'MLR',
        # assert 'maxLogLikelihood' in scores.columns, scores.columns
        scores.columns = ["position", "ongrid_shat", "ongrid_maxLogLikelihood", "MLR"]  # , "s2hat", "maxLogLikelihood"
        # neutral LL:
        neutLL = scores.maxLogLikelihood - (scores.MLR / 2)
        ongrid_MLR = (scores.ongrid_maxLogLikelihood - neutLL) * 2

        scores.rename(columns={'# locus': 'position'},  # 'MLR': 'offgrid_MLR', , 'ongrid_s2hat': 's2hat'
                      inplace=True)
        scores['MLR'] = ongrid_MLR

    elif on_off_grid == 'off':  # 'maxLogLikelihood': 'MLR',
        if scores.shape[1] == 6:
            scores.columns = ["position", "ongrid_shat", "ongrid_maxLogLikelihood", "shat", "maxLogLikelihood", "MLR"]
        elif scores.shape[1] == 7:
            h = trueH_regex.findall(mlrfilename)[0]
            assert float(h) == 0.5, mlrfilename
            scores.columns = ["position", "ongrid_shat", "ongrid_maxLogLikelihood", "shat", "maxLogLikelihood", "MLR",
                              "chi2_p"]
        else:
            raise ValueError(f"Invalid number of columns: {filename}:\n{scores.columns}")
    # take subset
    scores.position = scores.position.astype(int)
    # record sel_pos first
    if (sel_pos - 1) in pos_list:
        scores = scores[scores.position.isin(pos_list)]
        SelPos_pd = scores[scores.position == (sel_pos - 1)].reindex()
        # print(SelPos_pd)
        SelPos_pdRow = [int(SelPos_pd.position), float(SelPos_pd.ongrid_shat), float(SelPos_pd.shat),
                        float(SelPos_pd.MLR)]
        # SelPos_pdRow = [SelPos_pd.position[0], SelPos_pd.ongrid_shat[0], SelPos_pd.shat[0], SelPos_pd.MLR[0]]
        # print(SelPos_pdRow)
    elif sel_pos in pos_list:
        scores = scores[scores.position.isin(pos_list)]
        SelPos_pd = scores[scores.position == sel_pos].reindex()
        SelPos_pdRow = [int(SelPos_pd.position), float(SelPos_pd.ongrid_shat), float(SelPos_pd.shat),
                        float(SelPos_pd.MLR)]
        # SelPos_pdRow = [SelPos_pd.position[0], SelPos_pd.ongrid_shat[0], SelPos_pd.shat[0], SelPos_pd.MLR[0]]
        # print(SelPos_pdRow)
    else:
        SelPos_pdRow = None
        pdPos_peri = [pos for pos in scores.position if (pos < sel_pos + 500) and (pos > sel_pos - 500)]
        pos_list_peri = [pos for pos in pos_list if (pos < sel_pos + 500) and (pos > sel_pos - 500)]
        print(f'{sel_pos} or {sel_pos - 1} not in pos_list {pos_list_peri}')
        raise ValueError(f'{sel_pos} or {sel_pos - 1} not in the dataset. pos_list = {pdPos_peri}, {pos_list_peri}\n'
                         f'len(pos_list) = {len(pos_list)}, len(scores.position) = {len(scores.position)}\n'
                         f'type(pos_list[0])={type(pos_list[0])}; type(sel_pos)={type(sel_pos)}')
    # pick maxLR site next
    if cand_range is not None:
        if isinstance(cand_range, tuple):
            assert len(cand_range) == 2, cand_range
            scores = scores[(scores.position >= int(cand_range[0])) & (scores.position <= int(cand_range[1]))]
        else:
            assert seq_len is not None, 'Must provide sequence length if \"cand_range\" is set to be the radius (single value)'
            # left_bnd, right_bnd = min([sel_pos - cand_range, 0]), max([sel_pos + cand_range, scores.position.max()])
            # force the 180k to be the subset
            left_bnd, right_bnd = (sel_pos - cand_range), (sel_pos + cand_range)
            if left_bnd < 0 or right_bnd > seq_len:
                raise ValueError(f'target-centered window ({sel_pos} +- {cand_range})'
                                 f'fall outside of the {seq_len / 1e3:g}kb sequence:\n{mlrfilename}')
            scores = scores[(scores.position >= int(left_bnd)) & (scores.position <= int(right_bnd))]
    scores = scores.reindex()  # range(1, scores.shape[0]+1)
    maxPos_pd = scores.loc[scores['MLR'].idxmax()]
    # print(scores.loc[scores['MLR'].idxmax()])
    maxPos_pdRow = [maxPos_pd.position, maxPos_pd.ongrid_shat, maxPos_pd.shat, maxPos_pd.MLR]
    return SelPos_pdRow, maxPos_pdRow


def main():
    # wdir = 'diplolocus_manuscript_figs/figS8-S12_SLiM/'
    wdir = os.getcwd()
    # hack this in
    if (wdir[-1] != '/'):
        wdir = f"{wdir}/"

    repN = 200
    on_off_grid = "off"
    # llInit = sys.argv[1] # Unif or Freq
    llInit = "Unif"
    poolMAF = sys.argv[1]
    poolMAF = float(poolMAF)
    assert 0 <= poolMAF < 1
    filter_tag = f'_minMAF{str(poolMAF)[1:]}'
    name_prefix = sys.argv[2]

    if 's1' in sys.argv:
        which_s = 's1'
        s2_list = ['.0', '.0002', '.0004', '.0006', '.0008', '.001']
        # h_list = ['5', '.5']
        h_list = ['5']
    else:
        which_s = 's2'
        s2_list = ['.0', '.001', '.002', '.003', '.004', '.005']
        h_list = ['0', '.5', '1']

    # collect max scores for each scenario and write to file
    meta_maxFile = f'{wdir}{name_prefix}-var{which_s.upper()}_maxLR-sites{filter_tag}_200kb_t9x500gen_n40_{len(h_list)}h{len(s2_list)}s_Concatenated_{on_off_grid}Grid_Dist.tsv.gz'
    meta_selFile = f'{wdir}{name_prefix}-var{which_s.upper()}_selected-sites{filter_tag}_200kb_t9x500gen_n40_{len(h_list)}h{len(s2_list)}s_Concatenated_{on_off_grid}Grid_Dist.tsv.gz'
    if not os.path.exists(meta_maxFile) or not os.path.exists(meta_selFile) or ('recount' in sys.argv):
        # prep two containers, one for known selected site, another for maxMLR site
        Selected_sites = pd.DataFrame(columns=['Condition', 's2', 'h', 'rep', 'position', 'ongrid_shat', 'shat', 'MLR'])
        maxMLR_sites = pd.DataFrame(columns=['Condition', 's2', 'h', 'rep', 'position', 'ongrid_shat', 'shat', 'MLR'])
        counter = 0
        # read all data
        for s in s2_list:
            for h in h_list:
                if 's1' in sys.argv:
                    if (h != '5') and (f's{s}_h{h}' != 's.0_h.5'):
                        continue
                    # elif f's{s}_h{h}' == 's.0_h5':
                    #     cond = f's.0_h.5'  #
                    #     print(cond)
                    #     # batch_selFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/{name_prefix}_maxMLRs-selPos{filter_tag}_200kb_{cond}_fixH{h}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                    #     # batch_maxFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/{name_prefix}_maxMLRs-maxPos{filter_tag}_200kb_{cond}_fixH{h}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                    #     batch_selFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/{name_prefix}_maxMLRs-selPos{filter_tag}_200kb_{cond}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                    #     batch_maxFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/{name_prefix}_maxMLRs-maxPos{filter_tag}_200kb_{cond}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                    else:
                        cond = f's{s}_h{h}'  #
                        print(cond)
                        # batch_selFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/{name_prefix}_maxMLRs-selPos{filter_tag}_200kb_{cond}_fixH5_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                        # batch_maxFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/{name_prefix}_maxMLRs-maxPos{filter_tag}_200kb_{cond}_fixH5_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                        batch_selFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/{name_prefix}_maxMLRs-selPos{filter_tag}_200kb_{cond}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                        batch_maxFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/{name_prefix}_maxMLRs-maxPos{filter_tag}_200kb_{cond}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                # elif f's{s}_h{h}' not in ['s.0_h0', 's.0_h1']:
                else:
                    cond = f's{s}_h{h}'  # , 's.0_h5'
                    print(cond)
                    # batch_selFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/DL_maxMLRs-selPos{filter_tag}_200kb_{cond}_fixH{h}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                    # batch_maxFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/DL_maxMLRs-maxPos{filter_tag}_200kb_{cond}_fixH{h}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                    batch_selFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/DL_maxMLRs-selPos{filter_tag}_200kb_{cond}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                    batch_maxFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/DL_maxMLRs-maxPos{filter_tag}_200kb_{cond}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                # else:
                #     assert f's{s}_h{h}' in ['s.0_h0', 's.0_h1', 's.0_h5'], f's{s}_h{h}'
                #     cond = f's.0_h.5'  #
                #     print(cond)
                #     # batch_selFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/DL_maxMLRs-selPos{filter_tag}_200kb_{cond}_fixH{h}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                #     # batch_maxFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/DL_maxMLRs-maxPos{filter_tag}_200kb_{cond}_fixH{h}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                #     batch_selFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/DL_maxMLRs-selPos{filter_tag}_200kb_{cond}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                #     batch_maxFile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/DL_maxMLRs-maxPos{filter_tag}_200kb_{cond}_51xgeom75e-2_{llInit}_{repN}reps.tsv'
                # print(batch_maxFile, batch_selFile)
                batch_Selected_sites = pd.DataFrame(columns=['rep', 'position', 'ongrid_shat', 'shat', 'MLR'])
                batch_maxMLR_sites = pd.DataFrame(columns=['rep', 'position', 'ongrid_shat', 'shat', 'MLR'])
                batch_counter = 0
                for i in range(repN):
                    countfile = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/count/HC_{cond}_ConstN_t9x500gen_n40_rep{i}.count'
                    # `masks` be a list of strings representing different masks
                    # return: sel_pos, max-MLR pos, list of pos (to keep)
                    try:
                        pos_list, sel_pos = read_slim_counts_and_filter(countfile, masks=[" "],
                                                                        minMAF=poolMAF)  # 2morePoly
                    except FileNotFoundError:  # {countfile}
                        print(f'{cond} rep {i} sample count file doesn\'t exist! Skip.\n{countfile}')
                        # continue
                        sys.exit()
                    # sanity check
                    try:
                        assert ((sel_pos - 1) in pos_list) or (sel_pos in pos_list), \
                            f'{cond}, rep {i}, sel_pos={sel_pos} not in pos_list:\n' \
                            f'{[p for p in pos_list if (p > sel_pos - 500) and (p < sel_pos + 500)]}'
                        sel_target_exist = True
                    except AssertionError:
                        print(
                            f'{cond} rep {i} most probably need restarting. sel_pos not recorded in SLiM samples')
                        sel_target_exist = False

                    # if f's{s}_h{h}' in ['s.0_h0', 's.0_h1', 's.0_h5']:
                    #     # sanity check
                    #     assert cond == f's.0_h.5', cond
                    #     mlrfilename = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/likelihood/HC_{cond}_ConstN_t9x500gen_n40_rep{i}_DL_MAF{str(poolMAF)[1:]}_fixH{h}_51xgeom75e-2_{llInit}_off-grid_maxLLs.txt'
                    # else:
                    #     mlrfilename = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/likelihood/HC_{cond}_ConstN_t9x500gen_n40_rep{i}_DL_MAF{str(poolMAF)[1:]}_51xgeom75e-2_{llInit}_off-grid_maxLLs.txt'
                    mlrfilename = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/likelihood/HC_{cond}_ConstN_t9x500gen_n40_rep{i}_DL_MAF{str(poolMAF)[1:]}_51xgeom75e-2_{llInit}_off-grid_maxLLs.txt'

                    if not os.path.exists(mlrfilename):
                        # HC_s.0_h.5_newSV_ConstN_t9x500gen_n40_rep195_DL_MAF.05_fixH1_51xgeom75e-2_Unif_off-grid_maxLLs.txt
                        # HC_s.0_h.5_newSV_ConstN_t9x500gen_n40_rep150_DL_MAF.05_51xgeom75e-2_Unif_off-grid_maxLLs.txt
                        # mlrfilename = f'{wdir}200kb_HC_{cond}_ConstN_t9x500gen/likelihood/HC_{cond}_ConstN_t9x500gen_n40_rep{i}_DL_51xgeom75e-2_{llInit}_off-grid_maxLLs.txt'.split("/")[-1]
                        print(f'{cond} rep {i} maxLL file doesn\'t exist! Skip.\n\t{mlrfilename}')
                        # continue
                        sys.exit()
                    # only read LLs of sites in `pos_list`, return `max_pos`, `ongrid_shat`, `ongrid_MLR`, `shat`, `MLR` as pd.Series
                    # else:
                    try:
                        SelPos_pdRow, maxPos_pdRow = extract_MLRs_for_poslist(mlrfilename, pos_list, sel_pos,
                                                                              on_off_grid, cand_range=9e4, seq_len=2e5)
                    except Exception as e:
                        print(e)
                        print(f'Abnormal {cond} rep {i} maxLL file. Skip.\n\t{mlrfilename.split("/")[-1]}')
                        continue

                    # add to existing DF
                    if sel_target_exist:
                        assert SelPos_pdRow is not None
                        temp_new_row = {'rep': i, 'position': SelPos_pdRow[0],
                                        'ongrid_shat': SelPos_pdRow[1],
                                        'shat': SelPos_pdRow[2], 'MLR': SelPos_pdRow[3]}
                        batch_Selected_sites.loc[batch_Selected_sites.shape[0]] = temp_new_row

                    batch_maxMLR_sites.loc[batch_counter] = {'rep': i, 'position': maxPos_pdRow[0],
                                                             'ongrid_shat': maxPos_pdRow[1],
                                                             'shat': maxPos_pdRow[2],
                                                             'MLR': maxPos_pdRow[3]}
                    # increment
                    batch_counter += 1
                    counter += 1
                # save to file
                try:
                    batch_Selected_sites.to_csv(batch_selFile, sep="\t", index=False)
                    batch_maxMLR_sites.to_csv(batch_maxFile, sep="\t", index=False)
                except Exception as e:
                    print(e)
                    print(cond, batch_selFile, batch_maxFile)
                # integrate
                for DF in (batch_maxMLR_sites, batch_Selected_sites):
                    DF['Condition'] = cond
                    DF['s2'] = s
                    DF['h'] = h
                    # DF = DF.reindex(columns=['Condition', 's2', 'h', 'rep', 'position', 'ongrid_shat', 'shat', 'MLR'])
                # out of loop for Cond, save all
                Selected_sites = pd.concat([Selected_sites, batch_Selected_sites], ignore_index=True)
                maxMLR_sites = pd.concat([maxMLR_sites, batch_maxMLR_sites], ignore_index=True)
                # sanity check
                assert maxMLR_sites.shape[0] == counter, \
                    f'maxMLR_sites.shape[0]={maxMLR_sites.shape[0]}\ncounter={counter}'
                print("Done with ", cond, 'LL_H =', h, time.ctime())
        # End of looping through all s & h, save to file
        Selected_sites = Selected_sites.reindex(
            columns=['Condition', 's2', 'h', 'rep', 'position', 'ongrid_shat', 'shat', 'MLR'])
        Selected_sites.to_csv(meta_selFile, sep="\t", index=False)
        maxMLR_sites.reindex(columns=['Condition', 's2', 'h', 'rep', 'position', 'ongrid_shat', 'shat', 'MLR']).to_csv(
            meta_maxFile, sep="\t", index=False)
    else:
        assert os.path.isfile(meta_maxFile)
        assert os.path.isfile(meta_selFile)
        Selected_sites = pd.read_csv(meta_selFile, sep="\t")
        maxMLR_sites = pd.read_csv(meta_maxFile, sep="\t")
        if not isinstance(Selected_sites.position[0], int):
            print(Selected_sites.position[0], type(Selected_sites.position[0]))
            Selected_sites = pd.read_csv(meta_selFile, sep="\t",
                                         dtype={'Condition': str, 's2': float, 'h': float, 'rep': int,
                                                'position': np.int64, 'ongrid_shat': np.float64, 'shat': np.float64,
                                                'MLR': np.float64})
        if not isinstance(maxMLR_sites.position[0], int):
            maxMLR_sites = pd.read_csv(meta_maxFile, sep="\t",
                                       dtype={'Condition': str, 's2': float, 'h': float, 'rep': int,
                                              'position': np.int64, 'ongrid_shat': np.float64,
                                              'shat': np.float64, 'MLR': np.float64})
    print(Selected_sites.describe())
    print(maxMLR_sites.describe())


if __name__ == '__main__':
    main()
