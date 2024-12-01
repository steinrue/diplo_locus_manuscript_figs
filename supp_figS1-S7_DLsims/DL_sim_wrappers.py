"""Wrappers around `diplo_locus` functions to set up DL simulations for the manuscript

"""
# set up
import sys, time, os, pickle
import numpy as np
from datetime import datetime, timedelta
import diplo_locus.likelihood as likelihood
import diplo_locus.simulate as simulate
from diplo_locus.utility import _get_geom_grid, _find_1D_Max_perSite, get_onGrid_max_only


# pretty much the first part of `_find_1D_Max_perSite`
def _get_onGrid_Max_perSite(s_grid, ith_LLs_persite):
    idx, LLs_persite = ith_LLs_persite[0], ith_LLs_persite[1:]
    # get on-grid max:
    ongrid_max_idx = np.argmax(LLs_persite)
    s_0 = s_grid[ongrid_max_idx]
    ongrid_max = LLs_persite[ongrid_max_idx]
    return [s_0, ongrid_max, np.nan, np.nan]


def interpolate_offgrid_max(s1_list, s2_list, s_pairs, num_pairs, LLmatrix):
    # decide how many s to optimize & define the function
    s1_gridnum = len(set(s1_list))
    s2_gridnum = len(set(s2_list))
    if min(s1_gridnum, s2_gridnum) > 3 and (s1_gridnum * s2_gridnum) == num_pairs:
        # 2d-grid
        print(f'DiploLocus do not support 2D optimization for now. Please re-specify the parameter grid.')
        return False
    # decide that it's 1d
    elif num_pairs in {s1_gridnum, s2_gridnum}:
        # fix h
        if s1_gridnum == num_pairs and s2_gridnum == num_pairs:
            var_list = list(zip(*s_pairs))[1]
            header = '\t'.join(['locus', 'ongrid_s2hat', 'ongrid_maxLogLikelihood', 's2hat', 'maxLogLikelihood'])
        # fix s1
        elif s1_gridnum == 1:
            var_list = list(zip(*s_pairs))[1]
            header = '\t'.join(['locus', 'ongrid_s2hat', 'ongrid_maxLogLikelihood', 's2hat', 'maxLogLikelihood'])
        # fix s2
        elif s2_gridnum == 1:
            var_list = list(zip(*s_pairs))[0]
            header = '\t'.join(['locus', 'ongrid_s1hat', 'ongrid_maxLogLikelihood', 's1hat', 'maxLogLikelihood'])
        else:
            print(f's1_gridnum={s1_gridnum}; s2_gridnum={s2_gridnum}; num_pairs={num_pairs}')
            print(f's1_list={s1_list}\ns2_list={s2_list}\ns_pairs={s_pairs}.')
            sys.exit()

        var_list = np.array(list(map(float, var_list)))
        bnds = (min(var_list), max(var_list))

        def _get_max_persite(x):
            try:
                wrapper = _find_1D_Max_perSite(np.array(var_list), ith_LLs_persite=x, bounds=bnds)
            except Exception as e:
                print(e, var_list, bnds)
                print(header, num_pairs)
                print(s1_list, s2_list)
                # return on-grid max instead
                wrapper = _get_onGrid_Max_perSite(np.array(var_list), ith_LLs_persite=x)
            return wrapper

        # add row index as a col
        iLLmatrix = np.insert(LLmatrix, 0, np.arange(1, LLmatrix.shape[0] + 1), axis=1)
        # maxLLs = np.apply_along_axis(_get_max_persite, 1, iLLmatrix)
    else:
        print(f'''Insufficient grid points for effective interpolation. 
        At least 4 different values are needed for either `s1` or `s2`.
        s1_list={s1_list};\ns2_list={s2_list}\ns_pairs={s_pairs}.''')
        sys.exit()

    maxLLs = np.apply_along_axis(_get_max_persite, 1, iLLmatrix)

    return maxLLs, header


# ["init_freq", "standing_var"]
# values, probs = init_distn[:, 0], init_distn[:, 1]
def simulate_samples(Ne: int or float, u01: float, true_s_pair_list, sample_times, sampleSizes,
                     init: str, init_freq: float, init_distn, numRep: int, seed: int,
                     u10: float = None, minMAF: float = 0., sim_deltaT: float = 0.05):
    K = len(sample_times) - 1
    Scenarios_varS = {}
    thrownOut = {}
    # for s2 in true_s_list:
    #     s1 = H * s2
    for (s1, s2) in true_s_pair_list:
        samples = np.empty((numRep, K))
        trajs = np.empty((numRep, int((sample_times[-1] - sample_times[0]) / sim_deltaT) + 2))
        if init == 'init_freq':
            def _itermulate(numreps, itr_seed):
                thisSample, thisTraj = simulate.simulateSamples(Ne, float(s1), float(s2), u01, u10, sample_times,
                                                                sampleSizes, itr_seed, initCond='initFreq',
                                                                numReplicates=numreps,
                                                                initFreq=init_freq, deltaT=sim_deltaT)
                return thisSample, thisTraj

        elif init == 'standing_var':
            # assert distn_file is not None
            def _itermulate(numreps, itr_seed):
                values, probs = init_distn[:, 0], init_distn[:, 1]
                thisSample, thisTraj = simulate.simulateSamples(Ne, float(s1), float(s2), u01, u10, sample_times,
                                                                sampleSizes, itr_seed, initCond='choice',
                                                                initValues=values, initProbs=probs,
                                                                numReplicates=numreps,
                                                                deltaT=sim_deltaT)
                return thisSample, thisTraj
        else:
            raise ValueError(f"{init} is not accepted as an init")

        r = 0
        seed_thisItr = seed
        while r < numRep:
            thisSample, thisTraj = _itermulate(numRep - r, seed_thisItr)
            # check this here
            # ## more than two time points polymorphic
            # thisSampleFreq = thisSample / sampleSizes
            # polymorph_samples = (0 < thisSampleFreq) & (thisSampleFreq < 1)
            # toRemove = (np.array(polymorph_samples, dtype=int).sum(axis=1) < 2)
            ## last sample nonzero
            # toRemove = np.any(thisSample[:, -1] == np.zeros((1, 1)), axis=0)
            ## not all fixed (sample-wise)
            # toRemove |= (thisSample.sum(axis=1) == np.array(sampleSizes).sum())
            # not lost throughout
            toRemove = np.any(thisTraj[:, 2:] < 1 / (4 * Ne), axis=1)
            # print(toRemove.sum())

            # filter for pooled MAF regardless of minMAF value ## sampleSizes only have 1 dimension here
            ## sample freqs
            pooledMAF = np.sum(thisSample, axis=1) / np.sum(sampleSizes)
            ## fold
            pooledMAF = np.where(pooledMAF > 0.5, 1 - pooledMAF, pooledMAF)
            maf_filter = (pooledMAF <= minMAF)
            # print(maf_filter.sum())
            toRemove |= maf_filter
            # print(toRemove.sum())

            # get the mask
            # toKeep = np.array(1 - toRemove, dtype="bool")
            toKeep = np.array(~toRemove, dtype=bool)
            numToKeep = toKeep.sum()
            if numToKeep <= (numRep - r):
                assert sum(toRemove) == numRep - r - numToKeep
                if (str(s1), str(s2)) not in thrownOut: thrownOut[(str(s1), str(s2))] = 0
                thrownOut[(str(s1), str(s2))] += sum(toRemove)
                thisSample = thisSample[toKeep, :]
                thisTraj = thisTraj[toKeep, :]
            # store it
            try:
                samples[r:(r + numToKeep), :] = thisSample
                trajs[r:(r + numToKeep), :] = thisTraj
            except Exception as e:
                print(e)
                print(toRemove.shape, sum(toRemove), toKeep.shape, sum(toKeep))
                print(thisSample[toKeep, :].shape)
                print(r, samples[r:(r + numToKeep), :].shape, thisSample.shape)
                print(trajs[r:(r + numToKeep), :].shape, thisTraj.shape)
                assert False
            # update r: number of passed reps
            r += numToKeep
            # update seed
            seed_thisItr += 1
        # end of while loop, save
        Scenarios_varS[(str(s1), str(s2))] = (samples, trajs)
    # check out
    print(thrownOut)
    return Scenarios_varS, thrownOut


def compute_likelihood_1D(sample_times, sampleSizes, s_pairs, s_grid, Scenarios_varS, init: str, init_freq: float,
                          numRep: int or float, Ne: int or float, u01: float, u10: float, LL_deltaT: int or float = 1.):
    K = len(sample_times) - 1
    # obtain the wrapper function based on initCond args
    if init == 'Unif':
        def _computeLL(sampleSizes, samples, s1, s2):
            HMMcore = likelihood.SelHmm(Ne, s1, s2, u01, u10, initCond="uniform",
                                        sampleSizesSet=set((np.array(sampleSizes).reshape(numRep * K)).astype('int')),
                                        emissionType="integer", transitionType='constant', deltaT=LL_deltaT)
            LLmatrix = HMMcore.computeLogLikelihood(sample_times, sampleSizes.astype('int'), samples.astype('int'))
            return LLmatrix
    elif init == 'Freq':
        assert init_freq is not None

        def _computeLL(sampleSizes, samples, s1, s2):
            HMMcore = likelihood.SelHmm(Ne, s1, s2, u01, u10, initCond="initFreq", initFreq=init_freq,
                                        sampleSizesSet=set((np.array(sampleSizes).reshape(numRep * K)).astype('int')),
                                        emissionType="integer", transitionType='constant', deltaT=LL_deltaT)
            LLmatrix = HMMcore.computeLogLikelihood(sample_times, sampleSizes.astype('int'), samples.astype('int'))
            return LLmatrix
    else:
        return False

    LRs_varS = {}
    sHat_varS = {}
    LLsurfaces_varS = {}
    num_pairs = len(s_grid)
    # make sure it's 1D
    assert num_pairs == len(s_pairs)
    s1_grid, s2_grid = zip(*s_pairs)
    # make sure they're not strings
    s1_grid = np.array(list(map(float, s1_grid)))
    s2_grid = np.array(list(map(float, s2_grid)))

    # now loop through samples to get LRs
    print("Simulated scenarios:", Scenarios_varS.keys())
    for true_s_pair in Scenarios_varS.keys():
        print(time.ctime(), true_s_pair)
        samples_here = Scenarios_varS[true_s_pair][0]
        sampleSizes_here = np.ones((numRep, 1)) * np.array(sampleSizes).reshape(1, K)
        print(K, samples_here.shape, sampleSizes_here.shape)
        LLcontainer = np.empty((numRep, num_pairs))
        for i, s_pair in enumerate(s_pairs):
            s1, s2 = s_pair
            LL_array = _computeLL(sampleSizes_here, samples_here, float(s1), float(s2))
            LLcontainer[:, i] = LL_array
        LLsurfaces_varS[true_s_pair] = LLcontainer
        # get neut
        neut_idx = [i for i, s_pair in enumerate(s_pairs) if np.isclose(sum([float(s) for s in s_pair]), 0)]
        neut_idx = neut_idx[0]
        neutLL = LLcontainer[:, neut_idx]
        # neutLL = LLcontainer[:, s_pairs.index((0,0))]
        # now get max
        ## header (if s2): 'locus', 'ongrid_s2hat', 'ongrid_maxLogLikelihood', 's2hat', 'maxLogLikelihood'
        ## header (if s1): 'locus', 'ongrid_s1hat', 'ongrid_maxLogLikelihood', 's1hat', 'maxLogLikelihood'
        maxLLs, header = interpolate_offgrid_max(s1_grid, s2_grid, s_pairs, num_pairs, LLcontainer)
        MLRs = (maxLLs[:, -1] - neutLL) * 2
        MLRs = MLRs.reshape(numRep)
        # sanity check: MLR > 0:
        try:
            assert np.all(MLRs >= 0), f'## Negative MLR: {MLRs[np.where(MLRs < 0)]}. Replace with on-grid MLR instead.'
        except AssertionError:
            to_replace = np.where(MLRs < 0)
            # replace off-grid MLR with on-grid ones
            ongrid_MLRs = (maxLLs[:, 1] - neutLL) * 2
            ongrid_MLRs = ongrid_MLRs.reshape(numRep)
            MLRs[to_replace] = ongrid_MLRs[to_replace]
            # replace off-grid MLEs with on-grid MLEs
            maxLLs[to_replace, -2] = maxLLs[to_replace, 0]
            maxLLs[to_replace, -1] = maxLLs[to_replace, 1]
        LRs_varS[true_s_pair] = MLRs # <-- should already be on-grid ones if MLEs are on-grid 
        # sHat_varS[true_s_pair] = maxLLs[:, :2] #<-- these are on-grid ones!
        sHat_varS[true_s_pair] = maxLLs[:, (0, -2)]  # (on-grid, off-grid) <-- get the off-grid MLEs

    return LRs_varS, sHat_varS, LLsurfaces_varS


def print_time(duration):
    return time.strftime("%H hr %M min %S sec", time.gmtime(duration))


def print_timedelta(duration: timedelta):
    output_str = f'{duration.days} days {print_time(duration.total_seconds())}'
    return output_str


def get_neut_samp_file_name(out_prefix, seed, Ne, u01, sim_init, K, sampleSize, numRep,
                            init_freq: float = 0., minMAF: float = 0., filter_tag: str = 'notLost'):
    # filter_tag = 'lastNotLost-notFixed'
    if sim_init == "standing_var":
        sim_init_tag = '_simInitSV'
    else:
        assert sim_init == "init_freq", f'init_cond \"{sim_init}\" not supported.'
        assert init_freq > 0
        sim_init_tag = f'_simInitFreq{str(init_freq)[1:]}'

    if minMAF > 0:
        maf_tag = '_MAF' + str(minMAF)[1:]
    else:
        maf_tag = ''

    neut_samps_file = f'{out_prefix}_seed{seed}_N{Ne}-u{u01}{sim_init_tag}_Neut{maf_tag}_{filter_tag}_t{K}n{sampleSize}_{numRep}reps_samps-n-trajs.pkl'

    return neut_samps_file


def get_sample_file_names(out_prefix, seed, Ne, u01, sim_init, ture_s_coeff_tag, K, sampleSize, numRep,
                          init_freq: float = 0., minMAF: float = 0., filter_tag: str = ''):  # 2morePoly
    # filter_tag = 'lastNotLost-notFixed'
    if sim_init == "standing_var":
        sim_init_tag = '_simInitSV'
    else:
        assert sim_init == "init_freq", f'init_cond \"{sim_init}\" not supported.'
        assert init_freq > 0
        sim_init_tag = f'_simInitFreq{str(init_freq)[1:]}'

    if minMAF > 0:
        maf_tag = '_MAF' + str(minMAF)[1:]
    else:
        maf_tag = ''

    if filter_tag == '':
        sample_file = f'{out_prefix}_seed{seed}_N{Ne}-u{u01}{sim_init_tag}_{ture_s_coeff_tag}{maf_tag}_t{K}n{sampleSize}_{numRep}reps_samples-n-trajs.pkl'
        forLL_abbr = f'{out_prefix}_fromSeed{seed}{sim_init_tag}_{ture_s_coeff_tag}{maf_tag}_t{K}n{sampleSize}_{numRep}reps'
    else:
        sample_file = f'{out_prefix}_seed{seed}_N{Ne}-u{u01}{sim_init_tag}_{ture_s_coeff_tag}{maf_tag}_{filter_tag}_t{K}n{sampleSize}_{numRep}reps_samples-n-trajs.pkl'
        forLL_abbr = f'{out_prefix}_fromSeed{seed}{sim_init_tag}_{ture_s_coeff_tag}{maf_tag}_{filter_tag}_t{K}n{sampleSize}_{numRep}reps'

    return sample_file, forLL_abbr


def get_LL_file_names(prefix_from_sim, LL_init, grid_tag):
    # grid_tag = f'_{len(s2_grid)}x{gridType}Grid75e-2'
    LL_pkl_file = f'{prefix_from_sim}_LLinit{LL_init}{grid_tag}_offGridMLE.pkl'
    return LL_pkl_file


def simulate_and_calculate_by_s2(Ne: int or float, u01: float, u10: float, numRep: int, seed: int,
                                 true_H_list, true_s_list, sampleTimes, sampleSizes, out_prefix: str,
                                 sim_init_cond: str, sim_init_freq=None, sim_init_distn=None, minMAF: float = 0.,
                                 sim_deltaT: float = .5, neut_samps: str = '', LL_init_freq=None, LL_s_grid=None,
                                 LL_s_grid_tag: str = '', LL_H: float = 0.5,
                                 LL_deltaT: float = 1.):  # LL_init_cond: str = '' ,
    K = len(sampleTimes) - 1
    # decide whether to generate or read neutral:
    # decide neut filenames first
    neut_samps_default = get_neut_samp_file_name(out_prefix, seed, Ne, u01, sim_init=sim_init_cond,
                                                 K=K, sampleSize=sampleSizes[0], numRep=numRep,
                                                 init_freq=sim_init_freq, minMAF=minMAF)
    ## read neut if not simulated
    simulate_neut = True
    if neut_samps != '':
        if os.path.exists(neut_samps):
            simulate_neut = False
            neut_samps_file = neut_samps
        elif os.path.exists(neut_samps_default):
            simulate_neut = False
            neut_samps_file = neut_samps_default
        else:
            print(f'File {neut_samps} does not exist. Simulate a new batch with seed {seed}')
            neut_samps_file = neut_samps_default
    else:
        print('Simulate neutral replicates from scratch. Seed:', seed)
        neut_samps_file = neut_samps_default

    if not simulate_neut:
        assert os.path.isfile(neut_samps_file), f'Neutral file {neut_samps_file} invalid.'
        print('Loading pre-simulated neutral replicates from ', neut_samps_file)
        with open(neut_samps_file, 'rb') as npk:
            neutScenarios, neutThrownOuts = pickle.load(npk)
        npk.close()
    else:
        # simulate a fresh batch
        print(time.ctime(), 'Simulating neutral replicates...')
        last = datetime.now()
        last_time = str(time.ctime())
        neutScenarios, neutThrownOuts = simulate_samples(Ne, u01, true_s_pair_list=[('0', '0')],
                                                         sample_times=sampleTimes, sampleSizes=sampleSizes,
                                                         init=sim_init_cond, init_freq=sim_init_freq,
                                                         init_distn=sim_init_distn, numRep=numRep, seed=seed,
                                                         u10=u10, minMAF=minMAF, sim_deltaT=sim_deltaT)
        print("from ", last_time)
        print(time.ctime())
        print(f'Simulating {numRep} neutral reps took {print_timedelta(datetime.now() - last)}.')
        # save
        print(f'Saving neutral samples to {neut_samps_file}.')
        with open(neut_samps_file, 'wb') as npk:
            pickle.dump((neutScenarios, neutThrownOuts), npk)
        npk.close()

    # sanity check
    assert ('0', '0') in neutScenarios.keys(), neutScenarios.keys()
    # print(neutScenarios.keys())
    neutScenarios = neutScenarios[('0', '0')]
    neutThrownOuts = neutThrownOuts[('0', '0')]

    metaSamples, metaThrownOuts = {}, {}
    metaLRs, metaLLs, metaS2hats = {'Unif': {}, 'Freq': {}}, {'Unif': {}, 'Freq': {}}, {'Unif': {}, 'Freq': {}}
    for sim_H in true_H_list:
        h = float(sim_H)
        true_s_pair_list = [(f'{float(s) * h:g}', s) for s in true_s_list]
        this_sample_pkl, this_forLL_abbr = get_sample_file_names(out_prefix, seed, Ne, u01, sim_init=sim_init_cond,
                                                                 ture_s_coeff_tag=f'h{sim_H}x{len(true_s_list)}s',
                                                                 K=K, sampleSize=sampleSizes[0],
                                                                 numRep=numRep, init_freq=sim_init_freq, minMAF=minMAF)
        if not os.path.exists(this_sample_pkl):
            last = datetime.now()
            print(time.ctime(),
                  f'Simulating  h={sim_H}, s \\in {true_s_list}.Starting with standard Watternson neutral.')
            metaSamples[sim_H], metaThrownOuts[sim_H] = simulate_samples(Ne, u01, true_s_pair_list,
                                                                         sample_times=sampleTimes,
                                                                         sampleSizes=sampleSizes,
                                                                         init=sim_init_cond, init_freq=sim_init_freq,
                                                                         init_distn=sim_init_distn, numRep=numRep,
                                                                         seed=seed, u10=u10, minMAF=minMAF,
                                                                         sim_deltaT=sim_deltaT)
            print(
                f"Simulation completed in {print_timedelta(datetime.now() - last)} with seed {seed}.\n Saving to {this_sample_pkl}")
            with open(this_sample_pkl, "wb") as pk:
                pickle.dump((metaSamples[sim_H], metaThrownOuts[sim_H]), pk)
            # make sure it closes
            pk.close()
        # read it in if already exist
        else:
            print(f'Loading pre-computed replicates from {this_sample_pkl}')
            with open(this_sample_pkl, 'rb') as pk:
                metaSamples[sim_H], metaThrownOuts[sim_H] = pickle.load(pk)
            pk.close()

        if ('0', '0') not in metaSamples[sim_H].keys():
            metaSamples[sim_H][('0', '0')] = neutScenarios
            metaThrownOuts[sim_H][('0', '0')] = neutThrownOuts
        print("s_list for h=", h, metaSamples[sim_H].keys())

        # then compute LLs
        LL_s_pairs = [(f'{s * float(LL_H):g}', f'{s:g}') for s in LL_s_grid]
        for LL_initType in ('Unif', 'Freq'):
            LL_pkl = get_LL_file_names(f'{this_forLL_abbr}_LLfixH{LL_H}', LL_init=LL_initType, grid_tag=LL_s_grid_tag)
            if not os.path.exists(LL_pkl):
                print(f'\n{time.ctime()} compute likelihoods on {LL_s_grid_tag} grid with {LL_initType} initCond.')
                last = datetime.now()
                metaLRs[LL_initType][sim_H], metaS2hats[LL_initType][sim_H], metaLLs[LL_initType][sim_H] = \
                    compute_likelihood_1D(sampleTimes, sampleSizes, s_pairs=LL_s_pairs, s_grid=LL_s_grid,
                                          Scenarios_varS=metaSamples[sim_H], init=LL_initType, init_freq=LL_init_freq,
                                          numRep=numRep, Ne=Ne, u01=u01, u10=u10, LL_deltaT=LL_deltaT)
                print(f"Computing (LL_H={LL_H}) completed in {print_timedelta(datetime.now() - last)}.")
                print(f'Saving LLs of all true_H={h}, LL_H={LL_H} to {LL_pkl}')
                with open(LL_pkl, "wb") as pk:
                    pickle.dump(
                        (metaLRs[LL_initType][sim_H], metaS2hats[LL_initType][sim_H], metaLLs[LL_initType][sim_H]), pk)
                pk.close()
            else:
                metaLRs[LL_initType][sim_H], metaS2hats[LL_initType][sim_H], metaLLs[LL_initType][sim_H] = pickle.load(
                    open(LL_pkl, "rb"))
                # sanity check
                print('for true h=', sim_H)
                print('Samples[sim_H].keys()=', metaSamples[sim_H].keys())
                print('metaLRs[sim_H].keys()=', metaLRs[LL_initType][sim_H].keys())

    return metaSamples, metaThrownOuts, metaLRs, metaLLs, metaS2hats


def simulate_and_calculate_by_s1(Ne: int or float, u01: float, u10: float, numRep: int, seed: int,
                                 true_H_list, true_s_list, sampleTimes, sampleSizes, out_prefix: str,
                                 sim_init_cond: str, sim_init_freq=None, sim_init_distn=None, minMAF: float = 0.,
                                 sim_deltaT: float = .5, neut_samps: str = '', LL_init_freq=None, LL_s_grid=None,
                                 LL_s_grid_tag: str = '', LL_H: str = '5',
                                 LL_deltaT: float = 1.):  # LL_init_cond: str = '' ,
    K = len(sampleTimes) - 1
    # decide whether to generate or read neutral:
    # decide neut filenames first
    neut_samps_default = get_neut_samp_file_name(out_prefix, seed, Ne, u01, sim_init=sim_init_cond,
                                                 K=K, sampleSize=sampleSizes[0], numRep=numRep,
                                                 init_freq=sim_init_freq, minMAF=minMAF)
    ## read neut if not simulated
    simulate_neut = True
    if neut_samps != '':
        if os.path.exists(neut_samps):
            simulate_neut = False
            neut_samps_file = neut_samps
        elif os.path.exists(neut_samps_default):
            simulate_neut = False
            neut_samps_file = neut_samps_default
        else:
            print(f'File {neut_samps} does not exist. Simulate a new batch with seed {seed}')
            neut_samps_file = neut_samps_default
    else:
        print('Simulate neutral replicates from scratch. Seed:', seed)
        neut_samps_file = neut_samps_default

    if not simulate_neut:
        assert os.path.isfile(neut_samps_file), f'Neutral file {neut_samps_file} invalid.'
        print('Loading pre-simulated neutral replicates from ', neut_samps_file)
        with open(neut_samps_file, 'rb') as npk:
            neutScenarios, neutThrownOuts = pickle.load(npk)
        npk.close()
    else:
        # simulate a fresh batch of neut
        print(time.ctime(), 'Simulating neutral replicates...')
        last = datetime.now()
        last_time = str(time.ctime())
        neutScenarios, neutThrownOuts = simulate_samples(Ne, u01, true_s_pair_list=[('0', '0')],
                                                         sample_times=sampleTimes, sampleSizes=sampleSizes,
                                                         init=sim_init_cond, init_freq=sim_init_freq,
                                                         init_distn=sim_init_distn, numRep=numRep, seed=seed,
                                                         u10=u10, minMAF=minMAF, sim_deltaT=sim_deltaT)
        print("from ", last_time)
        print(time.ctime())
        print(f'Simulating {numRep} neutral reps took {print_timedelta(datetime.now() - last)}.')
        # save
        print(f'Saving neutral samples to {neut_samps_file}.')
        with open(neut_samps_file, 'wb') as npk:
            pickle.dump((neutScenarios, neutThrownOuts), npk)
        npk.close()

    # sanity check
    assert ('0', '0') in neutScenarios.keys(), neutScenarios.keys()
    # print(neutScenarios.keys())
    neutScenarios = neutScenarios[('0', '0')]
    neutThrownOuts = neutThrownOuts[('0', '0')]

    # containers
    metaSamples, metaThrownOuts = {}, {}
    metaLRs, metaS2hats, metaLLs = {'Unif': {}, 'Freq': {}}, {'Unif': {}, 'Freq': {}}, {'Unif': {}, 'Freq': {}}

    # before simulating, decide LL-computing parameters
    ## bc when len(s1_grid) == len(s2_grid), `interpolate_offgrid_max` assumes it's s2 x h, the shat reported be s2
    if LL_H.lower() == 'inf':
        LL_s_pairs = [(f'{s:g}', '0') for s in LL_s_grid]
        dom_tag = '_LLfixS2-0'
    else:
        try:
            LL_h = float(LL_H)
            LL_s_pairs = [(f'{s:g}', f'{s / LL_h:g}') for s in LL_s_grid]
            dom_tag = f'_LLfixH{LL_H}'
        except Exception as e:
            print(e)
            raise ValueError(f'Do not recognize LL_H: {LL_H}')

    for sim_H in true_H_list:
        # true_s_pair_list = [(f'{float(s) * h:g}', s) for s in true_s_list]
        # true_s_pair_list = [(s, f'{float(s) / sim_H:g}') for s in true_s_list]
        try:
            h = float(sim_H)
            true_s_pair_list = [(f'{float(s):g}', f'{float(s) / h:g}') for s in true_s_list]
        except ValueError:
            if sim_H.lower() == "inf":
                true_s_pair_list = [(f'{float(s):g}', '0') for s in true_s_list]
            else:
                raise ValueError(f'Do not recognize sim_H: {sim_H}')
        except Exception as e:
            print(e)
            raise ValueError(f'Do not recognize sim_H: {sim_H}')

        this_sample_pkl, this_forLL_abbr = get_sample_file_names(out_prefix, seed, Ne, u01, sim_init=sim_init_cond,
                                                                 ture_s_coeff_tag=f'simH{sim_H}x{len(true_s_list)}s1',
                                                                 K=K, sampleSize=sampleSizes[0],
                                                                 numRep=numRep, init_freq=sim_init_freq, minMAF=minMAF)
        if not os.path.exists(this_sample_pkl):
            last = datetime.now()
            print(time.ctime(),
                  f'Simulating  h={sim_H}, s1 \\in {true_s_list}.Starting with standard Watternson neutral.')
            metaSamples[sim_H], metaThrownOuts[sim_H] = simulate_samples(Ne, u01, true_s_pair_list,
                                                                         sampleTimes, sampleSizes,
                                                                         init=sim_init_cond, init_freq=sim_init_freq,
                                                                         init_distn=sim_init_distn, numRep=numRep,
                                                                         seed=seed, u10=u10, minMAF=minMAF,
                                                                         sim_deltaT=sim_deltaT)
            print(f"Simulation completed in {print_timedelta(datetime.now() - last)} with seed {seed}.\n"
                  f"Saving to {this_sample_pkl}")
            with open(this_sample_pkl, "wb") as pk:
                pickle.dump((metaSamples[sim_H], metaThrownOuts[sim_H]), pk)
            # make sure it closes
            pk.close()
        # read it in if already exist
        else:
            print(f'Loading pre-computed replicates from {this_sample_pkl}')
            with open(this_sample_pkl, 'rb') as pk:
                metaSamples[sim_H], metaThrownOuts[sim_H] = pickle.load(pk)
            pk.close()

        if ('0', '0') not in metaSamples[sim_H].keys():
            metaSamples[sim_H][('0', '0')] = neutScenarios
            metaThrownOuts[sim_H][('0', '0')] = neutThrownOuts
        print("s_list for h=", sim_H, metaSamples[sim_H].keys())

        # then compute LLs
        for LL_initType in ('Unif', 'Freq'):
            LL_pkl = get_LL_file_names(this_forLL_abbr + dom_tag, LL_init=LL_initType, grid_tag=LL_s_grid_tag)
            if not os.path.exists(LL_pkl):
                print(f'\n{time.ctime()} compute likelihoods on {LL_s_grid_tag} grid with {LL_initType} initCond.')
                last = datetime.now()
                last_time = str(time.ctime())
                res_placeholder = compute_likelihood_1D(sampleTimes, sampleSizes, s_pairs=LL_s_pairs,
                                                        s_grid=LL_s_grid, Scenarios_varS=metaSamples[sim_H],
                                                        init=LL_initType, init_freq=LL_init_freq, numRep=numRep,
                                                        Ne=Ne, u01=u01, u10=u10, LL_deltaT=LL_deltaT)
                tempLRs, tempShats, tempLLs = res_placeholder
                if 'H' in dom_tag:
                    assert LL_H.isdecimal(), f'LL_H={LL_H} is not decimal.'
                    # the reported s is s2, gotta convert it
                    H = float(LL_H)
                    try:
                        #tempShats *= H
                        for pair in tempShats.keys():
                            batch_shats = tempShats[pair]
                            batch_shats *= H
                            tempShats[pair] = batch_shats
                    except Exception as e:
                        print(e)
                        print(type(tempShats))
                        print(tempShats.keys())
                metaLRs[LL_initType][sim_H] = tempLRs
                metaS2hats[LL_initType][sim_H] = tempShats
                metaLLs[LL_initType][sim_H] = tempLLs
                print("from ", last_time)
                print("to ", time.ctime())
                print(f"Computing (sim_H={sim_H}) completed in {print_timedelta(datetime.now() - last)}.")
                print(f'Saving LLs of all true_H={sim_H} & LL_H={LL_H} to {LL_pkl}')
                with open(LL_pkl, "wb") as pk:
                    pickle.dump(res_placeholder, pk)
                pk.close()
            else:
                with open(LL_pkl, "rb") as pk:
                    res_placeholder = pickle.load(pk)
                pk.close()
                assert len(res_placeholder) == 3, f'len(res_placeholder) = {len(res_placeholder)}'
                tempLRs, tempShats, tempLLs = res_placeholder
                metaLRs[LL_initType][sim_H] = tempLRs
                metaS2hats[LL_initType][sim_H] = tempShats
                metaLLs[LL_initType][sim_H] = tempLLs
                # sanity check
                print('for true h=', sim_H)
                print('Samples[H].keys()=', metaSamples[sim_H].keys())
                print('metaLRs[H].keys()=', metaLRs[LL_initType][sim_H].keys())

    return metaSamples, metaThrownOuts, metaLRs, metaLLs, metaS2hats
