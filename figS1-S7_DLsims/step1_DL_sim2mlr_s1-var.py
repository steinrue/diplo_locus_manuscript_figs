"""Pipeline to run simulations with variable s_Aa with both initial conditions.
fix s2 value, vary s1
Usage:
    python %prog <outprefix> <neut_samples_pkl> <init_distn> <minMAF> [seed]
"""
# set up
import sys, os, random, pickle
import numpy as np
from diplo_locus.utility import _get_geom_grid
from DL_sim_wrappers import simulate_and_calculate_by_s1


def main():
    outprefix = sys.argv[1]
    neut_samps_file = sys.argv[2]
    distn_file = sys.argv[3]

    minMAF = float(sys.argv[4])

    if len(sys.argv) > 5:
        seed = int(sys.argv[5])
    else:
        seed = random.randrange(int(1e6))  # sys.maxsize
    # set seed
    print("Seed:", seed)
    random.seed(seed)
    np.random.seed(seed)

    # define parameters:
    ## basics
    # u01 = 1.25e-8
    # u10 = 1.25e-8
    u01 = 0.
    u10 = 0.
    Ne = 10000
    numRep = 500

    ## for simulations
    # key parameter
    true_H_list = ['5', 'Inf']
    # others are the same as s2
    true_s1_list = ('0.001', '0.002', '0.003', '0.004', '0.005', '0.01', '0.02')  # , 0.006, 0.007, 0.008, 0.0090,
    sampleTimes = [0, 0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]  #
    K = len(sampleTimes) - 1
    sampleSizes = [40] * K
    sim_deltaT = 0.5
    sim_init_freq = 0.01

    ## for computing LL
    gridPoints = 50
    LL_deltaT = 1
    LL_init_freq = 0.01
    LL_s1_grid = _get_geom_grid(-0.75, 0.75, gridPoints, Ne)

    # action
    for simInit in ('Freq.01', 'SV'): #
        # metaSamples, metaThrownOuts, metaLRs, metaLLs, metaShats = {}, {}, {}, {}, {}
        if simInit == 'Freq.01':
            # metaSamples['Freq.01'], metaThrownOuts['Freq.01'], metaLRs['Freq.01'], metaLLs['Freq.01'], metaShats['Freq.01'] = {}, {}, {}, {}, {}
            metaSamples, metaThrownOuts, metaLRs, metaLLs, metaShats = {}, {}, {}, {}, {}
            # LL_H as keys
            for LL_H in true_H_list:
                metaSamples[LL_H], metaThrownOuts[LL_H], metaLRs[LL_H], metaLLs[LL_H], metaShats[LL_H] = \
                    simulate_and_calculate_by_s1(Ne, u01, u10, numRep, seed,
                                                 true_H_list, true_s1_list, sampleTimes, sampleSizes,
                                                 out_prefix=outprefix, sim_init_cond="init_freq",
                                                 sim_init_freq=sim_init_freq, minMAF=minMAF, sim_deltaT=sim_deltaT,
                                                 neut_samps=neut_samps_file, LL_init_freq=LL_init_freq,
                                                 LL_s_grid=LL_s1_grid, LL_H=LL_H,
                                                 LL_s_grid_tag=f'_{len(LL_s1_grid)}xgeomGrid75e-2',
                                                 LL_deltaT=LL_deltaT)
            # big dump
            lump_pkl_file = f'{outprefix}_simInitFreq{str(sim_init_freq)[1:]}_seed{seed}_simH-{"-".join(true_H_list)}_{len(true_s1_list)}xS1_t{K}n{sampleSizes[0]}_{len(LL_s1_grid)}xgeomGrid75e-2_big-lump.pkl'
            with open(lump_pkl_file, 'wb') as lump_pkl:
                pickle.dump((metaSamples, metaThrownOuts, metaLRs, metaLLs, metaShats), lump_pkl)
            lump_pkl.close()
        elif simInit == 'SV':
            # metaSamples['SV'], metaThrownOuts['SV'], metaLRs['SV'], metaLLs['SV'], metaShats['SV'] = {}, {}, {}, {}, {}
            metaSamples, metaThrownOuts, metaLRs, metaLLs, metaShats = {}, {}, {}, {}, {}
            for LL_H in true_H_list:
                # get init distn for sims
                assert os.path.exists(distn_file), f'Initial distribution file {distn_file} invalid'
                init_distn = np.loadtxt(distn_file)
                metaSamples[LL_H], metaThrownOuts[LL_H], metaLRs[LL_H], metaLLs[LL_H], metaShats[LL_H] = \
                    simulate_and_calculate_by_s1(Ne, u01, u10, numRep, seed,
                                                 true_H_list, true_s1_list, sampleTimes, sampleSizes,
                                                 out_prefix=outprefix, sim_init_cond="standing_var",
                                                 sim_init_distn=init_distn, minMAF=minMAF, sim_deltaT=sim_deltaT,
                                                 neut_samps=neut_samps_file, LL_init_freq=LL_init_freq,
                                                 LL_s_grid=LL_s1_grid, LL_H=LL_H,
                                                 LL_s_grid_tag=f'_{len(LL_s1_grid)}xgeomGrid75e-2',
                                                 LL_deltaT=LL_deltaT)
            # big dump
            lump_pkl_file = f'{outprefix}_simInitSV_seed{seed}_simH-{"-".join(true_H_list)}_{len(true_s1_list)}xS1_t{K}n{sampleSizes[0]}_{len(LL_s1_grid)}xgeomGrid75e-2_big-lump.pkl'
            with open(lump_pkl_file, 'wb') as lump_pkl:
                pickle.dump((metaSamples, metaThrownOuts, metaLRs, metaLLs, metaShats), lump_pkl)
            lump_pkl.close()


if __name__ == '__main__':
    main()
