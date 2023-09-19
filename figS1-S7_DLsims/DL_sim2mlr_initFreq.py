"""Pipeline to run simulations initiated at the given frequency of 0.01
Usage:
    python %prog <outprefix> <neut_samples_pkl> <minMAF> [seed]
"""
# set up
import sys, random, pickle
import numpy as np
from diplo_locus.utility import _get_geom_grid
from DL_sim_wrappers import simulate_and_calculate_by_s2


def main():
    outprefix = sys.argv[1]
    neut_samps_file = sys.argv[2]

    minMAF = float(sys.argv[3])

    if len(sys.argv) > 4:
        seed = int(sys.argv[4])
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
    true_H_list = ('0', '.5', '1', '5')
    true_s2_list = ('0.001', '0.002', '0.003', '0.004', '0.005', '0.01', '0.02')  # , 0.006, 0.007, 0.008, 0.0090,
    sampleTimes = [0, 0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]  #
    K = len(sampleTimes) - 1
    sampleSizes = [40] * K
    sim_deltaT = 0.5
    sim_init_freq = 0.01

    ## for computing LL
    gridPoints = 50
    LL_deltaT = 1
    LL_init_freq = 0.01
    LL_s2_grid = _get_geom_grid(-0.75, 0.75, gridPoints, Ne)

    # action
    metaSamples, metaThrownOuts, metaLRs, metaLLs, metaS2hats = {}, {}, {}, {}, {}
    for LL_H in true_H_list:
        LL_h = float(LL_H)
        metaSamples[LL_H], metaThrownOuts[LL_H], metaLRs[LL_H], metaLLs[LL_H], metaS2hats[LL_H] = \
            simulate_and_calculate_by_s2(Ne, u01, u10, numRep, seed,
                                         true_H_list, true_s2_list, sampleTimes, sampleSizes,
                                         out_prefix=outprefix, sim_init_cond="init_freq",
                                         sim_init_freq=sim_init_freq, minMAF=minMAF, sim_deltaT=sim_deltaT,
                                         neut_samps=neut_samps_file, LL_init_freq=LL_init_freq,
                                         LL_s_grid=LL_s2_grid, LL_H=LL_h,
                                         LL_s_grid_tag=f'_{len(LL_s2_grid)}xgeomGrid75e-2',
                                         LL_deltaT=LL_deltaT)

    # big dump
    lump_pkl_file = f'{outprefix}_simInitFreq{str(sim_init_freq)[1:]}_seed{seed}_{len(true_H_list)}xH_{len(true_s2_list)}xS2_t{K}n{sampleSizes[0]}_MAF{str(minMAF)[1:]}_{len(LL_s2_grid)}xgeomGrid75e-2_big-lump.pkl'
    with open(lump_pkl_file, 'wb') as lump_pkl:
        pickle.dump((metaSamples, metaThrownOuts, metaLRs, metaLLs, metaS2hats), lump_pkl)
    lump_pkl.close()


if __name__ == '__main__':
    main()
