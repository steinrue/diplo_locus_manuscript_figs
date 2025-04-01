import numpy
import pandas
import pathlib
# custom modules
import execution_engine

metaRNG = numpy.random.default_rng (4711)




# how to execute
THE_EXECUTOR = execution_engine.FileExecutor ("SLURM_SIMULATE.txt", append=False)
# THE_EXECUTOR = execution_engine.SubprocessExecutor()
# where execultable
DIPLO_SIM_EXEC = "DiploLocus-simulate"
SIMULATION_DIR = pathlib.Path ("simulated_data")
# put the simulated data in the respective directory
RESULTS_DIR = pathlib.Path ("results_diplo_locus")
# fixed parameters for simulations
# sampling scheme
GEN_BETWEEN_SAMPLES = 500
NUM_SAMPLING_TIMES = 9
NUM_SAMPLES_EACH = 100
# diffusion
DIFF_M_ALPHA = 0
DIFF_M_BETA = DIFF_M_ALPHA
DIFF_NE = 10000
INIT_FREQ = 0.01
# selection
ALL_S = [f'{x:.4f}' for x in [0, 0.001, 0.002, 0.003, 0.004, 0.005]]
FIXED_H = "0.5"
SELECTION_ONSET = 2000
# simulation
DELTA_T = 0.1
# DELTA_T = 1
# NUM_REPLICATES = lambda s : 500 if numpy.isclose(s, 0) else 100
# NUM_REPLICATES = lambda s : 2000
NUM_REPLICATES = 500
# NUM_REPLICATES = 5
MAF_FILTER = 0.05




def simulate_all_data (thisRNG):

    # create directory for all the simulations
    # SIMULATION_DIR.mkdir()
    SIMULATION_DIR.mkdir (exist_ok=True)

    # onset or constant
    for onsetFlag in [True, False]:
        onsetName = "onset" if onsetFlag else "constant"

        # go though all selection
        for stringS in ALL_S:
            print (f'==== {stringS} ({onsetName})')

            # name for this batch
            batchName = f'{onsetName}_{stringS}'

            # also basename for output
            thisSimOutBase = pathlib.Path (SIMULATION_DIR, batchName)

            # set up selection coefficients
            floatS = float(stringS)
            (s1, s2) = (float(FIXED_H) * floatS, floatS)
            # onset or not
            selectionChangeTimes = None
            if (not onsetFlag):
                stringS1 = f"{s1:.6f}"
                stringS2 = f"{s2:.6f}"
            else:
                selectionChangeTimes = [SELECTION_ONSET]
                # and modify coefficients
                selInd = numpy.array([0,1])
                s1 = list(s1*selInd)
                s2 = list(s2*selInd)
                stringS1 = ",".join([f"{x:.6f}" for x in s1])
                stringS2 = ",".join([f"{x:.6f}" for x in s2])

            # get some sampling times and sizes here, because they depend on key
            thisSamplingTimes = numpy.linspace(0, (NUM_SAMPLING_TIMES-1)*GEN_BETWEEN_SAMPLES, NUM_SAMPLING_TIMES)
            thisSampleSizes = numpy.array ([NUM_SAMPLES_EACH]*len(thisSamplingTimes), dtype=int)

            # put together the command line
            cmdList = [
                f"{DIPLO_SIM_EXEC}",
                f"--u01 {DIFF_M_ALPHA} --u10 {DIFF_M_BETA}",
                f"--Ne {DIFF_NE}",
                f"--num_rep {NUM_REPLICATES}",
                f"--seed {thisRNG.integers (99999999)}",
                f"--init initFreq --initFreq {INIT_FREQ}",
                f"--sample_times {','.join([str(int(x)) for x in thisSamplingTimes])}",
                f"--sample_sizes {','.join([str(int(x)) for x in thisSampleSizes])}",
                "--not_lost",
                f"--minMAF {MAF_FILTER}",
                f"--deltaT {DELTA_T}",
                f"-o {thisSimOutBase}",
                f"--s1 {stringS1} --s2 {stringS2}",
            ]
            if (onsetFlag):
                cmdList.append (f"--selection_change_times {','.join([str(int(x)) for x in selectionChangeTimes])}")

            # and simulate the batch
            print ("[DIPLO_LOCUS]")
            simCmd = " ".join(cmdList)
            THE_EXECUTOR.runCmd (simCmd)
            print ("[DIPLO_LOCUS_END]")



def load_all_simulations ():

    assert (SIMULATION_DIR.is_dir()), str(SIMULATION_DIR)

    simulatedSampleSizes = None
    simulatedSamplingTimes = {}
    simulatedSamples = {}

    # onset or constant
    for onsetFlag in [True, False]:
        onsetName = "onset" if onsetFlag else "constant"

        for stringS in ALL_S:
            # name for this batch
            batchName = f'{onsetName}_{stringS}'
            # get the input file
            candFiles = [str(x) for x in SIMULATION_DIR.glob (f"{batchName}_*_samples.count")]
            # should only be one
            assert (len(candFiles) == 1), len(candFiles)
            theFile = candFiles[0]
            print (theFile)

            # get the sampling times out of the file
            thisSamplingTimes = None
            ifs = open (theFile, 'r')
            for line in ifs.readlines ():
                line = line.strip()
                if (line.startswith ('## Sampling times:')):
                    thisSamplingTimes = numpy.array ([float(x) for x in line.split(':')[-1].split(';')[0].split(',')])
                    break
            ifs.close()
            assert (thisSamplingTimes is not None)

            # what's the time?
            if (onsetName not in simulatedSamplingTimes):
                simulatedSamplingTimes[onsetName] = thisSamplingTimes
            else:
                testSamplingTimes = simulatedSamplingTimes[onsetName]
                assert (numpy.all (testSamplingTimes == thisSamplingTimes)), (testSamplingTimes, thisSamplingTimes)

            # get the actual data
            simFrame = pandas.read_csv (theFile, sep='\t', comment='#')

            # and store/compare everything properly
            simulatedSamples[batchName] = simFrame
            # print (simFrame.shape)
            # print (simFrame.columns)
            allSizes = simFrame.iloc[:,2::2]
            # make sure all the sizes in a column match
            assert (numpy.all (numpy.diff (allSizes, axis=0) == 0))
            thisSampleSizes = numpy.array (allSizes.iloc[0,:])
            if (simulatedSampleSizes is None):
                simulatedSampleSizes = thisSampleSizes
            else:
                assert (numpy.all (simulatedSampleSizes == thisSampleSizes)), (simulatedSampleSizes, thisSampleSizes)

    return (simulatedSamples, simulatedSamplingTimes, simulatedSampleSizes)


def write_data_diplo_locus (combinedData, simulatedSamplingTimes):

    # diplo-locus & lls -- store each scenario in one file
    # RESULTS_DIR.mkdir()
    RESULTS_DIR.mkdir (exist_ok=True)

    # format: ID    d1  n1  d2  n2  ... (tab-separated)
    thisCombinedFile = pathlib.Path (RESULTS_DIR, f"onset_samples.count")

    ofs = open (thisCombinedFile, 'w')
    ofs.write (f"## Sampling times: {', '.join([str(int(x)) for x in simulatedSamplingTimes])}" + ";\n")
    combinedData.to_csv (ofs, sep='\t', index=False, header=True)
    ofs.close()


def main():

    # simulate some data
    simulate_all_data (metaRNG)


    # load the data, so we can write it in different formats
    (simulatedSamples, simulatedSamplingTimes, simulatedSampleSizes) = load_all_simulations ()

    # make sure all good
    assert (len(simulatedSamplingTimes['onset']) == len(simulatedSamplingTimes['constant']))
    assert (numpy.all(simulatedSamplingTimes['onset'] == simulatedSamplingTimes['constant']))
    # unify
    simulatedSamplingTimes = simulatedSamplingTimes['onset']

    # for (k, v) in simulatedSamples.items():
    #     print (k, v.shape)
    # for (k, v) in simulatedSamplingTimes.items():
    #     print (k, v)
    # print (simulatedSampleSizes)


    ####################
    # now store the data in the correct input format for the respective methods
    ####################

    # putting the data into blocks is probably good regardless of method
    combinedData = None
    for (k, v) in simulatedSamples.items():
        allSplit = k.split('_')
        onsetName = allSplit[0]
        thisCoefficient = allSplit[1]
        print (onsetName, thisCoefficient)

        # add the onsetName and selection coefficient in the index column though
        v['ID'] = [f'{x}_{onsetName}_{thisCoefficient}' for x in v['ID']]

        if (combinedData is None):
            combinedData = v
        else:
            combinedData = pandas.concat ([combinedData, v], ignore_index=True)
            

    # shuffle rows, so when methods only get part-way, they still have estimates for all parameters
    # times 2 for onset or not
    numCases = len(ALL_S) * 2
    blockStarts = numpy.tile(numpy.arange(0,numCases*NUM_REPLICATES,NUM_REPLICATES),NUM_REPLICATES)
    idxWithinBlock = numpy.repeat (numpy.arange(NUM_REPLICATES), numCases)
    shuffledIdxs = blockStarts + idxWithinBlock
    combinedData = combinedData.reindex(shuffledIdxs).reset_index(drop=True)

    # print (combinedData.shape)
    # print (combinedData)


    # write it
    write_data_diplo_locus (combinedData, simulatedSamplingTimes)


if __name__ == "__main__":
    main()

