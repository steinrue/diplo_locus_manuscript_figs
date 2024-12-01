import numpy
import pandas
import pathlib
# custom modules
import execution_engine

metaRNG = numpy.random.default_rng (4714)




# how to execute
THE_EXECUTOR = execution_engine.FileExecutor ("SLURM_SIMULATE.txt", append=False)
# THE_EXECUTOR = execution_engine.SubprocessExecutor()
# where execultable
DIPLO_SIM_EXEC = "DiploLocus-simulate"
SIMULATION_DIR = pathlib.Path ("simulated_data")
# put the simulated data in the respective directory
# some need it, so just do it for all
RESULTS_DIR = {
    'diplo_locus' : pathlib.Path ("results_diplo_locus"),
    'lls' : pathlib.Path ("results_lls"),
    'wfabc' : pathlib.Path ("results_wfabc"),
    'bmws' : pathlib.Path ("results_bmws"),
}
# fixed parameters for simulations
# sampling scheme
GEN_BETWEEN_SAMPLES = {
    'g4000' : 500,
    'g160' : 20,
}
NUM_SAMPLING_TIMES = 9
NUM_SAMPLES_EACH = 40
# diffusion
DIFF_M_ALPHA = 0
DIFF_M_BETA = DIFF_M_ALPHA
DIFF_NE = {
    'g4000' : 10000,
    'g160' : 1000,
}
# selection
ALL_S = {
    'g4000' : [f'{x:.4f}' for x in [0, 0.001, 0.002, 0.003, 0.004, 0.005]],
    # 'g4000' : [f'{x:.4f}' for x in [0, 0.002, 0.005]],
    'g160' : [f'{x:.4f}' for x in [0, 0.01, 0.02, 0.03, 0.04, 0.05]],
    # 'g160' : [f'{x:.4f}' for x in [0, 0.02, 0.05]],
}
FIXED_H = ["0.5", "1.0"]
SELECTION_CHANGE_TIMES = {
    'g4000' : [1000,3000],
    'g160' : [40,120],
}
# simulation
DELTA_T = 0.1
# DELTA_T = 1
# NUM_REPLICATES = lambda s : 500 if numpy.isclose(s, 0) else 100
# NUM_REPLICATES = lambda s : 2000
NUM_REPLICATES = 2000
MAF_FILTER = 0.05




def simulate_all_data (scenarios, thisRNG):

    # create directory for all the simulations
    # SIMULATION_DIR.mkdir()
    SIMULATION_DIR.mkdir (exist_ok=True)


    # prepare file with initial distribution for standing variation
    # all popsizes needed
    for thisNe in DIFF_NE.values():
        numFocalAlleles = numpy.arange(1,2*thisNe)
        freqs = numpy.arange(1,2*thisNe) / (2*thisNe)
        weights = 1/numFocalAlleles
        weights = weights / numpy.sum(weights)
        stdVarFrame = pandas.DataFrame()
        stdVarFrame['freqs'] = freqs
        stdVarFrame['weights'] = weights
        stdVarDistFile = pathlib.Path (SIMULATION_DIR, f'stdVar_{int(thisNe)}.dist')
        stdVarFrame.to_csv (stdVarDistFile, sep='\t', header=False, index=False)


    # go through all scenarios
    for (thisScenario, thisParameters) in scenarios.items():

        print (f'++ {thisScenario}')
        totalGensKey = thisScenario.split('_')[0].strip()

        # go though all selection
        for stringS in ALL_S[totalGensKey]:
            print (f'==== {stringS}')

            # name for this batch
            batchName = f'{thisScenario}_{stringS}'

            # also basename for output
            thisSimOutBase = pathlib.Path (SIMULATION_DIR, batchName)

            # Ne always important
            thisNe = DIFF_NE[totalGensKey]
            stdVarDistFile = pathlib.Path (SIMULATION_DIR, f'stdVar_{int(thisNe)}.dist')

            # set up selection coefficients
            floatS = float(stringS)
            thisH = thisParameters['fixedH']
            (s1, s2) = (thisH * floatS, floatS)
            # changing or not
            selectionChangeTimes = None
            if ('selChangeTimes' not in thisParameters):
                stringS1 = f"{s1:.6f}"
                stringS2 = f"{s2:.6f}"
            else:
                selectionChangeTimes = thisParameters['selChangeTimes']
                # and modify coefficients
                selInd = numpy.array(thisParameters['selIndicator'])
                s1 = list(s1*selInd)
                s2 = list(s2*selInd)
                stringS1 = ",".join([f"{x:.6f}" for x in s1])
                stringS2 = ",".join([f"{x:.6f}" for x in s2])

            # get some sampling times and sizes here, because they depend on key
            thisBetweenSamples = GEN_BETWEEN_SAMPLES[totalGensKey]
            thisSamplingTimes = numpy.linspace(0, (NUM_SAMPLING_TIMES-1)*thisBetweenSamples, NUM_SAMPLING_TIMES)
            thisSampleSizes = numpy.array ([NUM_SAMPLES_EACH]*len(thisSamplingTimes), dtype=int)

            # put together the command line
            # DiploLocus-simulate --u01 0 --u10 0 --Ne 10000 --num_rep 100 --seed 4711 --s1 0.025 --s2 0.05 --init custom
            # --initDistn simulated_data_diplo/stdVar_10000.dist --sample_times 0,10,42,68,95,100 --sample_sizes 6,10,12,24,16,40
            # --not_lost --minMAF 0.05 -o diplo_sim_test --write_traj --deltaT 0.5
            cmdList = [
                f"{DIPLO_SIM_EXEC}",
                f"--u01 {DIFF_M_ALPHA} --u10 {DIFF_M_BETA}",
                f"--Ne {thisNe}",
                # f"--num_rep {NUM_REPLICATES(floatS)}",
                f"--num_rep {NUM_REPLICATES}",
                f"--seed {thisRNG.integers (99999999)}",
                f"--init custom --initDistn {str(stdVarDistFile)}",
                f"--sample_times {','.join([str(int(x)) for x in thisSamplingTimes])}",
                f"--sample_sizes {','.join([str(int(x)) for x in thisSampleSizes])}",
                "--not_lost",
                f"--minMAF {MAF_FILTER}",
                f"--deltaT {DELTA_T}",
                f"-o {thisSimOutBase}",
                # "--write_traj",
                f"--s1 {stringS1} --s2 {stringS2}",
            ]
            if (selectionChangeTimes is not None):
                cmdList.append (f"--selection_change_times {','.join([str(int(x)) for x in selectionChangeTimes])}")

            # and simulate the batch
            print ("[DIPLO_LOCUS]")
            simCmd = " ".join(cmdList)
            THE_EXECUTOR.runCmd (simCmd)
            print ("[DIPLO_LOCUS_END]")


def simulationScenarios ():

    # make sure stuff checks out
    assert (GEN_BETWEEN_SAMPLES.keys() == DIFF_NE.keys())
    assert (GEN_BETWEEN_SAMPLES.keys() == ALL_S.keys())
    assert (GEN_BETWEEN_SAMPLES.keys() == SELECTION_CHANGE_TIMES.keys())

    # put together scenarios
    scenarios = {}
    for totalGens in GEN_BETWEEN_SAMPLES.keys():
        for changing in [False, True]:
            for thisH in FIXED_H:
                thisParameters = {}
                thisParameters['fixedH'] = float(thisH)
                if (changing):
                    assert (len(SELECTION_CHANGE_TIMES[totalGens]) == 2), len(SELECTION_CHANGE_TIMES[totalGens])
                    thisParameters['selChangeTimes'] = SELECTION_CHANGE_TIMES[totalGens]
                    thisParameters['selIndicator'] = [0,1,0]
                thisName = f"{totalGens}_h{thisH}_{'var' if changing else 'const'}"
                scenarios[thisName] = thisParameters

    return scenarios


def load_all_simulations (scenarios):

    assert (SIMULATION_DIR.is_dir()), str(SIMULATION_DIR)

    simulatedSampleSizes = None
    simulatedSamplingTimes = {}
    simulatedSamples = {}

    for (thisScenario, thisParameters) in scenarios.items():

        # need the key again
        totalGensKey = thisScenario.split('_')[0].strip()

        for stringS in ALL_S[totalGensKey]:
            # name for this batch
            batchName = f'{thisScenario}_{stringS}'
            # get the input file
            candFiles = [str(x) for x in SIMULATION_DIR.glob (f"{batchName}_*_samples.count")]
            # should only be one
            assert (len(candFiles) == 1), len(candFiles)
            theFile = candFiles[0]
            # print (theFile)

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
            if (totalGensKey not in simulatedSamplingTimes):
                simulatedSamplingTimes[totalGensKey] = thisSamplingTimes
            else:
                testSamplingTimes = simulatedSamplingTimes[totalGensKey]
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


def write_data_diplo_locus_lls (combinedData, simulatedSamplingTimes):

    # diplo-locus & lls -- store each scenario in one file
    # RESULTS_DIR['diplo_locus'].mkdir()
    RESULTS_DIR['diplo_locus'].mkdir (exist_ok=True)
    # RESULTS_DIR['lls'].mkdir()
    RESULTS_DIR['lls'].mkdir (exist_ok=True)

    # format: ID    d1  n1  d2  n2  ... (tab-separated)
    for k in combinedData.keys():
        print (f"=== {k}")
        print (combinedData[k].shape)
        totalGensKey = k.split('_')[0].strip()

        for thisDir in [RESULTS_DIR['diplo_locus'], RESULTS_DIR['lls']]:
            # we already have a nice data-frame, so just the header and then the rest
            thisCombinedFile = pathlib.Path (thisDir, f"{k}_samples.count")
            ofs = open (thisCombinedFile, 'w')
            ofs.write (f"## Sampling times: {', '.join([str(int(x)) for x in simulatedSamplingTimes[totalGensKey]])}" + ";\n")
            combinedData[k].to_csv (ofs, sep='\t', index=False, header=True)
            ofs.close()


def write_data_wfabc (combinedData, simulatedSamplingTimes):

    # wfabc has its own format -- store each scenario in one file as well though
    # RESULTS_DIR['wfabc'].mkdir()
    RESULTS_DIR['wfabc'].mkdir (exist_ok=True)
            
    # format:
    # first line: two numbers: nr_loci num_time_points
    # second line: sampling times in generations
    # then alternating pair of lines for each locus (comma seperated)
    # one line sample sizes
    # one line sampled alleles
    for k in combinedData.keys():
        print (f"=== {k}")
        totalGensKey = k.split('_')[0].strip()
        thisFrame = combinedData[k]
        print (thisFrame.shape)

        # how many loci?
        numLoci = thisFrame.shape[0]
        numTimes = len(simulatedSamplingTimes[totalGensKey])

        # we need to reshape the frame a bit for the right output here
        wfabcMatrix = numpy.zeros ((2*numLoci, numTimes), dtype=int)
        # where do we find numDers and sizes? 
        numDerCols = [f'd{a}' for a in numpy.arange(1,numTimes+1)]
        sizeCols = [f'n{a}' for a in numpy.arange(1,numTimes+1)]
        # put them in the right places
        wfabcMatrix[0::2] = thisFrame[sizeCols]
        wfabcMatrix[1::2] = thisFrame[numDerCols]

        # write the file
        thisCombinedFile = pathlib.Path (RESULTS_DIR['wfabc'], f"{k}_loci.txt")
        ofs = open (thisCombinedFile, 'w')
        ofs.write (f"{numLoci} {numTimes}\n")
        ofs.write (f"{','.join([str(int(x)) for x in simulatedSamplingTimes[totalGensKey]])},\n")
        # print (wfabcMatrix)
        numpy.savetxt (ofs, wfabcMatrix.astype(int), delimiter=',', fmt='%d', newline=',\n')
        # combinedData[k].to_csv (ofs, sep='\t', index=False, header=True)
        ofs.close()

        # I also need the ids to remeber which one was simulated with which selection coefficient
        thisIdFile = pathlib.Path (RESULTS_DIR['wfabc'], f"{k}.ids")
        numpy.savetxt (thisIdFile, numpy.array(thisFrame['ID']).reshape(thisFrame.shape[0],1),fmt="%s")


def write_data_bmws (combinedData, simulatedSamplingTimes):

    # bmws -- store each scenario in one file
    # RESULTS_DIR['bmws'].mkdir()
    RESULTS_DIR['bmws'].mkdir (exist_ok=True)

    # format: vcf, but pseudohaploid, where 0/0 is one ancestral and 1/1 is one derived allele
    # also a second file that has details on when samples are taken
    for k in combinedData.keys():
        print (f"=== {k}")
        totalGensKey = k.split('_')[0].strip()
        print (combinedData[k].shape)
        thisFrame = combinedData[k]

        # put all the relevant things in the right places
        numLoci = thisFrame.shape[0]
        numTimes = len(simulatedSamplingTimes[totalGensKey])

        # get the sizes
        sizeCols = [f'n{a}' for a in numpy.arange(1,numTimes+1)]
        allSizes = numpy.array(thisFrame[sizeCols])
        assert (numpy.all (numpy.diff (allSizes, axis=0) == 0))
        sampleSizes = allSizes[0,:]

        # and the number of derived alleles
        numDerCols = [f'd{a}' for a in numpy.arange(1,numTimes+1)]
        numDerAlleles = numpy.array(thisFrame[numDerCols])

        # put together individuals and genotypes
        fakeIds = []
        individualSamplingTimes = []
        genotypes = []
        for (gIdx, thisGen) in enumerate([int(x) for x in simulatedSamplingTimes[totalGensKey]]):
            thisNumDer = numDerAlleles[:,gIdx]
            for iIdx in numpy.arange(sampleSizes[gIdx]):
                # get a fake id
                fakeIds.append (f"gen{thisGen}_ind{iIdx}")
                # get a real sampling time
                individualSamplingTimes.append (thisGen)

                # get some genotype for this ind (LD shouldn't matter)
                indGeno = numpy.array(['0/0']*numDerAlleles.shape[0])
                # num derived alleles decides how many ones to get at corresponding "locus"
                indGeno[thisNumDer > iIdx] = '1/1'

                # and remember it
                genotypes.append (indGeno)

        # now we have everything in place
        # fake an id file (do all columns, because who knows what we need)
        idFrame = pandas.DataFrame({
            '#ID' : fakeIds,
            'DateBP' : individualSamplingTimes,
            'Region' : ['Home']*len(fakeIds),
            'Latitude' : ['25.00']*len(fakeIds),
            'Longitude' : ['-70.00']*len(fakeIds),
        })
        idFile = pathlib.Path (RESULTS_DIR['bmws'], f"{k}.meta")
        idFrame.to_csv (idFile, index=False, sep='\t')

        # and fake a vcf file
        # info columns
        vcfDict = {
            '#CHROM' : 1,
            'POS' : numpy.arange(1,numLoci+1),
            'ID' : thisFrame['ID'],
            'REF' : ['A']*numLoci,
            'ALT' : ['T']*numLoci,
            'QUAL' : ['100']*numLoci,
            'FILTER' : ['PASS']*numLoci,
            'INFO' : ['.']*numLoci,
            'FORMAT' : ['GT']*numLoci
        }
        # individuals
        for (gIdx, thisGeno) in enumerate(genotypes):
            vcfDict[fakeIds[gIdx]] = thisGeno
        # and write the file
        vcfFile = pathlib.Path (RESULTS_DIR['bmws'], f"{k}.vcf")
        # header
        ofs = open (vcfFile, 'w')
        ofs.write ('##fileformat=VCFv4.0\n##source=faked\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        # body
        vcfFrame = pandas.DataFrame(vcfDict)
        # print (numDerAlleles[7,:])
        # print (numpy.array(vcfFrame.iloc[7,:]))
        vcfFrame.to_csv (ofs, index=False, sep='\t')
        ofs.close()


def main():

    # get the scenarios to simulate
    theScenarios = simulationScenarios()


    # simulate some data
    simulate_all_data (theScenarios, metaRNG)


    # load the data, so we can write it in different formats
    (simulatedSamples, simulatedSamplingTimes, simulatedSampleSizes) = load_all_simulations (theScenarios)

    # for (k, v) in simulatedSamples.items():
    #     print (k, v.shape)
    # for (k, v) in simulatedSamplingTimes.items():
    #     print (k, v)
    # print (simulatedSampleSizes)


    ####################
    # now store the data in the correct input format for the respective methods
    ####################

    # putting the data into blocks is probably good regardless of method
    combinedData = {}
    for (k, v) in simulatedSamples.items():
        allSplit = k.split('_')
        thisScenario = '_'.join(allSplit[:3])
        thisCoefficient = allSplit[3]
        print (thisScenario, thisCoefficient)

        # add the selection coefficient in the index column though
        v['ID'] = [f'{x}_{thisCoefficient}' for x in v['ID']]

        if (thisScenario not in combinedData):
            combinedData[thisScenario] = v
        else:
            combinedData[thisScenario] = pandas.concat ([combinedData[thisScenario], v], ignore_index=True)
            

    # shuffle rows, so when methods only get part-way, they still have estimates for all parameters
    shuffledData = {}
    for (thisScenario, thisFrame) in combinedData.items():
        thisGenKey = thisScenario.strip().split('_')[0]
        numS = len(ALL_S[thisGenKey])
        blockStarts = numpy.tile(numpy.arange(0,numS*NUM_REPLICATES,NUM_REPLICATES),NUM_REPLICATES)
        idxWithinBlock = numpy.repeat (numpy.arange(NUM_REPLICATES), numS)
        shuffledIdxs = blockStarts + idxWithinBlock
        shuffledData[thisScenario] = thisFrame.reindex(shuffledIdxs).reset_index(drop=True)
    combinedData = shuffledData

    # for (k,v) in combinedData.items():
    #     print (k, v.shape)
    #     print (v)


    write_data_diplo_locus_lls (combinedData, simulatedSamplingTimes)

    write_data_wfabc (combinedData, simulatedSamplingTimes)
            
    write_data_bmws (combinedData, simulatedSamplingTimes)


if __name__ == "__main__":
    main()

