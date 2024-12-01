import numpy
import pickle
import bz2
import time
import sys
import diplo_locus.likelihood
import scipy
import joblib
import pandas




# which replicates do we have?
SYNC_NUM_REPLICATES = 10
SYNC_NUM_TIMEPOINTS = 7
# filtering
MAF_FILTER = 0.05
# parameter for diplo-locus
M_ALPHA = 5e-9
M_BETA = M_ALPHA
NE = 300
GENIC_DIRECTION = numpy.array([0.5, 1])
PRE_COEFFS = numpy.geomspace (0.0075, 0.75, 31)
SEL_COEFFS = numpy.concatenate (( - numpy.flip(PRE_COEFFS), [0], PRE_COEFFS))
INIT_DICT = {'initCond' : 'uniform'}
NUM_STATES = 501
# some statistics & numerics
MIN_PVALUE = 1e-32
EPSILON = 1e-12




def loadReplicate (inputFile, replicateToAnalyze):

    # get the pickle out of the jar
    pickelData = pickle.load (bz2.open (inputFile, 'rb'))
    # print (pickelData.keys()) 
    thePositions = pickelData['thePositions']
    majorReadCounts = pickelData['majorReadCounts']
    minorReadCounts = pickelData['minorReadCounts']
    assert (majorReadCounts.shape == minorReadCounts.shape)
    assert (len(thePositions) == majorReadCounts.shape[0])
    assert (SYNC_NUM_REPLICATES*SYNC_NUM_TIMEPOINTS == majorReadCounts.shape[1])
    # print (len(thePositions), majorReadCounts.shape, minorReadCounts.shape)

    # get the replicate we want, and convert it to proper time-series
    # columns in original sync-file are (according to README_for_F0-F60SNP_CMH_FET_blockID.sync.docx):
    # Chromosome	position	base	Dsim_Fl_Base_1	...	Dsim_Fl_Base_10	Dsim_Fl_Hot_F10_1	...	Dsim_Fl_Hot_F10_10	...	Dsim_Fl_Hot_F60_1	...	Dsim_Fl_Hot_F60_10 \
    #               -log10(pvalue)_CMH  -log10(p-value)_FET_rep1    ...	-log10(p-value)_FET_rep10	blockID_0.75cor	blockID_0.35cor
    # also somewhat compatible with https://code.google.com/archive/p/popoolation2/wikis/Tutorial.wiki -- though the p-values are not allowed in a real *.sync file
    assert ((replicateToAnalyze >= 0) and (replicateToAnalyze < SYNC_NUM_REPLICATES)), f"Replicate with given index {replicateToAnalyze} does not exist. Only {SYNC_NUM_REPLICATES} in file."
    specifiedMajorCounts = majorReadCounts[:,replicateToAnalyze::SYNC_NUM_REPLICATES]
    specifiedMinorCounts = minorReadCounts[:,replicateToAnalyze::SYNC_NUM_REPLICATES]
    assert (specifiedMajorCounts.min() >= 0)
    assert (specifiedMinorCounts.min() >= 0)
    assert (specifiedMajorCounts.shape == specifiedMinorCounts.shape)
    assert (specifiedMajorCounts.shape[1] == SYNC_NUM_TIMEPOINTS)

    # make proper time-series
    sampleSizes = specifiedMajorCounts + specifiedMinorCounts
    samples = specifiedMinorCounts
    samplingTimes = numpy.arange(0,SYNC_NUM_TIMEPOINTS*10,10)
    totalSampleSizes = sampleSizes.sum (axis=1)
    print (f"min total size: {totalSampleSizes.min()}")
    assert (totalSampleSizes.min() > 0)
    print (f"max total size: {totalSampleSizes.max()}")
    sampleTimePoints = (sampleSizes > 0).sum (axis=1)
    print (f"min sampling times: {sampleTimePoints.min()}")
    assert (sampleTimePoints.min() >= 2)

    # apply MAF filter here
    assert (numpy.all(totalSampleSizes > 0))
    fs = samples.sum(axis=1)/totalSampleSizes
    mafs = numpy.minimum (fs, 1-fs)
    mafMask = (mafs > MAF_FILTER)
    print (f"pass maf-filter: {mafMask.sum()} out of {len(mafMask)}")
    thePositions = thePositions[mafMask]
    sampleSizes = sampleSizes[mafMask,:]
    samples = samples[mafMask,:]

    return (thePositions, samplingTimes, sampleSizes, samples)
    
    
# helpful function
def findMax (params, ll, neutralLL):

    # need to replace some potential -infs, before interpolating
    negInfMask = numpy.isneginf (ll)
    theMin = ll[~negInfMask].min()
    assert (theMin < 0)
    ll[negInfMask] = 10 * theMin
    assert (not numpy.any(numpy.isinf (ll))), ll

    # normalize it (put the 2 already here)
    normLL = 2 * (ll - neutralLL)

    # make an interp function
    f = scipy.interpolate.interp1d (params, normLL, kind="cubic")

    # the right bounds
    bounds = (min(params),max(params))

    # and find the optimum (minus, cause we can only minimize)
    theOpt = scipy.optimize.minimize_scalar (lambda x : - f(x), method='bounded', bounds=(bounds[0], bounds[1]))

    # return the results
    maxParam = float(theOpt.x)
    maxLL = f(maxParam)

    # if this is negative, we did not find a good off-grid maximum
    # in this case take on grid
    offGridInd = True
    if (maxLL < 0):
        # take on grid maximum
        offGridInd = False
        onGridMaxIdx = numpy.argmax (normLL)
        print (f'[DANGER] ({maxParam:.5f}, {maxLL:.5f}) vs. ({params[onGridMaxIdx]:.5f}, {normLL[onGridMaxIdx]:.5f})')
        maxParam = params[onGridMaxIdx]
        maxLL = normLL[onGridMaxIdx]
    else:
        # take off-grid, as is
        pass

    return (numpy.array ([maxParam, maxLL]), offGridInd)


def surfaceOneDParallel (numCpus, times, samplesSizes, samples, selDirection, selCoefficients, **kwargs):

    # some argument stuff
    assert ("emissionType" in kwargs)
    if ((kwargs["emissionType"] == "integer") and ("sampleSizesSet" not in kwargs)):
        assert (issubclass (numpy.array(samplesSizes).dtype.type, numpy.integer))
        assert (issubclass (numpy.array(samples).dtype.type, numpy.integer))
        kwargs["sampleSizesSet"] = set(samplesSizes.flatten())

    # now build a proper grid
    assert (len(selDirection.shape) == 1)
    assert (selDirection.shape[0] == 2)
    assert (len(selCoefficients.shape) == 1)
    selGrid = numpy.zeros ((len(selCoefficients),2))
    selGrid[:,0] = selCoefficients * selDirection[0]
    selGrid[:,1] = selCoefficients * selDirection[1]

    # and iterate over the coefficients
    assert (samples.shape == samplesSizes.shape)
    thisLL = numpy.zeros ((samples.shape[0], selGrid.shape[0]))

    # parallelized function
    def singleCoefficent (selIdx):
        # get them coefficents
        s1 = selGrid[selIdx,0]
        s2 = selGrid[selIdx,1]
        print (s1,s2)
        # set up hmm
        # hopefully someone figured out the initial conditions in initialKeywords
        assert ("initCond" in kwargs)
        # selection coefficients have to be given
        # just give all kwargs to the constructor and let it figure out whether the right things are there
        thisHMM = diplo_locus.likelihood.SelHmm (s1=s1, s2=s2, **kwargs)
        # compute the likelihood for this coefficient
        return  thisHMM.computeLogLikelihood (times, samplesSizes, samples)

    # iterate over selection (in parallel)
    thisLL = joblib.Parallel(n_jobs=int(numCpus))(joblib.delayed(singleCoefficent)(selIdx) for selIdx in range(selGrid.shape[0]))
    # need to convert it from list to numpy.array for later
    # it supposedly preserves order, which we absolutely need
    thisLL = numpy.array(thisLL).transpose()

    # return it
    return (thisLL, selGrid)


def getLikelihoodRatioParallel (numCpus, thisLL, selGrid):

    # how many times is grid-MLE on boundary?
    maxIdxs = thisLL.argmax(axis=1)
    # which ones are on the boundaries, cause for these we can just set the maxLL
    boundMask = (maxIdxs <= 0) | (maxIdxs >= (len(selGrid[:,1]) - 1))
    print (f"num mles on boundary: {boundMask.sum()}")

    # find max for each surface
    # but now in parallel
    neutralIdx = numpy.where (numpy.isclose (selGrid[:,1], 0))[0][0]

    # parallelized function
    def findMaxForIdx (idx):

        if (idx % 5000 == 0):
            print (f"it {idx}")

        # find max for this one
        if (boundMask[idx]):
            # is one the bound, so just take the value from the grid
            thisMaxIdx = thisLL[idx,:].argmax()
            thisLLR = 2*(thisLL[idx,thisMaxIdx] - thisLL[neutralIdx,thisMaxIdx])
            thisMLE = selGrid[thisMaxIdx,1]
        else:
            # not on the bound, so the interpolation stuff could work
            (mleAndLL, offGridInd) = findMax (selGrid[:,1], thisLL[idx,:].copy(), thisLL[idx,neutralIdx])
            thisLLR = mleAndLL[1]
            thisMLE = mleAndLL[0]

        return (thisMLE, thisLLR)

    # and run it
    thisStats = joblib.Parallel(n_jobs=int(numCpus))(joblib.delayed(findMaxForIdx)(idx) for idx in range (thisLL.shape[0]))
    thisStats = numpy.array(thisStats)

    # might as well turn it into a data frame already
    statsFrame = pandas.DataFrame({
        'mle' : thisStats[:,0],
        'llr' : thisStats[:,1],
    })

    # make sure we values everywhere that look ok
    assert (statsFrame['mle'].max() <= (selGrid[:,1].max() + EPSILON)), 'missing mle'
    assert (statsFrame['llr'].min() > -EPSILON), 'LLR too negative'

    # I guess we keep this vectorized, but don't run it fully parallel (yet)
    # get some p-values from a chi2
    rawPValues = 1 - scipy.stats.chi2.cdf (thisStats[:,1], 1)
    # just so we have stuff working downstream
    assert (rawPValues.max() <= 1.0)
    statsFrame['pValue'] = numpy.maximum (rawPValues, MIN_PVALUE)

    return statsFrame


def analyzeDiploLocus (numCpus, blockSize, samplingTimes, sampleSizes, samples):

    assert (sampleSizes.shape[1] == len(samplingTimes))
    assert (sampleSizes.shape == samples.shape)

    # compute likelihood on a genic surface using diplo-locus
    # need some modded sampling times
    modSamplingTimes = numpy.concatenate (([0], samplingTimes))

    # do it in blocks though
    blockStartIdx = 0
    blockEndIdx = blockSize
    statsFrame = None
    # until no more data there
    while (blockStartIdx < samples.shape[0]):

        print (f"Block: {blockStartIdx}:{blockEndIdx} out of {samples.shape[0]} total")

        # compute surface for this block
        startTime = time.time()
        blockedSampleSizes = sampleSizes[blockStartIdx:blockEndIdx,:]
        blockedSamples = samples[blockStartIdx:blockEndIdx,:]
        (thisLL, selGrid) = surfaceOneDParallel (numCpus, modSamplingTimes, blockedSampleSizes, blockedSamples, GENIC_DIRECTION, SEL_COEFFS,
                                            **INIT_DICT, emissionType="integer", sampleSizesSet=set(blockedSampleSizes.flatten()),
                                            transitionType="constant", mAlpha=M_ALPHA, mBeta=M_BETA, Ne=NE, numStates=NUM_STATES)
        print (f"Surface compute time: {int(int(time.time() - startTime))}")

        # and then get some likelihood-ratios and p-values
        startTime = time.time()
        thisStatsFrame = getLikelihoodRatioParallel (numCpus, thisLL, selGrid)
        print (f"Likelihood-ratio compute time: {int(int(time.time() - startTime))}")

        # append the results
        if (statsFrame is None):
            statsFrame = thisStatsFrame
        else:
            statsFrame = pandas.concat ([statsFrame, thisStatsFrame], ignore_index=True)

        # step to next block
        blockStartIdx = blockEndIdx
        blockEndIdx = blockStartIdx + blockSize

    # return hopefully all the stats
    return statsFrame


def saveResults (thePositions, statsFrame, outputFile):

    # add positions to statsFrame
    assert (len(thePositions) == statsFrame.shape[0]), (len(thePositions), statsFrame.shape[0])
    statsFrame['positions'] = thePositions

    # reorder things
    assert (statsFrame.shape[1] == 4)
    statsFrame = statsFrame.iloc[:,[3,0,1,2]]

    # pickle in a jar
    pickleData = {
        'statsFrame' : statsFrame,
    }
    pickle.dump (pickleData, bz2.open (outputFile, 'wb'))


def wholeAnalysis (numCpus, blockSize, inputFile, replicateToAnalyze, outputFile):

    # first, load some data
    startTime = time.time()
    print ('[LOAD]')
    (thePositions, samplingTimes, sampleSizes, samples) = loadReplicate (inputFile, replicateToAnalyze)
    print (f'[LOAD_DONE] {int(time.time() - startTime)}')

    # and then analyze it
    print ('[ANALYZE_DIPLO_LOCUS]')
    statsFrame = analyzeDiploLocus (numCpus, blockSize, samplingTimes, sampleSizes, samples)
    print ('[ANALYZE_DIPLO_LOCUS_DONE]')

    # then save it
    startTime = time.time()
    print ('[SAVE_RESULTS]')
    saveResults (thePositions, statsFrame, outputFile)
    print (f'[SAVE_RESULTS_DONE] {int(time.time() - startTime)}')


def main():

    if (len(sys.argv) != 6):
        print ("usage: python <script_name> <num_cpus> <block_size> <input_file.pickle.bz2> <replicate_to_analyze> <output_file.pickle.bz2>")
        sys.exit (1)

    numCpus = int(sys.argv[1])
    blockSize = int(sys.argv[2])
    inputFile = sys.argv[3]
    replicateToAnalyze = int(sys.argv[4])
    outputFile = sys.argv[5]
    print (f"[ANALYZE] Repl. {replicateToAnalyze} from file {inputFile}.")
    print (f"[ANALYZE] Using {numCpus} cores with blocks of size {blockSize}.")
    print (f"[ANALYZE] Saving output to file {outputFile}.")

    # do it
    wholeAnalysis (numCpus, blockSize, inputFile, int(replicateToAnalyze), outputFile)


if __name__ == "__main__":
    main()

