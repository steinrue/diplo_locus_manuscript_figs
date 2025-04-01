import numpy
import pathlib
import pandas
import pickle
import bz2




# some meta parameters for simulation and analysis
MIN_ONSET_GEN = 0
MAX_ONSET_GEN = 4000
NUM_ONSET_GRID = 51
ONSETS = [int(x) for x in numpy.linspace (MIN_ONSET_GEN+1, MAX_ONSET_GEN-1, NUM_ONSET_GRID)]
ONSETS_WITH_CONST = ONSETS + ['nan']
ONSET_FILE_SKELETON = "results_diplo_locus/onset_%s_results_LLmatrices.table"
NUM_SIM_ONSET_OR_CONST = 2
NUM_SIM_SEL_COEFF = 6
NUM_SIM_SCENARIOS = NUM_SIM_ONSET_OR_CONST * NUM_SIM_SEL_COEFF
NUM_SIM_REPS = 500
# file name for storing
PICKLE_FILE = "final_lls.pickle.bz2"




def load_files():

    # get all the names of the input files
    llFiles = {}
    for thisOnset in ONSETS_WITH_CONST:
        # llFiles[thisOnset] = pathlib.Path (f"results_diplo_locus/onset_{thisOnset}_results_LLmatrices.table")
        # get the ones were we set selection in the first epoch to 0, unless there is no epoch in 'nan'
        if (thisOnset == 'nan'):
            llFiles[thisOnset] = pathlib.Path (ONSET_FILE_SKELETON % f"{thisOnset}")
        else:
            llFiles[thisOnset] = pathlib.Path (ONSET_FILE_SKELETON % f"{thisOnset}_0.000000e+00")

    # load the files
    llFrames = {}
    theHeader = None
    theIDs = None
    for (k,v) in llFiles.items():
        print (str(v))
        llFrames[k] = pandas.read_csv (str(v), sep='\t', comment='#')
        # headers and IDs are consistent scross files
        if (theHeader is None):
            theHeader = llFrames[k].columns
        else:
            assert (numpy.all (theHeader == llFrames[k].columns))
        if (theIDs is None):
            theIDs = llFrames[k]["ID"]
        else:
            assert (numpy.all (theIDs == llFrames[k]["ID"]))
        
    # return what we read
    return (llFrames, theHeader, theIDs)


def transformToSurf (llFrames, theHeader, theIDs):

    # get the selection coeffcients (s2)
    assert (theHeader[-1] == "ID")
    preCoeffs = theHeader[:-1]
    # whatever
    selCoeffs = numpy.array ([float(x.split(',')[1].split(')')[0]) for x in preCoeffs])

    # prepare splitting up IDs by scenario and simulated coefficient
    splitArray = numpy.array([x.split("_") for x in theIDs])
    assert (splitArray.shape == (NUM_SIM_SCENARIOS * NUM_SIM_REPS, 3))
    onsetScenarios = sorted(set(splitArray[:,1]))
    simCoeffs = sorted(set(splitArray[:,2]))
    assert (len(onsetScenarios)*len(simCoeffs) == NUM_SIM_SCENARIOS)
    repMasks = {}
    for thisScenario in onsetScenarios:
        for thisSimCoeff in simCoeffs:
            thisMask = (splitArray[:,1] == thisScenario) & (splitArray[:,2] == thisSimCoeff)
            assert (numpy.sum(thisMask) == NUM_SIM_REPS)
            repMasks[(thisScenario, thisSimCoeff)] = thisMask

    # and split up the results by scenario x simCoeffs
    llSurfs = {}
    llConstant = {}
    # go though the different scenario x simCoeffs
    for (thisSimScenario, thisMask) in repMasks.items():
        print (thisSimScenario)

        # iterate through all onsets
        # good thing they are in order =), cause we need to have them in order
        assert (ONSETS_WITH_CONST[-1] == 'nan')
        assert (numpy.all(numpy.diff(ONSETS_WITH_CONST[:-1]) > 0))
        for thisOnset in ONSETS_WITH_CONST:
            # print (thisOnset)
            assert (thisOnset in llFrames)

            thisFrame = llFrames[thisOnset]
            preLL = numpy.array (thisFrame.iloc[:,:-1])
            thisLL = preLL[thisMask,:]

            # get combined ll
            combinedLL = None
            # ll computation under onset or constant?
            # only onset computation needs possible extension
            if (thisOnset != 'nan'):
                # ll onset
                if (thisSimScenario in llSurfs):
                    combinedLL = llSurfs[thisSimScenario]

            # start or extend
            # unfortunately, this needs some newxis-shenanigans
            if (combinedLL is None):
                combinedLL = thisLL[:,:,numpy.newaxis]
            else:
                # this results in reps x sel x onset
                combinedLL = numpy.concatenate ([combinedLL,thisLL[:,:,numpy.newaxis]],axis=2)

            # store in the right place
            if (thisOnset == 'nan'):
                # should be the only one
                assert (thisSimScenario not in llConstant)
                # undo newxis-shenanigans
                llConstant[thisSimScenario] = combinedLL[:,:,0]
            else:
                llSurfs[thisSimScenario] = combinedLL

    # make sure the dimesnion of the surfaces look good
    for (k,v) in llConstant.items():
        assert (v.shape == (NUM_SIM_REPS, len(selCoeffs))), f"{k}: {v.shape} == {(NUM_SIM_REPS, len(selCoeffs))}"
    for (k,v) in llSurfs.items():
        assert (v.shape == (NUM_SIM_REPS, len(selCoeffs), len(ONSETS))), f"{k}: {v.shape} == {(NUM_SIM_REPS, len(selCoeffs), len(ONSETS))}"

    # gite it all away
    return (llSurfs, llConstant, NUM_SIM_REPS, selCoeffs, ONSETS)


def getMLEs(llSurfs, llConstant, numSimReps, selCoeffs, onsets):

    # get the neutralIdx for selection
    neutralIdx = numpy.argmin(numpy.abs(selCoeffs - 0))
    # try to find some MLEs
    # first for the constant guys
    print ("+++++ constant")
    constMLEs = {}
    for (k,v) in llConstant.items():
        # make sure things are shaped the right way
        print (k, v.shape)
        assert (v.shape == (numSimReps, len(selCoeffs)))

        # get the MLEs (only s direction)
        maxIdxs = numpy.argmax (v, axis=1)
        # get ll at max
        maxLL = v[numpy.arange(len(maxIdxs)),maxIdxs]
        neutralLL = v[numpy.arange(len(maxIdxs)),neutralIdx]
        # get value of mle
        sMLEs = selCoeffs[maxIdxs]

        # stack and save
        llAndMle = numpy.stack ([neutralLL, maxLL, sMLEs], axis=1)
        constMLEs[k] = llAndMle

    # then the onset guys
    print ("+++++ onset")
    onsetMLEs = {}
    for (k,v) in llSurfs.items():
        # make sure things are shaped the right way
        print (k, v.shape)
        assert (v.shape == (numSimReps, len(selCoeffs), len(onsets)))

        # get the MLEs (2d: rep x s x onset)
        # flatten v in the surface direction
        flatV = numpy.reshape (v, (numSimReps, len(selCoeffs) * len(onsets)), order='C')
        # get max-flat-idx for each durfaces
        maxIdxs = numpy.argmax (flatV, axis=1)
        # get the ll at max
        maxLL = flatV[numpy.arange(len(maxIdxs)),maxIdxs]
        # get the corresponding baseLL
        baseLL = llConstant[k][numpy.arange(len(maxIdxs)),neutralIdx]
        # unflatten the max-idxs
        unravelledMaxIdxs = numpy.unravel_index (maxIdxs, (len(selCoeffs), len(onsets)), order='C')
        # and get the corresponding values
        sMLEs = selCoeffs[unravelledMaxIdxs[0]]
        oMLEs = numpy.array(onsets)[unravelledMaxIdxs[1]]

        # stack em
        mleAndLL = numpy.stack ([baseLL, maxLL, sMLEs ,oMLEs], axis=1)

        # and remember
        onsetMLEs[k] = mleAndLL

    return (constMLEs, onsetMLEs)


def pickleSurfaces (llSurfs, llConstant, numSimReps, selCoeffs, onsets, constMLEs, onsetMLEs):

    # pickle
    pickleData = {
        "llSurfs" : llSurfs,
        "llConstant" : llConstant,
        "numSimReps" : numSimReps,
        "selCoeffs" : selCoeffs,
        "onsets" : onsets,
        "constMLEs" : constMLEs,
        "onsetMLEs" : onsetMLEs,
    }
    pickle.dump (pickleData, bz2.open (PICKLE_FILE, "wb"), )


def collect_results():

    print ("[LOAD_FILES]")
    (llFrames, theHeader, theIDs) = load_files()
    print ("[LOAD_FILES_DONE]")

    print ("[TRANSFORM_TO_SURFACES]")
    (llSurfs, llConstant, numSimReps, selCoeffs, onsets) = transformToSurf (llFrames, theHeader, theIDs)
    print ("[TRANSFORM_TO_SURFACES_DONE]")
    
    print ("[GET_MLES]")
    (constMLEs, onsetMLEs) = getMLEs(llSurfs, llConstant, numSimReps, selCoeffs, onsets)
    print ("[GET_MLES_DONE]")

    print ("[PICKLE_RESULTS]")
    pickleSurfaces (llSurfs, llConstant, numSimReps, selCoeffs, onsets, constMLEs, onsetMLEs)
    print ("[PICKLE_RESULTS_DONE]")


def main():

    collect_results()


if __name__ == "__main__":
    main()