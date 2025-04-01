import numpy
import pathlib
import pickle
import bz2
import scipy


# files with p-values for each replicate
NUM_REPLICATES = 10
# CHR_NAMES = [f"chr{x}" for x in ['2L', '2R', '3L', '3R', '4', 'X']]
CHR_NAMES = [f"chr{x}" for x in ['2L']]
INPUT_FILES = {}
INPUT_DIR = "results/"
for thisChr in CHR_NAMES:
    repFileList = []
    for idx in numpy.arange(NUM_REPLICATES):
        repFileList.append (pathlib.Path (INPUT_DIR, f"analyzed_{thisChr}_F0-F60_rep{idx}.pickle.bz2"))
    INPUT_FILES[thisChr] = repFileList
OUTPUT_FILE = pathlib.Path(INPUT_DIR, "analyzed_F0-F60_fisher.pickle.bz2")



def loadPerReplicateResults():

    # load the data from all the input file
    allStatistics = {}
    numTests = 0
    for (thisChr, thisReplicateList) in INPUT_FILES.items():
        print (f"+++++ {thisChr}")
        statsList = []
        for repFile in thisReplicateList:
            print (repFile, pathlib.Path(repFile).is_file())
            # unpickle
            pickleData = pickle.load (bz2.open (repFile, 'rb'))
            statsList.append (pickleData['statsFrame'])
            numTests += statsList[-1].shape[0]

        allStatistics[thisChr] = statsList

    print (f"-- total number of loci: {numTests}")

    return allStatistics


def computeFisher (allStatistics):

    # get the subset of positions that have a p-value in every replicate
    commonPositions = {}
    for (thisChr, thisStatsList) in allStatistics.items():
        print (f"+++++ {thisChr}")
        thisCommonSet = None
        for thisStatsFrame in thisStatsList:
            thisPositions = thisStatsFrame['positions']
            print (f"-- {len(thisPositions)}")
            if (thisCommonSet is None):
                thisCommonSet = set(thisPositions)
            else:
                thisCommonSet = thisCommonSet.intersection(set(thisPositions))
        print (len(thisCommonSet))
        # sort the common positions
        commonPositions[thisChr] = numpy.array(sorted(thisCommonSet))

    # compute combned p-values using Fishers method at the common positions
    fisherPValues = {}
    for (thisChr, thisStatsList) in allStatistics.items():
        print (f"+++++ {thisChr}")
        thisFisherStatistic = None
        for thisStatsFrame in thisStatsList:
            # get the right positions
            thisPositions = thisStatsFrame['positions']
            # make current positions are sorted they are ordered
            assert (numpy.diff(thisPositions).min() > 0)
            allMask = numpy.isin (thisPositions, commonPositions[thisChr])
            # these are the correct values, since everything is ordered
            thisP = thisStatsFrame['pValue'][allMask].to_numpy(dtype=float)
            thisLogP = -2*numpy.log(thisP)
            # collect values for fishers statistic
            if (thisFisherStatistic is None):
                thisFisherStatistic = thisLogP
            else:
                thisFisherStatistic += thisLogP
                
        # do the chisq for fisher one
        df = 2*len(thisStatsList)
        # print (df)
        thisPValues =  1 - scipy.stats.chi2.cdf (thisFisherStatistic, df)
        thisPValues = numpy.minimum (thisPValues, 1-1e-10)
        thisPValues = numpy.maximum (thisPValues, 1e-16)
        # fisherPValues[thisChr] = thisFisherStatistic
        fisherPValues[thisChr] = thisPValues

    return (commonPositions, fisherPValues)


def saveFisher (commonPositions, fisherPValues):

    pickleDict = {
        "commonPositions" : commonPositions,
        "fisherPValues" : fisherPValues,
    }
    pickle.dump (pickleDict, bz2.open (OUTPUT_FILE, "wb"))


def collectResults():

    # load the per replicate results
    print("[LOAD_PER_REPLICATE]")    
    allStatistics = loadPerReplicateResults()
    print("[LOAD_PER_REPLICATE_DONE]")    

    # compute fishers p-values
    print ("[COMPUTE_FISHER]")
    (commonPositions, fisherPValues) = computeFisher (allStatistics)
    print ("[COMPUTE_FISHER_DONE]")

    # and save the values
    print ("[SAVE_FISHER]")
    saveFisher (commonPositions, fisherPValues)
    print ("[SAVE_FISHER_DONE]")


def main():

    collectResults()


if __name__ == "__main__":
    main()

